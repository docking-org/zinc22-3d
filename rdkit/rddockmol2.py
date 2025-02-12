
import ast
from cStringIO import StringIO
from collections import defaultdict
import itertools
from rdkit import Chem as C


_TYPED_MOL_PROPERTIES_POSSIBLE = all(
    hasattr(C.Mol, p) for p in (
        'SetIntProp', 'SetDoubleProp', 'SetBoolProp'
    )
)
_TYPED_ATOM_PROPERTIES_POSSIBLE = all(
    hasattr(C.Atom, p) for p in (
        'SetIntProp', 'SetDoubleProp', 'SetBoolProp'
    )
)


def monomer_field_to_setter(monomer, field):
    if callable(field):
        return lambda value: field(monomer, value)
    elif hasattr(monomer, field):
            return monomer.field
    else:
        parts = [part.capitalize() for part in field.split('_')]
        setter_field = ''.join(['Set'] + parts)
        if hasattr(monomer, setter_field):
            return getattr(monomer, setter_field)
    return None
    


def atom_props_to_monomer_info(mol, mappings, monomer_cls=C.AtomPDBMonomerInfo, prop_getter=C.Atom.GetProp):
    props_set = 0
    for atom in mol.GetAtoms():
        for prop, field in mappings.items():
            if atom.HasProp(prop):
                monomer, is_new = atom.GetMonomerInfo(), False
                if monomer is None:
                    monomer = monomer_cls()
                    is_new = True
                value = prop_getter(atom, prop)
                setter = monomer_field_to_setter(monomer, field)
                setter(value)
                if is_new:
                    atom.SetMonomerInfo(monomer)
                props_set += 1
    return props_set


def AssignPropsToPDBMonomerInfo(mol, 
                                tempProp=None, 
                                occupancyProp=None, 
                                getTypedProp=C.Atom.GetDoubleProp):
    mappings = {}
    if tempProp is not None:
        mappings[tempProp] = lambda monomer, value: monomer.SetTempFactor(float(value))
    if occupancyProp is not None:
        mappigns[occupancyProp] = lambda monomer, value: monomer.SetOccupancy(float(value))

    assigned = atom_props_to_monomer_info(mol, mappings, 
                                          monomer_cls=C.AtomPDBMonomerInfo, 
                                          prop_getter=getTypedProperty)
    return assigned
                

class ForwardDockMol2Loader(object):
    COMMENTS_PREFIX = '#'
    MOL_PROPERTY_PREFIX = '########## '
    MOL_PROPERTY_DELIMITER = ':'
    ATOM_PROPERTY_PREFIX = '# '
    ATOM_NUMBER_FIELD = 'num'
    MOL2_DELIMITER = '@<TRIPOS>MOLECULE'

    _EMPTY = object()

    def __init__(self, source, **kwargs):
        self._source = source
        self._done = False
        self._num_seen = 0
        self._peek = []
        self._comment_buffer = []
        self._mol2_block_buffer = []
        self._eval_mol_props = _TYPED_MOL_PROPERTIES_POSSIBLE
        self._eval_atom_props = _TYPED_ATOM_PROPERTIES_POSSIBLE
        self._atom_numbering_offset = 1
        self._atom_props_header = None
        self._strict_comments = True
        
        # removeHs after parsing to make property assignment work
        self._removeHs = kwargs.setdefault('removeHs', False)
        kwargs['removeHs'] = False
        self._kwargs = kwargs

    def __iter__(self):
        return self

    def __next__(self):
        return self._next_mol()

    def next(self):
        return self._next_mol()

    def tell(self):
        return self._num_seen

    def _get_peeked_source(self):
        source = self._source
        if self._peek:
            source = itertools.chain(self._peek, source)
            self._peek = []
        return source
    
    def _fill_buffers(self):
        in_mol2_block = False
        for line in self._get_peeked_source():
            if in_mol2_block:
                if line.startswith(self.MOL2_DELIMITER) \
                   or line.startswith(self.COMMENTS_PREFIX):
                    self._peek.append(line)
                    break
            elif line.startswith(self.MOL2_DELIMITER):
                in_mol2_block = True
                
            if in_mol2_block:
                self._mol2_block_buffer.append(line)
            else:
                self._comment_buffer.append(line)
        else:
            self._done = True

    def _reset_buffers(self):
        comments = self._comment_buffer
        mol2_block = self._mol2_block_buffer
        self._comment_buffer = []
        self._mol2_block_buffer = []
        return comments, mol2_block

    def _consume_buffers(self):
        comments, mol2_block = self._reset_buffers()
        mol2_block = ''.join(mol2_block)
        mol = C.MolFromMol2Block(mol2_block, **self._kwargs)
        mol = self._init_mol(mol)
        if mol is not None:
            self._num_seen += 1
            self._set_mol_props_from_comments(mol, comments)
        if self._removeHs:
            mol = C.RemoveHs(mol)
        return mol

    def _init_mol(self, mol):
        return mol

    def _next_mol(self, raise_stop=True):
        if not self._done:
            self._fill_buffers()
            result = self._consume_buffers()
            return result
        elif raise_stop:
            raise StopIteration
        else:
            return None

    def _set_mol_props_from_comments(self, mol, comments):
        atom_props_header = self._atom_props_header
        atom_index_offset = self._atom_numbering_offset
        for line in comments:
            try:
                if line.startswith(self.MOL_PROPERTY_PREFIX):
                    prop_data = line[len(self.MOL_PROPERTY_PREFIX):]
                    prop_name, prop_value = prop_data.split(self.MOL_PROPERTY_DELIMITER, 1)
                    prop_name = prop_name.strip()
                    prop_value = prop_value.strip()
                    if self._eval_mol_props:
                        prop_value = self._try_to_eval(prop_value)
                    self._set_typed_prop(mol, prop_name, prop_value)
                elif line.startswith(self.ATOM_PROPERTY_PREFIX):
                    prop_data = line[len(self.ATOM_PROPERTY_PREFIX):]
                    tokens = prop_data.split()
                    tokens = [t.strip() for t in tokens]
                    if atom_props_header is None:
                        atom_props_header = tokens
                    else:
                        if self._eval_atom_props:
                            values = [self._try_to_eval(token) for token in tokens]
                        else:
                            values = tokens
                        header = atom_props_header or map(str, range(len(values)))
                        props = dict(zip(header, values))
                        atom_idx = int(props.pop(self.ATOM_NUMBER_FIELD, header[0]))
                        atom = mol.GetAtomWithIdx(atom_idx - atom_index_offset)
                        for prop_name, prop_value in props.items():
                            self._set_typed_prop(atom, prop_name, prop_value)
            except Exception:
                if self._strict_comments:
                    raise

    def _set_typed_prop(self, obj, name, value, typed=True):
        if not typed or isinstance(value, basestring):
            obj.SetProp(name, str(value))
        elif isinstance(value, (int, long)):
            obj.SetIntProp(name, value)
        elif isinstance(value, (float,)):
            obj.SetDoubleProp(name, value)
        elif isinstance(value, (bool,)):
            obj.SetBoolProp(name, value)
        else:
            obj.SetProp(name, str(value))

                    
    def _try_to_eval(self, value, should_eval=True):
        if not should_eval:
            return value
        try:
            return ast.literal_eval(value)
        except (SyntaxError, ValueError):
            return value
        
