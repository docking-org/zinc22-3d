#!/usr/bin/env python



import argparse
import itertools
import sys
from rdkit import Chem as C

LIMIT = 4
TAGS = [
    C.ChiralType.CHI_TETRAHEDRAL_CCW, 
    C.ChiralType.CHI_TETRAHEDRAL_CW,
]



def is_new_chiral(atom):
    symbol = atom.GetSymbol()
    neighbors = atom.GetNeighbors()
    neighbor_symbols = [n.GetSymbol() for n in neighbors]
    if symbol == 'N':
        return len(neighbors) == 4 and neighbor_symbols.count('H') == 1
    elif symbol == 'P':
        return True
    else:
        return False


def create_stereoisomer(mol, *stereocenters):
    expanded = C.Mol(mol)
    for idx, tag in stereocenters:
        expanded.GetAtomWithIdx(idx).SetChiralTag(tag)
    return expanded
        

def expand_stereo_centers(mol, reference=None, limit=LIMIT, strict=True):
    molH = C.AddHs(mol, explicitOnly=True)
    if reference is not None:
        refH = C.AddHs(reference, explicitOnly=True)
        ref_match = [(ridx, midx) for midx, ridx in enumerate(molH.GetSubstructMatch(refH))]
        if not ref_match:
            ref_match = [(ridx, midx) for ridx, midx in enumerate(refH.GetSubstructMatch(molH))]
        ref_match = dict(ref_match)
    else:
        ref_match = {}

    if ref_match:
        ref_centers = C.FindMolChiralCenters(refH)
        ref_stereo = dict(ref_centers)
        ref_stereo = [(ref_match.get(idx), refH.GetAtomWithIdx(idx).GetChiralTag()) for idx in ref_stereo]
        ref_stereo = dict(ref_stereo)
    else:
        ref_stereo = {}

    centers = C.FindMolChiralCenters(molH, includeUnassigned=True)
    unassigned = [idx for idx, tag in centers if tag == '?']

    if ref_stereo:
        still_unassigned = []
        previously_assigned = []
        for idx in unassigned:
            if idx in ref_stereo:
                stereoisomer = (idx, ref_stereo[idx])
                previously_assigned.append(stereoisomer)
            else:
                still_unassigned.append(idx)
        unassigned = still_unassigned
    else:
        previously_assigned = []

    for stereocenter in previously_assigned:
        mol = create_stereoisomer(mol, stereocenter)

    if strict:
        unassigned = [idx for idx in unassigned 
                      if is_new_chiral(molH.GetAtomWithIdx(idx))]

    expansions = itertools.product(unassigned, TAGS)
    if limit is not None:
        expansions = itertools.islice(expansions, limit)

    stereocenter = None
    for stereocenter in expansions:
        expanded = create_stereoisomer(mol, stereocenter)
        yield expanded

    # If non were added, yield original molecule
    if stereocenter is None:
        yield mol


def expand_centers(mols, reference=None, limit=LIMIT, strict=True):
    if limit < 1:
        limit = None

    for mol in mols:
        expansions = expand_stereo_centers(mol, reference=None, limit=limit, strict=strict)
        for expanded in expansions:
            yield expanded
    

def main(args):
    parser = argparse.ArgumentParser()
    parser.add_argument('input', type=str, help="Where to read protomers from")
    parser.add_argument('output', type=str, nargs='?', default='-', help="Where to write protomers to")
    parser.add_argument('-r', '--reference', nargs='?', default=None, type=str, help="Reference SMILES string to re-assign any lost stereochemistry NOT IMPLEMENTED")
    parser.add_argument('-l', '--limit', type=int, nargs='?', default=LIMIT, help="Max expansions to accept")
    params = parser.parse_args(args)

    with open(params.input) as f:
        for char in next(f, ''):
            if char.isspace():
                sep = char
                break
        else:
            sep = ' '

    if params.reference is not None and False:
        reference = C.SmilesMolSupplier(params.reference, delimiter=sep, titleLine=False, sanitize=False)
    else:
        reference = None
    reader = C.SmilesMolSupplier(params.input, delimiter=sep, titleLine=False, sanitize=False)
    try:
        props = list(next(reader).GetPropNames())
    except AttributeError:
        props = []
    writer = C.SmilesWriter(params.output, delimiter=sep, includeHeader=False, isomericSmiles=True)
    writer.SetProps(props)
    expanded = expand_centers(reader,
                              reference=reference,
                              limit=params.limit, 
                              strict=True)

    for mol in expanded:
        writer.write(mol)
            
    return 0

if __name__ == '__main__':
    sys.exit(main(sys.argv[1:]))


