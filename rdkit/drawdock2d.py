

from rdkit import Chem as C
from rdkit.Chem import Draw as D
from rdkit.Chem.Draw import SimilarityMaps as DSM


def GetAtomicWeightsFromProp(mol, propName, default=0, getTypedProp=C.Atom.GetDoubleProp):
    weights = []
    for atomId in range(mol.GetNumAtoms()):
        atom = mol.GetAtomWithIdx(atomId)
        try:
            weight = getTypedProp(atom, propName)
        except Exception:
            if default is not None:
                weight = default
            else:
                raise
        weights.append(weight)
    return weights


def GetSimilarityMapForProp(mol, propName, default=0, getTypedProp=C.Atom.GetDoubleProp, **kwargs):
    weights = GetAtomicWeightsFromProp(mol, propName, default=default, getTypedProp=getTypedProp)
    weights, maxWeight = DSM.GetStandardizedWeights(weights)
    fig = DSM.GetSimilarityMapFromWeights(mol, weights, **kwargs)
    return fig, maxWeight


