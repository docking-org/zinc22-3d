#!/usr/bin/env python
from __future__ import print_function

# FROM
# http://rdkit.org/UGM/2012/Ebejer_20110926_RDKit_1stUGM.pdf
# http://pubs.acs.org/doi/abs/10.1021/ci2004658

import argparse
import logging
import multiprocessing
import operator
import sys


from rdkit.Chem import (
    MolFromSmiles,
    MolFromMol2File,
    SanitizeMol,
    AllChem,
    MolToMolBlock,  # TODO Implement real mol2
    AddHs,
)
from rdkit.Chem.rdDistGeom import EmbedMultipleConfs
from rdkit.Chem.ChemicalForceFields import (
    MMFFGetMoleculeForceField,
    MMFFGetMoleculeProperties,
    UFFGetMoleculeForceField,
)
from rdkit.Chem.rdMolAlign import AlignMol
from rdkit.Chem.rdMolDescriptors import CalcNumRotatableBonds


RD_NAME = '_Name'
CONF_ENERGY = 'Conformer Energy'
SDF_MODEL_END = '$$$$'

FORCEFIELDS = {
    'mmff94': lambda mol, *args, **kwargs: MMFFGetMoleculeForceField(mol, MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94'), *args, **kwargs),
    'mmff94s': lambda mol, *args, **kwargs: MMFFGetMoleculeForceField(mol, MMFFGetMoleculeProperties(mol, mmffVariant='MMFF94s'), *args, **kwargs),
    'uff': UFFGetMoleculeForceField,
}
DEFAULT_FORCEFIELD = 'mmff94s'


# Heuristic from original paper
# NOTE: increasing these numbers by x2 to x4 will occasioanlly find lower
# engery conformations without adding too much compute time.
# Something investigate in the future
def get_num_confs_for_mol(m):
    rb = CalcNumRotatableBonds(m)
    if rb <= 7:
        return 50
    elif 8 <= rb <= 12:
        return 200
    else:
        return 300


def optimize_single_conformer(mol, conf_id, interfragment, max_steps, forcefield=DEFAULT_FORCEFIELD):
    ignore_frag = not interfragment
    get_forcefield = FORCEFIELDS[forcefield]
    forcefield = get_forcefield(mol, confId=conf_id,
                                     ignoreInterfragInteractions=ignore_frag)
    for _ in range(max_steps):  # Run up to ten steps
        if forcefield.Minimize() == 0:   # If convergence
            break
    energy = forcefield.CalcEnergy()
    return energy


def optimize_single_conformer_star(args):
    return optimize_single_conformer(*args)


def optimize_conformers(mol, interfragment=True, max_steps=10, parallelism=None, forcefield=DEFAULT_FORCEFIELD):
    conf_energy = {}
    conf_ids = [conf.GetId() for conf in mol.GetConformers()]

    if parallelism is None:
        for conf_id in conf_ids:
            energy = optimize_single_conformer(mol, conf_id, forcefield=forcefield,
                                                             interfragment=interfragment,
                                                             max_steps=max_steps)
            conf_energy[conf_id] = energy
    else:
        args = [(mol, conf_id, interfragment, max_steps, forcefield) for conf_id in conf_ids] 
        pool = multiprocessing.Pool(processes=parallelism)
        energy = pool.map(optimize_single_conformer_star, args)
        conf_energy = dict(zip(conf_ids, energy))

    return conf_energy


def generate_conformers(mol, add_hydrogens=True,
                             rmsd_threshold=2.0,    # Arbitrarily selected
                             num_conformers=None,   # None means best guess
                             parallelism=None,
                             forcefield=DEFAULT_FORCEFIELD,
                             log=logging):
    if add_hydrogens:
        log.info("Adding implicit hydrogens")
        mol = AddHs(mol)

    if num_conformers is None:
        num_conformers = get_num_confs_for_mol(mol)

    log.info("Attempting to generate {0} conformations with min RMSD of {1:.4f}".format(num_conformers, rmsd_threshold))
    orig_conf_ids = EmbedMultipleConfs(mol, numConfs=num_conformers, 
                                            pruneRmsThresh=rmsd_threshold,
                                            ignoreSmoothingFailures=True)  # Prevents crashes in some situations
    log.info("Generated {0} initial conformations".format(len(orig_conf_ids)))

    log.info("Optimizing and calculating energies using {0}".format(forcefield))
    conf_energy = optimize_conformers(mol, interfragment=True, 
                                           parallelism=parallelism,
                                           forcefield=forcefield)
    sorted_by_energy = sorted(conf_energy.iteritems(), key=operator.itemgetter(1))

    log.info("Filtering similar conformers")
    selected = []
    min_rmsd, max_rmsd = float('inf'), float('-inf')
    for idx, id_energy in enumerate(sorted_by_energy):
        conf_id, energy = id_energy
        keep = True
        for comp_id, other_energy in sorted_by_energy[idx+1:]:
            rmsd = AlignMol(mol, mol, prbCid=comp_id, refCid=conf_id)
            if rmsd <= rmsd_threshold:
                mol.RemoveConformer(conf_id)
                keep = False
                break
            else:
                if rmsd < min_rmsd:
                    min_rmsd = rmsd
                if rmsd > max_rmsd:
                    max_rmsd = rmsd
        if keep:
            selected.append(id_energy)
    AllChem.AlignMolConformers(mol, atomIds=[3, 4, 5, 6, 7, 8, 9, 21, 24, 25, 26, 32], maxIters=200)
    log.debug("Removed {0} after post-optimization RMSD filtering".format(len(orig_conf_ids) - len(selected)))
    log.info("RMSD: min={0:.4f} max={1:.4f}".format(min_rmsd, max_rmsd))
    return mol, selected


def dump_conformers_sdf(mol, output, conf_ids=None, 
                                     energies=None, 
                                     renumber=True):
    if conf_ids is None:
        conformers = mol.GetConformers()
    else:
        conformers = (mol.GetConformer(conf_id) for conf_id in conf_ids)

    # Record state of properties that may be overwritten
    original_name = mol.GetProp(RD_NAME)
    if energies is not None and mol.HasProp(CONF_ENERGY):
        original_energy = mol.GetProp(CONF_ENERGY)
    else:
        original_energy = None

    conformer_names = []

    # Render conformers
    for idx, conf in enumerate(conformers):
        conf_id = conf.GetId()
        if renumber:
            conf_idx = idx
        else:
            conf_idx = conf_id

        if energies is not None and conf_id in energies:
            energy = energies[conf_id]
            mol.SetProp(CONF_ENERGY, "{0:0.4f}".format(energy))
   
        conf_name = "{0}_{1}".format(original_name, conf_idx)
        mol.SetProp(RD_NAME, conf_name)
        conformer_names.append(conf_name)
        block = MolToMolBlock(mol, includeStereo=True, confId=conf_id)
        print(block, file=output, end="")
        print(SDF_MODEL_END, file=output)
    # Reset changes to mol properties
    mol.SetProp(RD_NAME, original_name)
    if original_energy is not None:
        mol.SetProp(CONF_ENERGY, original_energy)
    else:
        mol.ClearProp(CONF_ENERGY)

    return conformer_names


def main(args, output=sys.stdout, log=logging):
    parser = argparse.ArgumentParser(
    """RDKit-based conformer generation proof-of-concept.
    This program accepts either a mol2 file or a SMILES string and produces an SD file
    """)
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-m', '--mol2', type=str, help="Mol2 file to gererate conformers for")
    input_group.add_argument('-s', '--smiles', type=str, help="SMILES string of molecule")

    parser.add_argument('-N', '--name', type=str, default=None, help="Molecule name")
    parser.add_argument('-H', '--no-hydrogens', action='store_true', 
                                                default=False, 
                                                help="Do NOT explicitly add implicit Hydrogens to conformers [default: %(default)s]")
    parser.add_argument('-r', '--rmsd-threshold', type=float,
                                                  default=2.0,
                                                  help="Only accept conformers that have an RMSD of at least this value from previously seen conformers [default: %(default)s")
    parser.add_argument('-n', '--num-conformers', type=int,
                                                  default=None,
                                                  help="Number of conformers to initially generate [default: auto]")
    parser.add_argument('-F', '--forcefield', type=str,
                                              default=DEFAULT_FORCEFIELD,
                                              choices=FORCEFIELDS.keys(),
                                              help="Forcefield to use for optimization [default: %(default)s]")
    parser.add_argument('-P', '--parallelism', type=int,
                                               default=None,
                                               help="Number of processes to use [default: 1]")
    params = parser.parse_args(args)

    # Load input molecule
    if hasattr(params, 'mol2') and params.mol2 is not None:
        mol = MolFromMol2File(params.mol2, sanitize=False)
    else:
        mol = MolFromSmiles(params.smiles, sanitize=False)

    try:
        SanitizeMol(mol)
    except ValueError as e:
        log.critical("Could not sanitize molecule: {0}:".format(str(e)))
        sys.exit(2)
    except Exception:  # This is `Boost.Python.ArgumentError`
        log.critical("Could not parse molecule!")
        sys.exit(2)
        

    # Assign user-provided name if applicable
    if params.name is not None:
        mol.SetProp(RD_NAME, params.name)
    elif not mol.HasProp(RD_NAME):
        mol.SetProp(RD_NAME, 'Ligand')

    # Generate 3D conformers
    embedded, selected = generate_conformers(mol, 
                                             add_hydrogens=not params.no_hydrogens,
                                             rmsd_threshold=params.rmsd_threshold,
                                             num_conformers=params.num_conformers,
                                             parallelism=params.parallelism,
                                             forcefield=params.forcefield,
                                             log=log)

    log.info("Conformers selected: {0}".format(len(selected)))
    log.info("Energy: min={0:.4f} kcal/mol max={1:.4f} kcal/mol".format(selected[0][1], selected[-1][1]))

    # Find lowest-energy conformers
    sorted_by_energy = [item[0] for item in selected]

    # Render SDF file
    names = dump_conformers_sdf(embedded, output, conf_ids=sorted_by_energy, 
                                                  renumber=True)

    for name, (conf_id, energy) in zip(names, selected):
        log.info("\t{0}: {1:0.4f} kcal/mol".format(name, energy))

    return 0


if __name__ == '__main__':
    log = logging.getLogger()
    log.addHandler(logging.StreamHandler(sys.stderr))
    log.setLevel(logging.DEBUG)
    main(sys.argv[1:], sys.stdout, log=log)




