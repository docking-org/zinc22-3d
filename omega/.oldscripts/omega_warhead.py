#!/usr/bin/env python
import os
import sys
from openeye.oechem import *
from openeye.oeomega import *

def read_parm(parm):
    if os.path.exists(parm):
        f = open(parm, 'r')
        for line in f:
            line = line.strip()
            if line and line[0] != '#':
                eval('omega.%s' % line)
        f.close()

def write_molecule(outfile, mol):
    ofs = oemolostream(outfile)
    OEWriteMolecule(ofs, mol)
    ofs.close()

def set_defaults(omega):
    # File Parameters
    omega.SetCommentEnergy(True)
    omega.SetIncludeInput(False)
    omega.SetRotorOffset(True)
    omega.SetSDEnergy(False)
    omega.SetWarts(True)
    # 3D Construction Parameters
    omega.SetBuildForceField('mmff94s_noEstat')
    omega.SetCanonOrder(True)
    omega.SetFixDeleteH(True)
    omega.SetDielectric(1.0)
    omega.SetExponent(1.0)
    omega.SetFixRMS(0.15)
    omega.SetFromCT(True)
    omega.SetFixMaxMatch(1)
    omega.SetFixUniqueMatch(True)
    # Structure Enumeration Parameters
    omega.SetEnumNitrogen(True)
    omega.SetEnumRing(True)
    # Torsion Driving Parameters
    omega.SetEnergyWindow(30.0)
    #omega.SetMaxConfGen(100000)
    #omega.SetMaxConfs(10000)
    #jklyu,2017.07.25
    omega.SetMaxConfs(50)
    #omega.SetMaxPoolSize(20000)
    omega.SetMaxRotors(-1)
    omega.SetMaxSearchTime(120.0)
    omega.SetRangeIncrement(5)
    omega.SetRMSThreshold(0.5) # 0.8 for large molecules.
    omega.SetSearchForceField('mmff94s_noEstat')
    omega.SetTorsionDrive(True)

# Parse arguments here
infile = sys.argv[1]

omega = OEOmega()
set_defaults(omega)

mol = OEMol()
inroot, inext = os.path.splitext(infile)
ifs = oemolistream(infile)
OEReadMolecule(ifs, mol)
ifs.close()

#initialize a substructure search with a smarts pattern
ss = OESubSearch("[Si]([H])([H])[H]")

#create an array of atom molecules to define the rigid part
phos = [0]*(OEMolBase.GetMaxAtomIdx(mol))

print ('warning *!*!* you are now running COVALENT omega script - expecting SiH3' )

#for every match to the smarts pattern (should be unique)
#flip the rigid atoms bit in phos
for count,match in enumerate(ss.Match(mol,"true")):
    for ma in match.GetAtoms():
        phos[ma.target.GetIdx()]=count+1

pred = OEPartPredAtom(phos)
for i in xrange(count+1):
    print ('i%d:',i+1)
    pred.SelectPart(i+1)
    outfile = "%s.%d.db2in.mol2" % (inroot, i+1) 
    molcopy = OEMol(mol)
    if omega(molcopy, pred):
        write_molecule(outfile, molcopy)

