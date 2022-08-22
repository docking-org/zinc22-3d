#!/usr/bin/env python
from __future__ import division
import os
import sys
import logging
from openeye.oechem import *
import copy
#from openeye.oeomega import *

def write_molecule(outfile, mol):
    ofs = oemolostream(outfile)
    OEWriteMolecule(ofs, mol)
    ofs.close()

#def set_defaults(omega, limitConfs=200, energyWindow=12):
#    # File Parameters
#    omega.SetCommentEnergy(True)
#    #omega.SetIncludeInput(False)
#    omega.SetIncludeInput(True)
#    omega.SetRotorOffset(False) #Enkhee
#    omega.SetSDEnergy(False)
#    omega.SetWarts(True)
#    # 3D Construction Parameters
#    omega.SetBuildForceField('mmff94s')
#    omega.SetCanonOrder(True)
#    omega.SetFixDeleteH(True)
#    omega.SetDielectric(1.0)
#    omega.SetExponent(1.0)
#    omega.SetFixRMS(0.15)
#    #omega.SetFixRMS(0.01)# jklyu modify this to keep more conformors
#    omega.SetFromCT(False)
#    omega.SetFixMaxMatch(1)
#    omega.SetFixUniqueMatch(True)
#    # Structure Enumeration Parameters
#    omega.SetEnumNitrogen(False)
#    omega.SetEnumRing(False)#Enkhee  # TS & MK 20160524 (from False) (Improves ring pucker sampling)
#    # Torsion Driving Parameters
#    omega.SetEnergyWindow(energyWindow)  # JJI 20160329  # TS & MK 20160524 (from 6)
#    omega.SetMaxConfs(limitConfs)
#    omega.SetMaxRotors(-1)
#    omega.SetMaxSearchTime(120.0)
#    omega.SetRangeIncrement(5)
#    omega.SetRMSThreshold(0.50)  # JJI 20160329 # TS & MK 20160524 (from .5)
#    #omega.SetRMSThreshold(0.01)  # jklyu modify this to keep more conformors
#    omega.SetSearchForceField('mmff94s')
#    
#    omega.SetTorsionDrive(True)
#    #omega.SetTorsionDrive(False)
#    # Stereochemsitry
#    #omega.SetStrictStereo(False)
#    omega.SetStrictAtomTypes(False)

# Parse arguments here
infile = sys.argv[1]

if os.environ.get('OMEGA_ENERGY_WINDOW', '').strip() != '':
    energyWindow = int(os.environ['OMEGA_ENERGY_WINDOW'])
else:
    energyWindow = 12

if os.environ.get('OMEGA_MAX_CONFS', '').strip() != '':
    limitConfs = int(os.environ['OMEGA_MAX_CONFS'])
else:
    limitConfs = 200

if len(sys.argv) > 2:
    if sys.argv[2].isdigit():
        numHs = int(sys.argv[2])
        if numHs > 5:
            print("Refusing to build conformations with > 5 rotatable hydrogens")
            sys.exit(-1)
        if numHs >= 4:
            logging.warning("4-5 Rotatable hydrogens  reported. Reducing confs by a factor of 30")
            limitConfs = limitConfs // 30
        elif numHs >= 2:
            logging.warning("2-3 Rotatable hydrogens  reported. Reducing confs by a factor of 3")
            limitConfs = limitConfs // 3


#logging.warn('Setting energy window to %d and max confs to %d' % (energyWindow, limitConfs))
logging.warning('Setting energy window to %d and max confs to %d' % (energyWindow, limitConfs))

#omega = OEOmega()
#set_defaults(omega, limitConfs=limitConfs, energyWindow=energyWindow)

mol = OEMol()
inroot, inext = os.path.splitext(infile)
ifs = oemolistream(infile)
OEReadMolecule(ifs, mol)
ifs.close()
write_molecule(inroot+'.proc.mol2', mol)
OEDetermineComponents(mol)
write_molecule(inroot+'.proc2.mol2', mol)
count, ringlist = OEDetermineRingSystems(mol)
rcf = open('ring_count', 'w')
rcf.write('%d\n' % count)
rcf.close()
outfiles = []
fail_count = 0

#print ("bonds")
#for bond in mol.GetBonds():
#    print (bond)
#    #print(bond.GetOrder())


#ofs = oemolostream('mol2')

for r in ringlist: 
    print(r)

for ringidx in range(1, count + 1):
    print(ringidx, ". ring system:", end=" ")
    for atom in mol.GetAtoms():
        if ringlist[atom.GetIdx()] == ringidx:
            print(atom.GetIdx(), end=" ")
    print()     

ring_p_bond_list = copy.copy(ringlist)

for ringidx in range(1, count + 1):
    #print("Neighbors of ring %d:"%ringidx, end=" ")
    strbonds = "Neighbors of ring %d:"%ringidx
    for atom in mol.GetAtoms():
       #print("Atom:", end=" ")
       #print(atom.GetIdx(), end=" ")
    
        if ringlist[atom.GetIdx()] == ringidx:
           for bond in atom.GetBonds():
              nbor = bond.GetNbr(atom)
              if ringlist[nbor.GetIdx()] == ringidx: # if bonded atom is part of the ring skip
                  continue
              print(bond)
              if ring_p_bond_list[nbor.GetIdx()] == 0: 
                  ring_p_bond_list[nbor.GetIdx()] = ringidx
              #if ring_p_bond_list[nbor.GetIdx()] != ringidx: 
              #    ring_p_bond_list[nbor.GetIdx()] = ringidx
              #print(nbor.GetIdx(), end=" ")
              strbonds = strbonds + " " + str(nbor.GetIdx())
    print(strbonds)     

#exit()

submol = OEGraphMol()
adjustHcount = True
#adjustHcount = False
#pred = OEPartPredAtom(ringlist)
pred = OEPartPredAtom(ring_p_bond_list)

for i in range(1, count+1):
     pred.SelectPart(i)
     print(pred)
     #print(d)
     outfile = "%s.%d.ring.mol2" % (inroot, i)
     #oechem.OESubsetMol(submol, mol, includeexo, adjustHcount)
     OESubsetMol(submol, mol, pred, adjustHcount)
     #molcopy = OEMol(mol, pred)
     #write_molecule(outfile, molcopy)
     write_molecule(outfile, submol)



#parms = "-ewindow %f  -maxconfs %f" % (energyWindow, limitConfs )
#parm = " -commentEnergy  true -in  failed/bi2852/0.failed/0.mol2 -includeInput  true -out  0.confs.mol2 -rotorOffsetCompress  false -sdEnergy  false -warts  true"
#parm = parm + " -buildff  mmff94s -canonOrder  true -deleteFixHydrogens  true -dielectric  1.0 -exponent  1 -fixfile  ring.1.mol2 -fixrms  0.15 -fromCT  false -maxmatch  1 -umatch  true"
#parm = parm + " -ewindow  10.0 -maxconfs  200 -maxtime  120.0 -rangeIncrement  5 -rms  0.5 -searchff  mmff94s"

parm = " -commentEnergy  true -includeInput  true -rotorOffsetCompress  false -sdEnergy  false -warts  true"
parm = parm + " -buildff  mmff94s -canonOrder  true -deleteFixHydrogens  true -dielectric  1.0 -exponent  1 -fixrms  0.15 -fromCT  false -maxmatch  1 -umatch  true"
parm = parm + " -enumNitrogen  false -enumRing  false"
parm = parm + " -ewindow %f -maxconfs  %d -maxtime  120.0 -rangeIncrement  5 -rms  0.5 -searchff  mmff94s" % (energyWindow, limitConfs)

if count == 0:
   print ("${OE_DIR}/bin/omega2 -in %s -out %s.%d.db2in.mol2 %s" % (infile, inroot, 0, parm))
   os.system("${OE_DIR}/bin/omega2 -in %s -out %s.%d.db2in.mol2 %s" % (infile, inroot, 0, parm))

for i in range(1, count + 1):
   print ("${OE_DIR}/bin/omega2 -fixfile %s.%d.ring.mol2 -in %s -out %s.%d.db2in.mol2 %s" % (inroot, i, infile, inroot, i, parm))
   os.system("${OE_DIR}/bin/omega2 -fixfile %s.%d.ring.mol2 -in %s -out %s.%d.db2in.mol2 %s" % (inroot, i, infile, inroot, i, parm))



#print ('energy: ', omega.GetEnergyWindow())
#print ('conf: ', omega.GetMaxConfs() )
#if count == 0:
#    omega.SetMaxConfs(30)
#    outfile = "%s.%d.db2in.mol2" % (inroot, count) 
#    outfiles.append(outfile)
#    if omega(mol):
#        write_molecule(outfile, mol)
#    else:
#        fail_count +=1 
#else:
#    pred = OEPartPredAtom(ringlist)
#    for i in range(1, count+1):
#        pred.SelectPart(i)
#        outfile = "%s.%d.db2in.mol2" % (inroot, i) 
#        outfiles.append(outfile)
#        molcopy = OEMol(mol)
#        if omega(molcopy, pred):
#            write_molecule(outfile, molcopy)
#        else:
#            fail_count +=1
#            
#print (outfiles)
#
#sys.exit(fail_count)
