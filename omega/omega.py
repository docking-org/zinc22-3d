#!/usr/bin/env python3.7
# written by Jiankun Lyu, 03/12/2020
from __future__ import division
import os
import sys
import logging
from openeye import oechem
from openeye import oeomega
from openeye import oeff

DOCKBASE='/'.join(__file__.split("/")[0:-3])
sys.path.append(DOCKBASE + "/ligand/mol2db2_py3_strain")
sys.path.append(DOCKBASE + "/ligand/generate")

from hydrogens import count_hydrogens
from mol2 import Mol2

# write out conformers to a mol2 file
def write_molecule(outfile, mol):
    ofs = oechem.oemolostream(outfile)
    oechem.OEWriteMolecule(ofs, mol)
    ofs.close()

# set default parameters for omega TorDriveOpts
# https://docs.eyesopen.com/toolkits/python/omegatk/OEConfGenClasses/OETorDriveOptions.html#OEConfGen::OETorDriveOptions
# set the omega torsion library: 1) Original; 2) GubaV21
def set_TorDriveOpts_defaults(TorDriveOpts, torLibType="GubaV21"):
    # Sets whether to add conformer energy in kcal/mol into the comments that appear if the output is written to a MOL2 file. If true, the comment energy is written. Default: False.
    TorDriveOpts.SetCommentEnergy(True)
    # Sets the limit of molecule with maximum number of rotors that would be handled by torsion driving. If a value of 9999 or higher is defined, this limit is ignored. Default: 9999.
    TorDriveOpts.SetMaxRotors(9999)
    # Sets the maximum time that should be spent for torsion driving a molecule, in units of seconds. Default: 120.0 sec.
    TorDriveOpts.SetMaxSearchTime(120.0)
    # Sets whether the rotor offset compression should be turned on. Turning on the rotor offset compression reduces the space required when saved in an OEB files. Default: False.
    TorDriveOpts.SetRotorOffset(False)
    # Sets the predicate that is used to decide whether a bond is rotatable.
    # it's an option, but we are not using this for db2 generation
    # TorDriveOpts.SetRotorPredicate(const OESystem::OEUnaryPredicate<OEChem::OEBondBase>&)
    # Sets whether the conformer energy in kcal/mol should be stored as SD data. Default: False.
    TorDriveOpts.SetSDEnergy(False)
    # Sets the torsion library that is used for torsion driving. Defaults: OETorLib with OETorLibType_GubaV21 mode.
    if torLibType=="GubaV21":
        TorDriveOpts.SetTorLib(oeomega.OETorLib(oeomega.OETorLibType_GubaV21))
    elif torLibType=="Original":
        TorDriveOpts.SetTorLib(oeomega.OETorLib(oeomega.OETorLibType_Original))
    else:
        print("Incorrect torType: %s" % torLibType)
        exit()

# set default parameters for omega SliceEnsembleOpts
# https://docs.eyesopen.com/toolkits/python/omegatk/OEConfGenClasses/OESliceEnsembleOptions.html#OEConfGen::OESliceEnsembleOptions
def set_SliceEnsembleOpts_defaults(SliceEnsembleOpts, limitConfs=200, energyWindow=12, rmsd=0.5, energyRangeFlag=False, limitConfsRangeFlag=False, rmsdRangeFlag=False):
    if energyRangeFlag == False:
        # Sets the maximum allowable energy difference between the lowest and the highest energy conformers, in units of kcal/mol. Default: 10.0 kcal/mol.
        SliceEnsembleOpts.SetEnergyWindow(energyWindow)
    else:
        # Sets the energy range that allows to define a varying EnergyWindow energy window in conjunction with the RangeIncrement. Defining the range overrides the singly defined energy widnow value.
        SliceEnsembleOpts.SetEnergyRange("5.0, 10.0, 15.0, 20.0")
    if limitConfsRangeFlag == False:
        # Sets the maximum number of conformers to be kept. Default: 200.
        SliceEnsembleOpts.SetMaxConfs(limitConfs)
    else:
        # Sets the maximum number of conformers range that allows to define a varying number of MaxConfs in conjunction with the RangeIncrement. Defining the range overrides the singly defined max confs value.
        SliceEnsembleOpts.SetMaxConfRange("600, 600, 600, 600")
    # Sets the maximum number of terminal heavy atoms that should be allowed in calculating RMSD between conformers. If the number of terminal heavy atoms in the molecule accedes the specified limit, all terminal heavy atoms are ignored in RMSD calculation. If the number is set to 0, this values is treated as not set. Default: 0.
    SliceEnsembleOpts.SetMaxTerminalHeavy(0)
    # Sets the number of rotatable bonds range to be used when the OEOmegaOptions.SetEnergyRange, OEOmegaOptions.SetRMSRange, and OEOmegaOptions.SetMaxConfRange are defines in terms of ranges. Default: 5.
    #SliceEnsembleOpts.SetRangeIncrement(5)
    SliceEnsembleOpts.SetRangeIncrement(3)
    if rmsdRangeFlag == False:
        # Sets the RMS threshold (Root Mean Square Cartesian distance) below which two conformers are treated as duplicates. Default: 0.5.
        SliceEnsembleOpts.SetRMSThreshold(rmsd)
    else:
        #Sets the RMS range that allows to define a varying RMSThreshold in conjunction with the RangeIncrement. Defining the range overrides the singly defined RMS threshold value.
        #SliceEnsembleOpts.SetRMSRange("0.2, 0.3, 0.4, 0.5")
        SliceEnsembleOpts.SetRMSRange("0.1, 0.2, 0.3, 0.4")
    # Sets flag if energy of the conformers should be used for slicing of conformers. This flag should be turned on if energies are optimized energies of the conformers. Default: False
    SliceEnsembleOpts.SetSliceByEnergy(False)
    # Sets the Energy threshold (Difference between energies) below which two conformers are treated as duplicates. This value is only meaningful when SetSliceByEnergy is set to true. Default: 0.01.
    SliceEnsembleOpts.SetEnergyThreshold(0.01)

#https://docs.eyesopen.com/toolkits/cpp/oefftk/OEFFClasses/OEMMFFSheffieldOptions.html#OEFF::OEMMFFSheffieldOptions
def set_ffOpts_defaults(ffOpts,ffType="MMFF94Smod"):
    # The Dielectric defines the solvent dielectric constant to be used for coulombic interactions with OEMMFF. Default: 1.0.
    ffOpts.SetDielectric(1.0)
    # The Exponent defines the coulombic exponent term to be used for coulombic interactions with OEMMFF. Default: 1.0.
    ffOpts.SetExponent(1.0)
    # The ForceFieldType defines the variation of forcefield to be used. Options are defined in the OEMMFFSheffieldFFType namespace. Default: OEMMFFSheffieldFFType::MMFF94Smod_NOESTAT.
    if ffType == "MMFF94Smod":
        ffOpts.SetForceFieldType(oeff.OEMMFFSheffieldFFType_MMFF94Smod)
    elif ffType == "MMFF94Smod_NOESTAT":
        ffOpts.SetForceFieldType(oeff.OEMMFFSheffieldFFType_MMFF94Smod_NOESTAT)
    elif ffType == "MMFF94Smod_TRUNC":
        ffOpts.SetForceFieldType(oeff.OEMMFFSheffieldFFType_MMFF94Smod_TRUNC)
    elif ffType == "MMFF94Smod_SHEFF":
        ffOpts.SetForceFieldType(oeff.OEMMFFSheffieldFFType_MMFF94Smod_SHEFF)
    else:
        print("Doesn't support %s ff" % ffType)
        exit()
    #ffOpts.SetForceFieldType(oeff.OEMMFFSheffieldFFType_MMFF94S)

def add_torison_rules(torLib):
    # add torsion rules
    # keep amide bonds plane 
    rule1 = "[O:1]=[C:2]-[N:3][!H:4] 0 180"
    torLib.AddTorsionRule(rule1)
    # keep amide bonds plane
    rule2 = "[O:1]=[C:2]-[N:3][H:4] 180"
    torLib.AddTorsionRule(rule2)
    # keep guanidinium groups plane
    rule3 = "[H:1][NHX3:2]-[CH0X3:3]=[NH2X3+,NHX2+0:4] 180"
    torLib.AddTorsionRule(rule3)
    # keep 2-Aminopyridine groups plane
    rule4 = "[H:1][NHX3:2]-[cX3:3][nHX3+,nX2+0:4] 0 180"
    torLib.AddTorsionRule(rule4)
    # keep dicarbonyl groups plane
    # ZINC000672399020
    rule5 = "[O:1]=[C:2]-[C:3]=[O:4] 180"
    torLib.AddTorsionRule(rule5)
    # ZINC000726467063
    rule6 = "[S:1]=[C:2]-[N:3][!H:4] 0 180"
    torLib.AddTorsionRule(rule6)
    rule7 = "[S:1]=[C:2]-[N:3][H:4] 180"
    torLib.AddTorsionRule(rule7)
    # phenol
    rule8 = "[H:1][O:2]-[c:3][c:4] 0 180"
    torLib.AddTorsionRule(rule8)

def generate_conformations(mol, h):
    
    # Parse arguments here
    #infile = sys.argv[1]

    # parameters related to omega
    # set omega energy window, if it equals 0, rotatable-bond-dependent window method.
    if os.environ.get('OMEGA_ENERGY_WINDOW', '').strip() != '':
        energyWindow = int(os.environ['OMEGA_ENERGY_WINDOW'])
        if energyWindow == 0:
            energyRangeFlag = True
        else:
            energyRangeFlag = False
    else:
        energyWindow = 12

    # set omega max number of confs, if it equals 0, rotatable-bond-dependent window method.
    if os.environ.get('OMEGA_MAX_CONFS', '').strip() != '':
        limitConfs = int(os.environ['OMEGA_MAX_CONFS'])
        if limitConfs == 0:
            limitConfsRangeFlag = True
        else:
            limitConfsRangeFlag = False
    else:
        limitConfs = 200

    # set the omega rmsd for clustering and filtering conformations, if it equals 0, rotatable-bond-dependent window method.
    if os.environ.get('OMEGA_RMSD', '').strip() != '':
        rmsd = float(os.environ['OMEGA_RMSD'])
        if rmsd == 0.0:
            rmsdRangeFlag = True
        else:
            rmsdRangeFlag = False
    else:
        rmsd = 0.5

    # set the omega torsion library: 1) Original; 2) GubaV21
    if os.environ.get('OMEGA_TORLIB', '').strip() != '':
        torLibType = str(os.environ['OMEGA_TORLIB'])
    else:
        torLibType = "GubaV21"

    # set the omega force field. Options are in the link below
    # https://docs.eyesopen.com/toolkits/cpp/oefftk/OEFFConstants/OEMMFFSheffieldFFType.html#OEFF::OEMMFFSheffieldFFType::MMFF94Smod
    if os.environ.get('OMEGA_FF', '').strip() != '':
        ffType = str(os.environ['OMEGA_FF'])
    else:
        ffType = "MMFF94Smod"

    # set the flag if omega uses hard-coded torsion patterns
    if os.environ.get('OMEGA_HARD_CODED_TOR_PATTERN', '').strip() != '':
        restrainTorPattern = bool(os.environ['OMEGA_HARD_CODED_TOR_PATTERN'])
    else:
        restrainTorPattern = True

    numHs = h
    if numHs > 5:
        raise ValueError("Refusing to build conformations with >5 rotatable hydrogens")
    if numHs >= 4:
        ##logging.warn("4-5 Rotatable hydrogens  reported. Reducing confs by a factor of 30")
        limitConfs = limitConfs // 30
    elif numHs >= 2:
        #logging.warn("2-3 Rotatable hydrogens  reported. Reducing confs by a factor of 3")
        limitConfs = limitConfs // 3

    TorDriveOpts = oeomega.OETorDriveOptions()
    set_TorDriveOpts_defaults(TorDriveOpts, torLibType)
    #logging.warning('Setting torLibType to %s' % (torLibType))

    # use restrainTorPatterns
    #logging.warning('Setting restrainTorPattern to %s' % str(restrainTorPattern))
    if restrainTorPattern == True:
        torLib = TorDriveOpts.GetTorLib()
        add_torison_rules(torLib)

    print(f'limitConfs={limitConfs}, energyWindow={energyWindow}, rmsd=${rmsd}, rangeflags=[energy={energyRangeFlag}, confs={limitConfsRangeFlag}, rmsd={rmsdRangeFlag}]')

    SliceEnsembleOpts = oeomega.OESliceEnsembleOptions()
    set_SliceEnsembleOpts_defaults(SliceEnsembleOpts, limitConfs, energyWindow, rmsd, energyRangeFlag, limitConfsRangeFlag, rmsdRangeFlag)
    #logging.warning('Setting energy window to %d; max confs to %d; rmsd to %.2f; energyRangeFlag to %s; limitConfsRangeFlag to %s; rmsdRangeFlag to %s.' % (energyWindow, limitConfs, rmsd, str(energyRangeFlag), str(limitConfsRangeFlag), str(rmsdRangeFlag)))

    tordriver = oeomega.OETorDriver(TorDriveOpts)
    tordriver.SetSliceEnsembleOptions(SliceEnsembleOpts)

    ffOpts = oeff.OEMMFFSheffieldOptions()
    set_ffOpts_defaults(ffOpts,ffType)
    #logging.warning('Setting ffType to %s' % (ffType))
    ff = oeff.OEMMFFSheffield(ffOpts)
    tordriver.SetForceField(ff)
    omegaFixOpts = oeomega.OEConfFixOptions()

    oechem.OEDetermineComponents(mol)
    count, ringlist = oechem.OEDetermineRingSystems(mol)
    
    print(ringlist)
    print(count)
    #rcf = open('ring_count', 'w')
    #rcf.write('%d\n' % count)
    #rcf.close()

    db2ins = []
    outfiles = []
    fail_count = 0

    if count == 0:
        SliceEnsembleOpts.SetMaxConfs(30)
        tordriver.SetSliceEnsembleOptions(SliceEnsembleOpts)
        ret_code = tordriver.GenerateConfs(mol)

        if ret_code == oeomega.OEOmegaReturnCode_Success:
            db2ins.append(mol)
        else:
            fail_count +=1 
    else:
        pred = oechem.OEPartPredAtom(ringlist)

        for i in range(1, count+1):

            pred.SelectPart(i)
            print(pred)

            #ringmol = oechem.OEMol()
            molcopy = oechem.OEMol(mol)
            #oechem.OESubsetMol(ringmol, molcopy, pred)
            ##print(i, "->", oechem.OEMolToSmiles(ringmol))
            omegaFixOpts.SetFixPredicate(pred)
            ret_code = tordriver.GenerateConfs(molcopy,omegaFixOpts)
            if ret_code == oeomega.OEOmegaReturnCode_Success:
                db2ins.append(molcopy)
            else:
                fail_count +=1
                
    return db2ins, fail_count

# class for handling the conversion between Mol2 classes needed by various stages of the pipeline
class MultiMol2:
    def __init__(self, data):
        self.data = data
        self.dockFormat = Mol2(mol2text=[line+'\n' for line in data.split('\n')])
        ifs = oechem.oemolistream()
        ifs.SetFormat(oechem.OEFormat_MOL2)
        ifs.openstring(data)
        ms = [oechem.OEMol(m) for m in ifs.GetOEMols()]
        self.oeFormat = ms[0]
        ifs.close() # even though we're just reading from memory, i'm still paranoid this needs to be called
    @staticmethod
    def oe2str(oemol):
        ofs = oechem.oemolostream()
        ofs.SetFormat(oechem.OEFormat_MOL2)
        ofs.openstring()
        oechem.OEWriteMolecule(ofs, oemol)
        s = ofs.GetString().decode('utf-8')
        ofs.close()
        return s    
    @staticmethod
    def oe2dock(oemol):
        s = MultiMol2.oe2str(oemol).split('\n')
        return Mol2(mol2text=[line+'\n' for line in s])
    @staticmethod
    def oe2supplier(oemol):
        s = MultiMol2.oe2str(oemol)
        return Mol2MolSupplier_noF(s) # new function added to TL_Functions.py, creates a mol supplier object without reading from a file

if __name__ == "__main__":

    infile = sys.argv[1]
    inroot, inext = os.path.splitext(infile)

    with open(infile) as mol2_data:
        mol = MultiMol2(mol2_data.read())

    h = count_hydrogens(mol.dockFormat)
    print(h)

    db2ins, fail_count = generate_conformations(mol.oeFormat, h)

    for i, db2in in enumerate(db2ins):
        outfile = "{}.{}.db2in.mol2".format(inroot, i+1)
        print(outfile)
        write_molecule(outfile, db2in)
