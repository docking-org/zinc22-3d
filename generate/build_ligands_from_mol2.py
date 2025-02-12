# author: Benjamin Tingle 8/17/2020
import sys
import os
import subprocess
import shutil
import tarfile
import io
import openbabel
import time
import signal
import hashlib
from openeye import oechem

DOCKBASE='/'.join(__file__.split("/")[0:-3])
sys.path.append(DOCKBASE + "/ligand/mol2db2_py3_strain")
sys.path.append(DOCKBASE + "/ligand/strain")
sys.path.append(DOCKBASE + "/ligand/omega")

from mol2 import Mol2
from hydrogens import count_hydrogens
from mol2db2 import mol2db2_quick
from Torsion_Strain import calc_strain
from TL_Functions import Mol2MolSupplier_noF, Mol2MolSupplier
from omega import generate_conformations

# initialize working directories
subprocess.call(["mkdir", "solv"])
subprocess.call(["mkdir", "3d"])

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
    def mol2smiles(mol2):
        ifs = oechem.oemolistream()
        ifs.SetFormat(oechem.OEFormat_MOL2)
        ofs = oechem.oemolostream() 
        ofs.SetFormat(oechem.OEFormat_CAN)
        ifs.openstring(mol2)
        m = list(ifs.GetOEMols())[0]
        ofs.openstring()
        oechem.OEWriteMolecule(ofs, m)
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

mol2s_in = sys.argv[1:]

mols_in = []
i = 0
for mol2 in mol2s_in:
    names, mols = Mol2MolSupplier(mol2, dontconvert=True)
    for name in names:
        mf = open(f"3d/{i}.mol2", 'w')
        mf.write(mols[name])
        mf.close()
        smiles = "none"
        mols_in.append((mf.name, name, smiles))
        i += 1

start = time.time()

# read mol2 data into memory
mol2_data = []
solv_data = []
for i, entry in enumerate(mols_in):
    mol2, name, smiles = entry
    prot_id = 0

    # calculate solvation as we go along
    solv_dir = "solv/" + str(i)
    if not os.path.isdir(solv_dir):
        subprocess.call(["mkdir", solv_dir])
        mol2_fullname = '.'.join([name, 'mol2'])
        subprocess.call(["ln", "-svfn", "../../" + mol2, solv_dir + "/" + mol2_fullname])
        os.chdir(solv_dir)
        subprocess.call(["/bin/csh", DOCKBASE + "/ligand/amsol/calc_solvation.csh", mol2_fullname])
        os.chdir("../..")

    if os.path.isfile(solv_dir + "/output.solv"):

        with open(solv_dir + "/output.solv") as solv:
            solv_data.append(solv.read())

        with open(solv_dir + "/output.mol2") as mol2_f:
            mol = MultiMol2(mol2_f.read())
            mol.idx = i
            mol.name = name
            mol.smiles = smiles
           
            mol.prot_id = prot_id
            mol.charge = int(float(solv_data[-1].split('\n')[0].split()[2]))
            mol2_data.append(mol)

t_solvation = (time.time() - start)
print("{} / {} solvation successful".format(len(mol2_data), len(mol2s_in)))

# hard-code this bad boy in here
mp_range = ["M500","M400","M300","M200","M100","M000","P000","P010","P020","P030","P040","P050","P060","P070","P080","P090","P100","P110","P120","P130","P140","P150","P160","P170","P180","P190","P200","P210","P220","P230","P240","P250","P260","P270","P280","P290","P300","P310","P320","P330","P340","P350","P360","P370","P380","P390","P400","P410","P420","P430","P440","P450","P460","P470","P480","P490","P500","P600","P700","P800","P900"]
def get_zinc_directory_hash(zinc_id):
    if not zinc_id.startswith("ZINC"):
        return "."
    digits  = "0123456789abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ"
    hp_b62  = zinc_id[4:6]
    si_b62  = zinc_id[6:16]
    hash_l1 = zinc_id[14:16]
    hash_l2 = zinc_id[12:14]
    h, p    = digits.index(hp_b62[0]), digits.index(hp_b62[1])
    tranche = "H{0:>02d}".format(h) + mp_range[p]
    return '/'.join([tranche[0:3], tranche, hash_l2, hash_l1])

def convert(data, inf, otf):
    mol = openbabel.OBMol()
    convert = openbabel.OBConversion()
    convert.SetInAndOutFormats(inf, otf)
    convert.ReadString(mol, data)
    return convert.WriteString(mol)

def write_to_tarball(ball, data, name):
    tar = tarfile.TarInfo(name=name)
    tar.size = len(data)
    ball.addfile(tar, io.BytesIO(data))

tar_name = lambda x, *y:((x+'/') if x else '')+'.'.join([*y])

if os.path.isfile("output.tar.gz"):
    print("found ")
    subprocess.call(["mv", "output.tar.gz", "restart.tar.gz"])

stop = False

with tarfile.open("output.tar.gz", mode='w:gz') as output:

    # make sure to handle the SIGUSR1 interrupt so we properly close the tarfile
    def timeup_handler(signum, frame):
        global stop
        print("received SIGUSR1 : build_ligands.py")
        stop = True
        #output.close()
        #sys.exit(1)

    signal.signal(10, timeup_handler)
    signal.signal(2, timeup_handler)

    # allow for restartability
    # really need to rewrite this at some point...
    if os.path.isfile("restart.tar.gz"):
        with tarfile.open("restart.tar.gz", mode='r:gz') as restart:
            found_mol2s = []
            dest_mol2_map = {
                tar_name(get_zinc_directory_hash(mol.name), '.'.join([mol.name, str(mol.prot_id), chr(mol.charge+78)]), 'mol2') : mol for mol in mol2_data
            }
            for name in restart.getnames():
                with restart.extractfile(name) as memberfile:
                    data = memberfile.read()
                    write_to_tarball(output, data, name)
                # remove any mol2s from the worklist that already exist in the restart tarball
                if name.endswith('.mol2'):
                    found_mol2s.append(name)
            # set() fixes a weird bug where distinct protomers are generated but given id 0, causing identical mol2s in the output
            # update from the future: this is fixed at the protomer generation step now, this fix is no longer necessary
            found_mol2s = list(set(found_mol2s))
            for mol2 in [dest_mol2_map[m] for m in found_mol2s]:
                mol2_data.remove(mol2)

    t_strain_tot = t_db2_tot = t_omega_tot = 0

    for i, mol in enumerate(mol2_data):
        zinc_hash = get_zinc_directory_hash(mol.name)
        mol_fullname = '.'.join([mol.name, chr(mol.charge+78)])
        print(zinc_hash + '/' + mol_fullname)

        # new wrapper function added to hydrogens.py, count_hydrogens
        h = count_hydrogens(mol.dockFormat)
        # generate_conformations is a new function in the new file omega.py
        start = time.time()
        db2ins, fails = generate_conformations(mol.oeFormat, h)
        t_omega_tot += (time.time() - start)
        assert(fails == 0)
        db2ins = list(map(lambda mo:(MultiMol2.oe2dock(mo), MultiMol2.oe2supplier(mo)), db2ins))

        db2ins.append((MultiMol2.oe2dock(mol.oeFormat), MultiMol2.oe2supplier(mol.oeFormat))) # include the original conformation as a "separate" db2

        db2_all_data = ""
        for j, db2in in enumerate(db2ins):
            if stop == True:
                output.close()
                sys.exit()

            db2in_standard, db2in_strain = db2in
            # new wrapper function added to Torsion_Strain.py, calc_strain
            start = time.time()
            tE, pE = calc_strain(*db2in_strain)
            t_strain_tot += (time.time() - start)
            # new utility function added to mol2.Mol2 class, addStrainInfo
            db2in_standard.addStrainInfo(tE, pE)
            # new wrapper function added to mol2db2.py, mol2db2_quick
            start = time.time()
            db2_data = mol2db2_quick(db2in_standard, solvfile="solv/" + str(mol.idx) + "/output.solv", clashfile=DOCKBASE + "/ligand/mol2db2/clashfile.txt")
            db2_all_data += db2_data
            t_db2_tot += (time.time() - start)

            print(j+1, '/', len(db2ins))

        pdbqt_data = convert(mol2_data[i].data, 'mol2', 'pdbqt')
        sdf_data   = convert(mol2_data[i].data, 'mol2', 'sdf')

        write_to_tarball(output, solv_data[i].encode('utf-8'),      name=tar_name(zinc_hash, mol_fullname, 'solv'))
        write_to_tarball(output, mol2_data[i].data.encode('utf-8'), name=tar_name(zinc_hash, mol_fullname, 'mol2'))
        write_to_tarball(output, sdf_data.encode('utf-8'),          name=tar_name(zinc_hash, mol_fullname, 'sdf'))
        write_to_tarball(output, pdbqt_data.encode('utf-8'),        name=tar_name(zinc_hash, mol_fullname, 'pdbqt'))
        write_to_tarball(output, db2_all_data.encode('utf-8'),      name=tar_name(zinc_hash, mol_fullname, 'db2.gz'))

    write_to_tarball(output, 'version={}'.format(os.environ.get("DOCK_VERSION")).encode('utf-8'), name=".dock_version")
    print("elapsed times:")
    print("solvation: {}".format(t_solvation))
    print("strain:    {}".format(t_strain_tot))
    print("db2:       {}".format(t_db2_tot))
    print("omega:     {}".format(t_omega_tot))
# that's all!
