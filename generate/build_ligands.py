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
from openeye import oechem

DOCKBASE='/'.join(__file__.split("/")[0:-3])
sys.path.append(DOCKBASE + "/ligand/mol2db2_py3_strain")
sys.path.append(DOCKBASE + "/ligand/strain")
sys.path.append(DOCKBASE + "/ligand/omega")

from mol2 import Mol2
from hydrogens import count_hydrogens
from mol2db2 import mol2db2_quick
from Torsion_Strain import calc_strain
from TL_Functions import Mol2MolSupplier_noF
from omega import generate_conformations

# read in protomers
protomers = sys.argv[1]
protonated_flat = []
with open(protomers) as prot_f:
    for line in prot_f:
        smiles, name, prot_id = line.split()
        protonated_flat.append((name + "." + str(prot_id), smiles, int(prot_id)))

# sort the protonated list for the next step
protonated_flat = sorted(protonated_flat, key=lambda x:x[0])

# when protomers get expanded, expanded protomers will not be assigned a unique id
# e.g if a new protomer was expanded from ZINC123.1, the new protomer will still be named ZINC123.1
# we fix this here
protomers_aug = open(os.path.dirname(protomers) + "/protomers_in_corina", 'w')
name_prev = None
protonated_flat_aug = []
for protomer in protonated_flat:
    name = protomer[0]
    if name_prev == name:
        protid = int(name.split('.')[-1])
        name = '.'.join(name.split('.')[:-1] + [str(protid + 1)])
        protomers_aug.write(protomer[1] + " " + name + "\n")
        protonated_flat_aug.append((name, protomer[1], protid + 1))
        name_prev = name
    else:
        protomers_aug.write(protomer[1] + " " + name + "\n")
        protonated_flat_aug.append(protomer)
        name_prev = name
        dup_protid_cnt = 0
protomers_aug.close()
protonated_flat = protonated_flat_aug

# initialize working directories
subprocess.call(["mkdir", "3d"])
subprocess.call(["mkdir", "solv"])

# read in number of corina conformations to (attempt to) generate per protomer
if not os.environ.get("CORINA_MAX_CONFS"):
    os.environ["CORINA_MAX_CONFS"] = "1"
corinaconfs = int(os.environ["CORINA_MAX_CONFS"])

# create 3d embeddings for protomers
subprocess.call([DOCKBASE + "/ligand/3D/embed3d_corina.sh", protomers_aug.name, "-o", "3d/"])
# remove size zero files created by spliton.py (not sure why it does this)
with subprocess.Popen(["find", "3d", "-size", "0"], stdout=subprocess.PIPE) as proc:
    subprocess.call(["xargs", "rm"], stdin=proc.stdout)

# find the names of every mol that succeeded corina 3d embedding
corina_success = []
for mol2fn in os.listdir("3d"):
    with subprocess.Popen(["awk", "{if (NR == 2) print $0}", "3d/" + mol2fn], stdout=subprocess.PIPE) as proc:
        mol_name = proc.stdout.read().strip().decode('utf-8')
        corina_success.append((mol_name, mol2fn))
# sort the list by mol number so that it lines up with the list of protonated mols
corina_success = sorted(corina_success, key=lambda x:int(x[1]))

# rewrites the name of a mol2, only used when CORINA_CONFS > 1 to denote different corina conformations
def rewrite_mol2_name(mol2fn, name):
    mol2_f = open(mol2fn, 'r')
    mol2_t = open(mol2fn + '.t', 'w')
    i = 0
    for line in mol2_f:
        i += 1
        if i == 2:
            mol2_t.write(name + '\n')
        else:
            mol2_t.write(line)
    mol2_f.close()
    mol2_t.close()
    os.rename(mol2fn + '.t', mol2fn)

# identify which molecules succeeded, assign corina confs numbers (if CORINA_CONFS > 1)
corina_success_names = [c[0] for c in corina_success]
protonated_success = []
total_success = len(protonated_flat)
total_protomers = len(protonated_flat)

for p in protonated_flat:
    try:
        ci = corina_success_names.index(p[0])
    except:
        total_success -= 1
        continue

    confnum = 0
    while ci >= 0:        
        mol_name, mol2fn = corina_success[ci]
        corina_success_names.pop(ci)
        corina_success.pop(ci)
        if corinaconfs > 1:
            mol_name = '.'.join([p[0], str(confnum)])
            p1 = (mol_name, p[1], p[2])
            rewrite_mol2_name("3d/" + mol2fn, mol_name)
            protonated_success.append((mol2fn, p1))
        else:
            protonated_success.append((mol2fn, p))
        confnum += 1
        try:
            ci = corina_success_names.index(p[0])
        except:
            ci = -1
protonated_flat = protonated_success

if corinaconfs > 1:
    print("{} / {} protomers built successfully, {} conformations generated (avg {} confs per protomer, mc={})".format(total_success, total_protomers, len(protonated_flat), len(protonated_flat) / total_success, corinaconfs))
else:
    print("{} / {} protomers built successfully".format(total_success, total_protomers))

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

start = time.time()

# read mol2 data into memory
mol2_data = []
solv_data = []
for mol2, protomer in protonated_flat:
    name, smiles, prot_id = protomer

    # calculate solvation as we go along
    if not os.path.isdir("solv/" + str(mol2)):
        subprocess.call(["mkdir", "solv/" + str(mol2)])
        mol2_fullname = '.'.join([name, 'mol2'])
        subprocess.call(["ln", "-svfn", os.getcwd() + "/3d/" + str(mol2), "solv/" + str(mol2) + "/" + mol2_fullname])
        os.chdir("solv/" + str(mol2))
        subprocess.call(["/bin/csh", "-f", DOCKBASE + "/ligand/amsol/calc_solvation.csh", mol2_fullname])
        os.chdir("../..")

    if os.path.isfile("solv/" + str(mol2) + "/output.solv"):

        with open("solv/" + str(mol2) + "/output.solv") as solv:
            solv_data.append(solv.read())

        with open("solv/" + str(mol2) + "/output.mol2") as mol2_f:
            mol = MultiMol2(mol2_f.read())
            mol.idx = mol2
            mol.name = name
            mol.smiles = smiles
            mol.prot_id = prot_id
            mol.charge = int(float(solv_data[-1].split('\n')[0].split()[2]))
            mol2_data.append(mol)

t_solvation = (time.time() - start)
print("{} / {} solvation successful".format(len(mol2_data), len(protonated_flat)))

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
    
    skip_omega = True if os.getenv("SKIP_OMEGA") else False

    for i, mol in enumerate(mol2_data):
        zinc_hash = get_zinc_directory_hash(mol.name)
        mol_fullname = '.'.join([mol.name, chr(mol.charge+78)])
        print(zinc_hash + '/' + mol_fullname)

        # new wrapper function added to hydrogens.py, count_hydrogens
        h = count_hydrogens(mol.dockFormat)
        # generate_conformations is a new function in the new file omega.py
        if not skip_omega:
            start = time.time()
            db2ins, fails = generate_conformations(mol.oeFormat, h)
            t_omega_tot += (time.time() - start)
            assert(fails == 0)
            db2ins = list(map(lambda mo:(MultiMol2.oe2dock(mo), MultiMol2.oe2supplier(mo)), db2ins))
        else:
            db2ins = [MultiMol2.oe2dock(mol.oeFormat), MultiMol2.oe2supplier(mol.oeFormat)]

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
