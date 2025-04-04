# author: Benjamin Tingle 8/17/2020
import time

import sys
import os

DOCKBASE=os.environ['DOCKBASE']
sys.path.append(DOCKBASE + "/ligand/mol2db2_py3_strain")
sys.path.append(DOCKBASE + "/ligand/strain")
sys.path.append(DOCKBASE + "/ligand/omega")

import subprocess
import shutil
import tarfile
import io
import openbabel
import signal
from openeye import oechem


from mol2 import Mol2
from hydrogens import count_hydrogens
from mol2db2 import mol2db2_quick
from Torsion_Strain import calc_strain
from TL_Functions import Mol2MolSupplier_noF
from omega import generate_conformations


start = time.time()

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

def write_to_tarball(ball, data, name):
    tar = tarfile.TarInfo(name=name)
    tar.size = len(data)
    ball.addfile(tar, io.BytesIO(data))


def load_neutral_smiles_dict(input_file='/data/input.smi'):
    smiles_dict = dict()
    with open(input_file) as f:
        for line1 in f:  # using line1, name1 because of global variable names... -.-
            smiles1, name1 = line1.split()
            if name1 in smiles_dict:
                raise ValueError(f"Name {name1} is not unique.")
            smiles_dict[name1] = smiles1
    return smiles_dict


neutral_smiles_dict = load_neutral_smiles_dict()

# TODO: add the neutral smiles to the "protonated flat" list
protonated_flat = []
prot_ids = {}
for line in sys.stdin:
    smiles, name = line.split()[:2]
    if name in prot_ids:
        prot_ids[name] += 1
    else:
        prot_ids[name] = 0
    prot_id = prot_ids[name]
    protonated_flat.append((name + "." + str(prot_id), smiles, int(prot_id)))

t_upstream = (time.time() - start)

start = time.time()

# sort the protonated list for the next step
protonated_flat = sorted(protonated_flat, key=lambda x:x[0])

# when protomers get expanded, expanded protomers will not be assigned a unique id
# e.g if a new protomer was expanded from ZINC123.1, the new protomer will still be named ZINC123.1
# we fix this here
# protomers_aug = open(os.getcwd() + "/protomers_in_corina", 'w')
corina_input = ""
name_prev = None
protonated_flat_aug = []
for protomer in protonated_flat:
    name = protomer[0]
    if name_prev == name:
        protid = int(name.split('.')[-1])
        name = '.'.join(name.split('.')[:-1] + [str(protid + 1)])
        corina_input += protomer[1] + " " + name + "\n"
        # protomers_aug.write(protomer[1] + " " + name + "\n")
        protonated_flat_aug.append((name, protomer[1], protid + 1))
        name_prev = name
    else:
        corina_input += protomer[1] + " " + name + "\n"
        # protomers_aug.write(protomer[1] + " " + name + "\n")
        protonated_flat_aug.append(protomer)
        name_prev = name
        dup_protid_cnt = 0
# protomers_aug.close()
protonated_flat = protonated_flat_aug

t_corina_aug = time.time() - start

start = time.time()

# initialize working directories
os.makedirs("3d", exist_ok=True)
os.makedirs("solv", exist_ok=True)

# read in number of corina conformations to (attempt to) generate per protomer
corinaconfs = int(os.environ["CORINA_MAX_CONFS"])

# create 3d embeddings for protomers
# subprocess.call([DOCKBASE + "/ligand/3D/embed3d_corina.sh", protomers_aug.name, "-o", "3d/"])
# Brendan rewrite in the function here to speed things up (avoid spawing sh script which i think can cause slow downs with lots of workers)
corina_result = subprocess.run([f'{os.environ["CORENAEXE"]}', "-i", "t=smiles", "-o", "t=mol2", "-d", "rc,flapn,de=6,mc=1,wh"], input=corina_input, text=True, capture_output=True)
subprocess.run([f'{os.environ["SPLITONEXE"]}', "#\tEnd of record", "3d/"], input=corina_result.stdout, text=True)

t_embed3d = time.time() - start

start = time.time()
# remove size zero files created by spliton.py (not sure why it does this)
# (The reason it does this is because when you get to the last delimeter it opens another file but then closes it because you got to end)
with subprocess.Popen(["find", "3d", "-size", "0"], stdout=subprocess.PIPE) as proc:
    subprocess.call(["xargs", "rm"], stdin=proc.stdout)

t_find = time.time() - start

start = time.time()
# find the names of every mol that succeeded corina 3d embedding
corina_success = []
for mol2fn in os.listdir("3d"):
    with subprocess.Popen(["awk", "{if (NR == 2) print $0}", "3d/" + mol2fn], stdout=subprocess.PIPE) as proc:
        mol_name = proc.stdout.read().strip().decode('utf-8')
        corina_success.append((mol_name, mol2fn))
# sort the list by mol number so that it lines up with the list of protonated mols
corina_success = sorted(corina_success, key=lambda x:int(x[1]))
t_find_success = time.time() - start 


start = time.time()
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
t_change_success = time.time() - start

# if corinaconfs > 1:
#     print("{} / {} protomers built successfully, {} conformations generated (avg {} confs per protomer, mc={})".format(total_success, total_protomers, len(protonated_flat), len(protonated_flat) / total_success, corinaconfs))
# else:
#     print("{} / {} protomers built successfully".format(total_success, total_protomers))

# t_corina = (time.time() - start)

start = time.time()

SOLV_TIMEOUT = 60

# read mol2 data into memory
mol2_data = []
solv_data = []
t_solvation = 0
failed_solvation_names = set()
for mol2, protomer in protonated_flat:
    name, smiles, prot_id = protomer
    if name in failed_solvation_names:
        continue
    # calculate solvation as we go along
    if not os.path.isdir("solv/" + str(mol2)):
        subprocess.call(["mkdir", "solv/" + str(mol2)])
        # mol2_fullname = '.'.join([name, 'mol2'])
        subprocess.call(["ln", "-svfn", os.getcwd() + "/3d/" + str(mol2), "solv/" + str(mol2) + "/temp.mol2"])
        os.chdir("solv/" + str(mol2))
        start = time.time()
        # This is replacing the old calc_solvation.csh call that I found had lots of overhead with many workers (should take ~15 sec but was getting ~30)
        # Run obabel
        obabel_result = subprocess.run([os.environ["OBABELEXE"], "-i", "mol2", "temp.mol2", "-o", "mopin"], capture_output=True, text=True)
        obabel_output = obabel_result.stdout
        # process input
        subprocess.run(["python", f'{os.environ["DOCKBASE"]}/ligand/amsol/make_amsol71_input_docker.py3.py', name], input=obabel_output, text=True)
        # run amsol
        with open('temp.in-wat', 'r') as infile, open('temp.o-wat', 'w') as outfile:
            # Some molecules break or make amsol take forever - calc_solvation.csh had a 1m timeout to stop this so I have one here too
            try:
                subprocess.run([os.environ["AMSOLEXE"]], stdin=infile, stdout=outfile, stderr=subprocess.PIPE, text=True, timeout=SOLV_TIMEOUT)
            except subprocess.TimeoutExpired:
                failed_solvation_names.add(name)
                os.chdir("../..")
                continue
        with open('temp.in-hex', 'r') as infile, open('temp.o-hex', 'w') as outfile:
            try:
                subprocess.run([os.environ["AMSOLEXE"]], stdin=infile, stdout=outfile, stderr=subprocess.PIPE, text=True, timeout=SOLV_TIMEOUT)
            except subprocess.TimeoutExpired:
                failed_solvation_names.add(name)
                os.chdir("../..")
                continue
        # process output
        try:
            subprocess.run(["python", f'{os.environ["DOCKBASE"]}/ligand/amsol/process_amsol_mol2.py3.py', "temp.o-wat", "temp.o-hex", "temp.mol2", "output"], check=True, capture_output=True, text=True)
        except Exception as e:
            print(f'processing amsol output failed with exception: {e}')
            failed_solvation_names.add(name)
            os.chdir("../..")
            continue

        os.chdir("../..")
        t_solvation += (time.time() - start)

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



# t_solvation = (time.time() - start)
# print("{} / {} solvation successful".format(len(mol2_data), len(protonated_flat)))

start = time.time()

stop = False

with tarfile.open("bundle.db2.tgz", mode='w:gz') as output:

    # make sure to handle the SIGUSR1 interrupt so we properly close the tarfile
    def timeup_handler(signum, frame):
        global stop
        # print("received SIGUSR1 : build_ligands.py")
        stop = True
        #output.close()
        #sys.exit(1)

    signal.signal(10, timeup_handler)
    signal.signal(2, timeup_handler)


    t_calc_strain = t_add_strain = t_db2_tot = t_omega_tot = t_count_h = t_db2ins = t_convert_write = 0
    
    skip_omega = True if os.getenv("SKIP_OMEGA") else False    

    t_get_to_loop = (time.time() - start)

    successfully_built = set()

    for i, mol in enumerate(mol2_data):
        if mol.name in failed_solvation_names:
            continue
        start = time.time()
        mol_fullname = '.'.join([mol.name, chr(mol.charge+78)])

        # new wrapper function added to hydrogens.py, count_hydrogens
        h = count_hydrogens(mol.dockFormat)
        t_count_h += (time.time() - start)
        # generate_conformations is a new function in the new file omega.py
        if not skip_omega:
            start = time.time()
            try:
                db2ins, fails = generate_conformations(mol.oeFormat, h)
            except Exception as e:
                print(f'Mol {mol.name} failed omega conformer generation with exception: {e}')
                continue
            t_omega_tot += (time.time() - start)
            assert(fails == 0)
            start = time.time()
            db2ins = list(map(lambda mo:(MultiMol2.oe2dock(mo), MultiMol2.oe2supplier(mo)), db2ins))
            t_db2ins += (time.time() - start)
        else:
            start = time.time()
            db2ins = [MultiMol2.oe2dock(mol.oeFormat), MultiMol2.oe2supplier(mol.oeFormat)]
            t_db2ins += (time.time() - start)

        db2_all_data = ""
        failed_strain_db2 = False
        for j, db2in in enumerate(db2ins):
            if stop == True:
                output.close()
                sys.exit()

            # TODO: Put the microspecies percentage in the db2

            db2in_standard, db2in_strain = db2in
            # new wrapper function added to Torsion_Strain.py, calc_strain
            start = time.time()
            try:
                tE, pE = calc_strain(*db2in_strain)
            except Exception as e:
                print(f"Mol {mol.name} failed strain calculation with exception: {e}")
                failed_strain_db2 = True
                break
            t_calc_strain += (time.time() - start)
            # new utility function added to mol2.Mol2 class, addStrainInfo
            start = time.time()
            db2in_standard.addStrainInfo(tE, pE)
            t_add_strain += (time.time() - start)
            # new wrapper function added to mol2db2.py, mol2db2_quick
            start = time.time()
            try:
                db2_data = mol2db2_quick(db2in_standard, solvfile="solv/" + str(mol.idx) + "/output.solv", clashfile=DOCKBASE + "/ligand/mol2db2/clashfile.txt")
            except Exception as e:
                print(f"Mol {mol.name} failed mol2db2_quick with exception: {e}")
                failed_strain_db2 = True
                break

            # Remove the protomer ID inside the DB2 file itself
            db2_data = db2_data[:2] + f"{mol.name.split('.')[0]:16}" + db2_data[18:]

            # Add the neutral smiles to the DB2
            shortname = mol.name.split('.')[0]
            neutral_smiles = neutral_smiles_dict[shortname] if shortname in neutral_smiles_dict else "N/A"
            if len(neutral_smiles) > 76:
                neutral_smiles = "smilestoolong"
            db2_data_list = db2_data.split('\n')
            db2_data_list[2] = f"M {neutral_smiles:76}"
            db2_data = '\n'.join(db2_data_list)

            db2_all_data += db2_data
            t_db2_tot += (time.time() - start)

            # print(j+1, '/', len(db2ins))

        if failed_strain_db2:
            continue

        start = time.time()
        write_to_tarball(output, db2_all_data.encode('utf-8'), name=mol_fullname + '.db2')
        successfully_built.add(shortname)

        t_convert_write += (time.time() - start)
        start = time.time()
        t_convert_write += (time.time() - start)

    start = time.time()
    t_convert_write += (time.time() - start)

    failed_out = open('/data/failed_to_build.smi', 'w')
    success_out = open('/data/successfully_built.smi', 'w')
    for name in neutral_smiles_dict:
        if name not in successfully_built:
            failed_out.write(f"{neutral_smiles_dict[name]} {name}\n")
        else:
            success_out.write(f"{neutral_smiles_dict[name]} {name}\n")
    failed_out.close()
    success_out.close()


# Everything is run in /tmp so put the final output to the correct places
start = time.time()
shutil.move('bundle.db2.tgz', '/data/bundle.db2.tgz')

t_cleanup = time.time() - start

print(
    "elapsed times:\n"
    "upstream: {}\n"
    "num_for_solvation: {}\n"
    "solvation: {}\n"
    "calc strain: {}\n"
    "add strain: {}\n"
    "db2: {}\n"
    "omega: {}\n"
    "corina_aug: {}\n"
    "corina_embed3d: {}\n"
    "corina_find: {}\n"
    "corina_find_success: {}\n"
    "corina_change_success: {}\n"
    "to_to_loop: {}\n"
    "count H: {}\n"
    "db2ins: {}\n"
    "con/write: {}\n"
    "cleanup_folders: {}".format(
        t_upstream,
        len(protonated_flat),
        t_solvation,
        t_calc_strain,
        t_add_strain,
        t_db2_tot,
        t_omega_tot,
        t_corina_aug,
        t_embed3d,
        t_find,
        t_find_success,
        t_change_success,
        t_get_to_loop,
        t_count_h,
        t_db2ins,
        t_convert_write,
        t_cleanup
    )
)
        

# that's all!
