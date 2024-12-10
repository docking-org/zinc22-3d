# 3D Building Pipeline

This code will generate db2 files from smiles/mol2.

# Required Software

- CORINA
- Openeye OMEGA (code is setup to use the openeye python toolkit, could be reworked to use command line version)
- Chemaxon JCHEM (cxcalc and molconvert)
- Openbabel
- Anaconda (Miniconda 3)

# Getting Started

1. Fork this repository
2. Clone directory, cd into created directory
3. Run install.sh
4. Install conda to either a shared nfs path or the '/soft' directory
5. Create the conda environment using 'conda env create -f environment.yml -n envname'
6. Move required software to '/soft' directory. Software must be .tar.gz format
7. Move required licenses to '/licenses' directory (Jchem/Omega). Software looks for jchem-license.cxl and oe-license.txt by default. These can be changed in the env.sh files
8. Change paths in env.sh file (required changes are lines that are commented out)
9. Modify the 'submit/build-3d.bash' file, source the conda environment on the first line

# Setting Up Software

In the soft folder, you must add the following folders

1. Openbabel - https://github.com/openbabel/openbabel/releases/tag/openbabel-3-1-1 (you must clone the repository and compile it from the source code)
2. JChem - https://download.chemaxon.com/jchem-engines
3. CORINA - containing only the corina executable
4. Omega - containing the
5. extralibs - a folder containing extra lib files that your system needs. This can vary from system to system. Common files include libg2c.so, libz.so, etc.

All software folders will be compressed through the scripts. In the case that these aren't automatically compressed, run;

- tar -czvf FOLDER_NAME.tar.gz FOLDER_NAME

for each folder

# Running The Pipeline

1. Source the conda environment
2. Set optional/required environment variables ([see here](https://wiki.docking.org/index.php/Building_The_3D_Pipeline_ZINC22))
3. Run 'bash submit/submit-all.bash' to start job
