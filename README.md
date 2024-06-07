3D Building Pipeline
====

This code will generate db2 files from smiles/mol2.


Required Software
===
- CORINA
- Openeye OMEGA (code is setup to use the openeye python toolkit, could be reworked to use command line version)
- Chemaxon JCHEM (cxcalc and molconvert)
- Openbabel
- Anaconda (Miniconda 3)

Getting Started
===
  1. Fork this repository
  2. Clone directory, cd into created directory
  4. Run install.sh 
  5. Install conda to either a shared nfs path or the '/soft' directory
  6. Create the conda environment using 'conda env create -f environment.yml -n envname'
  7. Move required software to '/soft' directory. Software must be .tar.gz format
  8. Move required licenses to '/licenses' directory (Jchem/Omega)
  9. Change paths in env.sh file (required changes are lines that are commented out)
  10. Modify the 'submit/build-3d.bash' file, source the conda environment on the first line

Running The Pipeline
===  
  1. Source the conda environment
  2. Set optional/required environment variables ([see here](https://wiki.docking.org/index.php/Building_The_3D_Pipeline_ZINC22))
  3. Run 'bash submit/submit-all.bash' to start job 

