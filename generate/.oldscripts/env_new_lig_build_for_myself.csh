setenv SCRATCH_DIR /scratch
if ( ! -d $SCRATCH_DIR ) then
    setenv SCRATCH_DIR /tmp
endif

#setenv DOCKBASE /nfs/soft/dock/versions/dock37/DOCK-3.7-trunk
# debugging stage
setenv DOCKBASE /mnt/nfs/home/jklyu/zzz.github/DOCK3/DOCK
setenv AMSOLEXE /nfs/soft/amsol/in-house/amsol7.1-colinear-fix/amsol7.1

# Experimental changes to DOCK ligand pipeline
setenv EMBED_PROTOMERS_3D_EXE $DOCKBASE/ligand/3D/embed3d_corina.sh

# parameters related to omega
# set omega energy window, if it equals 0, rotatable-bond-dependent window method.
setenv OMEGA_ENERGY_WINDOW 12
# set omega max number of confs, if it equals 0, rotatable-bond-dependent window method.
setenv OMEGA_MAX_CONFS 600
# set the omega torsion library: 1) Original; 2) GubaV21
#setenv OMEGA_TORLIB GubaV21
setenv OMEGA_TORLIB Original
# set the omega force field. Options are in the link below
# https://docs.eyesopen.com/toolkits/cpp/oefftk/OEFFConstants/OEMMFFSheffieldFFType.html#OEFF::OEMMFFSheffieldFFType::MMFF94Smod
setenv OMEGA_FF MMFF94Smod
# set the omega rmsd for clustering and filtering conformations, if it equals 0, rotatable-bond-dependent window method.
setenv OMEGA_RMSD 0.5
# set the flag if omega uses hard-coded torsion patterns
setenv OMEGA_HARD_CODED_TOR_PATTERN True

# Dependencies
#source /nfs/soft/python/envs/dock/python-2.7/env.csh
#source /nfs/soft/corina/current/env.csh
source /nfs/soft/openbabel/openbabel-2.3.2/env.csh
# activate the openeye license
setenv OE_LICENSE /nfs/soft/openeye/oe_license.txt
#source /nfs/soft/openeye/openeye-201505/env.csh
source /nfs/soft/jchem/jchem-19.15/env.csh
#source /nfs/soft/jchem/jchem-15.11.23.0/env.csh

set LIMIT_JAVA="${DOCKBASE}/common/java-thread-limiter/mock-num-cpus 2"
setenv CXCALCEXE "${LIMIT_JAVA} `which cxcalc `"
setenv MOLCONVERTEXE "${LIMIT_JAVA} `which molconvert`"

setenv PATH "${PATH}:${DOCKBASE}/bin"

# activate the conda env
#source /mnt/nfs/home/jklyu/anaconda3/etc/profile.d/conda.csh
conda activate lig_build_py3-3.7
