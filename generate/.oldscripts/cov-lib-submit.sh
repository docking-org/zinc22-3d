#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N cov-lib

#example usage: qsub -t 1-filenum $DOCKBASE/ligand/generate/cov-lib-submit.sh /abs/path/to/input/smiles/

PERSIST=$1
INPUT_FILES=$PERSIST/input
COMPLETE=$PERSIST/output
TASK_INPUT=$( ls $INPUT_FILES | sed -n ${SGE_TASK_ID}p )
TASK_NAME=$( basename $TASK_INPUT )
TASK_FILE=$INPUT_FILES/$TASK_INPUT
ZINC_STORE=/nfs/db/zincraw

SCRATCH_DIR=/scratch
if [ ! -d $SCRATCH_DIR ]; then
    SCRATCH_DIR=/tmp
fi

source /nfs/soft/openeye/current/env.sh
source /nfs/soft/jchem/current/env.sh
# Need in-development ZINC libraries
source /nfs/soft/python/envs/dock/latest/env.sh
export OBABELBASE=/nfs/soft/openbabel/current
export AMSOLEXE=/nfs/soft/amsol/in-house/amsol7.1-colinear-fix/amsol7.1
export DOCKBASE=/mnt/nfs/home/londonir/code/git_trunk/DOCK

# Limit chemaxon programs to a single CPU
LIMIT_JCHEM="$DOCKBASE/common/on-one-core"
export CXCALCEXE="$LIMIT_JCHEM /nfs/soft/jchem/current/bin/cxcalc"
export MOLCONVERTEXE="$LIMIT_JCHEM /nfs/soft/jchem/current/bin/molconvert"

TASK_DIR=$SCRATCH_DIR/$( whoami )/$JOB_ID/$TASK_NAME

echo "Source file:" $TASK_FILE 1>&2
echo "Build dir:" $( hostname ):$TASK_DIR 1>&2

mkdir -pv $TASK_DIR 1>&2
pushd $TASK_DIR 1>&2
$DOCKBASE/ligand/generate/build_smiles_ligand.sh $TASK_FILE --covalent
popd

mkdir -pv $COMPLETE 1>&2

if [ -e $TASK_DIR/*.gz ]; then
mkdir -pv $COMPLETE/$TASK_NAME 
mv -v $TASK_DIR/*.gz $COMPLETE/$TASK_NAME
else
mkdir -pv $COMPLETE/failed
mv -v $TASK_DIR $COMPLETE/failed/$TASK_NAME 1>&2
fi

rm -rvf $TASK_DIR 1>&2

