#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N build-protomers

set -e

if [ -z "${JOB_ID}" ]; then
	export JOB_ID=$$
fi

if [ -z "${SGE_TASK_ID}" ]; then
	export SGE_TASK_ID=1
fi

PERSIST=$1
INPUT_FILES=$PERSIST/input
COMPLETE=$PERSIST/output
TASK_INPUT=$( ls $INPUT_FILES | sed -n ${SGE_TASK_ID}p )
TASK_NAME=$( basename $TASK_INPUT )
TASK_FILE=$( readlink -f $INPUT_FILES/$TASK_INPUT )
ZINC_STORE=/nfs/work/teague/zinc/testing/rebuild

SCRATCH_DIR=/scratch
if [ ! -d $SCRATCH_DIR ]; then
    SCRATCH_DIR=/tmp
fi

source /nfs/soft/openeye/current/env.sh
source /nfs/soft/jchem/current/env.sh
# Need in-development ZINC libraries
source /nfs/soft/python/envs/dock/latest/env.sh
#source /nfs/work/teague/Projects/zinc/dev/env.sh
export ZINCSRC=/nfs/home/teague/Code/ZINC
export ZINCPY=/nfs/work/teague/Projects/zinc/dev/bin/python
export PYTHONPATH=$ZINCSRC
#export STOREPROTOMERSEXE="$ZINCPY $ZINCSRC/zinc/load/load_protomers.py --to=$ZINC_STORE --log-level=debug"
export OBABELBASE=/nfs/soft/openbabel/current
export AMSOLEXE=/nfs/soft/amsol/in-house/amsol7.1-colinear-fix/amsol7.1
export DOCKBASE=/nfs/home/teague/Code/DOCK

# Limit chemaxon programs to a single CPU
LIMIT_JCHEM="$DOCKBASE/ligand/protonate/limiter/limit_cpus 1"
export CXCALCEXE="$LIMIT_JCHEM /nfs/soft/jchem/current/bin/cxcalc"
export MOLCONVERTEXE="$LIMIT_JCHEM /nfs/soft/jchem/current/bin/molconvert"

TASK_DIR=$SCRATCH_DIR/$( whoami )/$JOB_ID/$TASK_NAME

echo "Source file:" $TASK_FILE 1>&2
echo "Build dir:" $( hostname ):$TASK_DIR 1>&2

mkdir -pv $TASK_DIR 1>&2
pushd $TASK_DIR 1>&2
$DOCKBASE/ligand/generate/build_database_ligand.sh $TASK_FILE
popd

mkdir -pv $COMPLETE 1>&2
mv -v $TASK_DIR $COMPLETE/$TASK_NAME 1>&2

rm -rvf $TASK_DIR 1>&2


