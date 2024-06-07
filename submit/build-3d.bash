#!/bin/bash

#source your env.sh here
# req: JOB_ID // SLURM_ARRAY_JOB_ID
# req: SGE_TASK_ID // SLURM_ARRAY_TASK_ID
# req: INPUT
# req: OUTPUT
# req: SHRTCACHE
# req: LONGCACHE

# opt: WORK_DIR
cwd=$PWD
export SOFT_HOME=${SOFT_HOME-$HOME}
export SHRTCACHE=${SHRTCACHE-/tmp}
export LONGCACHE=${LONGCACHE-/tmp}
export COMMON_DIR=$LONGCACHE/build_3d_common_$(whoami)
#export COMMON_DIR=$LONGCACHE/build_3d_common/

# don't set default version in here- do it in submit-all.bash
export DOCK_VERSION=${DOCK_VERSION}
export CORINA_VERSION=${CORINA_VERSION}
export PYENV_VERSION=${PYENV_VERSION}
export JCHEM_VERSION=${JCHEM_VERSION}
export OPENBABEL_VERSION=${OPENBABEL_VERSION}
export EXTRALIBS_VERSION=${EXTRALIBS_VERSION}

export PYTHONBASE=$COMMON_DIR/$PYENV_VERSION
export DOCKBASE=$COMMON_DIR/${DOCK_VERSION}

failed=
for required_var in INPUT OUTPUT DOCK_VERSION CORINA_VERSION PYENV_VERSION JCHEM_VERSION OPENBABEL_VERSION EXTRALIBS_VERSION; do
	if [ -z ${!required_var} ]; then
		echo "missing $required_var!" 1>&2
		failed=1
	fi
done
! [ -z $failed ] && exit 1

if [ -f $SOFT_HOME/soft/${DOCK_VERSION}.tar.gz ]; then
	SOFT_HOME=$SOFT_HOME/soft
fi

#if ! [ -d $COMMON_DIR ]; then
#        mkdir -p $COMMON_DIR
#        chmod 777 $COMMON_DIR
#fi

mkdir -p $COMMON_DIR
chmod 777 $COMMON_DIR

JOB_ID=${SLURM_ARRAY_JOB_ID-$JOB_ID}
TASK_ID=${SLURM_ARRAY_TASK_ID-$SGE_TASK_ID}

[ -z $TASK_ID ] && echo "missing TASK_ID!" 1>&2 && exit 1
function log {
    echo "[build-3d $(date +%X)]: $@"
}
function extract_cmd {
	log "extracting $1"
	tar -C $COMMON_DIR -xzf $SOFT_HOME/$1.tar.gz && echo > $COMMON_DIR/$1/.done
}

# added an additional check to make sure software dir isn't empty, since this seems to have happened before
# "lib" is a bandaid to fix some libraries that weren't found- can probably include most of it with openbabel-install
for software in $DOCK_VERSION $PYENV_VERSION $EXTRALIBS_VERSION $OPENBABEL_VERSION $JCHEM_VERSION $CORINA_VERSION $EXTRALIBS_VERSION; do
	(
		flock -x 9
		if ! [ -f $COMMON_DIR/$software/.done ] || [ $(ls $COMMON_DIR/$software | wc -l) -eq 0 ]; then
			extract_cmd $software
		fi
		flock -u 9
	)9>$COMMON_DIR/install_${software}.lock
done

log $(hostname)

WORK_BASE=$SHRTCACHE/$(whoami)_build3d/${JOB_ID}/${TASK_ID}
mkdir -p $WORK_BASE

if [ -z $RESUBMIT ]; then
	SPLIT_FILE=$INPUT/`ls $INPUT | tr '\n' ' ' | cut -d' ' -f$TASK_ID`
	TARGET_FILE=$WORK_BASE/$TASK_ID
else
	SPLIT_FILE=$INPUT/resubmit/`ls $INPUT/resubmit | tr '\n' ' ' | cut -d' ' -f$TASK_ID`
	TARGET_FILE=$WORK_BASE/$(basename $SPLIT_FILE)
fi

cp $SPLIT_FILE $TARGET_FILE
log SPLIT_FILE=$SPLIT_FILE
cd $WORK_BASE

log PATH=$PATH
log $(hostname)
log "starting build job on $TARGET_FILE"
log "len($TARGET_FILE)=$(cat $TARGET_FILE | wc -l)"

# move logs to their final destination & clean up the working directory
cleanup() {
        if [ -z $SKIP_DELETE ] && [ -d $WORK_BASE ]; then rm -r $WORK_BASE; fi
}

trap cleanup EXIT

# save our progress if we've reached the time limit. DOCK has not been modified to take advantage of this, so this doesn't do much as of yet
# but this will be useful for re-doing as little work as possible in the future
reached_time_limit() {
	pushd $WORK_BASE
        log "time limit reached! saving progress..."
        tar -czf $(basename $TARGET_FILE).save.tar.gz .
        mkdir -p $OUTPUT/save
        mv $(basename $TARGET_FILE).save.tar.gz $OUTPUT/save
        popd $WORK_BASE
	cleanup 99
}

# on sge, SIGUSR1 is sent once a job surpasses it's "soft" time limit (-l s_rt=XX:XX:XX), usually specified a minute or two before the hard time limit (-l h_rt=XX:XX:XX) where SIGTERM is sent
# on slurm, the same can be achieved by adding the --signal=USR1@60 option to your sbatch args to send the SIGUSR1 signal 1 minute before the job is terminated
# trap reached_time_limit SIGUSR1

# jobs that have output already shouldn't be resubmitted, but this is just in case that doesn't happen
if [ -f $OUTPUT/$(basename $TARGET_FILE).tar.gz ]; then
        log "results already present in $OUTPUT_BASE for this job, exiting..."
        cleanup 0
fi

# un-archive our saved progress (if any) into the current working directory
if [ -f $OUTPUT/save/$(basename $TARGET_FILE).save.tar.gz ]; then
        echo "saved progress found!"
	tar -xzf $OUTPUT/save/$(basename $TARGET_FILE).save.tar.gz .
        rm $OUTPUT/save/$(basename $TARGET_FILE).save.tar.gz
fi

##### start initialize environment. Moving this from env_new_lig_build.sh to here so everything fits in one file

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$COMMON_DIR/$OPENBABEL_VERSION/lib
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$COMMON_DIR/$EXTRALIBS_VERSION
log LD_LIBRARY_PATH=$LD_LIBRARY_PATH

#export AMSOLEXE=$SOFTBASE/amsol/in-house/amsol7.1-colinear-fix/amsol7.1

# Experimental changes to DOCK ligand pipeline
export EMBED_PROTOMERS_3D_EXE=$DOCKBASE/ligand/3D/embed3d_corina.sh
# parameters related to omega
# set omega energy window, if it equals 0, rotatable-bond-dependent window method.
export OMEGA_ENERGY_WINDOW=${OMEGA_ENERGY_WINDOW-12}
# set omega max number of confs, if it equals 0, rotatable-bond-dependent window method.
export OMEGA_MAX_CONFS=${OMEGA_MAX_CONFS-600}
# set the omega torsion library: 1) Original; 2) GubaV21
export OMEGA_TORLIB=${OMEGA_TORLIB-Original}
# set the omega force field. Options are in the link below
# https://docs.eyesopen.com/toolkits/cpp/oefftk/OEFFConstants/OEMMFFSheffieldFFType.html#OEFF::OEMMFFSheffieldFFType::MMFF94Smod
export OMEGA_FF=${OMEGA_FF-MMFF94Smod}
# set the omega rmsd for clustering and filtering conformations, if it equals 0, rotatable-bond-dependent window method.
export OMEGA_RMSD=${OMEGA_RMSD-0.5}

# Dependencies

# CORINA env.sh on wynton has an incorrect path specified
export PATH="$COMMON_DIR/$CORINA_VERSION:${PATH}"

# so does openbabel
export OBABELBASE=$COMMON_DIR/$OPENBABEL_VERSION
OB_VER=$(echo $OPENBABEL_VERSION | cut -d'-' -f2-)
export BABEL_LIBDIR=$COMMON_DIR/$OPENBABEL_VERSION/lib/openbabel/$OB_VER
export BABEL_DATADIR=$COMMON_DIR/$OPENBABEL_VERSION/share/openbabel/$OB_VER
export PATH="${PATH}:${OBABELBASE}/bin"



# aaand jchem too. all the software
# activate the openeye license
LICENSE_HOME=${LICENSE_HOME-$SOFT_HOME}
export OE_LICENSE=$LICENSE_HOME/.oe-license.txt
export CHEMAXON_PATH=$COMMON_DIR/$JCHEM_VERSION
export CHEMAXON_LICENSE_URL=$LICENSE_HOME/.jchem-license.cxl
export PATH="$PATH:$CHEMAXON_PATH/bin"

LIMIT_JAVA="${DOCKBASE}/common/java-thread-limiter/mock-num-cpus 2"
export CXCALCEXE="`which cxcalc `"
export MOLCONVERTEXE="`which molconvert`"
export PATH="${PATH}:${DOCKBASE}/bin"

# activate python environment
source $PYTHONBASE/bin/activate

log $(which python) ::: $(which python3) :::

##### end initialize environment 

#source $cwd/env_new_lig_build.sh
# export DEBUG=TRUE
# this env variable used for debugging old versions only
if ! [ $BUILD_MOL2 = "true" ]; then
	MAIN_SCRIPT_NAME=${MAIN_SCRIPT_NAME-build_database_ligand_strain_noH_btingle.sh}
	${DOCKBASE}/ligand/generate/$MAIN_SCRIPT_NAME -H 7.4 --no-db $(basename $TARGET_FILE) &
	genpid=$!
else
	mkdir working
	cd working
	MAIN_SCRIPT_NAME=${MAIN_SCRIPT_NAME-build_ligands_from_mol2.py}
	python ${DOCKBASE}/ligand/generate/$MAIN_SCRIPT_NAME $(cat $TARGET_FILE) &
	genpid=$!
	cd ..
fi

function signal_generate_ligands {
        received_sigusr=TRUE
        kill -10 $genpid
}
trap signal_generate_ligands SIGUSR1
trap signal_generate_ligands SIGINT

while [ -z "$(kill -0 $genpid 2>&1)" ]; do
	sleep 5
done

if ! [ -z $received_sigusr ]; then
        reached_time_limit
fi

log "finished build job on $TARGET_FILE"

mv working/output.tar.gz $OUTPUT/$(basename $TARGET_FILE).tar.gz
exit 0
