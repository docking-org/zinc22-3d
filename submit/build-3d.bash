#!/bin/bash

# req: JOB_ID // SLURM_ARRAY_JOB_ID
# req: SGE_TASK_ID // SLURM_ARRAY_TASK_ID
# req: INPUT
# req: OUTPUT
# req: LOGGING
# req: SHRTCACHE
# req: LONGCACHE

# opt: WORK_DIR

cwd=$PWD
export SOFT_HOME=${SOFT_HOME-$HOME}
export BUILD_DIR=3dbuild_$(whoami)
export SHRTCACHE=${SHRTCACHE-/scratch}
export LONGCACHE=${LONGCACHE-/scratch}
export WORK_DIR=$SHRTCACHE/$BUILD_DIR
export DOCK_VERSION=${DOCK_VERSION-DOCK.3.8.4.1.3d}
export PYENV_VERSION=${PYENV_VERSION-lig_build_py3-3.7.1}
export COMMON_DIR=$LONGCACHE/build_3d_common
export PYTHONBASE=$COMMON_DIR/$PYENV_VERSION
export DOCKBASE=$COMMON_DIR/${DOCK_VERSION}
export CORINA_VERSION=${CORINA_VERSION-corina-2025}

if ! [ -d $COMMON_DIR ]; then
        mkdir -p $COMMON_DIR
        chmod 777 $COMMON_DIR
fi

# wipe out legacy directories if they're here
#if   [ -d $TEMPDIR/build_3d_$(whoami) ]; then
#        rm -r $TEMPDIR/build_3d_$(whoami)
#elif [ -d $TEMPDIR/build_3d ]; then
#        rm -r $TEMPDIR/build_3d
#fi

JOB_ID=${SLURM_ARRAY_JOB_ID-$JOB_ID}
TASK_ID=${SLURM_ARRAY_TASK_ID-$SGE_TASK_ID}

mkdir -p $WORK_DIR

# deprecated (and jank as hell)
#function synchronize_all_but_first {
#        if [ -f /tmp/${1}.done ]; then 
#		if [ $(( (`date +%s` - `stat -L --format %Y /tmp/${1}.done`) > (10) )) ]; then
#			rm /tmp/${1}.done
#		else
#			return;
#		fi
#	fi # in the case of a particularly short running command, it might be done by the time another job even enters this function
#        flock -n /tmp/${1}.lock -c "printf ${1}; ${@:2} || echo command failed; echo > /tmp/${1}.done" && FIRST=TRUE
#        if [ -z $FIRST ]; then
#                printf "waiting ${1}"
#                n=0
#                while ! [ -f /tmp/${1}.done ]; do sleep 0.1; n=$((n+1)); if [ $n -eq 10 ]; then printf "."; n=0; fi; done
#        else
#                sleep 1 && rm /tmp/${1}.done
#        fi
#	echo
#}

# any jobs that were cancelled previously should be cleaned up
old_work=$(find $WORK_DIR -mindepth 1 -maxdepth 1 -user $(whoami) -mmin +180 -name '*.build-3d.d' | wc -l)
old_logs=$(find $SHRTCACHE -mindepth 1 -maxdepth 1 -user $(whoami) -mmin +180 -name '*batch_3d*.???' | wc -l)

# to properly synchronize work done by multiple threads we need to do a couple things
# 1. we need a long-term flag to mark the work as completed/not completed
#    a. In the case of removing old files, this "flag" might be as simple as there no longer being old files present
#    b. In the case of an install/extract operation, this would be a ".done" file in the install directory
# 2. local synchronization
#    a. It doesn't serve us well to have multiple threads trying to extract/remove the same file at the same time
#       To prevent this it helps to use a flock lock to ensure that only one process can work on the task at a time
#    b. In the case of an install operation, we might not want other threads to continue until the first thread is done
#       with it's work. In this case, we need to have a temporary "done" flag for the first thread to signal the
#       other waiting threads that it's okay to move on, which is what I do in the synchronize_all_but_first function

function extract_cmd {
	tar -C $COMMON_DIR -xzf $SOFT_HOME/soft/$1.tar.gz && echo > $COMMON_DIR/$1/.done
}

# added an additional check to make sure software dir isn't empty, since this seems to have happened before
for software in $DOCK_VERSION $PYENV_VERSION lib openbabel-install jchem-19.15 $CORINA_VERSION; do
	(
		flock -x 9
		if ! [ -f $COMMON_DIR/$software/.done ] || [ $(ls $COMMON_DIR/$software | wc -l) -eq 0 ]; then
			extract_cmd $software
		fi
		flock -u 9
	)9>/dev/shm/3d_install_${software}_${USER}.lock

done

if [ $old_work -gt 0 ] || [ $old_logs -gt 0 ]; then
	(
		flock -x 9
		find $SHRTCACHE -mindepth 1 -maxdepth 1 -user $(whoami) -mmin +180 -name '*batch_3d*.???' | xargs rm -r
		find $WORK_DIR -mindepth 1 -maxdepth 1 -user $(whoami) -mmin +120 -name '*.build-3d.d' | xargs rm -r
		flock -u 9
	)9>/dev/shm/rm_old_${USER}.lock
fi

function log {
    echo "[build-3d $(date +%X)]: $@"
}

log $(hostname)

WORK_BASE=$WORK_DIR/${JOB_ID}_${TASK_ID}.build-3d.d
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

log $(hostname)
log "starting build job on $TARGET_FILE"
log "len($TARGET_FILE)=$(cat $TARGET_FILE | wc -l)"

# move logs to their final destination & clean up the working directory
cleanup() {
        #mv $ERROR_LOG $LOGGING/$(basename $TARGET_FILE).err
        #mv $OUTPUT_LOG $LOGGING/$(basename $TARGET_FILE).out
        if [ -z $SKIP_DELETE ] && [ -d $WORK_BASE ]; then rm -r $WORK_BASE; fi
        exit $1
}

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

export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:$COMMON_DIR/lib

#export AMSOLEXE=$SOFTBASE/amsol/in-house/amsol7.1-colinear-fix/amsol7.1

# Experimental changes to DOCK ligand pipeline
export EMBED_PROTOMERS_3D_EXE=$DOCKBASE/ligand/3D/embed3d_corina.sh
# parameters related to omega
# set omega energy window, if it equals 0, rotatable-bond-dependent window method.
export OMEGA_ENERGY_WINDOW=12
# set omega max number of confs, if it equals 0, rotatable-bond-dependent window method.
export OMEGA_MAX_CONFS=600
# set the omega torsion library: 1) Original; 2) GubaV21
export OMEGA_TORLIB=Original
# set the omega force field. Options are in the link below
# https://docs.eyesopen.com/toolkits/cpp/oefftk/OEFFConstants/OEMMFFSheffieldFFType.html#OEFF::OEMMFFSheffieldFFType::MMFF94Smod
export OMEGA_FF=MMFF94Smod
# set the omega rmsd for clustering and filtering conformations, if it equals 0, rotatable-bond-dependent window method.
export OMEGA_RMSD=0.5

# Dependencies

# CORINA env.sh on wynton has an incorrect path specified
export PATH="$COMMON_DIR/$CORINA_VERSION:${PATH}"

# so does openbabel
export OBABELBASE=$COMMON_DIR/openbabel-install
export BABEL_LIBDIR=$COMMON_DIR/lib/openbabel/2.3.1
export BABEL_DATADIR=$COMMON_DIR/openbabel-install/share/openbabel/2.3.1
export PATH="${PATH}:${OBABELBASE}/bin"

# aaand jchem too. all the software
# activate the openeye license
export OE_LICENSE=$SOFT_HOME/.oe-license.txt
export CHEMAXON_PATH=$COMMON_DIR/jchem-19.15
export CHEMAXON_LICENSE_URL=$SOFT_HOME/.jchem-license.cxl
export PATH="$PATH:$CHEMAXON_PATH/bin"

LIMIT_JAVA="${DOCKBASE}/common/java-thread-limiter/mock-num-cpus 2"
export CXCALCEXE="`which cxcalc `"
export MOLCONVERTEXE="`which molconvert`"
export PATH="${PATH}:${DOCKBASE}/bin"

# activate python environment
source $PYTHONBASE/bin/activate

##### end initialize environment 

#source $cwd/env_new_lig_build.sh
# export DEBUG=TRUE
# this env variable used for debugging old versions only
MAIN_SCRIPT_NAME=${MAIN_SCRIPT_NAME-build_database_ligand_strain_noH_btingle.sh}
${DOCKBASE}/ligand/generate/$MAIN_SCRIPT_NAME -H 7.4 --no-db $(basename $TARGET_FILE) &
genpid=$!

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
cleanup 0
