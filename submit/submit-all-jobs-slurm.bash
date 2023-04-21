#!/bin/bash

BINDIR=$(dirname $0)
BINDIR=${BINDIR-.}
if ! [[ "$BINDIR" == /* ]]; then
	BINDIR=$PWD/$BINDIR
fi

export SUBMIT_MODE=SLURM
$BINDIR/submit-all.bash
