#!/bin/bash

PROCS=$( grep processor /proc/cpuinfo | wc -l )
PICKED=$(( $RANDOM % $PROCS ))

echo "Notice: Limiting $@ to one core" 1>&2
numactl -C $PICKED "$@"
