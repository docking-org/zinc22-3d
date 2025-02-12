#!/bin/sh

if [[ $# -eq 0 ]]; then
    MOLDIR='.'
else
    MOLDIR=$1
fi

rm -v $MOLDIR/output* \
      $MOLDIR/fort.* \
      $MOLDIR/temp.* \
      $MOLDIR/*.db2.gz \
      $MOLDIR/name.txt
