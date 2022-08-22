#!/bin/bash

# A faster, smarter version of prot_marvin.py

set -e

INPUTS="-"
COLS=2
FORMAT="smiles:a"
EMBEDDING_ARGS="-T name"
EXPORTING_ARGS=""
OUTPUT=""
DELIMITER=" "

INPUTS=( "$1" )
OUTPUT="$2"
MOLCONVERT_ARGS=( -g "smiles:a" $EMBEDDING_ARGS $EXPORTING_ARGS ) 

${MOLCONVERTEXE} "${MOLCONVERT_ARGS[@]}" "${INPUTS[@]}" \
	| tail -n +2 > "${OUTPUT}"
