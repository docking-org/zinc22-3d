#!/bin/bash

# A faster, smarter version of prot_marvin.py

set -e

if [ -z "$MOLCONVERTEXE" ]; then
	MOLCONVERTEXE=$( which molconvert )
fi


EXPORTING_ARGS="-m"
EMBEDDING_ARGS=""  # This changes the conf?
OUTPUT=""
DELIMITER=" "

INPUTS=( "$1" )
OUTPUT="$2"
MOLCONVERT_ARGS=( -g "mol2:+H" $EMBEDDING_ARGS $EXPORTING_ARGS -o "${OUTPUT}" ) 

${MOLCONVERTEXE} "${MOLCONVERT_ARGS[@]}" "${INPUTS[@]}" 
