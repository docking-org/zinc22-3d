#!/bin/bash

# A faster, smarter version of prot_marvin.py

set -e

INPUTS="-"
COLS=2
FORMAT="mol2"
EMBEDDING_ARGS="-3:[prehydrogenize][S]{fine}"
EXPORTING_ARGS="-m"
OUTPUT=""
DELIMITER=" "

if [ -z "$MOLCONVERTEXE" ]; then
	MOLCONVERTEXE="$( which molconvert )"
fi

if [ -z "$OBABELEXE" ] ; then
	OBABELEXE="$( which obabel )"
fi

while [[ "$#" -gt 0 ]]; do
	ARG="$1"
	case $ARG in 
		-c|--columns)
			COLS="$2"
			shift 2
			;;
		-d|--delimiter)
			DELIMITER="${2}"
			shift 2
			;;
                -e|--embedding)
			EMBEDDING_ARGS="${2}"
			shift 2
			;;
                -f|--format)
			FORMAT="${2}"
			shift 2
			;;
                -A|--allow-stupidity)
			OBABELEXE=""
			shift 1
			;;
		-s|--single-file)
			EXPORTING_ARGS=""
			shift 1
			;;
		-m|--multiple-files)
			EXPORTING_ARGS="-m"
			shift 1
			;;
		-o|--output)
			OUTPUT="-o ${2}"
			shift 2
			;;
		*)
			if [ "${INPUTS}" == '-' ]; then
				INPUTS="${ARG}"
			else
				INPUTS="${INPUTS} ${ARG}"
			fi
			shift 1
			;;
	esac
done

INPUTS=( ${INPUTS} )
MOLCONVERT_ARGS=( -g ${FORMAT} )
if [ ! -z "${EMBEDDING_ARGS}" ]; then
	MOLCONVERT_ARGS=( ${MOLCONVERT_ARGS[@]} ${EMBEDDING_ARGS[@]} )
fi
if [ ! -z "${EXPORTING_ARGS}" ]; then
	MOLCONVERT_ARGS=( ${MOLCONVERT_ARGS[@]} ${EXPORTING_ARGS[@]} )
	if [ -z "${OUTPUT}" ] ; then
		OUTPUT='-o ./'
	fi
fi
if [ ! -z "${OUTPUT}" ] ; then
	MOLCONVERT_ARGS=( ${MOLCONVERT_ARGS[@]} ${OUTPUT[@]} )
fi

cat "${INPUTS[@]}" | \
 	cut -d "${DELIMITER}" -f "1-${COLS}" | \
 	${MOLCONVERTEXE} "${MOLCONVERT_ARGS[@]}"

if [ ! -z "$OBABELEXE" ] ; then
    echo  $( find . -name "${OUTPUT[1]}[0-9]" )
	for F in $( find . -name "${OUTPUT[1]}[0-9]*" ); do
		$OBABELEXE -i mol2 "$F" -o mol2 -O "${F}.oechem"
		[ ! -s "${F}.oechem" ] && mv -v "${F}.oechem" "$F" || rm -v "${F}.oechem"
	done
fi
