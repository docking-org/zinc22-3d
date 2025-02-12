#!/bin/bash

# A faster, smarter version of prot_marvin.py

set -e

END_OF_RECORD="$( echo -e "#\tEnd of record" )"

INPUTS=()
COLS=2
FORMAT="mol2"
EMBEDDING_ARGS="-d rc,flapn,de=6,mc=${CORINA_MAX_CONFS-1},wh"
SPLIT_ON_DELIMITER="${END_OF_RECORD}"
OUTPUT=""
DELIMITER=" "
DEBUG=""

echo "We are using corina for 3D embeding"

if [ -z "$CORENAEXE" ] ; then
	CORENAEXE="$( which corina )"
fi

if [ -z "${SPLITONEXE}" ] ; then
    SPLITONEXE="$( dirname "${0}" )/spliton.py"
fi

while [[ "$#" -gt 0 ]]; do
	ARG="$1"
	case $ARG in
	    -d|--delimiter)
	        DELIMITER="${2}"
	        shift 2
	        ;;
	    -c|--columns)
	        COLS="${2}"
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
        -s|--single-file)
			SPLIT_ON_DELIMITER=""
			shift 1
			;;
		-m|--multiple-files)
			SPLIT_ON_DELIMITER="${END_OF_RECORD}"
			shift 1
			;;
		-o|--output)
			OUTPUT="${2}"
			shift 2
			;;
		--debug)
			DEBUG="yes"
			shift 1
			;;
		*)
			INPUTS+="${ARG}"
			shift 1
			;;
	esac
done

INPUTS=( "${INPUTS[@]}" )
CORINA_ARGS=( -i t=smiles -o t=${FORMAT} )
if [ ! -z "${EMBEDDING_ARGS}" ]; then
	CORINA_ARGS=( ${CORINA_ARGS[@]} ${EMBEDDING_ARGS[@]} )
fi

echo "debuging info::" ${INPUTS[@]} ${CORENAEXE} ${CORINA_ARGS[@]}

if [ -z "${SPLIT_ON_DELIMITER}" ] ; then
    cat "${INPUTS[@]}" | \
 	cut -d "${DELIMITER}" -f "1-${COLS}" | \
 	${CORENAEXE} "${CORINA_ARGS[@]}" | \
	grep -v '^#' \
 	    > "${OUTPUT}"
else
    cat "${INPUTS[@]}" | \
 	cut -d "${DELIMITER}" -f "1-${COLS}" | \
 	${CORENAEXE} "${CORINA_ARGS[@]}" | \
 	    ${SPLITONEXE} "${SPLIT_ON_DELIMITER}" "${OUTPUT}"
	for F in "${OUTPUT}"* ; do
		sed -i -e 's/^#.*$//' -e '/^$/d' "${F}"
	done
fi

if [ "${DEBUG}" != "yes" ] ; then
	rm -v corina.trc
fi
