#!/bin/bash

# A faster, smarter version of prot_marvin.py

set -e

INPUT="-"
PH=7.4
START=2
TAUT_PROT_CUTOFF="${TAUT_PROT_CUTOFF-10}"
TAUTOMER_LIMIT="${TAUTOMER_LIMIT-20}"
PROTOMER_LIMIT="${PROTOMER_LIMIT-20}"

if [ -z "$CXCALCEXE" ]; then
	CXCALCEXE=$( which cxcalc )
fi
if [ -z "$MOLCONVERTEXE" ]; then
	MOLCONVERTEXE=$( which molconvert )
fi

while [[ "$#" -gt 0 ]]; do
	ARG="$1"
	case $ARG in 
		-H|--pH)
			PH="$2"
			shift 2
			;;
                -c|--cutoff)
			TAUT_PROT_CUTOFF="$2"
			shift 2
			;;
		-t|--tautomer-limit)
			TAUTOMER_LIMIT="$2"
			shift 2
			;;
		-p|--protomer-limit)
			PROTOMOER_LIMIT="$2"
			shift 2
			;;
		-h|--header)
			START=2
			shift 1;
			;;
		*)
			if [ "$INPUT" == '-' ]; then
				INPUT="$1"
			else
				INPUT="$INPUT $1"
			fi
			shift 1
			;;
	esac
done

sed 's/\s\+/ /g' "${INPUT}" | \
	${CXCALCEXE} -g dominanttautomerdistribution -H "${PH}" -C false -t tautomer-dist | \
	${MOLCONVERTEXE} sdf -g -c "tautomer-dist>=${TAUTOMER_LIMIT}" | \
	${CXCALCEXE} -g microspeciesdistribution -H $PH -t protomer-dist | \
	${MOLCONVERTEXE} smiles -g -c "protomer-dist>=${PROTOMER_LIMIT}" -T name:tautomer-dist:protomer-dist | \
        awk -v "cutoff=${TAUT_PROT_CUTOFF}" -v "start=${START}" '{
                if (NR == 1 && start < 2) { 
                        print $0, "score" 
                } else { 
                        score = ($3 * $4)/100 ; 
                        if (score >= cutoff) {
                                print $0, score
                        }
                }
        }'

