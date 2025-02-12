#!/bin/sh

#THIS IS AN EXAMPLE STUB!!!!!!!!
# EXAMPLE USAGE
# mkdir /tmp/3302344
# cd /tmp/330234
# $DOCKBASE/ligand/generate/build_ligand.sh /db/zinc/sources/3302344.mol2 \
#                                           --name=ZINC003302344 \
#                                           --smiles='Cc1cccc(c1)c2nnc(n2c3ccc(cc3)Cl)S[C@H](C)C(=O)Nc4ccc(cc4)NC(=O)C'
# cp -v 3302344.db2.gz /db/zinc/dockable
#
set -e

PH=7.4
DOMAIN_CUTOFF=20

if [ ! -z $( which cxcalc ) ]; then
    export CXCALCEXE=$( which cxcalc )
fi
if [ -z "$CXCALCEXE" ]; then
    echo "This proof-of-concept script requires CXCALCEXE (cxcalc) to be set! Exiting" 1>&2
    exit -1
fi

if [ ! -z $( which molconvert ) ]; then
    export MOLCONVERTEXE=$( which molconvert )
fi
if [ -z "$MOLCONVERTEXE" ]; then
    echo "This proof-of-concept script requires MOLCONVERTEXE (molconvert) to be set! Exiting" 1>&2
    exit -1
fi


DBNAME=$1
DBNAME=`basename $DBNAME .gz`
DBNAME=`basename $DBNAME .smi`
DBNAME=`basename $DBNAME .ism`

COVALENT=$2

echo $DBNAME

while read LINE
do
    if [ -z "$LINE" ]; then
        continue
    fi
    TOKENS=(${LINE//\s\+/})
    SMILES=${TOKENS[0]}
    NAME=${TOKENS[1]}
    if [ -z "$NAME" ] || [ -z "$SMILES" ]; then
        echo "Skipping, missing information $NAME $SMILES"
        continue
    fi
    echo "Building $NAME" 1>&2

    mkdir -pv $NAME
    pushd $NAME 1>&2

    echo ""
    echo "Generating protomers using Marvin (ChemAxon):" 1>&2
    echo "$LINE" | $DOCKBASE/ligand/protonate/prot_marvin.py -H $PH -d $DOMAIN_CUTOFF -o $NAME-protonated.ism
    echo $( wc -l $NAME-protonated.ism ) "protomers generated" 1>&2
    echo ""

    PROT_NUM=0
    while read PROTOMER
    do
        PROT_TOKENS=(${PROTOMER//\s\+/})
        PROT_SMILES=${TOKENS[0]}
        PROT_ID=${TOKENS[1]}
        
        PROT_NAME="${PROT_ID}.${PROT_NUM}"
        mkdir -pv $PROT_NUM
        pushd $PROT_NUM 1>&2

        echo ""
        echo "Embedding SMILES in 3D and adding Hydrogens for $PROT_NAME" 1>&2
        echo "${PROTOMER}" | ${MOLCONVERTEXE} mol2 "-3:[prehydrogenize][S]{fine}" > $PROT_NAME.mol2
        echo "Preparing input files" 1>&2
        $DOCKBASE/ligand/generate/prepare.py -n $PROT_NAME -s $PROT_SMILES $PROT_NAME.mol2
        echo ""
        echo "Starting the preparation of the solvation calculations (AMSOL7.1)" 1>&2
        echo "  (SMILES: $PROT_SMILES)" 1>&2
        echo ""
        $DOCKBASE/ligand/amsol/calc_solvation.csh $PROT_NAME.mol2

	#"if covalent" statement here
	if [ "${2}" == "--covalent" ]; then 
          $DOCKBASE/ligand/omega/omega_warhead.py output.mol2
	else
          $DOCKBASE/ligand/omega/omega_db2.e.py output.mol2
	fi

        for onemol in `\ls -1 *.db2in.mol2`
        do
             $DOCKBASE/ligand/mol2db2/mol2db2.py -v $2 --mol2=$onemol --solv=output.solv --db=$onemol.db2.gz 
        done

        # Build DB files from extracted conformations
        ln -svfn output.solv mol.solv 1>&2
        if ! ( for onemol in $( find . -maxdepth 1 -name '*.db2in.mol2' ) ; do
             ln -sv $onemol db.mol2 1>&2
             $DOCKBASE/ligand/mol2db/mol2db $DOCKBASE/ligand/mol2db/data/inhier 1>&2
             mv -v db.db $onemol.db 1>&2
             rm -v db.mol2 1>&2
             gzip $onemol.db 1>&2
          done ); then
            echo "DB File generation failed! Skipping $NAME $PROT_NAME" 1>&2
            popd 1>&2
            mv -v $PROT_NUM $PROT_NUM.failed 1>&2
            PROT_NUM=`expr $PROT_NUM + 1`
            break 
        fi
        popd 1>&2
        PROT_NUM=`expr $PROT_NUM + 1`
    done < $NAME-protonated.ism

    popd 1>&2 
    echo ""
    echo "Writing to $DBNAME swap file" 1>&2
    zcat $NAME/**/*.db2.gz >> $DBNAME-swap.db2
    echo ""
    echo "Finished preparing $NAME"
done < $1
echo "Finalizing..."
mv -v $DBNAME-swap.db2 $DBNAME.db2
gzip -9 $DBNAME.db2

