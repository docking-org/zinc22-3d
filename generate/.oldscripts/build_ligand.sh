#!/bin/sh

#THIS IS AN EXAMPLE STUB!!!!!!!!
# EXAMPLE USAGE
# mkdir /tmp/3302344
# cd /tmp/330234
# $DOCKBASE/ligand/generate/build_ligand.sh /db/zinc/sources/3302344.mol2 \
#                                           --name=ZINC003302344 \
#                                           --smiles='Cc1cccc(c1)c2nnc(n2c3ccc(cc3)Cl)S[C@H](C)C(=O)Nc4ccc(cc4)NC(=O)C'
# cp -v 3302344.db2.gz /db/zinc/dockable


MOL=`basename $1 .gz`
MOL=`basename $MOL .mol2`

$DOCKBASE/ligand/generate/prepare.py $*
$DOCKBASE/ligand/amsol/calc_solvation.csh $1
$DOCKBASE/ligand/mol2db2/mol2db2.py --mol2=$1 --solv=output.solv --db=$MOL.db2.gz
