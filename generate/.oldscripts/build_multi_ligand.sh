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

echo "name.txt $MOL $MOL $MOL $MOL $MOL" > name.txt
$DOCKBASE/ligand/generate/prepare.py -f name.txt $*
$DOCKBASE/ligand/amsol/calc_solvation.csh $1
$DOCKBASE/ligand/omega/omega_db2.e.py $1
for onemol in `\ls -1 *.db2in.mol2`
do
  $DOCKBASE/ligand/mol2db2/mol2db2.py -v --mol2=$onemol --solv=output.solv --db=$onemol.db2.gz 
done
zcat *.db2.gz > final.db2
gzip -9 final.db2
