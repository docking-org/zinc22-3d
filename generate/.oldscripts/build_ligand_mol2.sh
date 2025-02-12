#!/bin/sh

#THIS IS AN EXAMPLE STUB!!!!!!!!
# EXAMPLE USAGE
# mkdir /tmp/3302344
# cd /tmp/330234
# wget "http://zinc11.docking.org/fget.pl?l=0&z=21995818&f=m"; mv "fget.pl?l=0&z=21995818&f=m" 3302344.mol2 
# $DOCKBASE/ligand/generate/build_ligand_mol2.sh 3302344.mol2 
# Modified by Trent Balius 


MOL=`basename $1 .gz`
MOL=`basename $MOL .mol2`

$SOFT/openbabel/current/bin/obabel -imol2 -osmi ${MOL}.mol2 -O ${MOL}.smi
head -1 ${MOL}.smi | awk '{if(NF==1){print "None 0 "$1}; if(NF==2){print $2" 0 "$1}; if(NF>=3){print $2" "$3" "$1}}' > name.txt

#$DOCKBASE/ligand/generate/prepare.py $*
$DOCKBASE/ligand/amsol/calc_solvation.csh $1
$DOCKBASE/ligand/omega/omega_db2.e.py output.mol2
#$DOCKBASE/ligand/mol2db2/mol2db2.py --mol2=output.1.db2in.mol2 --solv=output.solv --db=${MOL}.db2.gz

count=1
for file in $( ls output.*.db2in.mol2 ); do
  $DOCKBASE/ligand/mol2db2/mol2db2.py --mol2=$file --solv=output.solv --db=${MOL}.${count}.db2.gz
  count=$((count+1))
done

