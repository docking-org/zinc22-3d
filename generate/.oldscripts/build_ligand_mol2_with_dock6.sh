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

$OBABELBASE/bin/obabel -h -imol2 -osmi ${MOL}.mol2 -O ${MOL}_OB.smi

#cp MRTX849_0.smi ${MOL}.smi 

head -1 ${MOL}_OB.smi | awk '{if(NF==1){print "None 0 "$1}; if(NF==2){print $2" 0 "$1}; if(NF>=3){print $2" "$3" "$1}}' > name.txt

python $DOCKBASE/ligand/generate/mol2tosmi_rdkit.py ${MOL}.mol2 ${MOL}.smi 

#$DOCKBASE/ligand/generate/prepare.py $*
$DOCKBASE/ligand/amsol/calc_solvation.csh $1

#$SOFT/openbabel/current/bin/obabel -imol2 -osmi output.mol2 -O output.smi

csh $DOCKBASE/ligand/generate/dock6_confgen_db2.csh output.mol2 output.solv name.txt 

