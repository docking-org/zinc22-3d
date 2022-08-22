#!/bin/sh

#THIS IS AN EXAMPLE STUB!!!!!!!!
# EXAMPLE USAGE
# mkdir /tmp/3302344
# cd /tmp/330234
# wget "http://zinc11.docking.org/fget.pl?l=0&z=21995818&f=m"; mv "fget.pl?l=0&z=21995818&f=m" 3302344.mol2 
# $DOCKBASE/ligand/generate/build_ligand_tautprot_smiles.sh 3302344.mol2 
# Modified by Trent Balius 

#MOL=`basename $1 .gz`
#MOL=`basename $MOL .mol2`
MOL=`basename $1 .smi`

if [[ $(wc -l < $1) > 1 ]]; then
    echo "This file has more than 1 line."
    exit
fi

#exit

CORENAEXE="$( which corina )"

echo "$CORENAEXE -i t=smiles -o t=mol2 -d rc,flapn,de=10,mc=10,wh ${MOL}.smi ${MOL}.mol2"
$CORENAEXE -i t=smiles -o t=mol2 -d rc,flapn,de=10,mc=10,wh ${MOL}.smi ${MOL}.mol2

#exit
# break up the mol2 
 # replace this with a internal script
 TEB_scripts=/nfs/home/tbalius/zzz.github/teb_scripts_programs/zzz.scripts
 python ${TEB_scripts}/separate_mol2_more10000.py ${MOL}.mol2 ${MOL}_break
# loop over all mol2

countmol2=1

mountdir="$( pwd )"

for mol2file in $( ls ${MOL}_break*.mol2 ); do

  echo  "${mountdir}/${countmol2} :: ${mol2file} "   
  mkdir -p "${mountdir}/${countmol2}"
  cd ${mountdir}/${countmol2}/
  cp ${mountdir}/${mol2file} .

  #$SOFT/openbabel/current/bin/obabel -imol2 -osmi ${MOL}.mol2 -O ${MOL}.smi
  head -1 ${mountdir}/${MOL}.smi | awk '{if(NF==1){print "None 0 "$1}; if(NF==2){print $2" 0 "$1}; if(NF>=3){print $2" "$3" "$1}}' > name.txt
  
  $DOCKBASE/ligand/amsol/calc_solvation.csh ${mol2file}
  $DOCKBASE/ligand/omega/omega_db2.e.py output.mol2

  ls output.*.db2in.mol2
  # loop over rings   
  count=1
  for file in $( ls output.*.db2in.mol2 ); do
    $DOCKBASE/ligand/mol2db2/mol2db2.py --mol2=$file --solv=output.solv --db=${MOL}.${count}.db2.gz
    count=$((count+1))
  done # file

  countmol2=$((countmol2+1))
done # mol2 file
