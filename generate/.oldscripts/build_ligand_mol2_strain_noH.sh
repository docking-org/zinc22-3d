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
#MOL=`basename $MOL .db2`

/nfs/soft/openbabel/current/bin/obabel -imol2 -osmi $1 -O ${MOL}.smi
#/nfs/soft/openbabel/current/bin/obabel -imol2 -osmi ${MOL}.mol2 -O ${MOL}.smi
head -1 ${MOL}.smi | awk '{if(NF==1){print "None 0 "$1}; if(NF==2){print $2" 0 "$1}; if(NF>=3){print $2" "$3" "$1}}' > name.txt

#$DOCKBASE/ligand/generate/prepare.py $*
$DOCKBASE/ligand/amsol/calc_solvation.csh $1
#$DOCKBASE/ligand/omega/omega_db2.e.py output.mol2
#python $DOCKBASE/ligand/strain/db2tomol2.py ${MOL}.db2.gz
$DOCKBASE/ligand/omega/omega_db2.e.py3.py output.mol2
#$DOCKBASE/ligand/mol2db2/mol2db2.py --mol2=output.1.db2in.mol2 --solv=output.solv --db=${MOL}.db2.gz

count=1
for file in $( ls output.*.db2in.mol2 ); do
#for file in $( ls ${MOL}.*.mol2 ); do
  #run UCSF strain
  #pwd=`pwd`
  #cd /nfs/ex5/work/jklyu/strainfilter/
  #python Torsion_Strain.py ${pwd}/${file}
  #cd -
  prefix=`basename ${file} .mol2`
  #echo $prefix
  #awk -F"," '{print $1" "$2" "$6}' ${prefix}_Torsion_Strain.csv > ${prefix}.strain.details.txt
  #add total strains into output.mol2
  #python $DOCKBASE/ligand/strain/add_strain_info.py ${file} ${prefix}.strain.details.txt
  #python $DOCKBASE/ligand/strain/add_strain_info.py ${file}
  python $DOCKBASE/ligand/strain/add_strain_noH_info.py ${file}
  $DOCKBASE/ligand/mol2db2_py3_strain/mol2db2.py --mol2=${prefix}.strain.mol2 --solv=output.solv --db=${MOL}.${count}.db2.gz
  count=$((count+1))
done

