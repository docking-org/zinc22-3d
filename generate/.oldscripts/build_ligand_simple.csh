#!/bin/csh
# Written by Trent Balius at FNLCR on Feb 6, 2020.  
# This script is to smiplify the Database generation or at less make if easer for me to modify it.   
# This only generated a mol2 file with one conformation for each molecules in smiles.  
#

  if (-e db_build_working) then
      echo "Remove this directory: db_build_working"
      exit
  endif

  mkdir db_build_working
  cd db_build_working
  cp ../$1 .


  ## Step 1. process smiles with chemAxon to protonate and tautomerize.  
  source ~baliuste/zzz.programs/jchem/env.csh
  #echo "I AM HERE"
  echo "filename = $1"
#  set smilist = `cat $1 | awk '{print "\""$1"\""}'`   
#  echo "$smilist"
  set PH = 7.2
  set TAUTOMER_LIMIT = 30
  set PROTOMER_LIMIT = 30
  set TAUT_PROT_CUTOFF = 1
  set START = 1 

  set CXCALCEXE = `which cxcalc`  
  set MOLCONVERTEXE = `which molconvert` 

  sed 's/\s\+/ /g' "${1}" | \
        ${CXCALCEXE} -g dominanttautomerdistribution -H "${PH}" -C false -t tautomer-dist | \
        ${MOLCONVERTEXE} sdf -g -c "tautomer-dist>=${TAUTOMER_LIMIT}" | \
        ${CXCALCEXE} -g microspeciesdistribution -H $PH -t protomer-dist | \
        ${MOLCONVERTEXE} smiles -g -c "protomer-dist>=${PROTOMER_LIMIT}" -T name:tautomer-dist:protomer-dist | \
        awk -v "cutoff=${TAUT_PROT_CUTOFF}" -v "start=${START}" '{ if (NR == 1 && start < 2) { print $0, "score" } else { score = ($3 * $4)/100 ; if (score >= cutoff) { print $0, score } } }' > prot-taut.info

  
  # step 2. convert smiles to mol2 files.  
  source /home/baliuste/zzz.programs/corina/env.csh
  

  split -a 8 -l 1 prot-taut.info  prot-taut_split_

  mv prot-taut_split_aaaaaaaa header # move the header
  set count = 0
  #foreach file (`ls prot-taut_split_???????[bcdghijklmnopqrstuvwxyz]`)

  set mountdir = `pwd`
  touch dirlist
  foreach file (`ls prot-taut_split_????????`)
     cd $mountdir
     set name = `awk -F'\t' '{print $2}' $file`
     set newname =  "${name}_$count"
     set workdir = ${mountdir}/${name}
     if !(-e $workdir) then 
        mkdir $workdir
        echo $name >> dirlist # remember all of the dir that we make
     endif
     
     echo "$newname"
     if (-e "$workdir/$newname.smi") then
        "Error. "
        exit
     endif
     mv "$file" "$workdir/$newname.smi"
     cd $workdir
     /home/baliuste/zzz.programs/corina/corina -i t=smiles -o t=mol2 -d rc,flapn,de=6,mc=1,wh $newname.smi $newname.mol2
      
     @ count = $count + 1
  end
  #awk '{print $1 " " $2}' prot-taut.info > prot-taut.smi
  
  #/home/baliuste/zzz.programs/corina/corina -i t=smiles -o t=mol2 -d rc,flapn,de=6,mc=1,wh prot-taut.smi prot-taut.mol2

  # Step 3. Run Amsol. 
  source ~/.cshrc.python3
  
  foreach  dir (`cat ${mountdir}/dirlist`)
     cd ${mountdir}/${dir}/
     foreach mol2 (`ls *.mol2`)
        echo $mol2
        ~/zzz.github/DOCK/ligand/amsol/calc_solvation.csh $mol2 
        mv output.mol2 ${mol2:r}_output.mol2
        mv output.solv ${mol2:r}_output.solv
     end
  end

  # Step 4. Perform conformational expantion and db2 generation. 
  #
 


  # Step 5. 
