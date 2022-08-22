#!/bin/csh

# $1 is the smiles file with all compounds to be generated
# $2 is the how many subdirs do you want to end up with 
# $3 is the prefix for the dirnames
# usage: build_covalent_lib.csh aldehydes.ism 500 aldehyde

set chunk_size=$2
set prefix=$3
set suffix=`echo $1 | awk -F. '{print $NF}'`
set nmols=`cat $1 | wc -l`
set ndirs=`date | awk '{print int('$nmols'/'$chunk_size')+1}'`

foreach i (`seq 1 $chunk_size`)
  echo "Now preping chunck $i"
  set start=`date | awk '{print (('$i'-1)*'$ndirs')+1}'`
  set end=`date | awk '{print (('$i'-1)*'$ndirs')+'$ndirs'}'`
  #echo $start $end
  sed -n ''$start','$end''p $1 > $prefix.$i.$suffix
  mkdir $prefix.$i
  set ncompsdir = `wc -l $prefix.$i.$suffix | awk '{print $1}'` 
  mv $prefix.$i.$suffix $prefix.$i
  if ($ncompsdir == 0) then
    rm -rf $prefix.$i
  else 
    cd $prefix.$i
    set rundir=`pwd`
    mkdir logs
    mkdir input
    mkdir output
    foreach j (`seq 1 $ncompsdir`)
      cat $prefix.$i.$suffix | awk 'NR=='$j'' > input/$prefix.$i.$j.$suffix
    end
    cd logs
    qsub -t 1-$ncompsdir $DOCKBASE/ligand/generate/cov-lib-submit.sh $rundir
    cd ../../
  endif
end
