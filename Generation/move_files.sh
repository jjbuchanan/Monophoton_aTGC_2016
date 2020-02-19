#!/bin/bash

# ### ### ### ### ### ### ### ### ### #
# Move files from nfs_scratch to hdfs #
# ### ### ### ### ### ### ### ### ### #

VERTEX="ZZg"
JOBNAME="positive_1"

for h3 in 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 10e-4 11e-4 12e-4
do
for h4 in 0.0
do
  echo "Moving h3=${h3}, h4=${h4}"
  
  h33=$(echo ${h3} | sed -e 's/\./p/' -e 's/-/m/')
  h44=$(echo ${h4} | sed -e 's/\./p/' -e 's/-/m/')
  mkdir /hdfs/store/user/jjbuch/aTGC/${VERTEX}/h3${h33}_h4${h44}_${JOBNAME}
  
  for i in {0..49}
  do
    NUM=$i
    if [ $i -lt 10 ]
    then
      NUM=0$NUM
    fi
    echo "NUM=${NUM}"
    
    DIR1=/nfs_scratch/jjbuchanan/gensim_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}_${JOBNAME}/gensim_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}_${JOBNAME}-00${NUM}
    DIR2=/hdfs/store/user/jjbuch/aTGC/${VERTEX}/h3${h33}_h4${h44}_${JOBNAME}
    FNAME=sherpa_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}_MASTER_cff_py_GEN_SIM
    mv ${DIR1}/${FNAME}.root ${DIR2}/${FNAME}_00${NUM}.root
  done
done
done
