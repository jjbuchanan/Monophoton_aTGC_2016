#!/bin/bash

cd /afs/hep.wisc.edu/cms/jjbuchanan/aTGC_generation

scram project CMSSW_9_3_15
cd CMSSW_9_3_15/src
cmsenv   
export TOPDIR=$PWD

git cms-addpkg -q GeneratorInterface/SherpaInterface

git clone https://github.com/SiewYan/SherpaGeneration.git -b EXO-2.2.5
cp $TOPDIR/GeneratorInterface/SherpaInterface/data/*SherpaLibs.sh $TOPDIR/SherpaGeneration/Generator/test/
scram b -j8

cd $TOPDIR/SherpaGeneration/Generator/
mkdir sherpant
./fetchSherpa.sh
mv buildSherpant.sh SHERPA-MC-2.2.5
cd $TOPDIR/SherpaGeneration/Generator/SHERPA-MC-2.2.5
./buildSherpant.sh

make install -j4

cd $TOPDIR/SherpaGeneration/Generator/
source sherpant.sh
cd $TOPDIR
scram b -j4
