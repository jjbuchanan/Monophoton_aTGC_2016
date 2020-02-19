#!/bin/bash

export WORKINGDIR="/afs/hep.wisc.edu/cms/jjbuchanan/aTGC_generation"

cd $WORKINGDIR
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

# Set up environment for NTGC_UFO model installation
cd $TOPDIR
cmsenv
cd SherpaGeneration/Generator
source sherpant.sh
cd data/models

# Copy NTGC_UFO model from Bhawna's directory
cp -r /cms/gomber/Sherpa_13TeV/CMSSW_9_3_15/src/SherpaGeneration/Generator/data/models/NTGC_UFO .
# Install model
Sherpa-generate-model NTGC_UFO
