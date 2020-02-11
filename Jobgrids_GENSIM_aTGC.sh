#!/bin/bash

# Run this from inside SherpaGeneration/Generator/test
# Make sure to have already set up the local CMSSW environment with cmsenv
# and also set up the proxy service with
# voms-proxy-init --voms cms --valid 168:00
# Reminder: cmsRun only works from lzlogin02.hep.wisc.edu

# Set up environment
cd ..
source sherpant.sh
cd test

export NEVENTS=50000
export JOBNAME="positive_1"

# export VERTEX=Zgg

# # for h3 in 0.0
# # do
# # for h4 in 5e-7 8e-7 11e-7 14e-7 17e-7 20e-7 23e-7 26e-7 29e-7 32e-7 35e-7
# # do
# # source Run_GENSIM_aTGC.sh
# # done
# # done

# export VERTEX=ZZg

# for h3 in 0.0
# do
# for h4 in 5e-7 8e-7 11e-7 14e-7 17e-7 20e-7 23e-7 26e-7 29e-7 32e-7 35e-7
# do
# source Run_GENSIM_aTGC.sh
# done
# done


# export VERTEX=Zgg

# for h3 in 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 10e-4 11e-4 12e-4
# do
# for h4 in 0.0
# do
# source Run_GENSIM_aTGC.sh
# done
# done

export VERTEX=ZZg

for h3 in 3e-4 4e-4 5e-4 6e-4 7e-4 8e-4 9e-4 10e-4 11e-4 12e-4
do
for h4 in 0.0
do
source Run_GENSIM_aTGC.sh
done
done

