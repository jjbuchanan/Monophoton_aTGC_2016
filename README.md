# Monophoton_aTGC_2016

### Step 1: Install Sherpa and aTGC model
ssh to lzlogin02.hep.wisc.edu\
Set up your voms proxy\
Run the install instructions in sherpa_installation.sh:\
\ Change line 3 to your own working directory\
\ Run the instructions one block at a time, to monitor progress and ensure that everything installs correctly (this doesn't always happen!)

### Step 2: GEN-SIM event generation
Copy the included *custom* version of cmsRun.sh to your working directory
Change these lines in cmsRun.sh to specify your own directory:\
359, 360\
\ These lines are necessary to source sherpant in a way that is accessible to farmout
Copy Jobgrids_GENSIM_aTGC.sh and Run_GENSIM_aTGC.sh to CMSSW_9_3_15/src/SherpaGeneration/Generator/test\
Change these lines in Run_GENSIM_aTGC.sh to specify your own directories:\
144, 149, 151\
\ Line 151 should point to the custom version of cmsRun.sh

Configure Jobgrids_GENSIM_aTGC.sh with the desired NEVENTS per sample, JOBNAME, VERTEX type, and grid of aTGC parameter values to run on\
Be sure to be logged into lzlogin02, and have your cmsenv and voms proxy active\
Execute Jobgrids_GENSIM_aTGC.sh

By default, generated events will appear in the nfs_scratch directory\
\ hdfs is *inaccessible* from lzlogin02\
\ nfs_scratch has limited storage space - don't generate more than 100 GB of events at one time!\
\ You will need to move events from nfs_scratch to hdfs, so first you must ssh to a login machine that is connected to hdfs\
\ move_files.sh provides a template for how to move these files to hdfs

