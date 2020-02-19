# aTGC sample generation

### Part 1: Install Sherpa and aTGC model
`ssh` to `lzlogin02.hep.wisc.edu`\
Set up your voms proxy\
Run the install instructions in `sherpa_installation.sh`:\
\ Change line 3 to your own working directory\
\ Run the instructions one block at a time, to monitor progress and ensure that everything installs correctly (this doesn't always happen!)

### Part 2: GEN-SIM event generation
Copy the included *custom* version of `cmsRun.sh` to your working directory\
Change these lines in `cmsRun.sh` to specify your own directory:\
359, 360\
\ These lines are necessary to source sherpant in a way that is accessible to farmout\
Copy `Jobgrids_GENSIM_aTGC.sh` and `Run_GENSIM_aTGC.sh` to `CMSSW_9_3_15/src/SherpaGeneration/Generator/test`\
Change these lines in `Run_GENSIM_aTGC.sh` to specify your own directories:\
144, 149, 151\
\ Line 151 should point to the custom version of cmsRun.sh

Configure `Jobgrids_GENSIM_aTGC.sh` with the desired `NEVENTS` per sample, `JOBNAME`, `VERTEX` type, and grid of aTGC parameter values to run on

Be sure to be logged into lzlogin02, and have your `cmsenv` and voms proxy active\
Execute `Jobgrids_GENSIM_aTGC.sh`

By default, generated events will appear in the `nfs_scratch` directory\
\ hdfs is **inaccessible** from lzlogin02\
\ `nfs_scratch` has limited storage space - don't generate more than 100 GB of events at one time!\
\ You will need to move events from nfs_scratch to hdfs, so first you must ssh to a login machine that is connected to hdfs\
\ `move_files.sh` provides a template for how to move these files to hdfs

### Part 3: MINIAOD
Set up a `CMSSW_9_4_0_patch1` environment in your working directory\
Copy `edit_cfg_step<N>.py`, `Jobgrids_step<N>_aTGC.sh`, and `Run_step<N>_aTGC.sh` to `CMSSW_9_4_0_patch1/src`\
\ Here and below, `<N>` = 1, 2, and 3\
Copy the partial list of Neutrino_E-10_gun files (`Neutrino_E-10_gun_GSDR_files.txt`) to the same directory\
\ These are used for pileup simulation. If necessary, change these to whatever is convenient.\
Change line 8 of `edit_cfg_step1.py` to refer to your own copy of `Neutrino_E-10_gun_GSDR_files.txt`\
Change these lines of each `Run_step<N>_aTGC.sh` to specify your own directories:\
\ 13, 14

Configure each `Jobgrids_step<N>_aTGC.sh` with the desired `JOBNAME`, `VERTEX` type, and grid of aTGC parameter values to run on

`ssh` to an hdfs-accessible login machine\
Have your `cmsenv` and voms proxy active\
Execute `Jobgrids_step1_aTGC.sh`. When all of its jobs are fully complete, follow that with `step2`, and finally `step3`.
