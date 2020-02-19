# Expected inputs:
# NEVENTS (int), JOBNAME (string), VERTEX (string), h3 (float), h4 (float)

echo
echo "*"
echo "Assembling job ${VERTEX} h3=${h3} h4=${h4}"
echo "*"
echo "*"

# Assemble datacard
if [ "$VERTEX" = "Zgg" ]; then
  export h3Zgg=${h3}
  export h4Zgg=${h4}
  export h3ZZg=0.0
  export h4ZZg=0.0
fi
if [ "$VERTEX" = "ZZg" ]; then
  export h3Zgg=0.0
  export h4Zgg=0.0
  export h3ZZg=${h3}
  export h4ZZg=${h4}
fi

cat>Run.dat<<EOF
(run){
  EVENTS = 1; ERROR 0.1;
  FSF:=1.; RSF:=1.; QSF:=1.;

  SCALES STRICT_METS{FSF*MU_F2}{RSF*MU_R2}{QSF*MU_Q2};
  ME_SIGNAL_GENERATOR Comix;
  EVENT_GENERATION_MODE PartiallyUnweighted;

  MASSIVE[5] 1
  MASSIVE[4] 1

  BEAM_1 2212; BEAM_ENERGY_1 6500.;
  BEAM_2 2212; BEAM_ENERGY_2 6500.;


  PARTICLE_CONTAINER 900[m:-1] gammaZ 22 23;

  MODEL NTGC_UFO;
}(run)


(processes){
  Process 93 93 -> 900[a] 22;
  Decay 900[a] -> 91 91;
  Order (*,*,*);
  Print_Graphs Process;
  Integration_Error 0.02;
  End process

}(processes)


(selector){
  # phase space cuts for matrix elements
  DecayMass 900 30.0 E_CMS;
  IsolationCut 22 0.4 1 0.1;
  ET  22  150 1200;
  PseudoRapidity 22 -2.6 2.6
}(selector)


(ufo){
block yukawa
  1  0.00504     # ymdo
  2  0.00255     # ymup
  3  0.101       # yms
  4  1.27        # ymc
  5  4.7         # ymb
  6  172         # ymt
  11 0.000511    # yme
  13 0.10566     # ymm
  15 1.777       # ymtau

block sminputs
  1  127.9       # aEWM1
  2  1.16637e-05 # Gf
  3  0.1184      # aS

block antgc
  1  0.0          # f4a
  2  0.0          # f4Z
  3  0.0          # f5a
  4  0.0          # f5Z
  5  0.0          # h1a
  6  0.0          # h1Z
  7  ${h3Zgg}     # h3a
  8  ${h3ZZg}     # h3Z
  9  0.0          # h2a
  10 0.0          # h2Z
  11 ${h4Zgg}     # h4a
  12 ${h4ZZg}     # h4Z

block mass
  23 91.1876     # MZ
  11 0.    # Me
  13 0.     # MMU
  15 1.777       # MTA
  2  0.0     # MU
  4  1.27        # MC
  6  172         # MT
  1  0.0     # MD
  3  0.       # MS
  5  4.7         # MB
  25 125         # MH


decay   23 2.4952      # WZ
decay   24 2.085       # WW
decay   6  1.50833649  # WT
decay   25 0.00407     # WH

}(ufo)

EOF

h33=$(echo ${h3} | sed -e 's/\./p/' -e 's/-/m/')
h44=$(echo ${h4} | sed -e 's/\./p/' -e 's/-/m/')

cp Run.dat Run.dat_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}

# Generate Sherpack from datacard
./run_MakeSherpaLibs.sh Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}

./run_PrepareSherpaLibs.sh Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}
scp sherpa_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}_MASTER_cff.py ../python 
cd ..
scram b 

# Make config file for GEN-SIM generation
cd test
cmsDriver.py SherpaGeneration/Generator/python/sherpa_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}_MASTER_cff.py --mc --eventcontent RAWSIM --datatier GEN-SIM --conditions 93X_mc2017_realistic_v3 --beamspot Realistic25ns13TeVEarly2017Collision --step GEN,SIM --geometry DB:Extended --era Run2_2017 --no_exec
SEARCH="input = cms.untracked.int32(1)"
REPLACE="input = cms.untracked.int32(\$nEventsPerJob)"
sed -i "s/${SEARCH}/${REPLACE}/" sherpa_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}_MASTER_cff_py_GEN_SIM.py
SEARCH="process.load('Configuration.StandardSequences.Services_cff')"
REPLACE="process.load('Configuration.StandardSequences.Services_cff')\nprocess.RandomNumberGeneratorService.generator.initialSeed = \$randomNumber"
sed -i "s/${SEARCH}/${REPLACE}/" sherpa_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}_MASTER_cff_py_GEN_SIM.py

# Set up farmout job, but don't submit yet
DIR=/afs/hep.wisc.edu/cms/jjbuchanan/aTGC_generation/CMSSW_9_3_15/src/SherpaGeneration/Generator
farmoutRandomSeedJobs --no-submit --use-singularity=SL6 --extra-inputs=$DIR/test/sherpa_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}_MASTER.tgz,$DIR/test/sherpa_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}_MASTER.md5 gensim_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}_${JOBNAME} ${NEVENTS} 1000 /afs/hep.wisc.edu/cms/jjbuchanan/aTGC_generation/CMSSW_9_3_15 $DIR/test/sherpa_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}_MASTER_cff_py_GEN_SIM.py 'outputFile=$outputFileName' 'maxEvents=$nEventsPerJob' 'randomSeed=$randomNumber'

# Specify cmsRun.sh location in submit file,
# and finally submit
cd /nfs_scratch/jjbuchanan/gensim_Znunugamma_aTGC_${VERTEX}_h3${h33}_h4${h44}_${JOBNAME}
SEARCH=/afs/hep.wisc.edu/cms/sw/farmout/cmsRun.sh
REPLACE=/afs/hep.wisc.edu/cms/jjbuchanan/aTGC_generation/cmsRun.sh
sed -i "s|${SEARCH}|${REPLACE}|" submit
condor_submit submit

# Return to test directory to set up next job
cd -

echo "*"
echo "*"
echo "Submitted job ${VERTEX} h3=${h3} h4=${h4}"
echo "*"
echo
