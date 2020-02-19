echo
echo "*"
echo "Assembling step1 job ${VERTEX} h3=${h3} h4=${h4}"
echo "*"
echo "*"

h33=$(echo ${h3} | sed -e 's/\./p/' -e 's/-/m/')
h44=$(echo ${h4} | sed -e 's/\./p/' -e 's/-/m/')

cmsDriver.py step1 --filein "\$inputFileNames" --fileout "\$outputFileName" --pileup_input "PILEUPINPUT" --mc --eventcontent PREMIXRAW --datatier GEN-SIM-RAW --conditions 94X_mc2017_realistic_v10 --step DIGIPREMIX_S2,DATAMIX,L1,DIGI2RAW,HLT:2e34v40 --nThreads 8 --datamix PreMix --era Run2_2017 --python_filename ${VERTEX}_h3${h33}_h4${h44}_step1_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 1000
python edit_cfg_step1.py --filename=${VERTEX}_h3${h33}_h4${h44}_step1_cfg.py

CMSSWPATH=/afs/hep.wisc.edu/cms/jjbuchanan/aTGC_generation/CMSSW_9_4_0_patch1
HDFSDIR=/hdfs/store/user/jjbuch/aTGC
farmoutAnalysisJobs --skip-existing-output --use-singularity=SL6 --memory-requirements=5000 --input-dir=${HDFSDIR}/${VERTEX}/h3${h33}_h4${h44}_${JOBNAME} ${VERTEX}_h3${h33}_h4${h44}_step1_${JOBNAME} ${CMSSWPATH} ${CMSSWPATH}/src/${VERTEX}_h3${h33}_h4${h44}_step1_cfg.py 'outputFile=$outputFileName' 'inputFiles=$inputFileNames'


echo "*"
echo "*"
echo "Submitted step1 job ${VERTEX} h3=${h3} h4=${h4}"
echo "*"
echo
