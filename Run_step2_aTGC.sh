echo
echo "*"
echo "Assembling step2 jobs for ${VERTEX}, h3=${h3}, h4=${h4}"
echo "*"
echo "*"

h33=$(echo ${h3} | sed -e 's/\./p/' -e 's/-/m/')
h44=$(echo ${h4} | sed -e 's/\./p/' -e 's/-/m/')

cmsDriver.py step2 --filein "\$inputFileNames" --fileout "\$outputFileName" --mc --eventcontent AODSIM --runUnscheduled --datatier AODSIM --conditions 94X_mc2017_realistic_v10 --step RAW2DIGI,RECO,RECOSIM,EI --nThreads 8 --era Run2_2017 --python_filename ${VERTEX}_h3${h33}_h4${h44}_step2_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 1000
python edit_cfg_step2.py --filename=${VERTEX}_h3${h33}_h4${h44}_step2_cfg.py

CMSSWPATH=/afs/hep.wisc.edu/cms/jjbuchanan/aTGC_generation/CMSSW_9_4_0_patch1
HDFSDIR=/hdfs/store/user/jjbuchanan
STEP1NAME=${VERTEX}_h3${h33}_h4${h44}_step1_${JOBNAME}
farmoutAnalysisJobs --skip-existing-output --use-singularity=SL6 --memory-requirements=5000 --input-dir=${HDFSDIR}/${STEP1NAME}-${VERTEX}_h3${h33}_h4${h44}_step1_cfg ${VERTEX}_h3${h33}_h4${h44}_step2_${JOBNAME} ${CMSSWPATH} ${CMSSWPATH}/src/${VERTEX}_h3${h33}_h4${h44}_step2_cfg.py 'outputFile=$outputFileName' 'inputFiles=$inputFileNames'

echo "*"
echo "*"
echo "Submitted step2 jobs for ${VERTEX}, h3=${h3}, h4=${h4}"
echo "*"
echo
