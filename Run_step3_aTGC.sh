echo
echo "*"
echo "Assembling step3 jobs for ${VERTEX}, h3=${h3}, h4=${h4}"
echo "*"
echo "*"

h33=$(echo ${h3} | sed -e 's/\./p/' -e 's/-/m/')
h44=$(echo ${h4} | sed -e 's/\./p/' -e 's/-/m/')

cmsDriver.py step3 --filein "\$inputFileNames" --fileout "\$outputFileName" --mc --eventcontent MINIAODSIM --runUnscheduled --datatier MINIAODSIM --conditions 94X_mc2017_realistic_v10 --step PAT --nThreads 8 --era Run2_2017 --python_filename ${VERTEX}_h3${h33}_h4${h44}_step3_cfg.py --no_exec --customise Configuration/DataProcessing/Utils.addMonitoring -n 1000
python edit_cfg_step3.py --filename=${VERTEX}_h3${h33}_h4${h44}_step3_cfg.py

CMSSWPATH=/afs/hep.wisc.edu/cms/jjbuchanan/aTGC_generation/CMSSW_9_4_0_patch1
HDFSDIR=/hdfs/store/user/jjbuchanan
STEP2NAME=${VERTEX}_h3${h33}_h4${h44}_step2_${JOBNAME}
farmoutAnalysisJobs --skip-existing-output --use-singularity=SL6 --memory-requirements=5000 --input-dir=${HDFSDIR}/${STEP2NAME}-${VERTEX}_h3${h33}_h4${h44}_step2_cfg ${VERTEX}_h3${h33}_h4${h44}_step3_${JOBNAME} ${CMSSWPATH} ${CMSSWPATH}/src/${VERTEX}_h3${h33}_h4${h44}_step3_cfg.py 'outputFile=$outputFileName' 'inputFiles=$inputFileNames'

echo "*"
echo "*"
echo "Submitted step3 jobs for ${VERTEX}, h3=${h3}, h4=${h4}"
echo "*"
echo
