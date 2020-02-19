# First, go to desired base directory for installation.
# Then run these steps:
export SCRAM_ARCH=slc6_amd64_gcc530
cmsrel CMSSW_8_1_0
cd CMSSW_8_1_0/src
cmsenv
git clone https://github.com/cms-analysis/HiggsAnalysis-CombinedLimit.git HiggsAnalysis/CombinedLimit
cd HiggsAnalysis/CombinedLimit
git fetch origin
git checkout v7.0.7
scramv1 b clean
scramv1 b

#  Install the portion of CombineHarvester needed to run impacts
#   Starting from CMSSW_8_1_0/src/HiggsAnalysis/CombinedLimit:
cd ../..
bash <(curl -s https://raw.githubusercontent.com/cms-analysis/CombineHarvester/master/CombineTools/scripts/sparse-checkout-ssh.sh)
scram b -j 8
