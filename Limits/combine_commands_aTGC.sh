#!/bin/sh

# Set up working area inside CMSSW_8_1_0/src/HiggsCombine
cmsenv

# Be sure all the needed histograms and transfer factors are in your working area
./copy_histos.sh

## Make data card ##
combineCards.py -S SA=datacard_signal_above0p5_aTGC.txt SB=datacard_signal_below0p5_aTGC.txt CR_ME=datacard_monoele.txt CR_MM=datacard_monomu.txt > comb_card_aTGC.txt

## Make workspace(s) ##
#  Choose which workspace(s) to make inside the script
#  Set connectWZ to false,
#   since the Zgamma and Wgamma estimates should not be connected
#   with a transfer factor for aTGC limits
root -l -q -b createWorkspaces_Pt_HaloFit.C
#  Now that a workspace has been made,
#    open comb_card_aTGC.txt and replace workspace.root with a specific workspace name

## Get limits and plots for one signal sample ##
combine -M AsymptoticLimits comb_card_aTGC.txt


# ------------------------------ #
# NOT YET TESTED FOR aTGC LIMITS #
# ------------------------------ #
# -- Run fit diagnostics --
text2workspace.py comb_card_aTGC.txt
combine comb_card_aTGC.root -M FitDiagnostics --saveShapes --saveWithUncertainties  --saveNormalizations
# -- Get total pre- and post-fit yields --
python ../test/mlfitNormsToText.py fitDiagnostics.root --uncertainties >& pre_and_post_fit_yields.txt
# -- Get magnitude of nuisance shifts --
python ../test/diffNuisances.py fitDiagnostics.root -g absoluteNuisanceDifferences.root
# -- Get pulls --
python ../test/diffNuisances.py fitDiagnostics.root -g pulls.root --pullDef relDiffAsymErrs
# -- Make phoPt plots --
root -l -q -b prefit_znng_SA_plotter.C
root -l -q -b prefit_znng_SB_plotter.C
root -l -q -b prefit_weng_plotter.C
root -l -q -b prefit_wmng_plotter.C
root -l -q -b prefit_zeeg_plotter.C
root -l -q -b prefit_zmmg_plotter.C
#  * To make background-only postfit plots, make sure shapes_fit_b (and not shapes_fit_s) is used everywhere
#    and make sure "SA_b_" (and not "SA_s_") is added to plotTitle.
#    For signal+background fits, do the opposite.
root -l -q -b postfit_znng_SA_plotter.C
root -l -q -b postfit_znng_SB_plotter.C
root -l -q -b postfit_weng_plotter.C
root -l -q -b postfit_wmng_plotter.C
root -l -q -b postfit_zeeg_plotter.C
root -l -q -b postfit_zmmg_plotter.C
# -- Impacts --
#  * CombineHarvester needs to be installed, according to the instructions above
#  * Remove -t -1 for unblinded impacts
text2workspace.py comb_card_aTGC.txt -m 125 -o comb_card_aTGC.root
combineTool.py -M Impacts -d comb_card_aTGC.root -m 125 --robustFit 1 --expectSignal=0 -t -1 --doInitialFit
combineTool.py -M Impacts -d comb_card_aTGC.root -m 125 --robustFit 1 --expectSignal=0 -t -1 --doFits --parallel 24
combineTool.py -M Impacts -d comb_card_aTGC.root -m 125 -o impacts.json_comb_card
plotImpacts.py -i impacts.json_comb_card -o impacts_comb_card
