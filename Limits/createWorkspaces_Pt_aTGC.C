#include <vector>
#include <utility>
#include <algorithm>
#include "TFile.h"
#include "TMath.h"
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

void addTemplate(string procname, RooArgList& varlist, RooWorkspace& ws, TH1F* hist) {
  RooDataHist rhist(procname.c_str(), "", varlist, hist);
  ws.import(rhist);
}

void makeBinList(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* hist, RooArgList& binlist, bool setConst = false) {
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    stringstream binss;
    binss << procname << "_bin" << i;
    RooRealVar* binvar;
    if (!setConst) binvar = new RooRealVar(binss.str().c_str(), "", hist->GetBinContent(i), 0.0, 2000.0);
    else           binvar = new RooRealVar(binss.str().c_str(), "", hist->GetBinContent(i));
    binlist.add(*binvar);
  }

  stringstream normss;
  normss << procname << "_norm";

  RooParametricHist phist(procname.c_str(), "", var, binlist, *hist);
  RooAddition norm(normss.str().c_str(), "", binlist);

  // ws.import(phist,RooFit::RecycleConflictNodes());
  ws.import(phist);
  ws.import(norm, RooFit::RecycleConflictNodes());

}


void makeBinList_halo(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* hist, vector<pair<RooRealVar*, TH1*>> syst, RooArgList& binlist) {
  RooRealVar* mu_halo_monophHighPhi_scale = new RooRealVar("mu_halo_monophHighPhi_scale", "mu_halo_monophHighPhi_scale", 1.0, 0.0, 5.0);
  RooRealVar* shape_err = new RooRealVar("Halo_SR_above0p5_MIPTotalEnergy", "Halo_SR_above0p5_MIPTotalEnergy", 0.0, -5.0, 5.0);
  for (int i = 1; i <= hist->GetNbinsX(); i++) {
    stringstream rawmuss;
    rawmuss << procname << "_rawmu_bin" << i;
    RooRealVar* rawmu_bin = new RooRealVar(rawmuss.str().c_str(), "", hist->GetBinContent(i)); // Constant
    RooArgList fobinlist;
    fobinlist.add(*mu_halo_monophHighPhi_scale);
    fobinlist.add(*rawmu_bin);
    fobinlist.add(*shape_err);
    
    stringstream binss;
    binss << procname << "_bin" << i;
    stringstream formss;
    formss << "@0*@1";
    formss << "*(TMath::Exp(" << hist->GetBinError(i)/hist->GetBinContent(i) << "*@2))";
    for (int j = 0; j < syst.size(); j++) {
      stringstream systbinss;
      if (syst[j].first == NULL) {
        systbinss << procname << "_" << syst[j].second->GetName() << "_bin" << i;
        RooRealVar* systbinvar = new RooRealVar(systbinss.str().c_str(), "", 0., -5., 5.);
        fobinlist.add(*systbinvar);
      }
      else {
        fobinlist.add(*syst[j].first);
      }
      formss << "*(TMath::Exp(" << syst[j].second->GetBinContent(i) << "*@" << j+3 << "))";
    }
    RooFormulaVar* binvar = new RooFormulaVar(binss.str().c_str(), "", formss.str().c_str(), RooArgList(fobinlist));
    binlist.add(*binvar);
  }

  stringstream normss;
  normss << procname << "_norm";

  RooParametricHist phist(procname.c_str(), "", var, binlist, *hist);
  RooAddition norm(normss.str().c_str(), "", binlist);

  // ws.import(phist,RooFit::RecycleConflictNodes());
  ws.import(phist);
  ws.import(norm, RooFit::RecycleConflictNodes());

}
// ///////////////////////
// // RooProduct version
// ///////////////////////
// void makeBinList_halo(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* hist, RooArgList& binlist) {
//   RooRealVar* mu_halo_monophHighPhi_scale = new RooRealVar("mu_halo_monophHighPhi_scale", "mu_halo_monophHighPhi_scale", 1.0, 0.0, 5.0);
//   for (int i = 1; i <= hist->GetNbinsX(); i++) {
//     stringstream rawmuss;
//     rawmuss << procname << "_rawmu_bin" << i;
//     stringstream prodss;
//     prodss << procname << "_prod_bin" << i;
//     RooRealVar* rawmu_bin = new RooRealVar(rawmuss.str().c_str(), "", hist->GetBinContent(i)); // Constant
//     RooProduct* binvar = new RooProduct(prodss.str().c_str(), "", RooArgSet(*mu_halo_monophHighPhi_scale,*rawmu_bin));
//     binlist.add(*binvar);
//   }

//   stringstream normss;
//   normss << procname << "_norm";

//   RooParametricHist phist(procname.c_str(), "", var, binlist, *hist);
//   RooAddition norm(normss.str().c_str(), "", binlist);

//   // ws.import(phist,RooFit::RecycleConflictNodes());
//   ws.import(phist);
//   ws.import(norm, RooFit::RecycleConflictNodes());

// }

void makeConnectedBinList(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* rhist, vector<pair<RooRealVar*, TH1*>> syst, const RooArgList& srbinlist, RooArgList* crbinlist=NULL) {
  if (crbinlist == NULL) crbinlist = new RooArgList();

  for (int i = 1; i <= rhist->GetNbinsX(); i++) {
    stringstream rbinss;
    rbinss << "r_" << procname << "_bin" << i;
    RooRealVar* rbinvar = new RooRealVar(rbinss.str().c_str(), "", rhist->GetBinContent(i));

    stringstream rerrbinss;
    rerrbinss << procname << "_bin" << i << "_Runc";
    RooRealVar* rerrbinvar = new RooRealVar(rerrbinss.str().c_str(), "", 0., -5., 5.);

    stringstream binss;
    binss << procname << "_bin" << i;

    RooArgList fobinlist;
    fobinlist.add(srbinlist[i-1]);
    fobinlist.add(*rbinvar);
    fobinlist.add(*rerrbinvar);

    stringstream formss;
    formss << "@0/";
    formss << "(";
    formss << "@1";
    formss << "*(TMath::Exp(" << rhist->GetBinError(i)/rhist->GetBinContent(i) << "*@2))";
    for (int j = 0; j < syst.size(); j++) {
      stringstream systbinss;
      if (syst[j].first == NULL) {
  systbinss << procname << "_" << syst[j].second->GetName() << "_bin" << i;
  RooRealVar* systbinvar = new RooRealVar(systbinss.str().c_str(), "", 0., -5., 5.);
  fobinlist.add(*systbinvar);
      }
      else {
  fobinlist.add(*syst[j].first);
      }
      formss << "*(TMath::Exp(" << syst[j].second->GetBinContent(i) << "*@" << j+3 << "))";
    }
    formss << ")";

    RooFormulaVar* binvar = new RooFormulaVar(binss.str().c_str(), "", formss.str().c_str(), RooArgList(fobinlist));
    crbinlist->add(*binvar);
  }

  stringstream normss;
  normss << procname << "_norm";

  RooParametricHist phist(procname.c_str(), "", var, *crbinlist, *rhist);
  RooAddition norm(normss.str().c_str(),"", *crbinlist);

  ws.import(phist, RooFit::RecycleConflictNodes());
  ws.import(norm, RooFit::RecycleConflictNodes());
}

void makeConnectedBinList_noStatErr(string procname, RooRealVar& var, RooWorkspace& ws, TH1F* rhist, vector<pair<RooRealVar*, TH1*>> syst, const RooArgList& srbinlist, RooArgList* crbinlist=NULL) {
  if (crbinlist == NULL) crbinlist = new RooArgList();

  for (int i = 1; i <= rhist->GetNbinsX(); i++) {
    stringstream rbinss;
    rbinss << "r_" << procname << "_bin" << i;
    RooRealVar* rbinvar = new RooRealVar(rbinss.str().c_str(), "", rhist->GetBinContent(i));

    stringstream binss;
    binss << procname << "_bin" << i;

    RooArgList fobinlist;
    fobinlist.add(srbinlist[i-1]);
    fobinlist.add(*rbinvar);

    stringstream formss;
    formss << "@0/";
    formss << "(";
    formss << "@1";
    for (int j = 0; j < syst.size(); j++) {
      stringstream systbinss;
      if (syst[j].first == NULL) {
  systbinss << procname << "_" << syst[j].second->GetName() << "_bin" << i;
  RooRealVar* systbinvar = new RooRealVar(systbinss.str().c_str(), "", 0., -5., 5.);
  fobinlist.add(*systbinvar);
      }
      else {
  fobinlist.add(*syst[j].first);
      }
      formss << "*(TMath::Exp(" << syst[j].second->GetBinContent(i) << "*@" << j+3 << "))";
    }
    formss << ")";

    RooFormulaVar* binvar = new RooFormulaVar(binss.str().c_str(), "", formss.str().c_str(), RooArgList(fobinlist));
    crbinlist->add(*binvar);
  }

  stringstream normss;
  normss << procname << "_norm";

  RooParametricHist phist(procname.c_str(), "", var, *crbinlist, *rhist);
  RooAddition norm(normss.str().c_str(),"", *crbinlist);

  ws.import(phist, RooFit::RecycleConflictNodes());
  ws.import(norm, RooFit::RecycleConflictNodes());
}

void do_createWorkspace(string signal_samplename, bool connectWZ, vector<string> &signal_samplenames, vector<float> &signal_multipliers, string signal_histos_filename_stub, bool aTGC=false, float h3=0.0, float h4=0.0, vector<vector<vector<double>>> aTGC_2Dfits=vector<vector<vector<double>>>()){
  gSystem->Load("libHiggsAnalysisCombinedLimit.so");
  
  string outfile_stub = "workspace_Pt_";
  if(connectWZ) outfile_stub += "yesZW_";
  else outfile_stub += "noZW_";
  if(aTGC) outfile_stub += "aTGC_";
  outfile_stub += signal_samplename;
  if(aTGC) outfile_stub += signal_histos_filename_stub;
  
  TFile *outfile = new TFile(TString(outfile_stub+".root"),"RECREATE");
  RooWorkspace wspace("w","w");

  // RooRealVar mt("mt","M_{T}",0,1200);
  // RooRealVar pfmet("pfmet","pfMET",170,1000);
  RooRealVar phopt("phopt","E_{T}",175,1000);
  // RooArgList vars(mt);
  // RooArgList vars(pfmet);
  RooArgList vars(phopt);

  // Templates
  TFile* transfer_factors_file = new TFile("transfer_factors_Pt.root");
    
  TFile* signalabove0p5file;
  TFile* signalbelow0p5file;
  TFile* ZnnGNNLOratiofile;
  TFile* ZnnGNLOEWKratiofile;
  TFile* ZnnGSherToMadratiofile;
  if(!aTGC) {
    signalabove0p5file = new TFile(TString(signal_histos_filename_stub+"_histos_above0p5_Pt.root")); // DM, ADD
    signalbelow0p5file = new TFile(TString(signal_histos_filename_stub+"_histos_below0p5_Pt.root")); // DM, ADD
    // signalabove0p5file = new TFile("aTGC_histos_above0p5_Pt.root"); // aTGC TEST
    // signalbelow0p5file = new TFile("aTGC_histos_below0p5_Pt.root"); // aTGC TEST
  }
  else {
    signalabove0p5file = new TFile("aTGC_histos_above0p5_Pt.root");
    signalbelow0p5file = new TFile("aTGC_histos_below0p5_Pt.root");
    ZnnGNNLOratiofile = new TFile("ZnnG_NNLOratio_forATGC_Pt.root");
    ZnnGNLOEWKratiofile = new TFile("ZnnG_NLOEWKratio_forATGC_Pt.root");
    ZnnGSherToMadratiofile = new TFile(TString("ZnnG_SherToMadratio_forATGC_"+signal_samplename+"_Pt.root"));
  }
  
  TFile* ZnnGabove0p5file = new TFile("ZnnG_histos_above0p5_Pt.root");
  TFile* ZnnGbelow0p5file = new TFile("ZnnG_histos_below0p5_Pt.root");

  TFile* WenGfile = new TFile("WenG_histos_Pt.root");
  
  TFile* WmnGfile = new TFile("WmnG_histos_Pt.root");
  
  TFile* ZeeGfile = new TFile("ZeeG_histos_Pt.root");
  
  TFile* ZmmGfile = new TFile("ZmmG_histos_Pt.root");
  
  // ---------------------------- SIGNAL REGION (above0p5) -------------------------------------------------------------------//
  string procname_ZG_SA = "ZnunuG_SR_above0p5";
  string procname_WG_SA = "WG_SR_above0p5";
  string procname_halo_SA = "Halo_SR_above0p5";
  
  TH1F* data_obs_SR_above0p5 = (TH1F*)ZnnGabove0p5file->Get("data_obs");
  const int nBins = data_obs_SR_above0p5->GetNbinsX();
  
  // Data
  addTemplate("data_obs_SR_above0p5", vars, wspace, data_obs_SR_above0p5);
  
  TH1F* znng_NLOEWKratio;
  TH1F* znng_NLOEWKratio_straightUp;
  TH1F* znng_NLOEWKratio_straightDown;
  TH1F* znng_NLOEWKratio_twistedUp;
  TH1F* znng_NLOEWKratio_twistedDown;
  TH1F* znng_NLOEWKratio_gammaUp;
  TH1F* znng_NLOEWKratio_gammaDown;
  TH1F* znng_NNLOratio;
  TH1F* znng_NNLOratio_JESUp;
  TH1F* znng_NNLOratio_JESDown;
  TH1F* znng_NNLOratio_PESUp;
  TH1F* znng_NNLOratio_PESDown;
  TH1F* znng_SherToMadratio;
  TH1F* znng_SherToMadratio_JESUp;
  TH1F* znng_SherToMadratio_JESDown;
  TH1F* znng_SherToMadratio_PESUp;
  TH1F* znng_SherToMadratio_PESDown;
  if (aTGC) {
    znng_NLOEWKratio = (TH1F*)ZnnGNLOEWKratiofile->Get("histo_znng_NLOEWKratio");
    znng_NLOEWKratio_straightUp = (TH1F*)ZnnGNLOEWKratiofile->Get("histo_znng_NLOEWKratio_straightUp");
    znng_NLOEWKratio_straightDown = (TH1F*)ZnnGNLOEWKratiofile->Get("histo_znng_NLOEWKratio_straightDown");
    znng_NLOEWKratio_twistedUp = (TH1F*)ZnnGNLOEWKratiofile->Get("histo_znng_NLOEWKratio_twistedUp");
    znng_NLOEWKratio_twistedDown = (TH1F*)ZnnGNLOEWKratiofile->Get("histo_znng_NLOEWKratio_twistedDown");
    znng_NLOEWKratio_gammaUp = (TH1F*)ZnnGNLOEWKratiofile->Get("histo_znng_NLOEWKratio_gammaUp");
    znng_NLOEWKratio_gammaDown = (TH1F*)ZnnGNLOEWKratiofile->Get("histo_znng_NLOEWKratio_gammaDown");
    znng_NNLOratio = (TH1F*)ZnnGNNLOratiofile->Get("histo_znng_NNLOratio");
    znng_NNLOratio_JESUp = (TH1F*)ZnnGNNLOratiofile->Get("histo_znng_NNLOratio_JESUp");
    znng_NNLOratio_JESDown = (TH1F*)ZnnGNNLOratiofile->Get("histo_znng_NNLOratio_JESDown");
    znng_NNLOratio_PESUp = (TH1F*)ZnnGNNLOratiofile->Get("histo_znng_NNLOratio_PESUp");
    znng_NNLOratio_PESDown = (TH1F*)ZnnGNNLOratiofile->Get("histo_znng_NNLOratio_PESDown");
    znng_SherToMadratio = (TH1F*)ZnnGSherToMadratiofile->Get("histo_znng_SherToMadratio");
    znng_SherToMadratio_JESUp = (TH1F*)ZnnGSherToMadratiofile->Get("histo_znng_SherToMadratio_JESUp");
    znng_SherToMadratio_JESDown = (TH1F*)ZnnGSherToMadratiofile->Get("histo_znng_SherToMadratio_JESDown");
    znng_SherToMadratio_PESUp = (TH1F*)ZnnGSherToMadratiofile->Get("histo_znng_SherToMadratio_PESUp");
    znng_SherToMadratio_PESDown = (TH1F*)ZnnGSherToMadratiofile->Get("histo_znng_SherToMadratio_PESDown");
    // DEBUG
    for(int i = 1; i <= 6; i++){
      cout<<"i="<<i<<endl;
      cout<<"znng_NLOEWKratio: "<<znng_NLOEWKratio->GetBinContent(i)<<endl;
      cout<<"znng_NLOEWKratio_straightUp: "<<znng_NLOEWKratio_straightUp->GetBinContent(i)<<endl;
      cout<<"znng_NLOEWKratio_straightDown: "<<znng_NLOEWKratio_straightDown->GetBinContent(i)<<endl;
      cout<<"znng_NLOEWKratio_twistedUp: "<<znng_NLOEWKratio_twistedUp->GetBinContent(i)<<endl;
      cout<<"znng_NLOEWKratio_twistedDown: "<<znng_NLOEWKratio_twistedDown->GetBinContent(i)<<endl;
      cout<<"znng_NLOEWKratio_gammaUp: "<<znng_NLOEWKratio_gammaUp->GetBinContent(i)<<endl;
      cout<<"znng_NLOEWKratio_gammaDown: "<<znng_NLOEWKratio_gammaDown->GetBinContent(i)<<endl;
      cout<<"znng_NNLOratio: "<<znng_NNLOratio->GetBinContent(i)<<endl;
      cout<<"znng_NNLOratio +/- "<<znng_NNLOratio->GetBinError(i)<<endl;
      cout<<"znng_NNLOratio_JESUp: "<<znng_NNLOratio_JESUp->GetBinContent(i)<<endl;
      cout<<"znng_NNLOratio_JESDown: "<<znng_NNLOratio_JESDown->GetBinContent(i)<<endl;
      cout<<"znng_NNLOratio_PESUp: "<<znng_NNLOratio_PESUp->GetBinContent(i)<<endl;
      cout<<"znng_NNLOratio_PESDown: "<<znng_NNLOratio_PESDown->GetBinContent(i)<<endl;
      cout<<"znng_SherToMadratio: "<<znng_SherToMadratio->GetBinContent(i)<<endl;
      cout<<"znng_SherToMadratio +/- "<<znng_SherToMadratio->GetBinError(i)<<endl;
      cout<<"znng_SherToMadratio_JESUp: "<<znng_SherToMadratio_JESUp->GetBinContent(i)<<endl;
      cout<<"znng_SherToMadratio_JESDown: "<<znng_SherToMadratio_JESDown->GetBinContent(i)<<endl;
      cout<<"znng_SherToMadratio_PESUp: "<<znng_SherToMadratio_PESUp->GetBinContent(i)<<endl;
      cout<<"znng_SherToMadratio_PESDown: "<<znng_SherToMadratio_PESDown->GetBinContent(i)<<endl;
    }
  }
  
  // Signal shape
  TH1F* signal_SR_above0p5_hist;
  TH1F* signal_SR_above0p5_hist_statUp;
  TH1F* signal_SR_above0p5_hist_statDown;
  TH1F* signal_SR_above0p5_hist_qcdscaleUp;
  TH1F* signal_SR_above0p5_hist_qcdscaleDown;
  TH1F* signal_SR_above0p5_hist_straightUp;
  TH1F* signal_SR_above0p5_hist_straightDown;
  TH1F* signal_SR_above0p5_hist_twistedUp;
  TH1F* signal_SR_above0p5_hist_twistedDown;
  TH1F* signal_SR_above0p5_hist_gammaUp;
  TH1F* signal_SR_above0p5_hist_gammaDown;
  TH1F* signal_SR_above0p5_hist_JESUp;
  TH1F* signal_SR_above0p5_hist_JESDown;
  TH1F* signal_SR_above0p5_hist_PESUp;
  TH1F* signal_SR_above0p5_hist_PESDown;
  if (!aTGC) {
    signal_SR_above0p5_hist = (TH1F*)signalabove0p5file->Get(TString("histo_"+signal_samplename));
    signal_samplenames.push_back(signal_samplename); // For printout at the end
  }
  else {
    signal_SR_above0p5_hist = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename));
    signal_SR_above0p5_hist_statUp = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_statUp"));
    signal_SR_above0p5_hist_statDown = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_statDown"));
    signal_SR_above0p5_hist_qcdscaleUp = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_qcdscaleUp"));
    signal_SR_above0p5_hist_qcdscaleDown = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_qcdscaleDown"));
    signal_SR_above0p5_hist_straightUp = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_straightUp"));
    signal_SR_above0p5_hist_straightDown = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_straightDown"));
    signal_SR_above0p5_hist_twistedUp = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_twistedUp"));
    signal_SR_above0p5_hist_twistedDown = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_twistedDown"));
    signal_SR_above0p5_hist_gammaUp = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_gammaUp"));
    signal_SR_above0p5_hist_gammaDown = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_gammaDown"));
    signal_SR_above0p5_hist_JESUp = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_JESUp"));
    signal_SR_above0p5_hist_JESDown = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_JESDown"));
    signal_SR_above0p5_hist_PESUp = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_PESUp"));
    signal_SR_above0p5_hist_PESDown = (TH1F*)data_obs_SR_above0p5->Clone(TString("histo_"+signal_samplename+"_PESDown"));
    for(int bin_idx = 1; bin_idx <= nBins; bin_idx++) {
      vector<double> c = aTGC_2Dfits[0][bin_idx-1];
      vector<double> c_JESUp = aTGC_2Dfits[1][bin_idx-1];
      vector<double> c_JESDown = aTGC_2Dfits[2][bin_idx-1];
      vector<double> c_PESUp = aTGC_2Dfits[3][bin_idx-1];
      vector<double> c_PESDown = aTGC_2Dfits[4][bin_idx-1];
      signal_SR_above0p5_hist->SetBinContent(bin_idx, c[0]+c[1]*h3+c[2]*h4+c[3]*h3*h3+c[4]*h4*h4+c[5]*h3*h4);
      signal_SR_above0p5_hist_statUp->SetBinContent(bin_idx, c[0]+c[1]*h3+c[2]*h4+c[3]*h3*h3+c[4]*h4*h4+c[5]*h3*h4);
      signal_SR_above0p5_hist_statDown->SetBinContent(bin_idx, c[0]+c[1]*h3+c[2]*h4+c[3]*h3*h3+c[4]*h4*h4+c[5]*h3*h4);
      signal_SR_above0p5_hist_qcdscaleUp->SetBinContent(bin_idx, c[0]+c[1]*h3+c[2]*h4+c[3]*h3*h3+c[4]*h4*h4+c[5]*h3*h4);
      signal_SR_above0p5_hist_qcdscaleDown->SetBinContent(bin_idx, c[0]+c[1]*h3+c[2]*h4+c[3]*h3*h3+c[4]*h4*h4+c[5]*h3*h4);
      signal_SR_above0p5_hist_straightUp->SetBinContent(bin_idx, c[0]+c[1]*h3+c[2]*h4+c[3]*h3*h3+c[4]*h4*h4+c[5]*h3*h4);
      signal_SR_above0p5_hist_straightDown->SetBinContent(bin_idx, c[0]+c[1]*h3+c[2]*h4+c[3]*h3*h3+c[4]*h4*h4+c[5]*h3*h4);
      signal_SR_above0p5_hist_twistedUp->SetBinContent(bin_idx, c[0]+c[1]*h3+c[2]*h4+c[3]*h3*h3+c[4]*h4*h4+c[5]*h3*h4);
      signal_SR_above0p5_hist_twistedDown->SetBinContent(bin_idx, c[0]+c[1]*h3+c[2]*h4+c[3]*h3*h3+c[4]*h4*h4+c[5]*h3*h4);
      signal_SR_above0p5_hist_gammaUp->SetBinContent(bin_idx, c[0]+c[1]*h3+c[2]*h4+c[3]*h3*h3+c[4]*h4*h4+c[5]*h3*h4);
      signal_SR_above0p5_hist_gammaDown->SetBinContent(bin_idx, c[0]+c[1]*h3+c[2]*h4+c[3]*h3*h3+c[4]*h4*h4+c[5]*h3*h4);
      signal_SR_above0p5_hist_JESUp->SetBinContent(bin_idx, c_JESUp[0]+c_JESUp[1]*h3+c_JESUp[2]*h4+c_JESUp[3]*h3*h3+c_JESUp[4]*h4*h4+c_JESUp[5]*h3*h4);
      signal_SR_above0p5_hist_JESDown->SetBinContent(bin_idx, c_JESDown[0]+c_JESDown[1]*h3+c_JESDown[2]*h4+c_JESDown[3]*h3*h3+c_JESDown[4]*h4*h4+c_JESDown[5]*h3*h4);
      signal_SR_above0p5_hist_PESUp->SetBinContent(bin_idx, c_PESUp[0]+c_PESUp[1]*h3+c_PESUp[2]*h4+c_PESUp[3]*h3*h3+c_PESUp[4]*h4*h4+c_PESUp[5]*h3*h4);
      signal_SR_above0p5_hist_PESDown->SetBinContent(bin_idx, c_PESDown[0]+c_PESDown[1]*h3+c_PESDown[2]*h4+c_PESDown[3]*h3*h3+c_PESDown[4]*h4*h4+c_PESDown[5]*h3*h4);
      // Scale from Sherpa to MadGraph and from LO to NNLO QCD*NLO EWK
      // Sherpa to MadGraph
      float SherToMad_factor = znng_SherToMadratio->GetBinContent(bin_idx);
      float SherToMad_stat_error = znng_SherToMadratio->GetBinError(bin_idx);
      signal_SR_above0p5_hist->SetBinContent(bin_idx, signal_SR_above0p5_hist->GetBinContent(bin_idx)*SherToMad_factor);
      signal_SR_above0p5_hist_statUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_statUp->GetBinContent(bin_idx)*(SherToMad_factor+SherToMad_stat_error));
      signal_SR_above0p5_hist_statDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_statDown->GetBinContent(bin_idx)*(SherToMad_factor-SherToMad_stat_error));
      signal_SR_above0p5_hist_qcdscaleUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_qcdscaleUp->GetBinContent(bin_idx)*SherToMad_factor);
      signal_SR_above0p5_hist_qcdscaleDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_qcdscaleDown->GetBinContent(bin_idx)*SherToMad_factor);
      signal_SR_above0p5_hist_straightUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_straightUp->GetBinContent(bin_idx)*SherToMad_factor);
      signal_SR_above0p5_hist_straightDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_straightDown->GetBinContent(bin_idx)*SherToMad_factor);
      signal_SR_above0p5_hist_twistedUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_twistedUp->GetBinContent(bin_idx)*SherToMad_factor);
      signal_SR_above0p5_hist_twistedDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_twistedDown->GetBinContent(bin_idx)*SherToMad_factor);
      signal_SR_above0p5_hist_gammaUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_gammaUp->GetBinContent(bin_idx)*SherToMad_factor);
      signal_SR_above0p5_hist_gammaDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_gammaDown->GetBinContent(bin_idx)*SherToMad_factor);
      float SherToMad_factor_JESUp = znng_SherToMadratio_JESUp->GetBinContent(bin_idx);
      float SherToMad_factor_JESDown = znng_SherToMadratio_JESDown->GetBinContent(bin_idx);
      float SherToMad_factor_PESUp = znng_SherToMadratio_PESUp->GetBinContent(bin_idx);
      float SherToMad_factor_PESDown = znng_SherToMadratio_PESDown->GetBinContent(bin_idx);
      signal_SR_above0p5_hist_JESUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_JESUp->GetBinContent(bin_idx)*SherToMad_factor_JESUp);
      signal_SR_above0p5_hist_JESDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_JESDown->GetBinContent(bin_idx)*SherToMad_factor_JESDown);
      signal_SR_above0p5_hist_PESUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_PESUp->GetBinContent(bin_idx)*SherToMad_factor_PESUp);
      signal_SR_above0p5_hist_PESDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_PESDown->GetBinContent(bin_idx)*SherToMad_factor_PESDown);
      // NNLO QCD*NLO EWK
      float allCorr_factor = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio->GetBinContent(bin_idx);
      float allCorr_qcdscale_error = znng_NNLOratio->GetBinError(bin_idx)*znng_NLOEWKratio->GetBinContent(bin_idx);
      signal_SR_above0p5_hist->SetBinContent(bin_idx, signal_SR_above0p5_hist->GetBinContent(bin_idx)*allCorr_factor);
      signal_SR_above0p5_hist_statUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_statUp->GetBinContent(bin_idx)*allCorr_factor);
      signal_SR_above0p5_hist_statDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_statDown->GetBinContent(bin_idx)*allCorr_factor);
      signal_SR_above0p5_hist_qcdscaleUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_qcdscaleUp->GetBinContent(bin_idx)*(allCorr_factor+allCorr_qcdscale_error));
      signal_SR_above0p5_hist_qcdscaleDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_qcdscaleDown->GetBinContent(bin_idx)*(allCorr_factor-allCorr_qcdscale_error));
      float allCorr_factor_straightUp = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio_straightUp->GetBinContent(bin_idx);
      float allCorr_factor_straightDown = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio_straightDown->GetBinContent(bin_idx);
      float allCorr_factor_twistedUp = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio_twistedUp->GetBinContent(bin_idx);
      float allCorr_factor_twistedDown = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio_twistedDown->GetBinContent(bin_idx);
      float allCorr_factor_gammaUp = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio_gammaUp->GetBinContent(bin_idx);
      float allCorr_factor_gammaDown = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio_gammaDown->GetBinContent(bin_idx);
      float allCorr_factor_JESUp = znng_NNLOratio_JESUp->GetBinContent(bin_idx)*znng_NLOEWKratio->GetBinContent(bin_idx);
      float allCorr_factor_JESDown = znng_NNLOratio_JESDown->GetBinContent(bin_idx)*znng_NLOEWKratio->GetBinContent(bin_idx);
      float allCorr_factor_PESUp = znng_NNLOratio_PESUp->GetBinContent(bin_idx)*znng_NLOEWKratio->GetBinContent(bin_idx);
      float allCorr_factor_PESDown = znng_NNLOratio_PESDown->GetBinContent(bin_idx)*znng_NLOEWKratio->GetBinContent(bin_idx);
      signal_SR_above0p5_hist_straightUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_straightUp->GetBinContent(bin_idx)*allCorr_factor_straightUp);
      signal_SR_above0p5_hist_straightDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_straightDown->GetBinContent(bin_idx)*allCorr_factor_straightDown);
      signal_SR_above0p5_hist_twistedUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_twistedUp->GetBinContent(bin_idx)*allCorr_factor_twistedUp);
      signal_SR_above0p5_hist_twistedDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_twistedDown->GetBinContent(bin_idx)*allCorr_factor_twistedDown);
      signal_SR_above0p5_hist_gammaUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_gammaUp->GetBinContent(bin_idx)*allCorr_factor_gammaUp);
      signal_SR_above0p5_hist_gammaDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_gammaDown->GetBinContent(bin_idx)*allCorr_factor_gammaDown);
      signal_SR_above0p5_hist_JESUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_JESUp->GetBinContent(bin_idx)*allCorr_factor_JESUp);
      signal_SR_above0p5_hist_JESDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_JESDown->GetBinContent(bin_idx)*allCorr_factor_JESDown);
      signal_SR_above0p5_hist_PESUp->SetBinContent(bin_idx, signal_SR_above0p5_hist_PESUp->GetBinContent(bin_idx)*allCorr_factor_PESUp);
      signal_SR_above0p5_hist_PESDown->SetBinContent(bin_idx, signal_SR_above0p5_hist_PESDown->GetBinContent(bin_idx)*allCorr_factor_PESDown);
    }
    signal_samplenames.push_back("aTGC_"+signal_samplename+signal_histos_filename_stub);
  }
  float signal_sampleyield_above0p5 = signal_SR_above0p5_hist->Integral();
  float signal_multiplier = 1.0;
  // // Normalize DM yields so that limits are close to 1.
  // float signal_multiplier = 80.0 / signal_sampleyield_above0p5;
  // Save signal_multiplier so the limit can be properly scaled later.
  signal_multipliers.push_back(signal_multiplier);
  signal_SR_above0p5_hist->Scale(signal_multiplier);
  if (aTGC) {
    signal_SR_above0p5_hist_statUp->Scale(signal_multiplier);
    signal_SR_above0p5_hist_statDown->Scale(signal_multiplier);
    signal_SR_above0p5_hist_qcdscaleUp->Scale(signal_multiplier);
    signal_SR_above0p5_hist_qcdscaleDown->Scale(signal_multiplier);
    signal_SR_above0p5_hist_straightUp->Scale(signal_multiplier);
    signal_SR_above0p5_hist_straightDown->Scale(signal_multiplier);
    signal_SR_above0p5_hist_twistedUp->Scale(signal_multiplier);
    signal_SR_above0p5_hist_twistedDown->Scale(signal_multiplier);
    signal_SR_above0p5_hist_gammaUp->Scale(signal_multiplier);
    signal_SR_above0p5_hist_gammaDown->Scale(signal_multiplier);
    signal_SR_above0p5_hist_JESUp->Scale(signal_multiplier);
    signal_SR_above0p5_hist_JESDown->Scale(signal_multiplier);
    signal_SR_above0p5_hist_PESUp->Scale(signal_multiplier);
    signal_SR_above0p5_hist_PESDown->Scale(signal_multiplier);
  }

  // DEBUG
  if (aTGC) {
    cout<<"signal_SR_above0p5_hist"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_statUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_statUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_statDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_statDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_qcdscaleUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_qcdscaleUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_qcdscaleDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_qcdscaleDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_straightUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_straightUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_straightDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_straightDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_twistedUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_twistedUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_twistedDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_twistedDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_gammaUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_gammaUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_gammaDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_gammaDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_JESUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_JESUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_JESDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_JESDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_PESUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_PESUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_PESDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_PESDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<endl;
  }
    
  TH1F* znng_SR_above0p5_hist;
  TH1F* znng_SR_above0p5_hist_statUp;
  TH1F* znng_SR_above0p5_hist_statDown;
  TH1F* znng_SR_above0p5_hist_qcdscaleUp;
  TH1F* znng_SR_above0p5_hist_qcdscaleDown;
  TH1F* znng_SR_above0p5_hist_straightUp;
  TH1F* znng_SR_above0p5_hist_straightDown;
  TH1F* znng_SR_above0p5_hist_twistedUp;
  TH1F* znng_SR_above0p5_hist_twistedDown;
  TH1F* znng_SR_above0p5_hist_gammaUp;
  TH1F* znng_SR_above0p5_hist_gammaDown;
  TH1F* znng_SR_above0p5_hist_JESUp;
  TH1F* znng_SR_above0p5_hist_JESDown;
  TH1F* znng_SR_above0p5_hist_PESUp;
  TH1F* znng_SR_above0p5_hist_PESDown;
  RooArgList znng_SR_above0p5_bins;
  if(!aTGC) {
    znng_SR_above0p5_hist = (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG");
    makeBinList(procname_ZG_SA, phopt, wspace, znng_SR_above0p5_hist, znng_SR_above0p5_bins);
  }
  else{
    // Start with LO Sherpa yield
    znng_SR_above0p5_hist = (TH1F*)signalabove0p5file->Get(TString("histo_aTGC-"+signal_samplename+"_h3-0p0_h4-0p0"));
    znng_SR_above0p5_hist_statUp = (TH1F*)znng_SR_above0p5_hist->Clone("znng_SR_above0p5_hist_statUp");
    znng_SR_above0p5_hist_statDown = (TH1F*)znng_SR_above0p5_hist->Clone("znng_SR_above0p5_hist_statDown");
    znng_SR_above0p5_hist_qcdscaleUp = (TH1F*)znng_SR_above0p5_hist->Clone("znng_SR_above0p5_hist_qcdscaleUp");
    znng_SR_above0p5_hist_qcdscaleDown = (TH1F*)znng_SR_above0p5_hist->Clone("znng_SR_above0p5_hist_qcdscaleDown");
    znng_SR_above0p5_hist_straightUp = (TH1F*)znng_SR_above0p5_hist->Clone("znng_SR_above0p5_hist_straightUp");
    znng_SR_above0p5_hist_straightDown = (TH1F*)znng_SR_above0p5_hist->Clone("znng_SR_above0p5_hist_straightDown");
    znng_SR_above0p5_hist_twistedUp = (TH1F*)znng_SR_above0p5_hist->Clone("znng_SR_above0p5_hist_twistedUp");
    znng_SR_above0p5_hist_twistedDown = (TH1F*)znng_SR_above0p5_hist->Clone("znng_SR_above0p5_hist_twistedDown");
    znng_SR_above0p5_hist_gammaUp = (TH1F*)znng_SR_above0p5_hist->Clone("znng_SR_above0p5_hist_gammaUp");
    znng_SR_above0p5_hist_gammaDown = (TH1F*)znng_SR_above0p5_hist->Clone("znng_SR_above0p5_hist_gammaDown");
    znng_SR_above0p5_hist_JESUp = (TH1F*)signalabove0p5file->Get(TString("histo_aTGC-"+signal_samplename+"_h3-0p0_h4-0p0_JESUp"));
    znng_SR_above0p5_hist_JESDown = (TH1F*)signalabove0p5file->Get(TString("histo_aTGC-"+signal_samplename+"_h3-0p0_h4-0p0_JESDown"));
    znng_SR_above0p5_hist_PESUp = (TH1F*)signalabove0p5file->Get(TString("histo_aTGC-"+signal_samplename+"_h3-0p0_h4-0p0_PESUp"));
    znng_SR_above0p5_hist_PESDown = (TH1F*)signalabove0p5file->Get(TString("histo_aTGC-"+signal_samplename+"_h3-0p0_h4-0p0_PESDown"));
    // Scale from Sherpa to MadGraph and from LO to NNLO QCD*NLO EWK
    for (int bin_idx = 1; bin_idx <= nBins; bin_idx++) {
      // Sherpa to MadGraph
      float SherToMad_factor = znng_SherToMadratio->GetBinContent(bin_idx);
      float SherToMad_stat_error = znng_SherToMadratio->GetBinError(bin_idx);
      znng_SR_above0p5_hist->SetBinContent(bin_idx, znng_SR_above0p5_hist->GetBinContent(bin_idx)*SherToMad_factor);
      znng_SR_above0p5_hist_statUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_statUp->GetBinContent(bin_idx)*(SherToMad_factor+SherToMad_stat_error));
      znng_SR_above0p5_hist_statDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_statDown->GetBinContent(bin_idx)*(SherToMad_factor-SherToMad_stat_error));
      znng_SR_above0p5_hist_qcdscaleUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_qcdscaleUp->GetBinContent(bin_idx)*SherToMad_factor);
      znng_SR_above0p5_hist_qcdscaleDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_qcdscaleDown->GetBinContent(bin_idx)*SherToMad_factor);
      znng_SR_above0p5_hist_straightUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_straightUp->GetBinContent(bin_idx)*SherToMad_factor);
      znng_SR_above0p5_hist_straightDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_straightDown->GetBinContent(bin_idx)*SherToMad_factor);
      znng_SR_above0p5_hist_twistedUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_twistedUp->GetBinContent(bin_idx)*SherToMad_factor);
      znng_SR_above0p5_hist_twistedDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_twistedDown->GetBinContent(bin_idx)*SherToMad_factor);
      znng_SR_above0p5_hist_gammaUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_gammaUp->GetBinContent(bin_idx)*SherToMad_factor);
      znng_SR_above0p5_hist_gammaDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_gammaDown->GetBinContent(bin_idx)*SherToMad_factor);
      float SherToMad_factor_JESUp = znng_SherToMadratio_JESUp->GetBinContent(bin_idx);
      float SherToMad_factor_JESDown = znng_SherToMadratio_JESDown->GetBinContent(bin_idx);
      float SherToMad_factor_PESUp = znng_SherToMadratio_PESUp->GetBinContent(bin_idx);
      float SherToMad_factor_PESDown = znng_SherToMadratio_PESDown->GetBinContent(bin_idx);
      znng_SR_above0p5_hist_JESUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_JESUp->GetBinContent(bin_idx)*SherToMad_factor_JESUp);
      znng_SR_above0p5_hist_JESDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_JESDown->GetBinContent(bin_idx)*SherToMad_factor_JESDown);
      znng_SR_above0p5_hist_PESUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_PESUp->GetBinContent(bin_idx)*SherToMad_factor_PESUp);
      znng_SR_above0p5_hist_PESDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_PESDown->GetBinContent(bin_idx)*SherToMad_factor_PESDown);
      // NNLO QCD*NLO EWK
      float allCorr_factor = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio->GetBinContent(bin_idx);
      float allCorr_qcdscale_error = znng_NNLOratio->GetBinError(bin_idx)*znng_NLOEWKratio->GetBinContent(bin_idx);
      znng_SR_above0p5_hist->SetBinContent(bin_idx, znng_SR_above0p5_hist->GetBinContent(bin_idx)*allCorr_factor);
      znng_SR_above0p5_hist_statUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_statUp->GetBinContent(bin_idx)*allCorr_factor);
      znng_SR_above0p5_hist_statDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_statDown->GetBinContent(bin_idx)*allCorr_factor);
      znng_SR_above0p5_hist_qcdscaleUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_qcdscaleUp->GetBinContent(bin_idx)*(allCorr_factor+allCorr_qcdscale_error));
      znng_SR_above0p5_hist_qcdscaleDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_qcdscaleDown->GetBinContent(bin_idx)*(allCorr_factor-allCorr_qcdscale_error));
      float allCorr_factor_straightUp = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio_straightUp->GetBinContent(bin_idx);
      float allCorr_factor_straightDown = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio_straightDown->GetBinContent(bin_idx);
      float allCorr_factor_twistedUp = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio_twistedUp->GetBinContent(bin_idx);
      float allCorr_factor_twistedDown = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio_twistedDown->GetBinContent(bin_idx);
      float allCorr_factor_gammaUp = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio_gammaUp->GetBinContent(bin_idx);
      float allCorr_factor_gammaDown = znng_NNLOratio->GetBinContent(bin_idx)*znng_NLOEWKratio_gammaDown->GetBinContent(bin_idx);
      float allCorr_factor_JESUp = znng_NNLOratio_JESUp->GetBinContent(bin_idx)*znng_NLOEWKratio->GetBinContent(bin_idx);
      float allCorr_factor_JESDown = znng_NNLOratio_JESDown->GetBinContent(bin_idx)*znng_NLOEWKratio->GetBinContent(bin_idx);
      float allCorr_factor_PESUp = znng_NNLOratio_PESUp->GetBinContent(bin_idx)*znng_NLOEWKratio->GetBinContent(bin_idx);
      float allCorr_factor_PESDown = znng_NNLOratio_PESDown->GetBinContent(bin_idx)*znng_NLOEWKratio->GetBinContent(bin_idx);
      znng_SR_above0p5_hist_straightUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_straightUp->GetBinContent(bin_idx)*allCorr_factor_straightUp);
      znng_SR_above0p5_hist_straightDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_straightDown->GetBinContent(bin_idx)*allCorr_factor_straightDown);
      znng_SR_above0p5_hist_twistedUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_twistedUp->GetBinContent(bin_idx)*allCorr_factor_twistedUp);
      znng_SR_above0p5_hist_twistedDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_twistedDown->GetBinContent(bin_idx)*allCorr_factor_twistedDown);
      znng_SR_above0p5_hist_gammaUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_gammaUp->GetBinContent(bin_idx)*allCorr_factor_gammaUp);
      znng_SR_above0p5_hist_gammaDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_gammaDown->GetBinContent(bin_idx)*allCorr_factor_gammaDown);
      znng_SR_above0p5_hist_JESUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_JESUp->GetBinContent(bin_idx)*allCorr_factor_JESUp);
      znng_SR_above0p5_hist_JESDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_JESDown->GetBinContent(bin_idx)*allCorr_factor_JESDown);
      znng_SR_above0p5_hist_PESUp->SetBinContent(bin_idx, znng_SR_above0p5_hist_PESUp->GetBinContent(bin_idx)*allCorr_factor_PESUp);
      znng_SR_above0p5_hist_PESDown->SetBinContent(bin_idx, znng_SR_above0p5_hist_PESDown->GetBinContent(bin_idx)*allCorr_factor_PESDown);
    }
    addTemplate(procname_ZG_SA, vars, wspace, znng_SR_above0p5_hist);
    addTemplate(procname_ZG_SA+"_SherToMadStatUp", vars, wspace, znng_SR_above0p5_hist_statUp);
    addTemplate(procname_ZG_SA+"_SherToMadStatDown", vars, wspace, znng_SR_above0p5_hist_statDown);
    addTemplate(procname_ZG_SA+"_qcdscaleUp", vars, wspace, znng_SR_above0p5_hist_qcdscaleUp);
    addTemplate(procname_ZG_SA+"_qcdscaleDown", vars, wspace, znng_SR_above0p5_hist_qcdscaleDown);
    addTemplate(procname_ZG_SA+"_straightUp", vars, wspace, znng_SR_above0p5_hist_straightUp);
    addTemplate(procname_ZG_SA+"_straightDown", vars, wspace, znng_SR_above0p5_hist_straightDown);
    addTemplate(procname_ZG_SA+"_twistedUp", vars, wspace, znng_SR_above0p5_hist_twistedUp);
    addTemplate(procname_ZG_SA+"_twistedDown", vars, wspace, znng_SR_above0p5_hist_twistedDown);
    addTemplate(procname_ZG_SA+"_gammaUp", vars, wspace, znng_SR_above0p5_hist_gammaUp);
    addTemplate(procname_ZG_SA+"_gammaDown", vars, wspace, znng_SR_above0p5_hist_gammaDown);
    addTemplate(procname_ZG_SA+"_JESUp", vars, wspace, znng_SR_above0p5_hist_JESUp);
    addTemplate(procname_ZG_SA+"_JESDown", vars, wspace, znng_SR_above0p5_hist_JESDown);
    addTemplate(procname_ZG_SA+"_PESUp", vars, wspace, znng_SR_above0p5_hist_PESUp);
    addTemplate(procname_ZG_SA+"_PESDown", vars, wspace, znng_SR_above0p5_hist_PESDown);
  }

  // DEBUG
  if (aTGC) {
    cout<<"znng_SR_above0p5_hist"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_statUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_statUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_statDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_statDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_qcdscaleUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_qcdscaleUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_qcdscaleDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_qcdscaleDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_straightUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_straightUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_straightDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_straightDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_twistedUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_twistedUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_twistedDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_twistedDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_gammaUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_gammaUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_gammaDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_gammaDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_JESUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_JESUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_JESDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_JESDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_PESUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_PESUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"znng_SR_above0p5_hist_PESDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<znng_SR_above0p5_hist_PESDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<endl;
  }

  if(aTGC) {
    // Let the signal hists represent the difference between the expected ZnnG yields
    // with and without aTGCs
    signal_SR_above0p5_hist->Add(znng_SR_above0p5_hist, -1.0);
    signal_SR_above0p5_hist_statUp->Add(znng_SR_above0p5_hist_statUp, -1.0);
    signal_SR_above0p5_hist_statDown->Add(znng_SR_above0p5_hist_statDown, -1.0);
    signal_SR_above0p5_hist_qcdscaleUp->Add(znng_SR_above0p5_hist_qcdscaleUp, -1.0);
    signal_SR_above0p5_hist_qcdscaleDown->Add(znng_SR_above0p5_hist_qcdscaleDown, -1.0);
    signal_SR_above0p5_hist_straightUp->Add(znng_SR_above0p5_hist_straightUp, -1.0);
    signal_SR_above0p5_hist_straightDown->Add(znng_SR_above0p5_hist_straightDown, -1.0);
    signal_SR_above0p5_hist_twistedUp->Add(znng_SR_above0p5_hist_twistedUp, -1.0);
    signal_SR_above0p5_hist_twistedDown->Add(znng_SR_above0p5_hist_twistedDown, -1.0);
    signal_SR_above0p5_hist_gammaUp->Add(znng_SR_above0p5_hist_gammaUp, -1.0);
    signal_SR_above0p5_hist_gammaDown->Add(znng_SR_above0p5_hist_gammaDown, -1.0);
    signal_SR_above0p5_hist_JESUp->Add(znng_SR_above0p5_hist_JESUp, -1.0);
    signal_SR_above0p5_hist_JESDown->Add(znng_SR_above0p5_hist_JESDown, -1.0);
    signal_SR_above0p5_hist_PESUp->Add(znng_SR_above0p5_hist_PESUp, -1.0);
    signal_SR_above0p5_hist_PESDown->Add(znng_SR_above0p5_hist_PESDown, -1.0);
    addTemplate("Signal_SR_above0p5_SherToMadStatUp", vars, wspace, signal_SR_above0p5_hist_statUp);
    addTemplate("Signal_SR_above0p5_SherToMadStatDown", vars, wspace, signal_SR_above0p5_hist_statDown);
    addTemplate("Signal_SR_above0p5_qcdscaleUp", vars, wspace, signal_SR_above0p5_hist_qcdscaleUp);
    addTemplate("Signal_SR_above0p5_qcdscaleDown", vars, wspace, signal_SR_above0p5_hist_qcdscaleDown);
    addTemplate("Signal_SR_above0p5_straightUp", vars, wspace, signal_SR_above0p5_hist_straightUp);
    addTemplate("Signal_SR_above0p5_straightDown", vars, wspace, signal_SR_above0p5_hist_straightDown);
    addTemplate("Signal_SR_above0p5_twistedUp", vars, wspace, signal_SR_above0p5_hist_twistedUp);
    addTemplate("Signal_SR_above0p5_twistedDown", vars, wspace, signal_SR_above0p5_hist_twistedDown);
    addTemplate("Signal_SR_above0p5_gammaUp", vars, wspace, signal_SR_above0p5_hist_gammaUp);
    addTemplate("Signal_SR_above0p5_gammaDown", vars, wspace, signal_SR_above0p5_hist_gammaDown);
    addTemplate("Signal_SR_above0p5_JESUp", vars, wspace, signal_SR_above0p5_hist_JESUp);
    addTemplate("Signal_SR_above0p5_JESDown", vars, wspace, signal_SR_above0p5_hist_JESDown);
    addTemplate("Signal_SR_above0p5_PESUp", vars, wspace, signal_SR_above0p5_hist_PESUp);
    addTemplate("Signal_SR_above0p5_PESDown", vars, wspace, signal_SR_above0p5_hist_PESDown);
  }
  addTemplate("Signal_SR_above0p5", vars, wspace, signal_SR_above0p5_hist);
  
  // DEBUG
  if (aTGC) {
    cout<<"Now that SM has been subtracted:"<<endl;
    cout<<"--------------------------------"<<endl;
    cout<<"signal_SR_above0p5_hist"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_statUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_statUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_statDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_statDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_qcdscaleUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_qcdscaleUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_qcdscaleDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_qcdscaleDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_straightUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_straightUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_straightDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_straightDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_twistedUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_twistedUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_twistedDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_twistedDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_gammaUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_gammaUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_gammaDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_gammaDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_JESUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_JESUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_JESDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_JESDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_PESUp"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_PESUp->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<"signal_SR_above0p5_hist_PESDown"<<endl;
    for (int i = 1; i <= nBins; i++){ cout<<signal_SR_above0p5_hist_PESDown->GetBinContent(i)<<" "; }
    cout<<endl;
    cout<<endl;
  }

  TH1F* halo_SR_above0p5_hist = (TH1F*)ZnnGabove0p5file->Get("histo_bhalo");
  float halo_above0p5_int = halo_SR_above0p5_hist->Integral();
  TH1F* halo_SR_below0p5_hist = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_bhalo"))->Clone("histo_bhalo_below");
  halo_SR_above0p5_hist->Add(halo_SR_below0p5_hist);
  halo_SR_above0p5_hist->Scale(halo_above0p5_int/(halo_SR_above0p5_hist->Integral()));
  RooArgList halo_SR_above0p5_bins;
  TH1F* halo_SR_above0p5_MIPTotEnergy_shiftUp = (TH1F*)ZnnGabove0p5file->Get("histo_bhalo_MIPTotEnergyUp");
  float halo_above0p5_MIPTotEnergy_shiftUp_int = halo_SR_above0p5_MIPTotEnergy_shiftUp->Integral();
  TH1F* halo_SR_below0p5_MIPTotEnergy_shiftUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_bhalo_MIPTotEnergyUp"))->Clone("histo_bhalo_MIPTotEnergyUp_below");
  halo_SR_above0p5_MIPTotEnergy_shiftUp->Add(halo_SR_below0p5_MIPTotEnergy_shiftUp);
  halo_SR_above0p5_MIPTotEnergy_shiftUp->Scale(halo_above0p5_MIPTotEnergy_shiftUp_int/(halo_SR_above0p5_MIPTotEnergy_shiftUp->Integral()));
  TH1F* halo_SR_above0p5_MIPTotEnergy_shiftDown = (TH1F*)ZnnGabove0p5file->Get("histo_bhalo_MIPTotEnergyUp");
  float halo_above0p5_MIPTotEnergy_shiftDown_int = halo_SR_above0p5_MIPTotEnergy_shiftDown->Integral();
  TH1F* halo_SR_below0p5_MIPTotEnergy_shiftDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_bhalo_MIPTotEnergyUp"))->Clone("histo_bhalo_MIPTotEnergyUp_below");
  halo_SR_above0p5_MIPTotEnergy_shiftDown->Add(halo_SR_below0p5_MIPTotEnergy_shiftDown);
  halo_SR_above0p5_MIPTotEnergy_shiftDown->Scale(halo_above0p5_MIPTotEnergy_shiftDown_int/(halo_SR_above0p5_MIPTotEnergy_shiftDown->Integral()));
  TH1F* halo_SR_above0p5_MIPTotEnergy_fractionalShifts = (TH1F*)halo_SR_above0p5_MIPTotEnergy_shiftUp->Clone("halo_SR_above0p5_MIPTotEnergy_fractionalShifts");
  for (int i = 1; i <= halo_SR_above0p5_MIPTotEnergy_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = halo_SR_above0p5_MIPTotEnergy_shiftUp->GetBinContent(i)/halo_SR_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = halo_SR_above0p5_MIPTotEnergy_shiftDown->GetBinContent(i)/halo_SR_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    halo_SR_above0p5_MIPTotEnergy_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  vector<pair<RooRealVar*, TH1*>> halo_SR_above0p5_syst;
  RooRealVar* halo_SR_above0p5_MIPTotEnergy = new RooRealVar("Halo_SR_above0p5_MIPTotEnergy", "", 0., -5., 5.);
  halo_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(halo_SR_above0p5_MIPTotEnergy, halo_SR_above0p5_MIPTotEnergy_fractionalShifts));
  makeBinList_halo(procname_halo_SA, phopt, wspace, halo_SR_above0p5_hist, halo_SR_above0p5_syst, halo_SR_above0p5_bins);

  // Without Z/W link
  TH1F* wlng_SR_above0p5_hist = (TH1F*)ZnnGabove0p5file->Get("histo_WG");
  RooArgList wlng_SR_above0p5_bins;
  // With Z/W link
  TH1F* zoverw_SRr_above0p5_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_noShift");
  
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_qcdscale_shiftUp");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_qcdscale_shiftDown");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_shiftUp->Clone(TString(procname_WG_SA+"_ZNuNuGJets_qcdscale"));
  for (int i = 1; i <= zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_ewkscale_shiftUp");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_ewkscale_shiftDown");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_shiftUp->Clone(TString(procname_WG_SA+"_ZNuNuGJets_ewkscale"));
  for (int i = 1; i <= zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_ewkshape_shiftUp");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_ewkshape_shiftDown");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_shiftUp->Clone(TString(procname_WG_SA+"_ZNuNuGJets_ewkshape"));
  for (int i = 1; i <= zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_gamma_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_gamma_shiftUp");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_gamma_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_ZNuNuGJets_gamma_shiftDown");
  TH1F* zoverw_SRr_above0p5_ZNuNuGJets_gamma_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_ZNuNuGJets_gamma_shiftUp->Clone(TString(procname_WG_SA+"_ZNuNuGJets_gamma"));
  for (int i = 1; i <= zoverw_SRr_above0p5_ZNuNuGJets_gamma_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_ZNuNuGJets_gamma_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_ZNuNuGJets_gamma_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_ZNuNuGJets_gamma_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_WGJets_qcdscale_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_qcdscale_shiftUp");
  TH1F* zoverw_SRr_above0p5_WGJets_qcdscale_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_qcdscale_shiftDown");
  TH1F* zoverw_SRr_above0p5_WGJets_qcdscale_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_WGJets_qcdscale_shiftUp->Clone(TString(procname_WG_SA+"_WGJets_qcdscale"));
  for (int i = 1; i <= zoverw_SRr_above0p5_WGJets_qcdscale_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_WGJets_qcdscale_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_WGJets_qcdscale_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_WGJets_qcdscale_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_WGJets_ewkscale_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_ewkscale_shiftUp");
  TH1F* zoverw_SRr_above0p5_WGJets_ewkscale_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_ewkscale_shiftDown");
  TH1F* zoverw_SRr_above0p5_WGJets_ewkscale_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_WGJets_ewkscale_shiftUp->Clone(TString(procname_WG_SA+"_WGJets_ewkscale"));
  for (int i = 1; i <= zoverw_SRr_above0p5_WGJets_ewkscale_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_WGJets_ewkscale_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_WGJets_ewkscale_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_WGJets_ewkscale_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_WGJets_ewkshape_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_ewkshape_shiftUp");
  TH1F* zoverw_SRr_above0p5_WGJets_ewkshape_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_ewkshape_shiftDown");
  TH1F* zoverw_SRr_above0p5_WGJets_ewkshape_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_WGJets_ewkshape_shiftUp->Clone(TString(procname_WG_SA+"_WGJets_ewkshape"));
  for (int i = 1; i <= zoverw_SRr_above0p5_WGJets_ewkshape_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_WGJets_ewkshape_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_WGJets_ewkshape_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_WGJets_ewkshape_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  TH1F* zoverw_SRr_above0p5_WGJets_gamma_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_gamma_shiftUp");
  TH1F* zoverw_SRr_above0p5_WGJets_gamma_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_WGJets_gamma_shiftDown");
  TH1F* zoverw_SRr_above0p5_WGJets_gamma_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_WGJets_gamma_shiftUp->Clone(TString(procname_WG_SA+"_WGJets_gamma"));
  for (int i = 1; i <= zoverw_SRr_above0p5_WGJets_gamma_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_WGJets_gamma_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_WGJets_gamma_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_WGJets_gamma_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  // TH1F* zoverw_SRr_above0p5_anticorrelated_ewkscale_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_anticorrelated_ewkscale_shiftUp");
  // TH1F* zoverw_SRr_above0p5_anticorrelated_ewkscale_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_anticorrelated_ewkscale_shiftDown");
  // TH1F* zoverw_SRr_above0p5_anticorrelated_ewkscale_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_anticorrelated_ewkscale_shiftUp->Clone(TString(procname_WG_SA+"_anticorrelated_ewkscale"));
  // for (int i = 1; i <= zoverw_SRr_above0p5_anticorrelated_ewkscale_fractionalShifts->GetNbinsX(); i++) {
  //   Float_t upshift = zoverw_SRr_above0p5_anticorrelated_ewkscale_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
  //   Float_t downshift = zoverw_SRr_above0p5_anticorrelated_ewkscale_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
  //   Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
  //   zoverw_SRr_above0p5_anticorrelated_ewkscale_fractionalShifts->SetBinContent(i, shiftEnvelope);
  // }
  // TH1F* zoverw_SRr_above0p5_anticorrelated_ewkshape_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_anticorrelated_ewkshape_shiftUp");
  // TH1F* zoverw_SRr_above0p5_anticorrelated_ewkshape_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_anticorrelated_ewkshape_shiftDown");
  // TH1F* zoverw_SRr_above0p5_anticorrelated_ewkshape_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_anticorrelated_ewkshape_shiftUp->Clone(TString(procname_WG_SA+"_anticorrelated_ewkshape"));
  // for (int i = 1; i <= zoverw_SRr_above0p5_anticorrelated_ewkshape_fractionalShifts->GetNbinsX(); i++) {
  //   Float_t upshift = zoverw_SRr_above0p5_anticorrelated_ewkshape_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
  //   Float_t downshift = zoverw_SRr_above0p5_anticorrelated_ewkshape_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
  //   Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
  //   zoverw_SRr_above0p5_anticorrelated_ewkshape_fractionalShifts->SetBinContent(i, shiftEnvelope);
  // }
  // TH1F* zoverw_SRr_above0p5_anticorrelated_gamma_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_anticorrelated_gamma_shiftUp");
  // TH1F* zoverw_SRr_above0p5_anticorrelated_gamma_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_anticorrelated_gamma_shiftDown");
  // TH1F* zoverw_SRr_above0p5_anticorrelated_gamma_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_anticorrelated_gamma_shiftUp->Clone(TString(procname_WG_SA+"_anticorrelated_gamma"));
  // for (int i = 1; i <= zoverw_SRr_above0p5_anticorrelated_gamma_fractionalShifts->GetNbinsX(); i++) {
  //   Float_t upshift = zoverw_SRr_above0p5_anticorrelated_gamma_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
  //   Float_t downshift = zoverw_SRr_above0p5_anticorrelated_gamma_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
  //   Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
  //   zoverw_SRr_above0p5_anticorrelated_gamma_fractionalShifts->SetBinContent(i, shiftEnvelope);
  // }
  TH1F* zoverw_SRr_above0p5_pdf_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_pdf_shiftUp");
  TH1F* zoverw_SRr_above0p5_pdf_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinZnnG_to_ZnnGinZnnG_above0p5_pdf_shiftDown");
  TH1F* zoverw_SRr_above0p5_pdf_fractionalShifts = (TH1F*)zoverw_SRr_above0p5_pdf_shiftUp->Clone(TString(procname_WG_SA+"_pdf"));
  for (int i = 1; i <= zoverw_SRr_above0p5_pdf_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = zoverw_SRr_above0p5_pdf_shiftUp->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = zoverw_SRr_above0p5_pdf_shiftDown->GetBinContent(i)/zoverw_SRr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    zoverw_SRr_above0p5_pdf_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  vector<pair<RooRealVar*, TH1*>> zoverw_SR_above0p5_syst;
  RooRealVar* zoverw_SR_above0p5_pdf = new RooRealVar("ZNuNuGoverWG_SR_above0p5_pdf", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_pdf, zoverw_SRr_above0p5_pdf_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_ZNuNuGJets_qcdscale = new RooRealVar("ZNuNuGoverWG_SR_above0p5_ZNuNuGJets_qcdscale", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_ZNuNuGJets_qcdscale, zoverw_SRr_above0p5_ZNuNuGJets_qcdscale_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_ZNuNuGJets_ewkscale = new RooRealVar("ZNuNuGoverWG_SR_above0p5_ZNuNuGJets_ewkscale", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_ZNuNuGJets_ewkscale, zoverw_SRr_above0p5_ZNuNuGJets_ewkscale_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_ZNuNuGJets_ewkshape = new RooRealVar("ZNuNuGoverWG_SR_above0p5_ZNuNuGJets_ewkshape", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_ZNuNuGJets_ewkshape, zoverw_SRr_above0p5_ZNuNuGJets_ewkshape_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_ZNuNuGJets_gamma = new RooRealVar("ZNuNuGoverWG_SR_above0p5_ZNuNuGJets_gamma", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_ZNuNuGJets_gamma, zoverw_SRr_above0p5_ZNuNuGJets_gamma_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_WGJets_qcdscale = new RooRealVar("ZNuNuGoverWG_SR_above0p5_WGJets_qcdscale", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_WGJets_qcdscale, zoverw_SRr_above0p5_WGJets_qcdscale_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_WGJets_ewkscale = new RooRealVar("ZNuNuGoverWG_SR_above0p5_WGJets_ewkscale", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_WGJets_ewkscale, zoverw_SRr_above0p5_WGJets_ewkscale_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_WGJets_ewkshape = new RooRealVar("ZNuNuGoverWG_SR_above0p5_WGJets_ewkshape", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_WGJets_ewkshape, zoverw_SRr_above0p5_WGJets_ewkshape_fractionalShifts));
  RooRealVar* zoverw_SR_above0p5_WGJets_gamma = new RooRealVar("ZNuNuGoverWG_SR_above0p5_WGJets_gamma", "", 0., -5., 5.);
  zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_WGJets_gamma, zoverw_SRr_above0p5_WGJets_gamma_fractionalShifts));
  // RooRealVar* zoverw_SR_above0p5_anticorrelated_qcdscale = new RooRealVar("ZNuNuGoverWG_SR_above0p5_anticorrelated_qcdscale", "", 0., -5., 5.);
  // zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_anticorrelated_qcdscale, zoverw_SRr_above0p5_anticorrelated_qcdscale_fractionalShifts));
  // RooRealVar* zoverw_SR_above0p5_anticorrelated_ewkscale = new RooRealVar("ZNuNuGoverWG_SR_above0p5_anticorrelated_ewkscale", "", 0., -5., 5.);
  // zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_anticorrelated_ewkscale, zoverw_SRr_above0p5_anticorrelated_ewkscale_fractionalShifts));
  // RooRealVar* zoverw_SR_above0p5_anticorrelated_ewkshape = new RooRealVar("ZNuNuGoverWG_SR_above0p5_anticorrelated_ewkshape", "", 0., -5., 5.);
  // zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_anticorrelated_ewkshape, zoverw_SRr_above0p5_anticorrelated_ewkshape_fractionalShifts));
  // RooRealVar* zoverw_SR_above0p5_anticorrelated_gamma = new RooRealVar("ZNuNuGoverWG_SR_above0p5_anticorrelated_gamma", "", 0., -5., 5.);
  // zoverw_SR_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(zoverw_SR_above0p5_anticorrelated_gamma, zoverw_SRr_above0p5_anticorrelated_gamma_fractionalShifts));

  if (!connectWZ) makeBinList(procname_WG_SA, phopt, wspace, wlng_SR_above0p5_hist, wlng_SR_above0p5_bins);
  else   makeConnectedBinList(procname_WG_SA, phopt, wspace, zoverw_SRr_above0p5_hist, zoverw_SR_above0p5_syst, znng_SR_above0p5_bins, &wlng_SR_above0p5_bins);
  
  // TH1F* hist_ZNuNuG_SR_above0p5_corrected = (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_above0p5_uncorrected = (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_uncorrected"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_above0p5_qcdscale = (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_qcdscale"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_above0p5_EWKUp = (TH1F*)hist_ZNuNuG_SR_above0p5_corrected->Clone("hist_ZNuNuG_SR_above0p5_EWKUp"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_above0p5_EWKDown = (TH1F*)hist_ZNuNuG_SR_above0p5_corrected->Clone("hist_ZNuNuG_SR_above0p5_EWKDown"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_above0p5_qcdscaleUp = (TH1F*)hist_ZNuNuG_SR_above0p5_corrected->Clone("hist_ZNuNuG_SR_above0p5_qcdscaleUp"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_above0p5_qcdscaleDown = (TH1F*)hist_ZNuNuG_SR_above0p5_corrected->Clone("hist_ZNuNuG_SR_above0p5_qcdscaleDown"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_corrected = (TH1F*)ZnnGabove0p5file->Get("histo_WG"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_uncorrected = (TH1F*)ZnnGabove0p5file->Get("histo_WG_uncorrected"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_qcdscale = (TH1F*)ZnnGabove0p5file->Get("histo_WG_qcdscale"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_EWKUp = (TH1F*)hist_WG_SR_above0p5_corrected->Clone("hist_WG_SR_above0p5_EWKUp"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_EWKDown = (TH1F*)hist_WG_SR_above0p5_corrected->Clone("hist_WG_SR_above0p5_EWKDown"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_qcdscaleUp = (TH1F*)hist_WG_SR_above0p5_corrected->Clone("hist_WG_SR_above0p5_qcdscaleUp"); // Only if not fitting
  // TH1F* hist_WG_SR_above0p5_qcdscaleDown = (TH1F*)hist_WG_SR_above0p5_corrected->Clone("hist_WG_SR_above0p5_qcdscaleDown"); // Only if not fitting
  // for(int i = 1; i <= nBins; i++){ // Only if not fitting
  //   Float_t diff = fabs(hist_ZNuNuG_SR_above0p5_corrected->GetBinContent(i) - hist_ZNuNuG_SR_above0p5_uncorrected->GetBinContent(i)); // Only if not fitting
  //   hist_ZNuNuG_SR_above0p5_EWKUp->SetBinContent(i, hist_ZNuNuG_SR_above0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_ZNuNuG_SR_above0p5_EWKDown->SetBinContent(i, hist_ZNuNuG_SR_above0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  //   diff = fabs(hist_ZNuNuG_SR_above0p5_corrected->GetBinContent(i) - hist_ZNuNuG_SR_above0p5_qcdscale->GetBinContent(i)); // Only if not fitting
  //   hist_ZNuNuG_SR_above0p5_qcdscaleUp->SetBinContent(i, hist_ZNuNuG_SR_above0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_ZNuNuG_SR_above0p5_qcdscaleDown->SetBinContent(i, hist_ZNuNuG_SR_above0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  //   diff = fabs(hist_WG_SR_above0p5_corrected->GetBinContent(i) - hist_WG_SR_above0p5_uncorrected->GetBinContent(i)); // Only if not fitting
  //   hist_WG_SR_above0p5_EWKUp->SetBinContent(i, hist_WG_SR_above0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_WG_SR_above0p5_EWKDown->SetBinContent(i, hist_WG_SR_above0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  //   diff = fabs(hist_WG_SR_above0p5_corrected->GetBinContent(i) - hist_WG_SR_above0p5_qcdscale->GetBinContent(i)); // Only if not fitting
  //   hist_WG_SR_above0p5_qcdscaleUp->SetBinContent(i, hist_WG_SR_above0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_WG_SR_above0p5_qcdscaleDown->SetBinContent(i, hist_WG_SR_above0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  // } // Only if not fitting

  
  // Data driven backgrounds
  addTemplate("QCD_SR_above0p5"                , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_jetfake"));
  addTemplate("QCD_SR_above0p5_QCDrUp"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_jetfake_errUp"));
  addTemplate("QCD_SR_above0p5_QCDrDown"       , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_jetfake_errDown"));
  addTemplate("Elefake_SR_above0p5"            , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_elefake"));
  // addTemplate("BHalo_SR_above0p5"              , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_bhalo"));
  addTemplate("Spikes_SR_above0p5"             , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_spikes"));
  // MC backgrounds
  // addTemplate("ZNuNuG_SR_above0p5"              , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_JESUp")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_JESDown")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_PESUp")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_PESDown")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_phoSFUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_phoSFUp")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_phoSFDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG_phoSFDown")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_EWKUp"      , vars, wspace, hist_ZNuNuG_SR_above0p5_EWKUp); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_EWKDown"      , vars, wspace, hist_ZNuNuG_SR_above0p5_EWKDown); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_qcdscaleUp"      , vars, wspace, hist_ZNuNuG_SR_above0p5_qcdscaleUp); // Only if not fitting
  // addTemplate("ZNuNuG_SR_above0p5_qcdscaleDown"      , vars, wspace, hist_ZNuNuG_SR_above0p5_qcdscaleDown); // Only if not fitting
  // addTemplate("WG_SR_above0p5"              , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG_JESUp")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG_JESDown")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG_PESUp")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG_PESDown")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_phoSFUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG_phoSFUp")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_phoSFDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WG_phoSFDown")); // Only if not fitting
  // addTemplate("WG_SR_above0p5_EWKUp"      , vars, wspace, hist_WG_SR_above0p5_EWKUp); // Only if not fitting
  // addTemplate("WG_SR_above0p5_EWKDown"      , vars, wspace, hist_WG_SR_above0p5_EWKDown); // Only if not fitting
  // addTemplate("WG_SR_above0p5_qcdscaleUp"      , vars, wspace, hist_WG_SR_above0p5_qcdscaleUp); // Only if not fitting
  // addTemplate("WG_SR_above0p5_qcdscaleDown"      , vars, wspace, hist_WG_SR_above0p5_qcdscaleDown); // Only if not fitting
  addTemplate("GJets_SR_above0p5"              , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_GJets"));
  addTemplate("GJets_SR_above0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_GJets_JESUp"));
  addTemplate("GJets_SR_above0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_GJets_JESDown"));
  addTemplate("GJets_SR_above0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_GJets_PESUp"));
  addTemplate("GJets_SR_above0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_GJets_PESDown"));
  addTemplate("ZllG_SR_above0p5"               , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZllG_combined"));
  addTemplate("ZllG_SR_above0p5_JESUp"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZllG_JESUp_combined"));
  addTemplate("ZllG_SR_above0p5_JESDown"       , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZllG_JESDown_combined"));
  addTemplate("ZllG_SR_above0p5_PESUp"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZllG_PESUp_combined"));
  addTemplate("ZllG_SR_above0p5_PESDown"       , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZllG_PESDown_combined"));
  addTemplate("TTG_SR_above0p5"                , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TTG"));
  addTemplate("TTG_SR_above0p5_JESUp"          , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TTG_JESUp"));
  addTemplate("TTG_SR_above0p5_JESDown"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TTG_JESDown"));
  addTemplate("TTG_SR_above0p5_PESUp"          , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TTG_PESUp"));
  addTemplate("TTG_SR_above0p5_PESDown"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TTG_PESDown"));
  addTemplate("TG_SR_above0p5"                 , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TG"));
  addTemplate("TG_SR_above0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TG_JESUp"));
  addTemplate("TG_SR_above0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TG_JESDown"));
  addTemplate("TG_SR_above0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TG_PESUp"));
  addTemplate("TG_SR_above0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_TG_PESDown"));
  addTemplate("Diphoton_SR_above0p5"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_diphoton"));
  addTemplate("Diphoton_SR_above0p5_JESUp"     , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_diphoton_JESUp"));
  addTemplate("Diphoton_SR_above0p5_JESDown"   , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_diphoton_JESDown"));
  addTemplate("Diphoton_SR_above0p5_PESUp"     , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_diphoton_PESUp"));
  addTemplate("Diphoton_SR_above0p5_PESDown"   , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_diphoton_PESDown"));
  addTemplate("WZ_SR_above0p5"                 , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WZ"));
  addTemplate("WZ_SR_above0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WZ_JESUp"));
  addTemplate("WZ_SR_above0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WZ_JESDown"));
  addTemplate("WZ_SR_above0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WZ_PESUp"));
  addTemplate("WZ_SR_above0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WZ_PESDown"));
  addTemplate("ZZ_SR_above0p5"                 , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZZ"));
  addTemplate("ZZ_SR_above0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZZ_JESUp"));
  addTemplate("ZZ_SR_above0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZZ_JESDown"));
  addTemplate("ZZ_SR_above0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZZ_PESUp"));
  addTemplate("ZZ_SR_above0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_ZZ_PESDown"));
  addTemplate("WMuNu_SR_above0p5"              , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WMuNu"));
  addTemplate("WMuNu_SR_above0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WMuNu_JESUp"));
  addTemplate("WMuNu_SR_above0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WMuNu_JESDown"));
  addTemplate("WMuNu_SR_above0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WMuNu_PESUp"));
  addTemplate("WMuNu_SR_above0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WMuNu_PESDown"));
  addTemplate("WTauNu_SR_above0p5"             , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WTauNu"));
  addTemplate("WTauNu_SR_above0p5_JESUp"       , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WTauNu_JESUp"));
  addTemplate("WTauNu_SR_above0p5_JESDown"     , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WTauNu_JESDown"));
  addTemplate("WTauNu_SR_above0p5_PESUp"       , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WTauNu_PESUp"));
  addTemplate("WTauNu_SR_above0p5_PESDown"     , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WTauNu_PESDown"));
  addTemplate("WW_SR_above0p5"                 , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WW"));
  addTemplate("WW_SR_above0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WW_JESUp"));
  addTemplate("WW_SR_above0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WW_JESDown"));
  addTemplate("WW_SR_above0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WW_PESUp"));
  addTemplate("WW_SR_above0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGabove0p5file->Get("histo_WW_PESDown"));
  
  // ---------------------------- SIGNAL REGION (below0p5) -------------------------------------------------------------------//
  string procname_ZG_SB = "ZnunuG_SR_below0p5";
  string procname_WG_SB = "WG_SR_below0p5";
  string procname_halo_SB = "Halo_SR_below0p5";
  
  // Data
  addTemplate("data_obs_SR_below0p5", vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("data_obs"));
  
  // Signal shape
  TH1F* signal_SR_below0p5_hist;
  if(!aTGC) {
    signal_SR_below0p5_hist = (TH1F*)signalbelow0p5file->Get(TString("histo_"+signal_samplename));
    signal_SR_below0p5_hist->Scale(signal_multiplier);
  }
  else {
    signal_SR_below0p5_hist = (TH1F*)signal_SR_above0p5_hist->Clone(TString("histo_"+signal_samplename));
    TH1F* signal_SR_below0p5_hist_statUp = (TH1F*)signal_SR_above0p5_hist_statUp->Clone(TString("histo_"+signal_samplename+"_statUp"));
    TH1F* signal_SR_below0p5_hist_statDown = (TH1F*)signal_SR_above0p5_hist_statDown->Clone(TString("histo_"+signal_samplename+"_statDown"));
    TH1F* signal_SR_below0p5_hist_qcdscaleUp = (TH1F*)signal_SR_above0p5_hist_qcdscaleUp->Clone(TString("histo_"+signal_samplename+"_qcdscaleUp"));
    TH1F* signal_SR_below0p5_hist_qcdscaleDown = (TH1F*)signal_SR_above0p5_hist_qcdscaleDown->Clone(TString("histo_"+signal_samplename+"_qcdscaleDown"));
    TH1F* signal_SR_below0p5_hist_straightUp = (TH1F*)signal_SR_above0p5_hist_straightUp->Clone(TString("histo_"+signal_samplename+"_straightUp"));
    TH1F* signal_SR_below0p5_hist_straightDown = (TH1F*)signal_SR_above0p5_hist_straightDown->Clone(TString("histo_"+signal_samplename+"_straightDown"));
    TH1F* signal_SR_below0p5_hist_twistedUp = (TH1F*)signal_SR_above0p5_hist_twistedUp->Clone(TString("histo_"+signal_samplename+"_twistedUp"));
    TH1F* signal_SR_below0p5_hist_twistedDown = (TH1F*)signal_SR_above0p5_hist_twistedDown->Clone(TString("histo_"+signal_samplename+"_twistedDown"));
    TH1F* signal_SR_below0p5_hist_gammaUp = (TH1F*)signal_SR_above0p5_hist_gammaUp->Clone(TString("histo_"+signal_samplename+"_gammaUp"));
    TH1F* signal_SR_below0p5_hist_gammaDown = (TH1F*)signal_SR_above0p5_hist_gammaDown->Clone(TString("histo_"+signal_samplename+"_gammaDown"));
    TH1F* signal_SR_below0p5_hist_JESUp = (TH1F*)signal_SR_above0p5_hist_JESUp->Clone(TString("histo_"+signal_samplename+"_JESUp"));
    TH1F* signal_SR_below0p5_hist_JESDown = (TH1F*)signal_SR_above0p5_hist_JESDown->Clone(TString("histo_"+signal_samplename+"_JESDown"));
    TH1F* signal_SR_below0p5_hist_PESUp = (TH1F*)signal_SR_above0p5_hist_PESUp->Clone(TString("histo_"+signal_samplename+"_PESUp"));
    TH1F* signal_SR_below0p5_hist_PESDown = (TH1F*)signal_SR_above0p5_hist_PESDown->Clone(TString("histo_"+signal_samplename+"_PESDown"));
    signal_SR_below0p5_hist->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_statUp->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_statDown->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_qcdscaleUp->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_qcdscaleDown->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_straightUp->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_straightDown->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_twistedUp->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_twistedDown->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_gammaUp->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_gammaDown->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_JESUp->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_JESDown->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_PESUp->Scale(1.0/2.14159);
    signal_SR_below0p5_hist_PESDown->Scale(1.0/2.14159);
    addTemplate("Signal_SR_below0p5_SherToMadStatUp", vars, wspace, signal_SR_below0p5_hist_statUp);
    addTemplate("Signal_SR_below0p5_SherToMadStatDown", vars, wspace, signal_SR_below0p5_hist_statDown);
    addTemplate("Signal_SR_below0p5_qcdscaleUp", vars, wspace, signal_SR_below0p5_hist_qcdscaleUp);
    addTemplate("Signal_SR_below0p5_qcdscaleDown", vars, wspace, signal_SR_below0p5_hist_qcdscaleDown);
    addTemplate("Signal_SR_below0p5_straightUp", vars, wspace, signal_SR_below0p5_hist_straightUp);
    addTemplate("Signal_SR_below0p5_straightDown", vars, wspace, signal_SR_below0p5_hist_straightDown);
    addTemplate("Signal_SR_below0p5_twistedUp", vars, wspace, signal_SR_below0p5_hist_twistedUp);
    addTemplate("Signal_SR_below0p5_twistedDown", vars, wspace, signal_SR_below0p5_hist_twistedDown);
    addTemplate("Signal_SR_below0p5_gammaUp", vars, wspace, signal_SR_below0p5_hist_gammaUp);
    addTemplate("Signal_SR_below0p5_gammaDown", vars, wspace, signal_SR_below0p5_hist_gammaDown);
    addTemplate("Signal_SR_below0p5_JESUp", vars, wspace, signal_SR_below0p5_hist_JESUp);
    addTemplate("Signal_SR_below0p5_JESDown", vars, wspace, signal_SR_below0p5_hist_JESDown);
    addTemplate("Signal_SR_below0p5_PESUp", vars, wspace, signal_SR_below0p5_hist_PESUp);
    addTemplate("Signal_SR_below0p5_PESDown", vars, wspace, signal_SR_below0p5_hist_PESDown);
  }
  addTemplate("Signal_SR_below0p5", vars, wspace, signal_SR_below0p5_hist);
  
  TH1F* znng_SRr_aboveOverBelow_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_SBtoSA");
  vector<pair<RooRealVar*, TH1*>>   znng_SR_below0p5_syst;
  if(!aTGC) makeConnectedBinList_noStatErr(procname_ZG_SB, phopt, wspace, znng_SRr_aboveOverBelow_hist, znng_SR_below0p5_syst, znng_SR_above0p5_bins);
  else{
    TH1F* znng_SR_below0p5_hist = (TH1F*)znng_SR_above0p5_hist->Clone("histo_ZNuNuG_below");
    TH1F* znng_SR_below0p5_hist_statUp = (TH1F*)znng_SR_above0p5_hist_statUp->Clone("histo_ZNuNuG_below_statUp");
    TH1F* znng_SR_below0p5_hist_statDown = (TH1F*)znng_SR_above0p5_hist_statDown->Clone("histo_ZNuNuG_below_statDown");
    TH1F* znng_SR_below0p5_hist_qcdscaleUp = (TH1F*)znng_SR_above0p5_hist_qcdscaleUp->Clone("histo_ZNuNuG_below_qcdscaleUp");
    TH1F* znng_SR_below0p5_hist_qcdscaleDown = (TH1F*)znng_SR_above0p5_hist_qcdscaleDown->Clone("histo_ZNuNuG_below_qcdscaleDown");
    TH1F* znng_SR_below0p5_hist_straightUp = (TH1F*)znng_SR_above0p5_hist_straightUp->Clone("histo_ZNuNuG_below_straightUp");
    TH1F* znng_SR_below0p5_hist_straightDown = (TH1F*)znng_SR_above0p5_hist_straightDown->Clone("histo_ZNuNuG_below_straightDown");
    TH1F* znng_SR_below0p5_hist_twistedUp = (TH1F*)znng_SR_above0p5_hist_twistedUp->Clone("histo_ZNuNuG_below_twistedUp");
    TH1F* znng_SR_below0p5_hist_twistedDown = (TH1F*)znng_SR_above0p5_hist_twistedDown->Clone("histo_ZNuNuG_below_twistedDown");
    TH1F* znng_SR_below0p5_hist_gammaUp = (TH1F*)znng_SR_above0p5_hist_gammaUp->Clone("histo_ZNuNuG_below_gammaUp");
    TH1F* znng_SR_below0p5_hist_gammaDown = (TH1F*)znng_SR_above0p5_hist_gammaDown->Clone("histo_ZNuNuG_below_gammaDown");
    TH1F* znng_SR_below0p5_hist_JESUp = (TH1F*)znng_SR_above0p5_hist_JESUp->Clone("histo_ZNuNuG_below_JESUp");
    TH1F* znng_SR_below0p5_hist_JESDown = (TH1F*)znng_SR_above0p5_hist_JESDown->Clone("histo_ZNuNuG_below_JESDown");
    TH1F* znng_SR_below0p5_hist_PESUp = (TH1F*)znng_SR_above0p5_hist_PESUp->Clone("histo_ZNuNuG_below_PESUp");
    TH1F* znng_SR_below0p5_hist_PESDown = (TH1F*)znng_SR_above0p5_hist_PESDown->Clone("histo_ZNuNuG_below_PESDown");
    znng_SR_below0p5_hist->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_statUp->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_statDown->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_qcdscaleUp->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_qcdscaleDown->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_straightUp->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_straightDown->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_twistedUp->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_twistedDown->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_gammaUp->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_gammaDown->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_JESUp->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_JESDown->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_PESUp->Scale(1.0/2.14159);
    znng_SR_below0p5_hist_PESDown->Scale(1.0/2.14159);
    addTemplate(procname_ZG_SB, vars, wspace, znng_SR_below0p5_hist);
    addTemplate(procname_ZG_SB+"_SherToMadStatUp", vars, wspace, znng_SR_below0p5_hist_statUp);
    addTemplate(procname_ZG_SB+"_SherToMadStatDown", vars, wspace, znng_SR_below0p5_hist_statDown);
    addTemplate(procname_ZG_SB+"_qcdscaleUp", vars, wspace, znng_SR_below0p5_hist_qcdscaleUp);
    addTemplate(procname_ZG_SB+"_qcdscaleDown", vars, wspace, znng_SR_below0p5_hist_qcdscaleDown);
    addTemplate(procname_ZG_SB+"_straightUp", vars, wspace, znng_SR_below0p5_hist_straightUp);
    addTemplate(procname_ZG_SB+"_straightDown", vars, wspace, znng_SR_below0p5_hist_straightDown);
    addTemplate(procname_ZG_SB+"_twistedUp", vars, wspace, znng_SR_below0p5_hist_twistedUp);
    addTemplate(procname_ZG_SB+"_twistedDown", vars, wspace, znng_SR_below0p5_hist_twistedDown);
    addTemplate(procname_ZG_SB+"_gammaUp", vars, wspace, znng_SR_below0p5_hist_gammaUp);
    addTemplate(procname_ZG_SB+"_gammaDown", vars, wspace, znng_SR_below0p5_hist_gammaDown);
    addTemplate(procname_ZG_SB+"_JESUp", vars, wspace, znng_SR_below0p5_hist_JESUp);
    addTemplate(procname_ZG_SB+"_JESDown", vars, wspace, znng_SR_below0p5_hist_JESDown);
    addTemplate(procname_ZG_SB+"_PESUp", vars, wspace, znng_SR_below0p5_hist_PESUp);
    addTemplate(procname_ZG_SB+"_PESDown", vars, wspace, znng_SR_below0p5_hist_PESDown);
  }
  
  TH1F* wlng_SRr_aboveOverBelow_hist = (TH1F*)znng_SRr_aboveOverBelow_hist->Clone("wlng_SRr_aboveOverBelow_hist");
  vector<pair<RooRealVar*, TH1*>>   wlng_SR_below0p5_syst;
  makeConnectedBinList_noStatErr(procname_WG_SB, phopt, wspace, wlng_SRr_aboveOverBelow_hist, wlng_SR_below0p5_syst, wlng_SR_above0p5_bins);
  
  TH1F* halo_SRr_aboveOverBelow_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_SBtoSA_halo");
  vector<pair<RooRealVar*, TH1*>>   halo_SR_below0p5_syst;
  makeConnectedBinList(procname_halo_SB, phopt, wspace, halo_SRr_aboveOverBelow_hist, halo_SR_below0p5_syst, halo_SR_above0p5_bins);
  
  // TH1F* hist_ZNuNuG_SR_below0p5_corrected = (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_below0p5_uncorrected = (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_uncorrected"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_below0p5_qcdscale = (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_qcdscale"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_below0p5_EWKUp = (TH1F*)hist_ZNuNuG_SR_below0p5_corrected->Clone("hist_ZNuNuG_SR_below0p5_EWKUp"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_below0p5_EWKDown = (TH1F*)hist_ZNuNuG_SR_below0p5_corrected->Clone("hist_ZNuNuG_SR_below0p5_EWKDown"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_below0p5_qcdscaleUp = (TH1F*)hist_ZNuNuG_SR_below0p5_corrected->Clone("hist_ZNuNuG_SR_below0p5_qcdscaleUp"); // Only if not fitting
  // TH1F* hist_ZNuNuG_SR_below0p5_qcdscaleDown = (TH1F*)hist_ZNuNuG_SR_below0p5_corrected->Clone("hist_ZNuNuG_SR_below0p5_qcdscaleDown"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_corrected = (TH1F*)ZnnGbelow0p5file->Get("histo_WG"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_uncorrected = (TH1F*)ZnnGbelow0p5file->Get("histo_WG_uncorrected"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_qcdscale = (TH1F*)ZnnGbelow0p5file->Get("histo_WG_qcdscale"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_EWKUp = (TH1F*)hist_WG_SR_below0p5_corrected->Clone("hist_WG_SR_below0p5_EWKUp"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_EWKDown = (TH1F*)hist_WG_SR_below0p5_corrected->Clone("hist_WG_SR_below0p5_EWKDown"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_qcdscaleUp = (TH1F*)hist_WG_SR_below0p5_corrected->Clone("hist_WG_SR_below0p5_qcdscaleUp"); // Only if not fitting
  // TH1F* hist_WG_SR_below0p5_qcdscaleDown = (TH1F*)hist_WG_SR_below0p5_corrected->Clone("hist_WG_SR_below0p5_qcdscaleDown"); // Only if not fitting
  // for(int i = 1; i <= nBins; i++){ // Only if not fitting
  //   Float_t diff = fabs(hist_ZNuNuG_SR_below0p5_corrected->GetBinContent(i) - hist_ZNuNuG_SR_below0p5_uncorrected->GetBinContent(i)); // Only if not fitting
  //   hist_ZNuNuG_SR_below0p5_EWKUp->SetBinContent(i, hist_ZNuNuG_SR_below0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_ZNuNuG_SR_below0p5_EWKDown->SetBinContent(i, hist_ZNuNuG_SR_below0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  //   diff = fabs(hist_ZNuNuG_SR_below0p5_corrected->GetBinContent(i) - hist_ZNuNuG_SR_below0p5_qcdscale->GetBinContent(i)); // Only if not fitting
  //   hist_ZNuNuG_SR_below0p5_qcdscaleUp->SetBinContent(i, hist_ZNuNuG_SR_below0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_ZNuNuG_SR_below0p5_qcdscaleDown->SetBinContent(i, hist_ZNuNuG_SR_below0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  //   diff = fabs(hist_WG_SR_below0p5_corrected->GetBinContent(i) - hist_WG_SR_below0p5_uncorrected->GetBinContent(i)); // Only if not fitting
  //   hist_WG_SR_below0p5_EWKUp->SetBinContent(i, hist_WG_SR_below0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_WG_SR_below0p5_EWKDown->SetBinContent(i, hist_WG_SR_below0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  //   diff = fabs(hist_WG_SR_below0p5_corrected->GetBinContent(i) - hist_WG_SR_below0p5_qcdscale->GetBinContent(i)); // Only if not fitting
  //   hist_WG_SR_below0p5_qcdscaleUp->SetBinContent(i, hist_WG_SR_below0p5_corrected->GetBinContent(i) + diff); // Only if not fitting
  //   hist_WG_SR_below0p5_qcdscaleDown->SetBinContent(i, hist_WG_SR_below0p5_corrected->GetBinContent(i) - diff); // Only if not fitting
  // } // Only if not fitting
  
  // Data driven backgrounds
  addTemplate("QCD_SR_below0p5"                , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_jetfake"));
  addTemplate("QCD_SR_below0p5_QCDrUp"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_jetfake_errUp"));
  addTemplate("QCD_SR_below0p5_QCDrDown"       , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_jetfake_errDown"));
  addTemplate("Elefake_SR_below0p5"            , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_elefake"));
  // addTemplate("BHalo_SR_below0p5"              , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_bhalo"));
  addTemplate("Spikes_SR_below0p5"             , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_spikes"));
  // MC backgrounds
  // addTemplate("ZNuNuG_SR_below0p5"              , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_JESUp")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_JESDown")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_PESUp")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_PESDown")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_phoSFUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_phoSFUp")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_phoSFDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG_phoSFDown")); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_EWKUp"      , vars, wspace, hist_ZNuNuG_SR_below0p5_EWKUp); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_EWKDown"      , vars, wspace, hist_ZNuNuG_SR_below0p5_EWKDown); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_qcdscaleUp"      , vars, wspace, hist_ZNuNuG_SR_below0p5_qcdscaleUp); // Only if not fitting
  // addTemplate("ZNuNuG_SR_below0p5_qcdscaleDown"      , vars, wspace, hist_ZNuNuG_SR_below0p5_qcdscaleDown); // Only if not fitting
  // addTemplate("WG_SR_below0p5"              , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG_JESUp")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG_JESDown")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG_PESUp")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG_PESDown")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_phoSFUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG_phoSFUp")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_phoSFDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WG_phoSFDown")); // Only if not fitting
  // addTemplate("WG_SR_below0p5_EWKUp"      , vars, wspace, hist_WG_SR_below0p5_EWKUp); // Only if not fitting
  // addTemplate("WG_SR_below0p5_EWKDown"      , vars, wspace, hist_WG_SR_below0p5_EWKDown); // Only if not fitting
  // addTemplate("WG_SR_below0p5_qcdscaleUp"      , vars, wspace, hist_WG_SR_below0p5_qcdscaleUp); // Only if not fitting
  // addTemplate("WG_SR_below0p5_qcdscaleDown"      , vars, wspace, hist_WG_SR_below0p5_qcdscaleDown); // Only if not fitting
  addTemplate("GJets_SR_below0p5"              , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_GJets"));
  addTemplate("GJets_SR_below0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_GJets_JESUp"));
  addTemplate("GJets_SR_below0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_GJets_JESDown"));
  addTemplate("GJets_SR_below0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_GJets_PESUp"));
  addTemplate("GJets_SR_below0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_GJets_PESDown"));
  addTemplate("ZllG_SR_below0p5"               , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_combined"));
  addTemplate("ZllG_SR_below0p5_JESUp"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_JESUp_combined"));
  addTemplate("ZllG_SR_below0p5_JESDown"       , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_JESDown_combined"));
  addTemplate("ZllG_SR_below0p5_PESUp"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_PESUp_combined"));
  addTemplate("ZllG_SR_below0p5_PESDown"       , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_PESDown_combined"));
  addTemplate("TTG_SR_below0p5"                , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TTG"));
  addTemplate("TTG_SR_below0p5_JESUp"          , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TTG_JESUp"));
  addTemplate("TTG_SR_below0p5_JESDown"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TTG_JESDown"));
  addTemplate("TTG_SR_below0p5_PESUp"          , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TTG_PESUp"));
  addTemplate("TTG_SR_below0p5_PESDown"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TTG_PESDown"));
  addTemplate("TG_SR_below0p5"                 , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TG"));
  addTemplate("TG_SR_below0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TG_JESUp"));
  addTemplate("TG_SR_below0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TG_JESDown"));
  addTemplate("TG_SR_below0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TG_PESUp"));
  addTemplate("TG_SR_below0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_TG_PESDown"));
  addTemplate("Diphoton_SR_below0p5"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_diphoton"));
  addTemplate("Diphoton_SR_below0p5_JESUp"     , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_diphoton_JESUp"));
  addTemplate("Diphoton_SR_below0p5_JESDown"   , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_diphoton_JESDown"));
  addTemplate("Diphoton_SR_below0p5_PESUp"     , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_diphoton_PESUp"));
  addTemplate("Diphoton_SR_below0p5_PESDown"   , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_diphoton_PESDown"));
  addTemplate("WZ_SR_below0p5"                 , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WZ"));
  addTemplate("WZ_SR_below0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WZ_JESUp"));
  addTemplate("WZ_SR_below0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WZ_JESDown"));
  addTemplate("WZ_SR_below0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WZ_PESUp"));
  addTemplate("WZ_SR_below0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WZ_PESDown"));
  addTemplate("ZZ_SR_below0p5"                 , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZZ"));
  addTemplate("ZZ_SR_below0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZZ_JESUp"));
  addTemplate("ZZ_SR_below0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZZ_JESDown"));
  addTemplate("ZZ_SR_below0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZZ_PESUp"));
  addTemplate("ZZ_SR_below0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_ZZ_PESDown"));
  addTemplate("WMuNu_SR_below0p5"              , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu"));
  addTemplate("WMuNu_SR_below0p5_JESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu_JESUp"));
  addTemplate("WMuNu_SR_below0p5_JESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu_JESDown"));
  addTemplate("WMuNu_SR_below0p5_PESUp"        , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu_PESUp"));
  addTemplate("WMuNu_SR_below0p5_PESDown"      , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu_PESDown"));
  addTemplate("WTauNu_SR_below0p5"             , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu"));
  addTemplate("WTauNu_SR_below0p5_JESUp"       , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu_JESUp"));
  addTemplate("WTauNu_SR_below0p5_JESDown"     , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu_JESDown"));
  addTemplate("WTauNu_SR_below0p5_PESUp"       , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu_PESUp"));
  addTemplate("WTauNu_SR_below0p5_PESDown"     , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu_PESDown"));
  addTemplate("WW_SR_below0p5"                 , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WW"));
  addTemplate("WW_SR_below0p5_JESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WW_JESUp"));
  addTemplate("WW_SR_below0p5_JESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WW_JESDown"));
  addTemplate("WW_SR_below0p5_PESUp"           , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WW_PESUp"));
  addTemplate("WW_SR_below0p5_PESDown"         , vars, wspace, (TH1F*)ZnnGbelow0p5file->Get("histo_WW_PESDown"));

  // ---------------------------- CONTROL REGION (Dimuon) -----------------------------------------------------------------//
  if(!aTGC){
    string procname_ZM = "ZllG_ZM";
    
    addTemplate("data_obs_ZM", vars, wspace, (TH1F*)ZmmGfile->Get("data_obs"));
    
    TH1F* znng_ZMr_above0p5_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_ZllGinZmmG_to_ZnnGinZnnG_above0p5_noShift");
    TH1F* znng_ZMr_above0p5_muEff_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_ZllGinZmmG_to_ZnnGinZnnG_above0p5_muEff_shiftUp");
    TH1F* znng_ZMr_above0p5_muEff_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_ZllGinZmmG_to_ZnnGinZnnG_above0p5_muEff_shiftDown");
    TH1F* znng_ZMr_above0p5_muEff_fractionalShifts = (TH1F*)znng_ZMr_above0p5_muEff_shiftUp->Clone(TString(procname_ZM+"_muEff"));
    for (int i = 1; i <= znng_ZMr_above0p5_muEff_fractionalShifts->GetNbinsX(); i++) {
      Float_t upshift = znng_ZMr_above0p5_muEff_shiftUp->GetBinContent(i)/znng_ZMr_above0p5_hist->GetBinContent(i) - 1.0;
      Float_t downshift = znng_ZMr_above0p5_muEff_shiftDown->GetBinContent(i)/znng_ZMr_above0p5_hist->GetBinContent(i) - 1.0;
      Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
      znng_ZMr_above0p5_muEff_fractionalShifts->SetBinContent(i, shiftEnvelope);
    }
    vector<pair<RooRealVar*, TH1*>> znng_ZM_above0p5_syst;
    RooRealVar* znng_ZM_above0p5_muEff = new RooRealVar("ZNuNuGoverZLLG_ZM_above0p5_muEff", "", 0., -5., 5.);
    // znng_ZM_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(znng_ZM_above0p5_muEff, znng_ZMr_above0p5_muEff_fractionalShifts));
    if (!connectWZ || connectWZ) makeConnectedBinList(procname_ZM, phopt, wspace, znng_ZMr_above0p5_hist, znng_ZM_above0p5_syst, znng_SR_above0p5_bins);
    
    // Data driven backgrounds
    addTemplate("QCD_ZM"              , vars, wspace, (TH1F*)ZmmGfile->Get("histo_jetfake"));
    addTemplate("QCD_ZM_QCDrUp"       , vars, wspace, (TH1F*)ZmmGfile->Get("histo_jetfake_errUp"));
    addTemplate("QCD_ZM_QCDrDown"     , vars, wspace, (TH1F*)ZmmGfile->Get("histo_jetfake_errDown"));
    // MC backgrounds
    addTemplate("TTG_ZM"              , vars, wspace, (TH1F*)ZmmGfile->Get("histo_TTG"));
    addTemplate("TTG_ZM_JESUp"        , vars, wspace, (TH1F*)ZmmGfile->Get("histo_TTG_JESUp"));
    addTemplate("TTG_ZM_JESDown"      , vars, wspace, (TH1F*)ZmmGfile->Get("histo_TTG_JESDown"));
    addTemplate("TTG_ZM_PESUp"        , vars, wspace, (TH1F*)ZmmGfile->Get("histo_TTG_PESUp"));
    addTemplate("TTG_ZM_PESDown"      , vars, wspace, (TH1F*)ZmmGfile->Get("histo_TTG_PESDown"));
    addTemplate("WZ_ZM"               , vars, wspace, (TH1F*)ZmmGfile->Get("histo_WZ"));
    addTemplate("WZ_ZM_JESUp"         , vars, wspace, (TH1F*)ZmmGfile->Get("histo_WZ_JESUp"));
    addTemplate("WZ_ZM_JESDown"       , vars, wspace, (TH1F*)ZmmGfile->Get("histo_WZ_JESDown"));
    addTemplate("WZ_ZM_PESUp"         , vars, wspace, (TH1F*)ZmmGfile->Get("histo_WZ_PESUp"));
    addTemplate("WZ_ZM_PESDown"       , vars, wspace, (TH1F*)ZmmGfile->Get("histo_WZ_PESDown"));
  }

  // ---------------------------- CONTROL REGION (Dielectron) -----------------------------------------------------------------//
  if(!aTGC){
    string procname_ZE = "ZllG_ZE";
    
    addTemplate("data_obs_ZE", vars, wspace, (TH1F*)ZeeGfile->Get("data_obs"));
    
    TH1F* znng_ZEr_above0p5_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_ZllGinZeeG_to_ZnnGinZnnG_above0p5_noShift");
    TH1F* znng_ZEr_above0p5_eleEff_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_ZllGinZeeG_to_ZnnGinZnnG_above0p5_eleEff_shiftUp");
    TH1F* znng_ZEr_above0p5_eleEff_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_ZllGinZeeG_to_ZnnGinZnnG_above0p5_eleEff_shiftDown");
    TH1F* znng_ZEr_above0p5_eleEff_fractionalShifts = (TH1F*)znng_ZEr_above0p5_eleEff_shiftUp->Clone(TString(procname_ZE+"_eleEff"));
    for (int i = 1; i <= znng_ZEr_above0p5_eleEff_fractionalShifts->GetNbinsX(); i++) {
      Float_t upshift = znng_ZEr_above0p5_eleEff_shiftUp->GetBinContent(i)/znng_ZEr_above0p5_hist->GetBinContent(i) - 1.0;
      Float_t downshift = znng_ZEr_above0p5_eleEff_shiftDown->GetBinContent(i)/znng_ZEr_above0p5_hist->GetBinContent(i) - 1.0;
      Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
      znng_ZEr_above0p5_eleEff_fractionalShifts->SetBinContent(i, shiftEnvelope);
    }
    vector<pair<RooRealVar*, TH1*>> znng_ZE_above0p5_syst;
    RooRealVar* znng_ZE_above0p5_eleEff = new RooRealVar("ZNuNuGoverZLLG_ZE_above0p5_eleEff", "", 0., -5., 5.);
    // znng_ZE_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(znng_ZE_above0p5_eleEff, znng_ZEr_above0p5_eleEff_fractionalShifts));
    if (!connectWZ || connectWZ) makeConnectedBinList(procname_ZE, phopt, wspace, znng_ZEr_above0p5_hist, znng_ZE_above0p5_syst, znng_SR_above0p5_bins);
    
    // Data driven backgrounds
    addTemplate("QCD_ZE"              , vars, wspace, (TH1F*)ZeeGfile->Get("histo_jetfake"));
    addTemplate("QCD_ZE_QCDrUp"       , vars, wspace, (TH1F*)ZeeGfile->Get("histo_jetfake_errUp"));
    addTemplate("QCD_ZE_QCDrDown"     , vars, wspace, (TH1F*)ZeeGfile->Get("histo_jetfake_errDown"));
    // MC backgrounds
    addTemplate("TTG_ZE"              , vars, wspace, (TH1F*)ZeeGfile->Get("histo_TTG"));
    addTemplate("TTG_ZE_JESUp"        , vars, wspace, (TH1F*)ZeeGfile->Get("histo_TTG_JESUp"));
    addTemplate("TTG_ZE_JESDown"      , vars, wspace, (TH1F*)ZeeGfile->Get("histo_TTG_JESDown"));
    addTemplate("TTG_ZE_PESUp"        , vars, wspace, (TH1F*)ZeeGfile->Get("histo_TTG_PESUp"));
    addTemplate("TTG_ZE_PESDown"      , vars, wspace, (TH1F*)ZeeGfile->Get("histo_TTG_PESDown"));
    addTemplate("WZ_ZE"               , vars, wspace, (TH1F*)ZeeGfile->Get("histo_WZ"));
    addTemplate("WZ_ZE_JESUp"         , vars, wspace, (TH1F*)ZeeGfile->Get("histo_WZ_JESUp"));
    addTemplate("WZ_ZE_JESDown"       , vars, wspace, (TH1F*)ZeeGfile->Get("histo_WZ_JESDown"));
    addTemplate("WZ_ZE_PESUp"         , vars, wspace, (TH1F*)ZeeGfile->Get("histo_WZ_PESUp"));
    addTemplate("WZ_ZE_PESDown"       , vars, wspace, (TH1F*)ZeeGfile->Get("histo_WZ_PESDown"));
  }

  // ---------------------------- CONTROL REGION (Single muon) -----------------------------------------------------------------//
  string procname_WM = "WG_above0p5_WM";
  
  addTemplate("data_obs_WM"  , vars, wspace, (TH1F*)WmnGfile->Get("data_obs"));

  // Without Z/W link
  TH1F* wlng_WMr_above0p5_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinWmnG_to_WGinZnnG_above0p5_noShift");
  TH1F* wlng_WMr_above0p5_muEff_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinWmnG_to_WGinZnnG_above0p5_muEff_shiftUp");
  TH1F* wlng_WMr_above0p5_muEff_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinWmnG_to_WGinZnnG_above0p5_muEff_shiftDown");
  TH1F* wlng_WMr_above0p5_muEff_fractionalShifts = (TH1F*)wlng_WMr_above0p5_muEff_shiftUp->Clone(TString(procname_WM+"_muEff"));
  for (int i = 1; i <= wlng_WMr_above0p5_muEff_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = wlng_WMr_above0p5_muEff_shiftUp->GetBinContent(i)/wlng_WMr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = wlng_WMr_above0p5_muEff_shiftDown->GetBinContent(i)/wlng_WMr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    wlng_WMr_above0p5_muEff_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  vector<pair<RooRealVar*, TH1*>> wlng_WM_above0p5_syst;
  RooRealVar* wlng_WM_above0p5_muEff = new RooRealVar("WLNuG_WM_above0p5_muEff", "", 0., -5., 5.);
  // wlng_WM_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(wlng_WM_above0p5_muEff, wlng_WMr_above0p5_muEff_fractionalShifts));
  
  makeConnectedBinList(procname_WM, phopt, wspace, wlng_WMr_above0p5_hist, wlng_WM_above0p5_syst, wlng_SR_above0p5_bins);
  
  // Data driven backgrounds
  addTemplate("QCD_WM"              , vars, wspace, (TH1F*)WmnGfile->Get("histo_jetfake"));
  addTemplate("QCD_WM_QCDrUp"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_jetfake_errUp"));
  addTemplate("QCD_WM_QCDrDown"     , vars, wspace, (TH1F*)WmnGfile->Get("histo_jetfake_errDown"));
  // MC backgrounds
  addTemplate("ZllG_WM"             , vars, wspace, (TH1F*)WmnGfile->Get("histo_ZllG_combined"));
  addTemplate("ZllG_WM_JESUp"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_ZllG_JESUp_combined"));
  addTemplate("ZllG_WM_JESDown"     , vars, wspace, (TH1F*)WmnGfile->Get("histo_ZllG_JESDown_combined"));
  addTemplate("ZllG_WM_PESUp"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_ZllG_PESUp_combined"));
  addTemplate("ZllG_WM_PESDown"     , vars, wspace, (TH1F*)WmnGfile->Get("histo_ZllG_PESDown_combined"));
  addTemplate("TTG_WM"              , vars, wspace, (TH1F*)WmnGfile->Get("histo_TTG"));
  addTemplate("TTG_WM_JESUp"        , vars, wspace, (TH1F*)WmnGfile->Get("histo_TTG_JESUp"));
  addTemplate("TTG_WM_JESDown"      , vars, wspace, (TH1F*)WmnGfile->Get("histo_TTG_JESDown"));
  addTemplate("TTG_WM_PESUp"        , vars, wspace, (TH1F*)WmnGfile->Get("histo_TTG_PESUp"));
  addTemplate("TTG_WM_PESDown"      , vars, wspace, (TH1F*)WmnGfile->Get("histo_TTG_PESDown"));
  addTemplate("TG_WM"               , vars, wspace, (TH1F*)WmnGfile->Get("histo_TG"));
  addTemplate("TG_WM_JESUp"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_TG_JESUp"));
  addTemplate("TG_WM_JESDown"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_TG_JESDown"));
  addTemplate("TG_WM_PESUp"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_TG_PESUp"));
  addTemplate("TG_WM_PESDown"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_TG_PESDown"));
  addTemplate("Diphoton_WM"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_diphoton"));
  addTemplate("Diphoton_WM_JESUp"   , vars, wspace, (TH1F*)WmnGfile->Get("histo_diphoton_JESUp"));
  addTemplate("Diphoton_WM_JESDown" , vars, wspace, (TH1F*)WmnGfile->Get("histo_diphoton_JESDown"));
  addTemplate("Diphoton_WM_PESUp"   , vars, wspace, (TH1F*)WmnGfile->Get("histo_diphoton_PESUp"));
  addTemplate("Diphoton_WM_PESDown" , vars, wspace, (TH1F*)WmnGfile->Get("histo_diphoton_PESDown"));
  addTemplate("WZ_WM"               , vars, wspace, (TH1F*)WmnGfile->Get("histo_WZ"));
  addTemplate("WZ_WM_JESUp"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_WZ_JESUp"));
  addTemplate("WZ_WM_JESDown"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_WZ_JESDown"));
  addTemplate("WZ_WM_PESUp"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_WZ_PESUp"));
  addTemplate("WZ_WM_PESDown"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_WZ_PESDown"));
  addTemplate("WW_WM"               , vars, wspace, (TH1F*)WmnGfile->Get("histo_WW"));
  addTemplate("WW_WM_JESUp"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_WW_JESUp"));
  addTemplate("WW_WM_JESDown"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_WW_JESDown"));
  addTemplate("WW_WM_PESUp"         , vars, wspace, (TH1F*)WmnGfile->Get("histo_WW_PESUp"));
  addTemplate("WW_WM_PESDown"       , vars, wspace, (TH1F*)WmnGfile->Get("histo_WW_PESDown"));
  
  // ---------------------------- CONTROL REGION (Single electron) -----------------------------------------------------------------//
  string procname_WE = "WG_above0p5_WE";
  
  addTemplate("data_obs_WE"  , vars, wspace, (TH1F*)WenGfile->Get("data_obs"));

  // Without Z/W link
  TH1F* wlng_WEr_above0p5_hist = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinWenG_to_WGinZnnG_above0p5_noShift");
  TH1F* wlng_WEr_above0p5_eleEff_shiftUp = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinWenG_to_WGinZnnG_above0p5_eleEff_shiftUp");
  TH1F* wlng_WEr_above0p5_eleEff_shiftDown = (TH1F*)transfer_factors_file->Get("transfer_factor_WGinWenG_to_WGinZnnG_above0p5_eleEff_shiftDown");
  TH1F* wlng_WEr_above0p5_eleEff_fractionalShifts = (TH1F*)wlng_WEr_above0p5_eleEff_shiftUp->Clone(TString(procname_WE+"_eleEff"));
  for (int i = 1; i <= wlng_WEr_above0p5_eleEff_fractionalShifts->GetNbinsX(); i++) {
    Float_t upshift = wlng_WEr_above0p5_eleEff_shiftUp->GetBinContent(i)/wlng_WEr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t downshift = wlng_WEr_above0p5_eleEff_shiftDown->GetBinContent(i)/wlng_WEr_above0p5_hist->GetBinContent(i) - 1.0;
    Float_t shiftEnvelope = TMath::Max(fabs(upshift), fabs(downshift));
    wlng_WEr_above0p5_eleEff_fractionalShifts->SetBinContent(i, shiftEnvelope);
  }
  vector<pair<RooRealVar*, TH1*>> wlng_WE_above0p5_syst;
  RooRealVar* wlng_WE_above0p5_eleEff = new RooRealVar("WLNuG_WE_above0p5_eleEff", "", 0., -5., 5.);
  // wlng_WE_above0p5_syst.push_back(pair<RooRealVar*, TH1*>(wlng_WE_above0p5_eleEff, wlng_WEr_above0p5_eleEff_fractionalShifts));

  makeConnectedBinList(procname_WE, phopt, wspace, wlng_WEr_above0p5_hist, wlng_WE_above0p5_syst, wlng_SR_above0p5_bins);
  
  // Data driven backgrounds
  addTemplate("QCD_WE"              , vars, wspace, (TH1F*)WenGfile->Get("histo_jetfake"));
  addTemplate("QCD_WE_QCDrUp"       , vars, wspace, (TH1F*)WenGfile->Get("histo_jetfake_errUp"));
  addTemplate("QCD_WE_QCDrDown"     , vars, wspace, (TH1F*)WenGfile->Get("histo_jetfake_errDown"));
  addTemplate("Elefake_WE"          , vars, wspace, (TH1F*)WenGfile->Get("histo_elefake"));
  // MC backgrounds
  addTemplate("ZllG_WE"             , vars, wspace, (TH1F*)WenGfile->Get("histo_ZllG_combined"));
  addTemplate("ZllG_WE_JESUp"       , vars, wspace, (TH1F*)WenGfile->Get("histo_ZllG_JESUp_combined"));
  addTemplate("ZllG_WE_JESDown"     , vars, wspace, (TH1F*)WenGfile->Get("histo_ZllG_JESDown_combined"));
  addTemplate("ZllG_WE_PESUp"       , vars, wspace, (TH1F*)WenGfile->Get("histo_ZllG_PESUp_combined"));
  addTemplate("ZllG_WE_PESDown"     , vars, wspace, (TH1F*)WenGfile->Get("histo_ZllG_PESDown_combined"));
  addTemplate("TTG_WE"              , vars, wspace, (TH1F*)WenGfile->Get("histo_TTG"));
  addTemplate("TTG_WE_JESUp"        , vars, wspace, (TH1F*)WenGfile->Get("histo_TTG_JESUp"));
  addTemplate("TTG_WE_JESDown"      , vars, wspace, (TH1F*)WenGfile->Get("histo_TTG_JESDown"));
  addTemplate("TTG_WE_PESUp"        , vars, wspace, (TH1F*)WenGfile->Get("histo_TTG_PESUp"));
  addTemplate("TTG_WE_PESDown"      , vars, wspace, (TH1F*)WenGfile->Get("histo_TTG_PESDown"));
  addTemplate("TG_WE"               , vars, wspace, (TH1F*)WenGfile->Get("histo_TG"));
  addTemplate("TG_WE_JESUp"         , vars, wspace, (TH1F*)WenGfile->Get("histo_TG_JESUp"));
  addTemplate("TG_WE_JESDown"       , vars, wspace, (TH1F*)WenGfile->Get("histo_TG_JESDown"));
  addTemplate("TG_WE_PESUp"         , vars, wspace, (TH1F*)WenGfile->Get("histo_TG_PESUp"));
  addTemplate("TG_WE_PESDown"       , vars, wspace, (TH1F*)WenGfile->Get("histo_TG_PESDown"));
  addTemplate("Diphoton_WE"         , vars, wspace, (TH1F*)WenGfile->Get("histo_diphoton"));
  addTemplate("Diphoton_WE_JESUp"   , vars, wspace, (TH1F*)WenGfile->Get("histo_diphoton_JESUp"));
  addTemplate("Diphoton_WE_JESDown" , vars, wspace, (TH1F*)WenGfile->Get("histo_diphoton_JESDown"));
  addTemplate("Diphoton_WE_PESUp"   , vars, wspace, (TH1F*)WenGfile->Get("histo_diphoton_PESUp"));
  addTemplate("Diphoton_WE_PESDown" , vars, wspace, (TH1F*)WenGfile->Get("histo_diphoton_PESDown"));
  addTemplate("WZ_WE"               , vars, wspace, (TH1F*)WenGfile->Get("histo_WZ"));
  addTemplate("WZ_WE_JESUp"         , vars, wspace, (TH1F*)WenGfile->Get("histo_WZ_JESUp"));
  addTemplate("WZ_WE_JESDown"       , vars, wspace, (TH1F*)WenGfile->Get("histo_WZ_JESDown"));
  addTemplate("WZ_WE_PESUp"         , vars, wspace, (TH1F*)WenGfile->Get("histo_WZ_PESUp"));
  addTemplate("WZ_WE_PESDown"       , vars, wspace, (TH1F*)WenGfile->Get("histo_WZ_PESDown"));
  addTemplate("WW_WE"               , vars, wspace, (TH1F*)WenGfile->Get("histo_WW"));
  addTemplate("WW_WE_JESUp"         , vars, wspace, (TH1F*)WenGfile->Get("histo_WW_JESUp"));
  addTemplate("WW_WE_JESDown"       , vars, wspace, (TH1F*)WenGfile->Get("histo_WW_JESDown"));
  addTemplate("WW_WE_PESUp"         , vars, wspace, (TH1F*)WenGfile->Get("histo_WW_PESUp"));
  addTemplate("WW_WE_PESDown"       , vars, wspace, (TH1F*)WenGfile->Get("histo_WW_PESDown"));
  
  //Statistical errors
  for(int i = 1; i <= nBins; i++){
    char binChar[10];
    sprintf(binChar, "%d", i);
    std::string binNumber(binChar);
    
    
// Primary
    // TH1F* ZNuNuG_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG"))->Clone("ZNuNuG_SR_above0p5_histBinUp"); // Only if not fitting
    // TH1F* ZNuNuG_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_ZNuNuG"))->Clone("ZNuNuG_SR_above0p5_histBinDown"); // Only if not fitting
    // ZNuNuG_SR_above0p5_histBinUp->SetBinContent(i, ZNuNuG_SR_above0p5_histBinUp->GetBinContent(i)+ZNuNuG_SR_above0p5_histBinUp->GetBinError(i)); // Only if not fitting
    // ZNuNuG_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*ZNuNuG_SR_above0p5_histBinDown->GetBinContent(i), ZNuNuG_SR_above0p5_histBinDown->GetBinContent(i)-ZNuNuG_SR_above0p5_histBinDown->GetBinError(i))); // Only if not fitting
    // addTemplate("ZNuNuG_SR_above0p5_ZNuNuGSignalSBin"+binNumber+"Up", vars, wspace, ZNuNuG_SR_above0p5_histBinUp); // Only if not fitting
    // addTemplate("ZNuNuG_SR_above0p5_ZNuNuGSignalSBin"+binNumber+"Down", vars, wspace, ZNuNuG_SR_above0p5_histBinDown); // Only if not fitting
    // TH1F* WG_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WG"))->Clone("WG_SR_above0p5_histBinUp"); // Only if not fitting
    // TH1F* WG_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WG"))->Clone("WG_SR_above0p5_histBinDown"); // Only if not fitting
    // WG_SR_above0p5_histBinUp->SetBinContent(i, WG_SR_above0p5_histBinUp->GetBinContent(i)+WG_SR_above0p5_histBinUp->GetBinError(i)); // Only if not fitting
    // WG_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WG_SR_above0p5_histBinDown->GetBinContent(i), WG_SR_above0p5_histBinDown->GetBinContent(i)-WG_SR_above0p5_histBinDown->GetBinError(i))); // Only if not fitting
    // addTemplate("WG_SR_above0p5_WGSignalSBin"+binNumber+"Up", vars, wspace, WG_SR_above0p5_histBinUp); // Only if not fitting
    // addTemplate("WG_SR_above0p5_WGSignalSBin"+binNumber+"Down", vars, wspace, WG_SR_above0p5_histBinDown); // Only if not fitting
    // TH1F* ZNuNuG_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG"))->Clone("ZNuNuG_SR_below0p5_histBinUp"); // Only if not fitting
    // TH1F* ZNuNuG_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_ZNuNuG"))->Clone("ZNuNuG_SR_below0p5_histBinDown"); // Only if not fitting
    // ZNuNuG_SR_below0p5_histBinUp->SetBinContent(i, ZNuNuG_SR_below0p5_histBinUp->GetBinContent(i)+ZNuNuG_SR_below0p5_histBinUp->GetBinError(i)); // Only if not fitting
    // ZNuNuG_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*ZNuNuG_SR_below0p5_histBinDown->GetBinContent(i), ZNuNuG_SR_below0p5_histBinDown->GetBinContent(i)-ZNuNuG_SR_below0p5_histBinDown->GetBinError(i))); // Only if not fitting
    // addTemplate("ZNuNuG_SR_below0p5_ZNuNuGSignalSBin"+binNumber+"Up", vars, wspace, ZNuNuG_SR_below0p5_histBinUp); // Only if not fitting
    // addTemplate("ZNuNuG_SR_below0p5_ZNuNuGSignalSBin"+binNumber+"Down", vars, wspace, ZNuNuG_SR_below0p5_histBinDown); // Only if not fitting
    // TH1F* WG_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WG"))->Clone("WG_SR_below0p5_histBinUp"); // Only if not fitting
    // TH1F* WG_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WG"))->Clone("WG_SR_below0p5_histBinDown"); // Only if not fitting
    // WG_SR_below0p5_histBinUp->SetBinContent(i, WG_SR_below0p5_histBinUp->GetBinContent(i)+WG_SR_below0p5_histBinUp->GetBinError(i)); // Only if not fitting
    // WG_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WG_SR_below0p5_histBinDown->GetBinContent(i), WG_SR_below0p5_histBinDown->GetBinContent(i)-WG_SR_below0p5_histBinDown->GetBinError(i))); // Only if not fitting
    // addTemplate("WG_SR_below0p5_WGSignalSBin"+binNumber+"Up", vars, wspace, WG_SR_below0p5_histBinUp); // Only if not fitting
    // addTemplate("WG_SR_below0p5_WGSignalSBin"+binNumber+"Down", vars, wspace, WG_SR_below0p5_histBinDown); // Only if not fitting
    
    TH1F* GJets_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_GJets"))->Clone("GJets_SR_above0p5_histBinUp");
    TH1F* GJets_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_GJets"))->Clone("GJets_SR_above0p5_histBinDown");
    GJets_SR_above0p5_histBinUp->SetBinContent(i, GJets_SR_above0p5_histBinUp->GetBinContent(i)+GJets_SR_above0p5_histBinUp->GetBinError(i));
    GJets_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*GJets_SR_above0p5_histBinDown->GetBinContent(i), GJets_SR_above0p5_histBinDown->GetBinContent(i)-GJets_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("GJets_SR_above0p5_GJetsSignalSBin"+binNumber+"Up", vars, wspace, GJets_SR_above0p5_histBinUp);
    addTemplate("GJets_SR_above0p5_GJetsSignalSBin"+binNumber+"Down", vars, wspace, GJets_SR_above0p5_histBinDown);
    TH1F* WZ_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WZ"))->Clone("WZ_SR_above0p5_histBinUp");
    TH1F* WZ_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WZ"))->Clone("WZ_SR_above0p5_histBinDown");
    WZ_SR_above0p5_histBinUp->SetBinContent(i, WZ_SR_above0p5_histBinUp->GetBinContent(i)+WZ_SR_above0p5_histBinUp->GetBinError(i));
    WZ_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WZ_SR_above0p5_histBinDown->GetBinContent(i), WZ_SR_above0p5_histBinDown->GetBinContent(i)-WZ_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("WZ_SR_above0p5_WZSignalSBin"+binNumber+"Up", vars, wspace, WZ_SR_above0p5_histBinUp);
    addTemplate("WZ_SR_above0p5_WZSignalSBin"+binNumber+"Down", vars, wspace, WZ_SR_above0p5_histBinDown);
    TH1F* WMuNu_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WMuNu"))->Clone("WMuNu_SR_above0p5_histBinUp");
    TH1F* WMuNu_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WMuNu"))->Clone("WMuNu_SR_above0p5_histBinDown");
    WMuNu_SR_above0p5_histBinUp->SetBinContent(i, WMuNu_SR_above0p5_histBinUp->GetBinContent(i)+WMuNu_SR_above0p5_histBinUp->GetBinError(i));
    WMuNu_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WMuNu_SR_above0p5_histBinDown->GetBinContent(i), WMuNu_SR_above0p5_histBinDown->GetBinContent(i)-WMuNu_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("WMuNu_SR_above0p5_WMuNuSignalSBin"+binNumber+"Up", vars, wspace, WMuNu_SR_above0p5_histBinUp);
    addTemplate("WMuNu_SR_above0p5_WMuNuSignalSBin"+binNumber+"Down", vars, wspace, WMuNu_SR_above0p5_histBinDown);
    TH1F* WTauNu_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WTauNu"))->Clone("WTauNu_SR_above0p5_histBinUp");
    TH1F* WTauNu_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WTauNu"))->Clone("WTauNu_SR_above0p5_histBinDown");
    WTauNu_SR_above0p5_histBinUp->SetBinContent(i, WTauNu_SR_above0p5_histBinUp->GetBinContent(i)+WTauNu_SR_above0p5_histBinUp->GetBinError(i));
    WTauNu_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WTauNu_SR_above0p5_histBinDown->GetBinContent(i), WTauNu_SR_above0p5_histBinDown->GetBinContent(i)-WTauNu_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("WTauNu_SR_above0p5_WTauNuSignalSBin"+binNumber+"Up", vars, wspace, WTauNu_SR_above0p5_histBinUp);
    addTemplate("WTauNu_SR_above0p5_WTauNuSignalSBin"+binNumber+"Down", vars, wspace, WTauNu_SR_above0p5_histBinDown);
    TH1F* GJets_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_GJets"))->Clone("GJets_SR_below0p5_histBinUp");
    TH1F* GJets_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_GJets"))->Clone("GJets_SR_below0p5_histBinDown");
    GJets_SR_below0p5_histBinUp->SetBinContent(i, GJets_SR_below0p5_histBinUp->GetBinContent(i)+GJets_SR_below0p5_histBinUp->GetBinError(i));
    GJets_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*GJets_SR_below0p5_histBinDown->GetBinContent(i), GJets_SR_below0p5_histBinDown->GetBinContent(i)-GJets_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("GJets_SR_below0p5_GJetsSignalSBin"+binNumber+"Up", vars, wspace, GJets_SR_below0p5_histBinUp);
    addTemplate("GJets_SR_below0p5_GJetsSignalSBin"+binNumber+"Down", vars, wspace, GJets_SR_below0p5_histBinDown);
    TH1F* WZ_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WZ"))->Clone("WZ_SR_below0p5_histBinUp");
    TH1F* WZ_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WZ"))->Clone("WZ_SR_below0p5_histBinDown");
    WZ_SR_below0p5_histBinUp->SetBinContent(i, WZ_SR_below0p5_histBinUp->GetBinContent(i)+WZ_SR_below0p5_histBinUp->GetBinError(i));
    WZ_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WZ_SR_below0p5_histBinDown->GetBinContent(i), WZ_SR_below0p5_histBinDown->GetBinContent(i)-WZ_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("WZ_SR_below0p5_WZSignalSBin"+binNumber+"Up", vars, wspace, WZ_SR_below0p5_histBinUp);
    addTemplate("WZ_SR_below0p5_WZSignalSBin"+binNumber+"Down", vars, wspace, WZ_SR_below0p5_histBinDown);
    TH1F* WMuNu_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu"))->Clone("WMuNu_SR_below0p5_histBinUp");
    TH1F* WMuNu_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WMuNu"))->Clone("WMuNu_SR_below0p5_histBinDown");
    WMuNu_SR_below0p5_histBinUp->SetBinContent(i, WMuNu_SR_below0p5_histBinUp->GetBinContent(i)+WMuNu_SR_below0p5_histBinUp->GetBinError(i));
    WMuNu_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WMuNu_SR_below0p5_histBinDown->GetBinContent(i), WMuNu_SR_below0p5_histBinDown->GetBinContent(i)-WMuNu_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("WMuNu_SR_below0p5_WMuNuSignalSBin"+binNumber+"Up", vars, wspace, WMuNu_SR_below0p5_histBinUp);
    addTemplate("WMuNu_SR_below0p5_WMuNuSignalSBin"+binNumber+"Down", vars, wspace, WMuNu_SR_below0p5_histBinDown);
    TH1F* WTauNu_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu"))->Clone("WTauNu_SR_below0p5_histBinUp");
    TH1F* WTauNu_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WTauNu"))->Clone("WTauNu_SR_below0p5_histBinDown");
    WTauNu_SR_below0p5_histBinUp->SetBinContent(i, WTauNu_SR_below0p5_histBinUp->GetBinContent(i)+WTauNu_SR_below0p5_histBinUp->GetBinError(i));
    WTauNu_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WTauNu_SR_below0p5_histBinDown->GetBinContent(i), WTauNu_SR_below0p5_histBinDown->GetBinContent(i)-WTauNu_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("WTauNu_SR_below0p5_WTauNuSignalSBin"+binNumber+"Up", vars, wspace, WTauNu_SR_below0p5_histBinUp);
    addTemplate("WTauNu_SR_below0p5_WTauNuSignalSBin"+binNumber+"Down", vars, wspace, WTauNu_SR_below0p5_histBinDown);
    TH1F* WZ_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_WZ"))->Clone("WZ_WM_histBinUp");
    TH1F* WZ_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_WZ"))->Clone("WZ_WM_histBinDown");
    WZ_WM_histBinUp->SetBinContent(i, WZ_WM_histBinUp->GetBinContent(i)+WZ_WM_histBinUp->GetBinError(i));
    WZ_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*WZ_WM_histBinDown->GetBinContent(i), WZ_WM_histBinDown->GetBinContent(i)-WZ_WM_histBinDown->GetBinError(i)));
    addTemplate("WZ_WM_WZMonomuSBin"+binNumber+"Up", vars, wspace, WZ_WM_histBinUp);
    addTemplate("WZ_WM_WZMonomuSBin"+binNumber+"Down", vars, wspace, WZ_WM_histBinDown);
    TH1F* WW_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_WW"))->Clone("WW_WM_histBinUp");
    TH1F* WW_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_WW"))->Clone("WW_WM_histBinDown");
    WW_WM_histBinUp->SetBinContent(i, WW_WM_histBinUp->GetBinContent(i)+WW_WM_histBinUp->GetBinError(i));
    WW_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*WW_WM_histBinDown->GetBinContent(i), WW_WM_histBinDown->GetBinContent(i)-WW_WM_histBinDown->GetBinError(i)));
    addTemplate("WW_WM_WWMonomuSBin"+binNumber+"Up", vars, wspace, WW_WM_histBinUp);
    addTemplate("WW_WM_WWMonomuSBin"+binNumber+"Down", vars, wspace, WW_WM_histBinDown);
    TH1F* WZ_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_WZ"))->Clone("WZ_WE_histBinUp");
    TH1F* WZ_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_WZ"))->Clone("WZ_WE_histBinDown");
    WZ_WE_histBinUp->SetBinContent(i, WZ_WE_histBinUp->GetBinContent(i)+WZ_WE_histBinUp->GetBinError(i));
    WZ_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*WZ_WE_histBinDown->GetBinContent(i), WZ_WE_histBinDown->GetBinContent(i)-WZ_WE_histBinDown->GetBinError(i)));
    addTemplate("WZ_WE_WZMonoeleSBin"+binNumber+"Up", vars, wspace, WZ_WE_histBinUp);
    addTemplate("WZ_WE_WZMonoeleSBin"+binNumber+"Down", vars, wspace, WZ_WE_histBinDown);
    TH1F* WW_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_WW"))->Clone("WW_WE_histBinUp");
    TH1F* WW_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_WW"))->Clone("WW_WE_histBinDown");
    WW_WE_histBinUp->SetBinContent(i, WW_WE_histBinUp->GetBinContent(i)+WW_WE_histBinUp->GetBinError(i));
    WW_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*WW_WE_histBinDown->GetBinContent(i), WW_WE_histBinDown->GetBinContent(i)-WW_WE_histBinDown->GetBinError(i)));
    addTemplate("WW_WE_WWMonoeleSBin"+binNumber+"Up", vars, wspace, WW_WE_histBinUp);
    addTemplate("WW_WE_WWMonoeleSBin"+binNumber+"Down", vars, wspace, WW_WE_histBinDown);

// Others    
    // TH1F* signal_SR_above0p5_histBinUp = (TH1F*)signal_SR_above0p5_hist->Clone("signal_SR_above0p5_histBinUp");
    // TH1F* signal_SR_above0p5_histBinDown = (TH1F*)signal_SR_above0p5_hist->Clone("signal_SR_above0p5_histBinDown");
    // signal_SR_above0p5_histBinUp->SetBinContent(i, signal_SR_above0p5_hist->GetBinContent(i)+signal_SR_above0p5_hist->GetBinError(i));
    // signal_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*signal_SR_above0p5_hist->GetBinContent(i), signal_SR_above0p5_hist->GetBinContent(i)-signal_SR_above0p5_hist->GetBinError(i)));
    // addTemplate("Signal_SR_above0p5_DMSignalSBin"+binNumber+"Up", vars, wspace, signal_SR_above0p5_histBinUp);
    // addTemplate("Signal_SR_above0p5_DMSignalSBin"+binNumber+"Down", vars, wspace, signal_SR_above0p5_histBinDown);
    
    TH1F* QCD_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_jetfake"))->Clone("QCD_SR_above0p5_histBinUp");
    TH1F* QCD_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_jetfake"))->Clone("QCD_SR_above0p5_histBinDown");
    QCD_SR_above0p5_histBinUp->SetBinContent(i, QCD_SR_above0p5_histBinUp->GetBinContent(i)+QCD_SR_above0p5_histBinUp->GetBinError(i));
    QCD_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*QCD_SR_above0p5_histBinDown->GetBinContent(i), QCD_SR_above0p5_histBinDown->GetBinContent(i)-QCD_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("QCD_SR_above0p5_QCDHiPhiSignalSBin"+binNumber+"Up", vars, wspace, QCD_SR_above0p5_histBinUp);
    addTemplate("QCD_SR_above0p5_QCDHiPhiSignalSBin"+binNumber+"Down", vars, wspace, QCD_SR_above0p5_histBinDown);
    // TH1F* Elefake_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_elefake"))->Clone("Elefake_SR_above0p5_histBinUp");
    // TH1F* Elefake_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_elefake"))->Clone("Elefake_SR_above0p5_histBinDown");
    // Elefake_SR_above0p5_histBinUp->SetBinContent(i, Elefake_SR_above0p5_histBinUp->GetBinContent(i)+Elefake_SR_above0p5_histBinUp->GetBinError(i));
    // Elefake_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*Elefake_SR_above0p5_histBinDown->GetBinContent(i), Elefake_SR_above0p5_histBinDown->GetBinContent(i)-Elefake_SR_above0p5_histBinDown->GetBinError(i)));
    // addTemplate("Elefake_SR_above0p5_EleHiPhiSignalSBin"+binNumber+"Up", vars, wspace, Elefake_SR_above0p5_histBinUp);
    // addTemplate("Elefake_SR_above0p5_EleHiPhiSignalSBin"+binNumber+"Down", vars, wspace, Elefake_SR_above0p5_histBinDown);
    // TH1F* ZllG_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_ZllG_combined"))->Clone("ZllG_SR_above0p5_histBinUp");
    // TH1F* ZllG_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_ZllG_combined"))->Clone("ZllG_SR_above0p5_histBinDown");
    // ZllG_SR_above0p5_histBinUp->SetBinContent(i, ZllG_SR_above0p5_histBinUp->GetBinContent(i)+ZllG_SR_above0p5_histBinUp->GetBinError(i));
    // ZllG_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*ZllG_SR_above0p5_histBinDown->GetBinContent(i), ZllG_SR_above0p5_histBinDown->GetBinContent(i)-ZllG_SR_above0p5_histBinDown->GetBinError(i)));
    // addTemplate("ZllG_SR_above0p5_ZllGSignalSBin"+binNumber+"Up", vars, wspace, ZllG_SR_above0p5_histBinUp);
    // addTemplate("ZllG_SR_above0p5_ZllGSignalSBin"+binNumber+"Down", vars, wspace, ZllG_SR_above0p5_histBinDown);
    // TH1F* TTG_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_TTG"))->Clone("TTG_SR_above0p5_histBinUp");
    // TH1F* TTG_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_TTG"))->Clone("TTG_SR_above0p5_histBinDown");
    // TTG_SR_above0p5_histBinUp->SetBinContent(i, TTG_SR_above0p5_histBinUp->GetBinContent(i)+TTG_SR_above0p5_histBinUp->GetBinError(i));
    // TTG_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*TTG_SR_above0p5_histBinDown->GetBinContent(i), TTG_SR_above0p5_histBinDown->GetBinContent(i)-TTG_SR_above0p5_histBinDown->GetBinError(i)));
    // addTemplate("TTG_SR_above0p5_TTGSignalSBin"+binNumber+"Up", vars, wspace, TTG_SR_above0p5_histBinUp);
    // addTemplate("TTG_SR_above0p5_TTGSignalSBin"+binNumber+"Down", vars, wspace, TTG_SR_above0p5_histBinDown);
    // TH1F* TG_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_TG"))->Clone("TG_SR_above0p5_histBinUp");
    // TH1F* TG_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_TG"))->Clone("TG_SR_above0p5_histBinDown");
    // TG_SR_above0p5_histBinUp->SetBinContent(i, TG_SR_above0p5_histBinUp->GetBinContent(i)+TG_SR_above0p5_histBinUp->GetBinError(i));
    // TG_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*TG_SR_above0p5_histBinDown->GetBinContent(i), TG_SR_above0p5_histBinDown->GetBinContent(i)-TG_SR_above0p5_histBinDown->GetBinError(i)));
    // addTemplate("TG_SR_above0p5_TGSignalSBin"+binNumber+"Up", vars, wspace, TG_SR_above0p5_histBinUp);
    // addTemplate("TG_SR_above0p5_TGSignalSBin"+binNumber+"Down", vars, wspace, TG_SR_above0p5_histBinDown);
    // TH1F* Diphoton_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_diphoton"))->Clone("Diphoton_SR_above0p5_histBinUp");
    // TH1F* Diphoton_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_diphoton"))->Clone("Diphoton_SR_above0p5_histBinDown");
    // Diphoton_SR_above0p5_histBinUp->SetBinContent(i, Diphoton_SR_above0p5_histBinUp->GetBinContent(i)+Diphoton_SR_above0p5_histBinUp->GetBinError(i));
    // Diphoton_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*Diphoton_SR_above0p5_histBinDown->GetBinContent(i), Diphoton_SR_above0p5_histBinDown->GetBinContent(i)-Diphoton_SR_above0p5_histBinDown->GetBinError(i)));
    // addTemplate("Diphoton_SR_above0p5_DiphotonSignalSBin"+binNumber+"Up", vars, wspace, Diphoton_SR_above0p5_histBinUp);
    // addTemplate("Diphoton_SR_above0p5_DiphotonSignalSBin"+binNumber+"Down", vars, wspace, Diphoton_SR_above0p5_histBinDown);
    TH1F* ZZ_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_ZZ"))->Clone("ZZ_SR_above0p5_histBinUp");
    TH1F* ZZ_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_ZZ"))->Clone("ZZ_SR_above0p5_histBinDown");
    ZZ_SR_above0p5_histBinUp->SetBinContent(i, ZZ_SR_above0p5_histBinUp->GetBinContent(i)+ZZ_SR_above0p5_histBinUp->GetBinError(i));
    ZZ_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*ZZ_SR_above0p5_histBinDown->GetBinContent(i), ZZ_SR_above0p5_histBinDown->GetBinContent(i)-ZZ_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("ZZ_SR_above0p5_ZZSignalSBin"+binNumber+"Up", vars, wspace, ZZ_SR_above0p5_histBinUp);
    addTemplate("ZZ_SR_above0p5_ZZSignalSBin"+binNumber+"Down", vars, wspace, ZZ_SR_above0p5_histBinDown);
    TH1F* WW_SR_above0p5_histBinUp = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WW"))->Clone("WW_SR_above0p5_histBinUp");
    TH1F* WW_SR_above0p5_histBinDown = (TH1F*)((TH1F*)ZnnGabove0p5file->Get("histo_WW"))->Clone("WW_SR_above0p5_histBinDown");
    WW_SR_above0p5_histBinUp->SetBinContent(i, WW_SR_above0p5_histBinUp->GetBinContent(i)+WW_SR_above0p5_histBinUp->GetBinError(i));
    WW_SR_above0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WW_SR_above0p5_histBinDown->GetBinContent(i), WW_SR_above0p5_histBinDown->GetBinContent(i)-WW_SR_above0p5_histBinDown->GetBinError(i)));
    addTemplate("WW_SR_above0p5_WWSignalSBin"+binNumber+"Up", vars, wspace, WW_SR_above0p5_histBinUp);
    addTemplate("WW_SR_above0p5_WWSignalSBin"+binNumber+"Down", vars, wspace, WW_SR_above0p5_histBinDown);
    
    // TH1F* signal_SR_below0p5_histBinUp = (TH1F*)signal_SR_below0p5_hist->Clone("signal_SR_below0p5_histBinUp");
    // TH1F* signal_SR_below0p5_histBinDown = (TH1F*)signal_SR_below0p5_hist->Clone("signal_SR_below0p5_histBinDown");
    // signal_SR_below0p5_histBinUp->SetBinContent(i, signal_SR_below0p5_hist->GetBinContent(i)+signal_SR_below0p5_hist->GetBinError(i));
    // signal_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*signal_SR_below0p5_hist->GetBinContent(i), signal_SR_below0p5_hist->GetBinContent(i)-signal_SR_below0p5_hist->GetBinError(i)));
    // addTemplate("Signal_SR_below0p5_DMSignalSBin"+binNumber+"Up", vars, wspace, signal_SR_below0p5_histBinUp);
    // addTemplate("Signal_SR_below0p5_DMSignalSBin"+binNumber+"Down", vars, wspace, signal_SR_below0p5_histBinDown);
    
    TH1F* QCD_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_jetfake"))->Clone("QCD_SR_below0p5_histBinUp");
    TH1F* QCD_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_jetfake"))->Clone("QCD_SR_below0p5_histBinDown");
    QCD_SR_below0p5_histBinUp->SetBinContent(i, QCD_SR_below0p5_histBinUp->GetBinContent(i)+QCD_SR_below0p5_histBinUp->GetBinError(i));
    QCD_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*QCD_SR_below0p5_histBinDown->GetBinContent(i), QCD_SR_below0p5_histBinDown->GetBinContent(i)-QCD_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("QCD_SR_below0p5_QCDLoPhiSignalSBin"+binNumber+"Up", vars, wspace, QCD_SR_below0p5_histBinUp);
    addTemplate("QCD_SR_below0p5_QCDLoPhiSignalSBin"+binNumber+"Down", vars, wspace, QCD_SR_below0p5_histBinDown);
    // TH1F* Elefake_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_elefake"))->Clone("Elefake_SR_below0p5_histBinUp");
    // TH1F* Elefake_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_elefake"))->Clone("Elefake_SR_below0p5_histBinDown");
    // Elefake_SR_below0p5_histBinUp->SetBinContent(i, Elefake_SR_below0p5_histBinUp->GetBinContent(i)+Elefake_SR_below0p5_histBinUp->GetBinError(i));
    // Elefake_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*Elefake_SR_below0p5_histBinDown->GetBinContent(i), Elefake_SR_below0p5_histBinDown->GetBinContent(i)-Elefake_SR_below0p5_histBinDown->GetBinError(i)));
    // addTemplate("Elefake_SR_below0p5_EleLoPhiSignalSBin"+binNumber+"Up", vars, wspace, Elefake_SR_below0p5_histBinUp);
    // addTemplate("Elefake_SR_below0p5_EleLoPhiSignalSBin"+binNumber+"Down", vars, wspace, Elefake_SR_below0p5_histBinDown);
    // TH1F* ZllG_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_combined"))->Clone("ZllG_SR_below0p5_histBinUp");
    // TH1F* ZllG_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_ZllG_combined"))->Clone("ZllG_SR_below0p5_histBinDown");
    // ZllG_SR_below0p5_histBinUp->SetBinContent(i, ZllG_SR_below0p5_histBinUp->GetBinContent(i)+ZllG_SR_below0p5_histBinUp->GetBinError(i));
    // ZllG_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*ZllG_SR_below0p5_histBinDown->GetBinContent(i), ZllG_SR_below0p5_histBinDown->GetBinContent(i)-ZllG_SR_below0p5_histBinDown->GetBinError(i)));
    // addTemplate("ZllG_SR_below0p5_ZllGSignalSBin"+binNumber+"Up", vars, wspace, ZllG_SR_below0p5_histBinUp);
    // addTemplate("ZllG_SR_below0p5_ZllGSignalSBin"+binNumber+"Down", vars, wspace, ZllG_SR_below0p5_histBinDown);
    // TH1F* TTG_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_TTG"))->Clone("TTG_SR_below0p5_histBinUp");
    // TH1F* TTG_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_TTG"))->Clone("TTG_SR_below0p5_histBinDown");
    // TTG_SR_below0p5_histBinUp->SetBinContent(i, TTG_SR_below0p5_histBinUp->GetBinContent(i)+TTG_SR_below0p5_histBinUp->GetBinError(i));
    // TTG_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*TTG_SR_below0p5_histBinDown->GetBinContent(i), TTG_SR_below0p5_histBinDown->GetBinContent(i)-TTG_SR_below0p5_histBinDown->GetBinError(i)));
    // addTemplate("TTG_SR_below0p5_TTGSignalSBin"+binNumber+"Up", vars, wspace, TTG_SR_below0p5_histBinUp);
    // addTemplate("TTG_SR_below0p5_TTGSignalSBin"+binNumber+"Down", vars, wspace, TTG_SR_below0p5_histBinDown);
    // TH1F* TG_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_TG"))->Clone("TG_SR_below0p5_histBinUp");
    // TH1F* TG_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_TG"))->Clone("TG_SR_below0p5_histBinDown");
    // TG_SR_below0p5_histBinUp->SetBinContent(i, TG_SR_below0p5_histBinUp->GetBinContent(i)+TG_SR_below0p5_histBinUp->GetBinError(i));
    // TG_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*TG_SR_below0p5_histBinDown->GetBinContent(i), TG_SR_below0p5_histBinDown->GetBinContent(i)-TG_SR_below0p5_histBinDown->GetBinError(i)));
    // addTemplate("TG_SR_below0p5_TGSignalSBin"+binNumber+"Up", vars, wspace, TG_SR_below0p5_histBinUp);
    // addTemplate("TG_SR_below0p5_TGSignalSBin"+binNumber+"Down", vars, wspace, TG_SR_below0p5_histBinDown);
    // TH1F* Diphoton_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_diphoton"))->Clone("Diphoton_SR_below0p5_histBinUp");
    // TH1F* Diphoton_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_diphoton"))->Clone("Diphoton_SR_below0p5_histBinDown");
    // Diphoton_SR_below0p5_histBinUp->SetBinContent(i, Diphoton_SR_below0p5_histBinUp->GetBinContent(i)+Diphoton_SR_below0p5_histBinUp->GetBinError(i));
    // Diphoton_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*Diphoton_SR_below0p5_histBinDown->GetBinContent(i), Diphoton_SR_below0p5_histBinDown->GetBinContent(i)-Diphoton_SR_below0p5_histBinDown->GetBinError(i)));
    // addTemplate("Diphoton_SR_below0p5_DiphotonSignalSBin"+binNumber+"Up", vars, wspace, Diphoton_SR_below0p5_histBinUp);
    // addTemplate("Diphoton_SR_below0p5_DiphotonSignalSBin"+binNumber+"Down", vars, wspace, Diphoton_SR_below0p5_histBinDown);
    TH1F* ZZ_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_ZZ"))->Clone("ZZ_SR_below0p5_histBinUp");
    TH1F* ZZ_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_ZZ"))->Clone("ZZ_SR_below0p5_histBinDown");
    ZZ_SR_below0p5_histBinUp->SetBinContent(i, ZZ_SR_below0p5_histBinUp->GetBinContent(i)+ZZ_SR_below0p5_histBinUp->GetBinError(i));
    ZZ_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*ZZ_SR_below0p5_histBinDown->GetBinContent(i), ZZ_SR_below0p5_histBinDown->GetBinContent(i)-ZZ_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("ZZ_SR_below0p5_ZZSignalSBin"+binNumber+"Up", vars, wspace, ZZ_SR_below0p5_histBinUp);
    addTemplate("ZZ_SR_below0p5_ZZSignalSBin"+binNumber+"Down", vars, wspace, ZZ_SR_below0p5_histBinDown);
    TH1F* WW_SR_below0p5_histBinUp = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WW"))->Clone("WW_SR_below0p5_histBinUp");
    TH1F* WW_SR_below0p5_histBinDown = (TH1F*)((TH1F*)ZnnGbelow0p5file->Get("histo_WW"))->Clone("WW_SR_below0p5_histBinDown");
    WW_SR_below0p5_histBinUp->SetBinContent(i, WW_SR_below0p5_histBinUp->GetBinContent(i)+WW_SR_below0p5_histBinUp->GetBinError(i));
    WW_SR_below0p5_histBinDown->SetBinContent(i, TMath::Max(0.01*WW_SR_below0p5_histBinDown->GetBinContent(i), WW_SR_below0p5_histBinDown->GetBinContent(i)-WW_SR_below0p5_histBinDown->GetBinError(i)));
    addTemplate("WW_SR_below0p5_WWSignalSBin"+binNumber+"Up", vars, wspace, WW_SR_below0p5_histBinUp);
    addTemplate("WW_SR_below0p5_WWSignalSBin"+binNumber+"Down", vars, wspace, WW_SR_below0p5_histBinDown);
    
    if(!aTGC){
      TH1F* QCD_ZM_histBinUp = (TH1F*)((TH1F*)ZmmGfile->Get("histo_jetfake"))->Clone("QCD_ZM_histBinUp");
      TH1F* QCD_ZM_histBinDown = (TH1F*)((TH1F*)ZmmGfile->Get("histo_jetfake"))->Clone("QCD_ZM_histBinDown");
      QCD_ZM_histBinUp->SetBinContent(i, QCD_ZM_histBinUp->GetBinContent(i)+QCD_ZM_histBinUp->GetBinError(i));
      QCD_ZM_histBinDown->SetBinContent(i, TMath::Max(0.01*QCD_ZM_histBinDown->GetBinContent(i), QCD_ZM_histBinDown->GetBinContent(i)-QCD_ZM_histBinDown->GetBinError(i)));
      addTemplate("QCD_ZM_QCDDimuSBin"+binNumber+"Up", vars, wspace, QCD_ZM_histBinUp);
      addTemplate("QCD_ZM_QCDDimuSBin"+binNumber+"Down", vars, wspace, QCD_ZM_histBinDown);
      // TH1F* TTG_ZM_histBinUp = (TH1F*)((TH1F*)ZmmGfile->Get("histo_TTG"))->Clone("TTG_ZM_histBinUp");
      // TH1F* TTG_ZM_histBinDown = (TH1F*)((TH1F*)ZmmGfile->Get("histo_TTG"))->Clone("TTG_ZM_histBinDown");
      // TTG_ZM_histBinUp->SetBinContent(i, TTG_ZM_histBinUp->GetBinContent(i)+TTG_ZM_histBinUp->GetBinError(i));
      // TTG_ZM_histBinDown->SetBinContent(i, TMath::Max(0.01*TTG_ZM_histBinDown->GetBinContent(i), TTG_ZM_histBinDown->GetBinContent(i)-TTG_ZM_histBinDown->GetBinError(i)));
      // addTemplate("TTG_ZM_TTGDimuSBin"+binNumber+"Up", vars, wspace, TTG_ZM_histBinUp);
      // addTemplate("TTG_ZM_TTGDimuSBin"+binNumber+"Down", vars, wspace, TTG_ZM_histBinDown);
      // TH1F* WZ_ZM_histBinUp = (TH1F*)((TH1F*)ZmmGfile->Get("histo_WZ"))->Clone("WZ_ZM_histBinUp");
      // TH1F* WZ_ZM_histBinDown = (TH1F*)((TH1F*)ZmmGfile->Get("histo_WZ"))->Clone("WZ_ZM_histBinDown");
      // WZ_ZM_histBinUp->SetBinContent(i, WZ_ZM_histBinUp->GetBinContent(i)+WZ_ZM_histBinUp->GetBinError(i));
      // WZ_ZM_histBinDown->SetBinContent(i, TMath::Max(0.01*WZ_ZM_histBinDown->GetBinContent(i), WZ_ZM_histBinDown->GetBinContent(i)-WZ_ZM_histBinDown->GetBinError(i)));
      // addTemplate("WZ_ZM_WZDimuSBin"+binNumber+"Up", vars, wspace, WZ_ZM_histBinUp);
      // addTemplate("WZ_ZM_WZDimuSBin"+binNumber+"Down", vars, wspace, WZ_ZM_histBinDown);
    }
    
    if(!aTGC){
      TH1F* QCD_ZE_histBinUp = (TH1F*)((TH1F*)ZeeGfile->Get("histo_jetfake"))->Clone("QCD_ZE_histBinUp");
      TH1F* QCD_ZE_histBinDown = (TH1F*)((TH1F*)ZeeGfile->Get("histo_jetfake"))->Clone("QCD_ZE_histBinDown");
      QCD_ZE_histBinUp->SetBinContent(i, QCD_ZE_histBinUp->GetBinContent(i)+QCD_ZE_histBinUp->GetBinError(i));
      QCD_ZE_histBinDown->SetBinContent(i, TMath::Max(0.01*QCD_ZE_histBinDown->GetBinContent(i), QCD_ZE_histBinDown->GetBinContent(i)-QCD_ZE_histBinDown->GetBinError(i)));
      addTemplate("QCD_ZE_QCDDieleSBin"+binNumber+"Up", vars, wspace, QCD_ZE_histBinUp);
      addTemplate("QCD_ZE_QCDDieleSBin"+binNumber+"Down", vars, wspace, QCD_ZE_histBinDown);
      // TH1F* TTG_ZE_histBinUp = (TH1F*)((TH1F*)ZeeGfile->Get("histo_TTG"))->Clone("TTG_ZE_histBinUp");
      // TH1F* TTG_ZE_histBinDown = (TH1F*)((TH1F*)ZeeGfile->Get("histo_TTG"))->Clone("TTG_ZE_histBinDown");
      // TTG_ZE_histBinUp->SetBinContent(i, TTG_ZE_histBinUp->GetBinContent(i)+TTG_ZE_histBinUp->GetBinError(i));
      // TTG_ZE_histBinDown->SetBinContent(i, TMath::Max(0.01*TTG_ZE_histBinDown->GetBinContent(i), TTG_ZE_histBinDown->GetBinContent(i)-TTG_ZE_histBinDown->GetBinError(i)));
      // addTemplate("TTG_ZE_TTGDieleSBin"+binNumber+"Up", vars, wspace, TTG_ZE_histBinUp);
      // addTemplate("TTG_ZE_TTGDieleSBin"+binNumber+"Down", vars, wspace, TTG_ZE_histBinDown);
      // TH1F* WZ_ZE_histBinUp = (TH1F*)((TH1F*)ZeeGfile->Get("histo_WZ"))->Clone("WZ_ZE_histBinUp");
      // TH1F* WZ_ZE_histBinDown = (TH1F*)((TH1F*)ZeeGfile->Get("histo_WZ"))->Clone("WZ_ZE_histBinDown");
      // WZ_ZE_histBinUp->SetBinContent(i, WZ_ZE_histBinUp->GetBinContent(i)+WZ_ZE_histBinUp->GetBinError(i));
      // WZ_ZE_histBinDown->SetBinContent(i, TMath::Max(0.01*WZ_ZE_histBinDown->GetBinContent(i), WZ_ZE_histBinDown->GetBinContent(i)-WZ_ZE_histBinDown->GetBinError(i)));
      // addTemplate("WZ_ZE_WZDieleSBin"+binNumber+"Up", vars, wspace, WZ_ZE_histBinUp);
      // addTemplate("WZ_ZE_WZDieleSBin"+binNumber+"Down", vars, wspace, WZ_ZE_histBinDown);
    }

    TH1F* QCD_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_jetfake"))->Clone("QCD_WM_histBinUp");
    TH1F* QCD_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_jetfake"))->Clone("QCD_WM_histBinDown");
    QCD_WM_histBinUp->SetBinContent(i, QCD_WM_histBinUp->GetBinContent(i)+QCD_WM_histBinUp->GetBinError(i));
    QCD_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*QCD_WM_histBinDown->GetBinContent(i), QCD_WM_histBinDown->GetBinContent(i)-QCD_WM_histBinDown->GetBinError(i)));
    addTemplate("QCD_WM_QCDMonomuSBin"+binNumber+"Up", vars, wspace, QCD_WM_histBinUp);
    addTemplate("QCD_WM_QCDMonomuSBin"+binNumber+"Down", vars, wspace, QCD_WM_histBinDown);
    // TH1F* ZllG_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_ZllG_combined"))->Clone("ZllG_WM_histBinUp");
    // TH1F* ZllG_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_ZllG_combined"))->Clone("ZllG_WM_histBinDown");
    // ZllG_WM_histBinUp->SetBinContent(i, ZllG_WM_histBinUp->GetBinContent(i)+ZllG_WM_histBinUp->GetBinError(i));
    // ZllG_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*ZllG_WM_histBinDown->GetBinContent(i), ZllG_WM_histBinDown->GetBinContent(i)-ZllG_WM_histBinDown->GetBinError(i)));
    // addTemplate("ZllG_WM_ZllGMonomuSBin"+binNumber+"Up", vars, wspace, ZllG_WM_histBinUp);
    // addTemplate("ZllG_WM_ZllGMonomuSBin"+binNumber+"Down", vars, wspace, ZllG_WM_histBinDown);
    // TH1F* TTG_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_TTG"))->Clone("TTG_WM_histBinUp");
    // TH1F* TTG_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_TTG"))->Clone("TTG_WM_histBinDown");
    // TTG_WM_histBinUp->SetBinContent(i, TTG_WM_histBinUp->GetBinContent(i)+TTG_WM_histBinUp->GetBinError(i));
    // TTG_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*TTG_WM_histBinDown->GetBinContent(i), TTG_WM_histBinDown->GetBinContent(i)-TTG_WM_histBinDown->GetBinError(i)));
    // addTemplate("TTG_WM_TTGMonomuSBin"+binNumber+"Up", vars, wspace, TTG_WM_histBinUp);
    // addTemplate("TTG_WM_TTGMonomuSBin"+binNumber+"Down", vars, wspace, TTG_WM_histBinDown);
    // TH1F* TG_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_TG"))->Clone("TG_WM_histBinUp");
    // TH1F* TG_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_TG"))->Clone("TG_WM_histBinDown");
    // TG_WM_histBinUp->SetBinContent(i, TG_WM_histBinUp->GetBinContent(i)+TG_WM_histBinUp->GetBinError(i));
    // TG_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*TG_WM_histBinDown->GetBinContent(i), TG_WM_histBinDown->GetBinContent(i)-TG_WM_histBinDown->GetBinError(i)));
    // addTemplate("TG_WM_TGMonomuSBin"+binNumber+"Up", vars, wspace, TG_WM_histBinUp);
    // addTemplate("TG_WM_TGMonomuSBin"+binNumber+"Down", vars, wspace, TG_WM_histBinDown);
    // TH1F* Diphoton_WM_histBinUp = (TH1F*)((TH1F*)WmnGfile->Get("histo_diphoton"))->Clone("Diphoton_WM_histBinUp");
    // TH1F* Diphoton_WM_histBinDown = (TH1F*)((TH1F*)WmnGfile->Get("histo_diphoton"))->Clone("Diphoton_WM_histBinDown");
    // Diphoton_WM_histBinUp->SetBinContent(i, Diphoton_WM_histBinUp->GetBinContent(i)+Diphoton_WM_histBinUp->GetBinError(i));
    // Diphoton_WM_histBinDown->SetBinContent(i, TMath::Max(0.01*Diphoton_WM_histBinDown->GetBinContent(i), Diphoton_WM_histBinDown->GetBinContent(i)-Diphoton_WM_histBinDown->GetBinError(i)));
    // addTemplate("Diphoton_WM_DiphotonMonomuSBin"+binNumber+"Up", vars, wspace, Diphoton_WM_histBinUp);
    // addTemplate("Diphoton_WM_DiphotonMonomuSBin"+binNumber+"Down", vars, wspace, Diphoton_WM_histBinDown);

    TH1F* QCD_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_jetfake"))->Clone("QCD_WE_histBinUp");
    TH1F* QCD_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_jetfake"))->Clone("QCD_WE_histBinDown");
    QCD_WE_histBinUp->SetBinContent(i, QCD_WE_histBinUp->GetBinContent(i)+QCD_WE_histBinUp->GetBinError(i));
    QCD_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*QCD_WE_histBinDown->GetBinContent(i), QCD_WE_histBinDown->GetBinContent(i)-QCD_WE_histBinDown->GetBinError(i)));
    addTemplate("QCD_WE_QCDMonoeleSBin"+binNumber+"Up", vars, wspace, QCD_WE_histBinUp);
    addTemplate("QCD_WE_QCDMonoeleSBin"+binNumber+"Down", vars, wspace, QCD_WE_histBinDown);
    // TH1F* Elefake_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_elefake"))->Clone("Elefake_WE_histBinUp");
    // TH1F* Elefake_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_elefake"))->Clone("Elefake_WE_histBinDown");
    // Elefake_WE_histBinUp->SetBinContent(i, Elefake_WE_histBinUp->GetBinContent(i)+Elefake_WE_histBinUp->GetBinError(i));
    // Elefake_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*Elefake_WE_histBinDown->GetBinContent(i), Elefake_WE_histBinDown->GetBinContent(i)-Elefake_WE_histBinDown->GetBinError(i)));
    // addTemplate("Elefake_WE_EleMonoeleSBin"+binNumber+"Up", vars, wspace, Elefake_WE_histBinUp);
    // addTemplate("Elefake_WE_EleMonoeleSBin"+binNumber+"Down", vars, wspace, Elefake_WE_histBinDown);
    // TH1F* ZllG_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_ZllG_combined"))->Clone("ZllG_WE_histBinUp");
    // TH1F* ZllG_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_ZllG_combined"))->Clone("ZllG_WE_histBinDown");
    // ZllG_WE_histBinUp->SetBinContent(i, ZllG_WE_histBinUp->GetBinContent(i)+ZllG_WE_histBinUp->GetBinError(i));
    // ZllG_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*ZllG_WE_histBinDown->GetBinContent(i), ZllG_WE_histBinDown->GetBinContent(i)-ZllG_WE_histBinDown->GetBinError(i)));
    // addTemplate("ZllG_WE_ZllGMonoeleSBin"+binNumber+"Up", vars, wspace, ZllG_WE_histBinUp);
    // addTemplate("ZllG_WE_ZllGMonoeleSBin"+binNumber+"Down", vars, wspace, ZllG_WE_histBinDown);
    // TH1F* TTG_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_TTG"))->Clone("TTG_WE_histBinUp");
    // TH1F* TTG_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_TTG"))->Clone("TTG_WE_histBinDown");
    // TTG_WE_histBinUp->SetBinContent(i, TTG_WE_histBinUp->GetBinContent(i)+TTG_WE_histBinUp->GetBinError(i));
    // TTG_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*TTG_WE_histBinDown->GetBinContent(i), TTG_WE_histBinDown->GetBinContent(i)-TTG_WE_histBinDown->GetBinError(i)));
    // addTemplate("TTG_WE_TTGMonoeleSBin"+binNumber+"Up", vars, wspace, TTG_WE_histBinUp);
    // addTemplate("TTG_WE_TTGMonoeleSBin"+binNumber+"Down", vars, wspace, TTG_WE_histBinDown);
    // TH1F* TG_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_TG"))->Clone("TG_WE_histBinUp");
    // TH1F* TG_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_TG"))->Clone("TG_WE_histBinDown");
    // TG_WE_histBinUp->SetBinContent(i, TG_WE_histBinUp->GetBinContent(i)+TG_WE_histBinUp->GetBinError(i));
    // TG_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*TG_WE_histBinDown->GetBinContent(i), TG_WE_histBinDown->GetBinContent(i)-TG_WE_histBinDown->GetBinError(i)));
    // addTemplate("TG_WE_TGMonoeleSBin"+binNumber+"Up", vars, wspace, TG_WE_histBinUp);
    // addTemplate("TG_WE_TGMonoeleSBin"+binNumber+"Down", vars, wspace, TG_WE_histBinDown);
    // TH1F* Diphoton_WE_histBinUp = (TH1F*)((TH1F*)WenGfile->Get("histo_diphoton"))->Clone("Diphoton_WE_histBinUp");
    // TH1F* Diphoton_WE_histBinDown = (TH1F*)((TH1F*)WenGfile->Get("histo_diphoton"))->Clone("Diphoton_WE_histBinDown");
    // Diphoton_WE_histBinUp->SetBinContent(i, Diphoton_WE_histBinUp->GetBinContent(i)+Diphoton_WE_histBinUp->GetBinError(i));
    // Diphoton_WE_histBinDown->SetBinContent(i, TMath::Max(0.01*Diphoton_WE_histBinDown->GetBinContent(i), Diphoton_WE_histBinDown->GetBinContent(i)-Diphoton_WE_histBinDown->GetBinError(i)));
    // addTemplate("Diphoton_WE_DiphotonMonoeleSBin"+binNumber+"Up", vars, wspace, Diphoton_WE_histBinUp);
    // addTemplate("Diphoton_WE_DiphotonMonoeleSBin"+binNumber+"Down", vars, wspace, Diphoton_WE_histBinDown);
  }

  // ---------------------------- Write out the workspace -----------------------------------------------------------------//
  outfile->cd();
  wspace.Write();
  outfile->Close();
  transfer_factors_file->Close();
  signalabove0p5file->Close();
  signalbelow0p5file->Close();
  ZnnGabove0p5file->Close();
  ZnnGbelow0p5file->Close();
  WenGfile->Close();
  WmnGfile->Close();
  ZeeGfile->Close();
  ZmmGfile->Close();
}
 

void createWorkspaces_Pt_aTGC(){
  vector<string> sample_names_aTGC;
  sample_names_aTGC.clear();
  vector<float> signal_multipliers_aTGC;
  signal_multipliers_aTGC.clear();
  
  // aTGC model
  vector<vector<double>> Zgg_yield_2Dfit;
  vector<vector<double>> Zgg_JESUp_yield_2Dfit;
  vector<vector<double>> Zgg_JESDown_yield_2Dfit;
  vector<vector<double>> Zgg_PESUp_yield_2Dfit;
  vector<vector<double>> Zgg_PESDown_yield_2Dfit;
  // vector<vector<double>> ZZg_yield_2Dfit;
  // vector<vector<double>> ZZg_JESUp_yield_2Dfit;
  // vector<vector<double>> ZZg_JESDown_yield_2Dfit;
  // vector<vector<double>> ZZg_PESUp_yield_2Dfit;
  // vector<vector<double>> ZZg_PESDown_yield_2Dfit;
  Zgg_yield_2Dfit.clear();
  Zgg_JESUp_yield_2Dfit.clear();
  Zgg_JESDown_yield_2Dfit.clear();
  Zgg_PESUp_yield_2Dfit.clear();
  Zgg_PESDown_yield_2Dfit.clear();
  // ZZg_yield_2Dfit.clear();
  // ZZg_JESUp_yield_2Dfit.clear();
  // ZZg_JESDown_yield_2Dfit.clear();
  // ZZg_PESUp_yield_2Dfit.clear();
  // ZZg_PESDown_yield_2Dfit.clear();

  // Direct output from znng_plotter_aTGC.C
  Zgg_yield_2Dfit.push_back({102.499, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_JESUp_yield_2Dfit.push_back({102.296, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_JESDown_yield_2Dfit.push_back({102.566, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_PESUp_yield_2Dfit.push_back({106.422, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_PESDown_yield_2Dfit.push_back({100.67, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_yield_2Dfit.push_back({135.51, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_JESUp_yield_2Dfit.push_back({134.97, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_JESDown_yield_2Dfit.push_back({135.982, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_PESUp_yield_2Dfit.push_back({137.271, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_PESDown_yield_2Dfit.push_back({133.209, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_yield_2Dfit.push_back({61.8197, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_JESUp_yield_2Dfit.push_back({61.8197, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_JESDown_yield_2Dfit.push_back({61.887, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_PESUp_yield_2Dfit.push_back({62.0274, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_PESDown_yield_2Dfit.push_back({60.2005, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_yield_2Dfit.push_back({46.3898, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_JESUp_yield_2Dfit.push_back({46.3898, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_JESDown_yield_2Dfit.push_back({46.4568, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_PESUp_yield_2Dfit.push_back({48.3411, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_PESDown_yield_2Dfit.push_back({46.0504, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_yield_2Dfit.push_back({20.0128, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_JESUp_yield_2Dfit.push_back({20.0128, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_JESDown_yield_2Dfit.push_back({20.0128, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_PESUp_yield_2Dfit.push_back({20.5483, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_PESDown_yield_2Dfit.push_back({19.4094, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_yield_2Dfit.push_back({7.69449, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_JESUp_yield_2Dfit.push_back({7.69449, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_JESDown_yield_2Dfit.push_back({7.69449, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_PESUp_yield_2Dfit.push_back({7.69449, 0.0, 0.0, 0.0, 0.0, 0.0});
  Zgg_PESDown_yield_2Dfit.push_back({7.56186, 0.0, 0.0, 0.0, 0.0, 0.0});

  // Alternative method:
  // // Output from aTGC_2Dpoly_fits.py
  // // Each pushed back vector is a list of fit parameters for a 2D fit of the form
  // //   [0] + [1]*h3 + [2]*h4 + [3]*h3*h3 + [4]*h4*h4 + [5]*h3*h4
  // // representing the aTGC signal yield within one reconstructed phoET bin
  // Zgg_yield_2Dfit.push_back({69.2946075759, 385.774642091, 30128.8514956, -126506.681288, -5009043398.81, 11871308.7794});
  // Zgg_JESUp_yield_2Dfit.push_back({67.0169265115, 387.778515299, 30987.4887287, 254615.908606, -185962.793264, 4115315.15205});
  // Zgg_JESDown_yield_2Dfit.push_back({69.7092150243, 340.34067475, 30763.0158595, -158346.266038, -5368082911.66, 10274496.9364});
  // Zgg_PESUp_yield_2Dfit.push_back({67.0169265115, 387.778515299, 30987.4887287, 254615.908606, -185962.793264, 4115315.15205});
  // Zgg_PESDown_yield_2Dfit.push_back({69.7092150243, 340.34067475, 30763.0158595, -158346.266038, -5368082911.66, 10274496.9364});
  // Zgg_yield_2Dfit.push_back({88.7205452524, -96.113306941, 114544.529909, 397167.384582, 1445779770.07, -55763504.1152});
  // Zgg_JESUp_yield_2Dfit.push_back({88.1483300502, -119.084762559, 109562.234537, 493365.382694, 1414213001.29, -51367880.2351});
  // Zgg_JESDown_yield_2Dfit.push_back({88.8339892787, -109.330642298, 116762.602041, 442097.013292, 1499900870.94, -56386950.7642});
  // Zgg_PESUp_yield_2Dfit.push_back({88.1483300502, -119.084762559, 109562.234537, 493365.382694, 1414213001.29, -51367880.2351});
  // Zgg_PESDown_yield_2Dfit.push_back({88.8339892787, -109.330642298, 116762.602041, 442097.013292, 1499900870.94, -56386950.7642});
  // Zgg_yield_2Dfit.push_back({37.9852605868, -221.762898166, -40249.1281505, 663773.241151, 8630072163.24, -60565880.8678});
  // Zgg_JESUp_yield_2Dfit.push_back({37.9332598827, -234.983396787, -40594.9914679, 656879.490193, 8686603637.42, -61987797.1096});
  // Zgg_JESDown_yield_2Dfit.push_back({38.2043786444, -232.848827624, -39558.2840223, 623995.563811, 8512503645.17, -62318098.5043});
  // Zgg_PESUp_yield_2Dfit.push_back({37.9332598827, -234.983396787, -40594.9914679, 656879.490193, 8686603637.42, -61987797.1096});
  // Zgg_PESDown_yield_2Dfit.push_back({38.2043786444, -232.848827624, -39558.2840223, 623995.563811, 8512503645.17, -62318098.5043});
  // Zgg_yield_2Dfit.push_back({28.4206249008, -291.009979862, -100.86346948, 1156103.2795, 16950774523.6, -204658082.941});
  // Zgg_JESUp_yield_2Dfit.push_back({28.5089605233, -285.177424798, 369.877518934, 1137183.05106, 16901104873.5, -205938985.612});
  // Zgg_JESDown_yield_2Dfit.push_back({35.2689154058, -318.984766524, 835.494787254, -156857.479785, 1548798.38447, -155727241.458});
  // Zgg_PESUp_yield_2Dfit.push_back({28.5089605233, -285.177424798, 369.877518934, 1137183.05106, 16901104873.5, -205938985.612});
  // Zgg_PESDown_yield_2Dfit.push_back({35.2689154058, -318.984766524, 835.494787254, -156857.479785, 1548798.38447, -155727241.458});
  // Zgg_yield_2Dfit.push_back({10.9408485806, 103.563544868, -11042.2160819, 2868530.34613, 90658377550.1, -812261690.586});
  // Zgg_JESUp_yield_2Dfit.push_back({10.8146445447, 92.9711143126, -9644.77501838, 2901687.88965, 90804078775.0, -813877951.935});
  // Zgg_JESDown_yield_2Dfit.push_back({10.8339368157, 97.6058016024, -9665.18730173, 2895221.94732, 90757783292.7, -813246366.952});
  // Zgg_PESUp_yield_2Dfit.push_back({10.8146445447, 92.9711143126, -9644.77501838, 2901687.88965, 90804078775.0, -813877951.935});
  // Zgg_PESDown_yield_2Dfit.push_back({10.8339368157, 97.6058016024, -9665.18730173, 2895221.94732, 90757783292.7, -813246366.952});
  // Zgg_yield_2Dfit.push_back({18.2879260694, 304.998808486, 369372.67741, 8670016.36118, 6.15447938713e+12, -13927199363.1});
  // Zgg_JESUp_yield_2Dfit.push_back({18.5134959749, 190.225548524, 286650.141605, 8621693.63699, 6.15556779586e+12, -13927100043.6});
  // Zgg_JESDown_yield_2Dfit.push_back({18.4378826604, 197.791979153, 289385.786391, 8632657.95548, 6.15391840835e+12, -13923538483.5});
  // Zgg_PESUp_yield_2Dfit.push_back({18.5134959749, 190.225548524, 286650.141605, 8621693.63699, 6.15556779586e+12, -13927100043.6});
  // Zgg_PESDown_yield_2Dfit.push_back({18.4378826604, 197.791979153, 289385.786391, 8632657.95548, 6.15391840835e+12, -13923538483.5});

  vector<vector<vector<double>>> Zgg_yield_2Dfits;
  Zgg_yield_2Dfits.clear();
  Zgg_yield_2Dfits.push_back(Zgg_yield_2Dfit);
  Zgg_yield_2Dfits.push_back(Zgg_JESUp_yield_2Dfit);
  Zgg_yield_2Dfits.push_back(Zgg_JESDown_yield_2Dfit);
  Zgg_yield_2Dfits.push_back(Zgg_PESUp_yield_2Dfit);
  Zgg_yield_2Dfits.push_back(Zgg_PESDown_yield_2Dfit);
  // vector<vector<vector<double>>> ZZg_yield_2Dfits;
  // ZZg_yield_2Dfits.clear();
  // ZZg_yield_2Dfits.push_back(ZZg_yield_2Dfit);
  // ZZg_yield_2Dfits.push_back(ZZg_JESUp_yield_2Dfit);
  // ZZg_yield_2Dfits.push_back(ZZg_JESDown_yield_2Dfit);
  // ZZg_yield_2Dfits.push_back(ZZg_PESUp_yield_2Dfit);
  // ZZg_yield_2Dfits.push_back(ZZg_PESDown_yield_2Dfit);

  double h3 = 0.;
  double h4 = 1e-6;
  string h3string = boost::lexical_cast<string>(boost::format("%4.3e") % h3);
  string h4string = boost::lexical_cast<string>(boost::format("%4.3e") % h4);
  std::replace(h3string.begin(), h3string.end(), '.', 'p');
  std::replace(h3string.begin(), h3string.end(), '-', 'M');
  std::replace(h3string.begin(), h3string.end(), '+', 'P');
  std::replace(h4string.begin(), h4string.end(), '.', 'p');
  std::replace(h4string.begin(), h4string.end(), '-', 'M');
  std::replace(h4string.begin(), h4string.end(), '+', 'P');
  string samplename_suffix = "_h3"+h3string+"_h4"+h4string;
  // Do not connect Wgamma and Zgamma yields when aTGC is the signal: set connectWZ (the second argument) to false
  do_createWorkspace("Zgg", false, sample_names_aTGC, signal_multipliers_aTGC, samplename_suffix, true, h3, h4, Zgg_yield_2Dfits);
  
  // Alternative method:
  // // Examine h3,h4 values evenly spaced around an ellipse
  // // The spacing of points around the ellipse is governed by nPhiSteps_per_quadrant
  // // The semimajor axis of the ellipse has a length = rScale*h3_halfSpan,
  // //   and the semiminor axis has a length = rScale*h4_halfSpan
  // // rScale=1.0 corresponds to an ellipse that passes through
  // //   h3 = +/-h3_halfSpan and h4 = +/-h4_halfSpan
  // double h3_halfSpan = 0.0004;
  // double h4_halfSpan = 0.0000004;
  // vector<pair<double, int>> aTGC_rScale_and_nPhiSteps_per_quadrant;
  // aTGC_rScale_and_nPhiSteps_per_quadrant.clear();
  // // Some options to try out:
  // aTGC_rScale_and_nPhiSteps_per_quadrant.push_back(pair<double, int>(1.0, 1));
  // aTGC_rScale_and_nPhiSteps_per_quadrant.push_back(pair<double, int>(0.1, 1));
  // aTGC_rScale_and_nPhiSteps_per_quadrant.push_back(pair<double, int>(0.01, 1));

  // for(int idx = 0; idx < aTGC_rScale_and_nPhiSteps_per_quadrant.size(); idx++) {
  //   double rScale = aTGC_rScale_and_nPhiSteps_per_quadrant[idx].first;
  //   int nPhiSteps_per_quadrant = aTGC_rScale_and_nPhiSteps_per_quadrant[idx].second;
  //   for(int phiParam_idx = 0; phiParam_idx < 4*nPhiSteps_per_quadrant; phiParam_idx++) {
  //     double phi, h3, h4;
  //     phi = h3 = h4 = 0.0;
  //     if (phiParam_idx/nPhiSteps_per_quadrant < 1) {
  //       phi = 2*TMath::Pi() * phiParam_idx/nPhiSteps_per_quadrant;
  //       h3 = rScale*h3_halfSpan*TMath::Cos(phi);
  //       h4 = rScale*h4_halfSpan*TMath::Sin(phi);
  //     }
  //     else if (phiParam_idx/nPhiSteps_per_quadrant < 2) {
  //       phi = 2*TMath::Pi() * (phiParam_idx-nPhiSteps_per_quadrant)/nPhiSteps_per_quadrant;
  //       h3 = -rScale*h3_halfSpan*TMath::Sin(phi);
  //       h4 = rScale*h4_halfSpan*TMath::Cos(phi);
  //     }
  //     else if (phiParam_idx/nPhiSteps_per_quadrant < 3) {
  //       phi = 2*TMath::Pi() * (phiParam_idx-2*nPhiSteps_per_quadrant)/nPhiSteps_per_quadrant;
  //       h3 = -rScale*h3_halfSpan*TMath::Cos(phi);
  //       h4 = -rScale*h4_halfSpan*TMath::Sin(phi);
  //     }
  //     else {
  //       phi = 2*TMath::Pi() * (phiParam_idx-3*nPhiSteps_per_quadrant)/nPhiSteps_per_quadrant;
  //       h3 = rScale*h3_halfSpan*TMath::Sin(phi);
  //       h4 = -rScale*h4_halfSpan*TMath::Cos(phi);
  //     }
  //     // string h3string = boost::lexical_cast<string>(boost::format("%.5f") % h3);
  //     // string h4string = boost::lexical_cast<string>(boost::format("%.7f") % h4);
  //     string h3string = boost::lexical_cast<string>(boost::format("%4.3e") % h3);
  //     string h4string = boost::lexical_cast<string>(boost::format("%4.3e") % h4);
  //     std::replace(h3string.begin(), h3string.end(), '.', 'p');
  //     std::replace(h3string.begin(), h3string.end(), '-', 'M');
  //     std::replace(h3string.begin(), h3string.end(), '+', 'P');
  //     std::replace(h4string.begin(), h4string.end(), '.', 'p');
  //     std::replace(h4string.begin(), h4string.end(), '-', 'M');
  //     std::replace(h4string.begin(), h4string.end(), '+', 'P');
  //     string samplename_suffix = "_h3"+h3string+"_h4"+h4string;

  //     // Do not connect Wgamma and Zgamma yields when aTGC is the signal: set connectWZ (the second argument) to false
  //     do_createWorkspace("Zgg", false, sample_names_aTGC, signal_multipliers_aTGC, samplename_suffix, true, h3, h4, Zgg_yield_2Dfits);
  //     do_createWorkspace("ZZg", false, sample_names_aTGC, signal_multipliers_aTGC, samplename_suffix, true, h3, h4, ZZg_yield_2Dfits);
  //   }
  // }
  
}