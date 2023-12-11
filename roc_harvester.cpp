#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>

#include "TROOT.h"
#include "TMath.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TSystem.h"
#include "TH2F.h"
#include "TEllipse.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TString.h"
#include "TLorentzVector.h"
#include "TChain.h"
#include "TH1D.h"
#include "TObjString.h"
#include "TStyle.h"
#include "TLatex.h"
#include "TLegend.h"
#include "TGraph.h"

void draw_iso_efficiency(bool zeroiso);       // draw iso. efficiency and normalization as a function of rel. iso.
void draw_reliso_roc(bool zeroiso);           // draw ROC curves as a function of rel. iso.
void draw_pt_roc(bool zeroiso);               // draw ROC curves as a function of pT sum of charged tracks
void pTeff();                                 // draw efficiency plots as a function of pT of muons
void N_genMatched();                          // draw plots to check whether muon's track is matched with a Gen particle

using namespace std;


void draw_iso_efficiency(bool zeroiso) {
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  // PU200
  f_PU200_prompt = new TFile("data/231203/harvester_PU200_prompt_160_20.root");
  f_PU200_nonprompt = new TFile("data/231203/harvester_PU200_nonprompt_90.root");
  // noPU
  f_noPU_prompt = new TFile("data/231203/harvester_noPU_prompt_19.root");
  f_noPU_nonprompt = new TFile("data/231203/harvester_noPU_nonprompt_9.root");

  // Histograms for Barrel region
  TH1D *h_PU200_prompt_EB, *h_PU200_nonprompt_EB, *h_noPU_prompt_EB, *h_noPU_nonprompt_EB;
  TH1D *h_PU200_prompt_2sigma_EB, *h_PU200_nonprompt_2sigma_EB, *h_noPU_prompt_2sigma_EB, *h_noPU_2sigma_nonprompt_EB;
  TH1D *h_PU200_prompt_3sigma_EB, *h_PU200_nonprompt_3sigma_EB, *h_noPU_prompt_3sigma_EB, *h_noPU_3sigma_nonprompt_EB;
  TH1D *h_PU200_prompt_4sigma_EB, *h_PU200_nonprompt_4sigma_EB, *h_noPU_prompt_4sigma_EB, *h_noPU_4sigma_nonprompt_EB;
  TH1D *h_PU200_prompt_40_EB, *h_PU200_nonprompt_40_EB, *h_noPU_prompt_40_EB, *h_noPU_40_nonprompt_EB;
  TH1D *h_PU200_prompt_60_EB, *h_PU200_nonprompt_60_EB, *h_noPU_prompt_60_EB, *h_noPU_60_nonprompt_EB;
  TH1D *h_PU200_prompt_80_EB, *h_PU200_nonprompt_80_EB, *h_noPU_prompt_80_EB, *h_noPU_80_nonprompt_EB;
  TH1D *h_PU200_prompt_100_EB, *h_PU200_nonprompt_100_EB, *h_noPU_prompt_100_EB, *h_noPU_100_nonprompt_EB;
  TH1D *h_PU200_prompt_gen_EB, *h_PU200_nonprompt_gen_EB, *h_noPU_prompt_gen_EB, *h_noPU_nonprompt_gen_EB;

  // Histograms for Endcap region
  TH1D *h_PU200_prompt_EE, *h_PU200_nonprompt_EE, *h_noPU_prompt_EE, *h_noPU_nonprompt_EE;
  TH1D *h_PU200_prompt_2sigma_EE, *h_PU200_nonprompt_2sigma_EE, *h_noPU_prompt_2sigma_EE, *h_noPU_2sigma_nonprompt_EE;
  TH1D *h_PU200_prompt_3sigma_EE, *h_PU200_nonprompt_3sigma_EE, *h_noPU_prompt_3sigma_EE, *h_noPU_3sigma_nonprompt_EE;
  TH1D *h_PU200_prompt_4sigma_EE, *h_PU200_nonprompt_4sigma_EE, *h_noPU_prompt_4sigma_EE, *h_noPU_4sigma_nonprompt_EE;
  TH1D *h_PU200_prompt_40_EE, *h_PU200_nonprompt_40_EE, *h_noPU_prompt_40_EE, *h_noPU_40_nonprompt_EE;
  TH1D *h_PU200_prompt_60_EE, *h_PU200_nonprompt_60_EE, *h_noPU_prompt_60_EE, *h_noPU_60_nonprompt_EE;
  TH1D *h_PU200_prompt_80_EE, *h_PU200_nonprompt_80_EE, *h_noPU_prompt_80_EE, *h_noPU_80_nonprompt_EE;
  TH1D *h_PU200_prompt_100_EE, *h_PU200_nonprompt_100_EE, *h_noPU_prompt_100_EE, *h_noPU_100_nonprompt_EE;
  TH1D *h_PU200_prompt_gen_EE, *h_PU200_nonprompt_gen_EE, *h_noPU_prompt_gen_EE, *h_noPU_nonprompt_gen_EE;


  TString dir="/DQMData/Run 1/MTD/Run summary/MuonIso/";
  // Prompt
    // Barrel region
  h_PU200_prompt_EB 	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_Sig_EB");
  h_PU200_prompt_4sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EB");
  h_PU200_prompt_3sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EB");
  h_PU200_prompt_2sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EB");
  h_PU200_prompt_40_EB     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Sig_EB");
  h_PU200_prompt_60_EB     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Sig_EB");
  h_PU200_prompt_80_EB     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Sig_EB");
  h_PU200_prompt_100_EB    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Sig_EB");
  h_PU200_prompt_gen_EB	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_gen_Sig_EB");
  h_noPU_prompt_EB    	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_Sig_EB");
  h_noPU_prompt_gen_EB	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_gen_Sig_EB");
    // Endcap region
  h_PU200_prompt_EE 	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_Sig_EE");
  h_PU200_prompt_4sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EE");
  h_PU200_prompt_3sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EE");
  h_PU200_prompt_2sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EE");
  h_PU200_prompt_40_EE     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Sig_EE");
  h_PU200_prompt_60_EE     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Sig_EE");
  h_PU200_prompt_80_EE     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Sig_EE");
  h_PU200_prompt_100_EE    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Sig_EE");
  h_PU200_prompt_gen_EE	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_gen_Sig_EE");
  h_noPU_prompt_EE    	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_Sig_EE");
  h_noPU_prompt_gen_EE	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_gen_Sig_EE");
  // Nonprompt
    // Barrel region
  h_PU200_nonprompt_EB 	      = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_Bkg_EB");
  h_PU200_nonprompt_4sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EB");
  h_PU200_nonprompt_3sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EB");
  h_PU200_nonprompt_2sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EB");
  h_PU200_nonprompt_40_EB     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Bkg_EB");
  h_PU200_nonprompt_60_EB     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Bkg_EB");
  h_PU200_nonprompt_80_EB     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Bkg_EB");
  h_PU200_nonprompt_100_EB    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Bkg_EB");
  h_PU200_nonprompt_gen_EB    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_gen_Bkg_EB");
  h_noPU_nonprompt_EB         = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_Bkg_EB");
  h_noPU_nonprompt_gen_EB         = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_gen_Bkg_EB");
    // Endcap region
  h_PU200_nonprompt_EE 	      = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_Bkg_EE");
  h_PU200_nonprompt_4sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EE");
  h_PU200_nonprompt_3sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EE");
  h_PU200_nonprompt_2sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EE");
  h_PU200_nonprompt_40_EE     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Bkg_EE");
  h_PU200_nonprompt_60_EE     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Bkg_EE");
  h_PU200_nonprompt_80_EE     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Bkg_EE");
  h_PU200_nonprompt_100_EE    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Bkg_EE");
  h_PU200_nonprompt_gen_EE    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_gen_Bkg_EE");
  h_noPU_nonprompt_EE         = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_Bkg_EE");
  h_noPU_nonprompt_gen_EE     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_gen_Bkg_EE");

  int nbin = h_PU200_prompt_EB->GetNbinsX();

  cout << "total muon of PU200 prompt in Barrel   : " << h_PU200_prompt_EB->Integral(1,nbin) << endl;
  cout << "total muon of noPU prompt in Barrel    : " << h_noPU_prompt_EB->Integral(1,nbin) << endl;
  cout << "total muon of PU200 nonprompt in Barrel: " << h_PU200_nonprompt_EB->Integral(1,nbin) << endl;
  cout << "total muon of noPU nonprompt in Barrel : " << h_noPU_nonprompt_EB->Integral(1,nbin) << endl;
  cout << "total muon of PU200 prompt in Endcap   : " << h_PU200_prompt_EE->Integral(1,nbin) << endl;
  cout << "total muon of noPU prompt in Endcap    : " << h_noPU_prompt_EE->Integral(1,nbin) << endl;
  cout << "total muon of PU200 nonprompt in Endcap: " << h_PU200_nonprompt_EE->Integral(1,nbin) << endl;
  cout << "total muon of noPU nonprompt in Endcap : " << h_noPU_nonprompt_EE->Integral(1,nbin) << endl;

  // Define vectors to store the iso efficiency or normalization of muons
    // PU200
      // prompt
  vector<double> prompt_eff_PU200_EB={0}, prompt_norm_PU200_EB={0}, prompt_eff_PU200_EE={0}, prompt_norm_PU200_EE={0};
  vector<double> prompt_eff_PU200_2sigma_EB={0}, prompt_norm_PU200_2sigma_EB={0}, prompt_eff_PU200_2sigma_EE={0}, prompt_norm_PU200_2sigma_EE={0};
  vector<double> prompt_eff_PU200_3sigma_EB={0}, prompt_norm_PU200_3sigma_EB={0}, prompt_eff_PU200_3sigma_EE={0}, prompt_norm_PU200_3sigma_EE={0};
  vector<double> prompt_eff_PU200_4sigma_EB={0}, prompt_norm_PU200_4sigma_EB={0}, prompt_eff_PU200_4sigma_EE={0}, prompt_norm_PU200_4sigma_EE={0};
  vector<double> prompt_eff_PU200_40_EB={0}, prompt_norm_PU200_40_EB={0}, prompt_eff_PU200_40_EE={0}, prompt_norm_PU200_40_EE={0};
  vector<double> prompt_eff_PU200_60_EB={0}, prompt_norm_PU200_60_EB={0}, prompt_eff_PU200_60_EE={0}, prompt_norm_PU200_60_EE={0};
  vector<double> prompt_eff_PU200_80_EB={0}, prompt_norm_PU200_80_EB={0}, prompt_eff_PU200_80_EE={0}, prompt_norm_PU200_80_EE={0};
  vector<double> prompt_eff_PU200_100_EB={0}, prompt_norm_PU200_100_EB={0}, prompt_eff_PU200_100_EE={0}, prompt_norm_PU200_100_EE={0};
      // nonprompt
  vector<double> nonprompt_eff_PU200_EB={0}, nonprompt_norm_PU200_EB={0}, nonprompt_eff_PU200_EE={0}, nonprompt_norm_PU200_EE={0};
  vector<double> nonprompt_eff_PU200_2sigma_EB={0}, nonprompt_norm_PU200_2sigma_EB={0}, nonprompt_eff_PU200_2sigma_EE={0}, nonprompt_norm_PU200_2sigma_EE={0};
  vector<double> nonprompt_eff_PU200_3sigma_EB={0}, nonprompt_norm_PU200_3sigma_EB={0}, nonprompt_eff_PU200_3sigma_EE={0}, nonprompt_norm_PU200_3sigma_EE={0};
  vector<double> nonprompt_eff_PU200_4sigma_EB={0}, nonprompt_norm_PU200_4sigma_EB={0}, nonprompt_eff_PU200_4sigma_EE={0}, nonprompt_norm_PU200_4sigma_EE={0};
  vector<double> nonprompt_eff_PU200_40_EB={0}, nonprompt_norm_PU200_40_EB={0}, nonprompt_eff_PU200_40_EE={0}, nonprompt_norm_PU200_40_EE={0};
  vector<double> nonprompt_eff_PU200_60_EB={0}, nonprompt_norm_PU200_60_EB={0}, nonprompt_eff_PU200_60_EE={0}, nonprompt_norm_PU200_60_EE={0};
  vector<double> nonprompt_eff_PU200_80_EB={0}, nonprompt_norm_PU200_80_EB={0}, nonprompt_eff_PU200_80_EE={0}, nonprompt_norm_PU200_80_EE={0};
  vector<double> nonprompt_eff_PU200_100_EB={0}, nonprompt_norm_PU200_100_EB={0}, nonprompt_eff_PU200_100_EE={0}, nonprompt_norm_PU200_100_EE={0};
      // GEN case
  vector<double> prompt_eff_PU200_gen_EB={0}, prompt_norm_PU200_gen_EB={0}, prompt_eff_PU200_gen_EE={0}, prompt_norm_PU200_gen_EE={0};
  vector<double> nonprompt_eff_PU200_gen_EB={0}, nonprompt_norm_PU200_gen_EB={0}, nonprompt_eff_PU200_gen_EE={0}, nonprompt_norm_PU200_gen_EE={0};
    // noPU
      // prompt
  vector<double> prompt_eff_noPU_EB={0}, prompt_norm_noPU_EB={0}, prompt_eff_noPU_EE={0}, prompt_norm_noPU_EE={0};
  vector<double> prompt_eff_noPU_2sigma_EB={0}, prompt_norm_noPU_2sigma_EB={0}, prompt_eff_noPU_2sigma_EE={0}, prompt_norm_noPU_2sigma_EE={0};
  vector<double> prompt_eff_noPU_3sigma_EB={0}, prompt_norm_noPU_3sigma_EB={0}, prompt_eff_noPU_3sigma_EE={0}, prompt_norm_noPU_3sigma_EE={0};
  vector<double> prompt_eff_noPU_4sigma_EB={0}, prompt_norm_noPU_4sigma_EB={0}, prompt_eff_noPU_4sigma_EE={0}, prompt_norm_noPU_4sigma_EE={0};
  vector<double> prompt_eff_noPU_40_EB={0}, prompt_norm_noPU_40_EB={0}, prompt_eff_noPU_EE_40_EE={0}, prompt_norm_noPU_EE_40_EE={0};
  vector<double> prompt_eff_noPU_60_EB={0}, prompt_norm_noPU_60_EB={0}, prompt_eff_noPU_EE_60_EE={0}, prompt_norm_noPU_EE_60_EE={0};
  vector<double> prompt_eff_noPU_80_EB={0}, prompt_norm_noPU_80_EB={0}, prompt_eff_noPU_EE_80_EE={0}, prompt_norm_noPU_EE_80_EE={0};
  vector<double> prompt_eff_noPU_100_EB={0}, prompt_norm_noPU_100_EB={0}, prompt_eff_noPU_EE_100_EE={0}, prompt_norm_noPU_EE_100_EE={0};
      // nonprompt
  vector<double> nonprompt_eff_noPU_EB={0}, nonprompt_norm_noPU_EB={0}, nonprompt_eff_noPU_EE={0}, nonprompt_norm_noPU_EE={0};
  vector<double> nonprompt_eff_noPU_2sigma_EB={0}, nonprompt_norm_noPU_2sigma_EB={0}, nonprompt_eff_noPU_2sigma_EE={0}, nonprompt_norm_noPU_2sigma_EE={0};
  vector<double> nonprompt_eff_noPU_3sigma_EB={0}, nonprompt_norm_noPU_3sigma_EB={0}, nonprompt_eff_noPU_3sigma_EE={0}, nonprompt_norm_noPU_3sigma_EE={0};
  vector<double> nonprompt_eff_noPU_4sigma_EB={0}, nonprompt_norm_noPU_4sigma_EB={0}, nonprompt_eff_noPU_4sigma_EE={0}, nonprompt_norm_noPU_4sigma_EE={0};
  vector<double> nonprompt_eff_noPU_40_EB={0}, nonprompt_norm_noPU_40_EB={0}, nonprompt_eff_noPU_40_EE={0}, nonprompt_norm_noPU_40_EE={0};
  vector<double> nonprompt_eff_noPU_60_EB={0}, nonprompt_norm_noPU_60_EB={0}, nonprompt_eff_noPU_60_EE={0}, nonprompt_norm_noPU_60_EE={0};
  vector<double> nonprompt_eff_noPU_80_EB={0}, nonprompt_norm_noPU_80_EB={0}, nonprompt_eff_noPU_80_EE={0}, nonprompt_norm_noPU_80_EE={0};
  vector<double> nonprompt_eff_noPU_100_EB={0}, nonprompt_norm_noPU_100_EB={0}, nonprompt_eff_noPU_100_EE={0}, nonprompt_norm_noPU_100_EE={0};
      // GEN case
  vector<double> prompt_eff_noPU_gen_EB={0}, prompt_norm_noPU_gen_EB={0}, prompt_eff_noPU_gen_EE={0}, prompt_norm_noPU_gen_EE={0};
  vector<double> nonprompt_eff_noPU_gen_EB={0}, nonprompt_norm_noPU_gen_EB={0}, nonprompt_eff_noPU_gen_EE={0}, nonprompt_norm_noPU_gen_EE={0};


  //////////////////////////
  // Calculate efficiency //
  //////////////////////////
  for(unsigned int i=0; i<h_PU200_prompt_EB->GetNbinsX(); i++) {
    if(zeroiso==true) {        // Include the cases where the muon is already isolated
      // prompt
        // Barrel region
          // efficiency
      prompt_eff_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(1,i+1)/h_PU200_prompt_EB->Integral(1,nbin));
      prompt_eff_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(1,i+1)/h_noPU_prompt_EB->Integral(1,nbin));
      prompt_eff_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(1,i+1)/h_PU200_prompt_40_EB->Integral(1,nbin));
      prompt_eff_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(1,i+1)/h_PU200_prompt_60_EB->Integral(1,nbin));
      prompt_eff_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(1,i+1)/h_PU200_prompt_80_EB->Integral(1,nbin));
      prompt_eff_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(1,i+1)/h_PU200_prompt_100_EB->Integral(1,nbin));
      prompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(1,i+1)/h_PU200_prompt_2sigma_EB->Integral(1,nbin));
      prompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(1,i+1)/h_PU200_prompt_3sigma_EB->Integral(1,nbin));
      prompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(1,i+1)/h_PU200_prompt_4sigma_EB->Integral(1,nbin));
      prompt_eff_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(1,i+1)/h_PU200_prompt_gen_EB->Integral(1,nbin));
      prompt_eff_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(1,i+1)/h_noPU_prompt_gen_EB->Integral(1,nbin));
          // normalization
      prompt_norm_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(1,i+1));
      prompt_norm_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(1,i+1));
      prompt_norm_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(1,i+1));
      prompt_norm_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(1,i+1));
      prompt_norm_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(1,i+1));
      prompt_norm_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(1,i+1));
      prompt_norm_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(1,i+1));
      prompt_norm_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(1,i+1));
      prompt_norm_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(1,i+1));
      prompt_norm_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(1,i+1));
      prompt_norm_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(1,i+1));
        // Endcap region
          // efficiency
      prompt_eff_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(1,i+1)/h_PU200_prompt_EE->Integral(1,nbin));
      prompt_eff_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(1,i+1)/h_noPU_prompt_EE->Integral(1,nbin));
      prompt_eff_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(1,i+1)/h_PU200_prompt_40_EE->Integral(1,nbin));
      prompt_eff_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(1,i+1)/h_PU200_prompt_60_EE->Integral(1,nbin));
      prompt_eff_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(1,i+1)/h_PU200_prompt_80_EE->Integral(1,nbin));
      prompt_eff_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(1,i+1)/h_PU200_prompt_100_EE->Integral(1,nbin));
      prompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(1,i+1)/h_PU200_prompt_2sigma_EE->Integral(1,nbin));
      prompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(1,i+1)/h_PU200_prompt_3sigma_EE->Integral(1,nbin));
      prompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(1,i+1)/h_PU200_prompt_4sigma_EE->Integral(1,nbin));
      prompt_eff_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(1,i+1)/h_PU200_prompt_gen_EE->Integral(1,nbin));
      prompt_eff_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(1,i+1)/h_noPU_prompt_gen_EE->Integral(1,nbin));
          // normalization
      prompt_norm_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(1,i+1));
      prompt_norm_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(1,i+1));
      prompt_norm_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(1,i+1));
      prompt_norm_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(1,i+1));
      prompt_norm_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(1,i+1));
      prompt_norm_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(1,i+1));
      prompt_norm_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(1,i+1));
      prompt_norm_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(1,i+1));
      prompt_norm_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(1,i+1));
      prompt_norm_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(1,i+1));
      prompt_norm_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(1,i+1));
      // nonprompt
        // Barrel region
          // efficiency
      nonprompt_eff_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(1,i+1)/h_PU200_nonprompt_EB->Integral(1,nbin));
      nonprompt_eff_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(1,i+1)/h_noPU_nonprompt_EB->Integral(1,nbin));
      nonprompt_eff_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(1,i+1)/h_PU200_nonprompt_40_EB->Integral(1,nbin));
      nonprompt_eff_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(1,i+1)/h_PU200_nonprompt_60_EB->Integral(1,nbin));
      nonprompt_eff_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(1,i+1)/h_PU200_nonprompt_80_EB->Integral(1,nbin));
      nonprompt_eff_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(1,i+1)/h_PU200_nonprompt_100_EB->Integral(1,nbin));
      nonprompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EB->Integral(1,nbin));
      nonprompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EB->Integral(1,nbin));
      nonprompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EB->Integral(1,nbin));
      nonprompt_eff_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(1,i+1)/h_PU200_nonprompt_gen_EB->Integral(1,nbin));
      nonprompt_eff_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(1,i+1)/h_noPU_nonprompt_gen_EB->Integral(1,nbin));
          // normalization
      nonprompt_norm_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(1,i+1));
      nonprompt_norm_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(1,i+1));
      nonprompt_norm_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(1,i+1));
      nonprompt_norm_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(1,i+1));
      nonprompt_norm_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(1,i+1));
      nonprompt_norm_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(1,i+1));
      nonprompt_norm_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(1,i+1));
      nonprompt_norm_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(1,i+1));
      nonprompt_norm_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(1,i+1));
      nonprompt_norm_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(1,i+1));
      nonprompt_norm_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(1,i+1));
        // Endcap region
          // efficiency
      nonprompt_eff_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(1,i+1)/h_PU200_nonprompt_EE->Integral(1,nbin));
      nonprompt_eff_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(1,i+1)/h_noPU_nonprompt_EE->Integral(1,nbin));
      nonprompt_eff_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(1,i+1)/h_PU200_nonprompt_40_EE->Integral(1,nbin));
      nonprompt_eff_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(1,i+1)/h_PU200_nonprompt_60_EE->Integral(1,nbin));
      nonprompt_eff_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(1,i+1)/h_PU200_nonprompt_80_EE->Integral(1,nbin));
      nonprompt_eff_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(1,i+1)/h_PU200_nonprompt_100_EE->Integral(1,nbin));
      nonprompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EE->Integral(1,nbin));
      nonprompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EE->Integral(1,nbin));
      nonprompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EE->Integral(1,nbin));
      nonprompt_eff_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(1,i+1)/h_PU200_nonprompt_gen_EE->Integral(1,nbin));
      nonprompt_eff_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(1,i+1)/h_noPU_nonprompt_gen_EE->Integral(1,nbin));
          // normalization
      nonprompt_norm_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(1,i+1));
      nonprompt_norm_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(1,i+1));
      nonprompt_norm_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(1,i+1));
      nonprompt_norm_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(1,i+1));
      nonprompt_norm_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(1,i+1));
      nonprompt_norm_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(1,i+1));
      nonprompt_norm_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(1,i+1));
      nonprompt_norm_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(1,i+1));
      nonprompt_norm_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(1,i+1));
      nonprompt_norm_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(1,i+1));
      nonprompt_norm_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(1,i+1));

    }
    else {             // Not include the cases where the muon is already isolated
      // prompt
        // Barrel region
          // efficiency
      prompt_eff_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(2,2+i)/h_PU200_prompt_EB->Integral(2,nbin));
      prompt_eff_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(2,2+i)/h_noPU_prompt_EB->Integral(2,nbin));
      prompt_eff_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(2,2+i)/h_PU200_prompt_40_EB->Integral(2,nbin));
      prompt_eff_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(2,2+i)/h_PU200_prompt_60_EB->Integral(2,nbin));
      prompt_eff_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(2,2+i)/h_PU200_prompt_80_EB->Integral(2,nbin));
      prompt_eff_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(2,2+i)/h_PU200_prompt_100_EB->Integral(2,nbin));
      prompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(2,2+i)/h_PU200_prompt_2sigma_EB->Integral(2,nbin));
      prompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(2,2+i)/h_PU200_prompt_3sigma_EB->Integral(2,nbin));
      prompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(2,2+i)/h_PU200_prompt_4sigma_EB->Integral(2,nbin));
      prompt_eff_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(2,2+i)/h_PU200_prompt_gen_EB->Integral(2,nbin));
      prompt_eff_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(2,2+i)/h_noPU_prompt_gen_EB->Integral(2,nbin));
          // normalization
      prompt_norm_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(2,2+i));
      prompt_norm_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(2,2+i));
      prompt_norm_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(2,2+i));
      prompt_norm_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(2,2+i));
      prompt_norm_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(2,2+i));
      prompt_norm_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(2,2+i));
      prompt_norm_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(2,2+i));
      prompt_norm_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(2,2+i));
      prompt_norm_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(2,2+i));
      prompt_norm_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(2,2+i));
      prompt_norm_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(2,2+i));
        // Endcap region
          // efficiency
      prompt_eff_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(2,2+i)/h_PU200_prompt_EE->Integral(2,nbin));
      prompt_eff_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(2,2+i)/h_noPU_prompt_EE->Integral(2,nbin));
      prompt_eff_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(2,2+i)/h_PU200_prompt_40_EE->Integral(2,nbin));
      prompt_eff_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(2,2+i)/h_PU200_prompt_60_EE->Integral(2,nbin));
      prompt_eff_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(2,2+i)/h_PU200_prompt_80_EE->Integral(2,nbin));
      prompt_eff_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(2,2+i)/h_PU200_prompt_100_EE->Integral(2,nbin));
      prompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(2,2+i)/h_PU200_prompt_2sigma_EE->Integral(2,nbin));
      prompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(2,2+i)/h_PU200_prompt_3sigma_EE->Integral(2,nbin));
      prompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(2,2+i)/h_PU200_prompt_4sigma_EE->Integral(2,nbin));
      prompt_eff_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(2,2+i)/h_PU200_prompt_gen_EE->Integral(2,nbin));
      prompt_eff_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(2,2+i)/h_noPU_prompt_gen_EE->Integral(2,nbin));
          // normalization
      prompt_norm_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(2,2+i));
      prompt_norm_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(2,2+i));
      prompt_norm_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(2,2+i));
      prompt_norm_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(2,2+i));
      prompt_norm_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(2,2+i));
      prompt_norm_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(2,2+i));
      prompt_norm_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(2,2+i));
      prompt_norm_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(2,2+i));
      prompt_norm_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(2,2+i));
      prompt_norm_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(2,2+i));
      prompt_norm_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(2,2+i));
      // nonprompt
        // Barrel region
          // efficiency
      nonprompt_eff_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(2,2+i)/h_PU200_nonprompt_EB->Integral(2,nbin));
      nonprompt_eff_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(2,2+i)/h_noPU_nonprompt_EB->Integral(2,nbin));
      nonprompt_eff_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(2,2+i)/h_PU200_nonprompt_40_EB->Integral(2,nbin));
      nonprompt_eff_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(2,2+i)/h_PU200_nonprompt_60_EB->Integral(2,nbin));
      nonprompt_eff_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(2,2+i)/h_PU200_nonprompt_80_EB->Integral(2,nbin));
      nonprompt_eff_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(2,2+i)/h_PU200_nonprompt_100_EB->Integral(2,nbin));
      nonprompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(2,2+i)/h_PU200_nonprompt_2sigma_EB->Integral(2,nbin));
      nonprompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(2,2+i)/h_PU200_nonprompt_3sigma_EB->Integral(2,nbin));
      nonprompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(2,2+i)/h_PU200_nonprompt_4sigma_EB->Integral(2,nbin));
      nonprompt_eff_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(2,2+i)/h_PU200_nonprompt_gen_EB->Integral(2,nbin));
      nonprompt_eff_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(2,2+i)/h_noPU_nonprompt_gen_EB->Integral(2,nbin));
          // normalization
      nonprompt_norm_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(2,2+i));
      nonprompt_norm_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(2,2+i));
      nonprompt_norm_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(2,2+i));
      nonprompt_norm_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(2,2+i));
      nonprompt_norm_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(2,2+i));
      nonprompt_norm_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(2,2+i));
      nonprompt_norm_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(2,2+i));
      nonprompt_norm_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(2,2+i));
      nonprompt_norm_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(2,2+i));
      nonprompt_norm_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(2,2+i));
      nonprompt_norm_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(2,2+i));
        // Endcap region
          // efficiency
      nonprompt_eff_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(2,2+i)/h_PU200_nonprompt_EE->Integral(2,nbin));
      nonprompt_eff_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(2,2+i)/h_noPU_nonprompt_EE->Integral(2,nbin));
      nonprompt_eff_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(2,2+i)/h_PU200_nonprompt_40_EE->Integral(2,nbin));
      nonprompt_eff_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(2,2+i)/h_PU200_nonprompt_60_EE->Integral(2,nbin));
      nonprompt_eff_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(2,2+i)/h_PU200_nonprompt_80_EE->Integral(2,nbin));
      nonprompt_eff_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(2,2+i)/h_PU200_nonprompt_100_EE->Integral(2,nbin));
      nonprompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(2,2+i)/h_PU200_nonprompt_2sigma_EE->Integral(2,nbin));
      nonprompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(2,2+i)/h_PU200_nonprompt_3sigma_EE->Integral(2,nbin));
      nonprompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(2,2+i)/h_PU200_nonprompt_4sigma_EE->Integral(2,nbin));
      nonprompt_eff_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(2,2+i)/h_PU200_nonprompt_gen_EE->Integral(2,nbin));
      nonprompt_eff_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(2,2+i)/h_noPU_nonprompt_gen_EE->Integral(2,nbin));
          // normalization
      nonprompt_norm_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(2,2+i));
      nonprompt_norm_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(2,2+i));
      nonprompt_norm_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(2,2+i));
      nonprompt_norm_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(2,2+i));
      nonprompt_norm_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(2,2+i));
      nonprompt_norm_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(2,2+i));
      nonprompt_norm_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(2,2+i));
      nonprompt_norm_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(2,2+i));
      nonprompt_norm_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(2,2+i));
      nonprompt_norm_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(2,2+i));
      nonprompt_norm_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(2,2+i));
    }
  }

  // Define TGraph
  // Prompt
    // Barrel region
  TGraph* gr_eff_PU200_prompt_EB = new TGraph();              TGraph* gr_norm_PU200_prompt_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_40_EB = new TGraph();           TGraph* gr_norm_PU200_prompt_40_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_60_EB = new TGraph();           TGraph* gr_norm_PU200_prompt_60_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_80_EB = new TGraph();           TGraph* gr_norm_PU200_prompt_80_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_100_EB = new TGraph();          TGraph* gr_norm_PU200_prompt_100_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_2sigma_EB = new TGraph();       TGraph* gr_norm_PU200_prompt_2sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_3sigma_EB = new TGraph();       TGraph* gr_norm_PU200_prompt_3sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_4sigma_EB = new TGraph();       TGraph* gr_norm_PU200_prompt_4sigma_EB = new TGraph();
  TGraph* gr_eff_noPU_prompt_EB = new TGraph();               TGraph* gr_norm_noPU_prompt_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_gen_EB = new TGraph();          TGraph* gr_norm_PU200_prompt_gen_EB = new TGraph();
  TGraph* gr_eff_noPU_prompt_gen_EB = new TGraph();           TGraph* gr_norm_noPU_prompt_gen_EB = new TGraph();
    // Endcap region
  TGraph* gr_eff_PU200_prompt_EE = new TGraph();              TGraph* gr_norm_PU200_prompt_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_40_EE = new TGraph();           TGraph* gr_norm_PU200_prompt_40_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_60_EE = new TGraph();           TGraph* gr_norm_PU200_prompt_60_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_80_EE = new TGraph();           TGraph* gr_norm_PU200_prompt_80_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_100_EE = new TGraph();          TGraph* gr_norm_PU200_prompt_100_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_2sigma_EE = new TGraph();       TGraph* gr_norm_PU200_prompt_2sigma_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_3sigma_EE = new TGraph();       TGraph* gr_norm_PU200_prompt_3sigma_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_4sigma_EE = new TGraph();       TGraph* gr_norm_PU200_prompt_4sigma_EE = new TGraph();
  TGraph* gr_eff_noPU_prompt_EE = new TGraph();               TGraph* gr_norm_noPU_prompt_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_gen_EE = new TGraph();          TGraph* gr_norm_PU200_prompt_gen_EE = new TGraph();
  TGraph* gr_eff_noPU_prompt_gen_EE = new TGraph();           TGraph* gr_norm_noPU_prompt_gen_EE = new TGraph();
  // Nonprompt
    // Barrel region
  TGraph* gr_eff_PU200_nonprompt_EB = new TGraph();              TGraph* gr_norm_PU200_nonprompt_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_40_EB = new TGraph();           TGraph* gr_norm_PU200_nonprompt_40_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_60_EB = new TGraph();           TGraph* gr_norm_PU200_nonprompt_60_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_80_EB = new TGraph();           TGraph* gr_norm_PU200_nonprompt_80_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_100_EB = new TGraph();          TGraph* gr_norm_PU200_nonprompt_100_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_2sigma_EB = new TGraph();       TGraph* gr_norm_PU200_nonprompt_2sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_3sigma_EB = new TGraph();       TGraph* gr_norm_PU200_nonprompt_3sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_4sigma_EB = new TGraph();       TGraph* gr_norm_PU200_nonprompt_4sigma_EB = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_EB = new TGraph();               TGraph* gr_norm_noPU_nonprompt_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_gen_EB = new TGraph();          TGraph* gr_norm_PU200_nonprompt_gen_EB = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_gen_EB = new TGraph();           TGraph* gr_norm_noPU_nonprompt_gen_EB = new TGraph();
    // Endcap region
  TGraph* gr_eff_PU200_nonprompt_EE = new TGraph();              TGraph* gr_norm_PU200_nonprompt_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_40_EE = new TGraph();           TGraph* gr_norm_PU200_nonprompt_40_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_60_EE = new TGraph();           TGraph* gr_norm_PU200_nonprompt_60_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_80_EE = new TGraph();           TGraph* gr_norm_PU200_nonprompt_80_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_100_EE = new TGraph();          TGraph* gr_norm_PU200_nonprompt_100_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_2sigma_EE = new TGraph();       TGraph* gr_norm_PU200_nonprompt_2sigma_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_3sigma_EE = new TGraph();       TGraph* gr_norm_PU200_nonprompt_3sigma_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_4sigma_EE = new TGraph();       TGraph* gr_norm_PU200_nonprompt_4sigma_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_EE = new TGraph();               TGraph* gr_norm_noPU_nonprompt_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_gen_EE = new TGraph();          TGraph* gr_norm_PU200_nonprompt_gen_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_gen_EE = new TGraph();           TGraph* gr_norm_noPU_nonprompt_gen_EE = new TGraph();


  //for(unsigned int i=0; i<1000; i++) {  // Store efficiency up to rel. iso. cut==4
  for(unsigned int i=0; i<63; i++) {      // Store efficiency up to rel. iso. cut==0.25
  // Prompt
    // Barrel region
      // efficiency
    gr_eff_PU200_prompt_EB->SetPoint(gr_eff_PU200_prompt_EB->GetN(), 0.004*i, prompt_eff_PU200_EB.at(i+1));
    gr_eff_PU200_prompt_40_EB->SetPoint(gr_eff_PU200_prompt_40_EB->GetN(), 0.004*i, prompt_eff_PU200_40_EB.at(i+1));
    gr_eff_PU200_prompt_60_EB->SetPoint(gr_eff_PU200_prompt_60_EB->GetN(), 0.004*i, prompt_eff_PU200_60_EB.at(i+1));
    gr_eff_PU200_prompt_80_EB->SetPoint(gr_eff_PU200_prompt_80_EB->GetN(), 0.004*i, prompt_eff_PU200_80_EB.at(i+1));
    gr_eff_PU200_prompt_100_EB->SetPoint(gr_eff_PU200_prompt_100_EB->GetN(), 0.004*i, prompt_eff_PU200_100_EB.at(i+1));
    gr_eff_PU200_prompt_2sigma_EB->SetPoint(gr_eff_PU200_prompt_2sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EB.at(i+1));
    gr_eff_PU200_prompt_3sigma_EB->SetPoint(gr_eff_PU200_prompt_3sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EB.at(i+1));
    gr_eff_PU200_prompt_4sigma_EB->SetPoint(gr_eff_PU200_prompt_4sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EB.at(i+1));
    gr_eff_noPU_prompt_EB->SetPoint(gr_eff_noPU_prompt_EB->GetN(), 0.004*i, prompt_eff_noPU_EB.at(i+1));
    gr_eff_PU200_prompt_gen_EB->SetPoint(gr_eff_PU200_prompt_gen_EB->GetN(), 0.004*i, prompt_eff_PU200_gen_EB.at(i+1));
    gr_eff_noPU_prompt_gen_EB->SetPoint(gr_eff_noPU_prompt_gen_EB->GetN(), 0.004*i, prompt_eff_noPU_gen_EB.at(i+1));
      // normalization
    gr_norm_PU200_prompt_EB->SetPoint(gr_norm_PU200_prompt_EB->GetN(), 0.004*i, prompt_norm_PU200_EB.at(i+1));
    gr_norm_PU200_prompt_40_EB->SetPoint(gr_norm_PU200_prompt_40_EB->GetN(), 0.004*i, prompt_norm_PU200_40_EB.at(i+1));
    gr_norm_PU200_prompt_60_EB->SetPoint(gr_norm_PU200_prompt_60_EB->GetN(), 0.004*i, prompt_norm_PU200_60_EB.at(i+1));
    gr_norm_PU200_prompt_80_EB->SetPoint(gr_norm_PU200_prompt_80_EB->GetN(), 0.004*i, prompt_norm_PU200_80_EB.at(i+1));
    gr_norm_PU200_prompt_100_EB->SetPoint(gr_norm_PU200_prompt_100_EB->GetN(), 0.004*i, prompt_norm_PU200_100_EB.at(i+1));
    gr_norm_PU200_prompt_2sigma_EB->SetPoint(gr_norm_PU200_prompt_2sigma_EB->GetN(), 0.004*i, prompt_norm_PU200_2sigma_EB.at(i+1));
    gr_norm_PU200_prompt_3sigma_EB->SetPoint(gr_norm_PU200_prompt_3sigma_EB->GetN(), 0.004*i, prompt_norm_PU200_3sigma_EB.at(i+1));
    gr_norm_PU200_prompt_4sigma_EB->SetPoint(gr_norm_PU200_prompt_4sigma_EB->GetN(), 0.004*i, prompt_norm_PU200_4sigma_EB.at(i+1));
    gr_norm_noPU_prompt_EB->SetPoint(gr_norm_noPU_prompt_EB->GetN(), 0.004*i, prompt_norm_noPU_EB.at(i+1));
    gr_norm_PU200_prompt_gen_EB->SetPoint(gr_norm_PU200_prompt_gen_EB->GetN(), 0.004*i, prompt_norm_PU200_gen_EB.at(i+1));
    gr_norm_noPU_prompt_gen_EB->SetPoint(gr_norm_noPU_prompt_gen_EB->GetN(), 0.004*i, prompt_norm_noPU_gen_EB.at(i+1));
    // Endcap region
      // efficiency
    gr_eff_PU200_prompt_EE->SetPoint(gr_eff_PU200_prompt_EE->GetN(), 0.004*i, prompt_eff_PU200_EE.at(i+1));
    gr_eff_PU200_prompt_40_EE->SetPoint(gr_eff_PU200_prompt_40_EE->GetN(), 0.004*i, prompt_eff_PU200_40_EE.at(i+1));
    gr_eff_PU200_prompt_60_EE->SetPoint(gr_eff_PU200_prompt_60_EE->GetN(), 0.004*i, prompt_eff_PU200_60_EE.at(i+1));
    gr_eff_PU200_prompt_80_EE->SetPoint(gr_eff_PU200_prompt_80_EE->GetN(), 0.004*i, prompt_eff_PU200_80_EE.at(i+1));
    gr_eff_PU200_prompt_100_EE->SetPoint(gr_eff_PU200_prompt_100_EE->GetN(), 0.004*i, prompt_eff_PU200_100_EE.at(i+1));
    gr_eff_PU200_prompt_2sigma_EE->SetPoint(gr_eff_PU200_prompt_2sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EE.at(i+1));
    gr_eff_PU200_prompt_3sigma_EE->SetPoint(gr_eff_PU200_prompt_3sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EE.at(i+1));
    gr_eff_PU200_prompt_4sigma_EE->SetPoint(gr_eff_PU200_prompt_4sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EE.at(i+1));
    gr_eff_noPU_prompt_EE->SetPoint(gr_eff_noPU_prompt_EE->GetN(), 0.004*i, prompt_eff_noPU_EE.at(i+1));
    gr_eff_PU200_prompt_gen_EE->SetPoint(gr_eff_PU200_prompt_gen_EE->GetN(), 0.004*i, prompt_eff_PU200_gen_EE.at(i+1));
    gr_eff_noPU_prompt_gen_EE->SetPoint(gr_eff_noPU_prompt_gen_EE->GetN(), 0.004*i, prompt_eff_noPU_gen_EE.at(i+1));
      // normalization
    gr_norm_PU200_prompt_EE->SetPoint(gr_norm_PU200_prompt_EE->GetN(), 0.004*i, prompt_norm_PU200_EE.at(i+1));
    gr_norm_PU200_prompt_40_EE->SetPoint(gr_norm_PU200_prompt_40_EE->GetN(), 0.004*i, prompt_norm_PU200_40_EE.at(i+1));
    gr_norm_PU200_prompt_60_EE->SetPoint(gr_norm_PU200_prompt_60_EE->GetN(), 0.004*i, prompt_norm_PU200_60_EE.at(i+1));
    gr_norm_PU200_prompt_80_EE->SetPoint(gr_norm_PU200_prompt_80_EE->GetN(), 0.004*i, prompt_norm_PU200_80_EE.at(i+1));
    gr_norm_PU200_prompt_100_EE->SetPoint(gr_norm_PU200_prompt_100_EE->GetN(), 0.004*i, prompt_norm_PU200_100_EE.at(i+1));
    gr_norm_PU200_prompt_2sigma_EE->SetPoint(gr_norm_PU200_prompt_2sigma_EE->GetN(), 0.004*i, prompt_norm_PU200_2sigma_EE.at(i+1));
    gr_norm_PU200_prompt_3sigma_EE->SetPoint(gr_norm_PU200_prompt_3sigma_EE->GetN(), 0.004*i, prompt_norm_PU200_3sigma_EE.at(i+1));
    gr_norm_PU200_prompt_4sigma_EE->SetPoint(gr_norm_PU200_prompt_4sigma_EE->GetN(), 0.004*i, prompt_norm_PU200_4sigma_EE.at(i+1));
    gr_norm_noPU_prompt_EE->SetPoint(gr_norm_noPU_prompt_EE->GetN(), 0.004*i, prompt_norm_noPU_EE.at(i+1));
    gr_norm_PU200_prompt_gen_EE->SetPoint(gr_norm_PU200_prompt_gen_EE->GetN(), 0.004*i, prompt_norm_PU200_gen_EE.at(i+1));
    gr_norm_noPU_prompt_gen_EE->SetPoint(gr_norm_noPU_prompt_gen_EE->GetN(), 0.004*i, prompt_norm_noPU_gen_EE.at(i+1));
  }
  //for(unsigned int i=0; i<1000; i++) {  // Store efficiency up to rel. iso. cut==4
  for(unsigned int i=0; i<63; i++) {      // Store efficiency up to rel. iso. cut==0.25
  // Nonprompt
    // Barrel region
      // efficiency
    gr_eff_PU200_nonprompt_EB->SetPoint(gr_eff_PU200_nonprompt_EB->GetN(), 0.004*i, nonprompt_eff_PU200_EB.at(i+1));
    gr_eff_PU200_nonprompt_40_EB->SetPoint(gr_eff_PU200_nonprompt_40_EB->GetN(), 0.004*i, nonprompt_eff_PU200_40_EB.at(i+1));
    gr_eff_PU200_nonprompt_60_EB->SetPoint(gr_eff_PU200_nonprompt_60_EB->GetN(), 0.004*i, nonprompt_eff_PU200_60_EB.at(i+1));
    gr_eff_PU200_nonprompt_80_EB->SetPoint(gr_eff_PU200_nonprompt_80_EB->GetN(), 0.004*i, nonprompt_eff_PU200_80_EB.at(i+1));
    gr_eff_PU200_nonprompt_100_EB->SetPoint(gr_eff_PU200_nonprompt_100_EB->GetN(), 0.004*i, nonprompt_eff_PU200_100_EB.at(i+1));
    gr_eff_PU200_nonprompt_2sigma_EB->SetPoint(gr_eff_PU200_nonprompt_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EB.at(i+1));
    gr_eff_PU200_nonprompt_3sigma_EB->SetPoint(gr_eff_PU200_nonprompt_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EB.at(i+1));
    gr_eff_PU200_nonprompt_4sigma_EB->SetPoint(gr_eff_PU200_nonprompt_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EB.at(i+1));
    gr_eff_noPU_nonprompt_EB->SetPoint(gr_eff_noPU_nonprompt_EB->GetN(), 0.004*i, nonprompt_eff_noPU_EB.at(i+1));
    gr_eff_PU200_nonprompt_gen_EB->SetPoint(gr_eff_PU200_nonprompt_gen_EB->GetN(), 0.004*i, nonprompt_eff_PU200_gen_EB.at(i+1));
    gr_eff_noPU_nonprompt_gen_EB->SetPoint(gr_eff_noPU_nonprompt_gen_EB->GetN(), 0.004*i, nonprompt_eff_noPU_gen_EB.at(i+1));
      // normalization
    gr_norm_PU200_nonprompt_EB->SetPoint(gr_norm_PU200_nonprompt_EB->GetN(), 0.004*i, nonprompt_norm_PU200_EB.at(i+1));
    gr_norm_PU200_nonprompt_40_EB->SetPoint(gr_norm_PU200_nonprompt_40_EB->GetN(), 0.004*i, nonprompt_norm_PU200_40_EB.at(i+1));
    gr_norm_PU200_nonprompt_60_EB->SetPoint(gr_norm_PU200_nonprompt_60_EB->GetN(), 0.004*i, nonprompt_norm_PU200_60_EB.at(i+1));
    gr_norm_PU200_nonprompt_80_EB->SetPoint(gr_norm_PU200_nonprompt_80_EB->GetN(), 0.004*i, nonprompt_norm_PU200_80_EB.at(i+1));
    gr_norm_PU200_nonprompt_100_EB->SetPoint(gr_norm_PU200_nonprompt_100_EB->GetN(), 0.004*i, nonprompt_norm_PU200_100_EB.at(i+1));
    gr_norm_PU200_nonprompt_2sigma_EB->SetPoint(gr_norm_PU200_nonprompt_2sigma_EB->GetN(), 0.004*i, nonprompt_norm_PU200_2sigma_EB.at(i+1));
    gr_norm_PU200_nonprompt_3sigma_EB->SetPoint(gr_norm_PU200_nonprompt_3sigma_EB->GetN(), 0.004*i, nonprompt_norm_PU200_3sigma_EB.at(i+1));
    gr_norm_PU200_nonprompt_4sigma_EB->SetPoint(gr_norm_PU200_nonprompt_4sigma_EB->GetN(), 0.004*i, nonprompt_norm_PU200_4sigma_EB.at(i+1));
    gr_norm_noPU_nonprompt_EB->SetPoint(gr_norm_noPU_nonprompt_EB->GetN(), 0.004*i, nonprompt_norm_noPU_EB.at(i+1));
    gr_norm_PU200_nonprompt_gen_EB->SetPoint(gr_norm_PU200_nonprompt_gen_EB->GetN(), 0.004*i, nonprompt_norm_PU200_gen_EB.at(i+1));
    gr_norm_noPU_nonprompt_gen_EB->SetPoint(gr_norm_noPU_nonprompt_gen_EB->GetN(), 0.004*i, nonprompt_norm_noPU_gen_EB.at(i+1));
    // Endcap region
      // efficiency
    gr_eff_PU200_nonprompt_EE->SetPoint(gr_eff_PU200_nonprompt_EE->GetN(), 0.004*i, nonprompt_eff_PU200_EE.at(i+1));
    gr_eff_PU200_nonprompt_40_EE->SetPoint(gr_eff_PU200_nonprompt_40_EE->GetN(), 0.004*i, nonprompt_eff_PU200_40_EE.at(i+1));
    gr_eff_PU200_nonprompt_60_EE->SetPoint(gr_eff_PU200_nonprompt_60_EE->GetN(), 0.004*i, nonprompt_eff_PU200_60_EE.at(i+1));
    gr_eff_PU200_nonprompt_80_EE->SetPoint(gr_eff_PU200_nonprompt_80_EE->GetN(), 0.004*i, nonprompt_eff_PU200_80_EE.at(i+1));
    gr_eff_PU200_nonprompt_100_EE->SetPoint(gr_eff_PU200_nonprompt_100_EE->GetN(), 0.004*i, nonprompt_eff_PU200_100_EE.at(i+1));
    gr_eff_PU200_nonprompt_2sigma_EE->SetPoint(gr_eff_PU200_nonprompt_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EE.at(i+1));
    gr_eff_PU200_nonprompt_3sigma_EE->SetPoint(gr_eff_PU200_nonprompt_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EE.at(i+1));
    gr_eff_PU200_nonprompt_4sigma_EE->SetPoint(gr_eff_PU200_nonprompt_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EE.at(i+1));
    gr_eff_noPU_nonprompt_EE->SetPoint(gr_eff_noPU_nonprompt_EE->GetN(), 0.004*i, nonprompt_eff_noPU_EE.at(i+1));
    gr_eff_PU200_nonprompt_gen_EE->SetPoint(gr_eff_PU200_nonprompt_gen_EE->GetN(), 0.004*i, nonprompt_eff_PU200_gen_EE.at(i+1));
    gr_eff_noPU_nonprompt_gen_EE->SetPoint(gr_eff_noPU_nonprompt_gen_EE->GetN(), 0.004*i, nonprompt_eff_noPU_gen_EE.at(i+1));
      // normalization
    gr_norm_PU200_nonprompt_EE->SetPoint(gr_norm_PU200_nonprompt_EE->GetN(), 0.004*i, nonprompt_norm_PU200_EE.at(i+1));
    gr_norm_PU200_nonprompt_40_EE->SetPoint(gr_norm_PU200_nonprompt_40_EE->GetN(), 0.004*i, nonprompt_norm_PU200_40_EE.at(i+1));
    gr_norm_PU200_nonprompt_60_EE->SetPoint(gr_norm_PU200_nonprompt_60_EE->GetN(), 0.004*i, nonprompt_norm_PU200_60_EE.at(i+1));
    gr_norm_PU200_nonprompt_80_EE->SetPoint(gr_norm_PU200_nonprompt_80_EE->GetN(), 0.004*i, nonprompt_norm_PU200_80_EE.at(i+1));
    gr_norm_PU200_nonprompt_100_EE->SetPoint(gr_norm_PU200_nonprompt_100_EE->GetN(), 0.004*i, nonprompt_norm_PU200_100_EE.at(i+1));
    gr_norm_PU200_nonprompt_2sigma_EE->SetPoint(gr_norm_PU200_nonprompt_2sigma_EE->GetN(), 0.004*i, nonprompt_norm_PU200_2sigma_EE.at(i+1));
    gr_norm_PU200_nonprompt_3sigma_EE->SetPoint(gr_norm_PU200_nonprompt_3sigma_EE->GetN(), 0.004*i, nonprompt_norm_PU200_3sigma_EE.at(i+1));
    gr_norm_PU200_nonprompt_4sigma_EE->SetPoint(gr_norm_PU200_nonprompt_4sigma_EE->GetN(), 0.004*i, nonprompt_norm_PU200_4sigma_EE.at(i+1));
    gr_norm_noPU_nonprompt_EE->SetPoint(gr_norm_noPU_nonprompt_EE->GetN(), 0.004*i, nonprompt_norm_noPU_EE.at(i+1));
    gr_norm_PU200_nonprompt_gen_EE->SetPoint(gr_norm_PU200_nonprompt_gen_EE->GetN(), 0.004*i, nonprompt_norm_PU200_gen_EE.at(i+1));
    gr_norm_noPU_nonprompt_gen_EE->SetPoint(gr_norm_noPU_nonprompt_gen_EE->GetN(), 0.004*i, nonprompt_norm_noPU_gen_EE.at(i+1));

  }

  ///////////////
  // Cosmetics //
  ///////////////
  // Prompt
    // Barrel region
      // efficiency
  gr_eff_PU200_prompt_EB->SetTitle("Prompt muon in barrel region"); gr_eff_noPU_prompt_EB->SetTitle("Prompt muon in barrel region");
  gr_eff_PU200_prompt_EB->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_prompt_EB->GetXaxis()->SetTitle("Relative isolation");
  gr_eff_PU200_prompt_EB->GetYaxis()->SetTitle("Isolation efficiency"); gr_eff_noPU_prompt_EB->GetYaxis()->SetTitle("Isolation efficiency");
  gr_eff_PU200_prompt_EB->SetLineColor(kBlack); gr_eff_PU200_prompt_gen_EB->SetLineColor(kBlack);
  gr_eff_PU200_prompt_40_EB->SetLineColor(kRed); gr_eff_PU200_prompt_60_EB->SetLineColor(kGreen); gr_eff_PU200_prompt_80_EB->SetLineColor(kBlue); gr_eff_PU200_prompt_100_EB->SetLineColor(kMagenta);
  gr_eff_PU200_prompt_2sigma_EB->SetLineColor(kRed); gr_eff_PU200_prompt_3sigma_EB->SetLineColor(kGreen); gr_eff_PU200_prompt_4sigma_EB->SetLineColor(kBlue);
  gr_eff_noPU_prompt_EB->SetLineColor(kGray); gr_eff_noPU_prompt_gen_EB->SetLineColor(kGray);
  gr_eff_PU200_prompt_EB->SetLineWidth(2); gr_eff_PU200_prompt_gen_EB->SetLineWidth(2);
  gr_eff_PU200_prompt_40_EB->SetLineWidth(2); gr_eff_PU200_prompt_60_EB->SetLineWidth(2); gr_eff_PU200_prompt_80_EB->SetLineWidth(2); gr_eff_PU200_prompt_100_EB->SetLineWidth(2);
  gr_eff_PU200_prompt_2sigma_EB->SetLineWidth(2); gr_eff_PU200_prompt_3sigma_EB->SetLineWidth(2); gr_eff_PU200_prompt_4sigma_EB->SetLineWidth(2);
  gr_eff_noPU_prompt_EB->SetLineWidth(2); gr_eff_noPU_prompt_gen_EB->SetLineWidth(2);
  gr_eff_PU200_prompt_gen_EB->SetLineStyle(7); gr_eff_noPU_prompt_gen_EB->SetLineStyle(7);
      // normalization
  gr_norm_PU200_prompt_EB->SetTitle("Prompt muon in barrel region"); gr_norm_noPU_prompt_EB->SetTitle("Prompt muon in barrel region");
  gr_norm_PU200_prompt_EB->GetXaxis()->SetTitle("Relative isolation"); gr_norm_noPU_prompt_EB->GetXaxis()->SetTitle("Relative isolation");
  gr_norm_PU200_prompt_EB->GetYaxis()->SetTitle("Counts / 0.004"); gr_norm_noPU_prompt_EB->GetYaxis()->SetTitle("Counts / 0.004");
  gr_norm_PU200_prompt_EB->SetLineColor(kBlack); gr_norm_PU200_prompt_gen_EB->SetLineColor(kBlack);
  gr_norm_PU200_prompt_40_EB->SetLineColor(kRed); gr_norm_PU200_prompt_60_EB->SetLineColor(kGreen); gr_norm_PU200_prompt_80_EB->SetLineColor(kBlue); gr_norm_PU200_prompt_100_EB->SetLineColor(kMagenta);
  gr_norm_PU200_prompt_2sigma_EB->SetLineColor(kRed); gr_norm_PU200_prompt_3sigma_EB->SetLineColor(kGreen); gr_norm_PU200_prompt_4sigma_EB->SetLineColor(kBlue);
  gr_norm_noPU_prompt_EB->SetLineColor(kGray); gr_norm_noPU_prompt_gen_EB->SetLineColor(kGray);
  gr_norm_PU200_prompt_EB->SetLineWidth(2); gr_norm_PU200_prompt_gen_EB->SetLineWidth(2);
  gr_norm_PU200_prompt_40_EB->SetLineWidth(2); gr_norm_PU200_prompt_60_EB->SetLineWidth(2); gr_norm_PU200_prompt_80_EB->SetLineWidth(2); gr_norm_PU200_prompt_100_EB->SetLineWidth(2);
  gr_norm_PU200_prompt_2sigma_EB->SetLineWidth(2); gr_norm_PU200_prompt_3sigma_EB->SetLineWidth(2); gr_norm_PU200_prompt_4sigma_EB->SetLineWidth(2);
  gr_norm_noPU_prompt_EB->SetLineWidth(2); gr_norm_noPU_prompt_gen_EB->SetLineWidth(2);
  gr_norm_PU200_prompt_gen_EB->SetLineStyle(7); gr_norm_noPU_prompt_gen_EB->SetLineStyle(7);
  gr_norm_noPU_prompt_EB->SetMinimum(0);
    // Endcap region
      // efficiency
  gr_eff_PU200_prompt_EE->SetTitle("Prompt muon in endcap region"); gr_eff_noPU_prompt_EE->SetTitle("Prompt muon in endcap region");
  gr_eff_PU200_prompt_EE->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_prompt_EE->GetXaxis()->SetTitle("Relative isolation");
  gr_eff_PU200_prompt_EE->GetYaxis()->SetTitle("Isolation efficiency"); gr_eff_noPU_prompt_EE->GetYaxis()->SetTitle("Isolation efficiency");
  gr_eff_PU200_prompt_EE->SetLineColor(kBlack); gr_eff_PU200_prompt_gen_EE->SetLineColor(kBlack);
  gr_eff_PU200_prompt_40_EE->SetLineColor(kRed); gr_eff_PU200_prompt_60_EE->SetLineColor(kGreen); gr_eff_PU200_prompt_80_EE->SetLineColor(kBlue); gr_eff_PU200_prompt_100_EE->SetLineColor(kMagenta);
  gr_eff_PU200_prompt_2sigma_EE->SetLineColor(kRed); gr_eff_PU200_prompt_3sigma_EE->SetLineColor(kGreen); gr_eff_PU200_prompt_4sigma_EE->SetLineColor(kBlue);
  gr_eff_noPU_prompt_EE->SetLineColor(kGray); gr_eff_noPU_prompt_gen_EE->SetLineColor(kGray);
  gr_eff_PU200_prompt_EE->SetLineWidth(2); gr_eff_PU200_prompt_gen_EE->SetLineWidth(2);
  gr_eff_PU200_prompt_40_EE->SetLineWidth(2); gr_eff_PU200_prompt_60_EE->SetLineWidth(2); gr_eff_PU200_prompt_80_EE->SetLineWidth(2); gr_eff_PU200_prompt_100_EE->SetLineWidth(2);
  gr_eff_PU200_prompt_2sigma_EE->SetLineWidth(2); gr_eff_PU200_prompt_3sigma_EE->SetLineWidth(2); gr_eff_PU200_prompt_4sigma_EE->SetLineWidth(2);
  gr_eff_noPU_prompt_EE->SetLineWidth(2); gr_eff_noPU_prompt_gen_EE->SetLineWidth(2);
  gr_eff_PU200_prompt_gen_EE->SetLineStyle(7); gr_eff_noPU_prompt_gen_EE->SetLineStyle(7);
      // normalization
  gr_norm_PU200_prompt_EE->SetTitle("Prompt muon in endcap region"); gr_norm_noPU_prompt_EE->SetTitle("Prompt muon in endcap region");
  gr_norm_PU200_prompt_EE->GetXaxis()->SetTitle("Relative isolation"); gr_norm_noPU_prompt_EE->GetXaxis()->SetTitle("Relative isolation");
  gr_norm_PU200_prompt_EE->GetYaxis()->SetTitle("Counts / 0.004"); gr_norm_noPU_prompt_EE->GetYaxis()->SetTitle("Counts / 0.004");
  gr_norm_PU200_prompt_EE->SetLineColor(kBlack); gr_norm_PU200_prompt_gen_EE->SetLineColor(kBlack);
  gr_norm_PU200_prompt_40_EE->SetLineColor(kRed); gr_norm_PU200_prompt_60_EE->SetLineColor(kGreen); gr_norm_PU200_prompt_80_EE->SetLineColor(kBlue); gr_norm_PU200_prompt_100_EE->SetLineColor(kMagenta);
  gr_norm_PU200_prompt_2sigma_EE->SetLineColor(kRed); gr_norm_PU200_prompt_3sigma_EE->SetLineColor(kGreen); gr_norm_PU200_prompt_4sigma_EE->SetLineColor(kBlue);
  gr_norm_noPU_prompt_EE->SetLineColor(kGray); gr_norm_noPU_prompt_gen_EE->SetLineColor(kGray);
  gr_norm_PU200_prompt_EE->SetLineWidth(2); gr_norm_PU200_prompt_gen_EE->SetLineWidth(2);
  gr_norm_PU200_prompt_40_EE->SetLineWidth(2); gr_norm_PU200_prompt_60_EE->SetLineWidth(2); gr_norm_PU200_prompt_80_EE->SetLineWidth(2); gr_norm_PU200_prompt_100_EE->SetLineWidth(2);
  gr_norm_PU200_prompt_2sigma_EE->SetLineWidth(2); gr_norm_PU200_prompt_3sigma_EE->SetLineWidth(2); gr_norm_PU200_prompt_4sigma_EE->SetLineWidth(2);
  gr_norm_noPU_prompt_EE->SetLineWidth(2); gr_norm_noPU_prompt_gen_EE->SetLineWidth(2);
  gr_norm_PU200_prompt_gen_EE->SetLineStyle(7); gr_norm_noPU_prompt_gen_EE->SetLineStyle(7);
  gr_norm_noPU_prompt_EE->SetMinimum(0);
  // Nonprompt
    // Barrel region
      // efficiency
  gr_eff_PU200_nonprompt_EB->SetTitle("Non-prompt muon in barrel region"); gr_eff_noPU_nonprompt_EB->SetTitle("Non-prompt muon in barrel region");
  gr_eff_PU200_nonprompt_EB->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_nonprompt_EB->GetXaxis()->SetTitle("Relative isolation");
  gr_eff_PU200_nonprompt_EB->GetYaxis()->SetTitle("Isolation efficiency"); gr_eff_noPU_nonprompt_EB->GetYaxis()->SetTitle("Isolation efficiency");
  gr_eff_PU200_nonprompt_EB->SetLineColor(kBlack); gr_eff_PU200_nonprompt_gen_EB->SetLineColor(kBlack);
  gr_eff_PU200_nonprompt_40_EB->SetLineColor(kRed); gr_eff_PU200_nonprompt_60_EB->SetLineColor(kGreen); gr_eff_PU200_nonprompt_80_EB->SetLineColor(kBlue); gr_eff_PU200_nonprompt_100_EB->SetLineColor(kMagenta);
  gr_eff_PU200_nonprompt_2sigma_EB->SetLineColor(kRed); gr_eff_PU200_nonprompt_3sigma_EB->SetLineColor(kGreen); gr_eff_PU200_nonprompt_4sigma_EB->SetLineColor(kBlue);
  gr_eff_noPU_nonprompt_EB->SetLineColor(kGray); gr_eff_noPU_nonprompt_gen_EB->SetLineColor(kGray);
  gr_eff_PU200_nonprompt_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_gen_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_40_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_60_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_80_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_100_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_2sigma_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_3sigma_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_4sigma_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_EB->SetLineWidth(2); gr_eff_noPU_nonprompt_gen_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_EB->SetLineStyle(7); gr_eff_noPU_nonprompt_gen_EB->SetLineStyle(7);
      // normalization
  gr_norm_PU200_nonprompt_EB->SetTitle("Non-prompt muon in barrel region"); gr_norm_noPU_nonprompt_EB->SetTitle("Non-prompt muon in barrel region");
  gr_norm_PU200_nonprompt_EB->GetXaxis()->SetTitle("Relative isolation"); gr_norm_noPU_nonprompt_EB->GetXaxis()->SetTitle("Relative isolation");
  gr_norm_PU200_nonprompt_EB->GetYaxis()->SetTitle("Counts / 0.004"); gr_norm_noPU_nonprompt_EB->GetYaxis()->SetTitle("Counts / 0.004");
  gr_norm_PU200_nonprompt_EB->SetLineColor(kBlack); gr_norm_PU200_nonprompt_gen_EB->SetLineColor(kBlack);
  gr_norm_PU200_nonprompt_40_EB->SetLineColor(kRed); gr_norm_PU200_nonprompt_60_EB->SetLineColor(kGreen); gr_norm_PU200_nonprompt_80_EB->SetLineColor(kBlue); gr_norm_PU200_nonprompt_100_EB->SetLineColor(kMagenta);
  gr_norm_PU200_nonprompt_2sigma_EB->SetLineColor(kRed); gr_norm_PU200_nonprompt_3sigma_EB->SetLineColor(kGreen); gr_norm_PU200_nonprompt_4sigma_EB->SetLineColor(kBlue);
  gr_norm_noPU_nonprompt_EB->SetLineColor(kGray); gr_norm_noPU_nonprompt_gen_EB->SetLineColor(kGray);
  gr_norm_PU200_nonprompt_EB->SetLineWidth(2); gr_norm_PU200_nonprompt_gen_EB->SetLineWidth(2);
  gr_norm_PU200_nonprompt_40_EB->SetLineWidth(2); gr_norm_PU200_nonprompt_60_EB->SetLineWidth(2); gr_norm_PU200_nonprompt_80_EB->SetLineWidth(2); gr_norm_PU200_nonprompt_100_EB->SetLineWidth(2);
  gr_norm_PU200_nonprompt_2sigma_EB->SetLineWidth(2); gr_norm_PU200_nonprompt_3sigma_EB->SetLineWidth(2); gr_norm_PU200_nonprompt_4sigma_EB->SetLineWidth(2);
  gr_norm_PU200_nonprompt_gen_EB->SetLineWidth(2); gr_norm_noPU_nonprompt_gen_EB->SetLineWidth(2);
  gr_norm_PU200_nonprompt_gen_EB->SetLineStyle(7); gr_norm_noPU_nonprompt_gen_EB->SetLineStyle(7);
  gr_norm_noPU_nonprompt_EB->SetMinimum(0);
    // Endcap region
      // efficiency
  gr_eff_PU200_nonprompt_EE->SetTitle("Non-prompt muon in endcap region"); gr_eff_noPU_nonprompt_EE->SetTitle("Non-prompt muon in endcap region");
  gr_eff_PU200_nonprompt_EE->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_nonprompt_EE->GetXaxis()->SetTitle("Relative isolation");
  gr_eff_PU200_nonprompt_EE->GetYaxis()->SetTitle("Isolation efficiency"); gr_eff_noPU_nonprompt_EE->GetYaxis()->SetTitle("Isolation efficiency");
  gr_eff_PU200_nonprompt_EE->SetLineColor(kBlack); gr_eff_PU200_nonprompt_gen_EE->SetLineColor(kBlack);
  gr_eff_PU200_nonprompt_40_EE->SetLineColor(kRed); gr_eff_PU200_nonprompt_60_EE->SetLineColor(kGreen); gr_eff_PU200_nonprompt_80_EE->SetLineColor(kBlue); gr_eff_PU200_nonprompt_100_EE->SetLineColor(kMagenta);
  gr_eff_PU200_nonprompt_2sigma_EE->SetLineColor(kRed); gr_eff_PU200_nonprompt_3sigma_EE->SetLineColor(kGreen); gr_eff_PU200_nonprompt_4sigma_EE->SetLineColor(kBlue);
  gr_eff_noPU_nonprompt_EE->SetLineColor(kGray); gr_eff_noPU_nonprompt_gen_EE->SetLineColor(kGray);
  gr_eff_PU200_nonprompt_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_gen_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_40_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_60_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_80_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_100_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_2sigma_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_3sigma_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_4sigma_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_EE->SetLineWidth(2); gr_eff_noPU_nonprompt_gen_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_EE->SetLineStyle(7); gr_eff_noPU_nonprompt_gen_EE->SetLineStyle(7);
      // normalization
  gr_norm_PU200_nonprompt_EE->SetTitle("Non-prompt muon in endcap region"); gr_norm_noPU_nonprompt_EE->SetTitle("Non-prompt muon in endcap region");
  gr_norm_PU200_nonprompt_EE->GetXaxis()->SetTitle("Relative isolation"); gr_norm_noPU_nonprompt_EE->GetXaxis()->SetTitle("Relative isolation");
  gr_norm_PU200_nonprompt_EE->GetYaxis()->SetTitle("Counts / 0.004"); gr_norm_noPU_nonprompt_EE->GetYaxis()->SetTitle("Counts / 0.004");
  gr_norm_PU200_nonprompt_EE->SetLineColor(kBlack); gr_norm_PU200_nonprompt_gen_EE->SetLineColor(kBlack);
  gr_norm_PU200_nonprompt_40_EE->SetLineColor(kRed); gr_norm_PU200_nonprompt_60_EE->SetLineColor(kGreen); gr_norm_PU200_nonprompt_80_EE->SetLineColor(kBlue); gr_norm_PU200_nonprompt_100_EE->SetLineColor(kMagenta);
  gr_norm_PU200_nonprompt_2sigma_EE->SetLineColor(kRed); gr_norm_PU200_nonprompt_3sigma_EE->SetLineColor(kGreen); gr_norm_PU200_nonprompt_4sigma_EE->SetLineColor(kBlue);
  gr_norm_noPU_nonprompt_EE->SetLineColor(kGray); gr_norm_noPU_nonprompt_gen_EE->SetLineColor(kGray);
  gr_norm_PU200_nonprompt_EE->SetLineWidth(2); gr_norm_PU200_nonprompt_gen_EE->SetLineWidth(2);
  gr_norm_PU200_nonprompt_40_EE->SetLineWidth(2); gr_norm_PU200_nonprompt_60_EE->SetLineWidth(2); gr_norm_PU200_nonprompt_80_EE->SetLineWidth(2); gr_norm_PU200_nonprompt_100_EE->SetLineWidth(2);
  gr_norm_PU200_nonprompt_2sigma_EE->SetLineWidth(2); gr_norm_PU200_nonprompt_3sigma_EE->SetLineWidth(2); gr_norm_PU200_nonprompt_4sigma_EE->SetLineWidth(2);
  gr_norm_PU200_nonprompt_gen_EE->SetLineWidth(2); gr_norm_noPU_nonprompt_gen_EE->SetLineWidth(2);
  gr_norm_PU200_nonprompt_gen_EE->SetLineStyle(7); gr_norm_noPU_nonprompt_gen_EE->SetLineStyle(7);
  gr_norm_noPU_nonprompt_EE->SetMinimum(0);

  /////////////
  // Legends //
  /////////////
  // Prompt
    // Barrel region
      // efficiency
  TLegend* leg_eff_prompt_EB_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_gen_EB, "gen PU200");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_EB, "no MTD PU200");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_40_EB, "40ps PU200");
//  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_60_EB, "60ps PU200");
//  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_80_EB, "80ps PU200");
//  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_100_EB, "100ps PU200");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_noPU_prompt_gen_EB, "gen noPU");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_noPU_prompt_EB, "noPU");
  leg_eff_prompt_EB_dt->SetTextSize(0.03);
  TLegend* leg_eff_prompt_EB_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_gen_EB, "gen PU200");
  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_EB, "no MTD PU200");
  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_2sigma_EB, "2sigma PU200");
//  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_3sigma_EB, "3sigma PU200");
//  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_4sigma_EB, "4sigma PU200");
  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_noPU_prompt_gen_EB, "gen noPU");
  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_noPU_prompt_EB, "noPU");
  leg_eff_prompt_EB_sigma->SetTextSize(0.03);
      // normalization
  TLegend* leg_norm_prompt_EB_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_prompt_EB_dt->AddEntry(gr_norm_PU200_prompt_gen_EB, "gen PU200");
  leg_norm_prompt_EB_dt->AddEntry(gr_norm_PU200_prompt_EB, "no MTD PU200");
  leg_norm_prompt_EB_dt->AddEntry(gr_norm_PU200_prompt_40_EB, "40ps PU200");
//  leg_norm_prompt_EB_dt->AddEntry(gr_norm_PU200_prompt_60_EB, "60ps PU200");
//  leg_norm_prompt_EB_dt->AddEntry(gr_norm_PU200_prompt_80_EB, "80ps PU200");
//  leg_norm_prompt_EB_dt->AddEntry(gr_norm_PU200_prompt_100_EB, "100ps PU200");
  leg_norm_prompt_EB_dt->AddEntry(gr_norm_noPU_prompt_gen_EB, "gen noPU");
  leg_norm_prompt_EB_dt->AddEntry(gr_norm_noPU_prompt_EB, "noPU");
  leg_norm_prompt_EB_dt->SetTextSize(0.03);
  TLegend* leg_norm_prompt_EB_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_prompt_EB_sigma->AddEntry(gr_norm_PU200_prompt_gen_EB, "gen PU200");
  leg_norm_prompt_EB_sigma->AddEntry(gr_norm_PU200_prompt_EB, "no MTD PU200");
  leg_norm_prompt_EB_sigma->AddEntry(gr_norm_PU200_prompt_2sigma_EB, "2sigma PU200");
//  leg_norm_prompt_EB_sigma->AddEntry(gr_norm_PU200_prompt_3sigma_EB, "3sigma PU200");
//  leg_norm_prompt_EB_sigma->AddEntry(gr_norm_PU200_prompt_4sigma_EB, "4sigma PU200");
  leg_norm_prompt_EB_sigma->AddEntry(gr_norm_noPU_prompt_gen_EB, "gen noPU");
  leg_norm_prompt_EB_sigma->AddEntry(gr_norm_noPU_prompt_EB, "noPU");
  leg_norm_prompt_EB_sigma->SetTextSize(0.03);
    // Endcap region
      // efficiency
  TLegend* leg_eff_prompt_EE_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_gen_EE, "gen PU200");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_EE, "no MTD PU200");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_40_EE, "40ps PU200");
//  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_60_EE, "60ps PU200");
//  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_80_EE, "80ps PU200");
//  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_100_EE, "100ps PU200");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_noPU_prompt_gen_EE, "gen noPU");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_noPU_prompt_EE, "noPU");
  leg_eff_prompt_EE_dt->SetTextSize(0.03);
  TLegend* leg_eff_prompt_EE_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_gen_EE, "gen PU200");
  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_EE, "no MTD PU200");
  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_2sigma_EE, "2sigma PU200");
//  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_3sigma_EE, "3sigma PU200");
//  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_4sigma_EE, "4sigma PU200");
  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_noPU_prompt_gen_EE, "gen noPU");
  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_noPU_prompt_EE, "noPU");
  leg_eff_prompt_EE_sigma->SetTextSize(0.03);
      // normalization
  TLegend* leg_norm_prompt_EE_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_prompt_EE_dt->AddEntry(gr_norm_PU200_prompt_gen_EE, "gen PU200");
  leg_norm_prompt_EE_dt->AddEntry(gr_norm_PU200_prompt_EE, "no MTD PU200");
  leg_norm_prompt_EE_dt->AddEntry(gr_norm_PU200_prompt_40_EE, "40ps PU200");
//  leg_norm_prompt_EE_dt->AddEntry(gr_norm_PU200_prompt_60_EE, "60ps PU200");
//  leg_norm_prompt_EE_dt->AddEntry(gr_norm_PU200_prompt_80_EE, "80ps PU200");
//  leg_norm_prompt_EE_dt->AddEntry(gr_norm_PU200_prompt_100_EE, "100ps PU200");
  leg_norm_prompt_EE_dt->AddEntry(gr_norm_noPU_prompt_gen_EE, "gen noPU");
  leg_norm_prompt_EE_dt->AddEntry(gr_norm_noPU_prompt_EE, "noPU");
  leg_norm_prompt_EE_dt->SetTextSize(0.03);
  TLegend* leg_norm_prompt_EE_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_PU200_prompt_gen_EE, "gen PU200");
  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_PU200_prompt_EE, "no MTD PU200");
  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_PU200_prompt_2sigma_EE, "2sigma PU200");
//  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_PU200_prompt_3sigma_EE, "3sigma PU200");
//  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_PU200_prompt_4sigma_EE, "4sigma PU200");
  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_noPU_prompt_gen_EE, "gen noPU");
  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_noPU_prompt_EE, "noPU");
  leg_norm_prompt_EE_sigma->SetTextSize(0.03);

  // Nonprompt
    // Barrel region
      // efficiency
  TLegend* leg_eff_nonprompt_EB_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_gen_EB, "gen PU200");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_EB, "no MTD PU200");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_40_EB, "40ps PU200");
//  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_60_EB, "60ps PU200");
//  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_80_EB, "80ps PU200");
//  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_100_EB, "100ps PU200");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_noPU_nonprompt_gen_EB, "gen noPU");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_noPU_nonprompt_EB, "noPU");
  leg_eff_nonprompt_EB_dt->SetTextSize(0.03);
  TLegend* leg_eff_nonprompt_EB_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_gen_EB, "gen PU200");
  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_EB, "no MTD PU200");
  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_2sigma_EB, "2sigma PU200");
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_3sigma_EB, "3sigma PU200");
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_4sigma_EB, "4sigma PU200");
  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_noPU_nonprompt_gen_EB, "gen noPU");
  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_noPU_nonprompt_EB, "noPU");
  leg_eff_nonprompt_EB_sigma->SetTextSize(0.03);
      // normalization
  TLegend* leg_norm_nonprompt_EB_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_nonprompt_EB_dt->AddEntry(gr_norm_PU200_nonprompt_gen_EB, "gen PU200");
  leg_norm_nonprompt_EB_dt->AddEntry(gr_norm_PU200_nonprompt_EB, "no MTD PU200");
  leg_norm_nonprompt_EB_dt->AddEntry(gr_norm_PU200_nonprompt_40_EB, "40ps PU200");
//  leg_norm_nonprompt_EB_dt->AddEntry(gr_norm_PU200_nonprompt_60_EB, "60ps PU200");
//  leg_norm_nonprompt_EB_dt->AddEntry(gr_norm_PU200_nonprompt_80_EB, "80ps PU200");
//  leg_norm_nonprompt_EB_dt->AddEntry(gr_norm_PU200_nonprompt_100_EB, "100ps PU200");
  leg_norm_nonprompt_EB_dt->AddEntry(gr_norm_noPU_nonprompt_gen_EB, "gen noPU");
  leg_norm_nonprompt_EB_dt->AddEntry(gr_norm_noPU_nonprompt_EB, "noPU");
  leg_norm_nonprompt_EB_dt->SetTextSize(0.03);
  TLegend* leg_norm_nonprompt_EB_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_PU200_nonprompt_gen_EB, "gen PU200");
  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_PU200_nonprompt_EB, "no MTD PU200");
  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_PU200_nonprompt_2sigma_EB, "2sigma PU200");
//  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_PU200_nonprompt_3sigma_EB, "3sigma PU200");
//  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_PU200_nonprompt_4sigma_EB, "4sigma PU200");
  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_noPU_nonprompt_gen_EB, "gen noPU");
  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_noPU_nonprompt_EB, "noPU");
  leg_norm_nonprompt_EB_sigma->SetTextSize(0.03);
    // Endcap region
      // efficiency
  TLegend* leg_eff_nonprompt_EE_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_PU200_nonprompt_gen_EE, "gen PU200");
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_PU200_nonprompt_EE, "no MTD PU200");
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_PU200_nonprompt_40_EE, "40ps PU200");
//  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_PU200_nonprompt_60_EE, "60ps PU200");
//  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_PU200_nonprompt_80_EE, "80ps PU200");
//  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_PU200_nonprompt_100_EE, "100ps PU200");
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_noPU_nonprompt_gen_EE, "gen noPU");
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_noPU_nonprompt_EE, "noPU");
  leg_eff_nonprompt_EE_dt->SetTextSize(0.03);
  TLegend* leg_eff_nonprompt_EE_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_gen_EE, "gen PU200");
  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_EE, "no MTD PU200");
  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_2sigma_EE, "2sigma PU200");
//  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_3sigma_EE, "3sigma PU200");
//  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_4sigma_EE, "4sigma PU200");
  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_noPU_nonprompt_gen_EE, "gen noPU");
  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_noPU_nonprompt_EE, "noPU");
  leg_eff_nonprompt_EE_sigma->SetTextSize(0.03);
      // normalization
  TLegend* leg_norm_nonprompt_EE_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_nonprompt_EE_dt->AddEntry(gr_norm_PU200_nonprompt_gen_EE, "gen PU200");
  leg_norm_nonprompt_EE_dt->AddEntry(gr_norm_PU200_nonprompt_EE, "no MTD PU200");
  leg_norm_nonprompt_EE_dt->AddEntry(gr_norm_PU200_nonprompt_40_EE, "40ps PU200");
//  leg_norm_nonprompt_EE_dt->AddEntry(gr_norm_PU200_nonprompt_60_EE, "60ps PU200");
//  leg_norm_nonprompt_EE_dt->AddEntry(gr_norm_PU200_nonprompt_80_EE, "80ps PU200");
//  leg_norm_nonprompt_EE_dt->AddEntry(gr_norm_PU200_nonprompt_100_EE, "100ps PU200");
  leg_norm_nonprompt_EE_dt->AddEntry(gr_norm_noPU_nonprompt_gen_EE, "gen noPU");
  leg_norm_nonprompt_EE_dt->AddEntry(gr_norm_noPU_nonprompt_EE, "noPU");
  leg_norm_nonprompt_EE_dt->SetTextSize(0.03);
  TLegend* leg_norm_nonprompt_EE_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_PU200_nonprompt_gen_EE, "gen PU200");
  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_PU200_nonprompt_EE, "no MTD PU200");
  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_PU200_nonprompt_2sigma_EE, "2sigma PU200");
//  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_PU200_nonprompt_3sigma_EE, "3sigma PU200");
//  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_PU200_nonprompt_4sigma_EE, "4sigma PU200");
  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_noPU_nonprompt_gen_EE, "gen noPU");
  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_noPU_nonprompt_EE, "noPU");
  leg_norm_nonprompt_EE_sigma->SetTextSize(0.03);


  ////////////////
  // Draw plots //
  ////////////////
  // Prompt
    // Barrel region
      // efficiency
  TCanvas* c_eff_prompt_dt_EB = new TCanvas("c_eff_prompt_dt_EB", "c_eff_prompt_dt_EB", 1500, 1500);
  c_eff_prompt_dt_EB->cd();
  c_eff_prompt_dt_EB->SetGrid();
  c_eff_prompt_dt_EB->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EB->Draw("AL");
  gr_eff_PU200_prompt_40_EB->Draw("same");
//  gr_eff_PU200_prompt_60_EB->Draw("same");
//  gr_eff_PU200_prompt_80_EB->Draw("same");
//  gr_eff_PU200_prompt_100_EB->Draw("same");
  gr_eff_noPU_prompt_EB->Draw("same");
  gr_eff_PU200_prompt_gen_EB->Draw("same");
  gr_eff_noPU_prompt_gen_EB->Draw("same");
  leg_eff_prompt_EB_dt->Draw();
  c_eff_prompt_dt_EB->Print("plots/eff_prompt_dt_EB.pdf");

  TCanvas* c_eff_prompt_sigma_EB = new TCanvas("c_eff_prompt_sigma_EB", "c_eff_prompt_sigma_EB", 1500, 1500);
  c_eff_prompt_sigma_EB->cd();
  c_eff_prompt_sigma_EB->SetGrid();
  c_eff_prompt_sigma_EB->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EB->Draw("AL");
  gr_eff_PU200_prompt_2sigma_EB->Draw("same");
//  gr_eff_PU200_prompt_3sigma_EB->Draw("same");
//  gr_eff_PU200_prompt_4sigma_EB->Draw("same");
  gr_eff_noPU_prompt_EB->Draw("same");
  gr_eff_PU200_prompt_gen_EB->Draw("same");
  gr_eff_noPU_prompt_gen_EB->Draw("same");
  leg_eff_prompt_EB_sigma->Draw();
  c_eff_prompt_sigma_EB->Print("plots/eff_prompt_sigma_EB.pdf");

      // normalization
  TCanvas* c_norm_prompt_dt_EB = new TCanvas("c_norm_prompt_dt_EB", "c_norm_prompt_dt_EB", 1500, 1500);
  c_norm_prompt_dt_EB->cd();
  c_norm_prompt_dt_EB->SetGrid();
  c_norm_prompt_dt_EB->SetLeftMargin(0.15);
  gr_norm_noPU_prompt_EB->Draw("AL");
  gr_norm_PU200_prompt_EB->Draw("same");
  gr_norm_PU200_prompt_40_EB->Draw("same");
//  gr_norm_PU200_prompt_60_EB->Draw("same");
//  gr_norm_PU200_prompt_80_EB->Draw("same");
//  gr_norm_PU200_prompt_100_EB->Draw("same");
  gr_norm_noPU_prompt_gen_EB->Draw("same");
  gr_norm_PU200_prompt_gen_EB->Draw("same");
  leg_norm_prompt_EB_dt->Draw();
  c_norm_prompt_dt_EB->Print("plots/norm_prompt_dt_EB.pdf");

  TCanvas* c_norm_prompt_sigma_EB = new TCanvas("c_norm_prompt_sigma_EB", "c_norm_prompt_sigma_EB", 1500, 1500);
  c_norm_prompt_sigma_EB->cd();
  c_norm_prompt_sigma_EB->SetGrid();
  c_norm_prompt_sigma_EB->SetLeftMargin(0.15);
  gr_norm_noPU_prompt_EB->Draw("AL");
  gr_norm_PU200_prompt_EB->Draw("same");
  gr_norm_PU200_prompt_2sigma_EB->Draw("same");
//  gr_norm_PU200_prompt_3sigma_EB->Draw("same");
//  gr_norm_PU200_prompt_4sigma_EB->Draw("same");
  gr_norm_noPU_prompt_gen_EB->Draw("same");
  gr_norm_PU200_prompt_gen_EB->Draw("same");
  leg_norm_prompt_EB_sigma->Draw();
  c_norm_prompt_sigma_EB->Print("plots/norm_prompt_sigma_EB.pdf");

    // Endcap region
      // efficiency
  TCanvas* c_eff_prompt_dt_EE = new TCanvas("c_eff_prompt_dt_EE", "c_eff_prompt_dt_EE", 1500, 1500);
  c_eff_prompt_dt_EE->cd();
  c_eff_prompt_dt_EE->SetGrid();
  c_eff_prompt_dt_EE->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EE->Draw("AL");
  gr_eff_PU200_prompt_40_EE->Draw("same");
//  gr_eff_PU200_prompt_60_EE->Draw("same");
//  gr_eff_PU200_prompt_80_EE->Draw("same");
//  gr_eff_PU200_prompt_100_EE->Draw("same");
  gr_eff_noPU_prompt_EE->Draw("same");
  gr_eff_PU200_prompt_gen_EE->Draw("same");
  gr_eff_noPU_prompt_gen_EE->Draw("same");
  leg_eff_prompt_EE_dt->Draw();
  c_eff_prompt_dt_EE->Print("plots/eff_prompt_dt_EE.pdf");

  TCanvas* c_eff_prompt_sigma_EE = new TCanvas("c_eff_prompt_sigma_EE", "c_eff_prompt_sigma_EE", 1500, 1500);
  c_eff_prompt_sigma_EE->cd();
  c_eff_prompt_sigma_EE->SetGrid();
  c_eff_prompt_sigma_EE->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EE->Draw("AL");
  gr_eff_PU200_prompt_2sigma_EE->Draw("same");
//  gr_eff_PU200_prompt_3sigma_EE->Draw("same");
//  gr_eff_PU200_prompt_4sigma_EE->Draw("same");
  gr_eff_noPU_prompt_EE->Draw("same");
  gr_eff_PU200_prompt_gen_EE->Draw("same");
  gr_eff_noPU_prompt_gen_EE->Draw("same");
  leg_eff_prompt_EE_sigma->Draw();
  c_eff_prompt_sigma_EE->Print("plots/eff_prompt_sigma_EE.pdf");

      // normalization
  TCanvas* c_norm_prompt_dt_EE = new TCanvas("c_norm_prompt_dt_EE", "c_norm_prompt_dt_EE", 1500, 1500);
  c_norm_prompt_dt_EE->cd();
  c_norm_prompt_dt_EE->SetGrid();
  c_norm_prompt_dt_EE->SetLeftMargin(0.15);
  gr_norm_noPU_prompt_EE->Draw("AL");
  gr_norm_PU200_prompt_EE->Draw("same");
  gr_norm_PU200_prompt_40_EE->Draw("same");
//  gr_norm_PU200_prompt_60_EE->Draw("same");
//  gr_norm_PU200_prompt_80_EE->Draw("same");
//  gr_norm_PU200_prompt_100_EE->Draw("same");
  gr_norm_noPU_prompt_gen_EE->Draw("same");
  gr_norm_PU200_prompt_gen_EE->Draw("same");
  leg_norm_prompt_EE_dt->Draw();
  c_norm_prompt_dt_EE->Print("plots/norm_prompt_dt_EE.pdf");

  TCanvas* c_norm_prompt_sigma_EE = new TCanvas("c_norm_prompt_sigma_EE", "c_norm_prompt_sigma_EE", 1500, 1500);
  c_norm_prompt_sigma_EE->cd();
  c_norm_prompt_sigma_EE->SetGrid();
  c_norm_prompt_sigma_EE->SetLeftMargin(0.15);
  gr_norm_noPU_prompt_EE->Draw("AL");
  gr_norm_PU200_prompt_EE->Draw("same");
  gr_norm_PU200_prompt_2sigma_EE->Draw("same");
//  gr_norm_PU200_prompt_3sigma_EE->Draw("same");
//  gr_norm_PU200_prompt_4sigma_EE->Draw("same");
  gr_norm_noPU_prompt_gen_EE->Draw("same");
  gr_norm_PU200_prompt_gen_EE->Draw("same");
  leg_norm_prompt_EE_sigma->Draw();
  c_norm_prompt_sigma_EE->Print("plots/norm_prompt_sigma_EE.pdf");

  // Nonprompt
    // Barrel region
      // efficiency
  TCanvas* c_eff_nonprompt_dt_EB = new TCanvas("c_eff_nonprompt_dt_EB", "c_eff_nonprompt_dt_EB", 1500, 1500);
  c_eff_nonprompt_dt_EB->cd();
  c_eff_nonprompt_dt_EB->SetGrid();
  c_eff_nonprompt_dt_EB->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EB->Draw("AL");
  gr_eff_PU200_nonprompt_40_EB->Draw("same");
//  gr_eff_PU200_nonprompt_60_EB->Draw("same");
//  gr_eff_PU200_nonprompt_80_EB->Draw("same");
//  gr_eff_PU200_nonprompt_100_EB->Draw("same");
  gr_eff_noPU_nonprompt_EB->Draw("same");
  gr_eff_PU200_nonprompt_gen_EB->Draw("same");
  gr_eff_noPU_nonprompt_gen_EB->Draw("same");
  leg_eff_nonprompt_EB_dt->Draw();
  c_eff_nonprompt_dt_EB->Print("plots/eff_nonprompt_dt_EB.pdf");

  TCanvas* c_eff_nonprompt_sigma_EB = new TCanvas("c_eff_nonprompt_sigma_EB", "c_eff_nonprompt_sigma_EB", 1500, 1500);
  c_eff_nonprompt_sigma_EB->cd();
  c_eff_nonprompt_sigma_EB->SetGrid();
  c_eff_nonprompt_sigma_EB->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EB->Draw("AL");
  gr_eff_PU200_nonprompt_2sigma_EB->Draw("same");
//  gr_eff_PU200_nonprompt_3sigma_EB->Draw("same");
//  gr_eff_PU200_nonprompt_4sigma_EB->Draw("same");
  gr_eff_noPU_nonprompt_EB->Draw("same");
  gr_eff_PU200_nonprompt_gen_EB->Draw("same");
  gr_eff_noPU_nonprompt_gen_EB->Draw("same");
  leg_eff_nonprompt_EB_sigma->Draw();
  c_eff_nonprompt_sigma_EB->Print("plots/eff_nonprompt_sigma_EB.pdf");

      // normalization
  TCanvas* c_norm_nonprompt_dt_EB = new TCanvas("c_norm_nonprompt_dt_EB", "c_norm_nonprompt_dt_EB", 1500, 1500);
  c_norm_nonprompt_dt_EB->cd();
  c_norm_nonprompt_dt_EB->SetGrid();
  c_norm_nonprompt_dt_EB->SetLeftMargin(0.15);
  gr_norm_noPU_nonprompt_EB->Draw("AL");
  gr_norm_PU200_nonprompt_EB->Draw("same");
  gr_norm_PU200_nonprompt_40_EB->Draw("same");
//  gr_norm_PU200_nonprompt_60_EB->Draw("same");
//  gr_norm_PU200_nonprompt_80_EB->Draw("same");
//  gr_norm_PU200_nonprompt_100_EB->Draw("same");
  gr_norm_noPU_nonprompt_gen_EB->Draw("same");
  gr_norm_PU200_nonprompt_gen_EB->Draw("same");
  leg_norm_nonprompt_EB_dt->Draw();
  c_norm_nonprompt_dt_EB->Print("plots/norm_nonprompt_dt_EB.pdf");

  TCanvas* c_norm_nonprompt_sigma_EB = new TCanvas("c_norm_nonprompt_sigma_EB", "c_norm_nonprompt_sigma_EB", 1500, 1500);
  c_norm_nonprompt_sigma_EB->cd();
  c_norm_nonprompt_sigma_EB->SetGrid();
  c_norm_nonprompt_sigma_EB->SetLeftMargin(0.15);
  gr_norm_noPU_nonprompt_EB->Draw("AL");
  gr_norm_PU200_nonprompt_EB->Draw("same");
  gr_norm_PU200_nonprompt_2sigma_EB->Draw("same");
//  gr_norm_PU200_nonprompt_3sigma_EB->Draw("same");
//  gr_norm_PU200_nonprompt_4sigma_EB->Draw("same");
  gr_norm_noPU_nonprompt_gen_EB->Draw("same");
  gr_norm_PU200_nonprompt_gen_EB->Draw("same");
  leg_norm_nonprompt_EB_sigma->Draw();
  c_norm_nonprompt_sigma_EB->Print("plots/norm_nonprompt_sigma_EB.pdf");

    // Endcap region
      // efficiency
  TCanvas* c_eff_nonprompt_dt_EE = new TCanvas("c_eff_nonprompt_dt_EE", "c_eff_nonprompt_dt_EE", 1500, 1500);
  c_eff_nonprompt_dt_EE->cd();
  c_eff_nonprompt_dt_EE->SetGrid();
  c_eff_nonprompt_dt_EE->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EE->Draw("AL");
  gr_eff_PU200_nonprompt_40_EE->Draw("same");
//  gr_eff_PU200_nonprompt_60_EE->Draw("same");
//  gr_eff_PU200_nonprompt_80_EE->Draw("same");
//  gr_eff_PU200_nonprompt_100_EE->Draw("same");
  gr_eff_noPU_nonprompt_EE->Draw("same");
  gr_eff_PU200_nonprompt_gen_EE->Draw("same");
  gr_eff_noPU_nonprompt_gen_EE->Draw("same");
  leg_eff_nonprompt_EE_dt->Draw();
  c_eff_nonprompt_dt_EE->Print("plots/eff_nonprompt_dt_EE.pdf");

  TCanvas* c_eff_nonprompt_sigma_EE = new TCanvas("c_eff_nonprompt_sigma_EE", "c_eff_nonprompt_sigma_EE", 1500, 1500);
  c_eff_nonprompt_sigma_EE->cd();
  c_eff_nonprompt_sigma_EE->SetGrid();
  c_eff_nonprompt_sigma_EE->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EE->Draw("AL");
  gr_eff_PU200_nonprompt_2sigma_EE->Draw("same");
//  gr_eff_PU200_nonprompt_3sigma_EE->Draw("same");
//  gr_eff_PU200_nonprompt_4sigma_EE->Draw("same");
  gr_eff_noPU_nonprompt_EE->Draw("same");
  gr_eff_PU200_nonprompt_gen_EE->Draw("same");
  gr_eff_noPU_nonprompt_gen_EE->Draw("same");
  leg_eff_nonprompt_EE_sigma->Draw();
  c_eff_nonprompt_sigma_EE->Print("plots/eff_nonprompt_sigma_EE.pdf");

      // normalization
  TCanvas* c_norm_nonprompt_dt_EE = new TCanvas("c_norm_nonprompt_dt_EE", "c_norm_nonprompt_dt_EE", 1500, 1500);
  c_norm_nonprompt_dt_EE->cd();
  c_norm_nonprompt_dt_EE->SetGrid();
  c_norm_nonprompt_dt_EE->SetLeftMargin(0.15);
  gr_norm_noPU_nonprompt_EE->Draw("AL");
  gr_norm_PU200_nonprompt_EE->Draw("same");
  gr_norm_PU200_nonprompt_40_EE->Draw("same");
//  gr_norm_PU200_nonprompt_60_EE->Draw("same");
//  gr_norm_PU200_nonprompt_80_EE->Draw("same");
//  gr_norm_PU200_nonprompt_100_EE->Draw("same");
  gr_norm_noPU_nonprompt_gen_EE->Draw("same");
  gr_norm_PU200_nonprompt_gen_EE->Draw("same");
  leg_norm_nonprompt_EE_dt->Draw();
  c_norm_nonprompt_dt_EE->Print("plots/norm_nonprompt_dt_EE.pdf");

  TCanvas* c_norm_nonprompt_sigma_EE = new TCanvas("c_norm_nonprompt_sigma_EE", "c_norm_nonprompt_sigma_EE", 1500, 1500);
  c_norm_nonprompt_sigma_EE->cd();
  c_norm_nonprompt_sigma_EE->SetGrid();
  c_norm_nonprompt_sigma_EE->SetLeftMargin(0.15);
  gr_norm_noPU_nonprompt_EE->Draw("AL");
  gr_norm_PU200_nonprompt_EE->Draw("same");
  gr_norm_PU200_nonprompt_2sigma_EE->Draw("same");
//  gr_norm_PU200_nonprompt_3sigma_EE->Draw("same");
//  gr_norm_PU200_nonprompt_4sigma_EE->Draw("same");
  gr_norm_noPU_nonprompt_gen_EE->Draw("same");
  gr_norm_PU200_nonprompt_gen_EE->Draw("same");
  leg_norm_nonprompt_EE_sigma->Draw();
  c_norm_nonprompt_sigma_EE->Print("plots/norm_nonprompt_sigma_EE.pdf");
}


void draw_reliso_roc(bool zeroiso) {
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  // PU200
  f_PU200_prompt = new TFile("data/231203/harvester_PU200_prompt_160_20.root");
  f_PU200_nonprompt = new TFile("data/231203/harvester_PU200_nonprompt_90.root");
  // noPU
  f_noPU_prompt = new TFile("data/231203/harvester_noPU_prompt_19.root");
  f_noPU_nonprompt = new TFile("data/231203/harvester_noPU_nonprompt_9.root");

  // Histograms for Barrel region
  TH1D *h_PU200_prompt_EB, *h_PU200_nonprompt_EB, *h_noPU_prompt_EB, *h_noPU_nonprompt_EB;
  TH1D *h_PU200_prompt_2sigma_EB, *h_PU200_nonprompt_2sigma_EB, *h_noPU_prompt_2sigma_EB, *h_noPU_2sigma_nonprompt_EB;
  TH1D *h_PU200_prompt_3sigma_EB, *h_PU200_nonprompt_3sigma_EB, *h_noPU_prompt_3sigma_EB, *h_noPU_3sigma_nonprompt_EB;
  TH1D *h_PU200_prompt_4sigma_EB, *h_PU200_nonprompt_4sigma_EB, *h_noPU_prompt_4sigma_EB, *h_noPU_4sigma_nonprompt_EB;
  TH1D *h_PU200_prompt_40_EB, *h_PU200_nonprompt_40_EB, *h_noPU_prompt_40_EB, *h_noPU_40_nonprompt_EB;
  TH1D *h_PU200_prompt_60_EB, *h_PU200_nonprompt_60_EB, *h_noPU_prompt_60_EB, *h_noPU_60_nonprompt_EB;
  TH1D *h_PU200_prompt_80_EB, *h_PU200_nonprompt_80_EB, *h_noPU_prompt_80_EB, *h_noPU_80_nonprompt_EB;
  TH1D *h_PU200_prompt_100_EB, *h_PU200_nonprompt_100_EB, *h_noPU_prompt_100_EB, *h_noPU_100_nonprompt_EB;
  TH1D *h_PU200_prompt_gen_EB, *h_PU200_nonprompt_gen_EB, *h_noPU_prompt_gen_EB, *h_noPU_nonprompt_gen_EB;

  // Histograms for Endcap region
  TH1D *h_PU200_prompt_EE, *h_PU200_nonprompt_EE, *h_noPU_prompt_EE, *h_noPU_nonprompt_EE;
  TH1D *h_PU200_prompt_2sigma_EE, *h_PU200_nonprompt_2sigma_EE, *h_noPU_prompt_2sigma_EE, *h_noPU_2sigma_nonprompt_EE;
  TH1D *h_PU200_prompt_3sigma_EE, *h_PU200_nonprompt_3sigma_EE, *h_noPU_prompt_3sigma_EE, *h_noPU_3sigma_nonprompt_EE;
  TH1D *h_PU200_prompt_4sigma_EE, *h_PU200_nonprompt_4sigma_EE, *h_noPU_prompt_4sigma_EE, *h_noPU_4sigma_nonprompt_EE;
  TH1D *h_PU200_prompt_40_EE, *h_PU200_nonprompt_40_EE, *h_noPU_prompt_40_EE, *h_noPU_40_nonprompt_EE;
  TH1D *h_PU200_prompt_60_EE, *h_PU200_nonprompt_60_EE, *h_noPU_prompt_60_EE, *h_noPU_60_nonprompt_EE;
  TH1D *h_PU200_prompt_80_EE, *h_PU200_nonprompt_80_EE, *h_noPU_prompt_80_EE, *h_noPU_80_nonprompt_EE;
  TH1D *h_PU200_prompt_100_EE, *h_PU200_nonprompt_100_EE, *h_noPU_prompt_100_EE, *h_noPU_100_nonprompt_EE;
  TH1D *h_PU200_prompt_gen_EE, *h_PU200_nonprompt_gen_EE, *h_noPU_prompt_gen_EE, *h_noPU_nonprompt_gen_EE;


  TString dir="/DQMData/Run 1/MTD/Run summary/MuonIso/";
  // Prompt
    // Barrel region
  h_PU200_prompt_EB 	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_Sig_EB");
  h_PU200_prompt_4sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EB");
  h_PU200_prompt_3sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EB");
  h_PU200_prompt_2sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EB");
  h_PU200_prompt_40_EB     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Sig_EB");
  h_PU200_prompt_60_EB     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Sig_EB");
  h_PU200_prompt_80_EB     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Sig_EB");
  h_PU200_prompt_100_EB    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Sig_EB");
  h_PU200_prompt_gen_EB	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_gen_Sig_EB");
  h_noPU_prompt_EB    	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_Sig_EB");
  h_noPU_prompt_gen_EB 	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_gen_Sig_EB");
    // Endcap region
  h_PU200_prompt_EE 	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_Sig_EE");
  h_PU200_prompt_4sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EE");
  h_PU200_prompt_3sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EE");
  h_PU200_prompt_2sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EE");
  h_PU200_prompt_40_EE     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Sig_EE");
  h_PU200_prompt_60_EE     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Sig_EE");
  h_PU200_prompt_80_EE     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Sig_EE");
  h_PU200_prompt_100_EE    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Sig_EE");
  h_PU200_prompt_gen_EE	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_gen_Sig_EE");
  h_noPU_prompt_EE    	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_Sig_EE");
  h_noPU_prompt_gen_EE 	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_gen_Sig_EE");
  // Nonprompt
    // Barrel region
  h_PU200_nonprompt_EB 	      = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_Bkg_EB");
  h_PU200_nonprompt_4sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EB");
  h_PU200_nonprompt_3sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EB");
  h_PU200_nonprompt_2sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EB");
  h_PU200_nonprompt_40_EB     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Bkg_EB");
  h_PU200_nonprompt_60_EB     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Bkg_EB");
  h_PU200_nonprompt_80_EB     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Bkg_EB");
  h_PU200_nonprompt_100_EB    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Bkg_EB");
  h_PU200_nonprompt_gen_EB    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_gen_Bkg_EB");
  h_noPU_nonprompt_EB         = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_Bkg_EB");
  h_noPU_nonprompt_gen_EB     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_gen_Bkg_EB");
    // Endcap region
  h_PU200_nonprompt_EE 	      = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_Bkg_EE");
  h_PU200_nonprompt_4sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EE");
  h_PU200_nonprompt_3sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EE");
  h_PU200_nonprompt_2sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EE");
  h_PU200_nonprompt_40_EE     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Bkg_EE");
  h_PU200_nonprompt_60_EE     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Bkg_EE");
  h_PU200_nonprompt_80_EE     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Bkg_EE");
  h_PU200_nonprompt_100_EE    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Bkg_EE");
  h_PU200_nonprompt_gen_EE    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_gen_Bkg_EE");
  h_noPU_nonprompt_EE         = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_Bkg_EE");
  h_noPU_nonprompt_gen_EE     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_gen_Bkg_EE");

  int nbin = h_PU200_prompt_EE->GetNbinsX();

  // Define vectors to store the iso efficiency or normalization of muons
    // PU200
      // Prompt
  vector<double> prompt_eff_PU200_EB={0}, prompt_norm_PU200_EB={0}, prompt_eff_PU200_EE={0}, prompt_norm_PU200_EE={0};
  vector<double> prompt_eff_PU200_2sigma_EB={0}, prompt_norm_PU200_2sigma_EB={0}, prompt_eff_PU200_2sigma_EE={0}, prompt_norm_PU200_2sigma_EE={0};
  vector<double> prompt_eff_PU200_3sigma_EB={0}, prompt_norm_PU200_3sigma_EB={0}, prompt_eff_PU200_3sigma_EE={0}, prompt_norm_PU200_3sigma_EE={0};
  vector<double> prompt_eff_PU200_4sigma_EB={0}, prompt_norm_PU200_4sigma_EB={0}, prompt_eff_PU200_4sigma_EE={0}, prompt_norm_PU200_4sigma_EE={0};
  vector<double> prompt_eff_PU200_40_EB={0}, prompt_norm_PU200_40_EB={0}, prompt_eff_PU200_40_EE={0}, prompt_norm_PU200_40_EE={0};
  vector<double> prompt_eff_PU200_60_EB={0}, prompt_norm_PU200_60_EB={0}, prompt_eff_PU200_60_EE={0}, prompt_norm_PU200_60_EE={0};
  vector<double> prompt_eff_PU200_80_EB={0}, prompt_norm_PU200_80_EB={0}, prompt_eff_PU200_80_EE={0}, prompt_norm_PU200_80_EE={0};
  vector<double> prompt_eff_PU200_100_EB={0}, prompt_norm_PU200_100_EB={0}, prompt_eff_PU200_100_EE={0}, prompt_norm_PU200_100_EE={0};
      // Nonprompt
  vector<double> nonprompt_eff_PU200_EB={0}, nonprompt_norm_PU200_EB={0}, nonprompt_eff_PU200_EE={0}, nonprompt_norm_PU200_EE={0};
  vector<double> nonprompt_eff_PU200_2sigma_EB={0}, nonprompt_norm_PU200_2sigma_EB={0}, nonprompt_eff_PU200_2sigma_EE={0}, nonprompt_norm_PU200_2sigma_EE={0};
  vector<double> nonprompt_eff_PU200_3sigma_EB={0}, nonprompt_norm_PU200_3sigma_EB={0}, nonprompt_eff_PU200_3sigma_EE={0}, nonprompt_norm_PU200_3sigma_EE={0};
  vector<double> nonprompt_eff_PU200_4sigma_EB={0}, nonprompt_norm_PU200_4sigma_EB={0}, nonprompt_eff_PU200_4sigma_EE={0}, nonprompt_norm_PU200_4sigma_EE={0};
  vector<double> nonprompt_eff_PU200_40_EB={0}, nonprompt_norm_PU200_40_EB={0}, nonprompt_eff_PU200_40_EE={0}, nonprompt_norm_PU200_40_EE={0};
  vector<double> nonprompt_eff_PU200_60_EB={0}, nonprompt_norm_PU200_60_EB={0}, nonprompt_eff_PU200_60_EE={0}, nonprompt_norm_PU200_60_EE={0};
  vector<double> nonprompt_eff_PU200_80_EB={0}, nonprompt_norm_PU200_80_EB={0}, nonprompt_eff_PU200_80_EE={0}, nonprompt_norm_PU200_80_EE={0};
  vector<double> nonprompt_eff_PU200_100_EB={0}, nonprompt_norm_PU200_100_EB={0}, nonprompt_eff_PU200_100_EE={0}, nonprompt_norm_PU200_100_EE={0};
      // Gen case
  vector<double> prompt_eff_PU200_gen_EB={0}, prompt_norm_PU200_gen_EB={0}, prompt_eff_PU200_gen_EE={0}, prompt_norm_PU200_gen_EE={0};
  vector<double> nonprompt_eff_PU200_gen_EB={0}, nonprompt_norm_PU200_gen_EB={0}, nonprompt_eff_PU200_gen_EE={0}, nonprompt_norm_PU200_gen_EE={0};
    // noPU
      // Prompt
  vector<double> prompt_eff_noPU_EB={0}, prompt_norm_noPU_EB={0}, prompt_eff_noPU_EE={0}, prompt_norm_noPU_EE={0};
  vector<double> prompt_eff_noPU_2sigma_EB={0}, prompt_norm_noPU_2sigma_EB={0}, prompt_eff_noPU_2sigma_EE={0}, prompt_norm_noPU_2sigma_EE={0};
  vector<double> prompt_eff_noPU_3sigma_EB={0}, prompt_norm_noPU_3sigma_EB={0}, prompt_eff_noPU_3sigma_EE={0}, prompt_norm_noPU_3sigma_EE={0};
  vector<double> prompt_eff_noPU_4sigma_EB={0}, prompt_norm_noPU_4sigma_EB={0}, prompt_eff_noPU_4sigma_EE={0}, prompt_norm_noPU_4sigma_EE={0};
  vector<double> prompt_eff_noPU_40_EB={0}, prompt_norm_noPU_40_EB={0}, prompt_eff_noPU_EE_40_EE={0}, prompt_norm_noPU_EE_40_EE={0};
  vector<double> prompt_eff_noPU_60_EB={0}, prompt_norm_noPU_60_EB={0}, prompt_eff_noPU_EE_60_EE={0}, prompt_norm_noPU_EE_60_EE={0};
  vector<double> prompt_eff_noPU_80_EB={0}, prompt_norm_noPU_80_EB={0}, prompt_eff_noPU_EE_80_EE={0}, prompt_norm_noPU_EE_80_EE={0};
  vector<double> prompt_eff_noPU_100_EB={0}, prompt_norm_noPU_100_EB={0}, prompt_eff_noPU_EE_100_EE={0}, prompt_norm_noPU_EE_100_EE={0};
      // Nonprompt
  vector<double> nonprompt_eff_noPU_EB={0}, nonprompt_norm_noPU_EB={0}, nonprompt_eff_noPU_EE={0}, nonprompt_norm_noPU_EE={0};
  vector<double> nonprompt_eff_noPU_2sigma_EB={0}, nonprompt_norm_noPU_2sigma_EB={0}, nonprompt_eff_noPU_2sigma_EE={0}, nonprompt_norm_noPU_2sigma_EE={0};
  vector<double> nonprompt_eff_noPU_3sigma_EB={0}, nonprompt_norm_noPU_3sigma_EB={0}, nonprompt_eff_noPU_3sigma_EE={0}, nonprompt_norm_noPU_3sigma_EE={0};
  vector<double> nonprompt_eff_noPU_4sigma_EB={0}, nonprompt_norm_noPU_4sigma_EB={0}, nonprompt_eff_noPU_4sigma_EE={0}, nonprompt_norm_noPU_4sigma_EE={0};
  vector<double> nonprompt_eff_noPU_40_EB={0}, nonprompt_norm_noPU_40_EB={0}, nonprompt_eff_noPU_40_EE={0}, nonprompt_norm_noPU_40_EE={0};
  vector<double> nonprompt_eff_noPU_60_EB={0}, nonprompt_norm_noPU_60_EB={0}, nonprompt_eff_noPU_60_EE={0}, nonprompt_norm_noPU_60_EE={0};
  vector<double> nonprompt_eff_noPU_80_EB={0}, nonprompt_norm_noPU_80_EB={0}, nonprompt_eff_noPU_80_EE={0}, nonprompt_norm_noPU_80_EE={0};
  vector<double> nonprompt_eff_noPU_100_EB={0}, nonprompt_norm_noPU_100_EB={0}, nonprompt_eff_noPU_100_EE={0}, nonprompt_norm_noPU_100_EE={0};
      // Gen case
  vector<double> prompt_eff_noPU_gen_EB={0}, prompt_norm_noPU_gen_EB={0}, prompt_eff_noPU_gen_EE={0}, prompt_norm_noPU_gen_EE={0};
  vector<double> nonprompt_eff_noPU_gen_EB={0}, nonprompt_norm_noPU_gen_EB={0}, nonprompt_eff_noPU_gen_EE={0}, nonprompt_norm_noPU_gen_EE={0};


  //////////////////////////
  // Calculate efficiency //
  //////////////////////////
  for(unsigned int i=0; i<h_PU200_prompt_EB->GetNbinsX(); i++) {
    if(zeroiso==true) {        // Include the cases where the muon is already isolated
    // Prompt
      // Barrel region
        // efficiency
      prompt_eff_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(1,i+1)/h_PU200_prompt_EB->Integral(1,nbin));
      prompt_eff_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(1,i+1)/h_noPU_prompt_EB->Integral(1,nbin));
      prompt_eff_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(1,i+1)/h_PU200_prompt_40_EB->Integral(1,nbin));
      prompt_eff_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(1,i+1)/h_PU200_prompt_60_EB->Integral(1,nbin));
      prompt_eff_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(1,i+1)/h_PU200_prompt_80_EB->Integral(1,nbin));
      prompt_eff_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(1,i+1)/h_PU200_prompt_100_EB->Integral(1,nbin));
      prompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(1,i+1)/h_PU200_prompt_2sigma_EB->Integral(1,nbin));
      prompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(1,i+1)/h_PU200_prompt_3sigma_EB->Integral(1,nbin));
      prompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(1,i+1)/h_PU200_prompt_4sigma_EB->Integral(1,nbin));
      prompt_eff_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(1,i+1)/h_PU200_prompt_gen_EB->Integral(1,nbin));
      prompt_eff_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(1,i+1)/h_noPU_prompt_gen_EB->Integral(1,nbin));
        // normalization
      prompt_norm_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(1,i+1));
      prompt_norm_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(1,i+1));
      prompt_norm_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(1,i+1));
      prompt_norm_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(1,i+1));
      prompt_norm_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(1,i+1));
      prompt_norm_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(1,i+1));
      prompt_norm_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(1,i+1));
      prompt_norm_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(1,i+1));
      prompt_norm_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(1,i+1));
      prompt_norm_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(1,i+1));
      prompt_norm_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(1,i+1));
      // Endcap region
        // efficiency
      prompt_eff_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(1,i+1)/h_PU200_prompt_EE->Integral(1,nbin));
      prompt_eff_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(1,i+1)/h_noPU_prompt_EE->Integral(1,nbin));
      prompt_eff_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(1,i+1)/h_PU200_prompt_40_EE->Integral(1,nbin));
      prompt_eff_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(1,i+1)/h_PU200_prompt_60_EE->Integral(1,nbin));
      prompt_eff_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(1,i+1)/h_PU200_prompt_80_EE->Integral(1,nbin));
      prompt_eff_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(1,i+1)/h_PU200_prompt_100_EE->Integral(1,nbin));
      prompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(1,i+1)/h_PU200_prompt_2sigma_EE->Integral(1,nbin));
      prompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(1,i+1)/h_PU200_prompt_3sigma_EE->Integral(1,nbin));
      prompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(1,i+1)/h_PU200_prompt_4sigma_EE->Integral(1,nbin));
      prompt_eff_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(1,i+1)/h_PU200_prompt_gen_EE->Integral(1,nbin));
      prompt_eff_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(1,i+1)/h_noPU_prompt_gen_EE->Integral(1,nbin));
        // normalization
      prompt_norm_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(1,i+1));
      prompt_norm_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(1,i+1));
      prompt_norm_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(1,i+1));
      prompt_norm_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(1,i+1));
      prompt_norm_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(1,i+1));
      prompt_norm_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(1,i+1));
      prompt_norm_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(1,i+1));
      prompt_norm_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(1,i+1));
      prompt_norm_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(1,i+1));
      prompt_norm_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(1,i+1));
      prompt_norm_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(1,i+1));
    // Nonprompt
      // Barrel region
        // efficiency
      nonprompt_eff_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(1,i+1)/h_PU200_nonprompt_EB->Integral(1,nbin));
      nonprompt_eff_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(1,i+1)/h_noPU_nonprompt_EB->Integral(1,nbin));
      nonprompt_eff_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(1,i+1)/h_PU200_nonprompt_40_EB->Integral(1,nbin));
      nonprompt_eff_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(1,i+1)/h_PU200_nonprompt_60_EB->Integral(1,nbin));
      nonprompt_eff_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(1,i+1)/h_PU200_nonprompt_80_EB->Integral(1,nbin));
      nonprompt_eff_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(1,i+1)/h_PU200_nonprompt_100_EB->Integral(1,nbin));
      nonprompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EB->Integral(1,nbin));
      nonprompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EB->Integral(1,nbin));
      nonprompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EB->Integral(1,nbin));
      nonprompt_eff_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(1,i+1)/h_PU200_nonprompt_gen_EB->Integral(1,nbin));
      nonprompt_eff_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(1,i+1)/h_noPU_nonprompt_gen_EB->Integral(1,nbin));
        // normalization
      nonprompt_norm_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(1,i+1));
      nonprompt_norm_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(1,i+1));
      nonprompt_norm_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(1,i+1));
      nonprompt_norm_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(1,i+1));
      nonprompt_norm_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(1,i+1));
      nonprompt_norm_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(1,i+1));
      nonprompt_norm_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(1,i+1));
      nonprompt_norm_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(1,i+1));
      nonprompt_norm_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(1,i+1));
      nonprompt_norm_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(1,i+1));
      nonprompt_norm_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(1,i+1));
      // Endcap region
        // efficiency
      nonprompt_eff_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(1,i+1)/h_PU200_nonprompt_EE->Integral(1,nbin));
      nonprompt_eff_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(1,i+1)/h_noPU_nonprompt_EE->Integral(1,nbin));
      nonprompt_eff_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(1,i+1)/h_PU200_nonprompt_40_EE->Integral(1,nbin));
      nonprompt_eff_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(1,i+1)/h_PU200_nonprompt_60_EE->Integral(1,nbin));
      nonprompt_eff_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(1,i+1)/h_PU200_nonprompt_80_EE->Integral(1,nbin));
      nonprompt_eff_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(1,i+1)/h_PU200_nonprompt_100_EE->Integral(1,nbin));
      nonprompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EE->Integral(1,nbin));
      nonprompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EE->Integral(1,nbin));
      nonprompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EE->Integral(1,nbin));
      nonprompt_eff_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(1,i+1)/h_PU200_nonprompt_gen_EE->Integral(1,nbin));
      nonprompt_eff_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(1,i+1)/h_noPU_nonprompt_gen_EE->Integral(1,nbin));
        // normalization
      nonprompt_norm_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(1,i+1));
      nonprompt_norm_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(1,i+1));
      nonprompt_norm_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(1,i+1));
      nonprompt_norm_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(1,i+1));
      nonprompt_norm_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(1,i+1));
      nonprompt_norm_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(1,i+1));
      nonprompt_norm_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(1,i+1));
      nonprompt_norm_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(1,i+1));
      nonprompt_norm_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(1,i+1));
      nonprompt_norm_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(1,i+1));
      nonprompt_norm_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(1,i+1));
    }
    else {             // Not include the cases where the muon is already isolated
    // Prompt
      // Barrel region
        // efficiency
      prompt_eff_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(2,2+i)/h_PU200_prompt_EB->Integral(2,nbin));
      prompt_eff_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(2,2+i)/h_noPU_prompt_EB->Integral(2,nbin));
      prompt_eff_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(2,2+i)/h_PU200_prompt_40_EB->Integral(2,nbin));
      prompt_eff_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(2,2+i)/h_PU200_prompt_60_EB->Integral(2,nbin));
      prompt_eff_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(2,2+i)/h_PU200_prompt_80_EB->Integral(2,nbin));
      prompt_eff_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(2,2+i)/h_PU200_prompt_100_EB->Integral(2,nbin));
      prompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(2,2+i)/h_PU200_prompt_2sigma_EB->Integral(2,nbin));
      prompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(2,2+i)/h_PU200_prompt_3sigma_EB->Integral(2,nbin));
      prompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(2,2+i)/h_PU200_prompt_4sigma_EB->Integral(2,nbin));
      prompt_eff_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(2,2+i)/h_PU200_prompt_gen_EB->Integral(2,nbin));
      prompt_eff_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(2,2+i)/h_noPU_prompt_gen_EB->Integral(2,nbin));
        // normalization
      prompt_norm_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(2,2+i));
      prompt_norm_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(2,2+i));
      prompt_norm_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(2,2+i));
      prompt_norm_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(2,2+i));
      prompt_norm_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(2,2+i));
      prompt_norm_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(2,2+i));
      prompt_norm_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(2,2+i));
      prompt_norm_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(2,2+i));
      prompt_norm_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(2,2+i));
      prompt_norm_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(2,2+i));
      prompt_norm_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(2,2+i));
      // Endcap region
        // efficiency
      prompt_eff_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(2,2+i)/h_PU200_prompt_EE->Integral(2,nbin));
      prompt_eff_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(2,2+i)/h_noPU_prompt_EE->Integral(2,nbin));
      prompt_eff_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(2,2+i)/h_PU200_prompt_40_EE->Integral(2,nbin));
      prompt_eff_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(2,2+i)/h_PU200_prompt_60_EE->Integral(2,nbin));
      prompt_eff_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(2,2+i)/h_PU200_prompt_80_EE->Integral(2,nbin));
      prompt_eff_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(2,2+i)/h_PU200_prompt_100_EE->Integral(2,nbin));
      prompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(2,2+i)/h_PU200_prompt_2sigma_EE->Integral(2,nbin));
      prompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(2,2+i)/h_PU200_prompt_3sigma_EE->Integral(2,nbin));
      prompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(2,2+i)/h_PU200_prompt_4sigma_EE->Integral(2,nbin));
      prompt_eff_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(2,2+i)/h_PU200_prompt_gen_EE->Integral(2,nbin));
      prompt_eff_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(2,2+i)/h_noPU_prompt_gen_EE->Integral(2,nbin));
        // normalization
      prompt_norm_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(2,2+i));
      prompt_norm_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(2,2+i));
      prompt_norm_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(2,2+i));
      prompt_norm_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(2,2+i));
      prompt_norm_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(2,2+i));
      prompt_norm_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(2,2+i));
      prompt_norm_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(2,2+i));
      prompt_norm_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(2,2+i));
      prompt_norm_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(2,2+i));
      prompt_norm_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(2,2+i));
      prompt_norm_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(2,2+i));
    // Nonprompt
      // Barrel region
        // efficiency
      nonprompt_eff_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(2,2+i)/h_PU200_nonprompt_EB->Integral(2,nbin));
      nonprompt_eff_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(2,2+i)/h_noPU_nonprompt_EB->Integral(2,nbin));
      nonprompt_eff_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(2,2+i)/h_PU200_nonprompt_40_EB->Integral(2,nbin));
      nonprompt_eff_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(2,2+i)/h_PU200_nonprompt_60_EB->Integral(2,nbin));
      nonprompt_eff_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(2,2+i)/h_PU200_nonprompt_80_EB->Integral(2,nbin));
      nonprompt_eff_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(2,2+i)/h_PU200_nonprompt_100_EB->Integral(2,nbin));
      nonprompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(2,2+i)/h_PU200_nonprompt_2sigma_EB->Integral(2,nbin));
      nonprompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(2,2+i)/h_PU200_nonprompt_3sigma_EB->Integral(2,nbin));
      nonprompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(2,2+i)/h_PU200_nonprompt_4sigma_EB->Integral(2,nbin));
      nonprompt_eff_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(2,2+i)/h_PU200_nonprompt_gen_EB->Integral(2,nbin));
      nonprompt_eff_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(2,2+i)/h_noPU_nonprompt_gen_EB->Integral(2,nbin));
        // normalization
      nonprompt_norm_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(2,2+i));
      nonprompt_norm_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(2,2+i));
      nonprompt_norm_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(2,2+i));
      nonprompt_norm_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(2,2+i));
      nonprompt_norm_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(2,2+i));
      nonprompt_norm_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(2,2+i));
      nonprompt_norm_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(2,2+i));
      nonprompt_norm_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(2,2+i));
      nonprompt_norm_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(2,2+i));
      nonprompt_norm_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(2,2+i));
      nonprompt_norm_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(2,2+i));
      // Endcap region
        // efficiency
      nonprompt_eff_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(2,2+i)/h_PU200_nonprompt_EE->Integral(2,nbin));
      nonprompt_eff_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(2,2+i)/h_noPU_nonprompt_EE->Integral(2,nbin));
      nonprompt_eff_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(2,2+i)/h_PU200_nonprompt_40_EE->Integral(2,nbin));
      nonprompt_eff_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(2,2+i)/h_PU200_nonprompt_60_EE->Integral(2,nbin));
      nonprompt_eff_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(2,2+i)/h_PU200_nonprompt_80_EE->Integral(2,nbin));
      nonprompt_eff_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(2,2+i)/h_PU200_nonprompt_100_EE->Integral(2,nbin));
      nonprompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(2,2+i)/h_PU200_nonprompt_2sigma_EE->Integral(2,nbin));
      nonprompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(2,2+i)/h_PU200_nonprompt_3sigma_EE->Integral(2,nbin));
      nonprompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(2,2+i)/h_PU200_nonprompt_4sigma_EE->Integral(2,nbin));
      nonprompt_eff_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(2,2+i)/h_PU200_nonprompt_gen_EE->Integral(2,nbin));
      nonprompt_eff_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(2,2+i)/h_noPU_nonprompt_gen_EE->Integral(2,nbin));
        // normalization
      nonprompt_norm_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(2,2+i));
      nonprompt_norm_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(2,2+i));
      nonprompt_norm_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(2,2+i));
      nonprompt_norm_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(2,2+i));
      nonprompt_norm_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(2,2+i));
      nonprompt_norm_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(2,2+i));
      nonprompt_norm_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(2,2+i));
      nonprompt_norm_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(2,2+i));
      nonprompt_norm_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(2,2+i));
      nonprompt_norm_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(2,2+i));
      nonprompt_norm_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(2,2+i));
    }
  }

  // Define TGraph
  // Prompt
    // Barrel region
  TGraph* gr_PU200_EB = new TGraph(h_PU200_prompt_EB->GetNbinsX(), &prompt_eff_PU200_EB[0], &nonprompt_eff_PU200_EB[0]);
  TGraph* gr_PU200_40_EB = new TGraph(h_PU200_prompt_40_EB->GetNbinsX(), &prompt_eff_PU200_40_EB[0], &nonprompt_eff_PU200_40_EB[0]);
  TGraph* gr_PU200_60_EB = new TGraph(h_PU200_prompt_60_EB->GetNbinsX(), &prompt_eff_PU200_60_EB[0], &nonprompt_eff_PU200_60_EB[0]);
  TGraph* gr_PU200_80_EB = new TGraph(h_PU200_prompt_80_EB->GetNbinsX(), &prompt_eff_PU200_80_EB[0], &nonprompt_eff_PU200_80_EB[0]);
  TGraph* gr_PU200_100_EB = new TGraph(h_PU200_prompt_100_EB->GetNbinsX(), &prompt_eff_PU200_100_EB[0], &nonprompt_eff_PU200_100_EB[0]);
  TGraph* gr_PU200_2sigma_EB = new TGraph(h_PU200_prompt_2sigma_EB->GetNbinsX(), &prompt_eff_PU200_2sigma_EB[0], &nonprompt_eff_PU200_2sigma_EB[0]);
  TGraph* gr_PU200_3sigma_EB = new TGraph(h_PU200_prompt_3sigma_EB->GetNbinsX(), &prompt_eff_PU200_3sigma_EB[0], &nonprompt_eff_PU200_3sigma_EB[0]);
  TGraph* gr_PU200_4sigma_EB = new TGraph(h_PU200_prompt_4sigma_EB->GetNbinsX(), &prompt_eff_PU200_4sigma_EB[0], &nonprompt_eff_PU200_4sigma_EB[0]);
  TGraph* gr_noPU_EB = new TGraph(h_noPU_prompt_EB->GetNbinsX(), &prompt_eff_noPU_EB[0], &nonprompt_eff_noPU_EB[0]);
  TGraph* gr_PU200_gen_EB = new TGraph(h_PU200_prompt_gen_EB->GetNbinsX(), &prompt_eff_PU200_gen_EB[0], &nonprompt_eff_PU200_gen_EB[0]);
  TGraph* gr_noPU_gen_EB = new TGraph(h_noPU_prompt_gen_EB->GetNbinsX(), &prompt_eff_noPU_gen_EB[0], &nonprompt_eff_noPU_gen_EB[0]);

    // Endcap region
  TGraph* gr_PU200_EE = new TGraph(h_PU200_prompt_EE->GetNbinsX(), &prompt_eff_PU200_EE[0], &nonprompt_eff_PU200_EE[0]);
  TGraph* gr_PU200_40_EE = new TGraph(h_PU200_prompt_40_EE->GetNbinsX(), &prompt_eff_PU200_40_EE[0], &nonprompt_eff_PU200_40_EE[0]);
  TGraph* gr_PU200_60_EE = new TGraph(h_PU200_prompt_60_EE->GetNbinsX(), &prompt_eff_PU200_60_EE[0], &nonprompt_eff_PU200_60_EE[0]);
  TGraph* gr_PU200_80_EE = new TGraph(h_PU200_prompt_80_EE->GetNbinsX(), &prompt_eff_PU200_80_EE[0], &nonprompt_eff_PU200_80_EE[0]);
  TGraph* gr_PU200_100_EE = new TGraph(h_PU200_prompt_100_EE->GetNbinsX(), &prompt_eff_PU200_100_EE[0], &nonprompt_eff_PU200_100_EE[0]);
  TGraph* gr_PU200_2sigma_EE = new TGraph(h_PU200_prompt_2sigma_EE->GetNbinsX(), &prompt_eff_PU200_2sigma_EE[0], &nonprompt_eff_PU200_2sigma_EE[0]);
  TGraph* gr_PU200_3sigma_EE = new TGraph(h_PU200_prompt_3sigma_EE->GetNbinsX(), &prompt_eff_PU200_3sigma_EE[0], &nonprompt_eff_PU200_3sigma_EE[0]);
  TGraph* gr_PU200_4sigma_EE = new TGraph(h_PU200_prompt_4sigma_EE->GetNbinsX(), &prompt_eff_PU200_4sigma_EE[0], &nonprompt_eff_PU200_4sigma_EE[0]);
  TGraph* gr_noPU_EE = new TGraph(h_noPU_prompt_EE->GetNbinsX(), &prompt_eff_noPU_EE[0], &nonprompt_eff_noPU_EE[0]);
  TGraph* gr_PU200_gen_EE = new TGraph(h_PU200_prompt_gen_EE->GetNbinsX(), &prompt_eff_PU200_gen_EE[0], &nonprompt_eff_PU200_gen_EE[0]);
  TGraph* gr_noPU_gen_EE = new TGraph(h_noPU_prompt_gen_EE->GetNbinsX(), &prompt_eff_noPU_gen_EE[0], &nonprompt_eff_noPU_gen_EE[0]);

  // Remove the value at dump point (0,0)
  gr_PU200_EB->RemovePoint(0); gr_PU200_40_EB->RemovePoint(0); gr_PU200_60_EB->RemovePoint(0); gr_PU200_80_EB->RemovePoint(0); gr_PU200_100_EB->RemovePoint(0); gr_PU200_2sigma_EB->RemovePoint(0); gr_PU200_3sigma_EB->RemovePoint(0); gr_PU200_4sigma_EB->RemovePoint(0); gr_noPU_EB->RemovePoint(0); gr_PU200_gen_EB->RemovePoint(0); gr_noPU_gen_EB->RemovePoint(0);
  gr_PU200_EE->RemovePoint(0); gr_PU200_40_EE->RemovePoint(0); gr_PU200_60_EE->RemovePoint(0); gr_PU200_80_EE->RemovePoint(0); gr_PU200_100_EE->RemovePoint(0); gr_PU200_2sigma_EE->RemovePoint(0); gr_PU200_3sigma_EE->RemovePoint(0); gr_PU200_4sigma_EE->RemovePoint(0); gr_noPU_EE->RemovePoint(0);gr_PU200_gen_EE->RemovePoint(0); gr_noPU_gen_EE->RemovePoint(0);


  ///////////////
  // Cosmetics //
  ///////////////
    // Barrel region
  gr_PU200_EB->SetTitle("ROC curve in barrel region"); gr_noPU_EB->SetTitle("ROC  curve in barrel region");
  gr_PU200_EB->GetXaxis()->SetTitle("Prompt efficiency"); gr_noPU_EB->GetXaxis()->SetTitle("Prompt efficiency");
  gr_PU200_EB->GetYaxis()->SetTitle("Non-prompt efficiency"); gr_noPU_EB->GetYaxis()->SetTitle("Non-prompt efficiency");
  gr_PU200_EB->SetLineWidth(2); gr_PU200_40_EB->SetLineWidth(2); gr_PU200_60_EB->SetLineWidth(2); gr_PU200_80_EB->SetLineWidth(2); gr_PU200_100_EB->SetLineWidth(2);
  gr_PU200_2sigma_EB->SetLineWidth(2); gr_PU200_3sigma_EB->SetLineWidth(2); gr_PU200_4sigma_EB->SetLineWidth(2);
  gr_noPU_EB->SetLineWidth(2); gr_PU200_gen_EB->SetLineWidth(2); gr_noPU_gen_EB->SetLineWidth(2);
  gr_PU200_EB->SetLineColor(kBlack); gr_PU200_40_EB->SetLineColor(kRed); gr_PU200_60_EB->SetLineColor(kGreen); gr_PU200_80_EB->SetLineColor(kBlue); gr_PU200_100_EB->SetLineColor(kMagenta);
  gr_PU200_2sigma_EB->SetLineColor(kRed); gr_PU200_3sigma_EB->SetLineColor(kGreen); gr_PU200_4sigma_EB->SetLineColor(kBlue);
  gr_noPU_EB->SetLineColor(kGray); gr_PU200_gen_EB->SetLineColor(kBlack); gr_noPU_gen_EB->SetLineColor(kGray);
  gr_PU200_gen_EB->SetLineStyle(7); gr_noPU_gen_EB->SetLineStyle(7);

    // Endcap region
  gr_PU200_EE->SetTitle("ROC curve in endcap region"); gr_noPU_EE->SetTitle("ROC curve in endcap region");
  gr_PU200_EE->GetXaxis()->SetTitle("Prompt efficiency"); gr_noPU_EE->GetXaxis()->SetTitle("Prompt efficiency");
  gr_PU200_EE->GetYaxis()->SetTitle("Non-prompt efficiency"); gr_noPU_EE->GetYaxis()->SetTitle("Non-prompt efficiency");
  gr_PU200_EE->SetLineWidth(2); gr_PU200_40_EE->SetLineWidth(2); gr_PU200_60_EE->SetLineWidth(2); gr_PU200_80_EE->SetLineWidth(2); gr_PU200_100_EE->SetLineWidth(2);
  gr_PU200_2sigma_EE->SetLineWidth(2); gr_PU200_3sigma_EE->SetLineWidth(2); gr_PU200_4sigma_EE->SetLineWidth(2);
  gr_noPU_EE->SetLineWidth(2); gr_PU200_gen_EE->SetLineWidth(2); gr_noPU_gen_EE->SetLineWidth(2);
  gr_PU200_EE->SetLineColor(kBlack); gr_PU200_40_EE->SetLineColor(kRed); gr_PU200_60_EE->SetLineColor(kGreen); gr_PU200_80_EE->SetLineColor(kBlue); gr_PU200_100_EE->SetLineColor(kMagenta);
  gr_PU200_2sigma_EE->SetLineColor(kRed); gr_PU200_3sigma_EE->SetLineColor(kGreen); gr_PU200_4sigma_EE->SetLineColor(kBlue);
  gr_noPU_EE->SetLineColor(kGray); gr_PU200_gen_EE->SetLineColor(kBlack); gr_noPU_gen_EE->SetLineColor(kGray);
  gr_PU200_gen_EE->SetLineStyle(7); gr_noPU_gen_EE->SetLineStyle(7);


  /////////////
  // Legends //
  /////////////
  // Barrel region
  TLegend* leg_EB_dt = new TLegend(0.15, 0.65, 0.48, 0.85);
  leg_EB_dt->AddEntry(gr_PU200_gen_EB, "gen PU200");
  leg_EB_dt->AddEntry(gr_PU200_EB, "no MTD PU200");
  leg_EB_dt->AddEntry(gr_PU200_40_EB, "40ps PU200");
//  leg_EB_dt->AddEntry(gr_PU200_60_EB, "60ps PU200");
//  leg_EB_dt->AddEntry(gr_PU200_80_EB, "80ps PU200");
//  leg_EB_dt->AddEntry(gr_PU200_100_EB, "100ps PU200");
  leg_EB_dt->AddEntry(gr_noPU_gen_EB, "gen noPU");
  leg_EB_dt->AddEntry(gr_noPU_EB, "noPU");
  leg_EB_dt->SetTextSize(0.03);

  TLegend* leg_EB_sigma = new TLegend(0.15, 0.65, 0.48, 0.85);
  leg_EB_sigma->AddEntry(gr_PU200_gen_EB, "gen PU200");
  leg_EB_sigma->AddEntry(gr_PU200_EB, "no MTD PU200");
  leg_EB_sigma->AddEntry(gr_PU200_2sigma_EB, "2sigma PU200");
//  leg_EB_sigma->AddEntry(gr_PU200_3sigma_EB, "3sigma PU200");
//  leg_EB_sigma->AddEntry(gr_PU200_4sigma_EB, "4sigma PU200");
  leg_EB_sigma->AddEntry(gr_noPU_gen_EB, "gen noPU");
  leg_EB_sigma->AddEntry(gr_noPU_EB, "noPU");
  leg_EB_sigma->SetTextSize(0.03);

  // Endcap region
  TLegend* leg_EE_dt = new TLegend(0.15, 0.65, 0.48, 0.85);
  leg_EE_dt->AddEntry(gr_PU200_gen_EE, "gen PU200");
  leg_EE_dt->AddEntry(gr_PU200_EE, "no MTD PU200");
  leg_EE_dt->AddEntry(gr_PU200_40_EE, "40ps PU200");
//  leg_EE_dt->AddEntry(gr_PU200_60_EE, "60ps PU200");
//  leg_EE_dt->AddEntry(gr_PU200_80_EE, "80ps PU200");
//  leg_EE_dt->AddEntry(gr_PU200_100_EE, "100ps PU200");
  leg_EE_dt->AddEntry(gr_noPU_gen_EE, "gen noPU");
  leg_EE_dt->AddEntry(gr_noPU_EE, "noPU");
  leg_EE_dt->SetTextSize(0.03);

  TLegend* leg_EE_sigma = new TLegend(0.15, 0.65, 0.48, 0.85);
  leg_EE_sigma->AddEntry(gr_PU200_gen_EE, "gen PU200");
  leg_EE_sigma->AddEntry(gr_PU200_EE, "no MTD PU200");
  leg_EE_sigma->AddEntry(gr_PU200_2sigma_EE, "2sigma PU200");
//  leg_EE_sigma->AddEntry(gr_PU200_3sigma_EE, "3sigma PU200");
//  leg_EE_sigma->AddEntry(gr_PU200_4sigma_EE, "4sigma PU200");
  leg_EE_sigma->AddEntry(gr_noPU_gen_EE, "gen noPU");
  leg_EE_sigma->AddEntry(gr_noPU_EE, "noPU");
  leg_EE_sigma->SetTextSize(0.03);


  ////////////////
  // Draw plots //
  ////////////////
  // Barrel region
  TCanvas* c_dt_EB = new TCanvas("c_dt_EB", "c_dt_EB", 1500, 1500);
  c_dt_EB->cd();
  c_dt_EB->SetGrid();
  gr_PU200_EB->GetXaxis()->SetRangeUser(0.,1.1);
  gr_PU200_EB->Draw("AL");
  gr_PU200_40_EB->Draw("same");
//  gr_PU200_60_EB->Draw("same");
//  gr_PU200_80_EB->Draw("same");
//  gr_PU200_100_EB->Draw("same");
  gr_noPU_EB->Draw("same");
  gr_PU200_gen_EB->Draw("same");
  gr_noPU_gen_EB->Draw("same");
  leg_EB_dt->Draw();
  c_dt_EB->Print("plots/roc_dt_EB_reliso.pdf");

  TCanvas* c_sigma_EB = new TCanvas("c_sigma_EB", "c_sigma_EB", 1500, 1500);
  c_sigma_EB->cd();
  c_sigma_EB->SetGrid();
  gr_PU200_EB->GetXaxis()->SetRangeUser(0.,1.1);
  gr_PU200_EB->Draw("AL");
  gr_PU200_2sigma_EB->Draw("same");
//  gr_PU200_3sigma_EB->Draw("same");
//  gr_PU200_4sigma_EB->Draw("same");
  gr_noPU_EB->Draw("same");
  gr_PU200_gen_EB->Draw("same");
  gr_noPU_gen_EB->Draw("same");
  leg_EB_sigma->Draw();
  c_sigma_EB->Print("plots/roc_sigma_EB_reliso.pdf");

  // Endcap region
  TCanvas* c_dt_EE = new TCanvas("c_dt_EE", "c_dt_EE", 1500, 1500);
  c_dt_EE->cd();
  c_dt_EE->SetGrid();
  gr_PU200_EE->GetXaxis()->SetRangeUser(0.,1.1);
  gr_PU200_EE->Draw("AL");
  gr_PU200_40_EE->Draw("same");
//  gr_PU200_60_EE->Draw("same");
//  gr_PU200_80_EE->Draw("same");
//  gr_PU200_100_EE->Draw("same");
  gr_noPU_EE->Draw("same");
  gr_PU200_gen_EE->Draw("same");
  gr_noPU_gen_EE->Draw("same");
  leg_EE_dt->Draw();
  c_dt_EE->Print("plots/roc_dt_EE_reliso.pdf");

  TCanvas* c_sigma_EE = new TCanvas("c_sigma_EE", "c_sigma_EE", 1500, 1500);
  c_sigma_EE->cd();
  c_sigma_EE->SetGrid();
  gr_PU200_EE->GetXaxis()->SetRangeUser(0.,1.1);
  gr_PU200_EE->Draw("AL");
  gr_PU200_2sigma_EE->Draw("same");
//  gr_PU200_3sigma_EE->Draw("same");
//  gr_PU200_4sigma_EE->Draw("same");
  gr_noPU_EE->Draw("same");
  gr_PU200_gen_EE->Draw("same");
  gr_noPU_gen_EE->Draw("same");
  leg_EE_sigma->Draw();
  c_sigma_EE->Print("plots/roc_sigma_EE_reliso.pdf");
}


void draw_pt_roc(bool zeroiso) {
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  // PU200
  f_PU200_prompt = new TFile("data/231203/harvester_PU200_prompt_160_20.root");
  f_PU200_nonprompt = new TFile("data/231203/harvester_PU200_nonprompt_90.root");
  // noPU
  f_noPU_prompt = new TFile("data/231203/harvester_noPU_prompt_19.root");
  f_noPU_nonprompt = new TFile("data/231203/harvester_noPU_nonprompt_9.root");

  // Histograms for Barrel region
  TH1D *h_PU200_prompt_EB, *h_PU200_nonprompt_EB, *h_noPU_prompt_EB, *h_noPU_nonprompt_EB;
  TH1D *h_PU200_prompt_2sigma_EB, *h_PU200_nonprompt_2sigma_EB, *h_noPU_prompt_2sigma_EB, *h_noPU_2sigma_nonprompt_EB;
  TH1D *h_PU200_prompt_3sigma_EB, *h_PU200_nonprompt_3sigma_EB, *h_noPU_prompt_3sigma_EB, *h_noPU_3sigma_nonprompt_EB;
  TH1D *h_PU200_prompt_4sigma_EB, *h_PU200_nonprompt_4sigma_EB, *h_noPU_prompt_4sigma_EB, *h_noPU_4sigma_nonprompt_EB;
  TH1D *h_PU200_prompt_40_EB, *h_PU200_nonprompt_40_EB, *h_noPU_prompt_40_EB, *h_noPU_40_nonprompt_EB;
  TH1D *h_PU200_prompt_60_EB, *h_PU200_nonprompt_60_EB, *h_noPU_prompt_60_EB, *h_noPU_60_nonprompt_EB;
  TH1D *h_PU200_prompt_80_EB, *h_PU200_nonprompt_80_EB, *h_noPU_prompt_80_EB, *h_noPU_80_nonprompt_EB;
  TH1D *h_PU200_prompt_100_EB, *h_PU200_nonprompt_100_EB, *h_noPU_prompt_100_EB, *h_noPU_100_nonprompt_EB;
  TH1D *h_PU200_prompt_gen_EB, *h_PU200_nonprompt_gen_EB, *h_noPU_prompt_gen_EB, *h_noPU_nonprompt_gen_EB;

  // Histograms for Endcap region
  TH1D *h_PU200_prompt_EE, *h_PU200_nonprompt_EE, *h_noPU_prompt_EE, *h_noPU_nonprompt_EE;
  TH1D *h_PU200_prompt_2sigma_EE, *h_PU200_nonprompt_2sigma_EE, *h_noPU_prompt_2sigma_EE, *h_noPU_2sigma_nonprompt_EE;
  TH1D *h_PU200_prompt_3sigma_EE, *h_PU200_nonprompt_3sigma_EE, *h_noPU_prompt_3sigma_EE, *h_noPU_3sigma_nonprompt_EE;
  TH1D *h_PU200_prompt_4sigma_EE, *h_PU200_nonprompt_4sigma_EE, *h_noPU_prompt_4sigma_EE, *h_noPU_4sigma_nonprompt_EE;
  TH1D *h_PU200_prompt_40_EE, *h_PU200_nonprompt_40_EE, *h_noPU_prompt_40_EE, *h_noPU_40_nonprompt_EE;
  TH1D *h_PU200_prompt_60_EE, *h_PU200_nonprompt_60_EE, *h_noPU_prompt_60_EE, *h_noPU_60_nonprompt_EE;
  TH1D *h_PU200_prompt_80_EE, *h_PU200_nonprompt_80_EE, *h_noPU_prompt_80_EE, *h_noPU_80_nonprompt_EE;
  TH1D *h_PU200_prompt_100_EE, *h_PU200_nonprompt_100_EE, *h_noPU_prompt_100_EE, *h_noPU_100_nonprompt_EE;
  TH1D *h_PU200_prompt_gen_EE, *h_PU200_nonprompt_gen_EE, *h_noPU_prompt_gen_EE, *h_noPU_nonprompt_gen_EE;


  TString dir="/DQMData/Run 1/MTD/Run summary/MuonIso/";
  // Prompt
    // Barrel region
  h_PU200_prompt_EB 	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_Sig_EB");
  h_PU200_prompt_4sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_4sigma_Sig_EB");
  h_PU200_prompt_3sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_3sigma_Sig_EB");
  h_PU200_prompt_2sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_2sigma_Sig_EB");
  h_PU200_prompt_40_EB     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_7_Sig_EB");
  h_PU200_prompt_60_EB     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_5_Sig_EB");
  h_PU200_prompt_80_EB     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_3_Sig_EB");
  h_PU200_prompt_100_EB    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_1_Sig_EB");
  h_PU200_prompt_gen_EB	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_gen_Sig_EB");
  h_noPU_prompt_EB    	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_chIso_sum_Sig_EB");
  h_noPU_prompt_gen_EB 	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_chIso_sum_gen_Sig_EB");
    // Endcap region
  h_PU200_prompt_EE 	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_Sig_EE");
  h_PU200_prompt_4sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_4sigma_Sig_EE");
  h_PU200_prompt_3sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_3sigma_Sig_EE");
  h_PU200_prompt_2sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_2sigma_Sig_EE");
  h_PU200_prompt_40_EE     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_7_Sig_EE");
  h_PU200_prompt_60_EE     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_5_Sig_EE");
  h_PU200_prompt_80_EE     = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_3_Sig_EE");
  h_PU200_prompt_100_EE    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_MTD_1_Sig_EE");
  h_PU200_prompt_gen_EE    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_chIso_sum_gen_Sig_EE");
  h_noPU_prompt_EE    	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_chIso_sum_Sig_EE");
  h_noPU_prompt_gen_EE 	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_chIso_sum_gen_Sig_EE");
  // Nonprompt
    // Barrel region
  h_PU200_nonprompt_EB 	      = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_Bkg_EB");
  h_PU200_nonprompt_4sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_4sigma_Bkg_EB");
  h_PU200_nonprompt_3sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_3sigma_Bkg_EB");
  h_PU200_nonprompt_2sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_2sigma_Bkg_EB");
  h_PU200_nonprompt_40_EB     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_7_Bkg_EB");
  h_PU200_nonprompt_60_EB     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_5_Bkg_EB");
  h_PU200_nonprompt_80_EB     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_3_Bkg_EB");
  h_PU200_nonprompt_100_EB    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_1_Bkg_EB");
  h_PU200_nonprompt_gen_EB    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_gen_Bkg_EB");
  h_noPU_nonprompt_EB         = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_chIso_sum_Bkg_EB");
  h_noPU_nonprompt_gen_EB     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_chIso_sum_gen_Bkg_EB");
    // Endcap region
  h_PU200_nonprompt_EE 	      = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_Bkg_EE");
  h_PU200_nonprompt_4sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_4sigma_Bkg_EE");
  h_PU200_nonprompt_3sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_3sigma_Bkg_EE");
  h_PU200_nonprompt_2sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_2sigma_Bkg_EE");
  h_PU200_nonprompt_40_EE     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_7_Bkg_EE");
  h_PU200_nonprompt_60_EE     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_5_Bkg_EE");
  h_PU200_nonprompt_80_EE     = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_3_Bkg_EE");
  h_PU200_nonprompt_100_EE    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_MTD_1_Bkg_EE");
  h_PU200_nonprompt_gen_EE    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_chIso_sum_gen_Bkg_EE");
  h_noPU_nonprompt_EE         = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_chIso_sum_Bkg_EE");
  h_noPU_nonprompt_gen_EE     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_chIso_sum_gen_Bkg_EE");

  // PU200
    // Prompt
  vector<double> prompt_eff_PU200_EB={0}, prompt_norm_PU200_EB={0}, prompt_eff_PU200_EE={0}, prompt_norm_PU200_EE={0};
  vector<double> prompt_eff_PU200_2sigma_EB={0}, prompt_norm_PU200_2sigma_EB={0}, prompt_eff_PU200_2sigma_EE={0}, prompt_norm_PU200_2sigma_EE={0};
  vector<double> prompt_eff_PU200_3sigma_EB={0}, prompt_norm_PU200_3sigma_EB={0}, prompt_eff_PU200_3sigma_EE={0}, prompt_norm_PU200_3sigma_EE={0};
  vector<double> prompt_eff_PU200_4sigma_EB={0}, prompt_norm_PU200_4sigma_EB={0}, prompt_eff_PU200_4sigma_EE={0}, prompt_norm_PU200_4sigma_EE={0};
  vector<double> prompt_eff_PU200_40_EB={0}, prompt_norm_PU200_40_EB={0}, prompt_eff_PU200_40_EE={0}, prompt_norm_PU200_40_EE={0};
  vector<double> prompt_eff_PU200_60_EB={0}, prompt_norm_PU200_60_EB={0}, prompt_eff_PU200_60_EE={0}, prompt_norm_PU200_60_EE={0};
  vector<double> prompt_eff_PU200_80_EB={0}, prompt_norm_PU200_80_EB={0}, prompt_eff_PU200_80_EE={0}, prompt_norm_PU200_80_EE={0};
  vector<double> prompt_eff_PU200_100_EB={0}, prompt_norm_PU200_100_EB={0}, prompt_eff_PU200_100_EE={0}, prompt_norm_PU200_100_EE={0};
    // Nonprompt
  vector<double> nonprompt_eff_PU200_EB={0}, nonprompt_norm_PU200_EB={0}, nonprompt_eff_PU200_EE={0}, nonprompt_norm_PU200_EE={0};
  vector<double> nonprompt_eff_PU200_2sigma_EB={0}, nonprompt_norm_PU200_2sigma_EB={0}, nonprompt_eff_PU200_2sigma_EE={0}, nonprompt_norm_PU200_2sigma_EE={0};
  vector<double> nonprompt_eff_PU200_3sigma_EB={0}, nonprompt_norm_PU200_3sigma_EB={0}, nonprompt_eff_PU200_3sigma_EE={0}, nonprompt_norm_PU200_3sigma_EE={0};
  vector<double> nonprompt_eff_PU200_4sigma_EB={0}, nonprompt_norm_PU200_4sigma_EB={0}, nonprompt_eff_PU200_4sigma_EE={0}, nonprompt_norm_PU200_4sigma_EE={0};
  vector<double> nonprompt_eff_PU200_40_EB={0}, nonprompt_norm_PU200_40_EB={0}, nonprompt_eff_PU200_40_EE={0}, nonprompt_norm_PU200_40_EE={0};
  vector<double> nonprompt_eff_PU200_60_EB={0}, nonprompt_norm_PU200_60_EB={0}, nonprompt_eff_PU200_60_EE={0}, nonprompt_norm_PU200_60_EE={0};
  vector<double> nonprompt_eff_PU200_80_EB={0}, nonprompt_norm_PU200_80_EB={0}, nonprompt_eff_PU200_80_EE={0}, nonprompt_norm_PU200_80_EE={0};
  vector<double> nonprompt_eff_PU200_100_EB={0}, nonprompt_norm_PU200_100_EB={0}, nonprompt_eff_PU200_100_EE={0}, nonprompt_norm_PU200_100_EE={0};
      // Gen case
  vector<double> prompt_eff_PU200_gen_EB={0}, prompt_norm_PU200_gen_EB={0}, prompt_eff_PU200_gen_EE={0}, prompt_norm_PU200_gen_EE={0};
  vector<double> nonprompt_eff_PU200_gen_EB={0}, nonprompt_norm_PU200_gen_EB={0}, nonprompt_eff_PU200_gen_EE={0}, nonprompt_norm_PU200_gen_EE={0};
  // noPU
    // Prompt
  vector<double> prompt_eff_noPU_EB={0}, prompt_norm_noPU_EB={0}, prompt_eff_noPU_EE={0}, prompt_norm_noPU_EE={0};
  vector<double> prompt_eff_noPU_2sigma_EB={0}, prompt_norm_noPU_2sigma_EB={0}, prompt_eff_noPU_2sigma_EE={0}, prompt_norm_noPU_2sigma_EE={0};
  vector<double> prompt_eff_noPU_3sigma_EB={0}, prompt_norm_noPU_3sigma_EB={0}, prompt_eff_noPU_3sigma_EE={0}, prompt_norm_noPU_3sigma_EE={0};
  vector<double> prompt_eff_noPU_4sigma_EB={0}, prompt_norm_noPU_4sigma_EB={0}, prompt_eff_noPU_4sigma_EE={0}, prompt_norm_noPU_4sigma_EE={0};
  vector<double> prompt_eff_noPU_40_EB={0}, prompt_norm_noPU_40_EB={0}, prompt_eff_noPU_EE_40_EE={0}, prompt_norm_noPU_EE_40_EE={0};
  vector<double> prompt_eff_noPU_60_EB={0}, prompt_norm_noPU_60_EB={0}, prompt_eff_noPU_EE_60_EE={0}, prompt_norm_noPU_EE_60_EE={0};
  vector<double> prompt_eff_noPU_80_EB={0}, prompt_norm_noPU_80_EB={0}, prompt_eff_noPU_EE_80_EE={0}, prompt_norm_noPU_EE_80_EE={0};
  vector<double> prompt_eff_noPU_100_EB={0}, prompt_norm_noPU_100_EB={0}, prompt_eff_noPU_EE_100_EE={0}, prompt_norm_noPU_EE_100_EE={0};
    // Nonprompt
  vector<double> nonprompt_eff_noPU_EB={0}, nonprompt_norm_noPU_EB={0}, nonprompt_eff_noPU_EE={0}, nonprompt_norm_noPU_EE={0};
  vector<double> nonprompt_eff_noPU_2sigma_EB={0}, nonprompt_norm_noPU_2sigma_EB={0}, nonprompt_eff_noPU_2sigma_EE={0}, nonprompt_norm_noPU_2sigma_EE={0};
  vector<double> nonprompt_eff_noPU_3sigma_EB={0}, nonprompt_norm_noPU_3sigma_EB={0}, nonprompt_eff_noPU_3sigma_EE={0}, nonprompt_norm_noPU_3sigma_EE={0};
  vector<double> nonprompt_eff_noPU_4sigma_EB={0}, nonprompt_norm_noPU_4sigma_EB={0}, nonprompt_eff_noPU_4sigma_EE={0}, nonprompt_norm_noPU_4sigma_EE={0};
  vector<double> nonprompt_eff_noPU_40_EB={0}, nonprompt_norm_noPU_40_EB={0}, nonprompt_eff_noPU_40_EE={0}, nonprompt_norm_noPU_40_EE={0};
  vector<double> nonprompt_eff_noPU_60_EB={0}, nonprompt_norm_noPU_60_EB={0}, nonprompt_eff_noPU_60_EE={0}, nonprompt_norm_noPU_60_EE={0};
  vector<double> nonprompt_eff_noPU_80_EB={0}, nonprompt_norm_noPU_80_EB={0}, nonprompt_eff_noPU_80_EE={0}, nonprompt_norm_noPU_80_EE={0};
  vector<double> nonprompt_eff_noPU_100_EB={0}, nonprompt_norm_noPU_100_EB={0}, nonprompt_eff_noPU_100_EE={0}, nonprompt_norm_noPU_100_EE={0};
    // Gen case
  vector<double> prompt_eff_noPU_gen_EB={0}, prompt_norm_noPU_gen_EB={0}, prompt_eff_noPU_gen_EE={0}, prompt_norm_noPU_gen_EE={0};
  vector<double> nonprompt_eff_noPU_gen_EB={0}, nonprompt_norm_noPU_gen_EB={0}, nonprompt_eff_noPU_gen_EE={0}, nonprompt_norm_noPU_gen_EE={0};


  //////////////////////////
  // Calculate efficiency //
  //////////////////////////
  for(unsigned int i=0; i<h_PU200_prompt_EB->GetNbinsX(); i++) {
    if(zeroiso==true) {        // Include the cases where the muon is already isolated
    // Prompt
      // Barrel region
        // efficiency
      prompt_eff_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(1,i+1)/h_PU200_prompt_EB->Integral(1,-1));
      prompt_eff_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(1,i+1)/h_noPU_prompt_EB->Integral(1,-1));
      prompt_eff_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(1,i+1)/h_PU200_prompt_40_EB->Integral(1,-1));
      prompt_eff_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(1,i+1)/h_PU200_prompt_60_EB->Integral(1,-1));
      prompt_eff_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(1,i+1)/h_PU200_prompt_80_EB->Integral(1,-1));
      prompt_eff_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(1,i+1)/h_PU200_prompt_100_EB->Integral(1,-1));
      prompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(1,i+1)/h_PU200_prompt_2sigma_EB->Integral(1,-1));
      prompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(1,i+1)/h_PU200_prompt_3sigma_EB->Integral(1,-1));
      prompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(1,i+1)/h_PU200_prompt_4sigma_EB->Integral(1,-1));
      prompt_eff_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(1,i+1)/h_PU200_prompt_gen_EB->Integral(1,-1));
      prompt_eff_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(1,i+1)/h_noPU_prompt_gen_EB->Integral(1,-1));
        // normalization
      prompt_norm_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(1,i+1));
      prompt_norm_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(1,i+1));
      prompt_norm_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(1,i+1));
      prompt_norm_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(1,i+1));
      prompt_norm_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(1,i+1));
      prompt_norm_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(1,i+1));
      prompt_norm_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(1,i+1));
      prompt_norm_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(1,i+1));
      prompt_norm_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(1,i+1));
      prompt_norm_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(1,i+1));
      prompt_norm_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(1,i+1));
      // Endcap region
        // efficiency
      prompt_eff_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(1,i+1)/h_PU200_prompt_EE->Integral(1,-1));
      prompt_eff_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(1,i+1)/h_noPU_prompt_EE->Integral(1,-1));
      prompt_eff_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(1,i+1)/h_PU200_prompt_40_EE->Integral(1,-1));
      prompt_eff_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(1,i+1)/h_PU200_prompt_60_EE->Integral(1,-1));
      prompt_eff_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(1,i+1)/h_PU200_prompt_80_EE->Integral(1,-1));
      prompt_eff_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(1,i+1)/h_PU200_prompt_100_EE->Integral(1,-1));
      prompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(1,i+1)/h_PU200_prompt_2sigma_EE->Integral(1,-1));
      prompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(1,i+1)/h_PU200_prompt_3sigma_EE->Integral(1,-1));
      prompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(1,i+1)/h_PU200_prompt_4sigma_EE->Integral(1,-1));
      prompt_eff_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(1,i+1)/h_PU200_prompt_gen_EE->Integral(1,-1));
      prompt_eff_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(1,i+1)/h_noPU_prompt_gen_EE->Integral(1,-1));
        // normalization
      prompt_norm_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(1,i+1));
      prompt_norm_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(1,i+1));
      prompt_norm_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(1,i+1));
      prompt_norm_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(1,i+1));
      prompt_norm_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(1,i+1));
      prompt_norm_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(1,i+1));
      prompt_norm_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(1,i+1));
      prompt_norm_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(1,i+1));
      prompt_norm_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(1,i+1));
      prompt_norm_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(1,i+1));
      prompt_norm_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(1,i+1));
    // nonprompt
      // Barrel region
        // efficiency
      nonprompt_eff_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(1,i+1)/h_PU200_nonprompt_EB->Integral(1,-1));
      nonprompt_eff_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(1,i+1)/h_noPU_nonprompt_EB->Integral(1,-1));
      nonprompt_eff_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(1,i+1)/h_PU200_nonprompt_40_EB->Integral(1,-1));
      nonprompt_eff_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(1,i+1)/h_PU200_nonprompt_60_EB->Integral(1,-1));
      nonprompt_eff_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(1,i+1)/h_PU200_nonprompt_80_EB->Integral(1,-1));
      nonprompt_eff_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(1,i+1)/h_PU200_nonprompt_100_EB->Integral(1,-1));
      nonprompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EB->Integral(1,-1));
      nonprompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EB->Integral(1,-1));
      nonprompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EB->Integral(1,-1));
      nonprompt_eff_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(1,i+1)/h_PU200_nonprompt_gen_EB->Integral(1,-1));
      nonprompt_eff_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(1,i+1)/h_noPU_nonprompt_gen_EB->Integral(1,-1));
        // normalization
      nonprompt_norm_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(1,i+1));
      nonprompt_norm_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(1,i+1));
      nonprompt_norm_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(1,i+1));
      nonprompt_norm_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(1,i+1));
      nonprompt_norm_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(1,i+1));
      nonprompt_norm_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(1,i+1));
      nonprompt_norm_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(1,i+1));
      nonprompt_norm_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(1,i+1));
      nonprompt_norm_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(1,i+1));
      nonprompt_norm_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(1,i+1));
      nonprompt_norm_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(1,i+1));
      // Endcap region
        // efficiency
      nonprompt_eff_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(1,i+1)/h_PU200_nonprompt_EE->Integral(1,-1));
      nonprompt_eff_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(1,i+1)/h_noPU_nonprompt_EE->Integral(1,-1));
      nonprompt_eff_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(1,i+1)/h_PU200_nonprompt_40_EE->Integral(1,-1));
      nonprompt_eff_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(1,i+1)/h_PU200_nonprompt_60_EE->Integral(1,-1));
      nonprompt_eff_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(1,i+1)/h_PU200_nonprompt_80_EE->Integral(1,-1));
      nonprompt_eff_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(1,i+1)/h_PU200_nonprompt_100_EE->Integral(1,-1));
      nonprompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EE->Integral(1,-1));
      nonprompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EE->Integral(1,-1));
      nonprompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EE->Integral(1,-1));
      nonprompt_eff_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(1,i+1)/h_PU200_nonprompt_gen_EE->Integral(1,-1));
      nonprompt_eff_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(1,i+1)/h_noPU_nonprompt_gen_EE->Integral(1,-1));
        // normalization
      nonprompt_norm_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(1,i+1));
      nonprompt_norm_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(1,i+1));
      nonprompt_norm_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(1,i+1));
      nonprompt_norm_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(1,i+1));
      nonprompt_norm_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(1,i+1));
      nonprompt_norm_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(1,i+1));
      nonprompt_norm_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(1,i+1));
      nonprompt_norm_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(1,i+1));
      nonprompt_norm_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(1,i+1));
      nonprompt_norm_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(1,i+1));
      nonprompt_norm_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(1,i+1));
    }
    else {             // Not include the cases where the muon is already isolated
    // Prompt
      // Barrel region
        // efficiency
      prompt_eff_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(2,2+i)/h_PU200_prompt_EB->Integral(2,-1));
      prompt_eff_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(2,2+i)/h_noPU_prompt_EB->Integral(2,-1));
      prompt_eff_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(2,2+i)/h_PU200_prompt_40_EB->Integral(2,-1));
      prompt_eff_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(2,2+i)/h_PU200_prompt_60_EB->Integral(2,-1));
      prompt_eff_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(2,2+i)/h_PU200_prompt_80_EB->Integral(2,-1));
      prompt_eff_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(2,2+i)/h_PU200_prompt_100_EB->Integral(2,-1));
      prompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(2,2+i)/h_PU200_prompt_2sigma_EB->Integral(2,-1));
      prompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(2,2+i)/h_PU200_prompt_3sigma_EB->Integral(2,-1));
      prompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(2,2+i)/h_PU200_prompt_4sigma_EB->Integral(2,-1));
      prompt_eff_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(2,2+i)/h_PU200_prompt_gen_EB->Integral(2,-1));
      prompt_eff_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(2,2+i)/h_noPU_prompt_gen_EB->Integral(2,-1));
        // normalization
      prompt_norm_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(2,2+i));
      prompt_norm_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(2,2+i));
      prompt_norm_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(2,2+i));
      prompt_norm_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(2,2+i));
      prompt_norm_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(2,2+i));
      prompt_norm_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(2,2+i));
      prompt_norm_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(2,2+i));
      prompt_norm_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(2,2+i));
      prompt_norm_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(2,2+i));
      prompt_norm_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(2,2+i));
      prompt_norm_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(2,2+i));
      // Endcap region
        // efficiency
      prompt_eff_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(2,2+i)/h_PU200_prompt_EE->Integral(2,-1));
      prompt_eff_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(2,2+i)/h_noPU_prompt_EE->Integral(2,-1));
      prompt_eff_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(2,2+i)/h_PU200_prompt_40_EE->Integral(2,-1));
      prompt_eff_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(2,2+i)/h_PU200_prompt_60_EE->Integral(2,-1));
      prompt_eff_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(2,2+i)/h_PU200_prompt_80_EE->Integral(2,-1));
      prompt_eff_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(2,2+i)/h_PU200_prompt_100_EE->Integral(2,-1));
      prompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(2,2+i)/h_PU200_prompt_2sigma_EE->Integral(2,-1));
      prompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(2,2+i)/h_PU200_prompt_3sigma_EE->Integral(2,-1));
      prompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(2,2+i)/h_PU200_prompt_4sigma_EE->Integral(2,-1));
      prompt_eff_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(2,2+i)/h_PU200_prompt_gen_EE->Integral(2,-1));
      prompt_eff_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(2,2+i)/h_noPU_prompt_gen_EE->Integral(2,-1));
        // normalization
      prompt_norm_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(2,2+i));
      prompt_norm_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(2,2+i));
      prompt_norm_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(2,2+i));
      prompt_norm_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(2,2+i));
      prompt_norm_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(2,2+i));
      prompt_norm_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(2,2+i));
      prompt_norm_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(2,2+i));
      prompt_norm_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(2,2+i));
      prompt_norm_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(2,2+i));
      prompt_norm_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(2,2+i));
      prompt_norm_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(2,2+i));
    // Nonprompt
      // Barrel region
        // efficiency
      nonprompt_eff_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(2,2+i)/h_PU200_nonprompt_EB->Integral(2,-1));
      nonprompt_eff_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(2,2+i)/h_noPU_nonprompt_EB->Integral(2,-1));
      nonprompt_eff_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(2,2+i)/h_PU200_nonprompt_40_EB->Integral(2,-1));
      nonprompt_eff_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(2,2+i)/h_PU200_nonprompt_60_EB->Integral(2,-1));
      nonprompt_eff_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(2,2+i)/h_PU200_nonprompt_80_EB->Integral(2,-1));
      nonprompt_eff_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(2,2+i)/h_PU200_nonprompt_100_EB->Integral(2,-1));
      nonprompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(2,2+i)/h_PU200_nonprompt_2sigma_EB->Integral(2,-1));
      nonprompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(2,2+i)/h_PU200_nonprompt_3sigma_EB->Integral(2,-1));
      nonprompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(2,2+i)/h_PU200_nonprompt_4sigma_EB->Integral(2,-1));
      nonprompt_eff_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(2,2+i)/h_PU200_nonprompt_gen_EB->Integral(2,-1));
      nonprompt_eff_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(2,2+i)/h_noPU_nonprompt_gen_EB->Integral(2,-1));
        // normalization
      nonprompt_norm_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(2,2+i));
      nonprompt_norm_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(2,2+i));
      nonprompt_norm_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(2,2+i));
      nonprompt_norm_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(2,2+i));
      nonprompt_norm_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(2,2+i));
      nonprompt_norm_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(2,2+i));
      nonprompt_norm_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(2,2+i));
      nonprompt_norm_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(2,2+i));
      nonprompt_norm_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(2,2+i));
      nonprompt_norm_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(2,2+i));
      nonprompt_norm_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(2,2+i));
      // Endcap region
        // efficiency
      nonprompt_eff_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(2,2+i)/h_PU200_nonprompt_EE->Integral(2,-1));
      nonprompt_eff_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(2,2+i)/h_noPU_nonprompt_EE->Integral(2,-1));
      nonprompt_eff_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(2,2+i)/h_PU200_nonprompt_40_EE->Integral(2,-1));
      nonprompt_eff_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(2,2+i)/h_PU200_nonprompt_60_EE->Integral(2,-1));
      nonprompt_eff_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(2,2+i)/h_PU200_nonprompt_80_EE->Integral(2,-1));
      nonprompt_eff_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(2,2+i)/h_PU200_nonprompt_100_EE->Integral(2,-1));
      nonprompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(2,2+i)/h_PU200_nonprompt_2sigma_EE->Integral(2,-1));
      nonprompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(2,2+i)/h_PU200_nonprompt_3sigma_EE->Integral(2,-1));
      nonprompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(2,2+i)/h_PU200_nonprompt_4sigma_EE->Integral(2,-1));
      nonprompt_eff_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(2,2+i)/h_PU200_nonprompt_gen_EE->Integral(2,-1));
      nonprompt_eff_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(2,2+i)/h_noPU_nonprompt_gen_EE->Integral(2,-1));
        // normalization
      nonprompt_norm_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(2,2+i));
      nonprompt_norm_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(2,2+i));
      nonprompt_norm_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(2,2+i));
      nonprompt_norm_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(2,2+i));
      nonprompt_norm_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(2,2+i));
      nonprompt_norm_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(2,2+i));
      nonprompt_norm_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(2,2+i));
      nonprompt_norm_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(2,2+i));
      nonprompt_norm_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(2,2+i));
      nonprompt_norm_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(2,2+i));
      nonprompt_norm_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(2,2+i));
    }
  }

  // Define TGraph
    // Barrel region
  TGraph* gr_PU200_EB = new TGraph(h_PU200_prompt_EB->GetNbinsX(), &prompt_eff_PU200_EB[0], &nonprompt_eff_PU200_EB[0]);
  TGraph* gr_PU200_40_EB = new TGraph(h_PU200_prompt_40_EB->GetNbinsX(), &prompt_eff_PU200_40_EB[0], &nonprompt_eff_PU200_40_EB[0]);
  TGraph* gr_PU200_60_EB = new TGraph(h_PU200_prompt_60_EB->GetNbinsX(), &prompt_eff_PU200_60_EB[0], &nonprompt_eff_PU200_60_EB[0]);
  TGraph* gr_PU200_80_EB = new TGraph(h_PU200_prompt_80_EB->GetNbinsX(), &prompt_eff_PU200_80_EB[0], &nonprompt_eff_PU200_80_EB[0]);
  TGraph* gr_PU200_100_EB = new TGraph(h_PU200_prompt_100_EB->GetNbinsX(), &prompt_eff_PU200_100_EB[0], &nonprompt_eff_PU200_100_EB[0]);
  TGraph* gr_PU200_2sigma_EB = new TGraph(h_PU200_prompt_2sigma_EB->GetNbinsX(), &prompt_eff_PU200_2sigma_EB[0], &nonprompt_eff_PU200_2sigma_EB[0]);
  TGraph* gr_PU200_3sigma_EB = new TGraph(h_PU200_prompt_3sigma_EB->GetNbinsX(), &prompt_eff_PU200_3sigma_EB[0], &nonprompt_eff_PU200_3sigma_EB[0]);
  TGraph* gr_PU200_4sigma_EB = new TGraph(h_PU200_prompt_4sigma_EB->GetNbinsX(), &prompt_eff_PU200_4sigma_EB[0], &nonprompt_eff_PU200_4sigma_EB[0]);
  TGraph* gr_noPU_EB = new TGraph(h_noPU_prompt_EB->GetNbinsX(), &prompt_eff_noPU_EB[0], &nonprompt_eff_noPU_EB[0]);
  TGraph* gr_PU200_gen_EB = new TGraph(h_PU200_prompt_gen_EB->GetNbinsX(), &prompt_eff_PU200_gen_EB[0], &nonprompt_eff_PU200_gen_EB[0]);
  TGraph* gr_noPU_gen_EB = new TGraph(h_noPU_prompt_gen_EB->GetNbinsX(), &prompt_eff_noPU_gen_EB[0], &nonprompt_eff_noPU_gen_EB[0]);
    // Endcap region
  TGraph* gr_PU200_EE = new TGraph(h_PU200_prompt_EE->GetNbinsX(), &prompt_eff_PU200_EE[0], &nonprompt_eff_PU200_EE[0]);
  TGraph* gr_PU200_40_EE = new TGraph(h_PU200_prompt_40_EE->GetNbinsX(), &prompt_eff_PU200_40_EE[0], &nonprompt_eff_PU200_40_EE[0]);
  TGraph* gr_PU200_60_EE = new TGraph(h_PU200_prompt_60_EE->GetNbinsX(), &prompt_eff_PU200_60_EE[0], &nonprompt_eff_PU200_60_EE[0]);
  TGraph* gr_PU200_80_EE = new TGraph(h_PU200_prompt_80_EE->GetNbinsX(), &prompt_eff_PU200_80_EE[0], &nonprompt_eff_PU200_80_EE[0]);
  TGraph* gr_PU200_100_EE = new TGraph(h_PU200_prompt_100_EE->GetNbinsX(), &prompt_eff_PU200_100_EE[0], &nonprompt_eff_PU200_100_EE[0]);
  TGraph* gr_PU200_2sigma_EE = new TGraph(h_PU200_prompt_2sigma_EE->GetNbinsX(), &prompt_eff_PU200_2sigma_EE[0], &nonprompt_eff_PU200_2sigma_EE[0]);
  TGraph* gr_PU200_3sigma_EE = new TGraph(h_PU200_prompt_3sigma_EE->GetNbinsX(), &prompt_eff_PU200_3sigma_EE[0], &nonprompt_eff_PU200_3sigma_EE[0]);
  TGraph* gr_PU200_4sigma_EE = new TGraph(h_PU200_prompt_4sigma_EE->GetNbinsX(), &prompt_eff_PU200_4sigma_EE[0], &nonprompt_eff_PU200_4sigma_EE[0]);
  TGraph* gr_noPU_EE = new TGraph(h_noPU_prompt_EE->GetNbinsX(), &prompt_eff_noPU_EE[0], &nonprompt_eff_noPU_EE[0]);
  TGraph* gr_PU200_gen_EE = new TGraph(h_PU200_prompt_gen_EE->GetNbinsX(), &prompt_eff_PU200_gen_EE[0], &nonprompt_eff_PU200_gen_EE[0]);
  TGraph* gr_noPU_gen_EE = new TGraph(h_noPU_prompt_gen_EE->GetNbinsX(), &prompt_eff_noPU_gen_EE[0], &nonprompt_eff_noPU_gen_EE[0]);

  // Remove the value at dump point (0,0)
  gr_PU200_EB->RemovePoint(0); gr_PU200_40_EB->RemovePoint(0); gr_PU200_60_EB->RemovePoint(0); gr_PU200_80_EB->RemovePoint(0); gr_PU200_100_EB->RemovePoint(0); gr_PU200_2sigma_EB->RemovePoint(0); gr_PU200_3sigma_EB->RemovePoint(0); gr_PU200_4sigma_EB->RemovePoint(0); gr_noPU_EB->RemovePoint(0); gr_PU200_gen_EB->RemovePoint(0); gr_noPU_gen_EB->RemovePoint(0);
  gr_PU200_EE->RemovePoint(0); gr_PU200_40_EE->RemovePoint(0); gr_PU200_60_EE->RemovePoint(0); gr_PU200_80_EE->RemovePoint(0); gr_PU200_100_EE->RemovePoint(0); gr_PU200_2sigma_EE->RemovePoint(0); gr_PU200_3sigma_EE->RemovePoint(0); gr_PU200_4sigma_EE->RemovePoint(0); gr_noPU_EE->RemovePoint(0); gr_PU200_gen_EE->RemovePoint(0); gr_noPU_gen_EE->RemovePoint(0);


  ///////////////
  // Cosmetics //
  ///////////////
    // Barrel region
  gr_PU200_EB->SetTitle("ROC curve in barrel region"); gr_noPU_EB->SetTitle("ROC  curve in barrel region");
  gr_PU200_EB->GetXaxis()->SetTitle("Prompt efficiency"); gr_noPU_EB->GetXaxis()->SetTitle("Prompt efficiency");
  gr_PU200_EB->GetYaxis()->SetTitle("Non-prompt efficiency"); gr_noPU_EB->GetYaxis()->SetTitle("Non-prompt efficiency");
  gr_PU200_EB->SetLineWidth(2); gr_PU200_40_EB->SetLineWidth(2); gr_PU200_60_EB->SetLineWidth(2); gr_PU200_80_EB->SetLineWidth(2); gr_PU200_100_EB->SetLineWidth(2);
  gr_PU200_2sigma_EB->SetLineWidth(2); gr_PU200_3sigma_EB->SetLineWidth(2); gr_PU200_4sigma_EB->SetLineWidth(2);
  gr_noPU_EB->SetLineWidth(2); gr_PU200_gen_EB->SetLineWidth(2); gr_noPU_gen_EB->SetLineWidth(2);
  gr_PU200_EB->SetLineColor(kBlack); gr_PU200_40_EB->SetLineColor(kRed); gr_PU200_60_EB->SetLineColor(kGreen); gr_PU200_80_EB->SetLineColor(kBlue); gr_PU200_100_EB->SetLineColor(kMagenta);
  gr_PU200_2sigma_EB->SetLineColor(kRed); gr_PU200_3sigma_EB->SetLineColor(kGreen); gr_PU200_4sigma_EB->SetLineColor(kBlue);
  gr_noPU_EB->SetLineColor(kGray); gr_PU200_gen_EB->SetLineColor(kBlack); gr_noPU_gen_EB->SetLineColor(kGray);
  gr_PU200_gen_EB->SetLineStyle(7); gr_noPU_gen_EB->SetLineStyle(7);

    // Endcap region
  gr_PU200_EE->SetTitle("ROC curve in endcap region"); gr_noPU_EE->SetTitle("ROC curve in endcap region");
  gr_PU200_EE->GetXaxis()->SetTitle("Prompt efficiency"); gr_noPU_EE->GetXaxis()->SetTitle("Prompt efficiency");
  gr_PU200_EE->GetYaxis()->SetTitle("Non-prompt efficiency"); gr_noPU_EE->GetYaxis()->SetTitle("Non-prompt efficiency");
  gr_PU200_EE->SetLineWidth(2); gr_PU200_40_EE->SetLineWidth(2); gr_PU200_60_EE->SetLineWidth(2); gr_PU200_80_EE->SetLineWidth(2); gr_PU200_100_EE->SetLineWidth(2);
  gr_PU200_2sigma_EE->SetLineWidth(2); gr_PU200_3sigma_EE->SetLineWidth(2); gr_PU200_4sigma_EE->SetLineWidth(2);
  gr_noPU_EE->SetLineWidth(2); gr_PU200_gen_EE->SetLineWidth(2); gr_noPU_gen_EE->SetLineWidth(2);
  gr_PU200_EE->SetLineColor(kBlack); gr_PU200_40_EE->SetLineColor(kRed); gr_PU200_60_EE->SetLineColor(kGreen); gr_PU200_80_EE->SetLineColor(kBlue); gr_PU200_100_EE->SetLineColor(kMagenta);
  gr_PU200_2sigma_EE->SetLineColor(kRed); gr_PU200_3sigma_EE->SetLineColor(kGreen); gr_PU200_4sigma_EE->SetLineColor(kBlue);
  gr_noPU_EE->SetLineColor(kGray); gr_PU200_gen_EE->SetLineColor(kBlack); gr_noPU_gen_EE->SetLineColor(kGray);
  gr_PU200_gen_EE->SetLineStyle(7); gr_noPU_gen_EE->SetLineStyle(7);


  /////////////
  // Legends //
  /////////////
  // Barrel region
  TLegend* leg_EB_dt = new TLegend(0.15, 0.65, 0.48, 0.85);
  leg_EB_dt->AddEntry(gr_PU200_gen_EB, "gen PU200");
  leg_EB_dt->AddEntry(gr_PU200_EB, "no MTD PU200");
  leg_EB_dt->AddEntry(gr_PU200_40_EB, "40ps PU200");
//  leg_EB_dt->AddEntry(gr_PU200_60_EB, "60ps PU200");
//  leg_EB_dt->AddEntry(gr_PU200_80_EB, "80ps PU200");
//  leg_EB_dt->AddEntry(gr_PU200_100_EB, "100ps PU200");
  leg_EB_dt->AddEntry(gr_noPU_gen_EB, "gen noPU");
  leg_EB_dt->AddEntry(gr_noPU_EB, "noPU");
  leg_EB_dt->SetTextSize(0.03);

  TLegend* leg_EB_sigma = new TLegend(0.15, 0.65, 0.48, 0.85);
  leg_EB_sigma->AddEntry(gr_PU200_gen_EB, "gen PU200");
  leg_EB_sigma->AddEntry(gr_PU200_EB, "no MTD PU200");
  leg_EB_sigma->AddEntry(gr_PU200_2sigma_EB, "2sigma PU200");
//  leg_EB_sigma->AddEntry(gr_PU200_3sigma_EB, "3sigma PU200");
//  leg_EB_sigma->AddEntry(gr_PU200_4sigma_EB, "4sigma PU200");
  leg_EB_sigma->AddEntry(gr_noPU_gen_EB, "gen noPU");
  leg_EB_sigma->AddEntry(gr_noPU_EB, "noPU");
  leg_EB_sigma->SetTextSize(0.03);

  // Endcap region
  TLegend* leg_EE_dt = new TLegend(0.15, 0.65, 0.48, 0.85);
  leg_EE_dt->AddEntry(gr_PU200_gen_EE, "gen PU200");
  leg_EE_dt->AddEntry(gr_PU200_EE, "no MTD PU200");
  leg_EE_dt->AddEntry(gr_PU200_40_EE, "40ps PU200");
//  leg_EE_dt->AddEntry(gr_PU200_60_EE, "60ps PU200");
//  leg_EE_dt->AddEntry(gr_PU200_80_EE, "80ps PU200");
//  leg_EE_dt->AddEntry(gr_PU200_100_EE, "100ps PU200");
  leg_EE_dt->AddEntry(gr_noPU_gen_EE, "gen noPU");
  leg_EE_dt->AddEntry(gr_noPU_EE, "noPU");
  leg_EE_dt->SetTextSize(0.03);

  TLegend* leg_EE_sigma = new TLegend(0.15, 0.65, 0.48, 0.85);
  leg_EE_sigma->AddEntry(gr_PU200_gen_EE, "gen PU200");
  leg_EE_sigma->AddEntry(gr_PU200_EE, "no MTD PU200");
  leg_EE_sigma->AddEntry(gr_PU200_2sigma_EE, "2sigma PU200");
//  leg_EE_sigma->AddEntry(gr_PU200_3sigma_EE, "3sigma PU200");
//  leg_EE_sigma->AddEntry(gr_PU200_4sigma_EE, "4sigma PU200");
  leg_EE_sigma->AddEntry(gr_noPU_gen_EE, "gen noPU");
  leg_EE_sigma->AddEntry(gr_noPU_EE, " noPU");
  leg_EE_sigma->SetTextSize(0.03);


  ////////////////
  // Draw plots //
  ////////////////
  // Barrel region
  TCanvas* c_dt_EB = new TCanvas("c_dt_EB", "c_dt_EB", 1500, 1500);
  c_dt_EB->cd();
  c_dt_EB->SetGrid();
  gr_PU200_EB->GetXaxis()->SetRangeUser(0.,1.1);
  gr_PU200_EB->Draw("AL");
  gr_PU200_40_EB->Draw("same");
//  gr_PU200_60_EB->Draw("same");
//  gr_PU200_80_EB->Draw("same");
//  gr_PU200_100_EB->Draw("same");
  gr_noPU_EB->Draw("same");
  gr_PU200_gen_EB->Draw("same");
  gr_noPU_gen_EB->Draw("same");
  leg_EB_dt->Draw();
  c_dt_EB->Print("plots/roc_dt_EB_pt.pdf");

  TCanvas* c_sigma_EB = new TCanvas("c_sigma_EB", "c_sigma_EB", 1500, 1500);
  c_sigma_EB->cd();
  c_sigma_EB->SetGrid();
  gr_PU200_EB->GetXaxis()->SetRangeUser(0.,1.1);
  gr_PU200_EB->Draw("AL");
  gr_PU200_2sigma_EB->Draw("same");
//  gr_PU200_3sigma_EB->Draw("same");
//  gr_PU200_4sigma_EB->Draw("same");
  gr_noPU_EB->Draw("same");
  gr_PU200_gen_EB->Draw("same");
  gr_noPU_gen_EB->Draw("same");
  leg_EB_sigma->Draw();
  c_sigma_EB->Print("plots/roc_sigma_EB_pt.pdf");

  // Endcap region
  TCanvas* c_dt_EE = new TCanvas("c_dt_EE", "c_dt_EE", 1500, 1500);
  c_dt_EE->cd();
  c_dt_EE->SetGrid();
  gr_PU200_EE->GetXaxis()->SetRangeUser(0.,1.1);
  gr_PU200_EE->Draw("AL");
  gr_PU200_40_EE->Draw("same");
//  gr_PU200_60_EE->Draw("same");
//  gr_PU200_80_EE->Draw("same");
//  gr_PU200_100_EE->Draw("same");
  gr_noPU_EE->Draw("same");
  gr_PU200_gen_EE->Draw("same");
  gr_noPU_gen_EE->Draw("same");
  leg_EE_dt->Draw();
  c_dt_EE->Print("plots/roc_dt_EE_pt.pdf");

  TCanvas* c_sigma_EE = new TCanvas("c_sigma_EE", "c_sigma_EE", 1500, 1500);
  c_sigma_EE->cd();
  c_sigma_EE->SetGrid();
  gr_PU200_EE->GetXaxis()->SetRangeUser(0.,1.1);
  gr_PU200_EE->Draw("AL");
  gr_PU200_2sigma_EE->Draw("same");
//  gr_PU200_3sigma_EE->Draw("same");
//  gr_PU200_4sigma_EE->Draw("same");
  gr_noPU_EE->Draw("same");
  gr_PU200_gen_EE->Draw("same");
  gr_noPU_gen_EE->Draw("same");
  leg_EE_sigma->Draw();
  c_sigma_EE->Print("plots/roc_sigma_EE_pt.pdf");
}


void pTeff() {
  gStyle->SetOptStat(0);
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  // PU200
  f_PU200_prompt = new TFile("data/231203/harvester_PU200_prompt_160_20.root");
  f_PU200_nonprompt = new TFile("data/231203/harvester_PU200_nonprompt_90.root");
  // noPU
  f_noPU_prompt = new TFile("data/231203/harvester_noPU_prompt_19.root");
  f_noPU_nonprompt = new TFile("data/231203/harvester_noPU_nonprompt_9.root");

  // Histograms for Barrel region
  TH1D *h_PU200_prompt_pT_EB, *h_PU200_nonprompt_pT_EB, *h_noPU_prompt_pT_EB, *h_noPU_nonprompt_pT_EB;
  TH1D *h_PU200_prompt_2sigma_pT_EB, *h_PU200_nonprompt_2sigma_pT_EB, *h_noPU_prompt_2sigma_pT_EB, *h_noPU_2sigma_nonprompt_pT_EB;
  TH1D *h_PU200_prompt_3sigma_pT_EB, *h_PU200_nonprompt_3sigma_pT_EB, *h_noPU_prompt_3sigma_pT_EB, *h_noPU_3sigma_nonprompt_pT_EB;
  TH1D *h_PU200_prompt_4sigma_pT_EB, *h_PU200_nonprompt_4sigma_pT_EB, *h_noPU_prompt_4sigma_pT_EB, *h_noPU_4sigma_nonprompt_pT_EB;
  TH1D *h_PU200_prompt_40_pT_EB, *h_PU200_nonprompt_40_pT_EB, *h_noPU_prompt_40_pT_EB, *h_noPU_40_nonprompt_pT_EB;
  TH1D *h_PU200_prompt_60_pT_EB, *h_PU200_nonprompt_60_pT_EB, *h_noPU_prompt_60_pT_EB, *h_noPU_60_nonprompt_pT_EB;
  TH1D *h_PU200_prompt_80_pT_EB, *h_PU200_nonprompt_80_pT_EB, *h_noPU_prompt_80_pT_EB, *h_noPU_80_nonprompt_pT_EB;
  TH1D *h_PU200_prompt_100_pT_EB, *h_PU200_nonprompt_100_pT_EB, *h_noPU_prompt_100_pT_EB, *h_noPU_100_nonprompt_pT_EB;
  TH1D *h_PU200_prompt_gen_pT_EB, *h_PU200_nonprompt_gen_pT_EB, *h_noPU_prompt_gen_pT_EB, *h_noPU_nonprompt_gen_pT_EB;

  // Histograms for Endcap region
  TH1D *h_PU200_prompt_pT_EE, *h_PU200_nonprompt_pT_EE, *h_noPU_prompt_pT_EE, *h_noPU_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_2sigma_pT_EE, *h_PU200_nonprompt_2sigma_pT_EE, *h_noPU_prompt_2sigma_pT_EE, *h_noPU_2sigma_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_3sigma_pT_EE, *h_PU200_nonprompt_3sigma_pT_EE, *h_noPU_prompt_3sigma_pT_EE, *h_noPU_3sigma_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_4sigma_pT_EE, *h_PU200_nonprompt_4sigma_pT_EE, *h_noPU_prompt_4sigma_pT_EE, *h_noPU_4sigma_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_40_pT_EE, *h_PU200_nonprompt_40_pT_EE, *h_noPU_prompt_40_pT_EE, *h_noPU_40_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_60_pT_EE, *h_PU200_nonprompt_60_pT_EE, *h_noPU_prompt_60_pT_EE, *h_noPU_60_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_80_pT_EE, *h_PU200_nonprompt_80_pT_EE, *h_noPU_prompt_80_pT_EE, *h_noPU_80_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_100_pT_EE, *h_PU200_nonprompt_100_pT_EE, *h_noPU_prompt_100_pT_EE, *h_noPU_100_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_gen_pT_EE, *h_PU200_nonprompt_gen_pT_EE, *h_noPU_prompt_gen_pT_EE, *h_noPU_nonprompt_gen_pT_EE;

  TString dir="/DQMData/Run 1/MTD/Run summary/MuonIso/";
  // Prompt
    // Barrel region
  h_PU200_prompt_pT_EB 	   	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffnoMTD_Sig_EB");
  h_PU200_prompt_4sigma_pT_EB   = (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_4sigma_Sig_EB");
  h_PU200_prompt_3sigma_pT_EB 	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_3sigma_Sig_EB");
  h_PU200_prompt_2sigma_pT_EB 	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_2sigma_Sig_EB");
  h_PU200_prompt_40_pT_EB     	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_7_Sig_EB");
  h_PU200_prompt_60_pT_EB     	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_5_Sig_EB");
  h_PU200_prompt_80_pT_EB     	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_3_Sig_EB");
  h_PU200_prompt_100_pT_EB     	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_1_Sig_EB");
  h_PU200_prompt_gen_pT_EB	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_gen_Sig_EB");
  h_noPU_prompt_pT_EB 	   	= (TH1D*)f_noPU_prompt->Get(dir+"pTeffnoMTD_Sig_EB");            // noPU
  h_noPU_prompt_gen_pT_EB	= (TH1D*)f_noPU_prompt->Get(dir+"pTeffMTD_gen_Sig_EB");          // noPU
    // Endcap region
  h_PU200_prompt_pT_EE 	   	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffnoMTD_Sig_EE");
  h_PU200_prompt_4sigma_pT_EE   = (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_4sigma_Sig_EE");
  h_PU200_prompt_3sigma_pT_EE 	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_3sigma_Sig_EE");
  h_PU200_prompt_2sigma_pT_EE 	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_2sigma_Sig_EE");
  h_PU200_prompt_40_pT_EE     	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_7_Sig_EE");
  h_PU200_prompt_60_pT_EE     	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_5_Sig_EE");
  h_PU200_prompt_80_pT_EE     	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_3_Sig_EE");
  h_PU200_prompt_100_pT_EE     	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_1_Sig_EE");
  h_PU200_prompt_gen_pT_EE	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_gen_Sig_EE");
  h_noPU_prompt_pT_EE 	   	= (TH1D*)f_noPU_prompt->Get(dir+"pTeffnoMTD_Sig_EE");            // noPU
  h_noPU_prompt_gen_pT_EE	= (TH1D*)f_noPU_prompt->Get(dir+"pTeffMTD_gen_Sig_EE");          // noPU
  // Nonprompt
    // Barrel region
  h_PU200_nonprompt_pT_EB 	 = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffnoMTD_Bkg_EB");
  h_PU200_nonprompt_4sigma_pT_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_4sigma_Bkg_EB");
  h_PU200_nonprompt_3sigma_pT_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_3sigma_Bkg_EB");
  h_PU200_nonprompt_2sigma_pT_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_2sigma_Bkg_EB");
  h_PU200_nonprompt_40_pT_EB     = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_7_Bkg_EB");
  h_PU200_nonprompt_60_pT_EB     = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_5_Bkg_EB");
  h_PU200_nonprompt_80_pT_EB     = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_3_Bkg_EB");
  h_PU200_nonprompt_100_pT_EB    = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_1_Bkg_EB");
  h_PU200_nonprompt_gen_pT_EB	 = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_gen_Bkg_EB");
  h_noPU_nonprompt_pT_EB 	 = (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffnoMTD_Bkg_EB");         // noPU
  h_noPU_nonprompt_gen_pT_EB	 = (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffMTD_gen_Bkg_EB");       // noPU
    // Endcap region
  h_PU200_nonprompt_pT_EE 	 = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffnoMTD_Bkg_EE");
  h_PU200_nonprompt_4sigma_pT_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_4sigma_Bkg_EE");
  h_PU200_nonprompt_3sigma_pT_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_3sigma_Bkg_EE");
  h_PU200_nonprompt_2sigma_pT_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_2sigma_Bkg_EE");
  h_PU200_nonprompt_40_pT_EE     = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_7_Bkg_EE");
  h_PU200_nonprompt_60_pT_EE     = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_5_Bkg_EE");
  h_PU200_nonprompt_80_pT_EE     = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_3_Bkg_EE");
  h_PU200_nonprompt_100_pT_EE    = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_1_Bkg_EE");
  h_PU200_nonprompt_gen_pT_EE	 = (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_gen_Bkg_EE");
  h_noPU_nonprompt_pT_EE 	 = (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffnoMTD_Bkg_EE");         // noPU
  h_noPU_nonprompt_gen_pT_EE	 = (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffMTD_gen_Bkg_EE");       // noPU
 
  // Define the frame of histogram to change the range of x axis
  TH1D* h_PU200_prompt_incl_pT_EB = new TH1D("h_PU200_prompt_incl_pT_EB", "h_PU200_prompt_incl_pT_EB", 30, 1, 71);
  TH1D* h_PU200_prompt_incl_pT_EE = new TH1D("h_PU200_prompt_incl_pT_EE", "h_PU200_prompt_incl_pT_EE", 30, 1, 71);
  TH1D* h_PU200_nonprompt_incl_pT_EB = new TH1D("h_PU200_nonprompt_incl_pT_EB", "h_PU200_nonprompt_incl_pT_EB", 30, 1, 71);
  TH1D* h_PU200_nonprompt_incl_pT_EE = new TH1D("h_PU200_nonprompt_incl_pT_EE", "h_PU200_nonprompt_incl_pT_EE", 30, 1, 71);


  ///////////////
  // Cosmetics //
  ///////////////
  // Prompt
    // Barrel region
  h_PU200_prompt_incl_pT_EB->SetTitle("Prompt muon in barrel region");
  h_PU200_prompt_incl_pT_EB->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_prompt_incl_pT_EB->GetYaxis()->SetTitle("Efficiency");
  h_PU200_prompt_incl_pT_EB->GetYaxis()->SetRangeUser(0.6, 1.05);
  h_PU200_prompt_gen_pT_EB->SetLineColor(kBlack); h_PU200_prompt_pT_EB->SetLineColor(kBlack); h_PU200_prompt_40_pT_EB->SetLineColor(kRed); h_PU200_prompt_2sigma_pT_EB->SetLineColor(kRed); h_noPU_prompt_gen_pT_EB->SetLineColor(kGray); h_noPU_prompt_pT_EB->SetLineColor(kGray);
  h_PU200_prompt_gen_pT_EB->SetLineWidth(2); h_PU200_prompt_pT_EB->SetLineWidth(2); h_PU200_prompt_40_pT_EB->SetLineWidth(2); h_PU200_prompt_2sigma_pT_EB->SetLineWidth(2); h_noPU_prompt_gen_pT_EB->SetLineWidth(2); h_noPU_prompt_pT_EB->SetLineWidth(2);
  h_PU200_prompt_gen_pT_EB->SetLineStyle(2); h_noPU_prompt_gen_pT_EB->SetLineStyle(2);
    // Endcap region
  h_PU200_prompt_incl_pT_EE->SetTitle("Prompt muon in endcap region");
  h_PU200_prompt_incl_pT_EE->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_prompt_incl_pT_EE->GetYaxis()->SetTitle("Efficiency");
  h_PU200_prompt_incl_pT_EE->GetYaxis()->SetRangeUser(0.6, 1.05);
  h_PU200_prompt_gen_pT_EE->SetLineColor(kBlack); h_PU200_prompt_pT_EE->SetLineColor(kBlack); h_PU200_prompt_40_pT_EE->SetLineColor(kRed); h_PU200_prompt_2sigma_pT_EE->SetLineColor(kRed); h_noPU_prompt_gen_pT_EE->SetLineColor(kGray); h_noPU_prompt_pT_EE->SetLineColor(kGray);
  h_PU200_prompt_gen_pT_EE->SetLineWidth(2); h_PU200_prompt_pT_EE->SetLineWidth(2); h_PU200_prompt_40_pT_EE->SetLineWidth(2); h_PU200_prompt_2sigma_pT_EE->SetLineWidth(2); h_noPU_prompt_gen_pT_EE->SetLineWidth(2); h_noPU_prompt_pT_EE->SetLineWidth(2);
  h_PU200_prompt_gen_pT_EE->SetLineStyle(2); h_noPU_prompt_gen_pT_EE->SetLineStyle(2);

  // Nonprompt
    // Barrel region
  h_PU200_nonprompt_incl_pT_EB->SetTitle("Non-prompt muon in barrel region");
  h_PU200_nonprompt_incl_pT_EB->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_nonprompt_incl_pT_EB->GetYaxis()->SetTitle("Efficiency");
  h_PU200_nonprompt_incl_pT_EB->GetYaxis()->SetRangeUser(0., 1.05);
  h_PU200_nonprompt_gen_pT_EB->SetLineColor(kBlack); h_PU200_nonprompt_pT_EB->SetLineColor(kBlack); h_PU200_nonprompt_40_pT_EB->SetLineColor(kRed); h_PU200_nonprompt_2sigma_pT_EB->SetLineColor(kRed); h_noPU_nonprompt_gen_pT_EB->SetLineColor(kGray); h_noPU_nonprompt_pT_EB->SetLineColor(kGray);
  h_PU200_nonprompt_gen_pT_EB->SetLineWidth(2); h_PU200_nonprompt_pT_EB->SetLineWidth(2); h_PU200_nonprompt_40_pT_EB->SetLineWidth(2); h_PU200_nonprompt_2sigma_pT_EB->SetLineWidth(2); h_noPU_nonprompt_gen_pT_EB->SetLineWidth(2); h_noPU_nonprompt_pT_EB->SetLineWidth(2);
  h_PU200_nonprompt_gen_pT_EB->SetLineStyle(2); h_noPU_nonprompt_gen_pT_EB->SetLineStyle(2);
    // Endcap region
  h_PU200_nonprompt_incl_pT_EE->SetTitle("Non-prompt muon in endcap region");
  h_PU200_nonprompt_incl_pT_EE->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_nonprompt_incl_pT_EE->GetYaxis()->SetTitle("Efficiency");
  h_PU200_nonprompt_incl_pT_EE->GetYaxis()->SetRangeUser(0., 1.05);
  h_PU200_nonprompt_gen_pT_EE->SetLineColor(kBlack); h_PU200_nonprompt_pT_EE->SetLineColor(kBlack); h_PU200_nonprompt_40_pT_EE->SetLineColor(kRed); h_PU200_nonprompt_2sigma_pT_EE->SetLineColor(kRed); h_noPU_nonprompt_gen_pT_EE->SetLineColor(kGray); h_noPU_nonprompt_pT_EE->SetLineColor(kGray);
  h_PU200_nonprompt_gen_pT_EE->SetLineWidth(2); h_PU200_nonprompt_pT_EE->SetLineWidth(2); h_PU200_nonprompt_40_pT_EE->SetLineWidth(2); h_PU200_nonprompt_2sigma_pT_EE->SetLineWidth(2); h_noPU_nonprompt_gen_pT_EE->SetLineWidth(2); h_noPU_nonprompt_pT_EE->SetLineWidth(2);
  h_PU200_nonprompt_gen_pT_EE->SetLineStyle(2); h_noPU_nonprompt_gen_pT_EE->SetLineStyle(2);


  /////////////
  // Legends //
  /////////////
  // Prompt
    // Barrel region
  TLegend* leg_prompt_incl_pT_EB_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_prompt_incl_pT_EB_sigma->AddEntry(h_PU200_prompt_gen_pT_EB, "gen PU200");
  leg_prompt_incl_pT_EB_sigma->AddEntry(h_PU200_prompt_pT_EB, "no MTD PU200");
  leg_prompt_incl_pT_EB_sigma->AddEntry(h_PU200_prompt_2sigma_pT_EB, "2sigma PU200");
  leg_prompt_incl_pT_EB_sigma->AddEntry(h_noPU_prompt_gen_pT_EB, "gen noPU");
  leg_prompt_incl_pT_EB_sigma->AddEntry(h_noPU_prompt_pT_EB, "noPU");
  leg_prompt_incl_pT_EB_sigma->SetTextSize(0.03);
    // Endcap region
  TLegend* leg_prompt_incl_pT_EE_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_prompt_incl_pT_EE_sigma->AddEntry(h_PU200_prompt_gen_pT_EE, "gen PU200");
  leg_prompt_incl_pT_EE_sigma->AddEntry(h_PU200_prompt_pT_EE, "no MTD PU200");
  leg_prompt_incl_pT_EE_sigma->AddEntry(h_PU200_prompt_2sigma_pT_EE, "2sigma PU200");
  leg_prompt_incl_pT_EE_sigma->AddEntry(h_noPU_prompt_gen_pT_EE, "gen noPU");
  leg_prompt_incl_pT_EE_sigma->AddEntry(h_noPU_prompt_pT_EE, "noPU");
  leg_prompt_incl_pT_EE_sigma->SetTextSize(0.03);
  // Nonprompt
    // Barrel region
  TLegend* leg_nonprompt_incl_pT_EB_sigma = new TLegend(0.14, 0.62, 0.43, 0.88);
  leg_nonprompt_incl_pT_EB_sigma->AddEntry(h_PU200_nonprompt_gen_pT_EB, "gen PU200");
  leg_nonprompt_incl_pT_EB_sigma->AddEntry(h_PU200_nonprompt_pT_EB, "no MTD PU200");
  leg_nonprompt_incl_pT_EB_sigma->AddEntry(h_PU200_nonprompt_2sigma_pT_EB, "2sigma PU200");
  leg_nonprompt_incl_pT_EB_sigma->AddEntry(h_noPU_nonprompt_gen_pT_EB, "gen noPU");
  leg_nonprompt_incl_pT_EB_sigma->AddEntry(h_noPU_nonprompt_pT_EB, "noPU");
  leg_nonprompt_incl_pT_EB_sigma->SetTextSize(0.03);
    // Endcap region
  TLegend* leg_nonprompt_incl_pT_EE_sigma = new TLegend(0.14, 0.62, 0.43, 0.88);
  leg_nonprompt_incl_pT_EE_sigma->AddEntry(h_PU200_nonprompt_gen_pT_EE, "gen PU200");
  leg_nonprompt_incl_pT_EE_sigma->AddEntry(h_PU200_nonprompt_pT_EE, "no MTD PU200");
  leg_nonprompt_incl_pT_EE_sigma->AddEntry(h_PU200_nonprompt_2sigma_pT_EE, "2sigma PU200");
  leg_nonprompt_incl_pT_EE_sigma->AddEntry(h_noPU_nonprompt_gen_pT_EE, "gen noPU");
  leg_nonprompt_incl_pT_EE_sigma->AddEntry(h_noPU_nonprompt_pT_EE, "noPU");
  leg_nonprompt_incl_pT_EE_sigma->SetTextSize(0.03);


  ////////////////
  // Draw plots //
  ////////////////
  // Prompt
    // Barrel region
  TCanvas* c_PU200_prompt_pT_EB_sigma = new TCanvas("c_PU200_prompt_pT_EB_sigma", "c_PU200_prompt_pT_EB_sigma", 1500, 1500);
  c_PU200_prompt_pT_EB_sigma->cd();
  c_PU200_prompt_pT_EB_sigma->SetGrid();
  c_PU200_prompt_pT_EB_sigma->SetLeftMargin(0.12);
  h_PU200_prompt_incl_pT_EB->Draw("hist");
  h_PU200_prompt_gen_pT_EB->Draw("same e");
  h_PU200_prompt_pT_EB->Draw("same e");
  h_PU200_prompt_2sigma_pT_EB->Draw("same e");
  h_noPU_prompt_gen_pT_EB->Draw("same e");
  h_noPU_prompt_pT_EB->Draw("same e");
  leg_prompt_incl_pT_EB_sigma->Draw();
  c_PU200_prompt_pT_EB_sigma->Print("plots/pT_PU200_prompt_EB_sigma.pdf");
    // Endcap region
  TCanvas* c_PU200_prompt_pT_EE_sigma = new TCanvas("c_PU200_prompt_pT_EE_sigma", "c_PU200_prompt_pT_EE_sigma", 1500, 1500);
  c_PU200_prompt_pT_EE_sigma->cd();
  c_PU200_prompt_pT_EE_sigma->SetGrid();
  c_PU200_prompt_pT_EE_sigma->SetLeftMargin(0.12);
  h_PU200_prompt_incl_pT_EE->Draw("hist");
  h_PU200_prompt_gen_pT_EE->Draw("same e");
  h_PU200_prompt_pT_EE->Draw("same e");
  h_PU200_prompt_2sigma_pT_EE->Draw("same e");
  h_noPU_prompt_gen_pT_EE->Draw("same e");
  h_noPU_prompt_pT_EE->Draw("same e");
  leg_prompt_incl_pT_EE_sigma->Draw();
  c_PU200_prompt_pT_EE_sigma->Print("plots/pT_PU200_prompt_EE_sigma.pdf");
  // Nonprompt
    // Barrel region
  TCanvas* c_PU200_nonprompt_pT_EB_sigma = new TCanvas("c_PU200_nonprompt_pT_EB_sigma", "c_PU200_nonprompt_pT_EB_sigma", 1500, 1500);
  c_PU200_nonprompt_pT_EB_sigma->cd();
  c_PU200_nonprompt_pT_EB_sigma->SetGrid();
  c_PU200_nonprompt_pT_EB_sigma->SetLeftMargin(0.12);
  h_PU200_nonprompt_incl_pT_EB->Draw("hist");
  h_PU200_nonprompt_gen_pT_EB->Draw("same e");
  h_PU200_nonprompt_pT_EB->Draw("same e");
  h_PU200_nonprompt_2sigma_pT_EB->Draw("same e");
  h_noPU_nonprompt_gen_pT_EB->Draw("same e");
  h_noPU_nonprompt_pT_EB->Draw("same e");
  leg_nonprompt_incl_pT_EB_sigma->Draw();
  c_PU200_nonprompt_pT_EB_sigma->Print("plots/pT_PU200_nonprompt_EB_sigma.pdf");
    // Endcap region
  TCanvas* c_PU200_nonprompt_pT_EE_sigma = new TCanvas("c_PU200_nonprompt_pT_EE_sigma", "c_PU200_nonprompt_pT_EE_sigma", 1500, 1500);
  c_PU200_nonprompt_pT_EE_sigma->cd();
  c_PU200_nonprompt_pT_EE_sigma->SetGrid();
  c_PU200_nonprompt_pT_EE_sigma->SetLeftMargin(0.12);
  h_PU200_nonprompt_incl_pT_EE->Draw("hist");
  h_PU200_nonprompt_gen_pT_EE->Draw("same e");
  h_PU200_nonprompt_pT_EE->Draw("same e");
  h_PU200_nonprompt_2sigma_pT_EE->Draw("same e");
  h_noPU_nonprompt_gen_pT_EE->Draw("same e");
  h_noPU_nonprompt_pT_EE->Draw("same e");
  leg_nonprompt_incl_pT_EE_sigma->Draw();
  c_PU200_nonprompt_pT_EE_sigma->Print("plots/pT_PU200_nonprompt_EE_sigma.pdf");
}


void N_genMatched() {
  gStyle->SetOptStat(0);
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  // PU200
  f_PU200_prompt = new TFile("data/231203/harvester_PU200_prompt_160_20.root");
  f_PU200_nonprompt = new TFile("data/231203/harvester_PU200_nonprompt_90.root");
  // noPU
  f_noPU_prompt = new TFile("data/231203/harvester_noPU_prompt_19.root");
  f_noPU_nonprompt = new TFile("data/231203/harvester_noPU_nonprompt_9.root");


  TString dir="/DQMData/Run 1/MTD/Run summary/MuonIso/";
  TH1D *h_PU200_prompt, *h_PU200_nonprompt, *h_noPU_prompt, *h_noPU_nonprompt;

  h_PU200_prompt 	= (TH1D*)f_PU200_prompt->Get(dir+"Track_genMatch_info_check");
  h_PU200_nonprompt 	= (TH1D*)f_PU200_nonprompt->Get(dir+"Track_genMatch_info_check");
  h_noPU_prompt 	= (TH1D*)f_noPU_prompt->Get(dir+"Track_genMatch_info_check");
  h_noPU_nonprompt 	= (TH1D*)f_noPU_nonprompt->Get(dir+"Track_genMatch_info_check");

  
  ///////////////
  // Cosmetics //
  ///////////////
  // Barrel region
  h_PU200_prompt->SetTitle("Check whether a track is matched with a GenParticle"); h_PU200_nonprompt->SetTitle("Check whether a track is matched with a GenParticle"); h_noPU_prompt->SetTitle("Check whether a track is matched with a GenParticle"); h_noPU_nonprompt->SetTitle("Check whether a track is matched with a GenParticle");
  h_PU200_prompt->GetXaxis()->SetTitle(""); h_PU200_nonprompt->GetXaxis()->SetTitle(""); h_noPU_prompt->GetXaxis()->SetTitle(""); h_noPU_nonprompt->GetXaxis()->SetTitle("");
  h_PU200_prompt->GetYaxis()->SetTitle("Counts"); h_PU200_nonprompt->GetYaxis()->SetTitle("Counts"); h_noPU_prompt->GetYaxis()->SetTitle("Counts"); h_noPU_nonprompt->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt->SetLineWidth(2); h_PU200_nonprompt->SetLineWidth(2); h_noPU_prompt->SetLineWidth(2); h_noPU_nonprompt->SetLineWidth(2);
  h_PU200_prompt->SetLineColor(kBlack); h_PU200_nonprompt->SetLineColor(kBlack); h_noPU_prompt->SetLineColor(kBlack); h_noPU_nonprompt->SetLineColor(kBlack);
  h_PU200_prompt->SetNdivisions(3); h_PU200_nonprompt->SetNdivisions(3); h_noPU_prompt->SetNdivisions(3); h_noPU_nonprompt->SetNdivisions(3);
  h_PU200_prompt->SetMinimum(0.); h_PU200_nonprompt->SetMinimum(0.); h_noPU_prompt->SetMinimum(0.); h_noPU_nonprompt->SetMinimum(0.);


  ////////////////
  // Draw plots //
  ////////////////
  // PU200
  TCanvas* c_PU200_prompt = new TCanvas("c_PU200_prompt", "c_PU200_prompt", 1500, 1500);
  c_PU200_prompt->cd();
  h_PU200_prompt->Draw("hist");
  c_PU200_prompt->SetLeftMargin(0.12);
  c_PU200_prompt->Print("plots/N_genMatched_PU200_prompt.pdf");

  TCanvas* c_PU200_nonprompt = new TCanvas("c_PU200_nonprompt", "c_PU200_nonprompt", 1500, 1500);
  c_PU200_nonprompt->cd();
  h_PU200_nonprompt->Draw("hist");
  c_PU200_nonprompt->SetLeftMargin(0.12);
  c_PU200_nonprompt->Print("plots/N_genMatched_PU200_nonprompt.pdf");

  // noPU
  TCanvas* c_noPU_prompt = new TCanvas("c_noPU_prompt", "c_noPU_prompt", 1500, 1500);
  c_noPU_prompt->cd();
  h_noPU_prompt->Draw("hist");
  c_noPU_prompt->SetLeftMargin(0.12);
  c_noPU_prompt->Print("plots/N_genMatched_noPU_prompt.pdf");

  TCanvas* c_noPU_nonprompt = new TCanvas("c_noPU_nonprompt", "c_noPU_nonprompt", 1500, 1500);
  c_noPU_nonprompt->cd();
  h_noPU_nonprompt->Draw("hist");
  c_noPU_nonprompt->SetLeftMargin(0.12);
  c_noPU_nonprompt->Print("plots/N_genMatched_noPU_nonprompt.pdf");

}


int main(int argc, char **argv)
{
  bool zeroiso;
  zeroiso = true;

  draw_iso_efficiency(zeroiso);
  draw_reliso_roc(zeroiso);
  draw_pt_roc(zeroiso);
  pTeff();
  N_genMatched();


  return 0;

}
