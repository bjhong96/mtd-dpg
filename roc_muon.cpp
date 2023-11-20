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


using namespace std;

void roc() {
  //gStyle->SetOptStat(0);
  gROOT->ProcessLine("#include <vector>");

  // PU200
  TChain* ch_PU200_prompt    = new TChain("mtdMuonIsoValid/muonIso");
  TChain* ch_PU200_nonprompt = new TChain("mtdMuonIsoValid/muonIso");
  ch_PU200_prompt   ->Add("data/ntuple_PU200_Zall_prompt_180.root");
  ch_PU200_nonprompt->Add("data/ntuple_PU200_nonprompt_90.root");
  // noPU
  TChain* ch_noPU_prompt    = new TChain("mtdMuonIsoValid/muonIso");
  TChain* ch_noPU_nonprompt = new TChain("mtdMuonIsoValid/muonIso");
  ch_noPU_prompt    ->Add("data/ntuple_prompt_noPU_63.root");
  ch_noPU_nonprompt ->Add("data/ntuple_nonprompt_noPU_20.root");


  cout << "Entries of PU200_prompt: " << ch_PU200_prompt->GetEntries() << endl;
  cout << "Entries of PU200_nonprompt: " << ch_PU200_nonprompt->GetEntries() << endl;
  cout << "Entries of noPU_prompt: " << ch_noPU_prompt->GetEntries() << endl;
  cout << "Entries of noPU_nonprompt: " << ch_noPU_nonprompt->GetEntries() << endl;

  int Nbin = 1000;
  float bin_i = 0.01;
  float bin_f = 0.5;

  // PU200
    // no MTD case
  TH1D* h_prompt_EB = new TH1D("h_prompt_EB", "h_prompt_EB", Nbin, bin_i, bin_f);
  TH1D* h_prompt_EE = new TH1D("h_prompt_EE", "h_prompt_EE", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_EB = new TH1D("h_nonprompt_EB", "h_nonprompt_EB", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_EE = new TH1D("h_nonprompt_EE", "h_nonprompt_EE", Nbin, bin_i, bin_f);
    // MTD case
      // prompt
  TH1D* h_prompt_EB_40 = new TH1D("h_prompt_EB_40", "h_prompt_EB_40", Nbin, bin_i, bin_f);
  TH1D* h_prompt_EE_40 = new TH1D("h_prompt_EE_40", "h_prompt_EE_40", Nbin, bin_i, bin_f);
  TH1D* h_prompt_EB_60 = new TH1D("h_prompt_EB_60", "h_prompt_EB_60", Nbin, bin_i, bin_f);
  TH1D* h_prompt_EE_60 = new TH1D("h_prompt_EE_60", "h_prompt_EE_60", Nbin, bin_i, bin_f);
  TH1D* h_prompt_EB_80 = new TH1D("h_prompt_EB_80", "h_prompt_EB_80", Nbin, bin_i, bin_f);
  TH1D* h_prompt_EE_80 = new TH1D("h_prompt_EE_80", "h_prompt_EE_80", Nbin, bin_i, bin_f);
  TH1D* h_prompt_EB_100 = new TH1D("h_prompt_EB_100", "h_prompt_EB_100", Nbin, bin_i, bin_f);
  TH1D* h_prompt_EE_100 = new TH1D("h_prompt_EE_100", "h_prompt_EE_100", Nbin, bin_i, bin_f);
      // nonprompt
  TH1D* h_nonprompt_EB_40 = new TH1D("h_nonprompt_EB_40", "h_nonprompt_EB_40", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_EE_40 = new TH1D("h_nonprompt_EE_40", "h_nonprompt_EE_40", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_EB_60 = new TH1D("h_nonprompt_EB_60", "h_nonprompt_EB_60", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_EE_60 = new TH1D("h_nonprompt_EE_60", "h_nonprompt_EE_60", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_EB_80 = new TH1D("h_nonprompt_EB_80", "h_nonprompt_EB_80", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_EE_80 = new TH1D("h_nonprompt_EE_80", "h_nonprompt_EE_80", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_EB_100 = new TH1D("h_nonprompt_EB_100", "h_nonprompt_EB_100", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_EE_100 = new TH1D("h_nonprompt_EE_100", "h_nonprompt_EE_100", Nbin, bin_i, bin_f);

  // noPU
    // no MTD case
  TH1D* h_prompt_noPU_EB = new TH1D("h_prompt_noPU_EB", "h_prompt_noPU_EB", Nbin, bin_i, bin_f);
  TH1D* h_prompt_noPU_EE = new TH1D("h_prompt_noPU_EE", "h_prompt_noPU_EE", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_noPU_EB = new TH1D("h_nonprompt_noPU_EB", "h_nonprompt_noPU_EB", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_noPU_EE = new TH1D("h_nonprompt_noPU_EE", "h_nonprompt_noPU_EE", Nbin, bin_i, bin_f);
    // MTD case
      // prompt
  TH1D* h_prompt_noPU_EB_40 = new TH1D("h_prompt_noPU_EB_40", "h_prompt_noPU_EB_40", Nbin, bin_i, bin_f);
  TH1D* h_prompt_noPU_EE_40 = new TH1D("h_prompt_noPU_EE_40", "h_prompt_noPU_EE_40", Nbin, bin_i, bin_f);
  TH1D* h_prompt_noPU_EB_60 = new TH1D("h_prompt_noPU_EB_60", "h_prompt_noPU_EB_60", Nbin, bin_i, bin_f);
  TH1D* h_prompt_noPU_EE_60 = new TH1D("h_prompt_noPU_EE_60", "h_prompt_noPU_EE_60", Nbin, bin_i, bin_f);
  TH1D* h_prompt_noPU_EB_80 = new TH1D("h_prompt_noPU_EB_80", "h_prompt_noPU_EB_80", Nbin, bin_i, bin_f);
  TH1D* h_prompt_noPU_EE_80 = new TH1D("h_prompt_noPU_EE_80", "h_prompt_noPU_EE_80", Nbin, bin_i, bin_f);
  TH1D* h_prompt_noPU_EB_100 = new TH1D("h_prompt_noPU_EB_100", "h_prompt_noPU_EB_100", Nbin, bin_i, bin_f);
  TH1D* h_prompt_noPU_EE_100 = new TH1D("h_prompt_noPU_EE_100", "h_prompt_noPU_EE_100", Nbin, bin_i, bin_f);
      // nonprompt
  TH1D* h_nonprompt_noPU_EB_40 = new TH1D("h_nonprompt_noPU_EB_40", "h_nonprompt_noPU_EB_40", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_noPU_EE_40 = new TH1D("h_nonprompt_noPU_EE_40", "h_nonprompt_noPU_EE_40", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_noPU_EB_60 = new TH1D("h_nonprompt_noPU_EB_60", "h_nonprompt_noPU_EB_60", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_noPU_EE_60 = new TH1D("h_nonprompt_noPU_EE_60", "h_nonprompt_noPU_EE_60", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_noPU_EB_80 = new TH1D("h_nonprompt_noPU_EB_80", "h_nonprompt_noPU_EB_80", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_noPU_EE_80 = new TH1D("h_nonprompt_noPU_EE_80", "h_nonprompt_noPU_EE_80", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_noPU_EB_100 = new TH1D("h_nonprompt_noPU_EB_100", "h_nonprompt_noPU_EB_100", Nbin, bin_i, bin_f);
  TH1D* h_nonprompt_noPU_EE_100 = new TH1D("h_nonprompt_noPU_EE_100", "h_nonprompt_noPU_EE_100", Nbin, bin_i, bin_f);


  // PU200
    // no MTD case
    // To get rid of the event where track doesn't exist, the selection of Sum$(track_pt)>0 is used
  ch_PU200_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_EB", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)<1.5) && idx_muon==idx_track_muon_match", "goff");
  ch_PU200_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_EE", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)>1.5) && idx_muon==idx_track_muon_match", "goff");
  ch_PU200_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_EB", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)<1.5) && idx_muon==idx_track_muon_match", "goff");
  ch_PU200_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_EE", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)>1.5) && idx_muon==idx_track_muon_match", "goff");
    // MTD case
      // prompt
  ch_PU200_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_EB_40", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)<1.5) && dt_track_muon<0.12 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_EE_40", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)>1.5) && dt_track_muon<0.12 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_EB_60", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)<1.5) && dt_track_muon<0.18 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_EE_60", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)>1.5) && dt_track_muon<0.18 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_EB_80", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)<1.5) && dt_track_muon<0.24 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_EE_80", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)>1.5) && dt_track_muon<0.24 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_EB_100", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)<1.5) && dt_track_muon<0.3 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_EE_100", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)>1.5) && dt_track_muon<0.3 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
      // nonprompt
  ch_PU200_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_EB_40", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)<1.5) && dt_track_muon<0.12 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_EE_40", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)>1.5) && dt_track_muon<0.12 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_EB_60", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)<1.5) && dt_track_muon<0.18 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_EE_60", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)>1.5) && dt_track_muon<0.18 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_EB_80", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)<1.5) && dt_track_muon<0.24 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_EE_80", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)>1.5) && dt_track_muon<0.24 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_EB_100", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)<1.5) && dt_track_muon<0.3 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_PU200_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_EE_100", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)>1.5) && dt_track_muon<0.3 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");

  // noPU
    // no MTD case
  ch_noPU_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_noPU_EB", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)<1.5) && idx_muon==idx_track_muon_match", "goff");
  ch_noPU_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_noPU_EE", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)>1.5) && idx_muon==idx_track_muon_match", "goff");
  ch_noPU_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_noPU_EB", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)<1.5) && idx_muon==idx_track_muon_match", "goff");
  ch_noPU_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_noPU_EE", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)>1.5) && idx_muon==idx_track_muon_match", "goff");
    // MTD case
      // prompt
  ch_noPU_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_noPU_EB_40", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)<1.5) && dt_track_muon<0.12 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_noPU_EE_40", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)>1.5) && dt_track_muon<0.12 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_noPU_EB_60", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)<1.5) && dt_track_muon<0.18 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_noPU_EE_60", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)>1.5) && dt_track_muon<0.18 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_noPU_EB_80", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)<1.5) && dt_track_muon<0.24 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_noPU_EE_80", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)>1.5) && dt_track_muon<0.24 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_noPU_EB_100", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)<1.5) && dt_track_muon<0.30 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_prompt->Draw("(Sum$(track_pt)/muon_pt)>>h_prompt_noPU_EE_100", "Sum$(track_pt)>0 && muon_prompt==1 && (fabs(muon_eta)>1.5) && dt_track_muon<0.30 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
      // nonprompt
  ch_noPU_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_noPU_EB_40", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)<1.5) && dt_track_muon<0.12 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_noPU_EE_40", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)>1.5) && dt_track_muon<0.12 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_noPU_EB_60", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)<1.5) && dt_track_muon<0.18 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_noPU_EE_60", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)>1.5) && dt_track_muon<0.18 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_noPU_EB_80", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)<1.5) && dt_track_muon<0.24 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_noPU_EE_80", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)>1.5) && dt_track_muon<0.24 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_noPU_EB_100", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)<1.5) && dt_track_muon<0.30 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");
  ch_noPU_nonprompt->Draw("(Sum$(track_pt)/muon_pt)>>h_nonprompt_noPU_EE_100", "Sum$(track_pt)>0 && muon_prompt==0 && (fabs(muon_eta)>1.5) && dt_track_muon<0.30 && idx_muon==idx_track_muon_match && muon_time_err>0", "goff");


  // PU200
    // no MTD case
  vector<double> prompt_eff_EB={0}, prompt_eff_EE={0};
  vector<double> nonprompt_eff_EB={0}, nonprompt_eff_EE={0};
    // MTD case
      // prompt
  vector<double> prompt_eff_EB_40={0}, prompt_eff_EE_40={0};
  vector<double> prompt_eff_EB_60={0}, prompt_eff_EE_60={0};
  vector<double> prompt_eff_EB_80={0}, prompt_eff_EE_80={0};
  vector<double> prompt_eff_EB_100={0}, prompt_eff_EE_100={0};
      // nonprompt
  vector<double> nonprompt_eff_EB_40={0}, nonprompt_eff_EE_40={0};
  vector<double> nonprompt_eff_EB_60={0}, nonprompt_eff_EE_60={0};
  vector<double> nonprompt_eff_EB_80={0}, nonprompt_eff_EE_80={0};
  vector<double> nonprompt_eff_EB_100={0}, nonprompt_eff_EE_100={0};
  // noPU
    // no MTD case
  vector<double> prompt_eff_noPU_EB={0}, prompt_eff_noPU_EE={0};
  vector<double> nonprompt_eff_noPU_EB={0}, nonprompt_eff_noPU_EE={0};
    // MTD case
      // prompt
  vector<double> prompt_eff_noPU_EB_40={0}, prompt_eff_noPU_EE_40={0};
  vector<double> prompt_eff_noPU_EB_60={0}, prompt_eff_noPU_EE_60={0};
  vector<double> prompt_eff_noPU_EB_80={0}, prompt_eff_noPU_EE_80={0};
  vector<double> prompt_eff_noPU_EB_100={0}, prompt_eff_noPU_EE_100={0};
      // nonprompt
  vector<double> nonprompt_eff_noPU_EB_40={0}, nonprompt_eff_noPU_EE_40={0};
  vector<double> nonprompt_eff_noPU_EB_60={0}, nonprompt_eff_noPU_EE_60={0};
  vector<double> nonprompt_eff_noPU_EB_80={0}, nonprompt_eff_noPU_EE_80={0};
  vector<double> nonprompt_eff_noPU_EB_100={0}, nonprompt_eff_noPU_EE_100={0};

  for(int i=0; i<h_prompt_EB->GetNbinsX(); i++) {

    // PU200
      // no MTD case
      // To consider underflow of the histogram, the range starts at 0
      // To consider overflow of the histogram, the range of denominator is from 0 to -1
    prompt_eff_EB.emplace_back(h_prompt_EB->Integral(0,i+1)/h_prompt_EB->Integral(0,-1));
    prompt_eff_EE.emplace_back(h_prompt_EE->Integral(0,i+1)/h_prompt_EE->Integral(0,-1));
    nonprompt_eff_EB.emplace_back(h_nonprompt_EB->Integral(0,i+1)/h_nonprompt_EB->Integral(0,-1));
    nonprompt_eff_EE.emplace_back(h_nonprompt_EE->Integral(0,i+1)/h_nonprompt_EE->Integral(0,-1));
      // MTD case
        // prompt
    prompt_eff_EB_40.emplace_back(h_prompt_EB_40->Integral(0,i+1)/h_prompt_EB_40->Integral(0,-1));
    prompt_eff_EE_40.emplace_back(h_prompt_EE_40->Integral(0,i+1)/h_prompt_EE_40->Integral(0,-1));
    prompt_eff_EB_60.emplace_back(h_prompt_EB_60->Integral(0,i+1)/h_prompt_EB_60->Integral(0,-1));
    prompt_eff_EE_60.emplace_back(h_prompt_EE_60->Integral(0,i+1)/h_prompt_EE_60->Integral(0,-1));
    prompt_eff_EB_80.emplace_back(h_prompt_EB_80->Integral(0,i+1)/h_prompt_EB_80->Integral(0,-1));
    prompt_eff_EE_80.emplace_back(h_prompt_EE_80->Integral(0,i+1)/h_prompt_EE_80->Integral(0,-1));
    prompt_eff_EB_100.emplace_back(h_prompt_EB_100->Integral(0,i+1)/h_prompt_EB_100->Integral(0,-1));
    prompt_eff_EE_100.emplace_back(h_prompt_EE_100->Integral(0,i+1)/h_prompt_EE_100->Integral(0,-1));
        // nonprompt
    nonprompt_eff_EB_40.emplace_back(h_nonprompt_EB_40->Integral(0,i+1)/h_nonprompt_EB_40->Integral(0,-1));
    nonprompt_eff_EE_40.emplace_back(h_nonprompt_EE_40->Integral(0,i+1)/h_nonprompt_EE_40->Integral(0,-1));
    nonprompt_eff_EB_60.emplace_back(h_nonprompt_EB_60->Integral(0,i+1)/h_nonprompt_EB_60->Integral(0,-1));
    nonprompt_eff_EE_60.emplace_back(h_nonprompt_EE_60->Integral(0,i+1)/h_nonprompt_EE_60->Integral(0,-1));
    nonprompt_eff_EB_80.emplace_back(h_nonprompt_EB_80->Integral(0,i+1)/h_nonprompt_EB_80->Integral(0,-1));
    nonprompt_eff_EE_80.emplace_back(h_nonprompt_EE_80->Integral(0,i+1)/h_nonprompt_EE_80->Integral(0,-1));
    nonprompt_eff_EB_100.emplace_back(h_nonprompt_EB_100->Integral(0,i+1)/h_nonprompt_EB_100->Integral(0,-1));
    nonprompt_eff_EE_100.emplace_back(h_nonprompt_EE_100->Integral(0,i+1)/h_nonprompt_EE_100->Integral(0,-1));

    // noPU
      // no MTD case
    prompt_eff_noPU_EB.emplace_back(h_prompt_noPU_EB->Integral(0,i+1)/h_prompt_noPU_EB->Integral(0,-1));
    prompt_eff_noPU_EE.emplace_back(h_prompt_noPU_EE->Integral(0,i+1)/h_prompt_noPU_EE->Integral(0,-1));
    nonprompt_eff_noPU_EB.emplace_back(h_nonprompt_noPU_EB->Integral(0,i+1)/h_nonprompt_noPU_EB->Integral(0,-1));
    nonprompt_eff_noPU_EE.emplace_back(h_nonprompt_noPU_EE->Integral(0,i+1)/h_nonprompt_noPU_EE->Integral(0,-1));
      // MTD case
        // prompt
    prompt_eff_noPU_EB_40.emplace_back(h_prompt_noPU_EB_40->Integral(0,i+1)/h_prompt_noPU_EB_40->Integral(0,-1));
    prompt_eff_noPU_EE_40.emplace_back(h_prompt_noPU_EE_40->Integral(0,i+1)/h_prompt_noPU_EE_40->Integral(0,-1));
    prompt_eff_noPU_EB_60.emplace_back(h_prompt_noPU_EB_60->Integral(0,i+1)/h_prompt_noPU_EB_60->Integral(0,-1));
    prompt_eff_noPU_EE_60.emplace_back(h_prompt_noPU_EE_60->Integral(0,i+1)/h_prompt_noPU_EE_60->Integral(0,-1));
    prompt_eff_noPU_EB_80.emplace_back(h_prompt_noPU_EB_80->Integral(0,i+1)/h_prompt_noPU_EB_80->Integral(0,-1));
    prompt_eff_noPU_EE_80.emplace_back(h_prompt_noPU_EE_80->Integral(0,i+1)/h_prompt_noPU_EE_80->Integral(0,-1));
    prompt_eff_noPU_EB_100.emplace_back(h_prompt_noPU_EB_100->Integral(0,i+1)/h_prompt_noPU_EB_100->Integral(0,-1));
    prompt_eff_noPU_EE_100.emplace_back(h_prompt_noPU_EE_100->Integral(0,i+1)/h_prompt_noPU_EE_100->Integral(0,-1));
        // nonprompt
    nonprompt_eff_noPU_EB_40.emplace_back(h_nonprompt_noPU_EB_40->Integral(0,i+1)/h_nonprompt_noPU_EB_40->Integral(0,-1));
    nonprompt_eff_noPU_EE_40.emplace_back(h_nonprompt_noPU_EE_40->Integral(0,i+1)/h_nonprompt_noPU_EE_40->Integral(0,-1));
    nonprompt_eff_noPU_EB_60.emplace_back(h_nonprompt_noPU_EB_60->Integral(0,i+1)/h_nonprompt_noPU_EB_60->Integral(0,-1));
    nonprompt_eff_noPU_EE_60.emplace_back(h_nonprompt_noPU_EE_60->Integral(0,i+1)/h_nonprompt_noPU_EE_60->Integral(0,-1));
    nonprompt_eff_noPU_EB_80.emplace_back(h_nonprompt_noPU_EB_80->Integral(0,i+1)/h_nonprompt_noPU_EB_80->Integral(0,-1));
    nonprompt_eff_noPU_EE_80.emplace_back(h_nonprompt_noPU_EE_80->Integral(0,i+1)/h_nonprompt_noPU_EE_80->Integral(0,-1));
    nonprompt_eff_noPU_EB_100.emplace_back(h_nonprompt_noPU_EB_100->Integral(0,i+1)/h_nonprompt_noPU_EB_100->Integral(0,-1));
    nonprompt_eff_noPU_EE_100.emplace_back(h_nonprompt_noPU_EE_100->Integral(0,i+1)/h_nonprompt_noPU_EE_100->Integral(0,-1));
  }


  // PU200
    // no MTD Case
  TGraph* gr_PU200_EB = new TGraph(h_prompt_EB->GetNbinsX(), &prompt_eff_EB[0], &nonprompt_eff_EB[0]);
  TGraph* gr_PU200_EE = new TGraph(h_prompt_EE->GetNbinsX(), &prompt_eff_EE[0], &nonprompt_eff_EE[0]);
    // MTD case
  TGraph* gr_PU200_EB_40 = new TGraph(h_prompt_EB_40->GetNbinsX(), &prompt_eff_EB_40[0], &nonprompt_eff_EB_40[0]);
  TGraph* gr_PU200_EE_40 = new TGraph(h_prompt_EE_40->GetNbinsX(), &prompt_eff_EE_40[0], &nonprompt_eff_EE_40[0]);
  TGraph* gr_PU200_EB_60 = new TGraph(h_prompt_EB_60->GetNbinsX(), &prompt_eff_EB_60[0], &nonprompt_eff_EB_60[0]);
  TGraph* gr_PU200_EE_60 = new TGraph(h_prompt_EE_60->GetNbinsX(), &prompt_eff_EE_60[0], &nonprompt_eff_EE_60[0]);
  TGraph* gr_PU200_EB_80 = new TGraph(h_prompt_EB_80->GetNbinsX(), &prompt_eff_EB_80[0], &nonprompt_eff_EB_80[0]);
  TGraph* gr_PU200_EE_80 = new TGraph(h_prompt_EE_80->GetNbinsX(), &prompt_eff_EE_80[0], &nonprompt_eff_EE_80[0]);
  TGraph* gr_PU200_EB_100 = new TGraph(h_prompt_EB_100->GetNbinsX(), &prompt_eff_EB_100[0], &nonprompt_eff_EB_100[0]);
  TGraph* gr_PU200_EE_100 = new TGraph(h_prompt_EE_100->GetNbinsX(), &prompt_eff_EE_100[0], &nonprompt_eff_EE_100[0]);

  // noPU
  TGraph* gr_noPU_EB = new TGraph(h_prompt_noPU_EB->GetNbinsX(), &prompt_eff_noPU_EB[0], &nonprompt_eff_noPU_EB[0]);
  TGraph* gr_noPU_EE = new TGraph(h_prompt_noPU_EE->GetNbinsX(), &prompt_eff_noPU_EE[0], &nonprompt_eff_noPU_EE[0]);
    // MTD case
  TGraph* gr_noPU_EB_40 = new TGraph(h_prompt_noPU_EB_40->GetNbinsX(), &prompt_eff_noPU_EB_40[0], &nonprompt_eff_noPU_EB_40[0]);
  TGraph* gr_noPU_EE_40 = new TGraph(h_prompt_noPU_EE_40->GetNbinsX(), &prompt_eff_noPU_EE_40[0], &nonprompt_eff_noPU_EE_40[0]);
  TGraph* gr_noPU_EB_60 = new TGraph(h_prompt_noPU_EB_60->GetNbinsX(), &prompt_eff_noPU_EB_60[0], &nonprompt_eff_noPU_EB_60[0]);
  TGraph* gr_noPU_EE_60 = new TGraph(h_prompt_noPU_EE_60->GetNbinsX(), &prompt_eff_noPU_EE_60[0], &nonprompt_eff_noPU_EE_60[0]);
  TGraph* gr_noPU_EB_80 = new TGraph(h_prompt_noPU_EB_80->GetNbinsX(), &prompt_eff_noPU_EB_80[0], &nonprompt_eff_noPU_EB_80[0]);
  TGraph* gr_noPU_EE_80 = new TGraph(h_prompt_noPU_EE_80->GetNbinsX(), &prompt_eff_noPU_EE_80[0], &nonprompt_eff_noPU_EE_80[0]);
  TGraph* gr_noPU_EB_100 = new TGraph(h_prompt_noPU_EB_100->GetNbinsX(), &prompt_eff_noPU_EB_100[0], &nonprompt_eff_noPU_EB_100[0]);
  TGraph* gr_noPU_EE_100 = new TGraph(h_prompt_noPU_EE_100->GetNbinsX(), &prompt_eff_noPU_EE_100[0], &nonprompt_eff_noPU_EE_100[0]);


  // Cosmetics
    // Barrel
  gr_PU200_EB->SetTitle(""); gr_noPU_EB->SetTitle("");
  gr_PU200_EB->GetXaxis()->SetTitle("Prompt efficiency"); gr_noPU_EB->GetXaxis()->SetTitle("Prompt efficiency");
  gr_PU200_EB->GetYaxis()->SetTitle("Non-prompt efficiency"); gr_noPU_EB->GetYaxis()->SetTitle("Non-prompt efficiency");
  gr_PU200_EB->SetLineWidth(2); gr_PU200_EB_40->SetLineWidth(2); gr_PU200_EB_60->SetLineWidth(2); gr_PU200_EB_80->SetLineWidth(2); gr_PU200_EB_100->SetLineWidth(2);
  gr_noPU_EB->SetLineWidth(2); gr_noPU_EB_40->SetLineWidth(2); gr_noPU_EB_60->SetLineWidth(2); gr_noPU_EB_80->SetLineWidth(2); gr_noPU_EB_100->SetLineWidth(2);
  gr_PU200_EB->SetLineColor(kBlack); gr_PU200_EB_40->SetLineColor(kRed); gr_PU200_EB_60->SetLineColor(kGreen); gr_PU200_EB_80->SetLineColor(kBlue); gr_PU200_EB_100->SetLineColor(kMagenta);
  gr_noPU_EB->SetLineColor(kGray); gr_noPU_EB_40->SetLineColor(kRed); gr_noPU_EB_60->SetLineColor(kGreen); gr_noPU_EB_80->SetLineColor(kBlue); gr_noPU_EB_100->SetLineColor(kMagenta);
    // Endcap
  gr_PU200_EE->SetTitle(""); gr_noPU_EE->SetTitle("");
  gr_PU200_EE->GetXaxis()->SetTitle("Prompt efficiency"); gr_noPU_EE->GetXaxis()->SetTitle("Prompt efficiency");
  gr_PU200_EE->GetYaxis()->SetTitle("Non-prompt efficiency"); gr_noPU_EE->GetYaxis()->SetTitle("Non-prompt efficiency");
  gr_PU200_EE->SetLineWidth(2); gr_PU200_EE_40->SetLineWidth(2); gr_PU200_EE_60->SetLineWidth(2); gr_PU200_EE_80->SetLineWidth(2); gr_PU200_EE_100->SetLineWidth(2);
  gr_noPU_EE->SetLineWidth(2); gr_noPU_EE_40->SetLineWidth(2); gr_noPU_EE_60->SetLineWidth(2); gr_noPU_EE_80->SetLineWidth(2); gr_noPU_EE_100->SetLineWidth(2);
  gr_PU200_EE->SetLineColor(kBlack); gr_PU200_EE_40->SetLineColor(kRed); gr_PU200_EE_60->SetLineColor(kGreen); gr_PU200_EE_80->SetLineColor(kBlue); gr_PU200_EE_100->SetLineColor(kMagenta);
  gr_noPU_EE->SetLineColor(kGray); gr_noPU_EE_40->SetLineColor(kRed); gr_noPU_EE_60->SetLineColor(kGreen); gr_noPU_EE_80->SetLineColor(kBlue); gr_noPU_EE_100->SetLineColor(kMagenta);

  // Legend
    // Barrel
  TLegend* leg_1_EB = new TLegend(0.15, 0.65, 0.48, 0.85);
  leg_1_EB->AddEntry(gr_PU200_EB, "no MTD PU200");
  leg_1_EB->AddEntry(gr_PU200_EB_40, "40 ps PU200");
  leg_1_EB->AddEntry(gr_noPU_EB, "no MTD noPU");
  leg_1_EB->SetTextSize(0.03);
  TLegend* leg_2_EB = new TLegend(0.15, 0.60, 0.48, 0.85);
  leg_2_EB->AddEntry(gr_PU200_EB, "no MTD PU200");
  leg_2_EB->AddEntry(gr_PU200_EB_40, "40 ps PU200");
  leg_2_EB->AddEntry(gr_PU200_EB_60, "60 ps PU200");
  leg_2_EB->AddEntry(gr_PU200_EB_80, "80 ps PU200");
  leg_2_EB->AddEntry(gr_PU200_EB_100, "100 ps PU200");
  leg_2_EB->AddEntry(gr_noPU_EB, "no MTD noPU");
  leg_2_EB->SetTextSize(0.03);
  TLegend* leg_3_EB = new TLegend(0.15, 0.60, 0.48, 0.85);
  leg_3_EB->AddEntry(gr_noPU_EB, "no MTD no PU");
  leg_3_EB->AddEntry(gr_noPU_EB_40, "40 ps no PU");
  leg_3_EB->AddEntry(gr_noPU_EB_60, "60 ps no PU");
  leg_3_EB->AddEntry(gr_noPU_EB_80, "80 ps no PU");
  leg_3_EB->AddEntry(gr_noPU_EB_100, "100 ps no PU");
  leg_3_EB->SetTextSize(0.03);
    // Endcap
  TLegend* leg_1_EE = new TLegend(0.15, 0.65, 0.48, 0.85);
  leg_1_EE->AddEntry(gr_PU200_EE, "no MTD PU200");
  leg_1_EE->AddEntry(gr_PU200_EE_40, "40 ps PU200");
  leg_1_EE->AddEntry(gr_noPU_EE, "no MTD noPU");
  leg_1_EE->SetTextSize(0.03);
  TLegend* leg_2_EE = new TLegend(0.15, 0.60, 0.48, 0.85);
  leg_2_EE->AddEntry(gr_PU200_EE, "no MTD PU200");
  leg_2_EE->AddEntry(gr_PU200_EE_40, "40 ps PU200");
  leg_2_EE->AddEntry(gr_PU200_EE_60, "60 ps PU200");
  leg_2_EE->AddEntry(gr_PU200_EE_80, "80 ps PU200");
  leg_2_EE->AddEntry(gr_PU200_EE_100, "100 ps PU200");
  leg_2_EE->AddEntry(gr_noPU_EE, "no MTD noPU");
  leg_2_EE->SetTextSize(0.03);
  TLegend* leg_3_EE = new TLegend(0.15, 0.60, 0.48, 0.85);
  leg_3_EE->AddEntry(gr_noPU_EE, "no MTD no PU");
  leg_3_EE->AddEntry(gr_noPU_EE_40, "40 ps no PU");
  leg_3_EE->AddEntry(gr_noPU_EE_60, "60 ps no PU");
  leg_3_EE->AddEntry(gr_noPU_EE_80, "80 ps no PU");
  leg_3_EE->AddEntry(gr_noPU_EE_100, "100 ps no PU");
  leg_3_EE->SetTextSize(0.03);



  // Barrel
  TCanvas* c_noMTD_vs_40ps_EB = new TCanvas("c_noMTD_vs_40ps_EB","c_noMTD_vs_40ps_EB",1500,1500);
  c_noMTD_vs_40ps_EB->cd();
  c_noMTD_vs_40ps_EB->SetGrid();
  gr_PU200_EB->Draw("AL");
  gr_PU200_EB_40->Draw("same");
  gr_noPU_EB->Draw("same");
  leg_1_EB->Draw();
  c_noMTD_vs_40ps_EB->Print("roc_reliso_noMTD_vs_40ps_EB.pdf");

  TCanvas* c_noMTD_vs_totps_EB = new TCanvas("c_noMTD_vs_totps_EB","c_noMTD_vs_totps_EB",1500,1500);
  c_noMTD_vs_totps_EB->cd();
  c_noMTD_vs_totps_EB->SetGrid();
  gr_PU200_EB->Draw("AL");
  gr_PU200_EB_40->Draw("same");
  gr_PU200_EB_60->Draw("same");
  gr_PU200_EB_80->Draw("same");
  gr_PU200_EB_100->Draw("same");
  gr_noPU_EB->Draw("same");
  leg_2_EB->Draw();
  c_noMTD_vs_totps_EB->Print("roc_reliso_noMTD_vs_totps_EB.pdf");

  TCanvas* c_noMTD_vs_noMTD_EB = new TCanvas("c_noMTD_vs_noMTD_EB","c_noMTD_vs_noMTD_EB",1500,1500);
  c_noMTD_vs_noMTD_EB->cd();
  c_noMTD_vs_noMTD_EB->SetGrid();
  gr_noPU_EB->Draw("AL");
  gr_noPU_EB_40->Draw("same");
  gr_noPU_EB_60->Draw("same");
  gr_noPU_EB_80->Draw("same");
  gr_noPU_EB_100->Draw("same");
  leg_3_EB->Draw();
  c_noMTD_vs_noMTD_EB->Print("roc_reliso_noMTD_vs_noMTD_EB.pdf");

  TCanvas* c_noMTD_vs_40ps_EB_zoomed = new TCanvas("c_noMTD_vs_40ps_EB_zoomed","c_noMTD_vs_40ps_EB_zoomed",1500,1500);
  c_noMTD_vs_40ps_EB_zoomed->cd();
  c_noMTD_vs_40ps_EB_zoomed->SetGrid();
  gr_PU200_EB->GetXaxis()->SetRangeUser(0.8,1.);
  gr_PU200_EB->Draw("AL");
  gr_PU200_EB_40->Draw("same");
  gr_noPU_EB->Draw("same");
  leg_1_EB->Draw();
  c_noMTD_vs_40ps_EB_zoomed->Print("roc_reliso_noMTD_vs_40ps_EB_zoomed.pdf");

  TCanvas* c_noMTD_vs_totps_EB_zoomed = new TCanvas("c_noMTD_vs_totps_EB_zoomed","c_noMTD_vs_totps_EB_zoomed",1500,1500);
  c_noMTD_vs_totps_EB_zoomed->cd();
  c_noMTD_vs_totps_EB_zoomed->SetGrid();
  gr_PU200_EB->GetXaxis()->SetRangeUser(0.8,1.);
  gr_PU200_EB->Draw("AL");
  gr_PU200_EB_40->Draw("same");
  gr_PU200_EB_60->Draw("same");
  gr_PU200_EB_80->Draw("same");
  gr_PU200_EB_100->Draw("same");
  gr_noPU_EB->Draw("same");
  leg_2_EB->Draw();
  c_noMTD_vs_totps_EB_zoomed->Print("roc_reliso_noMTD_vs_totps_EB_zoomed.pdf");

  TCanvas* c_noMTD_vs_noMTD_EB_zoomed = new TCanvas("c_noMTD_vs_noMTD_EB_zoomed","c_noMTD_vs_noMTD_EB_zoomed",1500,1500);
  c_noMTD_vs_noMTD_EB_zoomed->cd();
  c_noMTD_vs_noMTD_EB_zoomed->SetGrid();
  gr_noPU_EB->GetXaxis()->SetRangeUser(0.8,1.);
  gr_noPU_EB->Draw("AL");
  gr_noPU_EB_40->Draw("same");
  gr_noPU_EB_60->Draw("same");
  gr_noPU_EB_80->Draw("same");
  gr_noPU_EB_100->Draw("same");
  leg_3_EB->Draw();
  c_noMTD_vs_noMTD_EB_zoomed->Print("roc_reliso_noMTD_vs_noMTD_EB_zoomed.pdf");

  // Endcap
  TCanvas* c_noMTD_vs_40ps_EE = new TCanvas("c_noMTD_vs_40ps_EE","c_noMTD_vs_40ps_EE",1500,1500);
  c_noMTD_vs_40ps_EE->cd();
  c_noMTD_vs_40ps_EE->SetGrid();
  gr_PU200_EE->Draw("AL");
  gr_PU200_EE_40->Draw("same");
  gr_noPU_EE->Draw("same");
  leg_1_EE->Draw();
  c_noMTD_vs_40ps_EE->Print("roc_reliso_noMTD_vs_40ps_EE.pdf");

  TCanvas* c_noMTD_vs_totps_EE = new TCanvas("c_noMTD_vs_totps_EE","c_noMTD_vs_totps_EE",1500,1500);
  c_noMTD_vs_totps_EE->cd();
  c_noMTD_vs_totps_EE->SetGrid();
  gr_PU200_EE->Draw("AL");
  gr_PU200_EE_40->Draw("same");
  gr_PU200_EE_60->Draw("same");
  gr_PU200_EE_80->Draw("same");
  gr_PU200_EE_100->Draw("same");
  gr_noPU_EE->Draw("same");
  leg_2_EE->Draw();
  c_noMTD_vs_totps_EE->Print("roc_reliso_noMTD_vs_totps_EE.pdf");

  TCanvas* c_noMTD_vs_noMTD_EE = new TCanvas("c_noMTD_vs_noMTD_EE","c_noMTD_vs_noMTD_EE",1500,1500);
  c_noMTD_vs_noMTD_EE->cd();
  c_noMTD_vs_noMTD_EE->SetGrid();
  gr_noPU_EE->Draw("AL");
  gr_noPU_EE_40->Draw("same");
  gr_noPU_EE_60->Draw("same");
  gr_noPU_EE_80->Draw("same");
  gr_noPU_EE_100->Draw("same");
  leg_3_EE->Draw();
  c_noMTD_vs_noMTD_EE->Print("roc_reliso_noMTD_vs_noMTD_EE.pdf");

  TCanvas* c_noMTD_vs_40ps_EE_zoomed = new TCanvas("c_noMTD_vs_40ps_EE_zoomed","c_noMTD_vs_40ps_EE_zoomed",1500,1500);
  c_noMTD_vs_40ps_EE_zoomed->cd();
  c_noMTD_vs_40ps_EE_zoomed->SetGrid();
  gr_PU200_EE->GetXaxis()->SetRangeUser(0.8,1.);
  gr_PU200_EE->Draw("AL");
  gr_PU200_EE_40->Draw("same");
  gr_noPU_EE->Draw("same");
  leg_1_EE->Draw();
  c_noMTD_vs_40ps_EE_zoomed->Print("roc_reliso_noMTD_vs_40ps_EE_zoomed.pdf");

  TCanvas* c_noMTD_vs_totps_EE_zoomed = new TCanvas("c_noMTD_vs_totps_EE_zoomed","c_noMTD_vs_totps_EE_zoomed",1500,1500);
  c_noMTD_vs_totps_EE_zoomed->cd();
  c_noMTD_vs_totps_EE_zoomed->SetGrid();
  gr_PU200_EE->GetXaxis()->SetRangeUser(0.8,1.);
  gr_PU200_EE->Draw("AL");
  gr_PU200_EE_40->Draw("same");
  gr_PU200_EE_60->Draw("same");
  gr_PU200_EE_80->Draw("same");
  gr_PU200_EE_100->Draw("same");
  gr_noPU_EE->Draw("same");
  leg_2_EE->Draw();
  c_noMTD_vs_totps_EE_zoomed->Print("roc_reliso_noMTD_vs_totps_EE_zoomed.pdf");

  TCanvas* c_noMTD_vs_noMTD_EE_zoomed = new TCanvas("c_noMTD_vs_noMTD_EE_zoomed","c_noMTD_vs_noMTD_EE_zoomed",1500,1500);
  c_noMTD_vs_noMTD_EE_zoomed->cd();
  c_noMTD_vs_noMTD_EE_zoomed->SetGrid();
  gr_noPU_EE->GetXaxis()->SetRangeUser(0.8,1.);
  gr_noPU_EE->Draw("AL");
  gr_noPU_EE_40->Draw("same");
  gr_noPU_EE_60->Draw("same");
  gr_noPU_EE_80->Draw("same");
  gr_noPU_EE_100->Draw("same");
  leg_3_EE->Draw();
  c_noMTD_vs_noMTD_EE_zoomed->Print("roc_reliso_noMTD_vs_noMTD_EE_zoomed.pdf");



}


int main(int argc, char **argv)
{


  roc();


  return 0;

}
