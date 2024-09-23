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
#include "THStack.h"

using namespace std;

void draw_iso_efficiency_ntuple();
void draw_track_type_sigma_ntuple();

TString path_PU200_prompt        = "240623/harvester_PU200_prompt_muon.root";
TString path_noPU_prompt         = "240623/harvester_noPU_prompt_muon.root";
TString path_PU200_prompt_vtx    = "240610/harvester_PU200_prompt_muon_vtx.root";
TString path_noPU_prompt_vtx     = "240610/harvester_noPU_prompt_muon_vtx.root";
TString path_PU200_nonprompt     = "240623/harvester_PU200_nonprompt_muon_qcd.root";
TString path_noPU_nonprompt      = "240623/harvester_noPU_nonprompt_muon_qcd.root";
TString path_PU200_nonprompt_vtx = "240610/harvester_PU200_nonprompt_muon_qcd_vtx.root";
TString path_noPU_nonprompt_vtx  = "240610/harvester_noPU_nonprompt_muon_qcd_vtx.root";
//TString path_PU200_nonprompt     = "240610/harvester_PU200_nonprompt_muon_ttbar.root";
//TString path_noPU_nonprompt      = "240610/harvester_noPU_nonprompt_muon_ttbar.root";
//TString path_PU200_nonprompt_vtx = "240610/harvester_PU200_nonprompt_muon_ttbar_vtx.root";
//TString path_noPU_nonprompt_vtx  = "240610/harvester_noPU_nonprompt_muon_ttbar_vtx.root";

TString ntuple_PU200_prompt        = "240623/ntuple_PU200_prompt_muon.root";
TString ntuple_noPU_prompt         = "240623/ntuple_noPU_prompt_muon.root";
TString ntuple_PU200_prompt_vtx    = "240610/ntuple_PU200_prompt_muon_vtx.root";
TString ntuple_noPU_prompt_vtx     = "240610/ntuple_noPU_prompt_muon_vtx.root";
TString ntuple_PU200_nonprompt     = "240623/ntuple_PU200_nonprompt_muon_qcd.root";
TString ntuple_noPU_nonprompt      = "240623/ntuple_noPU_nonprompt_muon_qcd.root";
TString ntuple_PU200_nonprompt_vtx = "240610/ntuple_PU200_nonprompt_muon_qcd_vtx.root";
TString ntuple_noPU_nonprompt_vtx  = "240610/ntuple_noPU_nonprompt_muon_qcd_vtx.root";
//TString ntuple_PU200_nonprompt     = "240610/ntuple_PU200_nonprompt_muon_ttbar.root";
//TString ntuple_noPU_nonprompt      = "240610/ntuple_noPU_nonprompt_muon_ttbar.root";
//TString ntuple_PU200_nonprompt_vtx = "240610/ntuple_PU200_nonprompt_muon_ttbar_vtx.root";
//TString ntuple_noPU_nonprompt_vtx  = "240610/ntuple_noPU_nonprompt_muon_ttbar_vtx.root";

 





void draw_iso_efficiency_ntuple() {

  // generate the dictionary for vector<vector<sth>> collection
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<bool> >", "vector");

  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  TFile *f_PU200_prompt_vtx, *f_PU200_nonprompt_vtx, *f_noPU_prompt_vtx, *f_noPU_nonprompt_vtx;


  // muon track
  // PU200
  TChain* ch_PU200_prompt = new TChain("mtdMuonIsoValid/muonIso");
  TChain* ch_PU200_nonprompt = new TChain("mtdMuonIsoValid/muonIso");
  ch_PU200_prompt->Add(Form("data/%s", ntuple_PU200_prompt.Data()));
  ch_PU200_nonprompt->Add(Form("data/%s", ntuple_PU200_nonprompt.Data()));
  // noPU
  TChain* ch_noPU_prompt = new TChain("mtdMuonIsoValid/muonIso");
  TChain* ch_noPU_nonprompt = new TChain("mtdMuonIsoValid/muonIso");
  ch_noPU_prompt->Add(Form("data/%s", ntuple_noPU_prompt.Data()));
  ch_noPU_nonprompt->Add(Form("data/%s", ntuple_noPU_nonprompt.Data()));

  // vertex
  // PU200
  TChain* ch_PU200_prompt_vtx = new TChain("mtdMuonIsoValid/muonIso");
  TChain* ch_PU200_nonprompt_vtx = new TChain("mtdMuonIsoValid/muonIso");
  ch_PU200_prompt_vtx->Add(Form("data/%s", ntuple_PU200_prompt_vtx.Data()));
  ch_PU200_nonprompt_vtx->Add(Form("data/%s", ntuple_PU200_nonprompt_vtx.Data()));
  // noPU
  TChain* ch_noPU_prompt_vtx = new TChain("mtdMuonIsoValid/muonIso");
  TChain* ch_noPU_nonprompt_vtx = new TChain("mtdMuonIsoValid/muonIso");
  ch_noPU_prompt_vtx->Add(Form("data/%s", ntuple_noPU_prompt_vtx.Data()));
  ch_noPU_nonprompt_vtx->Add(Form("data/%s", ntuple_noPU_nonprompt_vtx.Data()));

  int nbin=1000;

  // muon
  vector<float> *muon_pt_(0), *muon_time_(0), *muon_time_err_(0), *muon_prompt_(0), *muon_isBarrel_(0);
  vector<float> *muon_status_(0), *muon_pv_dz_(0), *muon_pv_dxy_(0), *muon_vz_(0);
  vector<float> *muon_PVweight_(0);
  vector<bool>  *muon_isMuon_(0), *muon_isPFMuon_(0), *muon_isGlobalMuon_(0), *muon_isTrackerMuon_(0), *muon_isStandAloneMuon_(0), *muon_isCutBasedIdLoose_(0), *muon_isLooseMuon_(0);
  // track
  std::vector<std::vector<float>> *track_pt_(0), *track_time_(0), *track_time_err_(0);
  std::vector<std::vector<float>> *track_pv_dz_(0), *track_vz_(0);
  std::vector<std::vector<float>> *track_PVweight_(0);
  vector<int>   *muon_index_(0);
  std::vector<std::vector<int>> *track_index_(0);
  std::vector<std::vector<bool>> *selectedLV_(0), *match_vtx_sim2reco_(0), *match_vtx_reco2sim_(0);
  int vtx_index_=999;
  int recovtx_sim_=999, simvtx_reco_=999;
  int simvtx_bx_=999, simvtx_evtId_=999;
  int recovtx_original_index_=999;



  // PU200
  ch_PU200_prompt->SetBranchAddress("muon_pt_",        &muon_pt_);
  ch_PU200_prompt->SetBranchAddress("muon_time_",      &muon_time_);
  ch_PU200_prompt->SetBranchAddress("muon_time_err_",  &muon_time_err_);
  ch_PU200_prompt->SetBranchAddress("muon_prompt_",    &muon_prompt_);
  ch_PU200_prompt->SetBranchAddress("muon_isBarrel_",  &muon_isBarrel_);
  ch_PU200_prompt->SetBranchAddress("muon_status_",    &muon_status_);
  ch_PU200_prompt->SetBranchAddress("muon_pv_dz_",     &muon_pv_dz_);
  ch_PU200_prompt->SetBranchAddress("muon_pv_dxy_",    &muon_pv_dxy_);
  ch_PU200_prompt->SetBranchAddress("muon_vz_",        &muon_vz_);
  ch_PU200_prompt->SetBranchAddress("muon_PVweight_",  &muon_PVweight_);
  ch_PU200_prompt->SetBranchAddress("track_pt_",       &track_pt_);
  ch_PU200_prompt->SetBranchAddress("track_time_",     &track_time_);
  ch_PU200_prompt->SetBranchAddress("track_time_err_", &track_time_err_);
  ch_PU200_prompt->SetBranchAddress("track_vz_",       &track_vz_);
  ch_PU200_prompt->SetBranchAddress("track_pv_dz_",    &track_pv_dz_);
  ch_PU200_prompt->SetBranchAddress("track_PVweight_", &track_PVweight_);
  ch_PU200_prompt->SetBranchAddress("selectedLV_",     &selectedLV_);
  ch_PU200_prompt->SetBranchAddress("match_vtx_sim2reco_",     &match_vtx_sim2reco_);
  ch_PU200_prompt->SetBranchAddress("match_vtx_reco2sim_",     &match_vtx_reco2sim_);
  ch_PU200_prompt->SetBranchAddress("vtx_index_",      &vtx_index_);
  ch_PU200_prompt->SetBranchAddress("recovtx_sim_",    &recovtx_sim_);
  ch_PU200_prompt->SetBranchAddress("simvtx_reco_",    &simvtx_reco_);
  ch_PU200_prompt->SetBranchAddress("simvtx_bx_",      &simvtx_bx_);
  ch_PU200_prompt->SetBranchAddress("simvtx_evtId_",   &simvtx_evtId_);
  ch_PU200_prompt->SetBranchAddress("recovtx_original_index_", &recovtx_original_index_);
  ch_PU200_prompt->SetBranchAddress("muon_isMuon_",    &muon_isMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isPFMuon_",    &muon_isPFMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isGlobalMuon_",    &muon_isGlobalMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isTrackerMuon_",    &muon_isTrackerMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isStandAloneMuon_",    &muon_isStandAloneMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isCutBasedIdLoose_",    &muon_isCutBasedIdLoose_);
  ch_PU200_prompt->SetBranchAddress("muon_isLooseMuon_",    &muon_isLooseMuon_);

  ch_PU200_nonprompt->SetBranchAddress("muon_pt_",        &muon_pt_);
  ch_PU200_nonprompt->SetBranchAddress("muon_time_",      &muon_time_);
  ch_PU200_nonprompt->SetBranchAddress("muon_time_err_",  &muon_time_err_);
  ch_PU200_nonprompt->SetBranchAddress("muon_prompt_",    &muon_prompt_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isBarrel_",  &muon_isBarrel_);
  ch_PU200_nonprompt->SetBranchAddress("muon_status_",    &muon_status_);
  ch_PU200_nonprompt->SetBranchAddress("muon_pv_dz_",     &muon_pv_dz_);
  ch_PU200_nonprompt->SetBranchAddress("muon_pv_dxy_",    &muon_pv_dxy_);
  ch_PU200_nonprompt->SetBranchAddress("track_pt_",       &track_pt_);
  ch_PU200_nonprompt->SetBranchAddress("track_time_",     &track_time_);
  ch_PU200_nonprompt->SetBranchAddress("track_time_err_", &track_time_err_);
  ch_PU200_nonprompt->SetBranchAddress("muon_vz_",        &muon_vz_);
  ch_PU200_nonprompt->SetBranchAddress("track_vz_",       &track_vz_);
  ch_PU200_nonprompt->SetBranchAddress("track_pv_dz_",    &track_pv_dz_);
  ch_PU200_nonprompt->SetBranchAddress("muon_vz_",        &muon_vz_);
  ch_PU200_nonprompt->SetBranchAddress("muon_PVweight_",  &muon_PVweight_);
  ch_PU200_nonprompt->SetBranchAddress("track_PVweight_", &track_PVweight_);
  ch_PU200_nonprompt->SetBranchAddress("selectedLV_",     &selectedLV_);
  ch_PU200_nonprompt->SetBranchAddress("match_vtx_sim2reco_",     &match_vtx_sim2reco_);
  ch_PU200_nonprompt->SetBranchAddress("match_vtx_reco2sim_",     &match_vtx_reco2sim_);
  ch_PU200_nonprompt->SetBranchAddress("vtx_index_",      &vtx_index_);
  ch_PU200_nonprompt->SetBranchAddress("recovtx_sim_",    &recovtx_sim_);
  ch_PU200_nonprompt->SetBranchAddress("simvtx_reco_",    &simvtx_reco_);
  ch_PU200_nonprompt->SetBranchAddress("simvtx_bx_",      &simvtx_bx_);
  ch_PU200_nonprompt->SetBranchAddress("simvtx_evtId_",   &simvtx_evtId_);
  ch_PU200_nonprompt->SetBranchAddress("recovtx_original_index_", &recovtx_original_index_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isMuon_",    &muon_isMuon_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isPFMuon_",    &muon_isPFMuon_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isGlobalMuon_",    &muon_isGlobalMuon_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isTrackerMuon_",    &muon_isTrackerMuon_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isStandAloneMuon_",    &muon_isStandAloneMuon_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isCutBasedIdLoose_",    &muon_isCutBasedIdLoose_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isLooseMuon_",    &muon_isLooseMuon_);
  // noPU
  ch_noPU_prompt->SetBranchAddress("muon_pt_",        &muon_pt_);
  ch_noPU_prompt->SetBranchAddress("muon_time_",      &muon_time_);
  ch_noPU_prompt->SetBranchAddress("muon_time_err_",  &muon_time_err_);
  ch_noPU_prompt->SetBranchAddress("muon_prompt_",    &muon_prompt_);
  ch_noPU_prompt->SetBranchAddress("muon_isBarrel_",  &muon_isBarrel_);
  ch_noPU_prompt->SetBranchAddress("muon_status_",    &muon_status_);
  ch_noPU_prompt->SetBranchAddress("muon_pv_dz_",     &muon_pv_dz_);
  ch_noPU_prompt->SetBranchAddress("muon_pv_dxy_",    &muon_pv_dxy_);
  ch_noPU_prompt->SetBranchAddress("track_pt_",       &track_pt_);
  ch_noPU_prompt->SetBranchAddress("track_time_",     &track_time_);
  ch_noPU_prompt->SetBranchAddress("track_time_err_", &track_time_err_);
  ch_noPU_prompt->SetBranchAddress("muon_vz_",        &muon_vz_);
  ch_noPU_prompt->SetBranchAddress("track_vz_",       &track_vz_);
  ch_noPU_prompt->SetBranchAddress("track_pv_dz_",    &track_pv_dz_);
  ch_noPU_prompt->SetBranchAddress("muon_PVweight_",  &muon_PVweight_);
  ch_noPU_prompt->SetBranchAddress("track_PVweight_", &track_PVweight_);
  ch_noPU_prompt->SetBranchAddress("selectedLV_",     &selectedLV_);
  ch_noPU_prompt->SetBranchAddress("match_vtx_sim2reco_",     &match_vtx_sim2reco_);
  ch_noPU_prompt->SetBranchAddress("match_vtx_reco2sim_",     &match_vtx_reco2sim_);
  ch_noPU_prompt->SetBranchAddress("vtx_index_",      &vtx_index_);
  ch_noPU_prompt->SetBranchAddress("recovtx_sim_",    &recovtx_sim_);
  ch_noPU_prompt->SetBranchAddress("simvtx_reco_",    &simvtx_reco_);
  ch_noPU_prompt->SetBranchAddress("simvtx_bx_",      &simvtx_bx_);
  ch_noPU_prompt->SetBranchAddress("simvtx_evtId_",   &simvtx_evtId_);
  ch_noPU_prompt->SetBranchAddress("recovtx_original_index_", &recovtx_original_index_);
  ch_noPU_prompt->SetBranchAddress("muon_isMuon_",    &muon_isMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isPFMuon_",    &muon_isPFMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isGlobalMuon_",    &muon_isGlobalMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isTrackerMuon_",    &muon_isTrackerMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isStandAloneMuon_",    &muon_isStandAloneMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isCutBasedIdLoose_",    &muon_isCutBasedIdLoose_);
  ch_noPU_prompt->SetBranchAddress("muon_isLooseMuon_",    &muon_isLooseMuon_);

  ch_noPU_nonprompt->SetBranchAddress("muon_pt_",        &muon_pt_);
  ch_noPU_nonprompt->SetBranchAddress("muon_time_",      &muon_time_);
  ch_noPU_nonprompt->SetBranchAddress("muon_time_err_",  &muon_time_err_);
  ch_noPU_nonprompt->SetBranchAddress("muon_prompt_",    &muon_prompt_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isBarrel_",  &muon_isBarrel_);
  ch_noPU_nonprompt->SetBranchAddress("muon_status_",    &muon_status_);
  ch_noPU_nonprompt->SetBranchAddress("muon_pv_dz_",     &muon_pv_dz_);
  ch_noPU_nonprompt->SetBranchAddress("muon_pv_dxy_",    &muon_pv_dxy_);
  ch_noPU_nonprompt->SetBranchAddress("track_pt_",       &track_pt_);
  ch_noPU_nonprompt->SetBranchAddress("track_time_",     &track_time_);
  ch_noPU_nonprompt->SetBranchAddress("track_time_err_", &track_time_err_);
  ch_noPU_nonprompt->SetBranchAddress("muon_vz_",        &muon_vz_);
  ch_noPU_nonprompt->SetBranchAddress("track_vz_",       &track_vz_);
  ch_noPU_nonprompt->SetBranchAddress("track_pv_dz_",    &track_pv_dz_);
  ch_noPU_nonprompt->SetBranchAddress("muon_PVweight_",  &muon_PVweight_);
  ch_noPU_nonprompt->SetBranchAddress("track_PVweight_", &track_PVweight_);
  ch_noPU_nonprompt->SetBranchAddress("selectedLV_",     &selectedLV_);
  ch_noPU_nonprompt->SetBranchAddress("match_vtx_sim2reco_",     &match_vtx_sim2reco_);
  ch_noPU_nonprompt->SetBranchAddress("match_vtx_reco2sim_",     &match_vtx_reco2sim_);
  ch_noPU_nonprompt->SetBranchAddress("vtx_index_",      &vtx_index_);
  ch_noPU_nonprompt->SetBranchAddress("recovtx_sim_",    &recovtx_sim_);
  ch_noPU_nonprompt->SetBranchAddress("simvtx_reco_",    &simvtx_reco_);
  ch_noPU_nonprompt->SetBranchAddress("simvtx_bx_",      &simvtx_bx_);
  ch_noPU_nonprompt->SetBranchAddress("simvtx_evtId_",   &simvtx_evtId_);
  ch_noPU_nonprompt->SetBranchAddress("recovtx_original_index_", &recovtx_original_index_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isMuon_",    &muon_isMuon_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isPFMuon_",    &muon_isPFMuon_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isGlobalMuon_",    &muon_isGlobalMuon_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isTrackerMuon_",    &muon_isTrackerMuon_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isStandAloneMuon_",    &muon_isStandAloneMuon_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isCutBasedIdLoose_",    &muon_isCutBasedIdLoose_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isLooseMuon_",    &muon_isLooseMuon_);


  ///////////////////////
  // Define histograms //
  ///////////////////////

  // PU200
    // prompt
  TH1D* h_PU200_prompt_EB         = new TH1D("h_PU200_prompt_EB",           "h_PU200_prompt_EB",           nbin, 0, 4);
  TH1D* h_PU200_prompt_1sigma_EB  = new TH1D("h_PU200_prompt_1sigma_EB",    "h_PU200_prompt_1sigma_EB",    nbin, 0, 4);
  TH1D* h_PU200_prompt_2sigma_EB  = new TH1D("h_PU200_prompt_2sigma_EB",    "h_PU200_prompt_2sigma_EB",    nbin, 0, 4);
  TH1D* h_PU200_prompt_3sigma_EB  = new TH1D("h_PU200_prompt_3sigma_EB",    "h_PU200_prompt_3sigma_EB",    nbin, 0, 4);
  TH1D* h_PU200_prompt_4sigma_EB  = new TH1D("h_PU200_prompt_4sigma_EB",    "h_PU200_prompt_4sigma_EB",    nbin, 0, 4);
  TH1D* h_PU200_prompt_40_EB      = new TH1D("h_PU200_prompt_40_EB",        "h_PU200_prompt_40_EB",        nbin, 0, 4);
  TH1D* h_PU200_prompt_60_EB      = new TH1D("h_PU200_prompt_60_EB",        "h_PU200_prompt_60_EB",        nbin, 0, 4);
  TH1D* h_PU200_prompt_80_EB      = new TH1D("h_PU200_prompt_80_EB",        "h_PU200_prompt_80_EB",        nbin, 0, 4);
  TH1D* h_PU200_prompt_100_EB     = new TH1D("h_PU200_prompt_100_EB",       "h_PU200_prompt_100_EB",       nbin, 0, 4);
  TH1D* h_PU200_prompt_EE         = new TH1D("h_PU200_prompt_EE",           "h_PU200_prompt_EE",           nbin, 0, 4);
  TH1D* h_PU200_prompt_1sigma_EE  = new TH1D("h_PU200_prompt_1sigma_EE",    "h_PU200_prompt_1sigma_EE",    nbin, 0, 4);
  TH1D* h_PU200_prompt_2sigma_EE  = new TH1D("h_PU200_prompt_2sigma_EE",    "h_PU200_prompt_2sigma_EE",    nbin, 0, 4);
  TH1D* h_PU200_prompt_3sigma_EE  = new TH1D("h_PU200_prompt_3sigma_EE",    "h_PU200_prompt_3sigma_EE",    nbin, 0, 4);
  TH1D* h_PU200_prompt_4sigma_EE  = new TH1D("h_PU200_prompt_4sigma_EE",    "h_PU200_prompt_4sigma_EE",    nbin, 0, 4);
  TH1D* h_PU200_prompt_40_EE      = new TH1D("h_PU200_prompt_40_EE",        "h_PU200_prompt_40_EE",        nbin, 0, 4);
  TH1D* h_PU200_prompt_60_EE      = new TH1D("h_PU200_prompt_60_EE",        "h_PU200_prompt_60_EE",        nbin, 0, 4);
  TH1D* h_PU200_prompt_80_EE      = new TH1D("h_PU200_prompt_80_EE",        "h_PU200_prompt_80_EE",        nbin, 0, 4);
  TH1D* h_PU200_prompt_100_EE     = new TH1D("h_PU200_prompt_100_EE",       "h_PU200_prompt_100_EE",       nbin, 0, 4);
  float sumiso_PU200_prompt_EB=0, sumiso_PU200_prompt_1sigma_EB=0, sumiso_PU200_prompt_2sigma_EB=0, sumiso_PU200_prompt_3sigma_EB=0, sumiso_PU200_prompt_4sigma_EB=0, sumiso_PU200_prompt_40_EB=0, sumiso_PU200_prompt_60_EB=0, sumiso_PU200_prompt_80_EB=0, sumiso_PU200_prompt_100_EB=0;
  float sumiso_PU200_prompt_EE=0, sumiso_PU200_prompt_1sigma_EE=0, sumiso_PU200_prompt_2sigma_EE=0, sumiso_PU200_prompt_3sigma_EE=0, sumiso_PU200_prompt_4sigma_EE=0, sumiso_PU200_prompt_40_EE=0, sumiso_PU200_prompt_60_EE=0, sumiso_PU200_prompt_80_EE=0, sumiso_PU200_prompt_100_EE=0;
  float reliso_PU200_prompt_EB=0, reliso_PU200_prompt_1sigma_EB=0, reliso_PU200_prompt_2sigma_EB=0, reliso_PU200_prompt_3sigma_EB=0, reliso_PU200_prompt_4sigma_EB=0, reliso_PU200_prompt_40_EB=0, reliso_PU200_prompt_60_EB=0, reliso_PU200_prompt_80_EB=0, reliso_PU200_prompt_100_EB=0;
  float reliso_PU200_prompt_EE=0, reliso_PU200_prompt_1sigma_EE=0, reliso_PU200_prompt_2sigma_EE=0, reliso_PU200_prompt_3sigma_EE=0, reliso_PU200_prompt_4sigma_EE=0, reliso_PU200_prompt_40_EE=0, reliso_PU200_prompt_60_EE=0, reliso_PU200_prompt_80_EE=0, reliso_PU200_prompt_100_EE=0;
    // nonprompt
  TH1D* h_PU200_nonprompt_EB         = new TH1D("h_PU200_nonprompt_EB",           "h_PU200_nonprompt_EB",           nbin, 0, 4);
  TH1D* h_PU200_nonprompt_1sigma_EB  = new TH1D("h_PU200_nonprompt_1sigma_EB",    "h_PU200_nonprompt_1sigma_EB",    nbin, 0, 4);
  TH1D* h_PU200_nonprompt_2sigma_EB  = new TH1D("h_PU200_nonprompt_2sigma_EB",    "h_PU200_nonprompt_2sigma_EB",    nbin, 0, 4);
  TH1D* h_PU200_nonprompt_3sigma_EB  = new TH1D("h_PU200_nonprompt_3sigma_EB",    "h_PU200_nonprompt_3sigma_EB",    nbin, 0, 4);
  TH1D* h_PU200_nonprompt_4sigma_EB  = new TH1D("h_PU200_nonprompt_4sigma_EB",    "h_PU200_nonprompt_4sigma_EB",    nbin, 0, 4);
  TH1D* h_PU200_nonprompt_40_EB      = new TH1D("h_PU200_nonprompt_40_EB",        "h_PU200_nonprompt_40_EB",        nbin, 0, 4);
  TH1D* h_PU200_nonprompt_60_EB      = new TH1D("h_PU200_nonprompt_60_EB",        "h_PU200_nonprompt_60_EB",        nbin, 0, 4);
  TH1D* h_PU200_nonprompt_80_EB      = new TH1D("h_PU200_nonprompt_80_EB",        "h_PU200_nonprompt_80_EB",        nbin, 0, 4);
  TH1D* h_PU200_nonprompt_100_EB     = new TH1D("h_PU200_nonprompt_100_EB",       "h_PU200_nonprompt_100_EB",       nbin, 0, 4);
  TH1D* h_PU200_nonprompt_EE         = new TH1D("h_PU200_nonprompt_EE",           "h_PU200_nonprompt_EE",           nbin, 0, 4);
  TH1D* h_PU200_nonprompt_1sigma_EE  = new TH1D("h_PU200_nonprompt_1sigma_EE",    "h_PU200_nonprompt_1sigma_EE",    nbin, 0, 4);
  TH1D* h_PU200_nonprompt_2sigma_EE  = new TH1D("h_PU200_nonprompt_2sigma_EE",    "h_PU200_nonprompt_2sigma_EE",    nbin, 0, 4);
  TH1D* h_PU200_nonprompt_3sigma_EE  = new TH1D("h_PU200_nonprompt_3sigma_EE",    "h_PU200_nonprompt_3sigma_EE",    nbin, 0, 4);
  TH1D* h_PU200_nonprompt_4sigma_EE  = new TH1D("h_PU200_nonprompt_4sigma_EE",    "h_PU200_nonprompt_4sigma_EE",    nbin, 0, 4);
  TH1D* h_PU200_nonprompt_40_EE      = new TH1D("h_PU200_nonprompt_40_EE",        "h_PU200_nonprompt_40_EE",        nbin, 0, 4);
  TH1D* h_PU200_nonprompt_60_EE      = new TH1D("h_PU200_nonprompt_60_EE",        "h_PU200_nonprompt_60_EE",        nbin, 0, 4);
  TH1D* h_PU200_nonprompt_80_EE      = new TH1D("h_PU200_nonprompt_80_EE",        "h_PU200_nonprompt_80_EE",        nbin, 0, 4);
  TH1D* h_PU200_nonprompt_100_EE     = new TH1D("h_PU200_nonprompt_100_EE",       "h_PU200_nonprompt_100_EE",       nbin, 0, 4);
  float sumiso_PU200_nonprompt_EB=0, sumiso_PU200_nonprompt_1sigma_EB=0, sumiso_PU200_nonprompt_2sigma_EB=0, sumiso_PU200_nonprompt_3sigma_EB=0, sumiso_PU200_nonprompt_4sigma_EB=0, sumiso_PU200_nonprompt_40_EB=0, sumiso_PU200_nonprompt_60_EB=0, sumiso_PU200_nonprompt_80_EB=0, sumiso_PU200_nonprompt_100_EB=0;
  float sumiso_PU200_nonprompt_EE=0, sumiso_PU200_nonprompt_1sigma_EE=0, sumiso_PU200_nonprompt_2sigma_EE=0, sumiso_PU200_nonprompt_3sigma_EE=0, sumiso_PU200_nonprompt_4sigma_EE=0, sumiso_PU200_nonprompt_40_EE=0, sumiso_PU200_nonprompt_60_EE=0, sumiso_PU200_nonprompt_80_EE=0, sumiso_PU200_nonprompt_100_EE=0;
  float reliso_PU200_nonprompt_EB=0, reliso_PU200_nonprompt_1sigma_EB=0, reliso_PU200_nonprompt_2sigma_EB=0, reliso_PU200_nonprompt_3sigma_EB=0, reliso_PU200_nonprompt_4sigma_EB=0, reliso_PU200_nonprompt_40_EB=0, reliso_PU200_nonprompt_60_EB=0, reliso_PU200_nonprompt_80_EB=0, reliso_PU200_nonprompt_100_EB=0;
  float reliso_PU200_nonprompt_EE=0, reliso_PU200_nonprompt_1sigma_EE=0, reliso_PU200_nonprompt_2sigma_EE=0, reliso_PU200_nonprompt_3sigma_EE=0, reliso_PU200_nonprompt_4sigma_EE=0, reliso_PU200_nonprompt_40_EE=0, reliso_PU200_nonprompt_60_EE=0, reliso_PU200_nonprompt_80_EE=0, reliso_PU200_nonprompt_100_EE=0;

  // noPU
    // prompt
  TH1D* h_noPU_prompt_EB         = new TH1D("h_noPU_prompt_EB",           "h_noPU_prompt_EB",           nbin, 0, 4);
  TH1D* h_noPU_prompt_1sigma_EB  = new TH1D("h_noPU_prompt_1sigma_EB",    "h_noPU_prompt_1sigma_EB",    nbin, 0, 4);
  TH1D* h_noPU_prompt_2sigma_EB  = new TH1D("h_noPU_prompt_2sigma_EB",    "h_noPU_prompt_2sigma_EB",    nbin, 0, 4);
  TH1D* h_noPU_prompt_3sigma_EB  = new TH1D("h_noPU_prompt_3sigma_EB",    "h_noPU_prompt_3sigma_EB",    nbin, 0, 4);
  TH1D* h_noPU_prompt_4sigma_EB  = new TH1D("h_noPU_prompt_4sigma_EB",    "h_noPU_prompt_4sigma_EB",    nbin, 0, 4);
  TH1D* h_noPU_prompt_40_EB      = new TH1D("h_noPU_prompt_40_EB",        "h_noPU_prompt_40_EB",        nbin, 0, 4);
  TH1D* h_noPU_prompt_60_EB      = new TH1D("h_noPU_prompt_60_EB",        "h_noPU_prompt_60_EB",        nbin, 0, 4);
  TH1D* h_noPU_prompt_80_EB      = new TH1D("h_noPU_prompt_80_EB",        "h_noPU_prompt_80_EB",        nbin, 0, 4);
  TH1D* h_noPU_prompt_100_EB     = new TH1D("h_noPU_prompt_100_EB",       "h_noPU_prompt_100_EB",       nbin, 0, 4);
  TH1D* h_noPU_prompt_EE         = new TH1D("h_noPU_prompt_EE",           "h_noPU_prompt_EE",           nbin, 0, 4);
  TH1D* h_noPU_prompt_1sigma_EE  = new TH1D("h_noPU_prompt_1sigma_EE",    "h_noPU_prompt_1sigma_EE",    nbin, 0, 4);
  TH1D* h_noPU_prompt_2sigma_EE  = new TH1D("h_noPU_prompt_2sigma_EE",    "h_noPU_prompt_2sigma_EE",    nbin, 0, 4);
  TH1D* h_noPU_prompt_3sigma_EE  = new TH1D("h_noPU_prompt_3sigma_EE",    "h_noPU_prompt_3sigma_EE",    nbin, 0, 4);
  TH1D* h_noPU_prompt_4sigma_EE  = new TH1D("h_noPU_prompt_4sigma_EE",    "h_noPU_prompt_4sigma_EE",    nbin, 0, 4);
  TH1D* h_noPU_prompt_40_EE      = new TH1D("h_noPU_prompt_40_EE",        "h_noPU_prompt_40_EE",        nbin, 0, 4);
  TH1D* h_noPU_prompt_60_EE      = new TH1D("h_noPU_prompt_60_EE",        "h_noPU_prompt_60_EE",        nbin, 0, 4);
  TH1D* h_noPU_prompt_80_EE      = new TH1D("h_noPU_prompt_80_EE",        "h_noPU_prompt_80_EE",        nbin, 0, 4);
  TH1D* h_noPU_prompt_100_EE     = new TH1D("h_noPU_prompt_100_EE",       "h_noPU_prompt_100_EE",       nbin, 0, 4);
  float sumiso_noPU_prompt_EB=0, sumiso_noPU_prompt_1sigma_EB=0, sumiso_noPU_prompt_2sigma_EB=0, sumiso_noPU_prompt_3sigma_EB=0, sumiso_noPU_prompt_4sigma_EB=0, sumiso_noPU_prompt_40_EB=0, sumiso_noPU_prompt_60_EB=0, sumiso_noPU_prompt_80_EB=0, sumiso_noPU_prompt_100_EB=0;
  float sumiso_noPU_prompt_EE=0, sumiso_noPU_prompt_1sigma_EE=0, sumiso_noPU_prompt_2sigma_EE=0, sumiso_noPU_prompt_3sigma_EE=0, sumiso_noPU_prompt_4sigma_EE=0, sumiso_noPU_prompt_40_EE=0, sumiso_noPU_prompt_60_EE=0, sumiso_noPU_prompt_80_EE=0, sumiso_noPU_prompt_100_EE=0;
  float reliso_noPU_prompt_EB=0, reliso_noPU_prompt_1sigma_EB=0, reliso_noPU_prompt_2sigma_EB=0, reliso_noPU_prompt_3sigma_EB=0, reliso_noPU_prompt_4sigma_EB=0, reliso_noPU_prompt_40_EB=0, reliso_noPU_prompt_60_EB=0, reliso_noPU_prompt_80_EB=0, reliso_noPU_prompt_100_EB=0;
  float reliso_noPU_prompt_EE=0, reliso_noPU_prompt_1sigma_EE=0, reliso_noPU_prompt_2sigma_EE=0, reliso_noPU_prompt_3sigma_EE=0, reliso_noPU_prompt_4sigma_EE=0, reliso_noPU_prompt_40_EE=0, reliso_noPU_prompt_60_EE=0, reliso_noPU_prompt_80_EE=0, reliso_noPU_prompt_100_EE=0;
    // nonprompt
  TH1D* h_noPU_nonprompt_EB         = new TH1D("h_noPU_nonprompt_EB",           "h_noPU_nonprompt_EB",           nbin, 0, 4);
  TH1D* h_noPU_nonprompt_1sigma_EB  = new TH1D("h_noPU_nonprompt_1sigma_EB",    "h_noPU_nonprompt_1sigma_EB",    nbin, 0, 4);
  TH1D* h_noPU_nonprompt_2sigma_EB  = new TH1D("h_noPU_nonprompt_2sigma_EB",    "h_noPU_nonprompt_2sigma_EB",    nbin, 0, 4);
  TH1D* h_noPU_nonprompt_3sigma_EB  = new TH1D("h_noPU_nonprompt_3sigma_EB",    "h_noPU_nonprompt_3sigma_EB",    nbin, 0, 4);
  TH1D* h_noPU_nonprompt_4sigma_EB  = new TH1D("h_noPU_nonprompt_4sigma_EB",    "h_noPU_nonprompt_4sigma_EB",    nbin, 0, 4);
  TH1D* h_noPU_nonprompt_40_EB      = new TH1D("h_noPU_nonprompt_40_EB",        "h_noPU_nonprompt_40_EB",        nbin, 0, 4);
  TH1D* h_noPU_nonprompt_60_EB      = new TH1D("h_noPU_nonprompt_60_EB",        "h_noPU_nonprompt_60_EB",        nbin, 0, 4);
  TH1D* h_noPU_nonprompt_80_EB      = new TH1D("h_noPU_nonprompt_80_EB",        "h_noPU_nonprompt_80_EB",        nbin, 0, 4);
  TH1D* h_noPU_nonprompt_100_EB     = new TH1D("h_noPU_nonprompt_100_EB",       "h_noPU_nonprompt_100_EB",       nbin, 0, 4);
  TH1D* h_noPU_nonprompt_EE         = new TH1D("h_noPU_nonprompt_EE",           "h_noPU_nonprompt_EE",           nbin, 0, 4);
  TH1D* h_noPU_nonprompt_1sigma_EE  = new TH1D("h_noPU_nonprompt_1sigma_EE",    "h_noPU_nonprompt_1sigma_EE",    nbin, 0, 4);
  TH1D* h_noPU_nonprompt_2sigma_EE  = new TH1D("h_noPU_nonprompt_2sigma_EE",    "h_noPU_nonprompt_2sigma_EE",    nbin, 0, 4);
  TH1D* h_noPU_nonprompt_3sigma_EE  = new TH1D("h_noPU_nonprompt_3sigma_EE",    "h_noPU_nonprompt_3sigma_EE",    nbin, 0, 4);
  TH1D* h_noPU_nonprompt_4sigma_EE  = new TH1D("h_noPU_nonprompt_4sigma_EE",    "h_noPU_nonprompt_4sigma_EE",    nbin, 0, 4);
  TH1D* h_noPU_nonprompt_40_EE      = new TH1D("h_noPU_nonprompt_40_EE",        "h_noPU_nonprompt_40_EE",        nbin, 0, 4);
  TH1D* h_noPU_nonprompt_60_EE      = new TH1D("h_noPU_nonprompt_60_EE",        "h_noPU_nonprompt_60_EE",        nbin, 0, 4);
  TH1D* h_noPU_nonprompt_80_EE      = new TH1D("h_noPU_nonprompt_80_EE",        "h_noPU_nonprompt_80_EE",        nbin, 0, 4);
  TH1D* h_noPU_nonprompt_100_EE     = new TH1D("h_noPU_nonprompt_100_EE",       "h_noPU_nonprompt_100_EE",       nbin, 0, 4);
  float sumiso_noPU_nonprompt_EB=0, sumiso_noPU_nonprompt_1sigma_EB=0, sumiso_noPU_nonprompt_2sigma_EB=0, sumiso_noPU_nonprompt_3sigma_EB=0, sumiso_noPU_nonprompt_4sigma_EB=0, sumiso_noPU_nonprompt_40_EB=0, sumiso_noPU_nonprompt_60_EB=0, sumiso_noPU_nonprompt_80_EB=0, sumiso_noPU_nonprompt_100_EB=0;
  float sumiso_noPU_nonprompt_EE=0, sumiso_noPU_nonprompt_1sigma_EE=0, sumiso_noPU_nonprompt_2sigma_EE=0, sumiso_noPU_nonprompt_3sigma_EE=0, sumiso_noPU_nonprompt_4sigma_EE=0, sumiso_noPU_nonprompt_40_EE=0, sumiso_noPU_nonprompt_60_EE=0, sumiso_noPU_nonprompt_80_EE=0, sumiso_noPU_nonprompt_100_EE=0;
  float reliso_noPU_nonprompt_EB=0, reliso_noPU_nonprompt_1sigma_EB=0, reliso_noPU_nonprompt_2sigma_EB=0, reliso_noPU_nonprompt_3sigma_EB=0, reliso_noPU_nonprompt_4sigma_EB=0, reliso_noPU_nonprompt_40_EB=0, reliso_noPU_nonprompt_60_EB=0, reliso_noPU_nonprompt_80_EB=0, reliso_noPU_nonprompt_100_EB=0;
  float reliso_noPU_nonprompt_EE=0, reliso_noPU_nonprompt_1sigma_EE=0, reliso_noPU_nonprompt_2sigma_EE=0, reliso_noPU_nonprompt_3sigma_EE=0, reliso_noPU_nonprompt_4sigma_EE=0, reliso_noPU_nonprompt_40_EE=0, reliso_noPU_nonprompt_60_EE=0, reliso_noPU_nonprompt_80_EE=0, reliso_noPU_nonprompt_100_EE=0;


  bool flag_isMuon            = false;
  bool flag_isPFMuon          = false;
  bool flag_isGlobalMuon      = false;
  bool flag_isTrackerMuon     = false;
  bool flag_isStandAloneMuon  = false;
  bool flag_isCutBasedIdLoose = false;
  bool flag_isLooseMuon       = true;

  bool flag_muon_status=true;
  bool flag_muon_pv_dz=false,     flag_muon_pv_dxy=false;
  bool flag_track_pv_dz=false;
  bool flag_vtx_matching = true,  flag_vtx_selectedLV=true;

  bool flag_muon_PVweight=false;
  bool flag_track_PVweight=false;

  float dtsig_cut=0, dt_cut=0;

  // PARAM
  float muon_PVweight_cut=0.000000000000001;
  float track_PVweight_cut=0.000000000000001;
  float muon_pv_dz_cut_EB=0.5;
  float muon_pv_dz_cut_EE=0.5;
  float muon_pv_dxy_cut_EB=0.2;
  float muon_pv_dxy_cut_EE=0.2;
  float track_pv_dz_cut_EB=0.1;
  float track_pv_dz_cut_EE=0.2;

  //////////////////////
  //// PU200 prompt ////
  //////////////////////

  int n_status_failed_PU200_prompt_EB=0, n_status_failed_PU200_prompt_EE=0;
  int n_muon_PU200_prompt_EB=0, n_muon_PU200_prompt_EE=0;

  // event loop
  for(int ievt=0; ievt<ch_PU200_prompt->GetEntries(); ievt++) {
	ch_PU200_prompt->GetEntry(ievt);

	// Make vtx selection tighter
	if(flag_vtx_matching) {
	  if(simvtx_reco_!=vtx_index_) continue; // it includes reco2sim matching
	}
	if(flag_vtx_selectedLV) {
	  if(!(simvtx_bx_==0 && simvtx_evtId_==0 /*&& recovtx_original_index_==0*/)) continue;
	}

	// muon loop
	for(int im=0; im<muon_pt_->size(); im++) {

	  // Muon selection
      if(flag_isMuon) {
        if(muon_isMuon_->at(im)==0) continue;
      }
      if(flag_isPFMuon) {
        if(muon_isPFMuon_->at(im)==0) continue;
      }
      if(flag_isGlobalMuon) {
        if(muon_isGlobalMuon_->at(im)==0) continue;
      }
      if(flag_isTrackerMuon) {
        if(muon_isTrackerMuon_->at(im)==0) continue;
      }
      if(flag_isStandAloneMuon) {
        if(muon_isStandAloneMuon_->at(im)==0) continue;
      }
      if(flag_isCutBasedIdLoose) {
        if(muon_isCutBasedIdLoose_->at(im)==0) continue;
      }
      if(flag_isLooseMuon) {
        if(muon_isLooseMuon_->at(im)==0) continue;
      }


	  //////////////////////////
	  // prompt muon - Barrel //
	  //////////////////////////
	  if(muon_prompt_->at(im)==1 && muon_isBarrel_->at(im)==1) {
		n_muon_PU200_prompt_EB++;

        if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
		    n_status_failed_PU200_prompt_EB++;
		    continue;
	      }
		}
	    if(flag_muon_pv_dz) {
		  if(muon_pv_dz_->at(im)>muon_pv_dz_cut_EB) continue;
	    }
	    if(flag_muon_pv_dxy) {
	      if(muon_pv_dxy_->at(im)>muon_pv_dxy_cut_EB) continue;
	    }
	    if(flag_muon_PVweight) {
		  if(muon_PVweight_->at(im)<muon_PVweight_cut) continue;
	    }


		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {

		  if(flag_track_pv_dz) {
			//if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) >= 0.1) continue;
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EB) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  // noMTD case
		  sumiso_PU200_prompt_EB += track_pt_->at(im).at(it);
		  // MTD case
		  if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_PU200_prompt_1sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_PU200_prompt_2sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_PU200_prompt_3sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_PU200_prompt_4sigma_EB += track_pt_->at(im).at(it);
		  if( dt_cut<0.12  &&   dt_cut>0 ) sumiso_PU200_prompt_40_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.18  &&   dt_cut>0 ) sumiso_PU200_prompt_60_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.24  &&   dt_cut>0 ) sumiso_PU200_prompt_80_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.30  &&   dt_cut>0 ) sumiso_PU200_prompt_100_EB    += track_pt_->at(im).at(it);
		  else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
            sumiso_PU200_prompt_1sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_2sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_3sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_4sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_40_EB     += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_60_EB     += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_80_EB     += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_100_EB    += track_pt_->at(im).at(it);
		  }
		  dtsig_cut=0;
		  dt_cut=0;
		}
		reliso_PU200_prompt_EB        = sumiso_PU200_prompt_EB/muon_pt_->at(im);
		reliso_PU200_prompt_1sigma_EB = sumiso_PU200_prompt_1sigma_EB/muon_pt_->at(im);
		reliso_PU200_prompt_2sigma_EB = sumiso_PU200_prompt_2sigma_EB/muon_pt_->at(im);
		reliso_PU200_prompt_3sigma_EB = sumiso_PU200_prompt_3sigma_EB/muon_pt_->at(im);
		reliso_PU200_prompt_4sigma_EB = sumiso_PU200_prompt_4sigma_EB/muon_pt_->at(im);
		reliso_PU200_prompt_40_EB     = sumiso_PU200_prompt_40_EB/muon_pt_->at(im);
		reliso_PU200_prompt_60_EB     = sumiso_PU200_prompt_60_EB/muon_pt_->at(im);
		reliso_PU200_prompt_80_EB     = sumiso_PU200_prompt_80_EB/muon_pt_->at(im);
		reliso_PU200_prompt_100_EB    = sumiso_PU200_prompt_100_EB/muon_pt_->at(im);

		// Store value of relative isolation
		h_PU200_prompt_EB->Fill(reliso_PU200_prompt_EB);
		h_PU200_prompt_1sigma_EB->Fill(reliso_PU200_prompt_1sigma_EB);
		h_PU200_prompt_2sigma_EB->Fill(reliso_PU200_prompt_2sigma_EB);
		h_PU200_prompt_3sigma_EB->Fill(reliso_PU200_prompt_3sigma_EB);
		h_PU200_prompt_4sigma_EB->Fill(reliso_PU200_prompt_4sigma_EB);
		h_PU200_prompt_40_EB->Fill(reliso_PU200_prompt_40_EB);
		h_PU200_prompt_60_EB->Fill(reliso_PU200_prompt_60_EB);
		h_PU200_prompt_80_EB->Fill(reliso_PU200_prompt_80_EB);
		h_PU200_prompt_100_EB->Fill(reliso_PU200_prompt_100_EB);

		// Initialize
        sumiso_PU200_prompt_EB=0, sumiso_PU200_prompt_1sigma_EB=0, sumiso_PU200_prompt_2sigma_EB=0, sumiso_PU200_prompt_3sigma_EB=0, sumiso_PU200_prompt_4sigma_EB=0, sumiso_PU200_prompt_40_EB=0, sumiso_PU200_prompt_60_EB=0, sumiso_PU200_prompt_80_EB=0, sumiso_PU200_prompt_100_EB=0;
        reliso_PU200_prompt_EB=0, reliso_PU200_prompt_1sigma_EB=0, reliso_PU200_prompt_2sigma_EB=0, reliso_PU200_prompt_3sigma_EB=0, reliso_PU200_prompt_4sigma_EB=0, reliso_PU200_prompt_40_EB=0, reliso_PU200_prompt_60_EB=0, reliso_PU200_prompt_80_EB=0, reliso_PU200_prompt_100_EB=0;
	  }

	  //////////////////////////
	  // prompt muon - Endcap //
	  //////////////////////////
	  else if(muon_prompt_->at(im)==1 && muon_isBarrel_->at(im)==0) {
		n_muon_PU200_prompt_EE++;

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
		    n_status_failed_PU200_prompt_EE++;
		    continue;
	      }
		}
	    if(flag_muon_pv_dz) {
		  if(muon_pv_dz_->at(im)>muon_pv_dz_cut_EE) continue;
	    }
	    if(flag_muon_pv_dxy) {
	      if(muon_pv_dxy_->at(im)>muon_pv_dxy_cut_EE) continue;
	    }
	    if(flag_muon_PVweight) {
		  if(muon_PVweight_->at(im)<muon_PVweight_cut) continue;
	    }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {

		  if(flag_track_pv_dz) {
			//if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) >= 0.2) continue;
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EE) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  // noMTD case
		  sumiso_PU200_prompt_EE += track_pt_->at(im).at(it);
		  // MTD case
		  if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_PU200_prompt_1sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_PU200_prompt_2sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_PU200_prompt_3sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_PU200_prompt_4sigma_EE += track_pt_->at(im).at(it);
		  if( dt_cut<0.12  &&   dt_cut>0 ) sumiso_PU200_prompt_40_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.18  &&   dt_cut>0 ) sumiso_PU200_prompt_60_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.24  &&   dt_cut>0 ) sumiso_PU200_prompt_80_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.30  &&   dt_cut>0 ) sumiso_PU200_prompt_100_EE    += track_pt_->at(im).at(it);
		  else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
            sumiso_PU200_prompt_1sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_2sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_3sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_4sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_40_EE     += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_60_EE     += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_80_EE     += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_100_EE    += track_pt_->at(im).at(it);
		  }
		  dtsig_cut=0;
		  dt_cut=0;
		}
		reliso_PU200_prompt_EE        = sumiso_PU200_prompt_EE/muon_pt_->at(im);
		reliso_PU200_prompt_1sigma_EE = sumiso_PU200_prompt_1sigma_EE/muon_pt_->at(im);
		reliso_PU200_prompt_2sigma_EE = sumiso_PU200_prompt_2sigma_EE/muon_pt_->at(im);
		reliso_PU200_prompt_3sigma_EE = sumiso_PU200_prompt_3sigma_EE/muon_pt_->at(im);
		reliso_PU200_prompt_4sigma_EE = sumiso_PU200_prompt_4sigma_EE/muon_pt_->at(im);
		reliso_PU200_prompt_40_EE     = sumiso_PU200_prompt_40_EE/muon_pt_->at(im);
		reliso_PU200_prompt_60_EE     = sumiso_PU200_prompt_60_EE/muon_pt_->at(im);
		reliso_PU200_prompt_80_EE     = sumiso_PU200_prompt_80_EE/muon_pt_->at(im);
		reliso_PU200_prompt_100_EE    = sumiso_PU200_prompt_100_EE/muon_pt_->at(im);

		// Store value of relative isolation
		h_PU200_prompt_EE->Fill(reliso_PU200_prompt_EE);
		h_PU200_prompt_1sigma_EE->Fill(reliso_PU200_prompt_1sigma_EE);
		h_PU200_prompt_2sigma_EE->Fill(reliso_PU200_prompt_2sigma_EE);
		h_PU200_prompt_3sigma_EE->Fill(reliso_PU200_prompt_3sigma_EE);
		h_PU200_prompt_4sigma_EE->Fill(reliso_PU200_prompt_4sigma_EE);
		h_PU200_prompt_40_EE->Fill(reliso_PU200_prompt_40_EE);
		h_PU200_prompt_60_EE->Fill(reliso_PU200_prompt_60_EE);
		h_PU200_prompt_80_EE->Fill(reliso_PU200_prompt_80_EE);
		h_PU200_prompt_100_EE->Fill(reliso_PU200_prompt_100_EE);

		// Initialize
        sumiso_PU200_prompt_EE=0, sumiso_PU200_prompt_1sigma_EE=0, sumiso_PU200_prompt_2sigma_EE=0, sumiso_PU200_prompt_3sigma_EE=0, sumiso_PU200_prompt_4sigma_EE=0, sumiso_PU200_prompt_40_EE=0, sumiso_PU200_prompt_60_EE=0, sumiso_PU200_prompt_80_EE=0, sumiso_PU200_prompt_100_EE=0;
        reliso_PU200_prompt_EE=0, reliso_PU200_prompt_1sigma_EE=0, reliso_PU200_prompt_2sigma_EE=0, reliso_PU200_prompt_3sigma_EE=0, reliso_PU200_prompt_4sigma_EE=0, reliso_PU200_prompt_40_EE=0, reliso_PU200_prompt_60_EE=0, reliso_PU200_prompt_80_EE=0, reliso_PU200_prompt_100_EE=0;
	  }
	} // End of muon loop
  } // End of event loop


  /////////////////////////
  //// PU200 nonprompt ////
  /////////////////////////
  
  int n_status_failed_PU200_nonprompt_EB=0, n_status_failed_PU200_nonprompt_EE=0;
  int n_muon_PU200_nonprompt_EB=0, n_muon_PU200_nonprompt_EE=0;

  // event loop
  for(int ievt=0; ievt<ch_PU200_nonprompt->GetEntries(); ievt++) {
	ch_PU200_nonprompt->GetEntry(ievt);

	// Make vtx selection tighter
	if(flag_vtx_matching) {
	  if(simvtx_reco_!=vtx_index_) continue; // it includes reco2sim matching
	}
	if(flag_vtx_selectedLV) {
	  if(!(simvtx_bx_==0 && simvtx_evtId_==0 /*&& recovtx_original_index_==0*/)) continue;
	}

	// muon loop
	for(int im=0; im<muon_pt_->size(); im++) {

	  // Muon selection
      if(flag_isMuon) {
        if(muon_isMuon_->at(im)==0) continue;
      }
      if(flag_isPFMuon) {
        if(muon_isPFMuon_->at(im)==0) continue;
      }
      if(flag_isGlobalMuon) {
        if(muon_isGlobalMuon_->at(im)==0) continue;
      }
      if(flag_isTrackerMuon) {
        if(muon_isTrackerMuon_->at(im)==0) continue;
      }
      if(flag_isStandAloneMuon) {
        if(muon_isStandAloneMuon_->at(im)==0) continue;
      }
      if(flag_isCutBasedIdLoose) {
        if(muon_isCutBasedIdLoose_->at(im)==0) continue;
      }
      if(flag_isLooseMuon) {
        if(muon_isLooseMuon_->at(im)==0) continue;
      }


	  /////////////////////////////
	  // nonprompt muon - Barrel //
	  /////////////////////////////
	  if(muon_prompt_->at(im)==0 && muon_isBarrel_->at(im)==1) {
		n_muon_PU200_nonprompt_EB++;

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
		    n_status_failed_PU200_nonprompt_EB++;
		    continue;
	      }
		}
	    if(flag_muon_pv_dz) {
		  if(muon_pv_dz_->at(im)>muon_pv_dz_cut_EB) continue;
	    }
	    if(flag_muon_pv_dxy) {
	      if(muon_pv_dxy_->at(im)>muon_pv_dxy_cut_EB) continue;
	    }
	    if(flag_muon_PVweight) {
		  if(muon_PVweight_->at(im)<muon_PVweight_cut) continue;
	    }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {

		  if(flag_track_pv_dz) {
			//if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) >= 0.1) continue;
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EB) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  // noMTD case
		  sumiso_PU200_nonprompt_EB += track_pt_->at(im).at(it);
		  // MTD case
		  if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_PU200_nonprompt_1sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_PU200_nonprompt_2sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_PU200_nonprompt_3sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_PU200_nonprompt_4sigma_EB += track_pt_->at(im).at(it);
		  if( dt_cut<0.12  &&   dt_cut>0 ) sumiso_PU200_nonprompt_40_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.18  &&   dt_cut>0 ) sumiso_PU200_nonprompt_60_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.24  &&   dt_cut>0 ) sumiso_PU200_nonprompt_80_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.30  &&   dt_cut>0 ) sumiso_PU200_nonprompt_100_EB    += track_pt_->at(im).at(it);
		  else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
            sumiso_PU200_nonprompt_1sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_2sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_3sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_4sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_40_EB     += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_60_EB     += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_80_EB     += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_100_EB    += track_pt_->at(im).at(it);
		  }
		  dtsig_cut=0;
		  dt_cut=0;
		}
		reliso_PU200_nonprompt_EB        = sumiso_PU200_nonprompt_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_1sigma_EB = sumiso_PU200_nonprompt_1sigma_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_2sigma_EB = sumiso_PU200_nonprompt_2sigma_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_3sigma_EB = sumiso_PU200_nonprompt_3sigma_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_4sigma_EB = sumiso_PU200_nonprompt_4sigma_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_40_EB     = sumiso_PU200_nonprompt_40_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_60_EB     = sumiso_PU200_nonprompt_60_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_80_EB     = sumiso_PU200_nonprompt_80_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_100_EB    = sumiso_PU200_nonprompt_100_EB/muon_pt_->at(im);

		// Store value of relative isolation
		h_PU200_nonprompt_EB->Fill(reliso_PU200_nonprompt_EB);
		h_PU200_nonprompt_1sigma_EB->Fill(reliso_PU200_nonprompt_1sigma_EB);
		h_PU200_nonprompt_2sigma_EB->Fill(reliso_PU200_nonprompt_2sigma_EB);
		h_PU200_nonprompt_3sigma_EB->Fill(reliso_PU200_nonprompt_3sigma_EB);
		h_PU200_nonprompt_4sigma_EB->Fill(reliso_PU200_nonprompt_4sigma_EB);
		h_PU200_nonprompt_40_EB->Fill(reliso_PU200_nonprompt_40_EB);
		h_PU200_nonprompt_60_EB->Fill(reliso_PU200_nonprompt_60_EB);
		h_PU200_nonprompt_80_EB->Fill(reliso_PU200_nonprompt_80_EB);
		h_PU200_nonprompt_100_EB->Fill(reliso_PU200_nonprompt_100_EB);

		// Initialize
        sumiso_PU200_nonprompt_EB=0, sumiso_PU200_nonprompt_1sigma_EB=0, sumiso_PU200_nonprompt_2sigma_EB=0, sumiso_PU200_nonprompt_3sigma_EB=0, sumiso_PU200_nonprompt_4sigma_EB=0, sumiso_PU200_nonprompt_40_EB=0, sumiso_PU200_nonprompt_60_EB=0, sumiso_PU200_nonprompt_80_EB=0, sumiso_PU200_nonprompt_100_EB=0;
        reliso_PU200_nonprompt_EB=0, reliso_PU200_nonprompt_1sigma_EB=0, reliso_PU200_nonprompt_2sigma_EB=0, reliso_PU200_nonprompt_3sigma_EB=0, reliso_PU200_nonprompt_4sigma_EB=0, reliso_PU200_nonprompt_40_EB=0, reliso_PU200_nonprompt_60_EB=0, reliso_PU200_nonprompt_80_EB=0, reliso_PU200_nonprompt_100_EB=0;
	  }

	  /////////////////////////////
	  // nonprompt muon - Endcap //
	  /////////////////////////////
	  else if(muon_prompt_->at(im)==0 && muon_isBarrel_->at(im)==0) {
		n_muon_PU200_nonprompt_EE++;

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
		    n_status_failed_PU200_nonprompt_EE++;
		    continue;
	      }
		}
	    if(flag_muon_pv_dz) {
		  if(muon_pv_dz_->at(im)>muon_pv_dz_cut_EE) continue;
	    }
	    if(flag_muon_pv_dxy) {
	      if(muon_pv_dxy_->at(im)>muon_pv_dxy_cut_EE) continue;
	    }
	    if(flag_muon_PVweight) {
		  if(muon_PVweight_->at(im)<muon_PVweight_cut) continue;
	    }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {

		  if(flag_track_pv_dz) {
			//if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) >= 0.2) continue;
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EE) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  // noMTD case
		  sumiso_PU200_nonprompt_EE += track_pt_->at(im).at(it);
		  // MTD case
		  if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_PU200_nonprompt_1sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_PU200_nonprompt_2sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_PU200_nonprompt_3sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_PU200_nonprompt_4sigma_EE += track_pt_->at(im).at(it);
		  if( dt_cut<0.12  &&   dt_cut>0 ) sumiso_PU200_nonprompt_40_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.18  &&   dt_cut>0 ) sumiso_PU200_nonprompt_60_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.24  &&   dt_cut>0 ) sumiso_PU200_nonprompt_80_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.30  &&   dt_cut>0 ) sumiso_PU200_nonprompt_100_EE    += track_pt_->at(im).at(it);
		  else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
            sumiso_PU200_nonprompt_1sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_2sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_3sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_4sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_40_EE     += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_60_EE     += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_80_EE     += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_100_EE    += track_pt_->at(im).at(it);
		  }
		  dtsig_cut=0;
		  dt_cut=0;
		}
		reliso_PU200_nonprompt_EE        = sumiso_PU200_nonprompt_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_1sigma_EE = sumiso_PU200_nonprompt_1sigma_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_2sigma_EE = sumiso_PU200_nonprompt_2sigma_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_3sigma_EE = sumiso_PU200_nonprompt_3sigma_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_4sigma_EE = sumiso_PU200_nonprompt_4sigma_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_40_EE     = sumiso_PU200_nonprompt_40_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_60_EE     = sumiso_PU200_nonprompt_60_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_80_EE     = sumiso_PU200_nonprompt_80_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_100_EE    = sumiso_PU200_nonprompt_100_EE/muon_pt_->at(im);

		// Store value of relative isolation
		h_PU200_nonprompt_EE->Fill(reliso_PU200_nonprompt_EE);
		h_PU200_nonprompt_1sigma_EE->Fill(reliso_PU200_nonprompt_1sigma_EE);
		h_PU200_nonprompt_2sigma_EE->Fill(reliso_PU200_nonprompt_2sigma_EE);
		h_PU200_nonprompt_3sigma_EE->Fill(reliso_PU200_nonprompt_3sigma_EE);
		h_PU200_nonprompt_4sigma_EE->Fill(reliso_PU200_nonprompt_4sigma_EE);
		h_PU200_nonprompt_40_EE->Fill(reliso_PU200_nonprompt_40_EE);
		h_PU200_nonprompt_60_EE->Fill(reliso_PU200_nonprompt_60_EE);
		h_PU200_nonprompt_80_EE->Fill(reliso_PU200_nonprompt_80_EE);
		h_PU200_nonprompt_100_EE->Fill(reliso_PU200_nonprompt_100_EE);

		// Initialize
        sumiso_PU200_nonprompt_EE=0, sumiso_PU200_nonprompt_1sigma_EE=0, sumiso_PU200_nonprompt_2sigma_EE=0, sumiso_PU200_nonprompt_3sigma_EE=0, sumiso_PU200_nonprompt_4sigma_EE=0, sumiso_PU200_nonprompt_40_EE=0, sumiso_PU200_nonprompt_60_EE=0, sumiso_PU200_nonprompt_80_EE=0, sumiso_PU200_nonprompt_100_EE=0;
        reliso_PU200_nonprompt_EE=0, reliso_PU200_nonprompt_1sigma_EE=0, reliso_PU200_nonprompt_2sigma_EE=0, reliso_PU200_nonprompt_3sigma_EE=0, reliso_PU200_nonprompt_4sigma_EE=0, reliso_PU200_nonprompt_40_EE=0, reliso_PU200_nonprompt_60_EE=0, reliso_PU200_nonprompt_80_EE=0, reliso_PU200_nonprompt_100_EE=0;

	  }
	} // End of muon loop
  } // End of event loop
	

  //////////////////////
  //// noPU prompt ////
  //////////////////////
  
  int n_status_failed_noPU_prompt_EB=0, n_status_failed_noPU_prompt_EE=0;
  int n_muon_noPU_prompt_EB=0, n_muon_noPU_prompt_EE=0;

  // event loop
  for(int ievt=0; ievt<ch_noPU_prompt->GetEntries(); ievt++) {
	ch_noPU_prompt->GetEntry(ievt);

	// Make vtx selection tighter
	if(flag_vtx_matching) {
	  if(simvtx_reco_!=vtx_index_) continue; // it includes reco2sim matching
	}
	if(flag_vtx_selectedLV) {
	  if(!(simvtx_bx_==0 && simvtx_evtId_==0 /*&& recovtx_original_index_==0*/)) continue;
	}

	// muon loop
	for(int im=0; im<muon_pt_->size(); im++) {

	  // Muon selection
      if(flag_isMuon) {
        if(muon_isMuon_->at(im)==0) continue;
      }
      if(flag_isPFMuon) {
        if(muon_isPFMuon_->at(im)==0) continue;
      }
      if(flag_isGlobalMuon) {
        if(muon_isGlobalMuon_->at(im)==0) continue;
      }
      if(flag_isTrackerMuon) {
        if(muon_isTrackerMuon_->at(im)==0) continue;
      }
      if(flag_isStandAloneMuon) {
        if(muon_isStandAloneMuon_->at(im)==0) continue;
      }
      if(flag_isCutBasedIdLoose) {
        if(muon_isCutBasedIdLoose_->at(im)==0) continue;
      }
      if(flag_isLooseMuon) {
        if(muon_isLooseMuon_->at(im)==0) continue;
      }


	  //////////////////////////
	  // prompt muon - Barrel //
	  //////////////////////////
	  if(muon_prompt_->at(im)==1 && muon_isBarrel_->at(im)==1) {
		n_muon_noPU_prompt_EB++;

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
		    n_status_failed_noPU_prompt_EB++;
		    continue;
	      }
		}
	    if(flag_muon_pv_dz) {
		  if(muon_pv_dz_->at(im)>muon_pv_dz_cut_EB) continue;
	    }
	    if(flag_muon_pv_dxy) {
	      if(muon_pv_dxy_->at(im)>muon_pv_dxy_cut_EB) continue;
	    }
	    if(flag_muon_PVweight) {
		  if(muon_PVweight_->at(im)<muon_PVweight_cut) continue;
	    }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {

		  if(flag_track_pv_dz) {
			//if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) >= 0.1) continue;
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EB) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  // noMTD case
		  sumiso_noPU_prompt_EB += track_pt_->at(im).at(it);
		  // MTD case
		  if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_noPU_prompt_1sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_noPU_prompt_2sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_noPU_prompt_3sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_noPU_prompt_4sigma_EB += track_pt_->at(im).at(it);
		  if( dt_cut<0.12  &&   dt_cut>0 ) sumiso_noPU_prompt_40_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.18  &&   dt_cut>0 ) sumiso_noPU_prompt_60_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.24  &&   dt_cut>0 ) sumiso_noPU_prompt_80_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.30  &&   dt_cut>0 ) sumiso_noPU_prompt_100_EB    += track_pt_->at(im).at(it);
		  else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
            sumiso_noPU_prompt_1sigma_EB += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_2sigma_EB += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_3sigma_EB += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_4sigma_EB += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_40_EB     += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_60_EB     += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_80_EB     += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_100_EB    += track_pt_->at(im).at(it);
		  }
		  dtsig_cut=0;
		  dt_cut=0;
		}
		reliso_noPU_prompt_EB        = sumiso_noPU_prompt_EB/muon_pt_->at(im);
		reliso_noPU_prompt_1sigma_EB = sumiso_noPU_prompt_1sigma_EB/muon_pt_->at(im);
		reliso_noPU_prompt_2sigma_EB = sumiso_noPU_prompt_2sigma_EB/muon_pt_->at(im);
		reliso_noPU_prompt_3sigma_EB = sumiso_noPU_prompt_3sigma_EB/muon_pt_->at(im);
		reliso_noPU_prompt_4sigma_EB = sumiso_noPU_prompt_4sigma_EB/muon_pt_->at(im);
		reliso_noPU_prompt_40_EB     = sumiso_noPU_prompt_40_EB/muon_pt_->at(im);
		reliso_noPU_prompt_60_EB     = sumiso_noPU_prompt_60_EB/muon_pt_->at(im);
		reliso_noPU_prompt_80_EB     = sumiso_noPU_prompt_80_EB/muon_pt_->at(im);
		reliso_noPU_prompt_100_EB    = sumiso_noPU_prompt_100_EB/muon_pt_->at(im);

		// Store value of relative isolation
		h_noPU_prompt_EB->Fill(reliso_noPU_prompt_EB);
		h_noPU_prompt_1sigma_EB->Fill(reliso_noPU_prompt_1sigma_EB);
		h_noPU_prompt_2sigma_EB->Fill(reliso_noPU_prompt_2sigma_EB);
		h_noPU_prompt_3sigma_EB->Fill(reliso_noPU_prompt_3sigma_EB);
		h_noPU_prompt_4sigma_EB->Fill(reliso_noPU_prompt_4sigma_EB);
		h_noPU_prompt_40_EB->Fill(reliso_noPU_prompt_40_EB);
		h_noPU_prompt_60_EB->Fill(reliso_noPU_prompt_60_EB);
		h_noPU_prompt_80_EB->Fill(reliso_noPU_prompt_80_EB);
		h_noPU_prompt_100_EB->Fill(reliso_noPU_prompt_100_EB);

		// Initialize
        sumiso_noPU_prompt_EB=0, sumiso_noPU_prompt_1sigma_EB=0, sumiso_noPU_prompt_2sigma_EB=0, sumiso_noPU_prompt_3sigma_EB=0, sumiso_noPU_prompt_4sigma_EB=0, sumiso_noPU_prompt_40_EB=0, sumiso_noPU_prompt_60_EB=0, sumiso_noPU_prompt_80_EB=0, sumiso_noPU_prompt_100_EB=0;
        reliso_noPU_prompt_EB=0, reliso_noPU_prompt_1sigma_EB=0, reliso_noPU_prompt_2sigma_EB=0, reliso_noPU_prompt_3sigma_EB=0, reliso_noPU_prompt_4sigma_EB=0, reliso_noPU_prompt_40_EB=0, reliso_noPU_prompt_60_EB=0, reliso_noPU_prompt_80_EB=0, reliso_noPU_prompt_100_EB=0;
	  }

	  //////////////////////////
	  // prompt muon - Endcap //
	  //////////////////////////
	  else if(muon_prompt_->at(im)==1 && muon_isBarrel_->at(im)==0) {
		n_muon_noPU_prompt_EE++;

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
		    n_status_failed_noPU_prompt_EE++;
		    continue;
	      }
		}
	    if(flag_muon_pv_dz) {
		  if(muon_pv_dz_->at(im)>muon_pv_dz_cut_EE) continue;
	    }
	    if(flag_muon_pv_dxy) {
	      if(muon_pv_dxy_->at(im)>muon_pv_dxy_cut_EE) continue;
	    }
	    if(flag_muon_PVweight) {
		  if(muon_PVweight_->at(im)<muon_PVweight_cut) continue;
	    }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {

		  if(flag_track_pv_dz) {
		    //if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) >= 0.2) continue;
		    if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EE) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  // noMTD case
		  sumiso_noPU_prompt_EE += track_pt_->at(im).at(it);
		  // MTD case
		  if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_noPU_prompt_1sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_noPU_prompt_2sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_noPU_prompt_3sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_noPU_prompt_4sigma_EE += track_pt_->at(im).at(it);
		  if( dt_cut<0.12  &&   dt_cut>0 ) sumiso_noPU_prompt_40_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.18  &&   dt_cut>0 ) sumiso_noPU_prompt_60_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.24  &&   dt_cut>0 ) sumiso_noPU_prompt_80_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.30  &&   dt_cut>0 ) sumiso_noPU_prompt_100_EE    += track_pt_->at(im).at(it);
		  else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
            sumiso_noPU_prompt_1sigma_EE += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_2sigma_EE += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_3sigma_EE += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_4sigma_EE += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_40_EE     += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_60_EE     += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_80_EE     += track_pt_->at(im).at(it);
            sumiso_noPU_prompt_100_EE    += track_pt_->at(im).at(it);
		  }
		  dtsig_cut=0;
		  dt_cut=0;
		}
		reliso_noPU_prompt_EE        = sumiso_noPU_prompt_EE/muon_pt_->at(im);
		reliso_noPU_prompt_1sigma_EE = sumiso_noPU_prompt_1sigma_EE/muon_pt_->at(im);
		reliso_noPU_prompt_2sigma_EE = sumiso_noPU_prompt_2sigma_EE/muon_pt_->at(im);
		reliso_noPU_prompt_3sigma_EE = sumiso_noPU_prompt_3sigma_EE/muon_pt_->at(im);
		reliso_noPU_prompt_4sigma_EE = sumiso_noPU_prompt_4sigma_EE/muon_pt_->at(im);
		reliso_noPU_prompt_40_EE     = sumiso_noPU_prompt_40_EE/muon_pt_->at(im);
		reliso_noPU_prompt_60_EE     = sumiso_noPU_prompt_60_EE/muon_pt_->at(im);
		reliso_noPU_prompt_80_EE     = sumiso_noPU_prompt_80_EE/muon_pt_->at(im);
		reliso_noPU_prompt_100_EE    = sumiso_noPU_prompt_100_EE/muon_pt_->at(im);

		// Store value of relative isolation
		h_noPU_prompt_EE->Fill(reliso_noPU_prompt_EE);
		h_noPU_prompt_1sigma_EE->Fill(reliso_noPU_prompt_1sigma_EE);
		h_noPU_prompt_2sigma_EE->Fill(reliso_noPU_prompt_2sigma_EE);
		h_noPU_prompt_3sigma_EE->Fill(reliso_noPU_prompt_3sigma_EE);
		h_noPU_prompt_4sigma_EE->Fill(reliso_noPU_prompt_4sigma_EE);
		h_noPU_prompt_40_EE->Fill(reliso_noPU_prompt_40_EE);
		h_noPU_prompt_60_EE->Fill(reliso_noPU_prompt_60_EE);
		h_noPU_prompt_80_EE->Fill(reliso_noPU_prompt_80_EE);
		h_noPU_prompt_100_EE->Fill(reliso_noPU_prompt_100_EE);

		// Initialize
        sumiso_noPU_prompt_EE=0, sumiso_noPU_prompt_1sigma_EE=0, sumiso_noPU_prompt_2sigma_EE=0, sumiso_noPU_prompt_3sigma_EE=0, sumiso_noPU_prompt_4sigma_EE=0, sumiso_noPU_prompt_40_EE=0, sumiso_noPU_prompt_60_EE=0, sumiso_noPU_prompt_80_EE=0, sumiso_noPU_prompt_100_EE=0;
        reliso_noPU_prompt_EE=0, reliso_noPU_prompt_1sigma_EE=0, reliso_noPU_prompt_2sigma_EE=0, reliso_noPU_prompt_3sigma_EE=0, reliso_noPU_prompt_4sigma_EE=0, reliso_noPU_prompt_40_EE=0, reliso_noPU_prompt_60_EE=0, reliso_noPU_prompt_80_EE=0, reliso_noPU_prompt_100_EE=0;
	  }
	} // End of muon loop
  } // End of event loop


  /////////////////////////
  //// noPU nonprompt ////
  /////////////////////////
  
  int n_status_failed_noPU_nonprompt_EB=0, n_status_failed_noPU_nonprompt_EE=0;
  int n_muon_noPU_nonprompt_EB=0, n_muon_noPU_nonprompt_EE=0;

  // event loop
  for(int ievt=0; ievt<ch_noPU_nonprompt->GetEntries(); ievt++) {
	ch_noPU_nonprompt->GetEntry(ievt);

	// Make vtx selection tighter
	if(flag_vtx_matching) {
	  if(simvtx_reco_!=vtx_index_) continue; // it includes reco2sim matching
	}
	if(flag_vtx_selectedLV) {
	  if(!(simvtx_bx_==0 && simvtx_evtId_==0 /*&& recovtx_original_index_==0*/)) continue;
	}

	// muon loop
	for(int im=0; im<muon_pt_->size(); im++) {

	  // Muon selection
      if(flag_isMuon) {
        if(muon_isMuon_->at(im)==0) continue;
      }
      if(flag_isPFMuon) {
        if(muon_isPFMuon_->at(im)==0) continue;
      }
      if(flag_isGlobalMuon) {
        if(muon_isGlobalMuon_->at(im)==0) continue;
      }
      if(flag_isTrackerMuon) {
        if(muon_isTrackerMuon_->at(im)==0) continue;
      }
      if(flag_isStandAloneMuon) {
        if(muon_isStandAloneMuon_->at(im)==0) continue;
      }
      if(flag_isCutBasedIdLoose) {
        if(muon_isCutBasedIdLoose_->at(im)==0) continue;
      }
      if(flag_isLooseMuon) {
        if(muon_isLooseMuon_->at(im)==0) continue;
      }


	  /////////////////////////////
	  // nonprompt muon - Barrel //
	  /////////////////////////////
	  if(muon_prompt_->at(im)==0 && muon_isBarrel_->at(im)==1) {
		n_muon_noPU_nonprompt_EB++;

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
		    n_status_failed_noPU_nonprompt_EB++;
		    continue;
	      }
		}
	    if(flag_muon_pv_dz) {
		  if(muon_pv_dz_->at(im)>muon_pv_dz_cut_EB) continue;
	    }
	    if(flag_muon_pv_dxy) {
	      if(muon_pv_dxy_->at(im)>muon_pv_dxy_cut_EB) continue;
	    }
	    if(flag_muon_PVweight) {
		  if(muon_PVweight_->at(im)<muon_PVweight_cut) continue;
	    }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {

		  if(flag_track_pv_dz) {
		    //if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) >= 0.1) continue;
		    if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EB) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  // noMTD case
		  sumiso_noPU_nonprompt_EB += track_pt_->at(im).at(it);
		  // MTD case
		  if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_noPU_nonprompt_1sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_noPU_nonprompt_2sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_noPU_nonprompt_3sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_noPU_nonprompt_4sigma_EB += track_pt_->at(im).at(it);
		  if( dt_cut<0.12  &&   dt_cut>0 ) sumiso_noPU_nonprompt_40_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.18  &&   dt_cut>0 ) sumiso_noPU_nonprompt_60_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.24  &&   dt_cut>0 ) sumiso_noPU_nonprompt_80_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.30  &&   dt_cut>0 ) sumiso_noPU_nonprompt_100_EB    += track_pt_->at(im).at(it);
		  else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
            sumiso_noPU_nonprompt_1sigma_EB += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_2sigma_EB += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_3sigma_EB += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_4sigma_EB += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_40_EB     += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_60_EB     += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_80_EB     += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_100_EB    += track_pt_->at(im).at(it);
		  }
		  dtsig_cut=0;
		  dt_cut=0;
		}
		reliso_noPU_nonprompt_EB        = sumiso_noPU_nonprompt_EB/muon_pt_->at(im);
		reliso_noPU_nonprompt_1sigma_EB = sumiso_noPU_nonprompt_1sigma_EB/muon_pt_->at(im);
		reliso_noPU_nonprompt_2sigma_EB = sumiso_noPU_nonprompt_2sigma_EB/muon_pt_->at(im);
		reliso_noPU_nonprompt_3sigma_EB = sumiso_noPU_nonprompt_3sigma_EB/muon_pt_->at(im);
		reliso_noPU_nonprompt_4sigma_EB = sumiso_noPU_nonprompt_4sigma_EB/muon_pt_->at(im);
		reliso_noPU_nonprompt_40_EB     = sumiso_noPU_nonprompt_40_EB/muon_pt_->at(im);
		reliso_noPU_nonprompt_60_EB     = sumiso_noPU_nonprompt_60_EB/muon_pt_->at(im);
		reliso_noPU_nonprompt_80_EB     = sumiso_noPU_nonprompt_80_EB/muon_pt_->at(im);
		reliso_noPU_nonprompt_100_EB    = sumiso_noPU_nonprompt_100_EB/muon_pt_->at(im);

		// Store value of relative isolation
		h_noPU_nonprompt_EB->Fill(reliso_noPU_nonprompt_EB);
		h_noPU_nonprompt_1sigma_EB->Fill(reliso_noPU_nonprompt_1sigma_EB);
		h_noPU_nonprompt_2sigma_EB->Fill(reliso_noPU_nonprompt_2sigma_EB);
		h_noPU_nonprompt_3sigma_EB->Fill(reliso_noPU_nonprompt_3sigma_EB);
		h_noPU_nonprompt_4sigma_EB->Fill(reliso_noPU_nonprompt_4sigma_EB);
		h_noPU_nonprompt_40_EB->Fill(reliso_noPU_nonprompt_40_EB);
		h_noPU_nonprompt_60_EB->Fill(reliso_noPU_nonprompt_60_EB);
		h_noPU_nonprompt_80_EB->Fill(reliso_noPU_nonprompt_80_EB);
		h_noPU_nonprompt_100_EB->Fill(reliso_noPU_nonprompt_100_EB);

		// Initialize
        sumiso_noPU_nonprompt_EB=0, sumiso_noPU_nonprompt_1sigma_EB=0, sumiso_noPU_nonprompt_2sigma_EB=0, sumiso_noPU_nonprompt_3sigma_EB=0, sumiso_noPU_nonprompt_4sigma_EB=0, sumiso_noPU_nonprompt_40_EB=0, sumiso_noPU_nonprompt_60_EB=0, sumiso_noPU_nonprompt_80_EB=0, sumiso_noPU_nonprompt_100_EB=0;
        reliso_noPU_nonprompt_EB=0, reliso_noPU_nonprompt_1sigma_EB=0, reliso_noPU_nonprompt_2sigma_EB=0, reliso_noPU_nonprompt_3sigma_EB=0, reliso_noPU_nonprompt_4sigma_EB=0, reliso_noPU_nonprompt_40_EB=0, reliso_noPU_nonprompt_60_EB=0, reliso_noPU_nonprompt_80_EB=0, reliso_noPU_nonprompt_100_EB=0;
	  }

	  /////////////////////////////
	  // nonprompt muon - Endcap //
	  /////////////////////////////
	  else if(muon_prompt_->at(im)==0 && muon_isBarrel_->at(im)==0) {
		n_muon_noPU_nonprompt_EE++;

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
		    n_status_failed_noPU_nonprompt_EE++;
		    continue;
	      }
		}
	    if(flag_muon_pv_dz) {
		  if(muon_pv_dz_->at(im)>muon_pv_dz_cut_EE) continue;
	    }
	    if(flag_muon_pv_dxy) {
	      if(muon_pv_dxy_->at(im)>muon_pv_dxy_cut_EE) continue;
	    }
	    if(flag_muon_PVweight) {
		  if(muon_PVweight_->at(im)<muon_PVweight_cut) continue;
	    }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {

		  if(flag_track_pv_dz) {
		    //if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) >= 0.2) continue;
		    if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EE) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  // noMTD case
		  sumiso_noPU_nonprompt_EE += track_pt_->at(im).at(it);
		  // MTD case
		  if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_noPU_nonprompt_1sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_noPU_nonprompt_2sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_noPU_nonprompt_3sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_noPU_nonprompt_4sigma_EE += track_pt_->at(im).at(it);
		  if( dt_cut<0.12  &&   dt_cut>0 ) sumiso_noPU_nonprompt_40_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.18  &&   dt_cut>0 ) sumiso_noPU_nonprompt_60_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.24  &&   dt_cut>0 ) sumiso_noPU_nonprompt_80_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.30  &&   dt_cut>0 ) sumiso_noPU_nonprompt_100_EE    += track_pt_->at(im).at(it);
		  else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
            sumiso_noPU_nonprompt_1sigma_EE += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_2sigma_EE += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_3sigma_EE += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_4sigma_EE += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_40_EE     += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_60_EE     += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_80_EE     += track_pt_->at(im).at(it);
            sumiso_noPU_nonprompt_100_EE    += track_pt_->at(im).at(it);
		  }
		  dtsig_cut=0;
		  dt_cut=0;
		}
		reliso_noPU_nonprompt_EE        = sumiso_noPU_nonprompt_EE/muon_pt_->at(im);
		reliso_noPU_nonprompt_1sigma_EE = sumiso_noPU_nonprompt_1sigma_EE/muon_pt_->at(im);
		reliso_noPU_nonprompt_2sigma_EE = sumiso_noPU_nonprompt_2sigma_EE/muon_pt_->at(im);
		reliso_noPU_nonprompt_3sigma_EE = sumiso_noPU_nonprompt_3sigma_EE/muon_pt_->at(im);
		reliso_noPU_nonprompt_4sigma_EE = sumiso_noPU_nonprompt_4sigma_EE/muon_pt_->at(im);
		reliso_noPU_nonprompt_40_EE     = sumiso_noPU_nonprompt_40_EE/muon_pt_->at(im);
		reliso_noPU_nonprompt_60_EE     = sumiso_noPU_nonprompt_60_EE/muon_pt_->at(im);
		reliso_noPU_nonprompt_80_EE     = sumiso_noPU_nonprompt_80_EE/muon_pt_->at(im);
		reliso_noPU_nonprompt_100_EE    = sumiso_noPU_nonprompt_100_EE/muon_pt_->at(im);

		// Store value of relative isolation
		h_noPU_nonprompt_EE->Fill(reliso_noPU_nonprompt_EE);
		h_noPU_nonprompt_1sigma_EE->Fill(reliso_noPU_nonprompt_1sigma_EE);
		h_noPU_nonprompt_2sigma_EE->Fill(reliso_noPU_nonprompt_2sigma_EE);
		h_noPU_nonprompt_3sigma_EE->Fill(reliso_noPU_nonprompt_3sigma_EE);
		h_noPU_nonprompt_4sigma_EE->Fill(reliso_noPU_nonprompt_4sigma_EE);
		h_noPU_nonprompt_40_EE->Fill(reliso_noPU_nonprompt_40_EE);
		h_noPU_nonprompt_60_EE->Fill(reliso_noPU_nonprompt_60_EE);
		h_noPU_nonprompt_80_EE->Fill(reliso_noPU_nonprompt_80_EE);
		h_noPU_nonprompt_100_EE->Fill(reliso_noPU_nonprompt_100_EE);

		// Initialize
        sumiso_noPU_nonprompt_EE=0, sumiso_noPU_nonprompt_1sigma_EE=0, sumiso_noPU_nonprompt_2sigma_EE=0, sumiso_noPU_nonprompt_3sigma_EE=0, sumiso_noPU_nonprompt_4sigma_EE=0, sumiso_noPU_nonprompt_40_EE=0, sumiso_noPU_nonprompt_60_EE=0, sumiso_noPU_nonprompt_80_EE=0, sumiso_noPU_nonprompt_100_EE=0;
        reliso_noPU_nonprompt_EE=0, reliso_noPU_nonprompt_1sigma_EE=0, reliso_noPU_nonprompt_2sigma_EE=0, reliso_noPU_nonprompt_3sigma_EE=0, reliso_noPU_nonprompt_4sigma_EE=0, reliso_noPU_nonprompt_40_EE=0, reliso_noPU_nonprompt_60_EE=0, reliso_noPU_nonprompt_80_EE=0, reliso_noPU_nonprompt_100_EE=0;

	  }
	} // End of muon loop
  } // End of event loop
	

  cout << "Fraction that status of genParticles of muon is -99 or muon does not have genParticles (Fake)" << endl;
  cout << "[PU200 prompt Barrel]    " << n_status_failed_PU200_prompt_EB      << "/" << n_muon_PU200_prompt_EB << " = "
	  								  << n_status_failed_PU200_prompt_EB*1./n_muon_PU200_prompt_EB << endl;
  cout << "[PU200 prompt Endcap]    " << n_status_failed_PU200_prompt_EE      << "/" << n_muon_PU200_prompt_EE << " = "
	  								  << n_status_failed_PU200_prompt_EE*1./n_muon_PU200_prompt_EE << endl;
  cout << "[PU200 nonprompt Barrel] " << n_status_failed_PU200_nonprompt_EB      << "/" << n_muon_PU200_nonprompt_EB << " = "
	  								  << n_status_failed_PU200_nonprompt_EB*1./n_muon_PU200_nonprompt_EB << endl;
  cout << "[PU200 nonprompt Endcap] " << n_status_failed_PU200_nonprompt_EE      << "/" << n_muon_PU200_nonprompt_EE << " = "
	  								  << n_status_failed_PU200_nonprompt_EE*1./n_muon_PU200_nonprompt_EE << endl;
  cout << "[noPU prompt Barrel]     " << n_status_failed_noPU_prompt_EB      << "/" << n_muon_noPU_prompt_EB << " = "
	  								  << n_status_failed_noPU_prompt_EB*1./n_muon_noPU_prompt_EB << endl;
  cout << "[noPU prompt Endcap]     " << n_status_failed_noPU_prompt_EE      << "/" << n_muon_noPU_prompt_EE << " = "
	  								  << n_status_failed_noPU_prompt_EE*1./n_muon_noPU_prompt_EE << endl;
  cout << "[noPU nonprompt Barrel]  " << n_status_failed_noPU_nonprompt_EB      << "/" << n_muon_noPU_nonprompt_EB << " = "
	  								  << n_status_failed_noPU_nonprompt_EB*1./n_muon_noPU_nonprompt_EB << endl;
  cout << "[noPU nonprompt Endcap]  " << n_status_failed_noPU_nonprompt_EE      << "/" << n_muon_noPU_nonprompt_EE << " = "
	  								  << n_status_failed_noPU_nonprompt_EE*1./n_muon_noPU_nonprompt_EE << endl;

  cout << endl;
  cout << "Number of muons" << endl;
  cout << "[PU200 prompt]    " << n_muon_PU200_prompt_EB-n_status_failed_PU200_prompt_EB+n_muon_PU200_prompt_EE-n_status_failed_PU200_prompt_EE << " = " << n_muon_PU200_prompt_EB-n_status_failed_PU200_prompt_EB << " + " << n_muon_PU200_prompt_EE-n_status_failed_PU200_prompt_EE << endl;
  cout << "[PU200 nonprompt] " << n_muon_PU200_nonprompt_EB-n_status_failed_PU200_nonprompt_EB+n_muon_PU200_nonprompt_EE-n_status_failed_PU200_nonprompt_EE << " = " << n_muon_PU200_nonprompt_EB-n_status_failed_PU200_nonprompt_EB << " + " << n_muon_PU200_nonprompt_EE-n_status_failed_PU200_nonprompt_EE << endl;
  cout << endl;
  cout << "[noPU prompt]     " << n_muon_noPU_prompt_EB-n_status_failed_noPU_prompt_EB+n_muon_noPU_prompt_EE-n_status_failed_noPU_prompt_EE << " = " << n_muon_noPU_prompt_EB-n_status_failed_noPU_prompt_EB << " + " << n_muon_noPU_prompt_EE-n_status_failed_noPU_prompt_EE << endl;
  cout << "[noPU nonprompt]  " << n_muon_noPU_nonprompt_EB-n_status_failed_noPU_nonprompt_EB+n_muon_noPU_nonprompt_EE-n_status_failed_noPU_nonprompt_EE << " = " << n_muon_noPU_nonprompt_EB-n_status_failed_noPU_nonprompt_EB << " + " << n_muon_noPU_nonprompt_EE-n_status_failed_noPU_nonprompt_EE << endl;



/*
  // PU200
    // prompt
//  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_EB",       "muon_prompt_==1 && muon_isBarrel_==1", "goff");
//  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_1sigma_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<1.0", "goff");
//  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_2sigma_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<2.0", "goff");
//  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_3sigma_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<3.0", "goff");
//  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_40_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.12", "goff");
//  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_60_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.18", "goff");
//  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_80_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.24", "goff");
//  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_100_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.30", "goff");
//  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_gen_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && genMatched==1", "goff");
  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_EE",       "muon_prompt_==1 && muon_isBarrel_==0", "goff");
  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_1sigma_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<1.0", "goff");
  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_2sigma_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<2.0", "goff");
  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_3sigma_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<3.0", "goff");
  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_40_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.12", "goff");
  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_60_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.18", "goff");
  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_80_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.24", "goff");
  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_100_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.30", "goff");
//  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_gen_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && genMatched==1", "goff");
    // nonprompt
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_EB",       "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_1sigma_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<1.0", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_2sigma_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<2.0", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_3sigma_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<3.0", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_40_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.12", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_60_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.18", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_80_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.24", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_100_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.30", "goff");
//  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_gen_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && genMatched==1", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_EE",       "muon_prompt_==0 && muon_isBarrel_==0", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_1sigma_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<1.0", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_2sigma_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<2.0", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_3sigma_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<3.0", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_40_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.12", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_60_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.18", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_80_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.24", "goff");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_100_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.30", "goff");
//  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_gen_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && genMatched==1", "goff");
  // noPU
    // prompt
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_EB",       "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_1sigma_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<1.0", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_2sigma_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<2.0", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_3sigma_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<3.0", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_40_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.12", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_60_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.18", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_80_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.24", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_100_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.30", "goff");
//  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_gen_EB",       "muon_prompt_==1 && muon_isBarrel_==1 && genMatched==1", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_EE",       "muon_prompt_==1 && muon_isBarrel_==0", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_1sigma_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<1.0", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_2sigma_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<2.0", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_3sigma_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<3.0", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_40_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.12", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_60_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.18", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_80_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.24", "goff");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_100_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.30", "goff");
//  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_gen_EE",       "muon_prompt_==1 && muon_isBarrel_==0 && genMatched==1", "goff");
    // nonprompt
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_EB",       "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_1sigma_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<1.0", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_2sigma_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<2.0", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_3sigma_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<3.0", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_40_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.12", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_60_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.18", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_80_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.24", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_100_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && (TMath::Abs(track_time_-muon_time_))<0.30", "goff");
//  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_gen_EB",       "muon_prompt_==0 && muon_isBarrel_==1 && genMatched==1", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_EE",       "muon_prompt_==0 && muon_isBarrel_==0", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_1sigma_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<1.0", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_2sigma_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<2.0", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_3sigma_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(muon_time_-track_time_)/TMath::Sqrt(track_time_err_*track_time_err_+muon_time_err_*muon_time_err_))<3.0", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_40_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.12", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_60_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.18", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_80_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.24", "goff");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_100_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && (TMath::Abs(track_time_-muon_time_))<0.30", "goff");
//  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_gen_EE",       "muon_prompt_==0 && muon_isBarrel_==0 && genMatched==1", "goff");


  ////////////////////////////////////////////////////////////////////////
  // Consider the case where there are no tracks within the muon's cone //
  //          Only the case when using timing information               //
  ////////////////////////////////////////////////////////////////////////

  TH1D* h_PU200_prompt_EB_notrack    = new TH1D("h_PU200_prompt_EB_notrack",    "h_PU200_prompt_EB_notrack",    1, 0, 1);
  TH1D* h_PU200_prompt_EE_notrack    = new TH1D("h_PU200_prompt_EE_notrack",    "h_PU200_prompt_EE_notrack",    1, 0, 1);
  TH1D* h_PU200_nonprompt_EB_notrack = new TH1D("h_PU200_nonprompt_EB_notrack", "h_PU200_nonprompt_EB_notrack", 1, 0, 1);
  TH1D* h_PU200_nonprompt_EE_notrack = new TH1D("h_PU200_nonprompt_EE_notrack", "h_PU200_nonprompt_EE_notrack", 1, 0, 1);
  TH1D* h_noPU_prompt_EB_notrack     = new TH1D("h_noPU_prompt_EB_notrack",     "h_noPU_prompt_EB_notrack",     1, 0, 1);
  TH1D* h_noPU_prompt_EE_notrack     = new TH1D("h_noPU_prompt_EE_notrack",     "h_noPU_prompt_EE_notrack",     1, 0, 1);
  TH1D* h_noPU_nonprompt_EB_notrack  = new TH1D("h_noPU_nonprompt_EB_notrack",  "h_noPU_nonprompt_EB_notrack",  1, 0, 1);
  TH1D* h_noPU_nonprompt_EE_notrack  = new TH1D("h_noPU_nonprompt_EE_notrack",  "h_noPU_nonprompt_EE_notrack",  1, 0, 1);

  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_EB_notrack",           "muon_prompt_==1 && muon_isBarrel_==1 && muon_rel_iso_==0");
  ch_PU200_prompt->Draw("muon_rel_iso_>>h_PU200_prompt_EE_notrack",           "muon_prompt_==1 && muon_isBarrel_==0 && muon_rel_iso_==0");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_EB_notrack",     "muon_prompt_==0 && muon_isBarrel_==1 && muon_rel_iso_==0");
  ch_PU200_nonprompt->Draw("muon_rel_iso_>>h_PU200_nonprompt_EE_notrack",     "muon_prompt_==0 && muon_isBarrel_==0 && muon_rel_iso_==0");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_EB_notrack",           "muon_prompt_==1 && muon_isBarrel_==1 && muon_rel_iso_==0");
  ch_noPU_prompt->Draw("muon_rel_iso_>>h_noPU_prompt_EE_notrack",           "muon_prompt_==1 && muon_isBarrel_==0 && muon_rel_iso_==0");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_EB_notrack",     "muon_prompt_==0 && muon_isBarrel_==1 && muon_rel_iso_==0");
  ch_noPU_nonprompt->Draw("muon_rel_iso_>>h_noPU_nonprompt_EE_notrack",     "muon_prompt_==0 && muon_isBarrel_==0 && muon_rel_iso_==0");

  // PU200
    // prompt
//  h_PU200_prompt_1sigma_EB->SetBinContent(1, h_PU200_prompt_1sigma_EB->GetBinContent(1)+h_PU200_prompt_EB_notrack->GetBinContent(1));
//  h_PU200_prompt_2sigma_EB->SetBinContent(1, h_PU200_prompt_2sigma_EB->GetBinContent(1)+h_PU200_prompt_EB_notrack->GetBinContent(1));
//  h_PU200_prompt_3sigma_EB->SetBinContent(1, h_PU200_prompt_3sigma_EB->GetBinContent(1)+h_PU200_prompt_EB_notrack->GetBinContent(1));
//  h_PU200_prompt_40_EB->SetBinContent(1,     h_PU200_prompt_40_EB->GetBinContent(1)+h_PU200_prompt_EB_notrack->GetBinContent(1));
//  h_PU200_prompt_60_EB->SetBinContent(1,     h_PU200_prompt_60_EB->GetBinContent(1)+h_PU200_prompt_EB_notrack->GetBinContent(1));
//  h_PU200_prompt_80_EB->SetBinContent(1,     h_PU200_prompt_80_EB->GetBinContent(1)+h_PU200_prompt_EB_notrack->GetBinContent(1));
//  h_PU200_prompt_100_EB->SetBinContent(1,    h_PU200_prompt_100_EB->GetBinContent(1)+h_PU200_prompt_EB_notrack->GetBinContent(1));
  h_PU200_prompt_1sigma_EE->SetBinContent(1, h_PU200_prompt_1sigma_EE->GetBinContent(1)+h_PU200_prompt_EE_notrack->GetBinContent(1));
  h_PU200_prompt_2sigma_EE->SetBinContent(1, h_PU200_prompt_2sigma_EE->GetBinContent(1)+h_PU200_prompt_EE_notrack->GetBinContent(1));
  h_PU200_prompt_3sigma_EE->SetBinContent(1, h_PU200_prompt_3sigma_EE->GetBinContent(1)+h_PU200_prompt_EE_notrack->GetBinContent(1));
  h_PU200_prompt_40_EE->SetBinContent(1,     h_PU200_prompt_40_EE->GetBinContent(1)+h_PU200_prompt_EE_notrack->GetBinContent(1));
  h_PU200_prompt_60_EE->SetBinContent(1,     h_PU200_prompt_60_EE->GetBinContent(1)+h_PU200_prompt_EE_notrack->GetBinContent(1));
  h_PU200_prompt_80_EE->SetBinContent(1,     h_PU200_prompt_80_EE->GetBinContent(1)+h_PU200_prompt_EE_notrack->GetBinContent(1));
  h_PU200_prompt_100_EE->SetBinContent(1,    h_PU200_prompt_100_EE->GetBinContent(1)+h_PU200_prompt_EE_notrack->GetBinContent(1));

    // nonprompt
  h_PU200_nonprompt_1sigma_EB->SetBinContent(1, h_PU200_nonprompt_1sigma_EB->GetBinContent(1)+h_PU200_nonprompt_EB_notrack->GetBinContent(1));
  h_PU200_nonprompt_2sigma_EB->SetBinContent(1, h_PU200_nonprompt_2sigma_EB->GetBinContent(1)+h_PU200_nonprompt_EB_notrack->GetBinContent(1));
  h_PU200_nonprompt_3sigma_EB->SetBinContent(1, h_PU200_nonprompt_3sigma_EB->GetBinContent(1)+h_PU200_nonprompt_EB_notrack->GetBinContent(1));
  h_PU200_nonprompt_40_EB->SetBinContent(1,     h_PU200_nonprompt_40_EB->GetBinContent(1)+h_PU200_nonprompt_EB_notrack->GetBinContent(1));
  h_PU200_nonprompt_60_EB->SetBinContent(1,     h_PU200_nonprompt_60_EB->GetBinContent(1)+h_PU200_nonprompt_EB_notrack->GetBinContent(1));
  h_PU200_nonprompt_80_EB->SetBinContent(1,     h_PU200_nonprompt_80_EB->GetBinContent(1)+h_PU200_nonprompt_EB_notrack->GetBinContent(1));
  h_PU200_nonprompt_100_EB->SetBinContent(1,    h_PU200_nonprompt_100_EB->GetBinContent(1)+h_PU200_nonprompt_EB_notrack->GetBinContent(1));
  h_PU200_nonprompt_1sigma_EE->SetBinContent(1, h_PU200_nonprompt_1sigma_EE->GetBinContent(1)+h_PU200_nonprompt_EE_notrack->GetBinContent(1));
  h_PU200_nonprompt_2sigma_EE->SetBinContent(1, h_PU200_nonprompt_2sigma_EE->GetBinContent(1)+h_PU200_nonprompt_EE_notrack->GetBinContent(1));
  h_PU200_nonprompt_3sigma_EE->SetBinContent(1, h_PU200_nonprompt_3sigma_EE->GetBinContent(1)+h_PU200_nonprompt_EE_notrack->GetBinContent(1));
  h_PU200_nonprompt_40_EE->SetBinContent(1,     h_PU200_nonprompt_40_EE->GetBinContent(1)+h_PU200_nonprompt_EE_notrack->GetBinContent(1));
  h_PU200_nonprompt_60_EE->SetBinContent(1,     h_PU200_nonprompt_60_EE->GetBinContent(1)+h_PU200_nonprompt_EE_notrack->GetBinContent(1));
  h_PU200_nonprompt_80_EE->SetBinContent(1,     h_PU200_nonprompt_80_EE->GetBinContent(1)+h_PU200_nonprompt_EE_notrack->GetBinContent(1));
  h_PU200_nonprompt_100_EE->SetBinContent(1,    h_PU200_nonprompt_100_EE->GetBinContent(1)+h_PU200_nonprompt_EE_notrack->GetBinContent(1));

*/

  ////////////////////////////////////////////
  // Define vectors to store values of eff. //
  ////////////////////////////////////////////

  // PU200
  vector<double> prompt_eff_PU200_EB,    prompt_eff_PU200_1sigma_EB,    prompt_eff_PU200_2sigma_EB,    prompt_eff_PU200_3sigma_EB,    prompt_eff_PU200_4sigma_EB,    prompt_eff_PU200_40_EB,    prompt_eff_PU200_60_EB,    prompt_eff_PU200_80_EB,    prompt_eff_PU200_100_EB,    prompt_eff_PU200_gen_EB;
  vector<double> prompt_eff_PU200_EE,    prompt_eff_PU200_1sigma_EE,    prompt_eff_PU200_2sigma_EE,    prompt_eff_PU200_3sigma_EE,    prompt_eff_PU200_4sigma_EE,    prompt_eff_PU200_40_EE,    prompt_eff_PU200_60_EE,    prompt_eff_PU200_80_EE,    prompt_eff_PU200_100_EE,    prompt_eff_PU200_gen_EE;
  vector<double> nonprompt_eff_PU200_EB, nonprompt_eff_PU200_1sigma_EB, nonprompt_eff_PU200_2sigma_EB, nonprompt_eff_PU200_3sigma_EB, nonprompt_eff_PU200_4sigma_EB, nonprompt_eff_PU200_40_EB, nonprompt_eff_PU200_60_EB, nonprompt_eff_PU200_80_EB, nonprompt_eff_PU200_100_EB, nonprompt_eff_PU200_gen_EB;
  vector<double> nonprompt_eff_PU200_EE, nonprompt_eff_PU200_1sigma_EE, nonprompt_eff_PU200_2sigma_EE, nonprompt_eff_PU200_3sigma_EE, nonprompt_eff_PU200_4sigma_EE, nonprompt_eff_PU200_40_EE, nonprompt_eff_PU200_60_EE, nonprompt_eff_PU200_80_EE, nonprompt_eff_PU200_100_EE, nonprompt_eff_PU200_gen_EE;
  // noPU
  vector<double> prompt_eff_noPU_EB,    prompt_eff_noPU_1sigma_EB,    prompt_eff_noPU_2sigma_EB,    prompt_eff_noPU_3sigma_EB,    prompt_eff_noPU_4sigma_EB,    prompt_eff_noPU_40_EB,    prompt_eff_noPU_60_EB,    prompt_eff_noPU_80_EB,    prompt_eff_noPU_100_EB,    prompt_eff_noPU_gen_EB;
  vector<double> prompt_eff_noPU_EE,    prompt_eff_noPU_1sigma_EE,    prompt_eff_noPU_2sigma_EE,    prompt_eff_noPU_3sigma_EE,    prompt_eff_noPU_4sigma_EE,    prompt_eff_noPU_40_EE,    prompt_eff_noPU_60_EE,    prompt_eff_noPU_80_EE,    prompt_eff_noPU_100_EE,    prompt_eff_noPU_gen_EE;
  vector<double> nonprompt_eff_noPU_EB, nonprompt_eff_noPU_1sigma_EB, nonprompt_eff_noPU_2sigma_EB, nonprompt_eff_noPU_3sigma_EB, nonprompt_eff_noPU_4sigma_EB, nonprompt_eff_noPU_40_EB, nonprompt_eff_noPU_60_EB, nonprompt_eff_noPU_80_EB, nonprompt_eff_noPU_100_EB, nonprompt_eff_noPU_gen_EB;
  vector<double> nonprompt_eff_noPU_EE, nonprompt_eff_noPU_1sigma_EE, nonprompt_eff_noPU_2sigma_EE, nonprompt_eff_noPU_3sigma_EE, nonprompt_eff_noPU_4sigma_EE, nonprompt_eff_noPU_40_EE, nonprompt_eff_noPU_60_EE, nonprompt_eff_noPU_80_EE, nonprompt_eff_noPU_100_EE, nonprompt_eff_noPU_gen_EE;

  for(int i=0; i<nbin+1; i++) {
	// PU200
	  // prompt
	prompt_eff_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(1,i+1)/h_PU200_prompt_EB->Integral(1,nbin+1));
	prompt_eff_PU200_1sigma_EB.emplace_back(h_PU200_prompt_1sigma_EB->Integral(1,i+1)/h_PU200_prompt_1sigma_EB->Integral(1,nbin+1));
	prompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(1,i+1)/h_PU200_prompt_2sigma_EB->Integral(1,nbin+1));
	prompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(1,i+1)/h_PU200_prompt_3sigma_EB->Integral(1,nbin+1));
	prompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(1,i+1)/h_PU200_prompt_4sigma_EB->Integral(1,nbin+1));
	prompt_eff_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(1,i+1)/h_PU200_prompt_40_EB->Integral(1,nbin+1));
	prompt_eff_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(1,i+1)/h_PU200_prompt_60_EB->Integral(1,nbin+1));
	prompt_eff_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(1,i+1)/h_PU200_prompt_80_EB->Integral(1,nbin+1));
	prompt_eff_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(1,i+1)/h_PU200_prompt_100_EB->Integral(1,nbin+1));
//	prompt_eff_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(1,i+1)/h_PU200_prompt_gen_EB->Integral(1,nbin+1));
	prompt_eff_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(1,i+1)/h_PU200_prompt_EE->Integral(1,nbin+1));
	prompt_eff_PU200_1sigma_EE.emplace_back(h_PU200_prompt_1sigma_EE->Integral(1,i+1)/h_PU200_prompt_1sigma_EE->Integral(1,nbin+1));
	prompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(1,i+1)/h_PU200_prompt_2sigma_EE->Integral(1,nbin+1));
	prompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(1,i+1)/h_PU200_prompt_3sigma_EE->Integral(1,nbin+1));
	prompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(1,i+1)/h_PU200_prompt_4sigma_EE->Integral(1,nbin+1));
	prompt_eff_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(1,i+1)/h_PU200_prompt_40_EE->Integral(1,nbin+1));
	prompt_eff_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(1,i+1)/h_PU200_prompt_60_EE->Integral(1,nbin+1));
	prompt_eff_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(1,i+1)/h_PU200_prompt_80_EE->Integral(1,nbin+1));
	prompt_eff_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(1,i+1)/h_PU200_prompt_100_EE->Integral(1,nbin+1));
//	prompt_eff_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(1,i+1)/h_PU200_prompt_gen_EE->Integral(1,nbin+1));
	  // nonprompt
	nonprompt_eff_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(1,i+1)/h_PU200_nonprompt_EB->Integral(1,nbin+1));
	nonprompt_eff_PU200_1sigma_EB.emplace_back(h_PU200_nonprompt_1sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_1sigma_EB->Integral(1,nbin+1));
	nonprompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EB->Integral(1,nbin+1));
	nonprompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EB->Integral(1,nbin+1));
	nonprompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EB->Integral(1,nbin+1));
	nonprompt_eff_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(1,i+1)/h_PU200_nonprompt_40_EB->Integral(1,nbin+1));
	nonprompt_eff_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(1,i+1)/h_PU200_nonprompt_60_EB->Integral(1,nbin+1));
	nonprompt_eff_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(1,i+1)/h_PU200_nonprompt_80_EB->Integral(1,nbin+1));
	nonprompt_eff_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(1,i+1)/h_PU200_nonprompt_100_EB->Integral(1,nbin+1));
//	nonprompt_eff_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(1,i+1)/h_PU200_nonprompt_gen_EB->Integral(1,nbin+1));
	nonprompt_eff_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(1,i+1)/h_PU200_nonprompt_EE->Integral(1,nbin+1));
	nonprompt_eff_PU200_1sigma_EE.emplace_back(h_PU200_nonprompt_1sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_1sigma_EE->Integral(1,nbin+1));
	nonprompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EE->Integral(1,nbin+1));
	nonprompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EE->Integral(1,nbin+1));
	nonprompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EE->Integral(1,nbin+1));
	nonprompt_eff_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(1,i+1)/h_PU200_nonprompt_40_EE->Integral(1,nbin+1));
	nonprompt_eff_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(1,i+1)/h_PU200_nonprompt_60_EE->Integral(1,nbin+1));
	nonprompt_eff_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(1,i+1)/h_PU200_nonprompt_80_EE->Integral(1,nbin+1));
	nonprompt_eff_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(1,i+1)/h_PU200_nonprompt_100_EE->Integral(1,nbin+1));
//	nonprompt_eff_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(1,i+1)/h_PU200_nonprompt_gen_EE->Integral(1,nbin+1));
	// noPU
	  // prompt
	prompt_eff_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(1,i+1)/h_noPU_prompt_EB->Integral(1,nbin+1));
	prompt_eff_noPU_1sigma_EB.emplace_back(h_noPU_prompt_1sigma_EB->Integral(1,i+1)/h_noPU_prompt_1sigma_EB->Integral(1,nbin+1));
	prompt_eff_noPU_2sigma_EB.emplace_back(h_noPU_prompt_2sigma_EB->Integral(1,i+1)/h_noPU_prompt_2sigma_EB->Integral(1,nbin+1));
	prompt_eff_noPU_3sigma_EB.emplace_back(h_noPU_prompt_3sigma_EB->Integral(1,i+1)/h_noPU_prompt_3sigma_EB->Integral(1,nbin+1));
	prompt_eff_noPU_4sigma_EB.emplace_back(h_noPU_prompt_4sigma_EB->Integral(1,i+1)/h_noPU_prompt_4sigma_EB->Integral(1,nbin+1));
	prompt_eff_noPU_40_EB.emplace_back(h_noPU_prompt_40_EB->Integral(1,i+1)/h_noPU_prompt_40_EB->Integral(1,nbin+1));
	prompt_eff_noPU_60_EB.emplace_back(h_noPU_prompt_60_EB->Integral(1,i+1)/h_noPU_prompt_60_EB->Integral(1,nbin+1));
	prompt_eff_noPU_80_EB.emplace_back(h_noPU_prompt_80_EB->Integral(1,i+1)/h_noPU_prompt_80_EB->Integral(1,nbin+1));
	prompt_eff_noPU_100_EB.emplace_back(h_noPU_prompt_100_EB->Integral(1,i+1)/h_noPU_prompt_100_EB->Integral(1,nbin+1));
//	prompt_eff_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(1,i+1)/h_noPU_prompt_gen_EB->Integral(1,nbin+1));
	prompt_eff_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(1,i+1)/h_noPU_prompt_EE->Integral(1,nbin+1));
	prompt_eff_noPU_1sigma_EE.emplace_back(h_noPU_prompt_1sigma_EE->Integral(1,i+1)/h_noPU_prompt_1sigma_EE->Integral(1,nbin+1));
	prompt_eff_noPU_2sigma_EE.emplace_back(h_noPU_prompt_2sigma_EE->Integral(1,i+1)/h_noPU_prompt_2sigma_EE->Integral(1,nbin+1));
	prompt_eff_noPU_3sigma_EE.emplace_back(h_noPU_prompt_3sigma_EE->Integral(1,i+1)/h_noPU_prompt_3sigma_EE->Integral(1,nbin+1));
	prompt_eff_noPU_4sigma_EE.emplace_back(h_noPU_prompt_4sigma_EE->Integral(1,i+1)/h_noPU_prompt_4sigma_EE->Integral(1,nbin+1));
	prompt_eff_noPU_40_EE.emplace_back(h_noPU_prompt_40_EE->Integral(1,i+1)/h_noPU_prompt_40_EE->Integral(1,nbin+1));
	prompt_eff_noPU_60_EE.emplace_back(h_noPU_prompt_60_EE->Integral(1,i+1)/h_noPU_prompt_60_EE->Integral(1,nbin+1));
	prompt_eff_noPU_80_EE.emplace_back(h_noPU_prompt_80_EE->Integral(1,i+1)/h_noPU_prompt_80_EE->Integral(1,nbin+1));
	prompt_eff_noPU_100_EE.emplace_back(h_noPU_prompt_100_EE->Integral(1,i+1)/h_noPU_prompt_100_EE->Integral(1,nbin+1));
//	prompt_eff_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(1,i+1)/h_noPU_prompt_gen_EE->Integral(1,nbin+1));
	  // nonprompt
	nonprompt_eff_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(1,i+1)/h_noPU_nonprompt_EB->Integral(1,nbin+1));
	nonprompt_eff_noPU_1sigma_EB.emplace_back(h_noPU_nonprompt_1sigma_EB->Integral(1,i+1)/h_noPU_nonprompt_1sigma_EB->Integral(1,nbin+1));
	nonprompt_eff_noPU_2sigma_EB.emplace_back(h_noPU_nonprompt_2sigma_EB->Integral(1,i+1)/h_noPU_nonprompt_2sigma_EB->Integral(1,nbin+1));
	nonprompt_eff_noPU_3sigma_EB.emplace_back(h_noPU_nonprompt_3sigma_EB->Integral(1,i+1)/h_noPU_nonprompt_3sigma_EB->Integral(1,nbin+1));
	nonprompt_eff_noPU_4sigma_EB.emplace_back(h_noPU_nonprompt_4sigma_EB->Integral(1,i+1)/h_noPU_nonprompt_4sigma_EB->Integral(1,nbin+1));
	nonprompt_eff_noPU_40_EB.emplace_back(h_noPU_nonprompt_40_EB->Integral(1,i+1)/h_noPU_nonprompt_40_EB->Integral(1,nbin+1));
	nonprompt_eff_noPU_60_EB.emplace_back(h_noPU_nonprompt_60_EB->Integral(1,i+1)/h_noPU_nonprompt_60_EB->Integral(1,nbin+1));
	nonprompt_eff_noPU_80_EB.emplace_back(h_noPU_nonprompt_80_EB->Integral(1,i+1)/h_noPU_nonprompt_80_EB->Integral(1,nbin+1));
	nonprompt_eff_noPU_100_EB.emplace_back(h_noPU_nonprompt_100_EB->Integral(1,i+1)/h_noPU_nonprompt_100_EB->Integral(1,nbin+1));
//	nonprompt_eff_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(1,i+1)/h_noPU_nonprompt_gen_EB->Integral(1,nbin+1));
	nonprompt_eff_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(1,i+1)/h_noPU_nonprompt_EE->Integral(1,nbin+1));
	nonprompt_eff_noPU_1sigma_EE.emplace_back(h_noPU_nonprompt_1sigma_EE->Integral(1,i+1)/h_noPU_nonprompt_1sigma_EE->Integral(1,nbin+1));
	nonprompt_eff_noPU_2sigma_EE.emplace_back(h_noPU_nonprompt_2sigma_EE->Integral(1,i+1)/h_noPU_nonprompt_2sigma_EE->Integral(1,nbin+1));
	nonprompt_eff_noPU_3sigma_EE.emplace_back(h_noPU_nonprompt_3sigma_EE->Integral(1,i+1)/h_noPU_nonprompt_3sigma_EE->Integral(1,nbin+1));
	nonprompt_eff_noPU_4sigma_EE.emplace_back(h_noPU_nonprompt_4sigma_EE->Integral(1,i+1)/h_noPU_nonprompt_4sigma_EE->Integral(1,nbin+1));
	nonprompt_eff_noPU_40_EE.emplace_back(h_noPU_nonprompt_40_EE->Integral(1,i+1)/h_noPU_nonprompt_40_EE->Integral(1,nbin+1));
	nonprompt_eff_noPU_60_EE.emplace_back(h_noPU_nonprompt_60_EE->Integral(1,i+1)/h_noPU_nonprompt_60_EE->Integral(1,nbin+1));
	nonprompt_eff_noPU_80_EE.emplace_back(h_noPU_nonprompt_80_EE->Integral(1,i+1)/h_noPU_nonprompt_80_EE->Integral(1,nbin+1));
	nonprompt_eff_noPU_100_EE.emplace_back(h_noPU_nonprompt_100_EE->Integral(1,i+1)/h_noPU_nonprompt_100_EE->Integral(1,nbin+1));
//	nonprompt_eff_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(1,i+1)/h_noPU_nonprompt_gen_EE->Integral(1,nbin+1));


  }

  // Define TGraph
    // PU200
	  // prompt
	TGraph* gr_eff_PU200_prompt_EB = new TGraph();                        TGraph* gr_eff_PU200_prompt_EE = new TGraph();
	TGraph* gr_eff_PU200_prompt_1sigma_EB = new TGraph();                 TGraph* gr_eff_PU200_prompt_1sigma_EE = new TGraph();
	TGraph* gr_eff_PU200_prompt_2sigma_EB = new TGraph();                 TGraph* gr_eff_PU200_prompt_2sigma_EE = new TGraph();
	TGraph* gr_eff_PU200_prompt_3sigma_EB = new TGraph();                 TGraph* gr_eff_PU200_prompt_3sigma_EE = new TGraph();
	TGraph* gr_eff_PU200_prompt_4sigma_EB = new TGraph();                 TGraph* gr_eff_PU200_prompt_4sigma_EE = new TGraph();
	TGraph* gr_eff_PU200_prompt_40_EB = new TGraph();                     TGraph* gr_eff_PU200_prompt_40_EE = new TGraph();
	TGraph* gr_eff_PU200_prompt_60_EB = new TGraph();                     TGraph* gr_eff_PU200_prompt_60_EE = new TGraph();
	TGraph* gr_eff_PU200_prompt_80_EB = new TGraph();                     TGraph* gr_eff_PU200_prompt_80_EE = new TGraph();
	TGraph* gr_eff_PU200_prompt_100_EB = new TGraph();                    TGraph* gr_eff_PU200_prompt_100_EE = new TGraph();
	TGraph* gr_eff_PU200_prompt_gen_EB = new TGraph();                    TGraph* gr_eff_PU200_prompt_gen_EE = new TGraph();
	  // nonprompt
	TGraph* gr_eff_PU200_nonprompt_EB = new TGraph();                        TGraph* gr_eff_PU200_nonprompt_EE = new TGraph();
	TGraph* gr_eff_PU200_nonprompt_1sigma_EB = new TGraph();                 TGraph* gr_eff_PU200_nonprompt_1sigma_EE = new TGraph();
	TGraph* gr_eff_PU200_nonprompt_2sigma_EB = new TGraph();                 TGraph* gr_eff_PU200_nonprompt_2sigma_EE = new TGraph();
	TGraph* gr_eff_PU200_nonprompt_3sigma_EB = new TGraph();                 TGraph* gr_eff_PU200_nonprompt_3sigma_EE = new TGraph();
	TGraph* gr_eff_PU200_nonprompt_4sigma_EB = new TGraph();                 TGraph* gr_eff_PU200_nonprompt_4sigma_EE = new TGraph();
	TGraph* gr_eff_PU200_nonprompt_40_EB = new TGraph();                     TGraph* gr_eff_PU200_nonprompt_40_EE = new TGraph();
	TGraph* gr_eff_PU200_nonprompt_60_EB = new TGraph();                     TGraph* gr_eff_PU200_nonprompt_60_EE = new TGraph();
	TGraph* gr_eff_PU200_nonprompt_80_EB = new TGraph();                     TGraph* gr_eff_PU200_nonprompt_80_EE = new TGraph();
	TGraph* gr_eff_PU200_nonprompt_100_EB = new TGraph();                    TGraph* gr_eff_PU200_nonprompt_100_EE = new TGraph();
	TGraph* gr_eff_PU200_nonprompt_gen_EB = new TGraph();                    TGraph* gr_eff_PU200_nonprompt_gen_EE = new TGraph();
    // noPU
	  // prompt
	TGraph* gr_eff_noPU_prompt_EB = new TGraph();                        TGraph* gr_eff_noPU_prompt_EE = new TGraph();
	TGraph* gr_eff_noPU_prompt_1sigma_EB = new TGraph();                 TGraph* gr_eff_noPU_prompt_1sigma_EE = new TGraph();
	TGraph* gr_eff_noPU_prompt_2sigma_EB = new TGraph();                 TGraph* gr_eff_noPU_prompt_2sigma_EE = new TGraph();
	TGraph* gr_eff_noPU_prompt_3sigma_EB = new TGraph();                 TGraph* gr_eff_noPU_prompt_3sigma_EE = new TGraph();
	TGraph* gr_eff_noPU_prompt_4sigma_EB = new TGraph();                 TGraph* gr_eff_noPU_prompt_4sigma_EE = new TGraph();
	TGraph* gr_eff_noPU_prompt_40_EB = new TGraph();                     TGraph* gr_eff_noPU_prompt_40_EE = new TGraph();
	TGraph* gr_eff_noPU_prompt_60_EB = new TGraph();                     TGraph* gr_eff_noPU_prompt_60_EE = new TGraph();
	TGraph* gr_eff_noPU_prompt_80_EB = new TGraph();                     TGraph* gr_eff_noPU_prompt_80_EE = new TGraph();
	TGraph* gr_eff_noPU_prompt_100_EB = new TGraph();                    TGraph* gr_eff_noPU_prompt_100_EE = new TGraph();
	TGraph* gr_eff_noPU_prompt_gen_EB = new TGraph();                    TGraph* gr_eff_noPU_prompt_gen_EE = new TGraph();
	  // nonprompt
	TGraph* gr_eff_noPU_nonprompt_EB = new TGraph();                        TGraph* gr_eff_noPU_nonprompt_EE = new TGraph();
	TGraph* gr_eff_noPU_nonprompt_1sigma_EB = new TGraph();                 TGraph* gr_eff_noPU_nonprompt_1sigma_EE = new TGraph();
	TGraph* gr_eff_noPU_nonprompt_2sigma_EB = new TGraph();                 TGraph* gr_eff_noPU_nonprompt_2sigma_EE = new TGraph();
	TGraph* gr_eff_noPU_nonprompt_3sigma_EB = new TGraph();                 TGraph* gr_eff_noPU_nonprompt_3sigma_EE = new TGraph();
	TGraph* gr_eff_noPU_nonprompt_4sigma_EB = new TGraph();                 TGraph* gr_eff_noPU_nonprompt_4sigma_EE = new TGraph();
	TGraph* gr_eff_noPU_nonprompt_40_EB = new TGraph();                     TGraph* gr_eff_noPU_nonprompt_40_EE = new TGraph();
	TGraph* gr_eff_noPU_nonprompt_60_EB = new TGraph();                     TGraph* gr_eff_noPU_nonprompt_60_EE = new TGraph();
	TGraph* gr_eff_noPU_nonprompt_80_EB = new TGraph();                     TGraph* gr_eff_noPU_nonprompt_80_EE = new TGraph();
	TGraph* gr_eff_noPU_nonprompt_100_EB = new TGraph();                    TGraph* gr_eff_noPU_nonprompt_100_EE = new TGraph();
	TGraph* gr_eff_noPU_nonprompt_gen_EB = new TGraph();                    TGraph* gr_eff_noPU_nonprompt_gen_EE = new TGraph();

  for(unsigned int i=0; i<nbin+1; i++) {
	// PU200
	  // prompt
	gr_eff_PU200_prompt_EB->SetPoint(gr_eff_PU200_prompt_EB->GetN(), 0.004*i, prompt_eff_PU200_EB.at(i));
	gr_eff_PU200_prompt_1sigma_EB->SetPoint(gr_eff_PU200_prompt_1sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_1sigma_EB.at(i));
	gr_eff_PU200_prompt_2sigma_EB->SetPoint(gr_eff_PU200_prompt_2sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EB.at(i));
	gr_eff_PU200_prompt_3sigma_EB->SetPoint(gr_eff_PU200_prompt_3sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EB.at(i));
	gr_eff_PU200_prompt_4sigma_EB->SetPoint(gr_eff_PU200_prompt_4sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EB.at(i));
	gr_eff_PU200_prompt_40_EB->SetPoint(gr_eff_PU200_prompt_40_EB->GetN(), 0.004*i, prompt_eff_PU200_40_EB.at(i));
	gr_eff_PU200_prompt_60_EB->SetPoint(gr_eff_PU200_prompt_60_EB->GetN(), 0.004*i, prompt_eff_PU200_60_EB.at(i));
	gr_eff_PU200_prompt_80_EB->SetPoint(gr_eff_PU200_prompt_80_EB->GetN(), 0.004*i, prompt_eff_PU200_80_EB.at(i));
	gr_eff_PU200_prompt_100_EB->SetPoint(gr_eff_PU200_prompt_100_EB->GetN(), 0.004*i, prompt_eff_PU200_100_EB.at(i));
//	gr_eff_PU200_prompt_gen_EB->SetPoint(gr_eff_PU200_prompt_gen_EB->GetN(), 0.004*i, prompt_eff_PU200_gen_EB.at(i));
	gr_eff_PU200_prompt_EE->SetPoint(gr_eff_PU200_prompt_EE->GetN(), 0.004*i, prompt_eff_PU200_EE.at(i));
	gr_eff_PU200_prompt_1sigma_EE->SetPoint(gr_eff_PU200_prompt_1sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_1sigma_EE.at(i));
	gr_eff_PU200_prompt_2sigma_EE->SetPoint(gr_eff_PU200_prompt_2sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EE.at(i));
	gr_eff_PU200_prompt_3sigma_EE->SetPoint(gr_eff_PU200_prompt_3sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EE.at(i));
	gr_eff_PU200_prompt_4sigma_EE->SetPoint(gr_eff_PU200_prompt_4sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EE.at(i));
	gr_eff_PU200_prompt_40_EE->SetPoint(gr_eff_PU200_prompt_40_EE->GetN(), 0.004*i, prompt_eff_PU200_40_EE.at(i));
	gr_eff_PU200_prompt_60_EE->SetPoint(gr_eff_PU200_prompt_60_EE->GetN(), 0.004*i, prompt_eff_PU200_60_EE.at(i));
	gr_eff_PU200_prompt_80_EE->SetPoint(gr_eff_PU200_prompt_80_EE->GetN(), 0.004*i, prompt_eff_PU200_80_EE.at(i));
	gr_eff_PU200_prompt_100_EE->SetPoint(gr_eff_PU200_prompt_100_EE->GetN(), 0.004*i, prompt_eff_PU200_100_EE.at(i));
//	gr_eff_PU200_prompt_gen_EE->SetPoint(gr_eff_PU200_prompt_gen_EE->GetN(), 0.004*i, prompt_eff_PU200_gen_EE.at(i));
      // nonprompt
	gr_eff_PU200_nonprompt_EB->SetPoint(gr_eff_PU200_nonprompt_EB->GetN(), 0.004*i, nonprompt_eff_PU200_EB.at(i));
	gr_eff_PU200_nonprompt_1sigma_EB->SetPoint(gr_eff_PU200_nonprompt_1sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_1sigma_EB.at(i));
	gr_eff_PU200_nonprompt_2sigma_EB->SetPoint(gr_eff_PU200_nonprompt_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EB.at(i));
	gr_eff_PU200_nonprompt_3sigma_EB->SetPoint(gr_eff_PU200_nonprompt_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EB.at(i));
	gr_eff_PU200_nonprompt_4sigma_EB->SetPoint(gr_eff_PU200_nonprompt_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EB.at(i));
	gr_eff_PU200_nonprompt_40_EB->SetPoint(gr_eff_PU200_nonprompt_40_EB->GetN(), 0.004*i, nonprompt_eff_PU200_40_EB.at(i));
	gr_eff_PU200_nonprompt_60_EB->SetPoint(gr_eff_PU200_nonprompt_60_EB->GetN(), 0.004*i, nonprompt_eff_PU200_60_EB.at(i));
	gr_eff_PU200_nonprompt_80_EB->SetPoint(gr_eff_PU200_nonprompt_80_EB->GetN(), 0.004*i, nonprompt_eff_PU200_80_EB.at(i));
	gr_eff_PU200_nonprompt_100_EB->SetPoint(gr_eff_PU200_nonprompt_100_EB->GetN(), 0.004*i, nonprompt_eff_PU200_100_EB.at(i));
//	gr_eff_PU200_nonprompt_gen_EB->SetPoint(gr_eff_PU200_nonprompt_gen_EB->GetN(), 0.004*i, nonprompt_eff_PU200_gen_EB.at(i));
	gr_eff_PU200_nonprompt_EE->SetPoint(gr_eff_PU200_nonprompt_EE->GetN(), 0.004*i, nonprompt_eff_PU200_EE.at(i));
	gr_eff_PU200_nonprompt_1sigma_EE->SetPoint(gr_eff_PU200_nonprompt_1sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_1sigma_EE.at(i));
	gr_eff_PU200_nonprompt_2sigma_EE->SetPoint(gr_eff_PU200_nonprompt_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EE.at(i));
	gr_eff_PU200_nonprompt_3sigma_EE->SetPoint(gr_eff_PU200_nonprompt_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EE.at(i));
	gr_eff_PU200_nonprompt_4sigma_EE->SetPoint(gr_eff_PU200_nonprompt_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EE.at(i));
	gr_eff_PU200_nonprompt_40_EE->SetPoint(gr_eff_PU200_nonprompt_40_EE->GetN(), 0.004*i, nonprompt_eff_PU200_40_EE.at(i));
	gr_eff_PU200_nonprompt_60_EE->SetPoint(gr_eff_PU200_nonprompt_60_EE->GetN(), 0.004*i, nonprompt_eff_PU200_60_EE.at(i));
	gr_eff_PU200_nonprompt_80_EE->SetPoint(gr_eff_PU200_nonprompt_80_EE->GetN(), 0.004*i, nonprompt_eff_PU200_80_EE.at(i));
	gr_eff_PU200_nonprompt_100_EE->SetPoint(gr_eff_PU200_nonprompt_100_EE->GetN(), 0.004*i, nonprompt_eff_PU200_100_EE.at(i));
//	gr_eff_PU200_nonprompt_gen_EE->SetPoint(gr_eff_PU200_nonprompt_gen_EE->GetN(), 0.004*i, nonprompt_eff_PU200_gen_EE.at(i));

	// noPU
	  // prompt
	gr_eff_noPU_prompt_EB->SetPoint(gr_eff_noPU_prompt_EB->GetN(), 0.004*i, prompt_eff_noPU_EB.at(i));
	gr_eff_noPU_prompt_1sigma_EB->SetPoint(gr_eff_noPU_prompt_1sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_1sigma_EB.at(i));
	gr_eff_noPU_prompt_2sigma_EB->SetPoint(gr_eff_noPU_prompt_2sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_2sigma_EB.at(i));
	gr_eff_noPU_prompt_3sigma_EB->SetPoint(gr_eff_noPU_prompt_3sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_3sigma_EB.at(i));
	gr_eff_noPU_prompt_4sigma_EB->SetPoint(gr_eff_noPU_prompt_4sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_4sigma_EB.at(i));
	gr_eff_noPU_prompt_40_EB->SetPoint(gr_eff_noPU_prompt_40_EB->GetN(), 0.004*i, prompt_eff_noPU_40_EB.at(i));
	gr_eff_noPU_prompt_60_EB->SetPoint(gr_eff_noPU_prompt_60_EB->GetN(), 0.004*i, prompt_eff_noPU_60_EB.at(i));
	gr_eff_noPU_prompt_80_EB->SetPoint(gr_eff_noPU_prompt_80_EB->GetN(), 0.004*i, prompt_eff_noPU_80_EB.at(i));
	gr_eff_noPU_prompt_100_EB->SetPoint(gr_eff_noPU_prompt_100_EB->GetN(), 0.004*i, prompt_eff_noPU_100_EB.at(i));
//	gr_eff_noPU_prompt_gen_EB->SetPoint(gr_eff_noPU_prompt_gen_EB->GetN(), 0.004*i, prompt_eff_noPU_gen_EB.at(i));
	gr_eff_noPU_prompt_EE->SetPoint(gr_eff_noPU_prompt_EE->GetN(), 0.004*i, prompt_eff_noPU_EE.at(i));
	gr_eff_noPU_prompt_1sigma_EE->SetPoint(gr_eff_noPU_prompt_1sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_1sigma_EE.at(i));
	gr_eff_noPU_prompt_2sigma_EE->SetPoint(gr_eff_noPU_prompt_2sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_2sigma_EE.at(i));
	gr_eff_noPU_prompt_3sigma_EE->SetPoint(gr_eff_noPU_prompt_3sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_3sigma_EE.at(i));
	gr_eff_noPU_prompt_4sigma_EE->SetPoint(gr_eff_noPU_prompt_4sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_4sigma_EE.at(i));
	gr_eff_noPU_prompt_40_EE->SetPoint(gr_eff_noPU_prompt_40_EE->GetN(), 0.004*i, prompt_eff_noPU_40_EE.at(i));
	gr_eff_noPU_prompt_60_EE->SetPoint(gr_eff_noPU_prompt_60_EE->GetN(), 0.004*i, prompt_eff_noPU_60_EE.at(i));
	gr_eff_noPU_prompt_80_EE->SetPoint(gr_eff_noPU_prompt_80_EE->GetN(), 0.004*i, prompt_eff_noPU_80_EE.at(i));
	gr_eff_noPU_prompt_100_EE->SetPoint(gr_eff_noPU_prompt_100_EE->GetN(), 0.004*i, prompt_eff_noPU_100_EE.at(i));
//	gr_eff_noPU_prompt_gen_EE->SetPoint(gr_eff_noPU_prompt_gen_EE->GetN(), 0.004*i, prompt_eff_noPU_gen_EE.at(i));
      // nonprompt
	gr_eff_noPU_nonprompt_EB->SetPoint(gr_eff_noPU_nonprompt_EB->GetN(), 0.004*i, nonprompt_eff_noPU_EB.at(i));
	gr_eff_noPU_nonprompt_1sigma_EB->SetPoint(gr_eff_noPU_nonprompt_1sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_1sigma_EB.at(i));
	gr_eff_noPU_nonprompt_2sigma_EB->SetPoint(gr_eff_noPU_nonprompt_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_2sigma_EB.at(i));
	gr_eff_noPU_nonprompt_3sigma_EB->SetPoint(gr_eff_noPU_nonprompt_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_3sigma_EB.at(i));
	gr_eff_noPU_nonprompt_4sigma_EB->SetPoint(gr_eff_noPU_nonprompt_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_4sigma_EB.at(i));
	gr_eff_noPU_nonprompt_40_EB->SetPoint(gr_eff_noPU_nonprompt_40_EB->GetN(), 0.004*i, nonprompt_eff_noPU_40_EB.at(i));
	gr_eff_noPU_nonprompt_60_EB->SetPoint(gr_eff_noPU_nonprompt_60_EB->GetN(), 0.004*i, nonprompt_eff_noPU_60_EB.at(i));
	gr_eff_noPU_nonprompt_80_EB->SetPoint(gr_eff_noPU_nonprompt_80_EB->GetN(), 0.004*i, nonprompt_eff_noPU_80_EB.at(i));
	gr_eff_noPU_nonprompt_100_EB->SetPoint(gr_eff_noPU_nonprompt_100_EB->GetN(), 0.004*i, nonprompt_eff_noPU_100_EB.at(i));
//	gr_eff_noPU_nonprompt_gen_EB->SetPoint(gr_eff_noPU_nonprompt_gen_EB->GetN(), 0.004*i, nonprompt_eff_noPU_gen_EB.at(i));
	gr_eff_noPU_nonprompt_EE->SetPoint(gr_eff_noPU_nonprompt_EE->GetN(), 0.004*i, nonprompt_eff_noPU_EE.at(i));
	gr_eff_noPU_nonprompt_1sigma_EE->SetPoint(gr_eff_noPU_nonprompt_1sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_1sigma_EE.at(i));
	gr_eff_noPU_nonprompt_2sigma_EE->SetPoint(gr_eff_noPU_nonprompt_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_2sigma_EE.at(i));
	gr_eff_noPU_nonprompt_3sigma_EE->SetPoint(gr_eff_noPU_nonprompt_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_3sigma_EE.at(i));
	gr_eff_noPU_nonprompt_4sigma_EE->SetPoint(gr_eff_noPU_nonprompt_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_4sigma_EE.at(i));
	gr_eff_noPU_nonprompt_40_EE->SetPoint(gr_eff_noPU_nonprompt_40_EE->GetN(), 0.004*i, nonprompt_eff_noPU_40_EE.at(i));
	gr_eff_noPU_nonprompt_60_EE->SetPoint(gr_eff_noPU_nonprompt_60_EE->GetN(), 0.004*i, nonprompt_eff_noPU_60_EE.at(i));
	gr_eff_noPU_nonprompt_80_EE->SetPoint(gr_eff_noPU_nonprompt_80_EE->GetN(), 0.004*i, nonprompt_eff_noPU_80_EE.at(i));
	gr_eff_noPU_nonprompt_100_EE->SetPoint(gr_eff_noPU_nonprompt_100_EE->GetN(), 0.004*i, nonprompt_eff_noPU_100_EE.at(i));
//	gr_eff_noPU_nonprompt_gen_EE->SetPoint(gr_eff_noPU_nonprompt_gen_EE->GetN(), 0.004*i, nonprompt_eff_noPU_gen_EE.at(i));

  }

  ///////////////
  // Cosmetics //
  ///////////////
  // PU200
    // prompt
	  // EB
  gr_eff_PU200_prompt_EB->GetXaxis()->SetTitle("Relative isolation"); gr_eff_PU200_prompt_EB->GetYaxis()->SetTitle("Isolation efficiency");
  gr_eff_PU200_prompt_EB->SetLineColor(kBlack);
  gr_eff_PU200_prompt_40_EB->SetLineColor(kRed); gr_eff_PU200_prompt_60_EB->SetLineColor(kGreen); gr_eff_PU200_prompt_80_EB->SetLineColor(kBlue); gr_eff_PU200_prompt_100_EB->SetLineColor(kMagenta);
  gr_eff_PU200_prompt_1sigma_EB->SetLineColor(kRed); gr_eff_PU200_prompt_2sigma_EB->SetLineColor(kGreen); gr_eff_PU200_prompt_3sigma_EB->SetLineColor(kBlue); gr_eff_PU200_prompt_4sigma_EB->SetLineColor(kMagenta);
  gr_eff_noPU_prompt_EB->SetLineColor(kGray+1);
  gr_eff_PU200_prompt_EB->SetLineWidth(2); gr_eff_PU200_prompt_40_EB->SetLineWidth(2); gr_eff_PU200_prompt_60_EB->SetLineWidth(2); gr_eff_PU200_prompt_80_EB->SetLineWidth(2); gr_eff_PU200_prompt_100_EB->SetLineWidth(2); gr_eff_PU200_prompt_2sigma_EB->SetLineWidth(2); gr_eff_PU200_prompt_3sigma_EB->SetLineWidth(2); gr_eff_PU200_prompt_1sigma_EB->SetLineWidth(2); gr_eff_noPU_prompt_EB->SetLineWidth(2);
      // EE
  gr_eff_PU200_prompt_EE->GetXaxis()->SetTitle("Relative isolation"); gr_eff_PU200_prompt_EE->GetYaxis()->SetTitle("Isolation efficiency");
  gr_eff_PU200_prompt_EE->SetLineColor(kBlack);
  gr_eff_PU200_prompt_40_EE->SetLineColor(kRed); gr_eff_PU200_prompt_60_EE->SetLineColor(kGreen); gr_eff_PU200_prompt_80_EE->SetLineColor(kBlue); gr_eff_PU200_prompt_100_EE->SetLineColor(kMagenta);
  gr_eff_PU200_prompt_1sigma_EE->SetLineColor(kRed); gr_eff_PU200_prompt_2sigma_EE->SetLineColor(kGreen); gr_eff_PU200_prompt_3sigma_EE->SetLineColor(kBlue); gr_eff_PU200_prompt_4sigma_EE->SetLineColor(kMagenta);
  gr_eff_noPU_prompt_EE->SetLineColor(kGray+1);
  gr_eff_PU200_prompt_EE->SetLineWidth(2); gr_eff_PU200_prompt_40_EE->SetLineWidth(2); gr_eff_PU200_prompt_60_EE->SetLineWidth(2); gr_eff_PU200_prompt_80_EE->SetLineWidth(2); gr_eff_PU200_prompt_100_EE->SetLineWidth(2); gr_eff_PU200_prompt_2sigma_EE->SetLineWidth(2); gr_eff_PU200_prompt_3sigma_EE->SetLineWidth(2); gr_eff_PU200_prompt_1sigma_EE->SetLineWidth(2); gr_eff_noPU_prompt_EE->SetLineWidth(2);
    // nonprompt
	  // EB
  gr_eff_PU200_nonprompt_EB->GetXaxis()->SetTitle("Relative isolation"); gr_eff_PU200_nonprompt_EB->GetYaxis()->SetTitle("Isolation efficiency");
  gr_eff_PU200_nonprompt_EB->SetLineColor(kBlack);
  gr_eff_PU200_nonprompt_40_EB->SetLineColor(kRed); gr_eff_PU200_nonprompt_60_EB->SetLineColor(kGreen); gr_eff_PU200_nonprompt_80_EB->SetLineColor(kBlue); gr_eff_PU200_nonprompt_100_EB->SetLineColor(kMagenta);
  gr_eff_PU200_nonprompt_1sigma_EB->SetLineColor(kRed); gr_eff_PU200_nonprompt_2sigma_EB->SetLineColor(kGreen); gr_eff_PU200_nonprompt_3sigma_EB->SetLineColor(kBlue); gr_eff_PU200_nonprompt_4sigma_EB->SetLineColor(kMagenta);
  gr_eff_noPU_nonprompt_EB->SetLineColor(kGray+1);
  gr_eff_PU200_nonprompt_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_40_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_60_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_80_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_100_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_2sigma_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_3sigma_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_1sigma_EB->SetLineWidth(2); gr_eff_noPU_nonprompt_EB->SetLineWidth(2);
      // EE
  gr_eff_PU200_nonprompt_EE->GetXaxis()->SetTitle("Relative isolation"); gr_eff_PU200_nonprompt_EE->GetYaxis()->SetTitle("Isolation efficiency");
  gr_eff_PU200_nonprompt_EE->SetLineColor(kBlack);
  gr_eff_PU200_nonprompt_40_EE->SetLineColor(kRed); gr_eff_PU200_nonprompt_60_EE->SetLineColor(kGreen); gr_eff_PU200_nonprompt_80_EE->SetLineColor(kBlue); gr_eff_PU200_nonprompt_100_EE->SetLineColor(kMagenta);
  gr_eff_PU200_nonprompt_1sigma_EE->SetLineColor(kRed); gr_eff_PU200_nonprompt_2sigma_EE->SetLineColor(kGreen); gr_eff_PU200_nonprompt_3sigma_EE->SetLineColor(kBlue); gr_eff_PU200_nonprompt_4sigma_EE->SetLineColor(kMagenta);
  gr_eff_noPU_nonprompt_EE->SetLineColor(kGray+1);
  gr_eff_PU200_nonprompt_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_40_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_60_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_80_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_100_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_2sigma_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_3sigma_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_1sigma_EE->SetLineWidth(2); gr_eff_noPU_nonprompt_EE->SetLineWidth(2);



  /////////////
  // Legends //
  /////////////
  // prompt
    // dt
      // EB
  TLegend* leg_eff_prompt_EB_dt = new TLegend(0.40, 0.13, 0.88, 0.38);
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_40_EB, "40ps PU200 (muon, track)");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_60_EB, "60ps PU200 (muon, track)");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_80_EB, "80ps PU200 (muon, track)");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_100_EB, "100ps PU200 (muon, track)");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_noPU_prompt_EB, "no MTD noPU");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_EB, "no MTD PU200");
  leg_eff_prompt_EB_dt->SetTextSize(0.03);
      // EE
  TLegend* leg_eff_prompt_EE_dt = new TLegend(0.40, 0.13, 0.88, 0.38);
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_40_EE, "40ps PU200 (muon, track)");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_60_EE, "60ps PU200 (muon, track)");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_80_EE, "80ps PU200 (muon, track)");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_100_EE, "100ps PU200 (muon, track)");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_noPU_prompt_EE, "no MTD noPU");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_EE, "no MTD PU200");
  leg_eff_prompt_EE_dt->SetTextSize(0.03);
    // dtsig
      // EB
  TLegend* leg_eff_prompt_EB_dtsig = new TLegend(0.40, 0.13, 0.88, 0.38);
  leg_eff_prompt_EB_dtsig->AddEntry(gr_eff_PU200_prompt_1sigma_EB, "1sigma PU200 (muon, track)");
  leg_eff_prompt_EB_dtsig->AddEntry(gr_eff_PU200_prompt_2sigma_EB, "2sigma PU200 (muon, track)");
  leg_eff_prompt_EB_dtsig->AddEntry(gr_eff_PU200_prompt_3sigma_EB, "3sigma PU200 (muon, track)");
  leg_eff_prompt_EB_dtsig->AddEntry(gr_eff_PU200_prompt_4sigma_EB, "4sigma PU200 (muon, track)");
  leg_eff_prompt_EB_dtsig->AddEntry(gr_eff_noPU_prompt_EB, "no MTD noPU");
  leg_eff_prompt_EB_dtsig->AddEntry(gr_eff_PU200_prompt_EB, "no MTD PU200");
  leg_eff_prompt_EB_dtsig->SetTextSize(0.03);
      // EE
  TLegend* leg_eff_prompt_EE_dtsig = new TLegend(0.40, 0.13, 0.88, 0.38);
  leg_eff_prompt_EE_dtsig->AddEntry(gr_eff_PU200_prompt_1sigma_EE, "1sigma PU200 (muon, track)");
  leg_eff_prompt_EE_dtsig->AddEntry(gr_eff_PU200_prompt_2sigma_EE, "2sigma PU200 (muon, track)");
  leg_eff_prompt_EE_dtsig->AddEntry(gr_eff_PU200_prompt_3sigma_EE, "3sigma PU200 (muon, track)");
  leg_eff_prompt_EE_dtsig->AddEntry(gr_eff_PU200_prompt_4sigma_EE, "4sigma PU200 (muon, track)");
  leg_eff_prompt_EE_dtsig->AddEntry(gr_eff_noPU_prompt_EE, "no MTD noPU");
  leg_eff_prompt_EE_dtsig->AddEntry(gr_eff_PU200_prompt_EE, "no MTD PU200");
  leg_eff_prompt_EE_dtsig->SetTextSize(0.03);

  // nonprompt
    // dt
      // EB
  TLegend* leg_eff_nonprompt_EB_dt = new TLegend(0.40, 0.13, 0.88, 0.38);
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_40_EB, "40ps PU200 (muon, track)");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_60_EB, "60ps PU200 (muon, track)");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_80_EB, "80ps PU200 (muon, track)");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_100_EB, "100ps PU200 (muon, track)");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_noPU_nonprompt_EB, "no MTD noPU");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_EB, "no MTD PU200");
  leg_eff_nonprompt_EB_dt->SetTextSize(0.03);
      // EE
  TLegend* leg_eff_nonprompt_EE_dt = new TLegend(0.40, 0.13, 0.88, 0.38);
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_PU200_nonprompt_40_EE, "40ps PU200 (muon, track)");
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_PU200_nonprompt_60_EE, "60ps PU200 (muon, track)");
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_PU200_nonprompt_80_EE, "80ps PU200 (muon, track)");
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_PU200_nonprompt_100_EE, "100ps PU200 (muon, track)");
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_noPU_nonprompt_EE, "no MTD noPU");
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_PU200_nonprompt_EE, "no MTD PU200");
  leg_eff_nonprompt_EE_dt->SetTextSize(0.03);
    // dtsig
      // EB
  TLegend* leg_eff_nonprompt_EB_dtsig = new TLegend(0.40, 0.13, 0.88, 0.38);
  leg_eff_nonprompt_EB_dtsig->AddEntry(gr_eff_PU200_nonprompt_1sigma_EB, "1sigma PU200 (muon, track)");
  leg_eff_nonprompt_EB_dtsig->AddEntry(gr_eff_PU200_nonprompt_2sigma_EB, "2sigma PU200 (muon, track)");
  leg_eff_nonprompt_EB_dtsig->AddEntry(gr_eff_PU200_nonprompt_3sigma_EB, "3sigma PU200 (muon, track)");
  leg_eff_nonprompt_EB_dtsig->AddEntry(gr_eff_PU200_nonprompt_4sigma_EB, "4sigma PU200 (muon, track)");
  leg_eff_nonprompt_EB_dtsig->AddEntry(gr_eff_noPU_nonprompt_EB, "no MTD noPU");
  leg_eff_nonprompt_EB_dtsig->AddEntry(gr_eff_PU200_nonprompt_EB, "no MTD PU200");
  leg_eff_nonprompt_EB_dtsig->SetTextSize(0.03);
      // EE
  TLegend* leg_eff_nonprompt_EE_dtsig = new TLegend(0.40, 0.13, 0.88, 0.38);
  leg_eff_nonprompt_EE_dtsig->AddEntry(gr_eff_PU200_nonprompt_1sigma_EE, "1sigma PU200 (muon, track)");
  leg_eff_nonprompt_EE_dtsig->AddEntry(gr_eff_PU200_nonprompt_2sigma_EE, "2sigma PU200 (muon, track)");
  leg_eff_nonprompt_EE_dtsig->AddEntry(gr_eff_PU200_nonprompt_3sigma_EE, "3sigma PU200 (muon, track)");
  leg_eff_nonprompt_EE_dtsig->AddEntry(gr_eff_PU200_nonprompt_4sigma_EE, "4sigma PU200 (muon, track)");
  leg_eff_nonprompt_EE_dtsig->AddEntry(gr_eff_noPU_nonprompt_EE, "no MTD noPU");
  leg_eff_nonprompt_EE_dtsig->AddEntry(gr_eff_PU200_nonprompt_EE, "no MTD PU200");
  leg_eff_nonprompt_EE_dtsig->SetTextSize(0.03);



  ////////////////
  // Draw plots //
  ////////////////
  
  // prompt
    // Barrel
  TCanvas* c_PU200_prompt_EB_reliso = new TCanvas("c_PU200_prompt_EB_reliso", "c_PU200_prompt_EB_reliso", 1200, 1200);
  c_PU200_prompt_EB_reliso->cd();
  c_PU200_prompt_EB_reliso->SetGrid();
  c_PU200_prompt_EB_reliso->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EB->GetXaxis()->SetRangeUser(0.,0.28);
  gr_eff_PU200_prompt_EB->SetTitle("");
  gr_eff_PU200_prompt_EB->Draw("AL");
  gr_eff_PU200_prompt_1sigma_EB->Draw("same");
  gr_eff_PU200_prompt_2sigma_EB->Draw("same");
  gr_eff_PU200_prompt_3sigma_EB->Draw("same");
  gr_eff_PU200_prompt_4sigma_EB->Draw("same");
  gr_eff_noPU_prompt_EB->Draw("same");
  leg_eff_prompt_EB_dtsig->Draw();
  c_PU200_prompt_EB_reliso->Print("plots/ntuple/isoeff_PU200_prompt_EB.pdf");
	// Endcap
  TCanvas* c_PU200_prompt_EE_reliso = new TCanvas("c_PU200_prompt_EE_reliso", "c_PU200_prompt_EE_reliso", 1200, 1200);
  c_PU200_prompt_EE_reliso->cd();
  c_PU200_prompt_EE_reliso->SetGrid();
  c_PU200_prompt_EE_reliso->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EE->GetXaxis()->SetRangeUser(0.,0.28);
  gr_eff_PU200_prompt_EE->SetTitle("");
  gr_eff_PU200_prompt_EE->Draw("AL");
  gr_eff_PU200_prompt_1sigma_EE->Draw("same");
  gr_eff_PU200_prompt_2sigma_EE->Draw("same");
  gr_eff_PU200_prompt_3sigma_EE->Draw("same");
  gr_eff_PU200_prompt_4sigma_EE->Draw("same");
  gr_eff_noPU_prompt_EE->Draw("same");
  leg_eff_prompt_EE_dtsig->Draw();
  c_PU200_prompt_EE_reliso->Print("plots/ntuple/isoeff_PU200_prompt_EE.pdf");
  // nonprompt
	// Barrel
  TCanvas* c_PU200_nonprompt_EB_reliso = new TCanvas("c_PU200_nonprompt_EB_reliso", "c_PU200_nonprompt_EB_reliso", 1200, 1200);
  c_PU200_nonprompt_EB_reliso->cd();
  c_PU200_nonprompt_EB_reliso->SetGrid();
  c_PU200_nonprompt_EB_reliso->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EB->GetXaxis()->SetRangeUser(0.,4.);
  gr_eff_PU200_nonprompt_EB->SetTitle("");
  gr_eff_PU200_nonprompt_EB->Draw("AL");
  gr_eff_PU200_nonprompt_1sigma_EB->Draw("same");
  gr_eff_PU200_nonprompt_2sigma_EB->Draw("same");
  gr_eff_PU200_nonprompt_3sigma_EB->Draw("same");
  gr_eff_PU200_nonprompt_4sigma_EB->Draw("same");
  gr_eff_noPU_nonprompt_EB->Draw("same");
  leg_eff_nonprompt_EB_dtsig->Draw();
  c_PU200_nonprompt_EB_reliso->Print("plots/ntuple/isoeff_PU200_nonprompt_EB.pdf");
	// Endcap
  TCanvas* c_PU200_nonprompt_EE_reliso = new TCanvas("c_PU200_nonprompt_EE_reliso", "c_PU200_nonprompt_EE_reliso", 1200, 1200);
  c_PU200_nonprompt_EE_reliso->cd();
  c_PU200_nonprompt_EE_reliso->SetGrid();
  c_PU200_nonprompt_EE_reliso->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EE->GetXaxis()->SetRangeUser(0.,4.);
  gr_eff_PU200_nonprompt_EE->SetTitle("");
  gr_eff_PU200_nonprompt_EE->Draw("AL");
  gr_eff_PU200_nonprompt_1sigma_EE->Draw("same");
  gr_eff_PU200_nonprompt_2sigma_EE->Draw("same");
  gr_eff_PU200_nonprompt_3sigma_EE->Draw("same");
  gr_eff_PU200_nonprompt_4sigma_EE->Draw("same");
  gr_eff_noPU_nonprompt_EE->Draw("same");
  leg_eff_nonprompt_EE_dtsig->Draw();
  c_PU200_nonprompt_EE_reliso->Print("plots/ntuple/isoeff_PU200_nonprompt_EE.pdf");





  /////////////////////
  ///// ROC Curve /////
  /////////////////////

  // Define TGraph
//  TGraph* gr_roc_PU200_EB = new TGraph(gr_roc_PU200_EB->GetN(), &prompt_eff_PU200_EB[0], &nonprompt_eff_PU200_EB[0]);
  // PU200
  TGraph* gr_roc_PU200_EB = new TGraph();				TGraph* gr_roc_PU200_EE = new TGraph();
  TGraph* gr_roc_PU200_1sigma_EB = new TGraph();		TGraph* gr_roc_PU200_1sigma_EE = new TGraph();
  TGraph* gr_roc_PU200_2sigma_EB = new TGraph();		TGraph* gr_roc_PU200_2sigma_EE = new TGraph();
  TGraph* gr_roc_PU200_3sigma_EB = new TGraph();		TGraph* gr_roc_PU200_3sigma_EE = new TGraph();
  TGraph* gr_roc_PU200_4sigma_EB = new TGraph();		TGraph* gr_roc_PU200_4sigma_EE = new TGraph();
  TGraph* gr_roc_PU200_40_EB = new TGraph();			TGraph* gr_roc_PU200_40_EE = new TGraph();
  TGraph* gr_roc_PU200_60_EB = new TGraph();			TGraph* gr_roc_PU200_60_EE = new TGraph();
  TGraph* gr_roc_PU200_80_EB = new TGraph();			TGraph* gr_roc_PU200_80_EE = new TGraph();
  TGraph* gr_roc_PU200_100_EB = new TGraph();			TGraph* gr_roc_PU200_100_EE = new TGraph();
  // noPU
  TGraph* gr_roc_noPU_EB = new TGraph();				TGraph* gr_roc_noPU_EE = new TGraph();
  TGraph* gr_roc_noPU_1sigma_EB = new TGraph();			TGraph* gr_roc_noPU_1sigma_EE = new TGraph();
  TGraph* gr_roc_noPU_2sigma_EB = new TGraph();			TGraph* gr_roc_noPU_2sigma_EE = new TGraph();
  TGraph* gr_roc_noPU_3sigma_EB = new TGraph();			TGraph* gr_roc_noPU_3sigma_EE = new TGraph();
  TGraph* gr_roc_noPU_4sigma_EB = new TGraph();			TGraph* gr_roc_noPU_4sigma_EE = new TGraph();
  TGraph* gr_roc_noPU_40_EB = new TGraph();				TGraph* gr_roc_noPU_40_EE = new TGraph();
  TGraph* gr_roc_noPU_60_EB = new TGraph();				TGraph* gr_roc_noPU_60_EE = new TGraph();
  TGraph* gr_roc_noPU_80_EB = new TGraph();				TGraph* gr_roc_noPU_80_EE = new TGraph();
  TGraph* gr_roc_noPU_100_EB = new TGraph();			TGraph* gr_roc_noPU_100_EE = new TGraph();


  // Draw ROC Curve
  for(int i=0; i<nbin+1; i++) {
	// PU200
	  // EB
	gr_roc_PU200_EB->SetPoint(i, 		prompt_eff_PU200_EB.at(i), 			nonprompt_eff_PU200_EB.at(i));
	gr_roc_PU200_1sigma_EB->SetPoint(i, prompt_eff_PU200_1sigma_EB.at(i), 	nonprompt_eff_PU200_1sigma_EB.at(i));
	gr_roc_PU200_2sigma_EB->SetPoint(i, prompt_eff_PU200_2sigma_EB.at(i), 	nonprompt_eff_PU200_2sigma_EB.at(i));
	gr_roc_PU200_3sigma_EB->SetPoint(i, prompt_eff_PU200_3sigma_EB.at(i), 	nonprompt_eff_PU200_3sigma_EB.at(i));
	gr_roc_PU200_4sigma_EB->SetPoint(i, prompt_eff_PU200_4sigma_EB.at(i), 	nonprompt_eff_PU200_4sigma_EB.at(i));
	gr_roc_PU200_40_EB->SetPoint(i, 	prompt_eff_PU200_40_EB.at(i), 		nonprompt_eff_PU200_40_EB.at(i));
	gr_roc_PU200_60_EB->SetPoint(i, 	prompt_eff_PU200_60_EB.at(i), 		nonprompt_eff_PU200_60_EB.at(i));
	gr_roc_PU200_80_EB->SetPoint(i, 	prompt_eff_PU200_80_EB.at(i), 		nonprompt_eff_PU200_80_EB.at(i));
	gr_roc_PU200_100_EB->SetPoint(i, 	prompt_eff_PU200_100_EB.at(i), 		nonprompt_eff_PU200_100_EB.at(i));
	  // EE
	gr_roc_PU200_EE->SetPoint(i, 		prompt_eff_PU200_EE.at(i), 			nonprompt_eff_PU200_EE.at(i));
	gr_roc_PU200_1sigma_EE->SetPoint(i, prompt_eff_PU200_1sigma_EE.at(i), 	nonprompt_eff_PU200_1sigma_EE.at(i));
	gr_roc_PU200_2sigma_EE->SetPoint(i, prompt_eff_PU200_2sigma_EE.at(i), 	nonprompt_eff_PU200_2sigma_EE.at(i));
	gr_roc_PU200_3sigma_EE->SetPoint(i, prompt_eff_PU200_3sigma_EE.at(i), 	nonprompt_eff_PU200_3sigma_EE.at(i));
	gr_roc_PU200_4sigma_EE->SetPoint(i, prompt_eff_PU200_4sigma_EE.at(i), 	nonprompt_eff_PU200_4sigma_EE.at(i));
	gr_roc_PU200_40_EE->SetPoint(i, 	prompt_eff_PU200_40_EE.at(i), 		nonprompt_eff_PU200_40_EE.at(i));
	gr_roc_PU200_60_EE->SetPoint(i, 	prompt_eff_PU200_60_EE.at(i), 		nonprompt_eff_PU200_60_EE.at(i));
	gr_roc_PU200_80_EE->SetPoint(i, 	prompt_eff_PU200_80_EE.at(i), 		nonprompt_eff_PU200_80_EE.at(i));
	gr_roc_PU200_100_EE->SetPoint(i, 	prompt_eff_PU200_100_EE.at(i), 		nonprompt_eff_PU200_100_EE.at(i));
	// noPU
	  // EB
	gr_roc_noPU_EB->SetPoint(i, 		prompt_eff_noPU_EB.at(i), 			nonprompt_eff_noPU_EB.at(i));
	gr_roc_noPU_1sigma_EB->SetPoint(i, 	prompt_eff_noPU_1sigma_EB.at(i), 	nonprompt_eff_noPU_1sigma_EB.at(i));
	gr_roc_noPU_2sigma_EB->SetPoint(i, 	prompt_eff_noPU_2sigma_EB.at(i), 	nonprompt_eff_noPU_2sigma_EB.at(i));
	gr_roc_noPU_3sigma_EB->SetPoint(i, 	prompt_eff_noPU_3sigma_EB.at(i), 	nonprompt_eff_noPU_3sigma_EB.at(i));
	gr_roc_noPU_4sigma_EB->SetPoint(i, 	prompt_eff_noPU_4sigma_EB.at(i), 	nonprompt_eff_noPU_4sigma_EB.at(i));
	gr_roc_noPU_40_EB->SetPoint(i, 		prompt_eff_noPU_40_EB.at(i), 		nonprompt_eff_noPU_40_EB.at(i));
	gr_roc_noPU_60_EB->SetPoint(i, 		prompt_eff_noPU_60_EB.at(i), 		nonprompt_eff_noPU_60_EB.at(i));
	gr_roc_noPU_80_EB->SetPoint(i, 		prompt_eff_noPU_80_EB.at(i), 		nonprompt_eff_noPU_80_EB.at(i));
	gr_roc_noPU_100_EB->SetPoint(i, 	prompt_eff_noPU_100_EB.at(i), 		nonprompt_eff_noPU_100_EB.at(i));
	  // EE
	gr_roc_noPU_EE->SetPoint(i, 		prompt_eff_noPU_EE.at(i), 			nonprompt_eff_noPU_EE.at(i));
	gr_roc_noPU_1sigma_EE->SetPoint(i, 	prompt_eff_noPU_1sigma_EE.at(i), 	nonprompt_eff_noPU_1sigma_EE.at(i));
	gr_roc_noPU_2sigma_EE->SetPoint(i, 	prompt_eff_noPU_2sigma_EE.at(i), 	nonprompt_eff_noPU_2sigma_EE.at(i));
	gr_roc_noPU_3sigma_EE->SetPoint(i, 	prompt_eff_noPU_3sigma_EE.at(i), 	nonprompt_eff_noPU_3sigma_EE.at(i));
	gr_roc_noPU_4sigma_EE->SetPoint(i, 	prompt_eff_noPU_4sigma_EE.at(i), 	nonprompt_eff_noPU_4sigma_EE.at(i));
	gr_roc_noPU_40_EE->SetPoint(i, 		prompt_eff_noPU_40_EE.at(i), 		nonprompt_eff_noPU_40_EE.at(i));
	gr_roc_noPU_60_EE->SetPoint(i, 		prompt_eff_noPU_60_EE.at(i), 		nonprompt_eff_noPU_60_EE.at(i));
	gr_roc_noPU_80_EE->SetPoint(i, 		prompt_eff_noPU_80_EE.at(i), 		nonprompt_eff_noPU_80_EE.at(i));
	gr_roc_noPU_100_EE->SetPoint(i, 	prompt_eff_noPU_100_EE.at(i), 		nonprompt_eff_noPU_100_EE.at(i));
  }

  ///////////////
  // Cosmetics //
  ///////////////
  // PU200
  gr_roc_PU200_EB->SetTitle(""); gr_roc_PU200_EB->GetXaxis()->SetTitle("Prompt efficiency"); gr_roc_PU200_EB->GetYaxis()->SetTitle("Non-prompt efficiency");
  gr_roc_PU200_EB->SetLineWidth(2); gr_roc_PU200_1sigma_EB->SetLineWidth(2); gr_roc_PU200_2sigma_EB->SetLineWidth(2); gr_roc_PU200_3sigma_EB->SetLineWidth(2); gr_roc_PU200_4sigma_EB->SetLineWidth(2); gr_roc_PU200_40_EB->SetLineWidth(2); gr_roc_PU200_60_EB->SetLineWidth(2); gr_roc_PU200_80_EB->SetLineWidth(2); gr_roc_PU200_100_EB->SetLineWidth(2);
  gr_roc_PU200_EB->SetLineColor(kBlack); gr_roc_PU200_40_EB->SetLineColor(kRed); gr_roc_PU200_60_EB->SetLineColor(kGreen); gr_roc_PU200_80_EB->SetLineColor(kBlue); gr_roc_PU200_100_EB->SetLineColor(kMagenta); gr_roc_PU200_1sigma_EB->SetLineColor(kRed); gr_roc_PU200_2sigma_EB->SetLineColor(kGreen); gr_roc_PU200_3sigma_EB->SetLineColor(kBlue); gr_roc_PU200_4sigma_EB->SetLineColor(kMagenta);
  gr_roc_PU200_EE->SetTitle(""); gr_roc_PU200_EE->GetXaxis()->SetTitle("Prompt efficiency"); gr_roc_PU200_EE->GetYaxis()->SetTitle("Non-prompt efficiency");
  gr_roc_PU200_EE->SetLineWidth(2); gr_roc_PU200_1sigma_EE->SetLineWidth(2); gr_roc_PU200_2sigma_EE->SetLineWidth(2); gr_roc_PU200_3sigma_EE->SetLineWidth(2); gr_roc_PU200_4sigma_EE->SetLineWidth(2); gr_roc_PU200_40_EE->SetLineWidth(2); gr_roc_PU200_60_EE->SetLineWidth(2); gr_roc_PU200_80_EE->SetLineWidth(2); gr_roc_PU200_100_EE->SetLineWidth(2);
  gr_roc_PU200_EE->SetLineColor(kBlack); gr_roc_PU200_40_EE->SetLineColor(kRed); gr_roc_PU200_60_EE->SetLineColor(kGreen); gr_roc_PU200_80_EE->SetLineColor(kBlue); gr_roc_PU200_100_EE->SetLineColor(kMagenta); gr_roc_PU200_1sigma_EE->SetLineColor(kRed); gr_roc_PU200_2sigma_EE->SetLineColor(kGreen); gr_roc_PU200_3sigma_EE->SetLineColor(kBlue); gr_roc_PU200_4sigma_EE->SetLineColor(kMagenta);
  // noPU
  gr_roc_noPU_EB->SetTitle(""); gr_roc_noPU_EB->GetXaxis()->SetTitle("Prompt efficiency"); gr_roc_noPU_EB->GetYaxis()->SetTitle("Non-prompt efficiency");
  gr_roc_noPU_EB->SetLineWidth(2); gr_roc_noPU_1sigma_EB->SetLineWidth(2); gr_roc_noPU_2sigma_EB->SetLineWidth(2); gr_roc_noPU_3sigma_EB->SetLineWidth(2); gr_roc_noPU_4sigma_EB->SetLineWidth(2); gr_roc_noPU_40_EB->SetLineWidth(2); gr_roc_noPU_60_EB->SetLineWidth(2); gr_roc_noPU_80_EB->SetLineWidth(2); gr_roc_noPU_100_EB->SetLineWidth(2);
  gr_roc_noPU_EB->SetLineColor(kGray); gr_roc_noPU_40_EB->SetLineColor(kRed); gr_roc_noPU_60_EB->SetLineColor(kGreen); gr_roc_noPU_80_EB->SetLineColor(kBlue); gr_roc_noPU_100_EB->SetLineColor(kMagenta); gr_roc_noPU_1sigma_EB->SetLineColor(kRed); gr_roc_noPU_2sigma_EB->SetLineColor(kGreen); gr_roc_noPU_3sigma_EB->SetLineColor(kBlue); gr_roc_noPU_4sigma_EB->SetLineColor(kMagenta);
  gr_roc_noPU_EE->SetTitle(""); gr_roc_noPU_EE->GetXaxis()->SetTitle("Prompt efficiency"); gr_roc_noPU_EE->GetYaxis()->SetTitle("Non-prompt efficiency");
  gr_roc_noPU_EE->SetLineWidth(2); gr_roc_noPU_1sigma_EE->SetLineWidth(2); gr_roc_noPU_2sigma_EE->SetLineWidth(2); gr_roc_noPU_3sigma_EE->SetLineWidth(2); gr_roc_noPU_4sigma_EE->SetLineWidth(2); gr_roc_noPU_40_EE->SetLineWidth(2); gr_roc_noPU_60_EE->SetLineWidth(2); gr_roc_noPU_80_EE->SetLineWidth(2); gr_roc_noPU_100_EE->SetLineWidth(2);
  gr_roc_noPU_EE->SetLineColor(kGray); gr_roc_noPU_40_EE->SetLineColor(kRed); gr_roc_noPU_60_EE->SetLineColor(kGreen); gr_roc_noPU_80_EE->SetLineColor(kBlue); gr_roc_noPU_100_EE->SetLineColor(kMagenta); gr_roc_noPU_1sigma_EE->SetLineColor(kRed); gr_roc_noPU_2sigma_EE->SetLineColor(kGreen); gr_roc_noPU_3sigma_EE->SetLineColor(kBlue); gr_roc_noPU_1sigma_EE->SetLineColor(kMagenta);


  ///////////////
  //  Legends  //
  ///////////////
  TLegend* leg_roc_EB_dt = new TLegend(0.15, 0.68, 0.41, 0.88);
  leg_roc_EB_dt->AddEntry(gr_roc_noPU_EB, "no MTD noPU");
  leg_roc_EB_dt->AddEntry(gr_roc_PU200_40_EB, "40 ps PU200");
  leg_roc_EB_dt->AddEntry(gr_roc_PU200_60_EB, "60 ps PU200");
  leg_roc_EB_dt->AddEntry(gr_roc_PU200_80_EB, "80 ps PU200");
  leg_roc_EB_dt->AddEntry(gr_roc_PU200_100_EB, "100 ps PU200");
  leg_roc_EB_dt->AddEntry(gr_roc_PU200_EB, "no MTD PU200");
  leg_roc_EB_dt->SetTextSize(0.03);

  TLegend* leg_roc_EE_dt = new TLegend(0.15, 0.68, 0.41, 0.88);
  leg_roc_EE_dt->AddEntry(gr_roc_noPU_EE, "no MTD noPU");
  leg_roc_EE_dt->AddEntry(gr_roc_PU200_40_EE, "40 ps PU200");
  leg_roc_EE_dt->AddEntry(gr_roc_PU200_60_EE, "60 ps PU200");
  leg_roc_EE_dt->AddEntry(gr_roc_PU200_80_EE, "80 ps PU200");
  leg_roc_EE_dt->AddEntry(gr_roc_PU200_100_EE, "100 ps PU200");
  leg_roc_EE_dt->AddEntry(gr_roc_PU200_EE, "no MTD PU200");
  leg_roc_EE_dt->SetTextSize(0.03);

  TLegend* leg_roc_EB_dtsig = new TLegend(0.15, 0.68, 0.41, 0.88);
  leg_roc_EB_dtsig->AddEntry(gr_roc_noPU_EB, "no MTD noPU");
  leg_roc_EB_dtsig->AddEntry(gr_roc_PU200_1sigma_EB, "1sigma PU200");
  leg_roc_EB_dtsig->AddEntry(gr_roc_PU200_2sigma_EB, "2sigma PU200");
  leg_roc_EB_dtsig->AddEntry(gr_roc_PU200_3sigma_EB, "3sigma PU200");
  leg_roc_EB_dtsig->AddEntry(gr_roc_PU200_4sigma_EB, "4sigma PU200");
  leg_roc_EB_dtsig->AddEntry(gr_roc_PU200_EB, "no MTD PU200");
  leg_roc_EB_dtsig->SetTextSize(0.03);

  TLegend* leg_roc_EE_dtsig = new TLegend(0.15, 0.68, 0.41, 0.88);
  leg_roc_EE_dtsig->AddEntry(gr_roc_noPU_EE, "no MTD noPU");
  leg_roc_EE_dtsig->AddEntry(gr_roc_PU200_1sigma_EE, "1sigma PU200");
  leg_roc_EE_dtsig->AddEntry(gr_roc_PU200_2sigma_EE, "2sigma PU200");
  leg_roc_EE_dtsig->AddEntry(gr_roc_PU200_3sigma_EE, "3sigma PU200");
  leg_roc_EE_dtsig->AddEntry(gr_roc_PU200_4sigma_EE, "4sigma PU200");
  leg_roc_EE_dtsig->AddEntry(gr_roc_PU200_EE, "no MTD PU200");
  leg_roc_EE_dtsig->SetTextSize(0.03);




  ////////////////
  //    Draw    //
  ////////////////
  // dt
  TCanvas* c_PU200_EB_reliso_roc_dt = new TCanvas("c_PU200_EB_reliso_roc_dt", "c_PU200_EB_reliso_roc_dt", 1200, 1200);
  c_PU200_EB_reliso_roc_dt->cd();
  c_PU200_EB_reliso_roc_dt->SetGrid();
  c_PU200_EB_reliso_roc_dt->SetLeftMargin(0.12);
  gr_roc_PU200_EB->GetXaxis()->SetRangeUser(0.85,1.);
  gr_roc_PU200_EB->GetYaxis()->SetRangeUser(0.,0.4);
  gr_roc_PU200_EB->Draw("AL");
  gr_roc_PU200_40_EB->Draw("same");
  gr_roc_PU200_60_EB->Draw("same");
  gr_roc_PU200_80_EB->Draw("same");
  gr_roc_PU200_100_EB->Draw("same");
  gr_roc_noPU_EB->Draw("same");
  leg_roc_EB_dt->Draw();
  c_PU200_EB_reliso_roc_dt->Print("plots/ntuple/isoroc_dt_PU200_EB.pdf");

  TCanvas* c_PU200_EE_reliso_roc_dt = new TCanvas("c_PU200_EE_reliso_roc_dt", "c_PU200_EE_reliso_roc_dt", 1200, 1200);
  c_PU200_EE_reliso_roc_dt->cd();
  c_PU200_EE_reliso_roc_dt->SetGrid();
  c_PU200_EE_reliso_roc_dt->SetLeftMargin(0.12);
  gr_roc_PU200_EE->GetXaxis()->SetRangeUser(0.85,1.);
  gr_roc_PU200_EE->GetYaxis()->SetRangeUser(0.,0.4);
  gr_roc_PU200_EE->Draw("AL");
  gr_roc_PU200_40_EE->Draw("same");
  gr_roc_PU200_60_EE->Draw("same");
  gr_roc_PU200_80_EE->Draw("same");
  gr_roc_PU200_100_EE->Draw("same");
  gr_roc_noPU_EE->Draw("same");
  leg_roc_EE_dt->Draw();
  c_PU200_EE_reliso_roc_dt->Print("plots/ntuple/isoroc_dt_PU200_EE.pdf");

  // dtsig
  TCanvas* c_PU200_EB_reliso_roc_dtsig = new TCanvas("c_PU200_EB_reliso_roc_dtsig", "c_PU200_EB_reliso_roc_dtsig", 1200, 1200);
  c_PU200_EB_reliso_roc_dtsig->cd();
  c_PU200_EB_reliso_roc_dtsig->SetGrid();
  c_PU200_EB_reliso_roc_dtsig->SetLeftMargin(0.12);
  gr_roc_PU200_EB->GetXaxis()->SetRangeUser(0.85,1.);
  gr_roc_PU200_EB->GetYaxis()->SetRangeUser(0.,0.4);
  gr_roc_PU200_EB->Draw("AL");
  gr_roc_PU200_1sigma_EB->Draw("same");
  gr_roc_PU200_2sigma_EB->Draw("same");
  gr_roc_PU200_3sigma_EB->Draw("same");
  gr_roc_PU200_4sigma_EB->Draw("same");
  gr_roc_noPU_EB->Draw("same");
  leg_roc_EB_dtsig->Draw();
  c_PU200_EB_reliso_roc_dtsig->Print("plots/ntuple/isoroc_dtsig_PU200_EB.pdf");

  TCanvas* c_PU200_EE_reliso_roc_dtsig = new TCanvas("c_PU200_EE_reliso_roc_dtsig", "c_PU200_EE_reliso_roc_dtsig", 1200, 1200);
  c_PU200_EE_reliso_roc_dtsig->cd();
  c_PU200_EE_reliso_roc_dtsig->SetGrid();
  c_PU200_EE_reliso_roc_dtsig->SetLeftMargin(0.12);
  gr_roc_PU200_EE->GetXaxis()->SetRangeUser(0.85,1.);
  gr_roc_PU200_EE->GetYaxis()->SetRangeUser(0.,0.4);
  gr_roc_PU200_EE->Draw("AL");
  gr_roc_PU200_1sigma_EE->Draw("same");
  gr_roc_PU200_2sigma_EE->Draw("same");
  gr_roc_PU200_3sigma_EE->Draw("same");
  gr_roc_PU200_4sigma_EE->Draw("same");
  gr_roc_noPU_EE->Draw("same");
  leg_roc_EE_dtsig->Draw();
  c_PU200_EE_reliso_roc_dtsig->Print("plots/ntuple/isoroc_dtsig_PU200_EE.pdf");



}
/*
void draw_track_type_sigma_ntuple() {

  // generate the dictionary for vector<vector<sth>> collection
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<bool> >", "vector");

  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  TFile *f_PU200_prompt_vtx, *f_PU200_nonprompt_vtx, *f_noPU_prompt_vtx, *f_noPU_nonprompt_vtx;


  // muon track
  // PU200
  TChain* ch_PU200_prompt = new TChain("mtdMuonIsoValid/muonIso");
  TChain* ch_PU200_nonprompt = new TChain("mtdMuonIsoValid/muonIso");
  ch_PU200_prompt->Add(Form("data/%s", ntuple_PU200_prompt.Data()));
  ch_PU200_nonprompt->Add(Form("data/%s", ntuple_PU200_nonprompt.Data()));
  // noPU
  TChain* ch_noPU_prompt = new TChain("mtdMuonIsoValid/muonIso");
  TChain* ch_noPU_nonprompt = new TChain("mtdMuonIsoValid/muonIso");
  ch_noPU_prompt->Add(Form("data/%s", ntuple_noPU_prompt.Data()));
  ch_noPU_nonprompt->Add(Form("data/%s", ntuple_noPU_nonprompt.Data()));

  // vertex
  // PU200
  TChain* ch_PU200_prompt_vtx = new TChain("mtdMuonIsoValid/muonIso");
  TChain* ch_PU200_nonprompt_vtx = new TChain("mtdMuonIsoValid/muonIso");
  ch_PU200_prompt_vtx->Add(Form("data/%s", ntuple_PU200_prompt_vtx.Data()));
  ch_PU200_nonprompt_vtx->Add(Form("data/%s", ntuple_PU200_nonprompt_vtx.Data()));
  // noPU
  TChain* ch_noPU_prompt_vtx = new TChain("mtdMuonIsoValid/muonIso");
  TChain* ch_noPU_nonprompt_vtx = new TChain("mtdMuonIsoValid/muonIso");
  ch_noPU_prompt_vtx->Add(Form("data/%s", ntuple_noPU_prompt_vtx.Data()));
  ch_noPU_nonprompt_vtx->Add(Form("data/%s", ntuple_noPU_nonprompt_vtx.Data()));

  // muon
  vector<float> *muon_pt_(0), *muon_time_(0), *muon_time_err_(0), *muon_prompt_(0), *muon_isBarrel_(0);
  vector<float> *muon_status_(0), *muon_pv_dz_(0), *muon_pv_dxy_(0), *muon_vz_(0);
  vector<float> *muon_PVweight_(0);
  // track
  std::vector<std::vector<float>> *track_pt_(0), *track_time_(0), *track_time_err_(0);
  std::vector<std::vector<float>> *track_pv_dz_(0), *track_vz_(0);
  std::vector<std::vector<float>> *track_PVweight_(0);
  vector<int>   *muon_index_(0);
  std::vector<std::vector<int>> *track_index_(0);
  std::vector<std::vector<int>> *track_type_(0);
  std::vector<std::vector<bool>> *selectedLV_(0), *match_vtx_sim2reco_(0), *match_vtx_reco2sim_(0);

  // PU200
  ch_PU200_prompt->SetBranchAddress("muon_pt_",        &muon_pt_);
  ch_PU200_prompt->SetBranchAddress("muon_time_",      &muon_time_);
  ch_PU200_prompt->SetBranchAddress("muon_time_err_",  &muon_time_err_);
  ch_PU200_prompt->SetBranchAddress("muon_prompt_",    &muon_prompt_);
  ch_PU200_prompt->SetBranchAddress("muon_isBarrel_",  &muon_isBarrel_);
  ch_PU200_prompt->SetBranchAddress("muon_status_",    &muon_status_);
  ch_PU200_prompt->SetBranchAddress("muon_pv_dz_",     &muon_pv_dz_);
  ch_PU200_prompt->SetBranchAddress("muon_pv_dxy_",    &muon_pv_dxy_);
  ch_PU200_prompt->SetBranchAddress("muon_vz_",        &muon_vz_);
  ch_PU200_prompt->SetBranchAddress("muon_PVweight_",  &muon_PVweight_);
  ch_PU200_prompt->SetBranchAddress("track_pt_",       &track_pt_);
  ch_PU200_prompt->SetBranchAddress("track_time_",     &track_time_);
  ch_PU200_prompt->SetBranchAddress("track_time_err_", &track_time_err_);
  ch_PU200_prompt->SetBranchAddress("track_vz_",       &track_vz_);
  ch_PU200_prompt->SetBranchAddress("track_pv_dz_",    &track_pv_dz_);
  ch_PU200_prompt->SetBranchAddress("track_PVweight_", &track_PVweight_);
  ch_PU200_prompt->SetBranchAddress("track_type_",     &track_type_);
  ch_PU200_prompt->SetBranchAddress("selectedLV_",     &selectedLV_);
  ch_PU200_prompt->SetBranchAddress("match_vtx_sim2reco_",     &match_vtx_sim2reco_);
  ch_PU200_prompt->SetBranchAddress("match_vtx_reco2sim_",     &match_vtx_reco2sim_);
  ch_PU200_nonprompt->SetBranchAddress("muon_pt_",        &muon_pt_);
  ch_PU200_nonprompt->SetBranchAddress("muon_time_",      &muon_time_);
  ch_PU200_nonprompt->SetBranchAddress("muon_time_err_",  &muon_time_err_);
  ch_PU200_nonprompt->SetBranchAddress("muon_prompt_",    &muon_prompt_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isBarrel_",  &muon_isBarrel_);
  ch_PU200_nonprompt->SetBranchAddress("muon_status_",    &muon_status_);
  ch_PU200_nonprompt->SetBranchAddress("muon_pv_dz_",     &muon_pv_dz_);
  ch_PU200_nonprompt->SetBranchAddress("muon_pv_dxy_",    &muon_pv_dxy_);
  ch_PU200_nonprompt->SetBranchAddress("muon_vz_",        &muon_vz_);
  ch_PU200_nonprompt->SetBranchAddress("muon_vz_",        &muon_vz_);
  ch_PU200_nonprompt->SetBranchAddress("muon_PVweight_",  &muon_PVweight_);
  ch_PU200_nonprompt->SetBranchAddress("track_pt_",       &track_pt_);
  ch_PU200_nonprompt->SetBranchAddress("track_time_",     &track_time_);
  ch_PU200_nonprompt->SetBranchAddress("track_time_err_", &track_time_err_);
  ch_PU200_nonprompt->SetBranchAddress("track_vz_",       &track_vz_);
  ch_PU200_nonprompt->SetBranchAddress("track_pv_dz_",    &track_pv_dz_);
  ch_PU200_nonprompt->SetBranchAddress("track_PVweight_", &track_PVweight_);
  ch_PU200_nonprompt->SetBranchAddress("track_type_",     &track_type_);
  ch_PU200_nonprompt->SetBranchAddress("selectedLV_",     &selectedLV_);
  ch_PU200_nonprompt->SetBranchAddress("match_vtx_sim2reco_",     &match_vtx_sim2reco_);
  ch_PU200_nonprompt->SetBranchAddress("match_vtx_reco2sim_",     &match_vtx_reco2sim_);
  // noPU
  ch_noPU_prompt->SetBranchAddress("muon_pt_",        &muon_pt_);
  ch_noPU_prompt->SetBranchAddress("muon_time_",      &muon_time_);
  ch_noPU_prompt->SetBranchAddress("muon_time_err_",  &muon_time_err_);
  ch_noPU_prompt->SetBranchAddress("muon_prompt_",    &muon_prompt_);
  ch_noPU_prompt->SetBranchAddress("muon_isBarrel_",  &muon_isBarrel_);
  ch_noPU_prompt->SetBranchAddress("muon_status_",    &muon_status_);
  ch_noPU_prompt->SetBranchAddress("muon_pv_dz_",     &muon_pv_dz_);
  ch_noPU_prompt->SetBranchAddress("muon_pv_dxy_",    &muon_pv_dxy_);
  ch_noPU_prompt->SetBranchAddress("muon_vz_",        &muon_vz_);
  ch_noPU_prompt->SetBranchAddress("muon_PVweight_",  &muon_PVweight_);
  ch_noPU_prompt->SetBranchAddress("track_pt_",       &track_pt_);
  ch_noPU_prompt->SetBranchAddress("track_time_",     &track_time_);
  ch_noPU_prompt->SetBranchAddress("track_time_err_", &track_time_err_);
  ch_noPU_prompt->SetBranchAddress("track_vz_",       &track_vz_);
  ch_noPU_prompt->SetBranchAddress("track_pv_dz_",    &track_pv_dz_);
  ch_noPU_prompt->SetBranchAddress("track_PVweight_", &track_PVweight_);
  ch_noPU_prompt->SetBranchAddress("track_type_",     &track_type_);
  ch_noPU_prompt->SetBranchAddress("selectedLV_",     &selectedLV_);
  ch_noPU_prompt->SetBranchAddress("match_vtx_sim2reco_",     &match_vtx_sim2reco_);
  ch_noPU_prompt->SetBranchAddress("match_vtx_reco2sim_",     &match_vtx_reco2sim_);
  ch_noPU_nonprompt->SetBranchAddress("muon_pt_",        &muon_pt_);
  ch_noPU_nonprompt->SetBranchAddress("muon_time_",      &muon_time_);
  ch_noPU_nonprompt->SetBranchAddress("muon_time_err_",  &muon_time_err_);
  ch_noPU_nonprompt->SetBranchAddress("muon_prompt_",    &muon_prompt_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isBarrel_",  &muon_isBarrel_);
  ch_noPU_nonprompt->SetBranchAddress("muon_status_",    &muon_status_);
  ch_noPU_nonprompt->SetBranchAddress("muon_pv_dz_",     &muon_pv_dz_);
  ch_noPU_nonprompt->SetBranchAddress("muon_pv_dxy_",    &muon_pv_dxy_);
  ch_noPU_nonprompt->SetBranchAddress("muon_vz_",        &muon_vz_);
  ch_noPU_nonprompt->SetBranchAddress("muon_PVweight_",  &muon_PVweight_);
  ch_noPU_nonprompt->SetBranchAddress("track_pt_",       &track_pt_);
  ch_noPU_nonprompt->SetBranchAddress("track_time_",     &track_time_);
  ch_noPU_nonprompt->SetBranchAddress("track_time_err_", &track_time_err_);
  ch_noPU_nonprompt->SetBranchAddress("track_vz_",       &track_vz_);
  ch_noPU_nonprompt->SetBranchAddress("track_pv_dz_",    &track_pv_dz_);
  ch_noPU_nonprompt->SetBranchAddress("track_PVweight_", &track_PVweight_);
  ch_noPU_nonprompt->SetBranchAddress("track_type_",     &track_type_);
  ch_noPU_nonprompt->SetBranchAddress("selectedLV_",     &selectedLV_);
  ch_noPU_nonprompt->SetBranchAddress("match_vtx_sim2reco_",     &match_vtx_sim2reco_);
  ch_noPU_nonprompt->SetBranchAddress("match_vtx_reco2sim_",     &match_vtx_reco2sim_);


  ///////////////////////
  // Define histograms //
  ///////////////////////

  int nbin=10;

  // PU200
    // prompt
  TH1D* h_PU200_prompt_PV_EB         = new TH1D("h_PU200_prompt_PV_EB",           "h_PU200_prompt_PV_EB",           nbin, 0, 5);
  TH1D* h_PU200_prompt_SV_EB         = new TH1D("h_PU200_prompt_SV_EB",           "h_PU200_prompt_SV_EB",           nbin, 0, 5);
  TH1D* h_PU200_prompt_PU_EB         = new TH1D("h_PU200_prompt_PU_EB",           "h_PU200_prompt_PU_EB",           nbin, 0, 5);
  TH1D* h_PU200_prompt_fake_EB       = new TH1D("h_PU200_prompt_fake_EB",         "h_PU200_prompt_fake_EB",         nbin, 0, 5);
  TH1D* h_PU200_prompt_no_tErr_EB    = new TH1D("h_PU200_prompt_no_tErr_EB",      "h_PU200_prompt_no_tErr_EB",      nbin, 0, 5);
  TH1D* h_PU200_prompt_EE            = new TH1D("h_PU200_prompt_EE",              "h_PU200_prompt_EE",              nbin, 0, 4);
  TH1D* h_PU200_prompt_SV_EE         = new TH1D("h_PU200_prompt_SV_EE",           "h_PU200_prompt_SV_EE",           nbin, 0, 5);
  TH1D* h_PU200_prompt_PU_EE         = new TH1D("h_PU200_prompt_PU_EE",           "h_PU200_prompt_PU_EE",           nbin, 0, 5);
  TH1D* h_PU200_prompt_fake_EE       = new TH1D("h_PU200_prompt_fake_EE",         "h_PU200_prompt_fake_EE",         nbin, 0, 5);
  TH1D* h_PU200_prompt_no_tErr_EE    = new TH1D("h_PU200_prompt_no_tErr_EE",      "h_PU200_prompt_no_tErr_EE",      nbin, 0, 5);
  float sumiso_PU200_prompt_EB=0;
  float sumiso_PU200_prompt_EE=0;
  float reliso_PU200_prompt_EB=0;
  float reliso_PU200_prompt_EE=0;
    // nonprompt
  TH1D* h_PU200_nonprompt_PV_EB         = new TH1D("h_PU200_nonprompt_PV_EB",           "h_PU200_nonprompt_PV_EB",           nbin, 0, 5);
  TH1D* h_PU200_nonprompt_SV_EB         = new TH1D("h_PU200_nonprompt_SV_EB",           "h_PU200_nonprompt_SV_EB",           nbin, 0, 5);
  TH1D* h_PU200_nonprompt_PU_EB         = new TH1D("h_PU200_nonprompt_PU_EB",           "h_PU200_nonprompt_PU_EB",           nbin, 0, 5);
  TH1D* h_PU200_nonprompt_fake_EB       = new TH1D("h_PU200_nonprompt_fake_EB",         "h_PU200_nonprompt_fake_EB",         nbin, 0, 5);
  TH1D* h_PU200_nonprompt_no_tErr_EB    = new TH1D("h_PU200_nonprompt_no_tErr_EB",      "h_PU200_nonprompt_no_tErr_EB",      nbin, 0, 5);
  TH1D* h_PU200_nonprompt_EE            = new TH1D("h_PU200_nonprompt_EE",              "h_PU200_nonprompt_EE",              nbin, 0, 4);
  TH1D* h_PU200_nonprompt_SV_EE         = new TH1D("h_PU200_nonprompt_SV_EE",           "h_PU200_nonprompt_SV_EE",           nbin, 0, 5);
  TH1D* h_PU200_nonprompt_PU_EE         = new TH1D("h_PU200_nonprompt_PU_EE",           "h_PU200_nonprompt_PU_EE",           nbin, 0, 5);
  TH1D* h_PU200_nonprompt_fake_EE       = new TH1D("h_PU200_nonprompt_fake_EE",         "h_PU200_nonprompt_fake_EE",         nbin, 0, 5);
  TH1D* h_PU200_nonprompt_no_tErr_EE    = new TH1D("h_PU200_nonprompt_no_tErr_EE",      "h_PU200_nonprompt_no_tErr_EE",      nbin, 0, 5);
  float sumiso_PU200_nonprompt_EB=0;
  float sumiso_PU200_nonprompt_EE=0;
  float reliso_PU200_nonprompt_EB=0;
  float reliso_PU200_nonprompt_EE=0;

	
  // noPU
    // prompt
  TH1D* h_noPU_prompt_PV_EB         = new TH1D("h_noPU_prompt_PV_EB",           "h_noPU_prompt_PV_EB",           nbin, 0, 5);
  TH1D* h_noPU_prompt_SV_EB         = new TH1D("h_noPU_prompt_SV_EB",           "h_noPU_prompt_SV_EB",           nbin, 0, 5);
  TH1D* h_noPU_prompt_PU_EB         = new TH1D("h_noPU_prompt_PU_EB",           "h_noPU_prompt_PU_EB",           nbin, 0, 5);
  TH1D* h_noPU_prompt_fake_EB       = new TH1D("h_noPU_prompt_fake_EB",         "h_noPU_prompt_fake_EB",         nbin, 0, 5);
  TH1D* h_noPU_prompt_no_tErr_EB    = new TH1D("h_noPU_prompt_no_tErr_EB",      "h_noPU_prompt_no_tErr_EB",      nbin, 0, 5);
  TH1D* h_noPU_prompt_EE            = new TH1D("h_noPU_prompt_EE",              "h_noPU_prompt_EE",              nbin, 0, 4);
  TH1D* h_noPU_prompt_SV_EE         = new TH1D("h_noPU_prompt_SV_EE",           "h_noPU_prompt_SV_EE",           nbin, 0, 5);
  TH1D* h_noPU_prompt_PU_EE         = new TH1D("h_noPU_prompt_PU_EE",           "h_noPU_prompt_PU_EE",           nbin, 0, 5);
  TH1D* h_noPU_prompt_fake_EE       = new TH1D("h_noPU_prompt_fake_EE",         "h_noPU_prompt_fake_EE",         nbin, 0, 5);
  TH1D* h_noPU_prompt_no_tErr_EE    = new TH1D("h_noPU_prompt_no_tErr_EE",      "h_noPU_prompt_no_tErr_EE",      nbin, 0, 5);
  float sumiso_noPU_prompt_EB=0;
  float sumiso_noPU_prompt_EE=0;
  float reliso_noPU_prompt_EB=0;
  float reliso_noPU_prompt_EE=0;
    // nonprompt
  TH1D* h_noPU_nonprompt_PV_EB         = new TH1D("h_noPU_nonprompt_PV_EB",           "h_noPU_nonprompt_PV_EB",           nbin, 0, 5);
  TH1D* h_noPU_nonprompt_SV_EB         = new TH1D("h_noPU_nonprompt_SV_EB",           "h_noPU_nonprompt_SV_EB",           nbin, 0, 5);
  TH1D* h_noPU_nonprompt_PU_EB         = new TH1D("h_noPU_nonprompt_PU_EB",           "h_noPU_nonprompt_PU_EB",           nbin, 0, 5);
  TH1D* h_noPU_nonprompt_fake_EB       = new TH1D("h_noPU_nonprompt_fake_EB",         "h_noPU_nonprompt_fake_EB",         nbin, 0, 5);
  TH1D* h_noPU_nonprompt_no_tErr_EB    = new TH1D("h_noPU_nonprompt_no_tErr_EB",      "h_noPU_nonprompt_no_tErr_EB",      nbin, 0, 5);
  TH1D* h_noPU_nonprompt_EE            = new TH1D("h_noPU_nonprompt_EE",              "h_noPU_nonprompt_EE",              nbin, 0, 4);
  TH1D* h_noPU_nonprompt_SV_EE         = new TH1D("h_noPU_nonprompt_SV_EE",           "h_noPU_nonprompt_SV_EE",           nbin, 0, 5);
  TH1D* h_noPU_nonprompt_PU_EE         = new TH1D("h_noPU_nonprompt_PU_EE",           "h_noPU_nonprompt_PU_EE",           nbin, 0, 5);
  TH1D* h_noPU_nonprompt_fake_EE       = new TH1D("h_noPU_nonprompt_fake_EE",         "h_noPU_nonprompt_fake_EE",         nbin, 0, 5);
  TH1D* h_noPU_nonprompt_no_tErr_EE    = new TH1D("h_noPU_nonprompt_no_tErr_EE",      "h_noPU_nonprompt_no_tErr_EE",      nbin, 0, 5);
  float sumiso_noPU_nonprompt_EB=0;
  float sumiso_noPU_nonprompt_EE=0;
  float reliso_noPU_nonprompt_EB=0;
  float reliso_noPU_nonprompt_EE=0;



  float dtsig_cut=0, dt_cut=0;

  //////////////////////
  //// PU200 prompt ////
  //////////////////////

  int n_status_failed_PU200_prompt_EB=0, n_status_failed_PU200_prompt_EE=0;
  int n_muon_PU200_prompt_EB=0, n_muon_PU200_prompt_EE=0;

  // event loop
  for(int ievt=0; ievt<ch_PU200_prompt->GetEntries(); ievt++) {
	ch_PU200_prompt->GetEntry(ievt);



	// muon loop
	for(int im=0; im<muon_pt_->size(); im++) {

	  if(muon_pv_dz_->at(im)>0.5 || muon_pv_dxy_->at(im)>0.2) continue;

	  //////////////////////////
	  // prompt muon - Barrel //
	  //////////////////////////
	  if(muon_prompt_->at(im)==1 && muon_isBarrel_->at(im)==1) {
		n_muon_PU200_prompt_EB++;
	    if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
		  n_status_failed_PU200_prompt_EB++;
		  continue;
	    }
//		if(muon_PVweight_->at(im)==0) continue;

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {

		  if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) >= 0.1) continue;

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }

		  // noMTD case
		  sumiso_PU200_prompt_EB += track_pt_->at(im).at(it);
		  // MTD case
//		  if(track_type_



		  if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_PU200_prompt_1sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_PU200_prompt_2sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_PU200_prompt_3sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_PU200_prompt_4sigma_EB += track_pt_->at(im).at(it);
		  if( dt_cut<0.12  &&   dt_cut>0 ) sumiso_PU200_prompt_40_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.18  &&   dt_cut>0 ) sumiso_PU200_prompt_60_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.24  &&   dt_cut>0 ) sumiso_PU200_prompt_80_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.30  &&   dt_cut>0 ) sumiso_PU200_prompt_100_EB    += track_pt_->at(im).at(it);
		  else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
            sumiso_PU200_prompt_1sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_2sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_3sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_4sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_40_EB     += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_60_EB     += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_80_EB     += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_100_EB    += track_pt_->at(im).at(it);
		  }
		  dtsig_cut=0;
		  dt_cut=0;
		}
		reliso_PU200_prompt_EB        = sumiso_PU200_prompt_EB/muon_pt_->at(im);
		reliso_PU200_prompt_1sigma_EB = sumiso_PU200_prompt_1sigma_EB/muon_pt_->at(im);
		reliso_PU200_prompt_2sigma_EB = sumiso_PU200_prompt_2sigma_EB/muon_pt_->at(im);
		reliso_PU200_prompt_3sigma_EB = sumiso_PU200_prompt_3sigma_EB/muon_pt_->at(im);
		reliso_PU200_prompt_4sigma_EB = sumiso_PU200_prompt_4sigma_EB/muon_pt_->at(im);
		reliso_PU200_prompt_40_EB     = sumiso_PU200_prompt_40_EB/muon_pt_->at(im);
		reliso_PU200_prompt_60_EB     = sumiso_PU200_prompt_60_EB/muon_pt_->at(im);
		reliso_PU200_prompt_80_EB     = sumiso_PU200_prompt_80_EB/muon_pt_->at(im);
		reliso_PU200_prompt_100_EB    = sumiso_PU200_prompt_100_EB/muon_pt_->at(im);

		// Store value of relative isolation
		h_PU200_prompt_EB->Fill(reliso_PU200_prompt_EB);
		h_PU200_prompt_1sigma_EB->Fill(reliso_PU200_prompt_1sigma_EB);
		h_PU200_prompt_2sigma_EB->Fill(reliso_PU200_prompt_2sigma_EB);
		h_PU200_prompt_3sigma_EB->Fill(reliso_PU200_prompt_3sigma_EB);
		h_PU200_prompt_4sigma_EB->Fill(reliso_PU200_prompt_4sigma_EB);
		h_PU200_prompt_40_EB->Fill(reliso_PU200_prompt_40_EB);
		h_PU200_prompt_60_EB->Fill(reliso_PU200_prompt_60_EB);
		h_PU200_prompt_80_EB->Fill(reliso_PU200_prompt_80_EB);
		h_PU200_prompt_100_EB->Fill(reliso_PU200_prompt_100_EB);

		// Initialize
        sumiso_PU200_prompt_EB=0, sumiso_PU200_prompt_1sigma_EB=0, sumiso_PU200_prompt_2sigma_EB=0, sumiso_PU200_prompt_3sigma_EB=0, sumiso_PU200_prompt_4sigma_EB=0, sumiso_PU200_prompt_40_EB=0, sumiso_PU200_prompt_60_EB=0, sumiso_PU200_prompt_80_EB=0, sumiso_PU200_prompt_100_EB=0;
        reliso_PU200_prompt_EB=0, reliso_PU200_prompt_1sigma_EB=0, reliso_PU200_prompt_2sigma_EB=0, reliso_PU200_prompt_3sigma_EB=0, reliso_PU200_prompt_4sigma_EB=0, reliso_PU200_prompt_40_EB=0, reliso_PU200_prompt_60_EB=0, reliso_PU200_prompt_80_EB=0, reliso_PU200_prompt_100_EB=0;
	  }

	  //////////////////////////
	  // prompt muon - Endcap //
	  //////////////////////////
	  else if(muon_prompt_->at(im)==1 && muon_isBarrel_->at(im)==0) {
		n_muon_PU200_prompt_EE++;
	    if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
		  n_status_failed_PU200_prompt_EE++;
		  continue;
	    }
//		if(muon_PVweight_->at(im)==0) continue;

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {

		  if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) >= 0.2) continue;

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  // noMTD case
		  sumiso_PU200_prompt_EE += track_pt_->at(im).at(it);
		  // MTD case
		  if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_PU200_prompt_1sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_PU200_prompt_2sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_PU200_prompt_3sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_PU200_prompt_4sigma_EE += track_pt_->at(im).at(it);
		  if( dt_cut<0.12  &&   dt_cut>0 ) sumiso_PU200_prompt_40_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.18  &&   dt_cut>0 ) sumiso_PU200_prompt_60_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.24  &&   dt_cut>0 ) sumiso_PU200_prompt_80_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.30  &&   dt_cut>0 ) sumiso_PU200_prompt_100_EE    += track_pt_->at(im).at(it);
		  else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
            sumiso_PU200_prompt_1sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_2sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_3sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_4sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_40_EE     += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_60_EE     += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_80_EE     += track_pt_->at(im).at(it);
            sumiso_PU200_prompt_100_EE    += track_pt_->at(im).at(it);
		  }
		  dtsig_cut=0;
		  dt_cut=0;
		}
		reliso_PU200_prompt_EE        = sumiso_PU200_prompt_EE/muon_pt_->at(im);
		reliso_PU200_prompt_1sigma_EE = sumiso_PU200_prompt_1sigma_EE/muon_pt_->at(im);
		reliso_PU200_prompt_2sigma_EE = sumiso_PU200_prompt_2sigma_EE/muon_pt_->at(im);
		reliso_PU200_prompt_3sigma_EE = sumiso_PU200_prompt_3sigma_EE/muon_pt_->at(im);
		reliso_PU200_prompt_4sigma_EE = sumiso_PU200_prompt_4sigma_EE/muon_pt_->at(im);
		reliso_PU200_prompt_40_EE     = sumiso_PU200_prompt_40_EE/muon_pt_->at(im);
		reliso_PU200_prompt_60_EE     = sumiso_PU200_prompt_60_EE/muon_pt_->at(im);
		reliso_PU200_prompt_80_EE     = sumiso_PU200_prompt_80_EE/muon_pt_->at(im);
		reliso_PU200_prompt_100_EE    = sumiso_PU200_prompt_100_EE/muon_pt_->at(im);

		// Store value of relative isolation
		h_PU200_prompt_EE->Fill(reliso_PU200_prompt_EE);
		h_PU200_prompt_1sigma_EE->Fill(reliso_PU200_prompt_1sigma_EE);
		h_PU200_prompt_2sigma_EE->Fill(reliso_PU200_prompt_2sigma_EE);
		h_PU200_prompt_3sigma_EE->Fill(reliso_PU200_prompt_3sigma_EE);
		h_PU200_prompt_4sigma_EE->Fill(reliso_PU200_prompt_4sigma_EE);
		h_PU200_prompt_40_EE->Fill(reliso_PU200_prompt_40_EE);
		h_PU200_prompt_60_EE->Fill(reliso_PU200_prompt_60_EE);
		h_PU200_prompt_80_EE->Fill(reliso_PU200_prompt_80_EE);
		h_PU200_prompt_100_EE->Fill(reliso_PU200_prompt_100_EE);

		// Initialize
        sumiso_PU200_prompt_EE=0, sumiso_PU200_prompt_1sigma_EE=0, sumiso_PU200_prompt_2sigma_EE=0, sumiso_PU200_prompt_3sigma_EE=0, sumiso_PU200_prompt_4sigma_EE=0, sumiso_PU200_prompt_40_EE=0, sumiso_PU200_prompt_60_EE=0, sumiso_PU200_prompt_80_EE=0, sumiso_PU200_prompt_100_EE=0;
        reliso_PU200_prompt_EE=0, reliso_PU200_prompt_1sigma_EE=0, reliso_PU200_prompt_2sigma_EE=0, reliso_PU200_prompt_3sigma_EE=0, reliso_PU200_prompt_4sigma_EE=0, reliso_PU200_prompt_40_EE=0, reliso_PU200_prompt_60_EE=0, reliso_PU200_prompt_80_EE=0, reliso_PU200_prompt_100_EE=0;
	  }
	} // End of muon loop
  } // End of event loop


  /////////////////////////
  //// PU200 nonprompt ////
  /////////////////////////
  
  int n_status_failed_PU200_nonprompt_EB=0, n_status_failed_PU200_nonprompt_EE=0;
  int n_muon_PU200_nonprompt_EB=0, n_muon_PU200_nonprompt_EE=0;

  // event loop
  for(int ievt=0; ievt<ch_PU200_nonprompt->GetEntries(); ievt++) {
	ch_PU200_nonprompt->GetEntry(ievt);

	// muon loop
	for(int im=0; im<muon_pt_->size(); im++) {

	  if(muon_pv_dz_->at(im)>0.5 || muon_pv_dxy_->at(im)>0.2) continue;

	  /////////////////////////////
	  // nonprompt muon - Barrel //
	  /////////////////////////////
	  if(muon_prompt_->at(im)==0 && muon_isBarrel_->at(im)==1) {
		n_muon_PU200_nonprompt_EB++;
	    if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
		  n_status_failed_PU200_nonprompt_EB++;
		  continue;
	    }
//		if(muon_PVweight_->at(im)==0) continue;

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {

		  if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) >= 0.1) continue;

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  // noMTD case
		  sumiso_PU200_nonprompt_EB += track_pt_->at(im).at(it);
		  // MTD case
		  if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_PU200_nonprompt_1sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_PU200_nonprompt_2sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_PU200_nonprompt_3sigma_EB += track_pt_->at(im).at(it);
		  if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_PU200_nonprompt_4sigma_EB += track_pt_->at(im).at(it);
		  if( dt_cut<0.12  &&   dt_cut>0 ) sumiso_PU200_nonprompt_40_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.18  &&   dt_cut>0 ) sumiso_PU200_nonprompt_60_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.24  &&   dt_cut>0 ) sumiso_PU200_nonprompt_80_EB     += track_pt_->at(im).at(it);
		  if( dt_cut<0.30  &&   dt_cut>0 ) sumiso_PU200_nonprompt_100_EB    += track_pt_->at(im).at(it);
		  else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
            sumiso_PU200_nonprompt_1sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_2sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_3sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_4sigma_EB += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_40_EB     += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_60_EB     += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_80_EB     += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_100_EB    += track_pt_->at(im).at(it);
		  }
		  dtsig_cut=0;
		  dt_cut=0;
		}
		reliso_PU200_nonprompt_EB        = sumiso_PU200_nonprompt_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_1sigma_EB = sumiso_PU200_nonprompt_1sigma_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_2sigma_EB = sumiso_PU200_nonprompt_2sigma_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_3sigma_EB = sumiso_PU200_nonprompt_3sigma_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_4sigma_EB = sumiso_PU200_nonprompt_4sigma_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_40_EB     = sumiso_PU200_nonprompt_40_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_60_EB     = sumiso_PU200_nonprompt_60_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_80_EB     = sumiso_PU200_nonprompt_80_EB/muon_pt_->at(im);
		reliso_PU200_nonprompt_100_EB    = sumiso_PU200_nonprompt_100_EB/muon_pt_->at(im);

		// Store value of relative isolation
		h_PU200_nonprompt_EB->Fill(reliso_PU200_nonprompt_EB);
		h_PU200_nonprompt_1sigma_EB->Fill(reliso_PU200_nonprompt_1sigma_EB);
		h_PU200_nonprompt_2sigma_EB->Fill(reliso_PU200_nonprompt_2sigma_EB);
		h_PU200_nonprompt_3sigma_EB->Fill(reliso_PU200_nonprompt_3sigma_EB);
		h_PU200_nonprompt_4sigma_EB->Fill(reliso_PU200_nonprompt_4sigma_EB);
		h_PU200_nonprompt_40_EB->Fill(reliso_PU200_nonprompt_40_EB);
		h_PU200_nonprompt_60_EB->Fill(reliso_PU200_nonprompt_60_EB);
		h_PU200_nonprompt_80_EB->Fill(reliso_PU200_nonprompt_80_EB);
		h_PU200_nonprompt_100_EB->Fill(reliso_PU200_nonprompt_100_EB);

		// Initialize
        sumiso_PU200_nonprompt_EB=0, sumiso_PU200_nonprompt_1sigma_EB=0, sumiso_PU200_nonprompt_2sigma_EB=0, sumiso_PU200_nonprompt_3sigma_EB=0, sumiso_PU200_nonprompt_4sigma_EB=0, sumiso_PU200_nonprompt_40_EB=0, sumiso_PU200_nonprompt_60_EB=0, sumiso_PU200_nonprompt_80_EB=0, sumiso_PU200_nonprompt_100_EB=0;
        reliso_PU200_nonprompt_EB=0, reliso_PU200_nonprompt_1sigma_EB=0, reliso_PU200_nonprompt_2sigma_EB=0, reliso_PU200_nonprompt_3sigma_EB=0, reliso_PU200_nonprompt_4sigma_EB=0, reliso_PU200_nonprompt_40_EB=0, reliso_PU200_nonprompt_60_EB=0, reliso_PU200_nonprompt_80_EB=0, reliso_PU200_nonprompt_100_EB=0;
	  }

	  /////////////////////////////
	  // nonprompt muon - Endcap //
	  /////////////////////////////
	  else if(muon_prompt_->at(im)==0 && muon_isBarrel_->at(im)==0) {
		n_muon_PU200_nonprompt_EE++;
	    if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
		  n_status_failed_PU200_nonprompt_EE++;
		  continue;
	    }
//		if(muon_PVweight_->at(im)==0) continue;

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {

		  if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) >= 0.2) continue;

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  // noMTD case
		  sumiso_PU200_nonprompt_EE += track_pt_->at(im).at(it);
		  // MTD case
		  if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_PU200_nonprompt_1sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_PU200_nonprompt_2sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_PU200_nonprompt_3sigma_EE += track_pt_->at(im).at(it);
		  if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_PU200_nonprompt_4sigma_EE += track_pt_->at(im).at(it);
		  if( dt_cut<0.12  &&   dt_cut>0 ) sumiso_PU200_nonprompt_40_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.18  &&   dt_cut>0 ) sumiso_PU200_nonprompt_60_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.24  &&   dt_cut>0 ) sumiso_PU200_nonprompt_80_EE     += track_pt_->at(im).at(it);
		  if( dt_cut<0.30  &&   dt_cut>0 ) sumiso_PU200_nonprompt_100_EE    += track_pt_->at(im).at(it);
		  else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
            sumiso_PU200_nonprompt_1sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_2sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_3sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_4sigma_EE += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_40_EE     += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_60_EE     += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_80_EE     += track_pt_->at(im).at(it);
            sumiso_PU200_nonprompt_100_EE    += track_pt_->at(im).at(it);
		  }
		  dtsig_cut=0;
		  dt_cut=0;
		}
		reliso_PU200_nonprompt_EE        = sumiso_PU200_nonprompt_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_1sigma_EE = sumiso_PU200_nonprompt_1sigma_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_2sigma_EE = sumiso_PU200_nonprompt_2sigma_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_3sigma_EE = sumiso_PU200_nonprompt_3sigma_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_4sigma_EE = sumiso_PU200_nonprompt_4sigma_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_40_EE     = sumiso_PU200_nonprompt_40_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_60_EE     = sumiso_PU200_nonprompt_60_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_80_EE     = sumiso_PU200_nonprompt_80_EE/muon_pt_->at(im);
		reliso_PU200_nonprompt_100_EE    = sumiso_PU200_nonprompt_100_EE/muon_pt_->at(im);

		// Store value of relative isolation
		h_PU200_nonprompt_EE->Fill(reliso_PU200_nonprompt_EE);
		h_PU200_nonprompt_1sigma_EE->Fill(reliso_PU200_nonprompt_1sigma_EE);
		h_PU200_nonprompt_2sigma_EE->Fill(reliso_PU200_nonprompt_2sigma_EE);
		h_PU200_nonprompt_3sigma_EE->Fill(reliso_PU200_nonprompt_3sigma_EE);
		h_PU200_nonprompt_4sigma_EE->Fill(reliso_PU200_nonprompt_4sigma_EE);
		h_PU200_nonprompt_40_EE->Fill(reliso_PU200_nonprompt_40_EE);
		h_PU200_nonprompt_60_EE->Fill(reliso_PU200_nonprompt_60_EE);
		h_PU200_nonprompt_80_EE->Fill(reliso_PU200_nonprompt_80_EE);
		h_PU200_nonprompt_100_EE->Fill(reliso_PU200_nonprompt_100_EE);

		// Initialize
        sumiso_PU200_nonprompt_EE=0, sumiso_PU200_nonprompt_1sigma_EE=0, sumiso_PU200_nonprompt_2sigma_EE=0, sumiso_PU200_nonprompt_3sigma_EE=0, sumiso_PU200_nonprompt_4sigma_EE=0, sumiso_PU200_nonprompt_40_EE=0, sumiso_PU200_nonprompt_60_EE=0, sumiso_PU200_nonprompt_80_EE=0, sumiso_PU200_nonprompt_100_EE=0;
        reliso_PU200_nonprompt_EE=0, reliso_PU200_nonprompt_1sigma_EE=0, reliso_PU200_nonprompt_2sigma_EE=0, reliso_PU200_nonprompt_3sigma_EE=0, reliso_PU200_nonprompt_4sigma_EE=0, reliso_PU200_nonprompt_40_EE=0, reliso_PU200_nonprompt_60_EE=0, reliso_PU200_nonprompt_80_EE=0, reliso_PU200_nonprompt_100_EE=0;

	  }
	} // End of muon loop
  } // End of event loop
	
}

*/


int main(int argc, char **argv)
{
  
  draw_iso_efficiency_ntuple();
//  draw_track_type_sigma_ntuple();

  return 0;
}



















































