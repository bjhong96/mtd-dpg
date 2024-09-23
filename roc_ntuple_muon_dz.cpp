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

void draw_iso_efficiency_ntuple_dz();
void draw_track_type_sigma_ntuple();

TString path_PU200_prompt        = "240801/harvester_PU200_prompt_muon.root";
TString path_noPU_prompt         = "240801/harvester_noPU_prompt_muon.root";
TString path_PU200_prompt_vtx    = "240801/harvester_PU200_prompt_muon_vtx.root";
TString path_noPU_prompt_vtx     = "240801/harvester_noPU_prompt_muon_vtx.root";
TString path_PU200_nonprompt     = "240801/harvester_PU200_nonprompt_muon_qcd.root";
TString path_noPU_nonprompt      = "240801/harvester_noPU_nonprompt_muon_qcd.root";
TString path_PU200_nonprompt_vtx = "240801/harvester_PU200_nonprompt_muon_qcd_vtx.root";
TString path_noPU_nonprompt_vtx  = "240801/harvester_noPU_nonprompt_muon_qcd_vtx.root";
//TString path_PU200_nonprompt     = "240801/harvester_PU200_nonprompt_muon_ttbar.root";
//TString path_noPU_nonprompt      = "240801/harvester_noPU_nonprompt_muon_ttbar.root";
//TString path_PU200_nonprompt_vtx = "240801/harvester_PU200_nonprompt_muon_ttbar_vtx.root";
//TString path_noPU_nonprompt_vtx  = "240801/harvester_noPU_nonprompt_muon_ttbar_vtx.root";

TString ntuple_PU200_prompt        = "240801/ntuple_PU200_prompt_muon.root";
TString ntuple_noPU_prompt         = "240801/ntuple_noPU_prompt_muon.root";
TString ntuple_PU200_prompt_vtx    = "240801/ntuple_PU200_prompt_muon_vtx.root";
TString ntuple_noPU_prompt_vtx     = "240801/ntuple_noPU_prompt_muon_vtx.root";
TString ntuple_PU200_nonprompt     = "240801/ntuple_PU200_nonprompt_muon_qcd.root";
TString ntuple_noPU_nonprompt      = "240801/ntuple_noPU_nonprompt_muon_qcd.root";
TString ntuple_PU200_nonprompt_vtx = "240801/ntuple_PU200_nonprompt_muon_qcd_vtx.root";
TString ntuple_noPU_nonprompt_vtx  = "240801/ntuple_noPU_nonprompt_muon_qcd_vtx.root";
//TString ntuple_PU200_nonprompt     = "240801/ntuple_PU200_nonprompt_muon_ttbar.root";
//TString ntuple_noPU_nonprompt      = "240801/ntuple_noPU_nonprompt_muon_ttbar.root";
//TString ntuple_PU200_nonprompt_vtx = "240801/ntuple_PU200_nonprompt_muon_ttbar_vtx.root";
//TString ntuple_noPU_nonprompt_vtx  = "240801/ntuple_noPU_nonprompt_muon_ttbar_vtx.root";

 





void draw_iso_efficiency_ntuple_dz() {

  // generate the dictionary for vector<vector<sth>> collection
  gInterpreter->GenerateDictionary("vector<vector<float> >",  "vector");
  gInterpreter->GenerateDictionary("vector<vector<double> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<int> >",    "vector");
  gInterpreter->GenerateDictionary("vector<vector<bool> >",   "vector");

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
  std::vector<std::vector<float>> *track_evtId_(0), *track_bx_(0);
  vector<int>   *muon_index_(0);
  std::vector<std::vector<int>> *track_index_(0);
  std::vector<std::vector<bool>> *track_genMatched_(0);
  std::vector<std::vector<bool>> *selectedLV_(0), *match_vtx_sim2reco_(0), *match_vtx_reco2sim_(0);
  int vtx_index_=999;
  int recovtx_sim_=999, simvtx_reco_=999;
  int simvtx_bx_=999, simvtx_evtId_=999;
  int recovtx_original_index_=999;
  vector<float> *vtx_time_(0), *vtx_time_err_(0);



  // PU200
  ch_PU200_prompt->SetBranchAddress("muon_pt_",        					&muon_pt_);
  ch_PU200_prompt->SetBranchAddress("muon_time_",      					&muon_time_);
  ch_PU200_prompt->SetBranchAddress("muon_time_err_",  					&muon_time_err_);
  ch_PU200_prompt->SetBranchAddress("muon_prompt_",    					&muon_prompt_);
  ch_PU200_prompt->SetBranchAddress("muon_isBarrel_",  					&muon_isBarrel_);
  ch_PU200_prompt->SetBranchAddress("muon_status_",    					&muon_status_);
  ch_PU200_prompt->SetBranchAddress("muon_pv_dz_",      				&muon_pv_dz_);
  ch_PU200_prompt->SetBranchAddress("muon_pv_dxy_",     				&muon_pv_dxy_);
  ch_PU200_prompt->SetBranchAddress("muon_vz_",       					&muon_vz_);
  ch_PU200_prompt->SetBranchAddress("muon_PVweight_", 					&muon_PVweight_);
  ch_PU200_prompt->SetBranchAddress("track_pt_",      					&track_pt_);
  ch_PU200_prompt->SetBranchAddress("track_time_",     					&track_time_);
  ch_PU200_prompt->SetBranchAddress("track_time_err_", 					&track_time_err_);
  ch_PU200_prompt->SetBranchAddress("track_vz_",       					&track_vz_);
  ch_PU200_prompt->SetBranchAddress("track_pv_dz_",    					&track_pv_dz_);
  ch_PU200_prompt->SetBranchAddress("track_PVweight_", 					&track_PVweight_);
  ch_PU200_prompt->SetBranchAddress("selectedLV_",     					&selectedLV_);
  ch_PU200_prompt->SetBranchAddress("match_vtx_sim2reco_",     			&match_vtx_sim2reco_);
  ch_PU200_prompt->SetBranchAddress("match_vtx_reco2sim_",     			&match_vtx_reco2sim_);
  ch_PU200_prompt->SetBranchAddress("vtx_index_",      					&vtx_index_);
  ch_PU200_prompt->SetBranchAddress("recovtx_sim_",    					&recovtx_sim_);
  ch_PU200_prompt->SetBranchAddress("simvtx_reco_",    					&simvtx_reco_);
  ch_PU200_prompt->SetBranchAddress("simvtx_bx_",      					&simvtx_bx_);
  ch_PU200_prompt->SetBranchAddress("simvtx_evtId_",   					&simvtx_evtId_);
  ch_PU200_prompt->SetBranchAddress("recovtx_original_index_", 			&recovtx_original_index_);
  ch_PU200_prompt->SetBranchAddress("muon_isMuon_",    					&muon_isMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isPFMuon_",    				&muon_isPFMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isGlobalMuon_",    			&muon_isGlobalMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isTrackerMuon_",    			&muon_isTrackerMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isStandAloneMuon_",    		&muon_isStandAloneMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isCutBasedIdLoose_",    		&muon_isCutBasedIdLoose_);
  ch_PU200_prompt->SetBranchAddress("muon_isLooseMuon_",    			&muon_isLooseMuon_);
  ch_PU200_prompt->SetBranchAddress("track_genMatched_",    			&track_genMatched_);
  ch_PU200_prompt->SetBranchAddress("track_evtId_",    					&track_evtId_);
  ch_PU200_prompt->SetBranchAddress("track_bx_",       					&track_bx_);
  ch_PU200_prompt->SetBranchAddress("vtx_time_",       					&vtx_time_);
  ch_PU200_prompt->SetBranchAddress("vtx_time_err_",       				&vtx_time_err_);

  ch_PU200_nonprompt->SetBranchAddress("muon_pt_",        					&muon_pt_);
  ch_PU200_nonprompt->SetBranchAddress("muon_time_",      					&muon_time_);
  ch_PU200_nonprompt->SetBranchAddress("muon_time_err_",  					&muon_time_err_);
  ch_PU200_nonprompt->SetBranchAddress("muon_prompt_",    					&muon_prompt_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isBarrel_",  					&muon_isBarrel_);
  ch_PU200_nonprompt->SetBranchAddress("muon_status_",    					&muon_status_);
  ch_PU200_nonprompt->SetBranchAddress("muon_pv_dz_",      					&muon_pv_dz_);
  ch_PU200_nonprompt->SetBranchAddress("muon_pv_dxy_",     					&muon_pv_dxy_);
  ch_PU200_nonprompt->SetBranchAddress("muon_vz_",       					&muon_vz_);
  ch_PU200_nonprompt->SetBranchAddress("muon_PVweight_", 					&muon_PVweight_);
  ch_PU200_nonprompt->SetBranchAddress("track_pt_",      					&track_pt_);
  ch_PU200_nonprompt->SetBranchAddress("track_time_",     					&track_time_);
  ch_PU200_nonprompt->SetBranchAddress("track_time_err_", 					&track_time_err_);
  ch_PU200_nonprompt->SetBranchAddress("track_vz_",       					&track_vz_);
  ch_PU200_nonprompt->SetBranchAddress("track_pv_dz_",    					&track_pv_dz_);
  ch_PU200_nonprompt->SetBranchAddress("track_PVweight_", 					&track_PVweight_);
  ch_PU200_nonprompt->SetBranchAddress("selectedLV_",     					&selectedLV_);
  ch_PU200_nonprompt->SetBranchAddress("match_vtx_sim2reco_",     			&match_vtx_sim2reco_);
  ch_PU200_nonprompt->SetBranchAddress("match_vtx_reco2sim_",     			&match_vtx_reco2sim_);
  ch_PU200_nonprompt->SetBranchAddress("vtx_index_",      					&vtx_index_);
  ch_PU200_nonprompt->SetBranchAddress("recovtx_sim_",    					&recovtx_sim_);
  ch_PU200_nonprompt->SetBranchAddress("simvtx_reco_",    					&simvtx_reco_);
  ch_PU200_nonprompt->SetBranchAddress("simvtx_bx_",      					&simvtx_bx_);
  ch_PU200_nonprompt->SetBranchAddress("simvtx_evtId_",   					&simvtx_evtId_);
  ch_PU200_nonprompt->SetBranchAddress("recovtx_original_index_", 			&recovtx_original_index_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isMuon_",    					&muon_isMuon_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isPFMuon_",    				&muon_isPFMuon_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isGlobalMuon_",    			&muon_isGlobalMuon_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isTrackerMuon_",    			&muon_isTrackerMuon_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isStandAloneMuon_",    		&muon_isStandAloneMuon_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isCutBasedIdLoose_",    		&muon_isCutBasedIdLoose_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isLooseMuon_",    				&muon_isLooseMuon_);
  ch_PU200_nonprompt->SetBranchAddress("track_genMatched_",    				&track_genMatched_);
  ch_PU200_nonprompt->SetBranchAddress("track_evtId_",    					&track_evtId_);
  ch_PU200_nonprompt->SetBranchAddress("track_bx_",       					&track_bx_);
  ch_PU200_nonprompt->SetBranchAddress("vtx_time_",       					&vtx_time_);
  ch_PU200_nonprompt->SetBranchAddress("vtx_time_err_",       				&vtx_time_err_);

  // noPU
  ch_noPU_prompt->SetBranchAddress("muon_pt_",        					&muon_pt_);
  ch_noPU_prompt->SetBranchAddress("muon_time_",      					&muon_time_);
  ch_noPU_prompt->SetBranchAddress("muon_time_err_",  					&muon_time_err_);
  ch_noPU_prompt->SetBranchAddress("muon_prompt_",    					&muon_prompt_);
  ch_noPU_prompt->SetBranchAddress("muon_isBarrel_",  					&muon_isBarrel_);
  ch_noPU_prompt->SetBranchAddress("muon_status_",    					&muon_status_);
  ch_noPU_prompt->SetBranchAddress("muon_pv_dz_",      					&muon_pv_dz_);
  ch_noPU_prompt->SetBranchAddress("muon_pv_dxy_",     					&muon_pv_dxy_);
  ch_noPU_prompt->SetBranchAddress("muon_vz_",       					&muon_vz_);
  ch_noPU_prompt->SetBranchAddress("muon_PVweight_", 					&muon_PVweight_);
  ch_noPU_prompt->SetBranchAddress("track_pt_",      					&track_pt_);
  ch_noPU_prompt->SetBranchAddress("track_time_",     					&track_time_);
  ch_noPU_prompt->SetBranchAddress("track_time_err_", 					&track_time_err_);
  ch_noPU_prompt->SetBranchAddress("track_vz_",       					&track_vz_);
  ch_noPU_prompt->SetBranchAddress("track_pv_dz_",    					&track_pv_dz_);
  ch_noPU_prompt->SetBranchAddress("track_PVweight_", 					&track_PVweight_);
  ch_noPU_prompt->SetBranchAddress("selectedLV_",     					&selectedLV_);
  ch_noPU_prompt->SetBranchAddress("match_vtx_sim2reco_",     			&match_vtx_sim2reco_);
  ch_noPU_prompt->SetBranchAddress("match_vtx_reco2sim_",     			&match_vtx_reco2sim_);
  ch_noPU_prompt->SetBranchAddress("vtx_index_",      					&vtx_index_);
  ch_noPU_prompt->SetBranchAddress("recovtx_sim_",    					&recovtx_sim_);
  ch_noPU_prompt->SetBranchAddress("simvtx_reco_",    					&simvtx_reco_);
  ch_noPU_prompt->SetBranchAddress("simvtx_bx_",      					&simvtx_bx_);
  ch_noPU_prompt->SetBranchAddress("simvtx_evtId_",   					&simvtx_evtId_);
  ch_noPU_prompt->SetBranchAddress("recovtx_original_index_", 			&recovtx_original_index_);
  ch_noPU_prompt->SetBranchAddress("muon_isMuon_",    					&muon_isMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isPFMuon_",    				&muon_isPFMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isGlobalMuon_",    			&muon_isGlobalMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isTrackerMuon_",    			&muon_isTrackerMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isStandAloneMuon_",    		&muon_isStandAloneMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isCutBasedIdLoose_",    		&muon_isCutBasedIdLoose_);
  ch_noPU_prompt->SetBranchAddress("muon_isLooseMuon_",    				&muon_isLooseMuon_);
  ch_noPU_prompt->SetBranchAddress("track_genMatched_",    				&track_genMatched_);
  ch_noPU_prompt->SetBranchAddress("track_evtId_",    					&track_evtId_);
  ch_noPU_prompt->SetBranchAddress("track_bx_",       					&track_bx_);
  ch_noPU_prompt->SetBranchAddress("vtx_time_",       					&vtx_time_);
  ch_noPU_prompt->SetBranchAddress("vtx_time_err_",       				&vtx_time_err_);

  ch_noPU_nonprompt->SetBranchAddress("muon_pt_",        					&muon_pt_);
  ch_noPU_nonprompt->SetBranchAddress("muon_time_",      					&muon_time_);
  ch_noPU_nonprompt->SetBranchAddress("muon_time_err_",  					&muon_time_err_);
  ch_noPU_nonprompt->SetBranchAddress("muon_prompt_",    					&muon_prompt_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isBarrel_",  					&muon_isBarrel_);
  ch_noPU_nonprompt->SetBranchAddress("muon_status_",    					&muon_status_);
  ch_noPU_nonprompt->SetBranchAddress("muon_pv_dz_",      					&muon_pv_dz_);
  ch_noPU_nonprompt->SetBranchAddress("muon_pv_dxy_",     					&muon_pv_dxy_);
  ch_noPU_nonprompt->SetBranchAddress("muon_vz_",       					&muon_vz_);
  ch_noPU_nonprompt->SetBranchAddress("muon_PVweight_", 					&muon_PVweight_);
  ch_noPU_nonprompt->SetBranchAddress("track_pt_",      					&track_pt_);
  ch_noPU_nonprompt->SetBranchAddress("track_time_",     					&track_time_);
  ch_noPU_nonprompt->SetBranchAddress("track_time_err_", 					&track_time_err_);
  ch_noPU_nonprompt->SetBranchAddress("track_vz_",       					&track_vz_);
  ch_noPU_nonprompt->SetBranchAddress("track_pv_dz_",    					&track_pv_dz_);
  ch_noPU_nonprompt->SetBranchAddress("track_PVweight_", 					&track_PVweight_);
  ch_noPU_nonprompt->SetBranchAddress("selectedLV_",     					&selectedLV_);
  ch_noPU_nonprompt->SetBranchAddress("match_vtx_sim2reco_",     			&match_vtx_sim2reco_);
  ch_noPU_nonprompt->SetBranchAddress("match_vtx_reco2sim_",     			&match_vtx_reco2sim_);
  ch_noPU_nonprompt->SetBranchAddress("vtx_index_",      					&vtx_index_);
  ch_noPU_nonprompt->SetBranchAddress("recovtx_sim_",    					&recovtx_sim_);
  ch_noPU_nonprompt->SetBranchAddress("simvtx_reco_",    					&simvtx_reco_);
  ch_noPU_nonprompt->SetBranchAddress("simvtx_bx_",      					&simvtx_bx_);
  ch_noPU_nonprompt->SetBranchAddress("simvtx_evtId_",   					&simvtx_evtId_);
  ch_noPU_nonprompt->SetBranchAddress("recovtx_original_index_", 			&recovtx_original_index_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isMuon_",    					&muon_isMuon_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isPFMuon_",    					&muon_isPFMuon_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isGlobalMuon_",    				&muon_isGlobalMuon_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isTrackerMuon_",    			&muon_isTrackerMuon_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isStandAloneMuon_",    			&muon_isStandAloneMuon_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isCutBasedIdLoose_",    		&muon_isCutBasedIdLoose_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isLooseMuon_",    				&muon_isLooseMuon_);
  ch_noPU_nonprompt->SetBranchAddress("track_genMatched_",    				&track_genMatched_);
  ch_noPU_nonprompt->SetBranchAddress("track_evtId_",    					&track_evtId_);
  ch_noPU_nonprompt->SetBranchAddress("track_bx_",       					&track_bx_);
  ch_noPU_nonprompt->SetBranchAddress("vtx_time_",       					&vtx_time_);
  ch_noPU_nonprompt->SetBranchAddress("vtx_time_err_",       				&vtx_time_err_);



  bool flag_isMuon            = false;
  bool flag_isPFMuon          = false;
  bool flag_isGlobalMuon      = false;
  bool flag_isTrackerMuon     = false;
  bool flag_isStandAloneMuon  = false;
  bool flag_isCutBasedIdLoose = false;
  bool flag_isLooseMuon       = true;

  bool flag_muon_status=true;
  bool flag_vtx_matching=true,  flag_vtx_selectedLV=true;
//  bool flag_track_pv_dz=false;

  bool flag_muon_pv_dz=false,     flag_muon_pv_dxy=false;

  bool flag_muon_PVweight=false;
  bool flag_track_PVweight=false;

  float dtsig_cut=0, dt_cut=0;
  float dtsig_cut_vtx=0, dt_cut_vtx=0;

  // PARAM
  float muon_PVweight_cut=0.000000000000001;
  float track_PVweight_cut=0.000000000000001;
  //float muon_pv_dz_cut_EB=0.5;
  //float muon_pv_dz_cut_EE=0.5;
  float muon_pv_dz_cut_EB=0.2;
  float muon_pv_dz_cut_EE=0.2;
  float muon_pv_dxy_cut_EB=0.2;
  float muon_pv_dxy_cut_EE=0.2;
  //vector<float> track_pv_dz_cut_EB={0.01, 0.05, 0.1, 0.5, 1, 5, 10000};
  //vector<float> track_pv_dz_cut_EE={0.01, 0.05, 0.1, 0.5, 1, 5, 10000};
  vector<float> track_pv_dz_cut_EB={0.1, 0.2, 0.3, 0.4, 0.5};
  vector<float> track_pv_dz_cut_EE={0.1, 0.2, 0.3, 0.4, 0.5};

  vector<int> num_gen_PU200_prompt_EB(track_pv_dz_cut_EB.size(), 0),    num_gen_notPV_PU200_prompt_EB(track_pv_dz_cut_EB.size(), 0),    num_gen_PU200_prompt_EE(track_pv_dz_cut_EE.size(), 0),    num_gen_notPV_PU200_prompt_EE(track_pv_dz_cut_EE.size(), 0);
  vector<int> num_gen_PU200_nonprompt_EB(track_pv_dz_cut_EB.size(), 0), num_gen_notPV_PU200_nonprompt_EB(track_pv_dz_cut_EB.size(), 0), num_gen_PU200_nonprompt_EE(track_pv_dz_cut_EE.size(), 0), num_gen_notPV_PU200_nonprompt_EE(track_pv_dz_cut_EE.size(), 0);
  vector<int> num_gen_noPU_prompt_EB(track_pv_dz_cut_EB.size(), 0),     num_gen_notPV_noPU_prompt_EB(track_pv_dz_cut_EB.size(), 0),     num_gen_noPU_prompt_EE(track_pv_dz_cut_EE.size(), 0),     num_gen_notPV_noPU_prompt_EE(track_pv_dz_cut_EE.size(), 0);
  vector<int> num_gen_noPU_nonprompt_EB(track_pv_dz_cut_EB.size(), 0),  num_gen_notPV_noPU_nonprompt_EB(track_pv_dz_cut_EB.size(), 0),  num_gen_noPU_nonprompt_EE(track_pv_dz_cut_EE.size(), 0),  num_gen_notPV_noPU_nonprompt_EE(track_pv_dz_cut_EE.size(), 0);

  ///////////////////////
  // Define histograms //
  ///////////////////////

  // PU200
    // prompt
	// (muon, track)
  vector<TH1D*> list_h_PU200_prompt_EB,            list_h_PU200_prompt_EE;
  vector<TH1D*> list_h_PU200_prompt_1sigma_EB,     list_h_PU200_prompt_1sigma_EE;
  vector<TH1D*> list_h_PU200_prompt_2sigma_EB,     list_h_PU200_prompt_2sigma_EE;
  vector<TH1D*> list_h_PU200_prompt_3sigma_EB,     list_h_PU200_prompt_3sigma_EE;
  vector<TH1D*> list_h_PU200_prompt_4sigma_EB,     list_h_PU200_prompt_4sigma_EE;
  vector<TH1D*> list_h_PU200_prompt_genMatched_EB, list_h_PU200_prompt_genMatched_EE;
	// (PV, track)
  vector<TH1D*> list_h_PU200_prompt_EB_vtx,            list_h_PU200_prompt_EE_vtx;
  vector<TH1D*> list_h_PU200_prompt_1sigma_EB_vtx,     list_h_PU200_prompt_1sigma_EE_vtx;
  vector<TH1D*> list_h_PU200_prompt_2sigma_EB_vtx,     list_h_PU200_prompt_2sigma_EE_vtx;
  vector<TH1D*> list_h_PU200_prompt_3sigma_EB_vtx,     list_h_PU200_prompt_3sigma_EE_vtx;
  vector<TH1D*> list_h_PU200_prompt_4sigma_EB_vtx,     list_h_PU200_prompt_4sigma_EE_vtx;
  vector<TH1D*> list_h_PU200_prompt_genMatched_EB_vtx, list_h_PU200_prompt_genMatched_EE_vtx;
  for (int ih=0; ih<track_pv_dz_cut_EB.size(); ih++) {
	// (muon, track)
	TH1D* h_PU200_prompt_EB = new TH1D(Form("h_PU200_prompt_EB_dz%d", ih), Form("h_PU200_prompt_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_prompt_EE = new TH1D(Form("h_PU200_prompt_EE_dz%d", ih), Form("h_PU200_prompt_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_prompt_1sigma_EB = new TH1D(Form("h_PU200_prompt_1sigma_EB_dz%d", ih), Form("h_PU200_prompt_1sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_prompt_1sigma_EE = new TH1D(Form("h_PU200_prompt_1sigma_EE_dz%d", ih), Form("h_PU200_prompt_1sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_prompt_2sigma_EB = new TH1D(Form("h_PU200_prompt_2sigma_EB_dz%d", ih), Form("h_PU200_prompt_2sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_prompt_2sigma_EE = new TH1D(Form("h_PU200_prompt_2sigma_EE_dz%d", ih), Form("h_PU200_prompt_2sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_prompt_3sigma_EB = new TH1D(Form("h_PU200_prompt_3sigma_EB_dz%d", ih), Form("h_PU200_prompt_3sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_prompt_3sigma_EE = new TH1D(Form("h_PU200_prompt_3sigma_EE_dz%d", ih), Form("h_PU200_prompt_3sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_prompt_4sigma_EB = new TH1D(Form("h_PU200_prompt_4sigma_EB_dz%d", ih), Form("h_PU200_prompt_4sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_prompt_4sigma_EE = new TH1D(Form("h_PU200_prompt_4sigma_EE_dz%d", ih), Form("h_PU200_prompt_4sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_prompt_genMatched_EB = new TH1D(Form("h_PU200_prompt_genMatched_EB_dz%d", ih), Form("h_PU200_prompt_genMatched_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_prompt_genMatched_EE = new TH1D(Form("h_PU200_prompt_genMatched_EE_dz%d", ih), Form("h_PU200_prompt_genMatched_EE_dz%d", ih), nbin, 0, 4);
	list_h_PU200_prompt_EB.emplace_back(h_PU200_prompt_EB); list_h_PU200_prompt_EE.emplace_back(h_PU200_prompt_EE);
	list_h_PU200_prompt_1sigma_EB.emplace_back(h_PU200_prompt_1sigma_EB); list_h_PU200_prompt_1sigma_EE.emplace_back(h_PU200_prompt_1sigma_EE);
	list_h_PU200_prompt_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB); list_h_PU200_prompt_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE);
	list_h_PU200_prompt_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB); list_h_PU200_prompt_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE);
	list_h_PU200_prompt_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB); list_h_PU200_prompt_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE);
	list_h_PU200_prompt_genMatched_EB.emplace_back(h_PU200_prompt_genMatched_EB); list_h_PU200_prompt_genMatched_EE.emplace_back(h_PU200_prompt_genMatched_EE);
	// (PV, track)
	TH1D* h_PU200_prompt_EB_vtx = new TH1D(Form("h_PU200_prompt_EB_vtx_dz%d", ih), Form("h_PU200_prompt_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_prompt_EE_vtx = new TH1D(Form("h_PU200_prompt_EE_vtx_dz%d", ih), Form("h_PU200_prompt_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_prompt_1sigma_EB_vtx = new TH1D(Form("h_PU200_prompt_1sigma_EB_vtx_dz%d", ih), Form("h_PU200_prompt_1sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_prompt_1sigma_EE_vtx = new TH1D(Form("h_PU200_prompt_1sigma_EE_vtx_dz%d", ih), Form("h_PU200_prompt_1sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_prompt_2sigma_EB_vtx = new TH1D(Form("h_PU200_prompt_2sigma_EB_vtx_dz%d", ih), Form("h_PU200_prompt_2sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_prompt_2sigma_EE_vtx = new TH1D(Form("h_PU200_prompt_2sigma_EE_vtx_dz%d", ih), Form("h_PU200_prompt_2sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_prompt_3sigma_EB_vtx = new TH1D(Form("h_PU200_prompt_3sigma_EB_vtx_dz%d", ih), Form("h_PU200_prompt_3sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_prompt_3sigma_EE_vtx = new TH1D(Form("h_PU200_prompt_3sigma_EE_vtx_dz%d", ih), Form("h_PU200_prompt_3sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_prompt_4sigma_EB_vtx = new TH1D(Form("h_PU200_prompt_4sigma_EB_vtx_dz%d", ih), Form("h_PU200_prompt_4sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_prompt_4sigma_EE_vtx = new TH1D(Form("h_PU200_prompt_4sigma_EE_vtx_dz%d", ih), Form("h_PU200_prompt_4sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_prompt_genMatched_EB_vtx = new TH1D(Form("h_PU200_prompt_genMatched_EB_vtx_dz%d", ih), Form("h_PU200_prompt_genMatched_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_prompt_genMatched_EE_vtx = new TH1D(Form("h_PU200_prompt_genMatched_EE_vtx_dz%d", ih), Form("h_PU200_prompt_genMatched_EE_vtx_dz%d", ih), nbin, 0, 4);
	list_h_PU200_prompt_EB_vtx.emplace_back(h_PU200_prompt_EB_vtx); list_h_PU200_prompt_EE_vtx.emplace_back(h_PU200_prompt_EE_vtx);
	list_h_PU200_prompt_1sigma_EB_vtx.emplace_back(h_PU200_prompt_1sigma_EB_vtx); list_h_PU200_prompt_1sigma_EE_vtx.emplace_back(h_PU200_prompt_1sigma_EE_vtx);
	list_h_PU200_prompt_2sigma_EB_vtx.emplace_back(h_PU200_prompt_2sigma_EB_vtx); list_h_PU200_prompt_2sigma_EE_vtx.emplace_back(h_PU200_prompt_2sigma_EE_vtx);
	list_h_PU200_prompt_3sigma_EB_vtx.emplace_back(h_PU200_prompt_3sigma_EB_vtx); list_h_PU200_prompt_3sigma_EE_vtx.emplace_back(h_PU200_prompt_3sigma_EE_vtx);
	list_h_PU200_prompt_4sigma_EB_vtx.emplace_back(h_PU200_prompt_4sigma_EB_vtx); list_h_PU200_prompt_4sigma_EE_vtx.emplace_back(h_PU200_prompt_4sigma_EE_vtx);
	list_h_PU200_prompt_genMatched_EB_vtx.emplace_back(h_PU200_prompt_genMatched_EB_vtx); list_h_PU200_prompt_genMatched_EE_vtx.emplace_back(h_PU200_prompt_genMatched_EE_vtx);
  }
  // (muon, track)
  vector<float> sumiso_PU200_prompt_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_prompt_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_prompt_1sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_prompt_1sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_prompt_2sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_prompt_2sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_prompt_3sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_prompt_3sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_prompt_4sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_prompt_4sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_prompt_genMatched_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_prompt_genMatched_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_prompt_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_prompt_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_prompt_1sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_prompt_1sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_prompt_2sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_prompt_2sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_prompt_3sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_prompt_3sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_prompt_4sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_prompt_4sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_prompt_genMatched_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_prompt_genMatched_EE(track_pv_dz_cut_EE.size(), 0);
  // (PV, track)
  vector<float> sumiso_PU200_prompt_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_prompt_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_prompt_1sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_prompt_1sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_prompt_2sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_prompt_2sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_prompt_3sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_prompt_3sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_prompt_4sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_prompt_4sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_prompt_genMatched_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_prompt_genMatched_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_prompt_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_prompt_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_prompt_1sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_prompt_1sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_prompt_2sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_prompt_2sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_prompt_3sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_prompt_3sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_prompt_4sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_prompt_4sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_prompt_genMatched_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_prompt_genMatched_EE_vtx(track_pv_dz_cut_EE.size(), 0);


    // nonprompt
	// (muon, track)
  vector<TH1D*> list_h_PU200_nonprompt_EB,            list_h_PU200_nonprompt_EE;
  vector<TH1D*> list_h_PU200_nonprompt_1sigma_EB,     list_h_PU200_nonprompt_1sigma_EE;
  vector<TH1D*> list_h_PU200_nonprompt_2sigma_EB,     list_h_PU200_nonprompt_2sigma_EE;
  vector<TH1D*> list_h_PU200_nonprompt_3sigma_EB,     list_h_PU200_nonprompt_3sigma_EE;
  vector<TH1D*> list_h_PU200_nonprompt_4sigma_EB,     list_h_PU200_nonprompt_4sigma_EE;
  vector<TH1D*> list_h_PU200_nonprompt_genMatched_EB, list_h_PU200_nonprompt_genMatched_EE;
	// (PV, track)
  vector<TH1D*> list_h_PU200_nonprompt_EB_vtx,            list_h_PU200_nonprompt_EE_vtx;
  vector<TH1D*> list_h_PU200_nonprompt_1sigma_EB_vtx,     list_h_PU200_nonprompt_1sigma_EE_vtx;
  vector<TH1D*> list_h_PU200_nonprompt_2sigma_EB_vtx,     list_h_PU200_nonprompt_2sigma_EE_vtx;
  vector<TH1D*> list_h_PU200_nonprompt_3sigma_EB_vtx,     list_h_PU200_nonprompt_3sigma_EE_vtx;
  vector<TH1D*> list_h_PU200_nonprompt_4sigma_EB_vtx,     list_h_PU200_nonprompt_4sigma_EE_vtx;
  vector<TH1D*> list_h_PU200_nonprompt_genMatched_EB_vtx, list_h_PU200_nonprompt_genMatched_EE_vtx;
  for (int ih=0; ih<track_pv_dz_cut_EB.size(); ih++) {
	// (muon,track)
	TH1D* h_PU200_nonprompt_EB = new TH1D(Form("h_PU200_nonprompt_EB_dz%d", ih), Form("h_PU200_nonprompt_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_nonprompt_EE = new TH1D(Form("h_PU200_nonprompt_EE_dz%d", ih), Form("h_PU200_nonprompt_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_nonprompt_1sigma_EB = new TH1D(Form("h_PU200_nonprompt_1sigma_EB_dz%d", ih), Form("h_PU200_nonprompt_1sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_nonprompt_1sigma_EE = new TH1D(Form("h_PU200_nonprompt_1sigma_EE_dz%d", ih), Form("h_PU200_nonprompt_1sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_nonprompt_2sigma_EB = new TH1D(Form("h_PU200_nonprompt_2sigma_EB_dz%d", ih), Form("h_PU200_nonprompt_2sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_nonprompt_2sigma_EE = new TH1D(Form("h_PU200_nonprompt_2sigma_EE_dz%d", ih), Form("h_PU200_nonprompt_2sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_nonprompt_3sigma_EB = new TH1D(Form("h_PU200_nonprompt_3sigma_EB_dz%d", ih), Form("h_PU200_nonprompt_3sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_nonprompt_3sigma_EE = new TH1D(Form("h_PU200_nonprompt_3sigma_EE_dz%d", ih), Form("h_PU200_nonprompt_3sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_nonprompt_4sigma_EB = new TH1D(Form("h_PU200_nonprompt_4sigma_EB_dz%d", ih), Form("h_PU200_nonprompt_4sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_nonprompt_4sigma_EE = new TH1D(Form("h_PU200_nonprompt_4sigma_EE_dz%d", ih), Form("h_PU200_nonprompt_4sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_nonprompt_genMatched_EB = new TH1D(Form("h_PU200_nonprompt_genMatched_EB_dz%d", ih), Form("h_PU200_nonprompt_genMatched_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_nonprompt_genMatched_EE = new TH1D(Form("h_PU200_nonprompt_genMatched_EE_dz%d", ih), Form("h_PU200_nonprompt_genMatched_EE_dz%d", ih), nbin, 0, 4);
	list_h_PU200_nonprompt_EB.emplace_back(h_PU200_nonprompt_EB); list_h_PU200_nonprompt_EE.emplace_back(h_PU200_nonprompt_EE);
	list_h_PU200_nonprompt_1sigma_EB.emplace_back(h_PU200_nonprompt_1sigma_EB); list_h_PU200_nonprompt_1sigma_EE.emplace_back(h_PU200_nonprompt_1sigma_EE);
	list_h_PU200_nonprompt_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB); list_h_PU200_nonprompt_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE);
	list_h_PU200_nonprompt_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB); list_h_PU200_nonprompt_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE);
	list_h_PU200_nonprompt_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB); list_h_PU200_nonprompt_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE);
	list_h_PU200_nonprompt_genMatched_EB.emplace_back(h_PU200_nonprompt_genMatched_EB); list_h_PU200_nonprompt_genMatched_EE.emplace_back(h_PU200_nonprompt_genMatched_EE);
	// (PV,track)
	TH1D* h_PU200_nonprompt_EB_vtx = new TH1D(Form("h_PU200_nonprompt_EB_vtx_dz%d", ih), Form("h_PU200_nonprompt_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_nonprompt_EE_vtx = new TH1D(Form("h_PU200_nonprompt_EE_vtx_dz%d", ih), Form("h_PU200_nonprompt_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_nonprompt_1sigma_EB_vtx = new TH1D(Form("h_PU200_nonprompt_1sigma_EB_vtx_dz%d", ih), Form("h_PU200_nonprompt_1sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_nonprompt_1sigma_EE_vtx = new TH1D(Form("h_PU200_nonprompt_1sigma_EE_vtx_dz%d", ih), Form("h_PU200_nonprompt_1sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_nonprompt_2sigma_EB_vtx = new TH1D(Form("h_PU200_nonprompt_2sigma_EB_vtx_dz%d", ih), Form("h_PU200_nonprompt_2sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_nonprompt_2sigma_EE_vtx = new TH1D(Form("h_PU200_nonprompt_2sigma_EE_vtx_dz%d", ih), Form("h_PU200_nonprompt_2sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_nonprompt_3sigma_EB_vtx = new TH1D(Form("h_PU200_nonprompt_3sigma_EB_vtx_dz%d", ih), Form("h_PU200_nonprompt_3sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_nonprompt_3sigma_EE_vtx = new TH1D(Form("h_PU200_nonprompt_3sigma_EE_vtx_dz%d", ih), Form("h_PU200_nonprompt_3sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_nonprompt_4sigma_EB_vtx = new TH1D(Form("h_PU200_nonprompt_4sigma_EB_vtx_dz%d", ih), Form("h_PU200_nonprompt_4sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_nonprompt_4sigma_EE_vtx = new TH1D(Form("h_PU200_nonprompt_4sigma_EE_vtx_dz%d", ih), Form("h_PU200_nonprompt_4sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_PU200_nonprompt_genMatched_EB_vtx = new TH1D(Form("h_PU200_nonprompt_genMatched_EB_vtx_dz%d", ih), Form("h_PU200_nonprompt_genMatched_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_PU200_nonprompt_genMatched_EE_vtx = new TH1D(Form("h_PU200_nonprompt_genMatched_EE_vtx_dz%d", ih), Form("h_PU200_nonprompt_genMatched_EE_vtx_dz%d", ih), nbin, 0, 4);
	list_h_PU200_nonprompt_EB_vtx.emplace_back(h_PU200_nonprompt_EB_vtx); list_h_PU200_nonprompt_EE_vtx.emplace_back(h_PU200_nonprompt_EE_vtx);
	list_h_PU200_nonprompt_1sigma_EB_vtx.emplace_back(h_PU200_nonprompt_1sigma_EB_vtx); list_h_PU200_nonprompt_1sigma_EE_vtx.emplace_back(h_PU200_nonprompt_1sigma_EE_vtx);
	list_h_PU200_nonprompt_2sigma_EB_vtx.emplace_back(h_PU200_nonprompt_2sigma_EB_vtx); list_h_PU200_nonprompt_2sigma_EE_vtx.emplace_back(h_PU200_nonprompt_2sigma_EE_vtx);
	list_h_PU200_nonprompt_3sigma_EB_vtx.emplace_back(h_PU200_nonprompt_3sigma_EB_vtx); list_h_PU200_nonprompt_3sigma_EE_vtx.emplace_back(h_PU200_nonprompt_3sigma_EE_vtx);
	list_h_PU200_nonprompt_4sigma_EB_vtx.emplace_back(h_PU200_nonprompt_4sigma_EB_vtx); list_h_PU200_nonprompt_4sigma_EE_vtx.emplace_back(h_PU200_nonprompt_4sigma_EE_vtx);
	list_h_PU200_nonprompt_genMatched_EB_vtx.emplace_back(h_PU200_nonprompt_genMatched_EB_vtx); list_h_PU200_nonprompt_genMatched_EE_vtx.emplace_back(h_PU200_nonprompt_genMatched_EE_vtx);

  }
  // (muon, track)
  vector<float> sumiso_PU200_nonprompt_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_nonprompt_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_nonprompt_1sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_nonprompt_1sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_nonprompt_2sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_nonprompt_2sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_nonprompt_3sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_nonprompt_3sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_nonprompt_4sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_nonprompt_4sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_nonprompt_genMatched_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_nonprompt_genMatched_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_nonprompt_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_nonprompt_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_nonprompt_1sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_nonprompt_1sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_nonprompt_2sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_nonprompt_2sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_nonprompt_3sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_nonprompt_3sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_nonprompt_4sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_nonprompt_4sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_nonprompt_genMatched_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_nonprompt_genMatched_EE(track_pv_dz_cut_EE.size(), 0);
  // (PV, track)
  vector<float> sumiso_PU200_nonprompt_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_nonprompt_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_nonprompt_1sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_nonprompt_1sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_nonprompt_2sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_nonprompt_2sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_nonprompt_3sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_nonprompt_3sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_nonprompt_4sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_nonprompt_4sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_PU200_nonprompt_genMatched_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_PU200_nonprompt_genMatched_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_nonprompt_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_nonprompt_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_nonprompt_1sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_nonprompt_1sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_nonprompt_2sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_nonprompt_2sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_nonprompt_3sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_nonprompt_3sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_nonprompt_4sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_nonprompt_4sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_PU200_nonprompt_genMatched_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_PU200_nonprompt_genMatched_EE_vtx(track_pv_dz_cut_EE.size(), 0);


  // noPU
    // prompt
	// (muon, track)
  vector<TH1D*> list_h_noPU_prompt_EB, list_h_noPU_prompt_EE;
  vector<TH1D*> list_h_noPU_prompt_1sigma_EB, list_h_noPU_prompt_1sigma_EE;
  vector<TH1D*> list_h_noPU_prompt_2sigma_EB, list_h_noPU_prompt_2sigma_EE;
  vector<TH1D*> list_h_noPU_prompt_3sigma_EB, list_h_noPU_prompt_3sigma_EE;
  vector<TH1D*> list_h_noPU_prompt_4sigma_EB, list_h_noPU_prompt_4sigma_EE;
  vector<TH1D*> list_h_noPU_prompt_genMatched_EB, list_h_noPU_prompt_genMatched_EE;
	// (PV, track)
  vector<TH1D*> list_h_noPU_prompt_EB_vtx, list_h_noPU_prompt_EE_vtx;
  vector<TH1D*> list_h_noPU_prompt_1sigma_EB_vtx, list_h_noPU_prompt_1sigma_EE_vtx;
  vector<TH1D*> list_h_noPU_prompt_2sigma_EB_vtx, list_h_noPU_prompt_2sigma_EE_vtx;
  vector<TH1D*> list_h_noPU_prompt_3sigma_EB_vtx, list_h_noPU_prompt_3sigma_EE_vtx;
  vector<TH1D*> list_h_noPU_prompt_4sigma_EB_vtx, list_h_noPU_prompt_4sigma_EE_vtx;
  vector<TH1D*> list_h_noPU_prompt_genMatched_EB_vtx, list_h_noPU_prompt_genMatched_EE_vtx;
  for (int ih=0; ih<track_pv_dz_cut_EB.size(); ih++) {
	// (muon, track)
	TH1D* h_noPU_prompt_EB = new TH1D(Form("h_noPU_prompt_EB_dz%d", ih), Form("h_noPU_prompt_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_prompt_EE = new TH1D(Form("h_noPU_prompt_EE_dz%d", ih), Form("h_noPU_prompt_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_prompt_1sigma_EB = new TH1D(Form("h_noPU_prompt_1sigma_EB_dz%d", ih), Form("h_noPU_prompt_1sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_prompt_1sigma_EE = new TH1D(Form("h_noPU_prompt_1sigma_EE_dz%d", ih), Form("h_noPU_prompt_1sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_prompt_2sigma_EB = new TH1D(Form("h_noPU_prompt_2sigma_EB_dz%d", ih), Form("h_noPU_prompt_2sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_prompt_2sigma_EE = new TH1D(Form("h_noPU_prompt_2sigma_EE_dz%d", ih), Form("h_noPU_prompt_2sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_prompt_3sigma_EB = new TH1D(Form("h_noPU_prompt_3sigma_EB_dz%d", ih), Form("h_noPU_prompt_3sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_prompt_3sigma_EE = new TH1D(Form("h_noPU_prompt_3sigma_EE_dz%d", ih), Form("h_noPU_prompt_3sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_prompt_4sigma_EB = new TH1D(Form("h_noPU_prompt_4sigma_EB_dz%d", ih), Form("h_noPU_prompt_4sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_prompt_4sigma_EE = new TH1D(Form("h_noPU_prompt_4sigma_EE_dz%d", ih), Form("h_noPU_prompt_4sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_prompt_genMatched_EB = new TH1D(Form("h_noPU_prompt_genMatched_EB_dz%d", ih), Form("h_noPU_prompt_genMatched_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_prompt_genMatched_EE = new TH1D(Form("h_noPU_prompt_genMatched_EE_dz%d", ih), Form("h_noPU_prompt_genMatched_EE_dz%d", ih), nbin, 0, 4);
	list_h_noPU_prompt_EB.emplace_back(h_noPU_prompt_EB); list_h_noPU_prompt_EE.emplace_back(h_noPU_prompt_EE);
	list_h_noPU_prompt_1sigma_EB.emplace_back(h_noPU_prompt_1sigma_EB); list_h_noPU_prompt_1sigma_EE.emplace_back(h_noPU_prompt_1sigma_EE);
	list_h_noPU_prompt_2sigma_EB.emplace_back(h_noPU_prompt_2sigma_EB); list_h_noPU_prompt_2sigma_EE.emplace_back(h_noPU_prompt_2sigma_EE);
	list_h_noPU_prompt_3sigma_EB.emplace_back(h_noPU_prompt_3sigma_EB); list_h_noPU_prompt_3sigma_EE.emplace_back(h_noPU_prompt_3sigma_EE);
	list_h_noPU_prompt_4sigma_EB.emplace_back(h_noPU_prompt_4sigma_EB); list_h_noPU_prompt_4sigma_EE.emplace_back(h_noPU_prompt_4sigma_EE);
	list_h_noPU_prompt_genMatched_EB.emplace_back(h_noPU_prompt_genMatched_EB); list_h_noPU_prompt_genMatched_EE.emplace_back(h_noPU_prompt_genMatched_EE);
	// (PV, track)
	TH1D* h_noPU_prompt_EB_vtx = new TH1D(Form("h_noPU_prompt_EB_vtx_dz%d", ih), Form("h_noPU_prompt_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_prompt_EE_vtx = new TH1D(Form("h_noPU_prompt_EE_vtx_dz%d", ih), Form("h_noPU_prompt_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_prompt_1sigma_EB_vtx = new TH1D(Form("h_noPU_prompt_1sigma_EB_vtx_dz%d", ih), Form("h_noPU_prompt_1sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_prompt_1sigma_EE_vtx = new TH1D(Form("h_noPU_prompt_1sigma_EE_vtx_dz%d", ih), Form("h_noPU_prompt_1sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_prompt_2sigma_EB_vtx = new TH1D(Form("h_noPU_prompt_2sigma_EB_vtx_dz%d", ih), Form("h_noPU_prompt_2sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_prompt_2sigma_EE_vtx = new TH1D(Form("h_noPU_prompt_2sigma_EE_vtx_dz%d", ih), Form("h_noPU_prompt_2sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_prompt_3sigma_EB_vtx = new TH1D(Form("h_noPU_prompt_3sigma_EB_vtx_dz%d", ih), Form("h_noPU_prompt_3sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_prompt_3sigma_EE_vtx = new TH1D(Form("h_noPU_prompt_3sigma_EE_vtx_dz%d", ih), Form("h_noPU_prompt_3sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_prompt_4sigma_EB_vtx = new TH1D(Form("h_noPU_prompt_4sigma_EB_vtx_dz%d", ih), Form("h_noPU_prompt_4sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_prompt_4sigma_EE_vtx = new TH1D(Form("h_noPU_prompt_4sigma_EE_vtx_dz%d", ih), Form("h_noPU_prompt_4sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_prompt_genMatched_EB_vtx = new TH1D(Form("h_noPU_prompt_genMatched_EB_vtx_dz%d", ih), Form("h_noPU_prompt_genMatched_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_prompt_genMatched_EE_vtx = new TH1D(Form("h_noPU_prompt_genMatched_EE_vtx_dz%d", ih), Form("h_noPU_prompt_genMatched_EE_vtx_dz%d", ih), nbin, 0, 4);
	list_h_noPU_prompt_EB_vtx.emplace_back(h_noPU_prompt_EB_vtx); list_h_noPU_prompt_EE_vtx.emplace_back(h_noPU_prompt_EE_vtx);
	list_h_noPU_prompt_1sigma_EB_vtx.emplace_back(h_noPU_prompt_1sigma_EB_vtx); list_h_noPU_prompt_1sigma_EE_vtx.emplace_back(h_noPU_prompt_1sigma_EE_vtx);
	list_h_noPU_prompt_2sigma_EB_vtx.emplace_back(h_noPU_prompt_2sigma_EB_vtx); list_h_noPU_prompt_2sigma_EE_vtx.emplace_back(h_noPU_prompt_2sigma_EE_vtx);
	list_h_noPU_prompt_3sigma_EB_vtx.emplace_back(h_noPU_prompt_3sigma_EB_vtx); list_h_noPU_prompt_3sigma_EE_vtx.emplace_back(h_noPU_prompt_3sigma_EE_vtx);
	list_h_noPU_prompt_4sigma_EB_vtx.emplace_back(h_noPU_prompt_4sigma_EB_vtx); list_h_noPU_prompt_4sigma_EE_vtx.emplace_back(h_noPU_prompt_4sigma_EE_vtx);
	list_h_noPU_prompt_genMatched_EB_vtx.emplace_back(h_noPU_prompt_genMatched_EB_vtx); list_h_noPU_prompt_genMatched_EE_vtx.emplace_back(h_noPU_prompt_genMatched_EE_vtx);

  }
  // (muon, track)
  vector<float> sumiso_noPU_prompt_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_prompt_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_prompt_1sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_prompt_1sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_prompt_2sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_prompt_2sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_prompt_3sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_prompt_3sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_prompt_4sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_prompt_4sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_prompt_genMatched_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_prompt_genMatched_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_prompt_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_prompt_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_prompt_1sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_prompt_1sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_prompt_2sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_prompt_2sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_prompt_3sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_prompt_3sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_prompt_4sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_prompt_4sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_prompt_genMatched_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_prompt_genMatched_EE(track_pv_dz_cut_EE.size(), 0);
  // (PV, track)
  vector<float> sumiso_noPU_prompt_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_prompt_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_prompt_1sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_prompt_1sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_prompt_2sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_prompt_2sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_prompt_3sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_prompt_3sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_prompt_4sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_prompt_4sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_prompt_genMatched_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_prompt_genMatched_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_prompt_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_prompt_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_prompt_1sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_prompt_1sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_prompt_2sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_prompt_2sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_prompt_3sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_prompt_3sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_prompt_4sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_prompt_4sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_prompt_genMatched_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_prompt_genMatched_EE_vtx(track_pv_dz_cut_EE.size(), 0);


    // nonprompt
	// (muon, track)
  vector<TH1D*> list_h_noPU_nonprompt_EB, list_h_noPU_nonprompt_EE;
  vector<TH1D*> list_h_noPU_nonprompt_1sigma_EB, list_h_noPU_nonprompt_1sigma_EE;
  vector<TH1D*> list_h_noPU_nonprompt_2sigma_EB, list_h_noPU_nonprompt_2sigma_EE;
  vector<TH1D*> list_h_noPU_nonprompt_3sigma_EB, list_h_noPU_nonprompt_3sigma_EE;
  vector<TH1D*> list_h_noPU_nonprompt_4sigma_EB, list_h_noPU_nonprompt_4sigma_EE;
  vector<TH1D*> list_h_noPU_nonprompt_genMatched_EB, list_h_noPU_nonprompt_genMatched_EE;
	// (PV, track)
  vector<TH1D*> list_h_noPU_nonprompt_EB_vtx, list_h_noPU_nonprompt_EE_vtx;
  vector<TH1D*> list_h_noPU_nonprompt_1sigma_EB_vtx, list_h_noPU_nonprompt_1sigma_EE_vtx;
  vector<TH1D*> list_h_noPU_nonprompt_2sigma_EB_vtx, list_h_noPU_nonprompt_2sigma_EE_vtx;
  vector<TH1D*> list_h_noPU_nonprompt_3sigma_EB_vtx, list_h_noPU_nonprompt_3sigma_EE_vtx;
  vector<TH1D*> list_h_noPU_nonprompt_4sigma_EB_vtx, list_h_noPU_nonprompt_4sigma_EE_vtx;
  vector<TH1D*> list_h_noPU_nonprompt_genMatched_EB_vtx, list_h_noPU_nonprompt_genMatched_EE_vtx;

  for (int ih=0; ih<track_pv_dz_cut_EB.size(); ih++) {
	// (muon, track)
	TH1D* h_noPU_nonprompt_EB = new TH1D(Form("h_noPU_nonprompt_EB_dz%d", ih), Form("h_noPU_nonprompt_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_nonprompt_EE = new TH1D(Form("h_noPU_nonprompt_EE_dz%d", ih), Form("h_noPU_nonprompt_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_nonprompt_1sigma_EB = new TH1D(Form("h_noPU_nonprompt_1sigma_EB_dz%d", ih), Form("h_noPU_nonprompt_1sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_nonprompt_1sigma_EE = new TH1D(Form("h_noPU_nonprompt_1sigma_EE_dz%d", ih), Form("h_noPU_nonprompt_1sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_nonprompt_2sigma_EB = new TH1D(Form("h_noPU_nonprompt_2sigma_EB_dz%d", ih), Form("h_noPU_nonprompt_2sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_nonprompt_2sigma_EE = new TH1D(Form("h_noPU_nonprompt_2sigma_EE_dz%d", ih), Form("h_noPU_nonprompt_2sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_nonprompt_3sigma_EB = new TH1D(Form("h_noPU_nonprompt_3sigma_EB_dz%d", ih), Form("h_noPU_nonprompt_3sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_nonprompt_3sigma_EE = new TH1D(Form("h_noPU_nonprompt_3sigma_EE_dz%d", ih), Form("h_noPU_nonprompt_3sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_nonprompt_4sigma_EB = new TH1D(Form("h_noPU_nonprompt_4sigma_EB_dz%d", ih), Form("h_noPU_nonprompt_4sigma_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_nonprompt_4sigma_EE = new TH1D(Form("h_noPU_nonprompt_4sigma_EE_dz%d", ih), Form("h_noPU_nonprompt_4sigma_EE_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_nonprompt_genMatched_EB = new TH1D(Form("h_noPU_nonprompt_genMatched_EB_dz%d", ih), Form("h_noPU_nonprompt_genMatched_EB_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_nonprompt_genMatched_EE = new TH1D(Form("h_noPU_nonprompt_genMatched_EE_dz%d", ih), Form("h_noPU_nonprompt_genMatched_EE_dz%d", ih), nbin, 0, 4);
	list_h_noPU_nonprompt_EB.emplace_back(h_noPU_nonprompt_EB); list_h_noPU_nonprompt_EE.emplace_back(h_noPU_nonprompt_EE);
	list_h_noPU_nonprompt_1sigma_EB.emplace_back(h_noPU_nonprompt_1sigma_EB); list_h_noPU_nonprompt_1sigma_EE.emplace_back(h_noPU_nonprompt_1sigma_EE);
	list_h_noPU_nonprompt_2sigma_EB.emplace_back(h_noPU_nonprompt_2sigma_EB); list_h_noPU_nonprompt_2sigma_EE.emplace_back(h_noPU_nonprompt_2sigma_EE);
	list_h_noPU_nonprompt_3sigma_EB.emplace_back(h_noPU_nonprompt_3sigma_EB); list_h_noPU_nonprompt_3sigma_EE.emplace_back(h_noPU_nonprompt_3sigma_EE);
	list_h_noPU_nonprompt_4sigma_EB.emplace_back(h_noPU_nonprompt_4sigma_EB); list_h_noPU_nonprompt_4sigma_EE.emplace_back(h_noPU_nonprompt_4sigma_EE);
	list_h_noPU_nonprompt_genMatched_EB.emplace_back(h_noPU_nonprompt_genMatched_EB); list_h_noPU_nonprompt_genMatched_EE.emplace_back(h_noPU_nonprompt_genMatched_EE);
	// (PV, track)
	TH1D* h_noPU_nonprompt_EB_vtx = new TH1D(Form("h_noPU_nonprompt_EB_vtx_dz%d", ih), Form("h_noPU_nonprompt_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_nonprompt_EE_vtx = new TH1D(Form("h_noPU_nonprompt_EE_vtx_dz%d", ih), Form("h_noPU_nonprompt_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_nonprompt_1sigma_EB_vtx = new TH1D(Form("h_noPU_nonprompt_1sigma_EB_vtx_dz%d", ih), Form("h_noPU_nonprompt_1sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_nonprompt_1sigma_EE_vtx = new TH1D(Form("h_noPU_nonprompt_1sigma_EE_vtx_dz%d", ih), Form("h_noPU_nonprompt_1sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_nonprompt_2sigma_EB_vtx = new TH1D(Form("h_noPU_nonprompt_2sigma_EB_vtx_dz%d", ih), Form("h_noPU_nonprompt_2sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_nonprompt_2sigma_EE_vtx = new TH1D(Form("h_noPU_nonprompt_2sigma_EE_vtx_dz%d", ih), Form("h_noPU_nonprompt_2sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_nonprompt_3sigma_EB_vtx = new TH1D(Form("h_noPU_nonprompt_3sigma_EB_vtx_dz%d", ih), Form("h_noPU_nonprompt_3sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_nonprompt_3sigma_EE_vtx = new TH1D(Form("h_noPU_nonprompt_3sigma_EE_vtx_dz%d", ih), Form("h_noPU_nonprompt_3sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_nonprompt_4sigma_EB_vtx = new TH1D(Form("h_noPU_nonprompt_4sigma_EB_vtx_dz%d", ih), Form("h_noPU_nonprompt_4sigma_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_nonprompt_4sigma_EE_vtx = new TH1D(Form("h_noPU_nonprompt_4sigma_EE_vtx_dz%d", ih), Form("h_noPU_nonprompt_4sigma_EE_vtx_dz%d", ih), nbin, 0, 4);
	TH1D* h_noPU_nonprompt_genMatched_EB_vtx = new TH1D(Form("h_noPU_nonprompt_genMatched_EB_vtx_dz%d", ih), Form("h_noPU_nonprompt_genMatched_EB_vtx_dz%d", ih), nbin, 0, 4);
    TH1D* h_noPU_nonprompt_genMatched_EE_vtx = new TH1D(Form("h_noPU_nonprompt_genMatched_EE_vtx_dz%d", ih), Form("h_noPU_nonprompt_genMatched_EE_vtx_dz%d", ih), nbin, 0, 4);
	list_h_noPU_nonprompt_EB_vtx.emplace_back(h_noPU_nonprompt_EB_vtx); list_h_noPU_nonprompt_EE_vtx.emplace_back(h_noPU_nonprompt_EE_vtx);
	list_h_noPU_nonprompt_1sigma_EB_vtx.emplace_back(h_noPU_nonprompt_1sigma_EB_vtx); list_h_noPU_nonprompt_1sigma_EE_vtx.emplace_back(h_noPU_nonprompt_1sigma_EE_vtx);
	list_h_noPU_nonprompt_2sigma_EB_vtx.emplace_back(h_noPU_nonprompt_2sigma_EB_vtx); list_h_noPU_nonprompt_2sigma_EE_vtx.emplace_back(h_noPU_nonprompt_2sigma_EE_vtx);
	list_h_noPU_nonprompt_3sigma_EB_vtx.emplace_back(h_noPU_nonprompt_3sigma_EB_vtx); list_h_noPU_nonprompt_3sigma_EE_vtx.emplace_back(h_noPU_nonprompt_3sigma_EE_vtx);
	list_h_noPU_nonprompt_4sigma_EB_vtx.emplace_back(h_noPU_nonprompt_4sigma_EB_vtx); list_h_noPU_nonprompt_4sigma_EE_vtx.emplace_back(h_noPU_nonprompt_4sigma_EE_vtx);
	list_h_noPU_nonprompt_genMatched_EB_vtx.emplace_back(h_noPU_nonprompt_genMatched_EB_vtx); list_h_noPU_nonprompt_genMatched_EE_vtx.emplace_back(h_noPU_nonprompt_genMatched_EE_vtx);

  }
  // (muon, track)
  vector<float> sumiso_noPU_nonprompt_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_nonprompt_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_nonprompt_1sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_nonprompt_1sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_nonprompt_2sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_nonprompt_2sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_nonprompt_3sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_nonprompt_3sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_nonprompt_4sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_nonprompt_4sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_nonprompt_genMatched_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_nonprompt_genMatched_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_nonprompt_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_nonprompt_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_nonprompt_1sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_nonprompt_1sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_nonprompt_2sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_nonprompt_2sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_nonprompt_3sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_nonprompt_3sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_nonprompt_4sigma_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_nonprompt_4sigma_EE(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_nonprompt_genMatched_EB(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_nonprompt_genMatched_EE(track_pv_dz_cut_EE.size(), 0);
  // (PV, track)
  vector<float> sumiso_noPU_nonprompt_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_nonprompt_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_nonprompt_1sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_nonprompt_1sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_nonprompt_2sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_nonprompt_2sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_nonprompt_3sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_nonprompt_3sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_nonprompt_4sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_nonprompt_4sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> sumiso_noPU_nonprompt_genMatched_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> sumiso_noPU_nonprompt_genMatched_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_nonprompt_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_nonprompt_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_nonprompt_1sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_nonprompt_1sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_nonprompt_2sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_nonprompt_2sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_nonprompt_3sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_nonprompt_3sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_nonprompt_4sigma_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_nonprompt_4sigma_EE_vtx(track_pv_dz_cut_EE.size(), 0);
  vector<float> reliso_noPU_nonprompt_genMatched_EB_vtx(track_pv_dz_cut_EB.size(), 0);
  vector<float> reliso_noPU_nonprompt_genMatched_EE_vtx(track_pv_dz_cut_EE.size(), 0);



  //////////////////////
  //// PU200 prompt ////
  //////////////////////

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

        if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
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
		  // dz loop
		  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EB.at(idz)) continue;
  		    if(flag_track_PVweight) {
		  	  if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		    }
		    // noMTD case
			sumiso_PU200_prompt_EB[idz] += track_pt_->at(im).at(it);

			// MTD case
			// (muon, track)
			if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
			}
			if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_PU200_prompt_1sigma_EB[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_PU200_prompt_2sigma_EB[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_PU200_prompt_3sigma_EB[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_PU200_prompt_4sigma_EB[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_PU200_prompt_1sigma_EB[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_prompt_2sigma_EB[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_prompt_3sigma_EB[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_prompt_4sigma_EB[idz]+= track_pt_->at(im).at(it);
            }
            // (PV, track)
			if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
			}
			if(dtsig_cut_vtx<1.0 && dtsig_cut_vtx>0) sumiso_PU200_prompt_1sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<2.0 && dtsig_cut_vtx>0) sumiso_PU200_prompt_2sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<3.0 && dtsig_cut_vtx>0) sumiso_PU200_prompt_3sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<4.0 && dtsig_cut_vtx>0) sumiso_PU200_prompt_4sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut_vtx==0 && (vtx_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_PU200_prompt_1sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_prompt_2sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_prompt_3sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_prompt_4sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
            }


			// GEN case
			if(track_genMatched_->at(im).at(it)==1) {
			  sumiso_PU200_prompt_genMatched_EB[idz] += track_pt_->at(im).at(it);
			  if(track_bx_->at(im).at(it)!=0 || track_evtId_->at(im).at(it)!=0) num_gen_notPV_PU200_prompt_EB[idz]++;
			  num_gen_PU200_prompt_EB[idz]++;
			}

            dtsig_cut=0, dtsig_cut_vtx=0;
            dt_cut=0, dt_cut_vtx=0;
		  }
		}
		for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
		  // (muon, track)
		  reliso_PU200_prompt_EB[idz]        = sumiso_PU200_prompt_EB.at(idz)/muon_pt_->at(im);
		  reliso_PU200_prompt_1sigma_EB[idz] = sumiso_PU200_prompt_1sigma_EB.at(idz)/muon_pt_->at(im);
		  reliso_PU200_prompt_2sigma_EB[idz] = sumiso_PU200_prompt_2sigma_EB.at(idz)/muon_pt_->at(im);
		  reliso_PU200_prompt_3sigma_EB[idz] = sumiso_PU200_prompt_3sigma_EB.at(idz)/muon_pt_->at(im);
		  reliso_PU200_prompt_4sigma_EB[idz] = sumiso_PU200_prompt_4sigma_EB.at(idz)/muon_pt_->at(im);
		  reliso_PU200_prompt_genMatched_EB[idz] = sumiso_PU200_prompt_genMatched_EB.at(idz)/muon_pt_->at(im);
          // (PV, track)
		  reliso_PU200_prompt_EB_vtx[idz]        = sumiso_PU200_prompt_EB_vtx.at(idz)/muon_pt_->at(im);
		  reliso_PU200_prompt_1sigma_EB_vtx[idz] = sumiso_PU200_prompt_1sigma_EB_vtx.at(idz)/muon_pt_->at(im);
		  reliso_PU200_prompt_2sigma_EB_vtx[idz] = sumiso_PU200_prompt_2sigma_EB_vtx.at(idz)/muon_pt_->at(im);
		  reliso_PU200_prompt_3sigma_EB_vtx[idz] = sumiso_PU200_prompt_3sigma_EB_vtx.at(idz)/muon_pt_->at(im);
		  reliso_PU200_prompt_4sigma_EB_vtx[idz] = sumiso_PU200_prompt_4sigma_EB_vtx.at(idz)/muon_pt_->at(im);
		  reliso_PU200_prompt_genMatched_EB_vtx[idz] = sumiso_PU200_prompt_genMatched_EB_vtx.at(idz)/muon_pt_->at(im);

		  // Store value of relative isolation
		  // (muon, track)
		  list_h_PU200_prompt_EB.at(idz)->Fill(reliso_PU200_prompt_EB.at(idz));
		  list_h_PU200_prompt_1sigma_EB.at(idz)->Fill(reliso_PU200_prompt_1sigma_EB.at(idz));
		  list_h_PU200_prompt_2sigma_EB.at(idz)->Fill(reliso_PU200_prompt_2sigma_EB.at(idz));
		  list_h_PU200_prompt_3sigma_EB.at(idz)->Fill(reliso_PU200_prompt_3sigma_EB.at(idz));
		  list_h_PU200_prompt_4sigma_EB.at(idz)->Fill(reliso_PU200_prompt_4sigma_EB.at(idz));
		  list_h_PU200_prompt_genMatched_EB.at(idz)->Fill(reliso_PU200_prompt_genMatched_EB.at(idz));
          // (PV, track)
		  list_h_PU200_prompt_EB_vtx.at(idz)->Fill(reliso_PU200_prompt_EB_vtx.at(idz));
		  list_h_PU200_prompt_1sigma_EB_vtx.at(idz)->Fill(reliso_PU200_prompt_1sigma_EB_vtx.at(idz));
		  list_h_PU200_prompt_2sigma_EB_vtx.at(idz)->Fill(reliso_PU200_prompt_2sigma_EB_vtx.at(idz));
		  list_h_PU200_prompt_3sigma_EB_vtx.at(idz)->Fill(reliso_PU200_prompt_3sigma_EB_vtx.at(idz));
		  list_h_PU200_prompt_4sigma_EB_vtx.at(idz)->Fill(reliso_PU200_prompt_4sigma_EB_vtx.at(idz));
		  list_h_PU200_prompt_genMatched_EB_vtx.at(idz)->Fill(reliso_PU200_prompt_genMatched_EB_vtx.at(idz));

		}

		// Initialize
		for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
		  // (muon, track)
          sumiso_PU200_prompt_EB[idz]=0;
          sumiso_PU200_prompt_1sigma_EB[idz]=0;
          sumiso_PU200_prompt_2sigma_EB[idz]=0;
          sumiso_PU200_prompt_3sigma_EB[idz]=0;
          sumiso_PU200_prompt_4sigma_EB[idz]=0;
          sumiso_PU200_prompt_genMatched_EB[idz]=0;
          reliso_PU200_prompt_EB[idz]=0;
          reliso_PU200_prompt_1sigma_EB[idz]=0;
          reliso_PU200_prompt_2sigma_EB[idz]=0;
          reliso_PU200_prompt_3sigma_EB[idz]=0;
          reliso_PU200_prompt_4sigma_EB[idz]=0;
          reliso_PU200_prompt_genMatched_EB[idz]=0;
          // (PV, track)
          sumiso_PU200_prompt_EB_vtx[idz]=0;
          sumiso_PU200_prompt_1sigma_EB_vtx[idz]=0;
          sumiso_PU200_prompt_2sigma_EB_vtx[idz]=0;
          sumiso_PU200_prompt_3sigma_EB_vtx[idz]=0;
          sumiso_PU200_prompt_4sigma_EB_vtx[idz]=0;
          sumiso_PU200_prompt_genMatched_EB_vtx[idz]=0;
          reliso_PU200_prompt_EB_vtx[idz]=0;
          reliso_PU200_prompt_1sigma_EB_vtx[idz]=0;
          reliso_PU200_prompt_2sigma_EB_vtx[idz]=0;
          reliso_PU200_prompt_3sigma_EB_vtx[idz]=0;
          reliso_PU200_prompt_4sigma_EB_vtx[idz]=0;
          reliso_PU200_prompt_genMatched_EB_vtx[idz]=0;

		}
	  }

	  //////////////////////////
	  // prompt muon - Endcap //
	  //////////////////////////
	  else if(muon_prompt_->at(im)==1 && muon_isBarrel_->at(im)==0) {

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
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
		  // dz loop
		  for(int idz=0; idz<track_pv_dz_cut_EE.size(); idz++) {
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EE.at(idz)) continue;
		    if(flag_track_PVweight) {
			  if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		    }
		    // noMTD case
		    sumiso_PU200_prompt_EE[idz] += track_pt_->at(im).at(it);

			// MTD case
			// (muon, track)
            if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
            }
            if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_PU200_prompt_1sigma_EE[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_PU200_prompt_2sigma_EE[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_PU200_prompt_3sigma_EE[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_PU200_prompt_4sigma_EE[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_PU200_prompt_1sigma_EE[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_prompt_2sigma_EE[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_prompt_3sigma_EE[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_prompt_4sigma_EE[idz]+= track_pt_->at(im).at(it);
            }
            // (PV, track)
			if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
			}
			if(dtsig_cut_vtx<1.0 && dtsig_cut_vtx>0) sumiso_PU200_prompt_1sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<2.0 && dtsig_cut_vtx>0) sumiso_PU200_prompt_2sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<3.0 && dtsig_cut_vtx>0) sumiso_PU200_prompt_3sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<4.0 && dtsig_cut_vtx>0) sumiso_PU200_prompt_4sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut_vtx==0 && (vtx_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_PU200_prompt_1sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_prompt_2sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_prompt_3sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_prompt_4sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
            }


			// GEN case
			if(track_genMatched_->at(im).at(it)==1) {
			  sumiso_PU200_prompt_genMatched_EE[idz] += track_pt_->at(im).at(it);
			  if(track_bx_->at(im).at(it)!=0 || track_evtId_->at(im).at(it)!=0) num_gen_notPV_PU200_prompt_EE[idz]++;
			  num_gen_PU200_prompt_EE[idz]++;
			}

            dtsig_cut=0, dtsig_cut_vtx=0;
            dt_cut=0, dt_cut_vtx=0;
		  }
		}
		for(int idz=0; idz<track_pv_dz_cut_EE.size(); idz++) {
		  // (muon, track)
		  reliso_PU200_prompt_EE[idz] = sumiso_PU200_prompt_EE.at(idz)/muon_pt_->at(im);
		  reliso_PU200_prompt_1sigma_EE[idz] = sumiso_PU200_prompt_1sigma_EE[idz]/muon_pt_->at(im);
          reliso_PU200_prompt_2sigma_EE[idz] = sumiso_PU200_prompt_2sigma_EE[idz]/muon_pt_->at(im);
          reliso_PU200_prompt_3sigma_EE[idz] = sumiso_PU200_prompt_3sigma_EE[idz]/muon_pt_->at(im);
          reliso_PU200_prompt_4sigma_EE[idz] = sumiso_PU200_prompt_4sigma_EE[idz]/muon_pt_->at(im);
          reliso_PU200_prompt_genMatched_EE[idz] = sumiso_PU200_prompt_genMatched_EE[idz]/muon_pt_->at(im);
          // (PV, track)
		  reliso_PU200_prompt_EE_vtx[idz] = sumiso_PU200_prompt_EE_vtx.at(idz)/muon_pt_->at(im);
		  reliso_PU200_prompt_1sigma_EE_vtx[idz] = sumiso_PU200_prompt_1sigma_EE_vtx[idz]/muon_pt_->at(im);
          reliso_PU200_prompt_2sigma_EE_vtx[idz] = sumiso_PU200_prompt_2sigma_EE_vtx[idz]/muon_pt_->at(im);
          reliso_PU200_prompt_3sigma_EE_vtx[idz] = sumiso_PU200_prompt_3sigma_EE_vtx[idz]/muon_pt_->at(im);
          reliso_PU200_prompt_4sigma_EE_vtx[idz] = sumiso_PU200_prompt_4sigma_EE_vtx[idz]/muon_pt_->at(im);
          reliso_PU200_prompt_genMatched_EE_vtx[idz] = sumiso_PU200_prompt_genMatched_EE_vtx[idz]/muon_pt_->at(im);

		  // Store value of relative isolation
		  // (muon, track)
		  list_h_PU200_prompt_EE.at(idz)->Fill(reliso_PU200_prompt_EE.at(idz));
		  list_h_PU200_prompt_1sigma_EE.at(idz)->Fill(reliso_PU200_prompt_1sigma_EE.at(idz));
          list_h_PU200_prompt_2sigma_EE.at(idz)->Fill(reliso_PU200_prompt_2sigma_EE.at(idz));
          list_h_PU200_prompt_3sigma_EE.at(idz)->Fill(reliso_PU200_prompt_3sigma_EE.at(idz));
          list_h_PU200_prompt_4sigma_EE.at(idz)->Fill(reliso_PU200_prompt_4sigma_EE.at(idz));
          list_h_PU200_prompt_genMatched_EE.at(idz)->Fill(reliso_PU200_prompt_genMatched_EE.at(idz));
          // (PV, track)
		  list_h_PU200_prompt_EE_vtx.at(idz)->Fill(reliso_PU200_prompt_EE_vtx.at(idz));
		  list_h_PU200_prompt_1sigma_EE_vtx.at(idz)->Fill(reliso_PU200_prompt_1sigma_EE_vtx.at(idz));
          list_h_PU200_prompt_2sigma_EE_vtx.at(idz)->Fill(reliso_PU200_prompt_2sigma_EE_vtx.at(idz));
          list_h_PU200_prompt_3sigma_EE_vtx.at(idz)->Fill(reliso_PU200_prompt_3sigma_EE_vtx.at(idz));
          list_h_PU200_prompt_4sigma_EE_vtx.at(idz)->Fill(reliso_PU200_prompt_4sigma_EE_vtx.at(idz));
          list_h_PU200_prompt_genMatched_EE_vtx.at(idz)->Fill(reliso_PU200_prompt_genMatched_EE_vtx.at(idz));
		}

		// Initialize
		for(int idz=0; idz<track_pv_dz_cut_EE.size(); idz++) {
		  // (muon, track)
		  sumiso_PU200_prompt_EE[idz]=0;
          sumiso_PU200_prompt_1sigma_EE[idz]=0;
          sumiso_PU200_prompt_2sigma_EE[idz]=0;
          sumiso_PU200_prompt_3sigma_EE[idz]=0;
          sumiso_PU200_prompt_4sigma_EE[idz]=0;
          sumiso_PU200_prompt_genMatched_EE[idz]=0;
          reliso_PU200_prompt_EE[idz]=0;
          reliso_PU200_prompt_1sigma_EE[idz]=0;
          reliso_PU200_prompt_2sigma_EE[idz]=0;
          reliso_PU200_prompt_3sigma_EE[idz]=0;
          reliso_PU200_prompt_4sigma_EE[idz]=0;
          reliso_PU200_prompt_genMatched_EE[idz]=0;
          // (PV, track)
		  sumiso_PU200_prompt_EE_vtx[idz]=0;
          sumiso_PU200_prompt_1sigma_EE_vtx[idz]=0;
          sumiso_PU200_prompt_2sigma_EE_vtx[idz]=0;
          sumiso_PU200_prompt_3sigma_EE_vtx[idz]=0;
          sumiso_PU200_prompt_4sigma_EE_vtx[idz]=0;
          sumiso_PU200_prompt_genMatched_EE_vtx[idz]=0;
          reliso_PU200_prompt_EE_vtx[idz]=0;
          reliso_PU200_prompt_1sigma_EE_vtx[idz]=0;
          reliso_PU200_prompt_2sigma_EE_vtx[idz]=0;
          reliso_PU200_prompt_3sigma_EE_vtx[idz]=0;
          reliso_PU200_prompt_4sigma_EE_vtx[idz]=0;
          reliso_PU200_prompt_genMatched_EE_vtx[idz]=0;

		}
	  }
	} // End of muon loop
  } // End of event loop


  /////////////////////////
  //// PU200 nonprompt ////
  /////////////////////////
  

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

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
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
		  // dz loop
		  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EB.at(idz)) continue;
		    if(flag_track_PVweight) {
			  if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		    }
		    // noMTD case
		    sumiso_PU200_nonprompt_EB[idz] += track_pt_->at(im).at(it);

			// MTD case
			// (muon,track)
            if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
            }
            if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_PU200_nonprompt_1sigma_EB[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_PU200_nonprompt_2sigma_EB[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_PU200_nonprompt_3sigma_EB[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_PU200_nonprompt_4sigma_EB[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_PU200_nonprompt_1sigma_EB[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_nonprompt_2sigma_EB[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_nonprompt_3sigma_EB[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_nonprompt_4sigma_EB[idz]+= track_pt_->at(im).at(it);
            }
            // (PV, track)
			if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
			}
			if(dtsig_cut_vtx<1.0 && dtsig_cut_vtx>0) sumiso_PU200_nonprompt_1sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<2.0 && dtsig_cut_vtx>0) sumiso_PU200_nonprompt_2sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<3.0 && dtsig_cut_vtx>0) sumiso_PU200_nonprompt_3sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<4.0 && dtsig_cut_vtx>0) sumiso_PU200_nonprompt_4sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut_vtx==0 && (vtx_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_PU200_nonprompt_1sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_nonprompt_2sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_nonprompt_3sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_nonprompt_4sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
            }


			// GEN case
			if(track_genMatched_->at(im).at(it)==1) {
			  sumiso_PU200_nonprompt_genMatched_EB[idz] += track_pt_->at(im).at(it);
			  if(track_bx_->at(im).at(it)!=0 || track_evtId_->at(im).at(it)!=0) num_gen_notPV_PU200_nonprompt_EB[idz]++;
			  num_gen_PU200_nonprompt_EB[idz]++;
			}

            dtsig_cut=0, dtsig_cut_vtx=0;
            dt_cut=0, dt_cut_vtx=0;
		  }
		}
		for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
		  // (muon, track)
		  reliso_PU200_nonprompt_EB[idz] = sumiso_PU200_nonprompt_EB.at(idz)/muon_pt_->at(im);
		  reliso_PU200_nonprompt_1sigma_EB[idz] = sumiso_PU200_nonprompt_1sigma_EB.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_2sigma_EB[idz] = sumiso_PU200_nonprompt_2sigma_EB.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_3sigma_EB[idz] = sumiso_PU200_nonprompt_3sigma_EB.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_4sigma_EB[idz] = sumiso_PU200_nonprompt_4sigma_EB.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_genMatched_EB[idz] = sumiso_PU200_nonprompt_genMatched_EB.at(idz)/muon_pt_->at(im);
		  // (PV, track)
		  reliso_PU200_nonprompt_EB_vtx[idz] = sumiso_PU200_nonprompt_EB_vtx.at(idz)/muon_pt_->at(im);
		  reliso_PU200_nonprompt_1sigma_EB_vtx[idz] = sumiso_PU200_nonprompt_1sigma_EB_vtx.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_2sigma_EB_vtx[idz] = sumiso_PU200_nonprompt_2sigma_EB_vtx.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_3sigma_EB_vtx[idz] = sumiso_PU200_nonprompt_3sigma_EB_vtx.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_4sigma_EB_vtx[idz] = sumiso_PU200_nonprompt_4sigma_EB_vtx.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_genMatched_EB_vtx[idz] = sumiso_PU200_nonprompt_genMatched_EB_vtx.at(idz)/muon_pt_->at(im);

		  // Store value of relative isolation
		  // (muon, track)
		  list_h_PU200_nonprompt_EB.at(idz)->Fill(reliso_PU200_nonprompt_EB.at(idz));
		  list_h_PU200_nonprompt_1sigma_EB.at(idz)->Fill(reliso_PU200_nonprompt_1sigma_EB.at(idz));
          list_h_PU200_nonprompt_2sigma_EB.at(idz)->Fill(reliso_PU200_nonprompt_2sigma_EB.at(idz));
          list_h_PU200_nonprompt_3sigma_EB.at(idz)->Fill(reliso_PU200_nonprompt_3sigma_EB.at(idz));
          list_h_PU200_nonprompt_4sigma_EB.at(idz)->Fill(reliso_PU200_nonprompt_4sigma_EB.at(idz));
          list_h_PU200_nonprompt_genMatched_EB.at(idz)->Fill(reliso_PU200_nonprompt_genMatched_EB.at(idz));
          // (PV, track)
		  list_h_PU200_nonprompt_EB_vtx.at(idz)->Fill(reliso_PU200_nonprompt_EB_vtx.at(idz));
		  list_h_PU200_nonprompt_1sigma_EB_vtx.at(idz)->Fill(reliso_PU200_nonprompt_1sigma_EB_vtx.at(idz));
          list_h_PU200_nonprompt_2sigma_EB_vtx.at(idz)->Fill(reliso_PU200_nonprompt_2sigma_EB_vtx.at(idz));
          list_h_PU200_nonprompt_3sigma_EB_vtx.at(idz)->Fill(reliso_PU200_nonprompt_3sigma_EB_vtx.at(idz));
          list_h_PU200_nonprompt_4sigma_EB_vtx.at(idz)->Fill(reliso_PU200_nonprompt_4sigma_EB_vtx.at(idz));
          list_h_PU200_nonprompt_genMatched_EB_vtx.at(idz)->Fill(reliso_PU200_nonprompt_genMatched_EB_vtx.at(idz));
		}

		// Initialize
		for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
		  // (muon, track)
		  sumiso_PU200_nonprompt_EB[idz]=0;
          sumiso_PU200_nonprompt_1sigma_EB[idz]=0;
          sumiso_PU200_nonprompt_2sigma_EB[idz]=0;
          sumiso_PU200_nonprompt_3sigma_EB[idz]=0;
          sumiso_PU200_nonprompt_4sigma_EB[idz]=0;
          sumiso_PU200_nonprompt_genMatched_EB[idz]=0;
          reliso_PU200_nonprompt_EB[idz]=0;
          reliso_PU200_nonprompt_1sigma_EB[idz]=0;
          reliso_PU200_nonprompt_2sigma_EB[idz]=0;
          reliso_PU200_nonprompt_3sigma_EB[idz]=0;
          reliso_PU200_nonprompt_4sigma_EB[idz]=0;
          reliso_PU200_nonprompt_genMatched_EB[idz]=0;
		  // (PV, track)
		  sumiso_PU200_nonprompt_EB_vtx[idz]=0;
          sumiso_PU200_nonprompt_1sigma_EB_vtx[idz]=0;
          sumiso_PU200_nonprompt_2sigma_EB_vtx[idz]=0;
          sumiso_PU200_nonprompt_3sigma_EB_vtx[idz]=0;
          sumiso_PU200_nonprompt_4sigma_EB_vtx[idz]=0;
          sumiso_PU200_nonprompt_genMatched_EB_vtx[idz]=0;
          reliso_PU200_nonprompt_EB_vtx[idz]=0;
          reliso_PU200_nonprompt_1sigma_EB_vtx[idz]=0;
          reliso_PU200_nonprompt_2sigma_EB_vtx[idz]=0;
          reliso_PU200_nonprompt_3sigma_EB_vtx[idz]=0;
          reliso_PU200_nonprompt_4sigma_EB_vtx[idz]=0;
          reliso_PU200_nonprompt_genMatched_EB_vtx[idz]=0;

		}
	  }

	  /////////////////////////////
	  // nonprompt muon - Endcap //
	  /////////////////////////////
	  else if(muon_prompt_->at(im)==0 && muon_isBarrel_->at(im)==0) {

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
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
		  // dz loop
		  for(int idz=0; idz<track_pv_dz_cut_EE.size(); idz++) {
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EE.at(idz)) continue;
		    if(flag_track_PVweight) {
		  	  if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		    }
		    // noMTD case
		    sumiso_PU200_nonprompt_EE[idz] += track_pt_->at(im).at(it);

			// MTD case
            if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
            }
            if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_PU200_nonprompt_1sigma_EE[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_PU200_nonprompt_2sigma_EE[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_PU200_nonprompt_3sigma_EE[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_PU200_nonprompt_4sigma_EE[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_PU200_nonprompt_1sigma_EE[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_nonprompt_2sigma_EE[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_nonprompt_3sigma_EE[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_nonprompt_4sigma_EE[idz]+= track_pt_->at(im).at(it);
            }
            // (PV, track)
			if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
			}
			if(dtsig_cut_vtx<1.0 && dtsig_cut_vtx>0) sumiso_PU200_nonprompt_1sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<2.0 && dtsig_cut_vtx>0) sumiso_PU200_nonprompt_2sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<3.0 && dtsig_cut_vtx>0) sumiso_PU200_nonprompt_3sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<4.0 && dtsig_cut_vtx>0) sumiso_PU200_nonprompt_4sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut_vtx==0 && (vtx_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_PU200_nonprompt_1sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_nonprompt_2sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_nonprompt_3sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_PU200_nonprompt_4sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
            }


			// GEN case
			if(track_genMatched_->at(im).at(it)==1) {
			  sumiso_PU200_nonprompt_genMatched_EE[idz] += track_pt_->at(im).at(it);
			  if(track_bx_->at(im).at(it)!=0 || track_evtId_->at(im).at(it)!=0) num_gen_notPV_PU200_nonprompt_EE[idz]++;
			  num_gen_PU200_nonprompt_EE[idz]++;
			}

            dtsig_cut=0, dtsig_cut_vtx=0;
            dt_cut=0, dt_cut=0;
		  }
		}
		for(int idz=0; idz<track_pv_dz_cut_EE.size(); idz++) {
		  // (muon, track)
		  reliso_PU200_nonprompt_EE[idz] = sumiso_PU200_nonprompt_EE.at(idz)/muon_pt_->at(im);
		  reliso_PU200_nonprompt_1sigma_EE[idz] = sumiso_PU200_nonprompt_1sigma_EE.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_2sigma_EE[idz] = sumiso_PU200_nonprompt_2sigma_EE.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_3sigma_EE[idz] = sumiso_PU200_nonprompt_3sigma_EE.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_4sigma_EE[idz] = sumiso_PU200_nonprompt_4sigma_EE.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_genMatched_EE[idz] = sumiso_PU200_nonprompt_genMatched_EE.at(idz)/muon_pt_->at(im);
		  // (PV, track)
		  reliso_PU200_nonprompt_EE_vtx[idz] = sumiso_PU200_nonprompt_EE_vtx.at(idz)/muon_pt_->at(im);
		  reliso_PU200_nonprompt_1sigma_EE_vtx[idz] = sumiso_PU200_nonprompt_1sigma_EE_vtx.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_2sigma_EE_vtx[idz] = sumiso_PU200_nonprompt_2sigma_EE_vtx.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_3sigma_EE_vtx[idz] = sumiso_PU200_nonprompt_3sigma_EE_vtx.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_4sigma_EE_vtx[idz] = sumiso_PU200_nonprompt_4sigma_EE_vtx.at(idz)/muon_pt_->at(im);
          reliso_PU200_nonprompt_genMatched_EE_vtx[idz] = sumiso_PU200_nonprompt_genMatched_EE_vtx.at(idz)/muon_pt_->at(im);

		  // Store value of relative isolation
		  // (muon, track)
		  list_h_PU200_nonprompt_EE.at(idz)->Fill(reliso_PU200_nonprompt_EE.at(idz));
		  list_h_PU200_nonprompt_1sigma_EE.at(idz)->Fill(reliso_PU200_nonprompt_1sigma_EE.at(idz));
          list_h_PU200_nonprompt_2sigma_EE.at(idz)->Fill(reliso_PU200_nonprompt_2sigma_EE.at(idz));
          list_h_PU200_nonprompt_3sigma_EE.at(idz)->Fill(reliso_PU200_nonprompt_3sigma_EE.at(idz));
          list_h_PU200_nonprompt_4sigma_EE.at(idz)->Fill(reliso_PU200_nonprompt_4sigma_EE.at(idz));
          list_h_PU200_nonprompt_genMatched_EE.at(idz)->Fill(reliso_PU200_nonprompt_genMatched_EE.at(idz));
          // (PV, track)
		  list_h_PU200_nonprompt_EE_vtx.at(idz)->Fill(reliso_PU200_nonprompt_EE_vtx.at(idz));
		  list_h_PU200_nonprompt_1sigma_EE_vtx.at(idz)->Fill(reliso_PU200_nonprompt_1sigma_EE_vtx.at(idz));
          list_h_PU200_nonprompt_2sigma_EE_vtx.at(idz)->Fill(reliso_PU200_nonprompt_2sigma_EE_vtx.at(idz));
          list_h_PU200_nonprompt_3sigma_EE_vtx.at(idz)->Fill(reliso_PU200_nonprompt_3sigma_EE_vtx.at(idz));
          list_h_PU200_nonprompt_4sigma_EE_vtx.at(idz)->Fill(reliso_PU200_nonprompt_4sigma_EE_vtx.at(idz));
          list_h_PU200_nonprompt_genMatched_EE_vtx.at(idz)->Fill(reliso_PU200_nonprompt_genMatched_EE_vtx.at(idz));

		}

		// Initialize
		for(int idz=0; idz<track_pv_dz_cut_EE.size(); idz++) {
		  // (muon, track)
		  sumiso_PU200_nonprompt_EE[idz]=0;
          sumiso_PU200_nonprompt_1sigma_EE[idz]=0;
          sumiso_PU200_nonprompt_2sigma_EE[idz]=0;
          sumiso_PU200_nonprompt_3sigma_EE[idz]=0;
          sumiso_PU200_nonprompt_4sigma_EE[idz]=0;
          sumiso_PU200_nonprompt_genMatched_EE[idz]=0;
          reliso_PU200_nonprompt_EE[idz]=0;
          reliso_PU200_nonprompt_1sigma_EE[idz]=0;
          reliso_PU200_nonprompt_2sigma_EE[idz]=0;
          reliso_PU200_nonprompt_3sigma_EE[idz]=0;
          reliso_PU200_nonprompt_4sigma_EE[idz]=0;
          reliso_PU200_nonprompt_genMatched_EE[idz]=0;
          // (PV, track)
		  sumiso_PU200_nonprompt_EE_vtx[idz]=0;
          sumiso_PU200_nonprompt_1sigma_EE_vtx[idz]=0;
          sumiso_PU200_nonprompt_2sigma_EE_vtx[idz]=0;
          sumiso_PU200_nonprompt_3sigma_EE_vtx[idz]=0;
          sumiso_PU200_nonprompt_4sigma_EE_vtx[idz]=0;
          sumiso_PU200_nonprompt_genMatched_EE_vtx[idz]=0;
          reliso_PU200_nonprompt_EE_vtx[idz]=0;
          reliso_PU200_nonprompt_1sigma_EE_vtx[idz]=0;
          reliso_PU200_nonprompt_2sigma_EE_vtx[idz]=0;
          reliso_PU200_nonprompt_3sigma_EE_vtx[idz]=0;
          reliso_PU200_nonprompt_4sigma_EE_vtx[idz]=0;
          reliso_PU200_nonprompt_genMatched_EE_vtx[idz]=0;

		}
	  }
	} // End of muon loop
  } // End of event loop
	

  //////////////////////
  //// noPU prompt ////
  //////////////////////
  

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

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
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
		  // dz loop
		  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EB.at(idz)) continue;
		    if(flag_track_PVweight) {
			  if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		    }
		    // noMTD case
		    sumiso_noPU_prompt_EB[idz] += track_pt_->at(im).at(it);

			// MTD case
			// (muon, track)
            if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
            }
            if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_noPU_prompt_1sigma_EB[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_noPU_prompt_2sigma_EB[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_noPU_prompt_3sigma_EB[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_noPU_prompt_4sigma_EB[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_noPU_prompt_1sigma_EB[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_prompt_2sigma_EB[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_prompt_3sigma_EB[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_prompt_4sigma_EB[idz]+= track_pt_->at(im).at(it);
            }
            // (PV, track)
			if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
			}
			if(dtsig_cut_vtx<1.0 && dtsig_cut_vtx>0) sumiso_noPU_prompt_1sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<2.0 && dtsig_cut_vtx>0) sumiso_noPU_prompt_2sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<3.0 && dtsig_cut_vtx>0) sumiso_noPU_prompt_3sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<4.0 && dtsig_cut_vtx>0) sumiso_noPU_prompt_4sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut_vtx==0 && (vtx_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_noPU_prompt_1sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_prompt_2sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_prompt_3sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_prompt_4sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
            }


			// GEN case
			if(track_genMatched_->at(im).at(it)==1) {
			  sumiso_noPU_prompt_genMatched_EB[idz] += track_pt_->at(im).at(it);
			  if(track_bx_->at(im).at(it)!=0 || track_evtId_->at(im).at(it)!=0) num_gen_notPV_noPU_prompt_EB[idz]++;
			  num_gen_noPU_prompt_EB[idz]++;
			}

            dtsig_cut=0, dtsig_cut_vtx=0;
            dt_cut=0, dt_cut_vtx=0;
		  }
		}
		for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
		  // (muon, track)
		  reliso_noPU_prompt_EB[idz] = sumiso_noPU_prompt_EB.at(idz)/muon_pt_->at(im);
		  reliso_noPU_prompt_1sigma_EB[idz] = sumiso_noPU_prompt_1sigma_EB.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_2sigma_EB[idz] = sumiso_noPU_prompt_2sigma_EB.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_3sigma_EB[idz] = sumiso_noPU_prompt_3sigma_EB.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_4sigma_EB[idz] = sumiso_noPU_prompt_4sigma_EB.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_genMatched_EB[idz] = sumiso_noPU_prompt_genMatched_EB.at(idz)/muon_pt_->at(im);
		  // (PV, track)
		  reliso_noPU_prompt_EB_vtx[idz] = sumiso_noPU_prompt_EB_vtx.at(idz)/muon_pt_->at(im);
		  reliso_noPU_prompt_1sigma_EB_vtx[idz] = sumiso_noPU_prompt_1sigma_EB_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_2sigma_EB_vtx[idz] = sumiso_noPU_prompt_2sigma_EB_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_3sigma_EB_vtx[idz] = sumiso_noPU_prompt_3sigma_EB_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_4sigma_EB_vtx[idz] = sumiso_noPU_prompt_4sigma_EB_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_genMatched_EB_vtx[idz] = sumiso_noPU_prompt_genMatched_EB_vtx.at(idz)/muon_pt_->at(im);

		  // Store value of relative isolation
		  // (muon, track)
		  list_h_noPU_prompt_EB.at(idz)->Fill(reliso_noPU_prompt_EB.at(idz));
		  list_h_noPU_prompt_1sigma_EB.at(idz)->Fill(reliso_noPU_prompt_1sigma_EB.at(idz));
          list_h_noPU_prompt_2sigma_EB.at(idz)->Fill(reliso_noPU_prompt_2sigma_EB.at(idz));
          list_h_noPU_prompt_3sigma_EB.at(idz)->Fill(reliso_noPU_prompt_3sigma_EB.at(idz));
          list_h_noPU_prompt_4sigma_EB.at(idz)->Fill(reliso_noPU_prompt_4sigma_EB.at(idz));
          list_h_noPU_prompt_genMatched_EB.at(idz)->Fill(reliso_noPU_prompt_genMatched_EB.at(idz));
		  // (PV, track)
		  list_h_noPU_prompt_EB_vtx.at(idz)->Fill(reliso_noPU_prompt_EB_vtx.at(idz));
		  list_h_noPU_prompt_1sigma_EB_vtx.at(idz)->Fill(reliso_noPU_prompt_1sigma_EB_vtx.at(idz));
          list_h_noPU_prompt_2sigma_EB_vtx.at(idz)->Fill(reliso_noPU_prompt_2sigma_EB_vtx.at(idz));
          list_h_noPU_prompt_3sigma_EB_vtx.at(idz)->Fill(reliso_noPU_prompt_3sigma_EB_vtx.at(idz));
          list_h_noPU_prompt_4sigma_EB_vtx.at(idz)->Fill(reliso_noPU_prompt_4sigma_EB_vtx.at(idz));
          list_h_noPU_prompt_genMatched_EB_vtx.at(idz)->Fill(reliso_noPU_prompt_genMatched_EB_vtx.at(idz));
		}

		// Initialize
		for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
		  // (muon, track)
		  sumiso_noPU_prompt_EB[idz]=0;
          sumiso_noPU_prompt_1sigma_EB[idz]=0;
          sumiso_noPU_prompt_2sigma_EB[idz]=0;
          sumiso_noPU_prompt_3sigma_EB[idz]=0;
          sumiso_noPU_prompt_4sigma_EB[idz]=0;
          sumiso_noPU_prompt_genMatched_EB[idz]=0;
          reliso_noPU_prompt_EB[idz]=0;
          reliso_noPU_prompt_1sigma_EB[idz]=0;
          reliso_noPU_prompt_2sigma_EB[idz]=0;
          reliso_noPU_prompt_3sigma_EB[idz]=0;
          reliso_noPU_prompt_4sigma_EB[idz]=0;
          reliso_noPU_prompt_genMatched_EB[idz]=0;
		  // (PV, track)
		  sumiso_noPU_prompt_EB_vtx[idz]=0;
          sumiso_noPU_prompt_1sigma_EB_vtx[idz]=0;
          sumiso_noPU_prompt_2sigma_EB_vtx[idz]=0;
          sumiso_noPU_prompt_3sigma_EB_vtx[idz]=0;
          sumiso_noPU_prompt_4sigma_EB_vtx[idz]=0;
          sumiso_noPU_prompt_genMatched_EB_vtx[idz]=0;
          reliso_noPU_prompt_EB_vtx[idz]=0;
          reliso_noPU_prompt_1sigma_EB_vtx[idz]=0;
          reliso_noPU_prompt_2sigma_EB_vtx[idz]=0;
          reliso_noPU_prompt_3sigma_EB_vtx[idz]=0;
          reliso_noPU_prompt_4sigma_EB_vtx[idz]=0;
          reliso_noPU_prompt_genMatched_EB_vtx[idz]=0;

		}
	  }

	  //////////////////////////
	  // prompt muon - Endcap //
	  //////////////////////////
	  else if(muon_prompt_->at(im)==1 && muon_isBarrel_->at(im)==0) {

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
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
		  // dz loop
		  for(int idz=0; idz<track_pv_dz_cut_EE.size(); idz++) {
		    if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EE.at(idz)) continue;
		    if(flag_track_PVweight) {
			  if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		    }
		    // noMTD case
		    sumiso_noPU_prompt_EE[idz] += track_pt_->at(im).at(it);

			// MTD case
            // (muon, track)
            if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
            }
            if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_noPU_prompt_1sigma_EE[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_noPU_prompt_2sigma_EE[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_noPU_prompt_3sigma_EE[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_noPU_prompt_4sigma_EE[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_noPU_prompt_1sigma_EE[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_prompt_2sigma_EE[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_prompt_3sigma_EE[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_prompt_4sigma_EE[idz]+= track_pt_->at(im).at(it);
            }
            // (PV, track)
			if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
			}
			if(dtsig_cut_vtx<1.0 && dtsig_cut_vtx>0) sumiso_noPU_prompt_1sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<2.0 && dtsig_cut_vtx>0) sumiso_noPU_prompt_2sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<3.0 && dtsig_cut_vtx>0) sumiso_noPU_prompt_3sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<4.0 && dtsig_cut_vtx>0) sumiso_noPU_prompt_4sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut_vtx==0 && (vtx_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_noPU_prompt_1sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_prompt_2sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_prompt_3sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_prompt_4sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
            }


			// GEN case
			if(track_genMatched_->at(im).at(it)==1) {
			  sumiso_noPU_prompt_genMatched_EE[idz] += track_pt_->at(im).at(it);
			  if(track_bx_->at(im).at(it)!=0 || track_evtId_->at(im).at(it)!=0) num_gen_notPV_noPU_prompt_EE[idz]++;
			  num_gen_noPU_prompt_EE[idz]++;
			}

            dtsig_cut=0, dtsig_cut_vtx=0;
            dt_cut=0, dt_cut_vtx=0;
		  }
		}
		for(int idz=0; idz<track_pv_dz_cut_EE.size(); idz++) {
		  // (muon, track)
		  reliso_noPU_prompt_EE[idz] = sumiso_noPU_prompt_EE.at(idz)/muon_pt_->at(im);
		  reliso_noPU_prompt_1sigma_EE[idz] = sumiso_noPU_prompt_1sigma_EE.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_2sigma_EE[idz] = sumiso_noPU_prompt_2sigma_EE.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_3sigma_EE[idz] = sumiso_noPU_prompt_3sigma_EE.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_4sigma_EE[idz] = sumiso_noPU_prompt_4sigma_EE.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_genMatched_EE[idz] = sumiso_noPU_prompt_genMatched_EE.at(idz)/muon_pt_->at(im);
		  // (PV, track)
		  reliso_noPU_prompt_EE_vtx[idz] = sumiso_noPU_prompt_EE_vtx.at(idz)/muon_pt_->at(im);
		  reliso_noPU_prompt_1sigma_EE_vtx[idz] = sumiso_noPU_prompt_1sigma_EE_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_2sigma_EE_vtx[idz] = sumiso_noPU_prompt_2sigma_EE_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_3sigma_EE_vtx[idz] = sumiso_noPU_prompt_3sigma_EE_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_4sigma_EE_vtx[idz] = sumiso_noPU_prompt_4sigma_EE_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_prompt_genMatched_EE_vtx[idz] = sumiso_noPU_prompt_genMatched_EE_vtx.at(idz)/muon_pt_->at(im);

		  // Store value of relative isolation
		  // (muon, track)
		  list_h_noPU_prompt_EE.at(idz)->Fill(reliso_noPU_prompt_EE.at(idz));
		  list_h_noPU_prompt_1sigma_EE.at(idz)->Fill(reliso_noPU_prompt_1sigma_EE.at(idz));
          list_h_noPU_prompt_2sigma_EE.at(idz)->Fill(reliso_noPU_prompt_2sigma_EE.at(idz));
          list_h_noPU_prompt_3sigma_EE.at(idz)->Fill(reliso_noPU_prompt_3sigma_EE.at(idz));
          list_h_noPU_prompt_4sigma_EE.at(idz)->Fill(reliso_noPU_prompt_4sigma_EE.at(idz));
          list_h_noPU_prompt_genMatched_EE.at(idz)->Fill(reliso_noPU_prompt_genMatched_EE.at(idz));
		  // (PV, track)
		  list_h_noPU_prompt_EE_vtx.at(idz)->Fill(reliso_noPU_prompt_EE_vtx.at(idz));
		  list_h_noPU_prompt_1sigma_EE_vtx.at(idz)->Fill(reliso_noPU_prompt_1sigma_EE_vtx.at(idz));
          list_h_noPU_prompt_2sigma_EE_vtx.at(idz)->Fill(reliso_noPU_prompt_2sigma_EE_vtx.at(idz));
          list_h_noPU_prompt_3sigma_EE_vtx.at(idz)->Fill(reliso_noPU_prompt_3sigma_EE_vtx.at(idz));
          list_h_noPU_prompt_4sigma_EE_vtx.at(idz)->Fill(reliso_noPU_prompt_4sigma_EE_vtx.at(idz));
          list_h_noPU_prompt_genMatched_EE_vtx.at(idz)->Fill(reliso_noPU_prompt_genMatched_EE_vtx.at(idz));
		}

		// Initialize
		for(int idz=0; idz<track_pv_dz_cut_EE.size(); idz++) {
		  // (muon, track)
		  sumiso_noPU_prompt_EE[idz]=0;
          sumiso_noPU_prompt_1sigma_EE[idz]=0;
          sumiso_noPU_prompt_2sigma_EE[idz]=0;
          sumiso_noPU_prompt_3sigma_EE[idz]=0;
          sumiso_noPU_prompt_4sigma_EE[idz]=0;
          sumiso_noPU_prompt_genMatched_EE[idz]=0;
          reliso_noPU_prompt_EE[idz]=0;
          reliso_noPU_prompt_1sigma_EE[idz]=0;
          reliso_noPU_prompt_2sigma_EE[idz]=0;
          reliso_noPU_prompt_3sigma_EE[idz]=0;
          reliso_noPU_prompt_4sigma_EE[idz]=0;
          reliso_noPU_prompt_genMatched_EE[idz]=0;
		  // (PV, track)
		  sumiso_noPU_prompt_EE_vtx[idz]=0;
          sumiso_noPU_prompt_1sigma_EE_vtx[idz]=0;
          sumiso_noPU_prompt_2sigma_EE_vtx[idz]=0;
          sumiso_noPU_prompt_3sigma_EE_vtx[idz]=0;
          sumiso_noPU_prompt_4sigma_EE_vtx[idz]=0;
          sumiso_noPU_prompt_genMatched_EE_vtx[idz]=0;
          reliso_noPU_prompt_EE_vtx[idz]=0;
          reliso_noPU_prompt_1sigma_EE_vtx[idz]=0;
          reliso_noPU_prompt_2sigma_EE_vtx[idz]=0;
          reliso_noPU_prompt_3sigma_EE_vtx[idz]=0;
          reliso_noPU_prompt_4sigma_EE_vtx[idz]=0;
          reliso_noPU_prompt_genMatched_EE_vtx[idz]=0;

		}
	  }
	} // End of muon loop
  } // End of event loop


  /////////////////////////
  //// noPU nonprompt ////
  /////////////////////////
  

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

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
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
		  // dz loop
		  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
		    if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EB.at(idz)) continue;
		    if(flag_track_PVweight) {
			  if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		    }
		    // noMTD case
		    sumiso_noPU_nonprompt_EB[idz] += track_pt_->at(im).at(it);

			// MTD case
			// (muon, track)
            if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
            }
            if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_noPU_nonprompt_1sigma_EB[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_noPU_nonprompt_2sigma_EB[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_noPU_nonprompt_3sigma_EB[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_noPU_nonprompt_4sigma_EB[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_noPU_nonprompt_1sigma_EB[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_nonprompt_2sigma_EB[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_nonprompt_3sigma_EB[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_nonprompt_4sigma_EB[idz]+= track_pt_->at(im).at(it);
            }
            // (PV, track)
			if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
			}
			if(dtsig_cut_vtx<1.0 && dtsig_cut_vtx>0) sumiso_noPU_nonprompt_1sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<2.0 && dtsig_cut_vtx>0) sumiso_noPU_nonprompt_2sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<3.0 && dtsig_cut_vtx>0) sumiso_noPU_nonprompt_3sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<4.0 && dtsig_cut_vtx>0) sumiso_noPU_nonprompt_4sigma_EB_vtx[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut_vtx==0 && (vtx_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_noPU_nonprompt_1sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_nonprompt_2sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_nonprompt_3sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_nonprompt_4sigma_EB_vtx[idz]+= track_pt_->at(im).at(it);
            }


			// GEN case
			if(track_genMatched_->at(im).at(it)==1) {
			  sumiso_noPU_nonprompt_genMatched_EB[idz] += track_pt_->at(im).at(it);
			  if(track_bx_->at(im).at(it)!=0 || track_evtId_->at(im).at(it)!=0) num_gen_notPV_noPU_nonprompt_EB[idz]++;
			  num_gen_noPU_nonprompt_EB[idz]++;
			}

            dtsig_cut=0, dtsig_cut_vtx=0;
            dt_cut=0, dt_cut_vtx=0;
		  }
		}
		for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
		  // (muon, track)
		  reliso_noPU_nonprompt_EB[idz] = sumiso_noPU_nonprompt_EB.at(idz)/muon_pt_->at(im);
		  reliso_noPU_nonprompt_1sigma_EB[idz] = sumiso_noPU_nonprompt_1sigma_EB.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_2sigma_EB[idz] = sumiso_noPU_nonprompt_2sigma_EB.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_3sigma_EB[idz] = sumiso_noPU_nonprompt_3sigma_EB.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_4sigma_EB[idz] = sumiso_noPU_nonprompt_4sigma_EB.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_genMatched_EB[idz] = sumiso_noPU_nonprompt_genMatched_EB.at(idz)/muon_pt_->at(im);
		  // (PV, track)
		  reliso_noPU_nonprompt_EB_vtx[idz] = sumiso_noPU_nonprompt_EB_vtx.at(idz)/muon_pt_->at(im);
		  reliso_noPU_nonprompt_1sigma_EB_vtx[idz] = sumiso_noPU_nonprompt_1sigma_EB_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_2sigma_EB_vtx[idz] = sumiso_noPU_nonprompt_2sigma_EB_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_3sigma_EB_vtx[idz] = sumiso_noPU_nonprompt_3sigma_EB_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_4sigma_EB_vtx[idz] = sumiso_noPU_nonprompt_4sigma_EB_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_genMatched_EB_vtx[idz] = sumiso_noPU_nonprompt_genMatched_EB_vtx.at(idz)/muon_pt_->at(im);

		  // Store value of relative isolation
		  // (muon, track)
		  list_h_noPU_nonprompt_EB.at(idz)->Fill(reliso_noPU_nonprompt_EB.at(idz));
		  list_h_noPU_nonprompt_1sigma_EB.at(idz)->Fill(reliso_noPU_nonprompt_1sigma_EB.at(idz));
          list_h_noPU_nonprompt_2sigma_EB.at(idz)->Fill(reliso_noPU_nonprompt_2sigma_EB.at(idz));
          list_h_noPU_nonprompt_3sigma_EB.at(idz)->Fill(reliso_noPU_nonprompt_3sigma_EB.at(idz));
          list_h_noPU_nonprompt_4sigma_EB.at(idz)->Fill(reliso_noPU_nonprompt_4sigma_EB.at(idz));
          list_h_noPU_nonprompt_genMatched_EB.at(idz)->Fill(reliso_noPU_nonprompt_genMatched_EB.at(idz));
		  // (PV, track)
		  list_h_noPU_nonprompt_EB_vtx.at(idz)->Fill(reliso_noPU_nonprompt_EB_vtx.at(idz));
		  list_h_noPU_nonprompt_1sigma_EB_vtx.at(idz)->Fill(reliso_noPU_nonprompt_1sigma_EB_vtx.at(idz));
          list_h_noPU_nonprompt_2sigma_EB_vtx.at(idz)->Fill(reliso_noPU_nonprompt_2sigma_EB_vtx.at(idz));
          list_h_noPU_nonprompt_3sigma_EB_vtx.at(idz)->Fill(reliso_noPU_nonprompt_3sigma_EB_vtx.at(idz));
          list_h_noPU_nonprompt_4sigma_EB_vtx.at(idz)->Fill(reliso_noPU_nonprompt_4sigma_EB_vtx.at(idz));
          list_h_noPU_nonprompt_genMatched_EB_vtx.at(idz)->Fill(reliso_noPU_nonprompt_genMatched_EB_vtx.at(idz));
		}

		// Initialize
		for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
		  // (muon, track)
		  sumiso_noPU_nonprompt_EB[idz]=0;
          sumiso_noPU_nonprompt_1sigma_EB[idz]=0;
          sumiso_noPU_nonprompt_2sigma_EB[idz]=0;
          sumiso_noPU_nonprompt_3sigma_EB[idz]=0;
          sumiso_noPU_nonprompt_4sigma_EB[idz]=0;
          sumiso_noPU_nonprompt_genMatched_EB[idz]=0;
          reliso_noPU_nonprompt_EB[idz]=0;
          reliso_noPU_nonprompt_1sigma_EB[idz]=0;
          reliso_noPU_nonprompt_2sigma_EB[idz]=0;
          reliso_noPU_nonprompt_3sigma_EB[idz]=0;
          reliso_noPU_nonprompt_4sigma_EB[idz]=0;
          reliso_noPU_nonprompt_genMatched_EB[idz]=0;
		  // (PV, track)
		  sumiso_noPU_nonprompt_EB_vtx[idz]=0;
          sumiso_noPU_nonprompt_1sigma_EB_vtx[idz]=0;
          sumiso_noPU_nonprompt_2sigma_EB_vtx[idz]=0;
          sumiso_noPU_nonprompt_3sigma_EB_vtx[idz]=0;
          sumiso_noPU_nonprompt_4sigma_EB_vtx[idz]=0;
          sumiso_noPU_nonprompt_genMatched_EB_vtx[idz]=0;
          reliso_noPU_nonprompt_EB_vtx[idz]=0;
          reliso_noPU_nonprompt_1sigma_EB_vtx[idz]=0;
          reliso_noPU_nonprompt_2sigma_EB_vtx[idz]=0;
          reliso_noPU_nonprompt_3sigma_EB_vtx[idz]=0;
          reliso_noPU_nonprompt_4sigma_EB_vtx[idz]=0;
          reliso_noPU_nonprompt_genMatched_EB_vtx[idz]=0;

		}
	  }

	  /////////////////////////////
	  // nonprompt muon - Endcap //
	  /////////////////////////////
	  else if(muon_prompt_->at(im)==0 && muon_isBarrel_->at(im)==0) {

		if(flag_muon_status) {
	      if(muon_status_->at(im)==-99 || muon_status_->at(im)==999999) {
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
		  // dz loop
		  for(int idz=0; idz<track_pv_dz_cut_EE.size(); idz++) {
		    if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EE.at(idz)) continue;
		    if(flag_track_PVweight) {
			  if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		    }
		    // noMTD case
		    sumiso_noPU_nonprompt_EE[idz] += track_pt_->at(im).at(it);

			// MTD case
            if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
            }
            if(dtsig_cut<1.0 && dtsig_cut>0) sumiso_noPU_nonprompt_1sigma_EE[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<2.0 && dtsig_cut>0) sumiso_noPU_nonprompt_2sigma_EE[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<3.0 && dtsig_cut>0) sumiso_noPU_nonprompt_3sigma_EE[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut<4.0 && dtsig_cut>0) sumiso_noPU_nonprompt_4sigma_EE[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut==0 && (muon_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_noPU_nonprompt_1sigma_EE[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_nonprompt_2sigma_EE[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_nonprompt_3sigma_EE[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_nonprompt_4sigma_EE[idz]+= track_pt_->at(im).at(it);
            }
            // (PV, track)
			if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
			}
			if(dtsig_cut_vtx<1.0 && dtsig_cut_vtx>0) sumiso_noPU_nonprompt_1sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<2.0 && dtsig_cut_vtx>0) sumiso_noPU_nonprompt_2sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<3.0 && dtsig_cut_vtx>0) sumiso_noPU_nonprompt_3sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            if(dtsig_cut_vtx<4.0 && dtsig_cut_vtx>0) sumiso_noPU_nonprompt_4sigma_EE_vtx[idz] += track_pt_->at(im).at(it);
            else if(dtsig_cut_vtx==0 && (vtx_time_err_->at(im)<0 || track_time_err_->at(im).at(it)<0)) {
              sumiso_noPU_nonprompt_1sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_nonprompt_2sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_nonprompt_3sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
              sumiso_noPU_nonprompt_4sigma_EE_vtx[idz]+= track_pt_->at(im).at(it);
            }


			// GEN case
			if(track_genMatched_->at(im).at(it)==1) {
			  sumiso_noPU_nonprompt_genMatched_EE[idz] += track_pt_->at(im).at(it);
			  if(track_bx_->at(im).at(it)!=0 || track_evtId_->at(im).at(it)!=0) num_gen_notPV_noPU_nonprompt_EE[idz]++;
			  num_gen_noPU_nonprompt_EE[idz]++;
			}

            dtsig_cut=0, dtsig_cut_vtx=0;
            dt_cut=0, dt_cut_vtx=0;
		  }
		}
		for(int idz=0; idz<track_pv_dz_cut_EE.size(); idz++) {
		  // (muon, track)
		  reliso_noPU_nonprompt_EE[idz]        = sumiso_noPU_nonprompt_EE.at(idz)/muon_pt_->at(im);
		  reliso_noPU_nonprompt_1sigma_EE[idz] = sumiso_noPU_nonprompt_1sigma_EE.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_2sigma_EE[idz] = sumiso_noPU_nonprompt_2sigma_EE.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_3sigma_EE[idz] = sumiso_noPU_nonprompt_3sigma_EE.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_4sigma_EE[idz] = sumiso_noPU_nonprompt_4sigma_EE.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_genMatched_EE[idz] = sumiso_noPU_nonprompt_genMatched_EE.at(idz)/muon_pt_->at(im);
		  // (PV, track)
		  reliso_noPU_nonprompt_EE_vtx[idz]        = sumiso_noPU_nonprompt_EE_vtx.at(idz)/muon_pt_->at(im);
		  reliso_noPU_nonprompt_1sigma_EE_vtx[idz] = sumiso_noPU_nonprompt_1sigma_EE_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_2sigma_EE_vtx[idz] = sumiso_noPU_nonprompt_2sigma_EE_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_3sigma_EE_vtx[idz] = sumiso_noPU_nonprompt_3sigma_EE_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_4sigma_EE_vtx[idz] = sumiso_noPU_nonprompt_4sigma_EE_vtx.at(idz)/muon_pt_->at(im);
          reliso_noPU_nonprompt_genMatched_EE_vtx[idz] = sumiso_noPU_nonprompt_genMatched_EE_vtx.at(idz)/muon_pt_->at(im);

		  // Store value of relative isolation
		  // (muon, track)
		  list_h_noPU_nonprompt_EE.at(idz)->Fill(reliso_noPU_nonprompt_EE.at(idz));
		  list_h_noPU_nonprompt_1sigma_EE.at(idz)->Fill(reliso_noPU_nonprompt_1sigma_EE.at(idz));
          list_h_noPU_nonprompt_2sigma_EE.at(idz)->Fill(reliso_noPU_nonprompt_2sigma_EE.at(idz));
          list_h_noPU_nonprompt_3sigma_EE.at(idz)->Fill(reliso_noPU_nonprompt_3sigma_EE.at(idz));
          list_h_noPU_nonprompt_4sigma_EE.at(idz)->Fill(reliso_noPU_nonprompt_4sigma_EE.at(idz));
          list_h_noPU_nonprompt_genMatched_EE.at(idz)->Fill(reliso_noPU_nonprompt_genMatched_EE.at(idz));
		  // (muon, track)
		  list_h_noPU_nonprompt_EE_vtx.at(idz)->Fill(reliso_noPU_nonprompt_EE_vtx.at(idz));
		  list_h_noPU_nonprompt_1sigma_EE_vtx.at(idz)->Fill(reliso_noPU_nonprompt_1sigma_EE_vtx.at(idz));
          list_h_noPU_nonprompt_2sigma_EE_vtx.at(idz)->Fill(reliso_noPU_nonprompt_2sigma_EE_vtx.at(idz));
          list_h_noPU_nonprompt_3sigma_EE_vtx.at(idz)->Fill(reliso_noPU_nonprompt_3sigma_EE_vtx.at(idz));
          list_h_noPU_nonprompt_4sigma_EE_vtx.at(idz)->Fill(reliso_noPU_nonprompt_4sigma_EE_vtx.at(idz));
          list_h_noPU_nonprompt_genMatched_EE_vtx.at(idz)->Fill(reliso_noPU_nonprompt_genMatched_EE_vtx.at(idz));
		}

		// Initialize
		for(int idz=0; idz<track_pv_dz_cut_EE.size(); idz++) {
		  // (muon, track)
		  sumiso_noPU_nonprompt_EE[idz]=0;
          sumiso_noPU_nonprompt_1sigma_EE[idz]=0;
          sumiso_noPU_nonprompt_2sigma_EE[idz]=0;
          sumiso_noPU_nonprompt_3sigma_EE[idz]=0;
          sumiso_noPU_nonprompt_4sigma_EE[idz]=0;
          sumiso_noPU_nonprompt_genMatched_EE[idz]=0;
          reliso_noPU_nonprompt_EE[idz]=0;
          reliso_noPU_nonprompt_1sigma_EE[idz]=0;
          reliso_noPU_nonprompt_2sigma_EE[idz]=0;
          reliso_noPU_nonprompt_3sigma_EE[idz]=0;
          reliso_noPU_nonprompt_4sigma_EE[idz]=0;
          reliso_noPU_nonprompt_genMatched_EE[idz]=0;
		  // (PV, track)
		  sumiso_noPU_nonprompt_EE_vtx[idz]=0;
          sumiso_noPU_nonprompt_1sigma_EE_vtx[idz]=0;
          sumiso_noPU_nonprompt_2sigma_EE_vtx[idz]=0;
          sumiso_noPU_nonprompt_3sigma_EE_vtx[idz]=0;
          sumiso_noPU_nonprompt_4sigma_EE_vtx[idz]=0;
          sumiso_noPU_nonprompt_genMatched_EE_vtx[idz]=0;
          reliso_noPU_nonprompt_EE_vtx[idz]=0;
          reliso_noPU_nonprompt_1sigma_EE_vtx[idz]=0;
          reliso_noPU_nonprompt_2sigma_EE_vtx[idz]=0;
          reliso_noPU_nonprompt_3sigma_EE_vtx[idz]=0;
          reliso_noPU_nonprompt_4sigma_EE_vtx[idz]=0;
          reliso_noPU_nonprompt_genMatched_EE_vtx[idz]=0;

		}
	  }
	} // End of muon loop
  } // End of event loop




  ////////////////////////////////////////////
  // Define vectors to store values of eff. //
  ////////////////////////////////////////////

  // (muon, track)
  // PU200
  vector<vector<double>> prompt_eff_PU200_EB,           prompt_eff_PU200_EE;
  vector<vector<double>> prompt_eff_PU200_1sigma_EB,    prompt_eff_PU200_1sigma_EE;
  vector<vector<double>> prompt_eff_PU200_2sigma_EB,    prompt_eff_PU200_2sigma_EE;
  vector<vector<double>> prompt_eff_PU200_3sigma_EB,    prompt_eff_PU200_3sigma_EE;
  vector<vector<double>> prompt_eff_PU200_4sigma_EB,    prompt_eff_PU200_4sigma_EE;
  vector<vector<double>> prompt_eff_PU200_genMatched_EB,    prompt_eff_PU200_genMatched_EE;
  vector<vector<double>> nonprompt_eff_PU200_EB,        nonprompt_eff_PU200_EE;
  vector<vector<double>> nonprompt_eff_PU200_1sigma_EB, nonprompt_eff_PU200_1sigma_EE;
  vector<vector<double>> nonprompt_eff_PU200_2sigma_EB, nonprompt_eff_PU200_2sigma_EE;
  vector<vector<double>> nonprompt_eff_PU200_3sigma_EB, nonprompt_eff_PU200_3sigma_EE;
  vector<vector<double>> nonprompt_eff_PU200_4sigma_EB, nonprompt_eff_PU200_4sigma_EE;
  vector<vector<double>> nonprompt_eff_PU200_genMatched_EB, nonprompt_eff_PU200_genMatched_EE;
  vector<double> prompt_eff_PU200_EB_comp,              prompt_eff_PU200_EE_comp;
  vector<double> prompt_eff_PU200_1sigma_EB_comp,       prompt_eff_PU200_1sigma_EE_comp;
  vector<double> prompt_eff_PU200_2sigma_EB_comp,       prompt_eff_PU200_2sigma_EE_comp;
  vector<double> prompt_eff_PU200_3sigma_EB_comp,       prompt_eff_PU200_3sigma_EE_comp;
  vector<double> prompt_eff_PU200_4sigma_EB_comp,       prompt_eff_PU200_4sigma_EE_comp;
  vector<double> prompt_eff_PU200_genMatched_EB_comp,       prompt_eff_PU200_genMatched_EE_comp;
  vector<double> nonprompt_eff_PU200_EB_comp,           nonprompt_eff_PU200_EE_comp;
  vector<double> nonprompt_eff_PU200_1sigma_EB_comp,    nonprompt_eff_PU200_1sigma_EE_comp;
  vector<double> nonprompt_eff_PU200_2sigma_EB_comp,    nonprompt_eff_PU200_2sigma_EE_comp;
  vector<double> nonprompt_eff_PU200_3sigma_EB_comp,    nonprompt_eff_PU200_3sigma_EE_comp;
  vector<double> nonprompt_eff_PU200_4sigma_EB_comp,    nonprompt_eff_PU200_4sigma_EE_comp;
  vector<double> nonprompt_eff_PU200_genMatched_EB_comp,    nonprompt_eff_PU200_genMatched_EE_comp;
  // noPU
  vector<vector<double>> prompt_eff_noPU_EB,           prompt_eff_noPU_EE;
  vector<vector<double>> prompt_eff_noPU_1sigma_EB,    prompt_eff_noPU_1sigma_EE;
  vector<vector<double>> prompt_eff_noPU_2sigma_EB,    prompt_eff_noPU_2sigma_EE;
  vector<vector<double>> prompt_eff_noPU_3sigma_EB,    prompt_eff_noPU_3sigma_EE;
  vector<vector<double>> prompt_eff_noPU_4sigma_EB,    prompt_eff_noPU_4sigma_EE;
  vector<vector<double>> prompt_eff_noPU_genMatched_EB,    prompt_eff_noPU_genMatched_EE;
  vector<vector<double>> nonprompt_eff_noPU_EB,        nonprompt_eff_noPU_EE;
  vector<vector<double>> nonprompt_eff_noPU_1sigma_EB, nonprompt_eff_noPU_1sigma_EE;
  vector<vector<double>> nonprompt_eff_noPU_2sigma_EB, nonprompt_eff_noPU_2sigma_EE;
  vector<vector<double>> nonprompt_eff_noPU_3sigma_EB, nonprompt_eff_noPU_3sigma_EE;
  vector<vector<double>> nonprompt_eff_noPU_4sigma_EB, nonprompt_eff_noPU_4sigma_EE;
  vector<vector<double>> nonprompt_eff_noPU_genMatched_EB, nonprompt_eff_noPU_genMatched_EE;
  vector<double> prompt_eff_noPU_EB_comp,              prompt_eff_noPU_EE_comp;
  vector<double> prompt_eff_noPU_1sigma_EB_comp,       prompt_eff_noPU_1sigma_EE_comp;
  vector<double> prompt_eff_noPU_2sigma_EB_comp,       prompt_eff_noPU_2sigma_EE_comp;
  vector<double> prompt_eff_noPU_3sigma_EB_comp,       prompt_eff_noPU_3sigma_EE_comp;
  vector<double> prompt_eff_noPU_4sigma_EB_comp,       prompt_eff_noPU_4sigma_EE_comp;
  vector<double> prompt_eff_noPU_genMatched_EB_comp,       prompt_eff_noPU_genMatched_EE_comp;
  vector<double> nonprompt_eff_noPU_EB_comp,           nonprompt_eff_noPU_EE_comp;
  vector<double> nonprompt_eff_noPU_1sigma_EB_comp,    nonprompt_eff_noPU_1sigma_EE_comp;
  vector<double> nonprompt_eff_noPU_2sigma_EB_comp,    nonprompt_eff_noPU_2sigma_EE_comp;
  vector<double> nonprompt_eff_noPU_3sigma_EB_comp,    nonprompt_eff_noPU_3sigma_EE_comp;
  vector<double> nonprompt_eff_noPU_4sigma_EB_comp,    nonprompt_eff_noPU_4sigma_EE_comp;
  vector<double> nonprompt_eff_noPU_genMatched_EB_comp,    nonprompt_eff_noPU_genMatched_EE_comp;

  // (PV, track)
  // PU200
  vector<vector<double>> prompt_eff_PU200_EB_vtx,           prompt_eff_PU200_EE_vtx;
  vector<vector<double>> prompt_eff_PU200_1sigma_EB_vtx,    prompt_eff_PU200_1sigma_EE_vtx;
  vector<vector<double>> prompt_eff_PU200_2sigma_EB_vtx,    prompt_eff_PU200_2sigma_EE_vtx;
  vector<vector<double>> prompt_eff_PU200_3sigma_EB_vtx,    prompt_eff_PU200_3sigma_EE_vtx;
  vector<vector<double>> prompt_eff_PU200_4sigma_EB_vtx,    prompt_eff_PU200_4sigma_EE_vtx;
  vector<vector<double>> prompt_eff_PU200_genMatched_EB_vtx,    prompt_eff_PU200_genMatched_EE_vtx;
  vector<vector<double>> nonprompt_eff_PU200_EB_vtx,        nonprompt_eff_PU200_EE_vtx;
  vector<vector<double>> nonprompt_eff_PU200_1sigma_EB_vtx, nonprompt_eff_PU200_1sigma_EE_vtx;
  vector<vector<double>> nonprompt_eff_PU200_2sigma_EB_vtx, nonprompt_eff_PU200_2sigma_EE_vtx;
  vector<vector<double>> nonprompt_eff_PU200_3sigma_EB_vtx, nonprompt_eff_PU200_3sigma_EE_vtx;
  vector<vector<double>> nonprompt_eff_PU200_4sigma_EB_vtx, nonprompt_eff_PU200_4sigma_EE_vtx;
  vector<vector<double>> nonprompt_eff_PU200_genMatched_EB_vtx, nonprompt_eff_PU200_genMatched_EE_vtx;
  vector<double> prompt_eff_PU200_EB_vtx_comp,              prompt_eff_PU200_EE_vtx_comp;
  vector<double> prompt_eff_PU200_1sigma_EB_vtx_comp,       prompt_eff_PU200_1sigma_EE_vtx_comp;
  vector<double> prompt_eff_PU200_2sigma_EB_vtx_comp,       prompt_eff_PU200_2sigma_EE_vtx_comp;
  vector<double> prompt_eff_PU200_3sigma_EB_vtx_comp,       prompt_eff_PU200_3sigma_EE_vtx_comp;
  vector<double> prompt_eff_PU200_4sigma_EB_vtx_comp,       prompt_eff_PU200_4sigma_EE_vtx_comp;
  vector<double> prompt_eff_PU200_genMatched_EB_vtx_comp,       prompt_eff_PU200_genMatched_EE_vtx_comp;
  vector<double> nonprompt_eff_PU200_EB_vtx_comp,           nonprompt_eff_PU200_EE_vtx_comp;
  vector<double> nonprompt_eff_PU200_1sigma_EB_vtx_comp,    nonprompt_eff_PU200_1sigma_EE_vtx_comp;
  vector<double> nonprompt_eff_PU200_2sigma_EB_vtx_comp,    nonprompt_eff_PU200_2sigma_EE_vtx_comp;
  vector<double> nonprompt_eff_PU200_3sigma_EB_vtx_comp,    nonprompt_eff_PU200_3sigma_EE_vtx_comp;
  vector<double> nonprompt_eff_PU200_4sigma_EB_vtx_comp,    nonprompt_eff_PU200_4sigma_EE_vtx_comp;
  vector<double> nonprompt_eff_PU200_genMatched_EB_vtx_comp,    nonprompt_eff_PU200_genMatched_EE_vtx_comp;
  // noPU
  vector<vector<double>> prompt_eff_noPU_EB_vtx,           prompt_eff_noPU_EE_vtx;
  vector<vector<double>> prompt_eff_noPU_1sigma_EB_vtx,    prompt_eff_noPU_1sigma_EE_vtx;
  vector<vector<double>> prompt_eff_noPU_2sigma_EB_vtx,    prompt_eff_noPU_2sigma_EE_vtx;
  vector<vector<double>> prompt_eff_noPU_3sigma_EB_vtx,    prompt_eff_noPU_3sigma_EE_vtx;
  vector<vector<double>> prompt_eff_noPU_4sigma_EB_vtx,    prompt_eff_noPU_4sigma_EE_vtx;
  vector<vector<double>> prompt_eff_noPU_genMatched_EB_vtx,    prompt_eff_noPU_genMatched_EE_vtx;
  vector<vector<double>> nonprompt_eff_noPU_EB_vtx,        nonprompt_eff_noPU_EE_vtx;
  vector<vector<double>> nonprompt_eff_noPU_1sigma_EB_vtx, nonprompt_eff_noPU_1sigma_EE_vtx;
  vector<vector<double>> nonprompt_eff_noPU_2sigma_EB_vtx, nonprompt_eff_noPU_2sigma_EE_vtx;
  vector<vector<double>> nonprompt_eff_noPU_3sigma_EB_vtx, nonprompt_eff_noPU_3sigma_EE_vtx;
  vector<vector<double>> nonprompt_eff_noPU_4sigma_EB_vtx, nonprompt_eff_noPU_4sigma_EE_vtx;
  vector<vector<double>> nonprompt_eff_noPU_genMatched_EB_vtx, nonprompt_eff_noPU_genMatched_EE_vtx;
  vector<double> prompt_eff_noPU_EB_vtx_comp,              prompt_eff_noPU_EE_vtx_comp;
  vector<double> prompt_eff_noPU_1sigma_EB_vtx_comp,       prompt_eff_noPU_1sigma_EE_vtx_comp;
  vector<double> prompt_eff_noPU_2sigma_EB_vtx_comp,       prompt_eff_noPU_2sigma_EE_vtx_comp;
  vector<double> prompt_eff_noPU_3sigma_EB_vtx_comp,       prompt_eff_noPU_3sigma_EE_vtx_comp;
  vector<double> prompt_eff_noPU_4sigma_EB_vtx_comp,       prompt_eff_noPU_4sigma_EE_vtx_comp;
  vector<double> prompt_eff_noPU_genMatched_EB_vtx_comp,       prompt_eff_noPU_genMatched_EE_vtx_comp;
  vector<double> nonprompt_eff_noPU_EB_vtx_comp,           nonprompt_eff_noPU_EE_vtx_comp;
  vector<double> nonprompt_eff_noPU_1sigma_EB_vtx_comp,    nonprompt_eff_noPU_1sigma_EE_vtx_comp;
  vector<double> nonprompt_eff_noPU_2sigma_EB_vtx_comp,    nonprompt_eff_noPU_2sigma_EE_vtx_comp;
  vector<double> nonprompt_eff_noPU_3sigma_EB_vtx_comp,    nonprompt_eff_noPU_3sigma_EE_vtx_comp;
  vector<double> nonprompt_eff_noPU_4sigma_EB_vtx_comp,    nonprompt_eff_noPU_4sigma_EE_vtx_comp;
  vector<double> nonprompt_eff_noPU_genMatched_EB_vtx_comp,    nonprompt_eff_noPU_genMatched_EE_vtx_comp;


  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
    for(int i=0; i<nbin+1; i++) {
	  // (muon, track)
	  // PU200
	    // prompt
	  prompt_eff_PU200_EB_comp.emplace_back(list_h_PU200_prompt_EB.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_EB.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_EE_comp.emplace_back(list_h_PU200_prompt_EE.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_EE.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_1sigma_EB_comp.emplace_back(list_h_PU200_prompt_1sigma_EB.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_1sigma_EB.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_1sigma_EE_comp.emplace_back(list_h_PU200_prompt_1sigma_EE.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_1sigma_EE.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_2sigma_EB_comp.emplace_back(list_h_PU200_prompt_2sigma_EB.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_2sigma_EB.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_2sigma_EE_comp.emplace_back(list_h_PU200_prompt_2sigma_EE.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_2sigma_EE.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_3sigma_EB_comp.emplace_back(list_h_PU200_prompt_3sigma_EB.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_3sigma_EB.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_3sigma_EE_comp.emplace_back(list_h_PU200_prompt_3sigma_EE.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_3sigma_EE.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_4sigma_EB_comp.emplace_back(list_h_PU200_prompt_4sigma_EB.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_4sigma_EB.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_4sigma_EE_comp.emplace_back(list_h_PU200_prompt_4sigma_EE.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_4sigma_EE.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_genMatched_EB_comp.emplace_back(list_h_PU200_prompt_genMatched_EB.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_genMatched_EB.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_genMatched_EE_comp.emplace_back(list_h_PU200_prompt_genMatched_EE.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_genMatched_EE.at(idz)->Integral(1,nbin+1));
	    // nonprompt
	  nonprompt_eff_PU200_EB_comp.emplace_back(list_h_PU200_nonprompt_EB.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_EB.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_EE_comp.emplace_back(list_h_PU200_nonprompt_EE.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_EE.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_1sigma_EB_comp.emplace_back(list_h_PU200_nonprompt_1sigma_EB.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_1sigma_EB.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_1sigma_EE_comp.emplace_back(list_h_PU200_nonprompt_1sigma_EE.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_1sigma_EE.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_2sigma_EB_comp.emplace_back(list_h_PU200_nonprompt_2sigma_EB.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_2sigma_EB.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_2sigma_EE_comp.emplace_back(list_h_PU200_nonprompt_2sigma_EE.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_2sigma_EE.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_3sigma_EB_comp.emplace_back(list_h_PU200_nonprompt_3sigma_EB.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_3sigma_EB.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_3sigma_EE_comp.emplace_back(list_h_PU200_nonprompt_3sigma_EE.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_3sigma_EE.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_4sigma_EB_comp.emplace_back(list_h_PU200_nonprompt_4sigma_EB.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_4sigma_EB.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_4sigma_EE_comp.emplace_back(list_h_PU200_nonprompt_4sigma_EE.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_4sigma_EE.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_genMatched_EB_comp.emplace_back(list_h_PU200_nonprompt_genMatched_EB.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_genMatched_EB.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_genMatched_EE_comp.emplace_back(list_h_PU200_nonprompt_genMatched_EE.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_genMatched_EE.at(idz)->Integral(1,nbin+1));
	  // noPU
	    // prompt
	  prompt_eff_noPU_EB_comp.emplace_back(list_h_noPU_prompt_EB.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_EB.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_EE_comp.emplace_back(list_h_noPU_prompt_EE.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_EE.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_1sigma_EB_comp.emplace_back(list_h_noPU_prompt_1sigma_EB.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_1sigma_EB.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_1sigma_EE_comp.emplace_back(list_h_noPU_prompt_1sigma_EE.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_1sigma_EE.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_2sigma_EB_comp.emplace_back(list_h_noPU_prompt_2sigma_EB.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_2sigma_EB.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_2sigma_EE_comp.emplace_back(list_h_noPU_prompt_2sigma_EE.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_2sigma_EE.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_3sigma_EB_comp.emplace_back(list_h_noPU_prompt_3sigma_EB.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_3sigma_EB.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_3sigma_EE_comp.emplace_back(list_h_noPU_prompt_3sigma_EE.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_3sigma_EE.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_4sigma_EB_comp.emplace_back(list_h_noPU_prompt_4sigma_EB.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_4sigma_EB.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_4sigma_EE_comp.emplace_back(list_h_noPU_prompt_4sigma_EE.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_4sigma_EE.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_genMatched_EB_comp.emplace_back(list_h_noPU_prompt_genMatched_EB.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_genMatched_EB.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_genMatched_EE_comp.emplace_back(list_h_noPU_prompt_genMatched_EE.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_genMatched_EE.at(idz)->Integral(1,nbin+1));
	    // nonprompt
	  nonprompt_eff_noPU_EB_comp.emplace_back(list_h_noPU_nonprompt_EB.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_EB.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_EE_comp.emplace_back(list_h_noPU_nonprompt_EE.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_EE.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_1sigma_EB_comp.emplace_back(list_h_noPU_nonprompt_1sigma_EB.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_1sigma_EB.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_1sigma_EE_comp.emplace_back(list_h_noPU_nonprompt_1sigma_EE.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_1sigma_EE.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_2sigma_EB_comp.emplace_back(list_h_noPU_nonprompt_2sigma_EB.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_2sigma_EB.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_2sigma_EE_comp.emplace_back(list_h_noPU_nonprompt_2sigma_EE.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_2sigma_EE.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_3sigma_EB_comp.emplace_back(list_h_noPU_nonprompt_3sigma_EB.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_3sigma_EB.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_3sigma_EE_comp.emplace_back(list_h_noPU_nonprompt_3sigma_EE.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_3sigma_EE.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_4sigma_EB_comp.emplace_back(list_h_noPU_nonprompt_4sigma_EB.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_4sigma_EB.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_4sigma_EE_comp.emplace_back(list_h_noPU_nonprompt_4sigma_EE.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_4sigma_EE.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_genMatched_EB_comp.emplace_back(list_h_noPU_nonprompt_genMatched_EB.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_genMatched_EB.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_genMatched_EE_comp.emplace_back(list_h_noPU_nonprompt_genMatched_EE.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_genMatched_EE.at(idz)->Integral(1,nbin+1));

	  // (PV, track)
	  // PU200
	    // prompt
	  prompt_eff_PU200_EB_vtx_comp.emplace_back(list_h_PU200_prompt_EB_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_EB_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_EE_vtx_comp.emplace_back(list_h_PU200_prompt_EE_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_EE_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_1sigma_EB_vtx_comp.emplace_back(list_h_PU200_prompt_1sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_1sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_1sigma_EE_vtx_comp.emplace_back(list_h_PU200_prompt_1sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_1sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_2sigma_EB_vtx_comp.emplace_back(list_h_PU200_prompt_2sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_2sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_2sigma_EE_vtx_comp.emplace_back(list_h_PU200_prompt_2sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_2sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_3sigma_EB_vtx_comp.emplace_back(list_h_PU200_prompt_3sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_3sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_3sigma_EE_vtx_comp.emplace_back(list_h_PU200_prompt_3sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_3sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_4sigma_EB_vtx_comp.emplace_back(list_h_PU200_prompt_4sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_4sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_4sigma_EE_vtx_comp.emplace_back(list_h_PU200_prompt_4sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_4sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_genMatched_EB_vtx_comp.emplace_back(list_h_PU200_prompt_genMatched_EB_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_genMatched_EB_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_PU200_genMatched_EE_vtx_comp.emplace_back(list_h_PU200_prompt_genMatched_EE_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_prompt_genMatched_EE_vtx.at(idz)->Integral(1,nbin+1));
	    // nonprompt
	  nonprompt_eff_PU200_EB_vtx_comp.emplace_back(list_h_PU200_nonprompt_EB_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_EB_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_EE_vtx_comp.emplace_back(list_h_PU200_nonprompt_EE_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_EE_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_1sigma_EB_vtx_comp.emplace_back(list_h_PU200_nonprompt_1sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_1sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_1sigma_EE_vtx_comp.emplace_back(list_h_PU200_nonprompt_1sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_1sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_2sigma_EB_vtx_comp.emplace_back(list_h_PU200_nonprompt_2sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_2sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_2sigma_EE_vtx_comp.emplace_back(list_h_PU200_nonprompt_2sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_2sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_3sigma_EB_vtx_comp.emplace_back(list_h_PU200_nonprompt_3sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_3sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_3sigma_EE_vtx_comp.emplace_back(list_h_PU200_nonprompt_3sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_3sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_4sigma_EB_vtx_comp.emplace_back(list_h_PU200_nonprompt_4sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_4sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_4sigma_EE_vtx_comp.emplace_back(list_h_PU200_nonprompt_4sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_4sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_genMatched_EB_vtx_comp.emplace_back(list_h_PU200_nonprompt_genMatched_EB_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_genMatched_EB_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_PU200_genMatched_EE_vtx_comp.emplace_back(list_h_PU200_nonprompt_genMatched_EE_vtx.at(idz)->Integral(1,i+1)/list_h_PU200_nonprompt_genMatched_EE_vtx.at(idz)->Integral(1,nbin+1));
	  // noPU
	    // prompt
	  prompt_eff_noPU_EB_vtx_comp.emplace_back(list_h_noPU_prompt_EB_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_EB_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_EE_vtx_comp.emplace_back(list_h_noPU_prompt_EE_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_EE_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_1sigma_EB_vtx_comp.emplace_back(list_h_noPU_prompt_1sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_1sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_1sigma_EE_vtx_comp.emplace_back(list_h_noPU_prompt_1sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_1sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_2sigma_EB_vtx_comp.emplace_back(list_h_noPU_prompt_2sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_2sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_2sigma_EE_vtx_comp.emplace_back(list_h_noPU_prompt_2sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_2sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_3sigma_EB_vtx_comp.emplace_back(list_h_noPU_prompt_3sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_3sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_3sigma_EE_vtx_comp.emplace_back(list_h_noPU_prompt_3sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_3sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_4sigma_EB_vtx_comp.emplace_back(list_h_noPU_prompt_4sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_4sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_4sigma_EE_vtx_comp.emplace_back(list_h_noPU_prompt_4sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_4sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_genMatched_EB_vtx_comp.emplace_back(list_h_noPU_prompt_genMatched_EB_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_genMatched_EB_vtx.at(idz)->Integral(1,nbin+1));
	  prompt_eff_noPU_genMatched_EE_vtx_comp.emplace_back(list_h_noPU_prompt_genMatched_EE_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_prompt_genMatched_EE_vtx.at(idz)->Integral(1,nbin+1));
	    // nonprompt
	  nonprompt_eff_noPU_EB_vtx_comp.emplace_back(list_h_noPU_nonprompt_EB_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_EB_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_EE_vtx_comp.emplace_back(list_h_noPU_nonprompt_EE_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_EE_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_1sigma_EB_vtx_comp.emplace_back(list_h_noPU_nonprompt_1sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_1sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_1sigma_EE_vtx_comp.emplace_back(list_h_noPU_nonprompt_1sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_1sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_2sigma_EB_vtx_comp.emplace_back(list_h_noPU_nonprompt_2sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_2sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_2sigma_EE_vtx_comp.emplace_back(list_h_noPU_nonprompt_2sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_2sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_3sigma_EB_vtx_comp.emplace_back(list_h_noPU_nonprompt_3sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_3sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_3sigma_EE_vtx_comp.emplace_back(list_h_noPU_nonprompt_3sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_3sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_4sigma_EB_vtx_comp.emplace_back(list_h_noPU_nonprompt_4sigma_EB_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_4sigma_EB_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_4sigma_EE_vtx_comp.emplace_back(list_h_noPU_nonprompt_4sigma_EE_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_4sigma_EE_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_genMatched_EB_vtx_comp.emplace_back(list_h_noPU_nonprompt_genMatched_EB_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_genMatched_EB_vtx.at(idz)->Integral(1,nbin+1));
	  nonprompt_eff_noPU_genMatched_EE_vtx_comp.emplace_back(list_h_noPU_nonprompt_genMatched_EE_vtx.at(idz)->Integral(1,i+1)/list_h_noPU_nonprompt_genMatched_EE_vtx.at(idz)->Integral(1,nbin+1));

	}

	// (muon, track)
	prompt_eff_PU200_EB.emplace_back(prompt_eff_PU200_EB_comp); prompt_eff_PU200_EE.emplace_back(prompt_eff_PU200_EE_comp);
	prompt_eff_PU200_1sigma_EB.emplace_back(prompt_eff_PU200_1sigma_EB_comp); prompt_eff_PU200_1sigma_EE.emplace_back(prompt_eff_PU200_1sigma_EE_comp);
	prompt_eff_PU200_2sigma_EB.emplace_back(prompt_eff_PU200_2sigma_EB_comp); prompt_eff_PU200_2sigma_EE.emplace_back(prompt_eff_PU200_2sigma_EE_comp);
	prompt_eff_PU200_3sigma_EB.emplace_back(prompt_eff_PU200_3sigma_EB_comp); prompt_eff_PU200_3sigma_EE.emplace_back(prompt_eff_PU200_3sigma_EE_comp);
	prompt_eff_PU200_4sigma_EB.emplace_back(prompt_eff_PU200_4sigma_EB_comp); prompt_eff_PU200_4sigma_EE.emplace_back(prompt_eff_PU200_4sigma_EE_comp);
	prompt_eff_PU200_genMatched_EB.emplace_back(prompt_eff_PU200_genMatched_EB_comp); prompt_eff_PU200_genMatched_EE.emplace_back(prompt_eff_PU200_genMatched_EE_comp);
	nonprompt_eff_PU200_EB.emplace_back(nonprompt_eff_PU200_EB_comp); nonprompt_eff_PU200_EE.emplace_back(nonprompt_eff_PU200_EE_comp);
	nonprompt_eff_PU200_1sigma_EB.emplace_back(nonprompt_eff_PU200_1sigma_EB_comp); nonprompt_eff_PU200_1sigma_EE.emplace_back(nonprompt_eff_PU200_1sigma_EE_comp);
	nonprompt_eff_PU200_2sigma_EB.emplace_back(nonprompt_eff_PU200_2sigma_EB_comp); nonprompt_eff_PU200_2sigma_EE.emplace_back(nonprompt_eff_PU200_2sigma_EE_comp);
	nonprompt_eff_PU200_3sigma_EB.emplace_back(nonprompt_eff_PU200_3sigma_EB_comp); nonprompt_eff_PU200_3sigma_EE.emplace_back(nonprompt_eff_PU200_3sigma_EE_comp);
	nonprompt_eff_PU200_4sigma_EB.emplace_back(nonprompt_eff_PU200_4sigma_EB_comp); nonprompt_eff_PU200_4sigma_EE.emplace_back(nonprompt_eff_PU200_4sigma_EE_comp);
	nonprompt_eff_PU200_genMatched_EB.emplace_back(nonprompt_eff_PU200_genMatched_EB_comp); nonprompt_eff_PU200_genMatched_EE.emplace_back(nonprompt_eff_PU200_genMatched_EE_comp);
	prompt_eff_noPU_EB.emplace_back(prompt_eff_noPU_EB_comp); prompt_eff_noPU_EE.emplace_back(prompt_eff_noPU_EE_comp);
	prompt_eff_noPU_1sigma_EB.emplace_back(prompt_eff_noPU_1sigma_EB_comp); prompt_eff_noPU_1sigma_EE.emplace_back(prompt_eff_noPU_1sigma_EE_comp);
	prompt_eff_noPU_2sigma_EB.emplace_back(prompt_eff_noPU_2sigma_EB_comp); prompt_eff_noPU_2sigma_EE.emplace_back(prompt_eff_noPU_2sigma_EE_comp);
	prompt_eff_noPU_3sigma_EB.emplace_back(prompt_eff_noPU_3sigma_EB_comp); prompt_eff_noPU_3sigma_EE.emplace_back(prompt_eff_noPU_3sigma_EE_comp);
	prompt_eff_noPU_4sigma_EB.emplace_back(prompt_eff_noPU_4sigma_EB_comp); prompt_eff_noPU_4sigma_EE.emplace_back(prompt_eff_noPU_4sigma_EE_comp);
	prompt_eff_noPU_genMatched_EB.emplace_back(prompt_eff_noPU_genMatched_EB_comp); prompt_eff_noPU_genMatched_EE.emplace_back(prompt_eff_noPU_genMatched_EE_comp);
	nonprompt_eff_noPU_EB.emplace_back(nonprompt_eff_noPU_EB_comp); nonprompt_eff_noPU_EE.emplace_back(nonprompt_eff_noPU_EE_comp);
	nonprompt_eff_noPU_1sigma_EB.emplace_back(nonprompt_eff_noPU_1sigma_EB_comp); nonprompt_eff_noPU_1sigma_EE.emplace_back(nonprompt_eff_noPU_1sigma_EE_comp);
	nonprompt_eff_noPU_2sigma_EB.emplace_back(nonprompt_eff_noPU_2sigma_EB_comp); nonprompt_eff_noPU_2sigma_EE.emplace_back(nonprompt_eff_noPU_2sigma_EE_comp);
	nonprompt_eff_noPU_3sigma_EB.emplace_back(nonprompt_eff_noPU_3sigma_EB_comp); nonprompt_eff_noPU_3sigma_EE.emplace_back(nonprompt_eff_noPU_3sigma_EE_comp);
	nonprompt_eff_noPU_4sigma_EB.emplace_back(nonprompt_eff_noPU_4sigma_EB_comp); nonprompt_eff_noPU_4sigma_EE.emplace_back(nonprompt_eff_noPU_4sigma_EE_comp);
	nonprompt_eff_noPU_genMatched_EB.emplace_back(nonprompt_eff_noPU_genMatched_EB_comp); nonprompt_eff_noPU_genMatched_EE.emplace_back(nonprompt_eff_noPU_genMatched_EE_comp);
	// (muon, track)
	prompt_eff_PU200_EB_vtx.emplace_back(prompt_eff_PU200_EB_vtx_comp); prompt_eff_PU200_EE_vtx.emplace_back(prompt_eff_PU200_EE_vtx_comp);
	prompt_eff_PU200_1sigma_EB_vtx.emplace_back(prompt_eff_PU200_1sigma_EB_vtx_comp); prompt_eff_PU200_1sigma_EE_vtx.emplace_back(prompt_eff_PU200_1sigma_EE_vtx_comp);
	prompt_eff_PU200_2sigma_EB_vtx.emplace_back(prompt_eff_PU200_2sigma_EB_vtx_comp); prompt_eff_PU200_2sigma_EE_vtx.emplace_back(prompt_eff_PU200_2sigma_EE_vtx_comp);
	prompt_eff_PU200_3sigma_EB_vtx.emplace_back(prompt_eff_PU200_3sigma_EB_vtx_comp); prompt_eff_PU200_3sigma_EE_vtx.emplace_back(prompt_eff_PU200_3sigma_EE_vtx_comp);
	prompt_eff_PU200_4sigma_EB_vtx.emplace_back(prompt_eff_PU200_4sigma_EB_vtx_comp); prompt_eff_PU200_4sigma_EE_vtx.emplace_back(prompt_eff_PU200_4sigma_EE_vtx_comp);
	prompt_eff_PU200_genMatched_EB_vtx.emplace_back(prompt_eff_PU200_genMatched_EB_vtx_comp); prompt_eff_PU200_genMatched_EE_vtx.emplace_back(prompt_eff_PU200_genMatched_EE_vtx_comp);
	nonprompt_eff_PU200_EB_vtx.emplace_back(nonprompt_eff_PU200_EB_vtx_comp); nonprompt_eff_PU200_EE_vtx.emplace_back(nonprompt_eff_PU200_EE_vtx_comp);
	nonprompt_eff_PU200_1sigma_EB_vtx.emplace_back(nonprompt_eff_PU200_1sigma_EB_vtx_comp); nonprompt_eff_PU200_1sigma_EE_vtx.emplace_back(nonprompt_eff_PU200_1sigma_EE_vtx_comp);
	nonprompt_eff_PU200_2sigma_EB_vtx.emplace_back(nonprompt_eff_PU200_2sigma_EB_vtx_comp); nonprompt_eff_PU200_2sigma_EE_vtx.emplace_back(nonprompt_eff_PU200_2sigma_EE_vtx_comp);
	nonprompt_eff_PU200_3sigma_EB_vtx.emplace_back(nonprompt_eff_PU200_3sigma_EB_vtx_comp); nonprompt_eff_PU200_3sigma_EE_vtx.emplace_back(nonprompt_eff_PU200_3sigma_EE_vtx_comp);
	nonprompt_eff_PU200_4sigma_EB_vtx.emplace_back(nonprompt_eff_PU200_4sigma_EB_vtx_comp); nonprompt_eff_PU200_4sigma_EE_vtx.emplace_back(nonprompt_eff_PU200_4sigma_EE_vtx_comp);
	nonprompt_eff_PU200_genMatched_EB_vtx.emplace_back(nonprompt_eff_PU200_genMatched_EB_vtx_comp); nonprompt_eff_PU200_genMatched_EE_vtx.emplace_back(nonprompt_eff_PU200_genMatched_EE_vtx_comp);
	prompt_eff_noPU_EB_vtx.emplace_back(prompt_eff_noPU_EB_vtx_comp); prompt_eff_noPU_EE_vtx.emplace_back(prompt_eff_noPU_EE_vtx_comp);
	prompt_eff_noPU_1sigma_EB_vtx.emplace_back(prompt_eff_noPU_1sigma_EB_vtx_comp); prompt_eff_noPU_1sigma_EE_vtx.emplace_back(prompt_eff_noPU_1sigma_EE_vtx_comp);
	prompt_eff_noPU_2sigma_EB_vtx.emplace_back(prompt_eff_noPU_2sigma_EB_vtx_comp); prompt_eff_noPU_2sigma_EE_vtx.emplace_back(prompt_eff_noPU_2sigma_EE_vtx_comp);
	prompt_eff_noPU_3sigma_EB_vtx.emplace_back(prompt_eff_noPU_3sigma_EB_vtx_comp); prompt_eff_noPU_3sigma_EE_vtx.emplace_back(prompt_eff_noPU_3sigma_EE_vtx_comp);
	prompt_eff_noPU_4sigma_EB_vtx.emplace_back(prompt_eff_noPU_4sigma_EB_vtx_comp); prompt_eff_noPU_4sigma_EE_vtx.emplace_back(prompt_eff_noPU_4sigma_EE_vtx_comp);
	prompt_eff_noPU_genMatched_EB_vtx.emplace_back(prompt_eff_noPU_genMatched_EB_vtx_comp); prompt_eff_noPU_genMatched_EE_vtx.emplace_back(prompt_eff_noPU_genMatched_EE_vtx_comp);
	nonprompt_eff_noPU_EB_vtx.emplace_back(nonprompt_eff_noPU_EB_vtx_comp); nonprompt_eff_noPU_EE_vtx.emplace_back(nonprompt_eff_noPU_EE_vtx_comp);
	nonprompt_eff_noPU_1sigma_EB_vtx.emplace_back(nonprompt_eff_noPU_1sigma_EB_vtx_comp); nonprompt_eff_noPU_1sigma_EE_vtx.emplace_back(nonprompt_eff_noPU_1sigma_EE_vtx_comp);
	nonprompt_eff_noPU_2sigma_EB_vtx.emplace_back(nonprompt_eff_noPU_2sigma_EB_vtx_comp); nonprompt_eff_noPU_2sigma_EE_vtx.emplace_back(nonprompt_eff_noPU_2sigma_EE_vtx_comp);
	nonprompt_eff_noPU_3sigma_EB_vtx.emplace_back(nonprompt_eff_noPU_3sigma_EB_vtx_comp); nonprompt_eff_noPU_3sigma_EE_vtx.emplace_back(nonprompt_eff_noPU_3sigma_EE_vtx_comp);
	nonprompt_eff_noPU_4sigma_EB_vtx.emplace_back(nonprompt_eff_noPU_4sigma_EB_vtx_comp); nonprompt_eff_noPU_4sigma_EE_vtx.emplace_back(nonprompt_eff_noPU_4sigma_EE_vtx_comp);
	nonprompt_eff_noPU_genMatched_EB_vtx.emplace_back(nonprompt_eff_noPU_genMatched_EB_vtx_comp); nonprompt_eff_noPU_genMatched_EE_vtx.emplace_back(nonprompt_eff_noPU_genMatched_EE_vtx_comp);

	// Initialize
	// (muon, track)
	prompt_eff_PU200_EB_comp={}, prompt_eff_PU200_EE_comp={}, nonprompt_eff_PU200_EB_comp={}, nonprompt_eff_PU200_EE_comp={};
	prompt_eff_PU200_1sigma_EB_comp={}, prompt_eff_PU200_1sigma_EE_comp={}, nonprompt_eff_PU200_1sigma_EB_comp={}, nonprompt_eff_PU200_1sigma_EE_comp={};
	prompt_eff_PU200_2sigma_EB_comp={}, prompt_eff_PU200_2sigma_EE_comp={}, nonprompt_eff_PU200_2sigma_EB_comp={}, nonprompt_eff_PU200_2sigma_EE_comp={};
	prompt_eff_PU200_3sigma_EB_comp={}, prompt_eff_PU200_3sigma_EE_comp={}, nonprompt_eff_PU200_3sigma_EB_comp={}, nonprompt_eff_PU200_3sigma_EE_comp={};
	prompt_eff_PU200_4sigma_EB_comp={}, prompt_eff_PU200_4sigma_EE_comp={}, nonprompt_eff_PU200_4sigma_EB_comp={}, nonprompt_eff_PU200_4sigma_EE_comp={};
	prompt_eff_PU200_genMatched_EB_comp={}, prompt_eff_PU200_genMatched_EE_comp={}, nonprompt_eff_PU200_genMatched_EB_comp={}, nonprompt_eff_PU200_genMatched_EE_comp={};
	prompt_eff_noPU_EB_comp={}, prompt_eff_noPU_EE_comp={}, nonprompt_eff_noPU_EB_comp={}, nonprompt_eff_noPU_EE_comp={};
	prompt_eff_noPU_1sigma_EB_comp={}, prompt_eff_noPU_1sigma_EE_comp={}, nonprompt_eff_noPU_1sigma_EB_comp={}, nonprompt_eff_noPU_1sigma_EE_comp={};
	prompt_eff_noPU_2sigma_EB_comp={}, prompt_eff_noPU_2sigma_EE_comp={}, nonprompt_eff_noPU_2sigma_EB_comp={}, nonprompt_eff_noPU_2sigma_EE_comp={};
	prompt_eff_noPU_3sigma_EB_comp={}, prompt_eff_noPU_3sigma_EE_comp={}, nonprompt_eff_noPU_3sigma_EB_comp={}, nonprompt_eff_noPU_3sigma_EE_comp={};
	prompt_eff_noPU_4sigma_EB_comp={}, prompt_eff_noPU_4sigma_EE_comp={}, nonprompt_eff_noPU_4sigma_EB_comp={}, nonprompt_eff_noPU_4sigma_EE_comp={};
	prompt_eff_noPU_genMatched_EB_comp={}, prompt_eff_noPU_genMatched_EE_comp={}, nonprompt_eff_noPU_genMatched_EB_comp={}, nonprompt_eff_noPU_genMatched_EE_comp={};
	// (muon, track)
	prompt_eff_PU200_EB_vtx_comp={}, prompt_eff_PU200_EE_vtx_comp={}, nonprompt_eff_PU200_EB_vtx_comp={}, nonprompt_eff_PU200_EE_vtx_comp={};
	prompt_eff_PU200_1sigma_EB_vtx_comp={}, prompt_eff_PU200_1sigma_EE_vtx_comp={}, nonprompt_eff_PU200_1sigma_EB_vtx_comp={}, nonprompt_eff_PU200_1sigma_EE_vtx_comp={};
	prompt_eff_PU200_2sigma_EB_vtx_comp={}, prompt_eff_PU200_2sigma_EE_vtx_comp={}, nonprompt_eff_PU200_2sigma_EB_vtx_comp={}, nonprompt_eff_PU200_2sigma_EE_vtx_comp={};
	prompt_eff_PU200_3sigma_EB_vtx_comp={}, prompt_eff_PU200_3sigma_EE_vtx_comp={}, nonprompt_eff_PU200_3sigma_EB_vtx_comp={}, nonprompt_eff_PU200_3sigma_EE_vtx_comp={};
	prompt_eff_PU200_4sigma_EB_vtx_comp={}, prompt_eff_PU200_4sigma_EE_vtx_comp={}, nonprompt_eff_PU200_4sigma_EB_vtx_comp={}, nonprompt_eff_PU200_4sigma_EE_vtx_comp={};
	prompt_eff_PU200_genMatched_EB_vtx_comp={}, prompt_eff_PU200_genMatched_EE_vtx_comp={}, nonprompt_eff_PU200_genMatched_EB_vtx_comp={}, nonprompt_eff_PU200_genMatched_EE_vtx_comp={};
	prompt_eff_noPU_EB_vtx_comp={}, prompt_eff_noPU_EE_vtx_comp={}, nonprompt_eff_noPU_EB_vtx_comp={}, nonprompt_eff_noPU_EE_vtx_comp={};
	prompt_eff_noPU_1sigma_EB_vtx_comp={}, prompt_eff_noPU_1sigma_EE_vtx_comp={}, nonprompt_eff_noPU_1sigma_EB_vtx_comp={}, nonprompt_eff_noPU_1sigma_EE_vtx_comp={};
	prompt_eff_noPU_2sigma_EB_vtx_comp={}, prompt_eff_noPU_2sigma_EE_vtx_comp={}, nonprompt_eff_noPU_2sigma_EB_vtx_comp={}, nonprompt_eff_noPU_2sigma_EE_vtx_comp={};
	prompt_eff_noPU_3sigma_EB_vtx_comp={}, prompt_eff_noPU_3sigma_EE_vtx_comp={}, nonprompt_eff_noPU_3sigma_EB_vtx_comp={}, nonprompt_eff_noPU_3sigma_EE_vtx_comp={};
	prompt_eff_noPU_4sigma_EB_vtx_comp={}, prompt_eff_noPU_4sigma_EE_vtx_comp={}, nonprompt_eff_noPU_4sigma_EB_vtx_comp={}, nonprompt_eff_noPU_4sigma_EE_vtx_comp={};
	prompt_eff_noPU_genMatched_EB_vtx_comp={}, prompt_eff_noPU_genMatched_EE_vtx_comp={}, nonprompt_eff_noPU_genMatched_EB_vtx_comp={}, nonprompt_eff_noPU_genMatched_EE_vtx_comp={};

  }


  // Define TGraph
  // (muon, track)
  vector<TGraph*> list_gr_eff_PU200_prompt_EB,            list_gr_eff_PU200_prompt_EE;
  vector<TGraph*> list_gr_eff_PU200_prompt_1sigma_EB,     list_gr_eff_PU200_prompt_1sigma_EE;
  vector<TGraph*> list_gr_eff_PU200_prompt_2sigma_EB,     list_gr_eff_PU200_prompt_2sigma_EE;
  vector<TGraph*> list_gr_eff_PU200_prompt_3sigma_EB,     list_gr_eff_PU200_prompt_3sigma_EE;
  vector<TGraph*> list_gr_eff_PU200_prompt_4sigma_EB,     list_gr_eff_PU200_prompt_4sigma_EE;
  vector<TGraph*> list_gr_eff_PU200_prompt_genMatched_EB,     list_gr_eff_PU200_prompt_genMatched_EE;
  vector<TGraph*> list_gr_eff_PU200_nonprompt_EB,         list_gr_eff_PU200_nonprompt_EE;
  vector<TGraph*> list_gr_eff_PU200_nonprompt_1sigma_EB,  list_gr_eff_PU200_nonprompt_1sigma_EE;
  vector<TGraph*> list_gr_eff_PU200_nonprompt_2sigma_EB,  list_gr_eff_PU200_nonprompt_2sigma_EE;
  vector<TGraph*> list_gr_eff_PU200_nonprompt_3sigma_EB,  list_gr_eff_PU200_nonprompt_3sigma_EE;
  vector<TGraph*> list_gr_eff_PU200_nonprompt_4sigma_EB,  list_gr_eff_PU200_nonprompt_4sigma_EE;
  vector<TGraph*> list_gr_eff_PU200_nonprompt_genMatched_EB,  list_gr_eff_PU200_nonprompt_genMatched_EE;
  vector<TGraph*> list_gr_eff_noPU_prompt_EB,             list_gr_eff_noPU_prompt_EE;
  vector<TGraph*> list_gr_eff_noPU_prompt_1sigma_EB,      list_gr_eff_noPU_prompt_1sigma_EE;
  vector<TGraph*> list_gr_eff_noPU_prompt_2sigma_EB,      list_gr_eff_noPU_prompt_2sigma_EE;
  vector<TGraph*> list_gr_eff_noPU_prompt_3sigma_EB,      list_gr_eff_noPU_prompt_3sigma_EE;
  vector<TGraph*> list_gr_eff_noPU_prompt_4sigma_EB,      list_gr_eff_noPU_prompt_4sigma_EE;
  vector<TGraph*> list_gr_eff_noPU_prompt_genMatched_EB,      list_gr_eff_noPU_prompt_genMatched_EE;
  vector<TGraph*> list_gr_eff_noPU_nonprompt_EB,          list_gr_eff_noPU_nonprompt_EE;
  vector<TGraph*> list_gr_eff_noPU_nonprompt_1sigma_EB,   list_gr_eff_noPU_nonprompt_1sigma_EE;
  vector<TGraph*> list_gr_eff_noPU_nonprompt_2sigma_EB,   list_gr_eff_noPU_nonprompt_2sigma_EE;
  vector<TGraph*> list_gr_eff_noPU_nonprompt_3sigma_EB,   list_gr_eff_noPU_nonprompt_3sigma_EE;
  vector<TGraph*> list_gr_eff_noPU_nonprompt_4sigma_EB,   list_gr_eff_noPU_nonprompt_4sigma_EE;
  vector<TGraph*> list_gr_eff_noPU_nonprompt_genMatched_EB,   list_gr_eff_noPU_nonprompt_genMatched_EE;
  // (PV, track)
  vector<TGraph*> list_gr_eff_PU200_prompt_1sigma_EB_vtx,     list_gr_eff_PU200_prompt_1sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_PU200_prompt_2sigma_EB_vtx,     list_gr_eff_PU200_prompt_2sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_PU200_prompt_3sigma_EB_vtx,     list_gr_eff_PU200_prompt_3sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_PU200_prompt_4sigma_EB_vtx,     list_gr_eff_PU200_prompt_4sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_PU200_nonprompt_1sigma_EB_vtx,  list_gr_eff_PU200_nonprompt_1sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_PU200_nonprompt_2sigma_EB_vtx,  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_PU200_nonprompt_3sigma_EB_vtx,  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_PU200_nonprompt_4sigma_EB_vtx,  list_gr_eff_PU200_nonprompt_4sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_noPU_prompt_1sigma_EB_vtx,      list_gr_eff_noPU_prompt_1sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_noPU_prompt_2sigma_EB_vtx,      list_gr_eff_noPU_prompt_2sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_noPU_prompt_3sigma_EB_vtx,      list_gr_eff_noPU_prompt_3sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_noPU_prompt_4sigma_EB_vtx,      list_gr_eff_noPU_prompt_4sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_noPU_nonprompt_1sigma_EB_vtx,   list_gr_eff_noPU_nonprompt_1sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_noPU_nonprompt_2sigma_EB_vtx,   list_gr_eff_noPU_nonprompt_2sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_noPU_nonprompt_3sigma_EB_vtx,   list_gr_eff_noPU_nonprompt_3sigma_EE_vtx;
  vector<TGraph*> list_gr_eff_noPU_nonprompt_4sigma_EB_vtx,   list_gr_eff_noPU_nonprompt_4sigma_EE_vtx;


  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
	// (muon, track)
    TGraph* gr_eff_PU200_prompt_EB = new TGraph();             TGraph* gr_eff_PU200_prompt_EE = new TGraph();
    TGraph* gr_eff_PU200_prompt_1sigma_EB = new TGraph();      TGraph* gr_eff_PU200_prompt_1sigma_EE = new TGraph();
    TGraph* gr_eff_PU200_prompt_2sigma_EB = new TGraph();      TGraph* gr_eff_PU200_prompt_2sigma_EE = new TGraph();
    TGraph* gr_eff_PU200_prompt_3sigma_EB = new TGraph();      TGraph* gr_eff_PU200_prompt_3sigma_EE = new TGraph();
    TGraph* gr_eff_PU200_prompt_4sigma_EB = new TGraph();      TGraph* gr_eff_PU200_prompt_4sigma_EE = new TGraph();
    TGraph* gr_eff_PU200_prompt_genMatched_EB = new TGraph();      TGraph* gr_eff_PU200_prompt_genMatched_EE = new TGraph();
    TGraph* gr_eff_PU200_nonprompt_EB = new TGraph();          TGraph* gr_eff_PU200_nonprompt_EE = new TGraph();
    TGraph* gr_eff_PU200_nonprompt_1sigma_EB = new TGraph();   TGraph* gr_eff_PU200_nonprompt_1sigma_EE = new TGraph();
    TGraph* gr_eff_PU200_nonprompt_2sigma_EB = new TGraph();   TGraph* gr_eff_PU200_nonprompt_2sigma_EE = new TGraph();
    TGraph* gr_eff_PU200_nonprompt_3sigma_EB = new TGraph();   TGraph* gr_eff_PU200_nonprompt_3sigma_EE = new TGraph();
    TGraph* gr_eff_PU200_nonprompt_4sigma_EB = new TGraph();   TGraph* gr_eff_PU200_nonprompt_4sigma_EE = new TGraph();
    TGraph* gr_eff_PU200_nonprompt_genMatched_EB = new TGraph();   TGraph* gr_eff_PU200_nonprompt_genMatched_EE = new TGraph();
    TGraph* gr_eff_noPU_prompt_EB = new TGraph();              TGraph* gr_eff_noPU_prompt_EE = new TGraph();
    TGraph* gr_eff_noPU_prompt_1sigma_EB = new TGraph();       TGraph* gr_eff_noPU_prompt_1sigma_EE = new TGraph();
    TGraph* gr_eff_noPU_prompt_2sigma_EB = new TGraph();       TGraph* gr_eff_noPU_prompt_2sigma_EE = new TGraph();
    TGraph* gr_eff_noPU_prompt_3sigma_EB = new TGraph();       TGraph* gr_eff_noPU_prompt_3sigma_EE = new TGraph();
    TGraph* gr_eff_noPU_prompt_4sigma_EB = new TGraph();       TGraph* gr_eff_noPU_prompt_4sigma_EE = new TGraph();
    TGraph* gr_eff_noPU_prompt_genMatched_EB = new TGraph();       TGraph* gr_eff_noPU_prompt_genMatched_EE = new TGraph();
    TGraph* gr_eff_noPU_nonprompt_EB = new TGraph();           TGraph* gr_eff_noPU_nonprompt_EE = new TGraph();
    TGraph* gr_eff_noPU_nonprompt_1sigma_EB = new TGraph();    TGraph* gr_eff_noPU_nonprompt_1sigma_EE = new TGraph();
    TGraph* gr_eff_noPU_nonprompt_2sigma_EB = new TGraph();    TGraph* gr_eff_noPU_nonprompt_2sigma_EE = new TGraph();
    TGraph* gr_eff_noPU_nonprompt_3sigma_EB = new TGraph();    TGraph* gr_eff_noPU_nonprompt_3sigma_EE = new TGraph();
    TGraph* gr_eff_noPU_nonprompt_4sigma_EB = new TGraph();    TGraph* gr_eff_noPU_nonprompt_4sigma_EE = new TGraph();
    TGraph* gr_eff_noPU_nonprompt_genMatched_EB = new TGraph();    TGraph* gr_eff_noPU_nonprompt_genMatched_EE = new TGraph();
	// (PV, track)
    TGraph* gr_eff_PU200_prompt_1sigma_EB_vtx = new TGraph();      TGraph* gr_eff_PU200_prompt_1sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_PU200_prompt_2sigma_EB_vtx = new TGraph();      TGraph* gr_eff_PU200_prompt_2sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_PU200_prompt_3sigma_EB_vtx = new TGraph();      TGraph* gr_eff_PU200_prompt_3sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_PU200_prompt_4sigma_EB_vtx = new TGraph();      TGraph* gr_eff_PU200_prompt_4sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_PU200_nonprompt_1sigma_EB_vtx = new TGraph();   TGraph* gr_eff_PU200_nonprompt_1sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_PU200_nonprompt_2sigma_EB_vtx = new TGraph();   TGraph* gr_eff_PU200_nonprompt_2sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_PU200_nonprompt_3sigma_EB_vtx = new TGraph();   TGraph* gr_eff_PU200_nonprompt_3sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_PU200_nonprompt_4sigma_EB_vtx = new TGraph();   TGraph* gr_eff_PU200_nonprompt_4sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_noPU_prompt_1sigma_EB_vtx = new TGraph();       TGraph* gr_eff_noPU_prompt_1sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_noPU_prompt_2sigma_EB_vtx = new TGraph();       TGraph* gr_eff_noPU_prompt_2sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_noPU_prompt_3sigma_EB_vtx = new TGraph();       TGraph* gr_eff_noPU_prompt_3sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_noPU_prompt_4sigma_EB_vtx = new TGraph();       TGraph* gr_eff_noPU_prompt_4sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_noPU_nonprompt_1sigma_EB_vtx = new TGraph();    TGraph* gr_eff_noPU_nonprompt_1sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_noPU_nonprompt_2sigma_EB_vtx = new TGraph();    TGraph* gr_eff_noPU_nonprompt_2sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_noPU_nonprompt_3sigma_EB_vtx = new TGraph();    TGraph* gr_eff_noPU_nonprompt_3sigma_EE_vtx = new TGraph();
    TGraph* gr_eff_noPU_nonprompt_4sigma_EB_vtx = new TGraph();    TGraph* gr_eff_noPU_nonprompt_4sigma_EE_vtx = new TGraph();

    for(unsigned int i=0; i<nbin+1; i++) {
	  // (muon, track)
	  // PU200
	    // prompt
	  gr_eff_PU200_prompt_EB->SetPoint(gr_eff_PU200_prompt_EB->GetN(), 0.004*i, prompt_eff_PU200_EB.at(idz).at(i));
	  gr_eff_PU200_prompt_EE->SetPoint(gr_eff_PU200_prompt_EE->GetN(), 0.004*i, prompt_eff_PU200_EE.at(idz).at(i));
	  gr_eff_PU200_prompt_1sigma_EB->SetPoint(gr_eff_PU200_prompt_1sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_1sigma_EB.at(idz).at(i));
	  gr_eff_PU200_prompt_1sigma_EE->SetPoint(gr_eff_PU200_prompt_1sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_1sigma_EE.at(idz).at(i));
	  gr_eff_PU200_prompt_2sigma_EB->SetPoint(gr_eff_PU200_prompt_2sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EB.at(idz).at(i));
	  gr_eff_PU200_prompt_2sigma_EE->SetPoint(gr_eff_PU200_prompt_2sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EE.at(idz).at(i));
	  gr_eff_PU200_prompt_3sigma_EB->SetPoint(gr_eff_PU200_prompt_3sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EB.at(idz).at(i));
	  gr_eff_PU200_prompt_3sigma_EE->SetPoint(gr_eff_PU200_prompt_3sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EE.at(idz).at(i));
	  gr_eff_PU200_prompt_4sigma_EB->SetPoint(gr_eff_PU200_prompt_4sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EB.at(idz).at(i));
	  gr_eff_PU200_prompt_4sigma_EE->SetPoint(gr_eff_PU200_prompt_4sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EE.at(idz).at(i));
	  gr_eff_PU200_prompt_genMatched_EB->SetPoint(gr_eff_PU200_prompt_genMatched_EB->GetN(), 0.004*i, prompt_eff_PU200_genMatched_EB.at(idz).at(i));
	  gr_eff_PU200_prompt_genMatched_EE->SetPoint(gr_eff_PU200_prompt_genMatched_EE->GetN(), 0.004*i, prompt_eff_PU200_genMatched_EE.at(idz).at(i));
        // nonprompt
	  gr_eff_PU200_nonprompt_EB->SetPoint(gr_eff_PU200_nonprompt_EB->GetN(), 0.004*i, nonprompt_eff_PU200_EB.at(idz).at(i));
	  gr_eff_PU200_nonprompt_EE->SetPoint(gr_eff_PU200_nonprompt_EE->GetN(), 0.004*i, nonprompt_eff_PU200_EE.at(idz).at(i));
	  gr_eff_PU200_nonprompt_1sigma_EB->SetPoint(gr_eff_PU200_nonprompt_1sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_1sigma_EB.at(idz).at(i));
	  gr_eff_PU200_nonprompt_1sigma_EE->SetPoint(gr_eff_PU200_nonprompt_1sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_1sigma_EE.at(idz).at(i));
	  gr_eff_PU200_nonprompt_2sigma_EB->SetPoint(gr_eff_PU200_nonprompt_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EB.at(idz).at(i));
	  gr_eff_PU200_nonprompt_2sigma_EE->SetPoint(gr_eff_PU200_nonprompt_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EE.at(idz).at(i));
	  gr_eff_PU200_nonprompt_3sigma_EB->SetPoint(gr_eff_PU200_nonprompt_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EB.at(idz).at(i));
	  gr_eff_PU200_nonprompt_3sigma_EE->SetPoint(gr_eff_PU200_nonprompt_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EE.at(idz).at(i));
	  gr_eff_PU200_nonprompt_4sigma_EB->SetPoint(gr_eff_PU200_nonprompt_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EB.at(idz).at(i));
	  gr_eff_PU200_nonprompt_4sigma_EE->SetPoint(gr_eff_PU200_nonprompt_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EE.at(idz).at(i));
	  gr_eff_PU200_nonprompt_genMatched_EB->SetPoint(gr_eff_PU200_nonprompt_genMatched_EB->GetN(), 0.004*i, nonprompt_eff_PU200_genMatched_EB.at(idz).at(i));
	  gr_eff_PU200_nonprompt_genMatched_EE->SetPoint(gr_eff_PU200_nonprompt_genMatched_EE->GetN(), 0.004*i, nonprompt_eff_PU200_genMatched_EE.at(idz).at(i));
	  // noPU
	    // prompt
	  gr_eff_noPU_prompt_EB->SetPoint(gr_eff_noPU_prompt_EB->GetN(), 0.004*i, prompt_eff_noPU_EB.at(idz).at(i));
	  gr_eff_noPU_prompt_EE->SetPoint(gr_eff_noPU_prompt_EE->GetN(), 0.004*i, prompt_eff_noPU_EE.at(idz).at(i));
	  gr_eff_noPU_prompt_1sigma_EB->SetPoint(gr_eff_noPU_prompt_1sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_1sigma_EB.at(idz).at(i));
	  gr_eff_noPU_prompt_1sigma_EE->SetPoint(gr_eff_noPU_prompt_1sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_1sigma_EE.at(idz).at(i));
	  gr_eff_noPU_prompt_2sigma_EB->SetPoint(gr_eff_noPU_prompt_2sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_2sigma_EB.at(idz).at(i));
	  gr_eff_noPU_prompt_2sigma_EE->SetPoint(gr_eff_noPU_prompt_2sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_2sigma_EE.at(idz).at(i));
	  gr_eff_noPU_prompt_3sigma_EB->SetPoint(gr_eff_noPU_prompt_3sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_3sigma_EB.at(idz).at(i));
	  gr_eff_noPU_prompt_3sigma_EE->SetPoint(gr_eff_noPU_prompt_3sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_3sigma_EE.at(idz).at(i));
	  gr_eff_noPU_prompt_4sigma_EB->SetPoint(gr_eff_noPU_prompt_4sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_4sigma_EB.at(idz).at(i));
	  gr_eff_noPU_prompt_4sigma_EE->SetPoint(gr_eff_noPU_prompt_4sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_4sigma_EE.at(idz).at(i));
	  gr_eff_noPU_prompt_genMatched_EB->SetPoint(gr_eff_noPU_prompt_genMatched_EB->GetN(), 0.004*i, prompt_eff_noPU_genMatched_EB.at(idz).at(i));
	  gr_eff_noPU_prompt_genMatched_EE->SetPoint(gr_eff_noPU_prompt_genMatched_EE->GetN(), 0.004*i, prompt_eff_noPU_genMatched_EE.at(idz).at(i));
        // nonprompt
	  gr_eff_noPU_nonprompt_EB->SetPoint(gr_eff_noPU_nonprompt_EB->GetN(), 0.004*i, nonprompt_eff_noPU_EB.at(idz).at(i));
	  gr_eff_noPU_nonprompt_EE->SetPoint(gr_eff_noPU_nonprompt_EE->GetN(), 0.004*i, nonprompt_eff_noPU_EE.at(idz).at(i));
	  gr_eff_noPU_nonprompt_1sigma_EB->SetPoint(gr_eff_noPU_nonprompt_1sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_1sigma_EB.at(idz).at(i));
	  gr_eff_noPU_nonprompt_1sigma_EE->SetPoint(gr_eff_noPU_nonprompt_1sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_1sigma_EE.at(idz).at(i));
	  gr_eff_noPU_nonprompt_2sigma_EB->SetPoint(gr_eff_noPU_nonprompt_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_2sigma_EB.at(idz).at(i));
	  gr_eff_noPU_nonprompt_2sigma_EE->SetPoint(gr_eff_noPU_nonprompt_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_2sigma_EE.at(idz).at(i));
	  gr_eff_noPU_nonprompt_3sigma_EB->SetPoint(gr_eff_noPU_nonprompt_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_3sigma_EB.at(idz).at(i));
	  gr_eff_noPU_nonprompt_3sigma_EE->SetPoint(gr_eff_noPU_nonprompt_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_3sigma_EE.at(idz).at(i));
	  gr_eff_noPU_nonprompt_4sigma_EB->SetPoint(gr_eff_noPU_nonprompt_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_4sigma_EB.at(idz).at(i));
	  gr_eff_noPU_nonprompt_4sigma_EE->SetPoint(gr_eff_noPU_nonprompt_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_4sigma_EE.at(idz).at(i));
	  gr_eff_noPU_nonprompt_genMatched_EB->SetPoint(gr_eff_noPU_nonprompt_genMatched_EB->GetN(), 0.004*i, nonprompt_eff_noPU_genMatched_EB.at(idz).at(i));
	  gr_eff_noPU_nonprompt_genMatched_EE->SetPoint(gr_eff_noPU_nonprompt_genMatched_EE->GetN(), 0.004*i, nonprompt_eff_noPU_genMatched_EE.at(idz).at(i));

	  // (PV, track)
	  // PU200
	    // prompt
	  gr_eff_PU200_prompt_1sigma_EB_vtx->SetPoint(gr_eff_PU200_prompt_1sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_PU200_1sigma_EB_vtx.at(idz).at(i));
	  gr_eff_PU200_prompt_1sigma_EE_vtx->SetPoint(gr_eff_PU200_prompt_1sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_PU200_1sigma_EE_vtx.at(idz).at(i));
	  gr_eff_PU200_prompt_2sigma_EB_vtx->SetPoint(gr_eff_PU200_prompt_2sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EB_vtx.at(idz).at(i));
	  gr_eff_PU200_prompt_2sigma_EE_vtx->SetPoint(gr_eff_PU200_prompt_2sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EE_vtx.at(idz).at(i));
	  gr_eff_PU200_prompt_3sigma_EB_vtx->SetPoint(gr_eff_PU200_prompt_3sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EB_vtx.at(idz).at(i));
	  gr_eff_PU200_prompt_3sigma_EE_vtx->SetPoint(gr_eff_PU200_prompt_3sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EE_vtx.at(idz).at(i));
	  gr_eff_PU200_prompt_4sigma_EB_vtx->SetPoint(gr_eff_PU200_prompt_4sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EB_vtx.at(idz).at(i));
	  gr_eff_PU200_prompt_4sigma_EE_vtx->SetPoint(gr_eff_PU200_prompt_4sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EE_vtx.at(idz).at(i));
        // nonprompt
	  gr_eff_PU200_nonprompt_1sigma_EB_vtx->SetPoint(gr_eff_PU200_nonprompt_1sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_1sigma_EB_vtx.at(idz).at(i));
	  gr_eff_PU200_nonprompt_1sigma_EE_vtx->SetPoint(gr_eff_PU200_nonprompt_1sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_1sigma_EE_vtx.at(idz).at(i));
	  gr_eff_PU200_nonprompt_2sigma_EB_vtx->SetPoint(gr_eff_PU200_nonprompt_2sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EB_vtx.at(idz).at(i));
	  gr_eff_PU200_nonprompt_2sigma_EE_vtx->SetPoint(gr_eff_PU200_nonprompt_2sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EE_vtx.at(idz).at(i));
	  gr_eff_PU200_nonprompt_3sigma_EB_vtx->SetPoint(gr_eff_PU200_nonprompt_3sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EB_vtx.at(idz).at(i));
	  gr_eff_PU200_nonprompt_3sigma_EE_vtx->SetPoint(gr_eff_PU200_nonprompt_3sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EE_vtx.at(idz).at(i));
	  gr_eff_PU200_nonprompt_4sigma_EB_vtx->SetPoint(gr_eff_PU200_nonprompt_4sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EB_vtx.at(idz).at(i));
	  gr_eff_PU200_nonprompt_4sigma_EE_vtx->SetPoint(gr_eff_PU200_nonprompt_4sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EE_vtx.at(idz).at(i));
	  // noPU
	    // prompt
	  gr_eff_noPU_prompt_1sigma_EB_vtx->SetPoint(gr_eff_noPU_prompt_1sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_noPU_1sigma_EB_vtx.at(idz).at(i));
	  gr_eff_noPU_prompt_1sigma_EE_vtx->SetPoint(gr_eff_noPU_prompt_1sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_noPU_1sigma_EE_vtx.at(idz).at(i));
	  gr_eff_noPU_prompt_2sigma_EB_vtx->SetPoint(gr_eff_noPU_prompt_2sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_noPU_2sigma_EB_vtx.at(idz).at(i));
	  gr_eff_noPU_prompt_2sigma_EE_vtx->SetPoint(gr_eff_noPU_prompt_2sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_noPU_2sigma_EE_vtx.at(idz).at(i));
	  gr_eff_noPU_prompt_3sigma_EB_vtx->SetPoint(gr_eff_noPU_prompt_3sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_noPU_3sigma_EB_vtx.at(idz).at(i));
	  gr_eff_noPU_prompt_3sigma_EE_vtx->SetPoint(gr_eff_noPU_prompt_3sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_noPU_3sigma_EE_vtx.at(idz).at(i));
	  gr_eff_noPU_prompt_4sigma_EB_vtx->SetPoint(gr_eff_noPU_prompt_4sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_noPU_4sigma_EB_vtx.at(idz).at(i));
	  gr_eff_noPU_prompt_4sigma_EE_vtx->SetPoint(gr_eff_noPU_prompt_4sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_noPU_4sigma_EE_vtx.at(idz).at(i));
        // nonprompt
	  gr_eff_noPU_nonprompt_1sigma_EB_vtx->SetPoint(gr_eff_noPU_nonprompt_1sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_noPU_1sigma_EB_vtx.at(idz).at(i));
	  gr_eff_noPU_nonprompt_1sigma_EE_vtx->SetPoint(gr_eff_noPU_nonprompt_1sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_noPU_1sigma_EE_vtx.at(idz).at(i));
	  gr_eff_noPU_nonprompt_2sigma_EB_vtx->SetPoint(gr_eff_noPU_nonprompt_2sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_noPU_2sigma_EB_vtx.at(idz).at(i));
	  gr_eff_noPU_nonprompt_2sigma_EE_vtx->SetPoint(gr_eff_noPU_nonprompt_2sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_noPU_2sigma_EE_vtx.at(idz).at(i));
	  gr_eff_noPU_nonprompt_3sigma_EB_vtx->SetPoint(gr_eff_noPU_nonprompt_3sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_noPU_3sigma_EB_vtx.at(idz).at(i));
	  gr_eff_noPU_nonprompt_3sigma_EE_vtx->SetPoint(gr_eff_noPU_nonprompt_3sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_noPU_3sigma_EE_vtx.at(idz).at(i));
	  gr_eff_noPU_nonprompt_4sigma_EB_vtx->SetPoint(gr_eff_noPU_nonprompt_4sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_noPU_4sigma_EB_vtx.at(idz).at(i));
	  gr_eff_noPU_nonprompt_4sigma_EE_vtx->SetPoint(gr_eff_noPU_nonprompt_4sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_noPU_4sigma_EE_vtx.at(idz).at(i));

    }

	// (muon, track)
	list_gr_eff_PU200_prompt_EB.emplace_back(gr_eff_PU200_prompt_EB);
	list_gr_eff_PU200_prompt_EE.emplace_back(gr_eff_PU200_prompt_EE);
	list_gr_eff_PU200_prompt_1sigma_EB.emplace_back(gr_eff_PU200_prompt_1sigma_EB);
	list_gr_eff_PU200_prompt_1sigma_EE.emplace_back(gr_eff_PU200_prompt_1sigma_EE);
	list_gr_eff_PU200_prompt_2sigma_EB.emplace_back(gr_eff_PU200_prompt_2sigma_EB);
	list_gr_eff_PU200_prompt_2sigma_EE.emplace_back(gr_eff_PU200_prompt_2sigma_EE);
	list_gr_eff_PU200_prompt_3sigma_EB.emplace_back(gr_eff_PU200_prompt_3sigma_EB);
	list_gr_eff_PU200_prompt_3sigma_EE.emplace_back(gr_eff_PU200_prompt_3sigma_EE);
	list_gr_eff_PU200_prompt_4sigma_EB.emplace_back(gr_eff_PU200_prompt_4sigma_EB);
	list_gr_eff_PU200_prompt_4sigma_EE.emplace_back(gr_eff_PU200_prompt_4sigma_EE);
	list_gr_eff_PU200_prompt_genMatched_EB.emplace_back(gr_eff_PU200_prompt_genMatched_EB);
	list_gr_eff_PU200_prompt_genMatched_EE.emplace_back(gr_eff_PU200_prompt_genMatched_EE);
	list_gr_eff_PU200_nonprompt_EB.emplace_back(gr_eff_PU200_nonprompt_EB);
	list_gr_eff_PU200_nonprompt_EE.emplace_back(gr_eff_PU200_nonprompt_EE);
	list_gr_eff_PU200_nonprompt_1sigma_EB.emplace_back(gr_eff_PU200_nonprompt_1sigma_EB);
	list_gr_eff_PU200_nonprompt_1sigma_EE.emplace_back(gr_eff_PU200_nonprompt_1sigma_EE);
	list_gr_eff_PU200_nonprompt_2sigma_EB.emplace_back(gr_eff_PU200_nonprompt_2sigma_EB);
	list_gr_eff_PU200_nonprompt_2sigma_EE.emplace_back(gr_eff_PU200_nonprompt_2sigma_EE);
	list_gr_eff_PU200_nonprompt_3sigma_EB.emplace_back(gr_eff_PU200_nonprompt_3sigma_EB);
	list_gr_eff_PU200_nonprompt_3sigma_EE.emplace_back(gr_eff_PU200_nonprompt_3sigma_EE);
	list_gr_eff_PU200_nonprompt_4sigma_EB.emplace_back(gr_eff_PU200_nonprompt_4sigma_EB);
	list_gr_eff_PU200_nonprompt_4sigma_EE.emplace_back(gr_eff_PU200_nonprompt_4sigma_EE);
	list_gr_eff_PU200_nonprompt_genMatched_EB.emplace_back(gr_eff_PU200_nonprompt_genMatched_EB);
	list_gr_eff_PU200_nonprompt_genMatched_EE.emplace_back(gr_eff_PU200_nonprompt_genMatched_EE);
	list_gr_eff_noPU_prompt_EB.emplace_back(gr_eff_noPU_prompt_EB);
	list_gr_eff_noPU_prompt_EE.emplace_back(gr_eff_noPU_prompt_EE);
	list_gr_eff_noPU_prompt_1sigma_EB.emplace_back(gr_eff_noPU_prompt_1sigma_EB);
	list_gr_eff_noPU_prompt_1sigma_EE.emplace_back(gr_eff_noPU_prompt_1sigma_EE);
	list_gr_eff_noPU_prompt_2sigma_EB.emplace_back(gr_eff_noPU_prompt_2sigma_EB);
	list_gr_eff_noPU_prompt_2sigma_EE.emplace_back(gr_eff_noPU_prompt_2sigma_EE);
	list_gr_eff_noPU_prompt_3sigma_EB.emplace_back(gr_eff_noPU_prompt_3sigma_EB);
	list_gr_eff_noPU_prompt_3sigma_EE.emplace_back(gr_eff_noPU_prompt_3sigma_EE);
	list_gr_eff_noPU_prompt_4sigma_EB.emplace_back(gr_eff_noPU_prompt_4sigma_EB);
	list_gr_eff_noPU_prompt_4sigma_EE.emplace_back(gr_eff_noPU_prompt_4sigma_EE);
	list_gr_eff_noPU_prompt_genMatched_EB.emplace_back(gr_eff_noPU_prompt_genMatched_EB);
	list_gr_eff_noPU_prompt_genMatched_EE.emplace_back(gr_eff_noPU_prompt_genMatched_EE);
	list_gr_eff_noPU_nonprompt_EB.emplace_back(gr_eff_noPU_nonprompt_EB);
	list_gr_eff_noPU_nonprompt_EE.emplace_back(gr_eff_noPU_nonprompt_EE);
	list_gr_eff_noPU_nonprompt_1sigma_EB.emplace_back(gr_eff_noPU_nonprompt_1sigma_EB);
	list_gr_eff_noPU_nonprompt_1sigma_EE.emplace_back(gr_eff_noPU_nonprompt_1sigma_EE);
	list_gr_eff_noPU_nonprompt_2sigma_EB.emplace_back(gr_eff_noPU_nonprompt_2sigma_EB);
	list_gr_eff_noPU_nonprompt_2sigma_EE.emplace_back(gr_eff_noPU_nonprompt_2sigma_EE);
	list_gr_eff_noPU_nonprompt_3sigma_EB.emplace_back(gr_eff_noPU_nonprompt_3sigma_EB);
	list_gr_eff_noPU_nonprompt_3sigma_EE.emplace_back(gr_eff_noPU_nonprompt_3sigma_EE);
	list_gr_eff_noPU_nonprompt_4sigma_EB.emplace_back(gr_eff_noPU_nonprompt_4sigma_EB);
	list_gr_eff_noPU_nonprompt_4sigma_EE.emplace_back(gr_eff_noPU_nonprompt_4sigma_EE);
	list_gr_eff_noPU_nonprompt_genMatched_EB.emplace_back(gr_eff_noPU_nonprompt_genMatched_EB);
	list_gr_eff_noPU_nonprompt_genMatched_EE.emplace_back(gr_eff_noPU_nonprompt_genMatched_EE);

	// (PV, track)
	list_gr_eff_PU200_prompt_1sigma_EB_vtx.emplace_back(gr_eff_PU200_prompt_1sigma_EB_vtx);
	list_gr_eff_PU200_prompt_1sigma_EE_vtx.emplace_back(gr_eff_PU200_prompt_1sigma_EE_vtx);
	list_gr_eff_PU200_prompt_2sigma_EB_vtx.emplace_back(gr_eff_PU200_prompt_2sigma_EB_vtx);
	list_gr_eff_PU200_prompt_2sigma_EE_vtx.emplace_back(gr_eff_PU200_prompt_2sigma_EE_vtx);
	list_gr_eff_PU200_prompt_3sigma_EB_vtx.emplace_back(gr_eff_PU200_prompt_3sigma_EB_vtx);
	list_gr_eff_PU200_prompt_3sigma_EE_vtx.emplace_back(gr_eff_PU200_prompt_3sigma_EE_vtx);
	list_gr_eff_PU200_prompt_4sigma_EB_vtx.emplace_back(gr_eff_PU200_prompt_4sigma_EB_vtx);
	list_gr_eff_PU200_prompt_4sigma_EE_vtx.emplace_back(gr_eff_PU200_prompt_4sigma_EE_vtx);
	list_gr_eff_PU200_nonprompt_1sigma_EB_vtx.emplace_back(gr_eff_PU200_nonprompt_1sigma_EB_vtx);
	list_gr_eff_PU200_nonprompt_1sigma_EE_vtx.emplace_back(gr_eff_PU200_nonprompt_1sigma_EE_vtx);
	list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.emplace_back(gr_eff_PU200_nonprompt_2sigma_EB_vtx);
	list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.emplace_back(gr_eff_PU200_nonprompt_2sigma_EE_vtx);
	list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.emplace_back(gr_eff_PU200_nonprompt_3sigma_EB_vtx);
	list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.emplace_back(gr_eff_PU200_nonprompt_3sigma_EE_vtx);
	list_gr_eff_PU200_nonprompt_4sigma_EB_vtx.emplace_back(gr_eff_PU200_nonprompt_4sigma_EB_vtx);
	list_gr_eff_PU200_nonprompt_4sigma_EE_vtx.emplace_back(gr_eff_PU200_nonprompt_4sigma_EE_vtx);
	list_gr_eff_noPU_prompt_1sigma_EB_vtx.emplace_back(gr_eff_noPU_prompt_1sigma_EB_vtx);
	list_gr_eff_noPU_prompt_1sigma_EE_vtx.emplace_back(gr_eff_noPU_prompt_1sigma_EE_vtx);
	list_gr_eff_noPU_prompt_2sigma_EB_vtx.emplace_back(gr_eff_noPU_prompt_2sigma_EB_vtx);
	list_gr_eff_noPU_prompt_2sigma_EE_vtx.emplace_back(gr_eff_noPU_prompt_2sigma_EE_vtx);
	list_gr_eff_noPU_prompt_3sigma_EB_vtx.emplace_back(gr_eff_noPU_prompt_3sigma_EB_vtx);
	list_gr_eff_noPU_prompt_3sigma_EE_vtx.emplace_back(gr_eff_noPU_prompt_3sigma_EE_vtx);
	list_gr_eff_noPU_prompt_4sigma_EB_vtx.emplace_back(gr_eff_noPU_prompt_4sigma_EB_vtx);
	list_gr_eff_noPU_prompt_4sigma_EE_vtx.emplace_back(gr_eff_noPU_prompt_4sigma_EE_vtx);
	list_gr_eff_noPU_nonprompt_1sigma_EB_vtx.emplace_back(gr_eff_noPU_nonprompt_1sigma_EB_vtx);
	list_gr_eff_noPU_nonprompt_1sigma_EE_vtx.emplace_back(gr_eff_noPU_nonprompt_1sigma_EE_vtx);
	list_gr_eff_noPU_nonprompt_2sigma_EB_vtx.emplace_back(gr_eff_noPU_nonprompt_2sigma_EB_vtx);
	list_gr_eff_noPU_nonprompt_2sigma_EE_vtx.emplace_back(gr_eff_noPU_nonprompt_2sigma_EE_vtx);
	list_gr_eff_noPU_nonprompt_3sigma_EB_vtx.emplace_back(gr_eff_noPU_nonprompt_3sigma_EB_vtx);
	list_gr_eff_noPU_nonprompt_3sigma_EE_vtx.emplace_back(gr_eff_noPU_nonprompt_3sigma_EE_vtx);
	list_gr_eff_noPU_nonprompt_4sigma_EB_vtx.emplace_back(gr_eff_noPU_nonprompt_4sigma_EB_vtx);
	list_gr_eff_noPU_nonprompt_4sigma_EE_vtx.emplace_back(gr_eff_noPU_nonprompt_4sigma_EE_vtx);


  }

  ///////////////
  // Cosmetics //
  ///////////////
  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
	// (muon, track)
	list_gr_eff_PU200_prompt_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_1sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_1sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_1sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_1sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_2sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_2sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_2sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_2sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_3sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_3sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_3sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_3sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_4sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_4sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_4sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_4sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_genMatched_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_genMatched_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_genMatched_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_genMatched_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_nonprompt_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_nonprompt_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_nonprompt_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_nonprompt_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_nonprompt_1sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_nonprompt_1sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_nonprompt_1sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_nonprompt_1sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_nonprompt_2sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_nonprompt_2sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_nonprompt_2sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_nonprompt_2sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_nonprompt_3sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_nonprompt_3sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_nonprompt_3sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_nonprompt_3sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_nonprompt_4sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_nonprompt_4sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_nonprompt_4sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_nonprompt_4sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_nonprompt_genMatched_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_nonprompt_genMatched_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_nonprompt_genMatched_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_nonprompt_genMatched_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_prompt_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_prompt_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_prompt_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_prompt_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_prompt_1sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_prompt_1sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_prompt_1sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_prompt_1sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_prompt_2sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_prompt_2sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_prompt_2sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_prompt_2sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_prompt_3sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_prompt_3sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_prompt_3sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_prompt_3sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_prompt_4sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_prompt_4sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_prompt_4sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_prompt_4sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_prompt_genMatched_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_prompt_genMatched_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_prompt_genMatched_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_prompt_genMatched_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_nonprompt_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_nonprompt_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_nonprompt_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_nonprompt_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_nonprompt_1sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_nonprompt_1sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_nonprompt_1sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_nonprompt_1sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_nonprompt_2sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_nonprompt_2sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_nonprompt_2sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_nonprompt_2sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_nonprompt_3sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_nonprompt_3sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_nonprompt_3sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_nonprompt_3sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_nonprompt_4sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_nonprompt_4sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_nonprompt_4sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_nonprompt_4sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_nonprompt_genMatched_EB.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_nonprompt_genMatched_EB.at(idz)->SetLineWidth(2);
	list_gr_eff_noPU_nonprompt_genMatched_EE.at(idz)->SetLineColor(idz+1); list_gr_eff_noPU_nonprompt_genMatched_EE.at(idz)->SetLineWidth(2);

	list_gr_eff_PU200_prompt_1sigma_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_1sigma_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_1sigma_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_1sigma_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_4sigma_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_4sigma_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_eff_PU200_prompt_4sigma_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_eff_PU200_prompt_4sigma_EE_vtx.at(idz)->SetLineWidth(2);
	if(idz==4) {
	  list_gr_eff_PU200_prompt_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_1sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_1sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_2sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_2sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_3sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_3sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_4sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_4sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_genMatched_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_genMatched_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_nonprompt_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_nonprompt_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_nonprompt_1sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_nonprompt_1sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_nonprompt_2sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_nonprompt_2sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_nonprompt_3sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_nonprompt_3sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_nonprompt_4sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_nonprompt_4sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_nonprompt_genMatched_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_nonprompt_genMatched_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_prompt_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_prompt_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_prompt_1sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_prompt_1sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_prompt_2sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_prompt_2sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_prompt_3sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_prompt_3sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_prompt_4sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_prompt_4sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_prompt_genMatched_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_prompt_genMatched_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_nonprompt_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_nonprompt_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_nonprompt_1sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_nonprompt_1sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_nonprompt_2sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_nonprompt_2sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_nonprompt_3sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_nonprompt_3sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_nonprompt_4sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_nonprompt_4sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_nonprompt_genMatched_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_noPU_nonprompt_genMatched_EE.at(idz)->SetLineColor(kMagenta);

	  list_gr_eff_PU200_prompt_1sigma_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_1sigma_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_4sigma_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_eff_PU200_prompt_4sigma_EE_vtx.at(idz)->SetLineColor(kMagenta);
	}
  }
  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
    list_gr_eff_PU200_prompt_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_1sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_1sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_1sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_1sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_2sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_2sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_2sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_2sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_3sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_3sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_3sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_3sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_4sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_4sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_4sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_4sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_genMatched_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_genMatched_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_genMatched_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_genMatched_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_nonprompt_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_nonprompt_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_nonprompt_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_nonprompt_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_nonprompt_1sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_nonprompt_1sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_nonprompt_1sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_nonprompt_1sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_nonprompt_2sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_nonprompt_2sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_nonprompt_2sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_nonprompt_2sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_nonprompt_3sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_nonprompt_3sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_nonprompt_3sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_nonprompt_3sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_nonprompt_4sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_nonprompt_4sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_nonprompt_4sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_nonprompt_4sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_nonprompt_genMatched_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_nonprompt_genMatched_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_nonprompt_genMatched_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_nonprompt_genMatched_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_prompt_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_prompt_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_prompt_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_prompt_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_prompt_1sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_prompt_1sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_prompt_1sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_prompt_1sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_prompt_2sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_prompt_2sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_prompt_2sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_prompt_2sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_prompt_3sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_prompt_3sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_prompt_3sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_prompt_3sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_prompt_4sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_prompt_4sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_prompt_4sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_prompt_4sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_prompt_genMatched_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_prompt_genMatched_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_prompt_genMatched_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_prompt_genMatched_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_nonprompt_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_nonprompt_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_nonprompt_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_nonprompt_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_nonprompt_1sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_nonprompt_1sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_nonprompt_1sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_nonprompt_1sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_nonprompt_2sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_nonprompt_2sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_nonprompt_2sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_nonprompt_2sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_nonprompt_3sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_nonprompt_3sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_nonprompt_3sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_nonprompt_3sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_nonprompt_4sigma_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_nonprompt_4sigma_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_nonprompt_4sigma_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_nonprompt_4sigma_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_nonprompt_genMatched_EB.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_nonprompt_genMatched_EB.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_noPU_nonprompt_genMatched_EE.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_noPU_nonprompt_genMatched_EE.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");


    list_gr_eff_PU200_prompt_1sigma_EB_vtx.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_1sigma_EB_vtx.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_1sigma_EE_vtx.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_1sigma_EE_vtx.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_4sigma_EB_vtx.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_4sigma_EB_vtx.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
    list_gr_eff_PU200_prompt_4sigma_EE_vtx.at(idz)->GetXaxis()->SetTitle("Relative isolation"); list_gr_eff_PU200_prompt_4sigma_EE_vtx.at(idz)->GetYaxis()->SetTitle("Isolation efficiency");
  }


  /////////////
  // Legends //
  /////////////
  TLegend* leg_eff_PU200_prompt_EB = new TLegend(0.40, 0.13, 0.88, 0.40); TLegend* leg_eff_PU200_prompt_EE = new TLegend(0.40, 0.13, 0.88, 0.40);
  TLegend* leg_eff_PU200_nonprompt_EB = new TLegend(0.40, 0.13, 0.88, 0.40); TLegend* leg_eff_PU200_nonprompt_EE = new TLegend(0.40, 0.13, 0.88, 0.40);
  TLegend* leg_eff_noPU_prompt_EB = new TLegend(0.40, 0.13, 0.88, 0.40); TLegend* leg_eff_noPU_prompt_EE = new TLegend(0.40, 0.13, 0.88, 0.40);
  TLegend* leg_eff_noPU_nonprompt_EB = new TLegend(0.40, 0.13, 0.88, 0.40); TLegend* leg_eff_noPU_nonprompt_EE = new TLegend(0.40, 0.13, 0.88, 0.40);
  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
	leg_eff_PU200_prompt_EB->AddEntry(list_gr_eff_PU200_prompt_EB.at(idz), Form("no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(idz)*10), "L");
	leg_eff_PU200_prompt_EB->SetTextSize(0.03);
	leg_eff_PU200_prompt_EE->AddEntry(list_gr_eff_PU200_prompt_EE.at(idz), Form("no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(idz)*10), "L");
	leg_eff_PU200_prompt_EE->SetTextSize(0.03);
	leg_eff_PU200_nonprompt_EB->AddEntry(list_gr_eff_PU200_nonprompt_EB.at(idz), Form("no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(idz)*10), "L");
	leg_eff_PU200_nonprompt_EB->SetTextSize(0.03);
	leg_eff_PU200_nonprompt_EE->AddEntry(list_gr_eff_PU200_nonprompt_EE.at(idz), Form("no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(idz)*10), "L");
	leg_eff_PU200_nonprompt_EE->SetTextSize(0.03);
	leg_eff_noPU_prompt_EB->AddEntry(list_gr_eff_noPU_prompt_EB.at(idz), Form("no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(idz)*10), "L");
	leg_eff_noPU_prompt_EB->SetTextSize(0.03);
	leg_eff_noPU_prompt_EE->AddEntry(list_gr_eff_noPU_prompt_EE.at(idz), Form("no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(idz)*10), "L");
	leg_eff_noPU_prompt_EE->SetTextSize(0.03);
	leg_eff_noPU_nonprompt_EB->AddEntry(list_gr_eff_noPU_nonprompt_EB.at(idz), Form("no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(idz)*10), "L");
	leg_eff_noPU_nonprompt_EB->SetTextSize(0.03);
	leg_eff_noPU_nonprompt_EE->AddEntry(list_gr_eff_noPU_nonprompt_EE.at(idz), Form("no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(idz)*10), "L");
	leg_eff_noPU_nonprompt_EE->SetTextSize(0.03);
  }

  /////////////////
  ////// eff //////
  /////////////////

  ////////////////
  // Draw plots //
  ////////////////

  // test
  for(int i=0; i<track_pv_dz_cut_EB.size(); i++) {
	cout << "dz = " << track_pv_dz_cut_EB.at(i)*10 << " mm" << endl;
    cout << "PU200 prompt EB:    " << num_gen_notPV_PU200_prompt_EB.at(i)    << "/" << num_gen_PU200_prompt_EB.at(i) << endl;
    cout << "PU200 prompt EE:    " << num_gen_notPV_PU200_prompt_EE.at(i)    << "/" << num_gen_PU200_prompt_EE.at(i) << endl;
    cout << "PU200 nonprompt EB: " << num_gen_notPV_PU200_nonprompt_EB.at(i) << "/" << num_gen_PU200_nonprompt_EB.at(i) << endl;
    cout << "PU200 nonprompt EE: " << num_gen_notPV_PU200_nonprompt_EE.at(i) << "/" << num_gen_PU200_nonprompt_EE.at(i) << endl;
    cout << "noPU prompt EB:     " << num_gen_notPV_noPU_prompt_EB.at(i)     << "/" << num_gen_noPU_prompt_EB.at(i) << endl;
    cout << "noPU prompt EE:     " << num_gen_notPV_noPU_prompt_EE.at(i)     << "/" << num_gen_noPU_prompt_EE.at(i) << endl;
    cout << "noPU nonprompt EB:  " << num_gen_notPV_noPU_nonprompt_EB.at(i)  << "/" << num_gen_noPU_nonprompt_EB.at(i) << endl;
    cout << "noPU nonprompt EE:  " << num_gen_notPV_noPU_nonprompt_EE.at(i)  << "/" << num_gen_noPU_nonprompt_EE.at(i) << endl;
	cout << endl;
  }

  // Iso distribution
  TCanvas* c_prompt_EB_reliso_dist = new TCanvas("c_prompt_EB_reliso_dist", "c_prompt_EB_reliso_dist", 1200, 1200);
  c_prompt_EB_reliso_dist->cd();
  c_prompt_EB_reliso_dist->SetGrid();
  c_prompt_EB_reliso_dist->SetLeftMargin(0.12);
  list_h_noPU_prompt_EB.at(1)->SetTitle("");
  list_h_noPU_prompt_EB.at(1)->SetLineColor(kBlack);
  list_h_PU200_prompt_genMatched_EB.at(1)->SetLineColor(kRed);
  list_h_noPU_prompt_EB.at(1)->GetXaxis()->SetRangeUser(0.,0.1);
//  list_h_PU200_prompt_EB.at(1)->GetYaxis()->SetRangeUser(0.,100.);
  list_h_noPU_prompt_EB.at(1)->Scale(1./list_h_noPU_prompt_EB.at(1)->Integral());
  list_h_PU200_prompt_genMatched_EB.at(1)->Scale(1./list_h_PU200_prompt_genMatched_EB.at(1)->Integral());
  list_h_noPU_prompt_EB.at(1)->Draw("hist");
  list_h_PU200_prompt_genMatched_EB.at(1)->SetLineStyle(7);
  list_h_PU200_prompt_genMatched_EB.at(1)->Draw("hist same");
  TLegend* leg_prompt_EB_reliso_dist = new TLegend(0.53, 0.73, 0.88, 0.88);
  leg_prompt_EB_reliso_dist->AddEntry(list_h_noPU_prompt_EB.at(1), "noPU no MTD");
  leg_prompt_EB_reliso_dist->AddEntry(list_h_PU200_prompt_genMatched_EB.at(1), "PU200 no MTD gen");
  leg_prompt_EB_reliso_dist->Draw();
  c_prompt_EB_reliso_dist->Print("plots/ntuple/dz_study/isoeff_prompt_EB_genMatched.pdf");

  
  // test end

  
  ///////////
  // PU200 //
  ///////////
  
  // prompt
    // Barrel
  TCanvas* c_PU200_prompt_EB_reliso = new TCanvas("c_PU200_prompt_EB_reliso", "c_PU200_prompt_EB_reliso", 1200, 1200);
  c_PU200_prompt_EB_reliso->cd();
  c_PU200_prompt_EB_reliso->SetGrid();
  c_PU200_prompt_EB_reliso->SetLeftMargin(0.12);
  //list_gr_eff_PU200_prompt_EB.at(0)->GetXaxis()->SetRangeUser(0.,2.);
  //list_gr_eff_PU200_prompt_EB.at(0)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_prompt_EB.at(0)->GetXaxis()->SetRangeUser(0.,0.4);
  list_gr_eff_PU200_prompt_EB.at(0)->GetYaxis()->SetRangeUser(0.35,1.);
  list_gr_eff_PU200_prompt_EB.at(0)->SetTitle("");
  list_gr_eff_PU200_prompt_EB.at(0)->Draw("AL");
  for(int idz=1; idz<track_pv_dz_cut_EB.size(); idz++) {
	list_gr_eff_PU200_prompt_EB.at(idz)->Draw("same");
  }
//  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
//	list_gr_eff_PU200_prompt_1sigma_EB.at(idz)->Draw("same");
//  }
  leg_eff_PU200_prompt_EB->Draw();
  c_PU200_prompt_EB_reliso->Print("plots/ntuple/dz_study/isoeff_PU200_prompt_EB_dzcut.pdf");
    // Endcap
  TCanvas* c_PU200_prompt_EE_reliso = new TCanvas("c_PU200_prompt_EE_reliso", "c_PU200_prompt_EE_reliso", 1200, 1200);
  c_PU200_prompt_EE_reliso->cd();
  c_PU200_prompt_EE_reliso->SetGrid();
  c_PU200_prompt_EE_reliso->SetLeftMargin(0.12);
  //list_gr_eff_PU200_prompt_EE.at(0)->GetXaxis()->SetRangeUser(0.,2.);
  //list_gr_eff_PU200_prompt_EE.at(0)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_prompt_EE.at(0)->GetXaxis()->SetRangeUser(0.,0.4);
  list_gr_eff_PU200_prompt_EE.at(0)->GetYaxis()->SetRangeUser(0.15,1.);
  list_gr_eff_PU200_prompt_EE.at(0)->SetTitle("");
  list_gr_eff_PU200_prompt_EE.at(0)->Draw("AL");
  for(int idz=1; idz<track_pv_dz_cut_EE.size(); idz++) {
	list_gr_eff_PU200_prompt_EE.at(idz)->Draw("same");
  }
  leg_eff_PU200_prompt_EE->Draw();
  c_PU200_prompt_EE_reliso->Print("plots/ntuple/dz_study/isoeff_PU200_prompt_EE_dzcut.pdf");

  // nonprompt
    // Barrel
  TCanvas* c_PU200_nonprompt_EB_reliso = new TCanvas("c_PU200_nonprompt_EB_reliso", "c_PU200_nonprompt_EB_reliso", 1200, 1200);
  c_PU200_nonprompt_EB_reliso->cd();
  c_PU200_nonprompt_EB_reliso->SetGrid();
  c_PU200_nonprompt_EB_reliso->SetLeftMargin(0.12);
  list_gr_eff_PU200_nonprompt_EB.at(0)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_PU200_nonprompt_EB.at(0)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_nonprompt_EB.at(0)->SetTitle("");
  list_gr_eff_PU200_nonprompt_EB.at(0)->Draw("AL");
  for(int idz=1; idz<track_pv_dz_cut_EB.size(); idz++) {
	list_gr_eff_PU200_nonprompt_EB.at(idz)->Draw("same");
  }
//  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
//	list_gr_eff_PU200_nonprompt_1sigma_EB.at(idz)->Draw("same");
//  }
  leg_eff_PU200_nonprompt_EB->Draw();
  c_PU200_nonprompt_EB_reliso->Print("plots/ntuple/dz_study/isoeff_PU200_nonprompt_EB_dzcut.pdf");
    // Endcap
  TCanvas* c_PU200_nonprompt_EE_reliso = new TCanvas("c_PU200_nonprompt_EE_reliso", "c_PU200_nonprompt_EE_reliso", 1200, 1200);
  c_PU200_nonprompt_EE_reliso->cd();
  c_PU200_nonprompt_EE_reliso->SetGrid();
  c_PU200_nonprompt_EE_reliso->SetLeftMargin(0.12);
  list_gr_eff_PU200_nonprompt_EE.at(0)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_PU200_nonprompt_EE.at(0)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_nonprompt_EE.at(0)->SetTitle("");
  list_gr_eff_PU200_nonprompt_EE.at(0)->Draw("AL");
  for(int idz=1; idz<track_pv_dz_cut_EE.size(); idz++) {
	list_gr_eff_PU200_nonprompt_EE.at(idz)->Draw("same");
  }
  leg_eff_PU200_nonprompt_EE->Draw();
  c_PU200_nonprompt_EE_reliso->Print("plots/ntuple/dz_study/isoeff_PU200_nonprompt_EE_dzcut.pdf");


  ////////////
  /// noPU ///
  ////////////

  // prompt
    // Barrel
  TCanvas* c_noPU_prompt_EB_reliso = new TCanvas("c_noPU_prompt_EB_reliso", "c_noPU_prompt_EB_reliso", 1200, 1200);
  c_noPU_prompt_EB_reliso->cd();
  c_noPU_prompt_EB_reliso->SetGrid();
  c_noPU_prompt_EB_reliso->SetLeftMargin(0.12);
  list_gr_eff_noPU_prompt_EB.at(0)->GetXaxis()->SetRangeUser(0.,0.15);
  list_gr_eff_noPU_prompt_EB.at(0)->GetYaxis()->SetRangeUser(0.8,1.);
  list_gr_eff_noPU_prompt_EB.at(0)->SetTitle("");
  list_gr_eff_noPU_prompt_EB.at(0)->Draw("AL");
  for(int idz=1; idz<track_pv_dz_cut_EB.size(); idz++) {
	list_gr_eff_noPU_prompt_EB.at(idz)->Draw("same");
  }
//  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
//	list_gr_eff_noPU_prompt_1sigma_EB.at(idz)->Draw("same");
//  }
  leg_eff_noPU_prompt_EB->Draw();
  c_noPU_prompt_EB_reliso->Print("plots/ntuple/dz_study/isoeff_noPU_prompt_EB_dzcut.pdf");

    // Endcap
  TCanvas* c_noPU_prompt_EE_reliso = new TCanvas("c_noPU_prompt_EE_reliso", "c_noPU_prompt_EE_reliso", 1200, 1200);
  c_noPU_prompt_EE_reliso->cd();
  c_noPU_prompt_EE_reliso->SetGrid();
  c_noPU_prompt_EE_reliso->SetLeftMargin(0.12);
  list_gr_eff_noPU_prompt_EE.at(0)->GetXaxis()->SetRangeUser(0.,0.15);
  list_gr_eff_noPU_prompt_EE.at(0)->GetYaxis()->SetRangeUser(0.65,1.);
  list_gr_eff_noPU_prompt_EE.at(0)->SetTitle("");
  list_gr_eff_noPU_prompt_EE.at(0)->Draw("AL");
  for(int idz=1; idz<track_pv_dz_cut_EE.size(); idz++) {
	list_gr_eff_noPU_prompt_EE.at(idz)->Draw("same");
  }
  leg_eff_noPU_prompt_EE->Draw();
  c_noPU_prompt_EE_reliso->Print("plots/ntuple/dz_study/isoeff_noPU_prompt_EE_dzcut.pdf");

  // nonprompt
    // Barrel
  TCanvas* c_noPU_nonprompt_EB_reliso = new TCanvas("c_noPU_nonprompt_EB_reliso", "c_noPU_nonprompt_EB_reliso", 1200, 1200);
  c_noPU_nonprompt_EB_reliso->cd();
  c_noPU_nonprompt_EB_reliso->SetGrid();
  c_noPU_nonprompt_EB_reliso->SetLeftMargin(0.12);
  list_gr_eff_noPU_nonprompt_EB.at(0)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_noPU_nonprompt_EB.at(0)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_noPU_nonprompt_EB.at(0)->SetTitle("");
  list_gr_eff_noPU_nonprompt_EB.at(0)->Draw("AL");
  for(int idz=1; idz<track_pv_dz_cut_EB.size(); idz++) {
	list_gr_eff_noPU_nonprompt_EB.at(idz)->Draw("same");
  }
//  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
//	list_gr_eff_noPU_nonprompt_1sigma_EB.at(idz)->Draw("same");
//  }
  leg_eff_noPU_nonprompt_EB->Draw();
  c_noPU_nonprompt_EB_reliso->Print("plots/ntuple/dz_study/isoeff_noPU_nonprompt_EB_dzcut.pdf");

    // Endcap
  TCanvas* c_noPU_nonprompt_EE_reliso = new TCanvas("c_noPU_nonprompt_EE_reliso", "c_noPU_nonprompt_EE_reliso", 1200, 1200);
  c_noPU_nonprompt_EE_reliso->cd();
  c_noPU_nonprompt_EE_reliso->SetGrid();
  c_noPU_nonprompt_EE_reliso->SetLeftMargin(0.12);
  list_gr_eff_noPU_nonprompt_EE.at(0)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_noPU_nonprompt_EE.at(0)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_noPU_nonprompt_EE.at(0)->SetTitle("");
  list_gr_eff_noPU_nonprompt_EE.at(0)->Draw("AL");
  for(int idz=1; idz<track_pv_dz_cut_EE.size(); idz++) {
	list_gr_eff_noPU_nonprompt_EE.at(idz)->Draw("same");
  }
  leg_eff_noPU_nonprompt_EE->Draw();
  c_noPU_nonprompt_EE_reliso->Print("plots/ntuple/dz_study/isoeff_noPU_nonprompt_EE_dzcut.pdf");

  // test
  // dz 1
  // prompt Barrel
  TCanvas* c_prompt_EB_reliso_dz1 = new TCanvas("c_prompt_EB_reliso_dz1", "c_prompt_EB_reliso_dz1", 1200, 1200);
  TLegend* leg_eff_prompt_EB_dz1 = new TLegend(0.40, 0.13, 0.88, 0.55);
  leg_eff_prompt_EB_dz1->AddEntry(list_gr_eff_noPU_prompt_genMatched_EB.at(0), "noPU no MTD gen", "L");
  leg_eff_prompt_EB_dz1->AddEntry(list_gr_eff_noPU_prompt_EB.at(0), "noPU no MTD", "L");
  leg_eff_prompt_EB_dz1->AddEntry(list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(0), "PU200 2sigma (PV, track)", "L");
  leg_eff_prompt_EB_dz1->AddEntry(list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(0), "PU200 3sigma (PV, track)", "L");
  leg_eff_prompt_EB_dz1->AddEntry(list_gr_eff_PU200_prompt_2sigma_EB.at(0), "PU200 2sigma (muon, track)", "L");
  leg_eff_prompt_EB_dz1->AddEntry(list_gr_eff_PU200_prompt_3sigma_EB.at(0), "PU200 3sigma (muon, track)", "L");
  leg_eff_prompt_EB_dz1->AddEntry(list_gr_eff_PU200_prompt_genMatched_EB.at(0), "PU200 no MTD gen", "L");
  leg_eff_prompt_EB_dz1->AddEntry(list_gr_eff_PU200_prompt_EB.at(0), "PU200 no MTD", "L");
//  leg_eff_prompt_EB_dz1->AddEntry(list_gr_eff_noPU_prompt_2sigma_EB.at(0), "noPU 2sigma", "L");
  list_gr_eff_PU200_prompt_EB.at(0)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_2sigma_EB.at(0)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(0)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(0)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_3sigma_EB.at(0)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(0)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(0)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_genMatched_EB.at(0)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_genMatched_EB.at(0)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_EB.at(0)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EB.at(0)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EB.at(0)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_2sigma_EB.at(0)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_2sigma_EB.at(0)->SetLineStyle(10);
  c_prompt_EB_reliso_dz1->cd();
  c_prompt_EB_reliso_dz1->SetGrid();
  c_prompt_EB_reliso_dz1->SetLeftMargin(0.12);
  list_gr_eff_PU200_prompt_EB.at(0)->GetXaxis()->SetRangeUser(0.,0.4);
  list_gr_eff_PU200_prompt_EB.at(0)->GetYaxis()->SetRangeUser(0.6,1.);
  list_gr_eff_PU200_prompt_EB.at(0)->Draw("AL");
  list_gr_eff_PU200_prompt_2sigma_EB.at(0)->Draw("same");
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(0)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EB.at(0)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(0)->Draw("same");
  list_gr_eff_PU200_prompt_genMatched_EB.at(0)->Draw("same");
  list_gr_eff_noPU_prompt_EB.at(0)->Draw("same");
  list_gr_eff_noPU_prompt_genMatched_EB.at(0)->Draw("same");
  leg_eff_prompt_EB_dz1->Draw();
  c_prompt_EB_reliso_dz1->Print("plots/ntuple/dz_study/isoeff_prompt_EB_dz1.pdf");

  // prompt Endcap
  TCanvas* c_prompt_EE_reliso_dz1 = new TCanvas("c_prompt_EE_reliso_dz1", "c_prompt_EE_reliso_dz1", 1200, 1200);
  TLegend* leg_eff_prompt_EE_dz1 = new TLegend(0.40, 0.13, 0.88, 0.55);
  leg_eff_prompt_EE_dz1->AddEntry(list_gr_eff_noPU_prompt_genMatched_EE.at(0), "noPU no MTD gen", "L");
  leg_eff_prompt_EE_dz1->AddEntry(list_gr_eff_noPU_prompt_EE.at(0), "noPU no MTD", "L");
  leg_eff_prompt_EE_dz1->AddEntry(list_gr_eff_PU200_prompt_genMatched_EE.at(0), "PU200 no MTD gen", "L");
  leg_eff_prompt_EE_dz1->AddEntry(list_gr_eff_PU200_prompt_EE.at(0), "PU200 no MTD", "L");
  leg_eff_prompt_EE_dz1->AddEntry(list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(0), "PU200 2sigma (PV, track)", "L");
  leg_eff_prompt_EE_dz1->AddEntry(list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(0), "PU200 3sigma (PV, track)", "L");
  leg_eff_prompt_EE_dz1->AddEntry(list_gr_eff_PU200_prompt_2sigma_EE.at(0), "PU200 2sigma (muon, track)", "L");
  leg_eff_prompt_EE_dz1->AddEntry(list_gr_eff_PU200_prompt_3sigma_EE.at(0), "PU200 3sigma (muon, track)", "L");
//  leg_eff_prompt_EE_dz1->AddEntry(list_gr_eff_noPU_prompt_2sigma_EE.at(0), "noPU 2sigma"), "L");
  list_gr_eff_PU200_prompt_EE.at(0)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_2sigma_EE.at(0)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(0)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(0)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_3sigma_EE.at(0)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(0)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(0)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_genMatched_EE.at(0)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_genMatched_EE.at(0)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_EE.at(0)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EE.at(0)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EE.at(0)->SetLineStyle(7);
  c_prompt_EE_reliso_dz1->cd();
  c_prompt_EE_reliso_dz1->SetGrid();
  c_prompt_EE_reliso_dz1->SetLeftMargin(0.12);
  list_gr_eff_PU200_prompt_EE.at(0)->GetXaxis()->SetRangeUser(0.,0.4);
  list_gr_eff_PU200_prompt_EE.at(0)->GetYaxis()->SetRangeUser(0.5,1.);
  list_gr_eff_PU200_prompt_EE.at(0)->Draw("AL");
  list_gr_eff_PU200_prompt_2sigma_EE.at(0)->Draw("same");
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(0)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EE.at(0)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(0)->Draw("same");
  list_gr_eff_PU200_prompt_genMatched_EE.at(0)->Draw("same");
  list_gr_eff_noPU_prompt_EE.at(0)->Draw("same");
  list_gr_eff_noPU_prompt_genMatched_EE.at(0)->Draw("same");
  leg_eff_prompt_EE_dz1->Draw();
  c_prompt_EE_reliso_dz1->Print("plots/ntuple/dz_study/isoeff_prompt_EE_dz1.pdf");

  // nonprompt Barrel
  TCanvas* c_nonprompt_EB_reliso_dz1 = new TCanvas("c_nonprompt_EB_reliso_dz1", "c_nonprompt_EB_reliso_dz1", 1200, 1200);
  TLegend* leg_eff_nonprompt_EB_dz1 = new TLegend(0.40, 0.13, 0.88, 0.40);
  leg_eff_nonprompt_EB_dz1->AddEntry(list_gr_eff_noPU_nonprompt_genMatched_EB.at(0), "noPU no MTD gen", "L");
  leg_eff_nonprompt_EB_dz1->AddEntry(list_gr_eff_noPU_nonprompt_EB.at(0), "noPU no MTD", "L");
  leg_eff_nonprompt_EB_dz1->AddEntry(list_gr_eff_PU200_nonprompt_genMatched_EB.at(0), "PU200 no MTD gen", "L");
  leg_eff_nonprompt_EB_dz1->AddEntry(list_gr_eff_PU200_nonprompt_EB.at(0), "PU200 no MTD", "L");
  leg_eff_nonprompt_EB_dz1->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(0), "PU200 2sigma (PV, track)", "L");
  leg_eff_nonprompt_EB_dz1->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(0), "PU200 3sigma (PV, track)", "L");
  leg_eff_nonprompt_EB_dz1->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EB.at(0), "PU200 2sigma (muon, track)", "L");
  leg_eff_nonprompt_EB_dz1->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EB.at(0), "PU200 3sigma (muon, track)", "L");
 // leg_eff_nonprompt_EB_dz1->AddEntry(list_gr_eff_noPU_nonprompt_2sigma_EB.at(0), "noPU 2sigma", "L");
  list_gr_eff_PU200_nonprompt_EB.at(0)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_2sigma_EB.at(0)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(0)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(0)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_3sigma_EB.at(0)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(0)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(0)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(0)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(0)->SetLineStyle(7);
  list_gr_eff_noPU_nonprompt_EB.at(0)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(0)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(0)->SetLineStyle(7);
  c_nonprompt_EB_reliso_dz1->cd();
  c_nonprompt_EB_reliso_dz1->SetGrid();
  c_nonprompt_EB_reliso_dz1->SetLeftMargin(0.12);
  list_gr_eff_PU200_nonprompt_EB.at(0)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_PU200_nonprompt_EB.at(0)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_nonprompt_EB.at(0)->Draw("AL");
  list_gr_eff_PU200_nonprompt_2sigma_EB.at(0)->Draw("same");
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(0)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EB.at(0)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(0)->Draw("same");
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(0)->Draw("same");
  list_gr_eff_noPU_nonprompt_EB.at(0)->Draw("same");
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(0)->Draw("same");
  leg_eff_nonprompt_EB_dz1->Draw();
  c_nonprompt_EB_reliso_dz1->Print("plots/ntuple/dz_study/isoeff_nonprompt_EB_dz1.pdf");

  // nonprompt Endcap
  TCanvas* c_nonprompt_EE_reliso_dz1 = new TCanvas("c_nonprompt_EE_reliso_dz1", "c_nonprompt_EE_reliso_dz1", 1200, 1200);
  TLegend* leg_eff_nonprompt_EE_dz1 = new TLegend(0.40, 0.13, 0.88, 0.40);
  leg_eff_nonprompt_EE_dz1->AddEntry(list_gr_eff_noPU_nonprompt_genMatched_EE.at(0), "noPU no MTD gen", "L");
  leg_eff_nonprompt_EE_dz1->AddEntry(list_gr_eff_noPU_nonprompt_EE.at(0), "noPU no MTD", "L");
  leg_eff_nonprompt_EE_dz1->AddEntry(list_gr_eff_PU200_nonprompt_genMatched_EE.at(0), "PU200 no MTD gen", "L");
  leg_eff_nonprompt_EE_dz1->AddEntry(list_gr_eff_PU200_nonprompt_EE.at(0), "PU200 no MTD", "L");
  leg_eff_nonprompt_EE_dz1->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(0), "PU200 2sigma (PV, track)", "L");
  leg_eff_nonprompt_EE_dz1->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(0), "PU200 3sigma (PV, track)", "L");
  leg_eff_nonprompt_EE_dz1->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EE.at(0), "PU200 2sigma (muon, track)", "L");
  leg_eff_nonprompt_EE_dz1->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EE.at(0), "PU200 3sigma (muon, track)", "L");
//  leg_eff_nonprompt_EE_dz1->AddEntry(list_gr_eff_noPU_nonprompt_2sigma_EE.at(0), "noPU 2sigma", "L");
  list_gr_eff_PU200_nonprompt_EE.at(0)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_2sigma_EE.at(0)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(0)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(0)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_3sigma_EE.at(0)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(0)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(0)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(0)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(0)->SetLineStyle(7);
  list_gr_eff_noPU_nonprompt_EE.at(0)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(0)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(0)->SetLineStyle(7);
  c_nonprompt_EE_reliso_dz1->cd();
  c_nonprompt_EE_reliso_dz1->SetGrid();
  c_nonprompt_EE_reliso_dz1->SetLeftMargin(0.12);
  list_gr_eff_PU200_nonprompt_EE.at(0)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_PU200_nonprompt_EE.at(0)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_nonprompt_EE.at(0)->Draw("AL");
  list_gr_eff_PU200_nonprompt_2sigma_EE.at(0)->Draw("same");
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(0)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EE.at(0)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(0)->Draw("same");
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(0)->Draw("same");
  list_gr_eff_noPU_nonprompt_EE.at(0)->Draw("same");
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(0)->Draw("same");
  leg_eff_nonprompt_EE_dz1->Draw();
  c_nonprompt_EE_reliso_dz1->Print("plots/ntuple/dz_study/isoeff_nonprompt_EE_dz1.pdf");

  // dz 2
  // prompt Barrel
  TCanvas* c_prompt_EB_reliso_dz2 = new TCanvas("c_prompt_EB_reliso_dz2", "c_prompt_EB_reliso_dz2", 1200, 1200);
  TLegend* leg_eff_prompt_EB_dz2 = new TLegend(0.40, 0.13, 0.88, 0.55);
  leg_eff_prompt_EB_dz2->AddEntry(list_gr_eff_noPU_prompt_genMatched_EB.at(1), "noPU no MTD gen", "L");
  leg_eff_prompt_EB_dz2->AddEntry(list_gr_eff_noPU_prompt_EB.at(1), "noPU no MTD", "L");
  leg_eff_prompt_EB_dz2->AddEntry(list_gr_eff_PU200_prompt_genMatched_EB.at(1), "PU200 no MTD gen", "L");
  leg_eff_prompt_EB_dz2->AddEntry(list_gr_eff_PU200_prompt_EB.at(1), "PU200 no MTD", "L");
  leg_eff_prompt_EB_dz2->AddEntry(list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(1), "PU200 2sigma (PV, track)", "L");
  leg_eff_prompt_EB_dz2->AddEntry(list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(1), "PU200 3sigma (PV, track)", "L");
  leg_eff_prompt_EB_dz2->AddEntry(list_gr_eff_PU200_prompt_2sigma_EB.at(1), "PU200 2sigma (muon, track)", "L");
  leg_eff_prompt_EB_dz2->AddEntry(list_gr_eff_PU200_prompt_3sigma_EB.at(1), "PU200 3sigma (muon, track)", "L");
//  leg_eff_prompt_EB_dz2->AddEntry(list_gr_eff_noPU_prompt_2sigma_EB.at(1), "noPU 2sigma", "L");
  list_gr_eff_PU200_prompt_EB.at(1)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_2sigma_EB.at(1)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(1)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(1)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_3sigma_EB.at(1)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(1)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(1)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_genMatched_EB.at(1)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_genMatched_EB.at(1)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_EB.at(1)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EB.at(1)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EB.at(1)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_2sigma_EB.at(1)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_2sigma_EB.at(1)->SetLineStyle(10);
  c_prompt_EB_reliso_dz2->cd();
  c_prompt_EB_reliso_dz2->SetGrid();
  c_prompt_EB_reliso_dz2->SetLeftMargin(0.12);
  list_gr_eff_PU200_prompt_EB.at(1)->GetXaxis()->SetRangeUser(0.,0.4);
  list_gr_eff_PU200_prompt_EB.at(1)->GetYaxis()->SetRangeUser(0.5,1.);
  list_gr_eff_PU200_prompt_EB.at(1)->Draw("AL");
  list_gr_eff_PU200_prompt_2sigma_EB.at(1)->Draw("same");
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(1)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EB.at(1)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(1)->Draw("same");
  list_gr_eff_PU200_prompt_genMatched_EB.at(1)->Draw("same");
  list_gr_eff_noPU_prompt_EB.at(1)->Draw("same");
  list_gr_eff_noPU_prompt_genMatched_EB.at(1)->Draw("same");
  leg_eff_prompt_EB_dz2->Draw();
  c_prompt_EB_reliso_dz2->Print("plots/ntuple/dz_study/isoeff_prompt_EB_dz2.pdf");

  // prompt Endcap
  TCanvas* c_prompt_EE_reliso_dz2 = new TCanvas("c_prompt_EE_reliso_dz2", "c_prompt_EE_reliso_dz2", 1200, 1200);
  TLegend* leg_eff_prompt_EE_dz2 = new TLegend(0.40, 0.13, 0.88, 0.55);
  leg_eff_prompt_EE_dz2->AddEntry(list_gr_eff_noPU_prompt_genMatched_EE.at(1), "noPU no MTD gen", "L");
  leg_eff_prompt_EE_dz2->AddEntry(list_gr_eff_noPU_prompt_EE.at(1), "noPU no MTD", "L");
  leg_eff_prompt_EE_dz2->AddEntry(list_gr_eff_PU200_prompt_genMatched_EE.at(1), "PU200 no MTD gen", "L");
  leg_eff_prompt_EE_dz2->AddEntry(list_gr_eff_PU200_prompt_EE.at(1), "PU200 no MTD", "L");
  leg_eff_prompt_EE_dz2->AddEntry(list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(1), "PU200 2sigma (PV, track)", "L");
  leg_eff_prompt_EE_dz2->AddEntry(list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(1), "PU200 3sigma (PV, track)", "L");
  leg_eff_prompt_EE_dz2->AddEntry(list_gr_eff_PU200_prompt_2sigma_EE.at(1), "PU200 2sigma (muon, track)", "L");
  leg_eff_prompt_EE_dz2->AddEntry(list_gr_eff_PU200_prompt_3sigma_EE.at(1), "PU200 3sigma (muon, track)", "L");
//  leg_eff_prompt_EE_dz2->AddEntry(list_gr_eff_noPU_prompt_2sigma_EE.at(1), "noPU 2sigma"), "L");
  list_gr_eff_PU200_prompt_EE.at(1)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_2sigma_EE.at(1)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(1)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(1)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_3sigma_EE.at(1)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(1)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(1)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_genMatched_EE.at(1)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_genMatched_EE.at(1)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_EE.at(1)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EE.at(1)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EE.at(1)->SetLineStyle(7);
  c_prompt_EE_reliso_dz2->cd();
  c_prompt_EE_reliso_dz2->SetGrid();
  c_prompt_EE_reliso_dz2->SetLeftMargin(0.12);
  list_gr_eff_PU200_prompt_EE.at(1)->GetXaxis()->SetRangeUser(0.,0.4);
  list_gr_eff_PU200_prompt_EE.at(1)->GetYaxis()->SetRangeUser(0.3,1.);
  list_gr_eff_PU200_prompt_EE.at(1)->Draw("AL");
  list_gr_eff_PU200_prompt_2sigma_EE.at(1)->Draw("same");
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(1)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EE.at(1)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(1)->Draw("same");
  list_gr_eff_PU200_prompt_genMatched_EE.at(1)->Draw("same");
  list_gr_eff_noPU_prompt_EE.at(1)->Draw("same");
  list_gr_eff_noPU_prompt_genMatched_EE.at(1)->Draw("same");
  leg_eff_prompt_EE_dz2->Draw();
  c_prompt_EE_reliso_dz2->Print("plots/ntuple/dz_study/isoeff_prompt_EE_dz2.pdf");

  // nonprompt Barrel
  TCanvas* c_nonprompt_EB_reliso_dz2 = new TCanvas("c_nonprompt_EB_reliso_dz2", "c_nonprompt_EB_reliso_dz2", 1200, 1200);
  TLegend* leg_eff_nonprompt_EB_dz2 = new TLegend(0.40, 0.13, 0.88, 0.40);
  leg_eff_nonprompt_EB_dz2->AddEntry(list_gr_eff_noPU_nonprompt_genMatched_EB.at(1), "noPU no MTD gen", "L");
  leg_eff_nonprompt_EB_dz2->AddEntry(list_gr_eff_noPU_nonprompt_EB.at(1), "noPU no MTD", "L");
  leg_eff_nonprompt_EB_dz2->AddEntry(list_gr_eff_PU200_nonprompt_genMatched_EB.at(1), "PU200 no MTD gen", "L");
  leg_eff_nonprompt_EB_dz2->AddEntry(list_gr_eff_PU200_nonprompt_EB.at(1), "PU200 no MTD", "L");
  leg_eff_nonprompt_EB_dz2->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(1), "PU200 2sigma (PV, track)", "L");
  leg_eff_nonprompt_EB_dz2->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(1), "PU200 3sigma (PV, track)", "L");
  leg_eff_nonprompt_EB_dz2->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EB.at(1), "PU200 2sigma (muon, track)", "L");
  leg_eff_nonprompt_EB_dz2->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EB.at(1), "PU200 3sigma (muon, track)", "L");
 // leg_eff_nonprompt_EB_dz2->AddEntry(list_gr_eff_noPU_nonprompt_2sigma_EB.at(1), "noPU 2sigma", "L");
  list_gr_eff_PU200_nonprompt_EB.at(1)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_2sigma_EB.at(1)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(1)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(1)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_3sigma_EB.at(1)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(1)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(1)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(1)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(1)->SetLineStyle(7);
  list_gr_eff_noPU_nonprompt_EB.at(1)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(1)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(1)->SetLineStyle(7);
  c_nonprompt_EB_reliso_dz2->cd();
  c_nonprompt_EB_reliso_dz2->SetGrid();
  c_nonprompt_EB_reliso_dz2->SetLeftMargin(0.12);
  list_gr_eff_PU200_nonprompt_EB.at(1)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_PU200_nonprompt_EB.at(1)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_nonprompt_EB.at(1)->Draw("AL");
  list_gr_eff_PU200_nonprompt_2sigma_EB.at(1)->Draw("same");
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(1)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EB.at(1)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(1)->Draw("same");
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(1)->Draw("same");
  list_gr_eff_noPU_nonprompt_EB.at(1)->Draw("same");
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(1)->Draw("same");
  leg_eff_nonprompt_EB_dz2->Draw();
  c_nonprompt_EB_reliso_dz2->Print("plots/ntuple/dz_study/isoeff_nonprompt_EB_dz2.pdf");

  // nonprompt Endcap
  TCanvas* c_nonprompt_EE_reliso_dz2 = new TCanvas("c_nonprompt_EE_reliso_dz2", "c_nonprompt_EE_reliso_dz2", 1200, 1200);
  TLegend* leg_eff_nonprompt_EE_dz2 = new TLegend(0.40, 0.13, 0.88, 0.40);
  leg_eff_nonprompt_EE_dz2->AddEntry(list_gr_eff_noPU_nonprompt_genMatched_EE.at(1), "noPU no MTD gen", "L");
  leg_eff_nonprompt_EE_dz2->AddEntry(list_gr_eff_noPU_nonprompt_EE.at(1), "noPU no MTD", "L");
  leg_eff_nonprompt_EE_dz2->AddEntry(list_gr_eff_PU200_nonprompt_genMatched_EE.at(1), "PU200 no MTD gen", "L");
  leg_eff_nonprompt_EE_dz2->AddEntry(list_gr_eff_PU200_nonprompt_EE.at(1), "PU200 no MTD", "L");
  leg_eff_nonprompt_EE_dz2->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(1), "PU200 2sigma (PV, track)", "L");
  leg_eff_nonprompt_EE_dz2->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(1), "PU200 3sigma (PV, track)", "L");
  leg_eff_nonprompt_EE_dz2->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EE.at(1), "PU200 2sigma (muon, track)", "L");
  leg_eff_nonprompt_EE_dz2->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EE.at(1), "PU200 3sigma (muon, track)", "L");
//  leg_eff_nonprompt_EE_dz2->AddEntry(list_gr_eff_noPU_nonprompt_2sigma_EE.at(1), "noPU 2sigma", "L");
  list_gr_eff_PU200_nonprompt_EE.at(1)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_2sigma_EE.at(1)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(1)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(1)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_3sigma_EE.at(1)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(1)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(1)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(1)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(1)->SetLineStyle(7);
  list_gr_eff_noPU_nonprompt_EE.at(1)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(1)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(1)->SetLineStyle(7);
  c_nonprompt_EE_reliso_dz2->cd();
  c_nonprompt_EE_reliso_dz2->SetGrid();
  c_nonprompt_EE_reliso_dz2->SetLeftMargin(0.12);
  list_gr_eff_PU200_nonprompt_EE.at(1)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_PU200_nonprompt_EE.at(1)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_nonprompt_EE.at(1)->Draw("AL");
  list_gr_eff_PU200_nonprompt_2sigma_EE.at(1)->Draw("same");
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(1)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EE.at(1)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(1)->Draw("same");
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(1)->Draw("same");
  list_gr_eff_noPU_nonprompt_EE.at(1)->Draw("same");
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(1)->Draw("same");
  leg_eff_nonprompt_EE_dz2->Draw();
  c_nonprompt_EE_reliso_dz2->Print("plots/ntuple/dz_study/isoeff_nonprompt_EE_dz2.pdf");


  // dz 3
  // prompt Barrel
  TCanvas* c_prompt_EB_reliso_dz3 = new TCanvas("c_prompt_EB_reliso_dz3", "c_prompt_EB_reliso_dz3", 1200, 1200);
  TLegend* leg_eff_prompt_EB_dz3 = new TLegend(0.40, 0.13, 0.88, 0.55);
  leg_eff_prompt_EB_dz3->AddEntry(list_gr_eff_noPU_prompt_genMatched_EB.at(2), "noPU no MTD gen", "L");
  leg_eff_prompt_EB_dz3->AddEntry(list_gr_eff_noPU_prompt_EB.at(2), "noPU no MTD", "L");
  leg_eff_prompt_EB_dz3->AddEntry(list_gr_eff_PU200_prompt_genMatched_EB.at(2), "PU200 no MTD gen", "L");
  leg_eff_prompt_EB_dz3->AddEntry(list_gr_eff_PU200_prompt_EB.at(2), "PU200 no MTD", "L");
  leg_eff_prompt_EB_dz3->AddEntry(list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(2), "PU200 2sigma (PV, track)", "L");
  leg_eff_prompt_EB_dz3->AddEntry(list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(2), "PU200 3sigma (PV, track)", "L");
  leg_eff_prompt_EB_dz3->AddEntry(list_gr_eff_PU200_prompt_2sigma_EB.at(2), "PU200 2sigma (muon, track)", "L");
  leg_eff_prompt_EB_dz3->AddEntry(list_gr_eff_PU200_prompt_3sigma_EB.at(2), "PU200 3sigma (muon, track)", "L");
//  leg_eff_prompt_EB_dz3->AddEntry(list_gr_eff_noPU_prompt_2sigma_EB.at(2), "noPU 2sigma", "L");
  list_gr_eff_PU200_prompt_EB.at(2)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_2sigma_EB.at(2)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(2)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(2)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_3sigma_EB.at(2)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(2)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(2)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_genMatched_EB.at(2)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_genMatched_EB.at(2)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_EB.at(2)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EB.at(2)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EB.at(2)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_2sigma_EB.at(2)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_2sigma_EB.at(2)->SetLineStyle(10);
  c_prompt_EB_reliso_dz3->cd();
  c_prompt_EB_reliso_dz3->SetGrid();
  c_prompt_EB_reliso_dz3->SetLeftMargin(0.12);
  list_gr_eff_PU200_prompt_EB.at(2)->GetXaxis()->SetRangeUser(0.,0.4);
  list_gr_eff_PU200_prompt_EB.at(2)->GetYaxis()->SetRangeUser(0.6,1.);
  list_gr_eff_PU200_prompt_EB.at(2)->Draw("AL");
  list_gr_eff_PU200_prompt_2sigma_EB.at(2)->Draw("same");
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(2)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EB.at(2)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(2)->Draw("same");
  list_gr_eff_PU200_prompt_genMatched_EB.at(2)->Draw("same");
  list_gr_eff_noPU_prompt_EB.at(2)->Draw("same");
  list_gr_eff_noPU_prompt_genMatched_EB.at(2)->Draw("same");
  leg_eff_prompt_EB_dz3->Draw();
  c_prompt_EB_reliso_dz3->Print("plots/ntuple/dz_study/isoeff_prompt_EB_dz3.pdf");

  // prompt Endcap
  TCanvas* c_prompt_EE_reliso_dz3 = new TCanvas("c_prompt_EE_reliso_dz3", "c_prompt_EE_reliso_dz3", 1200, 1200);
  TLegend* leg_eff_prompt_EE_dz3 = new TLegend(0.40, 0.13, 0.88, 0.55);
  leg_eff_prompt_EE_dz3->AddEntry(list_gr_eff_noPU_prompt_genMatched_EE.at(2), "noPU no MTD gen", "L");
  leg_eff_prompt_EE_dz3->AddEntry(list_gr_eff_noPU_prompt_EE.at(2), "noPU no MTD", "L");
  leg_eff_prompt_EE_dz3->AddEntry(list_gr_eff_PU200_prompt_genMatched_EE.at(2), "PU200 no MTD gen", "L");
  leg_eff_prompt_EE_dz3->AddEntry(list_gr_eff_PU200_prompt_EE.at(2), "PU200 no MTD", "L");
  leg_eff_prompt_EE_dz3->AddEntry(list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(2), "PU200 2sigma (PV, track)", "L");
  leg_eff_prompt_EE_dz3->AddEntry(list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(2), "PU200 3sigma (PV, track)", "L");
  leg_eff_prompt_EE_dz3->AddEntry(list_gr_eff_PU200_prompt_2sigma_EE.at(2), "PU200 2sigma (muon, track)", "L");
  leg_eff_prompt_EE_dz3->AddEntry(list_gr_eff_PU200_prompt_3sigma_EE.at(2), "PU200 3sigma (muon, track)", "L");
//  leg_eff_prompt_EE_dz3->AddEntry(list_gr_eff_noPU_prompt_2sigma_EE.at(2), "noPU 2sigma"), "L");
  list_gr_eff_PU200_prompt_EE.at(2)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_2sigma_EE.at(2)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(2)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(2)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_3sigma_EE.at(2)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(2)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(2)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_genMatched_EE.at(2)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_genMatched_EE.at(2)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_EE.at(2)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EE.at(2)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EE.at(2)->SetLineStyle(7);
  c_prompt_EE_reliso_dz3->cd();
  c_prompt_EE_reliso_dz3->SetGrid();
  c_prompt_EE_reliso_dz3->SetLeftMargin(0.12);
  list_gr_eff_PU200_prompt_EE.at(2)->GetXaxis()->SetRangeUser(0.,0.4);
  list_gr_eff_PU200_prompt_EE.at(2)->GetYaxis()->SetRangeUser(0.5,1.);
  list_gr_eff_PU200_prompt_EE.at(2)->Draw("AL");
  list_gr_eff_PU200_prompt_2sigma_EE.at(2)->Draw("same");
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(2)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EE.at(2)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(2)->Draw("same");
  list_gr_eff_PU200_prompt_genMatched_EE.at(2)->Draw("same");
  list_gr_eff_noPU_prompt_EE.at(2)->Draw("same");
  list_gr_eff_noPU_prompt_genMatched_EE.at(2)->Draw("same");
  leg_eff_prompt_EE_dz3->Draw();
  c_prompt_EE_reliso_dz3->Print("plots/ntuple/dz_study/isoeff_prompt_EE_dz3.pdf");

  // nonprompt Barrel
  TCanvas* c_nonprompt_EB_reliso_dz3 = new TCanvas("c_nonprompt_EB_reliso_dz3", "c_nonprompt_EB_reliso_dz3", 1200, 1200);
  TLegend* leg_eff_nonprompt_EB_dz3 = new TLegend(0.40, 0.13, 0.88, 0.40);
  leg_eff_nonprompt_EB_dz3->AddEntry(list_gr_eff_noPU_nonprompt_genMatched_EB.at(2), "noPU no MTD gen", "L");
  leg_eff_nonprompt_EB_dz3->AddEntry(list_gr_eff_noPU_nonprompt_EB.at(2), "noPU no MTD", "L");
  leg_eff_nonprompt_EB_dz3->AddEntry(list_gr_eff_PU200_nonprompt_genMatched_EB.at(2), "PU200 no MTD gen", "L");
  leg_eff_nonprompt_EB_dz3->AddEntry(list_gr_eff_PU200_nonprompt_EB.at(2), "PU200 no MTD", "L");
  leg_eff_nonprompt_EB_dz3->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(2), "PU200 2sigma (PV, track)", "L");
  leg_eff_nonprompt_EB_dz3->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(2), "PU200 3sigma (PV, track)", "L");
  leg_eff_nonprompt_EB_dz3->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EB.at(2), "PU200 2sigma (muon, track)", "L");
  leg_eff_nonprompt_EB_dz3->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EB.at(2), "PU200 3sigma (muon, track)", "L");
 // leg_eff_nonprompt_EB_dz3->AddEntry(list_gr_eff_noPU_nonprompt_2sigma_EB.at(2), "noPU 2sigma", "L");
  list_gr_eff_PU200_nonprompt_EB.at(2)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_2sigma_EB.at(2)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(2)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(2)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_3sigma_EB.at(2)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(2)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(2)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(2)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(2)->SetLineStyle(7);
  list_gr_eff_noPU_nonprompt_EB.at(2)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(2)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(2)->SetLineStyle(7);
  c_nonprompt_EB_reliso_dz3->cd();
  c_nonprompt_EB_reliso_dz3->SetGrid();
  c_nonprompt_EB_reliso_dz3->SetLeftMargin(0.12);
  list_gr_eff_PU200_nonprompt_EB.at(2)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_PU200_nonprompt_EB.at(2)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_nonprompt_EB.at(2)->Draw("AL");
  list_gr_eff_PU200_nonprompt_2sigma_EB.at(2)->Draw("same");
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(2)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EB.at(2)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(2)->Draw("same");
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(2)->Draw("same");
  list_gr_eff_noPU_nonprompt_EB.at(2)->Draw("same");
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(2)->Draw("same");
  leg_eff_nonprompt_EB_dz3->Draw();
  c_nonprompt_EB_reliso_dz3->Print("plots/ntuple/dz_study/isoeff_nonprompt_EB_dz3.pdf");

  // nonprompt Endcap
  TCanvas* c_nonprompt_EE_reliso_dz3 = new TCanvas("c_nonprompt_EE_reliso_dz3", "c_nonprompt_EE_reliso_dz3", 1200, 1200);
  TLegend* leg_eff_nonprompt_EE_dz3 = new TLegend(0.40, 0.13, 0.88, 0.40);
  leg_eff_nonprompt_EE_dz3->AddEntry(list_gr_eff_noPU_nonprompt_genMatched_EE.at(2), "noPU no MTD gen", "L");
  leg_eff_nonprompt_EE_dz3->AddEntry(list_gr_eff_noPU_nonprompt_EE.at(2), "noPU no MTD", "L");
  leg_eff_nonprompt_EE_dz3->AddEntry(list_gr_eff_PU200_nonprompt_genMatched_EE.at(2), "PU200 no MTD gen", "L");
  leg_eff_nonprompt_EE_dz3->AddEntry(list_gr_eff_PU200_nonprompt_EE.at(2), "PU200 no MTD", "L");
  leg_eff_nonprompt_EE_dz3->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(2), "PU200 2sigma (PV, track)", "L");
  leg_eff_nonprompt_EE_dz3->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(2), "PU200 3sigma (PV, track)", "L");
  leg_eff_nonprompt_EE_dz3->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EE.at(2), "PU200 2sigma (muon, track)", "L");
  leg_eff_nonprompt_EE_dz3->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EE.at(2), "PU200 3sigma (muon, track)", "L");
//  leg_eff_nonprompt_EE_dz3->AddEntry(list_gr_eff_noPU_nonprompt_2sigma_EE.at(2), "noPU 2sigma", "L");
  list_gr_eff_PU200_nonprompt_EE.at(2)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_2sigma_EE.at(2)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(2)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(2)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_3sigma_EE.at(2)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(2)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(2)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(2)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(2)->SetLineStyle(7);
  list_gr_eff_noPU_nonprompt_EE.at(2)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(2)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(2)->SetLineStyle(7);
  c_nonprompt_EE_reliso_dz3->cd();
  c_nonprompt_EE_reliso_dz3->SetGrid();
  c_nonprompt_EE_reliso_dz3->SetLeftMargin(0.12);
  list_gr_eff_PU200_nonprompt_EE.at(2)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_PU200_nonprompt_EE.at(2)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_nonprompt_EE.at(2)->Draw("AL");
  list_gr_eff_PU200_nonprompt_2sigma_EE.at(2)->Draw("same");
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(2)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EE.at(2)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(2)->Draw("same");
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(2)->Draw("same");
  list_gr_eff_noPU_nonprompt_EE.at(2)->Draw("same");
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(2)->Draw("same");
  leg_eff_nonprompt_EE_dz3->Draw();
  c_nonprompt_EE_reliso_dz3->Print("plots/ntuple/dz_study/isoeff_nonprompt_EE_dz3.pdf");


  // dz 4
  // prompt Barrel
  TCanvas* c_prompt_EB_reliso_dz4 = new TCanvas("c_prompt_EB_reliso_dz4", "c_prompt_EB_reliso_dz4", 1200, 1200);
  TLegend* leg_eff_prompt_EB_dz4 = new TLegend(0.40, 0.13, 0.88, 0.55);
  leg_eff_prompt_EB_dz4->AddEntry(list_gr_eff_noPU_prompt_genMatched_EB.at(3), "noPU no MTD gen", "L");
  leg_eff_prompt_EB_dz4->AddEntry(list_gr_eff_noPU_prompt_EB.at(3), "noPU no MTD", "L");
  leg_eff_prompt_EB_dz4->AddEntry(list_gr_eff_PU200_prompt_genMatched_EB.at(3), "PU200 no MTD gen", "L");
  leg_eff_prompt_EB_dz4->AddEntry(list_gr_eff_PU200_prompt_EB.at(3), "PU200 no MTD", "L");
  leg_eff_prompt_EB_dz4->AddEntry(list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(3), "PU200 2sigma (PV, track)", "L");
  leg_eff_prompt_EB_dz4->AddEntry(list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(3), "PU200 3sigma (PV, track)", "L");
  leg_eff_prompt_EB_dz4->AddEntry(list_gr_eff_PU200_prompt_2sigma_EB.at(3), "PU200 2sigma (muon, track)", "L");
  leg_eff_prompt_EB_dz4->AddEntry(list_gr_eff_PU200_prompt_3sigma_EB.at(3), "PU200 3sigma (muon, track)", "L");
//  leg_eff_prompt_EB_dz4->AddEntry(list_gr_eff_noPU_prompt_2sigma_EB.at(3), "noPU 2sigma", "L");
  list_gr_eff_PU200_prompt_EB.at(3)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_2sigma_EB.at(3)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(3)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(3)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_3sigma_EB.at(3)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(3)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(3)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_genMatched_EB.at(3)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_genMatched_EB.at(3)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_EB.at(3)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EB.at(3)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EB.at(3)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_2sigma_EB.at(3)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_2sigma_EB.at(3)->SetLineStyle(10);
  c_prompt_EB_reliso_dz4->cd();
  c_prompt_EB_reliso_dz4->SetGrid();
  c_prompt_EB_reliso_dz4->SetLeftMargin(0.12);
  list_gr_eff_PU200_prompt_EB.at(3)->GetXaxis()->SetRangeUser(0.,0.4);
  list_gr_eff_PU200_prompt_EB.at(3)->GetYaxis()->SetRangeUser(0.6,1.);
  list_gr_eff_PU200_prompt_EB.at(3)->Draw("AL");
  list_gr_eff_PU200_prompt_2sigma_EB.at(3)->Draw("same");
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(3)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EB.at(3)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(3)->Draw("same");
  list_gr_eff_PU200_prompt_genMatched_EB.at(3)->Draw("same");
  list_gr_eff_noPU_prompt_EB.at(3)->Draw("same");
  list_gr_eff_noPU_prompt_genMatched_EB.at(3)->Draw("same");
  leg_eff_prompt_EB_dz4->Draw();
  c_prompt_EB_reliso_dz4->Print("plots/ntuple/dz_study/isoeff_prompt_EB_dz4.pdf");

  // prompt Endcap
  TCanvas* c_prompt_EE_reliso_dz4 = new TCanvas("c_prompt_EE_reliso_dz4", "c_prompt_EE_reliso_dz4", 1200, 1200);
  TLegend* leg_eff_prompt_EE_dz4 = new TLegend(0.40, 0.13, 0.88, 0.55);
  leg_eff_prompt_EE_dz4->AddEntry(list_gr_eff_noPU_prompt_genMatched_EE.at(3), "noPU no MTD gen", "L");
  leg_eff_prompt_EE_dz4->AddEntry(list_gr_eff_noPU_prompt_EE.at(3), "noPU no MTD", "L");
  leg_eff_prompt_EE_dz4->AddEntry(list_gr_eff_PU200_prompt_genMatched_EE.at(3), "PU200 no MTD gen", "L");
  leg_eff_prompt_EE_dz4->AddEntry(list_gr_eff_PU200_prompt_EE.at(3), "PU200 no MTD", "L");
  leg_eff_prompt_EE_dz4->AddEntry(list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(3), "PU200 2sigma (PV, track)", "L");
  leg_eff_prompt_EE_dz4->AddEntry(list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(3), "PU200 3sigma (PV, track)", "L");
  leg_eff_prompt_EE_dz4->AddEntry(list_gr_eff_PU200_prompt_2sigma_EE.at(3), "PU200 2sigma (muon, track)", "L");
  leg_eff_prompt_EE_dz4->AddEntry(list_gr_eff_PU200_prompt_3sigma_EE.at(3), "PU200 3sigma (muon, track)", "L");
//  leg_eff_prompt_EE_dz4->AddEntry(list_gr_eff_noPU_prompt_2sigma_EE.at(3), "noPU 2sigma"), "L");
  list_gr_eff_PU200_prompt_EE.at(3)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_2sigma_EE.at(3)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(3)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(3)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_3sigma_EE.at(3)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(3)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(3)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_genMatched_EE.at(3)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_genMatched_EE.at(3)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_EE.at(3)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EE.at(3)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EE.at(3)->SetLineStyle(7);
  c_prompt_EE_reliso_dz4->cd();
  c_prompt_EE_reliso_dz4->SetGrid();
  c_prompt_EE_reliso_dz4->SetLeftMargin(0.12);
  list_gr_eff_PU200_prompt_EE.at(3)->GetXaxis()->SetRangeUser(0.,0.4);
  list_gr_eff_PU200_prompt_EE.at(3)->GetYaxis()->SetRangeUser(0.5,1.);
  list_gr_eff_PU200_prompt_EE.at(3)->Draw("AL");
  list_gr_eff_PU200_prompt_2sigma_EE.at(3)->Draw("same");
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(3)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EE.at(3)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(3)->Draw("same");
  list_gr_eff_PU200_prompt_genMatched_EE.at(3)->Draw("same");
  list_gr_eff_noPU_prompt_EE.at(3)->Draw("same");
  list_gr_eff_noPU_prompt_genMatched_EE.at(3)->Draw("same");
  leg_eff_prompt_EE_dz4->Draw();
  c_prompt_EE_reliso_dz4->Print("plots/ntuple/dz_study/isoeff_prompt_EE_dz4.pdf");

  // nonprompt Barrel
  TCanvas* c_nonprompt_EB_reliso_dz4 = new TCanvas("c_nonprompt_EB_reliso_dz4", "c_nonprompt_EB_reliso_dz4", 1200, 1200);
  TLegend* leg_eff_nonprompt_EB_dz4 = new TLegend(0.40, 0.13, 0.88, 0.40);
  leg_eff_nonprompt_EB_dz4->AddEntry(list_gr_eff_noPU_nonprompt_genMatched_EB.at(3), "noPU no MTD gen", "L");
  leg_eff_nonprompt_EB_dz4->AddEntry(list_gr_eff_noPU_nonprompt_EB.at(3), "noPU no MTD", "L");
  leg_eff_nonprompt_EB_dz4->AddEntry(list_gr_eff_PU200_nonprompt_genMatched_EB.at(3), "PU200 no MTD gen", "L");
  leg_eff_nonprompt_EB_dz4->AddEntry(list_gr_eff_PU200_nonprompt_EB.at(3), "PU200 no MTD", "L");
  leg_eff_nonprompt_EB_dz4->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(3), "PU200 2sigma (PV, track)", "L");
  leg_eff_nonprompt_EB_dz4->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(3), "PU200 3sigma (PV, track)", "L");
  leg_eff_nonprompt_EB_dz4->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EB.at(3), "PU200 2sigma (muon, track)", "L");
  leg_eff_nonprompt_EB_dz4->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EB.at(3), "PU200 3sigma (muon, track)", "L");
 // leg_eff_nonprompt_EB_dz4->AddEntry(list_gr_eff_noPU_nonprompt_2sigma_EB.at(3), "noPU 2sigma", "L");
  list_gr_eff_PU200_nonprompt_EB.at(3)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_2sigma_EB.at(3)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(3)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(3)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_3sigma_EB.at(3)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(3)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(3)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(3)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(3)->SetLineStyle(7);
  list_gr_eff_noPU_nonprompt_EB.at(3)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(3)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(3)->SetLineStyle(7);
  c_nonprompt_EB_reliso_dz4->cd();
  c_nonprompt_EB_reliso_dz4->SetGrid();
  c_nonprompt_EB_reliso_dz4->SetLeftMargin(0.12);
  list_gr_eff_PU200_nonprompt_EB.at(3)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_PU200_nonprompt_EB.at(3)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_nonprompt_EB.at(3)->Draw("AL");
  list_gr_eff_PU200_nonprompt_2sigma_EB.at(3)->Draw("same");
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(3)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EB.at(3)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(3)->Draw("same");
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(3)->Draw("same");
  list_gr_eff_noPU_nonprompt_EB.at(3)->Draw("same");
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(3)->Draw("same");
  leg_eff_nonprompt_EB_dz4->Draw();
  c_nonprompt_EB_reliso_dz4->Print("plots/ntuple/dz_study/isoeff_nonprompt_EB_dz4.pdf");

  // nonprompt Endcap
  TCanvas* c_nonprompt_EE_reliso_dz4 = new TCanvas("c_nonprompt_EE_reliso_dz4", "c_nonprompt_EE_reliso_dz4", 1200, 1200);
  TLegend* leg_eff_nonprompt_EE_dz4 = new TLegend(0.40, 0.13, 0.88, 0.40);
  leg_eff_nonprompt_EE_dz4->AddEntry(list_gr_eff_noPU_nonprompt_genMatched_EE.at(3), "noPU no MTD gen", "L");
  leg_eff_nonprompt_EE_dz4->AddEntry(list_gr_eff_noPU_nonprompt_EE.at(3), "noPU no MTD", "L");
  leg_eff_nonprompt_EE_dz4->AddEntry(list_gr_eff_PU200_nonprompt_genMatched_EE.at(3), "PU200 no MTD gen", "L");
  leg_eff_nonprompt_EE_dz4->AddEntry(list_gr_eff_PU200_nonprompt_EE.at(3), "PU200 no MTD", "L");
  leg_eff_nonprompt_EE_dz4->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(3), "PU200 2sigma (PV, track)", "L");
  leg_eff_nonprompt_EE_dz4->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(3), "PU200 3sigma (PV, track)", "L");
  leg_eff_nonprompt_EE_dz4->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EE.at(3), "PU200 2sigma (muon, track)", "L");
  leg_eff_nonprompt_EE_dz4->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EE.at(3), "PU200 3sigma (muon, track)", "L");
//  leg_eff_nonprompt_EE_dz4->AddEntry(list_gr_eff_noPU_nonprompt_2sigma_EE.at(3), "noPU 2sigma", "L");
  list_gr_eff_PU200_nonprompt_EE.at(3)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_2sigma_EE.at(3)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(3)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(3)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_3sigma_EE.at(3)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(3)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(3)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(3)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(3)->SetLineStyle(7);
  list_gr_eff_noPU_nonprompt_EE.at(3)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(3)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(3)->SetLineStyle(7);
  c_nonprompt_EE_reliso_dz4->cd();
  c_nonprompt_EE_reliso_dz4->SetGrid();
  c_nonprompt_EE_reliso_dz4->SetLeftMargin(0.12);
  list_gr_eff_PU200_nonprompt_EE.at(3)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_PU200_nonprompt_EE.at(3)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_nonprompt_EE.at(3)->Draw("AL");
  list_gr_eff_PU200_nonprompt_2sigma_EE.at(3)->Draw("same");
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(3)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EE.at(3)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(3)->Draw("same");
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(3)->Draw("same");
  list_gr_eff_noPU_nonprompt_EE.at(3)->Draw("same");
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(3)->Draw("same");
  leg_eff_nonprompt_EE_dz4->Draw();
  c_nonprompt_EE_reliso_dz4->Print("plots/ntuple/dz_study/isoeff_nonprompt_EE_dz4.pdf");


  // dz 5
  // prompt Barrel
  TCanvas* c_prompt_EB_reliso_dz5 = new TCanvas("c_prompt_EB_reliso_dz5", "c_prompt_EB_reliso_dz5", 1200, 1200);
  TLegend* leg_eff_prompt_EB_dz5 = new TLegend(0.40, 0.13, 0.88, 0.55);
  leg_eff_prompt_EB_dz5->AddEntry(list_gr_eff_noPU_prompt_genMatched_EB.at(4), "noPU no MTD gen", "L");
  leg_eff_prompt_EB_dz5->AddEntry(list_gr_eff_noPU_prompt_EB.at(4), "noPU no MTD", "L");
  leg_eff_prompt_EB_dz5->AddEntry(list_gr_eff_PU200_prompt_genMatched_EB.at(4), "PU200 no MTD gen", "L");
  leg_eff_prompt_EB_dz5->AddEntry(list_gr_eff_PU200_prompt_EB.at(4), "PU200 no MTD", "L");
  leg_eff_prompt_EB_dz5->AddEntry(list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(4), "PU200 2sigma (PV, track)", "L");
  leg_eff_prompt_EB_dz5->AddEntry(list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(4), "PU200 3sigma (PV, track)", "L");
  leg_eff_prompt_EB_dz5->AddEntry(list_gr_eff_PU200_prompt_2sigma_EB.at(4), "PU200 2sigma (muon, track)", "L");
  leg_eff_prompt_EB_dz5->AddEntry(list_gr_eff_PU200_prompt_3sigma_EB.at(4), "PU200 3sigma (muon, track)", "L");
//  leg_eff_prompt_EB_dz5->AddEntry(list_gr_eff_noPU_prompt_2sigma_EB.at(4), "noPU 2sigma", "L");
  list_gr_eff_PU200_prompt_EB.at(4)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_2sigma_EB.at(4)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(4)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(4)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_3sigma_EB.at(4)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(4)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(4)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_genMatched_EB.at(4)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_genMatched_EB.at(4)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_EB.at(4)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EB.at(4)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EB.at(4)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_2sigma_EB.at(4)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_2sigma_EB.at(4)->SetLineStyle(10);
  c_prompt_EB_reliso_dz5->cd();
  c_prompt_EB_reliso_dz5->SetGrid();
  c_prompt_EB_reliso_dz5->SetLeftMargin(0.12);
  list_gr_eff_PU200_prompt_EB.at(4)->GetXaxis()->SetRangeUser(0.,0.4);
  list_gr_eff_PU200_prompt_EB.at(4)->GetYaxis()->SetRangeUser(0.6,1.);
  list_gr_eff_PU200_prompt_EB.at(4)->Draw("AL");
  list_gr_eff_PU200_prompt_2sigma_EB.at(4)->Draw("same");
  list_gr_eff_PU200_prompt_2sigma_EB_vtx.at(4)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EB.at(4)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EB_vtx.at(4)->Draw("same");
  list_gr_eff_PU200_prompt_genMatched_EB.at(4)->Draw("same");
  list_gr_eff_noPU_prompt_EB.at(4)->Draw("same");
  list_gr_eff_noPU_prompt_genMatched_EB.at(4)->Draw("same");
  leg_eff_prompt_EB_dz5->Draw();
  c_prompt_EB_reliso_dz5->Print("plots/ntuple/dz_study/isoeff_prompt_EB_dz5.pdf");

  // prompt Endcap
  TCanvas* c_prompt_EE_reliso_dz5 = new TCanvas("c_prompt_EE_reliso_dz5", "c_prompt_EE_reliso_dz5", 1200, 1200);
  TLegend* leg_eff_prompt_EE_dz5 = new TLegend(0.40, 0.13, 0.88, 0.55);
  leg_eff_prompt_EE_dz5->AddEntry(list_gr_eff_noPU_prompt_genMatched_EE.at(4), "noPU no MTD gen", "L");
  leg_eff_prompt_EE_dz5->AddEntry(list_gr_eff_noPU_prompt_EE.at(4), "noPU no MTD", "L");
  leg_eff_prompt_EE_dz5->AddEntry(list_gr_eff_PU200_prompt_genMatched_EE.at(4), "PU200 no MTD gen", "L");
  leg_eff_prompt_EE_dz5->AddEntry(list_gr_eff_PU200_prompt_EE.at(4), "PU200 no MTD", "L");
  leg_eff_prompt_EE_dz5->AddEntry(list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(4), "PU200 2sigma (PV, track)", "L");
  leg_eff_prompt_EE_dz5->AddEntry(list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(4), "PU200 3sigma (PV, track)", "L");
  leg_eff_prompt_EE_dz5->AddEntry(list_gr_eff_PU200_prompt_2sigma_EE.at(4), "PU200 2sigma (muon, track)", "L");
  leg_eff_prompt_EE_dz5->AddEntry(list_gr_eff_PU200_prompt_3sigma_EE.at(4), "PU200 3sigma (muon, track)", "L");
//  leg_eff_prompt_EE_dz5->AddEntry(list_gr_eff_noPU_prompt_2sigma_EE.at(4), "noPU 2sigma"), "L");
  list_gr_eff_PU200_prompt_EE.at(4)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_2sigma_EE.at(4)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(4)->SetLineColor(kRed);
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(4)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_3sigma_EE.at(4)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(4)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(4)->SetLineStyle(7);
  list_gr_eff_PU200_prompt_genMatched_EE.at(4)->SetLineColor(kBlack);
  list_gr_eff_PU200_prompt_genMatched_EE.at(4)->SetLineStyle(7);
  list_gr_eff_noPU_prompt_EE.at(4)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EE.at(4)->SetLineColor(kGray+1);
  list_gr_eff_noPU_prompt_genMatched_EE.at(4)->SetLineStyle(7);
  c_prompt_EE_reliso_dz5->cd();
  c_prompt_EE_reliso_dz5->SetGrid();
  c_prompt_EE_reliso_dz5->SetLeftMargin(0.12);
  list_gr_eff_PU200_prompt_EE.at(4)->GetXaxis()->SetRangeUser(0.,0.4);
  list_gr_eff_PU200_prompt_EE.at(4)->GetYaxis()->SetRangeUser(0.5,1.);
  list_gr_eff_PU200_prompt_EE.at(4)->Draw("AL");
  list_gr_eff_PU200_prompt_2sigma_EE.at(4)->Draw("same");
  list_gr_eff_PU200_prompt_2sigma_EE_vtx.at(4)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EE.at(4)->Draw("same");
  list_gr_eff_PU200_prompt_3sigma_EE_vtx.at(4)->Draw("same");
  list_gr_eff_PU200_prompt_genMatched_EE.at(4)->Draw("same");
  list_gr_eff_noPU_prompt_EE.at(4)->Draw("same");
  list_gr_eff_noPU_prompt_genMatched_EE.at(4)->Draw("same");
  leg_eff_prompt_EE_dz5->Draw();
  c_prompt_EE_reliso_dz5->Print("plots/ntuple/dz_study/isoeff_prompt_EE_dz5.pdf");

  // nonprompt Barrel
  TCanvas* c_nonprompt_EB_reliso_dz5 = new TCanvas("c_nonprompt_EB_reliso_dz5", "c_nonprompt_EB_reliso_dz5", 1200, 1200);
  TLegend* leg_eff_nonprompt_EB_dz5 = new TLegend(0.40, 0.13, 0.88, 0.40);
  leg_eff_nonprompt_EB_dz5->AddEntry(list_gr_eff_noPU_nonprompt_genMatched_EB.at(4), "noPU no MTD gen", "L");
  leg_eff_nonprompt_EB_dz5->AddEntry(list_gr_eff_noPU_nonprompt_EB.at(4), "noPU no MTD", "L");
  leg_eff_nonprompt_EB_dz5->AddEntry(list_gr_eff_PU200_nonprompt_genMatched_EB.at(4), "PU200 no MTD gen", "L");
  leg_eff_nonprompt_EB_dz5->AddEntry(list_gr_eff_PU200_nonprompt_EB.at(4), "PU200 no MTD", "L");
  leg_eff_nonprompt_EB_dz5->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(4), "PU200 2sigma (PV, track)", "L");
  leg_eff_nonprompt_EB_dz5->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(4), "PU200 3sigma (PV, track)", "L");
  leg_eff_nonprompt_EB_dz5->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EB.at(4), "PU200 2sigma (muon, track)", "L");
  leg_eff_nonprompt_EB_dz5->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EB.at(4), "PU200 3sigma (muon, track)", "L");
 // leg_eff_nonprompt_EB_dz5->AddEntry(list_gr_eff_noPU_nonprompt_2sigma_EB.at(4), "noPU 2sigma", "L");
  list_gr_eff_PU200_nonprompt_EB.at(4)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_2sigma_EB.at(4)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(4)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(4)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_3sigma_EB.at(4)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(4)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(4)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(4)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(4)->SetLineStyle(7);
  list_gr_eff_noPU_nonprompt_EB.at(4)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(4)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(4)->SetLineStyle(7);
  c_nonprompt_EB_reliso_dz5->cd();
  c_nonprompt_EB_reliso_dz5->SetGrid();
  c_nonprompt_EB_reliso_dz5->SetLeftMargin(0.12);
  list_gr_eff_PU200_nonprompt_EB.at(4)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_PU200_nonprompt_EB.at(4)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_nonprompt_EB.at(4)->Draw("AL");
  list_gr_eff_PU200_nonprompt_2sigma_EB.at(4)->Draw("same");
  list_gr_eff_PU200_nonprompt_2sigma_EB_vtx.at(4)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EB.at(4)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EB_vtx.at(4)->Draw("same");
  list_gr_eff_PU200_nonprompt_genMatched_EB.at(4)->Draw("same");
  list_gr_eff_noPU_nonprompt_EB.at(4)->Draw("same");
  list_gr_eff_noPU_nonprompt_genMatched_EB.at(4)->Draw("same");
  leg_eff_nonprompt_EB_dz5->Draw();
  c_nonprompt_EB_reliso_dz5->Print("plots/ntuple/dz_study/isoeff_nonprompt_EB_dz5.pdf");

  // nonprompt Endcap
  TCanvas* c_nonprompt_EE_reliso_dz5 = new TCanvas("c_nonprompt_EE_reliso_dz5", "c_nonprompt_EE_reliso_dz5", 1200, 1200);
  TLegend* leg_eff_nonprompt_EE_dz5 = new TLegend(0.40, 0.13, 0.88, 0.40);
  leg_eff_nonprompt_EE_dz5->AddEntry(list_gr_eff_noPU_nonprompt_genMatched_EE.at(4), "noPU no MTD gen", "L");
  leg_eff_nonprompt_EE_dz5->AddEntry(list_gr_eff_noPU_nonprompt_EE.at(4), "noPU no MTD", "L");
  leg_eff_nonprompt_EE_dz5->AddEntry(list_gr_eff_PU200_nonprompt_genMatched_EE.at(4), "PU200 no MTD gen", "L");
  leg_eff_nonprompt_EE_dz5->AddEntry(list_gr_eff_PU200_nonprompt_EE.at(4), "PU200 no MTD", "L");
  leg_eff_nonprompt_EE_dz5->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(4), "PU200 2sigma (PV, track)", "L");
  leg_eff_nonprompt_EE_dz5->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(4), "PU200 3sigma (PV, track)", "L");
  leg_eff_nonprompt_EE_dz5->AddEntry(list_gr_eff_PU200_nonprompt_2sigma_EE.at(4), "PU200 2sigma (muon, track)", "L");
  leg_eff_nonprompt_EE_dz5->AddEntry(list_gr_eff_PU200_nonprompt_3sigma_EE.at(4), "PU200 3sigma (muon, track)", "L");
//  leg_eff_nonprompt_EE_dz5->AddEntry(list_gr_eff_noPU_nonprompt_2sigma_EE.at(4), "noPU 2sigma", "L");
  list_gr_eff_PU200_nonprompt_EE.at(4)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_2sigma_EE.at(4)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(4)->SetLineColor(kRed);
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(4)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_3sigma_EE.at(4)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(4)->SetLineColor(kAzure+1);
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(4)->SetLineStyle(7);
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(4)->SetLineColor(kBlack);
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(4)->SetLineStyle(7);
  list_gr_eff_noPU_nonprompt_EE.at(4)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(4)->SetLineColor(kGray+1);
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(4)->SetLineStyle(7);
  c_nonprompt_EE_reliso_dz5->cd();
  c_nonprompt_EE_reliso_dz5->SetGrid();
  c_nonprompt_EE_reliso_dz5->SetLeftMargin(0.12);
  list_gr_eff_PU200_nonprompt_EE.at(4)->GetXaxis()->SetRangeUser(0.,1.2);
  list_gr_eff_PU200_nonprompt_EE.at(4)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_eff_PU200_nonprompt_EE.at(4)->Draw("AL");
  list_gr_eff_PU200_nonprompt_2sigma_EE.at(4)->Draw("same");
  list_gr_eff_PU200_nonprompt_2sigma_EE_vtx.at(4)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EE.at(4)->Draw("same");
  list_gr_eff_PU200_nonprompt_3sigma_EE_vtx.at(4)->Draw("same");
  list_gr_eff_PU200_nonprompt_genMatched_EE.at(4)->Draw("same");
  list_gr_eff_noPU_nonprompt_EE.at(4)->Draw("same");
  list_gr_eff_noPU_nonprompt_genMatched_EE.at(4)->Draw("same");
  leg_eff_nonprompt_EE_dz5->Draw();
  c_nonprompt_EE_reliso_dz5->Print("plots/ntuple/dz_study/isoeff_nonprompt_EE_dz5.pdf");

  // test end



  /////////////////////
  ///// ROC Curve /////
  /////////////////////

  // Define TGraph
  // (muon, track)
  vector<TGraph*> list_gr_roc_PU200_EB,            list_gr_roc_PU200_EE;
  vector<TGraph*> list_gr_roc_PU200_1sigma_EB,     list_gr_roc_PU200_1sigma_EE;
  vector<TGraph*> list_gr_roc_PU200_2sigma_EB,     list_gr_roc_PU200_2sigma_EE;
  vector<TGraph*> list_gr_roc_PU200_3sigma_EB,     list_gr_roc_PU200_3sigma_EE;
  vector<TGraph*> list_gr_roc_PU200_4sigma_EB,     list_gr_roc_PU200_4sigma_EE;
  vector<TGraph*> list_gr_roc_PU200_genMatched_EB, list_gr_roc_PU200_genMatched_EE;
  vector<TGraph*> list_gr_roc_noPU_EB,             list_gr_roc_noPU_EE;
  vector<TGraph*> list_gr_roc_noPU_1sigma_EB,      list_gr_roc_noPU_1sigma_EE;
  vector<TGraph*> list_gr_roc_noPU_2sigma_EB,      list_gr_roc_noPU_2sigma_EE;
  vector<TGraph*> list_gr_roc_noPU_3sigma_EB,      list_gr_roc_noPU_3sigma_EE;
  vector<TGraph*> list_gr_roc_noPU_4sigma_EB,      list_gr_roc_noPU_4sigma_EE;
  vector<TGraph*> list_gr_roc_noPU_genMatched_EB,  list_gr_roc_noPU_genMatched_EE;

  // (PV, track)
  vector<TGraph*> list_gr_roc_PU200_EB_vtx,            list_gr_roc_PU200_EE_vtx;
  vector<TGraph*> list_gr_roc_PU200_1sigma_EB_vtx,     list_gr_roc_PU200_1sigma_EE_vtx;
  vector<TGraph*> list_gr_roc_PU200_2sigma_EB_vtx,     list_gr_roc_PU200_2sigma_EE_vtx;
  vector<TGraph*> list_gr_roc_PU200_3sigma_EB_vtx,     list_gr_roc_PU200_3sigma_EE_vtx;
  vector<TGraph*> list_gr_roc_PU200_4sigma_EB_vtx,     list_gr_roc_PU200_4sigma_EE_vtx;
  vector<TGraph*> list_gr_roc_PU200_genMatched_EB_vtx, list_gr_roc_PU200_genMatched_EE_vtx;
  vector<TGraph*> list_gr_roc_noPU_EB_vtx,             list_gr_roc_noPU_EE_vtx;
  vector<TGraph*> list_gr_roc_noPU_1sigma_EB_vtx,      list_gr_roc_noPU_1sigma_EE_vtx;
  vector<TGraph*> list_gr_roc_noPU_2sigma_EB_vtx,      list_gr_roc_noPU_2sigma_EE_vtx;
  vector<TGraph*> list_gr_roc_noPU_3sigma_EB_vtx,      list_gr_roc_noPU_3sigma_EE_vtx;
  vector<TGraph*> list_gr_roc_noPU_4sigma_EB_vtx,      list_gr_roc_noPU_4sigma_EE_vtx;
  vector<TGraph*> list_gr_roc_noPU_genMatched_EB_vtx,  list_gr_roc_noPU_genMatched_EE_vtx;


  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
	// (muon, track)
	TGraph* gr_roc_PU200_EB = new TGraph(); TGraph* gr_roc_PU200_EE = new TGraph();
	TGraph* gr_roc_PU200_1sigma_EB = new TGraph(); TGraph* gr_roc_PU200_1sigma_EE = new TGraph();
	TGraph* gr_roc_PU200_2sigma_EB = new TGraph(); TGraph* gr_roc_PU200_2sigma_EE = new TGraph();
	TGraph* gr_roc_PU200_3sigma_EB = new TGraph(); TGraph* gr_roc_PU200_3sigma_EE = new TGraph();
	TGraph* gr_roc_PU200_4sigma_EB = new TGraph(); TGraph* gr_roc_PU200_4sigma_EE = new TGraph();
	TGraph* gr_roc_PU200_genMatched_EB = new TGraph(); TGraph* gr_roc_PU200_genMatched_EE = new TGraph();
	TGraph* gr_roc_noPU_EB = new TGraph();  TGraph* gr_roc_noPU_EE = new TGraph();
	TGraph* gr_roc_noPU_1sigma_EB = new TGraph(); TGraph* gr_roc_noPU_1sigma_EE = new TGraph();
	TGraph* gr_roc_noPU_2sigma_EB = new TGraph(); TGraph* gr_roc_noPU_2sigma_EE = new TGraph();
	TGraph* gr_roc_noPU_3sigma_EB = new TGraph(); TGraph* gr_roc_noPU_3sigma_EE = new TGraph();
	TGraph* gr_roc_noPU_4sigma_EB = new TGraph(); TGraph* gr_roc_noPU_4sigma_EE = new TGraph();
	TGraph* gr_roc_noPU_genMatched_EB = new TGraph(); TGraph* gr_roc_noPU_genMatched_EE = new TGraph();
	// (PV, track)
	TGraph* gr_roc_PU200_EB_vtx = new TGraph(); TGraph* gr_roc_PU200_EE_vtx = new TGraph();
	TGraph* gr_roc_PU200_1sigma_EB_vtx = new TGraph(); TGraph* gr_roc_PU200_1sigma_EE_vtx = new TGraph();
	TGraph* gr_roc_PU200_2sigma_EB_vtx = new TGraph(); TGraph* gr_roc_PU200_2sigma_EE_vtx = new TGraph();
	TGraph* gr_roc_PU200_3sigma_EB_vtx = new TGraph(); TGraph* gr_roc_PU200_3sigma_EE_vtx = new TGraph();
	TGraph* gr_roc_PU200_4sigma_EB_vtx = new TGraph(); TGraph* gr_roc_PU200_4sigma_EE_vtx = new TGraph();
	TGraph* gr_roc_PU200_genMatched_EB_vtx = new TGraph(); TGraph* gr_roc_PU200_genMatched_EE_vtx = new TGraph();
	TGraph* gr_roc_noPU_EB_vtx = new TGraph();  TGraph* gr_roc_noPU_EE_vtx = new TGraph();
	TGraph* gr_roc_noPU_1sigma_EB_vtx = new TGraph(); TGraph* gr_roc_noPU_1sigma_EE_vtx = new TGraph();
	TGraph* gr_roc_noPU_2sigma_EB_vtx = new TGraph(); TGraph* gr_roc_noPU_2sigma_EE_vtx = new TGraph();
	TGraph* gr_roc_noPU_3sigma_EB_vtx = new TGraph(); TGraph* gr_roc_noPU_3sigma_EE_vtx = new TGraph();
	TGraph* gr_roc_noPU_4sigma_EB_vtx = new TGraph(); TGraph* gr_roc_noPU_4sigma_EE_vtx = new TGraph();
	TGraph* gr_roc_noPU_genMatched_EB_vtx = new TGraph(); TGraph* gr_roc_noPU_genMatched_EE_vtx = new TGraph();


    // Draw ROC Curve
	for(int i=0; i<nbin+1; i++) {
	  // (muon, track)
	  gr_roc_PU200_EB->SetPoint(i, 	      prompt_eff_PU200_EB.at(idz).at(i), 	      nonprompt_eff_PU200_EB.at(idz).at(i));
	  gr_roc_PU200_EE->SetPoint(i, 	      prompt_eff_PU200_EE.at(idz).at(i), 	      nonprompt_eff_PU200_EE.at(idz).at(i));
	  gr_roc_PU200_1sigma_EB->SetPoint(i, prompt_eff_PU200_1sigma_EB.at(idz).at(i), nonprompt_eff_PU200_1sigma_EB.at(idz).at(i));
	  gr_roc_PU200_1sigma_EE->SetPoint(i, prompt_eff_PU200_1sigma_EE.at(idz).at(i), nonprompt_eff_PU200_1sigma_EE.at(idz).at(i));
	  gr_roc_PU200_2sigma_EB->SetPoint(i, prompt_eff_PU200_2sigma_EB.at(idz).at(i), nonprompt_eff_PU200_2sigma_EB.at(idz).at(i));
	  gr_roc_PU200_2sigma_EE->SetPoint(i, prompt_eff_PU200_2sigma_EE.at(idz).at(i), nonprompt_eff_PU200_2sigma_EE.at(idz).at(i));
	  gr_roc_PU200_3sigma_EB->SetPoint(i, prompt_eff_PU200_3sigma_EB.at(idz).at(i), nonprompt_eff_PU200_3sigma_EB.at(idz).at(i));
	  gr_roc_PU200_3sigma_EE->SetPoint(i, prompt_eff_PU200_3sigma_EE.at(idz).at(i), nonprompt_eff_PU200_3sigma_EE.at(idz).at(i));
	  gr_roc_PU200_4sigma_EB->SetPoint(i, prompt_eff_PU200_4sigma_EB.at(idz).at(i), nonprompt_eff_PU200_4sigma_EB.at(idz).at(i));
	  gr_roc_PU200_4sigma_EE->SetPoint(i, prompt_eff_PU200_4sigma_EE.at(idz).at(i), nonprompt_eff_PU200_4sigma_EE.at(idz).at(i));
	  gr_roc_PU200_genMatched_EB->SetPoint(i, prompt_eff_PU200_genMatched_EB.at(idz).at(i), nonprompt_eff_PU200_genMatched_EB.at(idz).at(i));
	  gr_roc_PU200_genMatched_EE->SetPoint(i, prompt_eff_PU200_genMatched_EE.at(idz).at(i), nonprompt_eff_PU200_genMatched_EE.at(idz).at(i));
	  gr_roc_noPU_EB-> SetPoint(i, 	prompt_eff_noPU_EB.at(idz).at(i), 	nonprompt_eff_noPU_EB.at(idz).at(i));
	  gr_roc_noPU_EE-> SetPoint(i, 	prompt_eff_noPU_EE.at(idz).at(i), 	nonprompt_eff_noPU_EE.at(idz).at(i));
	  gr_roc_noPU_1sigma_EB->SetPoint(i, prompt_eff_noPU_1sigma_EB.at(idz).at(i), nonprompt_eff_noPU_1sigma_EB.at(idz).at(i));
	  gr_roc_noPU_1sigma_EE->SetPoint(i, prompt_eff_noPU_1sigma_EE.at(idz).at(i), nonprompt_eff_noPU_1sigma_EE.at(idz).at(i));
	  gr_roc_noPU_2sigma_EB->SetPoint(i, prompt_eff_noPU_2sigma_EB.at(idz).at(i), nonprompt_eff_noPU_2sigma_EB.at(idz).at(i));
	  gr_roc_noPU_2sigma_EE->SetPoint(i, prompt_eff_noPU_2sigma_EE.at(idz).at(i), nonprompt_eff_noPU_2sigma_EE.at(idz).at(i));
	  gr_roc_noPU_3sigma_EB->SetPoint(i, prompt_eff_noPU_3sigma_EB.at(idz).at(i), nonprompt_eff_noPU_3sigma_EB.at(idz).at(i));
	  gr_roc_noPU_3sigma_EE->SetPoint(i, prompt_eff_noPU_3sigma_EE.at(idz).at(i), nonprompt_eff_noPU_3sigma_EE.at(idz).at(i));
	  gr_roc_noPU_4sigma_EB->SetPoint(i, prompt_eff_noPU_4sigma_EB.at(idz).at(i), nonprompt_eff_noPU_4sigma_EB.at(idz).at(i));
	  gr_roc_noPU_4sigma_EE->SetPoint(i, prompt_eff_noPU_4sigma_EE.at(idz).at(i), nonprompt_eff_noPU_4sigma_EE.at(idz).at(i));
	  gr_roc_noPU_genMatched_EB->SetPoint(i, prompt_eff_noPU_genMatched_EB.at(idz).at(i), nonprompt_eff_noPU_genMatched_EB.at(idz).at(i));
	  gr_roc_noPU_genMatched_EE->SetPoint(i, prompt_eff_noPU_genMatched_EE.at(idz).at(i), nonprompt_eff_noPU_genMatched_EE.at(idz).at(i));

	  // (PV, track)
	  gr_roc_PU200_EB_vtx->SetPoint(i, 	      prompt_eff_PU200_EB_vtx.at(idz).at(i), 	      nonprompt_eff_PU200_EB_vtx.at(idz).at(i));
	  gr_roc_PU200_EE_vtx->SetPoint(i, 	      prompt_eff_PU200_EE_vtx.at(idz).at(i), 	      nonprompt_eff_PU200_EE_vtx.at(idz).at(i));
	  gr_roc_PU200_1sigma_EB_vtx->SetPoint(i, prompt_eff_PU200_1sigma_EB_vtx.at(idz).at(i), nonprompt_eff_PU200_1sigma_EB_vtx.at(idz).at(i));
	  gr_roc_PU200_1sigma_EE_vtx->SetPoint(i, prompt_eff_PU200_1sigma_EE_vtx.at(idz).at(i), nonprompt_eff_PU200_1sigma_EE_vtx.at(idz).at(i));
	  gr_roc_PU200_2sigma_EB_vtx->SetPoint(i, prompt_eff_PU200_2sigma_EB_vtx.at(idz).at(i), nonprompt_eff_PU200_2sigma_EB_vtx.at(idz).at(i));
	  gr_roc_PU200_2sigma_EE_vtx->SetPoint(i, prompt_eff_PU200_2sigma_EE_vtx.at(idz).at(i), nonprompt_eff_PU200_2sigma_EE_vtx.at(idz).at(i));
	  gr_roc_PU200_3sigma_EB_vtx->SetPoint(i, prompt_eff_PU200_3sigma_EB_vtx.at(idz).at(i), nonprompt_eff_PU200_3sigma_EB_vtx.at(idz).at(i));
	  gr_roc_PU200_3sigma_EE_vtx->SetPoint(i, prompt_eff_PU200_3sigma_EE_vtx.at(idz).at(i), nonprompt_eff_PU200_3sigma_EE_vtx.at(idz).at(i));
	  gr_roc_PU200_4sigma_EB_vtx->SetPoint(i, prompt_eff_PU200_4sigma_EB_vtx.at(idz).at(i), nonprompt_eff_PU200_4sigma_EB_vtx.at(idz).at(i));
	  gr_roc_PU200_4sigma_EE_vtx->SetPoint(i, prompt_eff_PU200_4sigma_EE_vtx.at(idz).at(i), nonprompt_eff_PU200_4sigma_EE_vtx.at(idz).at(i));
	  gr_roc_PU200_genMatched_EB_vtx->SetPoint(i, prompt_eff_PU200_genMatched_EB_vtx.at(idz).at(i), nonprompt_eff_PU200_genMatched_EB_vtx.at(idz).at(i));
	  gr_roc_PU200_genMatched_EE_vtx->SetPoint(i, prompt_eff_PU200_genMatched_EE_vtx.at(idz).at(i), nonprompt_eff_PU200_genMatched_EE_vtx.at(idz).at(i));
	  gr_roc_noPU_EB_vtx-> SetPoint(i, 	prompt_eff_noPU_EB_vtx.at(idz).at(i), 	nonprompt_eff_noPU_EB_vtx.at(idz).at(i));
	  gr_roc_noPU_EE_vtx-> SetPoint(i, 	prompt_eff_noPU_EE_vtx.at(idz).at(i), 	nonprompt_eff_noPU_EE_vtx.at(idz).at(i));
	  gr_roc_noPU_1sigma_EB_vtx->SetPoint(i, prompt_eff_noPU_1sigma_EB_vtx.at(idz).at(i), nonprompt_eff_noPU_1sigma_EB_vtx.at(idz).at(i));
	  gr_roc_noPU_1sigma_EE_vtx->SetPoint(i, prompt_eff_noPU_1sigma_EE_vtx.at(idz).at(i), nonprompt_eff_noPU_1sigma_EE_vtx.at(idz).at(i));
	  gr_roc_noPU_2sigma_EB_vtx->SetPoint(i, prompt_eff_noPU_2sigma_EB_vtx.at(idz).at(i), nonprompt_eff_noPU_2sigma_EB_vtx.at(idz).at(i));
	  gr_roc_noPU_2sigma_EE_vtx->SetPoint(i, prompt_eff_noPU_2sigma_EE_vtx.at(idz).at(i), nonprompt_eff_noPU_2sigma_EE_vtx.at(idz).at(i));
	  gr_roc_noPU_3sigma_EB_vtx->SetPoint(i, prompt_eff_noPU_3sigma_EB_vtx.at(idz).at(i), nonprompt_eff_noPU_3sigma_EB_vtx.at(idz).at(i));
	  gr_roc_noPU_3sigma_EE_vtx->SetPoint(i, prompt_eff_noPU_3sigma_EE_vtx.at(idz).at(i), nonprompt_eff_noPU_3sigma_EE_vtx.at(idz).at(i));
	  gr_roc_noPU_4sigma_EB_vtx->SetPoint(i, prompt_eff_noPU_4sigma_EB_vtx.at(idz).at(i), nonprompt_eff_noPU_4sigma_EB_vtx.at(idz).at(i));
	  gr_roc_noPU_4sigma_EE_vtx->SetPoint(i, prompt_eff_noPU_4sigma_EE_vtx.at(idz).at(i), nonprompt_eff_noPU_4sigma_EE_vtx.at(idz).at(i));
	  gr_roc_noPU_genMatched_EB_vtx->SetPoint(i, prompt_eff_noPU_genMatched_EB_vtx.at(idz).at(i), nonprompt_eff_noPU_genMatched_EB_vtx.at(idz).at(i));
	  gr_roc_noPU_genMatched_EE_vtx->SetPoint(i, prompt_eff_noPU_genMatched_EE_vtx.at(idz).at(i), nonprompt_eff_noPU_genMatched_EE_vtx.at(idz).at(i));

	}
	// (muon, track)
	list_gr_roc_PU200_EB.emplace_back(gr_roc_PU200_EB);
	list_gr_roc_PU200_EE.emplace_back(gr_roc_PU200_EE);
	list_gr_roc_PU200_1sigma_EB.emplace_back(gr_roc_PU200_1sigma_EB);
	list_gr_roc_PU200_1sigma_EE.emplace_back(gr_roc_PU200_1sigma_EE);
	list_gr_roc_PU200_2sigma_EB.emplace_back(gr_roc_PU200_2sigma_EB);
	list_gr_roc_PU200_2sigma_EE.emplace_back(gr_roc_PU200_2sigma_EE);
	list_gr_roc_PU200_3sigma_EB.emplace_back(gr_roc_PU200_3sigma_EB);
	list_gr_roc_PU200_3sigma_EE.emplace_back(gr_roc_PU200_3sigma_EE);
	list_gr_roc_PU200_4sigma_EB.emplace_back(gr_roc_PU200_4sigma_EB);
	list_gr_roc_PU200_4sigma_EE.emplace_back(gr_roc_PU200_4sigma_EE);
	list_gr_roc_PU200_genMatched_EB.emplace_back(gr_roc_PU200_genMatched_EB);
	list_gr_roc_PU200_genMatched_EE.emplace_back(gr_roc_PU200_genMatched_EE);
	list_gr_roc_noPU_EB.emplace_back(gr_roc_noPU_EB);
	list_gr_roc_noPU_EE.emplace_back(gr_roc_noPU_EE);
	list_gr_roc_noPU_1sigma_EB.emplace_back(gr_roc_noPU_1sigma_EB);
	list_gr_roc_noPU_1sigma_EE.emplace_back(gr_roc_noPU_1sigma_EE);
	list_gr_roc_noPU_2sigma_EB.emplace_back(gr_roc_noPU_2sigma_EB);
	list_gr_roc_noPU_2sigma_EE.emplace_back(gr_roc_noPU_2sigma_EE);
	list_gr_roc_noPU_3sigma_EB.emplace_back(gr_roc_noPU_3sigma_EB);
	list_gr_roc_noPU_3sigma_EE.emplace_back(gr_roc_noPU_3sigma_EE);
	list_gr_roc_noPU_4sigma_EB.emplace_back(gr_roc_noPU_4sigma_EB);
	list_gr_roc_noPU_4sigma_EE.emplace_back(gr_roc_noPU_4sigma_EE);
	list_gr_roc_noPU_genMatched_EB.emplace_back(gr_roc_noPU_genMatched_EB);
	list_gr_roc_noPU_genMatched_EE.emplace_back(gr_roc_noPU_genMatched_EE);

	// (PV, track)
	list_gr_roc_PU200_EB_vtx.emplace_back(gr_roc_PU200_EB_vtx);
	list_gr_roc_PU200_EE_vtx.emplace_back(gr_roc_PU200_EE_vtx);
	list_gr_roc_PU200_1sigma_EB_vtx.emplace_back(gr_roc_PU200_1sigma_EB_vtx);
	list_gr_roc_PU200_1sigma_EE_vtx.emplace_back(gr_roc_PU200_1sigma_EE_vtx);
	list_gr_roc_PU200_2sigma_EB_vtx.emplace_back(gr_roc_PU200_2sigma_EB_vtx);
	list_gr_roc_PU200_2sigma_EE_vtx.emplace_back(gr_roc_PU200_2sigma_EE_vtx);
	list_gr_roc_PU200_3sigma_EB_vtx.emplace_back(gr_roc_PU200_3sigma_EB_vtx);
	list_gr_roc_PU200_3sigma_EE_vtx.emplace_back(gr_roc_PU200_3sigma_EE_vtx);
	list_gr_roc_PU200_4sigma_EB_vtx.emplace_back(gr_roc_PU200_4sigma_EB_vtx);
	list_gr_roc_PU200_4sigma_EE_vtx.emplace_back(gr_roc_PU200_4sigma_EE_vtx);
	list_gr_roc_PU200_genMatched_EB_vtx.emplace_back(gr_roc_PU200_genMatched_EB_vtx);
	list_gr_roc_PU200_genMatched_EE_vtx.emplace_back(gr_roc_PU200_genMatched_EE_vtx);
	list_gr_roc_noPU_EB_vtx.emplace_back(gr_roc_noPU_EB_vtx);
	list_gr_roc_noPU_EE_vtx.emplace_back(gr_roc_noPU_EE_vtx);
	list_gr_roc_noPU_1sigma_EB_vtx.emplace_back(gr_roc_noPU_1sigma_EB_vtx);
	list_gr_roc_noPU_1sigma_EE_vtx.emplace_back(gr_roc_noPU_1sigma_EE_vtx);
	list_gr_roc_noPU_2sigma_EB_vtx.emplace_back(gr_roc_noPU_2sigma_EB_vtx);
	list_gr_roc_noPU_2sigma_EE_vtx.emplace_back(gr_roc_noPU_2sigma_EE_vtx);
	list_gr_roc_noPU_3sigma_EB_vtx.emplace_back(gr_roc_noPU_3sigma_EB_vtx);
	list_gr_roc_noPU_3sigma_EE_vtx.emplace_back(gr_roc_noPU_3sigma_EE_vtx);
	list_gr_roc_noPU_4sigma_EB_vtx.emplace_back(gr_roc_noPU_4sigma_EB_vtx);
	list_gr_roc_noPU_4sigma_EE_vtx.emplace_back(gr_roc_noPU_4sigma_EE_vtx);
	list_gr_roc_noPU_genMatched_EB_vtx.emplace_back(gr_roc_noPU_genMatched_EB_vtx);
	list_gr_roc_noPU_genMatched_EE_vtx.emplace_back(gr_roc_noPU_genMatched_EE_vtx);


  }


  ///////////////
  // Cosmetics //
  ///////////////
  // PU200
  
  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
	// (muon, track)
	list_gr_roc_PU200_EB.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_EB.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_EE.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_EE.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_1sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_1sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_1sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_1sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_2sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_2sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_2sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_2sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_3sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_3sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_3sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_3sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_4sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_4sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_4sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_4sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_genMatched_EB.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_genMatched_EB.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_genMatched_EE.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_genMatched_EE.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_EB.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_EB.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_EE.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_EE.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_1sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_1sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_1sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_1sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_2sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_2sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_2sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_2sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_3sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_3sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_3sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_3sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_4sigma_EB.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_4sigma_EB.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_4sigma_EE.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_4sigma_EE.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_genMatched_EB.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_genMatched_EB.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_genMatched_EE.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_genMatched_EE.at(idz)->SetLineWidth(2);
	// (PV, track)
	list_gr_roc_PU200_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_1sigma_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_1sigma_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_1sigma_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_1sigma_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_2sigma_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_2sigma_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_2sigma_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_2sigma_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_3sigma_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_3sigma_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_3sigma_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_3sigma_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_4sigma_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_4sigma_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_4sigma_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_4sigma_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_genMatched_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_genMatched_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_PU200_genMatched_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_PU200_genMatched_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_1sigma_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_1sigma_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_1sigma_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_1sigma_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_2sigma_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_2sigma_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_2sigma_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_2sigma_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_3sigma_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_3sigma_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_3sigma_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_3sigma_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_4sigma_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_4sigma_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_4sigma_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_4sigma_EE_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_genMatched_EB_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_genMatched_EB_vtx.at(idz)->SetLineWidth(2);
	list_gr_roc_noPU_genMatched_EE_vtx.at(idz)->SetLineColor(idz+1); list_gr_roc_noPU_genMatched_EE_vtx.at(idz)->SetLineWidth(2);

	if(idz==4) {
	  // (muon, track)
	  list_gr_roc_PU200_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_1sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_1sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_2sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_2sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_3sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_3sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_4sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_4sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_genMatched_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_genMatched_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_1sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_1sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_2sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_2sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_3sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_3sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_4sigma_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_4sigma_EE.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_genMatched_EB.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_genMatched_EE.at(idz)->SetLineColor(kMagenta);
	  // (PV, track)
	  list_gr_roc_PU200_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_1sigma_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_1sigma_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_2sigma_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_2sigma_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_3sigma_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_3sigma_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_4sigma_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_4sigma_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_genMatched_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_PU200_genMatched_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_1sigma_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_1sigma_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_2sigma_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_2sigma_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_3sigma_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_3sigma_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_4sigma_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_4sigma_EE_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_genMatched_EB_vtx.at(idz)->SetLineColor(kMagenta);
	  list_gr_roc_noPU_genMatched_EE_vtx.at(idz)->SetLineColor(kMagenta);

	}
  }
  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
    list_gr_roc_PU200_EB.at(idz)->GetXaxis()->SetTitle("Prompt efficiency"); list_gr_roc_PU200_EB.at(idz)->GetYaxis()->SetTitle("Non-prompt efficiency");
    list_gr_roc_PU200_EE.at(idz)->GetXaxis()->SetTitle("Prompt efficiency"); list_gr_roc_PU200_EE.at(idz)->GetYaxis()->SetTitle("Non-prompt efficiency");
    list_gr_roc_noPU_EB.at(idz)->GetXaxis()->SetTitle("Prompt efficiency"); list_gr_roc_noPU_EB.at(idz)->GetYaxis()->SetTitle("Non-prompt efficiency");
    list_gr_roc_noPU_EE.at(idz)->GetXaxis()->SetTitle("Prompt efficiency"); list_gr_roc_noPU_EE.at(idz)->GetYaxis()->SetTitle("Non-prompt efficiency");
  }


  ///////////////
  //  Legends  //
  ///////////////
  TLegend* leg_roc_PU200_EB = new TLegend(0.15, 0.61, 0.63, 0.88); TLegend* leg_roc_PU200_EE = new TLegend(0.15, 0.61, 0.63, 0.88);
  TLegend* leg_roc_noPU_EB = new TLegend(0.15, 0.61, 0.63, 0.88); TLegend* leg_roc_noPU_EE = new TLegend(0.15, 0.61, 0.63, 0.88);
  for(int idz=0; idz<track_pv_dz_cut_EB.size(); idz++) {
	leg_roc_PU200_EB->AddEntry(list_gr_roc_PU200_EB.at(idz), Form("no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(idz)*10), "L");
	leg_roc_PU200_EB->SetTextSize(0.03);
	leg_roc_PU200_EE->AddEntry(list_gr_roc_PU200_EE.at(idz), Form("no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(idz)*10), "L");
	leg_roc_PU200_EE->SetTextSize(0.03);
	leg_roc_noPU_EB->AddEntry(list_gr_roc_noPU_EB.at(idz), Form("no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(idz)*10), "L");
	leg_roc_noPU_EB->SetTextSize(0.03);
	leg_roc_noPU_EE->AddEntry(list_gr_roc_noPU_EE.at(idz), Form("no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(idz)*10), "L");
	leg_roc_noPU_EE->SetTextSize(0.03);
  }

  /////////////////
  ////// ROC //////
  /////////////////

  ////////////////
  //    Draw    //
  ////////////////
  TCanvas* c_PU200_EB_reliso_roc = new TCanvas("c_PU200_EB_reliso_roc", "c_PU200_EB_reliso_roc", 1200, 1200);
  c_PU200_EB_reliso_roc->cd();
  c_PU200_EB_reliso_roc->SetGrid();
  c_PU200_EB_reliso_roc->SetLeftMargin(0.12);
  //list_gr_roc_PU200_EB.at(0)->GetXaxis()->SetRangeUser(0.,1.);
  //list_gr_roc_PU200_EB.at(0)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_roc_PU200_EB.at(0)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EB.at(0)->GetYaxis()->SetRangeUser(0.,0.4);
  list_gr_roc_PU200_EB.at(0)->Draw("AL");
  for(int idz=1; idz<track_pv_dz_cut_EB.size(); idz++) {
	list_gr_roc_PU200_EB.at(idz)->Draw("same");
  }
//  gr_roc_noPU_EB->Draw("same");
  leg_roc_PU200_EB->Draw();
  c_PU200_EB_reliso_roc->Print("plots/ntuple/dz_study/isoroc_PU200_EB_dzcut.pdf");

  TCanvas* c_PU200_EE_reliso_roc = new TCanvas("c_PU200_EE_reliso_roc", "c_PU200_EE_reliso_roc", 1200, 1200);
  c_PU200_EE_reliso_roc->cd();
  c_PU200_EE_reliso_roc->SetGrid();
  c_PU200_EE_reliso_roc->SetLeftMargin(0.12);
  //list_gr_roc_PU200_EE.at(0)->GetXaxis()->SetRangeUser(0.,1.);
  //list_gr_roc_PU200_EE.at(0)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_roc_PU200_EE.at(0)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EE.at(0)->GetYaxis()->SetRangeUser(0.,0.4);
  list_gr_roc_PU200_EE.at(0)->Draw("AL");
  for(int idz=1; idz<track_pv_dz_cut_EE.size(); idz++) {
	list_gr_roc_PU200_EE.at(idz)->Draw("same");
  }
//  gr_roc_noPU_EE->Draw("same");
  leg_roc_PU200_EE->Draw();
  c_PU200_EE_reliso_roc->Print("plots/ntuple/dz_study/isoroc_PU200_EE_dzcut.pdf");

  TCanvas* c_noPU_EB_reliso_roc = new TCanvas("c_noPU_EB_reliso_roc", "c_noPU_EB_reliso_roc", 1200, 1200);
  c_noPU_EB_reliso_roc->cd();
  c_noPU_EB_reliso_roc->SetGrid();
  c_noPU_EB_reliso_roc->SetLeftMargin(0.12);
  //list_gr_roc_noPU_EB.at(0)->GetXaxis()->SetRangeUser(0.,1.);
  //list_gr_roc_noPU_EB.at(0)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_roc_noPU_EB.at(0)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_noPU_EB.at(0)->GetYaxis()->SetRangeUser(0.,0.4);
  list_gr_roc_noPU_EB.at(0)->Draw("AL");
  for(int idz=1; idz<track_pv_dz_cut_EB.size(); idz++) {
	list_gr_roc_noPU_EB.at(idz)->Draw("same");
  }
//  gr_roc_noPU_EB->Draw("same");
  leg_roc_noPU_EB->Draw();
  c_noPU_EB_reliso_roc->Print("plots/ntuple/dz_study/isoroc_noPU_EB_dzcut.pdf");

  TCanvas* c_noPU_EE_reliso_roc = new TCanvas("c_noPU_EE_reliso_roc", "c_noPU_EE_reliso_roc", 1200, 1200);
  c_noPU_EE_reliso_roc->cd();
  c_noPU_EE_reliso_roc->SetGrid();
  c_noPU_EE_reliso_roc->SetLeftMargin(0.12);
  //list_gr_roc_noPU_EE.at(0)->GetXaxis()->SetRangeUser(0.,1.);
  //list_gr_roc_noPU_EE.at(0)->GetYaxis()->SetRangeUser(0.,1.);
  list_gr_roc_noPU_EE.at(0)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_noPU_EE.at(0)->GetYaxis()->SetRangeUser(0.,0.4);
  list_gr_roc_noPU_EE.at(0)->Draw("AL");
  for(int idz=1; idz<track_pv_dz_cut_EE.size(); idz++) {
	list_gr_roc_noPU_EE.at(idz)->Draw("same");
  }
//  gr_roc_noPU_EE->Draw("same");
  leg_roc_noPU_EE->Draw();
  c_noPU_EE_reliso_roc->Print("plots/ntuple/dz_study/isoroc_noPU_EE_dzcut.pdf");



  // test
    // EB
  // dz1
  TCanvas* c_EB_reliso_roc_dz1 = new TCanvas("c_EB_reliso_roc_dz1", "c_EB_reliso_roc_dz1", 1200, 1200);
  TLegend* leg_roc_EB_dz1 = new TLegend(0.15, 0.52, 0.63, 0.88);
  leg_roc_EB_dz1->AddEntry(list_gr_roc_noPU_genMatched_EB.at(0), "noPU no MTD gen", "L");
  leg_roc_EB_dz1->AddEntry(list_gr_roc_noPU_EB.at(0), "noPU no MTD", "L");
  leg_roc_EB_dz1->AddEntry(list_gr_roc_PU200_genMatched_EB.at(0), "PU200 no MTD gen", "L");
  leg_roc_EB_dz1->AddEntry(list_gr_roc_PU200_EB.at(0), "PU200 no MTD", "L");
  leg_roc_EB_dz1->AddEntry(list_gr_roc_PU200_2sigma_EB_vtx.at(0), "PU200 2sigma (PV, track)", "L");
  leg_roc_EB_dz1->AddEntry(list_gr_roc_PU200_3sigma_EB_vtx.at(0), "PU200 3sigma (PV, track)", "L");
  leg_roc_EB_dz1->AddEntry(list_gr_roc_PU200_2sigma_EB.at(0), "PU200 2sigma (muon, track)", "L");
  leg_roc_EB_dz1->AddEntry(list_gr_roc_PU200_3sigma_EB.at(0), "PU200 3sigma (muon, track)", "L");
//  leg_roc_EB_dz1->AddEntry(list_gr_roc_noPU_2sigma_EB.at(0), "noPU 2sigma", "L");
  list_gr_roc_PU200_EB.at(0)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EB.at(0)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EB_vtx.at(0)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EB_vtx.at(0)->SetLineStyle(7);
  list_gr_roc_PU200_3sigma_EB.at(0)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EB_vtx.at(0)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EB_vtx.at(0)->SetLineStyle(7);
  list_gr_roc_PU200_genMatched_EB.at(0)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EB.at(0)->SetLineStyle(7);
  list_gr_roc_noPU_EB.at(0)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(0)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(0)->SetLineStyle(7);
  c_EB_reliso_roc_dz1->cd();
  c_EB_reliso_roc_dz1->SetGrid();
  c_EB_reliso_roc_dz1->SetLeftMargin(0.12);
  list_gr_roc_PU200_EB.at(0)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EB.at(0)->GetYaxis()->SetRangeUser(0., 0.4);
  list_gr_roc_PU200_EB.at(0)->Draw("AL");
  list_gr_roc_PU200_2sigma_EB.at(0)->Draw("same");
  list_gr_roc_PU200_2sigma_EB_vtx.at(0)->Draw("same");
  list_gr_roc_PU200_3sigma_EB.at(0)->Draw("same");
  list_gr_roc_PU200_3sigma_EB_vtx.at(0)->Draw("same");
  list_gr_roc_PU200_genMatched_EB.at(0)->Draw("same");
  list_gr_roc_noPU_EB.at(0)->Draw("same");
  list_gr_roc_noPU_genMatched_EB.at(0)->Draw("same");
  list_gr_roc_PU200_EB.at(0)->Draw("same");
  leg_roc_EB_dz1->Draw("same");
  c_EB_reliso_roc_dz1->Print("plots/ntuple/dz_study/isoroc_EB_dz1.pdf");

  // dz2
  TCanvas* c_EB_reliso_roc_dz2 = new TCanvas("c_EB_reliso_roc_dz2", "c_EB_reliso_roc_dz2", 1200, 1200);
  TLegend* leg_roc_EB_dz2 = new TLegend(0.15, 0.52, 0.63, 0.88);
  leg_roc_EB_dz2->AddEntry(list_gr_roc_noPU_genMatched_EB.at(1), "noPU no MTD gen", "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_noPU_EB.at(1), "noPU no MTD", "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_PU200_genMatched_EB.at(1), "PU200 no MTD gen", "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_PU200_EB.at(1), "PU200 no MTD", "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_PU200_2sigma_EB_vtx.at(1), "PU200 2sigma (PV, track)", "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_PU200_3sigma_EB_vtx.at(1), "PU200 3sigma (PV, track)", "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_PU200_2sigma_EB.at(1), "PU200 2sigma (muon, track)", "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_PU200_3sigma_EB.at(1), "PU200 3sigma (muon, track)", "L");
//  leg_roc_EB_dz2->AddEntry(list_gr_roc_noPU_2sigma_EB.at(1), "noPU 2sigma", "L");
  list_gr_roc_PU200_EB.at(1)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EB.at(1)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EB_vtx.at(1)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EB_vtx.at(1)->SetLineStyle(7);
  list_gr_roc_PU200_3sigma_EB.at(1)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EB_vtx.at(1)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EB_vtx.at(1)->SetLineStyle(7);
  list_gr_roc_PU200_genMatched_EB.at(1)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EB.at(1)->SetLineStyle(7);
  list_gr_roc_noPU_EB.at(1)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(1)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(1)->SetLineStyle(7);
  c_EB_reliso_roc_dz2->cd();
  c_EB_reliso_roc_dz2->SetGrid();
  c_EB_reliso_roc_dz2->SetLeftMargin(0.12);
  list_gr_roc_PU200_EB.at(1)->GetXaxis()->SetRangeUser(0.85,1.);
  //list_gr_roc_PU200_EB.at(1)->GetYaxis()->SetRangeUser(0., 0.4);
  list_gr_roc_PU200_EB.at(1)->GetYaxis()->SetRangeUser(0., 1.);
  list_gr_roc_PU200_EB.at(1)->Draw("AL");
  list_gr_roc_PU200_2sigma_EB.at(1)->Draw("same");
  list_gr_roc_PU200_2sigma_EB_vtx.at(1)->Draw("same");
  list_gr_roc_PU200_3sigma_EB.at(1)->Draw("same");
  list_gr_roc_PU200_3sigma_EB_vtx.at(1)->Draw("same");
  list_gr_roc_PU200_genMatched_EB.at(1)->Draw("same");
  list_gr_roc_noPU_EB.at(1)->Draw("same");
  list_gr_roc_noPU_genMatched_EB.at(1)->Draw("same");
  list_gr_roc_PU200_EB.at(1)->Draw("same");
  leg_roc_EB_dz2->Draw("same");
  c_EB_reliso_roc_dz2->Print("plots/ntuple/dz_study/isoroc_EB_dz2.pdf");

  // dz3
  TCanvas* c_EB_reliso_roc_dz3 = new TCanvas("c_EB_reliso_roc_dz3", "c_EB_reliso_roc_dz3", 1200, 1200);
  TLegend* leg_roc_EB_dz3 = new TLegend(0.15, 0.52, 0.63, 0.88);
  leg_roc_EB_dz3->AddEntry(list_gr_roc_noPU_genMatched_EB.at(2), "noPU no MTD gen", "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_noPU_EB.at(2), "noPU no MTD", "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_PU200_genMatched_EB.at(2), "PU200 no MTD gen", "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_PU200_EB.at(2), "PU200 no MTD", "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_PU200_2sigma_EB_vtx.at(2), "PU200 2sigma (PV, track)", "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_PU200_3sigma_EB_vtx.at(2), "PU200 3sigma (PV, track)", "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_PU200_2sigma_EB.at(2), "PU200 2sigma (muon, track)", "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_PU200_3sigma_EB.at(2), "PU200 3sigma (muon, track)", "L");
//  leg_roc_EB_dz3->AddEntry(list_gr_roc_noPU_2sigma_EB.at(2), "noPU 2sigma", "L");
  list_gr_roc_PU200_EB.at(2)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EB.at(2)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EB_vtx.at(2)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EB_vtx.at(2)->SetLineStyle(7);
  list_gr_roc_PU200_3sigma_EB.at(2)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EB_vtx.at(2)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EB_vtx.at(2)->SetLineStyle(7);
  list_gr_roc_PU200_genMatched_EB.at(2)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EB.at(2)->SetLineStyle(7);
  list_gr_roc_noPU_EB.at(2)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(2)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(2)->SetLineStyle(7);
  c_EB_reliso_roc_dz3->cd();
  c_EB_reliso_roc_dz3->SetGrid();
  c_EB_reliso_roc_dz3->SetLeftMargin(0.12);
  list_gr_roc_PU200_EB.at(2)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EB.at(2)->GetYaxis()->SetRangeUser(0., 0.4);
  list_gr_roc_PU200_EB.at(2)->Draw("AL");
  list_gr_roc_PU200_2sigma_EB.at(2)->Draw("same");
  list_gr_roc_PU200_2sigma_EB_vtx.at(2)->Draw("same");
  list_gr_roc_PU200_3sigma_EB.at(2)->Draw("same");
  list_gr_roc_PU200_3sigma_EB_vtx.at(2)->Draw("same");
  list_gr_roc_PU200_genMatched_EB.at(2)->Draw("same");
  list_gr_roc_noPU_EB.at(2)->Draw("same");
  list_gr_roc_noPU_genMatched_EB.at(2)->Draw("same");
  list_gr_roc_PU200_EB.at(2)->Draw("same");
  leg_roc_EB_dz3->Draw("same");
  c_EB_reliso_roc_dz3->Print("plots/ntuple/dz_study/isoroc_EB_dz3.pdf");

  // dz4
  TCanvas* c_EB_reliso_roc_dz4 = new TCanvas("c_EB_reliso_roc_dz4", "c_EB_reliso_roc_dz4", 1200, 1200);
  TLegend* leg_roc_EB_dz4 = new TLegend(0.15, 0.52, 0.63, 0.88);
  leg_roc_EB_dz4->AddEntry(list_gr_roc_noPU_genMatched_EB.at(3), "noPU no MTD gen", "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_noPU_EB.at(3), "noPU no MTD", "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_PU200_genMatched_EB.at(3), "PU200 no MTD gen", "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_PU200_EB.at(3), "PU200 no MTD", "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_PU200_2sigma_EB_vtx.at(3), "PU200 2sigma (PV, track)", "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_PU200_3sigma_EB_vtx.at(3), "PU200 3sigma (PV, track)", "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_PU200_2sigma_EB.at(3), "PU200 2sigma (muon, track)", "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_PU200_3sigma_EB.at(3), "PU200 3sigma (muon, track)", "L");
//  leg_roc_EB_dz4->AddEntry(list_gr_roc_noPU_2sigma_EB.at(3), "noPU 2sigma", "L");
  list_gr_roc_PU200_EB.at(3)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EB.at(3)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EB_vtx.at(3)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EB_vtx.at(3)->SetLineStyle(7);
  list_gr_roc_PU200_3sigma_EB.at(3)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EB_vtx.at(3)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EB_vtx.at(3)->SetLineStyle(7);
  list_gr_roc_PU200_genMatched_EB.at(3)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EB.at(3)->SetLineStyle(7);
  list_gr_roc_noPU_EB.at(3)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(3)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(3)->SetLineStyle(7);
  c_EB_reliso_roc_dz4->cd();
  c_EB_reliso_roc_dz4->SetGrid();
  c_EB_reliso_roc_dz4->SetLeftMargin(0.12);
  list_gr_roc_PU200_EB.at(3)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EB.at(3)->GetYaxis()->SetRangeUser(0., 0.4);
  list_gr_roc_PU200_EB.at(3)->Draw("AL");
  list_gr_roc_PU200_2sigma_EB.at(3)->Draw("same");
  list_gr_roc_PU200_2sigma_EB_vtx.at(3)->Draw("same");
  list_gr_roc_PU200_3sigma_EB.at(3)->Draw("same");
  list_gr_roc_PU200_3sigma_EB_vtx.at(3)->Draw("same");
  list_gr_roc_PU200_genMatched_EB.at(3)->Draw("same");
  list_gr_roc_noPU_EB.at(3)->Draw("same");
  list_gr_roc_noPU_genMatched_EB.at(3)->Draw("same");
  list_gr_roc_PU200_EB.at(3)->Draw("same");
  leg_roc_EB_dz4->Draw("same");
  c_EB_reliso_roc_dz4->Print("plots/ntuple/dz_study/isoroc_EB_dz4.pdf");

  // dz5
  TCanvas* c_EB_reliso_roc_dz5 = new TCanvas("c_EB_reliso_roc_dz5", "c_EB_reliso_roc_dz5", 1200, 1200);
  TLegend* leg_roc_EB_dz5 = new TLegend(0.15, 0.52, 0.63, 0.88);
  leg_roc_EB_dz5->AddEntry(list_gr_roc_noPU_genMatched_EB.at(4), "noPU no MTD gen", "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_noPU_EB.at(4), "noPU no MTD", "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_PU200_genMatched_EB.at(4), "PU200 no MTD gen", "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_PU200_EB.at(4), "PU200 no MTD", "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_PU200_2sigma_EB_vtx.at(4), "PU200 2sigma (PV, track)", "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_PU200_3sigma_EB_vtx.at(4), "PU200 3sigma (PV, track)", "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_PU200_2sigma_EB.at(4), "PU200 2sigma (muon, track)", "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_PU200_3sigma_EB.at(4), "PU200 3sigma (muon, track)", "L");
//  leg_roc_EB_dz5->AddEntry(list_gr_roc_noPU_2sigma_EB.at(4), "noPU 2sigma", "L");
  list_gr_roc_PU200_EB.at(4)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EB.at(4)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EB_vtx.at(4)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EB_vtx.at(4)->SetLineStyle(7);
  list_gr_roc_PU200_3sigma_EB.at(4)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EB_vtx.at(4)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EB_vtx.at(4)->SetLineStyle(7);
  list_gr_roc_PU200_genMatched_EB.at(4)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EB.at(4)->SetLineStyle(7);
  list_gr_roc_noPU_EB.at(4)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(4)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(4)->SetLineStyle(7);
  c_EB_reliso_roc_dz5->cd();
  c_EB_reliso_roc_dz5->SetGrid();
  c_EB_reliso_roc_dz5->SetLeftMargin(0.12);
  list_gr_roc_PU200_EB.at(4)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EB.at(4)->GetYaxis()->SetRangeUser(0., 0.4);
  list_gr_roc_PU200_EB.at(4)->Draw("AL");
  list_gr_roc_PU200_2sigma_EB.at(4)->Draw("same");
  list_gr_roc_PU200_2sigma_EB_vtx.at(4)->Draw("same");
  list_gr_roc_PU200_3sigma_EB.at(4)->Draw("same");
  list_gr_roc_PU200_3sigma_EB_vtx.at(4)->Draw("same");
  list_gr_roc_PU200_genMatched_EB.at(4)->Draw("same");
  list_gr_roc_noPU_EB.at(4)->Draw("same");
  list_gr_roc_noPU_genMatched_EB.at(4)->Draw("same");
  list_gr_roc_PU200_EB.at(4)->Draw("same");
  leg_roc_EB_dz5->Draw("same");
  c_EB_reliso_roc_dz5->Print("plots/ntuple/dz_study/isoroc_EB_dz5.pdf");

    // EE
  // dz1
  TCanvas* c_EE_reliso_roc_dz1 = new TCanvas("c_EE_reliso_roc_dz1", "c_EE_reliso_roc_dz1", 1200, 1200);
  TLegend* leg_roc_EE_dz1 = new TLegend(0.15, 0.52, 0.63, 0.88);
  leg_roc_EE_dz1->AddEntry(list_gr_roc_noPU_genMatched_EE.at(0), "noPU no MTD gen", "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_noPU_EE.at(0), "noPU no MTD", "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_PU200_genMatched_EE.at(0), "PU200 no MTD gen", "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_PU200_EE.at(0), "PU200 no MTD", "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_PU200_2sigma_EE_vtx.at(0), "PU200 2sigma (PV, track)", "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_PU200_3sigma_EE_vtx.at(0), "PU200 3sigma (PV, track)", "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_PU200_2sigma_EE.at(0), "PU200 2sigma (muon, track)", "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_PU200_3sigma_EE.at(0), "PU200 3sigma (muon, track)", "L");
//  leg_roc_EE_dz1->AddEntry(list_gr_roc_noPU_2sigma_EE.at(0), "noPU 2sigma", "L");
  list_gr_roc_PU200_EE.at(0)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EE.at(0)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EE_vtx.at(0)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EE_vtx.at(0)->SetLineStyle(7);
  list_gr_roc_PU200_3sigma_EE.at(0)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EE_vtx.at(0)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EE_vtx.at(0)->SetLineStyle(7);
  list_gr_roc_PU200_genMatched_EE.at(0)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EE.at(0)->SetLineStyle(7);
  list_gr_roc_noPU_EE.at(0)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(0)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(0)->SetLineStyle(7);
  c_EE_reliso_roc_dz1->cd();
  c_EE_reliso_roc_dz1->SetGrid();
  c_EE_reliso_roc_dz1->SetLeftMargin(0.12);
  list_gr_roc_PU200_EE.at(0)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EE.at(0)->GetYaxis()->SetRangeUser(0., 0.4);
  list_gr_roc_PU200_EE.at(0)->Draw("AL");
  list_gr_roc_PU200_2sigma_EE.at(0)->Draw("same");
  list_gr_roc_PU200_2sigma_EE_vtx.at(0)->Draw("same");
  list_gr_roc_PU200_3sigma_EE.at(0)->Draw("same");
  list_gr_roc_PU200_3sigma_EE_vtx.at(0)->Draw("same");
  list_gr_roc_PU200_genMatched_EE.at(0)->Draw("same");
  list_gr_roc_noPU_EE.at(0)->Draw("same");
  list_gr_roc_noPU_genMatched_EE.at(0)->Draw("same");
  list_gr_roc_PU200_EE.at(0)->Draw("same");
  leg_roc_EE_dz1->Draw("same");
  c_EE_reliso_roc_dz1->Print("plots/ntuple/dz_study/isoroc_EE_dz1.pdf");

  // dz2
  TCanvas* c_EE_reliso_roc_dz2 = new TCanvas("c_EE_reliso_roc_dz2", "c_EE_reliso_roc_dz2", 1200, 1200);
  TLegend* leg_roc_EE_dz2 = new TLegend(0.15, 0.52, 0.63, 0.88);
  leg_roc_EE_dz2->AddEntry(list_gr_roc_noPU_genMatched_EE.at(1), "noPU no MTD gen", "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_noPU_EE.at(1), "noPU no MTD", "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_PU200_genMatched_EE.at(1), "PU200 no MTD gen", "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_PU200_EE.at(1), "PU200 no MTD", "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_PU200_2sigma_EE_vtx.at(1), "PU200 2sigma (PV, track)", "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_PU200_3sigma_EE_vtx.at(1), "PU200 3sigma (PV, track)", "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_PU200_2sigma_EE.at(1), "PU200 2sigma (muon, track)", "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_PU200_3sigma_EE.at(1), "PU200 3sigma (muon, track)", "L");
//  leg_roc_EE_dz2->AddEntry(list_gr_roc_noPU_2sigma_EE.at(1), "noPU 2sigma", "L");
  list_gr_roc_PU200_EE.at(1)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EE.at(1)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EE_vtx.at(1)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EE_vtx.at(1)->SetLineStyle(7);
  list_gr_roc_PU200_3sigma_EE.at(1)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EE_vtx.at(1)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EE_vtx.at(1)->SetLineStyle(7);
  list_gr_roc_PU200_genMatched_EE.at(1)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EE.at(1)->SetLineStyle(7);
  list_gr_roc_noPU_EE.at(1)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(1)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(1)->SetLineStyle(7);
  c_EE_reliso_roc_dz2->cd();
  c_EE_reliso_roc_dz2->SetGrid();
  c_EE_reliso_roc_dz2->SetLeftMargin(0.12);
  list_gr_roc_PU200_EE.at(1)->GetXaxis()->SetRangeUser(0.85,1.);
  //list_gr_roc_PU200_EE.at(1)->GetYaxis()->SetRangeUser(0., 0.4);
  list_gr_roc_PU200_EE.at(1)->GetYaxis()->SetRangeUser(0., 1.);
  list_gr_roc_PU200_EE.at(1)->Draw("AL");
  list_gr_roc_PU200_2sigma_EE.at(1)->Draw("same");
  list_gr_roc_PU200_2sigma_EE_vtx.at(1)->Draw("same");
  list_gr_roc_PU200_3sigma_EE.at(1)->Draw("same");
  list_gr_roc_PU200_3sigma_EE_vtx.at(1)->Draw("same");
  list_gr_roc_PU200_genMatched_EE.at(1)->Draw("same");
  list_gr_roc_noPU_EE.at(1)->Draw("same");
  list_gr_roc_noPU_genMatched_EE.at(1)->Draw("same");
  list_gr_roc_PU200_EE.at(1)->Draw("same");
  leg_roc_EE_dz2->Draw("same");
  c_EE_reliso_roc_dz2->Print("plots/ntuple/dz_study/isoroc_EE_dz2.pdf");

  // dz3
  TCanvas* c_EE_reliso_roc_dz3 = new TCanvas("c_EE_reliso_roc_dz3", "c_EE_reliso_roc_dz3", 1200, 1200);
  TLegend* leg_roc_EE_dz3 = new TLegend(0.15, 0.52, 0.63, 0.88);
  leg_roc_EE_dz3->AddEntry(list_gr_roc_noPU_genMatched_EE.at(2), "noPU no MTD gen", "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_noPU_EE.at(2), "noPU no MTD", "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_PU200_genMatched_EE.at(2), "PU200 no MTD gen", "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_PU200_EE.at(2), "PU200 no MTD", "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_PU200_2sigma_EE_vtx.at(2), "PU200 2sigma (PV, track)", "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_PU200_3sigma_EE_vtx.at(2), "PU200 3sigma (PV, track)", "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_PU200_2sigma_EE.at(2), "PU200 2sigma (muon, track)", "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_PU200_3sigma_EE.at(2), "PU200 3sigma (muon, track)", "L");
//  leg_roc_EE_dz3->AddEntry(list_gr_roc_noPU_2sigma_EE.at(2), "noPU 2sigma", "L");
  list_gr_roc_PU200_EE.at(2)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EE.at(2)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EE_vtx.at(2)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EE_vtx.at(2)->SetLineStyle(7);
  list_gr_roc_PU200_3sigma_EE.at(2)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EE_vtx.at(2)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EE_vtx.at(2)->SetLineStyle(7);
  list_gr_roc_PU200_genMatched_EE.at(2)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EE.at(2)->SetLineStyle(7);
  list_gr_roc_noPU_EE.at(2)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(2)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(2)->SetLineStyle(7);
  c_EE_reliso_roc_dz3->cd();
  c_EE_reliso_roc_dz3->SetGrid();
  c_EE_reliso_roc_dz3->SetLeftMargin(0.12);
  list_gr_roc_PU200_EE.at(2)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EE.at(2)->GetYaxis()->SetRangeUser(0., 0.4);
  list_gr_roc_PU200_EE.at(2)->Draw("AL");
  list_gr_roc_PU200_2sigma_EE.at(2)->Draw("same");
  list_gr_roc_PU200_2sigma_EE_vtx.at(2)->Draw("same");
  list_gr_roc_PU200_3sigma_EE.at(2)->Draw("same");
  list_gr_roc_PU200_3sigma_EE_vtx.at(2)->Draw("same");
  list_gr_roc_PU200_genMatched_EE.at(2)->Draw("same");
  list_gr_roc_noPU_EE.at(2)->Draw("same");
  list_gr_roc_noPU_genMatched_EE.at(2)->Draw("same");
  list_gr_roc_PU200_EE.at(2)->Draw("same");
  leg_roc_EE_dz3->Draw("same");
  c_EE_reliso_roc_dz3->Print("plots/ntuple/dz_study/isoroc_EE_dz3.pdf");

  // dz4
  TCanvas* c_EE_reliso_roc_dz4 = new TCanvas("c_EE_reliso_roc_dz4", "c_EE_reliso_roc_dz4", 1200, 1200);
  TLegend* leg_roc_EE_dz4 = new TLegend(0.15, 0.52, 0.63, 0.88);
  leg_roc_EE_dz4->AddEntry(list_gr_roc_noPU_genMatched_EE.at(3), "noPU no MTD gen", "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_noPU_EE.at(3), "noPU no MTD", "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_PU200_genMatched_EE.at(3), "PU200 no MTD gen", "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_PU200_EE.at(3), "PU200 no MTD", "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_PU200_2sigma_EE_vtx.at(3), "PU200 2sigma (PV, track)", "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_PU200_3sigma_EE_vtx.at(3), "PU200 3sigma (PV, track)", "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_PU200_2sigma_EE.at(3), "PU200 2sigma (muon, track)", "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_PU200_3sigma_EE.at(3), "PU200 3sigma (muon, track)", "L");
//  leg_roc_EE_dz4->AddEntry(list_gr_roc_noPU_2sigma_EE.at(3), "noPU 2sigma", "L");
  list_gr_roc_PU200_EE.at(3)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EE.at(3)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EE_vtx.at(3)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EE_vtx.at(3)->SetLineStyle(7);
  list_gr_roc_PU200_3sigma_EE.at(3)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EE_vtx.at(3)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EE_vtx.at(3)->SetLineStyle(7);
  list_gr_roc_PU200_genMatched_EE.at(3)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EE.at(3)->SetLineStyle(7);
  list_gr_roc_noPU_EE.at(3)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(3)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(3)->SetLineStyle(7);
  c_EE_reliso_roc_dz4->cd();
  c_EE_reliso_roc_dz4->SetGrid();
  c_EE_reliso_roc_dz4->SetLeftMargin(0.12);
  list_gr_roc_PU200_EE.at(3)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EE.at(3)->GetYaxis()->SetRangeUser(0., 0.4);
  list_gr_roc_PU200_EE.at(3)->Draw("AL");
  list_gr_roc_PU200_2sigma_EE.at(3)->Draw("same");
  list_gr_roc_PU200_2sigma_EE_vtx.at(3)->Draw("same");
  list_gr_roc_PU200_3sigma_EE.at(3)->Draw("same");
  list_gr_roc_PU200_3sigma_EE_vtx.at(3)->Draw("same");
  list_gr_roc_PU200_genMatched_EE.at(3)->Draw("same");
  list_gr_roc_noPU_EE.at(3)->Draw("same");
  list_gr_roc_noPU_genMatched_EE.at(3)->Draw("same");
  list_gr_roc_PU200_EE.at(3)->Draw("same");
  leg_roc_EE_dz4->Draw("same");
  c_EE_reliso_roc_dz4->Print("plots/ntuple/dz_study/isoroc_EE_dz4.pdf");

  // dz5
  TCanvas* c_EE_reliso_roc_dz5 = new TCanvas("c_EE_reliso_roc_dz5", "c_EE_reliso_roc_dz5", 1200, 1200);
  TLegend* leg_roc_EE_dz5 = new TLegend(0.15, 0.52, 0.63, 0.88);
  leg_roc_EE_dz5->AddEntry(list_gr_roc_noPU_genMatched_EE.at(4), "noPU no MTD gen", "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_noPU_EE.at(4), "noPU no MTD", "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_PU200_genMatched_EE.at(4), "PU200 no MTD gen", "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_PU200_EE.at(4), "PU200 no MTD", "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_PU200_2sigma_EE_vtx.at(4), "PU200 2sigma (PV, track)", "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_PU200_3sigma_EE_vtx.at(4), "PU200 3sigma (PV, track)", "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_PU200_2sigma_EE.at(4), "PU200 2sigma (muon, track)", "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_PU200_3sigma_EE.at(4), "PU200 3sigma (muon, track)", "L");
//  leg_roc_EE_dz5->AddEntry(list_gr_roc_noPU_2sigma_EE.at(4), "noPU 2sigma", "L");
  list_gr_roc_PU200_EE.at(4)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EE.at(4)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EE_vtx.at(4)->SetLineColor(kRed);
  list_gr_roc_PU200_2sigma_EE_vtx.at(4)->SetLineStyle(7);
  list_gr_roc_PU200_3sigma_EE.at(4)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EE_vtx.at(4)->SetLineColor(kAzure+1);
  list_gr_roc_PU200_3sigma_EE_vtx.at(4)->SetLineStyle(7);
  list_gr_roc_PU200_genMatched_EE.at(4)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EE.at(4)->SetLineStyle(7);
  list_gr_roc_noPU_EE.at(4)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(4)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(4)->SetLineStyle(7);
  c_EE_reliso_roc_dz5->cd();
  c_EE_reliso_roc_dz5->SetGrid();
  c_EE_reliso_roc_dz5->SetLeftMargin(0.12);
  list_gr_roc_PU200_EE.at(4)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EE.at(4)->GetYaxis()->SetRangeUser(0., 0.4);
  list_gr_roc_PU200_EE.at(4)->Draw("AL");
  list_gr_roc_PU200_2sigma_EE.at(4)->Draw("same");
  list_gr_roc_PU200_2sigma_EE_vtx.at(4)->Draw("same");
  list_gr_roc_PU200_3sigma_EE.at(4)->Draw("same");
  list_gr_roc_PU200_3sigma_EE_vtx.at(4)->Draw("same");
  list_gr_roc_PU200_genMatched_EE.at(4)->Draw("same");
  list_gr_roc_noPU_EE.at(4)->Draw("same");
  list_gr_roc_noPU_genMatched_EE.at(4)->Draw("same");
  list_gr_roc_PU200_EE.at(4)->Draw("same");
  leg_roc_EE_dz5->Draw("same");
  c_EE_reliso_roc_dz5->Print("plots/ntuple/dz_study/isoroc_EE_dz5.pdf");














/*

  // dz2
  TCanvas* c_EB_reliso_roc_dz2 = new TCanvas("c_EB_reliso_roc_dz2", "c_EB_reliso_roc_dz2", 1200, 1200);
  TLegend* leg_roc_EB_dz2 = new TLegend(0.15, 0.61, 0.63, 0.88);
  leg_roc_EB_dz2->AddEntry(list_gr_roc_PU200_EB.at(1), Form("PU200 no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(1)*10), "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_PU200_genMatched_EB.at(1), Form("PU200 no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(1)*10), "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_PU200_2sigma_EB.at(1), Form("PU200 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(1)*10), "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_PU200_3sigma_EB.at(1), Form("PU200 3sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(1)*10), "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_PU200_4sigma_EB.at(1), Form("PU200 4sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(1)*10), "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_noPU_EB.at(1), Form("noPU no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(1)*10), "L");
  leg_roc_EB_dz2->AddEntry(list_gr_roc_noPU_genMatched_EB.at(1), Form("noPU no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(1)*10), "L");
 // leg_roc_EB_dz2->AddEntry(list_gr_roc_noPU_2sigma_EB.at(1), Form("noPU 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(1)*10), "L");
  list_gr_roc_PU200_EB.at(1)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EB.at(1)->SetLineColor(kGreen);
  list_gr_roc_PU200_3sigma_EB.at(1)->SetLineColor(kBlue);
  list_gr_roc_PU200_4sigma_EB.at(1)->SetLineColor(kMagenta);
  list_gr_roc_PU200_genMatched_EB.at(1)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EB.at(1)->SetLineStyle(7);
  list_gr_roc_noPU_EB.at(1)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(1)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(1)->SetLineStyle(7);
  list_gr_roc_noPU_2sigma_EB.at(1)->SetLineColor(kGray+1);
  list_gr_roc_noPU_2sigma_EB.at(1)->SetLineStyle(10);
  c_EB_reliso_roc_dz2->cd();
  c_EB_reliso_roc_dz2->SetGrid();
  c_EB_reliso_roc_dz2->SetLeftMargin(0.12);
  list_gr_roc_PU200_EB.at(1)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EB.at(1)->GetYaxis()->SetRangeUser(0., 1.);
  list_gr_roc_PU200_EB.at(1)->Draw("AL");
  list_gr_roc_PU200_2sigma_EB.at(1)->Draw("same");
  list_gr_roc_PU200_3sigma_EB.at(1)->Draw("same");
  list_gr_roc_PU200_4sigma_EB.at(1)->Draw("same");
  list_gr_roc_PU200_genMatched_EB.at(1)->Draw("same");
  list_gr_roc_noPU_EB.at(1)->Draw("same");
  list_gr_roc_noPU_genMatched_EB.at(1)->Draw("same");
  list_gr_roc_noPU_2sigma_EB.at(1)->Draw("same");
  leg_roc_EB_dz2->Draw("same");
  c_EB_reliso_roc_dz2->Print("plots/ntuple/dz_study/isoroc_EB_dz2.pdf");

  // dz3
  TCanvas* c_EB_reliso_roc_dz3 = new TCanvas("c_EB_reliso_roc_dz3", "c_EB_reliso_roc_dz3", 1200, 1200);
  TLegend* leg_roc_EB_dz3 = new TLegend(0.15, 0.61, 0.63, 0.88);
  leg_roc_EB_dz3->AddEntry(list_gr_roc_PU200_EB.at(2), Form("PU200 no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(2)*10), "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_PU200_genMatched_EB.at(2), Form("PU200 no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(2)*10), "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_PU200_2sigma_EB.at(2), Form("PU200 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(2)*10), "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_PU200_3sigma_EB.at(2), Form("PU200 3sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(2)*10), "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_PU200_4sigma_EB.at(2), Form("PU200 4sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(2)*10), "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_noPU_EB.at(2), Form("noPU no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(2)*10), "L");
  leg_roc_EB_dz3->AddEntry(list_gr_roc_noPU_genMatched_EB.at(2), Form("noPU no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(2)*10), "L");
//  leg_roc_EB_dz3->AddEntry(list_gr_roc_noPU_2sigma_EB.at(2), Form("noPU 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(2)*10), "L");
  list_gr_roc_PU200_EB.at(2)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EB.at(2)->SetLineColor(kGreen);
  list_gr_roc_PU200_3sigma_EB.at(2)->SetLineColor(kBlue);
  list_gr_roc_PU200_4sigma_EB.at(2)->SetLineColor(kMagenta);
  list_gr_roc_PU200_genMatched_EB.at(2)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EB.at(2)->SetLineStyle(7);
  list_gr_roc_noPU_EB.at(2)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(2)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(2)->SetLineStyle(7);
  list_gr_roc_noPU_2sigma_EB.at(2)->SetLineColor(kGray+1);
  list_gr_roc_noPU_2sigma_EB.at(2)->SetLineStyle(10);
  c_EB_reliso_roc_dz3->cd();
  c_EB_reliso_roc_dz3->SetGrid();
  c_EB_reliso_roc_dz3->SetLeftMargin(0.12);
  list_gr_roc_PU200_EB.at(2)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EB.at(2)->GetYaxis()->SetRangeUser(0., 1.);
  list_gr_roc_PU200_EB.at(2)->Draw("AL");
  list_gr_roc_PU200_2sigma_EB.at(2)->Draw("same");
  list_gr_roc_PU200_3sigma_EB.at(2)->Draw("same");
  list_gr_roc_PU200_4sigma_EB.at(2)->Draw("same");
  list_gr_roc_PU200_genMatched_EB.at(2)->Draw("same");
  list_gr_roc_noPU_EB.at(2)->Draw("same");
  list_gr_roc_noPU_genMatched_EB.at(2)->Draw("same");
  list_gr_roc_noPU_2sigma_EB.at(2)->Draw("same");
  leg_roc_EB_dz3->Draw("same");
  c_EB_reliso_roc_dz3->Print("plots/ntuple/dz_study/isoroc_EB_dz3.pdf");

  // dz4
  TCanvas* c_EB_reliso_roc_dz4 = new TCanvas("c_EB_reliso_roc_dz4", "c_EB_reliso_roc_dz4", 1200, 1200);
  TLegend* leg_roc_EB_dz4 = new TLegend(0.15, 0.61, 0.63, 0.88);
  leg_roc_EB_dz4->AddEntry(list_gr_roc_PU200_EB.at(3), Form("PU200 no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(3)*10), "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_PU200_genMatched_EB.at(3), Form("PU200 no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(3)*10), "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_PU200_2sigma_EB.at(3), Form("PU200 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(3)*10), "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_PU200_3sigma_EB.at(3), Form("PU200 3sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(3)*10), "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_PU200_4sigma_EB.at(3), Form("PU200 4sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(3)*10), "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_noPU_EB.at(3), Form("noPU no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(3)*10), "L");
  leg_roc_EB_dz4->AddEntry(list_gr_roc_noPU_genMatched_EB.at(3), Form("noPU no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(3)*10), "L");
//  leg_roc_EB_dz4->AddEntry(list_gr_roc_noPU_2sigma_EB.at(3), Form("noPU 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(3)*10), "L");
  list_gr_roc_PU200_EB.at(3)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EB.at(3)->SetLineColor(kGreen);
  list_gr_roc_PU200_3sigma_EB.at(3)->SetLineColor(kBlue);
  list_gr_roc_PU200_4sigma_EB.at(3)->SetLineColor(kMagenta);
  list_gr_roc_PU200_genMatched_EB.at(3)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EB.at(3)->SetLineStyle(7);
  list_gr_roc_noPU_EB.at(3)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(3)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(3)->SetLineStyle(7);
  list_gr_roc_noPU_2sigma_EB.at(3)->SetLineColor(kGray+1);
  list_gr_roc_noPU_2sigma_EB.at(3)->SetLineStyle(10);
  c_EB_reliso_roc_dz4->cd();
  c_EB_reliso_roc_dz4->SetGrid();
  c_EB_reliso_roc_dz4->SetLeftMargin(0.12);
  list_gr_roc_PU200_EB.at(3)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EB.at(3)->GetYaxis()->SetRangeUser(0., 1.);
  list_gr_roc_PU200_EB.at(3)->Draw("AL");
  list_gr_roc_PU200_2sigma_EB.at(3)->Draw("same");
  list_gr_roc_PU200_3sigma_EB.at(3)->Draw("same");
  list_gr_roc_PU200_4sigma_EB.at(3)->Draw("same");
  list_gr_roc_PU200_genMatched_EB.at(3)->Draw("same");
  list_gr_roc_noPU_EB.at(3)->Draw("same");
  list_gr_roc_noPU_genMatched_EB.at(3)->Draw("same");
  list_gr_roc_noPU_2sigma_EB.at(3)->Draw("same");
  leg_roc_EB_dz4->Draw("same");
  c_EB_reliso_roc_dz4->Print("plots/ntuple/dz_study/isoroc_EB_dz4.pdf");

  // dz5
  TCanvas* c_EB_reliso_roc_dz5 = new TCanvas("c_EB_reliso_roc_dz5", "c_EB_reliso_roc_dz5", 1200, 1200);
  TLegend* leg_roc_EB_dz5 = new TLegend(0.15, 0.61, 0.63, 0.88);
  leg_roc_EB_dz5->AddEntry(list_gr_roc_PU200_EB.at(4), Form("PU200 no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(4)*10), "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_PU200_genMatched_EB.at(4), Form("PU200 no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(4)*10), "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_PU200_2sigma_EB.at(4), Form("PU200 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(4)*10), "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_PU200_3sigma_EB.at(4), Form("PU200 3sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(4)*10), "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_PU200_4sigma_EB.at(4), Form("PU200 4sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(4)*10), "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_noPU_EB.at(4), Form("noPU no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(4)*10), "L");
  leg_roc_EB_dz5->AddEntry(list_gr_roc_noPU_genMatched_EB.at(4), Form("noPU no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(4)*10), "L");
//  leg_roc_EB_dz5->AddEntry(list_gr_roc_noPU_2sigma_EB.at(4), Form("noPU 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EB.at(4)*10), "L");
  list_gr_roc_PU200_EB.at(4)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EB.at(4)->SetLineColor(kGreen);
  list_gr_roc_PU200_3sigma_EB.at(4)->SetLineColor(kBlue);
  list_gr_roc_PU200_4sigma_EB.at(4)->SetLineColor(kMagenta);
  list_gr_roc_PU200_genMatched_EB.at(4)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EB.at(4)->SetLineStyle(7);
  list_gr_roc_noPU_EB.at(4)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(4)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EB.at(4)->SetLineStyle(7);
  list_gr_roc_noPU_2sigma_EB.at(4)->SetLineColor(kGray+1);
  list_gr_roc_noPU_2sigma_EB.at(4)->SetLineStyle(10);
  c_EB_reliso_roc_dz5->cd();
  c_EB_reliso_roc_dz5->SetGrid();
  c_EB_reliso_roc_dz5->SetLeftMargin(0.12);
  list_gr_roc_PU200_EB.at(4)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EB.at(4)->GetYaxis()->SetRangeUser(0., 1.);
  list_gr_roc_PU200_EB.at(4)->Draw("AL");
  list_gr_roc_PU200_2sigma_EB.at(4)->Draw("same");
  list_gr_roc_PU200_3sigma_EB.at(4)->Draw("same");
  list_gr_roc_PU200_4sigma_EB.at(4)->Draw("same");
  list_gr_roc_PU200_genMatched_EB.at(4)->Draw("same");
  list_gr_roc_noPU_EB.at(4)->Draw("same");
  list_gr_roc_noPU_genMatched_EB.at(4)->Draw("same");
  list_gr_roc_noPU_2sigma_EB.at(4)->Draw("same");
  leg_roc_EB_dz5->Draw("same");
  c_EB_reliso_roc_dz5->Print("plots/ntuple/dz_study/isoroc_EB_dz5.pdf");


    // EE
  // dz1
  TCanvas* c_EE_reliso_roc_dz1 = new TCanvas("c_EE_reliso_roc_dz1", "c_EE_reliso_roc_dz1", 1200, 1200);
  TLegend* leg_roc_EE_dz1 = new TLegend(0.15, 0.61, 0.63, 0.88);
  leg_roc_EE_dz1->AddEntry(list_gr_roc_PU200_EE.at(0), Form("PU200 no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(0)*10), "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_PU200_genMatched_EE.at(0), Form("PU200 no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(0)*10), "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_PU200_2sigma_EE.at(0), Form("PU200 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(0)*10), "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_PU200_3sigma_EE.at(0), Form("PU200 3sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(0)*10), "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_PU200_4sigma_EE.at(0), Form("PU200 4sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(0)*10), "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_noPU_EE.at(0), Form("noPU no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(0)*10), "L");
  leg_roc_EE_dz1->AddEntry(list_gr_roc_noPU_genMatched_EE.at(0), Form("noPU no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(0)*10), "L");
//  leg_roc_EE_dz1->AddEntry(list_gr_roc_noPU_2sigma_EE.at(0), Form("noPU 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(0)*10), "L");
  list_gr_roc_PU200_EE.at(0)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EE.at(0)->SetLineColor(kGreen);
  list_gr_roc_PU200_3sigma_EE.at(0)->SetLineColor(kBlue);
  list_gr_roc_PU200_4sigma_EE.at(0)->SetLineColor(kMagenta);
  list_gr_roc_PU200_genMatched_EE.at(0)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EE.at(0)->SetLineStyle(7);
  list_gr_roc_noPU_EE.at(0)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(0)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(0)->SetLineStyle(7);
  list_gr_roc_noPU_2sigma_EE.at(0)->SetLineColor(kGray+1);
  list_gr_roc_noPU_2sigma_EE.at(0)->SetLineStyle(10);
  c_EE_reliso_roc_dz1->cd();
  c_EE_reliso_roc_dz1->SetGrid();
  c_EE_reliso_roc_dz1->SetLeftMargin(0.12);
  list_gr_roc_PU200_EE.at(0)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EE.at(0)->GetYaxis()->SetRangeUser(0., 1.);
  list_gr_roc_PU200_EE.at(0)->Draw("AL");
  list_gr_roc_PU200_2sigma_EE.at(0)->Draw("same");
  list_gr_roc_PU200_3sigma_EE.at(0)->Draw("same");
  list_gr_roc_PU200_4sigma_EE.at(0)->Draw("same");
  list_gr_roc_PU200_genMatched_EE.at(0)->Draw("same");
  list_gr_roc_noPU_EE.at(0)->Draw("same");
  list_gr_roc_noPU_genMatched_EE.at(0)->Draw("same");
  list_gr_roc_noPU_2sigma_EE.at(0)->Draw("same");
  leg_roc_EE_dz1->Draw("same");
  c_EE_reliso_roc_dz1->Print("plots/ntuple/dz_study/isoroc_EE_dz1.pdf");

  // dz2
  TCanvas* c_EE_reliso_roc_dz2 = new TCanvas("c_EE_reliso_roc_dz2", "c_EE_reliso_roc_dz2", 1200, 1200);
  TLegend* leg_roc_EE_dz2 = new TLegend(0.15, 0.61, 0.63, 0.88);
  leg_roc_EE_dz2->AddEntry(list_gr_roc_PU200_EE.at(1), Form("PU200 no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(1)*10), "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_PU200_genMatched_EE.at(1), Form("PU200 no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(1)*10), "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_PU200_2sigma_EE.at(1), Form("PU200 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(1)*10), "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_PU200_3sigma_EE.at(1), Form("PU200 3sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(1)*10), "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_PU200_4sigma_EE.at(1), Form("PU200 4sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(1)*10), "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_noPU_EE.at(1), Form("noPU no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(1)*10), "L");
  leg_roc_EE_dz2->AddEntry(list_gr_roc_noPU_genMatched_EE.at(1), Form("noPU no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(1)*10), "L");
//  leg_roc_EE_dz2->AddEntry(list_gr_roc_noPU_2sigma_EE.at(1), Form("noPU 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(1)*10), "L");
  list_gr_roc_PU200_EE.at(1)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EE.at(1)->SetLineColor(kGreen);
  list_gr_roc_PU200_3sigma_EE.at(1)->SetLineColor(kBlue);
  list_gr_roc_PU200_4sigma_EE.at(1)->SetLineColor(kMagenta);
  list_gr_roc_PU200_genMatched_EE.at(1)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EE.at(1)->SetLineStyle(7);
  list_gr_roc_noPU_EE.at(1)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(1)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(1)->SetLineStyle(7);
  list_gr_roc_noPU_2sigma_EE.at(1)->SetLineColor(kGray+1);
  list_gr_roc_noPU_2sigma_EE.at(1)->SetLineStyle(10);
  c_EE_reliso_roc_dz2->cd();
  c_EE_reliso_roc_dz2->SetGrid();
  c_EE_reliso_roc_dz2->SetLeftMargin(0.12);
  list_gr_roc_PU200_EE.at(1)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EE.at(1)->GetYaxis()->SetRangeUser(0., 1.);
  list_gr_roc_PU200_EE.at(1)->Draw("AL");
  list_gr_roc_PU200_2sigma_EE.at(1)->Draw("same");
  list_gr_roc_PU200_3sigma_EE.at(1)->Draw("same");
  list_gr_roc_PU200_4sigma_EE.at(1)->Draw("same");
  list_gr_roc_PU200_genMatched_EE.at(1)->Draw("same");
  list_gr_roc_noPU_EE.at(1)->Draw("same");
  list_gr_roc_noPU_genMatched_EE.at(1)->Draw("same");
  list_gr_roc_noPU_2sigma_EE.at(1)->Draw("same");
  leg_roc_EE_dz2->Draw("same");
  c_EE_reliso_roc_dz2->Print("plots/ntuple/dz_study/isoroc_EE_dz2.pdf");

  // dz3
  TCanvas* c_EE_reliso_roc_dz3 = new TCanvas("c_EE_reliso_roc_dz3", "c_EE_reliso_roc_dz3", 1200, 1200);
  TLegend* leg_roc_EE_dz3 = new TLegend(0.15, 0.61, 0.63, 0.88);
  leg_roc_EE_dz3->AddEntry(list_gr_roc_PU200_EE.at(2), Form("PU200 no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(2)*10), "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_PU200_genMatched_EE.at(2), Form("PU200 no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(2)*10), "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_PU200_2sigma_EE.at(2), Form("PU200 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(2)*10), "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_PU200_3sigma_EE.at(2), Form("PU200 3sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(2)*10), "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_PU200_4sigma_EE.at(2), Form("PU200 4sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(2)*10), "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_noPU_EE.at(2), Form("noPU no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(2)*10), "L");
  leg_roc_EE_dz3->AddEntry(list_gr_roc_noPU_genMatched_EE.at(2), Form("noPU no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(2)*10), "L");
 // leg_roc_EE_dz3->AddEntry(list_gr_roc_noPU_2sigma_EE.at(2), Form("noPU 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(2)*10), "L");
  list_gr_roc_PU200_EE.at(2)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EE.at(2)->SetLineColor(kGreen);
  list_gr_roc_PU200_3sigma_EE.at(2)->SetLineColor(kBlue);
  list_gr_roc_PU200_4sigma_EE.at(2)->SetLineColor(kMagenta);
  list_gr_roc_PU200_genMatched_EE.at(2)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EE.at(2)->SetLineStyle(7);
  list_gr_roc_noPU_EE.at(2)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(2)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(2)->SetLineStyle(7);
  list_gr_roc_noPU_2sigma_EE.at(2)->SetLineColor(kGray+1);
  list_gr_roc_noPU_2sigma_EE.at(2)->SetLineStyle(10);
  c_EE_reliso_roc_dz3->cd();
  c_EE_reliso_roc_dz3->SetGrid();
  c_EE_reliso_roc_dz3->SetLeftMargin(0.12);
  list_gr_roc_PU200_EE.at(2)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EE.at(2)->GetYaxis()->SetRangeUser(0., 1.);
  list_gr_roc_PU200_EE.at(2)->Draw("AL");
  list_gr_roc_PU200_2sigma_EE.at(2)->Draw("same");
  list_gr_roc_PU200_3sigma_EE.at(2)->Draw("same");
  list_gr_roc_PU200_4sigma_EE.at(2)->Draw("same");
  list_gr_roc_PU200_genMatched_EE.at(2)->Draw("same");
  list_gr_roc_noPU_EE.at(2)->Draw("same");
  list_gr_roc_noPU_genMatched_EE.at(2)->Draw("same");
  list_gr_roc_noPU_2sigma_EE.at(2)->Draw("same");
  leg_roc_EE_dz3->Draw("same");
  c_EE_reliso_roc_dz3->Print("plots/ntuple/dz_study/isoroc_EE_dz3.pdf");

  // dz4
  TCanvas* c_EE_reliso_roc_dz4 = new TCanvas("c_EE_reliso_roc_dz4", "c_EE_reliso_roc_dz4", 1200, 1200);
  TLegend* leg_roc_EE_dz4 = new TLegend(0.15, 0.61, 0.63, 0.88);
  leg_roc_EE_dz4->AddEntry(list_gr_roc_PU200_EE.at(3), Form("PU200 no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(3)*10), "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_PU200_genMatched_EE.at(3), Form("PU200 no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(3)*10), "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_PU200_2sigma_EE.at(3), Form("PU200 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(3)*10), "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_PU200_3sigma_EE.at(3), Form("PU200 3sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(3)*10), "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_PU200_4sigma_EE.at(3), Form("PU200 4sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(3)*10), "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_noPU_EE.at(3), Form("noPU no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(3)*10), "L");
  leg_roc_EE_dz4->AddEntry(list_gr_roc_noPU_genMatched_EE.at(3), Form("noPU no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(3)*10), "L");
//  leg_roc_EE_dz4->AddEntry(list_gr_roc_noPU_2sigma_EE.at(3), Form("noPU 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(3)*10), "L");
  list_gr_roc_PU200_EE.at(3)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EE.at(3)->SetLineColor(kGreen);
  list_gr_roc_PU200_3sigma_EE.at(3)->SetLineColor(kBlue);
  list_gr_roc_PU200_4sigma_EE.at(3)->SetLineColor(kMagenta);
  list_gr_roc_PU200_genMatched_EE.at(3)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EE.at(3)->SetLineStyle(7);
  list_gr_roc_noPU_EE.at(3)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(3)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(3)->SetLineStyle(7);
  list_gr_roc_noPU_2sigma_EE.at(3)->SetLineColor(kGray+1);
  list_gr_roc_noPU_2sigma_EE.at(3)->SetLineStyle(10);
  c_EE_reliso_roc_dz4->cd();
  c_EE_reliso_roc_dz4->SetGrid();
  c_EE_reliso_roc_dz4->SetLeftMargin(0.12);
  list_gr_roc_PU200_EE.at(3)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EE.at(3)->GetYaxis()->SetRangeUser(0., 1.);
  list_gr_roc_PU200_EE.at(3)->Draw("AL");
  list_gr_roc_PU200_2sigma_EE.at(3)->Draw("same");
  list_gr_roc_PU200_3sigma_EE.at(3)->Draw("same");
  list_gr_roc_PU200_4sigma_EE.at(3)->Draw("same");
  list_gr_roc_PU200_genMatched_EE.at(3)->Draw("same");
  list_gr_roc_noPU_EE.at(3)->Draw("same");
  list_gr_roc_noPU_genMatched_EE.at(3)->Draw("same");
  list_gr_roc_noPU_2sigma_EE.at(3)->Draw("same");
  leg_roc_EE_dz4->Draw("same");
  c_EE_reliso_roc_dz4->Print("plots/ntuple/dz_study/isoroc_EE_dz4.pdf");

  // dz5
  TCanvas* c_EE_reliso_roc_dz5 = new TCanvas("c_EE_reliso_roc_dz5", "c_EE_reliso_roc_dz5", 1200, 1200);
  TLegend* leg_roc_EE_dz5 = new TLegend(0.15, 0.61, 0.63, 0.88);
  leg_roc_EE_dz5->AddEntry(list_gr_roc_PU200_EE.at(4), Form("PU200 no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(4)*10), "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_PU200_genMatched_EE.at(4), Form("PU200 no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(4)*10), "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_PU200_2sigma_EE.at(4), Form("PU200 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(4)*10), "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_PU200_3sigma_EE.at(4), Form("PU200 3sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(4)*10), "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_PU200_4sigma_EE.at(4), Form("PU200 4sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(4)*10), "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_noPU_EE.at(4), Form("noPU no MTD: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(4)*10), "L");
  leg_roc_EE_dz5->AddEntry(list_gr_roc_noPU_genMatched_EE.at(4), Form("noPU no MTD gen: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(4)*10), "L");
//  leg_roc_EE_dz5->AddEntry(list_gr_roc_noPU_2sigma_EE.at(4), Form("noPU 2sigma: |#Delta z| < %.1f mm", track_pv_dz_cut_EE.at(4)*10), "L");
  list_gr_roc_PU200_EE.at(4)->SetLineColor(kBlack);
  list_gr_roc_PU200_2sigma_EE.at(4)->SetLineColor(kGreen);
  list_gr_roc_PU200_3sigma_EE.at(4)->SetLineColor(kBlue);
  list_gr_roc_PU200_4sigma_EE.at(4)->SetLineColor(kMagenta);
  list_gr_roc_PU200_genMatched_EE.at(4)->SetLineColor(kBlack);
  list_gr_roc_PU200_genMatched_EE.at(4)->SetLineStyle(7);
  list_gr_roc_noPU_EE.at(4)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(4)->SetLineColor(kGray+1);
  list_gr_roc_noPU_genMatched_EE.at(4)->SetLineStyle(7);
  list_gr_roc_noPU_2sigma_EE.at(4)->SetLineColor(kGray+1);
  list_gr_roc_noPU_2sigma_EE.at(4)->SetLineStyle(10);
  c_EE_reliso_roc_dz5->cd();
  c_EE_reliso_roc_dz5->SetGrid();
  c_EE_reliso_roc_dz5->SetLeftMargin(0.12);
  list_gr_roc_PU200_EE.at(4)->GetXaxis()->SetRangeUser(0.85,1.);
  list_gr_roc_PU200_EE.at(4)->GetYaxis()->SetRangeUser(0., 1.);
  list_gr_roc_PU200_EE.at(4)->Draw("AL");
  list_gr_roc_PU200_2sigma_EE.at(4)->Draw("same");
  list_gr_roc_PU200_3sigma_EE.at(4)->Draw("same");
  list_gr_roc_PU200_4sigma_EE.at(4)->Draw("same");
  list_gr_roc_PU200_genMatched_EE.at(4)->Draw("same");
  list_gr_roc_noPU_EE.at(4)->Draw("same");
  list_gr_roc_noPU_genMatched_EE.at(4)->Draw("same");
  list_gr_roc_noPU_2sigma_EE.at(4)->Draw("same");
  leg_roc_EE_dz5->Draw("same");
  c_EE_reliso_roc_dz5->Print("plots/ntuple/dz_study/isoroc_EE_dz5.pdf");


*/



  // test end







}

int main(int argc, char **argv)
{
  
  draw_iso_efficiency_ntuple_dz();
//  draw_track_type_sigma_ntuple();

  return 0;
}



















































