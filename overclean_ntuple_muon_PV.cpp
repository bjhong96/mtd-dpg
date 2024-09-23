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
//TString ntuple_PU200_nonprompt_vtx = "240715/ntuple_PU200_nonprompt_muon_ttbar_vtx.root";
//TString ntuple_noPU_nonprompt_vtx  = "240715/ntuple_noPU_nonprompt_muon_ttbar_vtx.root";

 





void draw_track_type_sigma_ntuple() {

  // generate the dictionary for vector<vector<sth>> collection
  gInterpreter->GenerateDictionary("vector<vector<float> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<double> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<int> >", "vector");
  gInterpreter->GenerateDictionary("vector<vector<bool> >", "vector");

  gErrorIgnoreLevel = kWarning;
//  gStyle->SetOptStat(0);
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
  vector<float> *muon_pt_(0), *muon_time_(0), *muon_time_err_(0), *muon_time_err_i_(0), *muon_mva_(0), *muon_prompt_(0), *muon_isBarrel_(0);
  vector<float> *muon_status_(0), *muon_pv_dz_(0), *muon_pv_dxy_(0), *muon_vz_(0);
  vector<float> *muon_PVweight_(0);
  vector<bool>  *muon_isMuon_(0), *muon_isPFMuon_(0), *muon_isGlobalMuon_(0), *muon_isTrackerMuon_(0), *muon_isStandAloneMuon_(0), *muon_isCutBasedIdLoose_(0), *muon_isLooseMuon_(0);
  // track
  std::vector<std::vector<float>> *track_pt_(0), *track_time_(0), *track_time_err_(0), *track_time_err_i_(0), *track_mva_(0);
  std::vector<std::vector<float>> *track_pv_dz_(0), *track_vz_(0);
  std::vector<std::vector<float>> *track_PVweight_(0);
  std::vector<std::vector<int>>   *track_type_(0);
  std::vector<std::vector<int>>   *track_evtId_(0), *track_bx_(0);
  vector<int>   *muon_index_(0);
  std::vector<std::vector<int>> *track_index_(0);
  std::vector<std::vector<bool>> *track_genMatched_(0);
  std::vector<std::vector<bool>> *selectedLV_(0), *match_vtx_sim2reco_(0), *match_vtx_reco2sim_(0);
  int vtx_index_=999;
  int recovtx_sim_=999, simvtx_reco_=999;
  int simvtx_bx_=999, simvtx_evtId_=999;
  int recovtx_original_index_=999;
  vector<float> *vtx_time_(0), *vtx_time_err_(0);

  std::vector<float> *muon_time_sim_(0);
  std::vector<std::vector<float>> *track_time_sim_(0);
  float simvtx_time_=999;


  // PU200
  ch_PU200_prompt->SetBranchAddress("muon_pt_",        						&muon_pt_);
  ch_PU200_prompt->SetBranchAddress("muon_time_",      						&muon_time_);
  ch_PU200_prompt->SetBranchAddress("muon_time_err_",  						&muon_time_err_);
  ch_PU200_prompt->SetBranchAddress("muon_time_err_i_", 					&muon_time_err_i_);
  ch_PU200_prompt->SetBranchAddress("muon_mva_",     						&muon_mva_);
  ch_PU200_prompt->SetBranchAddress("muon_prompt_",   						&muon_prompt_);
  ch_PU200_prompt->SetBranchAddress("muon_isBarrel_",  						&muon_isBarrel_);
  ch_PU200_prompt->SetBranchAddress("muon_status_",    						&muon_status_);
  ch_PU200_prompt->SetBranchAddress("muon_pv_dz_",     						&muon_pv_dz_);
  ch_PU200_prompt->SetBranchAddress("muon_pv_dxy_",    						&muon_pv_dxy_);
  ch_PU200_prompt->SetBranchAddress("muon_vz_",        						&muon_vz_);
  ch_PU200_prompt->SetBranchAddress("muon_PVweight_",  						&muon_PVweight_);
  ch_PU200_prompt->SetBranchAddress("track_pt_",       						&track_pt_);
  ch_PU200_prompt->SetBranchAddress("track_time_",     						&track_time_);
  ch_PU200_prompt->SetBranchAddress("track_time_err_", 						&track_time_err_);
  ch_PU200_prompt->SetBranchAddress("track_time_err_i_", 					&track_time_err_i_);
  ch_PU200_prompt->SetBranchAddress("track_mva_",     						&track_mva_);
  ch_PU200_prompt->SetBranchAddress("track_vz_",       						&track_vz_);
  ch_PU200_prompt->SetBranchAddress("track_pv_dz_",    						&track_pv_dz_);
  ch_PU200_prompt->SetBranchAddress("track_PVweight_", 						&track_PVweight_);
  ch_PU200_prompt->SetBranchAddress("track_type_", 	   						&track_type_);
  ch_PU200_prompt->SetBranchAddress("selectedLV_",     						&selectedLV_);
  ch_PU200_prompt->SetBranchAddress("match_vtx_sim2reco_",     				&match_vtx_sim2reco_);
  ch_PU200_prompt->SetBranchAddress("match_vtx_reco2sim_",     				&match_vtx_reco2sim_);
  ch_PU200_prompt->SetBranchAddress("vtx_index_",      						&vtx_index_);
  ch_PU200_prompt->SetBranchAddress("recovtx_sim_",    						&recovtx_sim_);
  ch_PU200_prompt->SetBranchAddress("simvtx_reco_", 						&simvtx_reco_);
  ch_PU200_prompt->SetBranchAddress("simvtx_bx_",      						&simvtx_bx_);
  ch_PU200_prompt->SetBranchAddress("simvtx_evtId_",   						&simvtx_evtId_);
  ch_PU200_prompt->SetBranchAddress("recovtx_original_index_", 				&recovtx_original_index_);
  ch_PU200_prompt->SetBranchAddress("muon_isMuon_",    						&muon_isMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isPFMuon_",    					&muon_isPFMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isGlobalMuon_",    				&muon_isGlobalMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isTrackerMuon_",    				&muon_isTrackerMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isStandAloneMuon_",    			&muon_isStandAloneMuon_);
  ch_PU200_prompt->SetBranchAddress("muon_isCutBasedIdLoose_",    			&muon_isCutBasedIdLoose_);
  ch_PU200_prompt->SetBranchAddress("muon_isLooseMuon_",    				&muon_isLooseMuon_);
  ch_PU200_prompt->SetBranchAddress("track_genMatched_",    				&track_genMatched_);
  ch_PU200_prompt->SetBranchAddress("track_evtId_", 	   					&track_evtId_);
  ch_PU200_prompt->SetBranchAddress("track_bx_", 	       					&track_bx_);
  ch_PU200_prompt->SetBranchAddress("vtx_time_", 	       					&vtx_time_);
  ch_PU200_prompt->SetBranchAddress("vtx_time_err_", 	      				&vtx_time_err_);
  ch_PU200_prompt->SetBranchAddress("muon_time_sim_", 	       				&muon_time_sim_);
  ch_PU200_prompt->SetBranchAddress("track_time_sim_", 	       				&track_time_sim_);
  ch_PU200_prompt->SetBranchAddress("simvtx_time_", 	       				&simvtx_time_);

  ch_PU200_nonprompt->SetBranchAddress("muon_pt_",        					&muon_pt_);
  ch_PU200_nonprompt->SetBranchAddress("muon_time_",      					&muon_time_);
  ch_PU200_nonprompt->SetBranchAddress("muon_time_err_",  					&muon_time_err_);
  ch_PU200_nonprompt->SetBranchAddress("muon_time_err_i_", 					&muon_time_err_i_);
  ch_PU200_nonprompt->SetBranchAddress("muon_mva_",     					&muon_mva_);
  ch_PU200_nonprompt->SetBranchAddress("muon_prompt_",   					&muon_prompt_);
  ch_PU200_nonprompt->SetBranchAddress("muon_isBarrel_",  					&muon_isBarrel_);
  ch_PU200_nonprompt->SetBranchAddress("muon_status_",    					&muon_status_);
  ch_PU200_nonprompt->SetBranchAddress("muon_pv_dz_",     					&muon_pv_dz_);
  ch_PU200_nonprompt->SetBranchAddress("muon_pv_dxy_",    					&muon_pv_dxy_);
  ch_PU200_nonprompt->SetBranchAddress("muon_vz_",        					&muon_vz_);
  ch_PU200_nonprompt->SetBranchAddress("muon_PVweight_",  					&muon_PVweight_);
  ch_PU200_nonprompt->SetBranchAddress("track_pt_",       					&track_pt_);
  ch_PU200_nonprompt->SetBranchAddress("track_time_",     					&track_time_);
  ch_PU200_nonprompt->SetBranchAddress("track_time_err_", 					&track_time_err_);
  ch_PU200_nonprompt->SetBranchAddress("track_time_err_i_", 				&track_time_err_i_);
  ch_PU200_nonprompt->SetBranchAddress("track_mva_",     					&track_mva_);
  ch_PU200_nonprompt->SetBranchAddress("track_vz_",       					&track_vz_);
  ch_PU200_nonprompt->SetBranchAddress("track_pv_dz_",    					&track_pv_dz_);
  ch_PU200_nonprompt->SetBranchAddress("track_PVweight_", 					&track_PVweight_);
  ch_PU200_nonprompt->SetBranchAddress("track_type_", 	   					&track_type_);
  ch_PU200_nonprompt->SetBranchAddress("selectedLV_",     					&selectedLV_);
  ch_PU200_nonprompt->SetBranchAddress("match_vtx_sim2reco_",     			&match_vtx_sim2reco_);
  ch_PU200_nonprompt->SetBranchAddress("match_vtx_reco2sim_",     			&match_vtx_reco2sim_);
  ch_PU200_nonprompt->SetBranchAddress("vtx_index_",      					&vtx_index_);
  ch_PU200_nonprompt->SetBranchAddress("recovtx_sim_",    					&recovtx_sim_);
  ch_PU200_nonprompt->SetBranchAddress("simvtx_reco_", 						&simvtx_reco_);
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
  ch_PU200_nonprompt->SetBranchAddress("track_evtId_", 	   					&track_evtId_);
  ch_PU200_nonprompt->SetBranchAddress("track_bx_", 	       				&track_bx_);
  ch_PU200_nonprompt->SetBranchAddress("vtx_time_", 	       				&vtx_time_);
  ch_PU200_nonprompt->SetBranchAddress("vtx_time_err_", 	      			&vtx_time_err_);
  ch_PU200_nonprompt->SetBranchAddress("muon_time_sim_", 	       			&muon_time_sim_);
  ch_PU200_nonprompt->SetBranchAddress("track_time_sim_", 	       			&track_time_sim_);
  ch_PU200_nonprompt->SetBranchAddress("simvtx_time_", 	       				&simvtx_time_);

  // noPU
  ch_noPU_prompt->SetBranchAddress("muon_pt_",        						&muon_pt_);
  ch_noPU_prompt->SetBranchAddress("muon_time_",      						&muon_time_);
  ch_noPU_prompt->SetBranchAddress("muon_time_err_",  						&muon_time_err_);
  ch_noPU_prompt->SetBranchAddress("muon_time_err_i_", 						&muon_time_err_i_);
  ch_noPU_prompt->SetBranchAddress("muon_mva_",     						&muon_mva_);
  ch_noPU_prompt->SetBranchAddress("muon_prompt_",   						&muon_prompt_);
  ch_noPU_prompt->SetBranchAddress("muon_isBarrel_",  						&muon_isBarrel_);
  ch_noPU_prompt->SetBranchAddress("muon_status_",    						&muon_status_);
  ch_noPU_prompt->SetBranchAddress("muon_pv_dz_",     						&muon_pv_dz_);
  ch_noPU_prompt->SetBranchAddress("muon_pv_dxy_",    						&muon_pv_dxy_);
  ch_noPU_prompt->SetBranchAddress("muon_vz_",        						&muon_vz_);
  ch_noPU_prompt->SetBranchAddress("muon_PVweight_",  						&muon_PVweight_);
  ch_noPU_prompt->SetBranchAddress("track_pt_",       						&track_pt_);
  ch_noPU_prompt->SetBranchAddress("track_time_",     						&track_time_);
  ch_noPU_prompt->SetBranchAddress("track_time_err_", 						&track_time_err_);
  ch_noPU_prompt->SetBranchAddress("track_time_err_i_", 					&track_time_err_i_);
  ch_noPU_prompt->SetBranchAddress("track_mva_",     						&track_mva_);
  ch_noPU_prompt->SetBranchAddress("track_vz_",       						&track_vz_);
  ch_noPU_prompt->SetBranchAddress("track_pv_dz_",    						&track_pv_dz_);
  ch_noPU_prompt->SetBranchAddress("track_PVweight_", 						&track_PVweight_);
  ch_noPU_prompt->SetBranchAddress("track_type_", 	   						&track_type_);
  ch_noPU_prompt->SetBranchAddress("selectedLV_",     						&selectedLV_);
  ch_noPU_prompt->SetBranchAddress("match_vtx_sim2reco_",     				&match_vtx_sim2reco_);
  ch_noPU_prompt->SetBranchAddress("match_vtx_reco2sim_",     				&match_vtx_reco2sim_);
  ch_noPU_prompt->SetBranchAddress("vtx_index_",      						&vtx_index_);
  ch_noPU_prompt->SetBranchAddress("recovtx_sim_",    						&recovtx_sim_);
  ch_noPU_prompt->SetBranchAddress("simvtx_reco_", 							&simvtx_reco_);
  ch_noPU_prompt->SetBranchAddress("simvtx_bx_",      						&simvtx_bx_);
  ch_noPU_prompt->SetBranchAddress("simvtx_evtId_",   						&simvtx_evtId_);
  ch_noPU_prompt->SetBranchAddress("recovtx_original_index_", 				&recovtx_original_index_);
  ch_noPU_prompt->SetBranchAddress("muon_isMuon_",    						&muon_isMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isPFMuon_",    					&muon_isPFMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isGlobalMuon_",    				&muon_isGlobalMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isTrackerMuon_",    				&muon_isTrackerMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isStandAloneMuon_",    			&muon_isStandAloneMuon_);
  ch_noPU_prompt->SetBranchAddress("muon_isCutBasedIdLoose_",    			&muon_isCutBasedIdLoose_);
  ch_noPU_prompt->SetBranchAddress("muon_isLooseMuon_",    					&muon_isLooseMuon_);
  ch_noPU_prompt->SetBranchAddress("track_genMatched_",    					&track_genMatched_);
  ch_noPU_prompt->SetBranchAddress("track_evtId_", 	   						&track_evtId_);
  ch_noPU_prompt->SetBranchAddress("track_bx_", 	       					&track_bx_);
  ch_noPU_prompt->SetBranchAddress("vtx_time_", 	       					&vtx_time_);
  ch_noPU_prompt->SetBranchAddress("vtx_time_err_", 	      				&vtx_time_err_);
  ch_noPU_prompt->SetBranchAddress("muon_time_sim_", 	       				&muon_time_sim_);
  ch_noPU_prompt->SetBranchAddress("track_time_sim_", 	       				&track_time_sim_);
  ch_noPU_prompt->SetBranchAddress("simvtx_time_", 	       					&simvtx_time_);

  ch_noPU_nonprompt->SetBranchAddress("muon_pt_",        					&muon_pt_);
  ch_noPU_nonprompt->SetBranchAddress("muon_time_",      					&muon_time_);
  ch_noPU_nonprompt->SetBranchAddress("muon_time_err_",  					&muon_time_err_);
  ch_noPU_nonprompt->SetBranchAddress("muon_time_err_i_", 					&muon_time_err_i_);
  ch_noPU_nonprompt->SetBranchAddress("muon_mva_",     						&muon_mva_);
  ch_noPU_nonprompt->SetBranchAddress("muon_prompt_",   					&muon_prompt_);
  ch_noPU_nonprompt->SetBranchAddress("muon_isBarrel_",  					&muon_isBarrel_);
  ch_noPU_nonprompt->SetBranchAddress("muon_status_",    					&muon_status_);
  ch_noPU_nonprompt->SetBranchAddress("muon_pv_dz_",     					&muon_pv_dz_);
  ch_noPU_nonprompt->SetBranchAddress("muon_pv_dxy_",    					&muon_pv_dxy_);
  ch_noPU_nonprompt->SetBranchAddress("muon_vz_",        					&muon_vz_);
  ch_noPU_nonprompt->SetBranchAddress("muon_PVweight_",  					&muon_PVweight_);
  ch_noPU_nonprompt->SetBranchAddress("track_pt_",       					&track_pt_);
  ch_noPU_nonprompt->SetBranchAddress("track_time_",     					&track_time_);
  ch_noPU_nonprompt->SetBranchAddress("track_time_err_", 					&track_time_err_);
  ch_noPU_nonprompt->SetBranchAddress("track_time_err_i_", 					&track_time_err_i_);
  ch_noPU_nonprompt->SetBranchAddress("track_mva_",     					&track_mva_);
  ch_noPU_nonprompt->SetBranchAddress("track_vz_",       					&track_vz_);
  ch_noPU_nonprompt->SetBranchAddress("track_pv_dz_",    					&track_pv_dz_);
  ch_noPU_nonprompt->SetBranchAddress("track_PVweight_", 					&track_PVweight_);
  ch_noPU_nonprompt->SetBranchAddress("track_type_", 	   					&track_type_);
  ch_noPU_nonprompt->SetBranchAddress("selectedLV_",     					&selectedLV_);
  ch_noPU_nonprompt->SetBranchAddress("match_vtx_sim2reco_",     			&match_vtx_sim2reco_);
  ch_noPU_nonprompt->SetBranchAddress("match_vtx_reco2sim_",     			&match_vtx_reco2sim_);
  ch_noPU_nonprompt->SetBranchAddress("vtx_index_",      					&vtx_index_);
  ch_noPU_nonprompt->SetBranchAddress("recovtx_sim_",    					&recovtx_sim_);
  ch_noPU_nonprompt->SetBranchAddress("simvtx_reco_", 						&simvtx_reco_);
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
  ch_noPU_nonprompt->SetBranchAddress("track_evtId_", 	   					&track_evtId_);
  ch_noPU_nonprompt->SetBranchAddress("track_bx_", 	       					&track_bx_);
  ch_noPU_nonprompt->SetBranchAddress("vtx_time_", 	       					&vtx_time_);
  ch_noPU_nonprompt->SetBranchAddress("vtx_time_err_", 	      				&vtx_time_err_);
  ch_noPU_nonprompt->SetBranchAddress("muon_time_sim_", 	       			&muon_time_sim_);
  ch_noPU_nonprompt->SetBranchAddress("track_time_sim_", 	       			&track_time_sim_);
  ch_noPU_nonprompt->SetBranchAddress("simvtx_time_", 	       				&simvtx_time_);



  ///////////////////////
  // Define histograms //
  ///////////////////////

  int nbin=20;


      ///////////////////////
      //   (muon, track)   //
      ///////////////////////

  // PU200
    // prompt
  TH1D* h_PU200_prompt_PV_EB         = new TH1D("h_PU200_prompt_PV_EB",           "h_PU200_prompt_PV_EB",           nbin, 0, 10);
  TH1D* h_PU200_prompt_SV_EB         = new TH1D("h_PU200_prompt_SV_EB",           "h_PU200_prompt_SV_EB",           nbin, 0, 10);
  TH1D* h_PU200_prompt_PU_EB         = new TH1D("h_PU200_prompt_PU_EB",           "h_PU200_prompt_PU_EB",           nbin, 0, 10);
  TH1D* h_PU200_prompt_fake_EB       = new TH1D("h_PU200_prompt_fake_EB",         "h_PU200_prompt_fake_EB",         nbin, 0, 10);
  TH1D* h_PU200_prompt_no_tErr_EB     = new TH1D("h_PU200_prompt_no_tErr_EB",       "h_PU200_prompt_no_tErr_EB",       nbin, 0, 10);
  TH1D* h_PU200_prompt_no_tErr_PV_EB     = new TH1D("h_PU200_prompt_no_tErr_PV_EB",       "h_PU200_prompt_no_tErr_PV_EB",       nbin, 0, 10);
  TH1D* h_PU200_prompt_not_vtx_matched_EB     = new TH1D("h_PU200_prompt_not_vtx_matched_EB",       "h_PU200_prompt_not_vtx_matched_EB",       nbin, 0, 10);
  TH1D* h_PU200_prompt_bx_PV_EB         = new TH1D("h_PU200_prompt_bx_PV_EB",           "h_PU200_prompt_bx_PV_EB",           10, 0, 10);
  TH1D* h_PU200_prompt_evtId_PV_EB         = new TH1D("h_PU200_prompt_evtId_PV_EB",           "h_PU200_prompt_evtId_PV_EB",           250, 0, 250);
  TH1D* h_PU200_prompt_PVweight_PV_EB         = new TH1D("h_PU200_prompt_PVweight_PV_EB",           "h_PU200_prompt_PVweight_PV_EB",           100, 0, 1);
  TH1D* h_PU200_prompt_simvtx_bx_PV_EB         = new TH1D("h_PU200_prompt_simvtx_bx_PV_EB",           "h_PU200_prompt_simvtx_bx_PV_EB",           10, 0, 10);
  TH1D* h_PU200_prompt_simvtx_evtId_PV_EB         = new TH1D("h_PU200_prompt_simvtx_evtId_PV_EB",           "h_PU200_prompt_simvtx_evtId_PV_EB",           250, 0, 250);
  TH1D* h_PU200_prompt_PV_EE         = new TH1D("h_PU200_prompt_PV_EE",           "h_PU200_prompt_PV_EE",           nbin, 0, 10);
  TH1D* h_PU200_prompt_SV_EE         = new TH1D("h_PU200_prompt_SV_EE",           "h_PU200_prompt_SV_EE",           nbin, 0, 10);
  TH1D* h_PU200_prompt_PU_EE         = new TH1D("h_PU200_prompt_PU_EE",           "h_PU200_prompt_PU_EE",           nbin, 0, 10);
  TH1D* h_PU200_prompt_fake_EE       = new TH1D("h_PU200_prompt_fake_EE",         "h_PU200_prompt_fake_EE",         nbin, 0, 10);
  TH1D* h_PU200_prompt_no_tErr_EE     = new TH1D("h_PU200_prompt_no_tErr_EE",       "h_PU200_prompt_no_tErr_EE",       nbin, 0, 10);
  TH1D* h_PU200_prompt_not_vtx_matched_EE     = new TH1D("h_PU200_prompt_not_vtx_matched_EE",       "h_PU200_prompt_not_vtx_matched_EE",       nbin, 0, 10);
  TH1D* h_PU200_prompt_bx_PV_EE         = new TH1D("h_PU200_prompt_bx_PV_EE",           "h_PU200_prompt_bx_PV_EE",           10, 0, 10);
  TH1D* h_PU200_prompt_evtId_PV_EE         = new TH1D("h_PU200_prompt_evtId_PV_EE",           "h_PU200_prompt_evtId_PV_EE",           250, 0, 250);
  TH1D* h_PU200_prompt_PVweight_PV_EE         = new TH1D("h_PU200_prompt_PVweight_PV_EE",           "h_PU200_prompt_PVweight_PV_EE",           100, 0, 1);
  TH1D* h_PU200_prompt_simvtx_bx_PV_EE         = new TH1D("h_PU200_prompt_simvtx_bx_PV_EE",           "h_PU200_prompt_simvtx_bx_PV_EE",           10, 0, 10);
  TH1D* h_PU200_prompt_simvtx_evtId_PV_EE         = new TH1D("h_PU200_prompt_simvtx_evtId_PV_EE",           "h_PU200_prompt_simvtx_evtId_PV_EE",           250, 0, 250);
    // nonprompt
  TH1D* h_PU200_nonprompt_PV_EB         = new TH1D("h_PU200_nonprompt_PV_EB",           "h_PU200_nonprompt_PV_EB",           nbin, 0, 10);
  TH1D* h_PU200_nonprompt_SV_EB         = new TH1D("h_PU200_nonprompt_SV_EB",           "h_PU200_nonprompt_SV_EB",           nbin, 0, 10);
  TH1D* h_PU200_nonprompt_PU_EB         = new TH1D("h_PU200_nonprompt_PU_EB",           "h_PU200_nonprompt_PU_EB",           nbin, 0, 10);
  TH1D* h_PU200_nonprompt_fake_EB       = new TH1D("h_PU200_nonprompt_fake_EB",         "h_PU200_nonprompt_fake_EB",         nbin, 0, 10);
  TH1D* h_PU200_nonprompt_no_tErr_EB     = new TH1D("h_PU200_nonprompt_no_tErr_EB",       "h_PU200_nonprompt_no_tErr_EB",       nbin, 0, 10);
  TH1D* h_PU200_nonprompt_not_vtx_matched_EB     = new TH1D("h_PU200_nonprompt_not_vtx_matched_EB",       "h_PU200_nonprompt_not_vtx_matched_EB",       nbin, 0, 10);
  TH1D* h_PU200_nonprompt_bx_PV_EB         = new TH1D("h_PU200_nonprompt_bx_PV_EB",           "h_PU200_nonprompt_bx_PV_EB",           10, 0, 10);
  TH1D* h_PU200_nonprompt_evtId_PV_EB         = new TH1D("h_PU200_nonprompt_evtId_PV_EB",           "h_PU200_nonprompt_evtId_PV_EB",           250, 0, 250);
  TH1D* h_PU200_nonprompt_PVweight_PV_EB         = new TH1D("h_PU200_nonprompt_PVweight_PV_EB",           "h_PU200_nonprompt_PVweight_PV_EB",           100, 0, 1);
  TH1D* h_PU200_nonprompt_simvtx_bx_PV_EB         = new TH1D("h_PU200_nonprompt_simvtx_bx_PV_EB",           "h_PU200_nonprompt_simvtx_bx_PV_EB",           10, 0, 10);
  TH1D* h_PU200_nonprompt_simvtx_evtId_PV_EB         = new TH1D("h_PU200_nonprompt_simvtx_evtId_PV_EB",           "h_PU200_nonprompt_simvtx_evtId_PV_EB",           250, 0, 250);
  TH1D* h_PU200_nonprompt_PV_EE         = new TH1D("h_PU200_nonprompt_PV_EE",           "h_PU200_nonprompt_PV_EE",           nbin, 0, 10);
  TH1D* h_PU200_nonprompt_SV_EE         = new TH1D("h_PU200_nonprompt_SV_EE",           "h_PU200_nonprompt_SV_EE",           nbin, 0, 10);
  TH1D* h_PU200_nonprompt_PU_EE         = new TH1D("h_PU200_nonprompt_PU_EE",           "h_PU200_nonprompt_PU_EE",           nbin, 0, 10);
  TH1D* h_PU200_nonprompt_fake_EE       = new TH1D("h_PU200_nonprompt_fake_EE",         "h_PU200_nonprompt_fake_EE",         nbin, 0, 10);
  TH1D* h_PU200_nonprompt_no_tErr_EE     = new TH1D("h_PU200_nonprompt_no_tErr_EE",       "h_PU200_nonprompt_no_tErr_EE",       nbin, 0, 10);
  TH1D* h_PU200_nonprompt_not_vtx_matched_EE     = new TH1D("h_PU200_nonprompt_not_vtx_matched_EE",       "h_PU200_nonprompt_not_vtx_matched_EE",       nbin, 0, 10);
  TH1D* h_PU200_nonprompt_bx_PV_EE         = new TH1D("h_PU200_nonprompt_bx_PV_EE",           "h_PU200_nonprompt_bx_PV_EE",           10, 0, 10);
  TH1D* h_PU200_nonprompt_evtId_PV_EE         = new TH1D("h_PU200_nonprompt_evtId_PV_EE",           "h_PU200_nonprompt_evtId_PV_EE",           250, 0, 250);
  TH1D* h_PU200_nonprompt_PVweight_PV_EE         = new TH1D("h_PU200_nonprompt_PVweight_PV_EE",           "h_PU200_nonprompt_PVweight_PV_EE",           100, 0, 1);
  TH1D* h_PU200_nonprompt_simvtx_bx_PV_EE         = new TH1D("h_PU200_nonprompt_simvtx_bx_PV_EE",           "h_PU200_nonprompt_simvtx_bx_PV_EE",           10, 0, 10);
  TH1D* h_PU200_nonprompt_simvtx_evtId_PV_EE         = new TH1D("h_PU200_nonprompt_simvtx_evtId_PV_EE",           "h_PU200_nonprompt_simvtx_evtId_PV_EE",           250, 0, 250);

  // noPU
    // prompt
  TH1D* h_noPU_prompt_PV_EB         = new TH1D("h_noPU_prompt_PV_EB",           "h_noPU_prompt_PV_EB",           nbin, 0, 10);
  TH1D* h_noPU_prompt_SV_EB         = new TH1D("h_noPU_prompt_SV_EB",           "h_noPU_prompt_SV_EB",           nbin, 0, 10);
  TH1D* h_noPU_prompt_PU_EB         = new TH1D("h_noPU_prompt_PU_EB",           "h_noPU_prompt_PU_EB",           nbin, 0, 10);
  TH1D* h_noPU_prompt_fake_EB       = new TH1D("h_noPU_prompt_fake_EB",         "h_noPU_prompt_fake_EB",         nbin, 0, 10);
  TH1D* h_noPU_prompt_no_tErr_EB     = new TH1D("h_noPU_prompt_no_tErr_EB",       "h_noPU_prompt_no_tErr_EB",       nbin, 0, 10);
  TH1D* h_noPU_prompt_not_vtx_matched_EB     = new TH1D("h_noPU_prompt_not_vtx_matched_EB",       "h_noPU_prompt_not_vtx_matched_EB",       nbin, 0, 10);
  TH1D* h_noPU_prompt_bx_PV_EB         = new TH1D("h_noPU_prompt_bx_PV_EB",           "h_noPU_prompt_bx_PV_EB",           10, 0, 10);
  TH1D* h_noPU_prompt_evtId_PV_EB         = new TH1D("h_noPU_prompt_evtId_PV_EB",           "h_noPU_prompt_evtId_PV_EB",           250, 0, 250);
  TH1D* h_noPU_prompt_PVweight_PV_EB         = new TH1D("h_noPU_prompt_PVweight_PV_EB",           "h_noPU_prompt_PVweight_PV_EB",           100, 0, 1);
  TH1D* h_noPU_prompt_simvtx_bx_PV_EB         = new TH1D("h_noPU_prompt_simvtx_bx_PV_EB",           "h_noPU_prompt_simvtx_bx_PV_EB",           10, 0, 10);
  TH1D* h_noPU_prompt_simvtx_evtId_PV_EB         = new TH1D("h_noPU_prompt_simvtx_evtId_PV_EB",           "h_noPU_prompt_simvtx_evtId_PV_EB",           250, 0, 250);
  TH1D* h_noPU_prompt_PV_EE         = new TH1D("h_noPU_prompt_PV_EE",           "h_noPU_prompt_PV_EE",           nbin, 0, 10);
  TH1D* h_noPU_prompt_SV_EE         = new TH1D("h_noPU_prompt_SV_EE",           "h_noPU_prompt_SV_EE",           nbin, 0, 10);
  TH1D* h_noPU_prompt_PU_EE         = new TH1D("h_noPU_prompt_PU_EE",           "h_noPU_prompt_PU_EE",           nbin, 0, 10);
  TH1D* h_noPU_prompt_fake_EE       = new TH1D("h_noPU_prompt_fake_EE",         "h_noPU_prompt_fake_EE",         nbin, 0, 10);
  TH1D* h_noPU_prompt_no_tErr_EE     = new TH1D("h_noPU_prompt_no_tErr_EE",       "h_noPU_prompt_no_tErr_EE",       nbin, 0, 10);
  TH1D* h_noPU_prompt_not_vtx_matched_EE     = new TH1D("h_noPU_prompt_not_vtx_matched_EE",       "h_noPU_prompt_not_vtx_matched_EE",       nbin, 0, 10);
  TH1D* h_noPU_prompt_bx_PV_EE         = new TH1D("h_noPU_prompt_bx_PV_EE",           "h_noPU_prompt_bx_PV_EE",           10, 0, 10);
  TH1D* h_noPU_prompt_evtId_PV_EE         = new TH1D("h_noPU_prompt_evtId_PV_EE",           "h_noPU_prompt_evtId_PV_EE",           250, 0, 250);
  TH1D* h_noPU_prompt_PVweight_PV_EE         = new TH1D("h_noPU_prompt_PVweight_PV_EE",           "h_noPU_prompt_PVweight_PV_EE",           100, 0, 1);
  TH1D* h_noPU_prompt_simvtx_bx_PV_EE         = new TH1D("h_noPU_prompt_simvtx_bx_PV_EE",           "h_noPU_prompt_simvtx_bx_PV_EE",           10, 0, 10);
  TH1D* h_noPU_prompt_simvtx_evtId_PV_EE         = new TH1D("h_noPU_prompt_simvtx_evtId_PV_EE",           "h_noPU_prompt_simvtx_evtId_PV_EE",           250, 0, 250);
    // nonprompt
  TH1D* h_noPU_nonprompt_PV_EB         = new TH1D("h_noPU_nonprompt_PV_EB",           "h_noPU_nonprompt_PV_EB",           nbin, 0, 10);
  TH1D* h_noPU_nonprompt_SV_EB         = new TH1D("h_noPU_nonprompt_SV_EB",           "h_noPU_nonprompt_SV_EB",           nbin, 0, 10);
  TH1D* h_noPU_nonprompt_PU_EB         = new TH1D("h_noPU_nonprompt_PU_EB",           "h_noPU_nonprompt_PU_EB",           nbin, 0, 10);
  TH1D* h_noPU_nonprompt_fake_EB       = new TH1D("h_noPU_nonprompt_fake_EB",         "h_noPU_nonprompt_fake_EB",         nbin, 0, 10);
  TH1D* h_noPU_nonprompt_no_tErr_EB     = new TH1D("h_noPU_nonprompt_no_tErr_EB",       "h_noPU_nonprompt_no_tErr_EB",       nbin, 0, 10);
  TH1D* h_noPU_nonprompt_not_vtx_matched_EB     = new TH1D("h_noPU_nonprompt_not_vtx_matched_EB",       "h_noPU_nonprompt_not_vtx_matched_EB",       nbin, 0, 10);
  TH1D* h_noPU_nonprompt_bx_PV_EB         = new TH1D("h_noPU_nonprompt_bx_PV_EB",           "h_noPU_nonprompt_bx_PV_EB",           10, 0, 10);
  TH1D* h_noPU_nonprompt_evtId_PV_EB         = new TH1D("h_noPU_nonprompt_evtId_PV_EB",           "h_noPU_nonprompt_evtId_PV_EB",           250, 0, 250);
  TH1D* h_noPU_nonprompt_PVweight_PV_EB         = new TH1D("h_noPU_nonprompt_PVweight_PV_EB",           "h_noPU_nonprompt_PVweight_PV_EB",           100, 0, 1);
  TH1D* h_noPU_nonprompt_simvtx_bx_PV_EB         = new TH1D("h_noPU_nonprompt_simvtx_bx_PV_EB",           "h_noPU_nonprompt_simvtx_bx_PV_EB",           10, 0, 10);
  TH1D* h_noPU_nonprompt_simvtx_evtId_PV_EB         = new TH1D("h_noPU_nonprompt_simvtx_evtId_PV_EB",           "h_noPU_nonprompt_simvtx_evtId_PV_EB",           250, 0, 250);
  TH1D* h_noPU_nonprompt_PV_EE         = new TH1D("h_noPU_nonprompt_PV_EE",           "h_noPU_nonprompt_PV_EE",           nbin, 0, 10);
  TH1D* h_noPU_nonprompt_SV_EE         = new TH1D("h_noPU_nonprompt_SV_EE",           "h_noPU_nonprompt_SV_EE",           nbin, 0, 10);
  TH1D* h_noPU_nonprompt_PU_EE         = new TH1D("h_noPU_nonprompt_PU_EE",           "h_noPU_nonprompt_PU_EE",           nbin, 0, 10);
  TH1D* h_noPU_nonprompt_fake_EE       = new TH1D("h_noPU_nonprompt_fake_EE",         "h_noPU_nonprompt_fake_EE",         nbin, 0, 10);
  TH1D* h_noPU_nonprompt_no_tErr_EE     = new TH1D("h_noPU_nonprompt_no_tErr_EE",       "h_noPU_nonprompt_no_tErr_EE",       nbin, 0, 10);
  TH1D* h_noPU_nonprompt_not_vtx_matched_EE     = new TH1D("h_noPU_nonprompt_not_vtx_matched_EE",       "h_noPU_nonprompt_not_vtx_matched_EE",       nbin, 0, 10);
  TH1D* h_noPU_nonprompt_bx_PV_EE         = new TH1D("h_noPU_nonprompt_bx_PV_EE",           "h_noPU_nonprompt_bx_PV_EE",           10, 0, 10);
  TH1D* h_noPU_nonprompt_evtId_PV_EE         = new TH1D("h_noPU_nonprompt_evtId_PV_EE",           "h_noPU_nonprompt_evtId_PV_EE",           250, 0, 250);
  TH1D* h_noPU_nonprompt_PVweight_PV_EE         = new TH1D("h_noPU_nonprompt_PVweight_PV_EE",           "h_noPU_nonprompt_PVweight_PV_EE",           100, 0, 1);
  TH1D* h_noPU_nonprompt_simvtx_bx_PV_EE         = new TH1D("h_noPU_nonprompt_simvtx_bx_PV_EE",           "h_noPU_nonprompt_simvtx_bx_PV_EE",           10, 0, 10);
  TH1D* h_noPU_nonprompt_simvtx_evtId_PV_EE         = new TH1D("h_noPU_nonprompt_simvtx_evtId_PV_EE",           "h_noPU_nonprompt_simvtx_evtId_PV_EE",           250, 0, 250);


      ///////////////////////
      //    (PV, track)    //
      ///////////////////////

  // PU200
    // prompt
  TH1D* h_PU200_prompt_PV_EB_vtx         = new TH1D("h_PU200_prompt_PV_EB_vtx",           "h_PU200_prompt_PV_EB_vtx",           nbin, 0, 10);
  TH1D* h_PU200_prompt_SV_EB_vtx         = new TH1D("h_PU200_prompt_SV_EB_vtx",           "h_PU200_prompt_SV_EB_vtx",           nbin, 0, 10);
  TH1D* h_PU200_prompt_PU_EB_vtx         = new TH1D("h_PU200_prompt_PU_EB_vtx",           "h_PU200_prompt_PU_EB_vtx",           nbin, 0, 10);
  TH1D* h_PU200_prompt_fake_EB_vtx       = new TH1D("h_PU200_prompt_fake_EB_vtx",         "h_PU200_prompt_fake_EB_vtx",         nbin, 0, 10);
  TH1D* h_PU200_prompt_no_tErr_EB_vtx     = new TH1D("h_PU200_prompt_no_tErr_EB_vtx",       "h_PU200_prompt_no_tErr_EB_vtx",       nbin, 0, 10);
  TH1D* h_PU200_prompt_no_tErr_PV_EB_vtx     = new TH1D("h_PU200_prompt_no_tErr_PV_EB_vtx",       "h_PU200_prompt_no_tErr_PV_EB_vtx",       nbin, 0, 10);
  TH1D* h_PU200_prompt_not_vtx_matched_EB_vtx     = new TH1D("h_PU200_prompt_not_vtx_matched_EB_vtx",       "h_PU200_prompt_not_vtx_matched_EB_vtx",       nbin, 0, 10);
  TH1D* h_PU200_prompt_bx_PV_EB_vtx         = new TH1D("h_PU200_prompt_bx_PV_EB_vtx",           "h_PU200_prompt_bx_PV_EB_vtx",           10, 0, 10);
  TH1D* h_PU200_prompt_evtId_PV_EB_vtx         = new TH1D("h_PU200_prompt_evtId_PV_EB_vtx",           "h_PU200_prompt_evtId_PV_EB_vtx",           250, 0, 250);
  TH1D* h_PU200_prompt_PVweight_PV_EB_vtx         = new TH1D("h_PU200_prompt_PVweight_PV_EB_vtx",           "h_PU200_prompt_PVweight_PV_EB_vtx",           100, 0, 1);
  TH1D* h_PU200_prompt_simvtx_bx_PV_EB_vtx         = new TH1D("h_PU200_prompt_simvtx_bx_PV_EB_vtx",           "h_PU200_prompt_simvtx_bx_PV_EB_vtx",           10, 0, 10);
  TH1D* h_PU200_prompt_simvtx_evtId_PV_EB_vtx         = new TH1D("h_PU200_prompt_simvtx_evtId_PV_EB_vtx",           "h_PU200_prompt_simvtx_evtId_PV_EB_vtx",           250, 0, 250);
  TH1D* h_PU200_prompt_PV_EE_vtx         = new TH1D("h_PU200_prompt_PV_EE_vtx",           "h_PU200_prompt_PV_EE_vtx",           nbin, 0, 10);
  TH1D* h_PU200_prompt_SV_EE_vtx         = new TH1D("h_PU200_prompt_SV_EE_vtx",           "h_PU200_prompt_SV_EE_vtx",           nbin, 0, 10);
  TH1D* h_PU200_prompt_PU_EE_vtx         = new TH1D("h_PU200_prompt_PU_EE_vtx",           "h_PU200_prompt_PU_EE_vtx",           nbin, 0, 10);
  TH1D* h_PU200_prompt_fake_EE_vtx       = new TH1D("h_PU200_prompt_fake_EE_vtx",         "h_PU200_prompt_fake_EE_vtx",         nbin, 0, 10);
  TH1D* h_PU200_prompt_no_tErr_EE_vtx     = new TH1D("h_PU200_prompt_no_tErr_EE_vtx",       "h_PU200_prompt_no_tErr_EE_vtx",       nbin, 0, 10);
  TH1D* h_PU200_prompt_not_vtx_matched_EE_vtx     = new TH1D("h_PU200_prompt_not_vtx_matched_EE_vtx",       "h_PU200_prompt_not_vtx_matched_EE_vtx",       nbin, 0, 10);
  TH1D* h_PU200_prompt_bx_PV_EE_vtx         = new TH1D("h_PU200_prompt_bx_PV_EE_vtx",           "h_PU200_prompt_bx_PV_EE_vtx",           10, 0, 10);
  TH1D* h_PU200_prompt_evtId_PV_EE_vtx         = new TH1D("h_PU200_prompt_evtId_PV_EE_vtx",           "h_PU200_prompt_evtId_PV_EE_vtx",           250, 0, 250);
  TH1D* h_PU200_prompt_PVweight_PV_EE_vtx         = new TH1D("h_PU200_prompt_PVweight_PV_EE_vtx",           "h_PU200_prompt_PVweight_PV_EE_vtx",           100, 0, 1);
  TH1D* h_PU200_prompt_simvtx_bx_PV_EE_vtx         = new TH1D("h_PU200_prompt_simvtx_bx_PV_EE_vtx",           "h_PU200_prompt_simvtx_bx_PV_EE_vtx",           10, 0, 10);
  TH1D* h_PU200_prompt_simvtx_evtId_PV_EE_vtx         = new TH1D("h_PU200_prompt_simvtx_evtId_PV_EE_vtx",           "h_PU200_prompt_simvtx_evtId_PV_EE_vtx",           250, 0, 250);
    // nonprompt
  TH1D* h_PU200_nonprompt_PV_EB_vtx         = new TH1D("h_PU200_nonprompt_PV_EB_vtx",           "h_PU200_nonprompt_PV_EB_vtx",           nbin, 0, 10);
  TH1D* h_PU200_nonprompt_SV_EB_vtx         = new TH1D("h_PU200_nonprompt_SV_EB_vtx",           "h_PU200_nonprompt_SV_EB_vtx",           nbin, 0, 10);
  TH1D* h_PU200_nonprompt_PU_EB_vtx         = new TH1D("h_PU200_nonprompt_PU_EB_vtx",           "h_PU200_nonprompt_PU_EB_vtx",           nbin, 0, 10);
  TH1D* h_PU200_nonprompt_fake_EB_vtx       = new TH1D("h_PU200_nonprompt_fake_EB_vtx",         "h_PU200_nonprompt_fake_EB_vtx",         nbin, 0, 10);
  TH1D* h_PU200_nonprompt_no_tErr_EB_vtx     = new TH1D("h_PU200_nonprompt_no_tErr_EB_vtx",       "h_PU200_nonprompt_no_tErr_EB_vtx",       nbin, 0, 10);
  TH1D* h_PU200_nonprompt_not_vtx_matched_EB_vtx     = new TH1D("h_PU200_nonprompt_not_vtx_matched_EB_vtx",       "h_PU200_nonprompt_not_vtx_matched_EB_vtx",       nbin, 0, 10);
  TH1D* h_PU200_nonprompt_bx_PV_EB_vtx         = new TH1D("h_PU200_nonprompt_bx_PV_EB_vtx",           "h_PU200_nonprompt_bx_PV_EB_vtx",           10, 0, 10);
  TH1D* h_PU200_nonprompt_evtId_PV_EB_vtx         = new TH1D("h_PU200_nonprompt_evtId_PV_EB_vtx",           "h_PU200_nonprompt_evtId_PV_EB_vtx",           250, 0, 250);
  TH1D* h_PU200_nonprompt_PVweight_PV_EB_vtx         = new TH1D("h_PU200_nonprompt_PVweight_PV_EB_vtx",           "h_PU200_nonprompt_PVweight_PV_EB_vtx",           100, 0, 1);
  TH1D* h_PU200_nonprompt_simvtx_bx_PV_EB_vtx         = new TH1D("h_PU200_nonprompt_simvtx_bx_PV_EB_vtx",           "h_PU200_nonprompt_simvtx_bx_PV_EB_vtx",           10, 0, 10);
  TH1D* h_PU200_nonprompt_simvtx_evtId_PV_EB_vtx         = new TH1D("h_PU200_nonprompt_simvtx_evtId_PV_EB_vtx",           "h_PU200_nonprompt_simvtx_evtId_PV_EB_vtx",           250, 0, 250);
  TH1D* h_PU200_nonprompt_PV_EE_vtx         = new TH1D("h_PU200_nonprompt_PV_EE_vtx",           "h_PU200_nonprompt_PV_EE_vtx",           nbin, 0, 10);
  TH1D* h_PU200_nonprompt_SV_EE_vtx         = new TH1D("h_PU200_nonprompt_SV_EE_vtx",           "h_PU200_nonprompt_SV_EE_vtx",           nbin, 0, 10);
  TH1D* h_PU200_nonprompt_PU_EE_vtx         = new TH1D("h_PU200_nonprompt_PU_EE_vtx",           "h_PU200_nonprompt_PU_EE_vtx",           nbin, 0, 10);
  TH1D* h_PU200_nonprompt_fake_EE_vtx       = new TH1D("h_PU200_nonprompt_fake_EE_vtx",         "h_PU200_nonprompt_fake_EE_vtx",         nbin, 0, 10);
  TH1D* h_PU200_nonprompt_no_tErr_EE_vtx     = new TH1D("h_PU200_nonprompt_no_tErr_EE_vtx",       "h_PU200_nonprompt_no_tErr_EE_vtx",       nbin, 0, 10);
  TH1D* h_PU200_nonprompt_not_vtx_matched_EE_vtx     = new TH1D("h_PU200_nonprompt_not_vtx_matched_EE_vtx",       "h_PU200_nonprompt_not_vtx_matched_EE_vtx",       nbin, 0, 10);
  TH1D* h_PU200_nonprompt_bx_PV_EE_vtx         = new TH1D("h_PU200_nonprompt_bx_PV_EE_vtx",           "h_PU200_nonprompt_bx_PV_EE_vtx",           10, 0, 10);
  TH1D* h_PU200_nonprompt_evtId_PV_EE_vtx         = new TH1D("h_PU200_nonprompt_evtId_PV_EE_vtx",           "h_PU200_nonprompt_evtId_PV_EE_vtx",           250, 0, 250);
  TH1D* h_PU200_nonprompt_PVweight_PV_EE_vtx         = new TH1D("h_PU200_nonprompt_PVweight_PV_EE_vtx",           "h_PU200_nonprompt_PVweight_PV_EE_vtx",           100, 0, 1);
  TH1D* h_PU200_nonprompt_simvtx_bx_PV_EE_vtx         = new TH1D("h_PU200_nonprompt_simvtx_bx_PV_EE_vtx",           "h_PU200_nonprompt_simvtx_bx_PV_EE_vtx",           10, 0, 10);
  TH1D* h_PU200_nonprompt_simvtx_evtId_PV_EE_vtx         = new TH1D("h_PU200_nonprompt_simvtx_evtId_PV_EE_vtx",           "h_PU200_nonprompt_simvtx_evtId_PV_EE_vtx",           250, 0, 250);

  // noPU
    // prompt
  TH1D* h_noPU_prompt_PV_EB_vtx         = new TH1D("h_noPU_prompt_PV_EB_vtx",           "h_noPU_prompt_PV_EB_vtx",           nbin, 0, 10);
  TH1D* h_noPU_prompt_SV_EB_vtx         = new TH1D("h_noPU_prompt_SV_EB_vtx",           "h_noPU_prompt_SV_EB_vtx",           nbin, 0, 10);
  TH1D* h_noPU_prompt_PU_EB_vtx         = new TH1D("h_noPU_prompt_PU_EB_vtx",           "h_noPU_prompt_PU_EB_vtx",           nbin, 0, 10);
  TH1D* h_noPU_prompt_fake_EB_vtx       = new TH1D("h_noPU_prompt_fake_EB_vtx",         "h_noPU_prompt_fake_EB_vtx",         nbin, 0, 10);
  TH1D* h_noPU_prompt_no_tErr_EB_vtx     = new TH1D("h_noPU_prompt_no_tErr_EB_vtx",       "h_noPU_prompt_no_tErr_EB_vtx",       nbin, 0, 10);
  TH1D* h_noPU_prompt_not_vtx_matched_EB_vtx     = new TH1D("h_noPU_prompt_not_vtx_matched_EB_vtx",       "h_noPU_prompt_not_vtx_matched_EB_vtx",       nbin, 0, 10);
  TH1D* h_noPU_prompt_bx_PV_EB_vtx         = new TH1D("h_noPU_prompt_bx_PV_EB_vtx",           "h_noPU_prompt_bx_PV_EB_vtx",           10, 0, 10);
  TH1D* h_noPU_prompt_evtId_PV_EB_vtx         = new TH1D("h_noPU_prompt_evtId_PV_EB_vtx",           "h_noPU_prompt_evtId_PV_EB_vtx",           250, 0, 250);
  TH1D* h_noPU_prompt_PVweight_PV_EB_vtx         = new TH1D("h_noPU_prompt_PVweight_PV_EB_vtx",           "h_noPU_prompt_PVweight_PV_EB_vtx",           100, 0, 1);
  TH1D* h_noPU_prompt_simvtx_bx_PV_EB_vtx         = new TH1D("h_noPU_prompt_simvtx_bx_PV_EB_vtx",           "h_noPU_prompt_simvtx_bx_PV_EB_vtx",           10, 0, 10);
  TH1D* h_noPU_prompt_simvtx_evtId_PV_EB_vtx         = new TH1D("h_noPU_prompt_simvtx_evtId_PV_EB_vtx",           "h_noPU_prompt_simvtx_evtId_PV_EB_vtx",           250, 0, 250);
  TH1D* h_noPU_prompt_PV_EE_vtx         = new TH1D("h_noPU_prompt_PV_EE_vtx",           "h_noPU_prompt_PV_EE_vtx",           nbin, 0, 10);
  TH1D* h_noPU_prompt_SV_EE_vtx         = new TH1D("h_noPU_prompt_SV_EE_vtx",           "h_noPU_prompt_SV_EE_vtx",           nbin, 0, 10);
  TH1D* h_noPU_prompt_PU_EE_vtx         = new TH1D("h_noPU_prompt_PU_EE_vtx",           "h_noPU_prompt_PU_EE_vtx",           nbin, 0, 10);
  TH1D* h_noPU_prompt_fake_EE_vtx       = new TH1D("h_noPU_prompt_fake_EE_vtx",         "h_noPU_prompt_fake_EE_vtx",         nbin, 0, 10);
  TH1D* h_noPU_prompt_no_tErr_EE_vtx     = new TH1D("h_noPU_prompt_no_tErr_EE_vtx",       "h_noPU_prompt_no_tErr_EE_vtx",       nbin, 0, 10);
  TH1D* h_noPU_prompt_not_vtx_matched_EE_vtx     = new TH1D("h_noPU_prompt_not_vtx_matched_EE_vtx",       "h_noPU_prompt_not_vtx_matched_EE_vtx",       nbin, 0, 10);
  TH1D* h_noPU_prompt_bx_PV_EE_vtx         = new TH1D("h_noPU_prompt_bx_PV_EE_vtx",           "h_noPU_prompt_bx_PV_EE_vtx",           10, 0, 10);
  TH1D* h_noPU_prompt_evtId_PV_EE_vtx         = new TH1D("h_noPU_prompt_evtId_PV_EE_vtx",           "h_noPU_prompt_evtId_PV_EE_vtx",           250, 0, 250);
  TH1D* h_noPU_prompt_PVweight_PV_EE_vtx         = new TH1D("h_noPU_prompt_PVweight_PV_EE_vtx",           "h_noPU_prompt_PVweight_PV_EE_vtx",           100, 0, 1);
  TH1D* h_noPU_prompt_simvtx_bx_PV_EE_vtx         = new TH1D("h_noPU_prompt_simvtx_bx_PV_EE_vtx",           "h_noPU_prompt_simvtx_bx_PV_EE_vtx",           10, 0, 10);
  TH1D* h_noPU_prompt_simvtx_evtId_PV_EE_vtx         = new TH1D("h_noPU_prompt_simvtx_evtId_PV_EE_vtx",           "h_noPU_prompt_simvtx_evtId_PV_EE_vtx",           250, 0, 250);
    // nonprompt
  TH1D* h_noPU_nonprompt_PV_EB_vtx         = new TH1D("h_noPU_nonprompt_PV_EB_vtx",           "h_noPU_nonprompt_PV_EB_vtx",           nbin, 0, 10);
  TH1D* h_noPU_nonprompt_SV_EB_vtx         = new TH1D("h_noPU_nonprompt_SV_EB_vtx",           "h_noPU_nonprompt_SV_EB_vtx",           nbin, 0, 10);
  TH1D* h_noPU_nonprompt_PU_EB_vtx         = new TH1D("h_noPU_nonprompt_PU_EB_vtx",           "h_noPU_nonprompt_PU_EB_vtx",           nbin, 0, 10);
  TH1D* h_noPU_nonprompt_fake_EB_vtx       = new TH1D("h_noPU_nonprompt_fake_EB_vtx",         "h_noPU_nonprompt_fake_EB_vtx",         nbin, 0, 10);
  TH1D* h_noPU_nonprompt_no_tErr_EB_vtx     = new TH1D("h_noPU_nonprompt_no_tErr_EB_vtx",       "h_noPU_nonprompt_no_tErr_EB_vtx",       nbin, 0, 10);
  TH1D* h_noPU_nonprompt_not_vtx_matched_EB_vtx     = new TH1D("h_noPU_nonprompt_not_vtx_matched_EB_vtx",       "h_noPU_nonprompt_not_vtx_matched_EB_vtx",       nbin, 0, 10);
  TH1D* h_noPU_nonprompt_bx_PV_EB_vtx         = new TH1D("h_noPU_nonprompt_bx_PV_EB_vtx",           "h_noPU_nonprompt_bx_PV_EB_vtx",           10, 0, 10);
  TH1D* h_noPU_nonprompt_evtId_PV_EB_vtx         = new TH1D("h_noPU_nonprompt_evtId_PV_EB_vtx",           "h_noPU_nonprompt_evtId_PV_EB_vtx",           250, 0, 250);
  TH1D* h_noPU_nonprompt_PVweight_PV_EB_vtx         = new TH1D("h_noPU_nonprompt_PVweight_PV_EB_vtx",           "h_noPU_nonprompt_PVweight_PV_EB_vtx",           100, 0, 1);
  TH1D* h_noPU_nonprompt_simvtx_bx_PV_EB_vtx         = new TH1D("h_noPU_nonprompt_simvtx_bx_PV_EB_vtx",           "h_noPU_nonprompt_simvtx_bx_PV_EB_vtx",           10, 0, 10);
  TH1D* h_noPU_nonprompt_simvtx_evtId_PV_EB_vtx         = new TH1D("h_noPU_nonprompt_simvtx_evtId_PV_EB_vtx",           "h_noPU_nonprompt_simvtx_evtId_PV_EB_vtx",           250, 0, 250);
  TH1D* h_noPU_nonprompt_PV_EE_vtx         = new TH1D("h_noPU_nonprompt_PV_EE_vtx",           "h_noPU_nonprompt_PV_EE_vtx",           nbin, 0, 10);
  TH1D* h_noPU_nonprompt_SV_EE_vtx         = new TH1D("h_noPU_nonprompt_SV_EE_vtx",           "h_noPU_nonprompt_SV_EE_vtx",           nbin, 0, 10);
  TH1D* h_noPU_nonprompt_PU_EE_vtx         = new TH1D("h_noPU_nonprompt_PU_EE_vtx",           "h_noPU_nonprompt_PU_EE_vtx",           nbin, 0, 10);
  TH1D* h_noPU_nonprompt_fake_EE_vtx       = new TH1D("h_noPU_nonprompt_fake_EE_vtx",         "h_noPU_nonprompt_fake_EE_vtx",         nbin, 0, 10);
  TH1D* h_noPU_nonprompt_no_tErr_EE_vtx     = new TH1D("h_noPU_nonprompt_no_tErr_EE_vtx",       "h_noPU_nonprompt_no_tErr_EE_vtx",       nbin, 0, 10);
  TH1D* h_noPU_nonprompt_not_vtx_matched_EE_vtx     = new TH1D("h_noPU_nonprompt_not_vtx_matched_EE_vtx",       "h_noPU_nonprompt_not_vtx_matched_EE_vtx",       nbin, 0, 10);
  TH1D* h_noPU_nonprompt_bx_PV_EE_vtx         = new TH1D("h_noPU_nonprompt_bx_PV_EE_vtx",           "h_noPU_nonprompt_bx_PV_EE_vtx",           10, 0, 10);
  TH1D* h_noPU_nonprompt_evtId_PV_EE_vtx         = new TH1D("h_noPU_nonprompt_evtId_PV_EE_vtx",           "h_noPU_nonprompt_evtId_PV_EE_vtx",           250, 0, 250);
  TH1D* h_noPU_nonprompt_PVweight_PV_EE_vtx         = new TH1D("h_noPU_nonprompt_PVweight_PV_EE_vtx",           "h_noPU_nonprompt_PVweight_PV_EE_vtx",           100, 0, 1);
  TH1D* h_noPU_nonprompt_simvtx_bx_PV_EE_vtx         = new TH1D("h_noPU_nonprompt_simvtx_bx_PV_EE_vtx",           "h_noPU_nonprompt_simvtx_bx_PV_EE_vtx",           10, 0, 10);
  TH1D* h_noPU_nonprompt_simvtx_evtId_PV_EE_vtx         = new TH1D("h_noPU_nonprompt_simvtx_evtId_PV_EE_vtx",           "h_noPU_nonprompt_simvtx_evtId_PV_EE_vtx",           250, 0, 250);



      ////////////////////////
      //         2D         //
      ////////////////////////
  TH2D* h_PU200_prompt_PV_EB_2D = new TH2D("h_PU200_prompt_PV_EB_2D", "h_PU200_prompt_PV_EB_2D", nbin, 0, nbin, nbin, 0, nbin);
  TH2D* h_PU200_prompt_PV_EE_2D = new TH2D("h_PU200_prompt_PV_EE_2D", "h_PU200_prompt_PV_EE_2D", nbin, 0, nbin, nbin, 0, nbin);
  TH2D* h_PU200_nonprompt_PV_EB_2D = new TH2D("h_PU200_nonprompt_PV_EB_2D", "h_PU200_nonprompt_PV_EB_2D", nbin, 0, nbin, nbin, 0, nbin);
  TH2D* h_PU200_nonprompt_PV_EE_2D = new TH2D("h_PU200_nonprompt_PV_EE_2D", "h_PU200_nonprompt_PV_EE_2D", nbin, 0, nbin, nbin, 0, nbin);
  TH2D* h_noPU_prompt_PV_EB_2D = new TH2D("h_noPU_prompt_PV_EB_2D", "h_noPU_prompt_PV_EB_2D", nbin, 0, nbin, nbin, 0, nbin);
  TH2D* h_noPU_prompt_PV_EE_2D = new TH2D("h_noPU_prompt_PV_EE_2D", "h_noPU_prompt_PV_EE_2D", nbin, 0, nbin, nbin, 0, nbin);
  TH2D* h_noPU_nonprompt_PV_EB_2D = new TH2D("h_noPU_nonprompt_PV_EB_2D", "h_noPU_nonprompt_PV_EB_2D", nbin, 0, nbin, nbin, 0, nbin);
  TH2D* h_noPU_nonprompt_PV_EE_2D = new TH2D("h_noPU_nonprompt_PV_EE_2D", "h_noPU_nonprompt_PV_EE_2D", nbin, 0, nbin, nbin, 0, nbin);

      ////////////////////////
      //        tsim        //
      ////////////////////////
  TH1D* h_PU200_prompt_PV_vtx_time_EB = new TH1D("h_PU200_prompt_PV_vtx_time_EB", "h_PU200_prompt_PV_vtx_time_EB", 200,-1,1);
  TH1D* h_PU200_prompt_PV_muon_time_EB = new TH1D("h_PU200_prompt_PV_muon_time_EB", "h_PU200_prompt_PV_muon_time_EB", 200,-1,1);
  TH1D* h_PU200_prompt_PV_track_time_EB = new TH1D("h_PU200_prompt_PV_track_time_EB", "h_PU200_prompt_PV_track_time_EB", 200,-1,1);
  TH1D* h_PU200_prompt_PU_track_time_EB = new TH1D("h_PU200_prompt_PU_track_time_EB", "h_PU200_prompt_PU_track_time_EB", 200,-1,1);
  TH1D* h_PU200_prompt_SV_track_time_EB = new TH1D("h_PU200_prompt_SV_track_time_EB", "h_PU200_prompt_SV_track_time_EB", 200,-1,1);
  TH1D* h_PU200_nonprompt_PV_vtx_time_EB = new TH1D("h_PU200_nonprompt_PV_vtx_time_EB", "h_PU200_nonprompt_PV_vtx_time_EB", 200,-1,1);
  TH1D* h_PU200_nonprompt_PV_muon_time_EB = new TH1D("h_PU200_nonprompt_PV_muon_time_EB", "h_PU200_nonprompt_PV_muon_time_EB", 200,-1,1);
  TH1D* h_PU200_nonprompt_PV_track_time_EB = new TH1D("h_PU200_nonprompt_PV_track_time_EB", "h_PU200_nonprompt_PV_track_time_EB", 200,-1,1);
  TH1D* h_PU200_nonprompt_PU_track_time_EB = new TH1D("h_PU200_nonprompt_PU_track_time_EB", "h_PU200_nonprompt_PU_track_time_EB", 200,-1,1);
  TH1D* h_PU200_nonprompt_SV_track_time_EB = new TH1D("h_PU200_nonprompt_SV_track_time_EB", "h_PU200_nonprompt_SV_track_time_EB", 200,-1,1);
  TH1D* h_noPU_prompt_PV_vtx_time_EB = new TH1D("h_noPU_prompt_PV_vtx_time_EB", "h_noPU_prompt_PV_vtx_time_EB", 200,-1,1);
  TH1D* h_noPU_prompt_PV_muon_time_EB = new TH1D("h_noPU_prompt_PV_muon_time_EB", "h_noPU_prompt_PV_muon_time_EB", 200,-1,1);
  TH1D* h_noPU_prompt_PV_track_time_EB = new TH1D("h_noPU_prompt_PV_track_time_EB", "h_noPU_prompt_PV_track_time_EB", 200,-1,1);
  TH1D* h_noPU_prompt_PU_track_time_EB = new TH1D("h_noPU_prompt_PU_track_time_EB", "h_noPU_prompt_PU_track_time_EB", 200,-1,1);
  TH1D* h_noPU_prompt_SV_track_time_EB = new TH1D("h_noPU_prompt_SV_track_time_EB", "h_noPU_prompt_SV_track_time_EB", 200,-1,1);
  TH1D* h_noPU_nonprompt_PV_vtx_time_EB = new TH1D("h_noPU_nonprompt_PV_vtx_time_EB", "h_noPU_nonprompt_PV_vtx_time_EB", 200,-1,1);
  TH1D* h_noPU_nonprompt_PV_muon_time_EB = new TH1D("h_noPU_nonprompt_PV_muon_time_EB", "h_noPU_nonprompt_PV_muon_time_EB", 200,-1,1);
  TH1D* h_noPU_nonprompt_PV_track_time_EB = new TH1D("h_noPU_nonprompt_PV_track_time_EB", "h_noPU_nonprompt_PV_track_time_EB", 200,-1,1);
  TH1D* h_noPU_nonprompt_PU_track_time_EB = new TH1D("h_noPU_nonprompt_PU_track_time_EB", "h_noPU_nonprompt_PU_track_time_EB", 200,-1,1);
  TH1D* h_noPU_nonprompt_SV_track_time_EB = new TH1D("h_noPU_nonprompt_SV_track_time_EB", "h_noPU_nonprompt_SV_track_time_EB", 200,-1,1);
  
  TH1D* h_PU200_prompt_PV_vtx_time_EE = new TH1D("h_PU200_prompt_PV_vtx_time_EE", "h_PU200_prompt_PV_vtx_time_EE", 200,-1,1);
  TH1D* h_PU200_prompt_PV_muon_time_EE = new TH1D("h_PU200_prompt_PV_muon_time_EE", "h_PU200_prompt_PV_muon_time_EE", 200,-1,1);
  TH1D* h_PU200_prompt_PV_track_time_EE = new TH1D("h_PU200_prompt_PV_track_time_EE", "h_PU200_prompt_PV_track_time_EE", 200,-1,1);
  TH1D* h_PU200_prompt_PU_track_time_EE = new TH1D("h_PU200_prompt_PU_track_time_EE", "h_PU200_prompt_PU_track_time_EE", 200,-1,1);
  TH1D* h_PU200_prompt_SV_track_time_EE = new TH1D("h_PU200_prompt_SV_track_time_EE", "h_PU200_prompt_SV_track_time_EE", 200,-1,1);
  TH1D* h_PU200_nonprompt_PV_vtx_time_EE = new TH1D("h_PU200_nonprompt_PV_vtx_time_EE", "h_PU200_nonprompt_PV_vtx_time_EE", 200,-1,1);
  TH1D* h_PU200_nonprompt_PV_muon_time_EE = new TH1D("h_PU200_nonprompt_PV_muon_time_EE", "h_PU200_nonprompt_PV_muon_time_EE", 200,-1,1);
  TH1D* h_PU200_nonprompt_PV_track_time_EE = new TH1D("h_PU200_nonprompt_PV_track_time_EE", "h_PU200_nonprompt_PV_track_time_EE", 200,-1,1);
  TH1D* h_PU200_nonprompt_PU_track_time_EE = new TH1D("h_PU200_nonprompt_PU_track_time_EE", "h_PU200_nonprompt_PU_track_time_EE", 200,-1,1);
  TH1D* h_PU200_nonprompt_SV_track_time_EE = new TH1D("h_PU200_nonprompt_SV_track_time_EE", "h_PU200_nonprompt_SV_track_time_EE", 200,-1,1);
  TH1D* h_noPU_prompt_PV_vtx_time_EE = new TH1D("h_noPU_prompt_PV_vtx_time_EE", "h_noPU_prompt_PV_vtx_time_EE", 200,-1,1);
  TH1D* h_noPU_prompt_PV_muon_time_EE = new TH1D("h_noPU_prompt_PV_muon_time_EE", "h_noPU_prompt_PV_muon_time_EE", 200,-1,1);
  TH1D* h_noPU_prompt_PV_track_time_EE = new TH1D("h_noPU_prompt_PV_track_time_EE", "h_noPU_prompt_PV_track_time_EE", 200,-1,1);
  TH1D* h_noPU_prompt_PU_track_time_EE = new TH1D("h_noPU_prompt_PU_track_time_EE", "h_noPU_prompt_PU_track_time_EE", 200,-1,1);
  TH1D* h_noPU_prompt_SV_track_time_EE = new TH1D("h_noPU_prompt_SV_track_time_EE", "h_noPU_prompt_SV_track_time_EE", 200,-1,1);
  TH1D* h_noPU_nonprompt_PV_vtx_time_EE = new TH1D("h_noPU_nonprompt_PV_vtx_time_EE", "h_noPU_nonprompt_PV_vtx_time_EE", 200,-1,1);
  TH1D* h_noPU_nonprompt_PV_muon_time_EE = new TH1D("h_noPU_nonprompt_PV_muon_time_EE", "h_noPU_nonprompt_PV_muon_time_EE", 200,-1,1);
  TH1D* h_noPU_nonprompt_PV_track_time_EE = new TH1D("h_noPU_nonprompt_PV_track_time_EE", "h_noPU_nonprompt_PV_track_time_EE", 200,-1,1);
  TH1D* h_noPU_nonprompt_PU_track_time_EE = new TH1D("h_noPU_nonprompt_PU_track_time_EE", "h_noPU_nonprompt_PU_track_time_EE", 200,-1,1);
  TH1D* h_noPU_nonprompt_SV_track_time_EE = new TH1D("h_noPU_nonprompt_SV_track_time_EE", "h_noPU_nonprompt_SV_track_time_EE", 200,-1,1);




  bool flag_isMuon            = false;
  bool flag_isPFMuon          = false;
  bool flag_isGlobalMuon      = false;
  bool flag_isTrackerMuon     = false;
  bool flag_isStandAloneMuon  = false;
  bool flag_isCutBasedIdLoose = false;
  bool flag_isLooseMuon       = true;


  bool flag_muon_status=true;
  bool flag_muon_pv_dz=false,    flag_muon_pv_dxy=false;
  bool flag_track_pv_dz=true;
  bool flag_vtx_matching=true, flag_vtx_selectedLV=true;

  bool flag_muon_PVweight=false;
  bool flag_track_PVweight=false;

  bool flag_2sigma_PV=false;
  bool flag_2sigma_PV_vtx=false;

  float dtsig_cut=0, dt_cut=0;
  float dtsig_cut_vtx=0, dt_cut_vtx=0;

  // PARAM
  float muon_PVweight_cut=0.000000000000001;
  float track_PVweight_cut=0.000000000000001;
  float muon_pv_dz_cut_EB=0.5;
  float muon_pv_dz_cut_EE=0.5;
  float muon_pv_dxy_cut_EB=0.2;
  float muon_pv_dxy_cut_EE=0.2;
  float track_pv_dz_cut_EB=0.2;
  float track_pv_dz_cut_EE=0.2;

  //////////////////////
  //// PU200 prompt ////
  //////////////////////

  int n_muon_PU200_prompt_EB=0, n_muon_PU200_prompt_EE=0;
  int n_status_failed_PU200_prompt_EB=0, n_status_failed_PU200_prompt_EE=0;

  int num_muon_PU200_prompt_EB=0,  num_muon_mva_PU200_prompt_EB=0,  num_muon_mva_cut_PU200_prompt_EB=0,  num_muon_time_PU200_prompt_EB=0,  num_muon_time_err_PU200_prompt_EB=0,  num_muon_time_err_i_PU200_prompt_EB=0;
  int num_muon_PU200_prompt_EE=0,  num_muon_mva_PU200_prompt_EE=0,  num_muon_mva_cut_PU200_prompt_EE=0,  num_muon_time_PU200_prompt_EE=0,  num_muon_time_err_PU200_prompt_EE=0,  num_muon_time_err_i_PU200_prompt_EE=0;
  int num_track_PU200_prompt_EB=0, num_track_mva_PU200_prompt_EB=0, num_track_mva_cut_PU200_prompt_EB=0, num_track_time_PU200_prompt_EB=0, num_track_time_err_PU200_prompt_EB=0, num_track_time_err_i_PU200_prompt_EB=0;
  int num_track_PU200_prompt_EE=0, num_track_mva_PU200_prompt_EE=0, num_track_mva_cut_PU200_prompt_EE=0, num_track_time_PU200_prompt_EE=0, num_track_time_err_PU200_prompt_EE=0, num_track_time_err_i_PU200_prompt_EE=0;
  int num_vtx_PU200_prompt_EB=0,   num_vtx_time_PU200_prompt_EB=0,  num_vtx_time_err_PU200_prompt_EB=0;
  int num_vtx_PU200_prompt_EE=0,   num_vtx_time_PU200_prompt_EE=0,  num_vtx_time_err_PU200_prompt_EE=0;
  int num_both_muon_track_time_err_PU200_prompt_EB=0, num_both_muon_track_time_err_PU200_prompt_EE=0;
  int num_both_vtx_track_time_err_PU200_prompt_EB=0, num_both_vtx_track_time_err_PU200_prompt_EE=0;

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

		// calculate fraction of observables
		num_muon_PU200_prompt_EB++;
		if(muon_mva_->at(im)>0)          num_muon_mva_PU200_prompt_EB++;
		if(muon_mva_->at(im)>0.5)        num_muon_mva_cut_PU200_prompt_EB++;
		if(muon_time_->at(im)!=0)        num_muon_time_PU200_prompt_EB++;
		if(muon_time_err_->at(im)!=-1)   num_muon_time_err_PU200_prompt_EB++;
		if(muon_time_err_i_->at(im)!=-1) num_muon_time_err_i_PU200_prompt_EB++;
		if(im==0) {
		  num_vtx_PU200_prompt_EB++;
		  if(vtx_time_->at(0)!=0)     num_vtx_time_PU200_prompt_EB++;
		  if(vtx_time_err_->at(0)!=0) num_vtx_time_err_PU200_prompt_EB++;
		}

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {
		  int track_type = -999;
		  if(simvtx_reco_==-999) track_type = 999;
		  else track_type = track_type_->at(im).at(it);

		  if(flag_track_pv_dz) {
			//if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) > 0.1) continue;
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EB) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // calculate fraction of observables
		  num_track_PU200_prompt_EB++;
		  if(track_mva_->at(im).at(it)>0) num_track_mva_PU200_prompt_EB++;
		  if(track_mva_->at(im).at(it)>0.5) num_track_mva_cut_PU200_prompt_EB++;
		  if(track_time_->at(im).at(it)!=0) num_track_time_PU200_prompt_EB++;
		  if(track_time_err_->at(im).at(it)!=-1) num_track_time_err_PU200_prompt_EB++;
		  if(track_time_err_i_->at(im).at(it)!=-1) num_track_time_err_i_PU200_prompt_EB++;
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_muon_track_time_err_PU200_prompt_EB++;
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_vtx_track_time_err_PU200_prompt_EB++;

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }

		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it));
		  }

		  if(track_type==0) {
			//h_PU200_prompt_PV_EB->Fill(dtsig_cut);
			// (muon, track)
			if(dtsig_cut!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut>2.0) {
			      h_PU200_prompt_PV_EB->Fill(dtsig_cut);
			      h_PU200_prompt_bx_PV_EB->Fill(track_bx_->at(im).at(it));
			      h_PU200_prompt_evtId_PV_EB->Fill(track_evtId_->at(im).at(it));
			      h_PU200_prompt_PVweight_PV_EB->Fill(track_PVweight_->at(im).at(it));
			      h_PU200_prompt_simvtx_bx_PV_EB->Fill(simvtx_bx_);
			      h_PU200_prompt_simvtx_evtId_PV_EB->Fill(simvtx_evtId_);
//			      h_PU200_prompt_PV_EB_vtx->Fill(dtsig_cut_vtx);
				}
			  }
			  else {
			    h_PU200_prompt_PV_EB->Fill(dtsig_cut);
			    h_PU200_prompt_bx_PV_EB->Fill(track_bx_->at(im).at(it));
			    h_PU200_prompt_evtId_PV_EB->Fill(track_evtId_->at(im).at(it));
			    h_PU200_prompt_PVweight_PV_EB->Fill(track_PVweight_->at(im).at(it));
			    h_PU200_prompt_simvtx_bx_PV_EB->Fill(simvtx_bx_);
			    h_PU200_prompt_simvtx_evtId_PV_EB->Fill(simvtx_evtId_);
//			    h_PU200_prompt_PV_EB_vtx->Fill(dtsig_cut_vtx);
			  }
			}
			else {
			  h_PU200_prompt_no_tErr_EB->Fill(dtsig_cut);
			  h_PU200_prompt_no_tErr_PV_EB->Fill(dtsig_cut);
			}

			// (PV, track)
			if(dtsig_cut_vtx!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut_vtx>2.0) {
			      h_PU200_prompt_PV_EB_vtx->Fill(dtsig_cut_vtx);
			      h_PU200_prompt_bx_PV_EB_vtx->Fill(track_bx_->at(im).at(it));
			      h_PU200_prompt_evtId_PV_EB_vtx->Fill(track_evtId_->at(im).at(it));
			      h_PU200_prompt_PVweight_PV_EB_vtx->Fill(track_PVweight_->at(im).at(it));
			      h_PU200_prompt_simvtx_bx_PV_EB_vtx->Fill(simvtx_bx_);
			      h_PU200_prompt_simvtx_evtId_PV_EB_vtx->Fill(simvtx_evtId_);
				}
			  }
			  else {
			    h_PU200_prompt_PV_EB_vtx->Fill(dtsig_cut_vtx);
			    h_PU200_prompt_bx_PV_EB_vtx->Fill(track_bx_->at(im).at(it));
			    h_PU200_prompt_evtId_PV_EB_vtx->Fill(track_evtId_->at(im).at(it));
			    h_PU200_prompt_PVweight_PV_EB_vtx->Fill(track_PVweight_->at(im).at(it));
			    h_PU200_prompt_simvtx_bx_PV_EB_vtx->Fill(simvtx_bx_);
			    h_PU200_prompt_simvtx_evtId_PV_EB_vtx->Fill(simvtx_evtId_);
			  }
			}
			else {
			  h_PU200_prompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
			  h_PU200_prompt_no_tErr_PV_EB_vtx->Fill(dtsig_cut_vtx);
			}

			// 2D
			if(dtsig_cut!=0 && dtsig_cut_vtx!=0) h_PU200_prompt_PV_EB_2D->Fill(dtsig_cut, dtsig_cut_vtx);
		  }
		  else if(track_type==3) {
			//h_PU200_prompt_SV_EB->Fill(dtsig_cut);
			if(dtsig_cut!=0) h_PU200_prompt_SV_EB->Fill(dtsig_cut);
			else h_PU200_prompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_prompt_SV_EB_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_prompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==1) {
			//h_PU200_prompt_PU_EB->Fill(dtsig_cut);
			if(dtsig_cut!=0) h_PU200_prompt_PU_EB->Fill(dtsig_cut);
			else h_PU200_prompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_prompt_PU_EB_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_prompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==2) {
			//h_PU200_prompt_fake_EB->Fill(dtsig_cut);
			if(dtsig_cut!=0) h_PU200_prompt_fake_EB->Fill(dtsig_cut);
			else h_PU200_prompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_prompt_fake_EB_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_prompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  //else cout << "Something wrong with PU200 prompt EB" << endl;
		  else {
			if(dtsig_cut!=0) h_PU200_prompt_not_vtx_matched_EB->Fill(dtsig_cut);
			else h_PU200_prompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_prompt_not_vtx_matched_EB_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_prompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  

		  // tsim test
		    // track
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_prompt_PV_track_time_EB->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==1 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_prompt_PU_track_time_EB->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==3 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_prompt_SV_track_time_EB->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		    // muon
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_prompt_PV_muon_time_EB->Fill(muon_time_->at(im) - muon_time_sim_->at(im));
		  }
		    // vertex
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_prompt_PV_vtx_time_EB->Fill(vtx_time_->at(im) - simvtx_time_*1e9);
		  }

		  dtsig_cut=0, dtsig_cut_vtx=0;
		  dt_cut=0, dt_cut_vtx=0;

		}


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

		// calculate fraction of observables
        num_muon_PU200_prompt_EE++;
        if(muon_mva_->at(im)>0)          num_muon_mva_PU200_prompt_EE++;
        if(muon_mva_->at(im)>0.5)        num_muon_mva_cut_PU200_prompt_EE++;
        if(muon_time_->at(im)!=0)        num_muon_time_PU200_prompt_EE++;
        if(muon_time_err_->at(im)!=-1)   num_muon_time_err_PU200_prompt_EE++;
        if(muon_time_err_i_->at(im)!=-1) num_muon_time_err_i_PU200_prompt_EE++;
        if(im==0) {
          num_vtx_PU200_prompt_EE++;
          if(vtx_time_->at(0)!=0)     num_vtx_time_PU200_prompt_EE++;
          if(vtx_time_err_->at(0)!=0) num_vtx_time_err_PU200_prompt_EE++;
        }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {
		  int track_type = -999;
		  if(simvtx_reco_==-999) track_type = 999;
		  else track_type = track_type_->at(im).at(it);

		  if(flag_track_pv_dz) {
			//if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) > 0.2) continue;
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EE) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // calculate fraction of observables
		  num_track_PU200_prompt_EE++;
          if(track_mva_->at(im).at(it)>0) num_track_mva_PU200_prompt_EE++;
          if(track_mva_->at(im).at(it)>0.5) num_track_mva_cut_PU200_prompt_EE++;
          if(track_time_->at(im).at(it)!=0) num_track_time_PU200_prompt_EE++;
          if(track_time_err_->at(im).at(it)!=-1) num_track_time_err_PU200_prompt_EE++;
          if(track_time_err_i_->at(im).at(it)!=-1) num_track_time_err_i_PU200_prompt_EE++;
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_muon_track_time_err_PU200_prompt_EE++;
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_vtx_track_time_err_PU200_prompt_EE++;

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dt_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it));
          }

		  if(track_type==0) {
			//h_PU200_prompt_PV_EE->Fill(dtsig_cut);
			// (muon, track)
			if(dtsig_cut!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut>2.0) {
			      h_PU200_prompt_PV_EE->Fill(dtsig_cut);
			      h_PU200_prompt_bx_PV_EE->Fill(track_bx_->at(im).at(it));
			      h_PU200_prompt_evtId_PV_EE->Fill(track_evtId_->at(im).at(it));
			      h_PU200_prompt_PVweight_PV_EE->Fill(track_PVweight_->at(im).at(it));
			      h_PU200_prompt_simvtx_bx_PV_EE->Fill(simvtx_bx_);
			      h_PU200_prompt_simvtx_evtId_PV_EE->Fill(simvtx_evtId_);
				}    
			  }
			  else {
			    h_PU200_prompt_PV_EE->Fill(dtsig_cut);
			    h_PU200_prompt_bx_PV_EE->Fill(track_bx_->at(im).at(it));
			    h_PU200_prompt_evtId_PV_EE->Fill(track_evtId_->at(im).at(it));
			    h_PU200_prompt_PVweight_PV_EE->Fill(track_PVweight_->at(im).at(it));
			    h_PU200_prompt_simvtx_bx_PV_EE->Fill(simvtx_bx_);
			    h_PU200_prompt_simvtx_evtId_PV_EE->Fill(simvtx_evtId_);
			  }
			}
			else h_PU200_prompt_no_tErr_EE->Fill(dtsig_cut);

			// (PV, track)
			if(dtsig_cut_vtx!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut_vtx>2.0) {
			      h_PU200_prompt_PV_EE_vtx->Fill(dtsig_cut_vtx);
			      h_PU200_prompt_bx_PV_EE_vtx->Fill(track_bx_->at(im).at(it));
			      h_PU200_prompt_evtId_PV_EE_vtx->Fill(track_evtId_->at(im).at(it));
			      h_PU200_prompt_PVweight_PV_EE_vtx->Fill(track_PVweight_->at(im).at(it));
			      h_PU200_prompt_simvtx_bx_PV_EE_vtx->Fill(simvtx_bx_);
			      h_PU200_prompt_simvtx_evtId_PV_EE_vtx->Fill(simvtx_evtId_);
				}    
			  }
			  else {
			    h_PU200_prompt_PV_EE_vtx->Fill(dtsig_cut_vtx);
			    h_PU200_prompt_bx_PV_EE_vtx->Fill(track_bx_->at(im).at(it));
			    h_PU200_prompt_evtId_PV_EE_vtx->Fill(track_evtId_->at(im).at(it));
			    h_PU200_prompt_PVweight_PV_EE_vtx->Fill(track_PVweight_->at(im).at(it));
			    h_PU200_prompt_simvtx_bx_PV_EE_vtx->Fill(simvtx_bx_);
			    h_PU200_prompt_simvtx_evtId_PV_EE_vtx->Fill(simvtx_evtId_);
			  }
			}
			else h_PU200_prompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);

			// 2D
			if(dtsig_cut!=0 && dtsig_cut_vtx!=0) h_PU200_prompt_PV_EE_2D->Fill(dtsig_cut, dtsig_cut_vtx);
		  }
		  else if(track_type==3) {
			//h_PU200_prompt_SV_EE->Fill(dtsig_cut);
			if(dtsig_cut!=0) h_PU200_prompt_SV_EE->Fill(dtsig_cut);
			else h_PU200_prompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_prompt_SV_EE_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_prompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==1) {
			//h_PU200_prompt_PU_EE->Fill(dtsig_cut);
			if(dtsig_cut!=0) h_PU200_prompt_PU_EE->Fill(dtsig_cut);
			else h_PU200_prompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_prompt_PU_EE_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_prompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==2) {
			//h_PU200_prompt_fake_EE->Fill(dtsig_cut);
			if(dtsig_cut!=0) h_PU200_prompt_fake_EE->Fill(dtsig_cut);
			else h_PU200_prompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_prompt_fake_EE_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_prompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }
		  //else cout << "Something wrong with PU200 prompt EE" << endl;
		  else {
			if(dtsig_cut!=0) h_PU200_prompt_not_vtx_matched_EE->Fill(dtsig_cut);
			else h_PU200_prompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_prompt_not_vtx_matched_EE_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_prompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }


		  // tsim test
		    // track
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_prompt_PV_track_time_EE->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==1 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_prompt_PU_track_time_EE->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==3 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_prompt_SV_track_time_EE->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  // muon
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_prompt_PV_muon_time_EE->Fill(muon_time_->at(im) - muon_time_sim_->at(im));
		  }
		  // vertex
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_prompt_PV_vtx_time_EE->Fill(vtx_time_->at(im) - simvtx_time_*1e9);
		  }

		  dtsig_cut=0, dtsig_cut_vtx=0;
		  dt_cut=0, dt_cut_vtx=0;
		}
	  }
	} // End of muon loop
  } // End of event loop


  /////////////////////////
  //// PU200 nonprompt ////
  /////////////////////////
  
  int n_muon_PU200_nonprompt_EB=0, n_muon_PU200_nonprompt_EE=0;
  int n_status_failed_PU200_nonprompt_EB=0, n_status_failed_PU200_nonprompt_EE=0;

  int num_muon_PU200_nonprompt_EB=0,  num_muon_mva_PU200_nonprompt_EB=0,  num_muon_mva_cut_PU200_nonprompt_EB=0,  num_muon_time_PU200_nonprompt_EB=0,  num_muon_time_err_PU200_nonprompt_EB=0,  num_muon_time_err_i_PU200_nonprompt_EB=0;
  int num_muon_PU200_nonprompt_EE=0,  num_muon_mva_PU200_nonprompt_EE=0,  num_muon_mva_cut_PU200_nonprompt_EE=0,  num_muon_time_PU200_nonprompt_EE=0,  num_muon_time_err_PU200_nonprompt_EE=0,  num_muon_time_err_i_PU200_nonprompt_EE=0;
  int num_track_PU200_nonprompt_EB=0, num_track_mva_PU200_nonprompt_EB=0, num_track_mva_cut_PU200_nonprompt_EB=0, num_track_time_PU200_nonprompt_EB=0, num_track_time_err_PU200_nonprompt_EB=0, num_track_time_err_i_PU200_nonprompt_EB=0;
  int num_track_PU200_nonprompt_EE=0, num_track_mva_PU200_nonprompt_EE=0, num_track_mva_cut_PU200_nonprompt_EE=0, num_track_time_PU200_nonprompt_EE=0, num_track_time_err_PU200_nonprompt_EE=0, num_track_time_err_i_PU200_nonprompt_EE=0;
  int num_vtx_PU200_nonprompt_EB=0,   num_vtx_time_PU200_nonprompt_EB=0,  num_vtx_time_err_PU200_nonprompt_EB=0;
  int num_vtx_PU200_nonprompt_EE=0,   num_vtx_time_PU200_nonprompt_EE=0,  num_vtx_time_err_PU200_nonprompt_EE=0;
  int num_both_muon_track_time_err_PU200_nonprompt_EB=0, num_both_muon_track_time_err_PU200_nonprompt_EE=0;
  int num_both_vtx_track_time_err_PU200_nonprompt_EB=0, num_both_vtx_track_time_err_PU200_nonprompt_EE=0;

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

		// calculate fraction of observables
        num_muon_PU200_nonprompt_EB++;
        if(muon_mva_->at(im)>0)          num_muon_mva_PU200_nonprompt_EB++;
        if(muon_mva_->at(im)>0.5)        num_muon_mva_cut_PU200_nonprompt_EB++;
        if(muon_time_->at(im)!=0)        num_muon_time_PU200_nonprompt_EB++;
        if(muon_time_err_->at(im)!=-1)   num_muon_time_err_PU200_nonprompt_EB++;
        if(muon_time_err_i_->at(im)!=-1) num_muon_time_err_i_PU200_nonprompt_EB++;
        if(im==0) {
          num_vtx_PU200_nonprompt_EB++;
          if(vtx_time_->at(0)!=0)     num_vtx_time_PU200_nonprompt_EB++;
          if(vtx_time_err_->at(0)!=0) num_vtx_time_err_PU200_nonprompt_EB++;
        }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {
		  int track_type = -999;
		  if(simvtx_reco_==-999) track_type = 999;
		  else track_type = track_type_->at(im).at(it);

		  if(flag_track_pv_dz) {
			//if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) > 0.1) continue;
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EB) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // calculate fraction of observables
          num_track_PU200_nonprompt_EB++;
          if(track_mva_->at(im).at(it)>0) num_track_mva_PU200_nonprompt_EB++;
          if(track_mva_->at(im).at(it)>0.5) num_track_mva_cut_PU200_nonprompt_EB++;
          if(track_time_->at(im).at(it)!=0) num_track_time_PU200_nonprompt_EB++;
          if(track_time_err_->at(im).at(it)!=-1) num_track_time_err_PU200_nonprompt_EB++;
          if(track_time_err_i_->at(im).at(it)!=-1) num_track_time_err_i_PU200_nonprompt_EB++;
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_muon_track_time_err_PU200_nonprompt_EB++;
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_vtx_track_time_err_PU200_nonprompt_EB++;

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dt_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it));
          }

		  if(track_type==0) {
			// (muon, track)
			if(dtsig_cut!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut>2.0) {
			      h_PU200_nonprompt_PV_EB->Fill(dtsig_cut);
			      h_PU200_nonprompt_bx_PV_EB->Fill(track_bx_->at(im).at(it));
			      h_PU200_nonprompt_evtId_PV_EB->Fill(track_evtId_->at(im).at(it));
			      h_PU200_nonprompt_PVweight_PV_EB->Fill(track_PVweight_->at(im).at(it));
			      h_PU200_nonprompt_simvtx_bx_PV_EB->Fill(simvtx_bx_);
			      h_PU200_nonprompt_simvtx_evtId_PV_EB->Fill(simvtx_evtId_);
				}
			  }
			  else {
			    h_PU200_nonprompt_PV_EB->Fill(dtsig_cut);
			    h_PU200_nonprompt_bx_PV_EB->Fill(track_bx_->at(im).at(it));
			    h_PU200_nonprompt_evtId_PV_EB->Fill(track_evtId_->at(im).at(it));
			    h_PU200_nonprompt_PVweight_PV_EB->Fill(track_PVweight_->at(im).at(it));
			    h_PU200_nonprompt_simvtx_bx_PV_EB->Fill(simvtx_bx_);
			    h_PU200_nonprompt_simvtx_evtId_PV_EB->Fill(simvtx_evtId_);
			  }
			}
			else h_PU200_nonprompt_no_tErr_EB->Fill(dtsig_cut);

			// (PV, track)
			if(dtsig_cut_vtx!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut_vtx>2.0) {
			      h_PU200_nonprompt_PV_EB_vtx->Fill(dtsig_cut_vtx);
			      h_PU200_nonprompt_bx_PV_EB_vtx->Fill(track_bx_->at(im).at(it));
			      h_PU200_nonprompt_evtId_PV_EB_vtx->Fill(track_evtId_->at(im).at(it));
			      h_PU200_nonprompt_PVweight_PV_EB_vtx->Fill(track_PVweight_->at(im).at(it));
			      h_PU200_nonprompt_simvtx_bx_PV_EB_vtx->Fill(simvtx_bx_);
			      h_PU200_nonprompt_simvtx_evtId_PV_EB_vtx->Fill(simvtx_evtId_);
				}
			  }
			  else {
			    h_PU200_nonprompt_PV_EB_vtx->Fill(dtsig_cut_vtx);
			    h_PU200_nonprompt_bx_PV_EB_vtx->Fill(track_bx_->at(im).at(it));
			    h_PU200_nonprompt_evtId_PV_EB_vtx->Fill(track_evtId_->at(im).at(it));
			    h_PU200_nonprompt_PVweight_PV_EB_vtx->Fill(track_PVweight_->at(im).at(it));
			    h_PU200_nonprompt_simvtx_bx_PV_EB_vtx->Fill(simvtx_bx_);
			    h_PU200_nonprompt_simvtx_evtId_PV_EB_vtx->Fill(simvtx_evtId_);
			  }
			}
			else h_PU200_nonprompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);

			// 2D
			if(dtsig_cut!=0 && dtsig_cut_vtx!=0) h_PU200_nonprompt_PV_EB_2D->Fill(dtsig_cut, dtsig_cut_vtx);
		  }
		  else if(track_type==3) {
			if(dtsig_cut!=0) h_PU200_nonprompt_SV_EB->Fill(dtsig_cut);
			else h_PU200_nonprompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_nonprompt_SV_EB_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_nonprompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==1) {
			if(dtsig_cut!=0) h_PU200_nonprompt_PU_EB->Fill(dtsig_cut);
			else h_PU200_nonprompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_nonprompt_PU_EB_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_nonprompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==2) {
			if(dtsig_cut!=0) h_PU200_nonprompt_fake_EB->Fill(dtsig_cut);
			else h_PU200_nonprompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_nonprompt_fake_EB_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_nonprompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  //else cout << "Something wrong with PU200 nonprompt EB" << endl;
		  else {
			if(dtsig_cut!=0) h_PU200_nonprompt_not_vtx_matched_EB->Fill(dtsig_cut);
			else h_PU200_nonprompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_nonprompt_not_vtx_matched_EB_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_nonprompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }


		  // tsim test
		    // track
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_nonprompt_PV_track_time_EB->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==1 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_nonprompt_PU_track_time_EB->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==3 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_nonprompt_SV_track_time_EB->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		    // muon
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_nonprompt_PV_muon_time_EB->Fill(muon_time_->at(im) - muon_time_sim_->at(im));
		  }
		    // vertex
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_nonprompt_PV_vtx_time_EB->Fill(vtx_time_->at(im) - simvtx_time_*1e9);
		  }

		  dtsig_cut=0, dtsig_cut_vtx=0;
		  dt_cut=0, dt_cut_vtx=0;
		}
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

		// calculate fraction of observables
        num_muon_PU200_nonprompt_EE++;
        if(muon_mva_->at(im)>0)          num_muon_mva_PU200_nonprompt_EE++;
        if(muon_mva_->at(im)>0.5)        num_muon_mva_cut_PU200_nonprompt_EE++;
        if(muon_time_->at(im)!=0)        num_muon_time_PU200_nonprompt_EE++;
        if(muon_time_err_->at(im)!=-1)   num_muon_time_err_PU200_nonprompt_EE++;
        if(muon_time_err_i_->at(im)!=-1) num_muon_time_err_i_PU200_nonprompt_EE++;
        if(im==0) {
          num_vtx_PU200_nonprompt_EE++;
          if(vtx_time_->at(0)!=0)     num_vtx_time_PU200_nonprompt_EE++;
          if(vtx_time_err_->at(0)!=0) num_vtx_time_err_PU200_nonprompt_EE++;
        }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {
		  int track_type = -999;
		  if(simvtx_reco_==-999) track_type = 999;
		  else track_type = track_type_->at(im).at(it);

		  if(flag_track_pv_dz) {
			//if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) > 0.2) continue;
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EE) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // calculate fraction of observables
          num_track_PU200_nonprompt_EE++;
          if(track_mva_->at(im).at(it)>0) num_track_mva_PU200_nonprompt_EE++;
          if(track_mva_->at(im).at(it)>0.5) num_track_mva_cut_PU200_nonprompt_EE++;
          if(track_time_->at(im).at(it)!=0) num_track_time_PU200_nonprompt_EE++;
          if(track_time_err_->at(im).at(it)!=-1) num_track_time_err_PU200_nonprompt_EE++;
          if(track_time_err_i_->at(im).at(it)!=-1) num_track_time_err_i_PU200_nonprompt_EE++;
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_muon_track_time_err_PU200_nonprompt_EE++;
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_vtx_track_time_err_PU200_nonprompt_EE++;

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dt_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it));
          }

		  if(track_type==0) {
			// (muon, track)
			if(dtsig_cut!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut>2.0) {
			      h_PU200_nonprompt_PV_EE->Fill(dtsig_cut);
			      h_PU200_nonprompt_bx_PV_EE->Fill(track_bx_->at(im).at(it));
			      h_PU200_nonprompt_evtId_PV_EE->Fill(track_evtId_->at(im).at(it));
			      h_PU200_nonprompt_PVweight_PV_EE->Fill(track_PVweight_->at(im).at(it));
			      h_PU200_nonprompt_simvtx_bx_PV_EE->Fill(simvtx_bx_);
			      h_PU200_nonprompt_simvtx_evtId_PV_EE->Fill(simvtx_evtId_);
				}
			  }
			  else {
			    h_PU200_nonprompt_PV_EE->Fill(dtsig_cut);
			    h_PU200_nonprompt_bx_PV_EE->Fill(track_bx_->at(im).at(it));
			    h_PU200_nonprompt_evtId_PV_EE->Fill(track_evtId_->at(im).at(it));
			    h_PU200_nonprompt_PVweight_PV_EE->Fill(track_PVweight_->at(im).at(it));
			    h_PU200_nonprompt_simvtx_bx_PV_EE->Fill(simvtx_bx_);
			    h_PU200_nonprompt_simvtx_evtId_PV_EE->Fill(simvtx_evtId_);
			  }
			}
			else h_PU200_nonprompt_no_tErr_EE->Fill(dtsig_cut);

			// (PV, track)
			if(dtsig_cut_vtx!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut_vtx>2.0) {
			      h_PU200_nonprompt_PV_EE_vtx->Fill(dtsig_cut_vtx);
			      h_PU200_nonprompt_bx_PV_EE_vtx->Fill(track_bx_->at(im).at(it));
			      h_PU200_nonprompt_evtId_PV_EE_vtx->Fill(track_evtId_->at(im).at(it));
			      h_PU200_nonprompt_PVweight_PV_EE_vtx->Fill(track_PVweight_->at(im).at(it));
			      h_PU200_nonprompt_simvtx_bx_PV_EE_vtx->Fill(simvtx_bx_);
			      h_PU200_nonprompt_simvtx_evtId_PV_EE_vtx->Fill(simvtx_evtId_);
				}
			  }
			  else {
			    h_PU200_nonprompt_PV_EE_vtx->Fill(dtsig_cut_vtx);
			    h_PU200_nonprompt_bx_PV_EE_vtx->Fill(track_bx_->at(im).at(it));
			    h_PU200_nonprompt_evtId_PV_EE_vtx->Fill(track_evtId_->at(im).at(it));
			    h_PU200_nonprompt_PVweight_PV_EE_vtx->Fill(track_PVweight_->at(im).at(it));
			    h_PU200_nonprompt_simvtx_bx_PV_EE_vtx->Fill(simvtx_bx_);
			    h_PU200_nonprompt_simvtx_evtId_PV_EE_vtx->Fill(simvtx_evtId_);
			  }
			}
			else h_PU200_nonprompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);

			// 2D
			if(dtsig_cut!=0 && dtsig_cut_vtx!=0) h_PU200_nonprompt_PV_EE_2D->Fill(dtsig_cut, dtsig_cut_vtx);
		  }
		  else if(track_type==3) {
			if(dtsig_cut!=0) h_PU200_nonprompt_SV_EE->Fill(dtsig_cut);
			else h_PU200_nonprompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_nonprompt_SV_EE_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_nonprompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==1) {
			if(dtsig_cut!=0) h_PU200_nonprompt_PU_EE->Fill(dtsig_cut);
			else h_PU200_nonprompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_nonprompt_PU_EE_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_nonprompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==2) {
			if(dtsig_cut!=0) h_PU200_nonprompt_fake_EE->Fill(dtsig_cut);
			else h_PU200_nonprompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_nonprompt_fake_EE_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_nonprompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }
		  //else cout << "Something wrong with PU200 nonprompt EE" << endl;
		  else {
			if(dtsig_cut!=0) h_PU200_nonprompt_not_vtx_matched_EE->Fill(dtsig_cut);
			else h_PU200_nonprompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_PU200_nonprompt_not_vtx_matched_EE_vtx->Fill(dtsig_cut_vtx);
			else h_PU200_nonprompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }


		  // tsim test
		    // track
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_nonprompt_PV_track_time_EE->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==1 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_nonprompt_PU_track_time_EE->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==3 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_nonprompt_SV_track_time_EE->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		    // muon
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_nonprompt_PV_muon_time_EE->Fill(muon_time_->at(im) - muon_time_sim_->at(im));
		  }
		    // vertex
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_PU200_nonprompt_PV_vtx_time_EE->Fill(vtx_time_->at(im) - simvtx_time_*1e9);
		  }

		  dtsig_cut=0, dtsig_cut_vtx=0;
		  dt_cut=0, dt_cut_vtx=0;
		}
	  }
	} // End of muon loop
  } // End of event loop
	

  //////////////////////
  //// noPU prompt ////
  //////////////////////
  
  int n_muon_noPU_prompt_EB=0, n_muon_noPU_prompt_EE=0;
  int n_status_failed_noPU_prompt_EB=0, n_status_failed_noPU_prompt_EE=0;

  int num_muon_noPU_prompt_EB=0,  num_muon_mva_noPU_prompt_EB=0,  num_muon_mva_cut_noPU_prompt_EB=0,  num_muon_time_noPU_prompt_EB=0,  num_muon_time_err_noPU_prompt_EB=0,  num_muon_time_err_i_noPU_prompt_EB=0;
  int num_muon_noPU_prompt_EE=0,  num_muon_mva_noPU_prompt_EE=0,  num_muon_mva_cut_noPU_prompt_EE=0,  num_muon_time_noPU_prompt_EE=0,  num_muon_time_err_noPU_prompt_EE=0,  num_muon_time_err_i_noPU_prompt_EE=0;
  int num_track_noPU_prompt_EB=0, num_track_mva_noPU_prompt_EB=0, num_track_mva_cut_noPU_prompt_EB=0, num_track_time_noPU_prompt_EB=0, num_track_time_err_noPU_prompt_EB=0, num_track_time_err_i_noPU_prompt_EB=0;
  int num_track_noPU_prompt_EE=0, num_track_mva_noPU_prompt_EE=0, num_track_mva_cut_noPU_prompt_EE=0, num_track_time_noPU_prompt_EE=0, num_track_time_err_noPU_prompt_EE=0, num_track_time_err_i_noPU_prompt_EE=0;
  int num_vtx_noPU_prompt_EB=0,   num_vtx_time_noPU_prompt_EB=0,  num_vtx_time_err_noPU_prompt_EB=0;
  int num_vtx_noPU_prompt_EE=0,   num_vtx_time_noPU_prompt_EE=0,  num_vtx_time_err_noPU_prompt_EE=0;
  int num_both_muon_track_time_err_noPU_prompt_EB=0, num_both_muon_track_time_err_noPU_prompt_EE=0;
  int num_both_vtx_track_time_err_noPU_prompt_EB=0, num_both_vtx_track_time_err_noPU_prompt_EE=0;

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

		// calculate fraction of observables
        num_muon_noPU_prompt_EB++;
        if(muon_mva_->at(im)>0)          num_muon_mva_noPU_prompt_EB++;
        if(muon_mva_->at(im)>0.5)        num_muon_mva_cut_noPU_prompt_EB++;
        if(muon_time_->at(im)!=0)        num_muon_time_noPU_prompt_EB++;
        if(muon_time_err_->at(im)!=-1)   num_muon_time_err_noPU_prompt_EB++;
        if(muon_time_err_i_->at(im)!=-1) num_muon_time_err_i_noPU_prompt_EB++;
        if(im==0) {
          num_vtx_noPU_prompt_EB++;
          if(vtx_time_->at(0)!=0)     num_vtx_time_noPU_prompt_EB++;
          if(vtx_time_err_->at(0)!=0) num_vtx_time_err_noPU_prompt_EB++;
        }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {
		  int track_type = -999;
		  if(simvtx_reco_==-999) track_type = 999;
		  else track_type = track_type_->at(im).at(it);

		  if(flag_track_pv_dz) {
			//if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) > 0.1) continue;
			if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EB) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // calculate fraction of observables
          num_track_noPU_prompt_EB++;
          if(track_mva_->at(im).at(it)>0) num_track_mva_noPU_prompt_EB++;
          if(track_mva_->at(im).at(it)>0.5) num_track_mva_cut_noPU_prompt_EB++;
          if(track_time_->at(im).at(it)!=0) num_track_time_noPU_prompt_EB++;
          if(track_time_err_->at(im).at(it)!=-1) num_track_time_err_noPU_prompt_EB++;
          if(track_time_err_i_->at(im).at(it)!=-1) num_track_time_err_i_noPU_prompt_EB++;
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_muon_track_time_err_noPU_prompt_EB++;
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_vtx_track_time_err_noPU_prompt_EB++;

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dt_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it));
          }

		  if(track_type==0) {
			// (muon, track)
			if(dtsig_cut!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut>2.0) {
			      h_noPU_prompt_PV_EB->Fill(dtsig_cut);
			      h_noPU_prompt_bx_PV_EB->Fill(track_bx_->at(im).at(it));
			      h_noPU_prompt_evtId_PV_EB->Fill(track_evtId_->at(im).at(it));
			      h_noPU_prompt_PVweight_PV_EB->Fill(track_PVweight_->at(im).at(it));
			      h_noPU_prompt_simvtx_bx_PV_EB->Fill(simvtx_bx_);
			      h_noPU_prompt_simvtx_evtId_PV_EB->Fill(simvtx_evtId_);
				}    
			  }
			  else {
			    h_noPU_prompt_PV_EB->Fill(dtsig_cut);
			    h_noPU_prompt_bx_PV_EB->Fill(track_bx_->at(im).at(it));
			    h_noPU_prompt_evtId_PV_EB->Fill(track_evtId_->at(im).at(it));
			    h_noPU_prompt_PVweight_PV_EB->Fill(track_PVweight_->at(im).at(it));
			    h_noPU_prompt_simvtx_bx_PV_EB->Fill(simvtx_bx_);
			    h_noPU_prompt_simvtx_evtId_PV_EB->Fill(simvtx_evtId_);
			  }
			}
			else h_noPU_prompt_no_tErr_EB->Fill(dtsig_cut);

			// (PV, track)
			if(dtsig_cut_vtx!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut_vtx>2.0) {
			      h_noPU_prompt_PV_EB_vtx->Fill(dtsig_cut_vtx);
			      h_noPU_prompt_bx_PV_EB_vtx->Fill(track_bx_->at(im).at(it));
			      h_noPU_prompt_evtId_PV_EB_vtx->Fill(track_evtId_->at(im).at(it));
			      h_noPU_prompt_PVweight_PV_EB_vtx->Fill(track_PVweight_->at(im).at(it));
			      h_noPU_prompt_simvtx_bx_PV_EB_vtx->Fill(simvtx_bx_);
			      h_noPU_prompt_simvtx_evtId_PV_EB_vtx->Fill(simvtx_evtId_);
				}    
			  }
			  else {
			    h_noPU_prompt_PV_EB_vtx->Fill(dtsig_cut_vtx);
			    h_noPU_prompt_bx_PV_EB_vtx->Fill(track_bx_->at(im).at(it));
			    h_noPU_prompt_evtId_PV_EB_vtx->Fill(track_evtId_->at(im).at(it));
			    h_noPU_prompt_PVweight_PV_EB_vtx->Fill(track_PVweight_->at(im).at(it));
			    h_noPU_prompt_simvtx_bx_PV_EB_vtx->Fill(simvtx_bx_);
			    h_noPU_prompt_simvtx_evtId_PV_EB_vtx->Fill(simvtx_evtId_);
			  }
			}
			else h_noPU_prompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);

			// 2D
			if(dtsig_cut!=0 && dtsig_cut_vtx!=0) h_noPU_prompt_PV_EB_2D->Fill(dtsig_cut, dtsig_cut_vtx);
		  }
		  else if(track_type==3) {
			if(dtsig_cut!=0) h_noPU_prompt_SV_EB->Fill(dtsig_cut);
			else h_noPU_prompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_prompt_SV_EB_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_prompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==1) {
			if(dtsig_cut!=0) h_noPU_prompt_PU_EB->Fill(dtsig_cut);
			else h_noPU_prompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_prompt_PU_EB_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_prompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==2) {
			if(dtsig_cut!=0) h_noPU_prompt_fake_EB->Fill(dtsig_cut);
			else h_noPU_prompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_prompt_fake_EB_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_prompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  //else cout << "Something wrong with noPU prompt EB" << endl;
		  else {
			if(dtsig_cut!=0) h_noPU_prompt_not_vtx_matched_EB->Fill(dtsig_cut);
			else h_noPU_prompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_prompt_not_vtx_matched_EB_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_prompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }


		  // tsim test
		    // track
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_prompt_PV_track_time_EB->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==1 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_prompt_PU_track_time_EB->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==3 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_prompt_SV_track_time_EB->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		    // muon
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_prompt_PV_muon_time_EB->Fill(muon_time_->at(im) - muon_time_sim_->at(im));
		  }
		    // vertex
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_prompt_PV_vtx_time_EB->Fill(vtx_time_->at(im) - simvtx_time_*1e9);
		  }
    
		  dtsig_cut=0, dtsig_cut_vtx=0;
		  dt_cut=0, dt_cut_vtx=0;
		}
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

		// calculate fraction of observables
        num_muon_noPU_prompt_EE++;
        if(muon_mva_->at(im)>0)          num_muon_mva_noPU_prompt_EE++;
        if(muon_mva_->at(im)>0.5)        num_muon_mva_cut_noPU_prompt_EE++;
        if(muon_time_->at(im)!=0)        num_muon_time_noPU_prompt_EE++;
        if(muon_time_err_->at(im)!=-1)   num_muon_time_err_noPU_prompt_EE++;
        if(muon_time_err_i_->at(im)!=-1) num_muon_time_err_i_noPU_prompt_EE++;
        if(im==0) {
          num_vtx_noPU_prompt_EE++;
          if(vtx_time_->at(0)!=0)     num_vtx_time_noPU_prompt_EE++;
          if(vtx_time_err_->at(0)!=0) num_vtx_time_err_noPU_prompt_EE++;
        }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {
		  int track_type = -999;
		  if(simvtx_reco_==-999) track_type = 999;
		  else track_type = track_type_->at(im).at(it);

		  if(flag_track_pv_dz) {
		    //if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) > 0.2) continue;
		    if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EE) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // calculate fraction of observables
          num_track_noPU_prompt_EE++;
          if(track_mva_->at(im).at(it)>0) num_track_mva_noPU_prompt_EE++;
          if(track_mva_->at(im).at(it)>0.5) num_track_mva_cut_noPU_prompt_EE++;
          if(track_time_->at(im).at(it)!=0) num_track_time_noPU_prompt_EE++;
          if(track_time_err_->at(im).at(it)!=-1) num_track_time_err_noPU_prompt_EE++;
          if(track_time_err_i_->at(im).at(it)!=-1) num_track_time_err_i_noPU_prompt_EE++;
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_muon_track_time_err_noPU_prompt_EE++;
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_vtx_track_time_err_noPU_prompt_EE++;

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dt_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it));
          }

		  if(track_type==0) {
			// (muon, track)
			if(dtsig_cut!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut>2.0) {
			      h_noPU_prompt_PV_EE->Fill(dtsig_cut);
			      h_noPU_prompt_bx_PV_EE->Fill(track_bx_->at(im).at(it));
			      h_noPU_prompt_evtId_PV_EE->Fill(track_evtId_->at(im).at(it));
			      h_noPU_prompt_PVweight_PV_EE->Fill(track_PVweight_->at(im).at(it));
			      h_noPU_prompt_simvtx_bx_PV_EE->Fill(simvtx_bx_);
			      h_noPU_prompt_simvtx_evtId_PV_EE->Fill(simvtx_evtId_);
				}    
			  }
			  else {
			    h_noPU_prompt_PV_EE->Fill(dtsig_cut);
			    h_noPU_prompt_bx_PV_EE->Fill(track_bx_->at(im).at(it));
			    h_noPU_prompt_evtId_PV_EE->Fill(track_evtId_->at(im).at(it));
			    h_noPU_prompt_PVweight_PV_EE->Fill(track_PVweight_->at(im).at(it));
			    h_noPU_prompt_simvtx_bx_PV_EE->Fill(simvtx_bx_);
			    h_noPU_prompt_simvtx_evtId_PV_EE->Fill(simvtx_evtId_);
			  }
			}
			else h_noPU_prompt_no_tErr_EE->Fill(dtsig_cut);

			// (PV, track)
			if(dtsig_cut_vtx!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut_vtx>2.0) {
			      h_noPU_prompt_PV_EE_vtx->Fill(dtsig_cut_vtx);
			      h_noPU_prompt_bx_PV_EE_vtx->Fill(track_bx_->at(im).at(it));
			      h_noPU_prompt_evtId_PV_EE_vtx->Fill(track_evtId_->at(im).at(it));
			      h_noPU_prompt_PVweight_PV_EE_vtx->Fill(track_PVweight_->at(im).at(it));
			      h_noPU_prompt_simvtx_bx_PV_EE_vtx->Fill(simvtx_bx_);
			      h_noPU_prompt_simvtx_evtId_PV_EE_vtx->Fill(simvtx_evtId_);
				}    
			  }
			  else {
			    h_noPU_prompt_PV_EE_vtx->Fill(dtsig_cut_vtx);
			    h_noPU_prompt_bx_PV_EE_vtx->Fill(track_bx_->at(im).at(it));
			    h_noPU_prompt_evtId_PV_EE_vtx->Fill(track_evtId_->at(im).at(it));
			    h_noPU_prompt_PVweight_PV_EE_vtx->Fill(track_PVweight_->at(im).at(it));
			    h_noPU_prompt_simvtx_bx_PV_EE_vtx->Fill(simvtx_bx_);
			    h_noPU_prompt_simvtx_evtId_PV_EE_vtx->Fill(simvtx_evtId_);
			  }
			}
			else h_noPU_prompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);

			// 2D
			if(dtsig_cut!=0 && dtsig_cut_vtx!=0) h_noPU_prompt_PV_EE_2D->Fill(dtsig_cut, dtsig_cut_vtx);
		  }
		  else if(track_type==3) {
			if(dtsig_cut!=0) h_noPU_prompt_SV_EE->Fill(dtsig_cut);
			else h_noPU_prompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_prompt_SV_EE_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_prompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==1) {
			if(dtsig_cut!=0) h_noPU_prompt_PU_EE->Fill(dtsig_cut);
			else h_noPU_prompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_prompt_PU_EE_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_prompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==2) {
			if(dtsig_cut!=0) h_noPU_prompt_fake_EE->Fill(dtsig_cut);
			else h_noPU_prompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_prompt_fake_EE_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_prompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }
		  //else cout << "Something wrong with noPU prompt EE" << endl;
		  else {
			if(dtsig_cut!=0) h_noPU_prompt_not_vtx_matched_EE->Fill(dtsig_cut);
			else h_noPU_prompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_prompt_not_vtx_matched_EE_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_prompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }


		  // tsim test
		    // track
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_prompt_PV_track_time_EE->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==1 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_prompt_PU_track_time_EE->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==3 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_prompt_SV_track_time_EE->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		    // muon
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_prompt_PV_muon_time_EE->Fill(muon_time_->at(im) - muon_time_sim_->at(im));
		  }
		    // vertex
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_prompt_PV_vtx_time_EE->Fill(vtx_time_->at(im) - simvtx_time_*1e9);
		  }

		  dtsig_cut=0, dtsig_cut_vtx=0;
		  dt_cut=0, dt_cut_vtx=0;
		}
	  }
	} // End of muon loop
  } // End of event loop


  /////////////////////////
  //// noPU nonprompt ////
  /////////////////////////
  
  int n_muon_noPU_nonprompt_EB=0, n_muon_noPU_nonprompt_EE=0;
  int n_status_failed_noPU_nonprompt_EB=0, n_status_failed_noPU_nonprompt_EE=0;

  int num_muon_noPU_nonprompt_EB=0,  num_muon_mva_noPU_nonprompt_EB=0,  num_muon_mva_cut_noPU_nonprompt_EB=0,  num_muon_time_noPU_nonprompt_EB=0,  num_muon_time_err_noPU_nonprompt_EB=0,  num_muon_time_err_i_noPU_nonprompt_EB=0;
  int num_muon_noPU_nonprompt_EE=0,  num_muon_mva_noPU_nonprompt_EE=0,  num_muon_mva_cut_noPU_nonprompt_EE=0,  num_muon_time_noPU_nonprompt_EE=0,  num_muon_time_err_noPU_nonprompt_EE=0,  num_muon_time_err_i_noPU_nonprompt_EE=0;
  int num_track_noPU_nonprompt_EB=0, num_track_mva_noPU_nonprompt_EB=0, num_track_mva_cut_noPU_nonprompt_EB=0, num_track_time_noPU_nonprompt_EB=0, num_track_time_err_noPU_nonprompt_EB=0, num_track_time_err_i_noPU_nonprompt_EB=0;
  int num_track_noPU_nonprompt_EE=0, num_track_mva_noPU_nonprompt_EE=0, num_track_mva_cut_noPU_nonprompt_EE=0, num_track_time_noPU_nonprompt_EE=0, num_track_time_err_noPU_nonprompt_EE=0, num_track_time_err_i_noPU_nonprompt_EE=0;
  int num_vtx_noPU_nonprompt_EB=0,   num_vtx_time_noPU_nonprompt_EB=0,  num_vtx_time_err_noPU_nonprompt_EB=0;
  int num_vtx_noPU_nonprompt_EE=0,   num_vtx_time_noPU_nonprompt_EE=0,  num_vtx_time_err_noPU_nonprompt_EE=0;
  int num_both_muon_track_time_err_noPU_nonprompt_EB=0, num_both_muon_track_time_err_noPU_nonprompt_EE=0;
  int num_both_vtx_track_time_err_noPU_nonprompt_EB=0, num_both_vtx_track_time_err_noPU_nonprompt_EE=0;

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

		// calculate fraction of observables
        num_muon_noPU_nonprompt_EB++;
        if(muon_mva_->at(im)>0)          num_muon_mva_noPU_nonprompt_EB++;
        if(muon_mva_->at(im)>0.5)        num_muon_mva_cut_noPU_nonprompt_EB++;
        if(muon_time_->at(im)!=0)        num_muon_time_noPU_nonprompt_EB++;
        if(muon_time_err_->at(im)!=-1)   num_muon_time_err_noPU_nonprompt_EB++;
        if(muon_time_err_i_->at(im)!=-1) num_muon_time_err_i_noPU_nonprompt_EB++;
        if(im==0) {
          num_vtx_noPU_nonprompt_EB++;
          if(vtx_time_->at(0)!=0)     num_vtx_time_noPU_nonprompt_EB++;
          if(vtx_time_err_->at(0)!=0) num_vtx_time_err_noPU_nonprompt_EB++;
        }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {
		  int track_type = -999;
		  if(simvtx_reco_==-999) track_type = 999;
		  else track_type = track_type_->at(im).at(it);

		  if(flag_track_pv_dz) {
		    //if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) > 0.1) continue;
		    if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EB) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // calculate fraction of observables
          num_track_noPU_nonprompt_EB++;
          if(track_mva_->at(im).at(it)>0) num_track_mva_noPU_nonprompt_EB++;
          if(track_mva_->at(im).at(it)>0.5) num_track_mva_cut_noPU_nonprompt_EB++;
          if(track_time_->at(im).at(it)!=0) num_track_time_noPU_nonprompt_EB++;
          if(track_time_err_->at(im).at(it)!=-1) num_track_time_err_noPU_nonprompt_EB++;
          if(track_time_err_i_->at(im).at(it)!=-1) num_track_time_err_i_noPU_nonprompt_EB++;
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_muon_track_time_err_noPU_nonprompt_EB++;
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_vtx_track_time_err_noPU_nonprompt_EB++;

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dt_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it));
          }

		  if(track_type==0) {
			// (muon, track)
			if(dtsig_cut!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut>2.0) {
			      h_noPU_nonprompt_PV_EB->Fill(dtsig_cut);
			      h_noPU_nonprompt_bx_PV_EB->Fill(track_bx_->at(im).at(it));
			      h_noPU_nonprompt_evtId_PV_EB->Fill(track_evtId_->at(im).at(it));
			      h_noPU_nonprompt_PVweight_PV_EB->Fill(track_PVweight_->at(im).at(it));
			      h_noPU_nonprompt_simvtx_bx_PV_EB->Fill(simvtx_bx_);
			      h_noPU_nonprompt_simvtx_evtId_PV_EB->Fill(simvtx_evtId_);
				}    
			  }    
			  else {
			    h_noPU_nonprompt_PV_EB->Fill(dtsig_cut);
			    h_noPU_nonprompt_bx_PV_EB->Fill(track_bx_->at(im).at(it));
			    h_noPU_nonprompt_evtId_PV_EB->Fill(track_evtId_->at(im).at(it));
			    h_noPU_nonprompt_PVweight_PV_EB->Fill(track_PVweight_->at(im).at(it));
			    h_noPU_nonprompt_simvtx_bx_PV_EB->Fill(simvtx_bx_);
			    h_noPU_nonprompt_simvtx_evtId_PV_EB->Fill(simvtx_evtId_);
			  }
			}
			else h_noPU_nonprompt_no_tErr_EB->Fill(dtsig_cut);

			// (PV, track)
			if(dtsig_cut_vtx!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut_vtx>2.0) {
			      h_noPU_nonprompt_PV_EB_vtx->Fill(dtsig_cut_vtx);
			      h_noPU_nonprompt_bx_PV_EB_vtx->Fill(track_bx_->at(im).at(it));
			      h_noPU_nonprompt_evtId_PV_EB_vtx->Fill(track_evtId_->at(im).at(it));
			      h_noPU_nonprompt_PVweight_PV_EB_vtx->Fill(track_PVweight_->at(im).at(it));
			      h_noPU_nonprompt_simvtx_bx_PV_EB_vtx->Fill(simvtx_bx_);
			      h_noPU_nonprompt_simvtx_evtId_PV_EB_vtx->Fill(simvtx_evtId_);
				}    
			  }    
			  else {
			    h_noPU_nonprompt_PV_EB_vtx->Fill(dtsig_cut_vtx);
			    h_noPU_nonprompt_bx_PV_EB_vtx->Fill(track_bx_->at(im).at(it));
			    h_noPU_nonprompt_evtId_PV_EB_vtx->Fill(track_evtId_->at(im).at(it));
			    h_noPU_nonprompt_PVweight_PV_EB_vtx->Fill(track_PVweight_->at(im).at(it));
			    h_noPU_nonprompt_simvtx_bx_PV_EB_vtx->Fill(simvtx_bx_);
			    h_noPU_nonprompt_simvtx_evtId_PV_EB_vtx->Fill(simvtx_evtId_);
			  }
			}
			else h_noPU_nonprompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);


			if(dtsig_cut!=0 && dtsig_cut_vtx!=0) h_noPU_nonprompt_PV_EB_2D->Fill(dtsig_cut, dtsig_cut_vtx);
		  }
		  else if(track_type==3) {
			if(dtsig_cut!=0) h_noPU_nonprompt_SV_EB->Fill(dtsig_cut);
			else h_noPU_nonprompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_nonprompt_SV_EB_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_nonprompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==1) {
			if(dtsig_cut!=0) h_noPU_nonprompt_PU_EB->Fill(dtsig_cut);
			else h_noPU_nonprompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_nonprompt_PU_EB_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_nonprompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==2) {
			if(dtsig_cut!=0) h_noPU_nonprompt_fake_EB->Fill(dtsig_cut);
			else h_noPU_nonprompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_nonprompt_fake_EB_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_nonprompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }
		  //else cout << "Something wrong with noPU nonprompt EB" << endl;
		  else {
			if(dtsig_cut!=0) h_noPU_nonprompt_not_vtx_matched_EB->Fill(dtsig_cut);
			else h_noPU_nonprompt_no_tErr_EB->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_nonprompt_not_vtx_matched_EB_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_nonprompt_no_tErr_EB_vtx->Fill(dtsig_cut_vtx);
		  }


		  // tsim test
		    // track
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_nonprompt_PV_track_time_EB->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==1 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_nonprompt_PU_track_time_EB->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==3 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_nonprompt_SV_track_time_EB->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		    // muon
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_nonprompt_PV_muon_time_EB->Fill(muon_time_->at(im) - muon_time_sim_->at(im));
		  }
		    // vertex
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_nonprompt_PV_vtx_time_EB->Fill(vtx_time_->at(im) - simvtx_time_*1e9);
		  }

		  dtsig_cut=0, dtsig_cut_vtx=0;
		  dt_cut=0, dt_cut_vtx=0;
		}
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

		// calculate fraction of observables
        num_muon_noPU_nonprompt_EE++;
        if(muon_mva_->at(im)>0)          num_muon_mva_noPU_nonprompt_EE++;
        if(muon_mva_->at(im)>0.5)        num_muon_mva_cut_noPU_nonprompt_EE++;
        if(muon_time_->at(im)!=0)        num_muon_time_noPU_nonprompt_EE++;
        if(muon_time_err_->at(im)!=-1)   num_muon_time_err_noPU_nonprompt_EE++;
        if(muon_time_err_i_->at(im)!=-1) num_muon_time_err_i_noPU_nonprompt_EE++;
        if(im==0) {
          num_vtx_noPU_nonprompt_EE++;
          if(vtx_time_->at(0)!=0)     num_vtx_time_noPU_nonprompt_EE++;
          if(vtx_time_err_->at(0)!=0) num_vtx_time_err_noPU_nonprompt_EE++;
        }

		// track loop
		for(int it=0; it<track_pt_->at(im).size(); it++) {
		  int track_type = -999;
		  if(simvtx_reco_==-999) track_type = 999;
		  else track_type = track_type_->at(im).at(it);

		  if(flag_track_pv_dz) {
		    //if(TMath::Abs((muon_vz_->at(im)-track_vz_->at(im).at(it))) > 0.2) continue;
		    if(TMath::Abs(track_pv_dz_->at(im).at(it)) > track_pv_dz_cut_EE) continue;
		  }
		  if(flag_track_PVweight) {
			if(track_PVweight_->at(im).at(it)<track_PVweight_cut) continue;
		  }

		  // calculate fraction of observables
          num_track_noPU_nonprompt_EE++;
          if(track_mva_->at(im).at(it)>0) num_track_mva_noPU_nonprompt_EE++;
          if(track_mva_->at(im).at(it)>0.5) num_track_mva_cut_noPU_nonprompt_EE++;
          if(track_time_->at(im).at(it)!=0) num_track_time_noPU_nonprompt_EE++;
          if(track_time_err_->at(im).at(it)!=-1) num_track_time_err_noPU_nonprompt_EE++;
          if(track_time_err_i_->at(im).at(it)!=-1) num_track_time_err_i_noPU_nonprompt_EE++;
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_muon_track_time_err_noPU_nonprompt_EE++;
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) num_both_vtx_track_time_err_noPU_nonprompt_EE++;

		  // define dtsignif
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(muon_time_err_->at(im)*muon_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
		    dtsig_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it))/TMath::Sqrt(vtx_time_err_->at(im)*vtx_time_err_->at(im) + track_time_err_->at(im).at(it)*track_time_err_->at(im).at(it));
		  }
		  // define dt
		  if(muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			dt_cut = TMath::Abs(muon_time_->at(im)-track_time_->at(im).at(it));
		  }
		  if(vtx_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
            dt_cut_vtx = TMath::Abs(vtx_time_->at(im)-track_time_->at(im).at(it));
          }

		  if(track_type==0) {
			// (muon, track)
			if(dtsig_cut!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut>2.0) {
			      h_noPU_nonprompt_PV_EE->Fill(dtsig_cut);
			      h_noPU_nonprompt_bx_PV_EE->Fill(track_bx_->at(im).at(it));
			      h_noPU_nonprompt_evtId_PV_EE->Fill(track_evtId_->at(im).at(it));
			      h_noPU_nonprompt_PVweight_PV_EE->Fill(track_PVweight_->at(im).at(it));
			      h_noPU_nonprompt_simvtx_bx_PV_EE->Fill(simvtx_bx_);
			      h_noPU_nonprompt_simvtx_evtId_PV_EE->Fill(simvtx_evtId_);
				}
			  }
			  else {
			    h_noPU_nonprompt_PV_EE->Fill(dtsig_cut);
			    h_noPU_nonprompt_bx_PV_EE->Fill(track_bx_->at(im).at(it));
			    h_noPU_nonprompt_evtId_PV_EE->Fill(track_evtId_->at(im).at(it));
			    h_noPU_nonprompt_PVweight_PV_EE->Fill(track_PVweight_->at(im).at(it));
			    h_noPU_nonprompt_simvtx_bx_PV_EE->Fill(simvtx_bx_);
			    h_noPU_nonprompt_simvtx_evtId_PV_EE->Fill(simvtx_evtId_);
			  }
			}
			else h_noPU_nonprompt_no_tErr_EE->Fill(dtsig_cut);

			// (PV, track)
			if(dtsig_cut_vtx!=0) {
			  if(flag_2sigma_PV) {
				if(dtsig_cut_vtx>2.0) {
			      h_noPU_nonprompt_PV_EE_vtx->Fill(dtsig_cut_vtx);
			      h_noPU_nonprompt_bx_PV_EE_vtx->Fill(track_bx_->at(im).at(it));
			      h_noPU_nonprompt_evtId_PV_EE_vtx->Fill(track_evtId_->at(im).at(it));
			      h_noPU_nonprompt_PVweight_PV_EE_vtx->Fill(track_PVweight_->at(im).at(it));
			      h_noPU_nonprompt_simvtx_bx_PV_EE_vtx->Fill(simvtx_bx_);
			      h_noPU_nonprompt_simvtx_evtId_PV_EE_vtx->Fill(simvtx_evtId_);
				}
			  }
			  else {
			    h_noPU_nonprompt_PV_EE_vtx->Fill(dtsig_cut_vtx);
			    h_noPU_nonprompt_bx_PV_EE_vtx->Fill(track_bx_->at(im).at(it));
			    h_noPU_nonprompt_evtId_PV_EE_vtx->Fill(track_evtId_->at(im).at(it));
			    h_noPU_nonprompt_PVweight_PV_EE_vtx->Fill(track_PVweight_->at(im).at(it));
			    h_noPU_nonprompt_simvtx_bx_PV_EE_vtx->Fill(simvtx_bx_);
			    h_noPU_nonprompt_simvtx_evtId_PV_EE_vtx->Fill(simvtx_evtId_);
			  }
			}
			else h_noPU_nonprompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);

			// 2D
			if(dtsig_cut!=0 && dtsig_cut_vtx!=0) h_noPU_nonprompt_PV_EE_2D->Fill(dtsig_cut, dtsig_cut_vtx);
		  }
		  else if(track_type==3) {
			if(dtsig_cut!=0) h_noPU_nonprompt_SV_EE->Fill(dtsig_cut);
			else h_noPU_nonprompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_nonprompt_SV_EE_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_nonprompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==1) {
			if(dtsig_cut!=0) h_noPU_nonprompt_PU_EE->Fill(dtsig_cut);
			else h_noPU_nonprompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_nonprompt_PU_EE_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_nonprompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }
		  else if(track_type==2) {
			if(dtsig_cut!=0) h_noPU_nonprompt_fake_EE->Fill(dtsig_cut);
			else h_noPU_nonprompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_nonprompt_fake_EE_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_nonprompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }
		  //else cout << "Something wrong with noPU nonprompt EE" << endl;
		  else {
			if(dtsig_cut!=0) h_noPU_nonprompt_not_vtx_matched_EE->Fill(dtsig_cut);
			else h_noPU_nonprompt_no_tErr_EE->Fill(dtsig_cut);
			if(dtsig_cut_vtx!=0) h_noPU_nonprompt_not_vtx_matched_EE_vtx->Fill(dtsig_cut_vtx);
			else h_noPU_nonprompt_no_tErr_EE_vtx->Fill(dtsig_cut_vtx);
		  }


		  // tsim test
		    // track
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_nonprompt_PV_track_time_EE->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==1 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_nonprompt_PU_track_time_EE->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		  if(track_type_->at(im).at(it)==3 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_nonprompt_SV_track_time_EE->Fill(track_time_->at(im).at(it) - track_time_sim_->at(im).at(it));
		  }
		    // muon
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_nonprompt_PV_muon_time_EE->Fill(muon_time_->at(im) - muon_time_sim_->at(im));
		  }
		    // vertex
		  if(track_type_->at(im).at(it)==0 && dtsig_cut!=0 && track_time_sim_->at(im).at(it)!=-1 && muon_time_err_->at(im)>0 && track_time_err_->at(im).at(it)>0) {
			h_noPU_nonprompt_PV_vtx_time_EE->Fill(vtx_time_->at(im) - simvtx_time_*1e9);
		  }

		  dtsig_cut=0, dtsig_cut_vtx=0;
		  dt_cut=0, dt_cut_vtx=0;
		}

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


  // check fraction of observables
  cout << endl;
  cout << endl;
  cout << "////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "//// Check on whether muon/track/vtx have MVA/time/tErr or pass MVA cut ////" << endl;
  cout << "////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "[PU200 prompt]" << endl;
  cout << "               |           time          |          tErr           |          tErr_i         |           MVA           |    MVA > 0.5 / total    |" << endl;
  cout << "      muon     | "
       << num_muon_time_PU200_prompt_EB+num_muon_time_PU200_prompt_EE << "/" << num_muon_PU200_prompt_EB+num_muon_PU200_prompt_EE << " = " << 1.*(num_muon_time_PU200_prompt_EB+num_muon_time_PU200_prompt_EE)/(num_muon_PU200_prompt_EB+num_muon_PU200_prompt_EE) << " % | "
       << num_muon_time_err_PU200_prompt_EB+num_muon_time_err_PU200_prompt_EE << "/" << num_muon_PU200_prompt_EB+num_muon_PU200_prompt_EE << " = " << 1.*(num_muon_time_err_PU200_prompt_EB+num_muon_time_err_PU200_prompt_EE)/(num_muon_PU200_prompt_EB+num_muon_PU200_prompt_EE) << " % | "
       << num_muon_time_err_i_PU200_prompt_EB+num_muon_time_err_i_PU200_prompt_EE << "/" << num_muon_PU200_prompt_EB+num_muon_PU200_prompt_EE << " = " << 1.*(num_muon_time_err_i_PU200_prompt_EB+num_muon_time_err_i_PU200_prompt_EE)/(num_muon_PU200_prompt_EB+num_muon_PU200_prompt_EE) << " % | "
       << num_muon_mva_PU200_prompt_EB+num_muon_mva_PU200_prompt_EE << "/" << num_muon_PU200_prompt_EB+num_muon_PU200_prompt_EE << " = " << 1.*(num_muon_mva_PU200_prompt_EB+num_muon_mva_PU200_prompt_EE)/(num_muon_PU200_prompt_EB+num_muon_PU200_prompt_EE) << " % | "
       << num_muon_mva_cut_PU200_prompt_EB+num_muon_mva_cut_PU200_prompt_EE << "/" << num_muon_PU200_prompt_EB+num_muon_PU200_prompt_EE << " = " << 1.*(num_muon_mva_cut_PU200_prompt_EB+num_muon_mva_cut_PU200_prompt_EE)/(num_muon_PU200_prompt_EB+num_muon_PU200_prompt_EE) << " % | "
       << endl;
  cout << "     vertex    | "
       << num_vtx_time_PU200_prompt_EB+num_vtx_time_PU200_prompt_EE << "/" << num_vtx_PU200_prompt_EB+num_vtx_PU200_prompt_EE << " = " << 1.*(num_vtx_time_PU200_prompt_EB+num_vtx_time_PU200_prompt_EE)/(num_vtx_PU200_prompt_EB+num_vtx_PU200_prompt_EE) << " % | "
       << num_vtx_time_err_PU200_prompt_EB+num_vtx_time_err_PU200_prompt_EE << "/" << num_vtx_PU200_prompt_EB+num_vtx_PU200_prompt_EE << " = " << 1.*(num_vtx_time_err_PU200_prompt_EB+num_vtx_time_err_PU200_prompt_EE)/(num_vtx_PU200_prompt_EB+num_vtx_PU200_prompt_EE) << " % | "
	   << endl;
  cout << "      track    | "
       << num_track_time_PU200_prompt_EB+num_track_time_PU200_prompt_EE << "/" << num_track_PU200_prompt_EB+num_track_PU200_prompt_EE << " = " << 1.*(num_track_time_PU200_prompt_EB+num_track_time_PU200_prompt_EE)/(num_track_PU200_prompt_EB+num_track_PU200_prompt_EE) << " % | "
       << num_track_time_err_PU200_prompt_EB+num_track_time_err_PU200_prompt_EE << "/" << num_track_PU200_prompt_EB+num_track_PU200_prompt_EE << " = " << 1.*(num_track_time_err_PU200_prompt_EB+num_track_time_err_PU200_prompt_EE)/(num_track_PU200_prompt_EB+num_track_PU200_prompt_EE) << " % | "
       << num_track_time_err_i_PU200_prompt_EB+num_track_time_err_i_PU200_prompt_EE << "/" << num_track_PU200_prompt_EB+num_track_PU200_prompt_EE << " = " << 1.*(num_track_time_err_i_PU200_prompt_EB+num_track_time_err_i_PU200_prompt_EE)/(num_track_PU200_prompt_EB+num_track_PU200_prompt_EE) << " % | "
       << num_track_mva_PU200_prompt_EB+num_track_mva_PU200_prompt_EE << "/" << num_track_PU200_prompt_EB+num_track_PU200_prompt_EE << " = " << 1.*(num_track_mva_PU200_prompt_EB+num_track_mva_PU200_prompt_EE)/(num_track_PU200_prompt_EB+num_track_PU200_prompt_EE) << " % | "
       << num_track_mva_cut_PU200_prompt_EB+num_track_mva_cut_PU200_prompt_EE << "/" << num_track_PU200_prompt_EB+num_track_PU200_prompt_EE << " = " << 1.*(num_track_mva_cut_PU200_prompt_EB+num_track_mva_cut_PU200_prompt_EE)/(num_track_PU200_prompt_EB+num_track_PU200_prompt_EE) << " % | "
       << endl;
  cout << "  muon&&track  | "
       << num_both_muon_track_time_err_PU200_prompt_EB+num_both_muon_track_time_err_PU200_prompt_EE << "/" << num_track_PU200_prompt_EB+num_track_PU200_prompt_EE << " = " << 1.*(num_both_muon_track_time_err_PU200_prompt_EB+num_both_muon_track_time_err_PU200_prompt_EE)/(num_track_PU200_prompt_EB+num_track_PU200_prompt_EE) << " % | "
	   << endl;
  cout << "  vtx&&track   | "
       << num_both_vtx_track_time_err_PU200_prompt_EB+num_both_vtx_track_time_err_PU200_prompt_EE << "/" << num_track_PU200_prompt_EB+num_track_PU200_prompt_EE << " = " << 1.*(num_both_vtx_track_time_err_PU200_prompt_EB+num_both_vtx_track_time_err_PU200_prompt_EE)/(num_track_PU200_prompt_EB+num_track_PU200_prompt_EE) << " % | "
       << endl;
       cout << endl;
  cout << "[PU200 nonprompt]" << endl;
  cout << "               |           time          |          tErr           |          tErr_i         |           MVA           |    MVA > 0.5 / total    |" << endl;
  cout << "      muon     | "
       << num_muon_time_PU200_nonprompt_EB+num_muon_time_PU200_nonprompt_EE << "/" << num_muon_PU200_nonprompt_EB+num_muon_PU200_nonprompt_EE << " = " << 1.*(num_muon_time_PU200_nonprompt_EB+num_muon_time_PU200_nonprompt_EE)/(num_muon_PU200_nonprompt_EB+num_muon_PU200_nonprompt_EE) << " % | "
       << num_muon_time_err_PU200_nonprompt_EB+num_muon_time_err_PU200_nonprompt_EE << "/" << num_muon_PU200_nonprompt_EB+num_muon_PU200_nonprompt_EE << " = " << 1.*(num_muon_time_err_PU200_nonprompt_EB+num_muon_time_err_PU200_nonprompt_EE)/(num_muon_PU200_nonprompt_EB+num_muon_PU200_nonprompt_EE) << " % | "
       << num_muon_time_err_i_PU200_nonprompt_EB+num_muon_time_err_i_PU200_nonprompt_EE << "/" << num_muon_PU200_nonprompt_EB+num_muon_PU200_nonprompt_EE << " = " << 1.*(num_muon_time_err_i_PU200_nonprompt_EB+num_muon_time_err_i_PU200_nonprompt_EE)/(num_muon_PU200_nonprompt_EB+num_muon_PU200_nonprompt_EE) << " % | "
       << num_muon_mva_PU200_nonprompt_EB+num_muon_mva_PU200_nonprompt_EE << "/" << num_muon_PU200_nonprompt_EB+num_muon_PU200_nonprompt_EE << " = " << 1.*(num_muon_mva_PU200_nonprompt_EB+num_muon_mva_PU200_nonprompt_EE)/(num_muon_PU200_nonprompt_EB+num_muon_PU200_nonprompt_EE) << " % | "
       << num_muon_mva_cut_PU200_nonprompt_EB+num_muon_mva_cut_PU200_nonprompt_EE << "/" << num_muon_PU200_nonprompt_EB+num_muon_PU200_nonprompt_EE << " = " << 1.*(num_muon_mva_cut_PU200_nonprompt_EB+num_muon_mva_cut_PU200_nonprompt_EE)/(num_muon_PU200_nonprompt_EB+num_muon_PU200_nonprompt_EE) << " % | "
       << endl;
  cout << "     vertex    | "
       << num_vtx_time_PU200_nonprompt_EB+num_vtx_time_PU200_nonprompt_EE << "/" << num_vtx_PU200_nonprompt_EB+num_vtx_PU200_nonprompt_EE << " = " << 1.*(num_vtx_time_PU200_nonprompt_EB+num_vtx_time_PU200_nonprompt_EE)/(num_vtx_PU200_nonprompt_EB+num_vtx_PU200_nonprompt_EE) << " % | "
       << num_vtx_time_err_PU200_nonprompt_EB+num_vtx_time_err_PU200_nonprompt_EE << "/" << num_vtx_PU200_nonprompt_EB+num_vtx_PU200_nonprompt_EE << " = " << 1.*(num_vtx_time_err_PU200_nonprompt_EB+num_vtx_time_err_PU200_nonprompt_EE)/(num_vtx_PU200_nonprompt_EB+num_vtx_PU200_nonprompt_EE) << " % | "
	   << endl;
  cout << "      track    | "
       << num_track_time_PU200_nonprompt_EB+num_track_time_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE << " = " << 1.*(num_track_time_PU200_nonprompt_EB+num_track_time_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE) << " % | "
       << num_track_time_err_PU200_nonprompt_EB+num_track_time_err_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE << " = " << 1.*(num_track_time_err_PU200_nonprompt_EB+num_track_time_err_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE) << " % | "
       << num_track_time_err_i_PU200_nonprompt_EB+num_track_time_err_i_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE << " = " << 1.*(num_track_time_err_i_PU200_nonprompt_EB+num_track_time_err_i_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE) << " % | "
       << num_track_mva_PU200_nonprompt_EB+num_track_mva_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE << " = " << 1.*(num_track_mva_PU200_nonprompt_EB+num_track_mva_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE) << " % | "
       << num_track_mva_cut_PU200_nonprompt_EB+num_track_mva_cut_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE << " = " << 1.*(num_track_mva_cut_PU200_nonprompt_EB+num_track_mva_cut_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE) << " % | "
       << endl;
  cout << "  muon&&track  | "
       << num_both_muon_track_time_err_PU200_nonprompt_EB+num_both_muon_track_time_err_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE << " = " << 1.*(num_both_muon_track_time_err_PU200_nonprompt_EB+num_both_muon_track_time_err_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE) << " % | "
	   << endl;
  cout << "  vtx&&track   | "
       << num_both_vtx_track_time_err_PU200_nonprompt_EB+num_both_vtx_track_time_err_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE << " = " << 1.*(num_both_vtx_track_time_err_PU200_nonprompt_EB+num_both_vtx_track_time_err_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EB+num_track_PU200_nonprompt_EE) << " % | "
       << endl;
       cout << endl;
  cout << "[noPU prompt]" << endl;
  cout << "               |           time          |          tErr           |          tErr_i         |           MVA           |    MVA > 0.5 / total    |" << endl;
  cout << "      muon     | "
       << num_muon_time_noPU_prompt_EB+num_muon_time_noPU_prompt_EE << "/" << num_muon_noPU_prompt_EB+num_muon_noPU_prompt_EE << " = " << 1.*(num_muon_time_noPU_prompt_EB+num_muon_time_noPU_prompt_EE)/(num_muon_noPU_prompt_EB+num_muon_noPU_prompt_EE) << " % | "
       << num_muon_time_err_noPU_prompt_EB+num_muon_time_err_noPU_prompt_EE << "/" << num_muon_noPU_prompt_EB+num_muon_noPU_prompt_EE << " = " << 1.*(num_muon_time_err_noPU_prompt_EB+num_muon_time_err_noPU_prompt_EE)/(num_muon_noPU_prompt_EB+num_muon_noPU_prompt_EE) << " % | "
       << num_muon_time_err_i_noPU_prompt_EB+num_muon_time_err_i_noPU_prompt_EE << "/" << num_muon_noPU_prompt_EB+num_muon_noPU_prompt_EE << " = " << 1.*(num_muon_time_err_i_noPU_prompt_EB+num_muon_time_err_i_noPU_prompt_EE)/(num_muon_noPU_prompt_EB+num_muon_noPU_prompt_EE) << " % | "
       << num_muon_mva_noPU_prompt_EB+num_muon_mva_noPU_prompt_EE << "/" << num_muon_noPU_prompt_EB+num_muon_noPU_prompt_EE << " = " << 1.*(num_muon_mva_noPU_prompt_EB+num_muon_mva_noPU_prompt_EE)/(num_muon_noPU_prompt_EB+num_muon_noPU_prompt_EE) << " % | "
       << num_muon_mva_cut_noPU_prompt_EB+num_muon_mva_cut_noPU_prompt_EE << "/" << num_muon_noPU_prompt_EB+num_muon_noPU_prompt_EE << " = " << 1.*(num_muon_mva_cut_noPU_prompt_EB+num_muon_mva_cut_noPU_prompt_EE)/(num_muon_noPU_prompt_EB+num_muon_noPU_prompt_EE) << " % | "
       << endl;
  cout << "     vertex    | "
       << num_vtx_time_noPU_prompt_EB+num_vtx_time_noPU_prompt_EE << "/" << num_vtx_noPU_prompt_EB+num_vtx_noPU_prompt_EE << " = " << 1.*(num_vtx_time_noPU_prompt_EB+num_vtx_time_noPU_prompt_EE)/(num_vtx_noPU_prompt_EB+num_vtx_noPU_prompt_EE) << " % | "
       << num_vtx_time_err_noPU_prompt_EB+num_vtx_time_err_noPU_prompt_EE << "/" << num_vtx_noPU_prompt_EB+num_vtx_noPU_prompt_EE << " = " << 1.*(num_vtx_time_err_noPU_prompt_EB+num_vtx_time_err_noPU_prompt_EE)/(num_vtx_noPU_prompt_EB+num_vtx_noPU_prompt_EE) << " % | "
	   << endl;
  cout << "      track    | "
       << num_track_time_noPU_prompt_EB+num_track_time_noPU_prompt_EE << "/" << num_track_noPU_prompt_EB+num_track_noPU_prompt_EE << " = " << 1.*(num_track_time_noPU_prompt_EB+num_track_time_noPU_prompt_EE)/(num_track_noPU_prompt_EB+num_track_noPU_prompt_EE) << " % | "
       << num_track_time_err_noPU_prompt_EB+num_track_time_err_noPU_prompt_EE << "/" << num_track_noPU_prompt_EB+num_track_noPU_prompt_EE << " = " << 1.*(num_track_time_err_noPU_prompt_EB+num_track_time_err_noPU_prompt_EE)/(num_track_noPU_prompt_EB+num_track_noPU_prompt_EE) << " % | "
       << num_track_time_err_i_noPU_prompt_EB+num_track_time_err_i_noPU_prompt_EE << "/" << num_track_noPU_prompt_EB+num_track_noPU_prompt_EE << " = " << 1.*(num_track_time_err_i_noPU_prompt_EB+num_track_time_err_i_noPU_prompt_EE)/(num_track_noPU_prompt_EB+num_track_noPU_prompt_EE) << " % | "
       << num_track_mva_noPU_prompt_EB+num_track_mva_noPU_prompt_EE << "/" << num_track_noPU_prompt_EB+num_track_noPU_prompt_EE << " = " << 1.*(num_track_mva_noPU_prompt_EB+num_track_mva_noPU_prompt_EE)/(num_track_noPU_prompt_EB+num_track_noPU_prompt_EE) << " % | "
       << num_track_mva_cut_noPU_prompt_EB+num_track_mva_cut_noPU_prompt_EE << "/" << num_track_noPU_prompt_EB+num_track_noPU_prompt_EE << " = " << 1.*(num_track_mva_cut_noPU_prompt_EB+num_track_mva_cut_noPU_prompt_EE)/(num_track_noPU_prompt_EB+num_track_noPU_prompt_EE) << " % | "
       << endl;
  cout << "  muon&&track  | "
       << num_both_muon_track_time_err_noPU_prompt_EB+num_both_muon_track_time_err_noPU_prompt_EE << "/" << num_track_noPU_prompt_EB+num_track_noPU_prompt_EE << " = " << 1.*(num_both_muon_track_time_err_noPU_prompt_EB+num_both_muon_track_time_err_noPU_prompt_EE)/(num_track_noPU_prompt_EB+num_track_noPU_prompt_EE) << " % | "
	   << endl;
  cout << "  vtx&&track   | "
       << num_both_vtx_track_time_err_noPU_prompt_EB+num_both_vtx_track_time_err_noPU_prompt_EE << "/" << num_track_noPU_prompt_EB+num_track_noPU_prompt_EE << " = " << 1.*(num_both_vtx_track_time_err_noPU_prompt_EB+num_both_vtx_track_time_err_noPU_prompt_EE)/(num_track_noPU_prompt_EB+num_track_noPU_prompt_EE) << " % | "
       << endl;
       cout << endl;
  cout << "[noPU nonprompt]" << endl;
  cout << "               |           time          |          tErr           |          tErr_i         |           MVA           |    MVA > 0.5 / total    |" << endl;
  cout << "      muon     | "
       << num_muon_time_noPU_nonprompt_EB+num_muon_time_noPU_nonprompt_EE << "/" << num_muon_noPU_nonprompt_EB+num_muon_noPU_nonprompt_EE << " = " << 1.*(num_muon_time_noPU_nonprompt_EB+num_muon_time_noPU_nonprompt_EE)/(num_muon_noPU_nonprompt_EB+num_muon_noPU_nonprompt_EE) << " % | "
       << num_muon_time_err_noPU_nonprompt_EB+num_muon_time_err_noPU_nonprompt_EE << "/" << num_muon_noPU_nonprompt_EB+num_muon_noPU_nonprompt_EE << " = " << 1.*(num_muon_time_err_noPU_nonprompt_EB+num_muon_time_err_noPU_nonprompt_EE)/(num_muon_noPU_nonprompt_EB+num_muon_noPU_nonprompt_EE) << " % | "
       << num_muon_time_err_i_noPU_nonprompt_EB+num_muon_time_err_i_noPU_nonprompt_EE << "/" << num_muon_noPU_nonprompt_EB+num_muon_noPU_nonprompt_EE << " = " << 1.*(num_muon_time_err_i_noPU_nonprompt_EB+num_muon_time_err_i_noPU_nonprompt_EE)/(num_muon_noPU_nonprompt_EB+num_muon_noPU_nonprompt_EE) << " % | "
       << num_muon_mva_noPU_nonprompt_EB+num_muon_mva_noPU_nonprompt_EE << "/" << num_muon_noPU_nonprompt_EB+num_muon_noPU_nonprompt_EE << " = " << 1.*(num_muon_mva_noPU_nonprompt_EB+num_muon_mva_noPU_nonprompt_EE)/(num_muon_noPU_nonprompt_EB+num_muon_noPU_nonprompt_EE) << " % | "
       << num_muon_mva_cut_noPU_nonprompt_EB+num_muon_mva_cut_noPU_nonprompt_EE << "/" << num_muon_noPU_nonprompt_EB+num_muon_noPU_nonprompt_EE << " = " << 1.*(num_muon_mva_cut_noPU_nonprompt_EB+num_muon_mva_cut_noPU_nonprompt_EE)/(num_muon_noPU_nonprompt_EB+num_muon_noPU_nonprompt_EE) << " % | "
       << endl;
  cout << "     vertex    | "
       << num_vtx_time_noPU_nonprompt_EB+num_vtx_time_noPU_nonprompt_EE << "/" << num_vtx_noPU_nonprompt_EB+num_vtx_noPU_nonprompt_EE << " = " << 1.*(num_vtx_time_noPU_nonprompt_EB+num_vtx_time_noPU_nonprompt_EE)/(num_vtx_noPU_nonprompt_EB+num_vtx_noPU_nonprompt_EE) << " % | "
       << num_vtx_time_err_noPU_nonprompt_EB+num_vtx_time_err_noPU_nonprompt_EE << "/" << num_vtx_noPU_nonprompt_EB+num_vtx_noPU_nonprompt_EE << " = " << 1.*(num_vtx_time_err_noPU_nonprompt_EB+num_vtx_time_err_noPU_nonprompt_EE)/(num_vtx_noPU_nonprompt_EB+num_vtx_noPU_nonprompt_EE) << " % | "
	   << endl;
  cout << "      track    | "
       << num_track_time_noPU_nonprompt_EB+num_track_time_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE << " = " << 1.*(num_track_time_noPU_nonprompt_EB+num_track_time_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE) << " % | "
       << num_track_time_err_noPU_nonprompt_EB+num_track_time_err_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE << " = " << 1.*(num_track_time_err_noPU_nonprompt_EB+num_track_time_err_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE) << " % | "
       << num_track_time_err_i_noPU_nonprompt_EB+num_track_time_err_i_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE << " = " << 1.*(num_track_time_err_i_noPU_nonprompt_EB+num_track_time_err_i_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE) << " % | "
       << num_track_mva_noPU_nonprompt_EB+num_track_mva_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE << " = " << 1.*(num_track_mva_noPU_nonprompt_EB+num_track_mva_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE) << " % | "
       << num_track_mva_cut_noPU_nonprompt_EB+num_track_mva_cut_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE << " = " << 1.*(num_track_mva_cut_noPU_nonprompt_EB+num_track_mva_cut_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE) << " % | "
       << endl;
  cout << "  muon&&track  | "
       << num_both_muon_track_time_err_noPU_nonprompt_EB+num_both_muon_track_time_err_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE << " = " << 1.*(num_both_muon_track_time_err_noPU_nonprompt_EB+num_both_muon_track_time_err_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE) << " % | "
	   << endl;
  cout << "  vtx&&track   | "
       << num_both_vtx_track_time_err_noPU_nonprompt_EB+num_both_vtx_track_time_err_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE << " = " << 1.*(num_both_vtx_track_time_err_noPU_nonprompt_EB+num_both_vtx_track_time_err_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EB+num_track_noPU_nonprompt_EE) << " % | "
       << endl;

  cout << endl;
  cout << endl;
  // Barrel
  cout << "////////" << endl;
  cout << "[Barrel]" << endl;
  cout << "////////" << endl;
  cout << "[PU200 prompt]" << endl;
  cout << "               |           time          |          tErr           |          tErr_i         |           MVA           |    MVA > 0.5 / total    |" << endl;
  cout << "      muon     | "
       << num_muon_time_PU200_prompt_EB << "/" << num_muon_PU200_prompt_EB << " = " << 1.*(num_muon_time_PU200_prompt_EB)/(num_muon_PU200_prompt_EB) << " % | "
       << num_muon_time_err_PU200_prompt_EB << "/" << num_muon_PU200_prompt_EB << " = " << 1.*(num_muon_time_err_PU200_prompt_EB)/(num_muon_PU200_prompt_EB) << " % | "
       << num_muon_time_err_i_PU200_prompt_EB << "/" << num_muon_PU200_prompt_EB << " = " << 1.*(num_muon_time_err_i_PU200_prompt_EB)/(num_muon_PU200_prompt_EB) << " % | "
       << num_muon_mva_PU200_prompt_EB << "/" << num_muon_PU200_prompt_EB << " = " << 1.*(num_muon_mva_PU200_prompt_EB)/(num_muon_PU200_prompt_EB) << " % | "
       << num_muon_mva_cut_PU200_prompt_EB << "/" << num_muon_PU200_prompt_EB << " = " << 1.*(num_muon_mva_cut_PU200_prompt_EB)/(num_muon_PU200_prompt_EB) << " % | "
       << endl;
  cout << "     vertex    | "
       << num_vtx_time_PU200_prompt_EB << "/" << num_vtx_PU200_prompt_EB << " = " << 1.*(num_vtx_time_PU200_prompt_EB)/(num_vtx_PU200_prompt_EB) << " % | "
       << num_vtx_time_err_PU200_prompt_EB << "/" << num_vtx_PU200_prompt_EB << " = " << 1.*(num_vtx_time_err_PU200_prompt_EB)/(num_vtx_PU200_prompt_EB) << " % | "
	   << endl;
  cout << "      track    | "
       << num_track_time_PU200_prompt_EB << "/" << num_track_PU200_prompt_EB << " = " << 1.*(num_track_time_PU200_prompt_EB)/(num_track_PU200_prompt_EB) << " % | "
       << num_track_time_err_PU200_prompt_EB << "/" << num_track_PU200_prompt_EB << " = " << 1.*(num_track_time_err_PU200_prompt_EB)/(num_track_PU200_prompt_EB) << " % | "
       << num_track_time_err_i_PU200_prompt_EB << "/" << num_track_PU200_prompt_EB << " = " << 1.*(num_track_time_err_i_PU200_prompt_EB)/(num_track_PU200_prompt_EB) << " % | "
       << num_track_mva_PU200_prompt_EB << "/" << num_track_PU200_prompt_EB << " = " << 1.*(num_track_mva_PU200_prompt_EB)/(num_track_PU200_prompt_EB) << " % | "
       << num_track_mva_cut_PU200_prompt_EB << "/" << num_track_PU200_prompt_EB << " = " << 1.*(num_track_mva_cut_PU200_prompt_EB)/(num_track_PU200_prompt_EB) << " % | "
       << endl;
  cout << "  muon&&track  | "
       << num_both_muon_track_time_err_PU200_prompt_EB << "/" << num_track_PU200_prompt_EB << " = " << 1.*(num_both_muon_track_time_err_PU200_prompt_EB)/(num_track_PU200_prompt_EB) << " % | "
	   << endl;
  cout << "  vtx&&track   | "
       << num_both_vtx_track_time_err_PU200_prompt_EB << "/" << num_track_PU200_prompt_EB << " = " << 1.*(num_both_vtx_track_time_err_PU200_prompt_EB)/(num_track_PU200_prompt_EB) << " % | "
       << endl;
       cout << endl;
  cout << "[PU200 nonprompt]" << endl;
  cout << "               |           time          |          tErr           |          tErr_i         |           MVA           |    MVA > 0.5 / total    |" << endl;
  cout << "      muon     | "
       << num_muon_time_PU200_nonprompt_EB << "/" << num_muon_PU200_nonprompt_EB << " = " << 1.*(num_muon_time_PU200_nonprompt_EB)/(num_muon_PU200_nonprompt_EB) << " % | "
       << num_muon_time_err_PU200_nonprompt_EB << "/" << num_muon_PU200_nonprompt_EB << " = " << 1.*(num_muon_time_err_PU200_nonprompt_EB)/(num_muon_PU200_nonprompt_EB) << " % | "
       << num_muon_time_err_i_PU200_nonprompt_EB << "/" << num_muon_PU200_nonprompt_EB << " = " << 1.*(num_muon_time_err_i_PU200_nonprompt_EB)/(num_muon_PU200_nonprompt_EB) << " % | "
       << num_muon_mva_PU200_nonprompt_EB << "/" << num_muon_PU200_nonprompt_EB << " = " << 1.*(num_muon_mva_PU200_nonprompt_EB)/(num_muon_PU200_nonprompt_EB) << " % | "
       << num_muon_mva_cut_PU200_nonprompt_EB << "/" << num_muon_PU200_nonprompt_EB << " = " << 1.*(num_muon_mva_cut_PU200_nonprompt_EB)/(num_muon_PU200_nonprompt_EB) << " % | "
       << endl;
  cout << "     vertex    | "
       << num_vtx_time_PU200_nonprompt_EB << "/" << num_vtx_PU200_nonprompt_EB << " = " << 1.*(num_vtx_time_PU200_nonprompt_EB)/(num_vtx_PU200_nonprompt_EB) << " % | "
       << num_vtx_time_err_PU200_nonprompt_EB << "/" << num_vtx_PU200_nonprompt_EB << " = " << 1.*(num_vtx_time_err_PU200_nonprompt_EB)/(num_vtx_PU200_nonprompt_EB) << " % | "
	   << endl;
  cout << "      track    | "
       << num_track_time_PU200_nonprompt_EB << "/" << num_track_PU200_nonprompt_EB << " = " << 1.*(num_track_time_PU200_nonprompt_EB)/(num_track_PU200_nonprompt_EB) << " % | "
       << num_track_time_err_PU200_nonprompt_EB << "/" << num_track_PU200_nonprompt_EB << " = " << 1.*(num_track_time_err_PU200_nonprompt_EB)/(num_track_PU200_nonprompt_EB) << " % | "
       << num_track_time_err_i_PU200_nonprompt_EB << "/" << num_track_PU200_nonprompt_EB << " = " << 1.*(num_track_time_err_i_PU200_nonprompt_EB)/(num_track_PU200_nonprompt_EB) << " % | "
       << num_track_mva_PU200_nonprompt_EB << "/" << num_track_PU200_nonprompt_EB << " = " << 1.*(num_track_mva_PU200_nonprompt_EB)/(num_track_PU200_nonprompt_EB) << " % | "
       << num_track_mva_cut_PU200_nonprompt_EB << "/" << num_track_PU200_nonprompt_EB << " = " << 1.*(num_track_mva_cut_PU200_nonprompt_EB)/(num_track_PU200_nonprompt_EB) << " % | "
       << endl;
  cout << "  muon&&track  | "
       << num_both_muon_track_time_err_PU200_nonprompt_EB << "/" << num_track_PU200_nonprompt_EB << " = " << 1.*(num_both_muon_track_time_err_PU200_nonprompt_EB)/(num_track_PU200_nonprompt_EB) << " % | "
	   << endl;
  cout << "  vtx&&track   | "
       << num_both_vtx_track_time_err_PU200_nonprompt_EB << "/" << num_track_PU200_nonprompt_EB << " = " << 1.*(num_both_vtx_track_time_err_PU200_nonprompt_EB)/(num_track_PU200_nonprompt_EB) << " % | "
       << endl;
       cout << endl;
  cout << "[noPU prompt]" << endl;
  cout << "               |           time          |          tErr           |          tErr_i         |           MVA           |    MVA > 0.5 / total    |" << endl;
  cout << "      muon     | "
       << num_muon_time_noPU_prompt_EB << "/" << num_muon_noPU_prompt_EB << " = " << 1.*(num_muon_time_noPU_prompt_EB)/(num_muon_noPU_prompt_EB) << " % | "
       << num_muon_time_err_noPU_prompt_EB << "/" << num_muon_noPU_prompt_EB << " = " << 1.*(num_muon_time_err_noPU_prompt_EB)/(num_muon_noPU_prompt_EB) << " % | "
       << num_muon_time_err_i_noPU_prompt_EB << "/" << num_muon_noPU_prompt_EB << " = " << 1.*(num_muon_time_err_i_noPU_prompt_EB)/(num_muon_noPU_prompt_EB) << " % | "
       << num_muon_mva_noPU_prompt_EB << "/" << num_muon_noPU_prompt_EB << " = " << 1.*(num_muon_mva_noPU_prompt_EB)/(num_muon_noPU_prompt_EB) << " % | "
       << num_muon_mva_cut_noPU_prompt_EB << "/" << num_muon_noPU_prompt_EB << " = " << 1.*(num_muon_mva_cut_noPU_prompt_EB)/(num_muon_noPU_prompt_EB) << " % | "
       << endl;
  cout << "     vertex    | "
       << num_vtx_time_noPU_prompt_EB << "/" << num_vtx_noPU_prompt_EB << " = " << 1.*(num_vtx_time_noPU_prompt_EB)/(num_vtx_noPU_prompt_EB) << " % | "
       << num_vtx_time_err_noPU_prompt_EB << "/" << num_vtx_noPU_prompt_EB << " = " << 1.*(num_vtx_time_err_noPU_prompt_EB)/(num_vtx_noPU_prompt_EB) << " % | "
	   << endl;
  cout << "      track    | "
       << num_track_time_noPU_prompt_EB << "/" << num_track_noPU_prompt_EB << " = " << 1.*(num_track_time_noPU_prompt_EB)/(num_track_noPU_prompt_EB) << " % | "
       << num_track_time_err_noPU_prompt_EB << "/" << num_track_noPU_prompt_EB << " = " << 1.*(num_track_time_err_noPU_prompt_EB)/(num_track_noPU_prompt_EB) << " % | "
       << num_track_time_err_i_noPU_prompt_EB << "/" << num_track_noPU_prompt_EB << " = " << 1.*(num_track_time_err_i_noPU_prompt_EB)/(num_track_noPU_prompt_EB) << " % | "
       << num_track_mva_noPU_prompt_EB << "/" << num_track_noPU_prompt_EB << " = " << 1.*(num_track_mva_noPU_prompt_EB)/(num_track_noPU_prompt_EB) << " % | "
       << num_track_mva_cut_noPU_prompt_EB << "/" << num_track_noPU_prompt_EB << " = " << 1.*(num_track_mva_cut_noPU_prompt_EB)/(num_track_noPU_prompt_EB) << " % | "
       << endl;
  cout << "  muon&&track  | "
       << num_both_muon_track_time_err_noPU_prompt_EB << "/" << num_track_noPU_prompt_EB << " = " << 1.*(num_both_muon_track_time_err_noPU_prompt_EB)/(num_track_noPU_prompt_EB) << " % | "
	   << endl;
  cout << "  vtx&&track   | "
       << num_both_vtx_track_time_err_noPU_prompt_EB << "/" << num_track_noPU_prompt_EB << " = " << 1.*(num_both_vtx_track_time_err_noPU_prompt_EB)/(num_track_noPU_prompt_EB) << " % | "
       << endl;
       cout << endl;
  cout << "[noPU nonprompt]" << endl;
  cout << "               |           time          |          tErr           |          tErr_i         |           MVA           |    MVA > 0.5 / total    |" << endl;
  cout << "      muon     | "
       << num_muon_time_noPU_nonprompt_EB << "/" << num_muon_noPU_nonprompt_EB << " = " << 1.*(num_muon_time_noPU_nonprompt_EB)/(num_muon_noPU_nonprompt_EB) << " % | "
       << num_muon_time_err_noPU_nonprompt_EB << "/" << num_muon_noPU_nonprompt_EB << " = " << 1.*(num_muon_time_err_noPU_nonprompt_EB)/(num_muon_noPU_nonprompt_EB) << " % | "
       << num_muon_time_err_i_noPU_nonprompt_EB << "/" << num_muon_noPU_nonprompt_EB << " = " << 1.*(num_muon_time_err_i_noPU_nonprompt_EB)/(num_muon_noPU_nonprompt_EB) << " % | "
       << num_muon_mva_noPU_nonprompt_EB << "/" << num_muon_noPU_nonprompt_EB << " = " << 1.*(num_muon_mva_noPU_nonprompt_EB)/(num_muon_noPU_nonprompt_EB) << " % | "
       << num_muon_mva_cut_noPU_nonprompt_EB << "/" << num_muon_noPU_nonprompt_EB << " = " << 1.*(num_muon_mva_cut_noPU_nonprompt_EB)/(num_muon_noPU_nonprompt_EB) << " % | "
       << endl;
  cout << "     vertex    | "
       << num_vtx_time_noPU_nonprompt_EB << "/" << num_vtx_noPU_nonprompt_EB << " = " << 1.*(num_vtx_time_noPU_nonprompt_EB)/(num_vtx_noPU_nonprompt_EB) << " % | "
       << num_vtx_time_err_noPU_nonprompt_EB << "/" << num_vtx_noPU_nonprompt_EB << " = " << 1.*(num_vtx_time_err_noPU_nonprompt_EB)/(num_vtx_noPU_nonprompt_EB) << " % | "
	   << endl;
  cout << "      track    | "
       << num_track_time_noPU_nonprompt_EB << "/" << num_track_noPU_nonprompt_EB << " = " << 1.*(num_track_time_noPU_nonprompt_EB)/(num_track_noPU_nonprompt_EB) << " % | "
       << num_track_time_err_noPU_nonprompt_EB << "/" << num_track_noPU_nonprompt_EB << " = " << 1.*(num_track_time_err_noPU_nonprompt_EB)/(num_track_noPU_nonprompt_EB) << " % | "
       << num_track_time_err_i_noPU_nonprompt_EB << "/" << num_track_noPU_nonprompt_EB << " = " << 1.*(num_track_time_err_i_noPU_nonprompt_EB)/(num_track_noPU_nonprompt_EB) << " % | "
       << num_track_mva_noPU_nonprompt_EB << "/" << num_track_noPU_nonprompt_EB << " = " << 1.*(num_track_mva_noPU_nonprompt_EB)/(num_track_noPU_nonprompt_EB) << " % | "
       << num_track_mva_cut_noPU_nonprompt_EB << "/" << num_track_noPU_nonprompt_EB << " = " << 1.*(num_track_mva_cut_noPU_nonprompt_EB)/(num_track_noPU_nonprompt_EB) << " % | "
       << endl;
  cout << "  muon&&track  | "
       << num_both_muon_track_time_err_noPU_nonprompt_EB << "/" << num_track_noPU_nonprompt_EB << " = " << 1.*(num_both_muon_track_time_err_noPU_nonprompt_EB)/(num_track_noPU_nonprompt_EB) << " % | "
	   << endl;
  cout << "  vtx&&track   | "
       << num_both_vtx_track_time_err_noPU_nonprompt_EB << "/" << num_track_noPU_nonprompt_EB << " = " << 1.*(num_both_vtx_track_time_err_noPU_nonprompt_EB)/(num_track_noPU_nonprompt_EB) << " % | "
       << endl;
       cout << endl;

  cout << endl;
  cout << endl;
  // Endcap
  cout << "////////" << endl;
  cout << "[Endcap]" << endl;
  cout << "////////" << endl;
  cout << "[PU200 prompt]" << endl;
  cout << "               |           time          |          tErr           |          tErr_i         |           MVA           |    MVA > 0.5 / total    |" << endl;
  cout << "      muon     | "
       << num_muon_time_PU200_prompt_EE << "/" << num_muon_PU200_prompt_EE << " = " << 1.*(num_muon_time_PU200_prompt_EE)/(num_muon_PU200_prompt_EE) << " % | "
       << num_muon_time_err_PU200_prompt_EE << "/" << num_muon_PU200_prompt_EE << " = " << 1.*(num_muon_time_err_PU200_prompt_EE)/(num_muon_PU200_prompt_EE) << " % | "
       << num_muon_time_err_i_PU200_prompt_EE << "/" << num_muon_PU200_prompt_EE << " = " << 1.*(num_muon_time_err_i_PU200_prompt_EE)/(num_muon_PU200_prompt_EE) << " % | "
       << num_muon_mva_PU200_prompt_EE << "/" << num_muon_PU200_prompt_EE << " = " << 1.*(num_muon_mva_PU200_prompt_EE)/(num_muon_PU200_prompt_EE) << " % | "
       << num_muon_mva_cut_PU200_prompt_EE << "/" << num_muon_PU200_prompt_EE << " = " << 1.*(num_muon_mva_cut_PU200_prompt_EE)/(num_muon_PU200_prompt_EE) << " % | "
       << endl;
  cout << "     vertex    | "
       << num_vtx_time_PU200_prompt_EE << "/" << num_vtx_PU200_prompt_EE << " = " << 1.*(num_vtx_time_PU200_prompt_EE)/(num_vtx_PU200_prompt_EE) << " % | "
       << num_vtx_time_err_PU200_prompt_EE << "/" << num_vtx_PU200_prompt_EE << " = " << 1.*(num_vtx_time_err_PU200_prompt_EE)/(num_vtx_PU200_prompt_EE) << " % | "
	   << endl;
  cout << "      track    | "
       << num_track_time_PU200_prompt_EE << "/" << num_track_PU200_prompt_EE << " = " << 1.*(num_track_time_PU200_prompt_EE)/(num_track_PU200_prompt_EE) << " % | "
       << num_track_time_err_PU200_prompt_EE << "/" << num_track_PU200_prompt_EE << " = " << 1.*(num_track_time_err_PU200_prompt_EE)/(num_track_PU200_prompt_EE) << " % | "
       << num_track_time_err_i_PU200_prompt_EE << "/" << num_track_PU200_prompt_EE << " = " << 1.*(num_track_time_err_i_PU200_prompt_EE)/(num_track_PU200_prompt_EE) << " % | "
       << num_track_mva_PU200_prompt_EE << "/" << num_track_PU200_prompt_EE << " = " << 1.*(num_track_mva_PU200_prompt_EE)/(num_track_PU200_prompt_EE) << " % | "
       << num_track_mva_cut_PU200_prompt_EE << "/" << num_track_PU200_prompt_EE << " = " << 1.*(num_track_mva_cut_PU200_prompt_EE)/(num_track_PU200_prompt_EE) << " % | "
       << endl;
  cout << "  muon&&track  | "
       << num_both_muon_track_time_err_PU200_prompt_EE << "/" << num_track_PU200_prompt_EE << " = " << 1.*(num_both_muon_track_time_err_PU200_prompt_EE)/(num_track_PU200_prompt_EE) << " % | "
	   << endl;
  cout << "  vtx&&track   | "
       << num_both_vtx_track_time_err_PU200_prompt_EE << "/" << num_track_PU200_prompt_EE << " = " << 1.*(num_both_vtx_track_time_err_PU200_prompt_EE)/(num_track_PU200_prompt_EE) << " % | "
       << endl;
       cout << endl;
  cout << "[PU200 nonprompt]" << endl;
  cout << "               |           time          |          tErr           |          tErr_i         |           MVA           |    MVA > 0.5 / total    |" << endl;
  cout << "      muon     | "
       << num_muon_time_PU200_nonprompt_EE << "/" << num_muon_PU200_nonprompt_EE << " = " << 1.*(num_muon_time_PU200_nonprompt_EE)/(num_muon_PU200_nonprompt_EE) << " % | "
       << num_muon_time_err_PU200_nonprompt_EE << "/" << num_muon_PU200_nonprompt_EE << " = " << 1.*(num_muon_time_err_PU200_nonprompt_EE)/(num_muon_PU200_nonprompt_EE) << " % | "
       << num_muon_time_err_i_PU200_nonprompt_EE << "/" << num_muon_PU200_nonprompt_EE << " = " << 1.*(num_muon_time_err_i_PU200_nonprompt_EE)/(num_muon_PU200_nonprompt_EE) << " % | "
       << num_muon_mva_PU200_nonprompt_EE << "/" << num_muon_PU200_nonprompt_EE << " = " << 1.*(num_muon_mva_PU200_nonprompt_EE)/(num_muon_PU200_nonprompt_EE) << " % | "
       << num_muon_mva_cut_PU200_nonprompt_EE << "/" << num_muon_PU200_nonprompt_EE << " = " << 1.*(num_muon_mva_cut_PU200_nonprompt_EE)/(num_muon_PU200_nonprompt_EE) << " % | "
       << endl;
  cout << "     vertex    | "
       << num_vtx_time_PU200_nonprompt_EE << "/" << num_vtx_PU200_nonprompt_EE << " = " << 1.*(num_vtx_time_PU200_nonprompt_EE)/(num_vtx_PU200_nonprompt_EE) << " % | "
       << num_vtx_time_err_PU200_nonprompt_EE << "/" << num_vtx_PU200_nonprompt_EE << " = " << 1.*(num_vtx_time_err_PU200_nonprompt_EE)/(num_vtx_PU200_nonprompt_EE) << " % | "
	   << endl;
  cout << "      track    | "
       << num_track_time_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EE << " = " << 1.*(num_track_time_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EE) << " % | "
       << num_track_time_err_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EE << " = " << 1.*(num_track_time_err_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EE) << " % | "
       << num_track_time_err_i_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EE << " = " << 1.*(num_track_time_err_i_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EE) << " % | "
       << num_track_mva_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EE << " = " << 1.*(num_track_mva_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EE) << " % | "
       << num_track_mva_cut_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EE << " = " << 1.*(num_track_mva_cut_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EE) << " % | "
       << endl;
  cout << "  muon&&track  | "
       << num_both_muon_track_time_err_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EE << " = " << 1.*(num_both_muon_track_time_err_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EE) << " % | "
	   << endl;
  cout << "  vtx&&track   | "
       << num_both_vtx_track_time_err_PU200_nonprompt_EE << "/" << num_track_PU200_nonprompt_EE << " = " << 1.*(num_both_vtx_track_time_err_PU200_nonprompt_EE)/(num_track_PU200_nonprompt_EE) << " % | "
       << endl;
       cout << endl;
  cout << "[noPU prompt]" << endl;
  cout << "               |           time          |          tErr           |          tErr_i         |           MVA           |    MVA > 0.5 / total    |" << endl;
  cout << "      muon     | "
       << num_muon_time_noPU_prompt_EE << "/" << num_muon_noPU_prompt_EE << " = " << 1.*(num_muon_time_noPU_prompt_EE)/(num_muon_noPU_prompt_EE) << " % | "
       << num_muon_time_err_noPU_prompt_EE << "/" << num_muon_noPU_prompt_EE << " = " << 1.*(num_muon_time_err_noPU_prompt_EE)/(num_muon_noPU_prompt_EE) << " % | "
       << num_muon_time_err_i_noPU_prompt_EE << "/" << num_muon_noPU_prompt_EE << " = " << 1.*(num_muon_time_err_i_noPU_prompt_EE)/(num_muon_noPU_prompt_EE) << " % | "
       << num_muon_mva_noPU_prompt_EE << "/" << num_muon_noPU_prompt_EE << " = " << 1.*(num_muon_mva_noPU_prompt_EE)/(num_muon_noPU_prompt_EE) << " % | "
       << num_muon_mva_cut_noPU_prompt_EE << "/" << num_muon_noPU_prompt_EE << " = " << 1.*(num_muon_mva_cut_noPU_prompt_EE)/(num_muon_noPU_prompt_EE) << " % | "
       << endl;
  cout << "     vertex    | "
       << num_vtx_time_noPU_prompt_EE << "/" << num_vtx_noPU_prompt_EE << " = " << 1.*(num_vtx_time_noPU_prompt_EE)/(num_vtx_noPU_prompt_EE) << " % | "
       << num_vtx_time_err_noPU_prompt_EE << "/" << num_vtx_noPU_prompt_EE << " = " << 1.*(num_vtx_time_err_noPU_prompt_EE)/(num_vtx_noPU_prompt_EE) << " % | "
	   << endl;
  cout << "      track    | "
       << num_track_time_noPU_prompt_EE << "/" << num_track_noPU_prompt_EE << " = " << 1.*(num_track_time_noPU_prompt_EE)/(num_track_noPU_prompt_EE) << " % | "
       << num_track_time_err_noPU_prompt_EE << "/" << num_track_noPU_prompt_EE << " = " << 1.*(num_track_time_err_noPU_prompt_EE)/(num_track_noPU_prompt_EE) << " % | "
       << num_track_time_err_i_noPU_prompt_EE << "/" << num_track_noPU_prompt_EE << " = " << 1.*(num_track_time_err_i_noPU_prompt_EE)/(num_track_noPU_prompt_EE) << " % | "
       << num_track_mva_noPU_prompt_EE << "/" << num_track_noPU_prompt_EE << " = " << 1.*(num_track_mva_noPU_prompt_EE)/(num_track_noPU_prompt_EE) << " % | "
       << num_track_mva_cut_noPU_prompt_EE << "/" << num_track_noPU_prompt_EE << " = " << 1.*(num_track_mva_cut_noPU_prompt_EE)/(num_track_noPU_prompt_EE) << " % | "
       << endl;
  cout << "  muon&&track  | "
       << num_both_muon_track_time_err_noPU_prompt_EE << "/" << num_track_noPU_prompt_EE << " = " << 1.*(num_both_muon_track_time_err_noPU_prompt_EE)/(num_track_noPU_prompt_EE) << " % | "
	   << endl;
  cout << "  vtx&&track   | "
       << num_both_vtx_track_time_err_noPU_prompt_EE << "/" << num_track_noPU_prompt_EE << " = " << 1.*(num_both_vtx_track_time_err_noPU_prompt_EE)/(num_track_noPU_prompt_EE) << " % | "
       << endl;
       cout << endl;
  cout << "[noPU nonprompt]" << endl;
  cout << "               |           time          |          tErr           |          tErr_i         |           MVA           |    MVA > 0.5 / total    |" << endl;
  cout << "      muon     | "
       << num_muon_time_noPU_nonprompt_EE << "/" << num_muon_noPU_nonprompt_EE << " = " << 1.*(num_muon_time_noPU_nonprompt_EE)/(num_muon_noPU_nonprompt_EE) << " % | "
       << num_muon_time_err_noPU_nonprompt_EE << "/" << num_muon_noPU_nonprompt_EE << " = " << 1.*(num_muon_time_err_noPU_nonprompt_EE)/(num_muon_noPU_nonprompt_EE) << " % | "
       << num_muon_time_err_i_noPU_nonprompt_EE << "/" << num_muon_noPU_nonprompt_EE << " = " << 1.*(num_muon_time_err_i_noPU_nonprompt_EE)/(num_muon_noPU_nonprompt_EE) << " % | "
       << num_muon_mva_noPU_nonprompt_EE << "/" << num_muon_noPU_nonprompt_EE << " = " << 1.*(num_muon_mva_noPU_nonprompt_EE)/(num_muon_noPU_nonprompt_EE) << " % | "
       << num_muon_mva_cut_noPU_nonprompt_EE << "/" << num_muon_noPU_nonprompt_EE << " = " << 1.*(num_muon_mva_cut_noPU_nonprompt_EE)/(num_muon_noPU_nonprompt_EE) << " % | "
       << endl;
  cout << "     vertex    | "
       << num_vtx_time_noPU_nonprompt_EE << "/" << num_vtx_noPU_nonprompt_EE << " = " << 1.*(num_vtx_time_noPU_nonprompt_EE)/(num_vtx_noPU_nonprompt_EE) << " % | "
       << num_vtx_time_err_noPU_nonprompt_EE << "/" << num_vtx_noPU_nonprompt_EE << " = " << 1.*(num_vtx_time_err_noPU_nonprompt_EE)/(num_vtx_noPU_nonprompt_EE) << " % | "
	   << endl;
  cout << "      track    | "
       << num_track_time_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EE << " = " << 1.*(num_track_time_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EE) << " % | "
       << num_track_time_err_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EE << " = " << 1.*(num_track_time_err_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EE) << " % | "
       << num_track_time_err_i_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EE << " = " << 1.*(num_track_time_err_i_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EE) << " % | "
       << num_track_mva_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EE << " = " << 1.*(num_track_mva_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EE) << " % | "
       << num_track_mva_cut_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EE << " = " << 1.*(num_track_mva_cut_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EE) << " % | "
       << endl;
  cout << "  muon&&track  | "
       << num_both_muon_track_time_err_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EE << " = " << 1.*(num_both_muon_track_time_err_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EE) << " % | "
	   << endl;
  cout << "  vtx&&track   | "
       << num_both_vtx_track_time_err_noPU_nonprompt_EE << "/" << num_track_noPU_nonprompt_EE << " = " << 1.*(num_both_vtx_track_time_err_noPU_nonprompt_EE)/(num_track_noPU_nonprompt_EE) << " % | "
       << endl;
       cout << endl;


  
  ///////////////
  // Cosmetics //
  ///////////////
  // (muon, track)
  h_PU200_prompt_PV_EB->SetFillColor(kGray+1); h_PU200_prompt_PU_EB->SetFillColor(kAzure+7); h_PU200_prompt_fake_EB->SetFillColor(kGreen+2); h_PU200_prompt_no_tErr_EB->SetFillColor(kYellow-7); h_PU200_prompt_SV_EB->SetFillColor(kOrange+7); h_PU200_prompt_not_vtx_matched_EB->SetFillColor(kMagenta);
  h_PU200_prompt_PV_EE->SetFillColor(kGray+1); h_PU200_prompt_PU_EE->SetFillColor(kAzure+7); h_PU200_prompt_fake_EE->SetFillColor(kGreen+2); h_PU200_prompt_no_tErr_EE->SetFillColor(kYellow-7); h_PU200_prompt_SV_EE->SetFillColor(kOrange+7); h_PU200_prompt_not_vtx_matched_EE->SetFillColor(kMagenta);
  h_PU200_nonprompt_PV_EB->SetFillColor(kGray+1); h_PU200_nonprompt_PU_EB->SetFillColor(kAzure+7); h_PU200_nonprompt_fake_EB->SetFillColor(kGreen+2); h_PU200_nonprompt_no_tErr_EB->SetFillColor(kYellow-7); h_PU200_nonprompt_SV_EB->SetFillColor(kOrange+7); h_PU200_nonprompt_not_vtx_matched_EB->SetFillColor(kMagenta);
  h_PU200_nonprompt_PV_EE->SetFillColor(kGray+1); h_PU200_nonprompt_PU_EE->SetFillColor(kAzure+7); h_PU200_nonprompt_fake_EE->SetFillColor(kGreen+2); h_PU200_nonprompt_no_tErr_EE->SetFillColor(kYellow-7); h_PU200_nonprompt_SV_EE->SetFillColor(kOrange+7); h_PU200_nonprompt_not_vtx_matched_EE->SetFillColor(kMagenta);
  h_noPU_prompt_PV_EB->SetFillColor(kGray+1); h_noPU_prompt_PU_EB->SetFillColor(kAzure+7); h_noPU_prompt_fake_EB->SetFillColor(kGreen+2); h_noPU_prompt_no_tErr_EB->SetFillColor(kYellow-7); h_noPU_prompt_SV_EB->SetFillColor(kOrange+7); h_noPU_prompt_not_vtx_matched_EB->SetFillColor(kMagenta);
  h_noPU_prompt_PV_EE->SetFillColor(kGray+1); h_noPU_prompt_PU_EE->SetFillColor(kAzure+7); h_noPU_prompt_fake_EE->SetFillColor(kGreen+2); h_noPU_prompt_no_tErr_EE->SetFillColor(kYellow-7); h_noPU_prompt_SV_EE->SetFillColor(kOrange+7); h_noPU_prompt_not_vtx_matched_EE->SetFillColor(kMagenta);
  h_noPU_nonprompt_PV_EB->SetFillColor(kGray+1); h_noPU_nonprompt_PU_EB->SetFillColor(kAzure+7); h_noPU_nonprompt_fake_EB->SetFillColor(kGreen+2); h_noPU_nonprompt_no_tErr_EB->SetFillColor(kYellow-7); h_noPU_nonprompt_SV_EB->SetFillColor(kOrange+7); h_noPU_nonprompt_not_vtx_matched_EB->SetFillColor(kMagenta);
  h_noPU_nonprompt_PV_EE->SetFillColor(kGray+1); h_noPU_nonprompt_PU_EE->SetFillColor(kAzure+7); h_noPU_nonprompt_fake_EE->SetFillColor(kGreen+2); h_noPU_nonprompt_no_tErr_EE->SetFillColor(kYellow-7); h_noPU_nonprompt_SV_EE->SetFillColor(kOrange+7); h_noPU_nonprompt_not_vtx_matched_EE->SetFillColor(kMagenta);

  // (PV, track)
  h_PU200_prompt_PV_EB_vtx->SetFillColor(kGray+1); h_PU200_prompt_PU_EB_vtx->SetFillColor(kAzure+7); h_PU200_prompt_fake_EB_vtx->SetFillColor(kGreen+2); h_PU200_prompt_no_tErr_EB_vtx->SetFillColor(kYellow-7); h_PU200_prompt_SV_EB_vtx->SetFillColor(kOrange+7);
  h_PU200_prompt_PV_EE_vtx->SetFillColor(kGray+1); h_PU200_prompt_PU_EE_vtx->SetFillColor(kAzure+7); h_PU200_prompt_fake_EE_vtx->SetFillColor(kGreen+2); h_PU200_prompt_no_tErr_EE_vtx->SetFillColor(kYellow-7); h_PU200_prompt_SV_EE_vtx->SetFillColor(kOrange+7);
  h_PU200_nonprompt_PV_EB_vtx->SetFillColor(kGray+1); h_PU200_nonprompt_PU_EB_vtx->SetFillColor(kAzure+7); h_PU200_nonprompt_fake_EB_vtx->SetFillColor(kGreen+2); h_PU200_nonprompt_no_tErr_EB_vtx->SetFillColor(kYellow-7); h_PU200_nonprompt_SV_EB_vtx->SetFillColor(kOrange+7);
  h_PU200_nonprompt_PV_EE_vtx->SetFillColor(kGray+1); h_PU200_nonprompt_PU_EE_vtx->SetFillColor(kAzure+7); h_PU200_nonprompt_fake_EE_vtx->SetFillColor(kGreen+2); h_PU200_nonprompt_no_tErr_EE_vtx->SetFillColor(kYellow-7); h_PU200_nonprompt_SV_EE_vtx->SetFillColor(kOrange+7);
  h_noPU_prompt_PV_EB_vtx->SetFillColor(kGray+1); h_noPU_prompt_PU_EB_vtx->SetFillColor(kAzure+7); h_noPU_prompt_fake_EB_vtx->SetFillColor(kGreen+2); h_noPU_prompt_no_tErr_EB_vtx->SetFillColor(kYellow-7); h_noPU_prompt_SV_EB_vtx->SetFillColor(kOrange+7);
  h_noPU_prompt_PV_EE_vtx->SetFillColor(kGray+1); h_noPU_prompt_PU_EE_vtx->SetFillColor(kAzure+7); h_noPU_prompt_fake_EE_vtx->SetFillColor(kGreen+2); h_noPU_prompt_no_tErr_EE_vtx->SetFillColor(kYellow-7); h_noPU_prompt_SV_EE_vtx->SetFillColor(kOrange+7);
  h_noPU_nonprompt_PV_EB_vtx->SetFillColor(kGray+1); h_noPU_nonprompt_PU_EB_vtx->SetFillColor(kAzure+7); h_noPU_nonprompt_fake_EB_vtx->SetFillColor(kGreen+2); h_noPU_nonprompt_no_tErr_EB_vtx->SetFillColor(kYellow-7); h_noPU_nonprompt_SV_EB_vtx->SetFillColor(kOrange+7);
  h_noPU_nonprompt_PV_EE_vtx->SetFillColor(kGray+1); h_noPU_nonprompt_PU_EE_vtx->SetFillColor(kAzure+7); h_noPU_nonprompt_fake_EE_vtx->SetFillColor(kGreen+2); h_noPU_nonprompt_no_tErr_EE_vtx->SetFillColor(kYellow-7); h_noPU_nonprompt_SV_EE_vtx->SetFillColor(kOrange+7);


  /////////////////
  //// Legends ////
  /////////////////
  // (muon, track)
  TLegend *leg_PU200_prompt_EB = new TLegend(0.50, 0.62, 0.88, 0.88);
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_PV_EB, "Track from PV", "F");
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_SV_EB, "Track from SV", "F");
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_PU_EB, "PU Track", "F");
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_fake_EB, "Fake Track", "F");
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_not_vtx_matched_EB, "Track without gen-matched vertex", "F");
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_no_tErr_EB, "Track without tErr", "F");
  TLegend *leg_PU200_prompt_EE = new TLegend(0.50, 0.62, 0.88, 0.88);
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_PV_EE, "Track from PV", "F");
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_SV_EE, "Track from SV", "F");
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_PU_EE, "PU Track", "F");
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_fake_EE, "Fake Track", "F");
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_not_vtx_matched_EE, "Track without gen-matched vertex", "F");
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_no_tErr_EE, "Track without tErr", "F");
  TLegend *leg_PU200_nonprompt_EB = new TLegend(0.50, 0.62, 0.88, 0.88);
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_PV_EB, "Track from PV", "F");
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_SV_EB, "Track from SV", "F");
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_PU_EB, "PU Track", "F");
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_fake_EB, "Fake Track", "F");
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_not_vtx_matched_EB, "Track without gen-matched vertex", "F");
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_no_tErr_EB, "Track without tErr", "F");
  TLegend *leg_PU200_nonprompt_EE = new TLegend(0.50, 0.62, 0.88, 0.88);
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_PV_EE, "Track from PV", "F");
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_SV_EE, "Track from SV", "F");
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_PU_EE, "PU Track", "F");
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_fake_EE, "Fake Track", "F");
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_not_vtx_matched_EE, "Track without gen-matched vertex", "F");
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_no_tErr_EE, "Track without tErr", "F");
  TLegend *leg_noPU_prompt_EB = new TLegend(0.50, 0.62, 0.88, 0.88);
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_PV_EB, "Track from PV", "F");
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_SV_EB, "Track from SV", "F");
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_PU_EB, "PU Track", "F");
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_fake_EB, "Fake Track", "F");
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_not_vtx_matched_EB, "Track without gen-matched vertex", "F");
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_no_tErr_EB, "Track without tErr", "F");
  TLegend *leg_noPU_prompt_EE = new TLegend(0.50, 0.62, 0.88, 0.88);
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_PV_EE, "Track from PV", "F");
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_SV_EE, "Track from SV", "F");
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_PU_EE, "PU Track", "F");
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_fake_EE, "Fake Track", "F");
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_not_vtx_matched_EE, "Track without gen-matched vertex", "F");
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_no_tErr_EE, "Track without tErr", "F");
  TLegend *leg_noPU_nonprompt_EB = new TLegend(0.50, 0.62, 0.88, 0.88);
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_PV_EB, "Track from PV", "F");
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_SV_EB, "Track from SV", "F");
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_PU_EB, "PU Track", "F");
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_fake_EB, "Fake Track", "F");
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_not_vtx_matched_EB, "Track without gen-matched vertex", "F");
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_no_tErr_EB, "Track without tErr", "F");
  TLegend *leg_noPU_nonprompt_EE = new TLegend(0.50, 0.62, 0.88, 0.88);
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_PV_EE, "Track from PV", "F");
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_SV_EE, "Track from SV", "F");
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_PU_EE, "PU Track", "F");
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_fake_EE, "Fake Track", "F");
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_not_vtx_matched_EE, "Track without gen-matched vertex", "F");
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_no_tErr_EE, "Track without tErr", "F");
  /*
  // (PV, track)
  TLegend *leg_PU200_prompt_EB_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_prompt_EB_vtx->AddEntry(h_PU200_prompt_PV_EB_vtx, "Track from PV", "F");
  leg_PU200_prompt_EB_vtx->AddEntry(h_PU200_prompt_SV_EB_vtx, "Track from SV", "F");
  leg_PU200_prompt_EB_vtx->AddEntry(h_PU200_prompt_PU_EB_vtx, "PU Track", "F");
  leg_PU200_prompt_EB_vtx->AddEntry(h_PU200_prompt_fake_EB_vtx, "Fake Track", "F");
  leg_PU200_prompt_EB_vtx->AddEntry(h_PU200_prompt_no_tErr_EB_vtx, "Track without tErr", "F");
  TLegend *leg_PU200_prompt_EE_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_prompt_EE_vtx->AddEntry(h_PU200_prompt_PV_EE_vtx, "Track from PV", "F");
  leg_PU200_prompt_EE_vtx->AddEntry(h_PU200_prompt_SV_EE_vtx, "Track from SV", "F");
  leg_PU200_prompt_EE_vtx->AddEntry(h_PU200_prompt_PU_EE_vtx, "PU Track", "F");
  leg_PU200_prompt_EE_vtx->AddEntry(h_PU200_prompt_fake_EE_vtx, "Fake Track", "F");
  leg_PU200_prompt_EE_vtx->AddEntry(h_PU200_prompt_no_tErr_EE_vtx, "Track without tErr", "F");
  TLegend *leg_PU200_nonprompt_EB_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_nonprompt_EB_vtx->AddEntry(h_PU200_nonprompt_PV_EB_vtx, "Track from PV", "F");
  leg_PU200_nonprompt_EB_vtx->AddEntry(h_PU200_nonprompt_SV_EB_vtx, "Track from SV", "F");
  leg_PU200_nonprompt_EB_vtx->AddEntry(h_PU200_nonprompt_PU_EB_vtx, "PU Track", "F");
  leg_PU200_nonprompt_EB_vtx->AddEntry(h_PU200_nonprompt_fake_EB_vtx, "Fake Track", "F");
  leg_PU200_nonprompt_EB_vtx->AddEntry(h_PU200_nonprompt_no_tErr_EB_vtx, "Track without tErr", "F");
  TLegend *leg_PU200_nonprompt_EE_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_nonprompt_EE_vtx->AddEntry(h_PU200_nonprompt_PV_EE_vtx, "Track from PV", "F");
  leg_PU200_nonprompt_EE_vtx->AddEntry(h_PU200_nonprompt_SV_EE_vtx, "Track from SV", "F");
  leg_PU200_nonprompt_EE_vtx->AddEntry(h_PU200_nonprompt_PU_EE_vtx, "PU Track", "F");
  leg_PU200_nonprompt_EE_vtx->AddEntry(h_PU200_nonprompt_fake_EE_vtx, "Fake Track", "F");
  leg_PU200_nonprompt_EE_vtx->AddEntry(h_PU200_nonprompt_no_tErr_EE_vtx, "Track without tErr", "F");
  TLegend *leg_noPU_prompt_EB_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_prompt_EB_vtx->AddEntry(h_noPU_prompt_PV_EB_vtx, "Track from PV", "F");
  leg_noPU_prompt_EB_vtx->AddEntry(h_noPU_prompt_SV_EB_vtx, "Track from SV", "F");
  leg_noPU_prompt_EB_vtx->AddEntry(h_noPU_prompt_PU_EB_vtx, "PU Track", "F");
  leg_noPU_prompt_EB_vtx->AddEntry(h_noPU_prompt_fake_EB_vtx, "Fake Track", "F");
  leg_noPU_prompt_EB_vtx->AddEntry(h_noPU_prompt_no_tErr_EB_vtx, "Track without tErr", "F");
  TLegend *leg_noPU_prompt_EE_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_prompt_EE_vtx->AddEntry(h_noPU_prompt_PV_EE_vtx, "Track from PV", "F");
  leg_noPU_prompt_EE_vtx->AddEntry(h_noPU_prompt_SV_EE_vtx, "Track from SV", "F");
  leg_noPU_prompt_EE_vtx->AddEntry(h_noPU_prompt_PU_EE_vtx, "PU Track", "F");
  leg_noPU_prompt_EE_vtx->AddEntry(h_noPU_prompt_fake_EE_vtx, "Fake Track", "F");
  leg_noPU_prompt_EE_vtx->AddEntry(h_noPU_prompt_no_tErr_EE_vtx, "Track without tErr", "F");
  TLegend *leg_noPU_nonprompt_EB_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_nonprompt_EB_vtx->AddEntry(h_noPU_nonprompt_PV_EB_vtx, "Track from PV", "F");
  leg_noPU_nonprompt_EB_vtx->AddEntry(h_noPU_nonprompt_SV_EB_vtx, "Track from SV", "F");
  leg_noPU_nonprompt_EB_vtx->AddEntry(h_noPU_nonprompt_PU_EB_vtx, "PU Track", "F");
  leg_noPU_nonprompt_EB_vtx->AddEntry(h_noPU_nonprompt_fake_EB_vtx, "Fake Track", "F");
  leg_noPU_nonprompt_EB_vtx->AddEntry(h_noPU_nonprompt_no_tErr_EB_vtx, "Track without tErr", "F");
  TLegend *leg_noPU_nonprompt_EE_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_nonprompt_EE_vtx->AddEntry(h_noPU_nonprompt_PV_EE_vtx, "Track from PV", "F");
  leg_noPU_nonprompt_EE_vtx->AddEntry(h_noPU_nonprompt_SV_EE_vtx, "Track from SV", "F");
  leg_noPU_nonprompt_EE_vtx->AddEntry(h_noPU_nonprompt_PU_EE_vtx, "PU Track", "F");
  leg_noPU_nonprompt_EE_vtx->AddEntry(h_noPU_nonprompt_fake_EE_vtx, "Fake Track", "F");
  leg_noPU_nonprompt_EE_vtx->AddEntry(h_noPU_nonprompt_no_tErr_EE_vtx, "Track without tErr", "F");
  */


  /////////////////
  // Stack plots //
  /////////////////
  // (muon, track)
  THStack* st_PU200_prompt_EB    = new THStack("st_PU200_prompt_EB","st_PU200_prompt_EB");
  THStack* st_PU200_prompt_EE    = new THStack("st_PU200_prompt_EE","st_PU200_prompt_EE");
  THStack* st_PU200_nonprompt_EB = new THStack("st_PU200_nonprompt_EB","st_PU200_nonprompt_EB");
  THStack* st_PU200_nonprompt_EE = new THStack("st_PU200_nonprompt_EE","st_PU200_nonprompt_EE");
  THStack* st_noPU_prompt_EB     = new THStack("st_noPU_prompt_EB","st_noPU_prompt_EB");
  THStack* st_noPU_prompt_EE     = new THStack("st_noPU_prompt_EE","st_noPU_prompt_EE");
  THStack* st_noPU_nonprompt_EB  = new THStack("st_noPU_nonprompt_EB","st_noPU_nonprompt_EB");
  THStack* st_noPU_nonprompt_EE  = new THStack("st_noPU_nonprompt_EE","st_noPU_nonprompt_EE");

  st_PU200_prompt_EB->Add(h_PU200_prompt_PV_EB); st_PU200_prompt_EB->Add(h_PU200_prompt_SV_EB); st_PU200_prompt_EB->Add(h_PU200_prompt_PU_EB); st_PU200_prompt_EB->Add(h_PU200_prompt_fake_EB); st_PU200_prompt_EB->Add(h_PU200_prompt_not_vtx_matched_EB); st_PU200_prompt_EB->Add(h_PU200_prompt_no_tErr_EB); st_PU200_prompt_EB->SetTitle(";#sigma_{t} (muon, track);Counts");
  st_PU200_prompt_EE->Add(h_PU200_prompt_PV_EE); st_PU200_prompt_EE->Add(h_PU200_prompt_SV_EE); st_PU200_prompt_EE->Add(h_PU200_prompt_PU_EE); st_PU200_prompt_EE->Add(h_PU200_prompt_fake_EE); st_PU200_prompt_EE->Add(h_PU200_prompt_not_vtx_matched_EE); st_PU200_prompt_EE->Add(h_PU200_prompt_no_tErr_EE); st_PU200_prompt_EE->SetTitle(";#sigma_{t} (muon, track);Counts");
  st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_PV_EB); st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_SV_EB); st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_PU_EB); st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_fake_EB); st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_not_vtx_matched_EB); st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_no_tErr_EB); st_PU200_nonprompt_EB->SetTitle(";#sigma_{t} (muon, track);Counts");
  st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_PV_EE); st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_SV_EE); st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_PU_EE); st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_fake_EE); st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_not_vtx_matched_EE); st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_no_tErr_EE); st_PU200_nonprompt_EE->SetTitle(";#sigma_{t} (muon, track);Counts");
  st_noPU_prompt_EB->Add(h_noPU_prompt_PV_EB); st_noPU_prompt_EB->Add(h_noPU_prompt_SV_EB); st_noPU_prompt_EB->Add(h_noPU_prompt_PU_EB); st_noPU_prompt_EB->Add(h_noPU_prompt_fake_EB); st_noPU_prompt_EB->Add(h_noPU_prompt_not_vtx_matched_EB); st_noPU_prompt_EB->Add(h_noPU_prompt_no_tErr_EB); st_noPU_prompt_EB->SetTitle(";#sigma_{t} (muon, track);Counts");
  st_noPU_prompt_EE->Add(h_noPU_prompt_PV_EE); st_noPU_prompt_EE->Add(h_noPU_prompt_SV_EE); st_noPU_prompt_EE->Add(h_noPU_prompt_PU_EE); st_noPU_prompt_EE->Add(h_noPU_prompt_fake_EE); st_noPU_prompt_EE->Add(h_noPU_prompt_not_vtx_matched_EE); st_noPU_prompt_EE->Add(h_noPU_prompt_no_tErr_EE); st_noPU_prompt_EE->SetTitle(";#sigma_{t} (muon, track);Counts");
  st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_PV_EB); st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_SV_EB); st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_PU_EB); st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_fake_EB); st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_not_vtx_matched_EB); st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_no_tErr_EB); st_noPU_nonprompt_EB->SetTitle(";#sigma_{t} (muon, track);Counts");
  st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_PV_EE); st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_SV_EE); st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_PU_EE); st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_fake_EE); st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_not_vtx_matched_EE); st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_no_tErr_EE); st_noPU_nonprompt_EE->SetTitle(";#sigma_{t} (muon, track);Counts");


  // (PV, track)
  THStack* st_PU200_prompt_EB_vtx    = new THStack("st_PU200_prompt_EB_vtx","st_PU200_prompt_EB_vtx");
  THStack* st_PU200_prompt_EE_vtx    = new THStack("st_PU200_prompt_EE_vtx","st_PU200_prompt_EE_vtx");
  THStack* st_PU200_nonprompt_EB_vtx = new THStack("st_PU200_nonprompt_EB_vtx","st_PU200_nonprompt_EB_vtx");
  THStack* st_PU200_nonprompt_EE_vtx = new THStack("st_PU200_nonprompt_EE_vtx","st_PU200_nonprompt_EE_vtx");
  THStack* st_noPU_prompt_EB_vtx     = new THStack("st_noPU_prompt_EB_vtx","st_noPU_prompt_EB_vtx");
  THStack* st_noPU_prompt_EE_vtx     = new THStack("st_noPU_prompt_EE_vtx","st_noPU_prompt_EE_vtx");
  THStack* st_noPU_nonprompt_EB_vtx  = new THStack("st_noPU_nonprompt_EB_vtx","st_noPU_nonprompt_EB_vtx");
  THStack* st_noPU_nonprompt_EE_vtx  = new THStack("st_noPU_nonprompt_EE_vtx","st_noPU_nonprompt_EE_vtx");
  st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_PV_EB_vtx); st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_PU_EB_vtx); st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_fake_EB_vtx); st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_SV_EB_vtx); st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_no_tErr_EB_vtx); st_PU200_prompt_EB_vtx->SetTitle(";#sigma_{t} (PV, track);Counts");
  st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_PV_EE_vtx); st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_PU_EE_vtx); st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_fake_EE_vtx); st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_SV_EE_vtx); st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_no_tErr_EE_vtx); st_PU200_prompt_EE_vtx->SetTitle(";#sigma_{t} (PV, track);Counts");
  st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_PV_EB_vtx); st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_PU_EB_vtx); st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_fake_EB_vtx); st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_SV_EB_vtx); st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_no_tErr_EB_vtx); st_PU200_nonprompt_EB_vtx->SetTitle(";#sigma_{t} (PV, track);Counts");
  st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_PV_EE_vtx); st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_PU_EE_vtx); st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_fake_EE_vtx); st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_SV_EE_vtx); st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_no_tErr_EE_vtx); st_PU200_nonprompt_EE_vtx->SetTitle(";#sigma_{t} (PV, track);Counts");
  st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_PV_EB_vtx); st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_PU_EB_vtx); st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_fake_EB_vtx); st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_SV_EB_vtx); st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_no_tErr_EB_vtx); st_noPU_prompt_EB_vtx->SetTitle(";#sigma_{t} (PV, track);Counts");
  st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_PV_EE_vtx); st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_PU_EE_vtx); st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_fake_EE_vtx); st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_SV_EE_vtx); st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_no_tErr_EE_vtx); st_noPU_prompt_EE_vtx->SetTitle(";#sigma_{t} (PV, track);Counts");
  st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_PV_EB_vtx); st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_PU_EB_vtx); st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_fake_EB_vtx); st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_SV_EB_vtx); st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_no_tErr_EB_vtx); st_noPU_nonprompt_EB_vtx->SetTitle(";#sigma_{t} (PV, track);Counts");
  st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_PV_EE_vtx); st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_PU_EE_vtx); st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_fake_EE_vtx); st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_SV_EE_vtx); st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_no_tErr_EE_vtx); st_noPU_nonprompt_EE_vtx->SetTitle(";#sigma_{t} (PV, track);Counts");

  ////////////////
  // Draw plots //
  ////////////////
  // (muon, track)
  // PU200
    // prompt
      // Barrel
  TCanvas* c_st_PU200_prompt_EB_v2 = new TCanvas("c_st_PU200_prompt_EB_v2", "c_st_PU200_prompt_EB_v2", 1500, 1500);
  c_st_PU200_prompt_EB_v2->cd();
  c_st_PU200_prompt_EB_v2->SetLeftMargin(0.12);
  st_PU200_prompt_EB->Draw();
  leg_PU200_prompt_EB->Draw();
  c_st_PU200_prompt_EB_v2->Print("plots/ntuple/track_sigma_PU200_prompt_EB.pdf");
  cout << "///////////////////////" << endl;
  cout << "///// muon, track /////" << endl;
  cout << "///////////////////////" << endl;
  cout << "/////    2sigma   /////" << endl;
  cout << "///////////////////////" << endl;
  cout << "PU200 Barrel prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EB->Integral(1,-1) << " : " << h_PU200_prompt_PV_EB->Integral(1,2) << " : " << h_PU200_prompt_PV_EB->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EB->Integral(3,-1)/h_PU200_prompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EB->Integral(1,-1) << " : " << h_PU200_prompt_SV_EB->Integral(1,2) << " : " << h_PU200_prompt_SV_EB->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EB->Integral(3,-1)/h_PU200_prompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EB->Integral(1,-1) << " : " << h_PU200_prompt_PU_EB->Integral(1,2) << " : " << h_PU200_prompt_PU_EB->Integral(3,-1) << "  " << "(" << 100 -h_PU200_prompt_PU_EB->Integral(3,-1)/h_PU200_prompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EB->Integral(1,-1) << " : " << h_PU200_prompt_fake_EB->Integral(1,2) << " : " << h_PU200_prompt_fake_EB->Integral(3,-1) << "  " << "(" << 100 -h_PU200_prompt_fake_EB->Integral(3,-1)/h_PU200_prompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_prompt_not_vtx_matched_EB->Integral(1,-1) << " : " << h_PU200_prompt_not_vtx_matched_EB->Integral(1,2) << " : " << h_PU200_prompt_not_vtx_matched_EB->Integral(3,-1) << "  " << "(" << 100 -h_PU200_prompt_not_vtx_matched_EB->Integral(3,-1)/h_PU200_prompt_not_vtx_matched_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EB->Integral() << endl;
  cout << "[no tErr PV] " << h_PU200_prompt_no_tErr_PV_EB->GetBinContent(1) << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_PU200_prompt_EE_v2 = new TCanvas("c_st_PU200_prompt_EE_v2", "c_st_PU200_prompt_EE_v2", 1500, 1500);
  c_st_PU200_prompt_EE_v2->cd();
  c_st_PU200_prompt_EE_v2->SetLeftMargin(0.12);
  st_PU200_prompt_EE->Draw();
  leg_PU200_prompt_EE->Draw();
  c_st_PU200_prompt_EE_v2->Print("plots/ntuple/track_sigma_PU200_prompt_EE.pdf");
  cout << "PU200 Endcap prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EE->Integral(1,-1) << " : " << h_PU200_prompt_PV_EE->Integral(1,2) << " : " << h_PU200_prompt_PV_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EE->Integral(3,-1)/h_PU200_prompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EE->Integral(1,-1) << " : " << h_PU200_prompt_SV_EE->Integral(1,2) << " : " << h_PU200_prompt_SV_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EE->Integral(3,-1)/h_PU200_prompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EE->Integral(1,-1) << " : " << h_PU200_prompt_PU_EE->Integral(1,2) << " : " << h_PU200_prompt_PU_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EE->Integral(3,-1)/h_PU200_prompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EE->Integral(1,-1) << " : " << h_PU200_prompt_fake_EE->Integral(1,2) << " : " << h_PU200_prompt_fake_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EE->Integral(3,-1)/h_PU200_prompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_prompt_not_vtx_matched_EE->Integral(1,-1) << " : " << h_PU200_prompt_not_vtx_matched_EE->Integral(1,2) << " : " << h_PU200_prompt_not_vtx_matched_EE->Integral(3,-1) << "  " << "(" << 100 -h_PU200_prompt_not_vtx_matched_EE->Integral(3,-1)/h_PU200_prompt_not_vtx_matched_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EE->Integral() << endl;
  cout << endl;
    // nonprompt
      // Barrel
  TCanvas* c_st_PU200_nonprompt_EB_v2 = new TCanvas("c_st_PU200_nonprompt_EB_v2", "c_st_PU200_nonprompt_EB_v2", 1500, 1500);
  c_st_PU200_nonprompt_EB_v2->cd();
  c_st_PU200_nonprompt_EB_v2->SetLeftMargin(0.12);
  st_PU200_nonprompt_EB->Draw();
  leg_PU200_nonprompt_EB->Draw();
  c_st_PU200_nonprompt_EB_v2->Print("plots/ntuple/track_sigma_PU200_nonprompt_EB.pdf");
  cout << "PU200 Barrel nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EB->Integral(1,2) << " : " << h_PU200_nonprompt_PV_EB->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EB->Integral(3,-1)/h_PU200_nonprompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EB->Integral(1,2) << " : " << h_PU200_nonprompt_SV_EB->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EB->Integral(3,-1)/h_PU200_nonprompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EB->Integral(1,2) << " : " << h_PU200_nonprompt_PU_EB->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EB->Integral(3,-1)/h_PU200_nonprompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EB->Integral(1,2) << " : " << h_PU200_nonprompt_fake_EB->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EB->Integral(3,-1)/h_PU200_nonprompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_nonprompt_not_vtx_matched_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_not_vtx_matched_EB->Integral(1,2) << " : " << h_PU200_nonprompt_not_vtx_matched_EB->Integral(3,-1) << "  " << "(" << 100 -h_PU200_nonprompt_not_vtx_matched_EB->Integral(3,-1)/h_PU200_nonprompt_not_vtx_matched_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EB->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_PU200_nonprompt_EE_v2 = new TCanvas("c_st_PU200_nonprompt_EE_v2", "c_st_PU200_nonprompt_EE_v2", 1500, 1500);
  c_st_PU200_nonprompt_EE_v2->cd();
  c_st_PU200_nonprompt_EE_v2->SetLeftMargin(0.12);
  st_PU200_nonprompt_EE->Draw();
  leg_PU200_nonprompt_EE->Draw();
  c_st_PU200_nonprompt_EE_v2->Print("plots/ntuple/track_sigma_PU200_nonprompt_EE.pdf");
  cout << "PU200 Endcap nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EE->Integral(1,2) << " : " << h_PU200_nonprompt_PV_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EE->Integral(3,-1)/h_PU200_nonprompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EE->Integral(1,2) << " : " << h_PU200_nonprompt_SV_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EE->Integral(3,-1)/h_PU200_nonprompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EE->Integral(1,2) << " : " << h_PU200_nonprompt_PU_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EE->Integral(3,-1)/h_PU200_nonprompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EE->Integral(1,2) << " : " << h_PU200_nonprompt_fake_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EE->Integral(3,-1)/h_PU200_nonprompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_nonprompt_not_vtx_matched_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_not_vtx_matched_EE->Integral(1,2) << " : " << h_PU200_nonprompt_not_vtx_matched_EE->Integral(3,-1) << "  " << "(" << 100 -h_PU200_nonprompt_not_vtx_matched_EE->Integral(3,-1)/h_PU200_nonprompt_not_vtx_matched_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EE->Integral() << endl;
  cout << endl;
  
  // noPU
    // prompt
      // Barrel
  TCanvas* c_st_noPU_prompt_EB_v2 = new TCanvas("c_st_noPU_prompt_EB_v2", "c_st_noPU_prompt_EB_v2", 1500, 1500);
  c_st_noPU_prompt_EB_v2->cd();
  c_st_noPU_prompt_EB_v2->SetLeftMargin(0.12);
  st_noPU_prompt_EB->Draw();
  leg_noPU_prompt_EB->Draw();
  c_st_noPU_prompt_EB_v2->Print("plots/ntuple/track_sigma_noPU_prompt_EB.pdf");
  cout << "noPU Barrel prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EB->Integral(1,-1) << " : " << h_noPU_prompt_PV_EB->Integral(1,2) << " : " << h_noPU_prompt_PV_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EB->Integral(3,-1)/h_noPU_prompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EB->Integral(1,-1) << " : " << h_noPU_prompt_SV_EB->Integral(1,2) << " : " << h_noPU_prompt_SV_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EB->Integral(3,-1)/h_noPU_prompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EB->Integral(1,-1) << " : " << h_noPU_prompt_PU_EB->Integral(1,2) << " : " << h_noPU_prompt_PU_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EB->Integral(3,-1)/h_noPU_prompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EB->Integral(1,-1) << " : " << h_noPU_prompt_fake_EB->Integral(1,2) << " : " << h_noPU_prompt_fake_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EB->Integral(3,-1)/h_noPU_prompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_prompt_not_vtx_matched_EB->Integral(1,-1) << " : " << h_noPU_prompt_not_vtx_matched_EB->Integral(1,2) << " : " << h_noPU_prompt_not_vtx_matched_EB->Integral(3,-1) << "  " << "(" << 100 -h_noPU_prompt_not_vtx_matched_EB->Integral(3,-1)/h_noPU_prompt_not_vtx_matched_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EB->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_noPU_prompt_EE_v2 = new TCanvas("c_st_noPU_prompt_EE_v2", "c_st_noPU_prompt_EE_v2", 1500, 1500);
  c_st_noPU_prompt_EE_v2->cd();
  c_st_noPU_prompt_EE_v2->SetLeftMargin(0.12);
  st_noPU_prompt_EE->Draw();
  leg_noPU_prompt_EE->Draw();
  c_st_noPU_prompt_EE_v2->Print("plots/ntuple/track_sigma_noPU_prompt_EE.pdf");
  cout << "noPU Endcap prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EE->Integral(1,-1) << " : " << h_noPU_prompt_PV_EE->Integral(1,2) << " : " << h_noPU_prompt_PV_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EE->Integral(3,-1)/h_noPU_prompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EE->Integral(1,-1) << " : " << h_noPU_prompt_SV_EE->Integral(1,2) << " : " << h_noPU_prompt_SV_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EE->Integral(3,-1)/h_noPU_prompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EE->Integral(1,-1) << " : " << h_noPU_prompt_PU_EE->Integral(1,2) << " : " << h_noPU_prompt_PU_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EE->Integral(3,-1)/h_noPU_prompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EE->Integral(1,-1) << " : " << h_noPU_prompt_fake_EE->Integral(1,2) << " : " << h_noPU_prompt_fake_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EE->Integral(3,-1)/h_noPU_prompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_prompt_not_vtx_matched_EE->Integral(1,-1) << " : " << h_noPU_prompt_not_vtx_matched_EE->Integral(1,2) << " : " << h_noPU_prompt_not_vtx_matched_EE->Integral(3,-1) << "  " << "(" << 100 -h_noPU_prompt_not_vtx_matched_EE->Integral(3,-1)/h_noPU_prompt_not_vtx_matched_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EE->Integral() << endl;
  cout << endl;
    // nonprompt
      // Barrel
  TCanvas* c_st_noPU_nonprompt_EB_v2 = new TCanvas("c_st_noPU_nonprompt_EB_v2", "c_st_noPU_nonprompt_EB_v2", 1500, 1500);
  c_st_noPU_nonprompt_EB_v2->cd();
  c_st_noPU_nonprompt_EB_v2->SetLeftMargin(0.12);
  st_noPU_nonprompt_EB->Draw();
  leg_noPU_nonprompt_EB->Draw();
  c_st_noPU_nonprompt_EB_v2->Print("plots/ntuple/track_sigma_noPU_nonprompt_EB.pdf");
  cout << "noPU Barrel nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EB->Integral(1,2) << " : " << h_noPU_nonprompt_PV_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EB->Integral(3,-1)/h_noPU_nonprompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EB->Integral(1,2) << " : " << h_noPU_nonprompt_SV_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EB->Integral(3,-1)/h_noPU_nonprompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EB->Integral(1,2) << " : " << h_noPU_nonprompt_PU_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EB->Integral(3,-1)/h_noPU_nonprompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EB->Integral(1,2) << " : " << h_noPU_nonprompt_fake_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EB->Integral(3,-1)/h_noPU_nonprompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_nonprompt_not_vtx_matched_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_not_vtx_matched_EB->Integral(1,2) << " : " << h_noPU_nonprompt_not_vtx_matched_EB->Integral(3,-1) << "  " << "(" << 100 -h_noPU_nonprompt_not_vtx_matched_EB->Integral(3,-1)/h_noPU_nonprompt_not_vtx_matched_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EB->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_noPU_nonprompt_EE_v2 = new TCanvas("c_st_noPU_nonprompt_EE_v2", "c_st_noPU_nonprompt_EE_v2", 1500, 1500);
  c_st_noPU_nonprompt_EE_v2->cd();
  c_st_noPU_nonprompt_EE_v2->SetLeftMargin(0.12);
  st_noPU_nonprompt_EE->Draw();
  leg_noPU_nonprompt_EE->Draw();
  c_st_noPU_nonprompt_EE_v2->Print("plots/ntuple/track_sigma_noPU_nonprompt_EE.pdf");
  cout << "noPU Endcap nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EE->Integral(1,2) << " : " << h_noPU_nonprompt_PV_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EE->Integral(3,-1)/h_noPU_nonprompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EE->Integral(1,2) << " : " << h_noPU_nonprompt_SV_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EE->Integral(3,-1)/h_noPU_nonprompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EE->Integral(1,2) << " : " << h_noPU_nonprompt_PU_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EE->Integral(3,-1)/h_noPU_nonprompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EE->Integral(1,2) << " : " << h_noPU_nonprompt_fake_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EE->Integral(3,-1)/h_noPU_nonprompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_nonprompt_not_vtx_matched_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_not_vtx_matched_EE->Integral(1,2) << " : " << h_noPU_nonprompt_not_vtx_matched_EE->Integral(3,-1) << "  " << "(" << 100 -h_noPU_nonprompt_not_vtx_matched_EE->Integral(3,-1)/h_noPU_nonprompt_not_vtx_matched_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EE->Integral() << endl;
  cout << endl;

  cout << "///////////////////////" << endl;
  cout << "/////    3sigma   /////" << endl;
  cout << "///////////////////////" << endl;
  cout << "PU200 Barrel prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EB->Integral(1,-1) << " : " << h_PU200_prompt_PV_EB->Integral(1,3) << " : " << h_PU200_prompt_PV_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EB->Integral(4,-1)/h_PU200_prompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EB->Integral(1,-1) << " : " << h_PU200_prompt_SV_EB->Integral(1,3) << " : " << h_PU200_prompt_SV_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EB->Integral(4,-1)/h_PU200_prompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EB->Integral(1,-1) << " : " << h_PU200_prompt_PU_EB->Integral(1,3) << " : " << h_PU200_prompt_PU_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EB->Integral(4,-1)/h_PU200_prompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EB->Integral(1,-1) << " : " << h_PU200_prompt_fake_EB->Integral(1,3) << " : " << h_PU200_prompt_fake_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EB->Integral(4,-1)/h_PU200_prompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_prompt_not_vtx_matched_EB->Integral(1,-1) << " : " << h_PU200_prompt_not_vtx_matched_EB->Integral(1,3) << " : " << h_PU200_prompt_not_vtx_matched_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_not_vtx_matched_EB->Integral(4,-1)/h_PU200_prompt_not_vtx_matched_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EB->Integral() << endl;
  cout << endl;
      // Endcap
  cout << "PU200 Endcap prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EE->Integral(1,-1) << " : " << h_PU200_prompt_PV_EE->Integral(1,3) << " : " << h_PU200_prompt_PV_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EE->Integral(4,-1)/h_PU200_prompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EE->Integral(1,-1) << " : " << h_PU200_prompt_PU_EE->Integral(1,3) << " : " << h_PU200_prompt_PU_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EE->Integral(4,-1)/h_PU200_prompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EE->Integral(1,-1) << " : " << h_PU200_prompt_fake_EE->Integral(1,3) << " : " << h_PU200_prompt_fake_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EE->Integral(4,-1)/h_PU200_prompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EE->Integral(1,-1) << " : " << h_PU200_prompt_SV_EE->Integral(1,3) << " : " << h_PU200_prompt_SV_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EE->Integral(4,-1)/h_PU200_prompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_prompt_not_vtx_matched_EE->Integral(1,-1) << " : " << h_PU200_prompt_not_vtx_matched_EE->Integral(1,3) << " : " << h_PU200_prompt_not_vtx_matched_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_not_vtx_matched_EE->Integral(4,-1)/h_PU200_prompt_not_vtx_matched_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EE->Integral() << endl;
  cout << endl;
    // nonprompt
      // Barrel
  cout << "PU200 Barrel nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EB->Integral(1,3) << " : " << h_PU200_nonprompt_PV_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EB->Integral(4,-1)/h_PU200_nonprompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EB->Integral(1,3) << " : " << h_PU200_nonprompt_SV_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EB->Integral(4,-1)/h_PU200_nonprompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EB->Integral(1,3) << " : " << h_PU200_nonprompt_PU_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EB->Integral(4,-1)/h_PU200_nonprompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EB->Integral(1,3) << " : " << h_PU200_nonprompt_fake_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EB->Integral(4,-1)/h_PU200_nonprompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_nonprompt_not_vtx_matched_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_not_vtx_matched_EB->Integral(1,3) << " : " << h_PU200_nonprompt_not_vtx_matched_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_not_vtx_matched_EB->Integral(4,-1)/h_PU200_nonprompt_not_vtx_matched_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EB->Integral() << endl;
  cout << endl;
      // Endcap
  cout << "PU200 Endcap nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EE->Integral(1,3) << " : " << h_PU200_nonprompt_PV_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EE->Integral(4,-1)/h_PU200_nonprompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EE->Integral(1,3) << " : " << h_PU200_nonprompt_SV_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EE->Integral(4,-1)/h_PU200_nonprompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EE->Integral(1,3) << " : " << h_PU200_nonprompt_PU_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EE->Integral(4,-1)/h_PU200_nonprompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EE->Integral(1,3) << " : " << h_PU200_nonprompt_fake_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EE->Integral(4,-1)/h_PU200_nonprompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_nonprompt_not_vtx_matched_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_not_vtx_matched_EE->Integral(1,3) << " : " << h_PU200_nonprompt_not_vtx_matched_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_not_vtx_matched_EE->Integral(4,-1)/h_PU200_nonprompt_not_vtx_matched_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EE->Integral() << endl;
  cout << endl;

  // noPU
    // prompt
      // Barrel
  cout << "noPU Barrel prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EB->Integral(1,-1) << " : " << h_noPU_prompt_PV_EB->Integral(1,3) << " : " << h_noPU_prompt_PV_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EB->Integral(4,-1)/h_noPU_prompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EB->Integral(1,-1) << " : " << h_noPU_prompt_SV_EB->Integral(1,3) << " : " << h_noPU_prompt_SV_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EB->Integral(4,-1)/h_noPU_prompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EB->Integral(1,-1) << " : " << h_noPU_prompt_PU_EB->Integral(1,3) << " : " << h_noPU_prompt_PU_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EB->Integral(4,-1)/h_noPU_prompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EB->Integral(1,-1) << " : " << h_noPU_prompt_fake_EB->Integral(1,3) << " : " << h_noPU_prompt_fake_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EB->Integral(4,-1)/h_noPU_prompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_prompt_not_vtx_matched_EB->Integral(1,-1) << " : " << h_noPU_prompt_not_vtx_matched_EB->Integral(1,3) << " : " << h_noPU_prompt_not_vtx_matched_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_not_vtx_matched_EB->Integral(4,-1)/h_noPU_prompt_not_vtx_matched_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EB->Integral() << endl;
  cout << endl;
      // Endcap
  cout << "noPU Endcap prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EE->Integral(1,-1) << " : " << h_noPU_prompt_PV_EE->Integral(1,3) << " : " << h_noPU_prompt_PV_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EE->Integral(4,-1)/h_noPU_prompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EE->Integral(1,-1) << " : " << h_noPU_prompt_SV_EE->Integral(1,3) << " : " << h_noPU_prompt_SV_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EE->Integral(4,-1)/h_noPU_prompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EE->Integral(1,-1) << " : " << h_noPU_prompt_PU_EE->Integral(1,3) << " : " << h_noPU_prompt_PU_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EE->Integral(4,-1)/h_noPU_prompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EE->Integral(1,-1) << " : " << h_noPU_prompt_fake_EE->Integral(1,3) << " : " << h_noPU_prompt_fake_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EE->Integral(4,-1)/h_noPU_prompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_prompt_not_vtx_matched_EE->Integral(1,-1) << " : " << h_noPU_prompt_not_vtx_matched_EE->Integral(1,3) << " : " << h_noPU_prompt_not_vtx_matched_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_not_vtx_matched_EE->Integral(4,-1)/h_noPU_prompt_not_vtx_matched_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EE->Integral() << endl;
  cout << endl;
    // nonprompt
      // Barrel
  cout << "noPU Barrel nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EB->Integral(1,3) << " : " << h_noPU_nonprompt_PV_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EB->Integral(4,-1)/h_noPU_nonprompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EB->Integral(1,3) << " : " << h_noPU_nonprompt_SV_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EB->Integral(4,-1)/h_noPU_nonprompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EB->Integral(1,3) << " : " << h_noPU_nonprompt_PU_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EB->Integral(4,-1)/h_noPU_nonprompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EB->Integral(1,3) << " : " << h_noPU_nonprompt_fake_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EB->Integral(4,-1)/h_noPU_nonprompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_nonprompt_not_vtx_matched_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_not_vtx_matched_EB->Integral(1,3) << " : " << h_noPU_nonprompt_not_vtx_matched_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_not_vtx_matched_EB->Integral(4,-1)/h_noPU_nonprompt_not_vtx_matched_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EB->Integral() << endl;
  cout << endl;
      // Endcap
  cout << "noPU Endcap nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EE->Integral(1,3) << " : " << h_noPU_nonprompt_PV_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EE->Integral(4,-1)/h_noPU_nonprompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EE->Integral(1,3) << " : " << h_noPU_nonprompt_SV_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EE->Integral(4,-1)/h_noPU_nonprompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EE->Integral(1,3) << " : " << h_noPU_nonprompt_PU_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EE->Integral(4,-1)/h_noPU_nonprompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EE->Integral(1,3) << " : " << h_noPU_nonprompt_fake_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EE->Integral(4,-1)/h_noPU_nonprompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_nonprompt_not_vtx_matched_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_not_vtx_matched_EE->Integral(1,3) << " : " << h_noPU_nonprompt_not_vtx_matched_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_not_vtx_matched_EE->Integral(4,-1)/h_noPU_nonprompt_not_vtx_matched_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EE->Integral() << endl;
  cout << endl;




  /////////////////////
  /////// PVPVPV //////
  /////////////////////


  // prompt
  // Barrel
  TLegend *leg_PV_dtsig_PU200_prompt_EB_F = new TLegend(0.65, 0.75, 0.88, 0.85);
  leg_PV_dtsig_PU200_prompt_EB_F->AddEntry(h_PU200_prompt_PV_EB, "Track from PV", "F");

  TLegend *leg_PV_dtsig_PU200_prompt_EB_L = new TLegend(0.65, 0.75, 0.88, 0.85);
  leg_PV_dtsig_PU200_prompt_EB_L->AddEntry(h_PU200_prompt_PV_EB, "Track from PV", "L");
  TLegend *leg_simvtx_PV_dtsig_PU200_prompt_EB_L = new TLegend(0.65, 0.75, 0.88, 0.85);
  leg_simvtx_PV_dtsig_PU200_prompt_EB_L->AddEntry(h_PU200_prompt_PV_EB, "sim vertex", "L");
  TLegend *leg_PV_dtsig_PU200_prompt_EB_L_vtx = new TLegend(0.65, 0.75, 0.88, 0.85);
  leg_PV_dtsig_PU200_prompt_EB_L_vtx->AddEntry(h_PU200_prompt_PV_EB_vtx, "Track from PV", "L");

  TCanvas* c_PV_dtsig_PU200_prompt_EB = new TCanvas("c_PV_dtsig_PU200_prompt_EB", "c_PV_dtsig_PU200_prompt_EB", 1500, 1500);
  c_PV_dtsig_PU200_prompt_EB->cd();
  c_PV_dtsig_PU200_prompt_EB->SetLeftMargin(0.12);
  h_PU200_prompt_PV_EB->SetTitle(";#sigma_{t} (muon, track);Counts");
  h_PU200_prompt_PV_EB->GetYaxis()->SetRangeUser(0.,400.);
  h_PU200_prompt_PV_EB->Draw();
  leg_PV_dtsig_PU200_prompt_EB_F->Draw();
  c_PV_dtsig_PU200_prompt_EB->Print("plots/ntuple/PV_dtsig_PU200_prompt_EB.pdf");

  TCanvas* c_PV_bx_PU200_prompt_EB = new TCanvas("c_PV_bx_PU200_prompt_EB", "c_PV_bx_PU200_prompt_EB", 1500, 1500);
  c_PV_bx_PU200_prompt_EB->cd();
  c_PV_bx_PU200_prompt_EB->SetLeftMargin(0.12);
  h_PU200_prompt_bx_PV_EB->SetTitle(";Bunch crossing number of track from PV;Counts");
  h_PU200_prompt_bx_PV_EB->SetLineColor(kBlack);
  h_PU200_prompt_bx_PV_EB->SetLineWidth(2);
  h_PU200_prompt_bx_PV_EB->Draw();
  leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_bx_PU200_prompt_EB->Print("plots/ntuple/PV_bx_PU200_prompt_EB.pdf");
  cout << "bx = 0  : " << h_PU200_prompt_bx_PV_EB->GetBinContent(1) << endl;
  cout << "total bx: " << h_PU200_prompt_bx_PV_EB->Integral(0,-1) << endl;
  cout << endl;

  TCanvas* c_PV_evtId_PU200_prompt_EB = new TCanvas("c_PV_evtId_PU200_prompt_EB", "c_PV_evtId_PU200_prompt_EB", 1500, 1500);
  c_PV_evtId_PU200_prompt_EB->cd();
  c_PV_evtId_PU200_prompt_EB->SetLeftMargin(0.12);
  h_PU200_prompt_evtId_PV_EB->SetTitle(";Event ID of track from PV;Counts");
  h_PU200_prompt_evtId_PV_EB->SetLineColor(kBlack);
  h_PU200_prompt_evtId_PV_EB->SetLineWidth(2);
  h_PU200_prompt_evtId_PV_EB->Draw();
  c_PV_evtId_PU200_prompt_EB->SetLogy();
  leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_evtId_PU200_prompt_EB->Print("plots/ntuple/PV_evtId_PU200_prompt_EB.pdf");
  cout << "evtId = 0  : " << h_PU200_prompt_evtId_PV_EB->GetBinContent(1) << endl;
  cout << "total evtId: " << h_PU200_prompt_evtId_PV_EB->Integral(0,-1) << endl;
  cout << endl;

  TCanvas* c_PV_PVweight_PU200_prompt_EB = new TCanvas("c_PV_PVweight_PU200_prompt_EB", "c_PV_PVweight_PU200_prompt_EB", 1500, 1500);
  c_PV_PVweight_PU200_prompt_EB->cd();
  c_PV_PVweight_PU200_prompt_EB->SetLeftMargin(0.12);
  h_PU200_prompt_PVweight_PV_EB->SetTitle(";PV weight of track from PV;Counts");
  h_PU200_prompt_PVweight_PV_EB->SetLineColor(kBlack);
  h_PU200_prompt_PVweight_PV_EB->SetLineWidth(2);
  h_PU200_prompt_PVweight_PV_EB->Draw();
  leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_PVweight_PU200_prompt_EB->Print("plots/ntuple/PV_PVweight_PU200_prompt_EB.pdf");
  cout << "PV weight = 0  : " << h_PU200_prompt_PVweight_PV_EB->GetBinContent(1) << endl;
  cout << "less than 0.5  : " << h_PU200_prompt_PVweight_PV_EB->Integral(0,50) << endl;
  cout << "total PV weight: " << h_PU200_prompt_PVweight_PV_EB->Integral(0,-1) << endl;
  cout << endl;

  TCanvas* c_PV_simvtx_bx_PU200_prompt_EB = new TCanvas("c_PV_simvtx_bx_PU200_prompt_EB", "c_PV_simvtx_bx_PU200_prompt_EB", 1500, 1500);
  c_PV_simvtx_bx_PU200_prompt_EB->cd();
  c_PV_simvtx_bx_PU200_prompt_EB->SetLeftMargin(0.12);
  h_PU200_prompt_simvtx_bx_PV_EB->SetTitle(";Bunch crossing number of sim vertex;Counts");
  h_PU200_prompt_simvtx_bx_PV_EB->SetLineColor(kBlack);
  h_PU200_prompt_simvtx_bx_PV_EB->SetLineWidth(2);
  h_PU200_prompt_simvtx_bx_PV_EB->Draw();
  leg_simvtx_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_simvtx_bx_PU200_prompt_EB->Print("plots/ntuple/PV_simvtx_bx_PU200_prompt_EB.pdf");
  cout << "simvtx bx = 0  : " << h_PU200_prompt_simvtx_bx_PV_EB->GetBinContent(1) << endl;
  cout << "total simvtx bx: " << h_PU200_prompt_simvtx_bx_PV_EB->Integral(0,-1) << endl;
  cout << endl;


  TCanvas* c_PV_simvtx_evtId_PU200_prompt_EB = new TCanvas("c_PV_simvtx_evtId_PU200_prompt_EB", "c_PV_simvtx_evtId_PU200_prompt_EB", 1500, 1500);
  c_PV_simvtx_evtId_PU200_prompt_EB->cd();
  c_PV_simvtx_evtId_PU200_prompt_EB->SetLeftMargin(0.12);
  h_PU200_prompt_simvtx_evtId_PV_EB->SetTitle(";Event ID of sim vertex;Counts");
  h_PU200_prompt_simvtx_evtId_PV_EB->SetLineColor(kBlack);
  h_PU200_prompt_simvtx_evtId_PV_EB->SetLineWidth(2);
  h_PU200_prompt_simvtx_evtId_PV_EB->Draw();
  c_PV_simvtx_evtId_PU200_prompt_EB->SetLogy();
  leg_simvtx_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_simvtx_evtId_PU200_prompt_EB->Print("plots/ntuple/PV_simvtx_evtId_PU200_prompt_EB.pdf");
  cout << "simvtx evtId = 0  : " << h_PU200_prompt_simvtx_evtId_PV_EB->GetBinContent(1) << endl;
  cout << "total simvtx evtId: " << h_PU200_prompt_simvtx_evtId_PV_EB->Integral(0,-1) << endl;
  cout << endl;


  // nonprompt
  // Barrel
  TCanvas* c_PV_dtsig_PU200_nonprompt_EB = new TCanvas("c_PV_dtsig_PU200_nonprompt_EB", "c_PV_dtsig_PU200_nonprompt_EB", 1500, 1500);
  c_PV_dtsig_PU200_nonprompt_EB->cd();
  c_PV_dtsig_PU200_nonprompt_EB->SetLeftMargin(0.12);
  h_PU200_nonprompt_PV_EB->SetTitle(";#sigma_{t} (muon, track);Counts");
  h_PU200_nonprompt_PV_EB->GetYaxis()->SetRangeUser(0.,700.);
  h_PU200_nonprompt_PV_EB->Draw();
  leg_PV_dtsig_PU200_prompt_EB_F->Draw();
  c_PV_dtsig_PU200_nonprompt_EB->Print("plots/ntuple/PV_dtsig_PU200_nonprompt_EB.pdf");

  TCanvas* c_PV_bx_PU200_nonprompt_EB = new TCanvas("c_PV_bx_PU200_nonprompt_EB", "c_PV_bx_PU200_nonprompt_EB", 1500, 1500);
  c_PV_bx_PU200_nonprompt_EB->cd();
  c_PV_bx_PU200_nonprompt_EB->SetLeftMargin(0.12);
  h_PU200_nonprompt_bx_PV_EB->SetTitle(";Bunch crossing number of track from PV;Counts");
  h_PU200_nonprompt_bx_PV_EB->SetLineColor(kBlack);
  h_PU200_nonprompt_bx_PV_EB->SetLineWidth(2);
  h_PU200_nonprompt_bx_PV_EB->Draw();
  leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_bx_PU200_nonprompt_EB->Print("plots/ntuple/PV_bx_PU200_nonprompt_EB.pdf");
  cout << "bx = 0  : " << h_PU200_nonprompt_bx_PV_EB->GetBinContent(1) << endl;
  cout << "total bx: " << h_PU200_nonprompt_bx_PV_EB->Integral(0,-1) << endl;
  cout << endl;

  TCanvas* c_PV_evtId_PU200_nonprompt_EB = new TCanvas("c_PV_evtId_PU200_nonprompt_EB", "c_PV_evtId_PU200_nonprompt_EB", 1500, 1500);
  c_PV_evtId_PU200_nonprompt_EB->cd();
  c_PV_evtId_PU200_nonprompt_EB->SetLeftMargin(0.12);
  h_PU200_nonprompt_evtId_PV_EB->SetTitle(";Event ID of track from PV;Counts");
  h_PU200_nonprompt_evtId_PV_EB->SetLineColor(kBlack);
  h_PU200_nonprompt_evtId_PV_EB->SetLineWidth(2);
  h_PU200_nonprompt_evtId_PV_EB->Draw();
  c_PV_evtId_PU200_nonprompt_EB->SetLogy();
  leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_evtId_PU200_nonprompt_EB->Print("plots/ntuple/PV_evtId_PU200_nonprompt_EB.pdf");
  cout << "evtId = 0  : " << h_PU200_nonprompt_evtId_PV_EB->GetBinContent(1) << endl;
  cout << "total evtId: " << h_PU200_nonprompt_evtId_PV_EB->Integral(0,-1) << endl;
  cout << endl;

  TCanvas* c_PV_PVweight_PU200_nonprompt_EB = new TCanvas("c_PV_PVweight_PU200_nonprompt_EB", "c_PV_PVweight_PU200_nonprompt_EB", 1500, 1500);
  c_PV_PVweight_PU200_nonprompt_EB->cd();
  c_PV_PVweight_PU200_nonprompt_EB->SetLeftMargin(0.12);
  h_PU200_nonprompt_PVweight_PV_EB->SetTitle(";PV weight of track from PV;Counts");
  h_PU200_nonprompt_PVweight_PV_EB->SetLineColor(kBlack);
  h_PU200_nonprompt_PVweight_PV_EB->SetLineWidth(2);
  h_PU200_nonprompt_PVweight_PV_EB->Draw();
  leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_PVweight_PU200_nonprompt_EB->Print("plots/ntuple/PV_PVweight_PU200_nonprompt_EB.pdf");
  cout << "PV weight = 0  : " << h_PU200_nonprompt_PVweight_PV_EB->GetBinContent(1) << endl;
  cout << "less than 0.5  : " << h_PU200_nonprompt_PVweight_PV_EB->Integral(0,50) << endl;
  cout << "total PV weight: " << h_PU200_nonprompt_PVweight_PV_EB->Integral(0,-1) << endl;
  cout << endl;

  TCanvas* c_PV_simvtx_bx_PU200_nonprompt_EB = new TCanvas("c_PV_simvtx_bx_PU200_nonprompt_EB", "c_PV_simvtx_bx_PU200_nonprompt_EB", 1500, 1500);
  c_PV_simvtx_bx_PU200_nonprompt_EB->cd();
  c_PV_simvtx_bx_PU200_nonprompt_EB->SetLeftMargin(0.12);
  h_PU200_nonprompt_simvtx_bx_PV_EB->SetTitle(";Bunch crossing number of sim vertex;Counts");
  h_PU200_nonprompt_simvtx_bx_PV_EB->SetLineColor(kBlack);
  h_PU200_nonprompt_simvtx_bx_PV_EB->SetLineWidth(2);
  h_PU200_nonprompt_simvtx_bx_PV_EB->Draw();
  leg_simvtx_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_simvtx_bx_PU200_nonprompt_EB->Print("plots/ntuple/PV_simvtx_bx_PU200_nonprompt_EB.pdf");
  cout << "simvtx bx = 0  : " << h_PU200_nonprompt_simvtx_bx_PV_EB->GetBinContent(1) << endl;
  cout << "total simvtx bx: " << h_PU200_nonprompt_simvtx_bx_PV_EB->Integral(0,-1) << endl;
  cout << endl;


  TCanvas* c_PV_simvtx_evtId_PU200_nonprompt_EB = new TCanvas("c_PV_simvtx_evtId_PU200_nonprompt_EB", "c_PV_simvtx_evtId_PU200_nonprompt_EB", 1500, 1500);
  c_PV_simvtx_evtId_PU200_nonprompt_EB->cd();
  c_PV_simvtx_evtId_PU200_nonprompt_EB->SetLeftMargin(0.12);
  h_PU200_nonprompt_simvtx_evtId_PV_EB->SetTitle(";Event ID of sim vertex;Counts");
  h_PU200_nonprompt_simvtx_evtId_PV_EB->SetLineColor(kBlack);
  h_PU200_nonprompt_simvtx_evtId_PV_EB->SetLineWidth(2);
  h_PU200_nonprompt_simvtx_evtId_PV_EB->Draw();
  c_PV_simvtx_evtId_PU200_nonprompt_EB->SetLogy();
  leg_simvtx_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_simvtx_evtId_PU200_nonprompt_EB->Print("plots/ntuple/PV_simvtx_evtId_PU200_nonprompt_EB.pdf");
  cout << "simvtx evtId = 0  : " << h_PU200_nonprompt_simvtx_evtId_PV_EB->GetBinContent(1) << endl;
  cout << "total simvtx evtId: " << h_PU200_nonprompt_simvtx_evtId_PV_EB->Integral(0,-1) << endl;
  cout << endl;






  // test


  // 2D
  cout << endl;
  cout << "//////////////" << endl;
  cout << "//    2D    //" << endl;
  cout << "//////////////" << endl;
  TCanvas* c_PV_dtsig_PU200_prompt_EB_2D = new TCanvas("c_PV_dtsig_PU200_prompt_EB_2D", "c_PV_dtsig_PU200_prompt_EB_2D", 1500, 1500);
  c_PV_dtsig_PU200_prompt_EB_2D->cd();
  c_PV_dtsig_PU200_prompt_EB_2D->SetRightMargin(0.12);
  h_PU200_prompt_PV_EB_2D->SetTitle("");
  h_PU200_prompt_PV_EB_2D->Draw("colz text");
  h_PU200_prompt_PV_EB_2D->GetXaxis()->SetTitle("#sigma_{t} (muon, track)");
  h_PU200_prompt_PV_EB_2D->GetYaxis()->SetTitle("#sigma_{t} (PV, track)");
  c_PV_dtsig_PU200_prompt_EB_2D->Print("plots/ntuple/PV_dtsig_PU200_prompt_EB_2D.pdf");
  cout << "PU200 prompt Barrel" << endl;
  cout << h_PU200_prompt_PV_EB_2D->Integral(0,2,3,-1) << " / " << h_PU200_prompt_PV_EB_2D->Integral(0,2,0,-1) << endl;
  cout << h_PU200_prompt_PV_EB_2D->Integral(3,-1,0,2) << " / " << h_PU200_prompt_PV_EB_2D->Integral(0,-1,0,2) << endl;

  TCanvas* c_PV_dtsig_PU200_prompt_EE_2D = new TCanvas("c_PV_dtsig_PU200_prompt_EE_2D", "c_PV_dtsig_PU200_prompt_EE_2D", 1500, 1500);
  c_PV_dtsig_PU200_prompt_EE_2D->cd();
  c_PV_dtsig_PU200_prompt_EE_2D->SetRightMargin(0.12);
  h_PU200_prompt_PV_EE_2D->SetTitle("");
  h_PU200_prompt_PV_EE_2D->Draw("colz");
  h_PU200_prompt_PV_EE_2D->GetXaxis()->SetTitle("#sigma_{t} (muon, track)");
  h_PU200_prompt_PV_EE_2D->GetYaxis()->SetTitle("#sigma_{t} (PV, track)");
  c_PV_dtsig_PU200_prompt_EE_2D->Print("plots/ntuple/PV_dtsig_PU200_prompt_EE_2D.pdf");
  
  TCanvas* c_PV_dtsig_PU200_nonprompt_EB_2D = new TCanvas("c_PV_dtsig_PU200_nonprompt_EB_2D", "c_PV_dtsig_PU200_nonprompt_EB_2D", 1500, 1500);
  c_PV_dtsig_PU200_nonprompt_EB_2D->cd();
  c_PV_dtsig_PU200_nonprompt_EB_2D->SetRightMargin(0.12);
  h_PU200_nonprompt_PV_EB_2D->SetTitle("");
  h_PU200_nonprompt_PV_EB_2D->Draw("colz");
  h_PU200_nonprompt_PV_EB_2D->GetXaxis()->SetTitle("#sigma_{t} (muon, track)");
  h_PU200_nonprompt_PV_EB_2D->GetYaxis()->SetTitle("#sigma_{t} (PV, track)");
  c_PV_dtsig_PU200_nonprompt_EB_2D->Print("plots/ntuple/PV_dtsig_PU200_nonprompt_EB_2D.pdf");
  cout << "PU200 nonprompt Barrel" << endl;
  cout << h_PU200_nonprompt_PV_EB_2D->Integral(0,2,3,-1) << " / " << h_PU200_nonprompt_PV_EB_2D->Integral(0,2,0,-1) << endl;
  cout << h_PU200_nonprompt_PV_EB_2D->Integral(3,-1,0,2) << " / " << h_PU200_nonprompt_PV_EB_2D->Integral(0,-1,0,2) << endl;

  TCanvas* c_PV_dtsig_PU200_nonprompt_EE_2D = new TCanvas("c_PV_dtsig_PU200_nonprompt_EE_2D", "c_PV_dtsig_PU200_nonprompt_EE_2D", 1500, 1500);
  c_PV_dtsig_PU200_nonprompt_EE_2D->cd();
  c_PV_dtsig_PU200_nonprompt_EE_2D->SetRightMargin(0.12);
  h_PU200_nonprompt_PV_EE_2D->SetTitle("");
  h_PU200_nonprompt_PV_EE_2D->Draw("colz");
  h_PU200_nonprompt_PV_EE_2D->GetXaxis()->SetTitle("#sigma_{t} (muon, track)");
  h_PU200_nonprompt_PV_EE_2D->GetYaxis()->SetTitle("#sigma_{t} (PV, track)");
  c_PV_dtsig_PU200_nonprompt_EE_2D->Print("plots/ntuple/PV_dtsig_PU200_nonprompt_EE_2D.pdf");

  TCanvas* c_PV_dtsig_noPU_prompt_EB_2D = new TCanvas("c_PV_dtsig_noPU_prompt_EB_2D", "c_PV_dtsig_noPU_prompt_EB_2D", 1500, 1500);
  c_PV_dtsig_noPU_prompt_EB_2D->cd();
  c_PV_dtsig_noPU_prompt_EB_2D->SetRightMargin(0.12);
  h_noPU_prompt_PV_EB_2D->SetTitle("");
  h_noPU_prompt_PV_EB_2D->Draw("colz");
  h_noPU_prompt_PV_EB_2D->GetXaxis()->SetTitle("#sigma_{t} (muon, track)");
  h_noPU_prompt_PV_EB_2D->GetYaxis()->SetTitle("#sigma_{t} (PV, track)");
  c_PV_dtsig_noPU_prompt_EB_2D->Print("plots/ntuple/PV_dtsig_noPU_prompt_EB_2D.pdf");
  cout << "noPU prompt Barrel" << endl;
  cout << h_noPU_prompt_PV_EB_2D->Integral(0,2,3,-1) << " / " << h_noPU_prompt_PV_EB_2D->Integral(0,2,0,-1) << endl;
  cout << h_noPU_prompt_PV_EB_2D->Integral(3,-1,0,2) << " / " << h_noPU_prompt_PV_EB_2D->Integral(0,-1,0,2) << endl;

  TCanvas* c_PV_dtsig_noPU_prompt_EE_2D = new TCanvas("c_PV_dtsig_noPU_prompt_EE_2D", "c_PV_dtsig_noPU_prompt_EE_2D", 1500, 1500);
  c_PV_dtsig_noPU_prompt_EE_2D->cd();
  c_PV_dtsig_noPU_prompt_EE_2D->SetRightMargin(0.12);
  h_noPU_prompt_PV_EE_2D->SetTitle("");
  h_noPU_prompt_PV_EE_2D->Draw("colz");
  h_noPU_prompt_PV_EE_2D->GetXaxis()->SetTitle("#sigma_{t} (muon, track)");
  h_noPU_prompt_PV_EE_2D->GetYaxis()->SetTitle("#sigma_{t} (PV, track)");
  c_PV_dtsig_noPU_prompt_EE_2D->Print("plots/ntuple/PV_dtsig_noPU_prompt_EE_2D.pdf");
  
  TCanvas* c_PV_dtsig_noPU_nonprompt_EB_2D = new TCanvas("c_PV_dtsig_noPU_nonprompt_EB_2D", "c_PV_dtsig_noPU_nonprompt_EB_2D", 1500, 1500);
  c_PV_dtsig_noPU_nonprompt_EB_2D->cd();
  c_PV_dtsig_noPU_nonprompt_EB_2D->SetRightMargin(0.12);
  h_noPU_nonprompt_PV_EB_2D->SetTitle("");
  h_noPU_nonprompt_PV_EB_2D->Draw("colz");
  h_noPU_nonprompt_PV_EB_2D->GetXaxis()->SetTitle("#sigma_{t} (muon, track)");
  h_noPU_nonprompt_PV_EB_2D->GetYaxis()->SetTitle("#sigma_{t} (PV, track)");
  c_PV_dtsig_noPU_nonprompt_EB_2D->Print("plots/ntuple/PV_dtsig_noPU_nonprompt_EB_2D.pdf");
  cout << "noPU nonprompt Barrel" << endl;
  cout << h_noPU_nonprompt_PV_EB_2D->Integral(0,2,3,-1) << " / " << h_noPU_nonprompt_PV_EB_2D->Integral(0,2,0,-1) << endl;
  cout << h_noPU_nonprompt_PV_EB_2D->Integral(3,-1,0,2) << " / " << h_noPU_nonprompt_PV_EB_2D->Integral(0,-1,0,2) << endl;

  TCanvas* c_PV_dtsig_noPU_nonprompt_EE_2D = new TCanvas("c_PV_dtsig_noPU_nonprompt_EE_2D", "c_PV_dtsig_noPU_nonprompt_EE_2D", 1500, 1500);
  c_PV_dtsig_noPU_nonprompt_EE_2D->cd();
  c_PV_dtsig_noPU_nonprompt_EE_2D->SetRightMargin(0.12);
  h_noPU_nonprompt_PV_EE_2D->SetTitle("");
  h_noPU_nonprompt_PV_EE_2D->Draw("colz");
  h_noPU_nonprompt_PV_EE_2D->GetXaxis()->SetTitle("#sigma_{t} (muon, track)");
  h_noPU_nonprompt_PV_EE_2D->GetYaxis()->SetTitle("#sigma_{t} (PV, track)");
  c_PV_dtsig_noPU_nonprompt_EE_2D->Print("plots/ntuple/PV_dtsig_noPU_nonprompt_EE_2D.pdf");

  // 2D end
  cout << endl;





  // tsim
  // PU200
  // Barrel
    // track
  TCanvas* c_PV_track_time_PU200_prompt_EB = new TCanvas("c_PV_track_time_PU200_prompt_EB", "c_PV_track_time_PU200_prompt_EB", 1500, 1500);
  c_PV_track_time_PU200_prompt_EB->cd();
  c_PV_track_time_PU200_prompt_EB->SetLeftMargin(0.12);
  h_PU200_prompt_PV_track_time_EB->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_PU200_prompt_PV_track_time_EB->SetLineColor(kBlack);
  h_PU200_prompt_PV_track_time_EB->SetLineWidth(2);
  h_PU200_prompt_PV_track_time_EB->Draw();
  c_PV_track_time_PU200_prompt_EB->Print("plots/ntuple/PV_track_time_PU200_prompt_EB.pdf");
  cout << "PU200 track: " << h_PU200_prompt_PV_track_time_EB->Integral(0,-1) << endl;
  cout << "PU200 track 0th  : " << h_PU200_prompt_PV_track_time_EB->GetBinContent(0) << endl;
  cout << "PU200 track 201th: " << h_PU200_prompt_PV_track_time_EB->GetBinContent(201) << endl;

  TCanvas* c_PU_track_time_PU200_prompt_EB = new TCanvas("c_PU_track_time_PU200_prompt_EB", "c_PU_track_time_PU200_prompt_EB", 1500, 1500);
  c_PU_track_time_PU200_prompt_EB->cd();
  c_PU_track_time_PU200_prompt_EB->SetLeftMargin(0.12);
  h_PU200_prompt_PU_track_time_EB->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_PU200_prompt_PU_track_time_EB->SetLineColor(kBlack);
  h_PU200_prompt_PU_track_time_EB->SetLineWidth(2);
  h_PU200_prompt_PU_track_time_EB->Draw();
  c_PU_track_time_PU200_prompt_EB->Print("plots/ntuple/PU_track_time_PU200_prompt_EB.pdf");
  cout << "PU200 track: " << h_PU200_prompt_PU_track_time_EB->Integral(0,-1) << endl;
  cout << "PU200 track 0th  : " << h_PU200_prompt_PU_track_time_EB->GetBinContent(0) << endl;
  cout << "PU200 track 201th: " << h_PU200_prompt_PU_track_time_EB->GetBinContent(201) << endl;

  TCanvas* c_SV_track_time_PU200_prompt_EB = new TCanvas("c_SV_track_time_PU200_prompt_EB", "c_SV_track_time_PU200_prompt_EB", 1500, 1500);
  c_SV_track_time_PU200_prompt_EB->cd();
  c_SV_track_time_PU200_prompt_EB->SetLeftMargin(0.12);
  h_PU200_prompt_SV_track_time_EB->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_PU200_prompt_SV_track_time_EB->SetLineColor(kBlack);
  h_PU200_prompt_SV_track_time_EB->SetLineWidth(2);
  h_PU200_prompt_SV_track_time_EB->Draw();
  c_SV_track_time_PU200_prompt_EB->Print("plots/ntuple/SV_track_time_PU200_prompt_EB.pdf");
  cout << "PU200 track: " << h_PU200_prompt_SV_track_time_EB->Integral(0,-1) << endl;
  cout << "PU200 track 0th  : " << h_PU200_prompt_SV_track_time_EB->GetBinContent(0) << endl;
  cout << "PU200 track 201th: " << h_PU200_prompt_SV_track_time_EB->GetBinContent(201) << endl;

  TCanvas* c_PV_track_time_PU200_nonprompt_EB = new TCanvas("c_PV_track_time_PU200_nonprompt_EB", "c_PV_track_time_PU200_nonprompt_EB", 1500, 1500);
  c_PV_track_time_PU200_nonprompt_EB->cd();
  c_PV_track_time_PU200_nonprompt_EB->SetLeftMargin(0.12);
  h_PU200_nonprompt_PV_track_time_EB->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_PU200_nonprompt_PV_track_time_EB->SetLineColor(kBlack);
  h_PU200_nonprompt_PV_track_time_EB->SetLineWidth(2);
  h_PU200_nonprompt_PV_track_time_EB->Draw();
  c_PV_track_time_PU200_nonprompt_EB->Print("plots/ntuple/PV_track_time_PU200_nonprompt_EB.pdf");
  cout << "PU200 track: " << h_PU200_nonprompt_PV_track_time_EB->Integral(0,-1) << endl;
  cout << "PU200 track 0th  : " << h_PU200_nonprompt_PV_track_time_EB->GetBinContent(0) << endl;
  cout << "PU200 track 201th: " << h_PU200_nonprompt_PV_track_time_EB->GetBinContent(201) << endl;

  TCanvas* c_PU_track_time_PU200_nonprompt_EB = new TCanvas("c_PU_track_time_PU200_nonprompt_EB", "c_PU_track_time_PU200_nonprompt_EB", 1500, 1500);
  c_PU_track_time_PU200_nonprompt_EB->cd();
  c_PU_track_time_PU200_nonprompt_EB->SetLeftMargin(0.12);
  h_PU200_nonprompt_PU_track_time_EB->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_PU200_nonprompt_PU_track_time_EB->SetLineColor(kBlack);
  h_PU200_nonprompt_PU_track_time_EB->SetLineWidth(2);
  h_PU200_nonprompt_PU_track_time_EB->Draw();
  c_PU_track_time_PU200_nonprompt_EB->Print("plots/ntuple/PU_track_time_PU200_nonprompt_EB.pdf");
  cout << "PU200 track: " << h_PU200_nonprompt_PU_track_time_EB->Integral(0,-1) << endl;
  cout << "PU200 track 0th  : " << h_PU200_nonprompt_PU_track_time_EB->GetBinContent(0) << endl;
  cout << "PU200 track 201th: " << h_PU200_nonprompt_PU_track_time_EB->GetBinContent(201) << endl;

  TCanvas* c_SV_track_time_PU200_nonprompt_EB = new TCanvas("c_SV_track_time_PU200_nonprompt_EB", "c_SV_track_time_PU200_nonprompt_EB", 1500, 1500);
  c_SV_track_time_PU200_nonprompt_EB->cd();
  c_SV_track_time_PU200_nonprompt_EB->SetLeftMargin(0.12);
  h_PU200_nonprompt_SV_track_time_EB->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_PU200_nonprompt_SV_track_time_EB->SetLineColor(kBlack);
  h_PU200_nonprompt_SV_track_time_EB->SetLineWidth(2);
  h_PU200_nonprompt_SV_track_time_EB->Draw();
  c_SV_track_time_PU200_nonprompt_EB->Print("plots/ntuple/SV_track_time_PU200_nonprompt_EB.pdf");
  cout << "PU200 track: " << h_PU200_nonprompt_SV_track_time_EB->Integral(0,-1) << endl;
  cout << "PU200 track 0th  : " << h_PU200_nonprompt_SV_track_time_EB->GetBinContent(0) << endl;
  cout << "PU200 track 201th: " << h_PU200_nonprompt_SV_track_time_EB->GetBinContent(201) << endl;
    // muon
  TCanvas* c_PV_muon_time_PU200_prompt_EB = new TCanvas("c_PV_muon_time_PU200_prompt_EB", "c_PV_muon_time_PU200_prompt_EB", 1500, 1500);
  c_PV_muon_time_PU200_prompt_EB->cd();
  c_PV_muon_time_PU200_prompt_EB->SetLeftMargin(0.12);
  h_PU200_prompt_PV_muon_time_EB->SetTitle(";t^{reco}_{muon} - t^{sim}_{muon} (ns);Counts");
  h_PU200_prompt_PV_muon_time_EB->SetLineColor(kBlack);
  h_PU200_prompt_PV_muon_time_EB->SetLineWidth(2);
  h_PU200_prompt_PV_muon_time_EB->Draw();
  c_PV_muon_time_PU200_prompt_EB->Print("plots/ntuple/PV_muon_time_PU200_prompt_EB.pdf");
  cout << "PU200 muon: " << h_PU200_prompt_PV_muon_time_EB->Integral(0,-1) << endl;
  cout << "PU200 muon 0th  : " << h_PU200_prompt_PV_muon_time_EB->GetBinContent(0) << endl;
  cout << "PU200 muon 201th: " << h_PU200_prompt_PV_muon_time_EB->GetBinContent(201) << endl;
    // vertex
  TCanvas* c_PV_vtx_time_PU200_prompt_EB = new TCanvas("c_PV_vtx_time_PU200_prompt_EB", "c_PV_vtx_time_PU200_prompt_EB", 1500, 1500);
  c_PV_vtx_time_PU200_prompt_EB->cd();
  c_PV_vtx_time_PU200_prompt_EB->SetLeftMargin(0.12);
  h_PU200_prompt_PV_vtx_time_EB->SetTitle(";t^{reco}_{PV} - t^{sim}_{PV} (ns);Counts");
  h_PU200_prompt_PV_vtx_time_EB->SetLineColor(kBlack);
  h_PU200_prompt_PV_vtx_time_EB->SetLineWidth(2);
  h_PU200_prompt_PV_vtx_time_EB->Draw();
  c_PV_vtx_time_PU200_prompt_EB->Print("plots/ntuple/PV_vtx_time_PU200_prompt_EB.pdf");
  cout << "PU200 vtx: " << h_PU200_prompt_PV_vtx_time_EB->Integral(0,-1) << endl;
  cout << "PU200 vtx 0th  : " << h_PU200_prompt_PV_vtx_time_EB->GetBinContent(0) << endl;
  cout << "PU200 vtx 201th: " << h_PU200_prompt_PV_vtx_time_EB->GetBinContent(201) << endl;

  // Endcap
    // track
  TCanvas* c_PV_track_time_PU200_prompt_EE = new TCanvas("c_PV_track_time_PU200_prompt_EE", "c_PV_track_time_PU200_prompt_EE", 1500, 1500);
  c_PV_track_time_PU200_prompt_EE->cd();
  c_PV_track_time_PU200_prompt_EE->SetLeftMargin(0.12);
  h_PU200_prompt_PV_track_time_EE->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_PU200_prompt_PV_track_time_EE->SetLineColor(kBlack);
  h_PU200_prompt_PV_track_time_EE->SetLineWidth(2);
  h_PU200_prompt_PV_track_time_EE->Draw();
  c_PV_track_time_PU200_prompt_EE->Print("plots/ntuple/PV_track_time_PU200_prompt_EE.pdf");
  cout << "PU200 track: " << h_PU200_prompt_PV_track_time_EE->Integral(0,-1) << endl;
  cout << "PU200 track 0th  : " << h_PU200_prompt_PV_track_time_EE->GetBinContent(0) << endl;
  cout << "PU200 track 201th: " << h_PU200_prompt_PV_track_time_EE->GetBinContent(201) << endl;

  TCanvas* c_PU_track_time_PU200_prompt_EE = new TCanvas("c_PU_track_time_PU200_prompt_EE", "c_PU_track_time_PU200_prompt_EE", 1500, 1500);
  c_PU_track_time_PU200_prompt_EE->cd();
  c_PU_track_time_PU200_prompt_EE->SetLeftMargin(0.12);
  h_PU200_prompt_PU_track_time_EE->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_PU200_prompt_PU_track_time_EE->SetLineColor(kBlack);
  h_PU200_prompt_PU_track_time_EE->SetLineWidth(2);
  h_PU200_prompt_PU_track_time_EE->Draw();
  c_PU_track_time_PU200_prompt_EE->Print("plots/ntuple/PU_track_time_PU200_prompt_EE.pdf");
  cout << "PU200 track: " << h_PU200_prompt_PU_track_time_EE->Integral(0,-1) << endl;
  cout << "PU200 track 0th  : " << h_PU200_prompt_PU_track_time_EE->GetBinContent(0) << endl;
  cout << "PU200 track 201th: " << h_PU200_prompt_PU_track_time_EE->GetBinContent(201) << endl;

  TCanvas* c_SV_track_time_PU200_prompt_EE = new TCanvas("c_SV_track_time_PU200_prompt_EE", "c_SV_track_time_PU200_prompt_EE", 1500, 1500);
  c_SV_track_time_PU200_prompt_EE->cd();
  c_SV_track_time_PU200_prompt_EE->SetLeftMargin(0.12);
  h_PU200_prompt_SV_track_time_EE->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_PU200_prompt_SV_track_time_EE->SetLineColor(kBlack);
  h_PU200_prompt_SV_track_time_EE->SetLineWidth(2);
  h_PU200_prompt_SV_track_time_EE->Draw();
  c_SV_track_time_PU200_prompt_EE->Print("plots/ntuple/SV_track_time_PU200_prompt_EE.pdf");
  cout << "PU200 track: " << h_PU200_prompt_SV_track_time_EE->Integral(0,-1) << endl;
  cout << "PU200 track 0th  : " << h_PU200_prompt_SV_track_time_EE->GetBinContent(0) << endl;
  cout << "PU200 track 201th: " << h_PU200_prompt_SV_track_time_EE->GetBinContent(201) << endl;

  TCanvas* c_PV_track_time_PU200_nonprompt_EE = new TCanvas("c_PV_track_time_PU200_nonprompt_EE", "c_PV_track_time_PU200_nonprompt_EE", 1500, 1500);
  c_PV_track_time_PU200_nonprompt_EE->cd();
  c_PV_track_time_PU200_nonprompt_EE->SetLeftMargin(0.12);
  h_PU200_nonprompt_PV_track_time_EE->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_PU200_nonprompt_PV_track_time_EE->SetLineColor(kBlack);
  h_PU200_nonprompt_PV_track_time_EE->SetLineWidth(2);
  h_PU200_nonprompt_PV_track_time_EE->Draw();
  c_PV_track_time_PU200_nonprompt_EE->Print("plots/ntuple/PV_track_time_PU200_nonprompt_EE.pdf");
  cout << "PU200 track: " << h_PU200_nonprompt_PV_track_time_EE->Integral(0,-1) << endl;
  cout << "PU200 track 0th  : " << h_PU200_nonprompt_PV_track_time_EE->GetBinContent(0) << endl;
  cout << "PU200 track 201th: " << h_PU200_nonprompt_PV_track_time_EE->GetBinContent(201) << endl;

  TCanvas* c_PU_track_time_PU200_nonprompt_EE = new TCanvas("c_PU_track_time_PU200_nonprompt_EE", "c_PU_track_time_PU200_nonprompt_EE", 1500, 1500);
  c_PU_track_time_PU200_nonprompt_EE->cd();
  c_PU_track_time_PU200_nonprompt_EE->SetLeftMargin(0.12);
  h_PU200_nonprompt_PU_track_time_EE->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_PU200_nonprompt_PU_track_time_EE->SetLineColor(kBlack);
  h_PU200_nonprompt_PU_track_time_EE->SetLineWidth(2);
  h_PU200_nonprompt_PU_track_time_EE->Draw();
  c_PU_track_time_PU200_nonprompt_EE->Print("plots/ntuple/PU_track_time_PU200_nonprompt_EE.pdf");
  cout << "PU200 track: " << h_PU200_nonprompt_PU_track_time_EE->Integral(0,-1) << endl;
  cout << "PU200 track 0th  : " << h_PU200_nonprompt_PU_track_time_EE->GetBinContent(0) << endl;
  cout << "PU200 track 201th: " << h_PU200_nonprompt_PU_track_time_EE->GetBinContent(201) << endl;

  TCanvas* c_SV_track_time_PU200_nonprompt_EE = new TCanvas("c_SV_track_time_PU200_nonprompt_EE", "c_SV_track_time_PU200_nonprompt_EE", 1500, 1500);
  c_SV_track_time_PU200_nonprompt_EE->cd();
  c_SV_track_time_PU200_nonprompt_EE->SetLeftMargin(0.12);
  h_PU200_nonprompt_SV_track_time_EE->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_PU200_nonprompt_SV_track_time_EE->SetLineColor(kBlack);
  h_PU200_nonprompt_SV_track_time_EE->SetLineWidth(2);
  h_PU200_nonprompt_SV_track_time_EE->Draw();
  c_SV_track_time_PU200_nonprompt_EE->Print("plots/ntuple/SV_track_time_PU200_nonprompt_EE.pdf");
  cout << "PU200 track: " << h_PU200_nonprompt_SV_track_time_EE->Integral(0,-1) << endl;
  cout << "PU200 track 0th  : " << h_PU200_nonprompt_SV_track_time_EE->GetBinContent(0) << endl;
  cout << "PU200 track 201th: " << h_PU200_nonprompt_SV_track_time_EE->GetBinContent(201) << endl;
    // muon
  TCanvas* c_PV_muon_time_PU200_prompt_EE = new TCanvas("c_PV_muon_time_PU200_prompt_EE", "c_PV_muon_time_PU200_prompt_EE", 1500, 1500);
  c_PV_muon_time_PU200_prompt_EE->cd();
  c_PV_muon_time_PU200_prompt_EE->SetLeftMargin(0.12);
  h_PU200_prompt_PV_muon_time_EE->SetTitle(";t^{reco}_{muon} - t^{sim}_{muon} (ns);Counts");
  h_PU200_prompt_PV_muon_time_EE->SetLineColor(kBlack);
  h_PU200_prompt_PV_muon_time_EE->SetLineWidth(2);
  h_PU200_prompt_PV_muon_time_EE->Draw();
  c_PV_muon_time_PU200_prompt_EE->Print("plots/ntuple/PV_muon_time_PU200_prompt_EE.pdf");
  cout << "PU200 muon: " << h_PU200_prompt_PV_muon_time_EE->Integral(0,-1) << endl;
  cout << "PU200 muon 0th  : " << h_PU200_prompt_PV_muon_time_EE->GetBinContent(0) << endl;
  cout << "PU200 muon 201th: " << h_PU200_prompt_PV_muon_time_EE->GetBinContent(201) << endl;
    // vertex
  TCanvas* c_PV_vtx_time_PU200_prompt_EE = new TCanvas("c_PV_vtx_time_PU200_prompt_EE", "c_PV_vtx_time_PU200_prompt_EE", 1500, 1500);
  c_PV_vtx_time_PU200_prompt_EE->cd();
  c_PV_vtx_time_PU200_prompt_EE->SetLeftMargin(0.12);
  h_PU200_prompt_PV_vtx_time_EE->SetTitle(";t^{reco}_{PV} - t^{sim}_{PV} (ns);Counts");
  h_PU200_prompt_PV_vtx_time_EE->SetLineColor(kBlack);
  h_PU200_prompt_PV_vtx_time_EE->SetLineWidth(2);
  h_PU200_prompt_PV_vtx_time_EE->Draw();
  c_PV_vtx_time_PU200_prompt_EE->Print("plots/ntuple/PV_vtx_time_PU200_prompt_EE.pdf");
  cout << "PU200 vtx: " << h_PU200_prompt_PV_vtx_time_EE->Integral(0,-1) << endl;
  cout << "PU200 vtx 0th  : " << h_PU200_prompt_PV_vtx_time_EE->GetBinContent(0) << endl;
  cout << "PU200 vtx 201th: " << h_PU200_prompt_PV_vtx_time_EE->GetBinContent(201) << endl;

  // noPU
  // Barrel
    // track
  TCanvas* c_PV_track_time_noPU_prompt_EB = new TCanvas("c_PV_track_time_noPU_prompt_EB", "c_PV_track_time_noPU_prompt_EB", 1500, 1500);
  c_PV_track_time_noPU_prompt_EB->cd();
  c_PV_track_time_noPU_prompt_EB->SetLeftMargin(0.12);
  h_noPU_prompt_PV_track_time_EB->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_noPU_prompt_PV_track_time_EB->SetLineColor(kBlack);
  h_noPU_prompt_PV_track_time_EB->SetLineWidth(2);
  h_noPU_prompt_PV_track_time_EB->Draw();
  c_PV_track_time_noPU_prompt_EB->Print("plots/ntuple/PV_track_time_noPU_prompt_EB.pdf");
  cout << "noPU track: " << h_noPU_prompt_PV_track_time_EB->Integral(0,-1) << endl;
  cout << "noPU track 0th  : " << h_noPU_prompt_PV_track_time_EB->GetBinContent(0) << endl;
  cout << "noPU track 201th: " << h_noPU_prompt_PV_track_time_EB->GetBinContent(201) << endl;

  TCanvas* c_PU_track_time_noPU_prompt_EB = new TCanvas("c_PU_track_time_noPU_prompt_EB", "c_PU_track_time_noPU_prompt_EB", 1500, 1500);
  c_PU_track_time_noPU_prompt_EB->cd();
  c_PU_track_time_noPU_prompt_EB->SetLeftMargin(0.12);
  h_noPU_prompt_PU_track_time_EB->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_noPU_prompt_PU_track_time_EB->SetLineColor(kBlack);
  h_noPU_prompt_PU_track_time_EB->SetLineWidth(2);
  h_noPU_prompt_PU_track_time_EB->Draw();
  c_PU_track_time_noPU_prompt_EB->Print("plots/ntuple/PU_track_time_noPU_prompt_EB.pdf");
  cout << "noPU track: " << h_noPU_prompt_PU_track_time_EB->Integral(0,-1) << endl;
  cout << "noPU track 0th  : " << h_noPU_prompt_PU_track_time_EB->GetBinContent(0) << endl;
  cout << "noPU track 201th: " << h_noPU_prompt_PU_track_time_EB->GetBinContent(201) << endl;

  TCanvas* c_SV_track_time_noPU_prompt_EB = new TCanvas("c_SV_track_time_noPU_prompt_EB", "c_SV_track_time_noPU_prompt_EB", 1500, 1500);
  c_SV_track_time_noPU_prompt_EB->cd();
  c_SV_track_time_noPU_prompt_EB->SetLeftMargin(0.12);
  h_noPU_prompt_SV_track_time_EB->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_noPU_prompt_SV_track_time_EB->SetLineColor(kBlack);
  h_noPU_prompt_SV_track_time_EB->SetLineWidth(2);
  h_noPU_prompt_SV_track_time_EB->Draw();
  c_SV_track_time_noPU_prompt_EB->Print("plots/ntuple/SV_track_time_noPU_prompt_EB.pdf");
  cout << "noPU track: " << h_noPU_prompt_SV_track_time_EB->Integral(0,-1) << endl;
  cout << "noPU track 0th  : " << h_noPU_prompt_SV_track_time_EB->GetBinContent(0) << endl;
  cout << "noPU track 201th: " << h_noPU_prompt_SV_track_time_EB->GetBinContent(201) << endl;

  TCanvas* c_PV_track_time_noPU_nonprompt_EB = new TCanvas("c_PV_track_time_noPU_nonprompt_EB", "c_PV_track_time_noPU_nonprompt_EB", 1500, 1500);
  c_PV_track_time_noPU_nonprompt_EB->cd();
  c_PV_track_time_noPU_nonprompt_EB->SetLeftMargin(0.12);
  h_noPU_nonprompt_PV_track_time_EB->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_noPU_nonprompt_PV_track_time_EB->SetLineColor(kBlack);
  h_noPU_nonprompt_PV_track_time_EB->SetLineWidth(2);
  h_noPU_nonprompt_PV_track_time_EB->Draw();
  c_PV_track_time_noPU_nonprompt_EB->Print("plots/ntuple/PV_track_time_noPU_nonprompt_EB.pdf");
  cout << "noPU track: " << h_noPU_nonprompt_PV_track_time_EB->Integral(0,-1) << endl;
  cout << "noPU track 0th  : " << h_noPU_nonprompt_PV_track_time_EB->GetBinContent(0) << endl;
  cout << "noPU track 201th: " << h_noPU_nonprompt_PV_track_time_EB->GetBinContent(201) << endl;

  TCanvas* c_PU_track_time_noPU_nonprompt_EB = new TCanvas("c_PU_track_time_noPU_nonprompt_EB", "c_PU_track_time_noPU_nonprompt_EB", 1500, 1500);
  c_PU_track_time_noPU_nonprompt_EB->cd();
  c_PU_track_time_noPU_nonprompt_EB->SetLeftMargin(0.12);
  h_noPU_nonprompt_PU_track_time_EB->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_noPU_nonprompt_PU_track_time_EB->SetLineColor(kBlack);
  h_noPU_nonprompt_PU_track_time_EB->SetLineWidth(2);
  h_noPU_nonprompt_PU_track_time_EB->Draw();
  c_PU_track_time_noPU_nonprompt_EB->Print("plots/ntuple/PU_track_time_noPU_nonprompt_EB.pdf");
  cout << "noPU track: " << h_noPU_nonprompt_PU_track_time_EB->Integral(0,-1) << endl;
  cout << "noPU track 0th  : " << h_noPU_nonprompt_PU_track_time_EB->GetBinContent(0) << endl;
  cout << "noPU track 201th: " << h_noPU_nonprompt_PU_track_time_EB->GetBinContent(201) << endl;

  TCanvas* c_SV_track_time_noPU_nonprompt_EB = new TCanvas("c_SV_track_time_noPU_nonprompt_EB", "c_SV_track_time_noPU_nonprompt_EB", 1500, 1500);
  c_SV_track_time_noPU_nonprompt_EB->cd();
  c_SV_track_time_noPU_nonprompt_EB->SetLeftMargin(0.12);
  h_noPU_nonprompt_SV_track_time_EB->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_noPU_nonprompt_SV_track_time_EB->SetLineColor(kBlack);
  h_noPU_nonprompt_SV_track_time_EB->SetLineWidth(2);
  h_noPU_nonprompt_SV_track_time_EB->Draw();
  c_SV_track_time_noPU_nonprompt_EB->Print("plots/ntuple/SV_track_time_noPU_nonprompt_EB.pdf");
  cout << "noPU track: " << h_noPU_nonprompt_SV_track_time_EB->Integral(0,-1) << endl;
  cout << "noPU track 0th  : " << h_noPU_nonprompt_SV_track_time_EB->GetBinContent(0) << endl;
  cout << "noPU track 201th: " << h_noPU_nonprompt_SV_track_time_EB->GetBinContent(201) << endl;
    // muon
  TCanvas* c_PV_muon_time_noPU_prompt_EB = new TCanvas("c_PV_muon_time_noPU_prompt_EB", "c_PV_muon_time_noPU_prompt_EB", 1500, 1500);
  c_PV_muon_time_noPU_prompt_EB->cd();
  c_PV_muon_time_noPU_prompt_EB->SetLeftMargin(0.12);
  h_noPU_prompt_PV_muon_time_EB->SetTitle(";t^{reco}_{muon} - t^{sim}_{muon} (ns);Counts");
  h_noPU_prompt_PV_muon_time_EB->SetLineColor(kBlack);
  h_noPU_prompt_PV_muon_time_EB->SetLineWidth(2);
  h_noPU_prompt_PV_muon_time_EB->Draw();
  c_PV_muon_time_noPU_prompt_EB->Print("plots/ntuple/PV_muon_time_noPU_prompt_EB.pdf");
  cout << "noPU muon: " << h_noPU_prompt_PV_muon_time_EB->Integral(0,-1) << endl;
  cout << "noPU muon 0th  : " << h_noPU_prompt_PV_muon_time_EB->GetBinContent(0) << endl;
  cout << "noPU muon 201th: " << h_noPU_prompt_PV_muon_time_EB->GetBinContent(201) << endl;
    // vertex
  TCanvas* c_PV_vtx_time_noPU_prompt_EB = new TCanvas("c_PV_vtx_time_noPU_prompt_EB", "c_PV_vtx_time_noPU_prompt_EB", 1500, 1500);
  c_PV_vtx_time_noPU_prompt_EB->cd();
  c_PV_vtx_time_noPU_prompt_EB->SetLeftMargin(0.12);
  h_noPU_prompt_PV_vtx_time_EB->SetTitle(";t^{reco}_{PV} - t^{sim}_{PV} (ns);Counts");
  h_noPU_prompt_PV_vtx_time_EB->SetLineColor(kBlack);
  h_noPU_prompt_PV_vtx_time_EB->SetLineWidth(2);
  h_noPU_prompt_PV_vtx_time_EB->Draw();
  c_PV_vtx_time_noPU_prompt_EB->Print("plots/ntuple/PV_vtx_time_noPU_prompt_EB.pdf");
  cout << "noPU vtx: " << h_noPU_prompt_PV_vtx_time_EB->Integral(0,-1) << endl;
  cout << "noPU vtx 0th  : " << h_noPU_prompt_PV_vtx_time_EB->GetBinContent(0) << endl;
  cout << "noPU vtx 201th: " << h_noPU_prompt_PV_vtx_time_EB->GetBinContent(201) << endl;

  // Endcap
    // track
  TCanvas* c_PV_track_time_noPU_prompt_EE = new TCanvas("c_PV_track_time_noPU_prompt_EE", "c_PV_track_time_noPU_prompt_EE", 1500, 1500);
  c_PV_track_time_noPU_prompt_EE->cd();
  c_PV_track_time_noPU_prompt_EE->SetLeftMargin(0.12);
  h_noPU_prompt_PV_track_time_EE->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_noPU_prompt_PV_track_time_EE->SetLineColor(kBlack);
  h_noPU_prompt_PV_track_time_EE->SetLineWidth(2);
  h_noPU_prompt_PV_track_time_EE->Draw();
  c_PV_track_time_noPU_prompt_EE->Print("plots/ntuple/PV_track_time_noPU_prompt_EE.pdf");
  cout << "noPU track: " << h_noPU_prompt_PV_track_time_EE->Integral(0,-1) << endl;
  cout << "noPU track 0th  : " << h_noPU_prompt_PV_track_time_EE->GetBinContent(0) << endl;
  cout << "noPU track 201th: " << h_noPU_prompt_PV_track_time_EE->GetBinContent(201) << endl;

  TCanvas* c_PU_track_time_noPU_prompt_EE = new TCanvas("c_PU_track_time_noPU_prompt_EE", "c_PU_track_time_noPU_prompt_EE", 1500, 1500);
  c_PU_track_time_noPU_prompt_EE->cd();
  c_PU_track_time_noPU_prompt_EE->SetLeftMargin(0.12);
  h_noPU_prompt_PU_track_time_EE->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_noPU_prompt_PU_track_time_EE->SetLineColor(kBlack);
  h_noPU_prompt_PU_track_time_EE->SetLineWidth(2);
  h_noPU_prompt_PU_track_time_EE->Draw();
  c_PU_track_time_noPU_prompt_EE->Print("plots/ntuple/PU_track_time_noPU_prompt_EE.pdf");
  cout << "noPU track: " << h_noPU_prompt_PU_track_time_EE->Integral(0,-1) << endl;
  cout << "noPU track 0th  : " << h_noPU_prompt_PU_track_time_EE->GetBinContent(0) << endl;
  cout << "noPU track 201th: " << h_noPU_prompt_PU_track_time_EE->GetBinContent(201) << endl;

  TCanvas* c_SV_track_time_noPU_prompt_EE = new TCanvas("c_SV_track_time_noPU_prompt_EE", "c_SV_track_time_noPU_prompt_EE", 1500, 1500);
  c_SV_track_time_noPU_prompt_EE->cd();
  c_SV_track_time_noPU_prompt_EE->SetLeftMargin(0.12);
  h_noPU_prompt_SV_track_time_EE->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_noPU_prompt_SV_track_time_EE->SetLineColor(kBlack);
  h_noPU_prompt_SV_track_time_EE->SetLineWidth(2);
  h_noPU_prompt_SV_track_time_EE->Draw();
  c_SV_track_time_noPU_prompt_EE->Print("plots/ntuple/SV_track_time_noPU_prompt_EE.pdf");
  cout << "noPU track: " << h_noPU_prompt_SV_track_time_EE->Integral(0,-1) << endl;
  cout << "noPU track 0th  : " << h_noPU_prompt_SV_track_time_EE->GetBinContent(0) << endl;
  cout << "noPU track 201th: " << h_noPU_prompt_SV_track_time_EE->GetBinContent(201) << endl;

  TCanvas* c_PV_track_time_noPU_nonprompt_EE = new TCanvas("c_PV_track_time_noPU_nonprompt_EE", "c_PV_track_time_noPU_nonprompt_EE", 1500, 1500);
  c_PV_track_time_noPU_nonprompt_EE->cd();
  c_PV_track_time_noPU_nonprompt_EE->SetLeftMargin(0.12);
  h_noPU_nonprompt_PV_track_time_EE->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_noPU_nonprompt_PV_track_time_EE->SetLineColor(kBlack);
  h_noPU_nonprompt_PV_track_time_EE->SetLineWidth(2);
  h_noPU_nonprompt_PV_track_time_EE->Draw();
  c_PV_track_time_noPU_nonprompt_EE->Print("plots/ntuple/PV_track_time_noPU_nonprompt_EE.pdf");
  cout << "noPU track: " << h_noPU_nonprompt_PV_track_time_EE->Integral(0,-1) << endl;
  cout << "noPU track 0th  : " << h_noPU_nonprompt_PV_track_time_EE->GetBinContent(0) << endl;
  cout << "noPU track 201th: " << h_noPU_nonprompt_PV_track_time_EE->GetBinContent(201) << endl;

  TCanvas* c_PU_track_time_noPU_nonprompt_EE = new TCanvas("c_PU_track_time_noPU_nonprompt_EE", "c_PU_track_time_noPU_nonprompt_EE", 1500, 1500);
  c_PU_track_time_noPU_nonprompt_EE->cd();
  c_PU_track_time_noPU_nonprompt_EE->SetLeftMargin(0.12);
  h_noPU_nonprompt_PU_track_time_EE->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_noPU_nonprompt_PU_track_time_EE->SetLineColor(kBlack);
  h_noPU_nonprompt_PU_track_time_EE->SetLineWidth(2);
  h_noPU_nonprompt_PU_track_time_EE->Draw();
  c_PU_track_time_noPU_nonprompt_EE->Print("plots/ntuple/PU_track_time_noPU_nonprompt_EE.pdf");
  cout << "noPU track: " << h_noPU_nonprompt_PU_track_time_EE->Integral(0,-1) << endl;
  cout << "noPU track 0th  : " << h_noPU_nonprompt_PU_track_time_EE->GetBinContent(0) << endl;
  cout << "noPU track 201th: " << h_noPU_nonprompt_PU_track_time_EE->GetBinContent(201) << endl;

  TCanvas* c_SV_track_time_noPU_nonprompt_EE = new TCanvas("c_SV_track_time_noPU_nonprompt_EE", "c_SV_track_time_noPU_nonprompt_EE", 1500, 1500);
  c_SV_track_time_noPU_nonprompt_EE->cd();
  c_SV_track_time_noPU_nonprompt_EE->SetLeftMargin(0.12);
  h_noPU_nonprompt_SV_track_time_EE->SetTitle(";t^{reco}_{track} - t^{sim}_{track} (ns);Counts");
  h_noPU_nonprompt_SV_track_time_EE->SetLineColor(kBlack);
  h_noPU_nonprompt_SV_track_time_EE->SetLineWidth(2);
  h_noPU_nonprompt_SV_track_time_EE->Draw();
  c_SV_track_time_noPU_nonprompt_EE->Print("plots/ntuple/SV_track_time_noPU_nonprompt_EE.pdf");
  cout << "noPU track: " << h_noPU_nonprompt_SV_track_time_EE->Integral(0,-1) << endl;
  cout << "noPU track 0th  : " << h_noPU_nonprompt_SV_track_time_EE->GetBinContent(0) << endl;
  cout << "noPU track 201th: " << h_noPU_nonprompt_SV_track_time_EE->GetBinContent(201) << endl;
    // muon
  TCanvas* c_PV_muon_time_noPU_prompt_EE = new TCanvas("c_PV_muon_time_noPU_prompt_EE", "c_PV_muon_time_noPU_prompt_EE", 1500, 1500);
  c_PV_muon_time_noPU_prompt_EE->cd();
  c_PV_muon_time_noPU_prompt_EE->SetLeftMargin(0.12);
  h_noPU_prompt_PV_muon_time_EE->SetTitle(";t^{reco}_{muon} - t^{sim}_{muon} (ns);Counts");
  h_noPU_prompt_PV_muon_time_EE->SetLineColor(kBlack);
  h_noPU_prompt_PV_muon_time_EE->SetLineWidth(2);
  h_noPU_prompt_PV_muon_time_EE->Draw();
  c_PV_muon_time_noPU_prompt_EE->Print("plots/ntuple/PV_muon_time_noPU_prompt_EE.pdf");
  cout << "noPU muon: " << h_noPU_prompt_PV_muon_time_EE->Integral(0,-1) << endl;
  cout << "noPU muon 0th  : " << h_noPU_prompt_PV_muon_time_EE->GetBinContent(0) << endl;
  cout << "noPU muon 201th: " << h_noPU_prompt_PV_muon_time_EE->GetBinContent(201) << endl;
    // vertex
  TCanvas* c_PV_vtx_time_noPU_prompt_EE = new TCanvas("c_PV_vtx_time_noPU_prompt_EE", "c_PV_vtx_time_noPU_prompt_EE", 1500, 1500);
  c_PV_vtx_time_noPU_prompt_EE->cd();
  c_PV_vtx_time_noPU_prompt_EE->SetLeftMargin(0.12);
  h_noPU_prompt_PV_vtx_time_EE->SetTitle(";t^{reco}_{PV} - t^{sim}_{PV} (ns);Counts");
  h_noPU_prompt_PV_vtx_time_EE->SetLineColor(kBlack);
  h_noPU_prompt_PV_vtx_time_EE->SetLineWidth(2);
  h_noPU_prompt_PV_vtx_time_EE->Draw();
  c_PV_vtx_time_noPU_prompt_EE->Print("plots/ntuple/PV_vtx_time_noPU_prompt_EE.pdf");
  cout << "noPU vtx: " << h_noPU_prompt_PV_vtx_time_EE->Integral(0,-1) << endl;
  cout << "noPU vtx 0th  : " << h_noPU_prompt_PV_vtx_time_EE->GetBinContent(0) << endl;
  cout << "noPU vtx 201th: " << h_noPU_prompt_PV_vtx_time_EE->GetBinContent(201) << endl;



  // tsim end
  cout << endl;






  TCanvas* c_PV_dtsig_PU200_prompt_EE = new TCanvas("c_PV_dtsig_PU200_prompt_EE", "c_PV_dtsig_PU200_prompt_EE", 1500, 1500);
  c_PV_dtsig_PU200_prompt_EE->cd();
  c_PV_dtsig_PU200_prompt_EE->SetLeftMargin(0.12);
  h_PU200_prompt_PV_EE->SetTitle(";#sigma_{t};Counts");
  h_PU200_prompt_PV_EE->Draw();
  c_PV_dtsig_PU200_prompt_EE->Print("plots/ntuple/PV_dtsig_PU200_prompt_EE.pdf");


  TCanvas* c_PV_PVweight_noPU_prompt_EB = new TCanvas("c_PV_PVweight_noPU_prompt_EB", "c_PV_PVweight_noPU_prompt_EB", 1500, 1500);
  c_PV_PVweight_noPU_prompt_EB->cd();
  c_PV_PVweight_noPU_prompt_EB->SetLeftMargin(0.12);
  h_noPU_prompt_PVweight_PV_EB->SetTitle(";PV weight of track from PV;Counts");
  h_noPU_prompt_PVweight_PV_EB->SetLineColor(kBlack);
  h_noPU_prompt_PVweight_PV_EB->SetLineWidth(2);
  h_noPU_prompt_PVweight_PV_EB->Draw();
  //leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_PVweight_noPU_prompt_EB->Print("plots/ntuple/PV_PVweight_noPU_prompt_EB.pdf");
  cout << "PV weight = 0  : " << h_noPU_prompt_PVweight_PV_EB->GetBinContent(1) << endl;
  cout << "less than 0.5  : " << h_noPU_prompt_PVweight_PV_EB->Integral(0,50) << endl;
  cout << "total PV weight: " << h_noPU_prompt_PVweight_PV_EB->Integral(0,-1) << endl;
  cout << endl;




  // (PV, track)
  // PU200
    // prompt
      // Barrel
  TCanvas* c_st_PU200_prompt_EB_vtx_v2 = new TCanvas("c_st_PU200_prompt_EB_vtx_v2", "c_st_PU200_prompt_EB_vtx_v2", 1500, 1500);
  c_st_PU200_prompt_EB_vtx_v2->cd();
  c_st_PU200_prompt_EB_vtx_v2->SetLeftMargin(0.12);
  st_PU200_prompt_EB_vtx->Draw();
  leg_PU200_prompt_EB->Draw();
  c_st_PU200_prompt_EB_vtx_v2->Print("plots/ntuple/track_sigma_PU200_prompt_EB_vtx.pdf");
  cout << "///////////////////////" << endl;
  cout << "///// vtx, track //////" << endl;
  cout << "///////////////////////" << endl;
  cout << "/////   2sigma   //////" << endl;
  cout << "///////////////////////" << endl;
  cout << "PU200 Barrel prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PV_EB_vtx->Integral(1,2) << " : " << h_PU200_prompt_PV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EB_vtx->Integral(3,-1)/h_PU200_prompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_SV_EB_vtx->Integral(1,2) << " : " << h_PU200_prompt_SV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EB_vtx->Integral(3,-1)/h_PU200_prompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PU_EB_vtx->Integral(1,2) << " : " << h_PU200_prompt_PU_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EB_vtx->Integral(3,-1)/h_PU200_prompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_fake_EB_vtx->Integral(1,2) << " : " << h_PU200_prompt_fake_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EB_vtx->Integral(3,-1)/h_PU200_prompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_prompt_not_vtx_matched_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_not_vtx_matched_EB_vtx->Integral(1,2) << " : " << h_PU200_prompt_not_vtx_matched_EB_vtx->Integral(3,-1) << "  " << "(" << 100 -h_PU200_prompt_not_vtx_matched_EB_vtx->Integral(3,-1)/h_PU200_prompt_not_vtx_matched_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EB_vtx->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_PU200_prompt_EE_vtx_v2 = new TCanvas("c_st_PU200_prompt_EE_vtx_v2", "c_st_PU200_prompt_EE_vtx_v2", 1500, 1500);
  c_st_PU200_prompt_EE_vtx_v2->cd();
  c_st_PU200_prompt_EE_vtx_v2->SetLeftMargin(0.12);
  st_PU200_prompt_EE_vtx->Draw();
  leg_PU200_prompt_EE->Draw();
  c_st_PU200_prompt_EE_vtx_v2->Print("plots/ntuple/track_sigma_PU200_prompt_EE_vtx.pdf");
  cout << "PU200 Endcap prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PV_EE_vtx->Integral(1,2) << " : " << h_PU200_prompt_PV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EE_vtx->Integral(3,-1)/h_PU200_prompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_SV_EE_vtx->Integral(1,2) << " : " << h_PU200_prompt_SV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EE_vtx->Integral(3,-1)/h_PU200_prompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PU_EE_vtx->Integral(1,2) << " : " << h_PU200_prompt_PU_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EE_vtx->Integral(3,-1)/h_PU200_prompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_fake_EE_vtx->Integral(1,2) << " : " << h_PU200_prompt_fake_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EE_vtx->Integral(3,-1)/h_PU200_prompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_prompt_not_vtx_matched_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_not_vtx_matched_EE_vtx->Integral(1,2) << " : " << h_PU200_prompt_not_vtx_matched_EE_vtx->Integral(3,-1) << "  " << "(" << 100 -h_PU200_prompt_not_vtx_matched_EE_vtx->Integral(3,-1)/h_PU200_prompt_not_vtx_matched_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EE_vtx->Integral() << endl;
  cout << endl;
    // nonprompt
      // Barrel
  TCanvas* c_st_PU200_nonprompt_EB_vtx_v2 = new TCanvas("c_st_PU200_nonprompt_EB_vtx_v2", "c_st_PU200_nonprompt_EB_vtx_v2", 1500, 1500);
  c_st_PU200_nonprompt_EB_vtx_v2->cd();
  c_st_PU200_nonprompt_EB_vtx_v2->SetLeftMargin(0.12);
  st_PU200_nonprompt_EB_vtx->Draw();
  leg_PU200_nonprompt_EB->Draw();
  c_st_PU200_nonprompt_EB_vtx_v2->Print("plots/ntuple/track_sigma_PU200_nonprompt_EB_vtx.pdf");
  cout << "PU200 Barrel nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EB_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_PV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EB_vtx->Integral(3,-1)/h_PU200_nonprompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EB_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_SV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EB_vtx->Integral(3,-1)/h_PU200_nonprompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EB_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_PU_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EB_vtx->Integral(3,-1)/h_PU200_nonprompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EB_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_fake_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EB_vtx->Integral(3,-1)/h_PU200_nonprompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_nonprompt_not_vtx_matched_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_not_vtx_matched_EB_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_not_vtx_matched_EB_vtx->Integral(3,-1) << "  " << "(" << 100 -h_PU200_nonprompt_not_vtx_matched_EB_vtx->Integral(3,-1)/h_PU200_nonprompt_not_vtx_matched_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EB_vtx->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_PU200_nonprompt_EE_vtx_v2 = new TCanvas("c_st_PU200_nonprompt_EE_vtx_v2", "c_st_PU200_nonprompt_EE_vtx_v2", 1500, 1500);
  c_st_PU200_nonprompt_EE_vtx_v2->cd();
  c_st_PU200_nonprompt_EE_vtx_v2->SetLeftMargin(0.12);
  st_PU200_nonprompt_EE_vtx->Draw();
  leg_PU200_nonprompt_EE->Draw();
  c_st_PU200_nonprompt_EE_vtx_v2->Print("plots/ntuple/track_sigma_PU200_nonprompt_EE_vtx.pdf");
  cout << "PU200 Endcap nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EE_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_PV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EE_vtx->Integral(3,-1)/h_PU200_nonprompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EE_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_SV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EE_vtx->Integral(3,-1)/h_PU200_nonprompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EE_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_PU_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EE_vtx->Integral(3,-1)/h_PU200_nonprompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EE_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_fake_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EE_vtx->Integral(3,-1)/h_PU200_nonprompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_nonprompt_not_vtx_matched_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_not_vtx_matched_EE_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_not_vtx_matched_EE_vtx->Integral(3,-1) << "  " << "(" << 100 -h_PU200_nonprompt_not_vtx_matched_EE_vtx->Integral(3,-1)/h_PU200_nonprompt_not_vtx_matched_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EE_vtx->Integral() << endl;
  cout << endl;

  // noPU
    // prompt
      // Barrel
  TCanvas* c_st_noPU_prompt_EB_vtx_v2 = new TCanvas("c_st_noPU_prompt_EB_vtx_v2", "c_st_noPU_prompt_EB_vtx_v2", 1500, 1500);
  c_st_noPU_prompt_EB_vtx_v2->cd();
  c_st_noPU_prompt_EB_vtx_v2->SetLeftMargin(0.12);
  st_noPU_prompt_EB_vtx->Draw();
  leg_noPU_prompt_EB->Draw();
  c_st_noPU_prompt_EB_vtx_v2->Print("plots/ntuple/track_sigma_noPU_prompt_EB_vtx.pdf");
  cout << "noPU Barrel prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PV_EB_vtx->Integral(1,2) << " : " << h_noPU_prompt_PV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EB_vtx->Integral(3,-1)/h_noPU_prompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_SV_EB_vtx->Integral(1,2) << " : " << h_noPU_prompt_SV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EB_vtx->Integral(3,-1)/h_noPU_prompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PU_EB_vtx->Integral(1,2) << " : " << h_noPU_prompt_PU_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EB_vtx->Integral(3,-1)/h_noPU_prompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_fake_EB_vtx->Integral(1,2) << " : " << h_noPU_prompt_fake_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EB_vtx->Integral(3,-1)/h_noPU_prompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_prompt_not_vtx_matched_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_not_vtx_matched_EB_vtx->Integral(1,2) << " : " << h_noPU_prompt_not_vtx_matched_EB_vtx->Integral(3,-1) << "  " << "(" << 100 -h_noPU_prompt_not_vtx_matched_EB_vtx->Integral(3,-1)/h_noPU_prompt_not_vtx_matched_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EB_vtx->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_noPU_prompt_EE_vtx_v2 = new TCanvas("c_st_noPU_prompt_EE_vtx_v2", "c_st_noPU_prompt_EE_vtx_v2", 1500, 1500);
  c_st_noPU_prompt_EE_vtx_v2->cd();
  c_st_noPU_prompt_EE_vtx_v2->SetLeftMargin(0.12);
  st_noPU_prompt_EE_vtx->Draw();
  leg_noPU_prompt_EE->Draw();
  c_st_noPU_prompt_EE_vtx_v2->Print("plots/ntuple/track_sigma_noPU_prompt_EE_vtx.pdf");
  cout << "noPU Endcap prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PV_EE_vtx->Integral(1,2) << " : " << h_noPU_prompt_PV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EE_vtx->Integral(3,-1)/h_noPU_prompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_SV_EE_vtx->Integral(1,2) << " : " << h_noPU_prompt_SV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EE_vtx->Integral(3,-1)/h_noPU_prompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PU_EE_vtx->Integral(1,2) << " : " << h_noPU_prompt_PU_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EE_vtx->Integral(3,-1)/h_noPU_prompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_fake_EE_vtx->Integral(1,2) << " : " << h_noPU_prompt_fake_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EE_vtx->Integral(3,-1)/h_noPU_prompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_prompt_not_vtx_matched_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_not_vtx_matched_EE_vtx->Integral(1,2) << " : " << h_noPU_prompt_not_vtx_matched_EE_vtx->Integral(3,-1) << "  " << "(" << 100 -h_noPU_prompt_not_vtx_matched_EE_vtx->Integral(3,-1)/h_noPU_prompt_not_vtx_matched_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EE_vtx->Integral() << endl;
  cout << endl;
    // nonprompt
      // Barrel
  TCanvas* c_st_noPU_nonprompt_EB_vtx_v2 = new TCanvas("c_st_noPU_nonprompt_EB_vtx_v2", "c_st_noPU_nonprompt_EB_vtx_v2", 1500, 1500);
  c_st_noPU_nonprompt_EB_vtx_v2->cd();
  c_st_noPU_nonprompt_EB_vtx_v2->SetLeftMargin(0.12);
  st_noPU_nonprompt_EB_vtx->Draw();
  leg_noPU_nonprompt_EB->Draw();
  c_st_noPU_nonprompt_EB_vtx_v2->Print("plots/ntuple/track_sigma_noPU_nonprompt_EB_vtx.pdf");
  cout << "noPU Barrel nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EB_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_PV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EB_vtx->Integral(3,-1)/h_noPU_nonprompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EB_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_SV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EB_vtx->Integral(3,-1)/h_noPU_nonprompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EB_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_PU_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EB_vtx->Integral(3,-1)/h_noPU_nonprompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EB_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_fake_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EB_vtx->Integral(3,-1)/h_noPU_nonprompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_nonprompt_not_vtx_matched_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_not_vtx_matched_EB_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_not_vtx_matched_EB_vtx->Integral(3,-1) << "  " << "(" << 100 -h_noPU_nonprompt_not_vtx_matched_EB_vtx->Integral(3,-1)/h_noPU_nonprompt_not_vtx_matched_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EB_vtx->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_noPU_nonprompt_EE_vtx_v2 = new TCanvas("c_st_noPU_nonprompt_EE_vtx_v2", "c_st_noPU_nonprompt_EE_vtx_v2", 1500, 1500);
  c_st_noPU_nonprompt_EE_vtx_v2->cd();
  c_st_noPU_nonprompt_EE_vtx_v2->SetLeftMargin(0.12);
  st_noPU_nonprompt_EE_vtx->Draw();
  leg_noPU_nonprompt_EE->Draw();
  c_st_noPU_nonprompt_EE_vtx_v2->Print("plots/ntuple/track_sigma_noPU_nonprompt_EE_vtx.pdf");
  cout << "noPU Endcap nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EE_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_PV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EE_vtx->Integral(3,-1)/h_noPU_nonprompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EE_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_SV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EE_vtx->Integral(3,-1)/h_noPU_nonprompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EE_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_PU_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EE_vtx->Integral(3,-1)/h_noPU_nonprompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EE_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_fake_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EE_vtx->Integral(3,-1)/h_noPU_nonprompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_nonprompt_not_vtx_matched_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_not_vtx_matched_EE_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_not_vtx_matched_EE_vtx->Integral(3,-1) << "  " << "(" << 100 -h_noPU_nonprompt_not_vtx_matched_EE_vtx->Integral(3,-1)/h_noPU_nonprompt_not_vtx_matched_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EE_vtx->Integral() << endl;
  cout << endl;

  cout << "///////////////////////" << endl;
  cout << "/////   3sigma   //////" << endl;
  cout << "///////////////////////" << endl;
  cout << "PU200 Barrel prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PV_EB_vtx->Integral(1,3) << " : " << h_PU200_prompt_PV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EB_vtx->Integral(4,-1)/h_PU200_prompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_SV_EB_vtx->Integral(1,3) << " : " << h_PU200_prompt_SV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EB_vtx->Integral(4,-1)/h_PU200_prompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PU_EB_vtx->Integral(1,3) << " : " << h_PU200_prompt_PU_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EB_vtx->Integral(4,-1)/h_PU200_prompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_fake_EB_vtx->Integral(1,3) << " : " << h_PU200_prompt_fake_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EB_vtx->Integral(4,-1)/h_PU200_prompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_prompt_not_vtx_matched_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_not_vtx_matched_EB_vtx->Integral(1,3) << " : " << h_PU200_prompt_not_vtx_matched_EB_vtx->Integral(4,-1) << "  " << "(" << 100 -h_PU200_prompt_not_vtx_matched_EB_vtx->Integral(4,-1)/h_PU200_prompt_not_vtx_matched_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EB_vtx->Integral() << endl;
  cout << endl;
      // Endcap
  cout << "PU200 Endcap prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PV_EE_vtx->Integral(1,3) << " : " << h_PU200_prompt_PV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EE_vtx->Integral(4,-1)/h_PU200_prompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_SV_EE_vtx->Integral(1,3) << " : " << h_PU200_prompt_SV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EE_vtx->Integral(4,-1)/h_PU200_prompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PU_EE_vtx->Integral(1,3) << " : " << h_PU200_prompt_PU_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EE_vtx->Integral(4,-1)/h_PU200_prompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_fake_EE_vtx->Integral(1,3) << " : " << h_PU200_prompt_fake_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EE_vtx->Integral(4,-1)/h_PU200_prompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_prompt_not_vtx_matched_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_not_vtx_matched_EE_vtx->Integral(1,3) << " : " << h_PU200_prompt_not_vtx_matched_EE_vtx->Integral(4,-1) << "  " << "(" << 100 -h_PU200_prompt_not_vtx_matched_EE_vtx->Integral(4,-1)/h_PU200_prompt_not_vtx_matched_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EE_vtx->Integral() << endl;
  cout << endl;
    // nonprompt
      // Barrel
  cout << "PU200 Barrel nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EB_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_PV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EB_vtx->Integral(4,-1)/h_PU200_nonprompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EB_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_SV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EB_vtx->Integral(4,-1)/h_PU200_nonprompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EB_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_PU_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EB_vtx->Integral(4,-1)/h_PU200_nonprompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EB_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_fake_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EB_vtx->Integral(4,-1)/h_PU200_nonprompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_nonprompt_not_vtx_matched_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_not_vtx_matched_EB_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_not_vtx_matched_EB_vtx->Integral(4,-1) << "  " << "(" << 100 -h_PU200_nonprompt_not_vtx_matched_EB_vtx->Integral(4,-1)/h_PU200_nonprompt_not_vtx_matched_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EB_vtx->Integral() << endl;
  cout << endl;
      // Endcap
  cout << "PU200 Endcap nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EE_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_PV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EE_vtx->Integral(4,-1)/h_PU200_nonprompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EE_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_SV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EE_vtx->Integral(4,-1)/h_PU200_nonprompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EE_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_PU_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EE_vtx->Integral(4,-1)/h_PU200_nonprompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EE_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_fake_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EE_vtx->Integral(4,-1)/h_PU200_nonprompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_PU200_nonprompt_not_vtx_matched_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_not_vtx_matched_EE_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_not_vtx_matched_EE_vtx->Integral(4,-1) << "  " << "(" << 100 -h_PU200_nonprompt_not_vtx_matched_EE_vtx->Integral(4,-1)/h_PU200_nonprompt_not_vtx_matched_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EE_vtx->Integral() << endl;
  cout << endl;

  // noPU
    // prompt
      // Barrel
  cout << "noPU Barrel prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PV_EB_vtx->Integral(1,3) << " : " << h_noPU_prompt_PV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EB_vtx->Integral(4,-1)/h_noPU_prompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_SV_EB_vtx->Integral(1,3) << " : " << h_noPU_prompt_SV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EB_vtx->Integral(4,-1)/h_noPU_prompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PU_EB_vtx->Integral(1,3) << " : " << h_noPU_prompt_PU_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EB_vtx->Integral(4,-1)/h_noPU_prompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_fake_EB_vtx->Integral(1,3) << " : " << h_noPU_prompt_fake_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EB_vtx->Integral(4,-1)/h_noPU_prompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_prompt_not_vtx_matched_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_not_vtx_matched_EB_vtx->Integral(1,3) << " : " << h_noPU_prompt_not_vtx_matched_EB_vtx->Integral(4,-1) << "  " << "(" << 100 -h_noPU_prompt_not_vtx_matched_EB_vtx->Integral(4,-1)/h_noPU_prompt_not_vtx_matched_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EB_vtx->Integral() << endl;
  cout << endl;
      // Endcap
  cout << "noPU Endcap prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PV_EE_vtx->Integral(1,3) << " : " << h_noPU_prompt_PV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EE_vtx->Integral(4,-1)/h_noPU_prompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_SV_EE_vtx->Integral(1,3) << " : " << h_noPU_prompt_SV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EE_vtx->Integral(4,-1)/h_noPU_prompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PU_EE_vtx->Integral(1,3) << " : " << h_noPU_prompt_PU_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EE_vtx->Integral(4,-1)/h_noPU_prompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_fake_EE_vtx->Integral(1,3) << " : " << h_noPU_prompt_fake_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EE_vtx->Integral(4,-1)/h_noPU_prompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_prompt_not_vtx_matched_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_not_vtx_matched_EE_vtx->Integral(1,3) << " : " << h_noPU_prompt_not_vtx_matched_EE_vtx->Integral(4,-1) << "  " << "(" << 100 -h_noPU_prompt_not_vtx_matched_EE_vtx->Integral(4,-1)/h_noPU_prompt_not_vtx_matched_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EE_vtx->Integral() << endl;
  cout << endl;
    // nonprompt
      // Barrel
  cout << "noPU Barrel nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EB_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_PV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EB_vtx->Integral(4,-1)/h_noPU_nonprompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EB_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_SV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EB_vtx->Integral(4,-1)/h_noPU_nonprompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EB_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_PU_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EB_vtx->Integral(4,-1)/h_noPU_nonprompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EB_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_fake_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EB_vtx->Integral(4,-1)/h_noPU_nonprompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_nonprompt_not_vtx_matched_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_not_vtx_matched_EB_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_not_vtx_matched_EB_vtx->Integral(4,-1) << "  " << "(" << 100 -h_noPU_nonprompt_not_vtx_matched_EB_vtx->Integral(4,-1)/h_noPU_nonprompt_not_vtx_matched_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EB_vtx->Integral() << endl;
  cout << endl;
      // Endcap
  cout << "noPU Endcap nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EE_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_PV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EE_vtx->Integral(4,-1)/h_noPU_nonprompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EE_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_SV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EE_vtx->Integral(4,-1)/h_noPU_nonprompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EE_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_PU_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EE_vtx->Integral(4,-1)/h_noPU_nonprompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EE_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_fake_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EE_vtx->Integral(4,-1)/h_noPU_nonprompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[not vtx] " << h_noPU_nonprompt_not_vtx_matched_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_not_vtx_matched_EE_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_not_vtx_matched_EE_vtx->Integral(4,-1) << "  " << "(" << 100 -h_noPU_nonprompt_not_vtx_matched_EE_vtx->Integral(4,-1)/h_noPU_nonprompt_not_vtx_matched_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EE_vtx->Integral() << endl;
  cout << endl;


  TCanvas* c_PV_dtsig_PU200_prompt_EB_vtx = new TCanvas("c_PV_dtsig_PU200_prompt_EB_vtx", "c_PV_dtsig_PU200_prompt_EB_vtx", 1500, 1500);
  c_PV_dtsig_PU200_prompt_EB_vtx->cd();
  c_PV_dtsig_PU200_prompt_EB_vtx->SetLeftMargin(0.12);
  h_PU200_prompt_PV_EB_vtx->SetTitle(";#sigma_{t} (PV, track);Counts");
  h_PU200_prompt_PV_EB_vtx->GetYaxis()->SetRangeUser(0.,600.);
  h_PU200_prompt_PV_EB_vtx->Draw();
  leg_PV_dtsig_PU200_prompt_EB_F->Draw();
  c_PV_dtsig_PU200_prompt_EB_vtx->Print("plots/ntuple/PV_dtsig_PU200_prompt_EB_vtx.pdf");

  TCanvas* c_PV_bx_PU200_prompt_EB_vtx = new TCanvas("c_PV_bx_PU200_prompt_EB_vtx", "c_PV_bx_PU200_prompt_EB_vtx", 1500, 1500);
  c_PV_bx_PU200_prompt_EB_vtx->cd();
  c_PV_bx_PU200_prompt_EB_vtx->SetLeftMargin(0.12);
  h_PU200_prompt_bx_PV_EB_vtx->SetTitle(";Bunch crossing number of track from PV;Counts");
  h_PU200_prompt_bx_PV_EB_vtx->SetLineColor(kBlack);
  h_PU200_prompt_bx_PV_EB_vtx->SetLineWidth(2);
  h_PU200_prompt_bx_PV_EB_vtx->Draw();
  leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_bx_PU200_prompt_EB_vtx->Print("plots/ntuple/PV_bx_PU200_prompt_EB_vtx.pdf");
  cout << "bx = 0  : " << h_PU200_prompt_bx_PV_EB_vtx->GetBinContent(1) << endl;
  cout << "total bx: " << h_PU200_prompt_bx_PV_EB_vtx->Integral(0,-1) << endl;
  cout << endl;

  TCanvas* c_PV_evtId_PU200_prompt_EB_vtx = new TCanvas("c_PV_evtId_PU200_prompt_EB_vtx", "c_PV_evtId_PU200_prompt_EB_vtx", 1500, 1500);
  c_PV_evtId_PU200_prompt_EB_vtx->cd();
  c_PV_evtId_PU200_prompt_EB_vtx->SetLeftMargin(0.12);
  h_PU200_prompt_evtId_PV_EB_vtx->SetTitle(";Event ID of track from PV;Counts");
  h_PU200_prompt_evtId_PV_EB_vtx->SetLineColor(kBlack);
  h_PU200_prompt_evtId_PV_EB_vtx->SetLineWidth(2);
  h_PU200_prompt_evtId_PV_EB_vtx->Draw();
  c_PV_evtId_PU200_prompt_EB_vtx->SetLogy();
  leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_evtId_PU200_prompt_EB_vtx->Print("plots/ntuple/PV_evtId_PU200_prompt_EB_vtx.pdf");
  cout << "evtId = 0  : " << h_PU200_prompt_evtId_PV_EB_vtx->GetBinContent(1) << endl;
  cout << "total evtId: " << h_PU200_prompt_evtId_PV_EB_vtx->Integral(0,-1) << endl;
  cout << endl;

  TCanvas* c_PV_PVweight_PU200_prompt_EB_vtx = new TCanvas("c_PV_PVweight_PU200_prompt_EB_vtx", "c_PV_PVweight_PU200_prompt_EB_vtx", 1500, 1500);
  c_PV_PVweight_PU200_prompt_EB_vtx->cd();
  c_PV_PVweight_PU200_prompt_EB_vtx->SetLeftMargin(0.12);
  h_PU200_prompt_PVweight_PV_EB_vtx->SetTitle(";PV weight of track from PV;Counts");
  h_PU200_prompt_PVweight_PV_EB_vtx->SetLineColor(kBlack);
  h_PU200_prompt_PVweight_PV_EB_vtx->SetLineWidth(2);
  h_PU200_prompt_PVweight_PV_EB_vtx->Draw();
  leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_PVweight_PU200_prompt_EB_vtx->Print("plots/ntuple/PV_PVweight_PU200_prompt_EB_vtx.pdf");
  cout << "PV weight = 0  : " << h_PU200_prompt_PVweight_PV_EB_vtx->GetBinContent(1) << endl;
  cout << "less than 0.5  : " << h_PU200_prompt_PVweight_PV_EB_vtx->Integral(0,50) << endl;
  cout << "total PV weight: " << h_PU200_prompt_PVweight_PV_EB_vtx->Integral(0,-1) << endl;
  cout << endl;

  TCanvas* c_PV_simvtx_bx_PU200_prompt_EB_vtx = new TCanvas("c_PV_simvtx_bx_PU200_prompt_EB_vtx", "c_PV_simvtx_bx_PU200_prompt_EB_vtx", 1500, 1500);
  c_PV_simvtx_bx_PU200_prompt_EB_vtx->cd();
  c_PV_simvtx_bx_PU200_prompt_EB_vtx->SetLeftMargin(0.12);
  h_PU200_prompt_simvtx_bx_PV_EB_vtx->SetTitle(";Bunch crossing number of sim vertex;Counts");
  h_PU200_prompt_simvtx_bx_PV_EB_vtx->SetLineColor(kBlack);
  h_PU200_prompt_simvtx_bx_PV_EB_vtx->SetLineWidth(2);
  h_PU200_prompt_simvtx_bx_PV_EB_vtx->Draw();
  leg_simvtx_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_simvtx_bx_PU200_prompt_EB_vtx->Print("plots/ntuple/PV_simvtx_bx_PU200_prompt_EB_vtx.pdf");
  cout << "simvtx bx = 0  : " << h_PU200_prompt_simvtx_bx_PV_EB_vtx->GetBinContent(1) << endl;
  cout << "total simvtx bx: " << h_PU200_prompt_simvtx_bx_PV_EB_vtx->Integral(0,-1) << endl;
  cout << endl;


  TCanvas* c_PV_simvtx_evtId_PU200_prompt_EB_vtx = new TCanvas("c_PV_simvtx_evtId_PU200_prompt_EB_vtx", "c_PV_simvtx_evtId_PU200_prompt_EB_vtx", 1500, 1500);
  c_PV_simvtx_evtId_PU200_prompt_EB_vtx->cd();
  c_PV_simvtx_evtId_PU200_prompt_EB_vtx->SetLeftMargin(0.12);
  h_PU200_prompt_simvtx_evtId_PV_EB_vtx->SetTitle(";Event ID of sim vertex;Counts");
  h_PU200_prompt_simvtx_evtId_PV_EB_vtx->SetLineColor(kBlack);
  h_PU200_prompt_simvtx_evtId_PV_EB_vtx->SetLineWidth(2);
  h_PU200_prompt_simvtx_evtId_PV_EB_vtx->Draw();
  c_PV_simvtx_evtId_PU200_prompt_EB_vtx->SetLogy();
  leg_simvtx_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_simvtx_evtId_PU200_prompt_EB_vtx->Print("plots/ntuple/PV_simvtx_evtId_PU200_prompt_EB_vtx.pdf");
  cout << "simvtx evtId = 0  : " << h_PU200_prompt_simvtx_evtId_PV_EB_vtx->GetBinContent(1) << endl;
  cout << "total simvtx evtId: " << h_PU200_prompt_simvtx_evtId_PV_EB_vtx->Integral(0,-1) << endl;
  cout << endl;


  // nonprompt
  // Barrel
  TCanvas* c_PV_dtsig_PU200_nonprompt_EB_vtx = new TCanvas("c_PV_dtsig_PU200_nonprompt_EB_vtx", "c_PV_dtsig_PU200_nonprompt_EB_vtx", 1500, 1500);
  c_PV_dtsig_PU200_nonprompt_EB_vtx->cd();
  c_PV_dtsig_PU200_nonprompt_EB_vtx->SetLeftMargin(0.12);
  h_PU200_nonprompt_PV_EB_vtx->SetTitle(";#sigma_{t} (PV, muon);Counts");
  h_PU200_nonprompt_PV_EB_vtx->GetYaxis()->SetRangeUser(0.,1000.);
  h_PU200_nonprompt_PV_EB_vtx->Draw();
  leg_PV_dtsig_PU200_prompt_EB_F->Draw();
  c_PV_dtsig_PU200_nonprompt_EB_vtx->Print("plots/ntuple/PV_dtsig_PU200_nonprompt_EB_vtx.pdf");

  TCanvas* c_PV_bx_PU200_nonprompt_EB_vtx = new TCanvas("c_PV_bx_PU200_nonprompt_EB_vtx", "c_PV_bx_PU200_nonprompt_EB_vtx", 1500, 1500);
  c_PV_bx_PU200_nonprompt_EB_vtx->cd();
  c_PV_bx_PU200_nonprompt_EB_vtx->SetLeftMargin(0.12);
  h_PU200_nonprompt_bx_PV_EB_vtx->SetTitle(";Bunch crossing number of track from PV;Counts");
  h_PU200_nonprompt_bx_PV_EB_vtx->SetLineColor(kBlack);
  h_PU200_nonprompt_bx_PV_EB_vtx->SetLineWidth(2);
  h_PU200_nonprompt_bx_PV_EB_vtx->Draw();
  leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_bx_PU200_nonprompt_EB_vtx->Print("plots/ntuple/PV_bx_PU200_nonprompt_EB_vtx.pdf");
  cout << "bx = 0  : " << h_PU200_nonprompt_bx_PV_EB_vtx->GetBinContent(1) << endl;
  cout << "total bx: " << h_PU200_nonprompt_bx_PV_EB_vtx->Integral(0,-1) << endl;
  cout << endl;

  TCanvas* c_PV_evtId_PU200_nonprompt_EB_vtx = new TCanvas("c_PV_evtId_PU200_nonprompt_EB_vtx", "c_PV_evtId_PU200_nonprompt_EB_vtx", 1500, 1500);
  c_PV_evtId_PU200_nonprompt_EB_vtx->cd();
  c_PV_evtId_PU200_nonprompt_EB_vtx->SetLeftMargin(0.12);
  h_PU200_nonprompt_evtId_PV_EB_vtx->SetTitle(";Event ID of track from PV;Counts");
  h_PU200_nonprompt_evtId_PV_EB_vtx->SetLineColor(kBlack);
  h_PU200_nonprompt_evtId_PV_EB_vtx->SetLineWidth(2);
  h_PU200_nonprompt_evtId_PV_EB_vtx->Draw();
  c_PV_evtId_PU200_nonprompt_EB_vtx->SetLogy();
  leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_evtId_PU200_nonprompt_EB_vtx->Print("plots/ntuple/PV_evtId_PU200_nonprompt_EB_vtx.pdf");
  cout << "evtId = 0  : " << h_PU200_nonprompt_evtId_PV_EB_vtx->GetBinContent(1) << endl;
  cout << "total evtId: " << h_PU200_nonprompt_evtId_PV_EB_vtx->Integral(0,-1) << endl;
  cout << endl;

  TCanvas* c_PV_PVweight_PU200_nonprompt_EB_vtx = new TCanvas("c_PV_PVweight_PU200_nonprompt_EB_vtx", "c_PV_PVweight_PU200_nonprompt_EB_vtx", 1500, 1500);
  c_PV_PVweight_PU200_nonprompt_EB_vtx->cd();
  c_PV_PVweight_PU200_nonprompt_EB_vtx->SetLeftMargin(0.12);
  h_PU200_nonprompt_PVweight_PV_EB_vtx->SetTitle(";PV weight of track from PV;Counts");
  h_PU200_nonprompt_PVweight_PV_EB_vtx->SetLineColor(kBlack);
  h_PU200_nonprompt_PVweight_PV_EB_vtx->SetLineWidth(2);
  h_PU200_nonprompt_PVweight_PV_EB_vtx->Draw();
  leg_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_PVweight_PU200_nonprompt_EB_vtx->Print("plots/ntuple/PV_PVweight_PU200_nonprompt_EB_vtx.pdf");
  cout << "PV weight = 0  : " << h_PU200_nonprompt_PVweight_PV_EB_vtx->GetBinContent(1) << endl;
  cout << "less than 0.5  : " << h_PU200_nonprompt_PVweight_PV_EB_vtx->Integral(0,50) << endl;
  cout << "total PV weight: " << h_PU200_nonprompt_PVweight_PV_EB_vtx->Integral(0,-1) << endl;
  cout << endl;

  TCanvas* c_PV_simvtx_bx_PU200_nonprompt_EB_vtx = new TCanvas("c_PV_simvtx_bx_PU200_nonprompt_EB_vtx", "c_PV_simvtx_bx_PU200_nonprompt_EB_vtx", 1500, 1500);
  c_PV_simvtx_bx_PU200_nonprompt_EB_vtx->cd();
  c_PV_simvtx_bx_PU200_nonprompt_EB_vtx->SetLeftMargin(0.12);
  h_PU200_nonprompt_simvtx_bx_PV_EB_vtx->SetTitle(";Bunch crossing number of sim vertex;Counts");
  h_PU200_nonprompt_simvtx_bx_PV_EB_vtx->SetLineColor(kBlack);
  h_PU200_nonprompt_simvtx_bx_PV_EB_vtx->SetLineWidth(2);
  h_PU200_nonprompt_simvtx_bx_PV_EB_vtx->Draw();
  leg_simvtx_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_simvtx_bx_PU200_nonprompt_EB_vtx->Print("plots/ntuple/PV_simvtx_bx_PU200_nonprompt_EB_vtx.pdf");
  cout << "simvtx bx = 0  : " << h_PU200_nonprompt_simvtx_bx_PV_EB_vtx->GetBinContent(1) << endl;
  cout << "total simvtx bx: " << h_PU200_nonprompt_simvtx_bx_PV_EB_vtx->Integral(0,-1) << endl;
  cout << endl;


  TCanvas* c_PV_simvtx_evtId_PU200_nonprompt_EB_vtx = new TCanvas("c_PV_simvtx_evtId_PU200_nonprompt_EB_vtx", "c_PV_simvtx_evtId_PU200_nonprompt_EB_vtx", 1500, 1500);
  c_PV_simvtx_evtId_PU200_nonprompt_EB_vtx->cd();
  c_PV_simvtx_evtId_PU200_nonprompt_EB_vtx->SetLeftMargin(0.12);
  h_PU200_nonprompt_simvtx_evtId_PV_EB_vtx->SetTitle(";Event ID of sim vertex;Counts");
  h_PU200_nonprompt_simvtx_evtId_PV_EB_vtx->SetLineColor(kBlack);
  h_PU200_nonprompt_simvtx_evtId_PV_EB_vtx->SetLineWidth(2);
  h_PU200_nonprompt_simvtx_evtId_PV_EB_vtx->Draw();
  c_PV_simvtx_evtId_PU200_nonprompt_EB_vtx->SetLogy();
  leg_simvtx_PV_dtsig_PU200_prompt_EB_L->Draw();
  c_PV_simvtx_evtId_PU200_nonprompt_EB_vtx->Print("plots/ntuple/PV_simvtx_evtId_PU200_nonprompt_EB_vtx.pdf");
  cout << "simvtx evtId = 0  : " << h_PU200_nonprompt_simvtx_evtId_PV_EB_vtx->GetBinContent(1) << endl;
  cout << "total simvtx evtId: " << h_PU200_nonprompt_simvtx_evtId_PV_EB_vtx->Integral(0,-1) << endl;
  cout << endl;



}



int main(int argc, char **argv)
{
  
  draw_track_type_sigma_ntuple();

  return 0;
}



















































