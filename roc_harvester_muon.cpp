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

void N_muon();
void draw_iso_efficiency();       // draw iso. efficiency and normalization as a function of rel. iso.
void draw_iso_distribution();       // draw iso. efficiency and normalization as a function of rel. iso.
void draw_reliso_roc();           // draw ROC curves as a function of rel. iso.
void draw_pt_roc();               // draw ROC curves as a function of pT sum of charged tracks
void pTeff();                                 // draw efficiency plots as a function of pT of muons
void N_genMatched();                          // draw plots to check whether muon's track is matched with a Gen particle
void track_type_v1();
void track_sigma_type_v1();
void track_type_v2();
void track_sigma_type_v2();
void muon_mother_pdgId();
void fraction_mva_cut();
void evtId_pvtrk();
void pt_eta_phi();
string pdgId(int pdgid);


TString path_PU200_prompt        = "240610/harvester_PU200_prompt_muon.root";
TString path_noPU_prompt         = "240610/harvester_noPU_prompt_muon.root";
TString path_PU200_prompt_vtx    = "240610/harvester_PU200_prompt_muon_vtx.root";
TString path_noPU_prompt_vtx     = "240610/harvester_noPU_prompt_muon_vtx.root";
TString path_PU200_nonprompt     = "240610/harvester_PU200_nonprompt_muon_qcd.root";
TString path_noPU_nonprompt      = "240610/harvester_noPU_nonprompt_muon_qcd.root";
TString path_PU200_nonprompt_vtx = "240610/harvester_PU200_nonprompt_muon_qcd_vtx.root";
TString path_noPU_nonprompt_vtx  = "240610/harvester_noPU_nonprompt_muon_qcd_vtx.root";
//TString path_PU200_nonprompt     = "240610/harvester_PU200_nonprompt_muon_ttbar.root";
//TString path_noPU_nonprompt      = "240610/harvester_noPU_nonprompt_muon_ttbar.root";
//TString path_PU200_nonprompt_vtx = "240610/harvester_PU200_nonprompt_muon_ttbar_vtx.root";
//TString path_noPU_nonprompt_vtx  = "240610/harvester_noPU_nonprompt_muon_ttbar_vtx.root";

TString ntuple_PU200_prompt        = "240610/ntuple_PU200_prompt_muon.root";
TString ntuple_noPU_prompt         = "240610/ntuple_noPU_prompt_muon.root";
TString ntuple_PU200_prompt_vtx    = "240610/ntuple_PU200_prompt_muon_vtx.root";
TString ntuple_noPU_prompt_vtx     = "240610/ntuple_noPU_prompt_muon_vtx.root";
TString ntuple_PU200_nonprompt     = "240610/ntuple_PU200_nonprompt_muon_qcd.root";
TString ntuple_noPU_nonprompt      = "240610/ntuple_noPU_nonprompt_muon_qcd.root";
TString ntuple_PU200_nonprompt_vtx = "240610/ntuple_PU200_nonprompt_muon_qcd_vtx.root";
TString ntuple_noPU_nonprompt_vtx  = "240610/ntuple_noPU_nonprompt_muon_qcd_vtx.root";
//TString ntuple_PU200_nonprompt     = "240610/ntuple_PU200_nonprompt_muon_ttbar.root";
//TString ntuple_noPU_nonprompt      = "240610/ntuple_noPU_nonprompt_muon_ttbar.root";
//TString ntuple_PU200_nonprompt_vtx = "240610/ntuple_PU200_nonprompt_muon_ttbar_vtx.root";
//TString ntuple_noPU_nonprompt_vtx  = "240610/ntuple_noPU_nonprompt_muon_ttbar_vtx.root";


/*
TString path_PU200_prompt        = "240522/harvester_PU200_prompt_muon.root";
TString path_noPU_prompt         = "240522/harvester_noPU_prompt_muon.root";
TString path_PU200_prompt_vtx    = "240522/harvester_PU200_prompt_muon_vtx.root";
TString path_noPU_prompt_vtx     = "240522/harvester_noPU_prompt_muon_vtx.root";
TString path_PU200_nonprompt     = "240522/harvester_PU200_nonprompt_muon_qcd.root";
TString path_noPU_nonprompt      = "240522/harvester_noPU_nonprompt_muon_qcd.root";
TString path_PU200_nonprompt_vtx = "240522/harvester_PU200_nonprompt_muon_qcd_vtx.root";
TString path_noPU_nonprompt_vtx  = "240522/harvester_noPU_nonprompt_muon_qcd_vtx.root";
//TString path_PU200_nonprompt     = "240522/harvester_PU200_nonprompt_muon_ttbar.root";
//TString path_noPU_nonprompt      = "240522/harvester_noPU_nonprompt_muon_ttbar.root";
//TString path_PU200_nonprompt_vtx = "240522/harvester_PU200_nonprompt_muon_ttbar_vtx.root";
//TString path_noPU_nonprompt_vtx  = "240522/harvester_noPU_nonprompt_muon_ttbar_vtx.root";

TString ntuple_PU200_prompt        = "240523/ntuple_PU200_prompt_muon.root";
TString ntuple_noPU_prompt         = "240523/ntuple_noPU_prompt_muon.root";
TString ntuple_PU200_prompt_vtx    = "240522/ntuple_PU200_prompt_muon_vtx.root";
TString ntuple_noPU_prompt_vtx     = "240522/ntuple_noPU_prompt_muon_vtx.root";
TString ntuple_PU200_nonprompt     = "240523/ntuple_PU200_nonprompt_muon_qcd.root";
TString ntuple_noPU_nonprompt      = "240523/ntuple_noPU_nonprompt_muon_qcd.root";
TString ntuple_PU200_nonprompt_vtx = "240522/ntuple_PU200_nonprompt_muon_qcd_vtx.root";
TString ntuple_noPU_nonprompt_vtx  = "240522/ntuple_noPU_nonprompt_muon_qcd_vtx.root";
//TString ntuple_PU200_nonprompt     = "240523/ntuple_PU200_nonprompt_muon_ttbar.root";
//TString ntuple_noPU_nonprompt      = "240523/ntuple_noPU_nonprompt_muon_ttbar.root";
//TString ntuple_PU200_nonprompt_vtx = "240522/ntuple_PU200_nonprompt_muon_ttbar_vtx.root";
//TString ntuple_noPU_nonprompt_vtx  = "240522/ntuple_noPU_nonprompt_muon_ttbar_vtx.root";
*/


string pdgId(int pdgid) {
  if (pdgid == 1) return "d";
  if (pdgid == -1) return "d_bar";
  else if (pdgid == 2) return "u";
  else if (pdgid == -2) return "u_bar";
  else if (pdgid == 3) return "s";
  else if (pdgid == -3) return "s_bar";
  else if (pdgid == 4) return "c";
  else if (pdgid == -4) return "c_bar";
  else if (pdgid == 5) return "b";
  else if (pdgid == -5) return "b_bar";
  else if (pdgid == 6) return "t";
  else if (pdgid == -6) return "t_bar";
  else if (pdgid == 11) return "e-";
  else if (pdgid == -11) return "e+";
  else if (pdgid == 13) return "mu-";
  else if (pdgid == -13) return "mu+";
  else if (pdgid == 15) return "tau-";
  else if (pdgid == -15) return "tau+";
  else if (pdgid == 21) return "g";
  else if (pdgid == 22) return "gamma";
  else if (pdgid == 23) return "Z";
  else if (pdgid == 24) return "W+";
  else if (pdgid == -24) return "W-";
  else if (pdgid == 25) return "H";
  else if (pdgid == 111) return "pion";
  else if (pdgid == 211) return "pion+";
  else if (pdgid == -211) return "pion-";
  else if (pdgid == 113) return "rho";
  else if (pdgid == 213) return "rho+";
  else if (pdgid == -213) return "rho-";
  else if (pdgid == 130) return "K_L^0";
  else if (pdgid == 221) return "eta";
  else if (pdgid == 223) return "omega";
  else if (pdgid == 313) return "K*^0";
  else if (pdgid == -313) return "K*^0";
  else if (pdgid == -321) return "K^-";
  else if (pdgid == 321) return "K^+";
  else if (pdgid == -323) return "K*^-";
  else if (pdgid == 323) return "K*^+";
  else if (pdgid == 331) return "eta\'";
  else if (pdgid == 333) return "phi";
  else if (pdgid == -411) return "D^-";
  else if (pdgid == 411) return "D^+";
  else if (pdgid == -421) return "D^0";
  else if (pdgid == 421) return "D^0";
  else if (pdgid == -431) return "D_s^-";
  else if (pdgid == 431) return "D_s^+";
  else if (pdgid == 443) return "J/psi(1S)";
  else if (pdgid == -511) return "B0";
  else if (pdgid == 511) return "B0";
  else if (pdgid == -521) return "B-";
  else if (pdgid == 521) return "B+";
  else if (pdgid == -531) return "B_s^0";
  else if (pdgid == 531) return "B_s^0";
  else if (pdgid == -5232) return "Xi baryon (Xi_b^0)";
  else if (pdgid == 5232) return "Xi baryon (Xi_b^0)";
  else if (pdgid == -5132) return "Xi baryon (Xi_b^+)";
  else if (pdgid == 5132) return "Xi baryon (Xi_b^-)";
  else if (pdgid == -5122) return "Lambda baryon (Lambda_b^0)";
  else if (pdgid == 5122) return "Lambda baryon (Lambda_b^0)";
  else if (pdgid == -4122) return "Lambda baryon (Lambda_c^-)";
  else if (pdgid == 4122) return "Lambda baryon (Lambda_c^+)";
  else if (pdgid == 4132) return "Xi baryon (Xi_c^0)";
  else return Form("non-identified, %d", pdgid);
}

void N_muon() {
  
  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  TFile *f_PU200_prompt_vtx, *f_PU200_nonprompt_vtx, *f_noPU_prompt_vtx, *f_noPU_nonprompt_vtx;
  // PU200
  f_PU200_prompt = new TFile(Form("data/%s", path_PU200_prompt.Data()));
  f_PU200_nonprompt = new TFile(Form("data/%s", path_PU200_nonprompt.Data()));
  f_PU200_prompt_vtx = new TFile(Form("data/%s", path_PU200_prompt_vtx.Data()));
  f_PU200_nonprompt_vtx = new TFile(Form("data/%s", path_PU200_nonprompt_vtx.Data()));
  // noPU
  f_noPU_prompt = new TFile(Form("data/%s", path_noPU_prompt.Data()));
  f_noPU_nonprompt = new TFile(Form("data/%s", path_noPU_nonprompt.Data()));
  f_noPU_prompt_vtx = new TFile(Form("data/%s", path_noPU_prompt_vtx.Data()));
  f_noPU_nonprompt_vtx = new TFile(Form("data/%s", path_noPU_nonprompt_vtx.Data()));


  TString dir="/DQMData/Run 1/MTD/Run summary/MuonIso/";

  TH1D *h_PU200_prompt_EB, *h_PU200_nonprompt_EB, *h_noPU_prompt_EB, *h_noPU_nonprompt_EB;
  TH1D *h_PU200_prompt_EE, *h_PU200_nonprompt_EE, *h_noPU_prompt_EE, *h_noPU_nonprompt_EE;

  // PU200
  h_PU200_prompt_EB 	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_Nmuons_Sig_EB");
  h_PU200_prompt_EE 	   = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_Nmuons_Sig_EE");
  h_PU200_nonprompt_EB 	   = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_Nmuons_Bkg_EB");
  h_PU200_nonprompt_EE 	   = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_Nmuons_Bkg_EE");
  // noPU
  h_noPU_prompt_EB 	   	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_Nmuons_Sig_EB");
  h_noPU_prompt_EE 	       = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_Nmuons_Sig_EE");
  h_noPU_nonprompt_EB 	   = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_Nmuons_Bkg_EB");
  h_noPU_nonprompt_EE 	   = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_Nmuons_Bkg_EE");

  int nmus_PU200_prompt_EB=0, nmus_PU200_prompt_EE=0, nmus_PU200_nonprompt_EB=0, nmus_PU200_nonprompt_EE=0;
  int nmus_noPU_prompt_EB=0, nmus_noPU_prompt_EE=0, nmus_noPU_nonprompt_EB=0, nmus_noPU_nonprompt_EE=0;

  for(int i=1; i<10; i++) {
    nmus_PU200_prompt_EB    = nmus_PU200_prompt_EB    + h_PU200_prompt_EB->GetBinContent(i)*(i-1);
    nmus_PU200_prompt_EE    = nmus_PU200_prompt_EE    + h_PU200_prompt_EE->GetBinContent(i)*(i-1);
    nmus_PU200_nonprompt_EB = nmus_PU200_nonprompt_EB + h_PU200_nonprompt_EB->GetBinContent(i)*(i-1);
    nmus_PU200_nonprompt_EE = nmus_PU200_nonprompt_EE + h_PU200_nonprompt_EE->GetBinContent(i)*(i-1);
    nmus_noPU_prompt_EB     = nmus_noPU_prompt_EB     + h_noPU_prompt_EB->GetBinContent(i)*(i-1);
    nmus_noPU_prompt_EE     = nmus_noPU_prompt_EE     + h_noPU_prompt_EE->GetBinContent(i)*(i-1);
    nmus_noPU_nonprompt_EB  = nmus_noPU_nonprompt_EB  + h_noPU_nonprompt_EB->GetBinContent(i)*(i-1);
    nmus_noPU_nonprompt_EE  = nmus_noPU_nonprompt_EE  + h_noPU_nonprompt_EE->GetBinContent(i)*(i-1);
  }

  cout << "total muon of PU200 prompt in Barrel   : " << nmus_PU200_prompt_EB << endl;
  cout << "total muon of PU200 prompt in Endcap   : " << nmus_PU200_prompt_EE << endl;
  cout << "total muon of PU200 nonprompt in Barrel: " << nmus_PU200_nonprompt_EB << endl;
  cout << "total muon of PU200 nonprompt in Endcap: " << nmus_PU200_nonprompt_EE << endl;
  cout << "total muon of noPU prompt in Barrel    : " << nmus_noPU_prompt_EB << endl;
  cout << "total muon of noPU prompt in Endcap    : " << nmus_noPU_prompt_EE << endl;
  cout << "total muon of noPU nonprompt in Barrel : " << nmus_noPU_nonprompt_EB << endl;
  cout << "total muon of noPU nonprompt in Endcap : " << nmus_noPU_nonprompt_EE << endl;
  cout << endl;

}

void draw_iso_efficiency() {
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  TFile *f_PU200_prompt_vtx, *f_PU200_nonprompt_vtx, *f_noPU_prompt_vtx, *f_noPU_nonprompt_vtx;
  // PU200
  f_PU200_prompt = new TFile(Form("data/%s", path_PU200_prompt.Data()));
  f_PU200_nonprompt = new TFile(Form("data/%s", path_PU200_nonprompt.Data()));
  f_PU200_prompt_vtx = new TFile(Form("data/%s", path_PU200_prompt_vtx.Data()));
  f_PU200_nonprompt_vtx = new TFile(Form("data/%s", path_PU200_nonprompt_vtx.Data()));
  // noPU
  f_noPU_prompt = new TFile(Form("data/%s", path_noPU_prompt.Data()));
  f_noPU_nonprompt = new TFile(Form("data/%s", path_noPU_nonprompt.Data()));
  f_noPU_prompt_vtx = new TFile(Form("data/%s", path_noPU_prompt_vtx.Data()));
  f_noPU_nonprompt_vtx = new TFile(Form("data/%s", path_noPU_nonprompt_vtx.Data()));

  // Histograms for Barrel region
    // PU200
  TH1D *h_PU200_prompt_EB, *h_PU200_nonprompt_EB, *h_noPU_prompt_EB, *h_noPU_nonprompt_EB;
  TH1D *h_PU200_prompt_2sigma_EB, *h_PU200_nonprompt_2sigma_EB, *h_noPU_prompt_2sigma_EB, *h_noPU_nonprompt_2sigma_EB;
  TH1D *h_PU200_prompt_3sigma_EB, *h_PU200_nonprompt_3sigma_EB, *h_noPU_prompt_3sigma_EB, *h_noPU_nonprompt_3sigma_EB;
  TH1D *h_PU200_prompt_4sigma_EB, *h_PU200_nonprompt_4sigma_EB, *h_noPU_prompt_4sigma_EB, *h_noPU_nonprompt_4sigma_EB;
  TH1D *h_PU200_prompt_40_EB, *h_PU200_nonprompt_40_EB, *h_noPU_prompt_40_EB, *h_noPU_nonprompt_40_EB;
  TH1D *h_PU200_prompt_60_EB, *h_PU200_nonprompt_60_EB, *h_noPU_prompt_60_EB, *h_noPU_nonprompt_60_EB;
  TH1D *h_PU200_prompt_80_EB, *h_PU200_nonprompt_80_EB, *h_noPU_prompt_80_EB, *h_noPU_nonprompt_80_EB;
  TH1D *h_PU200_prompt_100_EB, *h_PU200_nonprompt_100_EB, *h_noPU_prompt_100_EB, *h_noPU_nonprompt_100_EB;
      // GEN case
  TH1D *h_PU200_prompt_gen_EB, *h_PU200_nonprompt_gen_EB, *h_noPU_prompt_gen_EB, *h_noPU_nonprompt_gen_EB;
  TH1D *h_PU200_prompt_gen_2sigma_EB, *h_PU200_nonprompt_gen_2sigma_EB, *h_noPU_prompt_gen_2sigma_EB, *h_noPU_nonprompt_gen_2sigma_EB;
  TH1D *h_PU200_prompt_gen_3sigma_EB, *h_PU200_nonprompt_gen_3sigma_EB, *h_noPU_prompt_gen_3sigma_EB, *h_noPU_nonprompt_gen_3sigma_EB;
  TH1D *h_PU200_prompt_gen_4sigma_EB, *h_PU200_nonprompt_gen_4sigma_EB, *h_noPU_prompt_gen_4sigma_EB, *h_noPU_nonprompt_gen_4sigma_EB;
      // vtx
  TH1D *h_PU200_prompt_2sigma_EB_vtx, *h_PU200_nonprompt_2sigma_EB_vtx, *h_noPU_prompt_2sigma_EB_vtx, *h_noPU_nonprompt_2sigma_EB_vtx;
  TH1D *h_PU200_prompt_3sigma_EB_vtx, *h_PU200_nonprompt_3sigma_EB_vtx, *h_noPU_prompt_3sigma_EB_vtx, *h_noPU_nonprompt_3sigma_EB_vtx;
  TH1D *h_PU200_prompt_4sigma_EB_vtx, *h_PU200_nonprompt_4sigma_EB_vtx, *h_noPU_prompt_4sigma_EB_vtx, *h_noPU_nonprompt_4sigma_EB_vtx;

  // Histograms for Endcap region
    // PU200
  TH1D *h_PU200_prompt_EE, *h_PU200_nonprompt_EE, *h_noPU_prompt_EE, *h_noPU_nonprompt_EE;
  TH1D *h_PU200_prompt_2sigma_EE, *h_PU200_nonprompt_2sigma_EE, *h_noPU_prompt_2sigma_EE, *h_noPU_nonprompt_2sigma_EE;
  TH1D *h_PU200_prompt_3sigma_EE, *h_PU200_nonprompt_3sigma_EE, *h_noPU_prompt_3sigma_EE, *h_noPU_nonprompt_3sigma_EE;
  TH1D *h_PU200_prompt_4sigma_EE, *h_PU200_nonprompt_4sigma_EE, *h_noPU_prompt_4sigma_EE, *h_noPU_nonprompt_4sigma_EE;
  TH1D *h_PU200_prompt_40_EE, *h_PU200_nonprompt_40_EE, *h_noPU_prompt_40_EE, *h_noPU_nonprompt_40_EE;
  TH1D *h_PU200_prompt_60_EE, *h_PU200_nonprompt_60_EE, *h_noPU_prompt_60_EE, *h_noPU_nonprompt_60_EE;
  TH1D *h_PU200_prompt_80_EE, *h_PU200_nonprompt_80_EE, *h_noPU_prompt_80_EE, *h_noPU_nonprompt_80_EE;
  TH1D *h_PU200_prompt_100_EE, *h_PU200_nonprompt_100_EE, *h_noPU_prompt_100_EE, *h_noPU_nonprompt_100_EE;
      // GEN case
  TH1D *h_PU200_prompt_gen_EE, *h_PU200_nonprompt_gen_EE, *h_noPU_prompt_gen_EE, *h_noPU_nonprompt_gen_EE;
  TH1D *h_PU200_prompt_gen_2sigma_EE, *h_PU200_nonprompt_gen_2sigma_EE, *h_noPU_prompt_gen_2sigma_EE, *h_noPU_nonprompt_gen_2sigma_EE;
  TH1D *h_PU200_prompt_gen_3sigma_EE, *h_PU200_nonprompt_gen_3sigma_EE, *h_noPU_prompt_gen_3sigma_EE, *h_noPU_nonprompt_gen_3sigma_EE;
  TH1D *h_PU200_prompt_gen_4sigma_EE, *h_PU200_nonprompt_gen_4sigma_EE, *h_noPU_prompt_gen_4sigma_EE, *h_noPU_nonprompt_gen_4sigma_EE;
      // vtx
  TH1D *h_PU200_prompt_2sigma_EE_vtx, *h_PU200_nonprompt_2sigma_EE_vtx, *h_noPU_prompt_2sigma_EE_vtx, *h_noPU_nonprompt_2sigma_EE_vtx;
  TH1D *h_PU200_prompt_3sigma_EE_vtx, *h_PU200_nonprompt_3sigma_EE_vtx, *h_noPU_prompt_3sigma_EE_vtx, *h_noPU_nonprompt_3sigma_EE_vtx;
  TH1D *h_PU200_prompt_4sigma_EE_vtx, *h_PU200_nonprompt_4sigma_EE_vtx, *h_noPU_prompt_4sigma_EE_vtx, *h_noPU_nonprompt_4sigma_EE_vtx;

  /*
  TH1D *h_PU200_prompt_EE, *h_PU200_nonprompt_EE, *h_noPU_prompt_EE, *h_noPU_nonprompt_EE;
  TH1D *h_PU200_prompt_2sigma_EE, *h_PU200_nonprompt_2sigma_EE, *h_noPU_prompt_2sigma_EE, *h_noPU_2sigma_nonprompt_EE;
  TH1D *h_PU200_prompt_3sigma_EE, *h_PU200_nonprompt_3sigma_EE, *h_noPU_prompt_3sigma_EE, *h_noPU_3sigma_nonprompt_EE;
  TH1D *h_PU200_prompt_4sigma_EE, *h_PU200_nonprompt_4sigma_EE, *h_noPU_prompt_4sigma_EE, *h_noPU_4sigma_nonprompt_EE;
  TH1D *h_PU200_prompt_40_EE, *h_PU200_nonprompt_40_EE, *h_noPU_prompt_40_EE, *h_noPU_40_nonprompt_EE;
  TH1D *h_PU200_prompt_60_EE, *h_PU200_nonprompt_60_EE, *h_noPU_prompt_60_EE, *h_noPU_60_nonprompt_EE;
  TH1D *h_PU200_prompt_80_EE, *h_PU200_nonprompt_80_EE, *h_noPU_prompt_80_EE, *h_noPU_80_nonprompt_EE;
  TH1D *h_PU200_prompt_100_EE, *h_PU200_nonprompt_100_EE, *h_noPU_prompt_100_EE, *h_noPU_100_nonprompt_EE;
  TH1D *h_PU200_prompt_gen_EE, *h_PU200_nonprompt_gen_EE, *h_noPU_prompt_gen_EE, *h_noPU_nonprompt_gen_EE;
  TH1D *h_PU200_prompt_gen_2sigma_EE, *h_PU200_nonprompt_gen_2sigma_EE, *h_noPU_prompt_gen_2sigma_EE, *h_noPU_nonprompt_gen_2sigma_EE;
  */


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
  h_PU200_prompt_gen_4sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Sig_EB");
  h_PU200_prompt_gen_3sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Sig_EB");
  h_PU200_prompt_gen_2sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Sig_EB");
    // vtx
  h_PU200_prompt_4sigma_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EB");
  h_PU200_prompt_3sigma_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EB");
  h_PU200_prompt_2sigma_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EB");

  h_noPU_prompt_EB 	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_Sig_EB");
  h_noPU_prompt_4sigma_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EB");
  h_noPU_prompt_3sigma_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EB");
  h_noPU_prompt_2sigma_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EB");
  h_noPU_prompt_40_EB     = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Sig_EB");
  h_noPU_prompt_60_EB     = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Sig_EB");
  h_noPU_prompt_80_EB     = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Sig_EB");
  h_noPU_prompt_100_EB    = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Sig_EB");
  h_noPU_prompt_gen_EB	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_gen_Sig_EB");
  h_noPU_prompt_gen_4sigma_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Sig_EB");
  h_noPU_prompt_gen_3sigma_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Sig_EB");
  h_noPU_prompt_gen_2sigma_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Sig_EB");
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
  h_PU200_prompt_gen_4sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Sig_EE");
  h_PU200_prompt_gen_3sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Sig_EE");
  h_PU200_prompt_gen_2sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Sig_EE");
    // vtx
  h_PU200_prompt_4sigma_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EE");
  h_PU200_prompt_3sigma_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EE");
  h_PU200_prompt_2sigma_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EE");

  h_noPU_prompt_EE 	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_Sig_EE");
  h_noPU_prompt_4sigma_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EE");
  h_noPU_prompt_3sigma_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EE");
  h_noPU_prompt_2sigma_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EE");
  h_noPU_prompt_40_EE     = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Sig_EE");
  h_noPU_prompt_60_EE     = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Sig_EE");
  h_noPU_prompt_80_EE     = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Sig_EE");
  h_noPU_prompt_100_EE    = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Sig_EE");
  h_noPU_prompt_gen_EE	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_gen_Sig_EE");
  h_noPU_prompt_gen_4sigma_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Sig_EE");
  h_noPU_prompt_gen_3sigma_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Sig_EE");
  h_noPU_prompt_gen_2sigma_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Sig_EE");
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
  h_PU200_nonprompt_gen_4sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Bkg_EB");
  h_PU200_nonprompt_gen_3sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Bkg_EB");
  h_PU200_nonprompt_gen_2sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Bkg_EB");
      // vtx
  h_PU200_nonprompt_4sigma_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EB");
  h_PU200_nonprompt_3sigma_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EB");
  h_PU200_nonprompt_2sigma_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EB");

  h_noPU_nonprompt_EB 	      = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_Bkg_EB");
  h_noPU_nonprompt_4sigma_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EB");
  h_noPU_nonprompt_3sigma_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EB");
  h_noPU_nonprompt_2sigma_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EB");
  h_noPU_nonprompt_40_EB     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Bkg_EB");
  h_noPU_nonprompt_60_EB     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Bkg_EB");
  h_noPU_nonprompt_80_EB     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Bkg_EB");
  h_noPU_nonprompt_100_EB    = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Bkg_EB");
  h_noPU_nonprompt_gen_EB    = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_gen_Bkg_EB");
  h_noPU_nonprompt_gen_4sigma_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Bkg_EB");
  h_noPU_nonprompt_gen_3sigma_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Bkg_EB");
  h_noPU_nonprompt_gen_2sigma_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Bkg_EB");
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
  h_PU200_nonprompt_gen_4sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Bkg_EE");
  h_PU200_nonprompt_gen_3sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Bkg_EE");
  h_PU200_nonprompt_gen_2sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Bkg_EE");
      // vtx
  h_PU200_nonprompt_4sigma_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EE");
  h_PU200_nonprompt_3sigma_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EE");
  h_PU200_nonprompt_2sigma_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EE");

  h_noPU_nonprompt_EE 	      = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_Bkg_EE");
  h_noPU_nonprompt_4sigma_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EE");
  h_noPU_nonprompt_3sigma_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EE");
  h_noPU_nonprompt_2sigma_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EE");
  h_noPU_nonprompt_40_EE     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Bkg_EE");
  h_noPU_nonprompt_60_EE     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Bkg_EE");
  h_noPU_nonprompt_80_EE     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Bkg_EE");
  h_noPU_nonprompt_100_EE    = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Bkg_EE");
  h_noPU_nonprompt_gen_EE    = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_gen_Bkg_EE");
  h_noPU_nonprompt_gen_4sigma_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Bkg_EE");
  h_noPU_nonprompt_gen_3sigma_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Bkg_EE");
  h_noPU_nonprompt_gen_2sigma_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Bkg_EE");
  int nbin = h_PU200_prompt_EB->GetNbinsX()+2;

  cout << "total muon of PU200 prompt in Barrel   : " << h_PU200_prompt_EB->Integral(1,nbin) << endl;
  cout << "total muon of PU200 prompt in Endcap   : " << h_PU200_prompt_EE->Integral(1,nbin) << endl;
  cout << "total muon of PU200 nonprompt in Barrel: " << h_PU200_nonprompt_EB->Integral(1,nbin) << endl;
  cout << "total muon of PU200 nonprompt in Endcap: " << h_PU200_nonprompt_EE->Integral(1,nbin) << endl;
  cout << "total muon of noPU prompt in Barrel    : " << h_noPU_prompt_EB->Integral(1,nbin) << endl;
  cout << "total muon of noPU prompt in Endcap    : " << h_noPU_prompt_EE->Integral(1,nbin) << endl;
  cout << "total muon of noPU nonprompt in Barrel : " << h_noPU_nonprompt_EB->Integral(1,nbin) << endl;
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
        // vtx
  vector<double> prompt_eff_PU200_2sigma_EB_vtx={0}, prompt_norm_PU200_2sigma_EB_vtx={0}, prompt_eff_PU200_2sigma_EE_vtx={0}, prompt_norm_PU200_2sigma_EE_vtx={0};
  vector<double> prompt_eff_PU200_3sigma_EB_vtx={0}, prompt_norm_PU200_3sigma_EB_vtx={0}, prompt_eff_PU200_3sigma_EE_vtx={0}, prompt_norm_PU200_3sigma_EE_vtx={0};
  vector<double> prompt_eff_PU200_4sigma_EB_vtx={0}, prompt_norm_PU200_4sigma_EB_vtx={0}, prompt_eff_PU200_4sigma_EE_vtx={0}, prompt_norm_PU200_4sigma_EE_vtx={0};
      // nonprompt
  vector<double> nonprompt_eff_PU200_EB={0}, nonprompt_norm_PU200_EB={0}, nonprompt_eff_PU200_EE={0}, nonprompt_norm_PU200_EE={0};
  vector<double> nonprompt_eff_PU200_2sigma_EB={0}, nonprompt_norm_PU200_2sigma_EB={0}, nonprompt_eff_PU200_2sigma_EE={0}, nonprompt_norm_PU200_2sigma_EE={0};
  vector<double> nonprompt_eff_PU200_3sigma_EB={0}, nonprompt_norm_PU200_3sigma_EB={0}, nonprompt_eff_PU200_3sigma_EE={0}, nonprompt_norm_PU200_3sigma_EE={0};
  vector<double> nonprompt_eff_PU200_4sigma_EB={0}, nonprompt_norm_PU200_4sigma_EB={0}, nonprompt_eff_PU200_4sigma_EE={0}, nonprompt_norm_PU200_4sigma_EE={0};
  vector<double> nonprompt_eff_PU200_40_EB={0}, nonprompt_norm_PU200_40_EB={0}, nonprompt_eff_PU200_40_EE={0}, nonprompt_norm_PU200_40_EE={0};
  vector<double> nonprompt_eff_PU200_60_EB={0}, nonprompt_norm_PU200_60_EB={0}, nonprompt_eff_PU200_60_EE={0}, nonprompt_norm_PU200_60_EE={0};
  vector<double> nonprompt_eff_PU200_80_EB={0}, nonprompt_norm_PU200_80_EB={0}, nonprompt_eff_PU200_80_EE={0}, nonprompt_norm_PU200_80_EE={0};
  vector<double> nonprompt_eff_PU200_100_EB={0}, nonprompt_norm_PU200_100_EB={0}, nonprompt_eff_PU200_100_EE={0}, nonprompt_norm_PU200_100_EE={0};
        // vtx
  vector<double> nonprompt_eff_PU200_2sigma_EB_vtx={0}, nonprompt_norm_PU200_2sigma_EB_vtx={0}, nonprompt_eff_PU200_2sigma_EE_vtx={0}, nonprompt_norm_PU200_2sigma_EE_vtx={0};
  vector<double> nonprompt_eff_PU200_3sigma_EB_vtx={0}, nonprompt_norm_PU200_3sigma_EB_vtx={0}, nonprompt_eff_PU200_3sigma_EE_vtx={0}, nonprompt_norm_PU200_3sigma_EE_vtx={0};
  vector<double> nonprompt_eff_PU200_4sigma_EB_vtx={0}, nonprompt_norm_PU200_4sigma_EB_vtx={0}, nonprompt_eff_PU200_4sigma_EE_vtx={0}, nonprompt_norm_PU200_4sigma_EE_vtx={0};
      // GEN case
        // prompt
  vector<double> prompt_eff_PU200_gen_EB={0}, prompt_norm_PU200_gen_EB={0}, prompt_eff_PU200_gen_EE={0}, prompt_norm_PU200_gen_EE={0};
  vector<double> prompt_eff_PU200_gen_2sigma_EB={0}, prompt_norm_PU200_gen_2sigma_EB={0}, prompt_eff_PU200_gen_2sigma_EE={0}, prompt_norm_PU200_gen_2sigma_EE={0};
  vector<double> prompt_eff_PU200_gen_3sigma_EB={0}, prompt_norm_PU200_gen_3sigma_EB={0}, prompt_eff_PU200_gen_3sigma_EE={0}, prompt_norm_PU200_gen_3sigma_EE={0};
  vector<double> prompt_eff_PU200_gen_4sigma_EB={0}, prompt_norm_PU200_gen_4sigma_EB={0}, prompt_eff_PU200_gen_4sigma_EE={0}, prompt_norm_PU200_gen_4sigma_EE={0};
        // nonprompt
  vector<double> nonprompt_eff_PU200_gen_EB={0}, nonprompt_norm_PU200_gen_EB={0}, nonprompt_eff_PU200_gen_EE={0}, nonprompt_norm_PU200_gen_EE={0};
  vector<double> nonprompt_eff_PU200_gen_2sigma_EB={0}, nonprompt_norm_PU200_gen_2sigma_EB={0}, nonprompt_eff_PU200_gen_2sigma_EE={0}, nonprompt_norm_PU200_gen_2sigma_EE={0};
  vector<double> nonprompt_eff_PU200_gen_3sigma_EB={0}, nonprompt_norm_PU200_gen_3sigma_EB={0}, nonprompt_eff_PU200_gen_3sigma_EE={0}, nonprompt_norm_PU200_gen_3sigma_EE={0};
  vector<double> nonprompt_eff_PU200_gen_4sigma_EB={0}, nonprompt_norm_PU200_gen_4sigma_EB={0}, nonprompt_eff_PU200_gen_4sigma_EE={0}, nonprompt_norm_PU200_gen_4sigma_EE={0};
    // noPU
      // prompt
  vector<double> prompt_eff_noPU_EB={0}, prompt_norm_noPU_EB={0}, prompt_eff_noPU_EE={0}, prompt_norm_noPU_EE={0};
  vector<double> prompt_eff_noPU_2sigma_EB={0}, prompt_norm_noPU_2sigma_EB={0}, prompt_eff_noPU_2sigma_EE={0}, prompt_norm_noPU_2sigma_EE={0};
  vector<double> prompt_eff_noPU_3sigma_EB={0}, prompt_norm_noPU_3sigma_EB={0}, prompt_eff_noPU_3sigma_EE={0}, prompt_norm_noPU_3sigma_EE={0};
  vector<double> prompt_eff_noPU_4sigma_EB={0}, prompt_norm_noPU_4sigma_EB={0}, prompt_eff_noPU_4sigma_EE={0}, prompt_norm_noPU_4sigma_EE={0};
  vector<double> prompt_eff_noPU_40_EB={0}, prompt_norm_noPU_40_EB={0}, prompt_eff_noPU_40_EE={0}, prompt_norm_noPU_40_EE={0};
  vector<double> prompt_eff_noPU_60_EB={0}, prompt_norm_noPU_60_EB={0}, prompt_eff_noPU_60_EE={0}, prompt_norm_noPU_60_EE={0};
  vector<double> prompt_eff_noPU_80_EB={0}, prompt_norm_noPU_80_EB={0}, prompt_eff_noPU_80_EE={0}, prompt_norm_noPU_80_EE={0};
  vector<double> prompt_eff_noPU_100_EB={0}, prompt_norm_noPU_100_EB={0}, prompt_eff_noPU_100_EE={0}, prompt_norm_noPU_100_EE={0};
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
      // prompt
  vector<double> prompt_eff_noPU_gen_EB={0}, prompt_norm_noPU_gen_EB={0}, prompt_eff_noPU_gen_EE={0}, prompt_norm_noPU_gen_EE={0};
  vector<double> prompt_eff_noPU_gen_2sigma_EB={0}, prompt_norm_noPU_gen_2sigma_EB={0}, prompt_eff_noPU_gen_2sigma_EE={0}, prompt_norm_noPU_gen_2sigma_EE={0};
  vector<double> prompt_eff_noPU_gen_3sigma_EB={0}, prompt_norm_noPU_gen_3sigma_EB={0}, prompt_eff_noPU_gen_3sigma_EE={0}, prompt_norm_noPU_gen_3sigma_EE={0};
  vector<double> prompt_eff_noPU_gen_4sigma_EB={0}, prompt_norm_noPU_gen_4sigma_EB={0}, prompt_eff_noPU_gen_4sigma_EE={0}, prompt_norm_noPU_gen_4sigma_EE={0};
      // nonprompt
  vector<double> nonprompt_eff_noPU_gen_EB={0}, nonprompt_norm_noPU_gen_EB={0}, nonprompt_eff_noPU_gen_EE={0}, nonprompt_norm_noPU_gen_EE={0};
  vector<double> nonprompt_eff_noPU_gen_2sigma_EB={0}, nonprompt_norm_noPU_gen_2sigma_EB={0}, nonprompt_eff_noPU_gen_2sigma_EE={0}, nonprompt_norm_noPU_gen_2sigma_EE={0};
  vector<double> nonprompt_eff_noPU_gen_3sigma_EB={0}, nonprompt_norm_noPU_gen_3sigma_EB={0}, nonprompt_eff_noPU_gen_3sigma_EE={0}, nonprompt_norm_noPU_gen_3sigma_EE={0};
  vector<double> nonprompt_eff_noPU_gen_4sigma_EB={0}, nonprompt_norm_noPU_gen_4sigma_EB={0}, nonprompt_eff_noPU_gen_4sigma_EE={0}, nonprompt_norm_noPU_gen_4sigma_EE={0};


  //////////////////////////
  // Calculate efficiency //
  //////////////////////////
  for(unsigned int i=0; i<h_PU200_prompt_EB->GetNbinsX()+2; i++) {
    // prompt
      // Barrel region
        // efficiency
	  // PU200
    prompt_eff_PU200_EB.emplace_back(h_PU200_prompt_EB->Integral(1,i+1)/h_PU200_prompt_EB->Integral(1,nbin));
    prompt_eff_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->Integral(1,i+1)/h_PU200_prompt_40_EB->Integral(1,nbin));
    prompt_eff_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->Integral(1,i+1)/h_PU200_prompt_60_EB->Integral(1,nbin));
    prompt_eff_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->Integral(1,i+1)/h_PU200_prompt_80_EB->Integral(1,nbin));
    prompt_eff_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->Integral(1,i+1)/h_PU200_prompt_100_EB->Integral(1,nbin));
    prompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(1,i+1)/h_PU200_prompt_2sigma_EB->Integral(1,nbin));
    prompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(1,i+1)/h_PU200_prompt_3sigma_EB->Integral(1,nbin));
    prompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(1,i+1)/h_PU200_prompt_4sigma_EB->Integral(1,nbin));
              // vtx
    prompt_eff_PU200_2sigma_EB_vtx.emplace_back(h_PU200_prompt_2sigma_EB_vtx->Integral(1,i+1)/h_PU200_prompt_2sigma_EB_vtx->Integral(1,nbin));
    prompt_eff_PU200_3sigma_EB_vtx.emplace_back(h_PU200_prompt_3sigma_EB_vtx->Integral(1,i+1)/h_PU200_prompt_3sigma_EB_vtx->Integral(1,nbin));
    prompt_eff_PU200_4sigma_EB_vtx.emplace_back(h_PU200_prompt_4sigma_EB_vtx->Integral(1,i+1)/h_PU200_prompt_4sigma_EB_vtx->Integral(1,nbin));
            // GEN case
    prompt_eff_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->Integral(1,i+1)/h_PU200_prompt_gen_EB->Integral(1,nbin));
    prompt_eff_PU200_gen_2sigma_EB.emplace_back(h_PU200_prompt_gen_2sigma_EB->Integral(1,i+1)/h_PU200_prompt_gen_2sigma_EB->Integral(1,nbin));
    prompt_eff_PU200_gen_3sigma_EB.emplace_back(h_PU200_prompt_gen_3sigma_EB->Integral(1,i+1)/h_PU200_prompt_gen_3sigma_EB->Integral(1,nbin));
    prompt_eff_PU200_gen_4sigma_EB.emplace_back(h_PU200_prompt_gen_4sigma_EB->Integral(1,i+1)/h_PU200_prompt_gen_4sigma_EB->Integral(1,nbin));
	  // noPU
    prompt_eff_noPU_EB.emplace_back(h_noPU_prompt_EB->Integral(1,i+1)/h_noPU_prompt_EB->Integral(1,nbin));
    prompt_eff_noPU_40_EB.emplace_back(h_noPU_prompt_40_EB->Integral(1,i+1)/h_noPU_prompt_40_EB->Integral(1,nbin));
    prompt_eff_noPU_60_EB.emplace_back(h_noPU_prompt_60_EB->Integral(1,i+1)/h_noPU_prompt_60_EB->Integral(1,nbin));
    prompt_eff_noPU_80_EB.emplace_back(h_noPU_prompt_80_EB->Integral(1,i+1)/h_noPU_prompt_80_EB->Integral(1,nbin));
    prompt_eff_noPU_100_EB.emplace_back(h_noPU_prompt_100_EB->Integral(1,i+1)/h_noPU_prompt_100_EB->Integral(1,nbin));
    prompt_eff_noPU_2sigma_EB.emplace_back(h_noPU_prompt_2sigma_EB->Integral(1,i+1)/h_noPU_prompt_2sigma_EB->Integral(1,nbin));
    prompt_eff_noPU_3sigma_EB.emplace_back(h_noPU_prompt_3sigma_EB->Integral(1,i+1)/h_noPU_prompt_3sigma_EB->Integral(1,nbin));
    prompt_eff_noPU_4sigma_EB.emplace_back(h_noPU_prompt_4sigma_EB->Integral(1,i+1)/h_noPU_prompt_4sigma_EB->Integral(1,nbin));
            // GEN case
    prompt_eff_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(1,i+1)/h_noPU_prompt_gen_EB->Integral(1,nbin));
    prompt_eff_noPU_gen_2sigma_EB.emplace_back(h_noPU_prompt_gen_2sigma_EB->Integral(1,i+1)/h_noPU_prompt_gen_2sigma_EB->Integral(1,nbin));
    prompt_eff_noPU_gen_3sigma_EB.emplace_back(h_noPU_prompt_gen_3sigma_EB->Integral(1,i+1)/h_noPU_prompt_gen_3sigma_EB->Integral(1,nbin));
    prompt_eff_noPU_gen_4sigma_EB.emplace_back(h_noPU_prompt_gen_4sigma_EB->Integral(1,i+1)/h_noPU_prompt_gen_4sigma_EB->Integral(1,nbin));

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
	  // PU200
    prompt_eff_PU200_EE.emplace_back(h_PU200_prompt_EE->Integral(1,i+1)/h_PU200_prompt_EE->Integral(1,nbin));
    prompt_eff_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->Integral(1,i+1)/h_PU200_prompt_40_EE->Integral(1,nbin));
    prompt_eff_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->Integral(1,i+1)/h_PU200_prompt_60_EE->Integral(1,nbin));
    prompt_eff_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->Integral(1,i+1)/h_PU200_prompt_80_EE->Integral(1,nbin));
    prompt_eff_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->Integral(1,i+1)/h_PU200_prompt_100_EE->Integral(1,nbin));
    prompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->Integral(1,i+1)/h_PU200_prompt_2sigma_EE->Integral(1,nbin));
    prompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->Integral(1,i+1)/h_PU200_prompt_3sigma_EE->Integral(1,nbin));
    prompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->Integral(1,i+1)/h_PU200_prompt_4sigma_EE->Integral(1,nbin));
              // vtx
    prompt_eff_PU200_2sigma_EE_vtx.emplace_back(h_PU200_prompt_2sigma_EE_vtx->Integral(1,i+1)/h_PU200_prompt_2sigma_EE_vtx->Integral(1,nbin));
    prompt_eff_PU200_3sigma_EE_vtx.emplace_back(h_PU200_prompt_3sigma_EE_vtx->Integral(1,i+1)/h_PU200_prompt_3sigma_EE_vtx->Integral(1,nbin));
    prompt_eff_PU200_4sigma_EE_vtx.emplace_back(h_PU200_prompt_4sigma_EE_vtx->Integral(1,i+1)/h_PU200_prompt_4sigma_EE_vtx->Integral(1,nbin));
            // GEN case
    prompt_eff_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->Integral(1,i+1)/h_PU200_prompt_gen_EE->Integral(1,nbin));
    prompt_eff_PU200_gen_2sigma_EE.emplace_back(h_PU200_prompt_gen_2sigma_EE->Integral(1,i+1)/h_PU200_prompt_gen_2sigma_EE->Integral(1,nbin));
    prompt_eff_PU200_gen_3sigma_EE.emplace_back(h_PU200_prompt_gen_3sigma_EE->Integral(1,i+1)/h_PU200_prompt_gen_3sigma_EE->Integral(1,nbin));
    prompt_eff_PU200_gen_4sigma_EE.emplace_back(h_PU200_prompt_gen_4sigma_EE->Integral(1,i+1)/h_PU200_prompt_gen_4sigma_EE->Integral(1,nbin));
	  // noPU
    prompt_eff_noPU_EE.emplace_back(h_noPU_prompt_EE->Integral(1,i+1)/h_noPU_prompt_EE->Integral(1,nbin));
    prompt_eff_noPU_40_EE.emplace_back(h_noPU_prompt_40_EE->Integral(1,i+1)/h_noPU_prompt_40_EE->Integral(1,nbin));
    prompt_eff_noPU_60_EE.emplace_back(h_noPU_prompt_60_EE->Integral(1,i+1)/h_noPU_prompt_60_EE->Integral(1,nbin));
    prompt_eff_noPU_80_EE.emplace_back(h_noPU_prompt_80_EE->Integral(1,i+1)/h_noPU_prompt_80_EE->Integral(1,nbin));
    prompt_eff_noPU_100_EE.emplace_back(h_noPU_prompt_100_EE->Integral(1,i+1)/h_noPU_prompt_100_EE->Integral(1,nbin));
    prompt_eff_noPU_2sigma_EE.emplace_back(h_noPU_prompt_2sigma_EE->Integral(1,i+1)/h_noPU_prompt_2sigma_EE->Integral(1,nbin));
    prompt_eff_noPU_3sigma_EE.emplace_back(h_noPU_prompt_3sigma_EE->Integral(1,i+1)/h_noPU_prompt_3sigma_EE->Integral(1,nbin));
    prompt_eff_noPU_4sigma_EE.emplace_back(h_noPU_prompt_4sigma_EE->Integral(1,i+1)/h_noPU_prompt_4sigma_EE->Integral(1,nbin));
            // GEN case
    prompt_eff_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->Integral(1,i+1)/h_noPU_prompt_gen_EE->Integral(1,nbin));
    prompt_eff_noPU_gen_2sigma_EE.emplace_back(h_noPU_prompt_gen_2sigma_EE->Integral(1,i+1)/h_noPU_prompt_gen_2sigma_EE->Integral(1,nbin));
    prompt_eff_noPU_gen_3sigma_EE.emplace_back(h_noPU_prompt_gen_3sigma_EE->Integral(1,i+1)/h_noPU_prompt_gen_3sigma_EE->Integral(1,nbin));
    prompt_eff_noPU_gen_4sigma_EE.emplace_back(h_noPU_prompt_gen_4sigma_EE->Integral(1,i+1)/h_noPU_prompt_gen_4sigma_EE->Integral(1,nbin));

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
	  // PU200
    nonprompt_eff_PU200_EB.emplace_back(h_PU200_nonprompt_EB->Integral(1,i+1)/h_PU200_nonprompt_EB->Integral(1,nbin));
    nonprompt_eff_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->Integral(1,i+1)/h_PU200_nonprompt_40_EB->Integral(1,nbin));
    nonprompt_eff_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->Integral(1,i+1)/h_PU200_nonprompt_60_EB->Integral(1,nbin));
    nonprompt_eff_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->Integral(1,i+1)/h_PU200_nonprompt_80_EB->Integral(1,nbin));
    nonprompt_eff_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->Integral(1,i+1)/h_PU200_nonprompt_100_EB->Integral(1,nbin));
    nonprompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EB->Integral(1,nbin));
    nonprompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EB->Integral(1,nbin));
    nonprompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EB->Integral(1,nbin));
              // vtx
    nonprompt_eff_PU200_2sigma_EB_vtx.emplace_back(h_PU200_nonprompt_2sigma_EB_vtx->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EB_vtx->Integral(1,nbin));
    nonprompt_eff_PU200_3sigma_EB_vtx.emplace_back(h_PU200_nonprompt_3sigma_EB_vtx->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EB_vtx->Integral(1,nbin));
    nonprompt_eff_PU200_4sigma_EB_vtx.emplace_back(h_PU200_nonprompt_4sigma_EB_vtx->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EB_vtx->Integral(1,nbin));
            // GEN case
    nonprompt_eff_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->Integral(1,i+1)/h_PU200_nonprompt_gen_EB->Integral(1,nbin));
    nonprompt_eff_PU200_gen_2sigma_EB.emplace_back(h_PU200_nonprompt_gen_2sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_gen_2sigma_EB->Integral(1,nbin));
    nonprompt_eff_PU200_gen_3sigma_EB.emplace_back(h_PU200_nonprompt_gen_3sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_gen_3sigma_EB->Integral(1,nbin));
    nonprompt_eff_PU200_gen_4sigma_EB.emplace_back(h_PU200_nonprompt_gen_4sigma_EB->Integral(1,i+1)/h_PU200_nonprompt_gen_4sigma_EB->Integral(1,nbin));
	  // noPU
    nonprompt_eff_noPU_EB.emplace_back(h_noPU_nonprompt_EB->Integral(1,i+1)/h_noPU_nonprompt_EB->Integral(1,nbin));
    nonprompt_eff_noPU_40_EB.emplace_back(h_noPU_nonprompt_40_EB->Integral(1,i+1)/h_noPU_nonprompt_40_EB->Integral(1,nbin));
    nonprompt_eff_noPU_60_EB.emplace_back(h_noPU_nonprompt_60_EB->Integral(1,i+1)/h_noPU_nonprompt_60_EB->Integral(1,nbin));
    nonprompt_eff_noPU_80_EB.emplace_back(h_noPU_nonprompt_80_EB->Integral(1,i+1)/h_noPU_nonprompt_80_EB->Integral(1,nbin));
    nonprompt_eff_noPU_100_EB.emplace_back(h_noPU_nonprompt_100_EB->Integral(1,i+1)/h_noPU_nonprompt_100_EB->Integral(1,nbin));
    nonprompt_eff_noPU_2sigma_EB.emplace_back(h_noPU_nonprompt_2sigma_EB->Integral(1,i+1)/h_noPU_nonprompt_2sigma_EB->Integral(1,nbin));
    nonprompt_eff_noPU_3sigma_EB.emplace_back(h_noPU_nonprompt_3sigma_EB->Integral(1,i+1)/h_noPU_nonprompt_3sigma_EB->Integral(1,nbin));
    nonprompt_eff_noPU_4sigma_EB.emplace_back(h_noPU_nonprompt_4sigma_EB->Integral(1,i+1)/h_noPU_nonprompt_4sigma_EB->Integral(1,nbin));
            // GEN case
    nonprompt_eff_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->Integral(1,i+1)/h_noPU_nonprompt_gen_EB->Integral(1,nbin));
    nonprompt_eff_noPU_gen_2sigma_EB.emplace_back(h_noPU_nonprompt_gen_2sigma_EB->Integral(1,i+1)/h_noPU_nonprompt_gen_2sigma_EB->Integral(1,nbin));
    nonprompt_eff_noPU_gen_3sigma_EB.emplace_back(h_noPU_nonprompt_gen_3sigma_EB->Integral(1,i+1)/h_noPU_nonprompt_gen_3sigma_EB->Integral(1,nbin));
    nonprompt_eff_noPU_gen_4sigma_EB.emplace_back(h_noPU_nonprompt_gen_4sigma_EB->Integral(1,i+1)/h_noPU_nonprompt_gen_4sigma_EB->Integral(1,nbin));

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
	  // PU200
    nonprompt_eff_PU200_EE.emplace_back(h_PU200_nonprompt_EE->Integral(1,i+1)/h_PU200_nonprompt_EE->Integral(1,nbin));
    nonprompt_eff_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->Integral(1,i+1)/h_PU200_nonprompt_40_EE->Integral(1,nbin));
    nonprompt_eff_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->Integral(1,i+1)/h_PU200_nonprompt_60_EE->Integral(1,nbin));
    nonprompt_eff_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->Integral(1,i+1)/h_PU200_nonprompt_80_EE->Integral(1,nbin));
    nonprompt_eff_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->Integral(1,i+1)/h_PU200_nonprompt_100_EE->Integral(1,nbin));
    nonprompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EE->Integral(1,nbin));
    nonprompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EE->Integral(1,nbin));
    nonprompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EE->Integral(1,nbin));
              // vtx
    nonprompt_eff_PU200_2sigma_EE_vtx.emplace_back(h_PU200_nonprompt_2sigma_EE_vtx->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EE_vtx->Integral(1,nbin));
    nonprompt_eff_PU200_3sigma_EE_vtx.emplace_back(h_PU200_nonprompt_3sigma_EE_vtx->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EE_vtx->Integral(1,nbin));
    nonprompt_eff_PU200_4sigma_EE_vtx.emplace_back(h_PU200_nonprompt_4sigma_EE_vtx->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EE_vtx->Integral(1,nbin));
            // GEN case
    nonprompt_eff_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->Integral(1,i+1)/h_PU200_nonprompt_gen_EE->Integral(1,nbin));
    nonprompt_eff_PU200_gen_2sigma_EE.emplace_back(h_PU200_nonprompt_gen_2sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_gen_2sigma_EE->Integral(1,nbin));
    nonprompt_eff_PU200_gen_3sigma_EE.emplace_back(h_PU200_nonprompt_gen_3sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_gen_3sigma_EE->Integral(1,nbin));
    nonprompt_eff_PU200_gen_4sigma_EE.emplace_back(h_PU200_nonprompt_gen_4sigma_EE->Integral(1,i+1)/h_PU200_nonprompt_gen_4sigma_EE->Integral(1,nbin));
	  // noPU
    nonprompt_eff_noPU_EE.emplace_back(h_noPU_nonprompt_EE->Integral(1,i+1)/h_noPU_nonprompt_EE->Integral(1,nbin));
    nonprompt_eff_noPU_40_EE.emplace_back(h_noPU_nonprompt_40_EE->Integral(1,i+1)/h_noPU_nonprompt_40_EE->Integral(1,nbin));
    nonprompt_eff_noPU_60_EE.emplace_back(h_noPU_nonprompt_60_EE->Integral(1,i+1)/h_noPU_nonprompt_60_EE->Integral(1,nbin));
    nonprompt_eff_noPU_80_EE.emplace_back(h_noPU_nonprompt_80_EE->Integral(1,i+1)/h_noPU_nonprompt_80_EE->Integral(1,nbin));
    nonprompt_eff_noPU_100_EE.emplace_back(h_noPU_nonprompt_100_EE->Integral(1,i+1)/h_noPU_nonprompt_100_EE->Integral(1,nbin));
    nonprompt_eff_noPU_2sigma_EE.emplace_back(h_noPU_nonprompt_2sigma_EE->Integral(1,i+1)/h_noPU_nonprompt_2sigma_EE->Integral(1,nbin));
    nonprompt_eff_noPU_3sigma_EE.emplace_back(h_noPU_nonprompt_3sigma_EE->Integral(1,i+1)/h_noPU_nonprompt_3sigma_EE->Integral(1,nbin));
    nonprompt_eff_noPU_4sigma_EE.emplace_back(h_noPU_nonprompt_4sigma_EE->Integral(1,i+1)/h_noPU_nonprompt_4sigma_EE->Integral(1,nbin));
            // GEN case
    nonprompt_eff_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->Integral(1,i+1)/h_noPU_nonprompt_gen_EE->Integral(1,nbin));
    nonprompt_eff_noPU_gen_2sigma_EE.emplace_back(h_noPU_nonprompt_gen_2sigma_EE->Integral(1,i+1)/h_noPU_nonprompt_gen_2sigma_EE->Integral(1,nbin));
    nonprompt_eff_noPU_gen_3sigma_EE.emplace_back(h_noPU_nonprompt_gen_3sigma_EE->Integral(1,i+1)/h_noPU_nonprompt_gen_3sigma_EE->Integral(1,nbin));
    nonprompt_eff_noPU_gen_4sigma_EE.emplace_back(h_noPU_nonprompt_gen_4sigma_EE->Integral(1,i+1)/h_noPU_nonprompt_gen_4sigma_EE->Integral(1,nbin));

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

  // Define TGraph
  // Prompt
    // Barrel region
      // PU200
  TGraph* gr_eff_PU200_prompt_EB = new TGraph();              TGraph* gr_norm_PU200_prompt_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_40_EB = new TGraph();           TGraph* gr_norm_PU200_prompt_40_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_60_EB = new TGraph();           TGraph* gr_norm_PU200_prompt_60_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_80_EB = new TGraph();           TGraph* gr_norm_PU200_prompt_80_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_100_EB = new TGraph();          TGraph* gr_norm_PU200_prompt_100_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_2sigma_EB = new TGraph();       TGraph* gr_norm_PU200_prompt_2sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_3sigma_EB = new TGraph();       TGraph* gr_norm_PU200_prompt_3sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_4sigma_EB = new TGraph();       TGraph* gr_norm_PU200_prompt_4sigma_EB = new TGraph();

  TGraph* gr_eff_PU200_prompt_gen_EB = new TGraph();              TGraph* gr_norm_PU200_prompt_gen_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_gen_2sigma_EB = new TGraph();       TGraph* gr_norm_PU200_prompt_gen_2sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_gen_3sigma_EB = new TGraph();       TGraph* gr_norm_PU200_prompt_gen_3sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_prompt_gen_4sigma_EB = new TGraph();       TGraph* gr_norm_PU200_prompt_gen_4sigma_EB = new TGraph();
        // vtx
  TGraph* gr_eff_PU200_prompt_2sigma_EB_vtx = new TGraph();       TGraph* gr_norm_PU200_prompt_2sigma_EB_vtx = new TGraph();
  TGraph* gr_eff_PU200_prompt_3sigma_EB_vtx = new TGraph();       TGraph* gr_norm_PU200_prompt_3sigma_EB_vtx = new TGraph();
  TGraph* gr_eff_PU200_prompt_4sigma_EB_vtx = new TGraph();       TGraph* gr_norm_PU200_prompt_4sigma_EB_vtx = new TGraph();
      // noPU
  TGraph* gr_eff_noPU_prompt_EB = new TGraph();              TGraph* gr_norm_noPU_prompt_EB = new TGraph();
  TGraph* gr_eff_noPU_prompt_40_EB = new TGraph();           TGraph* gr_norm_noPU_prompt_40_EB = new TGraph();
  TGraph* gr_eff_noPU_prompt_60_EB = new TGraph();           TGraph* gr_norm_noPU_prompt_60_EB = new TGraph();
  TGraph* gr_eff_noPU_prompt_80_EB = new TGraph();           TGraph* gr_norm_noPU_prompt_80_EB = new TGraph();
  TGraph* gr_eff_noPU_prompt_100_EB = new TGraph();          TGraph* gr_norm_noPU_prompt_100_EB = new TGraph();
  TGraph* gr_eff_noPU_prompt_2sigma_EB = new TGraph();       TGraph* gr_norm_noPU_prompt_2sigma_EB = new TGraph();
  TGraph* gr_eff_noPU_prompt_3sigma_EB = new TGraph();       TGraph* gr_norm_noPU_prompt_3sigma_EB = new TGraph();
  TGraph* gr_eff_noPU_prompt_4sigma_EB = new TGraph();       TGraph* gr_norm_noPU_prompt_4sigma_EB = new TGraph();

  TGraph* gr_eff_noPU_prompt_gen_EB = new TGraph();              TGraph* gr_norm_noPU_prompt_gen_EB = new TGraph();
  TGraph* gr_eff_noPU_prompt_gen_2sigma_EB = new TGraph();       TGraph* gr_norm_noPU_prompt_gen_2sigma_EB = new TGraph();
  TGraph* gr_eff_noPU_prompt_gen_3sigma_EB = new TGraph();       TGraph* gr_norm_noPU_prompt_gen_3sigma_EB = new TGraph();
  TGraph* gr_eff_noPU_prompt_gen_4sigma_EB = new TGraph();       TGraph* gr_norm_noPU_prompt_gen_4sigma_EB = new TGraph();

    // Endcap region
      // PU200
  TGraph* gr_eff_PU200_prompt_EE = new TGraph();              TGraph* gr_norm_PU200_prompt_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_40_EE = new TGraph();           TGraph* gr_norm_PU200_prompt_40_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_60_EE = new TGraph();           TGraph* gr_norm_PU200_prompt_60_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_80_EE = new TGraph();           TGraph* gr_norm_PU200_prompt_80_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_100_EE = new TGraph();          TGraph* gr_norm_PU200_prompt_100_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_2sigma_EE = new TGraph();       TGraph* gr_norm_PU200_prompt_2sigma_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_3sigma_EE = new TGraph();       TGraph* gr_norm_PU200_prompt_3sigma_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_4sigma_EE = new TGraph();       TGraph* gr_norm_PU200_prompt_4sigma_EE = new TGraph();

  TGraph* gr_eff_PU200_prompt_gen_EE = new TGraph();              TGraph* gr_norm_PU200_prompt_gen_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_gen_2sigma_EE = new TGraph();       TGraph* gr_norm_PU200_prompt_gen_2sigma_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_gen_3sigma_EE = new TGraph();       TGraph* gr_norm_PU200_prompt_gen_3sigma_EE = new TGraph();
  TGraph* gr_eff_PU200_prompt_gen_4sigma_EE = new TGraph();       TGraph* gr_norm_PU200_prompt_gen_4sigma_EE = new TGraph();
        // vtx
  TGraph* gr_eff_PU200_prompt_2sigma_EE_vtx = new TGraph();       TGraph* gr_norm_PU200_prompt_2sigma_EE_vtx = new TGraph();
  TGraph* gr_eff_PU200_prompt_3sigma_EE_vtx = new TGraph();       TGraph* gr_norm_PU200_prompt_3sigma_EE_vtx = new TGraph();
  TGraph* gr_eff_PU200_prompt_4sigma_EE_vtx = new TGraph();       TGraph* gr_norm_PU200_prompt_4sigma_EE_vtx = new TGraph();
      // noPU
  TGraph* gr_eff_noPU_prompt_EE = new TGraph();              TGraph* gr_norm_noPU_prompt_EE = new TGraph();
  TGraph* gr_eff_noPU_prompt_40_EE = new TGraph();           TGraph* gr_norm_noPU_prompt_40_EE = new TGraph();
  TGraph* gr_eff_noPU_prompt_60_EE = new TGraph();           TGraph* gr_norm_noPU_prompt_60_EE = new TGraph();
  TGraph* gr_eff_noPU_prompt_80_EE = new TGraph();           TGraph* gr_norm_noPU_prompt_80_EE = new TGraph();
  TGraph* gr_eff_noPU_prompt_100_EE = new TGraph();          TGraph* gr_norm_noPU_prompt_100_EE = new TGraph();
  TGraph* gr_eff_noPU_prompt_2sigma_EE = new TGraph();       TGraph* gr_norm_noPU_prompt_2sigma_EE = new TGraph();
  TGraph* gr_eff_noPU_prompt_3sigma_EE = new TGraph();       TGraph* gr_norm_noPU_prompt_3sigma_EE = new TGraph();
  TGraph* gr_eff_noPU_prompt_4sigma_EE = new TGraph();       TGraph* gr_norm_noPU_prompt_4sigma_EE = new TGraph();

  TGraph* gr_eff_noPU_prompt_gen_EE = new TGraph();              TGraph* gr_norm_noPU_prompt_gen_EE = new TGraph();
  TGraph* gr_eff_noPU_prompt_gen_2sigma_EE = new TGraph();       TGraph* gr_norm_noPU_prompt_gen_2sigma_EE = new TGraph();
  TGraph* gr_eff_noPU_prompt_gen_3sigma_EE = new TGraph();       TGraph* gr_norm_noPU_prompt_gen_3sigma_EE = new TGraph();
  TGraph* gr_eff_noPU_prompt_gen_4sigma_EE = new TGraph();       TGraph* gr_norm_noPU_prompt_gen_4sigma_EE = new TGraph();

  // Nonprompt
    // Barrel region
      // PU200
  TGraph* gr_eff_PU200_nonprompt_EB = new TGraph();              TGraph* gr_norm_PU200_nonprompt_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_40_EB = new TGraph();           TGraph* gr_norm_PU200_nonprompt_40_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_60_EB = new TGraph();           TGraph* gr_norm_PU200_nonprompt_60_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_80_EB = new TGraph();           TGraph* gr_norm_PU200_nonprompt_80_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_100_EB = new TGraph();          TGraph* gr_norm_PU200_nonprompt_100_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_2sigma_EB = new TGraph();       TGraph* gr_norm_PU200_nonprompt_2sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_3sigma_EB = new TGraph();       TGraph* gr_norm_PU200_nonprompt_3sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_4sigma_EB = new TGraph();       TGraph* gr_norm_PU200_nonprompt_4sigma_EB = new TGraph();

  TGraph* gr_eff_PU200_nonprompt_gen_EB = new TGraph();              TGraph* gr_norm_PU200_nonprompt_gen_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_gen_2sigma_EB = new TGraph();       TGraph* gr_norm_PU200_nonprompt_gen_2sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_gen_3sigma_EB = new TGraph();       TGraph* gr_norm_PU200_nonprompt_gen_3sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_gen_4sigma_EB = new TGraph();       TGraph* gr_norm_PU200_nonprompt_gen_4sigma_EB = new TGraph();
        // vtx
  TGraph* gr_eff_PU200_nonprompt_2sigma_EB_vtx = new TGraph();       TGraph* gr_norm_PU200_nonprompt_2sigma_EB_vtx = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_3sigma_EB_vtx = new TGraph();       TGraph* gr_norm_PU200_nonprompt_3sigma_EB_vtx = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_4sigma_EB_vtx = new TGraph();       TGraph* gr_norm_PU200_nonprompt_4sigma_EB_vtx = new TGraph();
      // noPU
  TGraph* gr_eff_noPU_nonprompt_EB = new TGraph();              TGraph* gr_norm_noPU_nonprompt_EB = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_40_EB = new TGraph();           TGraph* gr_norm_noPU_nonprompt_40_EB = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_60_EB = new TGraph();           TGraph* gr_norm_noPU_nonprompt_60_EB = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_80_EB = new TGraph();           TGraph* gr_norm_noPU_nonprompt_80_EB = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_100_EB = new TGraph();          TGraph* gr_norm_noPU_nonprompt_100_EB = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_2sigma_EB = new TGraph();       TGraph* gr_norm_noPU_nonprompt_2sigma_EB = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_3sigma_EB = new TGraph();       TGraph* gr_norm_noPU_nonprompt_3sigma_EB = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_4sigma_EB = new TGraph();       TGraph* gr_norm_noPU_nonprompt_4sigma_EB = new TGraph();

  TGraph* gr_eff_noPU_nonprompt_gen_EB = new TGraph();              TGraph* gr_norm_noPU_nonprompt_gen_EB = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_gen_2sigma_EB = new TGraph();       TGraph* gr_norm_noPU_nonprompt_gen_2sigma_EB = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_gen_3sigma_EB = new TGraph();       TGraph* gr_norm_noPU_nonprompt_gen_3sigma_EB = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_gen_4sigma_EB = new TGraph();       TGraph* gr_norm_noPU_nonprompt_gen_4sigma_EB = new TGraph();

    // Endcap region
      // PU200
  TGraph* gr_eff_PU200_nonprompt_EE = new TGraph();              TGraph* gr_norm_PU200_nonprompt_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_40_EE = new TGraph();           TGraph* gr_norm_PU200_nonprompt_40_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_60_EE = new TGraph();           TGraph* gr_norm_PU200_nonprompt_60_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_80_EE = new TGraph();           TGraph* gr_norm_PU200_nonprompt_80_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_100_EE = new TGraph();          TGraph* gr_norm_PU200_nonprompt_100_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_2sigma_EE = new TGraph();       TGraph* gr_norm_PU200_nonprompt_2sigma_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_3sigma_EE = new TGraph();       TGraph* gr_norm_PU200_nonprompt_3sigma_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_4sigma_EE = new TGraph();       TGraph* gr_norm_PU200_nonprompt_4sigma_EE = new TGraph();

  TGraph* gr_eff_PU200_nonprompt_gen_EE = new TGraph();              TGraph* gr_norm_PU200_nonprompt_gen_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_gen_2sigma_EE = new TGraph();       TGraph* gr_norm_PU200_nonprompt_gen_2sigma_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_gen_3sigma_EE = new TGraph();       TGraph* gr_norm_PU200_nonprompt_gen_3sigma_EE = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_gen_4sigma_EE = new TGraph();       TGraph* gr_norm_PU200_nonprompt_gen_4sigma_EE = new TGraph();
        // vtx
  TGraph* gr_eff_PU200_nonprompt_2sigma_EE_vtx = new TGraph();       TGraph* gr_norm_PU200_nonprompt_2sigma_EE_vtx = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_3sigma_EE_vtx = new TGraph();       TGraph* gr_norm_PU200_nonprompt_3sigma_EE_vtx = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_4sigma_EE_vtx = new TGraph();       TGraph* gr_norm_PU200_nonprompt_4sigma_EE_vtx = new TGraph();
      // noPU
  TGraph* gr_eff_noPU_nonprompt_EE = new TGraph();              TGraph* gr_norm_noPU_nonprompt_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_40_EE = new TGraph();           TGraph* gr_norm_noPU_nonprompt_40_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_60_EE = new TGraph();           TGraph* gr_norm_noPU_nonprompt_60_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_80_EE = new TGraph();           TGraph* gr_norm_noPU_nonprompt_80_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_100_EE = new TGraph();          TGraph* gr_norm_noPU_nonprompt_100_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_2sigma_EE = new TGraph();       TGraph* gr_norm_noPU_nonprompt_2sigma_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_3sigma_EE = new TGraph();       TGraph* gr_norm_noPU_nonprompt_3sigma_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_4sigma_EE = new TGraph();       TGraph* gr_norm_noPU_nonprompt_4sigma_EE = new TGraph();

  TGraph* gr_eff_noPU_nonprompt_gen_EE = new TGraph();              TGraph* gr_norm_noPU_nonprompt_gen_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_gen_2sigma_EE = new TGraph();       TGraph* gr_norm_noPU_nonprompt_gen_2sigma_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_gen_3sigma_EE = new TGraph();       TGraph* gr_norm_noPU_nonprompt_gen_3sigma_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_gen_4sigma_EE = new TGraph();       TGraph* gr_norm_noPU_nonprompt_gen_4sigma_EE = new TGraph();


  for(unsigned int i=0; i<1001; i++) {  // Store efficiency up to rel. iso. cut==4
  //for(unsigned int i=0; i<63; i++) {      // Store efficiency up to rel. iso. cut==0.25
  //for(unsigned int i=0; i<4; i++) {      // 
  // Prompt
    // Barrel region
      // efficiency
        // PU200
    gr_eff_PU200_prompt_EB->SetPoint(gr_eff_PU200_prompt_EB->GetN(), 0.004*i, prompt_eff_PU200_EB.at(i+1));
    gr_eff_PU200_prompt_40_EB->SetPoint(gr_eff_PU200_prompt_40_EB->GetN(), 0.004*i, prompt_eff_PU200_40_EB.at(i+1));
    gr_eff_PU200_prompt_60_EB->SetPoint(gr_eff_PU200_prompt_60_EB->GetN(), 0.004*i, prompt_eff_PU200_60_EB.at(i+1));
    gr_eff_PU200_prompt_80_EB->SetPoint(gr_eff_PU200_prompt_80_EB->GetN(), 0.004*i, prompt_eff_PU200_80_EB.at(i+1));
    gr_eff_PU200_prompt_100_EB->SetPoint(gr_eff_PU200_prompt_100_EB->GetN(), 0.004*i, prompt_eff_PU200_100_EB.at(i+1));
    gr_eff_PU200_prompt_2sigma_EB->SetPoint(gr_eff_PU200_prompt_2sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EB.at(i+1));
    gr_eff_PU200_prompt_3sigma_EB->SetPoint(gr_eff_PU200_prompt_3sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EB.at(i+1));
    gr_eff_PU200_prompt_4sigma_EB->SetPoint(gr_eff_PU200_prompt_4sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EB.at(i+1));

    gr_eff_PU200_prompt_gen_EB->SetPoint(gr_eff_PU200_prompt_gen_EB->GetN(), 0.004*i, prompt_eff_PU200_gen_EB.at(i+1));
    gr_eff_PU200_prompt_gen_2sigma_EB->SetPoint(gr_eff_PU200_prompt_gen_2sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_gen_2sigma_EB.at(i+1));
    gr_eff_PU200_prompt_gen_3sigma_EB->SetPoint(gr_eff_PU200_prompt_gen_3sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_gen_3sigma_EB.at(i+1));
    gr_eff_PU200_prompt_gen_4sigma_EB->SetPoint(gr_eff_PU200_prompt_gen_4sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_gen_4sigma_EB.at(i+1));
          // vtx
    gr_eff_PU200_prompt_2sigma_EB_vtx->SetPoint(gr_eff_PU200_prompt_2sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EB_vtx.at(i+1));
    gr_eff_PU200_prompt_3sigma_EB_vtx->SetPoint(gr_eff_PU200_prompt_3sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EB_vtx.at(i+1));
    gr_eff_PU200_prompt_4sigma_EB_vtx->SetPoint(gr_eff_PU200_prompt_4sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EB_vtx.at(i+1));
        // noPU
    gr_eff_noPU_prompt_EB->SetPoint(gr_eff_noPU_prompt_EB->GetN(), 0.004*i, prompt_eff_noPU_EB.at(i+1));
    gr_eff_noPU_prompt_40_EB->SetPoint(gr_eff_noPU_prompt_40_EB->GetN(), 0.004*i, prompt_eff_noPU_40_EB.at(i+1));
    gr_eff_noPU_prompt_60_EB->SetPoint(gr_eff_noPU_prompt_60_EB->GetN(), 0.004*i, prompt_eff_noPU_60_EB.at(i+1));
    gr_eff_noPU_prompt_80_EB->SetPoint(gr_eff_noPU_prompt_80_EB->GetN(), 0.004*i, prompt_eff_noPU_80_EB.at(i+1));
    gr_eff_noPU_prompt_100_EB->SetPoint(gr_eff_noPU_prompt_100_EB->GetN(), 0.004*i, prompt_eff_noPU_100_EB.at(i+1));
    gr_eff_noPU_prompt_2sigma_EB->SetPoint(gr_eff_noPU_prompt_2sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_2sigma_EB.at(i+1));
    gr_eff_noPU_prompt_3sigma_EB->SetPoint(gr_eff_noPU_prompt_3sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_3sigma_EB.at(i+1));
    gr_eff_noPU_prompt_4sigma_EB->SetPoint(gr_eff_noPU_prompt_4sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_4sigma_EB.at(i+1));

    gr_eff_noPU_prompt_gen_EB->SetPoint(gr_eff_noPU_prompt_gen_EB->GetN(), 0.004*i, prompt_eff_noPU_gen_EB.at(i+1));
    gr_eff_noPU_prompt_gen_2sigma_EB->SetPoint(gr_eff_noPU_prompt_gen_2sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_gen_2sigma_EB.at(i+1));
    gr_eff_noPU_prompt_gen_3sigma_EB->SetPoint(gr_eff_noPU_prompt_gen_3sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_gen_3sigma_EB.at(i+1));
    gr_eff_noPU_prompt_gen_4sigma_EB->SetPoint(gr_eff_noPU_prompt_gen_4sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_gen_4sigma_EB.at(i+1));

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
        // PU200
    gr_eff_PU200_prompt_EE->SetPoint(gr_eff_PU200_prompt_EE->GetN(), 0.004*i, prompt_eff_PU200_EE.at(i+1));
    gr_eff_PU200_prompt_40_EE->SetPoint(gr_eff_PU200_prompt_40_EE->GetN(), 0.004*i, prompt_eff_PU200_40_EE.at(i+1));
    gr_eff_PU200_prompt_60_EE->SetPoint(gr_eff_PU200_prompt_60_EE->GetN(), 0.004*i, prompt_eff_PU200_60_EE.at(i+1));
    gr_eff_PU200_prompt_80_EE->SetPoint(gr_eff_PU200_prompt_80_EE->GetN(), 0.004*i, prompt_eff_PU200_80_EE.at(i+1));
    gr_eff_PU200_prompt_100_EE->SetPoint(gr_eff_PU200_prompt_100_EE->GetN(), 0.004*i, prompt_eff_PU200_100_EE.at(i+1));
    gr_eff_PU200_prompt_2sigma_EE->SetPoint(gr_eff_PU200_prompt_2sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EE.at(i+1));
    gr_eff_PU200_prompt_3sigma_EE->SetPoint(gr_eff_PU200_prompt_3sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EE.at(i+1));
    gr_eff_PU200_prompt_4sigma_EE->SetPoint(gr_eff_PU200_prompt_4sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EE.at(i+1));

    gr_eff_PU200_prompt_gen_EE->SetPoint(gr_eff_PU200_prompt_gen_EE->GetN(), 0.004*i, prompt_eff_PU200_gen_EE.at(i+1));
    gr_eff_PU200_prompt_gen_2sigma_EE->SetPoint(gr_eff_PU200_prompt_gen_2sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_gen_2sigma_EE.at(i+1));
    gr_eff_PU200_prompt_gen_3sigma_EE->SetPoint(gr_eff_PU200_prompt_gen_3sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_gen_3sigma_EE.at(i+1));
    gr_eff_PU200_prompt_gen_4sigma_EE->SetPoint(gr_eff_PU200_prompt_gen_4sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_gen_4sigma_EE.at(i+1));
          // vtx
    gr_eff_PU200_prompt_2sigma_EE_vtx->SetPoint(gr_eff_PU200_prompt_2sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EE_vtx.at(i+1));
    gr_eff_PU200_prompt_3sigma_EE_vtx->SetPoint(gr_eff_PU200_prompt_3sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EE_vtx.at(i+1));
    gr_eff_PU200_prompt_4sigma_EE_vtx->SetPoint(gr_eff_PU200_prompt_4sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EE_vtx.at(i+1));
        // noPU
    gr_eff_noPU_prompt_EE->SetPoint(gr_eff_noPU_prompt_EE->GetN(), 0.004*i, prompt_eff_noPU_EE.at(i+1));
    gr_eff_noPU_prompt_40_EE->SetPoint(gr_eff_noPU_prompt_40_EE->GetN(), 0.004*i, prompt_eff_noPU_40_EE.at(i+1));
    gr_eff_noPU_prompt_60_EE->SetPoint(gr_eff_noPU_prompt_60_EE->GetN(), 0.004*i, prompt_eff_noPU_60_EE.at(i+1));
    gr_eff_noPU_prompt_80_EE->SetPoint(gr_eff_noPU_prompt_80_EE->GetN(), 0.004*i, prompt_eff_noPU_80_EE.at(i+1));
    gr_eff_noPU_prompt_100_EE->SetPoint(gr_eff_noPU_prompt_100_EE->GetN(), 0.004*i, prompt_eff_noPU_100_EE.at(i+1));
    gr_eff_noPU_prompt_2sigma_EE->SetPoint(gr_eff_noPU_prompt_2sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_2sigma_EE.at(i+1));
    gr_eff_noPU_prompt_3sigma_EE->SetPoint(gr_eff_noPU_prompt_3sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_3sigma_EE.at(i+1));
    gr_eff_noPU_prompt_4sigma_EE->SetPoint(gr_eff_noPU_prompt_4sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_4sigma_EE.at(i+1));

    gr_eff_noPU_prompt_gen_EE->SetPoint(gr_eff_noPU_prompt_gen_EE->GetN(), 0.004*i, prompt_eff_noPU_gen_EE.at(i+1));
    gr_eff_noPU_prompt_gen_2sigma_EE->SetPoint(gr_eff_noPU_prompt_gen_2sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_gen_2sigma_EE.at(i+1));
    gr_eff_noPU_prompt_gen_3sigma_EE->SetPoint(gr_eff_noPU_prompt_gen_3sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_gen_3sigma_EE.at(i+1));
    gr_eff_noPU_prompt_gen_4sigma_EE->SetPoint(gr_eff_noPU_prompt_gen_4sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_gen_4sigma_EE.at(i+1));

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
  for(unsigned int i=0; i<1001; i++) {      //
  //for(unsigned int i=0; i<63; i++) {      // Store efficiency up to rel. iso. cut==0.25
  // Nonprompt
    // Barrel region
      // efficiency
        // PU200
    gr_eff_PU200_nonprompt_EB->SetPoint(gr_eff_PU200_nonprompt_EB->GetN(), 0.004*i, nonprompt_eff_PU200_EB.at(i+1));
    gr_eff_PU200_nonprompt_40_EB->SetPoint(gr_eff_PU200_nonprompt_40_EB->GetN(), 0.004*i, nonprompt_eff_PU200_40_EB.at(i+1));
    gr_eff_PU200_nonprompt_60_EB->SetPoint(gr_eff_PU200_nonprompt_60_EB->GetN(), 0.004*i, nonprompt_eff_PU200_60_EB.at(i+1));
    gr_eff_PU200_nonprompt_80_EB->SetPoint(gr_eff_PU200_nonprompt_80_EB->GetN(), 0.004*i, nonprompt_eff_PU200_80_EB.at(i+1));
    gr_eff_PU200_nonprompt_100_EB->SetPoint(gr_eff_PU200_nonprompt_100_EB->GetN(), 0.004*i, nonprompt_eff_PU200_100_EB.at(i+1));
    gr_eff_PU200_nonprompt_2sigma_EB->SetPoint(gr_eff_PU200_nonprompt_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EB.at(i+1));
    gr_eff_PU200_nonprompt_3sigma_EB->SetPoint(gr_eff_PU200_nonprompt_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EB.at(i+1));
    gr_eff_PU200_nonprompt_4sigma_EB->SetPoint(gr_eff_PU200_nonprompt_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EB.at(i+1));

    gr_eff_PU200_nonprompt_gen_EB->SetPoint(gr_eff_PU200_nonprompt_gen_EB->GetN(), 0.004*i, nonprompt_eff_PU200_gen_EB.at(i+1));
    gr_eff_PU200_nonprompt_gen_2sigma_EB->SetPoint(gr_eff_PU200_nonprompt_gen_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_gen_2sigma_EB.at(i+1));
    gr_eff_PU200_nonprompt_gen_3sigma_EB->SetPoint(gr_eff_PU200_nonprompt_gen_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_gen_3sigma_EB.at(i+1));
    gr_eff_PU200_nonprompt_gen_4sigma_EB->SetPoint(gr_eff_PU200_nonprompt_gen_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_gen_4sigma_EB.at(i+1));
          // vtx
    gr_eff_PU200_nonprompt_2sigma_EB_vtx->SetPoint(gr_eff_PU200_nonprompt_2sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EB_vtx.at(i+1));
    gr_eff_PU200_nonprompt_3sigma_EB_vtx->SetPoint(gr_eff_PU200_nonprompt_3sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EB_vtx.at(i+1));
    gr_eff_PU200_nonprompt_4sigma_EB_vtx->SetPoint(gr_eff_PU200_nonprompt_4sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EB_vtx.at(i+1));
        // noPU
    gr_eff_noPU_nonprompt_EB->SetPoint(gr_eff_noPU_nonprompt_EB->GetN(), 0.004*i, nonprompt_eff_noPU_EB.at(i+1));
    gr_eff_noPU_nonprompt_40_EB->SetPoint(gr_eff_noPU_nonprompt_40_EB->GetN(), 0.004*i, nonprompt_eff_noPU_40_EB.at(i+1));
    gr_eff_noPU_nonprompt_60_EB->SetPoint(gr_eff_noPU_nonprompt_60_EB->GetN(), 0.004*i, nonprompt_eff_noPU_60_EB.at(i+1));
    gr_eff_noPU_nonprompt_80_EB->SetPoint(gr_eff_noPU_nonprompt_80_EB->GetN(), 0.004*i, nonprompt_eff_noPU_80_EB.at(i+1));
    gr_eff_noPU_nonprompt_100_EB->SetPoint(gr_eff_noPU_nonprompt_100_EB->GetN(), 0.004*i, nonprompt_eff_noPU_100_EB.at(i+1));
    gr_eff_noPU_nonprompt_2sigma_EB->SetPoint(gr_eff_noPU_nonprompt_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_2sigma_EB.at(i+1));
    gr_eff_noPU_nonprompt_3sigma_EB->SetPoint(gr_eff_noPU_nonprompt_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_3sigma_EB.at(i+1));
    gr_eff_noPU_nonprompt_4sigma_EB->SetPoint(gr_eff_noPU_nonprompt_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_4sigma_EB.at(i+1));

    gr_eff_noPU_nonprompt_gen_EB->SetPoint(gr_eff_noPU_nonprompt_gen_EB->GetN(), 0.004*i, nonprompt_eff_noPU_gen_EB.at(i+1));
    gr_eff_noPU_nonprompt_gen_2sigma_EB->SetPoint(gr_eff_noPU_nonprompt_gen_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_gen_2sigma_EB.at(i+1));
    gr_eff_noPU_nonprompt_gen_3sigma_EB->SetPoint(gr_eff_noPU_nonprompt_gen_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_gen_3sigma_EB.at(i+1));
    gr_eff_noPU_nonprompt_gen_4sigma_EB->SetPoint(gr_eff_noPU_nonprompt_gen_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_gen_4sigma_EB.at(i+1));

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
        // PU200
    gr_eff_PU200_nonprompt_EE->SetPoint(gr_eff_PU200_nonprompt_EE->GetN(), 0.004*i, nonprompt_eff_PU200_EE.at(i+1));
    gr_eff_PU200_nonprompt_40_EE->SetPoint(gr_eff_PU200_nonprompt_40_EE->GetN(), 0.004*i, nonprompt_eff_PU200_40_EE.at(i+1));
    gr_eff_PU200_nonprompt_60_EE->SetPoint(gr_eff_PU200_nonprompt_60_EE->GetN(), 0.004*i, nonprompt_eff_PU200_60_EE.at(i+1));
    gr_eff_PU200_nonprompt_80_EE->SetPoint(gr_eff_PU200_nonprompt_80_EE->GetN(), 0.004*i, nonprompt_eff_PU200_80_EE.at(i+1));
    gr_eff_PU200_nonprompt_100_EE->SetPoint(gr_eff_PU200_nonprompt_100_EE->GetN(), 0.004*i, nonprompt_eff_PU200_100_EE.at(i+1));
    gr_eff_PU200_nonprompt_2sigma_EE->SetPoint(gr_eff_PU200_nonprompt_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EE.at(i+1));
    gr_eff_PU200_nonprompt_3sigma_EE->SetPoint(gr_eff_PU200_nonprompt_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EE.at(i+1));
    gr_eff_PU200_nonprompt_4sigma_EE->SetPoint(gr_eff_PU200_nonprompt_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EE.at(i+1));

    gr_eff_PU200_nonprompt_gen_EE->SetPoint(gr_eff_PU200_nonprompt_gen_EE->GetN(), 0.004*i, nonprompt_eff_PU200_gen_EE.at(i+1));
    gr_eff_PU200_nonprompt_gen_2sigma_EE->SetPoint(gr_eff_PU200_nonprompt_gen_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_gen_2sigma_EE.at(i+1));
    gr_eff_PU200_nonprompt_gen_3sigma_EE->SetPoint(gr_eff_PU200_nonprompt_gen_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_gen_3sigma_EE.at(i+1));
    gr_eff_PU200_nonprompt_gen_4sigma_EE->SetPoint(gr_eff_PU200_nonprompt_gen_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_gen_4sigma_EE.at(i+1));
          // vtx
    gr_eff_PU200_nonprompt_2sigma_EE_vtx->SetPoint(gr_eff_PU200_nonprompt_2sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EE_vtx.at(i+1));
    gr_eff_PU200_nonprompt_3sigma_EE_vtx->SetPoint(gr_eff_PU200_nonprompt_3sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EE_vtx.at(i+1));
    gr_eff_PU200_nonprompt_4sigma_EE_vtx->SetPoint(gr_eff_PU200_nonprompt_4sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EE_vtx.at(i+1));
        // noPU
    gr_eff_noPU_nonprompt_EE->SetPoint(gr_eff_noPU_nonprompt_EE->GetN(), 0.004*i, nonprompt_eff_noPU_EE.at(i+1));
    gr_eff_noPU_nonprompt_40_EE->SetPoint(gr_eff_noPU_nonprompt_40_EE->GetN(), 0.004*i, nonprompt_eff_noPU_40_EE.at(i+1));
    gr_eff_noPU_nonprompt_60_EE->SetPoint(gr_eff_noPU_nonprompt_60_EE->GetN(), 0.004*i, nonprompt_eff_noPU_60_EE.at(i+1));
    gr_eff_noPU_nonprompt_80_EE->SetPoint(gr_eff_noPU_nonprompt_80_EE->GetN(), 0.004*i, nonprompt_eff_noPU_80_EE.at(i+1));
    gr_eff_noPU_nonprompt_100_EE->SetPoint(gr_eff_noPU_nonprompt_100_EE->GetN(), 0.004*i, nonprompt_eff_noPU_100_EE.at(i+1));
    gr_eff_noPU_nonprompt_2sigma_EE->SetPoint(gr_eff_noPU_nonprompt_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_2sigma_EE.at(i+1));
    gr_eff_noPU_nonprompt_3sigma_EE->SetPoint(gr_eff_noPU_nonprompt_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_3sigma_EE.at(i+1));
    gr_eff_noPU_nonprompt_4sigma_EE->SetPoint(gr_eff_noPU_nonprompt_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_4sigma_EE.at(i+1));

    gr_eff_noPU_nonprompt_gen_EE->SetPoint(gr_eff_noPU_nonprompt_gen_EE->GetN(), 0.004*i, nonprompt_eff_noPU_gen_EE.at(i+1));
    gr_eff_noPU_nonprompt_gen_2sigma_EE->SetPoint(gr_eff_noPU_nonprompt_gen_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_gen_2sigma_EE.at(i+1));
    gr_eff_noPU_nonprompt_gen_3sigma_EE->SetPoint(gr_eff_noPU_nonprompt_gen_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_gen_3sigma_EE.at(i+1));
    gr_eff_noPU_nonprompt_gen_4sigma_EE->SetPoint(gr_eff_noPU_nonprompt_gen_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_gen_4sigma_EE.at(i+1));

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
  gr_eff_noPU_prompt_EB->SetLineColor(kGray+1); gr_eff_noPU_prompt_gen_EB->SetLineColor(kGray+1);
  gr_eff_PU200_prompt_EB->SetLineWidth(2); gr_eff_PU200_prompt_gen_EB->SetLineWidth(2);
  gr_eff_PU200_prompt_40_EB->SetLineWidth(2); gr_eff_PU200_prompt_60_EB->SetLineWidth(2); gr_eff_PU200_prompt_80_EB->SetLineWidth(2); gr_eff_PU200_prompt_100_EB->SetLineWidth(2);
  gr_eff_PU200_prompt_2sigma_EB->SetLineWidth(2); gr_eff_PU200_prompt_3sigma_EB->SetLineWidth(2); gr_eff_PU200_prompt_4sigma_EB->SetLineWidth(2);
  gr_eff_noPU_prompt_EB->SetLineWidth(2); gr_eff_noPU_prompt_gen_EB->SetLineWidth(2);
  gr_eff_PU200_prompt_gen_EB->SetLineStyle(7); gr_eff_noPU_prompt_gen_EB->SetLineStyle(7);
  gr_eff_PU200_prompt_gen_2sigma_EB->SetLineWidth(2); gr_eff_noPU_prompt_gen_2sigma_EB->SetLineWidth(2);
  gr_eff_PU200_prompt_gen_2sigma_EB->SetLineStyle(7); gr_eff_noPU_prompt_gen_2sigma_EB->SetLineStyle(7);
  gr_eff_PU200_prompt_gen_2sigma_EB->SetLineColor(kBlack); gr_eff_noPU_prompt_gen_2sigma_EB->SetLineColor(kGray+1);
  gr_eff_PU200_prompt_2sigma_EB_vtx->SetLineColor(kGreen+2); gr_eff_PU200_prompt_2sigma_EB_vtx->SetLineWidth(2);

  gr_eff_PU200_prompt_EB->SetMarkerStyle(8); gr_eff_PU200_prompt_EB->SetMarkerSize(1.5);
  gr_eff_noPU_prompt_EB->SetMarkerStyle(8); gr_eff_noPU_prompt_EB->SetMarkerSize(1.5); gr_eff_noPU_prompt_EB->SetMarkerColor(kGray+1);

      // normalization
  gr_norm_PU200_prompt_EB->SetTitle("Prompt muon in barrel region"); gr_norm_noPU_prompt_EB->SetTitle("Prompt muon in barrel region");
  gr_norm_PU200_prompt_EB->GetXaxis()->SetTitle("Relative isolation"); gr_norm_noPU_prompt_EB->GetXaxis()->SetTitle("Relative isolation");
  gr_norm_PU200_prompt_EB->GetYaxis()->SetTitle("Counts / 0.004"); gr_norm_noPU_prompt_EB->GetYaxis()->SetTitle("Counts / 0.004");
  gr_norm_PU200_prompt_EB->SetLineColor(kBlack); gr_norm_PU200_prompt_gen_EB->SetLineColor(kBlack);
  gr_norm_PU200_prompt_40_EB->SetLineColor(kRed); gr_norm_PU200_prompt_60_EB->SetLineColor(kGreen); gr_norm_PU200_prompt_80_EB->SetLineColor(kBlue); gr_norm_PU200_prompt_100_EB->SetLineColor(kMagenta);
  gr_norm_PU200_prompt_2sigma_EB->SetLineColor(kRed); gr_norm_PU200_prompt_3sigma_EB->SetLineColor(kGreen); gr_norm_PU200_prompt_4sigma_EB->SetLineColor(kBlue);
  gr_norm_noPU_prompt_EB->SetLineColor(kGray+1); gr_norm_noPU_prompt_gen_EB->SetLineColor(kGray+1);
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
  gr_eff_noPU_prompt_EE->SetLineColor(kGray+1); gr_eff_noPU_prompt_gen_EE->SetLineColor(kGray+1);
  gr_eff_PU200_prompt_EE->SetLineWidth(2); gr_eff_PU200_prompt_gen_EE->SetLineWidth(2);
  gr_eff_PU200_prompt_40_EE->SetLineWidth(2); gr_eff_PU200_prompt_60_EE->SetLineWidth(2); gr_eff_PU200_prompt_80_EE->SetLineWidth(2); gr_eff_PU200_prompt_100_EE->SetLineWidth(2);
  gr_eff_PU200_prompt_2sigma_EE->SetLineWidth(2); gr_eff_PU200_prompt_3sigma_EE->SetLineWidth(2); gr_eff_PU200_prompt_4sigma_EE->SetLineWidth(2);
  gr_eff_noPU_prompt_EE->SetLineWidth(2); gr_eff_noPU_prompt_gen_EE->SetLineWidth(2);
  gr_eff_PU200_prompt_gen_EE->SetLineStyle(7); gr_eff_noPU_prompt_gen_EE->SetLineStyle(7);
  gr_eff_PU200_prompt_gen_2sigma_EE->SetLineWidth(2); gr_eff_noPU_prompt_gen_2sigma_EE->SetLineWidth(2);
  gr_eff_PU200_prompt_gen_2sigma_EE->SetLineStyle(7); gr_eff_noPU_prompt_gen_2sigma_EE->SetLineStyle(7);
  gr_eff_PU200_prompt_gen_2sigma_EE->SetLineColor(kBlack); gr_eff_noPU_prompt_gen_2sigma_EE->SetLineColor(kGray+1);
  gr_eff_PU200_prompt_2sigma_EE_vtx->SetLineColor(kGreen+2); gr_eff_PU200_prompt_2sigma_EE_vtx->SetLineWidth(2);

  gr_eff_PU200_prompt_EE->SetMarkerStyle(8); gr_eff_PU200_prompt_EE->SetMarkerSize(1.5);
  gr_eff_noPU_prompt_EE->SetMarkerStyle(8); gr_eff_noPU_prompt_EE->SetMarkerSize(1.5); gr_eff_noPU_prompt_EE->SetMarkerColor(kGray+1);

      // normalization
  gr_norm_PU200_prompt_EE->SetTitle("Prompt muon in endcap region"); gr_norm_noPU_prompt_EE->SetTitle("Prompt muon in endcap region");
  gr_norm_PU200_prompt_EE->GetXaxis()->SetTitle("Relative isolation"); gr_norm_noPU_prompt_EE->GetXaxis()->SetTitle("Relative isolation");
  gr_norm_PU200_prompt_EE->GetYaxis()->SetTitle("Counts / 0.004"); gr_norm_noPU_prompt_EE->GetYaxis()->SetTitle("Counts / 0.004");
  gr_norm_PU200_prompt_EE->SetLineColor(kBlack); gr_norm_PU200_prompt_gen_EE->SetLineColor(kBlack);
  gr_norm_PU200_prompt_40_EE->SetLineColor(kRed); gr_norm_PU200_prompt_60_EE->SetLineColor(kGreen); gr_norm_PU200_prompt_80_EE->SetLineColor(kBlue); gr_norm_PU200_prompt_100_EE->SetLineColor(kMagenta);
  gr_norm_PU200_prompt_2sigma_EE->SetLineColor(kRed); gr_norm_PU200_prompt_3sigma_EE->SetLineColor(kGreen); gr_norm_PU200_prompt_4sigma_EE->SetLineColor(kBlue);
  gr_norm_noPU_prompt_EE->SetLineColor(kGray+1); gr_norm_noPU_prompt_gen_EE->SetLineColor(kGray+1);
  gr_norm_PU200_prompt_EE->SetLineWidth(2); gr_norm_PU200_prompt_gen_EE->SetLineWidth(2);
  gr_norm_PU200_prompt_40_EE->SetLineWidth(2); gr_norm_PU200_prompt_60_EE->SetLineWidth(2); gr_norm_PU200_prompt_80_EE->SetLineWidth(2); gr_norm_PU200_prompt_100_EE->SetLineWidth(2);
  gr_norm_PU200_prompt_2sigma_EE->SetLineWidth(2); gr_norm_PU200_prompt_3sigma_EE->SetLineWidth(2); gr_norm_PU200_prompt_4sigma_EE->SetLineWidth(2);
  gr_norm_noPU_prompt_EE->SetLineWidth(2); gr_norm_noPU_prompt_gen_EE->SetLineWidth(2);
  gr_norm_PU200_prompt_gen_EE->SetLineStyle(7); gr_norm_noPU_prompt_gen_EE->SetLineStyle(7);
  gr_norm_noPU_prompt_EE->SetMinimum(0);
  // Nonprompt
    // Barrel region
      // efficiency
  gr_eff_noPU_nonprompt_gen_2sigma_EB->SetTitle("Non-prompt muon in barrel region"); gr_eff_noPU_nonprompt_gen_2sigma_EB->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_nonprompt_gen_2sigma_EB->GetYaxis()->SetTitle("Isolation efficiency");
  gr_eff_PU200_nonprompt_EB->SetTitle("Non-prompt muon in barrel region"); gr_eff_noPU_nonprompt_EB->SetTitle("Non-prompt muon in barrel region");
  gr_eff_PU200_nonprompt_EB->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_nonprompt_EB->GetXaxis()->SetTitle("Relative isolation");
  gr_eff_PU200_nonprompt_EB->GetYaxis()->SetTitle("Isolation efficiency"); gr_eff_noPU_nonprompt_EB->GetYaxis()->SetTitle("Isolation efficiency");
  gr_eff_PU200_nonprompt_EB->SetLineColor(kBlack); gr_eff_PU200_nonprompt_gen_EB->SetLineColor(kBlack);
  gr_eff_PU200_nonprompt_40_EB->SetLineColor(kRed); gr_eff_PU200_nonprompt_60_EB->SetLineColor(kGreen); gr_eff_PU200_nonprompt_80_EB->SetLineColor(kBlue); gr_eff_PU200_nonprompt_100_EB->SetLineColor(kMagenta);
  gr_eff_PU200_nonprompt_2sigma_EB->SetLineColor(kRed); gr_eff_PU200_nonprompt_3sigma_EB->SetLineColor(kGreen); gr_eff_PU200_nonprompt_4sigma_EB->SetLineColor(kBlue);
  gr_eff_noPU_nonprompt_EB->SetLineColor(kGray+1); gr_eff_noPU_nonprompt_gen_EB->SetLineColor(kGray+1);
  gr_eff_noPU_nonprompt_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_gen_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_40_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_60_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_80_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_100_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_2sigma_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_3sigma_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_4sigma_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_EB->SetLineWidth(2); gr_eff_noPU_nonprompt_gen_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_EB->SetLineStyle(7); gr_eff_noPU_nonprompt_gen_EB->SetLineStyle(7);
  gr_eff_noPU_nonprompt_2sigma_EB->SetLineColor(kBlue); gr_eff_noPU_nonprompt_2sigma_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_2sigma_EB->SetLineWidth(2); gr_eff_noPU_nonprompt_gen_2sigma_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_2sigma_EB->SetLineStyle(7); gr_eff_noPU_nonprompt_gen_2sigma_EB->SetLineStyle(7);
  gr_eff_PU200_nonprompt_gen_2sigma_EB->SetLineColor(kBlack); gr_eff_noPU_nonprompt_gen_2sigma_EB->SetLineColor(kGray+1);
  gr_eff_PU200_nonprompt_2sigma_EB_vtx->SetLineColor(kGreen+2); gr_eff_PU200_nonprompt_2sigma_EB_vtx->SetLineWidth(2);

  gr_eff_PU200_nonprompt_EB->SetMarkerStyle(8); gr_eff_PU200_nonprompt_EB->SetMarkerSize(1.5);
  gr_eff_noPU_nonprompt_EB->SetMarkerStyle(8); gr_eff_noPU_nonprompt_EB->SetMarkerSize(1.5); gr_eff_noPU_nonprompt_EB->SetMarkerColor(kGray+1);

      // normalization
  gr_norm_PU200_nonprompt_EB->SetTitle("Non-prompt muon in barrel region"); gr_norm_noPU_nonprompt_EB->SetTitle("Non-prompt muon in barrel region");
  gr_norm_PU200_nonprompt_EB->GetXaxis()->SetTitle("Relative isolation"); gr_norm_noPU_nonprompt_EB->GetXaxis()->SetTitle("Relative isolation");
  gr_norm_PU200_nonprompt_EB->GetYaxis()->SetTitle("Counts / 0.004"); gr_norm_noPU_nonprompt_EB->GetYaxis()->SetTitle("Counts / 0.004");
  gr_norm_PU200_nonprompt_EB->SetLineColor(kBlack); gr_norm_PU200_nonprompt_gen_EB->SetLineColor(kBlack);
  gr_norm_PU200_nonprompt_40_EB->SetLineColor(kRed); gr_norm_PU200_nonprompt_60_EB->SetLineColor(kGreen); gr_norm_PU200_nonprompt_80_EB->SetLineColor(kBlue); gr_norm_PU200_nonprompt_100_EB->SetLineColor(kMagenta);
  gr_norm_PU200_nonprompt_2sigma_EB->SetLineColor(kRed); gr_norm_PU200_nonprompt_3sigma_EB->SetLineColor(kGreen); gr_norm_PU200_nonprompt_4sigma_EB->SetLineColor(kBlue);
  gr_norm_noPU_nonprompt_EB->SetLineColor(kGray+1); gr_norm_noPU_nonprompt_gen_EB->SetLineColor(kGray+1);
  gr_norm_PU200_nonprompt_EB->SetLineWidth(2); gr_norm_PU200_nonprompt_gen_EB->SetLineWidth(2);
  gr_norm_PU200_nonprompt_40_EB->SetLineWidth(2); gr_norm_PU200_nonprompt_60_EB->SetLineWidth(2); gr_norm_PU200_nonprompt_80_EB->SetLineWidth(2); gr_norm_PU200_nonprompt_100_EB->SetLineWidth(2);
  gr_norm_PU200_nonprompt_2sigma_EB->SetLineWidth(2); gr_norm_PU200_nonprompt_3sigma_EB->SetLineWidth(2); gr_norm_PU200_nonprompt_4sigma_EB->SetLineWidth(2);
  gr_norm_PU200_nonprompt_gen_EB->SetLineWidth(2); gr_norm_noPU_nonprompt_gen_EB->SetLineWidth(2);
  gr_norm_PU200_nonprompt_gen_EB->SetLineStyle(7); gr_norm_noPU_nonprompt_gen_EB->SetLineStyle(7);
  gr_norm_noPU_nonprompt_EB->SetMinimum(0);
    // Endcap region
      // efficiency
  gr_eff_noPU_nonprompt_gen_2sigma_EE->SetTitle("Non-prompt muon in endcap region"); gr_eff_noPU_nonprompt_gen_2sigma_EE->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_nonprompt_gen_2sigma_EE->GetYaxis()->SetTitle("Isolation efficiency");
  gr_eff_PU200_nonprompt_EE->SetTitle("Non-prompt muon in endcap region"); gr_eff_noPU_nonprompt_EE->SetTitle("Non-prompt muon in endcap region");
  gr_eff_PU200_nonprompt_EE->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_nonprompt_EE->GetXaxis()->SetTitle("Relative isolation");
  gr_eff_PU200_nonprompt_EE->GetYaxis()->SetTitle("Isolation efficiency"); gr_eff_noPU_nonprompt_EE->GetYaxis()->SetTitle("Isolation efficiency");
  gr_eff_PU200_nonprompt_EE->SetLineColor(kBlack); gr_eff_PU200_nonprompt_gen_EE->SetLineColor(kBlack);
  gr_eff_PU200_nonprompt_40_EE->SetLineColor(kRed); gr_eff_PU200_nonprompt_60_EE->SetLineColor(kGreen); gr_eff_PU200_nonprompt_80_EE->SetLineColor(kBlue); gr_eff_PU200_nonprompt_100_EE->SetLineColor(kMagenta);
  gr_eff_PU200_nonprompt_2sigma_EE->SetLineColor(kRed); gr_eff_PU200_nonprompt_3sigma_EE->SetLineColor(kGreen); gr_eff_PU200_nonprompt_4sigma_EE->SetLineColor(kBlue);
  gr_eff_noPU_nonprompt_EE->SetLineColor(kGray+1); gr_eff_noPU_nonprompt_gen_EE->SetLineColor(kGray+1);
  gr_eff_noPU_nonprompt_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_gen_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_40_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_60_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_80_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_100_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_2sigma_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_3sigma_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_4sigma_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_EE->SetLineWidth(2); gr_eff_noPU_nonprompt_gen_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_EE->SetLineStyle(7); gr_eff_noPU_nonprompt_gen_EE->SetLineStyle(7);
  gr_eff_PU200_nonprompt_gen_2sigma_EE->SetLineWidth(2); gr_eff_noPU_nonprompt_gen_2sigma_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_2sigma_EE->SetLineStyle(7); gr_eff_noPU_nonprompt_gen_2sigma_EE->SetLineStyle(7);
  gr_eff_PU200_nonprompt_gen_2sigma_EE->SetLineColor(kBlack); gr_eff_noPU_nonprompt_gen_2sigma_EE->SetLineColor(kGray+1);
  gr_eff_PU200_nonprompt_2sigma_EE_vtx->SetLineColor(kGreen+2); gr_eff_PU200_nonprompt_2sigma_EE_vtx->SetLineWidth(2);

  gr_eff_PU200_nonprompt_EE->SetMarkerStyle(8); gr_eff_PU200_nonprompt_EE->SetMarkerSize(1.5);
  gr_eff_noPU_nonprompt_EE->SetMarkerStyle(8); gr_eff_noPU_nonprompt_EE->SetMarkerSize(1.5); gr_eff_noPU_nonprompt_EE->SetMarkerColor(kGray+1);

      // normalization
  gr_norm_PU200_nonprompt_EE->SetTitle("Non-prompt muon in endcap region"); gr_norm_noPU_nonprompt_EE->SetTitle("Non-prompt muon in endcap region");
  gr_norm_PU200_nonprompt_EE->GetXaxis()->SetTitle("Relative isolation"); gr_norm_noPU_nonprompt_EE->GetXaxis()->SetTitle("Relative isolation");
  gr_norm_PU200_nonprompt_EE->GetYaxis()->SetTitle("Counts / 0.004"); gr_norm_noPU_nonprompt_EE->GetYaxis()->SetTitle("Counts / 0.004");
  gr_norm_PU200_nonprompt_EE->SetLineColor(kBlack); gr_norm_PU200_nonprompt_gen_EE->SetLineColor(kBlack);
  gr_norm_PU200_nonprompt_40_EE->SetLineColor(kRed); gr_norm_PU200_nonprompt_60_EE->SetLineColor(kGreen); gr_norm_PU200_nonprompt_80_EE->SetLineColor(kBlue); gr_norm_PU200_nonprompt_100_EE->SetLineColor(kMagenta);
  gr_norm_PU200_nonprompt_2sigma_EE->SetLineColor(kRed); gr_norm_PU200_nonprompt_3sigma_EE->SetLineColor(kGreen); gr_norm_PU200_nonprompt_4sigma_EE->SetLineColor(kBlue);
  gr_norm_noPU_nonprompt_EE->SetLineColor(kGray+1); gr_norm_noPU_nonprompt_gen_EE->SetLineColor(kGray+1);
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
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_noPU_prompt_gen_EB, "gen noPU");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_40_EB, "40ps PU200");
//  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_60_EB, "60ps PU200");
//  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_80_EB, "80ps PU200");
//  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_100_EB, "100ps PU200");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_noPU_prompt_EB, "no MTD noPU");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_EB, "no MTD PU200");
  leg_eff_prompt_EB_dt->SetTextSize(0.03);
  TLegend* leg_eff_prompt_EB_sigma = new TLegend(0.62, 0.13, 0.88, 0.33);
  leg_eff_prompt_EB_sigma->SetMargin(0.2);
//  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_gen_EB, "gen PU200");
//  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_noPU_prompt_gen_EB, "gen noPU");
//  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_2sigma_EB_vtx, "2sigma PU200 (PV, track)");
//  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_2sigma_EB, "2sigma PU200 (muon, track)");
  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_EB, "no MTD PU200");
  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_noPU_prompt_EB, "no MTD noPU");
  leg_eff_prompt_EB_sigma->SetTextSize(0.03);
      // normalization
  TLegend* leg_norm_prompt_EB_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_prompt_EB_dt->AddEntry(gr_norm_PU200_prompt_gen_EB, "gen PU200");
  leg_norm_prompt_EB_dt->AddEntry(gr_norm_PU200_prompt_40_EB, "40ps PU200");
  leg_norm_prompt_EB_dt->AddEntry(gr_norm_PU200_prompt_EB, "no MTD PU200");
  leg_norm_prompt_EB_dt->AddEntry(gr_norm_noPU_prompt_EB, "no MTD noPU");
  leg_norm_prompt_EB_dt->AddEntry(gr_norm_noPU_prompt_gen_EB, "gen noPU");
  leg_norm_prompt_EB_dt->SetTextSize(0.03);
  TLegend* leg_norm_prompt_EB_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_prompt_EB_sigma->AddEntry(gr_norm_PU200_prompt_gen_EB, "gen PU200");
  leg_norm_prompt_EB_sigma->AddEntry(gr_norm_PU200_prompt_2sigma_EB, "2sigma PU200");
  leg_norm_prompt_EB_sigma->AddEntry(gr_norm_PU200_prompt_EB, "no MTD PU200");
  leg_norm_prompt_EB_sigma->AddEntry(gr_norm_noPU_prompt_gen_EB, "gen noPU");
  leg_norm_prompt_EB_sigma->AddEntry(gr_norm_noPU_prompt_EB, "no MTD noPU");
  leg_norm_prompt_EB_sigma->SetTextSize(0.03);
    // Endcap region
      // efficiency
  TLegend* leg_eff_prompt_EE_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_gen_EE, "gen PU200");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_noPU_prompt_gen_EE, "gen noPU");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_40_EE, "40ps PU200");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_EE, "no MTD PU200");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_noPU_prompt_EE, "no MTD noPU");
  leg_eff_prompt_EE_dt->SetTextSize(0.03);
  TLegend* leg_eff_prompt_EE_sigma = new TLegend(0.62, 0.13, 0.88, 0.33);
  leg_eff_prompt_EE_sigma->SetMargin(0.2);
//  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_gen_EE, "gen PU200");
//  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_noPU_prompt_gen_EE, "gen noPU");
//  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_2sigma_EE_vtx, "2sigma PU200 (PV, track)");
//  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_2sigma_EE, "2sigma PU200 (muon, track)");
  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_EE, "no MTD PU200");
  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_noPU_prompt_EE, "no MTD noPU");
  leg_eff_prompt_EE_sigma->SetTextSize(0.03);
  // test
  TLegend* leg_eff_prompt_EB_sigma_test = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_prompt_EB_sigma_test->AddEntry(gr_eff_PU200_prompt_gen_2sigma_EB, "gen 2sigma PU200");
  leg_eff_prompt_EB_sigma_test->AddEntry(gr_eff_noPU_prompt_gen_2sigma_EB, "gen 2sigma noPU");
  leg_eff_prompt_EB_sigma_test->AddEntry(gr_eff_PU200_prompt_2sigma_EB, "2sigma PU200");
  leg_eff_prompt_EB_sigma_test->AddEntry(gr_eff_PU200_prompt_EB, "no MTD PU200");
  leg_eff_prompt_EB_sigma_test->AddEntry(gr_eff_noPU_prompt_EB, "no MTD noPU");
  TLegend* leg_eff_prompt_EE_sigma_test = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_prompt_EE_sigma_test->AddEntry(gr_eff_PU200_prompt_gen_2sigma_EE, "gen 2sigma PU200");
  leg_eff_prompt_EE_sigma_test->AddEntry(gr_eff_PU200_prompt_2sigma_EE, "2sigma PU200");
  leg_eff_prompt_EE_sigma_test->AddEntry(gr_eff_PU200_prompt_EE, "no MTD PU200");
  leg_eff_prompt_EE_sigma_test->AddEntry(gr_eff_noPU_prompt_gen_2sigma_EE, "gen 2sigma noPU");
  leg_eff_prompt_EE_sigma_test->AddEntry(gr_eff_noPU_prompt_EE, "no MTD noPU");
  TLegend* leg_eff_nonprompt_EB_sigma_test = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_nonprompt_EB_sigma_test->AddEntry(gr_eff_PU200_nonprompt_gen_2sigma_EB, "gen 2sigma PU200");
  leg_eff_nonprompt_EB_sigma_test->AddEntry(gr_eff_PU200_nonprompt_2sigma_EB, "2sigma PU200");
  leg_eff_nonprompt_EB_sigma_test->AddEntry(gr_eff_PU200_nonprompt_EB, "no MTD PU200");
  leg_eff_nonprompt_EB_sigma_test->AddEntry(gr_eff_noPU_nonprompt_gen_2sigma_EB, "gen 2sigma noPU");
  leg_eff_nonprompt_EB_sigma_test->AddEntry(gr_eff_noPU_nonprompt_EB, "no MTD noPU");
  TLegend* leg_eff_nonprompt_EE_sigma_test = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_nonprompt_EE_sigma_test->AddEntry(gr_eff_PU200_nonprompt_gen_2sigma_EE, "gen 2sigma PU200");
  leg_eff_nonprompt_EE_sigma_test->AddEntry(gr_eff_PU200_nonprompt_2sigma_EE, "2sigma PU200");
  leg_eff_nonprompt_EE_sigma_test->AddEntry(gr_eff_PU200_nonprompt_EE, "no MTD PU200");
  leg_eff_nonprompt_EE_sigma_test->AddEntry(gr_eff_noPU_nonprompt_gen_2sigma_EE, "gen 2sigma noPU");
  leg_eff_nonprompt_EE_sigma_test->AddEntry(gr_eff_noPU_nonprompt_EE, "no MTD noPU");
    // Martina
  TLegend* leg_eff_prompt_EB_sigma_region = new TLegend(0.6, 0.13, 0.82, 0.30);
  leg_eff_prompt_EB_sigma_region->AddEntry(gr_eff_noPU_prompt_gen_EB, "gen noPU");
  leg_eff_prompt_EB_sigma_region->AddEntry(gr_eff_noPU_prompt_EB, "no MTD noPU");
  leg_eff_prompt_EB_sigma_region->SetTextSize(0.03);
  TLegend* leg_eff_prompt_EE_sigma_region = new TLegend(0.6, 0.13, 0.82, 0.30);
  leg_eff_prompt_EE_sigma_region->AddEntry(gr_eff_noPU_prompt_gen_EE, "gen noPU");
  leg_eff_prompt_EE_sigma_region->AddEntry(gr_eff_noPU_prompt_EE, "no MTD noPU");
  leg_eff_prompt_EE_sigma_region->SetTextSize(0.03);
  TLegend* leg_eff_nonprompt_EB_sigma_region = new TLegend(0.6, 0.13, 0.82, 0.30);
  leg_eff_nonprompt_EB_sigma_region->AddEntry(gr_eff_noPU_nonprompt_gen_EB, "gen noPU");
  leg_eff_nonprompt_EB_sigma_region->AddEntry(gr_eff_noPU_nonprompt_EB, "no MTD noPU");
  leg_eff_nonprompt_EB_sigma_region->SetTextSize(0.03);
  TLegend* leg_eff_nonprompt_EE_sigma_region = new TLegend(0.6, 0.13, 0.82, 0.30);
  leg_eff_nonprompt_EE_sigma_region->AddEntry(gr_eff_noPU_nonprompt_gen_EE, "gen noPU");
  leg_eff_nonprompt_EE_sigma_region->AddEntry(gr_eff_noPU_nonprompt_EE, "no MTD noPU");
  leg_eff_nonprompt_EE_sigma_region->SetTextSize(0.03);

  // test end

      // normalization
  TLegend* leg_norm_prompt_EE_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_prompt_EE_dt->AddEntry(gr_norm_PU200_prompt_gen_EE, "gen PU200");
  leg_norm_prompt_EE_dt->AddEntry(gr_norm_PU200_prompt_EE, "no MTD PU200");
  leg_norm_prompt_EE_dt->AddEntry(gr_norm_PU200_prompt_40_EE, "40ps PU200");
//  leg_norm_prompt_EE_dt->AddEntry(gr_norm_PU200_prompt_60_EE, "60ps PU200");
//  leg_norm_prompt_EE_dt->AddEntry(gr_norm_PU200_prompt_80_EE, "80ps PU200");
//  leg_norm_prompt_EE_dt->AddEntry(gr_norm_PU200_prompt_100_EE, "100ps PU200");
  leg_norm_prompt_EE_dt->AddEntry(gr_norm_noPU_prompt_gen_EE, "gen noPU");
  leg_norm_prompt_EE_dt->AddEntry(gr_norm_noPU_prompt_EE, "no MTD noPU");
  leg_norm_prompt_EE_dt->SetTextSize(0.03);
  TLegend* leg_norm_prompt_EE_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_PU200_prompt_gen_EE, "gen PU200");
  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_noPU_prompt_gen_EE, "gen noPU");
  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_PU200_prompt_2sigma_EE, "2sigma PU200");
//  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_PU200_prompt_3sigma_EE, "3sigma PU200");
//  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_PU200_prompt_4sigma_EE, "4sigma PU200");
  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_PU200_prompt_EE, "no MTD PU200");
  leg_norm_prompt_EE_sigma->AddEntry(gr_norm_noPU_prompt_EE, "no MTD noPU");
  leg_norm_prompt_EE_sigma->SetTextSize(0.03);

  // Nonprompt
    // Barrel region
      // efficiency
  TLegend* leg_eff_nonprompt_EB_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_gen_EB, "gen PU200");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_noPU_nonprompt_gen_EB, "gen noPU");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_40_EB, "40ps PU200");
//  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_60_EB, "60ps PU200");
//  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_80_EB, "80ps PU200");
//  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_100_EB, "100ps PU200");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_EB, "no MTD PU200");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_noPU_nonprompt_EB, "no MTD noPU");
  leg_eff_nonprompt_EB_dt->SetTextSize(0.03);
  TLegend* leg_eff_nonprompt_EB_sigma = new TLegend(0.62, 0.63, 0.88, 0.83);
  leg_eff_nonprompt_EB_sigma->SetMargin(0.2);
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_gen_EB, "gen PU200");
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_noPU_nonprompt_gen_EB, "gen noPU");
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_2sigma_EB_vtx, "2sigma PU200 (PV, track)");
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_2sigma_EB, "2sigma PU200 (muon, track)");
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_noPU_nonprompt_2sigma_EB, "2sigma noPU");
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_3sigma_EB, "3sigma PU200");
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_4sigma_EB, "4sigma PU200");
  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_EB, "no MTD PU200");
  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_noPU_nonprompt_EB, "no MTD noPU");
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
  leg_norm_nonprompt_EB_dt->AddEntry(gr_norm_noPU_nonprompt_EB, "no MTD noPU");
  leg_norm_nonprompt_EB_dt->SetTextSize(0.03);
  TLegend* leg_norm_nonprompt_EB_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_PU200_nonprompt_gen_EB, "gen PU200");
  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_PU200_nonprompt_EB, "no MTD PU200");
  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_PU200_nonprompt_2sigma_EB, "2sigma PU200");
//  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_PU200_nonprompt_3sigma_EB, "3sigma PU200");
//  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_PU200_nonprompt_4sigma_EB, "4sigma PU200");
  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_noPU_nonprompt_gen_EB, "gen noPU");
  leg_norm_nonprompt_EB_sigma->AddEntry(gr_norm_noPU_nonprompt_EB, "no MTD noPU");
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
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_noPU_nonprompt_EE, "no MTD noPU");
  leg_eff_nonprompt_EE_dt->SetTextSize(0.03);
  TLegend* leg_eff_nonprompt_EE_sigma = new TLegend(0.62, 0.63, 0.88, 0.83);
  leg_eff_nonprompt_EE_sigma->SetMargin(0.2);
//  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_gen_EE, "gen PU200");
//  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_noPU_nonprompt_gen_EE, "gen noPU");
//  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_2sigma_EE_vtx, "2sigma PU200 (PV, track)");
//  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_2sigma_EE, "2sigma PU200 (muon, track)");
//  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_3sigma_EE, "3sigma PU200");
//  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_4sigma_EE, "4sigma PU200");
  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_EE, "no MTD PU200");
  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_noPU_nonprompt_EE, "no MTD noPU");
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
  leg_norm_nonprompt_EE_dt->AddEntry(gr_norm_noPU_nonprompt_EE, "no MTD noPU");
  leg_norm_nonprompt_EE_dt->SetTextSize(0.03);
  TLegend* leg_norm_nonprompt_EE_sigma = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_PU200_nonprompt_gen_EE, "gen PU200");
  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_PU200_nonprompt_EE, "no MTD PU200");
  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_PU200_nonprompt_2sigma_EE, "2sigma PU200");
//  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_PU200_nonprompt_3sigma_EE, "3sigma PU200");
//  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_PU200_nonprompt_4sigma_EE, "4sigma PU200");
  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_noPU_nonprompt_gen_EE, "gen noPU");
  leg_norm_nonprompt_EE_sigma->AddEntry(gr_norm_noPU_nonprompt_EE, "no MTD noPU");
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
  gr_eff_PU200_prompt_EB->SetTitle("");
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
//  gr_eff_PU200_prompt_EB->GetXaxis()->SetRangeUser(0., 0.08);
  gr_eff_PU200_prompt_EB->SetTitle("");
  gr_eff_PU200_prompt_EB->Draw("ALP");
//  gr_eff_PU200_prompt_2sigma_EB->Draw("same");
//  gr_eff_PU200_prompt_2sigma_EB_vtx->Draw("same");
  gr_eff_noPU_prompt_EB->Draw("same LP");
//  gr_eff_PU200_prompt_gen_EB->Draw("same");
//  gr_eff_noPU_prompt_gen_EB->Draw("same");
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
  gr_eff_PU200_prompt_EE->GetXaxis()->SetRangeUser(0., 0.28);
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
//  gr_eff_PU200_prompt_EE->GetXaxis()->SetRangeUser(0., 0.08);
  gr_eff_PU200_prompt_EE->SetTitle("");
  gr_eff_PU200_prompt_EE->Draw("ALP");
//  gr_eff_PU200_prompt_2sigma_EE->Draw("same");
//  gr_eff_PU200_prompt_2sigma_EE_vtx->Draw("same");
//  gr_eff_PU200_prompt_3sigma_EE->Draw("same");
//  gr_eff_PU200_prompt_4sigma_EE->Draw("same");
  gr_eff_noPU_prompt_EE->Draw("same LP");
//  gr_eff_PU200_prompt_gen_EE->Draw("same");
//  gr_eff_noPU_prompt_gen_EE->Draw("same");
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
  gr_norm_noPU_prompt_EE->SetTitle("");
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
  gr_eff_PU200_nonprompt_EB->GetXaxis()->SetRangeUser(0., 0.08);
  gr_eff_PU200_nonprompt_EB->SetTitle("");
  gr_eff_PU200_nonprompt_EB->Draw("ALP");
//  gr_eff_PU200_nonprompt_2sigma_EB->Draw("same");
//  gr_eff_PU200_nonprompt_2sigma_EB_vtx->Draw("same");
//  gr_eff_noPU_nonprompt_2sigma_EB->Draw("same");
//  gr_eff_PU200_nonprompt_3sigma_EB->Draw("same");
//  gr_eff_PU200_nonprompt_4sigma_EB->Draw("same");
  gr_eff_noPU_nonprompt_EB->Draw("same LP");
//  gr_eff_PU200_nonprompt_gen_EB->Draw("same");
//  gr_eff_noPU_nonprompt_gen_EB->Draw("same");
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
  gr_norm_noPU_nonprompt_EB->SetTitle("");
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
  gr_eff_PU200_nonprompt_EE->SetTitle("");
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
  gr_eff_PU200_nonprompt_EE->GetXaxis()->SetRangeUser(0., 0.08);
  gr_eff_PU200_nonprompt_EE->SetTitle("");
  gr_eff_PU200_nonprompt_EE->Draw("ALP");
//  gr_eff_PU200_nonprompt_2sigma_EE->Draw("same");
//  gr_eff_PU200_nonprompt_2sigma_EE_vtx->Draw("same");
//  gr_eff_PU200_nonprompt_3sigma_EE->Draw("same");
//  gr_eff_PU200_nonprompt_4sigma_EE->Draw("same");
  gr_eff_noPU_nonprompt_EE->Draw("same LP");
//  gr_eff_PU200_nonprompt_gen_EE->Draw("same");
//  gr_eff_noPU_nonprompt_gen_EE->Draw("same");
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

  // test
  TCanvas* c_eff_prompt_sigma_EB_test = new TCanvas("c_eff_prompt_sigma_EB_test", "c_eff_prompt_sigma_EB_test", 1500, 1500);
  c_eff_prompt_sigma_EB_test->cd();
  c_eff_prompt_sigma_EB_test->SetGrid();
  c_eff_prompt_sigma_EB_test->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EB->Draw("AL");
  gr_eff_PU200_prompt_2sigma_EB->Draw("same");
  gr_eff_noPU_prompt_EB->Draw("same");
  gr_eff_PU200_prompt_gen_2sigma_EB->Draw("same");
  gr_eff_noPU_prompt_gen_2sigma_EB->Draw("same");
  leg_eff_prompt_EB_sigma_test->Draw();
  c_eff_prompt_sigma_EB_test->Print("plots/eff_prompt_sigma_EB_gen_test.pdf");

  TCanvas* c_eff_prompt_sigma_EE_test = new TCanvas("c_eff_prompt_sigma_EE_test", "c_eff_prompt_sigma_EE_test", 1500, 1500);
  c_eff_prompt_sigma_EE_test->cd();
  c_eff_prompt_sigma_EE_test->SetGrid();
  c_eff_prompt_sigma_EE_test->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EE->Draw("AL");
  gr_eff_PU200_prompt_2sigma_EE->Draw("same");
  gr_eff_noPU_prompt_EE->Draw("same");
  gr_eff_PU200_prompt_gen_2sigma_EE->Draw("same");
  gr_eff_noPU_prompt_gen_2sigma_EE->Draw("same");
  leg_eff_prompt_EE_sigma_test->Draw();
  c_eff_prompt_sigma_EE_test->Print("plots/eff_prompt_sigma_EE_gen_test.pdf");

  TCanvas* c_eff_nonprompt_sigma_EB_test = new TCanvas("c_eff_nonprompt_sigma_EB_test", "c_eff_nonprompt_sigma_EB_test", 1500, 1500);
  c_eff_nonprompt_sigma_EB_test->cd();
  c_eff_nonprompt_sigma_EB_test->SetGrid();
  c_eff_nonprompt_sigma_EB_test->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EB->GetYaxis()->SetRangeUser(0.,0.23);
  gr_eff_PU200_nonprompt_EB->Draw("AL");
  gr_eff_PU200_nonprompt_2sigma_EB->Draw("same");
  gr_eff_noPU_nonprompt_EB->Draw("same");
  gr_eff_PU200_nonprompt_gen_2sigma_EB->Draw("same");
  gr_eff_noPU_nonprompt_gen_2sigma_EB->Draw("same");
  leg_eff_nonprompt_EB_sigma_test->Draw();
  c_eff_nonprompt_sigma_EB_test->Print("plots/eff_nonprompt_sigma_EB_gen_test.pdf");

  TCanvas* c_eff_nonprompt_sigma_EE_test = new TCanvas("c_eff_nonprompt_sigma_EE_test", "c_eff_nonprompt_sigma_EE_test", 1500, 1500);
  c_eff_nonprompt_sigma_EE_test->cd();
  c_eff_nonprompt_sigma_EE_test->SetGrid();
  c_eff_nonprompt_sigma_EE_test->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EE->GetYaxis()->SetRangeUser(0.,0.33);
  gr_eff_PU200_nonprompt_EE->Draw("AL");
  gr_eff_PU200_nonprompt_2sigma_EE->Draw("same");
  gr_eff_noPU_nonprompt_EE->Draw("same");
  gr_eff_PU200_nonprompt_gen_2sigma_EE->Draw("same");
  gr_eff_noPU_nonprompt_gen_2sigma_EE->Draw("same");
  leg_eff_nonprompt_EE_sigma_test->Draw();
  c_eff_nonprompt_sigma_EE_test->Print("plots/eff_nonprompt_sigma_EE_gen_test.pdf");

  // Martina
  TCanvas* c_eff_prompt_sigma_EB_region = new TCanvas("c_eff_prompt_sigma_EB_region", "c_eff_prompt_sigma_EB_region", 1500, 1500);
  c_eff_prompt_sigma_EB_region->cd();
  c_eff_prompt_sigma_EB_region->SetGrid();
  c_eff_prompt_sigma_EB_region->SetLeftMargin(0.12);

  //gr_eff_noPU_prompt_EB->GetXaxis()->SetRangeUser(0., 0.07);
  
  //gr_eff_noPU_prompt_EB->GetXaxis()->SetRangeUser(0.07, 0.14);
  //gr_eff_noPU_prompt_EB->GetYaxis()->SetRangeUser(0.94, 1.);
  
  //gr_eff_noPU_prompt_EB->GetXaxis()->SetRangeUser(0.14, 0.21);
  //gr_eff_noPU_prompt_EB->GetYaxis()->SetRangeUser(0.98, 1.);

  gr_eff_noPU_prompt_EB->GetXaxis()->SetRangeUser(0.21, 0.28);
  gr_eff_noPU_prompt_EB->GetYaxis()->SetRangeUser(0.994, 1.);

  gr_eff_noPU_prompt_EB->Draw("AL");
  gr_eff_noPU_prompt_gen_EB->Draw("same");
  leg_eff_prompt_EB_sigma_region->Draw();
  c_eff_prompt_sigma_EB_region->Print("plots/eff_prompt_sigma_EB_region.pdf");

  TCanvas* c_eff_prompt_sigma_EE_region = new TCanvas("c_eff_prompt_sigma_EE_region", "c_eff_prompt_sigma_EE_region", 1500, 1500);
  c_eff_prompt_sigma_EE_region->cd();
  c_eff_prompt_sigma_EE_region->SetGrid();
  c_eff_prompt_sigma_EE_region->SetLeftMargin(0.12);

  //gr_eff_noPU_prompt_EE->GetXaxis()->SetRangeUser(0., 0.07);
  
  //gr_eff_noPU_prompt_EE->GetXaxis()->SetRangeUser(0.07, 0.14);
  //gr_eff_noPU_prompt_EE->GetYaxis()->SetRangeUser(0.94, 1.);
  
  //gr_eff_noPU_prompt_EE->GetXaxis()->SetRangeUser(0.14, 0.21);
  //gr_eff_noPU_prompt_EE->GetYaxis()->SetRangeUser(0.98, 1.);

  gr_eff_noPU_prompt_EE->GetXaxis()->SetRangeUser(0.21, 0.28);
  gr_eff_noPU_prompt_EE->GetYaxis()->SetRangeUser(0.993, 1.);

  gr_eff_noPU_prompt_EE->Draw("AL");
  gr_eff_noPU_prompt_gen_EE->Draw("same");
  leg_eff_prompt_EE_sigma_region->Draw();
  c_eff_prompt_sigma_EE_region->Print("plots/eff_prompt_sigma_EE_region.pdf");

  TCanvas* c_eff_nonprompt_sigma_EB_region = new TCanvas("c_eff_nonprompt_sigma_EB_region", "c_eff_nonprompt_sigma_EB_region", 1500, 1500);
  c_eff_nonprompt_sigma_EB_region->cd();
  c_eff_nonprompt_sigma_EB_region->SetGrid();
  c_eff_nonprompt_sigma_EB_region->SetLeftMargin(0.12);

  //gr_eff_noPU_nonprompt_EB->GetXaxis()->SetRangeUser(0., 1.);
  
  //gr_eff_noPU_nonprompt_EB->GetXaxis()->SetRangeUser(1., 2.);
  //gr_eff_noPU_nonprompt_EB->GetYaxis()->SetRangeUser(0.4, 0.8);
  
  //gr_eff_noPU_nonprompt_EB->GetXaxis()->SetRangeUser(2., 3.);
  //gr_eff_noPU_nonprompt_EB->GetYaxis()->SetRangeUser(0.72, 0.9);

  gr_eff_noPU_nonprompt_EB->GetXaxis()->SetRangeUser(3., 4.);
  gr_eff_noPU_nonprompt_EB->GetYaxis()->SetRangeUser(0.85, 1.);

  gr_eff_noPU_nonprompt_EB->Draw("AL");
  gr_eff_noPU_nonprompt_gen_EB->Draw("same");
  leg_eff_nonprompt_EB_sigma_region->Draw();
  c_eff_nonprompt_sigma_EB_region->Print("plots/eff_nonprompt_sigma_EB_region.pdf");

  TCanvas* c_eff_nonprompt_sigma_EE_region = new TCanvas("c_eff_nonprompt_sigma_EE_region", "c_eff_nonprompt_sigma_EE_region", 1500, 1500);
  c_eff_nonprompt_sigma_EE_region->cd();
  c_eff_nonprompt_sigma_EE_region->SetGrid();
  c_eff_nonprompt_sigma_EE_region->SetLeftMargin(0.12);

  //gr_eff_noPU_nonprompt_EE->GetXaxis()->SetRangeUser(0., 1.);
  
  //gr_eff_noPU_nonprompt_EE->GetXaxis()->SetRangeUser(1., 2.);
  //gr_eff_noPU_nonprompt_EE->GetYaxis()->SetRangeUser(0.58, 0.85);

  //gr_eff_noPU_nonprompt_EE->GetXaxis()->SetRangeUser(2., 3.);
  //gr_eff_noPU_nonprompt_EE->GetYaxis()->SetRangeUser(0.8, 0.94);

  gr_eff_noPU_nonprompt_EE->GetXaxis()->SetRangeUser(3., 4.);
  gr_eff_noPU_nonprompt_EE->GetYaxis()->SetRangeUser(0.9, 1.);

  gr_eff_noPU_nonprompt_EE->Draw("AL");
  gr_eff_noPU_nonprompt_gen_EE->Draw("same");
  leg_eff_nonprompt_EE_sigma_region->Draw();
  c_eff_nonprompt_sigma_EE_region->Print("plots/eff_nonprompt_sigma_EE_region.pdf");



  TCanvas* c_eff_prompt_sigma_EB_test2 = new TCanvas("c_eff_prompt_sigma_EB_test2", "c_eff_prompt_sigma_EB_test2", 1500, 1500);
  c_eff_prompt_sigma_EB_test2->cd();
  c_eff_prompt_sigma_EB_test2->SetGrid();
  c_eff_prompt_sigma_EB_test2->SetLeftMargin(0.12);
  gr_eff_noPU_prompt_EB->SetTitle("");
  gr_eff_noPU_prompt_EB->GetXaxis()->SetRangeUser(0., 0.15);
  gr_eff_noPU_prompt_EB->GetYaxis()->SetRangeUser(0.8, 1.);
  gr_eff_noPU_prompt_EB->Draw("AL");
  gr_eff_noPU_prompt_2sigma_EB->Draw("same");
  gr_eff_noPU_prompt_3sigma_EB->Draw("same");
  gr_eff_noPU_prompt_4sigma_EB->Draw("same");
  gr_eff_noPU_prompt_gen_EB->Draw("same");
  c_eff_prompt_sigma_EB_test2->Print("test_prompt_EB.pdf");

  TCanvas* c_eff_nonprompt_sigma_EB_test2 = new TCanvas("c_eff_nonprompt_sigma_EB_test2", "c_eff_nonprompt_sigma_EB_test2", 1500, 1500);
  c_eff_nonprompt_sigma_EB_test2->cd();
  c_eff_nonprompt_sigma_EB_test2->SetGrid();
  c_eff_nonprompt_sigma_EB_test2->SetLeftMargin(0.12);
  gr_eff_noPU_nonprompt_EB->SetTitle("");
  gr_eff_noPU_nonprompt_EB->GetXaxis()->SetRangeUser(0., 0.8);
  gr_eff_noPU_nonprompt_EB->GetYaxis()->SetRangeUser(0., 1.);
  gr_eff_noPU_nonprompt_EB->Draw("AL");
  gr_eff_noPU_nonprompt_2sigma_EB->Draw("same");
  gr_eff_noPU_nonprompt_3sigma_EB->Draw("same");
  gr_eff_noPU_nonprompt_4sigma_EB->Draw("same");
  gr_eff_noPU_nonprompt_gen_EB->Draw("same");
  c_eff_nonprompt_sigma_EB_test2->Print("test_nonprompt_EB.pdf");



  // test end

}

void draw_iso_distribution() {
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  TFile *f_PU200_prompt_vtx, *f_PU200_nonprompt_vtx, *f_noPU_prompt_vtx, *f_noPU_nonprompt_vtx;
  // PU200
  f_PU200_prompt = new TFile(Form("data/%s", path_PU200_prompt.Data()));
  f_PU200_nonprompt = new TFile(Form("data/%s", path_PU200_nonprompt.Data()));
  f_PU200_prompt_vtx = new TFile(Form("data/%s", path_PU200_prompt_vtx.Data()));
  f_PU200_nonprompt_vtx = new TFile(Form("data/%s", path_PU200_nonprompt_vtx.Data()));
  // noPU
  f_noPU_prompt = new TFile(Form("data/%s", path_noPU_prompt.Data()));
  f_noPU_nonprompt = new TFile(Form("data/%s", path_noPU_nonprompt.Data()));
  f_noPU_prompt_vtx = new TFile(Form("data/%s", path_noPU_prompt_vtx.Data()));
  f_noPU_nonprompt_vtx = new TFile(Form("data/%s", path_noPU_nonprompt_vtx.Data()));

  // Histograms for Barrel region
    // PU200
  TH1D *h_PU200_prompt_EB, *h_PU200_nonprompt_EB, *h_noPU_prompt_EB, *h_noPU_nonprompt_EB;
  TH1D *h_PU200_prompt_2sigma_EB, *h_PU200_nonprompt_2sigma_EB, *h_noPU_prompt_2sigma_EB, *h_noPU_nonprompt_2sigma_EB;
  TH1D *h_PU200_prompt_3sigma_EB, *h_PU200_nonprompt_3sigma_EB, *h_noPU_prompt_3sigma_EB, *h_noPU_nonprompt_3sigma_EB;
  TH1D *h_PU200_prompt_4sigma_EB, *h_PU200_nonprompt_4sigma_EB, *h_noPU_prompt_4sigma_EB, *h_noPU_nonprompt_4sigma_EB;
  TH1D *h_PU200_prompt_40_EB, *h_PU200_nonprompt_40_EB, *h_noPU_prompt_40_EB, *h_noPU_nonprompt_40_EB;
  TH1D *h_PU200_prompt_60_EB, *h_PU200_nonprompt_60_EB, *h_noPU_prompt_60_EB, *h_noPU_nonprompt_60_EB;
  TH1D *h_PU200_prompt_80_EB, *h_PU200_nonprompt_80_EB, *h_noPU_prompt_80_EB, *h_noPU_nonprompt_80_EB;
  TH1D *h_PU200_prompt_100_EB, *h_PU200_nonprompt_100_EB, *h_noPU_prompt_100_EB, *h_noPU_nonprompt_100_EB;
      // GEN case
  TH1D *h_PU200_prompt_gen_EB, *h_PU200_nonprompt_gen_EB, *h_noPU_prompt_gen_EB, *h_noPU_nonprompt_gen_EB;
  TH1D *h_PU200_prompt_gen_2sigma_EB, *h_PU200_nonprompt_gen_2sigma_EB, *h_noPU_prompt_gen_2sigma_EB, *h_noPU_nonprompt_gen_2sigma_EB;
  TH1D *h_PU200_prompt_gen_3sigma_EB, *h_PU200_nonprompt_gen_3sigma_EB, *h_noPU_prompt_gen_3sigma_EB, *h_noPU_nonprompt_gen_3sigma_EB;
  TH1D *h_PU200_prompt_gen_4sigma_EB, *h_PU200_nonprompt_gen_4sigma_EB, *h_noPU_prompt_gen_4sigma_EB, *h_noPU_nonprompt_gen_4sigma_EB;
      // vtx
  TH1D *h_PU200_prompt_2sigma_EB_vtx, *h_PU200_nonprompt_2sigma_EB_vtx, *h_noPU_prompt_2sigma_EB_vtx, *h_noPU_nonprompt_2sigma_EB_vtx;
  TH1D *h_PU200_prompt_3sigma_EB_vtx, *h_PU200_nonprompt_3sigma_EB_vtx, *h_noPU_prompt_3sigma_EB_vtx, *h_noPU_nonprompt_3sigma_EB_vtx;
  TH1D *h_PU200_prompt_4sigma_EB_vtx, *h_PU200_nonprompt_4sigma_EB_vtx, *h_noPU_prompt_4sigma_EB_vtx, *h_noPU_nonprompt_4sigma_EB_vtx;

  // Histograms for Endcap region
    // PU200
  TH1D *h_PU200_prompt_EE, *h_PU200_nonprompt_EE, *h_noPU_prompt_EE, *h_noPU_nonprompt_EE;
  TH1D *h_PU200_prompt_2sigma_EE, *h_PU200_nonprompt_2sigma_EE, *h_noPU_prompt_2sigma_EE, *h_noPU_nonprompt_2sigma_EE;
  TH1D *h_PU200_prompt_3sigma_EE, *h_PU200_nonprompt_3sigma_EE, *h_noPU_prompt_3sigma_EE, *h_noPU_nonprompt_3sigma_EE;
  TH1D *h_PU200_prompt_4sigma_EE, *h_PU200_nonprompt_4sigma_EE, *h_noPU_prompt_4sigma_EE, *h_noPU_nonprompt_4sigma_EE;
  TH1D *h_PU200_prompt_40_EE, *h_PU200_nonprompt_40_EE, *h_noPU_prompt_40_EE, *h_noPU_nonprompt_40_EE;
  TH1D *h_PU200_prompt_60_EE, *h_PU200_nonprompt_60_EE, *h_noPU_prompt_60_EE, *h_noPU_nonprompt_60_EE;
  TH1D *h_PU200_prompt_80_EE, *h_PU200_nonprompt_80_EE, *h_noPU_prompt_80_EE, *h_noPU_nonprompt_80_EE;
  TH1D *h_PU200_prompt_100_EE, *h_PU200_nonprompt_100_EE, *h_noPU_prompt_100_EE, *h_noPU_nonprompt_100_EE;
      // GEN case
  TH1D *h_PU200_prompt_gen_EE, *h_PU200_nonprompt_gen_EE, *h_noPU_prompt_gen_EE, *h_noPU_nonprompt_gen_EE;
  TH1D *h_PU200_prompt_gen_2sigma_EE, *h_PU200_nonprompt_gen_2sigma_EE, *h_noPU_prompt_gen_2sigma_EE, *h_noPU_nonprompt_gen_2sigma_EE;
  TH1D *h_PU200_prompt_gen_3sigma_EE, *h_PU200_nonprompt_gen_3sigma_EE, *h_noPU_prompt_gen_3sigma_EE, *h_noPU_nonprompt_gen_3sigma_EE;
  TH1D *h_PU200_prompt_gen_4sigma_EE, *h_PU200_nonprompt_gen_4sigma_EE, *h_noPU_prompt_gen_4sigma_EE, *h_noPU_nonprompt_gen_4sigma_EE;
      // vtx
  TH1D *h_PU200_prompt_2sigma_EE_vtx, *h_PU200_nonprompt_2sigma_EE_vtx, *h_noPU_prompt_2sigma_EE_vtx, *h_noPU_nonprompt_2sigma_EE_vtx;
  TH1D *h_PU200_prompt_3sigma_EE_vtx, *h_PU200_nonprompt_3sigma_EE_vtx, *h_noPU_prompt_3sigma_EE_vtx, *h_noPU_nonprompt_3sigma_EE_vtx;
  TH1D *h_PU200_prompt_4sigma_EE_vtx, *h_PU200_nonprompt_4sigma_EE_vtx, *h_noPU_prompt_4sigma_EE_vtx, *h_noPU_nonprompt_4sigma_EE_vtx;


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
  h_PU200_prompt_gen_4sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Sig_EB");
  h_PU200_prompt_gen_3sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Sig_EB");
  h_PU200_prompt_gen_2sigma_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Sig_EB");
    // vtx
  h_PU200_prompt_4sigma_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EB");
  h_PU200_prompt_3sigma_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EB");
  h_PU200_prompt_2sigma_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EB");

  h_noPU_prompt_EB 	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_Sig_EB");
  h_noPU_prompt_4sigma_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EB");
  h_noPU_prompt_3sigma_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EB");
  h_noPU_prompt_2sigma_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EB");
  h_noPU_prompt_40_EB     = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Sig_EB");
  h_noPU_prompt_60_EB     = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Sig_EB");
  h_noPU_prompt_80_EB     = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Sig_EB");
  h_noPU_prompt_100_EB    = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Sig_EB");
  h_noPU_prompt_gen_EB	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_gen_Sig_EB");
  h_noPU_prompt_gen_4sigma_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Sig_EB");
  h_noPU_prompt_gen_3sigma_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Sig_EB");
  h_noPU_prompt_gen_2sigma_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Sig_EB");
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
  h_PU200_prompt_gen_4sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Sig_EE");
  h_PU200_prompt_gen_3sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Sig_EE");
  h_PU200_prompt_gen_2sigma_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Sig_EE");
    // vtx
  h_PU200_prompt_4sigma_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EE");
  h_PU200_prompt_3sigma_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EE");
  h_PU200_prompt_2sigma_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EE");

  h_noPU_prompt_EE 	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_Sig_EE");
  h_noPU_prompt_4sigma_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EE");
  h_noPU_prompt_3sigma_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EE");
  h_noPU_prompt_2sigma_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EE");
  h_noPU_prompt_40_EE     = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Sig_EE");
  h_noPU_prompt_60_EE     = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Sig_EE");
  h_noPU_prompt_80_EE     = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Sig_EE");
  h_noPU_prompt_100_EE    = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Sig_EE");
  h_noPU_prompt_gen_EE	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_gen_Sig_EE");
  h_noPU_prompt_gen_4sigma_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Sig_EE");
  h_noPU_prompt_gen_3sigma_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Sig_EE");
  h_noPU_prompt_gen_2sigma_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Sig_EE");
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
  h_PU200_nonprompt_gen_4sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Bkg_EB");
  h_PU200_nonprompt_gen_3sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Bkg_EB");
  h_PU200_nonprompt_gen_2sigma_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Bkg_EB");
      // vtx
  h_PU200_nonprompt_4sigma_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EB");
  h_PU200_nonprompt_3sigma_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EB");
  h_PU200_nonprompt_2sigma_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EB");

  h_noPU_nonprompt_EB 	      = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_Bkg_EB");
  h_noPU_nonprompt_4sigma_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EB");
  h_noPU_nonprompt_3sigma_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EB");
  h_noPU_nonprompt_2sigma_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EB");
  h_noPU_nonprompt_40_EB     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Bkg_EB");
  h_noPU_nonprompt_60_EB     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Bkg_EB");
  h_noPU_nonprompt_80_EB     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Bkg_EB");
  h_noPU_nonprompt_100_EB    = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Bkg_EB");
  h_noPU_nonprompt_gen_EB    = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_gen_Bkg_EB");
  h_noPU_nonprompt_gen_4sigma_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Bkg_EB");
  h_noPU_nonprompt_gen_3sigma_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Bkg_EB");
  h_noPU_nonprompt_gen_2sigma_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Bkg_EB");
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
  h_PU200_nonprompt_gen_4sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Bkg_EE");
  h_PU200_nonprompt_gen_3sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Bkg_EE");
  h_PU200_nonprompt_gen_2sigma_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Bkg_EE");
      // vtx
  h_PU200_nonprompt_4sigma_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EE");
  h_PU200_nonprompt_3sigma_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EE");
  h_PU200_nonprompt_2sigma_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EE");

  h_noPU_nonprompt_EE 	      = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_Bkg_EE");
  h_noPU_nonprompt_4sigma_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EE");
  h_noPU_nonprompt_3sigma_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EE");
  h_noPU_nonprompt_2sigma_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EE");
  h_noPU_nonprompt_40_EE     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_7_Bkg_EE");
  h_noPU_nonprompt_60_EE     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_5_Bkg_EE");
  h_noPU_nonprompt_80_EE     = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_3_Bkg_EE");
  h_noPU_nonprompt_100_EE    = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_1_Bkg_EE");
  h_noPU_nonprompt_gen_EE    = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_gen_Bkg_EE");
  h_noPU_nonprompt_gen_4sigma_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_4sigma_Bkg_EE");
  h_noPU_nonprompt_gen_3sigma_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_3sigma_Bkg_EE");
  h_noPU_nonprompt_gen_2sigma_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_rel_chIso_sum_MTD_gen_2sigma_Bkg_EE");

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
        // vtx
  vector<double> prompt_eff_PU200_2sigma_EB_vtx={0}, prompt_norm_PU200_2sigma_EB_vtx={0}, prompt_eff_PU200_2sigma_EE_vtx={0}, prompt_norm_PU200_2sigma_EE_vtx={0};
  vector<double> prompt_eff_PU200_3sigma_EB_vtx={0}, prompt_norm_PU200_3sigma_EB_vtx={0}, prompt_eff_PU200_3sigma_EE_vtx={0}, prompt_norm_PU200_3sigma_EE_vtx={0};
  vector<double> prompt_eff_PU200_4sigma_EB_vtx={0}, prompt_norm_PU200_4sigma_EB_vtx={0}, prompt_eff_PU200_4sigma_EE_vtx={0}, prompt_norm_PU200_4sigma_EE_vtx={0};
      // nonprompt
  vector<double> nonprompt_eff_PU200_EB={0}, nonprompt_norm_PU200_EB={0}, nonprompt_eff_PU200_EE={0}, nonprompt_norm_PU200_EE={0};
  vector<double> nonprompt_eff_PU200_2sigma_EB={0}, nonprompt_norm_PU200_2sigma_EB={0}, nonprompt_eff_PU200_2sigma_EE={0}, nonprompt_norm_PU200_2sigma_EE={0};
  vector<double> nonprompt_eff_PU200_3sigma_EB={0}, nonprompt_norm_PU200_3sigma_EB={0}, nonprompt_eff_PU200_3sigma_EE={0}, nonprompt_norm_PU200_3sigma_EE={0};
  vector<double> nonprompt_eff_PU200_4sigma_EB={0}, nonprompt_norm_PU200_4sigma_EB={0}, nonprompt_eff_PU200_4sigma_EE={0}, nonprompt_norm_PU200_4sigma_EE={0};
  vector<double> nonprompt_eff_PU200_40_EB={0}, nonprompt_norm_PU200_40_EB={0}, nonprompt_eff_PU200_40_EE={0}, nonprompt_norm_PU200_40_EE={0};
  vector<double> nonprompt_eff_PU200_60_EB={0}, nonprompt_norm_PU200_60_EB={0}, nonprompt_eff_PU200_60_EE={0}, nonprompt_norm_PU200_60_EE={0};
  vector<double> nonprompt_eff_PU200_80_EB={0}, nonprompt_norm_PU200_80_EB={0}, nonprompt_eff_PU200_80_EE={0}, nonprompt_norm_PU200_80_EE={0};
  vector<double> nonprompt_eff_PU200_100_EB={0}, nonprompt_norm_PU200_100_EB={0}, nonprompt_eff_PU200_100_EE={0}, nonprompt_norm_PU200_100_EE={0};
        // vtx
  vector<double> nonprompt_eff_PU200_2sigma_EB_vtx={0}, nonprompt_norm_PU200_2sigma_EB_vtx={0}, nonprompt_eff_PU200_2sigma_EE_vtx={0}, nonprompt_norm_PU200_2sigma_EE_vtx={0};
  vector<double> nonprompt_eff_PU200_3sigma_EB_vtx={0}, nonprompt_norm_PU200_3sigma_EB_vtx={0}, nonprompt_eff_PU200_3sigma_EE_vtx={0}, nonprompt_norm_PU200_3sigma_EE_vtx={0};
  vector<double> nonprompt_eff_PU200_4sigma_EB_vtx={0}, nonprompt_norm_PU200_4sigma_EB_vtx={0}, nonprompt_eff_PU200_4sigma_EE_vtx={0}, nonprompt_norm_PU200_4sigma_EE_vtx={0};
      // GEN case
        // prompt
  vector<double> prompt_eff_PU200_gen_EB={0}, prompt_norm_PU200_gen_EB={0}, prompt_eff_PU200_gen_EE={0}, prompt_norm_PU200_gen_EE={0};
  vector<double> prompt_eff_PU200_gen_2sigma_EB={0}, prompt_norm_PU200_gen_2sigma_EB={0}, prompt_eff_PU200_gen_2sigma_EE={0}, prompt_norm_PU200_gen_2sigma_EE={0};
  vector<double> prompt_eff_PU200_gen_3sigma_EB={0}, prompt_norm_PU200_gen_3sigma_EB={0}, prompt_eff_PU200_gen_3sigma_EE={0}, prompt_norm_PU200_gen_3sigma_EE={0};
  vector<double> prompt_eff_PU200_gen_4sigma_EB={0}, prompt_norm_PU200_gen_4sigma_EB={0}, prompt_eff_PU200_gen_4sigma_EE={0}, prompt_norm_PU200_gen_4sigma_EE={0};
        // nonprompt
  vector<double> nonprompt_eff_PU200_gen_EB={0}, nonprompt_norm_PU200_gen_EB={0}, nonprompt_eff_PU200_gen_EE={0}, nonprompt_norm_PU200_gen_EE={0};
  vector<double> nonprompt_eff_PU200_gen_2sigma_EB={0}, nonprompt_norm_PU200_gen_2sigma_EB={0}, nonprompt_eff_PU200_gen_2sigma_EE={0}, nonprompt_norm_PU200_gen_2sigma_EE={0};
  vector<double> nonprompt_eff_PU200_gen_3sigma_EB={0}, nonprompt_norm_PU200_gen_3sigma_EB={0}, nonprompt_eff_PU200_gen_3sigma_EE={0}, nonprompt_norm_PU200_gen_3sigma_EE={0};
  vector<double> nonprompt_eff_PU200_gen_4sigma_EB={0}, nonprompt_norm_PU200_gen_4sigma_EB={0}, nonprompt_eff_PU200_gen_4sigma_EE={0}, nonprompt_norm_PU200_gen_4sigma_EE={0};
    // noPU
      // prompt
  vector<double> prompt_eff_noPU_EB={0}, prompt_norm_noPU_EB={0}, prompt_eff_noPU_EE={0}, prompt_norm_noPU_EE={0};
  vector<double> prompt_eff_noPU_2sigma_EB={0}, prompt_norm_noPU_2sigma_EB={0}, prompt_eff_noPU_2sigma_EE={0}, prompt_norm_noPU_2sigma_EE={0};
  vector<double> prompt_eff_noPU_3sigma_EB={0}, prompt_norm_noPU_3sigma_EB={0}, prompt_eff_noPU_3sigma_EE={0}, prompt_norm_noPU_3sigma_EE={0};
  vector<double> prompt_eff_noPU_4sigma_EB={0}, prompt_norm_noPU_4sigma_EB={0}, prompt_eff_noPU_4sigma_EE={0}, prompt_norm_noPU_4sigma_EE={0};
  vector<double> prompt_eff_noPU_40_EB={0}, prompt_norm_noPU_40_EB={0}, prompt_eff_noPU_40_EE={0}, prompt_norm_noPU_40_EE={0};
  vector<double> prompt_eff_noPU_60_EB={0}, prompt_norm_noPU_60_EB={0}, prompt_eff_noPU_60_EE={0}, prompt_norm_noPU_60_EE={0};
  vector<double> prompt_eff_noPU_80_EB={0}, prompt_norm_noPU_80_EB={0}, prompt_eff_noPU_80_EE={0}, prompt_norm_noPU_80_EE={0};
  vector<double> prompt_eff_noPU_100_EB={0}, prompt_norm_noPU_100_EB={0}, prompt_eff_noPU_100_EE={0}, prompt_norm_noPU_100_EE={0};
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
      // prompt
  vector<double> prompt_eff_noPU_gen_EB={0}, prompt_norm_noPU_gen_EB={0}, prompt_eff_noPU_gen_EE={0}, prompt_norm_noPU_gen_EE={0};
  vector<double> prompt_eff_noPU_gen_2sigma_EB={0}, prompt_norm_noPU_gen_2sigma_EB={0}, prompt_eff_noPU_gen_2sigma_EE={0}, prompt_norm_noPU_gen_2sigma_EE={0};
  vector<double> prompt_eff_noPU_gen_3sigma_EB={0}, prompt_norm_noPU_gen_3sigma_EB={0}, prompt_eff_noPU_gen_3sigma_EE={0}, prompt_norm_noPU_gen_3sigma_EE={0};
  vector<double> prompt_eff_noPU_gen_4sigma_EB={0}, prompt_norm_noPU_gen_4sigma_EB={0}, prompt_eff_noPU_gen_4sigma_EE={0}, prompt_norm_noPU_gen_4sigma_EE={0};
      // nonprompt
  vector<double> nonprompt_eff_noPU_gen_EB={0}, nonprompt_norm_noPU_gen_EB={0}, nonprompt_eff_noPU_gen_EE={0}, nonprompt_norm_noPU_gen_EE={0};
  vector<double> nonprompt_eff_noPU_gen_2sigma_EB={0}, nonprompt_norm_noPU_gen_2sigma_EB={0}, nonprompt_eff_noPU_gen_2sigma_EE={0}, nonprompt_norm_noPU_gen_2sigma_EE={0};
  vector<double> nonprompt_eff_noPU_gen_3sigma_EB={0}, nonprompt_norm_noPU_gen_3sigma_EB={0}, nonprompt_eff_noPU_gen_3sigma_EE={0}, nonprompt_norm_noPU_gen_3sigma_EE={0};
  vector<double> nonprompt_eff_noPU_gen_4sigma_EB={0}, nonprompt_norm_noPU_gen_4sigma_EB={0}, nonprompt_eff_noPU_gen_4sigma_EE={0}, nonprompt_norm_noPU_gen_4sigma_EE={0};


  //////////////////////////
  // Calculate efficiency //
  //////////////////////////
  
  int nbin = h_PU200_prompt_EB->GetNbinsX()+1;
  for(unsigned int i=0; i<h_PU200_prompt_EB->GetNbinsX()+1; i++) {
  //for(unsigned int i=0; i<h_PU200_prompt_EB->GetNbinsX()/2+1; i++) {
    // prompt
      // Barrel region
        // efficiency
	  // PU200
	
    prompt_eff_PU200_EB.emplace_back(h_PU200_prompt_EB->GetBinContent(i)/h_PU200_prompt_EB->Integral(1,nbin));
    prompt_eff_PU200_40_EB.emplace_back(h_PU200_prompt_40_EB->GetBinContent(i)/h_PU200_prompt_40_EB->Integral(1,nbin));
    prompt_eff_PU200_60_EB.emplace_back(h_PU200_prompt_60_EB->GetBinContent(i)/h_PU200_prompt_60_EB->Integral(1,nbin));
    prompt_eff_PU200_80_EB.emplace_back(h_PU200_prompt_80_EB->GetBinContent(i)/h_PU200_prompt_80_EB->Integral(1,nbin));
    prompt_eff_PU200_100_EB.emplace_back(h_PU200_prompt_100_EB->GetBinContent(i)/h_PU200_prompt_100_EB->Integral(1,nbin));
    prompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->GetBinContent(i)/h_PU200_prompt_2sigma_EB->Integral(1,nbin));
    prompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->GetBinContent(i)/h_PU200_prompt_3sigma_EB->Integral(1,nbin));
    prompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->GetBinContent(i)/h_PU200_prompt_4sigma_EB->Integral(1,nbin));
              // vtx
    prompt_eff_PU200_2sigma_EB_vtx.emplace_back(h_PU200_prompt_2sigma_EB_vtx->GetBinContent(i)/h_PU200_prompt_2sigma_EB_vtx->Integral(1,nbin));
    prompt_eff_PU200_3sigma_EB_vtx.emplace_back(h_PU200_prompt_3sigma_EB_vtx->GetBinContent(i)/h_PU200_prompt_3sigma_EB_vtx->Integral(1,nbin));
    prompt_eff_PU200_4sigma_EB_vtx.emplace_back(h_PU200_prompt_4sigma_EB_vtx->GetBinContent(i)/h_PU200_prompt_4sigma_EB_vtx->Integral(1,nbin));
            // GEN case
    prompt_eff_PU200_gen_EB.emplace_back(h_PU200_prompt_gen_EB->GetBinContent(i)/h_PU200_prompt_gen_EB->Integral(1,nbin));
    prompt_eff_PU200_gen_2sigma_EB.emplace_back(h_PU200_prompt_gen_2sigma_EB->GetBinContent(i)/h_PU200_prompt_gen_2sigma_EB->Integral(1,nbin));
    prompt_eff_PU200_gen_3sigma_EB.emplace_back(h_PU200_prompt_gen_3sigma_EB->GetBinContent(i)/h_PU200_prompt_gen_3sigma_EB->Integral(1,nbin));
    prompt_eff_PU200_gen_4sigma_EB.emplace_back(h_PU200_prompt_gen_4sigma_EB->GetBinContent(i)/h_PU200_prompt_gen_4sigma_EB->Integral(1,nbin));
	  // noPU
    prompt_eff_noPU_EB.emplace_back(h_noPU_prompt_EB->GetBinContent(i)/h_noPU_prompt_EB->Integral(1,nbin));
    prompt_eff_noPU_40_EB.emplace_back(h_noPU_prompt_40_EB->GetBinContent(i)/h_noPU_prompt_40_EB->Integral(1,nbin));
    prompt_eff_noPU_60_EB.emplace_back(h_noPU_prompt_60_EB->GetBinContent(i)/h_noPU_prompt_60_EB->Integral(1,nbin));
    prompt_eff_noPU_80_EB.emplace_back(h_noPU_prompt_80_EB->GetBinContent(i)/h_noPU_prompt_80_EB->Integral(1,nbin));
    prompt_eff_noPU_100_EB.emplace_back(h_noPU_prompt_100_EB->GetBinContent(i)/h_noPU_prompt_100_EB->Integral(1,nbin));
    prompt_eff_noPU_2sigma_EB.emplace_back(h_noPU_prompt_2sigma_EB->GetBinContent(i)/h_noPU_prompt_2sigma_EB->Integral(1,nbin));
    prompt_eff_noPU_3sigma_EB.emplace_back(h_noPU_prompt_3sigma_EB->GetBinContent(i)/h_noPU_prompt_3sigma_EB->Integral(1,nbin));
    prompt_eff_noPU_4sigma_EB.emplace_back(h_noPU_prompt_4sigma_EB->GetBinContent(i)/h_noPU_prompt_4sigma_EB->Integral(1,nbin));
            // GEN case
    prompt_eff_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->GetBinContent(i)/h_noPU_prompt_gen_EB->Integral(1,nbin));
    prompt_eff_noPU_gen_2sigma_EB.emplace_back(h_noPU_prompt_gen_2sigma_EB->GetBinContent(i)/h_noPU_prompt_gen_2sigma_EB->Integral(1,nbin));
    prompt_eff_noPU_gen_3sigma_EB.emplace_back(h_noPU_prompt_gen_3sigma_EB->GetBinContent(i)/h_noPU_prompt_gen_3sigma_EB->Integral(1,nbin));
    prompt_eff_noPU_gen_4sigma_EB.emplace_back(h_noPU_prompt_gen_4sigma_EB->GetBinContent(i)/h_noPU_prompt_gen_4sigma_EB->Integral(1,nbin));

      // Endcap region
        // efficiency
	  // PU200
    prompt_eff_PU200_EE.emplace_back(h_PU200_prompt_EE->GetBinContent(i)/h_PU200_prompt_EE->Integral(1,nbin));
    prompt_eff_PU200_40_EE.emplace_back(h_PU200_prompt_40_EE->GetBinContent(i)/h_PU200_prompt_40_EE->Integral(1,nbin));
    prompt_eff_PU200_60_EE.emplace_back(h_PU200_prompt_60_EE->GetBinContent(i)/h_PU200_prompt_60_EE->Integral(1,nbin));
    prompt_eff_PU200_80_EE.emplace_back(h_PU200_prompt_80_EE->GetBinContent(i)/h_PU200_prompt_80_EE->Integral(1,nbin));
    prompt_eff_PU200_100_EE.emplace_back(h_PU200_prompt_100_EE->GetBinContent(i)/h_PU200_prompt_100_EE->Integral(1,nbin));
    prompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_prompt_2sigma_EE->GetBinContent(i)/h_PU200_prompt_2sigma_EE->Integral(1,nbin));
    prompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_prompt_3sigma_EE->GetBinContent(i)/h_PU200_prompt_3sigma_EE->Integral(1,nbin));
    prompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_prompt_4sigma_EE->GetBinContent(i)/h_PU200_prompt_4sigma_EE->Integral(1,nbin));
              // vtx
    prompt_eff_PU200_2sigma_EE_vtx.emplace_back(h_PU200_prompt_2sigma_EE_vtx->GetBinContent(i)/h_PU200_prompt_2sigma_EE_vtx->Integral(1,nbin));
    prompt_eff_PU200_3sigma_EE_vtx.emplace_back(h_PU200_prompt_3sigma_EE_vtx->GetBinContent(i)/h_PU200_prompt_3sigma_EE_vtx->Integral(1,nbin));
    prompt_eff_PU200_4sigma_EE_vtx.emplace_back(h_PU200_prompt_4sigma_EE_vtx->GetBinContent(i)/h_PU200_prompt_4sigma_EE_vtx->Integral(1,nbin));
            // GEN case
    prompt_eff_PU200_gen_EE.emplace_back(h_PU200_prompt_gen_EE->GetBinContent(i)/h_PU200_prompt_gen_EE->Integral(1,nbin));
    prompt_eff_PU200_gen_2sigma_EE.emplace_back(h_PU200_prompt_gen_2sigma_EE->GetBinContent(i)/h_PU200_prompt_gen_2sigma_EE->Integral(1,nbin));
    prompt_eff_PU200_gen_3sigma_EE.emplace_back(h_PU200_prompt_gen_3sigma_EE->GetBinContent(i)/h_PU200_prompt_gen_3sigma_EE->Integral(1,nbin));
    prompt_eff_PU200_gen_4sigma_EE.emplace_back(h_PU200_prompt_gen_4sigma_EE->GetBinContent(i)/h_PU200_prompt_gen_4sigma_EE->Integral(1,nbin));
	  // noPU
    prompt_eff_noPU_EE.emplace_back(h_noPU_prompt_EE->GetBinContent(i)/h_noPU_prompt_EE->Integral(1,nbin));
    prompt_eff_noPU_40_EE.emplace_back(h_noPU_prompt_40_EE->GetBinContent(i)/h_noPU_prompt_40_EE->Integral(1,nbin));
    prompt_eff_noPU_60_EE.emplace_back(h_noPU_prompt_60_EE->GetBinContent(i)/h_noPU_prompt_60_EE->Integral(1,nbin));
    prompt_eff_noPU_80_EE.emplace_back(h_noPU_prompt_80_EE->GetBinContent(i)/h_noPU_prompt_80_EE->Integral(1,nbin));
    prompt_eff_noPU_100_EE.emplace_back(h_noPU_prompt_100_EE->GetBinContent(i)/h_noPU_prompt_100_EE->Integral(1,nbin));
    prompt_eff_noPU_2sigma_EE.emplace_back(h_noPU_prompt_2sigma_EE->GetBinContent(i)/h_noPU_prompt_2sigma_EE->Integral(1,nbin));
    prompt_eff_noPU_3sigma_EE.emplace_back(h_noPU_prompt_3sigma_EE->GetBinContent(i)/h_noPU_prompt_3sigma_EE->Integral(1,nbin));
    prompt_eff_noPU_4sigma_EE.emplace_back(h_noPU_prompt_4sigma_EE->GetBinContent(i)/h_noPU_prompt_4sigma_EE->Integral(1,nbin));
            // GEN case
    prompt_eff_noPU_gen_EE.emplace_back(h_noPU_prompt_gen_EE->GetBinContent(i)/h_noPU_prompt_gen_EE->Integral(1,nbin));
    prompt_eff_noPU_gen_2sigma_EE.emplace_back(h_noPU_prompt_gen_2sigma_EE->GetBinContent(i)/h_noPU_prompt_gen_2sigma_EE->Integral(1,nbin));
    prompt_eff_noPU_gen_3sigma_EE.emplace_back(h_noPU_prompt_gen_3sigma_EE->GetBinContent(i)/h_noPU_prompt_gen_3sigma_EE->Integral(1,nbin));
    prompt_eff_noPU_gen_4sigma_EE.emplace_back(h_noPU_prompt_gen_4sigma_EE->GetBinContent(i)/h_noPU_prompt_gen_4sigma_EE->Integral(1,nbin));

    // nonprompt
      // Barrel region
        // efficiency
	  // PU200
    nonprompt_eff_PU200_EB.emplace_back(h_PU200_nonprompt_EB->GetBinContent(i)/h_PU200_nonprompt_EB->Integral(1,nbin));
    nonprompt_eff_PU200_40_EB.emplace_back(h_PU200_nonprompt_40_EB->GetBinContent(i)/h_PU200_nonprompt_40_EB->Integral(1,nbin));
    nonprompt_eff_PU200_60_EB.emplace_back(h_PU200_nonprompt_60_EB->GetBinContent(i)/h_PU200_nonprompt_60_EB->Integral(1,nbin));
    nonprompt_eff_PU200_80_EB.emplace_back(h_PU200_nonprompt_80_EB->GetBinContent(i)/h_PU200_nonprompt_80_EB->Integral(1,nbin));
    nonprompt_eff_PU200_100_EB.emplace_back(h_PU200_nonprompt_100_EB->GetBinContent(i)/h_PU200_nonprompt_100_EB->Integral(1,nbin));
    nonprompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_nonprompt_2sigma_EB->GetBinContent(i)/h_PU200_nonprompt_2sigma_EB->Integral(1,nbin));
    nonprompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_nonprompt_3sigma_EB->GetBinContent(i)/h_PU200_nonprompt_3sigma_EB->Integral(1,nbin));
    nonprompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_nonprompt_4sigma_EB->GetBinContent(i)/h_PU200_nonprompt_4sigma_EB->Integral(1,nbin));
              // vtx
    nonprompt_eff_PU200_2sigma_EB_vtx.emplace_back(h_PU200_nonprompt_2sigma_EB_vtx->GetBinContent(i)/h_PU200_nonprompt_2sigma_EB_vtx->Integral(1,nbin));
    nonprompt_eff_PU200_3sigma_EB_vtx.emplace_back(h_PU200_nonprompt_3sigma_EB_vtx->GetBinContent(i)/h_PU200_nonprompt_3sigma_EB_vtx->Integral(1,nbin));
    nonprompt_eff_PU200_4sigma_EB_vtx.emplace_back(h_PU200_nonprompt_4sigma_EB_vtx->GetBinContent(i)/h_PU200_nonprompt_4sigma_EB_vtx->Integral(1,nbin));
            // GEN case
    nonprompt_eff_PU200_gen_EB.emplace_back(h_PU200_nonprompt_gen_EB->GetBinContent(i)/h_PU200_nonprompt_gen_EB->Integral(1,nbin));
    nonprompt_eff_PU200_gen_2sigma_EB.emplace_back(h_PU200_nonprompt_gen_2sigma_EB->GetBinContent(i)/h_PU200_nonprompt_gen_2sigma_EB->Integral(1,nbin));
    nonprompt_eff_PU200_gen_3sigma_EB.emplace_back(h_PU200_nonprompt_gen_3sigma_EB->GetBinContent(i)/h_PU200_nonprompt_gen_3sigma_EB->Integral(1,nbin));
    nonprompt_eff_PU200_gen_4sigma_EB.emplace_back(h_PU200_nonprompt_gen_4sigma_EB->GetBinContent(i)/h_PU200_nonprompt_gen_4sigma_EB->Integral(1,nbin));
	  // noPU
    nonprompt_eff_noPU_EB.emplace_back(h_noPU_nonprompt_EB->GetBinContent(i)/h_noPU_nonprompt_EB->Integral(1,nbin));
    nonprompt_eff_noPU_40_EB.emplace_back(h_noPU_nonprompt_40_EB->GetBinContent(i)/h_noPU_nonprompt_40_EB->Integral(1,nbin));
    nonprompt_eff_noPU_60_EB.emplace_back(h_noPU_nonprompt_60_EB->GetBinContent(i)/h_noPU_nonprompt_60_EB->Integral(1,nbin));
    nonprompt_eff_noPU_80_EB.emplace_back(h_noPU_nonprompt_80_EB->GetBinContent(i)/h_noPU_nonprompt_80_EB->Integral(1,nbin));
    nonprompt_eff_noPU_100_EB.emplace_back(h_noPU_nonprompt_100_EB->GetBinContent(i)/h_noPU_nonprompt_100_EB->Integral(1,nbin));
    nonprompt_eff_noPU_2sigma_EB.emplace_back(h_noPU_nonprompt_2sigma_EB->GetBinContent(i)/h_noPU_nonprompt_2sigma_EB->Integral(1,nbin));
    nonprompt_eff_noPU_3sigma_EB.emplace_back(h_noPU_nonprompt_3sigma_EB->GetBinContent(i)/h_noPU_nonprompt_3sigma_EB->Integral(1,nbin));
    nonprompt_eff_noPU_4sigma_EB.emplace_back(h_noPU_nonprompt_4sigma_EB->GetBinContent(i)/h_noPU_nonprompt_4sigma_EB->Integral(1,nbin));
            // GEN case
    nonprompt_eff_noPU_gen_EB.emplace_back(h_noPU_nonprompt_gen_EB->GetBinContent(i)/h_noPU_nonprompt_gen_EB->Integral(1,nbin));
    nonprompt_eff_noPU_gen_2sigma_EB.emplace_back(h_noPU_nonprompt_gen_2sigma_EB->GetBinContent(i)/h_noPU_nonprompt_gen_2sigma_EB->Integral(1,nbin));
    nonprompt_eff_noPU_gen_3sigma_EB.emplace_back(h_noPU_nonprompt_gen_3sigma_EB->GetBinContent(i)/h_noPU_nonprompt_gen_3sigma_EB->Integral(1,nbin));
    nonprompt_eff_noPU_gen_4sigma_EB.emplace_back(h_noPU_nonprompt_gen_4sigma_EB->GetBinContent(i)/h_noPU_nonprompt_gen_4sigma_EB->Integral(1,nbin));

      // Endcap region
        // efficiency
	  // PU200
    nonprompt_eff_PU200_EE.emplace_back(h_PU200_nonprompt_EE->GetBinContent(i)/h_PU200_nonprompt_EE->Integral(1,nbin));
    nonprompt_eff_PU200_40_EE.emplace_back(h_PU200_nonprompt_40_EE->GetBinContent(i)/h_PU200_nonprompt_40_EE->Integral(1,nbin));
    nonprompt_eff_PU200_60_EE.emplace_back(h_PU200_nonprompt_60_EE->GetBinContent(i)/h_PU200_nonprompt_60_EE->Integral(1,nbin));
    nonprompt_eff_PU200_80_EE.emplace_back(h_PU200_nonprompt_80_EE->GetBinContent(i)/h_PU200_nonprompt_80_EE->Integral(1,nbin));
    nonprompt_eff_PU200_100_EE.emplace_back(h_PU200_nonprompt_100_EE->GetBinContent(i)/h_PU200_nonprompt_100_EE->Integral(1,nbin));
    nonprompt_eff_PU200_2sigma_EE.emplace_back(h_PU200_nonprompt_2sigma_EE->GetBinContent(i)/h_PU200_nonprompt_2sigma_EE->Integral(1,nbin));
    nonprompt_eff_PU200_3sigma_EE.emplace_back(h_PU200_nonprompt_3sigma_EE->GetBinContent(i)/h_PU200_nonprompt_3sigma_EE->Integral(1,nbin));
    nonprompt_eff_PU200_4sigma_EE.emplace_back(h_PU200_nonprompt_4sigma_EE->GetBinContent(i)/h_PU200_nonprompt_4sigma_EE->Integral(1,nbin));
              // vtx
    nonprompt_eff_PU200_2sigma_EE_vtx.emplace_back(h_PU200_nonprompt_2sigma_EE_vtx->GetBinContent(i)/h_PU200_nonprompt_2sigma_EE_vtx->Integral(1,nbin));
    nonprompt_eff_PU200_3sigma_EE_vtx.emplace_back(h_PU200_nonprompt_3sigma_EE_vtx->GetBinContent(i)/h_PU200_nonprompt_3sigma_EE_vtx->Integral(1,nbin));
    nonprompt_eff_PU200_4sigma_EE_vtx.emplace_back(h_PU200_nonprompt_4sigma_EE_vtx->GetBinContent(i)/h_PU200_nonprompt_4sigma_EE_vtx->Integral(1,nbin));
            // GEN case
    nonprompt_eff_PU200_gen_EE.emplace_back(h_PU200_nonprompt_gen_EE->GetBinContent(i)/h_PU200_nonprompt_gen_EE->Integral(1,nbin));
    nonprompt_eff_PU200_gen_2sigma_EE.emplace_back(h_PU200_nonprompt_gen_2sigma_EE->GetBinContent(i)/h_PU200_nonprompt_gen_2sigma_EE->Integral(1,nbin));
    nonprompt_eff_PU200_gen_3sigma_EE.emplace_back(h_PU200_nonprompt_gen_3sigma_EE->GetBinContent(i)/h_PU200_nonprompt_gen_3sigma_EE->Integral(1,nbin));
    nonprompt_eff_PU200_gen_4sigma_EE.emplace_back(h_PU200_nonprompt_gen_4sigma_EE->GetBinContent(i)/h_PU200_nonprompt_gen_4sigma_EE->Integral(1,nbin));
	  // noPU
    nonprompt_eff_noPU_EE.emplace_back(h_noPU_nonprompt_EE->GetBinContent(i)/h_noPU_nonprompt_EE->Integral(1,nbin));
    nonprompt_eff_noPU_40_EE.emplace_back(h_noPU_nonprompt_40_EE->GetBinContent(i)/h_noPU_nonprompt_40_EE->Integral(1,nbin));
    nonprompt_eff_noPU_60_EE.emplace_back(h_noPU_nonprompt_60_EE->GetBinContent(i)/h_noPU_nonprompt_60_EE->Integral(1,nbin));
    nonprompt_eff_noPU_80_EE.emplace_back(h_noPU_nonprompt_80_EE->GetBinContent(i)/h_noPU_nonprompt_80_EE->Integral(1,nbin));
    nonprompt_eff_noPU_100_EE.emplace_back(h_noPU_nonprompt_100_EE->GetBinContent(i)/h_noPU_nonprompt_100_EE->Integral(1,nbin));
    nonprompt_eff_noPU_2sigma_EE.emplace_back(h_noPU_nonprompt_2sigma_EE->GetBinContent(i)/h_noPU_nonprompt_2sigma_EE->Integral(1,nbin));
    nonprompt_eff_noPU_3sigma_EE.emplace_back(h_noPU_nonprompt_3sigma_EE->GetBinContent(i)/h_noPU_nonprompt_3sigma_EE->Integral(1,nbin));
    nonprompt_eff_noPU_4sigma_EE.emplace_back(h_noPU_nonprompt_4sigma_EE->GetBinContent(i)/h_noPU_nonprompt_4sigma_EE->Integral(1,nbin));
            // GEN case
    nonprompt_eff_noPU_gen_EE.emplace_back(h_noPU_nonprompt_gen_EE->GetBinContent(i)/h_noPU_nonprompt_gen_EE->Integral(1,nbin));
    nonprompt_eff_noPU_gen_2sigma_EE.emplace_back(h_noPU_nonprompt_gen_2sigma_EE->GetBinContent(i)/h_noPU_nonprompt_gen_2sigma_EE->Integral(1,nbin));
    nonprompt_eff_noPU_gen_3sigma_EE.emplace_back(h_noPU_nonprompt_gen_3sigma_EE->GetBinContent(i)/h_noPU_nonprompt_gen_3sigma_EE->Integral(1,nbin));
    nonprompt_eff_noPU_gen_4sigma_EE.emplace_back(h_noPU_nonprompt_gen_4sigma_EE->GetBinContent(i)/h_noPU_nonprompt_gen_4sigma_EE->Integral(1,nbin));

  }

  // Define TGraph
  // Prompt
    // Barrel region
      // PU200
  TGraph* gr_eff_PU200_prompt_EB = new TGraph();            
  TGraph* gr_eff_PU200_prompt_40_EB = new TGraph();         
  TGraph* gr_eff_PU200_prompt_60_EB = new TGraph();         
  TGraph* gr_eff_PU200_prompt_80_EB = new TGraph();         
  TGraph* gr_eff_PU200_prompt_100_EB = new TGraph();        
  TGraph* gr_eff_PU200_prompt_2sigma_EB = new TGraph();     
  TGraph* gr_eff_PU200_prompt_3sigma_EB = new TGraph();     
  TGraph* gr_eff_PU200_prompt_4sigma_EB = new TGraph();     

  TGraph* gr_eff_PU200_prompt_gen_EB = new TGraph();         
  TGraph* gr_eff_PU200_prompt_gen_2sigma_EB = new TGraph();  
  TGraph* gr_eff_PU200_prompt_gen_3sigma_EB = new TGraph();  
  TGraph* gr_eff_PU200_prompt_gen_4sigma_EB = new TGraph();  
        // vtx
  TGraph* gr_eff_PU200_prompt_2sigma_EB_vtx = new TGraph();  
  TGraph* gr_eff_PU200_prompt_3sigma_EB_vtx = new TGraph();  
  TGraph* gr_eff_PU200_prompt_4sigma_EB_vtx = new TGraph();  
      // noPU
  TGraph* gr_eff_noPU_prompt_EB = new TGraph();             
  TGraph* gr_eff_noPU_prompt_40_EB = new TGraph();          
  TGraph* gr_eff_noPU_prompt_60_EB = new TGraph();          
  TGraph* gr_eff_noPU_prompt_80_EB = new TGraph();          
  TGraph* gr_eff_noPU_prompt_100_EB = new TGraph();         
  TGraph* gr_eff_noPU_prompt_2sigma_EB = new TGraph();      
  TGraph* gr_eff_noPU_prompt_3sigma_EB = new TGraph();      
  TGraph* gr_eff_noPU_prompt_4sigma_EB = new TGraph();      

  TGraph* gr_eff_noPU_prompt_gen_EB = new TGraph();         
  TGraph* gr_eff_noPU_prompt_gen_2sigma_EB = new TGraph();  
  TGraph* gr_eff_noPU_prompt_gen_3sigma_EB = new TGraph();  
  TGraph* gr_eff_noPU_prompt_gen_4sigma_EB = new TGraph();  

    // Endcap region
      // PU200
  TGraph* gr_eff_PU200_prompt_EE = new TGraph();           
  TGraph* gr_eff_PU200_prompt_40_EE = new TGraph();        
  TGraph* gr_eff_PU200_prompt_60_EE = new TGraph();        
  TGraph* gr_eff_PU200_prompt_80_EE = new TGraph();        
  TGraph* gr_eff_PU200_prompt_100_EE = new TGraph();       
  TGraph* gr_eff_PU200_prompt_2sigma_EE = new TGraph();    
  TGraph* gr_eff_PU200_prompt_3sigma_EE = new TGraph();    
  TGraph* gr_eff_PU200_prompt_4sigma_EE = new TGraph();    

  TGraph* gr_eff_PU200_prompt_gen_EE = new TGraph();          
  TGraph* gr_eff_PU200_prompt_gen_2sigma_EE = new TGraph();   
  TGraph* gr_eff_PU200_prompt_gen_3sigma_EE = new TGraph();   
  TGraph* gr_eff_PU200_prompt_gen_4sigma_EE = new TGraph();   
        // vtx
  TGraph* gr_eff_PU200_prompt_2sigma_EE_vtx = new TGraph();  
  TGraph* gr_eff_PU200_prompt_3sigma_EE_vtx = new TGraph();  
  TGraph* gr_eff_PU200_prompt_4sigma_EE_vtx = new TGraph();  
      // noPU
  TGraph* gr_eff_noPU_prompt_EE = new TGraph();             
  TGraph* gr_eff_noPU_prompt_40_EE = new TGraph();          
  TGraph* gr_eff_noPU_prompt_60_EE = new TGraph();          
  TGraph* gr_eff_noPU_prompt_80_EE = new TGraph();          
  TGraph* gr_eff_noPU_prompt_100_EE = new TGraph();         
  TGraph* gr_eff_noPU_prompt_2sigma_EE = new TGraph();      
  TGraph* gr_eff_noPU_prompt_3sigma_EE = new TGraph();      
  TGraph* gr_eff_noPU_prompt_4sigma_EE = new TGraph();      

  TGraph* gr_eff_noPU_prompt_gen_EE = new TGraph();        
  TGraph* gr_eff_noPU_prompt_gen_2sigma_EE = new TGraph(); 
  TGraph* gr_eff_noPU_prompt_gen_3sigma_EE = new TGraph(); 
  TGraph* gr_eff_noPU_prompt_gen_4sigma_EE = new TGraph(); 

  // Nonprompt
    // Barrel region
      // PU200
  TGraph* gr_eff_PU200_nonprompt_EB = new TGraph();       
  TGraph* gr_eff_PU200_nonprompt_40_EB = new TGraph();    
  TGraph* gr_eff_PU200_nonprompt_60_EB = new TGraph();    
  TGraph* gr_eff_PU200_nonprompt_80_EB = new TGraph();    
  TGraph* gr_eff_PU200_nonprompt_100_EB = new TGraph();   
  TGraph* gr_eff_PU200_nonprompt_2sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_3sigma_EB = new TGraph();
  TGraph* gr_eff_PU200_nonprompt_4sigma_EB = new TGraph();

  TGraph* gr_eff_PU200_nonprompt_gen_EB = new TGraph();          
  TGraph* gr_eff_PU200_nonprompt_gen_2sigma_EB = new TGraph();   
  TGraph* gr_eff_PU200_nonprompt_gen_3sigma_EB = new TGraph();   
  TGraph* gr_eff_PU200_nonprompt_gen_4sigma_EB = new TGraph();   
        // vtx
  TGraph* gr_eff_PU200_nonprompt_2sigma_EB_vtx = new TGraph();   
  TGraph* gr_eff_PU200_nonprompt_3sigma_EB_vtx = new TGraph();   
  TGraph* gr_eff_PU200_nonprompt_4sigma_EB_vtx = new TGraph();   
      // noPU
  TGraph* gr_eff_noPU_nonprompt_EB = new TGraph();        
  TGraph* gr_eff_noPU_nonprompt_40_EB = new TGraph();     
  TGraph* gr_eff_noPU_nonprompt_60_EB = new TGraph();     
  TGraph* gr_eff_noPU_nonprompt_80_EB = new TGraph();     
  TGraph* gr_eff_noPU_nonprompt_100_EB = new TGraph();    
  TGraph* gr_eff_noPU_nonprompt_2sigma_EB = new TGraph(); 
  TGraph* gr_eff_noPU_nonprompt_3sigma_EB = new TGraph(); 
  TGraph* gr_eff_noPU_nonprompt_4sigma_EB = new TGraph(); 

  TGraph* gr_eff_noPU_nonprompt_gen_EB = new TGraph();          
  TGraph* gr_eff_noPU_nonprompt_gen_2sigma_EB = new TGraph();   
  TGraph* gr_eff_noPU_nonprompt_gen_3sigma_EB = new TGraph();   
  TGraph* gr_eff_noPU_nonprompt_gen_4sigma_EB = new TGraph();   

    // Endcap region
      // PU200
  TGraph* gr_eff_PU200_nonprompt_EE = new TGraph();            
  TGraph* gr_eff_PU200_nonprompt_40_EE = new TGraph();         
  TGraph* gr_eff_PU200_nonprompt_60_EE = new TGraph();         
  TGraph* gr_eff_PU200_nonprompt_80_EE = new TGraph();         
  TGraph* gr_eff_PU200_nonprompt_100_EE = new TGraph();        
  TGraph* gr_eff_PU200_nonprompt_2sigma_EE = new TGraph();     
  TGraph* gr_eff_PU200_nonprompt_3sigma_EE = new TGraph();     
  TGraph* gr_eff_PU200_nonprompt_4sigma_EE = new TGraph();     

  TGraph* gr_eff_PU200_nonprompt_gen_EE = new TGraph();         
  TGraph* gr_eff_PU200_nonprompt_gen_2sigma_EE = new TGraph();  
  TGraph* gr_eff_PU200_nonprompt_gen_3sigma_EE = new TGraph();  
  TGraph* gr_eff_PU200_nonprompt_gen_4sigma_EE = new TGraph();  
        // vtx
  TGraph* gr_eff_PU200_nonprompt_2sigma_EE_vtx = new TGraph();  
  TGraph* gr_eff_PU200_nonprompt_3sigma_EE_vtx = new TGraph();  
  TGraph* gr_eff_PU200_nonprompt_4sigma_EE_vtx = new TGraph();  
      // noPU
  TGraph* gr_eff_noPU_nonprompt_EE = new TGraph();            
  TGraph* gr_eff_noPU_nonprompt_40_EE = new TGraph();         
  TGraph* gr_eff_noPU_nonprompt_60_EE = new TGraph();         
  TGraph* gr_eff_noPU_nonprompt_80_EE = new TGraph();         
  TGraph* gr_eff_noPU_nonprompt_100_EE = new TGraph();        
  TGraph* gr_eff_noPU_nonprompt_2sigma_EE = new TGraph();     
  TGraph* gr_eff_noPU_nonprompt_3sigma_EE = new TGraph();     
  TGraph* gr_eff_noPU_nonprompt_4sigma_EE = new TGraph();     

  TGraph* gr_eff_noPU_nonprompt_gen_EE = new TGraph();       
  TGraph* gr_eff_noPU_nonprompt_gen_2sigma_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_gen_3sigma_EE = new TGraph();
  TGraph* gr_eff_noPU_nonprompt_gen_4sigma_EE = new TGraph();


  for(unsigned int i=0; i<1000; i++) {  // Store efficiency up to rel. iso. cut==4
  //for(unsigned int i=0; i<500; i++) {  // Store efficiency up to rel. iso. cut==4
  //for(unsigned int i=0; i<63; i++) {      // Store efficiency up to rel. iso. cut==0.25
  //for(unsigned int i=0; i<4; i++) {      // 
  // Prompt
    // Barrel region
      // efficiency
        // PU200
    gr_eff_PU200_prompt_EB->SetPoint(gr_eff_PU200_prompt_EB->GetN(), 0.004*i, prompt_eff_PU200_EB.at(i+2));
    gr_eff_PU200_prompt_40_EB->SetPoint(gr_eff_PU200_prompt_40_EB->GetN(), 0.004*i, prompt_eff_PU200_40_EB.at(i+2));
    gr_eff_PU200_prompt_60_EB->SetPoint(gr_eff_PU200_prompt_60_EB->GetN(), 0.004*i, prompt_eff_PU200_60_EB.at(i+2));
    gr_eff_PU200_prompt_80_EB->SetPoint(gr_eff_PU200_prompt_80_EB->GetN(), 0.004*i, prompt_eff_PU200_80_EB.at(i+2));
    gr_eff_PU200_prompt_100_EB->SetPoint(gr_eff_PU200_prompt_100_EB->GetN(), 0.004*i, prompt_eff_PU200_100_EB.at(i+2));
    gr_eff_PU200_prompt_2sigma_EB->SetPoint(gr_eff_PU200_prompt_2sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EB.at(i+2));
    gr_eff_PU200_prompt_3sigma_EB->SetPoint(gr_eff_PU200_prompt_3sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EB.at(i+2));
    gr_eff_PU200_prompt_4sigma_EB->SetPoint(gr_eff_PU200_prompt_4sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EB.at(i+2));

    gr_eff_PU200_prompt_gen_EB->SetPoint(gr_eff_PU200_prompt_gen_EB->GetN(), 0.004*i, prompt_eff_PU200_gen_EB.at(i+2));
    gr_eff_PU200_prompt_gen_2sigma_EB->SetPoint(gr_eff_PU200_prompt_gen_2sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_gen_2sigma_EB.at(i+2));
    gr_eff_PU200_prompt_gen_3sigma_EB->SetPoint(gr_eff_PU200_prompt_gen_3sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_gen_3sigma_EB.at(i+2));
    gr_eff_PU200_prompt_gen_4sigma_EB->SetPoint(gr_eff_PU200_prompt_gen_4sigma_EB->GetN(), 0.004*i, prompt_eff_PU200_gen_4sigma_EB.at(i+2));
          // vtx
    gr_eff_PU200_prompt_2sigma_EB_vtx->SetPoint(gr_eff_PU200_prompt_2sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EB_vtx.at(i+2));
    gr_eff_PU200_prompt_3sigma_EB_vtx->SetPoint(gr_eff_PU200_prompt_3sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EB_vtx.at(i+2));
    gr_eff_PU200_prompt_4sigma_EB_vtx->SetPoint(gr_eff_PU200_prompt_4sigma_EB_vtx->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EB_vtx.at(i+2));
        // noPU
    gr_eff_noPU_prompt_EB->SetPoint(gr_eff_noPU_prompt_EB->GetN(), 0.004*i, prompt_eff_noPU_EB.at(i+2));
    gr_eff_noPU_prompt_40_EB->SetPoint(gr_eff_noPU_prompt_40_EB->GetN(), 0.004*i, prompt_eff_noPU_40_EB.at(i+2));
    gr_eff_noPU_prompt_60_EB->SetPoint(gr_eff_noPU_prompt_60_EB->GetN(), 0.004*i, prompt_eff_noPU_60_EB.at(i+2));
    gr_eff_noPU_prompt_80_EB->SetPoint(gr_eff_noPU_prompt_80_EB->GetN(), 0.004*i, prompt_eff_noPU_80_EB.at(i+2));
    gr_eff_noPU_prompt_100_EB->SetPoint(gr_eff_noPU_prompt_100_EB->GetN(), 0.004*i, prompt_eff_noPU_100_EB.at(i+2));
    gr_eff_noPU_prompt_2sigma_EB->SetPoint(gr_eff_noPU_prompt_2sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_2sigma_EB.at(i+2));
    gr_eff_noPU_prompt_3sigma_EB->SetPoint(gr_eff_noPU_prompt_3sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_3sigma_EB.at(i+2));
    gr_eff_noPU_prompt_4sigma_EB->SetPoint(gr_eff_noPU_prompt_4sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_4sigma_EB.at(i+2));

    gr_eff_noPU_prompt_gen_EB->SetPoint(gr_eff_noPU_prompt_gen_EB->GetN(), 0.004*i, prompt_eff_noPU_gen_EB.at(i+2));
    gr_eff_noPU_prompt_gen_2sigma_EB->SetPoint(gr_eff_noPU_prompt_gen_2sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_gen_2sigma_EB.at(i+2));
    gr_eff_noPU_prompt_gen_3sigma_EB->SetPoint(gr_eff_noPU_prompt_gen_3sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_gen_3sigma_EB.at(i+2));
    gr_eff_noPU_prompt_gen_4sigma_EB->SetPoint(gr_eff_noPU_prompt_gen_4sigma_EB->GetN(), 0.004*i, prompt_eff_noPU_gen_4sigma_EB.at(i+2));

    // Endcap region
      // efficiency
        // PU200
    gr_eff_PU200_prompt_EE->SetPoint(gr_eff_PU200_prompt_EE->GetN(), 0.004*i, prompt_eff_PU200_EE.at(i+2));
    gr_eff_PU200_prompt_40_EE->SetPoint(gr_eff_PU200_prompt_40_EE->GetN(), 0.004*i, prompt_eff_PU200_40_EE.at(i+2));
    gr_eff_PU200_prompt_60_EE->SetPoint(gr_eff_PU200_prompt_60_EE->GetN(), 0.004*i, prompt_eff_PU200_60_EE.at(i+2));
    gr_eff_PU200_prompt_80_EE->SetPoint(gr_eff_PU200_prompt_80_EE->GetN(), 0.004*i, prompt_eff_PU200_80_EE.at(i+2));
    gr_eff_PU200_prompt_100_EE->SetPoint(gr_eff_PU200_prompt_100_EE->GetN(), 0.004*i, prompt_eff_PU200_100_EE.at(i+2));
    gr_eff_PU200_prompt_2sigma_EE->SetPoint(gr_eff_PU200_prompt_2sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EE.at(i+2));
    gr_eff_PU200_prompt_3sigma_EE->SetPoint(gr_eff_PU200_prompt_3sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EE.at(i+2));
    gr_eff_PU200_prompt_4sigma_EE->SetPoint(gr_eff_PU200_prompt_4sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EE.at(i+2));

    gr_eff_PU200_prompt_gen_EE->SetPoint(gr_eff_PU200_prompt_gen_EE->GetN(), 0.004*i, prompt_eff_PU200_gen_EE.at(i+2));
    gr_eff_PU200_prompt_gen_2sigma_EE->SetPoint(gr_eff_PU200_prompt_gen_2sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_gen_2sigma_EE.at(i+2));
    gr_eff_PU200_prompt_gen_3sigma_EE->SetPoint(gr_eff_PU200_prompt_gen_3sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_gen_3sigma_EE.at(i+2));
    gr_eff_PU200_prompt_gen_4sigma_EE->SetPoint(gr_eff_PU200_prompt_gen_4sigma_EE->GetN(), 0.004*i, prompt_eff_PU200_gen_4sigma_EE.at(i+2));
          // vtx
    gr_eff_PU200_prompt_2sigma_EE_vtx->SetPoint(gr_eff_PU200_prompt_2sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_PU200_2sigma_EE_vtx.at(i+2));
    gr_eff_PU200_prompt_3sigma_EE_vtx->SetPoint(gr_eff_PU200_prompt_3sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_PU200_3sigma_EE_vtx.at(i+2));
    gr_eff_PU200_prompt_4sigma_EE_vtx->SetPoint(gr_eff_PU200_prompt_4sigma_EE_vtx->GetN(), 0.004*i, prompt_eff_PU200_4sigma_EE_vtx.at(i+2));
        // noPU
    gr_eff_noPU_prompt_EE->SetPoint(gr_eff_noPU_prompt_EE->GetN(), 0.004*i, prompt_eff_noPU_EE.at(i+2));
    gr_eff_noPU_prompt_40_EE->SetPoint(gr_eff_noPU_prompt_40_EE->GetN(), 0.004*i, prompt_eff_noPU_40_EE.at(i+2));
    gr_eff_noPU_prompt_60_EE->SetPoint(gr_eff_noPU_prompt_60_EE->GetN(), 0.004*i, prompt_eff_noPU_60_EE.at(i+2));
    gr_eff_noPU_prompt_80_EE->SetPoint(gr_eff_noPU_prompt_80_EE->GetN(), 0.004*i, prompt_eff_noPU_80_EE.at(i+2));
    gr_eff_noPU_prompt_100_EE->SetPoint(gr_eff_noPU_prompt_100_EE->GetN(), 0.004*i, prompt_eff_noPU_100_EE.at(i+2));
    gr_eff_noPU_prompt_2sigma_EE->SetPoint(gr_eff_noPU_prompt_2sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_2sigma_EE.at(i+2));
    gr_eff_noPU_prompt_3sigma_EE->SetPoint(gr_eff_noPU_prompt_3sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_3sigma_EE.at(i+2));
    gr_eff_noPU_prompt_4sigma_EE->SetPoint(gr_eff_noPU_prompt_4sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_4sigma_EE.at(i+2));

    gr_eff_noPU_prompt_gen_EE->SetPoint(gr_eff_noPU_prompt_gen_EE->GetN(), 0.004*i, prompt_eff_noPU_gen_EE.at(i+2));
    gr_eff_noPU_prompt_gen_2sigma_EE->SetPoint(gr_eff_noPU_prompt_gen_2sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_gen_2sigma_EE.at(i+2));
    gr_eff_noPU_prompt_gen_3sigma_EE->SetPoint(gr_eff_noPU_prompt_gen_3sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_gen_3sigma_EE.at(i+2));
    gr_eff_noPU_prompt_gen_4sigma_EE->SetPoint(gr_eff_noPU_prompt_gen_4sigma_EE->GetN(), 0.004*i, prompt_eff_noPU_gen_4sigma_EE.at(i+2));
  }

  for(unsigned int i=0; i<1000; i++) {      //
  //for(unsigned int i=0; i<500; i++) {      //
  //for(unsigned int i=0; i<63; i++) {      // Store efficiency up to rel. iso. cut==0.25
  // Nonprompt
    // Barrel region
      // efficiency
        // PU200
    gr_eff_PU200_nonprompt_EB->SetPoint(gr_eff_PU200_nonprompt_EB->GetN(), 0.004*i, nonprompt_eff_PU200_EB.at(i+2));
    gr_eff_PU200_nonprompt_40_EB->SetPoint(gr_eff_PU200_nonprompt_40_EB->GetN(), 0.004*i, nonprompt_eff_PU200_40_EB.at(i+2));
    gr_eff_PU200_nonprompt_60_EB->SetPoint(gr_eff_PU200_nonprompt_60_EB->GetN(), 0.004*i, nonprompt_eff_PU200_60_EB.at(i+2));
    gr_eff_PU200_nonprompt_80_EB->SetPoint(gr_eff_PU200_nonprompt_80_EB->GetN(), 0.004*i, nonprompt_eff_PU200_80_EB.at(i+2));
    gr_eff_PU200_nonprompt_100_EB->SetPoint(gr_eff_PU200_nonprompt_100_EB->GetN(), 0.004*i, nonprompt_eff_PU200_100_EB.at(i+2));
    gr_eff_PU200_nonprompt_2sigma_EB->SetPoint(gr_eff_PU200_nonprompt_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EB.at(i+2));
    gr_eff_PU200_nonprompt_3sigma_EB->SetPoint(gr_eff_PU200_nonprompt_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EB.at(i+2));
    gr_eff_PU200_nonprompt_4sigma_EB->SetPoint(gr_eff_PU200_nonprompt_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EB.at(i+2));

    gr_eff_PU200_nonprompt_gen_EB->SetPoint(gr_eff_PU200_nonprompt_gen_EB->GetN(), 0.004*i, nonprompt_eff_PU200_gen_EB.at(i+2));
    gr_eff_PU200_nonprompt_gen_2sigma_EB->SetPoint(gr_eff_PU200_nonprompt_gen_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_gen_2sigma_EB.at(i+2));
    gr_eff_PU200_nonprompt_gen_3sigma_EB->SetPoint(gr_eff_PU200_nonprompt_gen_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_gen_3sigma_EB.at(i+2));
    gr_eff_PU200_nonprompt_gen_4sigma_EB->SetPoint(gr_eff_PU200_nonprompt_gen_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_PU200_gen_4sigma_EB.at(i+2));
          // vtx
    gr_eff_PU200_nonprompt_2sigma_EB_vtx->SetPoint(gr_eff_PU200_nonprompt_2sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EB_vtx.at(i+2));
    gr_eff_PU200_nonprompt_3sigma_EB_vtx->SetPoint(gr_eff_PU200_nonprompt_3sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EB_vtx.at(i+2));
    gr_eff_PU200_nonprompt_4sigma_EB_vtx->SetPoint(gr_eff_PU200_nonprompt_4sigma_EB_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EB_vtx.at(i+2));
        // noPU
    gr_eff_noPU_nonprompt_EB->SetPoint(gr_eff_noPU_nonprompt_EB->GetN(), 0.004*i, nonprompt_eff_noPU_EB.at(i+2));
    gr_eff_noPU_nonprompt_40_EB->SetPoint(gr_eff_noPU_nonprompt_40_EB->GetN(), 0.004*i, nonprompt_eff_noPU_40_EB.at(i+2));
    gr_eff_noPU_nonprompt_60_EB->SetPoint(gr_eff_noPU_nonprompt_60_EB->GetN(), 0.004*i, nonprompt_eff_noPU_60_EB.at(i+2));
    gr_eff_noPU_nonprompt_80_EB->SetPoint(gr_eff_noPU_nonprompt_80_EB->GetN(), 0.004*i, nonprompt_eff_noPU_80_EB.at(i+2));
    gr_eff_noPU_nonprompt_100_EB->SetPoint(gr_eff_noPU_nonprompt_100_EB->GetN(), 0.004*i, nonprompt_eff_noPU_100_EB.at(i+2));
    gr_eff_noPU_nonprompt_2sigma_EB->SetPoint(gr_eff_noPU_nonprompt_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_2sigma_EB.at(i+2));
    gr_eff_noPU_nonprompt_3sigma_EB->SetPoint(gr_eff_noPU_nonprompt_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_3sigma_EB.at(i+2));
    gr_eff_noPU_nonprompt_4sigma_EB->SetPoint(gr_eff_noPU_nonprompt_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_4sigma_EB.at(i+2));

    gr_eff_noPU_nonprompt_gen_EB->SetPoint(gr_eff_noPU_nonprompt_gen_EB->GetN(), 0.004*i, nonprompt_eff_noPU_gen_EB.at(i+2));
    gr_eff_noPU_nonprompt_gen_2sigma_EB->SetPoint(gr_eff_noPU_nonprompt_gen_2sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_gen_2sigma_EB.at(i+2));
    gr_eff_noPU_nonprompt_gen_3sigma_EB->SetPoint(gr_eff_noPU_nonprompt_gen_3sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_gen_3sigma_EB.at(i+2));
    gr_eff_noPU_nonprompt_gen_4sigma_EB->SetPoint(gr_eff_noPU_nonprompt_gen_4sigma_EB->GetN(), 0.004*i, nonprompt_eff_noPU_gen_4sigma_EB.at(i+2));

    // Endcap region
      // efficiency
        // PU200
    gr_eff_PU200_nonprompt_EE->SetPoint(gr_eff_PU200_nonprompt_EE->GetN(), 0.004*i, nonprompt_eff_PU200_EE.at(i+2));
    gr_eff_PU200_nonprompt_40_EE->SetPoint(gr_eff_PU200_nonprompt_40_EE->GetN(), 0.004*i, nonprompt_eff_PU200_40_EE.at(i+2));
    gr_eff_PU200_nonprompt_60_EE->SetPoint(gr_eff_PU200_nonprompt_60_EE->GetN(), 0.004*i, nonprompt_eff_PU200_60_EE.at(i+2));
    gr_eff_PU200_nonprompt_80_EE->SetPoint(gr_eff_PU200_nonprompt_80_EE->GetN(), 0.004*i, nonprompt_eff_PU200_80_EE.at(i+2));
    gr_eff_PU200_nonprompt_100_EE->SetPoint(gr_eff_PU200_nonprompt_100_EE->GetN(), 0.004*i, nonprompt_eff_PU200_100_EE.at(i+2));
    gr_eff_PU200_nonprompt_2sigma_EE->SetPoint(gr_eff_PU200_nonprompt_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EE.at(i+2));
    gr_eff_PU200_nonprompt_3sigma_EE->SetPoint(gr_eff_PU200_nonprompt_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EE.at(i+2));
    gr_eff_PU200_nonprompt_4sigma_EE->SetPoint(gr_eff_PU200_nonprompt_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EE.at(i+2));

    gr_eff_PU200_nonprompt_gen_EE->SetPoint(gr_eff_PU200_nonprompt_gen_EE->GetN(), 0.004*i, nonprompt_eff_PU200_gen_EE.at(i+2));
    gr_eff_PU200_nonprompt_gen_2sigma_EE->SetPoint(gr_eff_PU200_nonprompt_gen_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_gen_2sigma_EE.at(i+2));
    gr_eff_PU200_nonprompt_gen_3sigma_EE->SetPoint(gr_eff_PU200_nonprompt_gen_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_gen_3sigma_EE.at(i+2));
    gr_eff_PU200_nonprompt_gen_4sigma_EE->SetPoint(gr_eff_PU200_nonprompt_gen_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_PU200_gen_4sigma_EE.at(i+2));
          // vtx
    gr_eff_PU200_nonprompt_2sigma_EE_vtx->SetPoint(gr_eff_PU200_nonprompt_2sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_2sigma_EE_vtx.at(i+2));
    gr_eff_PU200_nonprompt_3sigma_EE_vtx->SetPoint(gr_eff_PU200_nonprompt_3sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_3sigma_EE_vtx.at(i+2));
    gr_eff_PU200_nonprompt_4sigma_EE_vtx->SetPoint(gr_eff_PU200_nonprompt_4sigma_EE_vtx->GetN(), 0.004*i, nonprompt_eff_PU200_4sigma_EE_vtx.at(i+2));
        // noPU
    gr_eff_noPU_nonprompt_EE->SetPoint(gr_eff_noPU_nonprompt_EE->GetN(), 0.004*i, nonprompt_eff_noPU_EE.at(i+2));
    gr_eff_noPU_nonprompt_40_EE->SetPoint(gr_eff_noPU_nonprompt_40_EE->GetN(), 0.004*i, nonprompt_eff_noPU_40_EE.at(i+2));
    gr_eff_noPU_nonprompt_60_EE->SetPoint(gr_eff_noPU_nonprompt_60_EE->GetN(), 0.004*i, nonprompt_eff_noPU_60_EE.at(i+2));
    gr_eff_noPU_nonprompt_80_EE->SetPoint(gr_eff_noPU_nonprompt_80_EE->GetN(), 0.004*i, nonprompt_eff_noPU_80_EE.at(i+2));
    gr_eff_noPU_nonprompt_100_EE->SetPoint(gr_eff_noPU_nonprompt_100_EE->GetN(), 0.004*i, nonprompt_eff_noPU_100_EE.at(i+2));
    gr_eff_noPU_nonprompt_2sigma_EE->SetPoint(gr_eff_noPU_nonprompt_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_2sigma_EE.at(i+2));
    gr_eff_noPU_nonprompt_3sigma_EE->SetPoint(gr_eff_noPU_nonprompt_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_3sigma_EE.at(i+2));
    gr_eff_noPU_nonprompt_4sigma_EE->SetPoint(gr_eff_noPU_nonprompt_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_4sigma_EE.at(i+2));

    gr_eff_noPU_nonprompt_gen_EE->SetPoint(gr_eff_noPU_nonprompt_gen_EE->GetN(), 0.004*i, nonprompt_eff_noPU_gen_EE.at(i+2));
    gr_eff_noPU_nonprompt_gen_2sigma_EE->SetPoint(gr_eff_noPU_nonprompt_gen_2sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_gen_2sigma_EE.at(i+2));
    gr_eff_noPU_nonprompt_gen_3sigma_EE->SetPoint(gr_eff_noPU_nonprompt_gen_3sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_gen_3sigma_EE.at(i+2));
    gr_eff_noPU_nonprompt_gen_4sigma_EE->SetPoint(gr_eff_noPU_nonprompt_gen_4sigma_EE->GetN(), 0.004*i, nonprompt_eff_noPU_gen_4sigma_EE.at(i+2));
  }

  ///////////////
  // Cosmetics //
  ///////////////
  // Prompt
    // Barrel region
      // efficiency
  gr_eff_PU200_prompt_EB->SetTitle("Prompt muon in barrel region"); gr_eff_noPU_prompt_EB->SetTitle("Prompt muon in barrel region");
  gr_eff_PU200_prompt_EB->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_prompt_EB->GetXaxis()->SetTitle("Relative isolation");
  gr_eff_PU200_prompt_EB->GetYaxis()->SetTitle("% Counts"); gr_eff_noPU_prompt_EB->GetYaxis()->SetTitle("% Counts");
  gr_eff_PU200_prompt_EB->SetLineColor(kBlack); gr_eff_PU200_prompt_gen_EB->SetLineColor(kBlack);
  gr_eff_PU200_prompt_40_EB->SetLineColor(kRed); gr_eff_PU200_prompt_60_EB->SetLineColor(kGreen); gr_eff_PU200_prompt_80_EB->SetLineColor(kBlue); gr_eff_PU200_prompt_100_EB->SetLineColor(kMagenta);
  gr_eff_PU200_prompt_2sigma_EB->SetLineColor(kRed); gr_eff_PU200_prompt_3sigma_EB->SetLineColor(kGreen); gr_eff_PU200_prompt_4sigma_EB->SetLineColor(kBlue);
  gr_eff_noPU_prompt_EB->SetLineColor(kGray+1); gr_eff_noPU_prompt_gen_EB->SetLineColor(kGray+1);
  gr_eff_PU200_prompt_EB->SetLineWidth(2); gr_eff_PU200_prompt_gen_EB->SetLineWidth(2);
  gr_eff_PU200_prompt_40_EB->SetLineWidth(2); gr_eff_PU200_prompt_60_EB->SetLineWidth(2); gr_eff_PU200_prompt_80_EB->SetLineWidth(2); gr_eff_PU200_prompt_100_EB->SetLineWidth(2);
  gr_eff_PU200_prompt_2sigma_EB->SetLineWidth(2); gr_eff_PU200_prompt_3sigma_EB->SetLineWidth(2); gr_eff_PU200_prompt_4sigma_EB->SetLineWidth(2);
  gr_eff_noPU_prompt_EB->SetLineWidth(2); gr_eff_noPU_prompt_gen_EB->SetLineWidth(2);
  gr_eff_PU200_prompt_gen_EB->SetLineStyle(7); gr_eff_noPU_prompt_gen_EB->SetLineStyle(7);
  gr_eff_PU200_prompt_gen_2sigma_EB->SetLineWidth(2); gr_eff_noPU_prompt_gen_2sigma_EB->SetLineWidth(2);
  gr_eff_PU200_prompt_gen_2sigma_EB->SetLineStyle(7); gr_eff_noPU_prompt_gen_2sigma_EB->SetLineStyle(7);
  gr_eff_PU200_prompt_gen_2sigma_EB->SetLineColor(kBlack); gr_eff_noPU_prompt_gen_2sigma_EB->SetLineColor(kGray+1);
  gr_eff_PU200_prompt_2sigma_EB_vtx->SetLineColor(kGreen+2); gr_eff_PU200_prompt_2sigma_EB_vtx->SetLineWidth(2);

  gr_eff_PU200_prompt_EB->SetMarkerStyle(8); gr_eff_PU200_prompt_EB->SetMarkerSize(1.5);
  gr_eff_noPU_prompt_EB->SetMarkerStyle(8); gr_eff_noPU_prompt_EB->SetMarkerSize(1.5); gr_eff_noPU_prompt_EB->SetMarkerColor(kGray+1);

    // Endcap region
      // efficiency
  gr_eff_PU200_prompt_EE->SetTitle("Prompt muon in endcap region"); gr_eff_noPU_prompt_EE->SetTitle("Prompt muon in endcap region");
  gr_eff_PU200_prompt_EE->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_prompt_EE->GetXaxis()->SetTitle("Relative isolation");
  gr_eff_PU200_prompt_EE->GetYaxis()->SetTitle("% Counts"); gr_eff_noPU_prompt_EE->GetYaxis()->SetTitle("% Counts");
  gr_eff_PU200_prompt_EE->SetLineColor(kBlack); gr_eff_PU200_prompt_gen_EE->SetLineColor(kBlack);
  gr_eff_PU200_prompt_40_EE->SetLineColor(kRed); gr_eff_PU200_prompt_60_EE->SetLineColor(kGreen); gr_eff_PU200_prompt_80_EE->SetLineColor(kBlue); gr_eff_PU200_prompt_100_EE->SetLineColor(kMagenta);
  gr_eff_PU200_prompt_2sigma_EE->SetLineColor(kRed); gr_eff_PU200_prompt_3sigma_EE->SetLineColor(kGreen); gr_eff_PU200_prompt_4sigma_EE->SetLineColor(kBlue);
  gr_eff_noPU_prompt_EE->SetLineColor(kGray+1); gr_eff_noPU_prompt_gen_EE->SetLineColor(kGray+1);
  gr_eff_PU200_prompt_EE->SetLineWidth(2); gr_eff_PU200_prompt_gen_EE->SetLineWidth(2);
  gr_eff_PU200_prompt_40_EE->SetLineWidth(2); gr_eff_PU200_prompt_60_EE->SetLineWidth(2); gr_eff_PU200_prompt_80_EE->SetLineWidth(2); gr_eff_PU200_prompt_100_EE->SetLineWidth(2);
  gr_eff_PU200_prompt_2sigma_EE->SetLineWidth(2); gr_eff_PU200_prompt_3sigma_EE->SetLineWidth(2); gr_eff_PU200_prompt_4sigma_EE->SetLineWidth(2);
  gr_eff_noPU_prompt_EE->SetLineWidth(2); gr_eff_noPU_prompt_gen_EE->SetLineWidth(2);
  gr_eff_PU200_prompt_gen_EE->SetLineStyle(7); gr_eff_noPU_prompt_gen_EE->SetLineStyle(7);
  gr_eff_PU200_prompt_gen_2sigma_EE->SetLineWidth(2); gr_eff_noPU_prompt_gen_2sigma_EE->SetLineWidth(2);
  gr_eff_PU200_prompt_gen_2sigma_EE->SetLineStyle(7); gr_eff_noPU_prompt_gen_2sigma_EE->SetLineStyle(7);
  gr_eff_PU200_prompt_gen_2sigma_EE->SetLineColor(kBlack); gr_eff_noPU_prompt_gen_2sigma_EE->SetLineColor(kGray+1);
  gr_eff_PU200_prompt_2sigma_EE_vtx->SetLineColor(kGreen+2); gr_eff_PU200_prompt_2sigma_EE_vtx->SetLineWidth(2);

  gr_eff_PU200_prompt_EE->SetMarkerStyle(8); gr_eff_PU200_prompt_EE->SetMarkerSize(1.5);
  gr_eff_noPU_prompt_EE->SetMarkerStyle(8); gr_eff_noPU_prompt_EE->SetMarkerSize(1.5); gr_eff_noPU_prompt_EE->SetMarkerColor(kGray+1);

  // Nonprompt
    // Barrel region
      // efficiency
  gr_eff_noPU_nonprompt_gen_2sigma_EB->SetTitle("Non-prompt muon in barrel region"); gr_eff_noPU_nonprompt_gen_2sigma_EB->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_nonprompt_gen_2sigma_EB->GetYaxis()->SetTitle("% Counts");
  gr_eff_PU200_nonprompt_EB->SetTitle("Non-prompt muon in barrel region"); gr_eff_noPU_nonprompt_EB->SetTitle("Non-prompt muon in barrel region");
  gr_eff_PU200_nonprompt_EB->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_nonprompt_EB->GetXaxis()->SetTitle("Relative isolation");
  gr_eff_PU200_nonprompt_EB->GetYaxis()->SetTitle("% Counts"); gr_eff_noPU_nonprompt_EB->GetYaxis()->SetTitle("% Counts");
  gr_eff_PU200_nonprompt_EB->SetLineColor(kBlack); gr_eff_PU200_nonprompt_gen_EB->SetLineColor(kBlack);
  gr_eff_PU200_nonprompt_40_EB->SetLineColor(kRed); gr_eff_PU200_nonprompt_60_EB->SetLineColor(kGreen); gr_eff_PU200_nonprompt_80_EB->SetLineColor(kBlue); gr_eff_PU200_nonprompt_100_EB->SetLineColor(kMagenta);
  gr_eff_PU200_nonprompt_2sigma_EB->SetLineColor(kRed); gr_eff_PU200_nonprompt_3sigma_EB->SetLineColor(kGreen); gr_eff_PU200_nonprompt_4sigma_EB->SetLineColor(kBlue);
  gr_eff_noPU_nonprompt_EB->SetLineColor(kGray+1); gr_eff_noPU_nonprompt_gen_EB->SetLineColor(kGray+1);
  gr_eff_noPU_nonprompt_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_gen_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_40_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_60_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_80_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_100_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_2sigma_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_3sigma_EB->SetLineWidth(2); gr_eff_PU200_nonprompt_4sigma_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_EB->SetLineWidth(2); gr_eff_noPU_nonprompt_gen_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_EB->SetLineStyle(7); gr_eff_noPU_nonprompt_gen_EB->SetLineStyle(7);
  gr_eff_noPU_nonprompt_2sigma_EB->SetLineColor(kBlue); gr_eff_noPU_nonprompt_2sigma_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_2sigma_EB->SetLineWidth(2); gr_eff_noPU_nonprompt_gen_2sigma_EB->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_2sigma_EB->SetLineStyle(7); gr_eff_noPU_nonprompt_gen_2sigma_EB->SetLineStyle(7);
  gr_eff_PU200_nonprompt_gen_2sigma_EB->SetLineColor(kBlack); gr_eff_noPU_nonprompt_gen_2sigma_EB->SetLineColor(kGray+1);
  gr_eff_PU200_nonprompt_2sigma_EB_vtx->SetLineColor(kGreen+2); gr_eff_PU200_nonprompt_2sigma_EB_vtx->SetLineWidth(2);

  gr_eff_PU200_nonprompt_EB->SetMarkerStyle(8); gr_eff_PU200_nonprompt_EB->SetMarkerSize(1.5);
  gr_eff_noPU_nonprompt_EB->SetMarkerStyle(8); gr_eff_noPU_nonprompt_EB->SetMarkerSize(1.5); gr_eff_noPU_nonprompt_EB->SetMarkerColor(kGray+1);

    // Endcap region
      // efficiency
  gr_eff_noPU_nonprompt_gen_2sigma_EE->SetTitle("Non-prompt muon in endcap region"); gr_eff_noPU_nonprompt_gen_2sigma_EE->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_nonprompt_gen_2sigma_EE->GetYaxis()->SetTitle("% Counts");
  gr_eff_PU200_nonprompt_EE->SetTitle("Non-prompt muon in endcap region"); gr_eff_noPU_nonprompt_EE->SetTitle("Non-prompt muon in endcap region");
  gr_eff_PU200_nonprompt_EE->GetXaxis()->SetTitle("Relative isolation"); gr_eff_noPU_nonprompt_EE->GetXaxis()->SetTitle("Relative isolation");
  gr_eff_PU200_nonprompt_EE->GetYaxis()->SetTitle("% Counts"); gr_eff_noPU_nonprompt_EE->GetYaxis()->SetTitle("% Counts");
  gr_eff_PU200_nonprompt_EE->SetLineColor(kBlack); gr_eff_PU200_nonprompt_gen_EE->SetLineColor(kBlack);
  gr_eff_PU200_nonprompt_40_EE->SetLineColor(kRed); gr_eff_PU200_nonprompt_60_EE->SetLineColor(kGreen); gr_eff_PU200_nonprompt_80_EE->SetLineColor(kBlue); gr_eff_PU200_nonprompt_100_EE->SetLineColor(kMagenta);
  gr_eff_PU200_nonprompt_2sigma_EE->SetLineColor(kRed); gr_eff_PU200_nonprompt_3sigma_EE->SetLineColor(kGreen); gr_eff_PU200_nonprompt_4sigma_EE->SetLineColor(kBlue);
  gr_eff_noPU_nonprompt_EE->SetLineColor(kGray+1); gr_eff_noPU_nonprompt_gen_EE->SetLineColor(kGray+1);
  gr_eff_noPU_nonprompt_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_gen_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_40_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_60_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_80_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_100_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_2sigma_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_3sigma_EE->SetLineWidth(2); gr_eff_PU200_nonprompt_4sigma_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_EE->SetLineWidth(2); gr_eff_noPU_nonprompt_gen_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_EE->SetLineStyle(7); gr_eff_noPU_nonprompt_gen_EE->SetLineStyle(7);
  gr_eff_PU200_nonprompt_gen_2sigma_EE->SetLineWidth(2); gr_eff_noPU_nonprompt_gen_2sigma_EE->SetLineWidth(2);
  gr_eff_PU200_nonprompt_gen_2sigma_EE->SetLineStyle(7); gr_eff_noPU_nonprompt_gen_2sigma_EE->SetLineStyle(7);
  gr_eff_PU200_nonprompt_gen_2sigma_EE->SetLineColor(kBlack); gr_eff_noPU_nonprompt_gen_2sigma_EE->SetLineColor(kGray+1);
  gr_eff_PU200_nonprompt_2sigma_EE_vtx->SetLineColor(kGreen+2); gr_eff_PU200_nonprompt_2sigma_EE_vtx->SetLineWidth(2);

  gr_eff_PU200_nonprompt_EE->SetMarkerStyle(8); gr_eff_PU200_nonprompt_EE->SetMarkerSize(1.5);
  gr_eff_noPU_nonprompt_EE->SetMarkerStyle(8); gr_eff_noPU_nonprompt_EE->SetMarkerSize(1.5); gr_eff_noPU_nonprompt_EE->SetMarkerColor(kGray+1);

  /////////////
  // Legends //
  /////////////
  // Prompt
    // Barrel region
      // efficiency
  TLegend* leg_eff_prompt_EB_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_gen_EB, "gen PU200");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_noPU_prompt_gen_EB, "gen noPU");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_40_EB, "40ps PU200");
//  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_60_EB, "60ps PU200");
//  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_80_EB, "80ps PU200");
//  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_100_EB, "100ps PU200");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_noPU_prompt_EB, "no MTD noPU");
  leg_eff_prompt_EB_dt->AddEntry(gr_eff_PU200_prompt_EB, "no MTD PU200");
  leg_eff_prompt_EB_dt->SetTextSize(0.03);
  TLegend* leg_eff_prompt_EB_sigma = new TLegend(0.62, 0.63, 0.88, 0.83);
  leg_eff_prompt_EB_sigma->SetMargin(0.2);
//  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_gen_EB, "gen PU200");
//  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_noPU_prompt_gen_EB, "gen noPU");
//  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_2sigma_EB_vtx, "2sigma PU200 (PV, track)");
//  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_2sigma_EB, "2sigma PU200 (muon, track)");
  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_PU200_prompt_EB, "no MTD PU200");
  leg_eff_prompt_EB_sigma->AddEntry(gr_eff_noPU_prompt_EB, "no MTD noPU");
  leg_eff_prompt_EB_sigma->SetTextSize(0.03);
    // Endcap region
      // efficiency
  TLegend* leg_eff_prompt_EE_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_gen_EE, "gen PU200");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_noPU_prompt_gen_EE, "gen noPU");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_40_EE, "40ps PU200");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_PU200_prompt_EE, "no MTD PU200");
  leg_eff_prompt_EE_dt->AddEntry(gr_eff_noPU_prompt_EE, "no MTD noPU");
  leg_eff_prompt_EE_dt->SetTextSize(0.03);
  TLegend* leg_eff_prompt_EE_sigma = new TLegend(0.62, 0.63, 0.88, 0.83);
  leg_eff_prompt_EE_sigma->SetMargin(0.2);
//  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_gen_EE, "gen PU200");
//  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_noPU_prompt_gen_EE, "gen noPU");
//  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_2sigma_EE_vtx, "2sigma PU200 (PV, track)");
//  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_2sigma_EE, "2sigma PU200 (muon, track)");
  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_PU200_prompt_EE, "no MTD PU200");
  leg_eff_prompt_EE_sigma->AddEntry(gr_eff_noPU_prompt_EE, "no MTD noPU");
  leg_eff_prompt_EE_sigma->SetTextSize(0.03);

  // Nonprompt
    // Barrel region
      // efficiency
  TLegend* leg_eff_nonprompt_EB_dt = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_gen_EB, "gen PU200");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_noPU_nonprompt_gen_EB, "gen noPU");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_40_EB, "40ps PU200");
//  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_60_EB, "60ps PU200");
//  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_80_EB, "80ps PU200");
//  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_100_EB, "100ps PU200");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_PU200_nonprompt_EB, "no MTD PU200");
  leg_eff_nonprompt_EB_dt->AddEntry(gr_eff_noPU_nonprompt_EB, "no MTD noPU");
  leg_eff_nonprompt_EB_dt->SetTextSize(0.03);
  TLegend* leg_eff_nonprompt_EB_sigma = new TLegend(0.62, 0.63, 0.88, 0.83);
  leg_eff_nonprompt_EB_sigma->SetMargin(0.2);
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_gen_EB, "gen PU200");
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_noPU_nonprompt_gen_EB, "gen noPU");
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_2sigma_EB_vtx, "2sigma PU200 (PV, track)");
//  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_2sigma_EB, "2sigma PU200 (muon, track)");
  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_PU200_nonprompt_EB, "no MTD PU200");
  leg_eff_nonprompt_EB_sigma->AddEntry(gr_eff_noPU_nonprompt_EB, "no MTD noPU");
  leg_eff_nonprompt_EB_sigma->SetTextSize(0.03);
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
  leg_eff_nonprompt_EE_dt->AddEntry(gr_eff_noPU_nonprompt_EE, "no MTD noPU");
  leg_eff_nonprompt_EE_dt->SetTextSize(0.03);
  TLegend* leg_eff_nonprompt_EE_sigma = new TLegend(0.62, 0.63, 0.88, 0.83);
  leg_eff_nonprompt_EE_sigma->SetMargin(0.2);
//  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_gen_EE, "gen PU200");
//  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_noPU_nonprompt_gen_EE, "gen noPU");
//  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_2sigma_EE_vtx, "2sigma PU200 (PV, track)");
//  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_2sigma_EE, "2sigma PU200 (muon, track)");
  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_PU200_nonprompt_EE, "no MTD PU200");
  leg_eff_nonprompt_EE_sigma->AddEntry(gr_eff_noPU_nonprompt_EE, "no MTD noPU");
  leg_eff_nonprompt_EE_sigma->SetTextSize(0.03);


  ////////////////
  // Draw plots //
  ////////////////
  // Prompt
    // Barrel region
      // efficiency
/*
  TCanvas* c_iso_distribution_prompt_dt_EB = new TCanvas("c_iso_distribution_prompt_dt_EB", "c_iso_distribution_prompt_dt_EB", 1500, 1500);
  c_iso_distribution_prompt_dt_EB->cd();
  c_iso_distribution_prompt_dt_EB->SetGrid();
  c_iso_distribution_prompt_dt_EB->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EB->SetTitle("");
  gr_eff_PU200_prompt_EB->Draw("AL");
  gr_eff_PU200_prompt_40_EB->Draw("same");
//  gr_eff_PU200_prompt_60_EB->Draw("same");
//  gr_eff_PU200_prompt_80_EB->Draw("same");
//  gr_eff_PU200_prompt_100_EB->Draw("same");
  gr_eff_noPU_prompt_EB->Draw("same");
  gr_eff_PU200_prompt_gen_EB->Draw("same");
  gr_eff_noPU_prompt_gen_EB->Draw("same");
  leg_eff_prompt_EB_dt->Draw();
  c_iso_distribution_prompt_dt_EB->Print("plots/eff_prompt_dt_EB.pdf");
*/

  TCanvas* c_iso_distribution_prompt_sigma_EB = new TCanvas("c_iso_distribution_prompt_sigma_EB", "c_iso_distribution_prompt_sigma_EB", 1500, 1500);
  c_iso_distribution_prompt_sigma_EB->cd();
  c_iso_distribution_prompt_sigma_EB->SetGrid();
  c_iso_distribution_prompt_sigma_EB->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EB->GetXaxis()->SetRangeUser(0., 0.08);
  gr_eff_PU200_prompt_EB->GetYaxis()->SetRangeUser(0., 0.85);
  gr_eff_PU200_prompt_EB->SetTitle("");
  gr_eff_PU200_prompt_EB->Draw("ALP");
//  gr_eff_PU200_prompt_2sigma_EB->Draw("same");
//  gr_eff_PU200_prompt_2sigma_EB_vtx->Draw("same");
  gr_eff_noPU_prompt_EB->Draw("same LP");
//  gr_eff_PU200_prompt_gen_EB->Draw("same");
//  gr_eff_noPU_prompt_gen_EB->Draw("same");
  leg_eff_prompt_EB_sigma->Draw();
  c_iso_distribution_prompt_sigma_EB->Print("plots/iso_distrib_prompt_sigma_EB.pdf");

    // Endcap region
      // efficiency
/*
  TCanvas* c_iso_distribution_prompt_dt_EE = new TCanvas("c_iso_distribution_prompt_dt_EE", "c_iso_distribution_prompt_dt_EE", 1500, 1500);
  c_iso_distribution_prompt_dt_EE->cd();
  c_iso_distribution_prompt_dt_EE->SetGrid();
  c_iso_distribution_prompt_dt_EE->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EE->GetXaxis()->SetRangeUser(0., 0.28);
  gr_eff_PU200_prompt_EE->Draw("AL");
  gr_eff_PU200_prompt_40_EE->Draw("same");
//  gr_eff_PU200_prompt_60_EE->Draw("same");
//  gr_eff_PU200_prompt_80_EE->Draw("same");
//  gr_eff_PU200_prompt_100_EE->Draw("same");
  gr_eff_noPU_prompt_EE->Draw("same");
  gr_eff_PU200_prompt_gen_EE->Draw("same");
  gr_eff_noPU_prompt_gen_EE->Draw("same");
  leg_eff_prompt_EE_dt->Draw();
  c_iso_distribution_prompt_dt_EE->Print("plots/eff_prompt_dt_EE.pdf");
*/

  TCanvas* c_iso_distribution_prompt_sigma_EE = new TCanvas("c_iso_distribution_prompt_sigma_EE", "c_iso_distribution_prompt_sigma_EE", 1500, 1500);
  c_iso_distribution_prompt_sigma_EE->cd();
  c_iso_distribution_prompt_sigma_EE->SetGrid();
  c_iso_distribution_prompt_sigma_EE->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EE->GetXaxis()->SetRangeUser(0., 0.08);
  gr_eff_PU200_prompt_EE->GetYaxis()->SetRangeUser(0., 0.75);
  gr_eff_PU200_prompt_EE->SetTitle("");
  gr_eff_PU200_prompt_EE->Draw("ALP");
//  gr_eff_PU200_prompt_2sigma_EE->Draw("same");
//  gr_eff_PU200_prompt_2sigma_EE_vtx->Draw("same");
  gr_eff_noPU_prompt_EE->Draw("same LP");
//  gr_eff_PU200_prompt_gen_EE->Draw("same");
//  gr_eff_noPU_prompt_gen_EE->Draw("same");
  leg_eff_prompt_EE_sigma->Draw();
  c_iso_distribution_prompt_sigma_EE->Print("plots/iso_distrib_prompt_sigma_EE.pdf");

  // Nonprompt
    // Barrel region
      // efficiency
/*
  TCanvas* c_iso_distribution_nonprompt_dt_EB = new TCanvas("c_iso_distribution_nonprompt_dt_EB", "c_iso_distribution_nonprompt_dt_EB", 1500, 1500);
  c_iso_distribution_nonprompt_dt_EB->cd();
  c_iso_distribution_nonprompt_dt_EB->SetGrid();
  c_iso_distribution_nonprompt_dt_EB->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EB->Draw("AL");
  gr_eff_PU200_nonprompt_40_EB->Draw("same");
//  gr_eff_PU200_nonprompt_60_EB->Draw("same");
//  gr_eff_PU200_nonprompt_80_EB->Draw("same");
//  gr_eff_PU200_nonprompt_100_EB->Draw("same");
  gr_eff_noPU_nonprompt_EB->Draw("same");
  gr_eff_PU200_nonprompt_gen_EB->Draw("same");
  gr_eff_noPU_nonprompt_gen_EB->Draw("same");
  leg_eff_nonprompt_EB_dt->Draw();
  c_iso_distribution_nonprompt_dt_EB->Print("plots/eff_nonprompt_dt_EB.pdf");
*/

  TCanvas* c_iso_distribution_nonprompt_sigma_EB = new TCanvas("c_iso_distribution_nonprompt_sigma_EB", "c_iso_distribution_nonprompt_sigma_EB", 1500, 1500);
  c_iso_distribution_nonprompt_sigma_EB->cd();
  c_iso_distribution_nonprompt_sigma_EB->SetGrid();
  c_iso_distribution_nonprompt_sigma_EB->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EB->GetXaxis()->SetRangeUser(0., 0.08);
  gr_eff_PU200_nonprompt_EB->GetYaxis()->SetRangeUser(0., 0.15);
  gr_eff_PU200_nonprompt_EB->SetTitle("");
  gr_eff_PU200_nonprompt_EB->Draw("ALP");
//  gr_eff_PU200_nonprompt_2sigma_EB->Draw("same");
//  gr_eff_PU200_nonprompt_2sigma_EB_vtx->Draw("same");
  gr_eff_noPU_nonprompt_EB->Draw("same LP");
//  gr_eff_PU200_nonprompt_gen_EB->Draw("same");
//  gr_eff_noPU_nonprompt_gen_EB->Draw("same");
  leg_eff_nonprompt_EB_sigma->Draw();
  c_iso_distribution_nonprompt_sigma_EB->Print("plots/iso_distrib_nonprompt_sigma_EB.pdf");

    // Endcap region
      // efficiency
/*
  TCanvas* c_iso_distribution_nonprompt_dt_EE = new TCanvas("c_iso_distribution_nonprompt_dt_EE", "c_iso_distribution_nonprompt_dt_EE", 1500, 1500);
  c_iso_distribution_nonprompt_dt_EE->cd();
  c_iso_distribution_nonprompt_dt_EE->SetGrid();
  c_iso_distribution_nonprompt_dt_EE->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EE->SetTitle("");
  gr_eff_PU200_nonprompt_EE->Draw("AL");
  gr_eff_PU200_nonprompt_40_EE->Draw("same");
//  gr_eff_PU200_nonprompt_60_EE->Draw("same");
//  gr_eff_PU200_nonprompt_80_EE->Draw("same");
//  gr_eff_PU200_nonprompt_100_EE->Draw("same");
  gr_eff_noPU_nonprompt_EE->Draw("same");
  gr_eff_PU200_nonprompt_gen_EE->Draw("same");
  gr_eff_noPU_nonprompt_gen_EE->Draw("same");
  leg_eff_nonprompt_EE_dt->Draw();
  c_iso_distribution_nonprompt_dt_EE->Print("plots/eff_nonprompt_dt_EE.pdf");
*/

  TCanvas* c_iso_distribution_nonprompt_sigma_EE = new TCanvas("c_iso_distribution_nonprompt_sigma_EE", "c_iso_distribution_nonprompt_sigma_EE", 1500, 1500);
  c_iso_distribution_nonprompt_sigma_EE->cd();
  c_iso_distribution_nonprompt_sigma_EE->SetGrid();
  c_iso_distribution_nonprompt_sigma_EE->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EE->GetXaxis()->SetRangeUser(0., 0.08);
  gr_eff_PU200_nonprompt_EE->GetYaxis()->SetRangeUser(0., 0.15);
  gr_eff_PU200_nonprompt_EE->SetTitle("");
  gr_eff_PU200_nonprompt_EE->Draw("ALP");
//  gr_eff_PU200_nonprompt_2sigma_EE->Draw("same");
//  gr_eff_PU200_nonprompt_2sigma_EE_vtx->Draw("same");
  gr_eff_noPU_nonprompt_EE->Draw("same LP");
//  gr_eff_PU200_nonprompt_gen_EE->Draw("same");
//  gr_eff_noPU_nonprompt_gen_EE->Draw("same");
  leg_eff_nonprompt_EE_sigma->Draw();
  c_iso_distribution_nonprompt_sigma_EE->Print("plots/iso_distrib_nonprompt_sigma_EE.pdf");

/*
  // test
  TCanvas* c_iso_distribution_prompt_sigma_EB_test = new TCanvas("c_iso_distribution_prompt_sigma_EB_test", "c_iso_distribution_prompt_sigma_EB_test", 1500, 1500);
  c_iso_distribution_prompt_sigma_EB_test->cd();
  c_iso_distribution_prompt_sigma_EB_test->SetGrid();
  c_iso_distribution_prompt_sigma_EB_test->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EB->Draw("AL");
  gr_eff_PU200_prompt_2sigma_EB->Draw("same");
  gr_eff_noPU_prompt_EB->Draw("same");
  gr_eff_PU200_prompt_gen_2sigma_EB->Draw("same");
  gr_eff_noPU_prompt_gen_2sigma_EB->Draw("same");
  leg_eff_prompt_EB_sigma_test->Draw();
  c_iso_distribution_prompt_sigma_EB_test->Print("plots/eff_prompt_sigma_EB_gen_test.pdf");

  TCanvas* c_iso_distribution_prompt_sigma_EE_test = new TCanvas("c_iso_distribution_prompt_sigma_EE_test", "c_iso_distribution_prompt_sigma_EE_test", 1500, 1500);
  c_iso_distribution_prompt_sigma_EE_test->cd();
  c_iso_distribution_prompt_sigma_EE_test->SetGrid();
  c_iso_distribution_prompt_sigma_EE_test->SetLeftMargin(0.12);
  gr_eff_PU200_prompt_EE->Draw("AL");
  gr_eff_PU200_prompt_2sigma_EE->Draw("same");
  gr_eff_noPU_prompt_EE->Draw("same");
  gr_eff_PU200_prompt_gen_2sigma_EE->Draw("same");
  gr_eff_noPU_prompt_gen_2sigma_EE->Draw("same");
  leg_eff_prompt_EE_sigma_test->Draw();
  c_iso_distribution_prompt_sigma_EE_test->Print("plots/eff_prompt_sigma_EE_gen_test.pdf");

  TCanvas* c_iso_distribution_nonprompt_sigma_EB_test = new TCanvas("c_iso_distribution_nonprompt_sigma_EB_test", "c_iso_distribution_nonprompt_sigma_EB_test", 1500, 1500);
  c_iso_distribution_nonprompt_sigma_EB_test->cd();
  c_iso_distribution_nonprompt_sigma_EB_test->SetGrid();
  c_iso_distribution_nonprompt_sigma_EB_test->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EB->GetYaxis()->SetRangeUser(0.,0.23);
  gr_eff_PU200_nonprompt_EB->Draw("AL");
  gr_eff_PU200_nonprompt_2sigma_EB->Draw("same");
  gr_eff_noPU_nonprompt_EB->Draw("same");
  gr_eff_PU200_nonprompt_gen_2sigma_EB->Draw("same");
  gr_eff_noPU_nonprompt_gen_2sigma_EB->Draw("same");
  leg_eff_nonprompt_EB_sigma_test->Draw();
  c_iso_distribution_nonprompt_sigma_EB_test->Print("plots/eff_nonprompt_sigma_EB_gen_test.pdf");

  TCanvas* c_iso_distribution_nonprompt_sigma_EE_test = new TCanvas("c_iso_distribution_nonprompt_sigma_EE_test", "c_iso_distribution_nonprompt_sigma_EE_test", 1500, 1500);
  c_iso_distribution_nonprompt_sigma_EE_test->cd();
  c_iso_distribution_nonprompt_sigma_EE_test->SetGrid();
  c_iso_distribution_nonprompt_sigma_EE_test->SetLeftMargin(0.12);
  gr_eff_PU200_nonprompt_EE->GetYaxis()->SetRangeUser(0.,0.33);
  gr_eff_PU200_nonprompt_EE->Draw("AL");
  gr_eff_PU200_nonprompt_2sigma_EE->Draw("same");
  gr_eff_noPU_nonprompt_EE->Draw("same");
  gr_eff_PU200_nonprompt_gen_2sigma_EE->Draw("same");
  gr_eff_noPU_nonprompt_gen_2sigma_EE->Draw("same");
  leg_eff_nonprompt_EE_sigma_test->Draw();
  c_iso_distribution_nonprompt_sigma_EE_test->Print("plots/eff_nonprompt_sigma_EE_gen_test.pdf");

  // Martina
  TCanvas* c_iso_distribution_prompt_sigma_EB_region = new TCanvas("c_iso_distribution_prompt_sigma_EB_region", "c_iso_distribution_prompt_sigma_EB_region", 1500, 1500);
  c_iso_distribution_prompt_sigma_EB_region->cd();
  c_iso_distribution_prompt_sigma_EB_region->SetGrid();
  c_iso_distribution_prompt_sigma_EB_region->SetLeftMargin(0.12);

  //gr_eff_noPU_prompt_EB->GetXaxis()->SetRangeUser(0., 0.07);
  
  //gr_eff_noPU_prompt_EB->GetXaxis()->SetRangeUser(0.07, 0.14);
  //gr_eff_noPU_prompt_EB->GetYaxis()->SetRangeUser(0.94, 1.);
  
  //gr_eff_noPU_prompt_EB->GetXaxis()->SetRangeUser(0.14, 0.21);
  //gr_eff_noPU_prompt_EB->GetYaxis()->SetRangeUser(0.98, 1.);

  gr_eff_noPU_prompt_EB->GetXaxis()->SetRangeUser(0.21, 0.28);
  gr_eff_noPU_prompt_EB->GetYaxis()->SetRangeUser(0.994, 1.);

  gr_eff_noPU_prompt_EB->Draw("AL");
  gr_eff_noPU_prompt_gen_EB->Draw("same");
  leg_eff_prompt_EB_sigma_region->Draw();
  c_iso_distribution_prompt_sigma_EB_region->Print("plots/eff_prompt_sigma_EB_region.pdf");

  TCanvas* c_iso_distribution_prompt_sigma_EE_region = new TCanvas("c_iso_distribution_prompt_sigma_EE_region", "c_iso_distribution_prompt_sigma_EE_region", 1500, 1500);
  c_iso_distribution_prompt_sigma_EE_region->cd();
  c_iso_distribution_prompt_sigma_EE_region->SetGrid();
  c_iso_distribution_prompt_sigma_EE_region->SetLeftMargin(0.12);

  //gr_eff_noPU_prompt_EE->GetXaxis()->SetRangeUser(0., 0.07);
  
  //gr_eff_noPU_prompt_EE->GetXaxis()->SetRangeUser(0.07, 0.14);
  //gr_eff_noPU_prompt_EE->GetYaxis()->SetRangeUser(0.94, 1.);
  
  //gr_eff_noPU_prompt_EE->GetXaxis()->SetRangeUser(0.14, 0.21);
  //gr_eff_noPU_prompt_EE->GetYaxis()->SetRangeUser(0.98, 1.);

  gr_eff_noPU_prompt_EE->GetXaxis()->SetRangeUser(0.21, 0.28);
  gr_eff_noPU_prompt_EE->GetYaxis()->SetRangeUser(0.993, 1.);

  gr_eff_noPU_prompt_EE->Draw("AL");
  gr_eff_noPU_prompt_gen_EE->Draw("same");
  leg_eff_prompt_EE_sigma_region->Draw();
  c_iso_distribution_prompt_sigma_EE_region->Print("plots/eff_prompt_sigma_EE_region.pdf");

  TCanvas* c_iso_distribution_nonprompt_sigma_EB_region = new TCanvas("c_iso_distribution_nonprompt_sigma_EB_region", "c_iso_distribution_nonprompt_sigma_EB_region", 1500, 1500);
  c_iso_distribution_nonprompt_sigma_EB_region->cd();
  c_iso_distribution_nonprompt_sigma_EB_region->SetGrid();
  c_iso_distribution_nonprompt_sigma_EB_region->SetLeftMargin(0.12);

  //gr_eff_noPU_nonprompt_EB->GetXaxis()->SetRangeUser(0., 1.);
  
  //gr_eff_noPU_nonprompt_EB->GetXaxis()->SetRangeUser(1., 2.);
  //gr_eff_noPU_nonprompt_EB->GetYaxis()->SetRangeUser(0.4, 0.8);
  
  //gr_eff_noPU_nonprompt_EB->GetXaxis()->SetRangeUser(2., 3.);
  //gr_eff_noPU_nonprompt_EB->GetYaxis()->SetRangeUser(0.72, 0.9);

  gr_eff_noPU_nonprompt_EB->GetXaxis()->SetRangeUser(3., 4.);
  gr_eff_noPU_nonprompt_EB->GetYaxis()->SetRangeUser(0.85, 1.);

  gr_eff_noPU_nonprompt_EB->Draw("AL");
  gr_eff_noPU_nonprompt_gen_EB->Draw("same");
  leg_eff_nonprompt_EB_sigma_region->Draw();
  c_iso_distribution_nonprompt_sigma_EB_region->Print("plots/eff_nonprompt_sigma_EB_region.pdf");

  TCanvas* c_iso_distribution_nonprompt_sigma_EE_region = new TCanvas("c_iso_distribution_nonprompt_sigma_EE_region", "c_iso_distribution_nonprompt_sigma_EE_region", 1500, 1500);
  c_iso_distribution_nonprompt_sigma_EE_region->cd();
  c_iso_distribution_nonprompt_sigma_EE_region->SetGrid();
  c_iso_distribution_nonprompt_sigma_EE_region->SetLeftMargin(0.12);

  //gr_eff_noPU_nonprompt_EE->GetXaxis()->SetRangeUser(0., 1.);
  
  //gr_eff_noPU_nonprompt_EE->GetXaxis()->SetRangeUser(1., 2.);
  //gr_eff_noPU_nonprompt_EE->GetYaxis()->SetRangeUser(0.58, 0.85);

  //gr_eff_noPU_nonprompt_EE->GetXaxis()->SetRangeUser(2., 3.);
  //gr_eff_noPU_nonprompt_EE->GetYaxis()->SetRangeUser(0.8, 0.94);

  gr_eff_noPU_nonprompt_EE->GetXaxis()->SetRangeUser(3., 4.);
  gr_eff_noPU_nonprompt_EE->GetYaxis()->SetRangeUser(0.9, 1.);

  gr_eff_noPU_nonprompt_EE->Draw("AL");
  gr_eff_noPU_nonprompt_gen_EE->Draw("same");
  leg_eff_nonprompt_EE_sigma_region->Draw();
  c_iso_distribution_nonprompt_sigma_EE_region->Print("plots/eff_nonprompt_sigma_EE_region.pdf");

  // test end
*/

}



void draw_reliso_roc() {
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  TFile *f_PU200_prompt_vtx, *f_PU200_nonprompt_vtx, *f_noPU_prompt_vtx, *f_noPU_nonprompt_vtx;
  // PU200
  f_PU200_prompt = new TFile(Form("data/%s", path_PU200_prompt.Data()));
  f_PU200_nonprompt = new TFile(Form("data/%s", path_PU200_nonprompt.Data()));
  f_PU200_prompt_vtx = new TFile(Form("data/%s", path_PU200_prompt_vtx.Data()));
  f_PU200_nonprompt_vtx = new TFile(Form("data/%s", path_PU200_nonprompt_vtx.Data()));
  // noPU
  f_noPU_prompt = new TFile(Form("data/%s", path_noPU_prompt.Data()));
  f_noPU_nonprompt = new TFile(Form("data/%s", path_noPU_nonprompt.Data()));
  f_noPU_prompt_vtx = new TFile(Form("data/%s", path_noPU_prompt_vtx.Data()));
  f_noPU_nonprompt_vtx = new TFile(Form("data/%s", path_noPU_nonprompt_vtx.Data()));

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
    // vtx
  TH1D *h_PU200_prompt_2sigma_EB_vtx, *h_PU200_nonprompt_2sigma_EB_vtx, *h_noPU_prompt_2sigma_EB_vtx, *h_noPU_2sigma_nonprompt_EB_vtx;
  TH1D *h_PU200_prompt_3sigma_EB_vtx, *h_PU200_nonprompt_3sigma_EB_vtx, *h_noPU_prompt_3sigma_EB_vtx, *h_noPU_3sigma_nonprompt_EB_vtx;
  TH1D *h_PU200_prompt_4sigma_EB_vtx, *h_PU200_nonprompt_4sigma_EB_vtx, *h_noPU_prompt_4sigma_EB_vtx, *h_noPU_4sigma_nonprompt_EB_vtx;

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
    // vtx
  TH1D *h_PU200_prompt_2sigma_EE_vtx, *h_PU200_nonprompt_2sigma_EE_vtx, *h_noPU_prompt_2sigma_EE_vtx, *h_noPU_2sigma_nonprompt_EE_vtx;
  TH1D *h_PU200_prompt_3sigma_EE_vtx, *h_PU200_nonprompt_3sigma_EE_vtx, *h_noPU_prompt_3sigma_EE_vtx, *h_noPU_3sigma_nonprompt_EE_vtx;
  TH1D *h_PU200_prompt_4sigma_EE_vtx, *h_PU200_nonprompt_4sigma_EE_vtx, *h_noPU_prompt_4sigma_EE_vtx, *h_noPU_4sigma_nonprompt_EE_vtx;


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
      // vtx
  h_PU200_prompt_4sigma_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EB");
  h_PU200_prompt_3sigma_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EB");
  h_PU200_prompt_2sigma_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EB");
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
      // vtx
  h_PU200_prompt_4sigma_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Sig_EE");
  h_PU200_prompt_3sigma_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Sig_EE");
  h_PU200_prompt_2sigma_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Sig_EE");
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
      // vtx
  h_PU200_nonprompt_4sigma_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EB");
  h_PU200_nonprompt_3sigma_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EB");
  h_PU200_nonprompt_2sigma_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EB");
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
      // vtx
  h_PU200_nonprompt_4sigma_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_4sigma_Bkg_EE");
  h_PU200_nonprompt_3sigma_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_3sigma_Bkg_EE");
  h_PU200_nonprompt_2sigma_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_rel_chIso_sum_MTD_2sigma_Bkg_EE");

  int nbin = h_PU200_prompt_EE->GetNbinsX()+1;

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
        // vtx
  vector<double> prompt_eff_PU200_2sigma_EB_vtx={0}, prompt_norm_PU200_2sigma_EB_vtx={0}, prompt_eff_PU200_2sigma_EE_vtx={0}, prompt_norm_PU200_2sigma_EE_vtx={0};
  vector<double> prompt_eff_PU200_3sigma_EB_vtx={0}, prompt_norm_PU200_3sigma_EB_vtx={0}, prompt_eff_PU200_3sigma_EE_vtx={0}, prompt_norm_PU200_3sigma_EE_vtx={0};
  vector<double> prompt_eff_PU200_4sigma_EB_vtx={0}, prompt_norm_PU200_4sigma_EB_vtx={0}, prompt_eff_PU200_4sigma_EE_vtx={0}, prompt_norm_PU200_4sigma_EE_vtx={0};
      // Nonprompt
  vector<double> nonprompt_eff_PU200_EB={0}, nonprompt_norm_PU200_EB={0}, nonprompt_eff_PU200_EE={0}, nonprompt_norm_PU200_EE={0};
  vector<double> nonprompt_eff_PU200_2sigma_EB={0}, nonprompt_norm_PU200_2sigma_EB={0}, nonprompt_eff_PU200_2sigma_EE={0}, nonprompt_norm_PU200_2sigma_EE={0};
  vector<double> nonprompt_eff_PU200_3sigma_EB={0}, nonprompt_norm_PU200_3sigma_EB={0}, nonprompt_eff_PU200_3sigma_EE={0}, nonprompt_norm_PU200_3sigma_EE={0};
  vector<double> nonprompt_eff_PU200_4sigma_EB={0}, nonprompt_norm_PU200_4sigma_EB={0}, nonprompt_eff_PU200_4sigma_EE={0}, nonprompt_norm_PU200_4sigma_EE={0};
  vector<double> nonprompt_eff_PU200_40_EB={0}, nonprompt_norm_PU200_40_EB={0}, nonprompt_eff_PU200_40_EE={0}, nonprompt_norm_PU200_40_EE={0};
  vector<double> nonprompt_eff_PU200_60_EB={0}, nonprompt_norm_PU200_60_EB={0}, nonprompt_eff_PU200_60_EE={0}, nonprompt_norm_PU200_60_EE={0};
  vector<double> nonprompt_eff_PU200_80_EB={0}, nonprompt_norm_PU200_80_EB={0}, nonprompt_eff_PU200_80_EE={0}, nonprompt_norm_PU200_80_EE={0};
  vector<double> nonprompt_eff_PU200_100_EB={0}, nonprompt_norm_PU200_100_EB={0}, nonprompt_eff_PU200_100_EE={0}, nonprompt_norm_PU200_100_EE={0};
        // vtx
  vector<double> nonprompt_eff_PU200_2sigma_EB_vtx={0}, nonprompt_norm_PU200_2sigma_EB_vtx={0}, nonprompt_eff_PU200_2sigma_EE_vtx={0}, nonprompt_norm_PU200_2sigma_EE_vtx={0};
  vector<double> nonprompt_eff_PU200_3sigma_EB_vtx={0}, nonprompt_norm_PU200_3sigma_EB_vtx={0}, nonprompt_eff_PU200_3sigma_EE_vtx={0}, nonprompt_norm_PU200_3sigma_EE_vtx={0};
  vector<double> nonprompt_eff_PU200_4sigma_EB_vtx={0}, nonprompt_norm_PU200_4sigma_EB_vtx={0}, nonprompt_eff_PU200_4sigma_EE_vtx={0}, nonprompt_norm_PU200_4sigma_EE_vtx={0};
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
  //for(unsigned int i=0; i<h_PU200_prompt_EB->GetNbinsX()+2; i++) {
  for(unsigned int i=0; i<20; i++) {
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
    prompt_eff_PU200_2sigma_EB.emplace_back(h_PU200_prompt_2sigma_EB->Integral(1,i+1)/h_PU200_prompt_2sigma_EB->Integral(1,nbin));
    prompt_eff_PU200_3sigma_EB.emplace_back(h_PU200_prompt_3sigma_EB->Integral(1,i+1)/h_PU200_prompt_3sigma_EB->Integral(1,nbin));
    prompt_eff_PU200_4sigma_EB.emplace_back(h_PU200_prompt_4sigma_EB->Integral(1,i+1)/h_PU200_prompt_4sigma_EB->Integral(1,nbin));
    prompt_eff_noPU_gen_EB.emplace_back(h_noPU_prompt_gen_EB->Integral(1,i+1)/h_noPU_prompt_gen_EB->Integral(1,nbin));
        // vtx
    prompt_eff_PU200_2sigma_EB_vtx.emplace_back(h_PU200_prompt_2sigma_EB_vtx->Integral(1,i+1)/h_PU200_prompt_2sigma_EB_vtx->Integral(1,nbin));
    prompt_eff_PU200_3sigma_EB_vtx.emplace_back(h_PU200_prompt_3sigma_EB_vtx->Integral(1,i+1)/h_PU200_prompt_3sigma_EB_vtx->Integral(1,nbin));
    prompt_eff_PU200_4sigma_EB_vtx.emplace_back(h_PU200_prompt_4sigma_EB_vtx->Integral(1,i+1)/h_PU200_prompt_4sigma_EB_vtx->Integral(1,nbin));
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
        // vtx
    prompt_eff_PU200_2sigma_EE_vtx.emplace_back(h_PU200_prompt_2sigma_EE_vtx->Integral(1,i+1)/h_PU200_prompt_2sigma_EE_vtx->Integral(1,nbin));
    prompt_eff_PU200_3sigma_EE_vtx.emplace_back(h_PU200_prompt_3sigma_EE_vtx->Integral(1,i+1)/h_PU200_prompt_3sigma_EE_vtx->Integral(1,nbin));
    prompt_eff_PU200_4sigma_EE_vtx.emplace_back(h_PU200_prompt_4sigma_EE_vtx->Integral(1,i+1)/h_PU200_prompt_4sigma_EE_vtx->Integral(1,nbin));
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
        // vtx
    nonprompt_eff_PU200_2sigma_EB_vtx.emplace_back(h_PU200_nonprompt_2sigma_EB_vtx->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EB_vtx->Integral(1,nbin));
    nonprompt_eff_PU200_3sigma_EB_vtx.emplace_back(h_PU200_nonprompt_3sigma_EB_vtx->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EB_vtx->Integral(1,nbin));
    nonprompt_eff_PU200_4sigma_EB_vtx.emplace_back(h_PU200_nonprompt_4sigma_EB_vtx->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EB_vtx->Integral(1,nbin));
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
        // vtx
    nonprompt_eff_PU200_2sigma_EE_vtx.emplace_back(h_PU200_nonprompt_2sigma_EE_vtx->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EE_vtx->Integral(1,nbin));
    nonprompt_eff_PU200_3sigma_EE_vtx.emplace_back(h_PU200_nonprompt_3sigma_EE_vtx->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EE_vtx->Integral(1,nbin));
    nonprompt_eff_PU200_4sigma_EE_vtx.emplace_back(h_PU200_nonprompt_4sigma_EE_vtx->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EE_vtx->Integral(1,nbin));
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
  /*
  cout << endl;
  cout << endl;
  cout << endl;
  for(int i=0; i<prompt_eff_PU200_EB.size(); i++) {
	cout << "bin " << i << ": (" << nonprompt_eff_PU200_EB.at(i) << ", " << nonprompt_eff_noPU_EB.at(i) << ")" << endl;
  }
  cout << endl;
  cout << endl;
  cout << endl;
  for(int i=0; i<prompt_eff_PU200_EB.size(); i++) {
	cout << "bin " << i << ": (" << prompt_eff_PU200_EB.at(i) << ", " << nonprompt_eff_PU200_EB.at(i) << ")  /  (" << 
		                            prompt_eff_noPU_EB.at(i)  << ", " << nonprompt_eff_noPU_EB.at(i)  << ")" << endl;
  }
  */

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
  TGraph* gr_noPU_2sigma_EB = new TGraph(h_noPU_prompt_2sigma_EB->GetNbinsX(), &prompt_eff_noPU_2sigma_EB[0], &nonprompt_eff_noPU_2sigma_EB[0]);
  TGraph* gr_noPU_3sigma_EB = new TGraph(h_noPU_prompt_3sigma_EB->GetNbinsX(), &prompt_eff_noPU_3sigma_EB[0], &nonprompt_eff_noPU_3sigma_EB[0]);
  TGraph* gr_noPU_4sigma_EB = new TGraph(h_noPU_prompt_4sigma_EB->GetNbinsX(), &prompt_eff_noPU_4sigma_EB[0], &nonprompt_eff_noPU_4sigma_EB[0]);
  TGraph* gr_PU200_gen_EB = new TGraph(h_PU200_prompt_gen_EB->GetNbinsX(), &prompt_eff_PU200_gen_EB[0], &nonprompt_eff_PU200_gen_EB[0]);
  TGraph* gr_noPU_gen_EB = new TGraph(h_noPU_prompt_gen_EB->GetNbinsX(), &prompt_eff_noPU_gen_EB[0], &nonprompt_eff_noPU_gen_EB[0]);
      // vtx
  TGraph* gr_PU200_2sigma_EB_vtx = new TGraph(h_PU200_prompt_2sigma_EB_vtx->GetNbinsX(), &prompt_eff_PU200_2sigma_EB_vtx[0], &nonprompt_eff_PU200_2sigma_EB_vtx[0]);
  TGraph* gr_PU200_3sigma_EB_vtx = new TGraph(h_PU200_prompt_3sigma_EB_vtx->GetNbinsX(), &prompt_eff_PU200_3sigma_EB_vtx[0], &nonprompt_eff_PU200_3sigma_EB_vtx[0]);
  TGraph* gr_PU200_4sigma_EB_vtx = new TGraph(h_PU200_prompt_4sigma_EB_vtx->GetNbinsX(), &prompt_eff_PU200_4sigma_EB_vtx[0], &nonprompt_eff_PU200_4sigma_EB_vtx[0]);

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
      // vtx
  TGraph* gr_PU200_2sigma_EE_vtx = new TGraph(h_PU200_prompt_2sigma_EE_vtx->GetNbinsX(), &prompt_eff_PU200_2sigma_EE_vtx[0], &nonprompt_eff_PU200_2sigma_EE_vtx[0]);
  TGraph* gr_PU200_3sigma_EE_vtx = new TGraph(h_PU200_prompt_3sigma_EE_vtx->GetNbinsX(), &prompt_eff_PU200_3sigma_EE_vtx[0], &nonprompt_eff_PU200_3sigma_EE_vtx[0]);
  TGraph* gr_PU200_4sigma_EE_vtx = new TGraph(h_PU200_prompt_4sigma_EE_vtx->GetNbinsX(), &prompt_eff_PU200_4sigma_EE_vtx[0], &nonprompt_eff_PU200_4sigma_EE_vtx[0]);

  /*
  // test 
  // Define TGraph
  // Prompt
    // Barrel region
  TGraph* gr_PU200_EB = new TGraph(20, &prompt_eff_PU200_EB[0], &nonprompt_eff_PU200_EB[0]);
  TGraph* gr_PU200_40_EB = new TGraph(20, &prompt_eff_PU200_40_EB[0], &nonprompt_eff_PU200_40_EB[0]);
  TGraph* gr_PU200_60_EB = new TGraph(20, &prompt_eff_PU200_60_EB[0], &nonprompt_eff_PU200_60_EB[0]);
  TGraph* gr_PU200_80_EB = new TGraph(20, &prompt_eff_PU200_80_EB[0], &nonprompt_eff_PU200_80_EB[0]);
  TGraph* gr_PU200_100_EB = new TGraph(20, &prompt_eff_PU200_100_EB[0], &nonprompt_eff_PU200_100_EB[0]);
  TGraph* gr_PU200_2sigma_EB = new TGraph(20, &prompt_eff_PU200_2sigma_EB[0], &nonprompt_eff_PU200_2sigma_EB[0]);
  TGraph* gr_PU200_3sigma_EB = new TGraph(20, &prompt_eff_PU200_3sigma_EB[0], &nonprompt_eff_PU200_3sigma_EB[0]);
  TGraph* gr_PU200_4sigma_EB = new TGraph(20, &prompt_eff_PU200_4sigma_EB[0], &nonprompt_eff_PU200_4sigma_EB[0]);
  TGraph* gr_noPU_EB = new TGraph(20, &prompt_eff_noPU_EB[0], &nonprompt_eff_noPU_EB[0]);
  TGraph* gr_PU200_gen_EB = new TGraph(20, &prompt_eff_PU200_gen_EB[0], &nonprompt_eff_PU200_gen_EB[0]);
  TGraph* gr_noPU_gen_EB = new TGraph(20, &prompt_eff_noPU_gen_EB[0], &nonprompt_eff_noPU_gen_EB[0]);
      // vtx
  TGraph* gr_PU200_2sigma_EB_vtx = new TGraph(20, &prompt_eff_PU200_2sigma_EB_vtx[0], &nonprompt_eff_PU200_2sigma_EB_vtx[0]);
  TGraph* gr_PU200_3sigma_EB_vtx = new TGraph(20, &prompt_eff_PU200_3sigma_EB_vtx[0], &nonprompt_eff_PU200_3sigma_EB_vtx[0]);
  TGraph* gr_PU200_4sigma_EB_vtx = new TGraph(20, &prompt_eff_PU200_4sigma_EB_vtx[0], &nonprompt_eff_PU200_4sigma_EB_vtx[0]);

    // Endcap region
  TGraph* gr_PU200_EE = new TGraph(20, &prompt_eff_PU200_EE[0], &nonprompt_eff_PU200_EE[0]);
  TGraph* gr_PU200_40_EE = new TGraph(20, &prompt_eff_PU200_40_EE[0], &nonprompt_eff_PU200_40_EE[0]);
  TGraph* gr_PU200_60_EE = new TGraph(20, &prompt_eff_PU200_60_EE[0], &nonprompt_eff_PU200_60_EE[0]);
  TGraph* gr_PU200_80_EE = new TGraph(20, &prompt_eff_PU200_80_EE[0], &nonprompt_eff_PU200_80_EE[0]);
  TGraph* gr_PU200_100_EE = new TGraph(20, &prompt_eff_PU200_100_EE[0], &nonprompt_eff_PU200_100_EE[0]);
  TGraph* gr_PU200_2sigma_EE = new TGraph(20, &prompt_eff_PU200_2sigma_EE[0], &nonprompt_eff_PU200_2sigma_EE[0]);
  TGraph* gr_PU200_3sigma_EE = new TGraph(20, &prompt_eff_PU200_3sigma_EE[0], &nonprompt_eff_PU200_3sigma_EE[0]);
  TGraph* gr_PU200_4sigma_EE = new TGraph(20, &prompt_eff_PU200_4sigma_EE[0], &nonprompt_eff_PU200_4sigma_EE[0]);
  TGraph* gr_noPU_EE = new TGraph(20, &prompt_eff_noPU_EE[0], &nonprompt_eff_noPU_EE[0]);
  TGraph* gr_PU200_gen_EE = new TGraph(20, &prompt_eff_PU200_gen_EE[0], &nonprompt_eff_PU200_gen_EE[0]);
  TGraph* gr_noPU_gen_EE = new TGraph(20, &prompt_eff_noPU_gen_EE[0], &nonprompt_eff_noPU_gen_EE[0]);
      // vtx
  TGraph* gr_PU200_2sigma_EE_vtx = new TGraph(20, &prompt_eff_PU200_2sigma_EE_vtx[0], &nonprompt_eff_PU200_2sigma_EE_vtx[0]);
  TGraph* gr_PU200_3sigma_EE_vtx = new TGraph(20, &prompt_eff_PU200_3sigma_EE_vtx[0], &nonprompt_eff_PU200_3sigma_EE_vtx[0]);
  TGraph* gr_PU200_4sigma_EE_vtx = new TGraph(20, &prompt_eff_PU200_4sigma_EE_vtx[0], &nonprompt_eff_PU200_4sigma_EE_vtx[0]);

  // Remove the value at dump point (0,0)
  gr_PU200_EB->RemovePoint(0); gr_PU200_40_EB->RemovePoint(0); gr_PU200_60_EB->RemovePoint(0); gr_PU200_80_EB->RemovePoint(0); gr_PU200_100_EB->RemovePoint(0); gr_PU200_2sigma_EB->RemovePoint(0); gr_PU200_3sigma_EB->RemovePoint(0); gr_PU200_4sigma_EB->RemovePoint(0); gr_noPU_EB->RemovePoint(0); gr_PU200_gen_EB->RemovePoint(0); gr_noPU_gen_EB->RemovePoint(0);
  gr_PU200_EE->RemovePoint(0); gr_PU200_40_EE->RemovePoint(0); gr_PU200_60_EE->RemovePoint(0); gr_PU200_80_EE->RemovePoint(0); gr_PU200_100_EE->RemovePoint(0); gr_PU200_2sigma_EE->RemovePoint(0); gr_PU200_3sigma_EE->RemovePoint(0); gr_PU200_4sigma_EE->RemovePoint(0); gr_noPU_EE->RemovePoint(0);gr_PU200_gen_EE->RemovePoint(0); gr_noPU_gen_EE->RemovePoint(0);
  */

  /*
  ///////
  // test
  TGraph* gr_PU200_EB_test = new TGraph();
  TGraph* gr_PU200_2sigma_EB_test = new TGraph();
  for(int i=0; i<h_PU200_prompt_EB->GetNbinsX(); i++) {
    gr_PU200_EB_test->SetPoint(i, prompt_eff_PU200_EB.at(i), nonprompt_eff_PU200_EB.at(i));
    gr_PU200_2sigma_EB_test->SetPoint(i, prompt_eff_PU200_2sigma_EB.at(i), nonprompt_eff_PU200_2sigma_EB.at(i));
  }
  gr_PU200_EB_test->SetLineWidth(2); gr_PU200_2sigma_EB_test->SetLineWidth(2);
  gr_PU200_EB_test->SetLineColor(kBlack); gr_PU200_2sigma_EB_test->SetLineColor(kRed);
  TLegend* leg_EB_sigma_test = new TLegend(0.15, 0.65, 0.48, 0.85);
  leg_EB_sigma_test->AddEntry(gr_PU200_EB_test, "no MTD PU200");
  leg_EB_sigma_test->AddEntry(gr_PU200_2sigma_EB_test, "2sigma PU200");
  TCanvas* c_sigma_EB_test = new TCanvas("c_sigma_EB_test", "c_sigma_EB_test", 1500, 1500);
  c_sigma_EB_test->cd();
  c_sigma_EB_test->SetGrid();
  gr_PU200_EB_test->GetXaxis()->SetRangeUser(0.75,1.01);
  gr_PU200_EB_test->Draw("AL");
  gr_PU200_2sigma_EB_test->Draw("same");
  leg_EB_sigma_test->Draw();
  c_sigma_EB_test->Print("test.pdf");

// for(int i=0; i<prompt_eff_PU200_EB.size(); i++) {
//    cout << "isolation cut2sigma: prompt_eff_PU200_EB.at(i)" << endl;
//  }

  // test
  ///////
  */


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
  gr_PU200_2sigma_EB_vtx->SetLineWidth(2); gr_PU200_2sigma_EB_vtx->SetLineColor(kGreen+2);

  gr_PU200_EB->SetMarkerStyle(8); gr_PU200_EB->SetMarkerSize(1.5);
  gr_noPU_EB->SetMarkerStyle(8); gr_noPU_EB->SetMarkerSize(1.5); gr_noPU_EB->SetMarkerColor(kGray+1);

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
  gr_PU200_2sigma_EE_vtx->SetLineWidth(2); gr_PU200_2sigma_EE_vtx->SetLineColor(kGreen+2);

  gr_PU200_EE->SetMarkerStyle(8); gr_PU200_EE->SetMarkerSize(1.5);
  gr_noPU_EE->SetMarkerStyle(8); gr_noPU_EE->SetMarkerSize(1.5); gr_noPU_EE->SetMarkerColor(kGray+1);


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
  leg_EB_dt->AddEntry(gr_noPU_EB, "no MTD noPU");
  leg_EB_dt->SetTextSize(0.03);

  TLegend* leg_EB_sigma = new TLegend(0.15, 0.68, 0.41, 0.88);
  leg_EB_sigma->SetMargin(0.2);
//  leg_EB_sigma->AddEntry(gr_PU200_gen_EB, "gen PU200");
//  leg_EB_sigma->AddEntry(gr_noPU_gen_EB, "gen noPU");
//  leg_EB_sigma->AddEntry(gr_PU200_2sigma_EB_vtx, "2sigma PU200 (PV, track)");
//  leg_EB_sigma->AddEntry(gr_PU200_2sigma_EB, "2sigma PU200 (muon, track)");
  leg_EB_sigma->AddEntry(gr_PU200_EB, "no MTD PU200");
  leg_EB_sigma->AddEntry(gr_noPU_EB, "no MTD noPU");
  leg_EB_sigma->SetTextSize(0.03);

  // Endcap region
  TLegend* leg_EE_dt = new TLegend(0.15, 0.65, 0.48, 0.85);
  leg_EE_dt->AddEntry(gr_PU200_gen_EE, "gen PU200");
  leg_EE_dt->AddEntry(gr_noPU_gen_EE, "gen noPU");
  leg_EE_dt->AddEntry(gr_PU200_40_EE, "40ps PU200");
//  leg_EE_dt->AddEntry(gr_PU200_60_EE, "60ps PU200");
//  leg_EE_dt->AddEntry(gr_PU200_80_EE, "80ps PU200");
//  leg_EE_dt->AddEntry(gr_PU200_100_EE, "100ps PU200");
  leg_EE_dt->AddEntry(gr_PU200_EE, "no MTD PU200");
  leg_EE_dt->AddEntry(gr_noPU_EE, "no MTD noPU");
  leg_EE_dt->SetTextSize(0.03);

  TLegend* leg_EE_sigma = new TLegend(0.15, 0.68, 0.41, 0.88);
  leg_EE_sigma->SetMargin(0.2);
//  leg_EE_sigma->AddEntry(gr_PU200_gen_EE, "gen PU200");
//  leg_EE_sigma->AddEntry(gr_noPU_gen_EE, "gen noPU");
//  leg_EE_sigma->AddEntry(gr_PU200_2sigma_EE_vtx, "2sigma PU200 (PV, track)");
//  leg_EE_sigma->AddEntry(gr_PU200_2sigma_EE, "2sigma PU200 (muon, track)");
  leg_EE_sigma->AddEntry(gr_PU200_EE, "no MTD PU200");
  leg_EE_sigma->AddEntry(gr_noPU_EE, "no MTD noPU");
  leg_EE_sigma->SetTextSize(0.03);


  ////////////////
  // Draw plots //
  ////////////////
  // Barrel region
  TCanvas* c_dt_EB = new TCanvas("c_dt_EB", "c_dt_EB", 1500, 1500);
  c_dt_EB->cd();
  c_dt_EB->SetGrid();
  gr_PU200_EB->GetXaxis()->SetRangeUser(0.85,1.);
  gr_PU200_EB->GetYaxis()->SetRangeUser(0.,0.4);
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
  c_sigma_EB->SetLeftMargin(0.12);
  gr_PU200_EB->GetXaxis()->SetRangeUser(0.,1.1);
  //gr_PU200_EB->GetXaxis()->SetRangeUser(0.80,1.);
  gr_PU200_EB->GetYaxis()->SetRangeUser(0.,0.2);
  gr_PU200_EB->SetTitle("");
  gr_PU200_EB->Draw("ALP");
//  gr_PU200_2sigma_EB->Draw("same");
//  gr_PU200_2sigma_EB_vtx->Draw("same");
//  gr_PU200_3sigma_EB->Draw("same");
//  gr_PU200_4sigma_EB->Draw("same");
  gr_noPU_EB->Draw("same LP");
//  gr_PU200_gen_EB->Draw("same");
//  gr_noPU_gen_EB->Draw("same");
  leg_EB_sigma->Draw();
  c_sigma_EB->Print("plots/roc_sigma_EB_reliso.pdf");

  // Endcap region
  TCanvas* c_dt_EE = new TCanvas("c_dt_EE", "c_dt_EE", 1500, 1500);
  c_dt_EE->cd();
  c_dt_EE->SetGrid();
  gr_PU200_EE->GetXaxis()->SetRangeUser(0.85,1.);
  gr_PU200_EE->GetYaxis()->SetRangeUser(0.,0.4);
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
  c_sigma_EE->SetLeftMargin(0.12);
  //gr_PU200_EE->GetXaxis()->SetRangeUser(0.,1.1);
  gr_PU200_EE->GetXaxis()->SetRangeUser(0.70,1.);
  gr_PU200_EE->GetYaxis()->SetRangeUser(0.,0.2);
  gr_PU200_EE->SetTitle("");
  gr_PU200_EE->Draw("ALP");
//  gr_PU200_2sigma_EE->Draw("same");
//  gr_PU200_2sigma_EE_vtx->Draw("same");
//  gr_PU200_3sigma_EE->Draw("same");
//  gr_PU200_4sigma_EE->Draw("same");
  gr_noPU_EE->Draw("same LP");
//  gr_PU200_gen_EE->Draw("same");
//  gr_noPU_gen_EE->Draw("same");
  leg_EE_sigma->Draw();
  c_sigma_EE->Print("plots/roc_sigma_EE_reliso.pdf");
}


void draw_pt_roc() {
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  TFile *f_PU200_prompt_vtx, *f_PU200_nonprompt_vtx, *f_noPU_prompt_vtx, *f_noPU_nonprompt_vtx;
  // PU200
  f_PU200_prompt = new TFile(Form("data/%s", path_PU200_prompt.Data()));
  f_PU200_nonprompt = new TFile(Form("data/%s", path_PU200_nonprompt.Data()));
  f_PU200_prompt_vtx = new TFile(Form("data/%s", path_PU200_prompt_vtx.Data()));
  f_PU200_nonprompt_vtx = new TFile(Form("data/%s", path_PU200_nonprompt_vtx.Data()));
  // noPU
  f_noPU_prompt = new TFile(Form("data/%s", path_noPU_prompt.Data()));
  f_noPU_nonprompt = new TFile(Form("data/%s", path_noPU_nonprompt.Data()));
  f_noPU_prompt_vtx = new TFile(Form("data/%s", path_noPU_prompt_vtx.Data()));
  f_noPU_nonprompt_vtx = new TFile(Form("data/%s", path_noPU_nonprompt_vtx.Data()));

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
    // vtx
  TH1D *h_PU200_prompt_2sigma_EB_vtx, *h_PU200_nonprompt_2sigma_EB_vtx, *h_noPU_prompt_2sigma_EB_vtx, *h_noPU_2sigma_nonprompt_EB_vtx;
  TH1D *h_PU200_prompt_3sigma_EB_vtx, *h_PU200_nonprompt_3sigma_EB_vtx, *h_noPU_prompt_3sigma_EB_vtx, *h_noPU_3sigma_nonprompt_EB_vtx;
  TH1D *h_PU200_prompt_4sigma_EB_vtx, *h_PU200_nonprompt_4sigma_EB_vtx, *h_noPU_prompt_4sigma_EB_vtx, *h_noPU_4sigma_nonprompt_EB_vtx;

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
    // vtx
  TH1D *h_PU200_prompt_2sigma_EE_vtx, *h_PU200_nonprompt_2sigma_EE_vtx, *h_noPU_prompt_2sigma_EE_vtx, *h_noPU_2sigma_nonprompt_EE_vtx;
  TH1D *h_PU200_prompt_3sigma_EE_vtx, *h_PU200_nonprompt_3sigma_EE_vtx, *h_noPU_prompt_3sigma_EE_vtx, *h_noPU_3sigma_nonprompt_EE_vtx;
  TH1D *h_PU200_prompt_4sigma_EE_vtx, *h_PU200_nonprompt_4sigma_EE_vtx, *h_noPU_prompt_4sigma_EE_vtx, *h_noPU_4sigma_nonprompt_EE_vtx;


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
      // vtx
  h_PU200_prompt_4sigma_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_chIso_sum_MTD_4sigma_Sig_EB");
  h_PU200_prompt_3sigma_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_chIso_sum_MTD_3sigma_Sig_EB");
  h_PU200_prompt_2sigma_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_chIso_sum_MTD_2sigma_Sig_EB");
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
      // vtx
  h_PU200_prompt_4sigma_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_chIso_sum_MTD_4sigma_Sig_EE");
  h_PU200_prompt_3sigma_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_chIso_sum_MTD_3sigma_Sig_EE");
  h_PU200_prompt_2sigma_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_chIso_sum_MTD_2sigma_Sig_EE");
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
      // vtx
  h_PU200_nonprompt_4sigma_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_chIso_sum_MTD_4sigma_Bkg_EB");
  h_PU200_nonprompt_3sigma_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_chIso_sum_MTD_3sigma_Bkg_EB");
  h_PU200_nonprompt_2sigma_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_chIso_sum_MTD_2sigma_Bkg_EB");
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
      // vtx
  h_PU200_nonprompt_4sigma_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_chIso_sum_MTD_4sigma_Bkg_EE");
  h_PU200_nonprompt_3sigma_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_chIso_sum_MTD_3sigma_Bkg_EE");
  h_PU200_nonprompt_2sigma_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_chIso_sum_MTD_2sigma_Bkg_EE");

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
        // vtx
  vector<double> prompt_eff_PU200_2sigma_EB_vtx={0}, prompt_norm_PU200_2sigma_EB_vtx={0}, prompt_eff_PU200_2sigma_EE_vtx={0}, prompt_norm_PU200_2sigma_EE_vtx={0};
  vector<double> prompt_eff_PU200_3sigma_EB_vtx={0}, prompt_norm_PU200_3sigma_EB_vtx={0}, prompt_eff_PU200_3sigma_EE_vtx={0}, prompt_norm_PU200_3sigma_EE_vtx={0};
  vector<double> prompt_eff_PU200_4sigma_EB_vtx={0}, prompt_norm_PU200_4sigma_EB_vtx={0}, prompt_eff_PU200_4sigma_EE_vtx={0}, prompt_norm_PU200_4sigma_EE_vtx={0};
    // Nonprompt
  vector<double> nonprompt_eff_PU200_EB={0}, nonprompt_norm_PU200_EB={0}, nonprompt_eff_PU200_EE={0}, nonprompt_norm_PU200_EE={0};
  vector<double> nonprompt_eff_PU200_2sigma_EB={0}, nonprompt_norm_PU200_2sigma_EB={0}, nonprompt_eff_PU200_2sigma_EE={0}, nonprompt_norm_PU200_2sigma_EE={0};
  vector<double> nonprompt_eff_PU200_3sigma_EB={0}, nonprompt_norm_PU200_3sigma_EB={0}, nonprompt_eff_PU200_3sigma_EE={0}, nonprompt_norm_PU200_3sigma_EE={0};
  vector<double> nonprompt_eff_PU200_4sigma_EB={0}, nonprompt_norm_PU200_4sigma_EB={0}, nonprompt_eff_PU200_4sigma_EE={0}, nonprompt_norm_PU200_4sigma_EE={0};
  vector<double> nonprompt_eff_PU200_40_EB={0}, nonprompt_norm_PU200_40_EB={0}, nonprompt_eff_PU200_40_EE={0}, nonprompt_norm_PU200_40_EE={0};
  vector<double> nonprompt_eff_PU200_60_EB={0}, nonprompt_norm_PU200_60_EB={0}, nonprompt_eff_PU200_60_EE={0}, nonprompt_norm_PU200_60_EE={0};
  vector<double> nonprompt_eff_PU200_80_EB={0}, nonprompt_norm_PU200_80_EB={0}, nonprompt_eff_PU200_80_EE={0}, nonprompt_norm_PU200_80_EE={0};
  vector<double> nonprompt_eff_PU200_100_EB={0}, nonprompt_norm_PU200_100_EB={0}, nonprompt_eff_PU200_100_EE={0}, nonprompt_norm_PU200_100_EE={0};
        // vtx
  vector<double> nonprompt_eff_PU200_2sigma_EB_vtx={0}, nonprompt_norm_PU200_2sigma_EB_vtx={0}, nonprompt_eff_PU200_2sigma_EE_vtx={0}, nonprompt_norm_PU200_2sigma_EE_vtx={0};
  vector<double> nonprompt_eff_PU200_3sigma_EB_vtx={0}, nonprompt_norm_PU200_3sigma_EB_vtx={0}, nonprompt_eff_PU200_3sigma_EE_vtx={0}, nonprompt_norm_PU200_3sigma_EE_vtx={0};
  vector<double> nonprompt_eff_PU200_4sigma_EB_vtx={0}, nonprompt_norm_PU200_4sigma_EB_vtx={0}, nonprompt_eff_PU200_4sigma_EE_vtx={0}, nonprompt_norm_PU200_4sigma_EE_vtx={0};
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
  for(unsigned int i=0; i<h_PU200_prompt_EB->GetNbinsX()+2; i++) {
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
        // vtx
    prompt_eff_PU200_2sigma_EB_vtx.emplace_back(h_PU200_prompt_2sigma_EB_vtx->Integral(1,i+1)/h_PU200_prompt_2sigma_EB_vtx->Integral(1,-1));
    prompt_eff_PU200_3sigma_EB_vtx.emplace_back(h_PU200_prompt_3sigma_EB_vtx->Integral(1,i+1)/h_PU200_prompt_3sigma_EB_vtx->Integral(1,-1));
    prompt_eff_PU200_4sigma_EB_vtx.emplace_back(h_PU200_prompt_4sigma_EB_vtx->Integral(1,i+1)/h_PU200_prompt_4sigma_EB_vtx->Integral(1,-1));
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
        // vtx
    prompt_eff_PU200_2sigma_EE_vtx.emplace_back(h_PU200_prompt_2sigma_EE_vtx->Integral(1,i+1)/h_PU200_prompt_2sigma_EE_vtx->Integral(1,-1));
    prompt_eff_PU200_3sigma_EE_vtx.emplace_back(h_PU200_prompt_3sigma_EE_vtx->Integral(1,i+1)/h_PU200_prompt_3sigma_EE_vtx->Integral(1,-1));
    prompt_eff_PU200_4sigma_EE_vtx.emplace_back(h_PU200_prompt_4sigma_EE_vtx->Integral(1,i+1)/h_PU200_prompt_4sigma_EE_vtx->Integral(1,-1));
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
        // vtx
    nonprompt_eff_PU200_2sigma_EB_vtx.emplace_back(h_PU200_nonprompt_2sigma_EB_vtx->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EB_vtx->Integral(1,-1));
    nonprompt_eff_PU200_3sigma_EB_vtx.emplace_back(h_PU200_nonprompt_3sigma_EB_vtx->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EB_vtx->Integral(1,-1));
    nonprompt_eff_PU200_4sigma_EB_vtx.emplace_back(h_PU200_nonprompt_4sigma_EB_vtx->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EB_vtx->Integral(1,-1));
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
        // vtx
    nonprompt_eff_PU200_2sigma_EE_vtx.emplace_back(h_PU200_nonprompt_2sigma_EE_vtx->Integral(1,i+1)/h_PU200_nonprompt_2sigma_EE_vtx->Integral(1,-1));
    nonprompt_eff_PU200_3sigma_EE_vtx.emplace_back(h_PU200_nonprompt_3sigma_EE_vtx->Integral(1,i+1)/h_PU200_nonprompt_3sigma_EE_vtx->Integral(1,-1));
    nonprompt_eff_PU200_4sigma_EE_vtx.emplace_back(h_PU200_nonprompt_4sigma_EE_vtx->Integral(1,i+1)/h_PU200_nonprompt_4sigma_EE_vtx->Integral(1,-1));
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
      // vtx
  TGraph* gr_PU200_2sigma_EB_vtx = new TGraph(h_PU200_prompt_2sigma_EB_vtx->GetNbinsX(), &prompt_eff_PU200_2sigma_EB_vtx[0], &nonprompt_eff_PU200_2sigma_EB_vtx[0]);
  TGraph* gr_PU200_3sigma_EB_vtx = new TGraph(h_PU200_prompt_3sigma_EB_vtx->GetNbinsX(), &prompt_eff_PU200_3sigma_EB_vtx[0], &nonprompt_eff_PU200_3sigma_EB_vtx[0]);
  TGraph* gr_PU200_4sigma_EB_vtx = new TGraph(h_PU200_prompt_4sigma_EB_vtx->GetNbinsX(), &prompt_eff_PU200_4sigma_EB_vtx[0], &nonprompt_eff_PU200_4sigma_EB_vtx[0]);
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
      // vtx
  TGraph* gr_PU200_2sigma_EE_vtx = new TGraph(h_PU200_prompt_2sigma_EE_vtx->GetNbinsX(), &prompt_eff_PU200_2sigma_EE_vtx[0], &nonprompt_eff_PU200_2sigma_EE_vtx[0]);
  TGraph* gr_PU200_3sigma_EE_vtx = new TGraph(h_PU200_prompt_3sigma_EE_vtx->GetNbinsX(), &prompt_eff_PU200_3sigma_EE_vtx[0], &nonprompt_eff_PU200_3sigma_EE_vtx[0]);
  TGraph* gr_PU200_4sigma_EE_vtx = new TGraph(h_PU200_prompt_4sigma_EE_vtx->GetNbinsX(), &prompt_eff_PU200_4sigma_EE_vtx[0], &nonprompt_eff_PU200_4sigma_EE_vtx[0]);

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
  gr_PU200_2sigma_EB_vtx->SetLineWidth(2); gr_PU200_2sigma_EB_vtx->SetLineColor(kGreen+2);

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
  gr_PU200_2sigma_EE_vtx->SetLineWidth(2); gr_PU200_2sigma_EE_vtx->SetLineColor(kGreen+2);


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
  leg_EB_dt->AddEntry(gr_noPU_EB, "no MTD noPU");
  leg_EB_dt->SetTextSize(0.03);

  TLegend* leg_EB_sigma = new TLegend(0.15, 0.58, 0.60, 0.88);
  leg_EB_sigma->SetMargin(0.2);
  leg_EB_sigma->AddEntry(gr_PU200_gen_EB, "gen PU200");
  leg_EB_sigma->AddEntry(gr_noPU_gen_EB, "gen noPU");
  leg_EB_sigma->AddEntry(gr_PU200_2sigma_EB_vtx, "2sigma PU200 (PV, track)");
  leg_EB_sigma->AddEntry(gr_PU200_2sigma_EB, "2sigma PU200 (muon, track)");
//  leg_EB_sigma->AddEntry(gr_PU200_3sigma_EB, "3sigma PU200");
//  leg_EB_sigma->AddEntry(gr_PU200_4sigma_EB, "4sigma PU200");
  leg_EB_sigma->AddEntry(gr_PU200_EB, "no MTD PU200");
  leg_EB_sigma->AddEntry(gr_noPU_EB, "no MTD noPU");
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
  leg_EE_dt->AddEntry(gr_noPU_EE, "no MTD noPU");
  leg_EE_dt->SetTextSize(0.03);

  TLegend* leg_EE_sigma = new TLegend(0.15, 0.58, 0.60, 0.88);
  leg_EE_sigma->SetMargin(0.2);
  leg_EE_sigma->AddEntry(gr_PU200_gen_EE, "gen PU200");
  leg_EE_sigma->AddEntry(gr_noPU_gen_EE, "gen noPU");
  leg_EE_sigma->AddEntry(gr_PU200_2sigma_EE_vtx, "2sigma PU200 (PV, track)");
  leg_EE_sigma->AddEntry(gr_PU200_2sigma_EE, "2sigma PU200 (muon, track)");
//  leg_EE_sigma->AddEntry(gr_PU200_3sigma_EE, "3sigma PU200");
//  leg_EE_sigma->AddEntry(gr_PU200_4sigma_EE, "4sigma PU200");
  leg_EE_sigma->AddEntry(gr_PU200_EE, "no MTD PU200");
  leg_EE_sigma->AddEntry(gr_noPU_EE, "no MTD noPU");
  leg_EE_sigma->SetTextSize(0.03);


  ////////////////
  // Draw plots //
  ////////////////
  // Barrel region
  TCanvas* c_dt_EB = new TCanvas("c_dt_EB", "c_dt_EB", 1500, 1500);
  c_dt_EB->cd();
  c_dt_EB->SetGrid();
  gr_PU200_EB->GetXaxis()->SetRangeUser(0.85,1.);
  gr_PU200_EB->GetYaxis()->SetRangeUser(0.,0.4);
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
  c_sigma_EB->SetLeftMargin(0.12);
  gr_PU200_EB->GetXaxis()->SetRangeUser(0.85,1.);
  gr_PU200_EB->GetYaxis()->SetRangeUser(0.,0.4);
  gr_PU200_EB->SetTitle("");
  gr_PU200_EB->Draw("AL");
  gr_PU200_2sigma_EB->Draw("same");
  gr_PU200_2sigma_EB_vtx->Draw("same");
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
  gr_PU200_EE->GetXaxis()->SetRangeUser(0.85,1.);
  gr_PU200_EE->GetYaxis()->SetRangeUser(0.,0.4);
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
  c_sigma_EE->SetLeftMargin(0.12);
  gr_PU200_EE->GetXaxis()->SetRangeUser(0.85,1.);
  gr_PU200_EE->GetYaxis()->SetRangeUser(0.,0.4);
  gr_PU200_EE->SetTitle("");
  gr_PU200_EE->Draw("AL");
  gr_PU200_2sigma_EE->Draw("same");
  gr_PU200_2sigma_EE_vtx->Draw("same");
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
  TFile *f_PU200_prompt_vtx, *f_PU200_nonprompt_vtx, *f_noPU_prompt_vtx, *f_noPU_nonprompt_vtx;
  // PU200
  f_PU200_prompt = new TFile(Form("data/%s", path_PU200_prompt.Data()));
  f_PU200_nonprompt = new TFile(Form("data/%s", path_PU200_nonprompt.Data()));
  f_PU200_prompt_vtx = new TFile(Form("data/%s", path_PU200_prompt_vtx.Data()));
  f_PU200_nonprompt_vtx = new TFile(Form("data/%s", path_PU200_nonprompt_vtx.Data()));
  // noPU
  f_noPU_prompt = new TFile(Form("data/%s", path_noPU_prompt.Data()));
  f_noPU_nonprompt = new TFile(Form("data/%s", path_noPU_nonprompt.Data()));
  f_noPU_prompt_vtx = new TFile(Form("data/%s", path_noPU_prompt_vtx.Data()));
  f_noPU_nonprompt_vtx = new TFile(Form("data/%s", path_noPU_nonprompt_vtx.Data()));

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
  TH1D *h_PU200_prompt_gen_2sigma_pT_EB, *h_PU200_nonprompt_gen_2sigma_pT_EB, *h_noPU_prompt_gen_2sigma_pT_EB, *h_noPU_nonprompt_gen_2sigma_pT_EB;
  TH1D *h_PU200_prompt_gen_3sigma_pT_EB, *h_PU200_nonprompt_gen_3sigma_pT_EB, *h_noPU_prompt_gen_3sigma_pT_EB, *h_noPU_nonprompt_gen_3sigma_pT_EB;
  TH1D *h_PU200_prompt_gen_4sigma_pT_EB, *h_PU200_nonprompt_gen_4sigma_pT_EB, *h_noPU_prompt_gen_4sigma_pT_EB, *h_noPU_nonprompt_gen_4sigma_pT_EB;
    // vtx
  TH1D *h_PU200_prompt_2sigma_pT_EB_vtx, *h_PU200_nonprompt_2sigma_pT_EB_vtx, *h_noPU_prompt_2sigma_pT_EB_vtx, *h_noPU_2sigma_nonprompt_pT_EB_vtx;
  TH1D *h_PU200_prompt_3sigma_pT_EB_vtx, *h_PU200_nonprompt_3sigma_pT_EB_vtx, *h_noPU_prompt_3sigma_pT_EB_vtx, *h_noPU_3sigma_nonprompt_pT_EB_vtx;
  TH1D *h_PU200_prompt_4sigma_pT_EB_vtx, *h_PU200_nonprompt_4sigma_pT_EB_vtx, *h_noPU_prompt_4sigma_pT_EB_vtx, *h_noPU_4sigma_nonprompt_pT_EB_vtx;
  TH1D *h_PU200_prompt_gen_2sigma_pT_EB_vtx, *h_PU200_nonprompt_gen_2sigma_pT_EB_vtx, *h_noPU_prompt_gen_2sigma_pT_EB_vtx, *h_noPU_nonprompt_gen_2sigma_pT_EB_vtx;
  TH1D *h_PU200_prompt_gen_3sigma_pT_EB_vtx, *h_PU200_nonprompt_gen_3sigma_pT_EB_vtx, *h_noPU_prompt_gen_3sigma_pT_EB_vtx, *h_noPU_nonprompt_gen_3sigma_pT_EB_vtx;
  TH1D *h_PU200_prompt_gen_4sigma_pT_EB_vtx, *h_PU200_nonprompt_gen_4sigma_pT_EB_vtx, *h_noPU_prompt_gen_4sigma_pT_EB_vtx, *h_noPU_nonprompt_gen_4sigma_pT_EB_vtx;

  // Histograms for Endcap region
  TH1D *h_PU200_prompt_pT_EE, *h_PU200_nonprompt_pT_EE, *h_noPU_prompt_pT_EE, *h_noPU_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_2sigma_pT_EE, *h_PU200_nonprompt_2sigma_pT_EE, *h_noPU_prompt_2sigma_pT_EE, *h_noPU_nonprompt_2sigma_pT_EE;
  TH1D *h_PU200_prompt_3sigma_pT_EE, *h_PU200_nonprompt_3sigma_pT_EE, *h_noPU_prompt_3sigma_pT_EE, *h_noPU_nonprompt_3sigma_pT_EE;
  TH1D *h_PU200_prompt_4sigma_pT_EE, *h_PU200_nonprompt_4sigma_pT_EE, *h_noPU_prompt_4sigma_pT_EE, *h_noPU_nonprompt_4sigma_pT_EE;
  TH1D *h_PU200_prompt_40_pT_EE, *h_PU200_nonprompt_40_pT_EE, *h_noPU_prompt_40_pT_EE, *h_noPU_40_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_60_pT_EE, *h_PU200_nonprompt_60_pT_EE, *h_noPU_prompt_60_pT_EE, *h_noPU_60_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_80_pT_EE, *h_PU200_nonprompt_80_pT_EE, *h_noPU_prompt_80_pT_EE, *h_noPU_80_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_100_pT_EE, *h_PU200_nonprompt_100_pT_EE, *h_noPU_prompt_100_pT_EE, *h_noPU_100_nonprompt_pT_EE;
  TH1D *h_PU200_prompt_gen_pT_EE, *h_PU200_nonprompt_gen_pT_EE, *h_noPU_prompt_gen_pT_EE, *h_noPU_nonprompt_gen_pT_EE;
  TH1D *h_PU200_prompt_gen_2sigma_pT_EE, *h_PU200_nonprompt_gen_2sigma_pT_EE, *h_noPU_prompt_gen_2sigma_pT_EE, *h_noPU_nonprompt_gen_2sigma_pT_EE;
  TH1D *h_PU200_prompt_gen_3sigma_pT_EE, *h_PU200_nonprompt_gen_3sigma_pT_EE, *h_noPU_prompt_gen_3sigma_pT_EE, *h_noPU_nonprompt_gen_3sigma_pT_EE;
  TH1D *h_PU200_prompt_gen_4sigma_pT_EE, *h_PU200_nonprompt_gen_4sigma_pT_EE, *h_noPU_prompt_gen_4sigma_pT_EE, *h_noPU_nonprompt_gen_4sigma_pT_EE;
    // vtx
  TH1D *h_PU200_prompt_2sigma_pT_EE_vtx, *h_PU200_nonprompt_2sigma_pT_EE_vtx, *h_noPU_prompt_2sigma_pT_EE_vtx, *h_noPU_nonprompt_2sigma_pT_EE_vtx;
  TH1D *h_PU200_prompt_3sigma_pT_EE_vtx, *h_PU200_nonprompt_3sigma_pT_EE_vtx, *h_noPU_prompt_3sigma_pT_EE_vtx, *h_noPU_nonprompt_3sigma_pT_EE_vtx;
  TH1D *h_PU200_prompt_4sigma_pT_EE_vtx, *h_PU200_nonprompt_4sigma_pT_EE_vtx, *h_noPU_prompt_4sigma_pT_EE_vtx, *h_noPU_nonprompt_4sigma_pT_EE_vtx;
  TH1D *h_PU200_prompt_gen_2sigma_pT_EE_vtx, *h_PU200_nonprompt_gen_2sigma_pT_EE_vtx, *h_noPU_prompt_gen_2sigma_pT_EE_vtx, *h_noPU_nonprompt_gen_2sigma_pT_EE_vtx;
  TH1D *h_PU200_prompt_gen_3sigma_pT_EE_vtx, *h_PU200_nonprompt_gen_3sigma_pT_EE_vtx, *h_noPU_prompt_gen_3sigma_pT_EE_vtx, *h_noPU_nonprompt_gen_3sigma_pT_EE_vtx;
  TH1D *h_PU200_prompt_gen_4sigma_pT_EE_vtx, *h_PU200_nonprompt_gen_4sigma_pT_EE_vtx, *h_noPU_prompt_gen_4sigma_pT_EE_vtx, *h_noPU_nonprompt_gen_4sigma_pT_EE_vtx;

  TString dir="/DQMData/Run 1/MTD/Run summary/MuonIso/";
  // Prompt
    // Barrel region
  h_PU200_prompt_pT_EB 	   		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffnoMTD_Sig_EB");
  h_PU200_prompt_4sigma_pT_EB   	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_4sigma_Sig_EB");
  h_PU200_prompt_3sigma_pT_EB 		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_3sigma_Sig_EB");
  h_PU200_prompt_2sigma_pT_EB 		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_2sigma_Sig_EB");
  h_PU200_prompt_40_pT_EB     		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_7_Sig_EB");
  h_PU200_prompt_60_pT_EB     		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_5_Sig_EB");
  h_PU200_prompt_80_pT_EB     		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_3_Sig_EB");
  h_PU200_prompt_100_pT_EB     		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_1_Sig_EB");
  h_PU200_prompt_gen_pT_EB		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_gen_Sig_EB");
  h_noPU_prompt_pT_EB 	   		= (TH1D*)f_noPU_prompt->Get(dir+"pTeffnoMTD_Sig_EB");            // noPU
  h_noPU_prompt_gen_pT_EB		= (TH1D*)f_noPU_prompt->Get(dir+"pTeffMTD_gen_Sig_EB");          // noPU
  h_PU200_prompt_gen_4sigma_pT_EB   	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_gen_4sigma_Sig_EB");
  h_PU200_prompt_gen_3sigma_pT_EB 	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_gen_3sigma_Sig_EB");
  h_PU200_prompt_gen_2sigma_pT_EB 	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_gen_2sigma_Sig_EB");
  h_noPU_prompt_gen_4sigma_pT_EB   	= (TH1D*)f_noPU_prompt->Get(dir+"pTeffMTD_gen_4sigma_Sig_EB");
  h_noPU_prompt_gen_3sigma_pT_EB 	= (TH1D*)f_noPU_prompt->Get(dir+"pTeffMTD_gen_3sigma_Sig_EB");
  h_noPU_prompt_gen_2sigma_pT_EB 	= (TH1D*)f_noPU_prompt->Get(dir+"pTeffMTD_gen_2sigma_Sig_EB");
      // vtx
  h_PU200_prompt_4sigma_pT_EB_vtx   	= (TH1D*)f_PU200_prompt_vtx->Get(dir+"pTeffMTD_4sigma_Sig_EB");
  h_PU200_prompt_3sigma_pT_EB_vtx	= (TH1D*)f_PU200_prompt_vtx->Get(dir+"pTeffMTD_3sigma_Sig_EB");
  h_PU200_prompt_2sigma_pT_EB_vtx 	= (TH1D*)f_PU200_prompt_vtx->Get(dir+"pTeffMTD_2sigma_Sig_EB");
  h_PU200_prompt_gen_4sigma_pT_EB_vtx  	= (TH1D*)f_PU200_prompt_vtx->Get(dir+"pTeffMTD_gen_4sigma_Sig_EB");
  h_PU200_prompt_gen_3sigma_pT_EB_vtx 	= (TH1D*)f_PU200_prompt_vtx->Get(dir+"pTeffMTD_gen_3sigma_Sig_EB");
  h_PU200_prompt_gen_2sigma_pT_EB_vtx 	= (TH1D*)f_PU200_prompt_vtx->Get(dir+"pTeffMTD_gen_2sigma_Sig_EB");
  h_noPU_prompt_gen_4sigma_pT_EB_vtx   	= (TH1D*)f_noPU_prompt_vtx->Get(dir+"pTeffMTD_gen_4sigma_Sig_EB");
  h_noPU_prompt_gen_3sigma_pT_EB_vtx 	= (TH1D*)f_noPU_prompt_vtx->Get(dir+"pTeffMTD_gen_3sigma_Sig_EB");
  h_noPU_prompt_gen_2sigma_pT_EB_vtx 	= (TH1D*)f_noPU_prompt_vtx->Get(dir+"pTeffMTD_gen_2sigma_Sig_EB");
    // Endcap region
  h_PU200_prompt_pT_EE 	   		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffnoMTD_Sig_EE");
  h_PU200_prompt_4sigma_pT_EE   	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_4sigma_Sig_EE");
  h_PU200_prompt_3sigma_pT_EE 		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_3sigma_Sig_EE");
  h_PU200_prompt_2sigma_pT_EE 		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_2sigma_Sig_EE");
  h_PU200_prompt_40_pT_EE     		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_7_Sig_EE");
  h_PU200_prompt_60_pT_EE     		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_5_Sig_EE");
  h_PU200_prompt_80_pT_EE     		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_3_Sig_EE");
  h_PU200_prompt_100_pT_EE     		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_1_Sig_EE");
  h_PU200_prompt_gen_pT_EE		= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_gen_Sig_EE");
  h_noPU_prompt_pT_EE 	   		= (TH1D*)f_noPU_prompt->Get(dir+"pTeffnoMTD_Sig_EE");            // noPU
  h_noPU_prompt_gen_pT_EE		= (TH1D*)f_noPU_prompt->Get(dir+"pTeffMTD_gen_Sig_EE");          // noPU
  h_PU200_prompt_gen_4sigma_pT_EE   	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_gen_4sigma_Sig_EE");
  h_PU200_prompt_gen_3sigma_pT_EE 	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_gen_3sigma_Sig_EE");
  h_PU200_prompt_gen_2sigma_pT_EE 	= (TH1D*)f_PU200_prompt->Get(dir+"pTeffMTD_gen_2sigma_Sig_EE");
  h_noPU_prompt_gen_4sigma_pT_EE   	= (TH1D*)f_noPU_prompt->Get(dir+"pTeffMTD_gen_4sigma_Sig_EE");
  h_noPU_prompt_gen_3sigma_pT_EE 	= (TH1D*)f_noPU_prompt->Get(dir+"pTeffMTD_gen_3sigma_Sig_EE");
  h_noPU_prompt_gen_2sigma_pT_EE 	= (TH1D*)f_noPU_prompt->Get(dir+"pTeffMTD_gen_2sigma_Sig_EE");
      // vtx
  h_PU200_prompt_4sigma_pT_EE_vtx   	= (TH1D*)f_PU200_prompt_vtx->Get(dir+"pTeffMTD_4sigma_Sig_EE");
  h_PU200_prompt_3sigma_pT_EE_vtx 	= (TH1D*)f_PU200_prompt_vtx->Get(dir+"pTeffMTD_3sigma_Sig_EE");
  h_PU200_prompt_2sigma_pT_EE_vtx 	= (TH1D*)f_PU200_prompt_vtx->Get(dir+"pTeffMTD_2sigma_Sig_EE");
  h_PU200_prompt_gen_4sigma_pT_EE_vtx  	= (TH1D*)f_PU200_prompt_vtx->Get(dir+"pTeffMTD_gen_4sigma_Sig_EE");
  h_PU200_prompt_gen_3sigma_pT_EE_vtx 	= (TH1D*)f_PU200_prompt_vtx->Get(dir+"pTeffMTD_gen_3sigma_Sig_EE");
  h_PU200_prompt_gen_2sigma_pT_EE_vtx 	= (TH1D*)f_PU200_prompt_vtx->Get(dir+"pTeffMTD_gen_2sigma_Sig_EE");
  h_noPU_prompt_gen_4sigma_pT_EE_vtx   	= (TH1D*)f_noPU_prompt_vtx->Get(dir+"pTeffMTD_gen_4sigma_Sig_EE");
  h_noPU_prompt_gen_3sigma_pT_EE_vtx 	= (TH1D*)f_noPU_prompt_vtx->Get(dir+"pTeffMTD_gen_3sigma_Sig_EE");
  h_noPU_prompt_gen_2sigma_pT_EE_vtx 	= (TH1D*)f_noPU_prompt_vtx->Get(dir+"pTeffMTD_gen_2sigma_Sig_EE");
  // Nonprompt
    // Barrel region
  h_PU200_nonprompt_pT_EB 	 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffnoMTD_Bkg_EB");
  h_PU200_nonprompt_4sigma_pT_EB 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_4sigma_Bkg_EB");
  h_PU200_nonprompt_3sigma_pT_EB 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_3sigma_Bkg_EB");
  h_PU200_nonprompt_2sigma_pT_EB 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_2sigma_Bkg_EB");
  h_PU200_nonprompt_40_pT_EB     	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_7_Bkg_EB");
  h_PU200_nonprompt_60_pT_EB     	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_5_Bkg_EB");
  h_PU200_nonprompt_80_pT_EB     	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_3_Bkg_EB");
  h_PU200_nonprompt_100_pT_EB    	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_1_Bkg_EB");
  h_PU200_nonprompt_gen_pT_EB	 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_gen_Bkg_EB");
  h_noPU_nonprompt_pT_EB 	 	= (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffnoMTD_Bkg_EB");         // noPU
  h_noPU_nonprompt_gen_pT_EB	 	= (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffMTD_gen_Bkg_EB");       // noPU
  h_PU200_nonprompt_gen_4sigma_pT_EB   	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_gen_4sigma_Bkg_EB");
  h_PU200_nonprompt_gen_3sigma_pT_EB 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_gen_3sigma_Bkg_EB");
  h_PU200_nonprompt_gen_2sigma_pT_EB 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_gen_2sigma_Bkg_EB");
  h_noPU_nonprompt_gen_4sigma_pT_EB   	= (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffMTD_gen_4sigma_Bkg_EB");
  h_noPU_nonprompt_gen_3sigma_pT_EB 	= (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffMTD_gen_3sigma_Bkg_EB");
  h_noPU_nonprompt_gen_2sigma_pT_EB 	= (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffMTD_gen_2sigma_Bkg_EB");
      // vtx
  h_PU200_nonprompt_4sigma_pT_EB_vtx 	= (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"pTeffMTD_4sigma_Bkg_EB");
  h_PU200_nonprompt_3sigma_pT_EB_vtx 	= (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"pTeffMTD_3sigma_Bkg_EB");
  h_PU200_nonprompt_2sigma_pT_EB_vtx 	= (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"pTeffMTD_2sigma_Bkg_EB");
  h_PU200_nonprompt_gen_4sigma_pT_EB_vtx   	= (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"pTeffMTD_gen_4sigma_Bkg_EB");
  h_PU200_nonprompt_gen_3sigma_pT_EB_vtx 	= (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"pTeffMTD_gen_3sigma_Bkg_EB");
  h_PU200_nonprompt_gen_2sigma_pT_EB_vtx 	= (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"pTeffMTD_gen_2sigma_Bkg_EB");
  h_noPU_nonprompt_gen_4sigma_pT_EB_vtx   	= (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"pTeffMTD_gen_4sigma_Bkg_EB");
  h_noPU_nonprompt_gen_3sigma_pT_EB_vtx 	= (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"pTeffMTD_gen_3sigma_Bkg_EB");
  h_noPU_nonprompt_gen_2sigma_pT_EB_vtx 	= (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"pTeffMTD_gen_2sigma_Bkg_EB");
    // Endcap region
  h_PU200_nonprompt_pT_EE 	 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffnoMTD_Bkg_EE");
  h_PU200_nonprompt_4sigma_pT_EE 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_4sigma_Bkg_EE");
  h_PU200_nonprompt_3sigma_pT_EE 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_3sigma_Bkg_EE");
  h_PU200_nonprompt_2sigma_pT_EE 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_2sigma_Bkg_EE");
  h_PU200_nonprompt_40_pT_EE     	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_7_Bkg_EE");
  h_PU200_nonprompt_60_pT_EE     	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_5_Bkg_EE");
  h_PU200_nonprompt_80_pT_EE     	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_3_Bkg_EE");
  h_PU200_nonprompt_100_pT_EE    	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_1_Bkg_EE");
  h_PU200_nonprompt_gen_pT_EE	 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_gen_Bkg_EE");
  h_noPU_nonprompt_pT_EE 	 	= (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffnoMTD_Bkg_EE");         // noPU
  h_noPU_nonprompt_gen_pT_EE	 	= (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffMTD_gen_Bkg_EE");       // noPU
  h_PU200_nonprompt_gen_4sigma_pT_EE   	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_gen_4sigma_Bkg_EE");
  h_PU200_nonprompt_gen_3sigma_pT_EE 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_gen_3sigma_Bkg_EE");
  h_PU200_nonprompt_gen_2sigma_pT_EE 	= (TH1D*)f_PU200_nonprompt->Get(dir+"pTeffMTD_gen_2sigma_Bkg_EE");
  h_noPU_nonprompt_gen_4sigma_pT_EE   	= (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffMTD_gen_4sigma_Bkg_EE");
  h_noPU_nonprompt_gen_3sigma_pT_EE 	= (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffMTD_gen_3sigma_Bkg_EE");
  h_noPU_nonprompt_gen_2sigma_pT_EE 	= (TH1D*)f_noPU_nonprompt->Get(dir+"pTeffMTD_gen_2sigma_Bkg_EE");
      // vtx
  h_PU200_nonprompt_4sigma_pT_EE_vtx 	= (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"pTeffMTD_4sigma_Bkg_EE");
  h_PU200_nonprompt_3sigma_pT_EE_vtx 	= (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"pTeffMTD_3sigma_Bkg_EE");
  h_PU200_nonprompt_2sigma_pT_EE_vtx 	= (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"pTeffMTD_2sigma_Bkg_EE");
  h_PU200_nonprompt_gen_4sigma_pT_EE_vtx   	= (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"pTeffMTD_gen_4sigma_Bkg_EE");
  h_PU200_nonprompt_gen_3sigma_pT_EE_vtx 	= (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"pTeffMTD_gen_3sigma_Bkg_EE");
  h_PU200_nonprompt_gen_2sigma_pT_EE_vtx 	= (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"pTeffMTD_gen_2sigma_Bkg_EE");
  h_noPU_nonprompt_gen_4sigma_pT_EE_vtx   	= (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"pTeffMTD_gen_4sigma_Bkg_EE");
  h_noPU_nonprompt_gen_3sigma_pT_EE_vtx 	= (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"pTeffMTD_gen_3sigma_Bkg_EE");
  h_noPU_nonprompt_gen_2sigma_pT_EE_vtx 	= (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"pTeffMTD_gen_2sigma_Bkg_EE");
 
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
  //h_PU200_prompt_incl_pT_EB->SetTitle("Prompt muon in barrel region");
  h_PU200_prompt_incl_pT_EB->SetTitle("");
  h_PU200_prompt_incl_pT_EB->GetXaxis()->SetTitle("Muon pT (GeV)");
  h_PU200_prompt_incl_pT_EB->GetYaxis()->SetTitle("Efficiency");
//  h_PU200_prompt_incl_pT_EB->GetYaxis()->SetRangeUser(0.45, 1.05);
  h_PU200_prompt_incl_pT_EB->GetYaxis()->SetRangeUser(0., 1.05);
  h_PU200_prompt_gen_pT_EB->SetLineColor(kBlack); h_PU200_prompt_pT_EB->SetLineColor(kBlack); h_PU200_prompt_40_pT_EB->SetLineColor(kRed); h_PU200_prompt_2sigma_pT_EB->SetLineColor(kRed); h_noPU_prompt_gen_pT_EB->SetLineColor(kGray); h_noPU_prompt_pT_EB->SetLineColor(kGray);
  h_PU200_prompt_gen_pT_EB->SetLineWidth(2); h_PU200_prompt_pT_EB->SetLineWidth(2); h_PU200_prompt_40_pT_EB->SetLineWidth(2); h_PU200_prompt_2sigma_pT_EB->SetLineWidth(2); h_noPU_prompt_gen_pT_EB->SetLineWidth(2); h_noPU_prompt_pT_EB->SetLineWidth(2);
  h_PU200_prompt_gen_pT_EB->SetLineStyle(2); h_noPU_prompt_gen_pT_EB->SetLineStyle(2);
  h_PU200_prompt_gen_2sigma_pT_EB->SetLineColor(kBlack); h_noPU_prompt_gen_2sigma_pT_EB->SetLineColor(kGray);
  h_PU200_prompt_gen_2sigma_pT_EB->SetLineWidth(2); h_noPU_prompt_gen_2sigma_pT_EB->SetLineWidth(2);
  h_PU200_prompt_gen_2sigma_pT_EB->SetLineStyle(2); h_noPU_prompt_gen_2sigma_pT_EB->SetLineStyle(2);
  h_PU200_prompt_2sigma_pT_EB_vtx->SetLineWidth(2); h_PU200_prompt_2sigma_pT_EB_vtx->SetLineColor(kGreen+2);
  h_PU200_prompt_gen_2sigma_pT_EB_vtx->SetLineColor(kBlack); h_noPU_prompt_gen_2sigma_pT_EB_vtx->SetLineColor(kGray);
  h_PU200_prompt_gen_2sigma_pT_EB_vtx->SetLineWidth(2); h_noPU_prompt_gen_2sigma_pT_EB_vtx->SetLineWidth(2);
  h_PU200_prompt_gen_2sigma_pT_EB_vtx->SetLineStyle(2); h_noPU_prompt_gen_2sigma_pT_EB_vtx->SetLineStyle(2);
    // Endcap region
  //h_PU200_prompt_incl_pT_EE->SetTitle("Prompt muon in endcap region");
  h_PU200_prompt_incl_pT_EE->SetTitle("");
  h_PU200_prompt_incl_pT_EE->GetXaxis()->SetTitle("Muon pT (GeV)");
  h_PU200_prompt_incl_pT_EE->GetYaxis()->SetTitle("Efficiency");
  h_PU200_prompt_incl_pT_EE->GetYaxis()->SetRangeUser(0., 1.05);
  h_PU200_prompt_gen_pT_EE->SetLineColor(kBlack); h_PU200_prompt_pT_EE->SetLineColor(kBlack); h_PU200_prompt_40_pT_EE->SetLineColor(kRed); h_PU200_prompt_2sigma_pT_EE->SetLineColor(kRed); h_noPU_prompt_gen_pT_EE->SetLineColor(kGray); h_noPU_prompt_pT_EE->SetLineColor(kGray);
  h_PU200_prompt_gen_pT_EE->SetLineWidth(2); h_PU200_prompt_pT_EE->SetLineWidth(2); h_PU200_prompt_40_pT_EE->SetLineWidth(2); h_PU200_prompt_2sigma_pT_EE->SetLineWidth(2); h_noPU_prompt_gen_pT_EE->SetLineWidth(2); h_noPU_prompt_pT_EE->SetLineWidth(2);
  h_PU200_prompt_gen_pT_EE->SetLineStyle(2); h_noPU_prompt_gen_pT_EE->SetLineStyle(2);
  h_PU200_prompt_gen_2sigma_pT_EE->SetLineColor(kBlack); h_noPU_prompt_gen_2sigma_pT_EE->SetLineColor(kGray);
  h_PU200_prompt_gen_2sigma_pT_EE->SetLineWidth(2); h_noPU_prompt_gen_2sigma_pT_EE->SetLineWidth(2);
  h_PU200_prompt_gen_2sigma_pT_EE->SetLineStyle(2); h_noPU_prompt_gen_2sigma_pT_EE->SetLineStyle(2);
  h_PU200_prompt_2sigma_pT_EE_vtx->SetLineWidth(2); h_PU200_prompt_2sigma_pT_EE_vtx->SetLineColor(kGreen+2);
  h_PU200_prompt_gen_2sigma_pT_EE_vtx->SetLineColor(kBlack); h_noPU_prompt_gen_2sigma_pT_EE_vtx->SetLineColor(kGray);
  h_PU200_prompt_gen_2sigma_pT_EE_vtx->SetLineWidth(2); h_noPU_prompt_gen_2sigma_pT_EE_vtx->SetLineWidth(2);
  h_PU200_prompt_gen_2sigma_pT_EE_vtx->SetLineStyle(2); h_noPU_prompt_gen_2sigma_pT_EE_vtx->SetLineStyle(2);

  // Nonprompt
    // Barrel region
  //h_PU200_nonprompt_incl_pT_EB->SetTitle("Non-prompt muon in barrel region");
  h_PU200_nonprompt_incl_pT_EB->SetTitle("");
  h_PU200_nonprompt_incl_pT_EB->GetXaxis()->SetTitle("Muon pT (GeV)");
  h_PU200_nonprompt_incl_pT_EB->GetYaxis()->SetTitle("Efficiency");
  h_PU200_nonprompt_incl_pT_EB->GetYaxis()->SetRangeUser(0., 1.05);
  h_PU200_nonprompt_gen_pT_EB->SetLineColor(kBlack); h_PU200_nonprompt_pT_EB->SetLineColor(kBlack); h_PU200_nonprompt_40_pT_EB->SetLineColor(kRed); h_PU200_nonprompt_2sigma_pT_EB->SetLineColor(kRed); h_noPU_nonprompt_gen_pT_EB->SetLineColor(kGray); h_noPU_nonprompt_pT_EB->SetLineColor(kGray);
  h_PU200_nonprompt_gen_pT_EB->SetLineWidth(2); h_PU200_nonprompt_pT_EB->SetLineWidth(2); h_PU200_nonprompt_40_pT_EB->SetLineWidth(2); h_PU200_nonprompt_2sigma_pT_EB->SetLineWidth(2); h_noPU_nonprompt_gen_pT_EB->SetLineWidth(2); h_noPU_nonprompt_pT_EB->SetLineWidth(2);
  h_PU200_nonprompt_gen_pT_EB->SetLineStyle(2); h_noPU_nonprompt_gen_pT_EB->SetLineStyle(2);
  h_PU200_nonprompt_gen_2sigma_pT_EB->SetLineColor(kBlack); h_noPU_nonprompt_gen_2sigma_pT_EB->SetLineColor(kGray);
  h_PU200_nonprompt_gen_2sigma_pT_EB->SetLineWidth(2); h_noPU_nonprompt_gen_2sigma_pT_EB->SetLineWidth(2);
  h_PU200_nonprompt_gen_2sigma_pT_EB->SetLineStyle(2); h_noPU_nonprompt_gen_2sigma_pT_EB->SetLineStyle(2);
  h_PU200_nonprompt_2sigma_pT_EB_vtx->SetLineWidth(2); h_PU200_nonprompt_2sigma_pT_EB_vtx->SetLineColor(kGreen+2);
  h_PU200_nonprompt_gen_2sigma_pT_EB_vtx->SetLineColor(kBlack); h_noPU_nonprompt_gen_2sigma_pT_EB_vtx->SetLineColor(kGray);
  h_PU200_nonprompt_gen_2sigma_pT_EB_vtx->SetLineWidth(2); h_noPU_nonprompt_gen_2sigma_pT_EB_vtx->SetLineWidth(2);
  h_PU200_nonprompt_gen_2sigma_pT_EB_vtx->SetLineStyle(2); h_noPU_nonprompt_gen_2sigma_pT_EB_vtx->SetLineStyle(2);
    // Endcap region
  //h_PU200_nonprompt_incl_pT_EE->SetTitle("Non-prompt muon in endcap region");
  h_PU200_nonprompt_incl_pT_EE->SetTitle("");
  h_PU200_nonprompt_incl_pT_EE->GetXaxis()->SetTitle("Muon pT (GeV)");
  h_PU200_nonprompt_incl_pT_EE->GetYaxis()->SetTitle("Efficiency");
  h_PU200_nonprompt_incl_pT_EE->GetYaxis()->SetRangeUser(0., 1.05);
  h_PU200_nonprompt_gen_pT_EE->SetLineColor(kBlack); h_PU200_nonprompt_pT_EE->SetLineColor(kBlack); h_PU200_nonprompt_40_pT_EE->SetLineColor(kRed); h_PU200_nonprompt_2sigma_pT_EE->SetLineColor(kRed); h_noPU_nonprompt_gen_pT_EE->SetLineColor(kGray); h_noPU_nonprompt_pT_EE->SetLineColor(kGray);
  h_PU200_nonprompt_gen_pT_EE->SetLineWidth(2); h_PU200_nonprompt_pT_EE->SetLineWidth(2); h_PU200_nonprompt_40_pT_EE->SetLineWidth(2); h_PU200_nonprompt_2sigma_pT_EE->SetLineWidth(2); h_noPU_nonprompt_gen_pT_EE->SetLineWidth(2); h_noPU_nonprompt_pT_EE->SetLineWidth(2);
  h_PU200_nonprompt_gen_pT_EE->SetLineStyle(2); h_noPU_nonprompt_gen_pT_EE->SetLineStyle(2);
  h_PU200_nonprompt_gen_2sigma_pT_EE->SetLineColor(kBlack); h_noPU_nonprompt_gen_2sigma_pT_EE->SetLineColor(kGray);
  h_PU200_nonprompt_gen_2sigma_pT_EE->SetLineWidth(2); h_noPU_nonprompt_gen_2sigma_pT_EE->SetLineWidth(2);
  h_PU200_nonprompt_gen_2sigma_pT_EE->SetLineStyle(2); h_noPU_nonprompt_gen_2sigma_pT_EE->SetLineStyle(2);
  h_PU200_nonprompt_2sigma_pT_EE_vtx->SetLineWidth(2); h_PU200_nonprompt_2sigma_pT_EE_vtx->SetLineColor(kGreen+2);
  h_PU200_nonprompt_gen_2sigma_pT_EE_vtx->SetLineColor(kBlack); h_noPU_nonprompt_gen_2sigma_pT_EE_vtx->SetLineColor(kGray);
  h_PU200_nonprompt_gen_2sigma_pT_EE_vtx->SetLineWidth(2); h_noPU_nonprompt_gen_2sigma_pT_EE_vtx->SetLineWidth(2);
  h_PU200_nonprompt_gen_2sigma_pT_EE_vtx->SetLineStyle(2); h_noPU_nonprompt_gen_2sigma_pT_EE_vtx->SetLineStyle(2);


  /////////////
  // Legends //
  /////////////
  // Prompt
    // Barrel region
  TLegend* leg_prompt_incl_pT_EB_sigma = new TLegend(0.43, 0.13, 0.88, 0.43);
  leg_prompt_incl_pT_EB_sigma->SetMargin(0.2);
  leg_prompt_incl_pT_EB_sigma->AddEntry(h_PU200_prompt_gen_pT_EB, "gen PU200");
  leg_prompt_incl_pT_EB_sigma->AddEntry(h_noPU_prompt_gen_pT_EB, "gen noPU");
  leg_prompt_incl_pT_EB_sigma->AddEntry(h_PU200_prompt_2sigma_pT_EB_vtx, "2sigma PU200 (PV, track)");
  leg_prompt_incl_pT_EB_sigma->AddEntry(h_PU200_prompt_2sigma_pT_EB, "2sigma PU200 (muon, track)");
  leg_prompt_incl_pT_EB_sigma->AddEntry(h_PU200_prompt_pT_EB, "no MTD PU200");
  leg_prompt_incl_pT_EB_sigma->AddEntry(h_noPU_prompt_pT_EB, "no MTD noPU");
  leg_prompt_incl_pT_EB_sigma->SetTextSize(0.03);
    // Endcap region
  TLegend* leg_prompt_incl_pT_EE_sigma = new TLegend(0.43, 0.13, 0.88, 0.43);
  leg_prompt_incl_pT_EE_sigma->SetMargin(0.2);
  leg_prompt_incl_pT_EE_sigma->AddEntry(h_PU200_prompt_gen_pT_EE, "gen PU200");
  leg_prompt_incl_pT_EE_sigma->AddEntry(h_noPU_prompt_gen_pT_EE, "gen noPU");
  leg_prompt_incl_pT_EE_sigma->AddEntry(h_PU200_prompt_2sigma_pT_EE_vtx, "2sigma PU200 (PV, track)");
  leg_prompt_incl_pT_EE_sigma->AddEntry(h_PU200_prompt_2sigma_pT_EE, "2sigma PU200 (muon, track)");
  leg_prompt_incl_pT_EE_sigma->AddEntry(h_PU200_prompt_pT_EE, "no MTD PU200");
  leg_prompt_incl_pT_EE_sigma->AddEntry(h_noPU_prompt_pT_EE, "no MTD noPU");
  leg_prompt_incl_pT_EE_sigma->SetTextSize(0.03);
  // Nonprompt
    // Barrel region
  TLegend* leg_nonprompt_incl_pT_EB_sigma = new TLegend(0.15, 0.58, 0.60, 0.88);
  leg_nonprompt_incl_pT_EB_sigma->SetMargin(0.2);
  leg_nonprompt_incl_pT_EB_sigma->AddEntry(h_PU200_nonprompt_gen_pT_EB, "gen PU200");
  leg_nonprompt_incl_pT_EB_sigma->AddEntry(h_noPU_nonprompt_gen_pT_EB, "gen noPU");
  leg_nonprompt_incl_pT_EB_sigma->AddEntry(h_PU200_nonprompt_2sigma_pT_EB_vtx, "2sigma PU200 (PV, track)");
  leg_nonprompt_incl_pT_EB_sigma->AddEntry(h_PU200_nonprompt_2sigma_pT_EB, "2sigma PU200 (muon, track)");
  leg_nonprompt_incl_pT_EB_sigma->AddEntry(h_PU200_nonprompt_pT_EB, "no MTD PU200");
  leg_nonprompt_incl_pT_EB_sigma->AddEntry(h_noPU_nonprompt_pT_EB, "no MTD noPU");
  leg_nonprompt_incl_pT_EB_sigma->SetTextSize(0.03);
    // Endcap region
  TLegend* leg_nonprompt_incl_pT_EE_sigma = new TLegend(0.15, 0.58, 0.60, 0.88);
  leg_nonprompt_incl_pT_EE_sigma->SetMargin(0.2);
  leg_nonprompt_incl_pT_EE_sigma->AddEntry(h_PU200_nonprompt_gen_pT_EE, "gen PU200");
  leg_nonprompt_incl_pT_EE_sigma->AddEntry(h_noPU_nonprompt_gen_pT_EE, "gen noPU");
  leg_nonprompt_incl_pT_EE_sigma->AddEntry(h_PU200_nonprompt_2sigma_pT_EE_vtx, "2sigma PU200 (PV, track)");
  leg_nonprompt_incl_pT_EE_sigma->AddEntry(h_PU200_nonprompt_2sigma_pT_EE, "2sigma PU200 (muon, track)");
  leg_nonprompt_incl_pT_EE_sigma->AddEntry(h_PU200_nonprompt_pT_EE, "no MTD PU200");
  leg_nonprompt_incl_pT_EE_sigma->AddEntry(h_noPU_nonprompt_pT_EE, "no MTD noPU");
  leg_nonprompt_incl_pT_EE_sigma->SetTextSize(0.03);
  // test
  // prompt_EB
  TLegend* leg_prompt_test_EB = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_prompt_test_EB->SetMargin(0.2);
  leg_prompt_test_EB->AddEntry(h_PU200_prompt_gen_2sigma_pT_EB, "gen 2sigma PU200");
  leg_prompt_test_EB->AddEntry(h_PU200_prompt_pT_EB, "no MTD PU200");
  leg_prompt_test_EB->AddEntry(h_PU200_prompt_2sigma_pT_EB, "2sigma PU200");
  leg_prompt_test_EB->AddEntry(h_noPU_prompt_gen_2sigma_pT_EB, "gen 2sigma noPU");
  leg_prompt_test_EB->AddEntry(h_noPU_prompt_pT_EB, "no MTD noPU");
  leg_prompt_test_EB->SetTextSize(0.03);
  // prompt_EE
  TLegend* leg_prompt_test_EE = new TLegend(0.55, 0.13, 0.88, 0.38);
  leg_prompt_test_EE->SetMargin(0.2);
  leg_prompt_test_EE->AddEntry(h_PU200_prompt_gen_2sigma_pT_EE, "gen 2sigma PU200");
  leg_prompt_test_EE->AddEntry(h_PU200_prompt_pT_EE, "no MTD PU200");
  leg_prompt_test_EE->AddEntry(h_PU200_prompt_2sigma_pT_EE, "2sigma PU200");
  leg_prompt_test_EE->AddEntry(h_noPU_prompt_gen_2sigma_pT_EE, "gen 2sigma noPU");
  leg_prompt_test_EE->AddEntry(h_noPU_prompt_pT_EE, "no MTD noPU");
  leg_prompt_test_EE->SetTextSize(0.03);
  // nonprompt_EB
  TLegend* leg_nonprompt_test_EB = new TLegend(0.15, 0.58, 0.63, 0.88);
  leg_nonprompt_test_EB->SetMargin(0.2);
  leg_nonprompt_test_EB->AddEntry(h_PU200_nonprompt_gen_2sigma_pT_EB_vtx, "gen 2sigma PU200 (PV, track)");
  leg_nonprompt_test_EB->AddEntry(h_noPU_nonprompt_gen_2sigma_pT_EB_vtx, "gen 2sigma noPU (PV, track)");
  leg_nonprompt_test_EB->AddEntry(h_PU200_nonprompt_2sigma_pT_EB_vtx, "2sigma PU200 (PV, track)");
  leg_nonprompt_test_EB->AddEntry(h_PU200_nonprompt_2sigma_pT_EB, "2sigma PU200 (muon, track)");
  leg_nonprompt_test_EB->AddEntry(h_PU200_nonprompt_pT_EB, "no MTD PU200");
  leg_nonprompt_test_EB->AddEntry(h_noPU_nonprompt_pT_EB, "no MTD noPU");
  leg_nonprompt_test_EB->SetTextSize(0.03);
  // nonprompt_EE
  TLegend* leg_nonprompt_test_EE = new TLegend(0.15, 0.58, 0.63, 0.88);
  leg_nonprompt_test_EE->SetMargin(0.2);
  leg_nonprompt_test_EE->AddEntry(h_PU200_nonprompt_gen_2sigma_pT_EE_vtx, "gen 2sigma PU200 (PV, track)");
  leg_nonprompt_test_EE->AddEntry(h_noPU_nonprompt_gen_2sigma_pT_EE_vtx, "gen 2sigma noPU (PV, track)");
  leg_nonprompt_test_EE->AddEntry(h_PU200_nonprompt_2sigma_pT_EE_vtx, "2sigma PU200 (PV, track)");
  leg_nonprompt_test_EE->AddEntry(h_PU200_nonprompt_2sigma_pT_EE, "2sigma PU200 (muon, track)");
  leg_nonprompt_test_EE->AddEntry(h_PU200_nonprompt_pT_EE, "no MTD PU200");
  leg_nonprompt_test_EE->AddEntry(h_noPU_nonprompt_pT_EE, "no MTD noPU");
  leg_nonprompt_test_EE->SetTextSize(0.03);


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
  h_PU200_prompt_2sigma_pT_EB_vtx->Draw("same e");
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
  h_PU200_prompt_2sigma_pT_EE_vtx->Draw("same e");
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
  h_PU200_nonprompt_2sigma_pT_EB_vtx->Draw("same e");
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
  h_PU200_nonprompt_2sigma_pT_EE_vtx->Draw("same e");
  h_noPU_nonprompt_gen_pT_EE->Draw("same e");
  h_noPU_nonprompt_pT_EE->Draw("same e");
  leg_nonprompt_incl_pT_EE_sigma->Draw();
  c_PU200_nonprompt_pT_EE_sigma->Print("plots/pT_PU200_nonprompt_EE_sigma.pdf");

  TCanvas* c_test_prompt_EB = new TCanvas("c_test_prompt_EB", "c_test_prompt_EB", 1500, 1500);
  c_test_prompt_EB->cd();
  c_test_prompt_EB->SetGrid();
  c_test_prompt_EB->SetLeftMargin(0.12);
  h_PU200_prompt_incl_pT_EB->Draw("hist");
  h_PU200_prompt_pT_EB->Draw("same e");
  h_PU200_prompt_2sigma_pT_EB->Draw("same e");
  h_PU200_prompt_gen_2sigma_pT_EB->Draw("same e");
  h_noPU_prompt_pT_EB->Draw("same e");
  h_noPU_prompt_gen_2sigma_pT_EB->Draw("same e");
  leg_prompt_test_EB->Draw();
  c_test_prompt_EB->Print("plots/pT_PU200_prompt_EB_2sigma_test.pdf");

  TCanvas* c_test_prompt_EE = new TCanvas("c_test_prompt_EE", "c_test_prompt_EE", 1500, 1500);
  c_test_prompt_EE->cd();
  c_test_prompt_EE->SetGrid();
  c_test_prompt_EE->SetLeftMargin(0.12);
  h_PU200_prompt_incl_pT_EE->Draw("hist");
  h_PU200_prompt_pT_EE->Draw("same e");
  h_PU200_prompt_2sigma_pT_EE->Draw("same e");
  h_PU200_prompt_gen_2sigma_pT_EE->Draw("same e");
  h_noPU_prompt_pT_EE->Draw("same e");
  h_noPU_prompt_gen_2sigma_pT_EE->Draw("same e");
  leg_prompt_test_EE->Draw();
  c_test_prompt_EE->Print("plots/pT_PU200_prompt_EE_2sigma_test.pdf");

  TCanvas* c_test_nonprompt_EB = new TCanvas("c_test_nonprompt_EB", "c_test_nonprompt_EB", 1500, 1500);
  c_test_nonprompt_EB->cd();
  c_test_nonprompt_EB->SetGrid();
  c_test_nonprompt_EB->SetLeftMargin(0.12);
  h_PU200_nonprompt_incl_pT_EB->Draw("hist");
  h_PU200_nonprompt_pT_EB->Draw("same e");
  h_PU200_nonprompt_2sigma_pT_EB->Draw("same e");
  h_PU200_nonprompt_2sigma_pT_EB_vtx->Draw("same e");
  h_PU200_nonprompt_gen_2sigma_pT_EB_vtx->Draw("same e");
  h_noPU_nonprompt_pT_EB->Draw("same e");
  h_noPU_nonprompt_gen_2sigma_pT_EB_vtx->Draw("same e");
  leg_nonprompt_test_EB->Draw();
  c_test_nonprompt_EB->Print("plots/pT_PU200_nonprompt_EB_2sigma_test.pdf");

  TCanvas* c_test_nonprompt_EE = new TCanvas("c_test_nonprompt_EE", "c_test_nonprompt_EE", 1500, 1500);
  c_test_nonprompt_EE->cd();
  c_test_nonprompt_EE->SetGrid();
  c_test_nonprompt_EE->SetLeftMargin(0.12);
  h_PU200_nonprompt_incl_pT_EE->Draw("hist");
  h_PU200_nonprompt_pT_EE->Draw("same e");
  h_PU200_nonprompt_2sigma_pT_EE->Draw("same e");
  h_PU200_nonprompt_2sigma_pT_EE_vtx->Draw("same e");
  h_PU200_nonprompt_gen_2sigma_pT_EE_vtx->Draw("same e");
  h_noPU_nonprompt_pT_EE->Draw("same e");
  h_noPU_nonprompt_gen_2sigma_pT_EE_vtx->Draw("same e");
  leg_nonprompt_test_EE->Draw();
  c_test_nonprompt_EE->Print("plots/pT_PU200_nonprompt_EE_2sigma_test.pdf");


}


void N_genMatched() {
  gStyle->SetOptStat(0);
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  // PU200
  f_PU200_prompt = new TFile(Form("data/%s", path_PU200_prompt.Data()));
  f_PU200_nonprompt = new TFile(Form("data/%s", path_PU200_nonprompt.Data()));
  // noPU
  f_noPU_prompt = new TFile(Form("data/%s", path_noPU_prompt.Data()));
  f_noPU_nonprompt = new TFile(Form("data/%s", path_noPU_nonprompt.Data()));

  TString dir="/DQMData/Run 1/MTD/Run summary/MuonIso/";
  TH1D *h_PU200_prompt, *h_PU200_nonprompt, *h_noPU_prompt, *h_noPU_nonprompt;
  TH1D *h_PU200_prompt_EB, *h_PU200_nonprompt_EB, *h_noPU_prompt_EB, *h_noPU_nonprompt_EB;
  TH1D *h_PU200_prompt_EE, *h_PU200_nonprompt_EE, *h_noPU_prompt_EE, *h_noPU_nonprompt_EE;
/*
  h_PU200_prompt 	= (TH1D*)f_PU200_prompt->Get(dir+"Track_genMatch_info_check");
  h_PU200_nonprompt 	= (TH1D*)f_PU200_nonprompt->Get(dir+"Track_genMatch_info_check");
  h_noPU_prompt 	= (TH1D*)f_noPU_prompt->Get(dir+"Track_genMatch_info_check");
  h_noPU_nonprompt 	= (TH1D*)f_noPU_nonprompt->Get(dir+"Track_genMatch_info_check");
*/
  h_PU200_prompt 	= (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_track_genMatched_Sig");
  h_PU200_nonprompt 	= (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_track_genMatched_Bkg");
  h_noPU_prompt 	= (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_track_genMatched_Sig");
  h_noPU_nonprompt 	= (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_track_genMatched_Bkg");

  h_PU200_prompt_EB 	= (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_track_genMatched_Sig_EB");
  h_PU200_nonprompt_EB 	= (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_track_genMatched_Bkg_EB");
  h_noPU_prompt_EB 	= (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_track_genMatched_Sig_EB");
  h_noPU_nonprompt_EB 	= (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_track_genMatched_Bkg_EB");
  h_PU200_prompt_EE 	= (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_track_genMatched_Sig_EE");
  h_PU200_nonprompt_EE 	= (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_track_genMatched_Bkg_EE");
  h_noPU_prompt_EE 	= (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_track_genMatched_Sig_EE");
  h_noPU_nonprompt_EE 	= (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_track_genMatched_Bkg_EE");
  
  ///////////////
  // Cosmetics //
  ///////////////
  // Barrel region
  h_PU200_prompt->SetTitle("Check whether tracks are matched with GenParticles"); h_PU200_nonprompt->SetTitle("Check whether tracks are matched with GenParticles"); h_noPU_prompt->SetTitle("Check whether track are matched with GenParticles"); h_noPU_nonprompt->SetTitle("Check whether tracks are matched with a GenParticles");
  h_PU200_prompt->GetXaxis()->SetTitle(""); h_PU200_nonprompt->GetXaxis()->SetTitle(""); h_noPU_prompt->GetXaxis()->SetTitle(""); h_noPU_nonprompt->GetXaxis()->SetTitle("");
  h_PU200_prompt->GetYaxis()->SetTitle("Counts"); h_PU200_nonprompt->GetYaxis()->SetTitle("Counts"); h_noPU_prompt->GetYaxis()->SetTitle("Counts"); h_noPU_nonprompt->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt->SetLineWidth(2); h_PU200_nonprompt->SetLineWidth(2); h_noPU_prompt->SetLineWidth(2); h_noPU_nonprompt->SetLineWidth(2);
  h_PU200_prompt->SetLineColor(kBlack); h_PU200_nonprompt->SetLineColor(kBlack); h_noPU_prompt->SetLineColor(kBlack); h_noPU_nonprompt->SetLineColor(kBlack);
  h_PU200_prompt->SetNdivisions(3); h_PU200_nonprompt->SetNdivisions(3); h_noPU_prompt->SetNdivisions(3); h_noPU_nonprompt->SetNdivisions(3);
  h_PU200_prompt->SetMinimum(0.); h_PU200_nonprompt->SetMinimum(0.); h_noPU_prompt->SetMinimum(0.); h_noPU_nonprompt->SetMinimum(0.);
  // Barrel region
  h_PU200_prompt_EB->SetTitle("Check whether tracks are matched with GenParticles"); h_PU200_nonprompt_EB->SetTitle("Check whether tracks are matched with GenParticles"); h_noPU_prompt_EB->SetTitle("Check whether track are matched with GenParticles"); h_noPU_nonprompt_EB->SetTitle("Check whether tracks are matched with a GenParticles");
  h_PU200_prompt_EB->GetXaxis()->SetTitle(""); h_PU200_nonprompt_EB->GetXaxis()->SetTitle(""); h_noPU_prompt_EB->GetXaxis()->SetTitle(""); h_noPU_nonprompt_EB->GetXaxis()->SetTitle("");
  h_PU200_prompt_EB->GetYaxis()->SetTitle("Counts"); h_PU200_nonprompt_EB->GetYaxis()->SetTitle("Counts"); h_noPU_prompt_EB->GetYaxis()->SetTitle("Counts"); h_noPU_nonprompt_EB->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt_EB->SetLineWidth(2); h_PU200_nonprompt_EB->SetLineWidth(2); h_noPU_prompt_EB->SetLineWidth(2); h_noPU_nonprompt_EB->SetLineWidth(2);
  h_PU200_prompt_EB->SetLineColor(kBlack); h_PU200_nonprompt_EB->SetLineColor(kBlack); h_noPU_prompt_EB->SetLineColor(kBlack); h_noPU_nonprompt_EB->SetLineColor(kBlack);
  h_PU200_prompt_EB->SetNdivisions(3); h_PU200_nonprompt_EB->SetNdivisions(3); h_noPU_prompt_EB->SetNdivisions(3); h_noPU_nonprompt_EB->SetNdivisions(3);
  h_PU200_prompt_EB->SetMinimum(0.); h_PU200_nonprompt_EB->SetMinimum(0.); h_noPU_prompt_EB->SetMinimum(0.); h_noPU_nonprompt_EB->SetMinimum(0.);
  // Endcap region
  h_PU200_prompt_EE->SetTitle("Check whether tracks are matched with GenParticles"); h_PU200_nonprompt_EE->SetTitle("Check whether tracks are matched with GenParticles"); h_noPU_prompt_EE->SetTitle("Check whether track are matched with GenParticles"); h_noPU_nonprompt_EE->SetTitle("Check whether tracks are matched with a GenParticles");
  h_PU200_prompt_EE->GetXaxis()->SetTitle(""); h_PU200_nonprompt_EE->GetXaxis()->SetTitle(""); h_noPU_prompt_EE->GetXaxis()->SetTitle(""); h_noPU_nonprompt_EE->GetXaxis()->SetTitle("");
  h_PU200_prompt_EE->GetYaxis()->SetTitle("Counts"); h_PU200_nonprompt_EE->GetYaxis()->SetTitle("Counts"); h_noPU_prompt_EE->GetYaxis()->SetTitle("Counts"); h_noPU_nonprompt_EE->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt_EE->SetLineWidth(2); h_PU200_nonprompt_EE->SetLineWidth(2); h_noPU_prompt_EE->SetLineWidth(2); h_noPU_nonprompt_EE->SetLineWidth(2);
  h_PU200_prompt_EE->SetLineColor(kBlack); h_PU200_nonprompt_EE->SetLineColor(kBlack); h_noPU_prompt_EE->SetLineColor(kBlack); h_noPU_nonprompt_EE->SetLineColor(kBlack);
  h_PU200_prompt_EE->SetNdivisions(3); h_PU200_nonprompt_EE->SetNdivisions(3); h_noPU_prompt_EE->SetNdivisions(3); h_noPU_nonprompt_EE->SetNdivisions(3);
  h_PU200_prompt_EE->SetMinimum(0.); h_PU200_nonprompt_EE->SetMinimum(0.); h_noPU_prompt_EE->SetMinimum(0.); h_noPU_nonprompt_EE->SetMinimum(0.);


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
  // EB
  TCanvas* c_PU200_prompt_EB = new TCanvas("c_PU200_prompt_EB", "c_PU200_prompt_EB", 1500, 1500);
  c_PU200_prompt_EB->cd();
  h_PU200_prompt_EB->Draw("hist");
  c_PU200_prompt_EB->SetLeftMargin(0.12);
  c_PU200_prompt_EB->Print("plots/N_genMatched_PU200_prompt_EB.pdf");

  TCanvas* c_PU200_nonprompt_EB = new TCanvas("c_PU200_nonprompt_EB", "c_PU200_nonprompt_EB", 1500, 1500);
  c_PU200_nonprompt_EB->cd();
  h_PU200_nonprompt_EB->Draw("hist");
  c_PU200_nonprompt_EB->SetLeftMargin(0.12);
  c_PU200_nonprompt_EB->Print("plots/N_genMatched_PU200_nonprompt_EB.pdf");
  // EE
  TCanvas* c_PU200_prompt_EE = new TCanvas("c_PU200_prompt_EE", "c_PU200_prompt_EE", 1500, 1500);
  c_PU200_prompt_EE->cd();
  h_PU200_prompt_EE->Draw("hist");
  c_PU200_prompt_EE->SetLeftMargin(0.12);
  c_PU200_prompt_EE->Print("plots/N_genMatched_PU200_prompt_EE.pdf");

  TCanvas* c_PU200_nonprompt_EE = new TCanvas("c_PU200_nonprompt_EE", "c_PU200_nonprompt_EE", 1500, 1500);
  c_PU200_nonprompt_EE->cd();
  h_PU200_nonprompt_EE->Draw("hist");
  c_PU200_nonprompt_EE->SetLeftMargin(0.12);
  c_PU200_nonprompt_EE->Print("plots/N_genMatched_PU200_nonprompt_EE.pdf");


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
   // EB
  TCanvas* c_noPU_prompt_EB = new TCanvas("c_noPU_prompt_EB", "c_noPU_prompt_EB", 1500, 1500);
  c_noPU_prompt_EB->cd();
  h_noPU_prompt_EB->Draw("hist");
  c_noPU_prompt_EB->SetLeftMargin(0.12);
  c_noPU_prompt_EB->Print("plots/N_genMatched_noPU_prompt_EB.pdf");

  TCanvas* c_noPU_nonprompt_EB = new TCanvas("c_noPU_nonprompt_EB", "c_noPU_nonprompt_EB", 1500, 1500);
  c_noPU_nonprompt_EB->cd();
  h_noPU_nonprompt_EB->Draw("hist");
  c_noPU_nonprompt_EB->SetLeftMargin(0.12);
  c_noPU_nonprompt_EB->Print("plots/N_genMatched_noPU_nonprompt_EB.pdf");
  // EE
  TCanvas* c_noPU_prompt_EE = new TCanvas("c_noPU_prompt_EE", "c_noPU_prompt_EE", 1500, 1500);
  c_noPU_prompt_EE->cd();
  h_noPU_prompt_EE->Draw("hist");
  c_noPU_prompt_EE->SetLeftMargin(0.12);
  c_noPU_prompt_EE->Print("plots/N_genMatched_noPU_prompt_EE.pdf");

  TCanvas* c_noPU_nonprompt_EE = new TCanvas("c_noPU_nonprompt_EE", "c_noPU_nonprompt_EE", 1500, 1500);
  c_noPU_nonprompt_EE->cd();
  h_noPU_nonprompt_EE->Draw("hist");
  c_noPU_nonprompt_EE->SetLeftMargin(0.12);
  c_noPU_nonprompt_EE->Print("plots/N_genMatched_noPU_nonprompt_EE.pdf");
  c_noPU_nonprompt->Print("plots/N_genMatched_noPU_nonprompt.pdf");

}

void track_type_v1() {
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  TFile *f_PU200_prompt_vtx, *f_PU200_nonprompt_vtx, *f_noPU_prompt_vtx, *f_noPU_nonprompt_vtx;
  // PU200
  f_PU200_prompt = new TFile(Form("data/%s", path_PU200_prompt.Data()));
  f_PU200_nonprompt = new TFile(Form("data/%s", path_PU200_nonprompt.Data()));
  f_PU200_prompt_vtx = new TFile(Form("data/%s", path_PU200_prompt_vtx.Data()));
  f_PU200_nonprompt_vtx = new TFile(Form("data/%s", path_PU200_nonprompt_vtx.Data()));
  // noPU
  f_noPU_prompt = new TFile(Form("data/%s", path_noPU_prompt.Data()));
  f_noPU_nonprompt = new TFile(Form("data/%s", path_noPU_nonprompt.Data()));
  f_noPU_prompt_vtx = new TFile(Form("data/%s", path_noPU_prompt_vtx.Data()));
  f_noPU_nonprompt_vtx = new TFile(Form("data/%s", path_noPU_nonprompt_vtx.Data()));

  TString dir="/DQMData/Run 1/MTD/Run summary/MuonIso/";
  TH1D *h_PU200_prompt_type_v1_EB, *h_PU200_prompt_type_v1_EE, *h_PU200_nonprompt_type_v1_EB, *h_PU200_nonprompt_type_v1_EE;
  TH1D *h_noPU_prompt_type_v1_EB, *h_noPU_prompt_type_v1_EE, *h_noPU_nonprompt_type_v1_EB, *h_noPU_nonprompt_type_v1_EE;
    //vtx
  TH1D *h_PU200_prompt_type_v1_EB_vtx, *h_PU200_prompt_type_v1_EE_vtx, *h_PU200_nonprompt_type_v1_EB_vtx, *h_PU200_nonprompt_type_v1_EE_vtx;
  TH1D *h_noPU_prompt_type_v1_EB_vtx, *h_noPU_prompt_type_v1_EE_vtx, *h_noPU_nonprompt_type_v1_EB_vtx, *h_noPU_nonprompt_type_v1_EE_vtx;

  // Barrel
    // PU200
  h_PU200_prompt_type_v1_EB    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_track_type_Sig_EB");
  h_PU200_nonprompt_type_v1_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_track_type_Bkg_EB");
  h_PU200_prompt_type_v1_EB_vtx    = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_Iso_track_type_Sig_EB");
  h_PU200_nonprompt_type_v1_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_Iso_track_type_Bkg_EB");
    // noPU
  h_noPU_prompt_type_v1_EB    = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_track_type_Sig_EB");
  h_noPU_nonprompt_type_v1_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_track_type_Bkg_EB");
  h_noPU_prompt_type_v1_EB_vtx    = (TH1D*)f_noPU_prompt_vtx->Get(dir+"Muon_Iso_track_type_Sig_EB");
  h_noPU_nonprompt_type_v1_EB_vtx = (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"Muon_Iso_track_type_Bkg_EB");

  // Endcap
    // PU200
  h_PU200_prompt_type_v1_EE    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_track_type_Sig_EE");
  h_PU200_nonprompt_type_v1_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_track_type_Bkg_EE");
  h_PU200_prompt_type_v1_EE_vtx    = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_Iso_track_type_Sig_EE");
  h_PU200_nonprompt_type_v1_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_Iso_track_type_Bkg_EE");
    // noPU
  h_noPU_prompt_type_v1_EE    = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_track_type_Sig_EE");
  h_noPU_nonprompt_type_v1_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_track_type_Bkg_EE");
  h_noPU_prompt_type_v1_EE_vtx    = (TH1D*)f_noPU_prompt_vtx->Get(dir+"Muon_Iso_track_type_Sig_EE");
  h_noPU_nonprompt_type_v1_EE_vtx = (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"Muon_Iso_track_type_Bkg_EE");


  ////////////////
  // Draw plots //
  ////////////////
  // PU200
    // prompt
      // Barrel
  TCanvas* c_PU200_prompt_type_v1_EB = new TCanvas("c_PU200_prompt_type_v1_EB", "c_PU200_prompt_type_v1_EB", 1500, 1500);
  c_PU200_prompt_type_v1_EB->cd();
  c_PU200_prompt_type_v1_EB->SetLeftMargin(0.12);
  TH1D* h_PU200_prompt_type_v1_EB_incl = new TH1D("h_PU200_prompt_type_v1_EB_incl","h_PU200_prompt_type_v1_EB_incl",4,0,4);
  h_PU200_prompt_type_v1_EB_incl->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt_type_v1_EB_incl->SetTitle("");
  h_PU200_prompt_type_v1_EB_incl->SetMinimum(0.);
  h_PU200_prompt_type_v1_EB_incl->SetBinContent(1,h_PU200_prompt_type_v1_EB->GetBinContent(1));
  h_PU200_prompt_type_v1_EB_incl->SetBinContent(2,h_PU200_prompt_type_v1_EB->GetBinContent(2));
  h_PU200_prompt_type_v1_EB_incl->SetBinContent(3,h_PU200_prompt_type_v1_EB->GetBinContent(3));
  h_PU200_prompt_type_v1_EB_incl->SetBinContent(4,h_PU200_prompt_type_v1_EB->GetBinContent(4));
  h_PU200_prompt_type_v1_EB_incl->Draw("hist");
  c_PU200_prompt_type_v1_EB->Print("plots/trktype_v1_PU200_prompt_EB.pdf");
  cout << "PU200" << endl;
  cout << "Barrel prompt    |" << h_PU200_prompt_type_v1_EB->GetBinContent(1) << " : " << h_PU200_prompt_type_v1_EB->GetBinContent(2) << " : " << h_PU200_prompt_type_v1_EB->GetBinContent(3) << " : " << h_PU200_prompt_type_v1_EB->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_PU200_prompt_type_v1_EE = new TCanvas("c_PU200_prompt_type_v1_EE", "c_PU200_prompt_type_v1_EE", 1500, 1500);
  c_PU200_prompt_type_v1_EE->cd();
  c_PU200_prompt_type_v1_EE->SetLeftMargin(0.12);
  TH1D* h_PU200_prompt_type_v1_EE_incl = new TH1D("h_PU200_prompt_type_v1_EE_incl","h_PU200_prompt_type_v1_EE_incl",4,0,4);
  h_PU200_prompt_type_v1_EE_incl->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt_type_v1_EE_incl->SetTitle("");
  h_PU200_prompt_type_v1_EE_incl->SetMinimum(0.);
  h_PU200_prompt_type_v1_EE_incl->SetBinContent(1,h_PU200_prompt_type_v1_EE->GetBinContent(1));
  h_PU200_prompt_type_v1_EE_incl->SetBinContent(2,h_PU200_prompt_type_v1_EE->GetBinContent(2));
  h_PU200_prompt_type_v1_EE_incl->SetBinContent(3,h_PU200_prompt_type_v1_EE->GetBinContent(3));
  h_PU200_prompt_type_v1_EE_incl->SetBinContent(4,h_PU200_prompt_type_v1_EE->GetBinContent(4));
  h_PU200_prompt_type_v1_EE_incl->Draw("hist");
  c_PU200_prompt_type_v1_EE->Print("plots/trktype_v1_PU200_prompt_EE.pdf");
  cout << "Endcap prompt    |" << h_PU200_prompt_type_v1_EE->GetBinContent(1) << " : " << h_PU200_prompt_type_v1_EE->GetBinContent(2) << " : " << h_PU200_prompt_type_v1_EE->GetBinContent(3) << " : " << h_PU200_prompt_type_v1_EE->GetBinContent(4) << endl;
    // nonprompt
      // Barrel
  TCanvas* c_PU200_nonprompt_type_v1_EB = new TCanvas("c_PU200_nonprompt_type_v1_EB", "c_PU200_nonprompt_type_v1_EB", 1500, 1500);
  c_PU200_nonprompt_type_v1_EB->cd();
  c_PU200_nonprompt_type_v1_EB->SetLeftMargin(0.12);
  TH1D* h_PU200_nonprompt_type_v1_EB_incl = new TH1D("h_PU200_nonprompt_type_v1_EB_incl","h_PU200_nonprompt_type_v1_EB_incl",4,0,4);
  h_PU200_nonprompt_type_v1_EB_incl->GetYaxis()->SetTitle("Counts");
  h_PU200_nonprompt_type_v1_EB_incl->SetTitle("");
  h_PU200_nonprompt_type_v1_EB_incl->SetMinimum(0.);
  h_PU200_nonprompt_type_v1_EB_incl->SetBinContent(1,h_PU200_nonprompt_type_v1_EB->GetBinContent(1));
  h_PU200_nonprompt_type_v1_EB_incl->SetBinContent(2,h_PU200_nonprompt_type_v1_EB->GetBinContent(2));
  h_PU200_nonprompt_type_v1_EB_incl->SetBinContent(3,h_PU200_nonprompt_type_v1_EB->GetBinContent(3));
  h_PU200_nonprompt_type_v1_EB_incl->SetBinContent(4,h_PU200_nonprompt_type_v1_EB->GetBinContent(4));
  h_PU200_nonprompt_type_v1_EB_incl->Draw("hist");
  c_PU200_nonprompt_type_v1_EB->Print("plots/trktype_v1_PU200_nonprompt_EB.pdf");
  cout << "Barrel nonprompt |" << h_PU200_nonprompt_type_v1_EB->GetBinContent(1) << " : " << h_PU200_nonprompt_type_v1_EB->GetBinContent(2) << " : " << h_PU200_nonprompt_type_v1_EB->GetBinContent(3) << " : " << h_PU200_nonprompt_type_v1_EB->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_PU200_nonprompt_type_v1_EE = new TCanvas("c_PU200_nonprompt_type_v1_EE", "c_PU200_nonprompt_type_v1_EE", 1500, 1500);
  c_PU200_nonprompt_type_v1_EE->cd();
  c_PU200_nonprompt_type_v1_EE->SetLeftMargin(0.12);
  TH1D* h_PU200_nonprompt_type_v1_EE_incl = new TH1D("h_PU200_nonprompt_type_v1_EE_incl","h_PU200_nonprompt_type_v1_EE_incl",4,0,4);
  h_PU200_nonprompt_type_v1_EE_incl->GetYaxis()->SetTitle("Counts");
  h_PU200_nonprompt_type_v1_EE_incl->SetTitle("");
  h_PU200_nonprompt_type_v1_EE_incl->SetMinimum(0.);
  h_PU200_nonprompt_type_v1_EE_incl->SetBinContent(1,h_PU200_nonprompt_type_v1_EE->GetBinContent(1));
  h_PU200_nonprompt_type_v1_EE_incl->SetBinContent(2,h_PU200_nonprompt_type_v1_EE->GetBinContent(2));
  h_PU200_nonprompt_type_v1_EE_incl->SetBinContent(3,h_PU200_nonprompt_type_v1_EE->GetBinContent(3));
  h_PU200_nonprompt_type_v1_EE_incl->SetBinContent(4,h_PU200_nonprompt_type_v1_EE->GetBinContent(4));
  h_PU200_nonprompt_type_v1_EE_incl->Draw("hist");
  c_PU200_nonprompt_type_v1_EE->Print("plots/trktype_v1_PU200_nonprompt_EE.pdf");
  cout << "Endcap nonprompt |" << h_PU200_nonprompt_type_v1_EE->GetBinContent(1) << " : " << h_PU200_nonprompt_type_v1_EE->GetBinContent(2) << " : " << h_PU200_nonprompt_type_v1_EE->GetBinContent(3)  << " : " << h_PU200_nonprompt_type_v1_EE->GetBinContent(4)<< endl;
  // noPU
    // prompt
      // Barrel
  TCanvas* c_noPU_prompt_type_v1_EB = new TCanvas("c_noPU_prompt_type_v1_EB", "c_noPU_prompt_type_v1_EB", 1500, 1500);
  c_noPU_prompt_type_v1_EB->cd();
  c_noPU_prompt_type_v1_EB->SetLeftMargin(0.12);
  TH1D* h_noPU_prompt_type_v1_EB_incl = new TH1D("h_noPU_prompt_type_v1_EB_incl","h_noPU_prompt_type_v1_EB_incl",4,0,4);
  h_noPU_prompt_type_v1_EB_incl->GetYaxis()->SetTitle("Counts");
  h_noPU_prompt_type_v1_EB_incl->SetTitle("");
  h_noPU_prompt_type_v1_EB_incl->SetMinimum(0.);
  h_noPU_prompt_type_v1_EB_incl->SetBinContent(1,h_noPU_prompt_type_v1_EB->GetBinContent(1));
  h_noPU_prompt_type_v1_EB_incl->SetBinContent(2,h_noPU_prompt_type_v1_EB->GetBinContent(2));
  h_noPU_prompt_type_v1_EB_incl->SetBinContent(3,h_noPU_prompt_type_v1_EB->GetBinContent(3));
  h_noPU_prompt_type_v1_EB_incl->SetBinContent(4,h_noPU_prompt_type_v1_EB->GetBinContent(4));
  h_noPU_prompt_type_v1_EB_incl->Draw("hist");
  c_noPU_prompt_type_v1_EB->Print("plots/trktype_v1_noPU_prompt_EB.pdf");
  cout << endl;
  cout << "noPU" << endl;
  cout << "Barrel prompt    |" << h_noPU_prompt_type_v1_EB->GetBinContent(1) << " : " << h_noPU_prompt_type_v1_EB->GetBinContent(2) << " : " << h_noPU_prompt_type_v1_EB->GetBinContent(3) << " : " << h_noPU_prompt_type_v1_EB->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_noPU_prompt_type_v1_EE = new TCanvas("c_noPU_prompt_type_v1_EE", "c_noPU_prompt_type_v1_EE", 1500, 1500);
  c_noPU_prompt_type_v1_EE->cd();
  c_noPU_prompt_type_v1_EE->SetLeftMargin(0.12);
  TH1D* h_noPU_prompt_type_v1_EE_incl = new TH1D("h_noPU_prompt_type_v1_EE_incl","h_noPU_prompt_type_v1_EE_incl",4,0,4);
  h_noPU_prompt_type_v1_EE_incl->GetYaxis()->SetTitle("Counts");
  h_noPU_prompt_type_v1_EE_incl->SetTitle("");
  h_noPU_prompt_type_v1_EE_incl->SetMinimum(0.);
  h_noPU_prompt_type_v1_EE_incl->SetBinContent(1,h_noPU_prompt_type_v1_EE->GetBinContent(1));
  h_noPU_prompt_type_v1_EE_incl->SetBinContent(2,h_noPU_prompt_type_v1_EE->GetBinContent(2));
  h_noPU_prompt_type_v1_EE_incl->SetBinContent(3,h_noPU_prompt_type_v1_EE->GetBinContent(3));
  h_noPU_prompt_type_v1_EE_incl->SetBinContent(4,h_noPU_prompt_type_v1_EE->GetBinContent(4));
  h_noPU_prompt_type_v1_EE_incl->Draw("hist");
  c_noPU_prompt_type_v1_EE->Print("plots/trktype_v1_noPU_prompt_EE.pdf");
  cout << "Endcap prompt    |" << h_noPU_prompt_type_v1_EE->GetBinContent(1) << " : " << h_noPU_prompt_type_v1_EE->GetBinContent(2) << " : " << h_noPU_prompt_type_v1_EE->GetBinContent(3) << " : " << h_noPU_prompt_type_v1_EE->GetBinContent(4) << endl;
    // nonprompt
      // Barrel
  TCanvas* c_noPU_nonprompt_type_v1_EB = new TCanvas("c_noPU_nonprompt_type_v1_EB", "c_noPU_nonprompt_type_v1_EB", 1500, 1500);
  c_noPU_nonprompt_type_v1_EB->cd();
  c_noPU_nonprompt_type_v1_EB->SetLeftMargin(0.12);
  TH1D* h_noPU_nonprompt_type_v1_EB_incl = new TH1D("h_noPU_nonprompt_type_v1_EB_incl","h_noPU_nonprompt_type_v1_EB_incl",4,0,4);
  h_noPU_nonprompt_type_v1_EB_incl->GetYaxis()->SetTitle("Counts");
  h_noPU_nonprompt_type_v1_EB_incl->SetTitle("");
  h_noPU_nonprompt_type_v1_EB_incl->SetMinimum(0.);
  h_noPU_nonprompt_type_v1_EB_incl->SetBinContent(1,h_noPU_nonprompt_type_v1_EB->GetBinContent(1));
  h_noPU_nonprompt_type_v1_EB_incl->SetBinContent(2,h_noPU_nonprompt_type_v1_EB->GetBinContent(2));
  h_noPU_nonprompt_type_v1_EB_incl->SetBinContent(3,h_noPU_nonprompt_type_v1_EB->GetBinContent(3));
  h_noPU_nonprompt_type_v1_EB_incl->SetBinContent(4,h_noPU_nonprompt_type_v1_EB->GetBinContent(4));
  h_noPU_nonprompt_type_v1_EB_incl->Draw("hist");
  c_noPU_nonprompt_type_v1_EB->Print("plots/trktype_v1_noPU_nonprompt_EB.pdf");
  cout << "Barrel nonprompt |" << h_noPU_nonprompt_type_v1_EB->GetBinContent(1) << " : " << h_noPU_nonprompt_type_v1_EB->GetBinContent(2) << " : " << h_noPU_nonprompt_type_v1_EB->GetBinContent(3) << " : " << h_noPU_nonprompt_type_v1_EB->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_noPU_nonprompt_type_v1_EE = new TCanvas("c_noPU_nonprompt_type_v1_EE", "c_noPU_nonprompt_type_v1_EE", 1500, 1500);
  c_noPU_nonprompt_type_v1_EE->cd();
  c_noPU_nonprompt_type_v1_EE->SetLeftMargin(0.12);
  TH1D* h_noPU_nonprompt_type_v1_EE_incl = new TH1D("h_noPU_nonprompt_type_v1_EE_incl","h_noPU_nonprompt_type_v1_EE_incl",4,0,4);
  h_noPU_nonprompt_type_v1_EE_incl->GetYaxis()->SetTitle("Counts");
  h_noPU_nonprompt_type_v1_EE_incl->SetTitle("");
  h_noPU_nonprompt_type_v1_EE_incl->SetMinimum(0.);
  h_noPU_nonprompt_type_v1_EE_incl->SetBinContent(1,h_noPU_nonprompt_type_v1_EE->GetBinContent(1));
  h_noPU_nonprompt_type_v1_EE_incl->SetBinContent(2,h_noPU_nonprompt_type_v1_EE->GetBinContent(2));
  h_noPU_nonprompt_type_v1_EE_incl->SetBinContent(3,h_noPU_nonprompt_type_v1_EE->GetBinContent(3));
  h_noPU_nonprompt_type_v1_EE_incl->SetBinContent(4,h_noPU_nonprompt_type_v1_EE->GetBinContent(4));
  h_noPU_nonprompt_type_v1_EE_incl->Draw("hist");
  c_noPU_nonprompt_type_v1_EE->Print("plots/trktype_v1_noPU_nonprompt_EE.pdf");
  cout << "Endcap nonprompt |" << h_noPU_nonprompt_type_v1_EE->GetBinContent(1) << " : " << h_noPU_nonprompt_type_v1_EE->GetBinContent(2) << " : " << h_noPU_nonprompt_type_v1_EE->GetBinContent(3) << " : " << h_noPU_nonprompt_type_v1_EE->GetBinContent(4) << endl;
  cout << endl;


  ///////////
  /// vtx ///
  ///////////
  // PU200
    // prompt
      // Barrel
  TCanvas* c_PU200_prompt_type_v1_EB_vtx = new TCanvas("c_PU200_prompt_type_v1_EB_vtx", "c_PU200_prompt_type_v1_EB_vtx", 1500, 1500);
  c_PU200_prompt_type_v1_EB_vtx->cd();
  c_PU200_prompt_type_v1_EB_vtx->SetLeftMargin(0.12);
  TH1D* h_PU200_prompt_type_v1_EB_incl_vtx = new TH1D("h_PU200_prompt_type_v1_EB_incl_vtx","h_PU200_prompt_type_v1_EB_incl_vtx",4,0,4);
  h_PU200_prompt_type_v1_EB_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt_type_v1_EB_incl_vtx->SetTitle("");
  h_PU200_prompt_type_v1_EB_incl_vtx->SetMinimum(0.);
  h_PU200_prompt_type_v1_EB_incl_vtx->SetBinContent(1,h_PU200_prompt_type_v1_EB_vtx->GetBinContent(1));
  h_PU200_prompt_type_v1_EB_incl_vtx->SetBinContent(2,h_PU200_prompt_type_v1_EB_vtx->GetBinContent(2));
  h_PU200_prompt_type_v1_EB_incl_vtx->SetBinContent(3,h_PU200_prompt_type_v1_EB_vtx->GetBinContent(3));
  h_PU200_prompt_type_v1_EB_incl_vtx->SetBinContent(4,h_PU200_prompt_type_v1_EB_vtx->GetBinContent(4));
  h_PU200_prompt_type_v1_EB_incl_vtx->Draw("hist");
  c_PU200_prompt_type_v1_EB_vtx->Print("plots/trktype_v1_PU200_prompt_EB_vtx.pdf");
  cout << "PU200" << endl;
  cout << "Barrel prompt    |" << h_PU200_prompt_type_v1_EB_vtx->GetBinContent(1) << " : " << h_PU200_prompt_type_v1_EB_vtx->GetBinContent(2) << " : " << h_PU200_prompt_type_v1_EB_vtx->GetBinContent(3) << " : " << h_PU200_prompt_type_v1_EB_vtx->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_PU200_prompt_type_v1_EE_vtx = new TCanvas("c_PU200_prompt_type_v1_EE_vtx", "c_PU200_prompt_type_v1_EE_vtx", 1500, 1500);
  c_PU200_prompt_type_v1_EE_vtx->cd();
  c_PU200_prompt_type_v1_EE_vtx->SetLeftMargin(0.12);
  TH1D* h_PU200_prompt_type_v1_EE_incl_vtx = new TH1D("h_PU200_prompt_type_v1_EE_incl_vtx","h_PU200_prompt_type_v1_EE_incl_vtx",4,0,4);
  h_PU200_prompt_type_v1_EE_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt_type_v1_EE_incl_vtx->SetTitle("");
  h_PU200_prompt_type_v1_EE_incl_vtx->SetMinimum(0.);
  h_PU200_prompt_type_v1_EE_incl_vtx->SetBinContent(1,h_PU200_prompt_type_v1_EE_vtx->GetBinContent(1));
  h_PU200_prompt_type_v1_EE_incl_vtx->SetBinContent(2,h_PU200_prompt_type_v1_EE_vtx->GetBinContent(2));
  h_PU200_prompt_type_v1_EE_incl_vtx->SetBinContent(3,h_PU200_prompt_type_v1_EE_vtx->GetBinContent(3));
  h_PU200_prompt_type_v1_EE_incl_vtx->SetBinContent(4,h_PU200_prompt_type_v1_EE_vtx->GetBinContent(4));
  h_PU200_prompt_type_v1_EE_incl_vtx->Draw("hist");
  c_PU200_prompt_type_v1_EE_vtx->Print("plots/trktype_v1_PU200_prompt_EE_vtx.pdf");
  cout << "Endcap prompt    |" << h_PU200_prompt_type_v1_EE_vtx->GetBinContent(1) << " : " << h_PU200_prompt_type_v1_EE_vtx->GetBinContent(2) << " : " << h_PU200_prompt_type_v1_EE_vtx->GetBinContent(3) << " : " << h_PU200_prompt_type_v1_EE_vtx->GetBinContent(4) << endl;
    // nonprompt
      // Barrel
  TCanvas* c_PU200_nonprompt_type_v1_EB_vtx = new TCanvas("c_PU200_nonprompt_type_v1_EB_vtx", "c_PU200_nonprompt_type_v1_EB_vtx", 1500, 1500);
  c_PU200_nonprompt_type_v1_EB_vtx->cd();
  c_PU200_nonprompt_type_v1_EB_vtx->SetLeftMargin(0.12);
  TH1D* h_PU200_nonprompt_type_v1_EB_incl_vtx = new TH1D("h_PU200_nonprompt_type_v1_EB_incl_vtx","h_PU200_nonprompt_type_v1_EB_incl_vtx",4,0,4);
  h_PU200_nonprompt_type_v1_EB_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_PU200_nonprompt_type_v1_EB_incl_vtx->SetTitle("");
  h_PU200_nonprompt_type_v1_EB_incl_vtx->SetMinimum(0.);
  h_PU200_nonprompt_type_v1_EB_incl_vtx->SetBinContent(1,h_PU200_nonprompt_type_v1_EB_vtx->GetBinContent(1));
  h_PU200_nonprompt_type_v1_EB_incl_vtx->SetBinContent(2,h_PU200_nonprompt_type_v1_EB_vtx->GetBinContent(2));
  h_PU200_nonprompt_type_v1_EB_incl_vtx->SetBinContent(3,h_PU200_nonprompt_type_v1_EB_vtx->GetBinContent(3));
  h_PU200_nonprompt_type_v1_EB_incl_vtx->SetBinContent(4,h_PU200_nonprompt_type_v1_EB_vtx->GetBinContent(4));
  h_PU200_nonprompt_type_v1_EB_incl_vtx->Draw("hist");
  c_PU200_nonprompt_type_v1_EB_vtx->Print("plots/trktype_v1_PU200_nonprompt_EB_vtx.pdf");
  cout << "Barrel nonprompt |" << h_PU200_nonprompt_type_v1_EB_vtx->GetBinContent(1) << " : " << h_PU200_nonprompt_type_v1_EB_vtx->GetBinContent(2) << " : " << h_PU200_nonprompt_type_v1_EB_vtx->GetBinContent(3) << " : " << h_PU200_nonprompt_type_v1_EB_vtx->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_PU200_nonprompt_type_v1_EE_vtx = new TCanvas("c_PU200_nonprompt_type_v1_EE_vtx", "c_PU200_nonprompt_type_v1_EE_vtx", 1500, 1500);
  c_PU200_nonprompt_type_v1_EE_vtx->cd();
  c_PU200_nonprompt_type_v1_EE_vtx->SetLeftMargin(0.12);
  TH1D* h_PU200_nonprompt_type_v1_EE_incl_vtx = new TH1D("h_PU200_nonprompt_type_v1_EE_incl_vtx","h_PU200_nonprompt_type_v1_EE_incl_vtx",4,0,4);
  h_PU200_nonprompt_type_v1_EE_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_PU200_nonprompt_type_v1_EE_incl_vtx->SetTitle("");
  h_PU200_nonprompt_type_v1_EE_incl_vtx->SetMinimum(0.);
  h_PU200_nonprompt_type_v1_EE_incl_vtx->SetBinContent(1,h_PU200_nonprompt_type_v1_EE_vtx->GetBinContent(1));
  h_PU200_nonprompt_type_v1_EE_incl_vtx->SetBinContent(2,h_PU200_nonprompt_type_v1_EE_vtx->GetBinContent(2));
  h_PU200_nonprompt_type_v1_EE_incl_vtx->SetBinContent(3,h_PU200_nonprompt_type_v1_EE_vtx->GetBinContent(3));
  h_PU200_nonprompt_type_v1_EE_incl_vtx->SetBinContent(4,h_PU200_nonprompt_type_v1_EE_vtx->GetBinContent(4));
  h_PU200_nonprompt_type_v1_EE_incl_vtx->Draw("hist");
  c_PU200_nonprompt_type_v1_EE_vtx->Print("plots/trktype_v1_PU200_nonprompt_EE_vtx.pdf");
  cout << "Endcap nonprompt |" << h_PU200_nonprompt_type_v1_EE_vtx->GetBinContent(1) << " : " << h_PU200_nonprompt_type_v1_EE_vtx->GetBinContent(2) << " : " << h_PU200_nonprompt_type_v1_EE_vtx->GetBinContent(3) << " : " << h_PU200_nonprompt_type_v1_EE_vtx->GetBinContent(4) << endl;
  // noPU
    // prompt
      // Barrel
  TCanvas* c_noPU_prompt_type_v1_EB_vtx = new TCanvas("c_noPU_prompt_type_v1_EB_vtx", "c_noPU_prompt_type_v1_EB_vtx", 1500, 1500);
  c_noPU_prompt_type_v1_EB_vtx->cd();
  c_noPU_prompt_type_v1_EB_vtx->SetLeftMargin(0.12);
  TH1D* h_noPU_prompt_type_v1_EB_incl_vtx = new TH1D("h_noPU_prompt_type_v1_EB_incl_vtx","h_noPU_prompt_type_v1_EB_incl_vtx",4,0,4);
  h_noPU_prompt_type_v1_EB_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_noPU_prompt_type_v1_EB_incl_vtx->SetTitle("");
  h_noPU_prompt_type_v1_EB_incl_vtx->SetMinimum(0.);
  h_noPU_prompt_type_v1_EB_incl_vtx->SetBinContent(1,h_noPU_prompt_type_v1_EB_vtx->GetBinContent(1));
  h_noPU_prompt_type_v1_EB_incl_vtx->SetBinContent(2,h_noPU_prompt_type_v1_EB_vtx->GetBinContent(2));
  h_noPU_prompt_type_v1_EB_incl_vtx->SetBinContent(3,h_noPU_prompt_type_v1_EB_vtx->GetBinContent(3));
  h_noPU_prompt_type_v1_EB_incl_vtx->SetBinContent(4,h_noPU_prompt_type_v1_EB_vtx->GetBinContent(4));
  h_noPU_prompt_type_v1_EB_incl_vtx->Draw("hist");
  c_noPU_prompt_type_v1_EB_vtx->Print("plots/trktype_v1_noPU_prompt_EB_vtx.pdf");
  cout << endl;
  cout << "noPU" << endl;
  cout << "Barrel prompt    |" << h_noPU_prompt_type_v1_EB_vtx->GetBinContent(1) << " : " << h_noPU_prompt_type_v1_EB_vtx->GetBinContent(2) << " : " << h_noPU_prompt_type_v1_EB_vtx->GetBinContent(3) << " : " << h_noPU_prompt_type_v1_EB_vtx->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_noPU_prompt_type_v1_EE_vtx = new TCanvas("c_noPU_prompt_type_v1_EE_vtx", "c_noPU_prompt_type_v1_EE_vtx", 1500, 1500);
  c_noPU_prompt_type_v1_EE_vtx->cd();
  c_noPU_prompt_type_v1_EE_vtx->SetLeftMargin(0.12);
  TH1D* h_noPU_prompt_type_v1_EE_incl_vtx = new TH1D("h_noPU_prompt_type_v1_EE_incl_vtx","h_noPU_prompt_type_v1_EE_incl_vtx",4,0,4);
  h_noPU_prompt_type_v1_EE_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_noPU_prompt_type_v1_EE_incl_vtx->SetTitle("");
  h_noPU_prompt_type_v1_EE_incl_vtx->SetMinimum(0.);
  h_noPU_prompt_type_v1_EE_incl_vtx->SetBinContent(1,h_noPU_prompt_type_v1_EE_vtx->GetBinContent(1));
  h_noPU_prompt_type_v1_EE_incl_vtx->SetBinContent(2,h_noPU_prompt_type_v1_EE_vtx->GetBinContent(2));
  h_noPU_prompt_type_v1_EE_incl_vtx->SetBinContent(3,h_noPU_prompt_type_v1_EE_vtx->GetBinContent(3));
  h_noPU_prompt_type_v1_EE_incl_vtx->SetBinContent(4,h_noPU_prompt_type_v1_EE_vtx->GetBinContent(4));
  h_noPU_prompt_type_v1_EE_incl_vtx->Draw("hist");
  c_noPU_prompt_type_v1_EE_vtx->Print("plots/trktype_v1_noPU_prompt_EE_vtx.pdf");
  cout << "Endcap prompt    |" << h_noPU_prompt_type_v1_EE_vtx->GetBinContent(1) << " : " << h_noPU_prompt_type_v1_EE_vtx->GetBinContent(2) << " : " << h_noPU_prompt_type_v1_EE_vtx->GetBinContent(3) << " : " << h_noPU_prompt_type_v1_EE_vtx->GetBinContent(4) << endl;
    // nonprompt
      // Barrel
  TCanvas* c_noPU_nonprompt_type_v1_EB_vtx = new TCanvas("c_noPU_nonprompt_type_v1_EB_vtx", "c_noPU_nonprompt_type_v1_EB_vtx", 1500, 1500);
  c_noPU_nonprompt_type_v1_EB_vtx->cd();
  c_noPU_nonprompt_type_v1_EB_vtx->SetLeftMargin(0.12);
  TH1D* h_noPU_nonprompt_type_v1_EB_incl_vtx = new TH1D("h_noPU_nonprompt_type_v1_EB_incl_vtx","h_noPU_nonprompt_type_v1_EB_incl_vtx",4,0,4);
  h_noPU_nonprompt_type_v1_EB_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_noPU_nonprompt_type_v1_EB_incl_vtx->SetTitle("");
  h_noPU_nonprompt_type_v1_EB_incl_vtx->SetMinimum(0.);
  h_noPU_nonprompt_type_v1_EB_incl_vtx->SetBinContent(1,h_noPU_nonprompt_type_v1_EB_vtx->GetBinContent(1));
  h_noPU_nonprompt_type_v1_EB_incl_vtx->SetBinContent(2,h_noPU_nonprompt_type_v1_EB_vtx->GetBinContent(2));
  h_noPU_nonprompt_type_v1_EB_incl_vtx->SetBinContent(3,h_noPU_nonprompt_type_v1_EB_vtx->GetBinContent(3));
  h_noPU_nonprompt_type_v1_EB_incl_vtx->SetBinContent(4,h_noPU_nonprompt_type_v1_EB_vtx->GetBinContent(4));
  h_noPU_nonprompt_type_v1_EB_incl_vtx->Draw("hist");
  c_noPU_nonprompt_type_v1_EB_vtx->Print("plots/trktype_v1_noPU_nonprompt_EB_vtx.pdf");
  cout << "Barrel nonprompt |" << h_noPU_nonprompt_type_v1_EB_vtx->GetBinContent(1) << " : " << h_noPU_nonprompt_type_v1_EB_vtx->GetBinContent(2) << " : " << h_noPU_nonprompt_type_v1_EB_vtx->GetBinContent(3) << " : " << h_noPU_nonprompt_type_v1_EB_vtx->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_noPU_nonprompt_type_v1_EE_vtx = new TCanvas("c_noPU_nonprompt_type_v1_EE_vtx", "c_noPU_nonprompt_type_v1_EE_vtx", 1500, 1500);
  c_noPU_nonprompt_type_v1_EE_vtx->cd();
  c_noPU_nonprompt_type_v1_EE_vtx->SetLeftMargin(0.12);
  TH1D* h_noPU_nonprompt_type_v1_EE_incl_vtx = new TH1D("h_noPU_nonprompt_type_v1_EE_incl_vtx","h_noPU_nonprompt_type_v1_EE_incl_vtx",4,0,4);
  h_noPU_nonprompt_type_v1_EE_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_noPU_nonprompt_type_v1_EE_incl_vtx->SetTitle("");
  h_noPU_nonprompt_type_v1_EE_incl_vtx->SetMinimum(0.);
  h_noPU_nonprompt_type_v1_EE_incl_vtx->SetBinContent(1,h_noPU_nonprompt_type_v1_EE_vtx->GetBinContent(1));
  h_noPU_nonprompt_type_v1_EE_incl_vtx->SetBinContent(2,h_noPU_nonprompt_type_v1_EE_vtx->GetBinContent(2));
  h_noPU_nonprompt_type_v1_EE_incl_vtx->SetBinContent(3,h_noPU_nonprompt_type_v1_EE_vtx->GetBinContent(3));
  h_noPU_nonprompt_type_v1_EE_incl_vtx->SetBinContent(4,h_noPU_nonprompt_type_v1_EE_vtx->GetBinContent(4));
  h_noPU_nonprompt_type_v1_EE_incl_vtx->Draw("hist");
  c_noPU_nonprompt_type_v1_EE_vtx->Print("plots/trktype_v1_noPU_nonprompt_EE_vtx.pdf");
  cout << "Endcap nonprompt |" << h_noPU_nonprompt_type_v1_EE_vtx->GetBinContent(1) << " : " << h_noPU_nonprompt_type_v1_EE_vtx->GetBinContent(2) << " : " << h_noPU_nonprompt_type_v1_EE_vtx->GetBinContent(3) << " : " << h_noPU_nonprompt_type_v1_EE_vtx->GetBinContent(4) << endl;
  cout << endl;



}

void track_sigma_type_v1() {
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  TFile *f_PU200_prompt_vtx, *f_PU200_nonprompt_vtx, *f_noPU_prompt_vtx, *f_noPU_nonprompt_vtx;

  // PU200
  f_PU200_prompt = new TFile(Form("data/%s", path_PU200_prompt.Data()));
  f_PU200_nonprompt = new TFile(Form("data/%s", path_PU200_nonprompt.Data()));
  f_PU200_prompt_vtx = new TFile(Form("data/%s", path_PU200_prompt_vtx.Data()));
  f_PU200_nonprompt_vtx = new TFile(Form("data/%s", path_PU200_nonprompt_vtx.Data()));
  // noPU
  f_noPU_prompt = new TFile(Form("data/%s", path_noPU_prompt.Data()));
  f_noPU_nonprompt = new TFile(Form("data/%s", path_noPU_nonprompt.Data()));
  f_noPU_prompt_vtx = new TFile(Form("data/%s", path_noPU_prompt_vtx.Data()));
  f_noPU_nonprompt_vtx = new TFile(Form("data/%s", path_noPU_nonprompt_vtx.Data()));


  TString dir="/DQMData/Run 1/MTD/Run summary/MuonIso/";
    // (muon, track)
  TH1D *h_PU200_prompt_PV_EB, *h_PU200_prompt_PU_EB, *h_PU200_prompt_fake_EB, *h_PU200_prompt_no_tErr_EB;
  TH1D *h_PU200_prompt_PV_EE, *h_PU200_prompt_PU_EE, *h_PU200_prompt_fake_EE, *h_PU200_prompt_no_tErr_EE;
  TH1D *h_PU200_nonprompt_PV_EB, *h_PU200_nonprompt_PU_EB, *h_PU200_nonprompt_fake_EB, *h_PU200_nonprompt_no_tErr_EB;
  TH1D *h_PU200_nonprompt_PV_EE, *h_PU200_nonprompt_PU_EE, *h_PU200_nonprompt_fake_EE, *h_PU200_nonprompt_no_tErr_EE;
  TH1D *h_noPU_prompt_PV_EB, *h_noPU_prompt_PU_EB, *h_noPU_prompt_fake_EB, *h_noPU_prompt_no_tErr_EB;
  TH1D *h_noPU_prompt_PV_EE, *h_noPU_prompt_PU_EE, *h_noPU_prompt_fake_EE, *h_noPU_prompt_no_tErr_EE;
  TH1D *h_noPU_nonprompt_PV_EB, *h_noPU_nonprompt_PU_EB, *h_noPU_nonprompt_fake_EB, *h_noPU_nonprompt_no_tErr_EB;
  TH1D *h_noPU_nonprompt_PV_EE, *h_noPU_nonprompt_PU_EE, *h_noPU_nonprompt_fake_EE, *h_noPU_nonprompt_no_tErr_EE;
    // (PV, track)
  TH1D *h_PU200_prompt_PV_EB_vtx, *h_PU200_prompt_PU_EB_vtx, *h_PU200_prompt_fake_EB_vtx, *h_PU200_prompt_no_tErr_EB_vtx;
  TH1D *h_PU200_prompt_PV_EE_vtx, *h_PU200_prompt_PU_EE_vtx, *h_PU200_prompt_fake_EE_vtx, *h_PU200_prompt_no_tErr_EE_vtx;
  TH1D *h_PU200_nonprompt_PV_EB_vtx, *h_PU200_nonprompt_PU_EB_vtx, *h_PU200_nonprompt_fake_EB_vtx, *h_PU200_nonprompt_no_tErr_EB_vtx;
  TH1D *h_PU200_nonprompt_PV_EE_vtx, *h_PU200_nonprompt_PU_EE_vtx, *h_PU200_nonprompt_fake_EE_vtx, *h_PU200_nonprompt_no_tErr_EE_vtx;
  TH1D *h_noPU_prompt_PV_EB_vtx, *h_noPU_prompt_PU_EB_vtx, *h_noPU_prompt_fake_EB_vtx, *h_noPU_prompt_no_tErr_EB_vtx;
  TH1D *h_noPU_prompt_PV_EE_vtx, *h_noPU_prompt_PU_EE_vtx, *h_noPU_prompt_fake_EE_vtx, *h_noPU_prompt_no_tErr_EE_vtx;
  TH1D *h_noPU_nonprompt_PV_EB_vtx, *h_noPU_nonprompt_PU_EB_vtx, *h_noPU_nonprompt_fake_EB_vtx, *h_noPU_nonprompt_no_tErr_EB_vtx;
  TH1D *h_noPU_nonprompt_PV_EE_vtx, *h_noPU_nonprompt_PU_EE_vtx, *h_noPU_nonprompt_fake_EE_vtx, *h_noPU_nonprompt_no_tErr_EE_vtx;

  ///////////////////
  // (muon, track) //
  ///////////////////
  // Barrel
    // PU200
  h_PU200_prompt_PV_EB 	    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PVtrk_Sig_EB");
  h_PU200_prompt_PU_EB 	    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PUtrk_Sig_EB");
  h_PU200_prompt_fake_EB    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_faketrk_Sig_EB");
  h_PU200_prompt_no_tErr_EB = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_without_tErr_Sig_EB");
  h_PU200_nonprompt_PV_EB      = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PVtrk_Bkg_EB");
  h_PU200_nonprompt_PU_EB      = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PUtrk_Bkg_EB");
  h_PU200_nonprompt_fake_EB    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_faketrk_Bkg_EB");
  h_PU200_nonprompt_no_tErr_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_without_tErr_Bkg_EB");
    // noPU
  h_noPU_prompt_PV_EB 	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PVtrk_Sig_EB");
  h_noPU_prompt_PU_EB 	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PUtrk_Sig_EB");
  h_noPU_prompt_fake_EB    = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_faketrk_Sig_EB");
  h_noPU_prompt_no_tErr_EB = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_without_tErr_Sig_EB");
  h_noPU_nonprompt_PV_EB      = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PVtrk_Bkg_EB");
  h_noPU_nonprompt_PU_EB      = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PUtrk_Bkg_EB");
  h_noPU_nonprompt_fake_EB    = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_faketrk_Bkg_EB");
  h_noPU_nonprompt_no_tErr_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_without_tErr_Bkg_EB");

  // endcap
    // PU200
  h_PU200_prompt_PV_EE 	    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PVtrk_Sig_EE");
  h_PU200_prompt_PU_EE 	    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PUtrk_Sig_EE");
  h_PU200_prompt_fake_EE    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_faketrk_Sig_EE");
  h_PU200_prompt_no_tErr_EE = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_without_tErr_Sig_EE");
  h_PU200_nonprompt_PV_EE      = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PVtrk_Bkg_EE");
  h_PU200_nonprompt_PU_EE      = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PUtrk_Bkg_EE");
  h_PU200_nonprompt_fake_EE    = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_faketrk_Bkg_EE");
  h_PU200_nonprompt_no_tErr_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_without_tErr_Bkg_EE");
    // noPU
  h_noPU_prompt_PV_EE 	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PVtrk_Sig_EE");
  h_noPU_prompt_PU_EE 	   = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PUtrk_Sig_EE");
  h_noPU_prompt_fake_EE    = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_faketrk_Sig_EE");
  h_noPU_prompt_no_tErr_EE = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_without_tErr_Sig_EE");
  h_noPU_nonprompt_PV_EE      = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PVtrk_Bkg_EE");
  h_noPU_nonprompt_PU_EE      = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_PUtrk_Bkg_EE");
  h_noPU_nonprompt_fake_EE    = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_faketrk_Bkg_EE");
  h_noPU_nonprompt_no_tErr_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_dtSig_muon_track_reco_without_tErr_Bkg_EE");
  
  ///////////////////
  /// (PV, track) ///
  ///////////////////
  // Barrel
    // PU200
  h_PU200_prompt_PV_EB_vtx 	    = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PVtrk_Sig_EB");
  h_PU200_prompt_PU_EB_vtx 	    = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PUtrk_Sig_EB");
  h_PU200_prompt_fake_EB_vtx    = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_faketrk_Sig_EB");
  h_PU200_prompt_no_tErr_EB_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_without_tErr_Sig_EB");
  h_PU200_nonprompt_PV_EB_vtx      = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PVtrk_Bkg_EB");
  h_PU200_nonprompt_PU_EB_vtx      = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PUtrk_Bkg_EB");
  h_PU200_nonprompt_fake_EB_vtx    = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_faketrk_Bkg_EB");
  h_PU200_nonprompt_no_tErr_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_without_tErr_Bkg_EB");
    // noPU
  h_noPU_prompt_PV_EB_vtx 	   = (TH1D*)f_noPU_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PVtrk_Sig_EB");
  h_noPU_prompt_PU_EB_vtx 	   = (TH1D*)f_noPU_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PUtrk_Sig_EB");
  h_noPU_prompt_fake_EB_vtx    = (TH1D*)f_noPU_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_faketrk_Sig_EB");
  h_noPU_prompt_no_tErr_EB_vtx = (TH1D*)f_noPU_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_without_tErr_Sig_EB");
  h_noPU_nonprompt_PV_EB_vtx      = (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PVtrk_Bkg_EB");
  h_noPU_nonprompt_PU_EB_vtx      = (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PUtrk_Bkg_EB");
  h_noPU_nonprompt_fake_EB_vtx    = (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_faketrk_Bkg_EB");
  h_noPU_nonprompt_no_tErr_EB_vtx = (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_without_tErr_Bkg_EB");

  // endcap
    // PU200
  h_PU200_prompt_PV_EE_vtx 	    = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PVtrk_Sig_EE");
  h_PU200_prompt_PU_EE_vtx 	    = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PUtrk_Sig_EE");
  h_PU200_prompt_fake_EE_vtx    = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_faketrk_Sig_EE");
  h_PU200_prompt_no_tErr_EE_vtx = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_without_tErr_Sig_EE");
  h_PU200_nonprompt_PV_EE_vtx      = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PVtrk_Bkg_EE");
  h_PU200_nonprompt_PU_EE_vtx      = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PUtrk_Bkg_EE");
  h_PU200_nonprompt_fake_EE_vtx    = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_faketrk_Bkg_EE");
  h_PU200_nonprompt_no_tErr_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_without_tErr_Bkg_EE");
    // noPU
  h_noPU_prompt_PV_EE_vtx 	   = (TH1D*)f_noPU_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PVtrk_Sig_EE");
  h_noPU_prompt_PU_EE_vtx 	   = (TH1D*)f_noPU_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PUtrk_Sig_EE");
  h_noPU_prompt_fake_EE_vtx    = (TH1D*)f_noPU_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_faketrk_Sig_EE");
  h_noPU_prompt_no_tErr_EE_vtx = (TH1D*)f_noPU_prompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_without_tErr_Sig_EE");
  h_noPU_nonprompt_PV_EE_vtx      = (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PVtrk_Bkg_EE");
  h_noPU_nonprompt_PU_EE_vtx      = (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_PUtrk_Bkg_EE");
  h_noPU_nonprompt_fake_EE_vtx    = (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_faketrk_Bkg_EE");
  h_noPU_nonprompt_no_tErr_EE_vtx = (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"Muon_Iso_dtSig_PV_track_reco_without_tErr_Bkg_EE");


  ////////////////
  // Re-binning // for fake tracks (fake tracks do not have dtSig)
  ////////////////
  // (muon,track)
  TH1D* h_PU200_prompt_fake_EB_rebin = new TH1D("h_PU200_prompt_fake_EB_rebin","h_PU200_prompt_fake_EB_rebin",1,0,1); h_PU200_prompt_fake_EB_rebin->SetBinContent(1, h_PU200_prompt_fake_EB->GetBinContent(2));
  TH1D* h_PU200_prompt_fake_EE_rebin = new TH1D("h_PU200_prompt_fake_EE_rebin","h_PU200_prompt_fake_EE_rebin",1,0,1); h_PU200_prompt_fake_EE_rebin->SetBinContent(1, h_PU200_prompt_fake_EE->GetBinContent(2));
  TH1D* h_PU200_nonprompt_fake_EB_rebin = new TH1D("h_PU200_nonprompt_fake_EB_rebin","h_PU200_nonprompt_fake_EB_rebin",1,0,1); h_PU200_nonprompt_fake_EB_rebin->SetBinContent(1, h_PU200_nonprompt_fake_EB->GetBinContent(2));
  TH1D* h_PU200_nonprompt_fake_EE_rebin = new TH1D("h_PU200_nonprompt_fake_EE_rebin","h_PU200_nonprompt_fake_EE_rebin",1,0,1); h_PU200_nonprompt_fake_EE_rebin->SetBinContent(1, h_PU200_nonprompt_fake_EE->GetBinContent(2));
  TH1D* h_noPU_prompt_fake_EB_rebin = new TH1D("h_noPU_prompt_fake_EB_rebin","h_noPU_prompt_fake_EB_rebin",1,0,1); h_noPU_prompt_fake_EB_rebin->SetBinContent(1, h_noPU_prompt_fake_EB->GetBinContent(2));
  TH1D* h_noPU_prompt_fake_EE_rebin = new TH1D("h_noPU_prompt_fake_EE_rebin","h_noPU_prompt_fake_EE_rebin",1,0,1); h_noPU_prompt_fake_EE_rebin->SetBinContent(1, h_noPU_prompt_fake_EE->GetBinContent(2));
  TH1D* h_noPU_nonprompt_fake_EB_rebin = new TH1D("h_noPU_nonprompt_fake_EB_rebin","h_noPU_nonprompt_fake_EB_rebin",1,0,1); h_noPU_nonprompt_fake_EB_rebin->SetBinContent(1, h_noPU_nonprompt_fake_EB->GetBinContent(2));
  TH1D* h_noPU_nonprompt_fake_EE_rebin = new TH1D("h_noPU_nonprompt_fake_EE_rebin","h_noPU_nonprompt_fake_EE_rebin",1,0,1); h_noPU_nonprompt_fake_EE_rebin->SetBinContent(1, h_noPU_nonprompt_fake_EE->GetBinContent(2));

  TH1D* h_PU200_prompt_no_tErr_EB_rebin = new TH1D("h_PU200_prompt_no_tErr_EB_rebin","h_PU200_prompt_no_tErr_EB_rebin",1,0,1); h_PU200_prompt_no_tErr_EB_rebin->SetBinContent(1, h_PU200_prompt_no_tErr_EB->GetBinContent(2));
  TH1D* h_PU200_prompt_no_tErr_EE_rebin = new TH1D("h_PU200_prompt_no_tErr_EE_rebin","h_PU200_prompt_no_tErr_EE_rebin",1,0,1); h_PU200_prompt_no_tErr_EE_rebin->SetBinContent(1, h_PU200_prompt_no_tErr_EE->GetBinContent(2));
  TH1D* h_PU200_nonprompt_no_tErr_EB_rebin = new TH1D("h_PU200_nonprompt_no_tErr_EB_rebin","h_PU200_nonprompt_no_tErr_EB_rebin",1,0,1); h_PU200_nonprompt_no_tErr_EB_rebin->SetBinContent(1, h_PU200_nonprompt_no_tErr_EB->GetBinContent(2));
  TH1D* h_PU200_nonprompt_no_tErr_EE_rebin = new TH1D("h_PU200_nonprompt_no_tErr_EE_rebin","h_PU200_nonprompt_no_tErr_EE_rebin",1,0,1); h_PU200_nonprompt_no_tErr_EE_rebin->SetBinContent(1, h_PU200_nonprompt_no_tErr_EE->GetBinContent(2));
  TH1D* h_noPU_prompt_no_tErr_EB_rebin = new TH1D("h_noPU_prompt_no_tErr_EB_rebin","h_noPU_prompt_no_tErr_EB_rebin",1,0,1); h_noPU_prompt_no_tErr_EB_rebin->SetBinContent(1, h_noPU_prompt_no_tErr_EB->GetBinContent(2));
  TH1D* h_noPU_prompt_no_tErr_EE_rebin = new TH1D("h_noPU_prompt_no_tErr_EE_rebin","h_noPU_prompt_no_tErr_EE_rebin",1,0,1); h_noPU_prompt_no_tErr_EE_rebin->SetBinContent(1, h_noPU_prompt_no_tErr_EE->GetBinContent(2));
  TH1D* h_noPU_nonprompt_no_tErr_EB_rebin = new TH1D("h_noPU_nonprompt_no_tErr_EB_rebin","h_noPU_nonprompt_no_tErr_EB_rebin",1,0,1); h_noPU_nonprompt_no_tErr_EB_rebin->SetBinContent(1, h_noPU_nonprompt_no_tErr_EB->GetBinContent(2));
  TH1D* h_noPU_nonprompt_no_tErr_EE_rebin = new TH1D("h_noPU_nonprompt_no_tErr_EE_rebin","h_noPU_nonprompt_no_tErr_EE_rebin",1,0,1); h_noPU_nonprompt_no_tErr_EE_rebin->SetBinContent(1, h_noPU_nonprompt_no_tErr_EE->GetBinContent(2));
  // (PV,track)
  TH1D* h_PU200_prompt_fake_EB_rebin_vtx = new TH1D("h_PU200_prompt_fake_EB_rebin_vtx","h_PU200_prompt_fake_EB_rebin_vtx",1,0,1); h_PU200_prompt_fake_EB_rebin_vtx->SetBinContent(1, h_PU200_prompt_fake_EB_vtx->GetBinContent(2));
  TH1D* h_PU200_prompt_fake_EE_rebin_vtx = new TH1D("h_PU200_prompt_fake_EE_rebin_vtx","h_PU200_prompt_fake_EE_rebin_vtx",1,0,1); h_PU200_prompt_fake_EE_rebin_vtx->SetBinContent(1, h_PU200_prompt_fake_EE_vtx->GetBinContent(2));
  TH1D* h_PU200_nonprompt_fake_EB_rebin_vtx = new TH1D("h_PU200_nonprompt_fake_EB_rebin_vtx","h_PU200_nonprompt_fake_EB_rebin_vtx",1,0,1); h_PU200_nonprompt_fake_EB_rebin_vtx->SetBinContent(1, h_PU200_nonprompt_fake_EB_vtx->GetBinContent(2));
  TH1D* h_PU200_nonprompt_fake_EE_rebin_vtx = new TH1D("h_PU200_nonprompt_fake_EE_rebin_vtx","h_PU200_nonprompt_fake_EE_rebin_vtx",1,0,1); h_PU200_nonprompt_fake_EE_rebin_vtx->SetBinContent(1, h_PU200_nonprompt_fake_EE_vtx->GetBinContent(2));
  TH1D* h_noPU_prompt_fake_EB_rebin_vtx = new TH1D("h_noPU_prompt_fake_EB_rebin_vtx","h_noPU_prompt_fake_EB_rebin_vtx",1,0,1); h_noPU_prompt_fake_EB_rebin_vtx->SetBinContent(1, h_noPU_prompt_fake_EB_vtx->GetBinContent(2));
  TH1D* h_noPU_prompt_fake_EE_rebin_vtx = new TH1D("h_noPU_prompt_fake_EE_rebin_vtx","h_noPU_prompt_fake_EE_rebin_vtx",1,0,1); h_noPU_prompt_fake_EE_rebin_vtx->SetBinContent(1, h_noPU_prompt_fake_EE_vtx->GetBinContent(2));
  TH1D* h_noPU_nonprompt_fake_EB_rebin_vtx = new TH1D("h_noPU_nonprompt_fake_EB_rebin_vtx","h_noPU_nonprompt_fake_EB_rebin_vtx",1,0,1); h_noPU_nonprompt_fake_EB_rebin_vtx->SetBinContent(1, h_noPU_nonprompt_fake_EB_vtx->GetBinContent(2));
  TH1D* h_noPU_nonprompt_fake_EE_rebin_vtx = new TH1D("h_noPU_nonprompt_fake_EE_rebin_vtx","h_noPU_nonprompt_fake_EE_rebin_vtx",1,0,1); h_noPU_nonprompt_fake_EE_rebin_vtx->SetBinContent(1, h_noPU_nonprompt_fake_EE_vtx->GetBinContent(2));

  TH1D* h_PU200_prompt_no_tErr_EB_rebin_vtx = new TH1D("h_PU200_prompt_no_tErr_EB_rebin_vtx","h_PU200_prompt_no_tErr_EB_rebin_vtx",1,0,1); h_PU200_prompt_no_tErr_EB_rebin_vtx->SetBinContent(1, h_PU200_prompt_no_tErr_EB_vtx->GetBinContent(2));
  TH1D* h_PU200_prompt_no_tErr_EE_rebin_vtx = new TH1D("h_PU200_prompt_no_tErr_EE_rebin_vtx","h_PU200_prompt_no_tErr_EE_rebin_vtx",1,0,1); h_PU200_prompt_no_tErr_EE_rebin_vtx->SetBinContent(1, h_PU200_prompt_no_tErr_EE_vtx->GetBinContent(2));
  TH1D* h_PU200_nonprompt_no_tErr_EB_rebin_vtx = new TH1D("h_PU200_nonprompt_no_tErr_EB_rebin_vtx","h_PU200_nonprompt_no_tErr_EB_rebin_vtx",1,0,1); h_PU200_nonprompt_no_tErr_EB_rebin_vtx->SetBinContent(1, h_PU200_nonprompt_no_tErr_EB_vtx->GetBinContent(2));
  TH1D* h_PU200_nonprompt_no_tErr_EE_rebin_vtx = new TH1D("h_PU200_nonprompt_no_tErr_EE_rebin_vtx","h_PU200_nonprompt_no_tErr_EE_rebin_vtx",1,0,1); h_PU200_nonprompt_no_tErr_EE_rebin_vtx->SetBinContent(1, h_PU200_nonprompt_no_tErr_EE_vtx->GetBinContent(2));
  TH1D* h_noPU_prompt_no_tErr_EB_rebin_vtx = new TH1D("h_noPU_prompt_no_tErr_EB_rebin_vtx","h_noPU_prompt_no_tErr_EB_rebin_vtx",1,0,1); h_noPU_prompt_no_tErr_EB_rebin_vtx->SetBinContent(1, h_noPU_prompt_no_tErr_EB_vtx->GetBinContent(2));
  TH1D* h_noPU_prompt_no_tErr_EE_rebin_vtx = new TH1D("h_noPU_prompt_no_tErr_EE_rebin_vtx","h_noPU_prompt_no_tErr_EE_rebin_vtx",1,0,1); h_noPU_prompt_no_tErr_EE_rebin_vtx->SetBinContent(1, h_noPU_prompt_no_tErr_EE_vtx->GetBinContent(2));
  TH1D* h_noPU_nonprompt_no_tErr_EB_rebin_vtx = new TH1D("h_noPU_nonprompt_no_tErr_EB_rebin_vtx","h_noPU_nonprompt_no_tErr_EB_rebin_vtx",1,0,1); h_noPU_nonprompt_no_tErr_EB_rebin_vtx->SetBinContent(1, h_noPU_nonprompt_no_tErr_EB_vtx->GetBinContent(2));
  TH1D* h_noPU_nonprompt_no_tErr_EE_rebin_vtx = new TH1D("h_noPU_nonprompt_no_tErr_EE_rebin_vtx","h_noPU_nonprompt_no_tErr_EE_rebin_vtx",1,0,1); h_noPU_nonprompt_no_tErr_EE_rebin_vtx->SetBinContent(1, h_noPU_nonprompt_no_tErr_EE_vtx->GetBinContent(2));


  ///////////////
  // Cosmetics //
  ///////////////
  // (muon, track)
  h_PU200_prompt_PV_EB->SetFillColor(kGray+1); h_PU200_prompt_PU_EB->SetFillColor(kAzure+7); h_PU200_prompt_fake_EB_rebin->SetFillColor(kGreen+2); h_PU200_prompt_no_tErr_EB_rebin->SetFillColor(kYellow-7);
  h_PU200_prompt_PV_EE->SetFillColor(kGray+1); h_PU200_prompt_PU_EE->SetFillColor(kAzure+7); h_PU200_prompt_fake_EE_rebin->SetFillColor(kGreen+2); h_PU200_prompt_no_tErr_EE_rebin->SetFillColor(kYellow-7);
  h_PU200_nonprompt_PV_EB->SetFillColor(kGray+1); h_PU200_nonprompt_PU_EB->SetFillColor(kAzure+7); h_PU200_nonprompt_fake_EB_rebin->SetFillColor(kGreen+2); h_PU200_nonprompt_no_tErr_EB_rebin->SetFillColor(kYellow-7);
  h_PU200_nonprompt_PV_EE->SetFillColor(kGray+1); h_PU200_nonprompt_PU_EE->SetFillColor(kAzure+7); h_PU200_nonprompt_fake_EE_rebin->SetFillColor(kGreen+2); h_PU200_nonprompt_no_tErr_EE_rebin->SetFillColor(kYellow-7);
  h_noPU_prompt_PV_EB->SetFillColor(kGray+1); h_noPU_prompt_PU_EB->SetFillColor(kAzure+7); h_noPU_prompt_fake_EB_rebin->SetFillColor(kGreen+2); h_noPU_prompt_no_tErr_EB_rebin->SetFillColor(kYellow-7);
  h_noPU_prompt_PV_EE->SetFillColor(kGray+1); h_noPU_prompt_PU_EE->SetFillColor(kAzure+7); h_noPU_prompt_fake_EE_rebin->SetFillColor(kGreen+2); h_noPU_prompt_no_tErr_EE_rebin->SetFillColor(kYellow-7);
  h_noPU_nonprompt_PV_EB->SetFillColor(kGray+1); h_noPU_nonprompt_PU_EB->SetFillColor(kAzure+7); h_noPU_nonprompt_fake_EB_rebin->SetFillColor(kGreen+2); h_noPU_nonprompt_no_tErr_EB_rebin->SetFillColor(kYellow-7);
  h_noPU_nonprompt_PV_EE->SetFillColor(kGray+1); h_noPU_nonprompt_PU_EE->SetFillColor(kAzure+7); h_noPU_nonprompt_fake_EE_rebin->SetFillColor(kGreen+2); h_noPU_nonprompt_no_tErr_EE_rebin->SetFillColor(kYellow-7);
  // (PV, track)
  h_PU200_prompt_PV_EB_vtx->SetFillColor(kGray+1); h_PU200_prompt_PU_EB_vtx->SetFillColor(kAzure+7); h_PU200_prompt_fake_EB_rebin_vtx->SetFillColor(kGreen+2); h_PU200_prompt_no_tErr_EB_rebin_vtx->SetFillColor(kYellow-7);
  h_PU200_prompt_PV_EE_vtx->SetFillColor(kGray+1); h_PU200_prompt_PU_EE_vtx->SetFillColor(kAzure+7); h_PU200_prompt_fake_EE_rebin_vtx->SetFillColor(kGreen+2); h_PU200_prompt_no_tErr_EE_rebin_vtx->SetFillColor(kYellow-7);
  h_PU200_nonprompt_PV_EB_vtx->SetFillColor(kGray+1); h_PU200_nonprompt_PU_EB_vtx->SetFillColor(kAzure+7); h_PU200_nonprompt_fake_EB_rebin_vtx->SetFillColor(kGreen+2); h_PU200_nonprompt_no_tErr_EB_rebin_vtx->SetFillColor(kYellow-7);
  h_PU200_nonprompt_PV_EE_vtx->SetFillColor(kGray+1); h_PU200_nonprompt_PU_EE_vtx->SetFillColor(kAzure+7); h_PU200_nonprompt_fake_EE_rebin_vtx->SetFillColor(kGreen+2); h_PU200_nonprompt_no_tErr_EE_rebin_vtx->SetFillColor(kYellow-7);
  h_noPU_prompt_PV_EB_vtx->SetFillColor(kGray+1); h_noPU_prompt_PU_EB_vtx->SetFillColor(kAzure+7); h_noPU_prompt_fake_EB_rebin_vtx->SetFillColor(kGreen+2); h_noPU_prompt_no_tErr_EB_rebin_vtx->SetFillColor(kYellow-7);
  h_noPU_prompt_PV_EE_vtx->SetFillColor(kGray+1); h_noPU_prompt_PU_EE_vtx->SetFillColor(kAzure+7); h_noPU_prompt_fake_EE_rebin_vtx->SetFillColor(kGreen+2); h_noPU_prompt_no_tErr_EE_rebin_vtx->SetFillColor(kYellow-7);
  h_noPU_nonprompt_PV_EB_vtx->SetFillColor(kGray+1); h_noPU_nonprompt_PU_EB_vtx->SetFillColor(kAzure+7); h_noPU_nonprompt_fake_EB_rebin_vtx->SetFillColor(kGreen+2); h_noPU_nonprompt_no_tErr_EB_rebin_vtx->SetFillColor(kYellow-7);
  h_noPU_nonprompt_PV_EE_vtx->SetFillColor(kGray+1); h_noPU_nonprompt_PU_EE_vtx->SetFillColor(kAzure+7); h_noPU_nonprompt_fake_EE_rebin_vtx->SetFillColor(kGreen+2); h_noPU_nonprompt_no_tErr_EE_rebin_vtx->SetFillColor(kYellow-7);


  
  /////////////////
  //// Legends ////
  /////////////////
  // (muon, track)
  TLegend *leg_PU200_prompt_EB = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_PV_EB, "Track from PV", "F");
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_PU_EB, "PU Track", "F");
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_fake_EB_rebin, "Fake Track", "F");
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_no_tErr_EB_rebin, "Track without tErr", "F");
  TLegend *leg_PU200_prompt_EE = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_PV_EE, "Track from PV", "F");
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_PU_EE, "PU Track", "F");
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_fake_EE_rebin, "Fake Track", "F");
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_no_tErr_EE_rebin, "Track without tErr", "F");
  TLegend *leg_PU200_nonprompt_EB = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_PV_EB, "Track from PV", "F");
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_PU_EB, "PU Track", "F");
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_fake_EB_rebin, "Fake Track", "F");
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_no_tErr_EB_rebin, "Track without tErr", "F");
  TLegend *leg_PU200_nonprompt_EE = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_PV_EE, "Track from PV", "F");
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_PU_EE, "PU Track", "F");
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_fake_EE_rebin, "Fake Track", "F");
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_no_tErr_EE_rebin, "Track without tErr", "F");
  TLegend *leg_noPU_prompt_EB = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_PV_EB, "Track from PV", "F");
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_PU_EB, "PU Track", "F");
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_fake_EB_rebin, "Fake Track", "F");
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_no_tErr_EB_rebin, "Track without tErr", "F");
  TLegend *leg_noPU_prompt_EE = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_PV_EE, "Track from PV", "F");
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_PU_EE, "PU Track", "F");
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_fake_EE_rebin, "Fake Track", "F");
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_no_tErr_EE_rebin, "Track without tErr", "F");
  TLegend *leg_noPU_nonprompt_EB = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_PV_EB, "Track from PV", "F");
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_PU_EB, "PU Track", "F");
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_fake_EB_rebin, "Fake Track", "F");
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_no_tErr_EB_rebin, "Track without tErr", "F");
  TLegend *leg_noPU_nonprompt_EE = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_PV_EE, "Track from PV", "F");
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_PU_EE, "PU Track", "F");
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_fake_EE_rebin, "Fake Track", "F");
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_no_tErr_EE_rebin, "Track without tErr", "F");
  // (PV, track)
  TLegend *leg_PU200_prompt_EB_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_prompt_EB_vtx->AddEntry(h_PU200_prompt_PV_EB_vtx, "Track from PV", "F");
  leg_PU200_prompt_EB_vtx->AddEntry(h_PU200_prompt_PU_EB_vtx, "PU Track", "F");
  leg_PU200_prompt_EB_vtx->AddEntry(h_PU200_prompt_fake_EB_rebin_vtx, "Fake Track", "F");
  leg_PU200_prompt_EB_vtx->AddEntry(h_PU200_prompt_no_tErr_EB_rebin_vtx, "Track without tErr", "F");
  TLegend *leg_PU200_prompt_EE_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_prompt_EE_vtx->AddEntry(h_PU200_prompt_PV_EE_vtx, "Track from PV", "F");
  leg_PU200_prompt_EE_vtx->AddEntry(h_PU200_prompt_PU_EE_vtx, "PU Track", "F");
  leg_PU200_prompt_EE_vtx->AddEntry(h_PU200_prompt_fake_EE_rebin_vtx, "Fake Track", "F");
  leg_PU200_prompt_EE_vtx->AddEntry(h_PU200_prompt_no_tErr_EE_rebin_vtx, "Track without tErr", "F");
  TLegend *leg_PU200_nonprompt_EB_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_nonprompt_EB_vtx->AddEntry(h_PU200_nonprompt_PV_EB_vtx, "Track from PV", "F");
  leg_PU200_nonprompt_EB_vtx->AddEntry(h_PU200_nonprompt_PU_EB_vtx, "PU Track", "F");
  leg_PU200_nonprompt_EB_vtx->AddEntry(h_PU200_nonprompt_fake_EB_rebin_vtx, "Fake Track", "F");
  leg_PU200_nonprompt_EB_vtx->AddEntry(h_PU200_nonprompt_no_tErr_EB_rebin_vtx, "Track without tErr", "F");
  TLegend *leg_PU200_nonprompt_EE_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_nonprompt_EE_vtx->AddEntry(h_PU200_nonprompt_PV_EE_vtx, "Track from PV", "F");
  leg_PU200_nonprompt_EE_vtx->AddEntry(h_PU200_nonprompt_PU_EE_vtx, "PU Track", "F");
  leg_PU200_nonprompt_EE_vtx->AddEntry(h_PU200_nonprompt_fake_EE_rebin_vtx, "Fake Track", "F");
  leg_PU200_nonprompt_EE_vtx->AddEntry(h_PU200_nonprompt_no_tErr_EE_rebin_vtx, "Track without tErr", "F");
  TLegend *leg_noPU_prompt_EB_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_prompt_EB_vtx->AddEntry(h_noPU_prompt_PV_EB_vtx, "Track from PV", "F");
  leg_noPU_prompt_EB_vtx->AddEntry(h_noPU_prompt_PU_EB_vtx, "PU Track", "F");
  leg_noPU_prompt_EB_vtx->AddEntry(h_noPU_prompt_fake_EB_rebin_vtx, "Fake Track", "F");
  leg_noPU_prompt_EB_vtx->AddEntry(h_noPU_prompt_no_tErr_EB_rebin_vtx, "Track without tErr", "F");
  TLegend *leg_noPU_prompt_EE_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_prompt_EE_vtx->AddEntry(h_noPU_prompt_PV_EE_vtx, "Track from PV", "F");
  leg_noPU_prompt_EE_vtx->AddEntry(h_noPU_prompt_PU_EE_vtx, "PU Track", "F");
  leg_noPU_prompt_EE_vtx->AddEntry(h_noPU_prompt_fake_EE_rebin_vtx, "Fake Track", "F");
  leg_noPU_prompt_EE_vtx->AddEntry(h_noPU_prompt_no_tErr_EE_rebin_vtx, "Track without tErr", "F");
  TLegend *leg_noPU_nonprompt_EB_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_nonprompt_EB_vtx->AddEntry(h_noPU_nonprompt_PV_EB_vtx, "Track from PV", "F");
  leg_noPU_nonprompt_EB_vtx->AddEntry(h_noPU_nonprompt_PU_EB_vtx, "PU Track", "F");
  leg_noPU_nonprompt_EB_vtx->AddEntry(h_noPU_nonprompt_fake_EB_rebin_vtx, "Fake Track", "F");
  leg_noPU_nonprompt_EB_vtx->AddEntry(h_noPU_nonprompt_no_tErr_EB_rebin_vtx, "Track without tErr", "F");
  TLegend *leg_noPU_nonprompt_EE_vtx = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_nonprompt_EE_vtx->AddEntry(h_noPU_nonprompt_PV_EE_vtx, "Track from PV", "F");
  leg_noPU_nonprompt_EE_vtx->AddEntry(h_noPU_nonprompt_PU_EE_vtx, "PU Track", "F");
  leg_noPU_nonprompt_EE_vtx->AddEntry(h_noPU_nonprompt_fake_EE_rebin_vtx, "Fake Track", "F");
  leg_noPU_nonprompt_EE_vtx->AddEntry(h_noPU_nonprompt_no_tErr_EE_rebin_vtx, "Track without tErr", "F");


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

  st_PU200_prompt_EB->Add(h_PU200_prompt_PV_EB); st_PU200_prompt_EB->Add(h_PU200_prompt_PU_EB); st_PU200_prompt_EB->Add(h_PU200_prompt_fake_EB_rebin); st_PU200_prompt_EB->Add(h_PU200_prompt_no_tErr_EB_rebin); st_PU200_prompt_EB->SetTitle(";#sigma_{t};Counts");
  st_PU200_prompt_EE->Add(h_PU200_prompt_PV_EE); st_PU200_prompt_EE->Add(h_PU200_prompt_PU_EE); st_PU200_prompt_EE->Add(h_PU200_prompt_fake_EE_rebin); st_PU200_prompt_EE->Add(h_PU200_prompt_no_tErr_EE_rebin); st_PU200_prompt_EE->SetTitle(";#sigma_{t};Counts");
  st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_PV_EB); st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_PU_EB); st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_fake_EB_rebin); st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_no_tErr_EB_rebin); st_PU200_nonprompt_EB->SetTitle(";#sigma_{t};Counts");
  st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_PV_EE); st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_PU_EE); st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_fake_EE_rebin); st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_no_tErr_EE_rebin); st_PU200_nonprompt_EE->SetTitle(";#sigma_{t};Counts");
  st_noPU_prompt_EB->Add(h_noPU_prompt_PV_EB); st_noPU_prompt_EB->Add(h_noPU_prompt_PU_EB); st_noPU_prompt_EB->Add(h_noPU_prompt_fake_EB_rebin); st_noPU_prompt_EB->Add(h_noPU_prompt_no_tErr_EB_rebin); st_noPU_prompt_EB->SetTitle(";#sigma_{t};Counts");
  st_noPU_prompt_EE->Add(h_noPU_prompt_PV_EE); st_noPU_prompt_EE->Add(h_noPU_prompt_PU_EE); st_noPU_prompt_EE->Add(h_noPU_prompt_fake_EE_rebin); st_noPU_prompt_EE->Add(h_noPU_prompt_no_tErr_EE_rebin); st_noPU_prompt_EE->SetTitle(";#sigma_{t};Counts");
  st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_PV_EB); st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_PU_EB); st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_fake_EB_rebin); st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_no_tErr_EB_rebin); st_noPU_nonprompt_EB->SetTitle(";#sigma_{t};Counts");
  st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_PV_EE); st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_PU_EE); st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_fake_EE_rebin); st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_no_tErr_EE_rebin); st_noPU_nonprompt_EE->SetTitle(";#sigma_{t};Counts");
  // (PV, track)
  THStack* st_PU200_prompt_EB_vtx    = new THStack("st_PU200_prompt_EB_vtx","st_PU200_prompt_EB_vtx");
  THStack* st_PU200_prompt_EE_vtx    = new THStack("st_PU200_prompt_EE_vtx","st_PU200_prompt_EE_vtx");
  THStack* st_PU200_nonprompt_EB_vtx = new THStack("st_PU200_nonprompt_EB_vtx","st_PU200_nonprompt_EB_vtx");
  THStack* st_PU200_nonprompt_EE_vtx = new THStack("st_PU200_nonprompt_EE_vtx","st_PU200_nonprompt_EE_vtx");
  THStack* st_noPU_prompt_EB_vtx     = new THStack("st_noPU_prompt_EB_vtx","st_noPU_prompt_EB_vtx");
  THStack* st_noPU_prompt_EE_vtx     = new THStack("st_noPU_prompt_EE_vtx","st_noPU_prompt_EE_vtx");
  THStack* st_noPU_nonprompt_EB_vtx  = new THStack("st_noPU_nonprompt_EB_vtx","st_noPU_nonprompt_EB_vtx");
  THStack* st_noPU_nonprompt_EE_vtx  = new THStack("st_noPU_nonprompt_EE_vtx","st_noPU_nonprompt_EE_vtx");

  st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_PV_EB_vtx); st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_PU_EB_vtx); st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_fake_EB_rebin_vtx); st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_no_tErr_EB_rebin_vtx); st_PU200_prompt_EB_vtx->SetTitle(";#sigma_{t};Counts");
  st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_PV_EE_vtx); st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_PU_EE_vtx); st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_fake_EE_rebin_vtx); st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_no_tErr_EE_rebin_vtx); st_PU200_prompt_EE_vtx->SetTitle(";#sigma_{t};Counts");
  st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_PV_EB_vtx); st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_PU_EB_vtx); st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_fake_EB_rebin_vtx); st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_no_tErr_EB_rebin_vtx); st_PU200_nonprompt_EB_vtx->SetTitle(";#sigma_{t};Counts");
  st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_PV_EE_vtx); st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_PU_EE_vtx); st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_fake_EE_rebin_vtx); st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_no_tErr_EE_rebin_vtx); st_PU200_nonprompt_EE_vtx->SetTitle(";#sigma_{t};Counts");
  st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_PV_EB_vtx); st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_PU_EB_vtx); st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_fake_EB_rebin_vtx); st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_no_tErr_EB_rebin_vtx); st_noPU_prompt_EB_vtx->SetTitle(";#sigma_{t};Counts");
  st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_PV_EE_vtx); st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_PU_EE_vtx); st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_fake_EE_rebin_vtx); st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_no_tErr_EE_rebin_vtx); st_noPU_prompt_EE_vtx->SetTitle(";#sigma_{t};Counts");
  st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_PV_EB_vtx); st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_PU_EB_vtx); st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_fake_EB_rebin_vtx); st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_no_tErr_EB_rebin_vtx); st_noPU_nonprompt_EB_vtx->SetTitle(";#sigma_{t};Counts");
  st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_PV_EE_vtx); st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_PU_EE_vtx); st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_fake_EE_rebin_vtx); st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_no_tErr_EE_rebin_vtx); st_noPU_nonprompt_EE_vtx->SetTitle(";#sigma_{t};Counts");

  ////////////////
  // Draw plots //
  ////////////////
  // (muon, track)
  // PU200
    // prompt
      // Barrel
  TCanvas* c_st_PU200_prompt_EB = new TCanvas("c_st_PU200_prompt_EB", "c_st_PU200_prompt_EB", 1500, 1500);
  c_st_PU200_prompt_EB->cd();
  c_st_PU200_prompt_EB->SetLeftMargin(0.12);
  st_PU200_prompt_EB->Draw();
  leg_PU200_prompt_EB->Draw();
  c_st_PU200_prompt_EB->Print("plots/track_sigma_PU200_prompt_EB_v1.pdf");
  cout << "PU200 Barrel prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EB->Integral(1,-1) << " : " << h_PU200_prompt_PV_EB->Integral(1,2) << " : " << h_PU200_prompt_PV_EB->Integral(3,-1) << "  " << "(" << h_PU200_prompt_PV_EB->Integral(3,-1)/h_PU200_prompt_PV_EB->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EB->Integral(1,-1) << " : " << h_PU200_prompt_PU_EB->Integral(1,2) << " : " << h_PU200_prompt_PU_EB->Integral(3,-1) << "  " << "(" << h_PU200_prompt_PU_EB->Integral(3,-1)/h_PU200_prompt_PU_EB->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EB_rebin->Integral() << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EB_rebin->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_PU200_prompt_EE = new TCanvas("c_st_PU200_prompt_EE", "c_st_PU200_prompt_EE", 1500, 1500);
  c_st_PU200_prompt_EE->cd();
  c_st_PU200_prompt_EE->SetLeftMargin(0.12);
  st_PU200_prompt_EE->Draw();
  leg_PU200_prompt_EE->Draw();
  c_st_PU200_prompt_EE->Print("plots/track_sigma_PU200_prompt_EE_v1.pdf");
  cout << "PU200 Endcap prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EE->Integral(1,-1) << " : " << h_PU200_prompt_PV_EE->Integral(1,2) << " : " << h_PU200_prompt_PV_EE->Integral(3,-1) << "  " << "(" << h_PU200_prompt_PV_EE->Integral(3,-1)/h_PU200_prompt_PV_EE->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EE->Integral(1,-1) << " : " << h_PU200_prompt_PU_EE->Integral(1,2) << " : " << h_PU200_prompt_PU_EE->Integral(3,-1) << "  " << "(" << h_PU200_prompt_PU_EE->Integral(3,-1)/h_PU200_prompt_PU_EE->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EE_rebin->Integral() << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EE_rebin->Integral() << endl;
  cout << endl;
    // nonprompt
      // Barrel
  TCanvas* c_st_PU200_nonprompt_EB = new TCanvas("c_st_PU200_nonprompt_EB", "c_st_PU200_nonprompt_EB", 1500, 1500);
  c_st_PU200_nonprompt_EB->cd();
  c_st_PU200_nonprompt_EB->SetLeftMargin(0.12);
  st_PU200_nonprompt_EB->Draw();
  leg_PU200_nonprompt_EB->Draw();
  c_st_PU200_nonprompt_EB->Print("plots/track_sigma_PU200_nonprompt_EB_v1.pdf");
  cout << "PU200 Barrel nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EB->Integral(1,2) << " : " << h_PU200_nonprompt_PV_EB->Integral(3,-1) << "  " << "(" << h_PU200_nonprompt_PV_EB->Integral(3,-1)/h_PU200_nonprompt_PV_EB->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EB->Integral(1,2) << " : " << h_PU200_nonprompt_PU_EB->Integral(3,-1) << "  " << "(" << h_PU200_nonprompt_PU_EB->Integral(3,-1)/h_PU200_nonprompt_PU_EB->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EB_rebin->Integral() << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EB_rebin->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_PU200_nonprompt_EE = new TCanvas("c_st_PU200_nonprompt_EE", "c_st_PU200_nonprompt_EE", 1500, 1500);
  c_st_PU200_nonprompt_EE->cd();
  c_st_PU200_nonprompt_EE->SetLeftMargin(0.12);
  st_PU200_nonprompt_EE->Draw();
  leg_PU200_nonprompt_EE->Draw();
  c_st_PU200_nonprompt_EE->Print("plots/track_sigma_PU200_nonprompt_EE_v1.pdf");
  cout << "PU200 Endcap nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EE->Integral(1,2) << " : " << h_PU200_nonprompt_PV_EE->Integral(3,-1) << "  " << "(" << h_PU200_nonprompt_PV_EE->Integral(3,-1)/h_PU200_nonprompt_PV_EE->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EE->Integral(1,2) << " : " << h_PU200_nonprompt_PU_EE->Integral(3,-1) << "  " << "(" << h_PU200_nonprompt_PU_EE->Integral(3,-1)/h_PU200_nonprompt_PU_EE->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EE_rebin->Integral() << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EE_rebin->Integral() << endl;
  cout << endl;

  // noPU
    // prompt
      // Barrel
  TCanvas* c_st_noPU_prompt_EB = new TCanvas("c_st_noPU_prompt_EB", "c_st_noPU_prompt_EB", 1500, 1500);
  c_st_noPU_prompt_EB->cd();
  c_st_noPU_prompt_EB->SetLeftMargin(0.12);
  st_noPU_prompt_EB->Draw();
  leg_noPU_prompt_EB->Draw();
  c_st_noPU_prompt_EB->Print("plots/track_sigma_noPU_prompt_EB_v1.pdf");
  cout << "noPU Barrel prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EB->Integral(1,-1) << " : " << h_noPU_prompt_PV_EB->Integral(1,2) << " : " << h_noPU_prompt_PV_EB->Integral(3,-1) << "  " << "(" << h_noPU_prompt_PV_EB->Integral(3,-1)/h_noPU_prompt_PV_EB->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EB->Integral(1,-1) << " : " << h_noPU_prompt_PU_EB->Integral(1,2) << " : " << h_noPU_prompt_PU_EB->Integral(3,-1) << "  " << "(" << h_noPU_prompt_PU_EB->Integral(3,-1)/h_noPU_prompt_PU_EB->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EB_rebin->Integral() << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EB_rebin->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_noPU_prompt_EE = new TCanvas("c_st_noPU_prompt_EE", "c_st_noPU_prompt_EE", 1500, 1500);
  c_st_noPU_prompt_EE->cd();
  c_st_noPU_prompt_EE->SetLeftMargin(0.12);
  st_noPU_prompt_EE->Draw();
  leg_noPU_prompt_EE->Draw();
  c_st_noPU_prompt_EE->Print("plots/track_sigma_noPU_prompt_EE_v1.pdf");
  cout << "noPU Endcap prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EE->Integral(1,-1) << " : " << h_noPU_prompt_PV_EE->Integral(1,2) << " : " << h_noPU_prompt_PV_EE->Integral(3,-1) << "  " << "(" << h_noPU_prompt_PV_EE->Integral(3,-1)/h_noPU_prompt_PV_EE->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EE->Integral(1,-1) << " : " << h_noPU_prompt_PU_EE->Integral(1,2) << " : " << h_noPU_prompt_PU_EE->Integral(3,-1) << "  " << "(" << h_noPU_prompt_PU_EE->Integral(3,-1)/h_noPU_prompt_PU_EE->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EE_rebin->Integral() << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EE_rebin->Integral() << endl;
  cout << endl;
    // nonprompt
      // Barrel
  TCanvas* c_st_noPU_nonprompt_EB = new TCanvas("c_st_noPU_nonprompt_EB", "c_st_noPU_nonprompt_EB", 1500, 1500);
  c_st_noPU_nonprompt_EB->cd();
  c_st_noPU_nonprompt_EB->SetLeftMargin(0.12);
  st_noPU_nonprompt_EB->Draw();
  leg_noPU_nonprompt_EB->Draw();
  c_st_noPU_nonprompt_EB->Print("plots/track_sigma_noPU_nonprompt_EB_v1.pdf");
  cout << "noPU Barrel nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EB->Integral(1,2) << " : " << h_noPU_nonprompt_PV_EB->Integral(3,-1) << "  " << "(" << h_noPU_nonprompt_PV_EB->Integral(3,-1)/h_noPU_nonprompt_PV_EB->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EB->Integral(1,2) << " : " << h_noPU_nonprompt_PU_EB->Integral(3,-1) << "  " << "(" << h_noPU_nonprompt_PU_EB->Integral(3,-1)/h_noPU_nonprompt_PU_EB->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EB_rebin->Integral() << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EB_rebin->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_noPU_nonprompt_EE = new TCanvas("c_st_noPU_nonprompt_EE", "c_st_noPU_nonprompt_EE", 1500, 1500);
  c_st_noPU_nonprompt_EE->cd();
  c_st_noPU_nonprompt_EE->SetLeftMargin(0.12);
  st_noPU_nonprompt_EE->Draw();
  leg_noPU_nonprompt_EE->Draw();
  c_st_noPU_nonprompt_EE->Print("plots/track_sigma_noPU_nonprompt_EE.pdf");
  cout << "noPU Endcap nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EE->Integral(1,2) << " : " << h_noPU_nonprompt_PV_EE->Integral(3,-1) << "  " << "(" << h_noPU_nonprompt_PV_EE->Integral(3,-1)/h_noPU_nonprompt_PV_EE->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EE->Integral(1,2) << " : " << h_noPU_nonprompt_PU_EE->Integral(3,-1) << "  " << "(" << h_noPU_nonprompt_PU_EE->Integral(3,-1)/h_noPU_nonprompt_PU_EE->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EE_rebin->Integral() << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EE_rebin->Integral() << endl;
  cout << endl;


  // (PV, track)
  // PU200
    // prompt
      // Barrel
  TCanvas* c_st_PU200_prompt_EB_vtx = new TCanvas("c_st_PU200_prompt_EB_vtx", "c_st_PU200_prompt_EB_vtx", 1500, 1500);
  c_st_PU200_prompt_EB_vtx->cd();
  c_st_PU200_prompt_EB_vtx->SetLeftMargin(0.12);
  st_PU200_prompt_EB_vtx->Draw();
  leg_PU200_prompt_EB_vtx->Draw();
  c_st_PU200_prompt_EB_vtx->Print("plots/track_sigma_PU200_prompt_EB_vtx.pdf");
  cout << "PU200 Barrel prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PV_EB_vtx->Integral(1,2) << " : " << h_PU200_prompt_PV_EB_vtx->Integral(3,-1) << "  " << "(" << h_PU200_prompt_PV_EB_vtx->Integral(3,-1)/h_PU200_prompt_PV_EB_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PU_EB_vtx->Integral(1,2) << " : " << h_PU200_prompt_PU_EB_vtx->Integral(3,-1) << "  " << "(" << h_PU200_prompt_PU_EB_vtx->Integral(3,-1)/h_PU200_prompt_PU_EB_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EB_rebin_vtx->Integral() << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EB_rebin_vtx->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_PU200_prompt_EE_vtx = new TCanvas("c_st_PU200_prompt_EE_vtx", "c_st_PU200_prompt_EE_vtx", 1500, 1500);
  c_st_PU200_prompt_EE_vtx->cd();
  c_st_PU200_prompt_EE_vtx->SetLeftMargin(0.12);
  st_PU200_prompt_EE_vtx->Draw();
  leg_PU200_prompt_EE_vtx->Draw();
  c_st_PU200_prompt_EE_vtx->Print("plots/track_sigma_PU200_prompt_EE_vtx.pdf");
  cout << "PU200 Endcap prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PV_EE_vtx->Integral(1,2) << " : " << h_PU200_prompt_PV_EE_vtx->Integral(3,-1) << "  " << "(" << h_PU200_prompt_PV_EE_vtx->Integral(3,-1)/h_PU200_prompt_PV_EE_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PU_EE_vtx->Integral(1,2) << " : " << h_PU200_prompt_PU_EE_vtx->Integral(3,-1) << "  " << "(" << h_PU200_prompt_PU_EE_vtx->Integral(3,-1)/h_PU200_prompt_PU_EE_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EE_rebin_vtx->Integral() << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EE_rebin_vtx->Integral() << endl;
  cout << endl;
    // nonprompt
      // Barrel
  TCanvas* c_st_PU200_nonprompt_EB_vtx = new TCanvas("c_st_PU200_nonprompt_EB_vtx", "c_st_PU200_nonprompt_EB_vtx", 1500, 1500);
  c_st_PU200_nonprompt_EB_vtx->cd();
  c_st_PU200_nonprompt_EB_vtx->SetLeftMargin(0.12);
  st_PU200_nonprompt_EB_vtx->Draw();
  leg_PU200_nonprompt_EB_vtx->Draw();
  c_st_PU200_nonprompt_EB_vtx->Print("plots/track_sigma_PU200_nonprompt_EB_vtx.pdf");
  cout << "PU200 Barrel nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EB_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_PV_EB_vtx->Integral(3,-1) << "  " << "(" << h_PU200_nonprompt_PV_EB_vtx->Integral(3,-1)/h_PU200_nonprompt_PV_EB_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EB_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_PU_EB_vtx->Integral(3,-1) << "  " << "(" << h_PU200_nonprompt_PU_EB_vtx->Integral(3,-1)/h_PU200_nonprompt_PU_EB_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EB_rebin_vtx->Integral() << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EB_rebin_vtx->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_PU200_nonprompt_EE_vtx = new TCanvas("c_st_PU200_nonprompt_EE_vtx", "c_st_PU200_nonprompt_EE_vtx", 1500, 1500);
  c_st_PU200_nonprompt_EE_vtx->cd();
  c_st_PU200_nonprompt_EE_vtx->SetLeftMargin(0.12);
  st_PU200_nonprompt_EE_vtx->Draw();
  leg_PU200_nonprompt_EE_vtx->Draw();
  c_st_PU200_nonprompt_EE_vtx->Print("plots/track_sigma_PU200_nonprompt_EE_vtx.pdf");
  cout << "PU200 Endcap nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EE_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_PV_EE_vtx->Integral(3,-1) << "  " << "(" << h_PU200_nonprompt_PV_EE_vtx->Integral(3,-1)/h_PU200_nonprompt_PV_EE_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EE_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_PU_EE_vtx->Integral(3,-1) << "  " << "(" << h_PU200_nonprompt_PU_EE_vtx->Integral(3,-1)/h_PU200_nonprompt_PU_EE_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EE_rebin_vtx->Integral() << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EE_rebin_vtx->Integral() << endl;
  cout << endl;

  // noPU
    // prompt
      // Barrel
  TCanvas* c_st_noPU_prompt_EB_vtx = new TCanvas("c_st_noPU_prompt_EB_vtx", "c_st_noPU_prompt_EB_vtx", 1500, 1500);
  c_st_noPU_prompt_EB_vtx->cd();
  c_st_noPU_prompt_EB_vtx->SetLeftMargin(0.12);
  st_noPU_prompt_EB_vtx->Draw();
  leg_noPU_prompt_EB_vtx->Draw();
  c_st_noPU_prompt_EB_vtx->Print("plots/track_sigma_noPU_prompt_EB_vtx.pdf");
  cout << "noPU Barrel prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PV_EB_vtx->Integral(1,2) << " : " << h_noPU_prompt_PV_EB_vtx->Integral(3,-1) << "  " << "(" << h_noPU_prompt_PV_EB_vtx->Integral(3,-1)/h_noPU_prompt_PV_EB_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PU_EB_vtx->Integral(1,2) << " : " << h_noPU_prompt_PU_EB_vtx->Integral(3,-1) << "  " << "(" << h_noPU_prompt_PU_EB_vtx->Integral(3,-1)/h_noPU_prompt_PU_EB_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EB_rebin_vtx->Integral() << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EB_rebin_vtx->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_noPU_prompt_EE_vtx = new TCanvas("c_st_noPU_prompt_EE_vtx", "c_st_noPU_prompt_EE_vtx", 1500, 1500);
  c_st_noPU_prompt_EE_vtx->cd();
  c_st_noPU_prompt_EE_vtx->SetLeftMargin(0.12);
  st_noPU_prompt_EE_vtx->Draw();
  leg_noPU_prompt_EE_vtx->Draw();
  c_st_noPU_prompt_EE_vtx->Print("plots/track_sigma_noPU_prompt_EE_vtx.pdf");
  cout << "noPU Endcap prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PV_EE_vtx->Integral(1,2) << " : " << h_noPU_prompt_PV_EE_vtx->Integral(3,-1) << "  " << "(" << h_noPU_prompt_PV_EE_vtx->Integral(3,-1)/h_noPU_prompt_PV_EE_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PU_EE_vtx->Integral(1,2) << " : " << h_noPU_prompt_PU_EE_vtx->Integral(3,-1) << "  " << "(" << h_noPU_prompt_PU_EE_vtx->Integral(3,-1)/h_noPU_prompt_PU_EE_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EE_rebin_vtx->Integral() << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EE_rebin_vtx->Integral() << endl;
  cout << endl;
    // nonprompt
      // Barrel
  TCanvas* c_st_noPU_nonprompt_EB_vtx = new TCanvas("c_st_noPU_nonprompt_EB_vtx", "c_st_noPU_nonprompt_EB_vtx", 1500, 1500);
  c_st_noPU_nonprompt_EB_vtx->cd();
  c_st_noPU_nonprompt_EB_vtx->SetLeftMargin(0.12);
  st_noPU_nonprompt_EB_vtx->Draw();
  leg_noPU_nonprompt_EB_vtx->Draw();
  c_st_noPU_nonprompt_EB_vtx->Print("plots/track_sigma_noPU_nonprompt_EB_vtx.pdf");
  cout << "noPU Barrel nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EB_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_PV_EB_vtx->Integral(3,-1) << "  " << "(" << h_noPU_nonprompt_PV_EB_vtx->Integral(3,-1)/h_noPU_nonprompt_PV_EB_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EB_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_PU_EB_vtx->Integral(3,-1) << "  " << "(" << h_noPU_nonprompt_PU_EB_vtx->Integral(3,-1)/h_noPU_nonprompt_PU_EB_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EB_rebin_vtx->Integral() << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EB_rebin_vtx->Integral() << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_noPU_nonprompt_EE_vtx = new TCanvas("c_st_noPU_nonprompt_EE_vtx", "c_st_noPU_nonprompt_EE_vtx", 1500, 1500);
  c_st_noPU_nonprompt_EE_vtx->cd();
  c_st_noPU_nonprompt_EE_vtx->SetLeftMargin(0.12);
  st_noPU_nonprompt_EE_vtx->Draw();
  leg_noPU_nonprompt_EE_vtx->Draw();
  c_st_noPU_nonprompt_EE_vtx->Print("plots/track_sigma_noPU_nonprompt_EE_vtx.pdf");
  cout << "noPU Endcap nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EE_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_PV_EE_vtx->Integral(3,-1) << "  " << "(" << h_noPU_nonprompt_PV_EE_vtx->Integral(3,-1)/h_noPU_nonprompt_PV_EE_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EE_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_PU_EE_vtx->Integral(3,-1) << "  " << "(" << h_noPU_nonprompt_PU_EE_vtx->Integral(3,-1)/h_noPU_nonprompt_PU_EE_vtx->Integral(1,-1)*100 << " % is rejected)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EE_rebin_vtx->Integral() << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EE_rebin_vtx->Integral() << endl;
  cout << endl;

}

void track_type_v2() {
  gErrorIgnoreLevel = kWarning;
  gStyle->SetOptStat(0);
  gROOT->ProcessLine("#include <vector>");

  TFile *f_PU200_prompt, *f_PU200_nonprompt, *f_noPU_prompt, *f_noPU_nonprompt;
  TFile *f_PU200_prompt_vtx, *f_PU200_nonprompt_vtx, *f_noPU_prompt_vtx, *f_noPU_nonprompt_vtx;
  // PU200
  f_PU200_prompt = new TFile(Form("data/%s", path_PU200_prompt.Data()));
  f_PU200_nonprompt = new TFile(Form("data/%s", path_PU200_nonprompt.Data()));
  f_PU200_prompt_vtx = new TFile(Form("data/%s", path_PU200_prompt_vtx.Data()));
  f_PU200_nonprompt_vtx = new TFile(Form("data/%s", path_PU200_nonprompt_vtx.Data()));
  // noPU
  f_noPU_prompt = new TFile(Form("data/%s", path_noPU_prompt.Data()));
  f_noPU_nonprompt = new TFile(Form("data/%s", path_noPU_nonprompt.Data()));
  f_noPU_prompt_vtx = new TFile(Form("data/%s", path_noPU_prompt_vtx.Data()));
  f_noPU_nonprompt_vtx = new TFile(Form("data/%s", path_noPU_nonprompt_vtx.Data()));

  TString dir="/DQMData/Run 1/MTD/Run summary/MuonIso/";
  TH1D *h_PU200_prompt_type_v2_EB, *h_PU200_prompt_type_v2_EE, *h_PU200_nonprompt_type_v2_EB, *h_PU200_nonprompt_type_v2_EE;
  TH1D *h_noPU_prompt_type_v2_EB, *h_noPU_prompt_type_v2_EE, *h_noPU_nonprompt_type_v2_EB, *h_noPU_nonprompt_type_v2_EE;
    //vtx
  TH1D *h_PU200_prompt_type_v2_EB_vtx, *h_PU200_prompt_type_v2_EE_vtx, *h_PU200_nonprompt_type_v2_EB_vtx, *h_PU200_nonprompt_type_v2_EE_vtx;
  TH1D *h_noPU_prompt_type_v2_EB_vtx, *h_noPU_prompt_type_v2_EE_vtx, *h_noPU_nonprompt_type_v2_EB_vtx, *h_noPU_nonprompt_type_v2_EE_vtx;


  // Barrel
    // PU200
  h_PU200_prompt_type_v2_EB    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_track_type_Sig_EB");
  h_PU200_nonprompt_type_v2_EB = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_track_type_Bkg_EB");
  h_PU200_prompt_type_v2_EB_vtx    = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_Iso_track_type_Sig_EB");
  h_PU200_nonprompt_type_v2_EB_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_Iso_track_type_Bkg_EB");
    // noPU
  h_noPU_prompt_type_v2_EB    = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_track_type_Sig_EB");
  h_noPU_nonprompt_type_v2_EB = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_track_type_Bkg_EB");
  h_noPU_prompt_type_v2_EB_vtx    = (TH1D*)f_noPU_prompt_vtx->Get(dir+"Muon_Iso_track_type_Sig_EB");
  h_noPU_nonprompt_type_v2_EB_vtx = (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"Muon_Iso_track_type_Bkg_EB");

  // Endcap
    // PU200
  h_PU200_prompt_type_v2_EE    = (TH1D*)f_PU200_prompt->Get(dir+"Muon_Iso_track_type_Sig_EE");
  h_PU200_nonprompt_type_v2_EE = (TH1D*)f_PU200_nonprompt->Get(dir+"Muon_Iso_track_type_Bkg_EE");
  h_PU200_prompt_type_v2_EE_vtx    = (TH1D*)f_PU200_prompt_vtx->Get(dir+"Muon_Iso_track_type_Sig_EE");
  h_PU200_nonprompt_type_v2_EE_vtx = (TH1D*)f_PU200_nonprompt_vtx->Get(dir+"Muon_Iso_track_type_Bkg_EE");
    // noPU
  h_noPU_prompt_type_v2_EE    = (TH1D*)f_noPU_prompt->Get(dir+"Muon_Iso_track_type_Sig_EE");
  h_noPU_nonprompt_type_v2_EE = (TH1D*)f_noPU_nonprompt->Get(dir+"Muon_Iso_track_type_Bkg_EE");
  h_noPU_prompt_type_v2_EE_vtx    = (TH1D*)f_noPU_prompt_vtx->Get(dir+"Muon_Iso_track_type_Sig_EE");
  h_noPU_nonprompt_type_v2_EE_vtx = (TH1D*)f_noPU_nonprompt_vtx->Get(dir+"Muon_Iso_track_type_Bkg_EE");



  ////////////////
  // Draw plots //
  ////////////////
  // PU200
    // prompt
      // Barrel
  TCanvas* c_PU200_prompt_type_v2_EB = new TCanvas("c_PU200_prompt_type_v2_EB", "c_PU200_prompt_type_v2_EB", 1500, 1500);
  c_PU200_prompt_type_v2_EB->cd();
  c_PU200_prompt_type_v2_EB->SetLeftMargin(0.12);
  TH1D* h_PU200_prompt_type_v2_EB_incl = new TH1D("h_PU200_prompt_type_v2_EB_incl","h_PU200_prompt_type_v2_EB_incl",4,0,4);
  h_PU200_prompt_type_v2_EB_incl->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt_type_v2_EB_incl->SetTitle("");
  h_PU200_prompt_type_v2_EB_incl->SetMinimum(0.);
  h_PU200_prompt_type_v2_EB_incl->SetBinContent(1,h_PU200_prompt_type_v2_EB->GetBinContent(1));
  h_PU200_prompt_type_v2_EB_incl->SetBinContent(2,h_PU200_prompt_type_v2_EB->GetBinContent(2));
  h_PU200_prompt_type_v2_EB_incl->SetBinContent(3,h_PU200_prompt_type_v2_EB->GetBinContent(3));
  h_PU200_prompt_type_v2_EB_incl->SetBinContent(4,h_PU200_prompt_type_v2_EB->GetBinContent(4));
  h_PU200_prompt_type_v2_EB_incl->Draw("hist");
  c_PU200_prompt_type_v2_EB->Print("plots/trktype_v2_PU200_prompt_EB.pdf");
  cout << "PU200" << endl;
  cout << "Barrel prompt    |" << h_PU200_prompt_type_v2_EB->GetBinContent(1) << " : " << h_PU200_prompt_type_v2_EB->GetBinContent(2) << " : " << h_PU200_prompt_type_v2_EB->GetBinContent(3) << " : " << h_PU200_prompt_type_v2_EB->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_PU200_prompt_type_v2_EE = new TCanvas("c_PU200_prompt_type_v2_EE", "c_PU200_prompt_type_v2_EE", 1500, 1500);
  c_PU200_prompt_type_v2_EE->cd();
  c_PU200_prompt_type_v2_EE->SetLeftMargin(0.12);
  TH1D* h_PU200_prompt_type_v2_EE_incl = new TH1D("h_PU200_prompt_type_v2_EE_incl","h_PU200_prompt_type_v2_EE_incl",4,0,4);
  h_PU200_prompt_type_v2_EE_incl->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt_type_v2_EE_incl->SetTitle("");
  h_PU200_prompt_type_v2_EE_incl->SetMinimum(0.);
  h_PU200_prompt_type_v2_EE_incl->SetBinContent(1,h_PU200_prompt_type_v2_EE->GetBinContent(1));
  h_PU200_prompt_type_v2_EE_incl->SetBinContent(2,h_PU200_prompt_type_v2_EE->GetBinContent(2));
  h_PU200_prompt_type_v2_EE_incl->SetBinContent(3,h_PU200_prompt_type_v2_EE->GetBinContent(3));
  h_PU200_prompt_type_v2_EE_incl->SetBinContent(4,h_PU200_prompt_type_v2_EE->GetBinContent(4));
  h_PU200_prompt_type_v2_EE_incl->Draw("hist");
  c_PU200_prompt_type_v2_EE->Print("plots/trktype_v2_PU200_prompt_EE.pdf");
  cout << "Endcap prompt    |" << h_PU200_prompt_type_v2_EE->GetBinContent(1) << " : " << h_PU200_prompt_type_v2_EE->GetBinContent(2) << " : " << h_PU200_prompt_type_v2_EE->GetBinContent(3) << " : " << h_PU200_prompt_type_v2_EE->GetBinContent(4) << endl;
    // nonprompt
      // Barrel
  TCanvas* c_PU200_nonprompt_type_v2_EB = new TCanvas("c_PU200_nonprompt_type_v2_EB", "c_PU200_nonprompt_type_v2_EB", 1500, 1500);
  c_PU200_nonprompt_type_v2_EB->cd();
  c_PU200_nonprompt_type_v2_EB->SetLeftMargin(0.12);
  TH1D* h_PU200_nonprompt_type_v2_EB_incl = new TH1D("h_PU200_nonprompt_type_v2_EB_incl","h_PU200_nonprompt_type_v2_EB_incl",4,0,4);
  h_PU200_nonprompt_type_v2_EB_incl->GetYaxis()->SetTitle("Counts");
  h_PU200_nonprompt_type_v2_EB_incl->SetTitle("");
  h_PU200_nonprompt_type_v2_EB_incl->SetMinimum(0.);
  h_PU200_nonprompt_type_v2_EB_incl->SetBinContent(1,h_PU200_nonprompt_type_v2_EB->GetBinContent(1));
  h_PU200_nonprompt_type_v2_EB_incl->SetBinContent(2,h_PU200_nonprompt_type_v2_EB->GetBinContent(2));
  h_PU200_nonprompt_type_v2_EB_incl->SetBinContent(3,h_PU200_nonprompt_type_v2_EB->GetBinContent(3));
  h_PU200_nonprompt_type_v2_EB_incl->SetBinContent(4,h_PU200_nonprompt_type_v2_EB->GetBinContent(4));
  h_PU200_nonprompt_type_v2_EB_incl->Draw("hist");
  c_PU200_nonprompt_type_v2_EB->Print("plots/trktype_v2_PU200_nonprompt_EB.pdf");
  cout << "Barrel nonprompt |" << h_PU200_nonprompt_type_v2_EB->GetBinContent(1) << " : " << h_PU200_nonprompt_type_v2_EB->GetBinContent(2) << " : " << h_PU200_nonprompt_type_v2_EB->GetBinContent(3) << " : " << h_PU200_nonprompt_type_v2_EB->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_PU200_nonprompt_type_v2_EE = new TCanvas("c_PU200_nonprompt_type_v2_EE", "c_PU200_nonprompt_type_v2_EE", 1500, 1500);
  c_PU200_nonprompt_type_v2_EE->cd();
  c_PU200_nonprompt_type_v2_EE->SetLeftMargin(0.12);
  TH1D* h_PU200_nonprompt_type_v2_EE_incl = new TH1D("h_PU200_nonprompt_type_v2_EE_incl","h_PU200_nonprompt_type_v2_EE_incl",4,0,4);
  h_PU200_nonprompt_type_v2_EE_incl->GetYaxis()->SetTitle("Counts");
  h_PU200_nonprompt_type_v2_EE_incl->SetTitle("");
  h_PU200_nonprompt_type_v2_EE_incl->SetMinimum(0.);
  h_PU200_nonprompt_type_v2_EE_incl->SetBinContent(1,h_PU200_nonprompt_type_v2_EE->GetBinContent(1));
  h_PU200_nonprompt_type_v2_EE_incl->SetBinContent(2,h_PU200_nonprompt_type_v2_EE->GetBinContent(2));
  h_PU200_nonprompt_type_v2_EE_incl->SetBinContent(3,h_PU200_nonprompt_type_v2_EE->GetBinContent(3));
  h_PU200_nonprompt_type_v2_EE_incl->SetBinContent(4,h_PU200_nonprompt_type_v2_EE->GetBinContent(4));
  h_PU200_nonprompt_type_v2_EE_incl->Draw("hist");
  c_PU200_nonprompt_type_v2_EE->Print("plots/trktype_v2_PU200_nonprompt_EE.pdf");
  cout << "Endcap nonprompt |" << h_PU200_nonprompt_type_v2_EE->GetBinContent(1) << " : " << h_PU200_nonprompt_type_v2_EE->GetBinContent(2) << " : " << h_PU200_nonprompt_type_v2_EE->GetBinContent(3)  << " : " << h_PU200_nonprompt_type_v2_EE->GetBinContent(4)<< endl;
  // noPU
    // prompt
      // Barrel
  TCanvas* c_noPU_prompt_type_v2_EB = new TCanvas("c_noPU_prompt_type_v2_EB", "c_noPU_prompt_type_v2_EB", 1500, 1500);
  c_noPU_prompt_type_v2_EB->cd();
  c_noPU_prompt_type_v2_EB->SetLeftMargin(0.12);
  TH1D* h_noPU_prompt_type_v2_EB_incl = new TH1D("h_noPU_prompt_type_v2_EB_incl","h_noPU_prompt_type_v2_EB_incl",4,0,4);
  h_noPU_prompt_type_v2_EB_incl->GetYaxis()->SetTitle("Counts");
  h_noPU_prompt_type_v2_EB_incl->SetTitle("");
  h_noPU_prompt_type_v2_EB_incl->SetMinimum(0.);
  h_noPU_prompt_type_v2_EB_incl->SetBinContent(1,h_noPU_prompt_type_v2_EB->GetBinContent(1));
  h_noPU_prompt_type_v2_EB_incl->SetBinContent(2,h_noPU_prompt_type_v2_EB->GetBinContent(2));
  h_noPU_prompt_type_v2_EB_incl->SetBinContent(3,h_noPU_prompt_type_v2_EB->GetBinContent(3));
  h_noPU_prompt_type_v2_EB_incl->SetBinContent(4,h_noPU_prompt_type_v2_EB->GetBinContent(4));
  h_noPU_prompt_type_v2_EB_incl->Draw("hist");
  c_noPU_prompt_type_v2_EB->Print("plots/trktype_v2_noPU_prompt_EB.pdf");
  cout << endl;
  cout << "noPU" << endl;
  cout << "Barrel prompt    |" << h_noPU_prompt_type_v2_EB->GetBinContent(1) << " : " << h_noPU_prompt_type_v2_EB->GetBinContent(2) << " : " << h_noPU_prompt_type_v2_EB->GetBinContent(3) << " : " << h_noPU_prompt_type_v2_EB->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_noPU_prompt_type_v2_EE = new TCanvas("c_noPU_prompt_type_v2_EE", "c_noPU_prompt_type_v2_EE", 1500, 1500);
  c_noPU_prompt_type_v2_EE->cd();
  c_noPU_prompt_type_v2_EE->SetLeftMargin(0.12);
  TH1D* h_noPU_prompt_type_v2_EE_incl = new TH1D("h_noPU_prompt_type_v2_EE_incl","h_noPU_prompt_type_v2_EE_incl",4,0,4);
  h_noPU_prompt_type_v2_EE_incl->GetYaxis()->SetTitle("Counts");
  h_noPU_prompt_type_v2_EE_incl->SetTitle("");
  h_noPU_prompt_type_v2_EE_incl->SetMinimum(0.);
  h_noPU_prompt_type_v2_EE_incl->SetBinContent(1,h_noPU_prompt_type_v2_EE->GetBinContent(1));
  h_noPU_prompt_type_v2_EE_incl->SetBinContent(2,h_noPU_prompt_type_v2_EE->GetBinContent(2));
  h_noPU_prompt_type_v2_EE_incl->SetBinContent(3,h_noPU_prompt_type_v2_EE->GetBinContent(3));
  h_noPU_prompt_type_v2_EE_incl->SetBinContent(4,h_noPU_prompt_type_v2_EE->GetBinContent(4));
  h_noPU_prompt_type_v2_EE_incl->Draw("hist");
  c_noPU_prompt_type_v2_EE->Print("plots/trktype_v2_noPU_prompt_EE.pdf");
  cout << "Endcap prompt    |" << h_noPU_prompt_type_v2_EE->GetBinContent(1) << " : " << h_noPU_prompt_type_v2_EE->GetBinContent(2) << " : " << h_noPU_prompt_type_v2_EE->GetBinContent(3) << " : " << h_noPU_prompt_type_v2_EE->GetBinContent(4) << endl;
    // nonprompt
      // Barrel
  TCanvas* c_noPU_nonprompt_type_v2_EB = new TCanvas("c_noPU_nonprompt_type_v2_EB", "c_noPU_nonprompt_type_v2_EB", 1500, 1500);
  c_noPU_nonprompt_type_v2_EB->cd();
  c_noPU_nonprompt_type_v2_EB->SetLeftMargin(0.12);
  TH1D* h_noPU_nonprompt_type_v2_EB_incl = new TH1D("h_noPU_nonprompt_type_v2_EB_incl","h_noPU_nonprompt_type_v2_EB_incl",4,0,4);
  h_noPU_nonprompt_type_v2_EB_incl->GetYaxis()->SetTitle("Counts");
  h_noPU_nonprompt_type_v2_EB_incl->SetTitle("");
  h_noPU_nonprompt_type_v2_EB_incl->SetMinimum(0.);
  h_noPU_nonprompt_type_v2_EB_incl->SetBinContent(1,h_noPU_nonprompt_type_v2_EB->GetBinContent(1));
  h_noPU_nonprompt_type_v2_EB_incl->SetBinContent(2,h_noPU_nonprompt_type_v2_EB->GetBinContent(2));
  h_noPU_nonprompt_type_v2_EB_incl->SetBinContent(3,h_noPU_nonprompt_type_v2_EB->GetBinContent(3));
  h_noPU_nonprompt_type_v2_EB_incl->SetBinContent(4,h_noPU_nonprompt_type_v2_EB->GetBinContent(4));
  h_noPU_nonprompt_type_v2_EB_incl->Draw("hist");
  c_noPU_nonprompt_type_v2_EB->Print("plots/trktype_v2_noPU_nonprompt_EB.pdf");
  cout << "Barrel nonprompt |" << h_noPU_nonprompt_type_v2_EB->GetBinContent(1) << " : " << h_noPU_nonprompt_type_v2_EB->GetBinContent(2) << " : " << h_noPU_nonprompt_type_v2_EB->GetBinContent(3) << " : " << h_noPU_nonprompt_type_v2_EB->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_noPU_nonprompt_type_v2_EE = new TCanvas("c_noPU_nonprompt_type_v2_EE", "c_noPU_nonprompt_type_v2_EE", 1500, 1500);
  c_noPU_nonprompt_type_v2_EE->cd();
  c_noPU_nonprompt_type_v2_EE->SetLeftMargin(0.12);
  TH1D* h_noPU_nonprompt_type_v2_EE_incl = new TH1D("h_noPU_nonprompt_type_v2_EE_incl","h_noPU_nonprompt_type_v2_EE_incl",4,0,4);
  h_noPU_nonprompt_type_v2_EE_incl->GetYaxis()->SetTitle("Counts");
  h_noPU_nonprompt_type_v2_EE_incl->SetTitle("");
  h_noPU_nonprompt_type_v2_EE_incl->SetMinimum(0.);
  h_noPU_nonprompt_type_v2_EE_incl->SetBinContent(1,h_noPU_nonprompt_type_v2_EE->GetBinContent(1));
  h_noPU_nonprompt_type_v2_EE_incl->SetBinContent(2,h_noPU_nonprompt_type_v2_EE->GetBinContent(2));
  h_noPU_nonprompt_type_v2_EE_incl->SetBinContent(3,h_noPU_nonprompt_type_v2_EE->GetBinContent(3));
  h_noPU_nonprompt_type_v2_EE_incl->SetBinContent(4,h_noPU_nonprompt_type_v2_EE->GetBinContent(4));
  h_noPU_nonprompt_type_v2_EE_incl->Draw("hist");
  c_noPU_nonprompt_type_v2_EE->Print("plots/trktype_v2_noPU_nonprompt_EE.pdf");
  cout << "Endcap nonprompt |" << h_noPU_nonprompt_type_v2_EE->GetBinContent(1) << " : " << h_noPU_nonprompt_type_v2_EE->GetBinContent(2) << " : " << h_noPU_nonprompt_type_v2_EE->GetBinContent(3) << " : " << h_noPU_nonprompt_type_v2_EE->GetBinContent(4) << endl;
  cout << endl;


  ///////////
  /// vtx ///
  ///////////
  // PU200
    // prompt
      // Barrel
  TCanvas* c_PU200_prompt_type_v2_EB_vtx = new TCanvas("c_PU200_prompt_type_v2_EB_vtx", "c_PU200_prompt_type_v2_EB_vtx", 1500, 1500);
  c_PU200_prompt_type_v2_EB_vtx->cd();
  c_PU200_prompt_type_v2_EB_vtx->SetLeftMargin(0.12);
  TH1D* h_PU200_prompt_type_v2_EB_incl_vtx = new TH1D("h_PU200_prompt_type_v2_EB_incl_vtx","h_PU200_prompt_type_v2_EB_incl_vtx",4,0,4);
  h_PU200_prompt_type_v2_EB_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt_type_v2_EB_incl_vtx->SetTitle("");
  h_PU200_prompt_type_v2_EB_incl_vtx->SetMinimum(0.);
  h_PU200_prompt_type_v2_EB_incl_vtx->SetBinContent(1,h_PU200_prompt_type_v2_EB_vtx->GetBinContent(1));
  h_PU200_prompt_type_v2_EB_incl_vtx->SetBinContent(2,h_PU200_prompt_type_v2_EB_vtx->GetBinContent(2));
  h_PU200_prompt_type_v2_EB_incl_vtx->SetBinContent(3,h_PU200_prompt_type_v2_EB_vtx->GetBinContent(3));
  h_PU200_prompt_type_v2_EB_incl_vtx->SetBinContent(4,h_PU200_prompt_type_v2_EB_vtx->GetBinContent(4));
  h_PU200_prompt_type_v2_EB_incl_vtx->Draw("hist");
  c_PU200_prompt_type_v2_EB_vtx->Print("plots/trktype_v2_PU200_prompt_EB_vtx.pdf");
  cout << "PU200" << endl;
  cout << "Barrel prompt    |" << h_PU200_prompt_type_v2_EB_vtx->GetBinContent(1) << " : " << h_PU200_prompt_type_v2_EB_vtx->GetBinContent(2) << " : " << h_PU200_prompt_type_v2_EB_vtx->GetBinContent(3) << " : " << h_PU200_prompt_type_v2_EB_vtx->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_PU200_prompt_type_v2_EE_vtx = new TCanvas("c_PU200_prompt_type_v2_EE_vtx", "c_PU200_prompt_type_v2_EE_vtx", 1500, 1500);
  c_PU200_prompt_type_v2_EE_vtx->cd();
  c_PU200_prompt_type_v2_EE_vtx->SetLeftMargin(0.12);
  TH1D* h_PU200_prompt_type_v2_EE_incl_vtx = new TH1D("h_PU200_prompt_type_v2_EE_incl_vtx","h_PU200_prompt_type_v2_EE_incl_vtx",4,0,4);
  h_PU200_prompt_type_v2_EE_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt_type_v2_EE_incl_vtx->SetTitle("");
  h_PU200_prompt_type_v2_EE_incl_vtx->SetMinimum(0.);
  h_PU200_prompt_type_v2_EE_incl_vtx->SetBinContent(1,h_PU200_prompt_type_v2_EE_vtx->GetBinContent(1));
  h_PU200_prompt_type_v2_EE_incl_vtx->SetBinContent(2,h_PU200_prompt_type_v2_EE_vtx->GetBinContent(2));
  h_PU200_prompt_type_v2_EE_incl_vtx->SetBinContent(3,h_PU200_prompt_type_v2_EE_vtx->GetBinContent(3));
  h_PU200_prompt_type_v2_EE_incl_vtx->SetBinContent(4,h_PU200_prompt_type_v2_EE_vtx->GetBinContent(4));
  h_PU200_prompt_type_v2_EE_incl_vtx->Draw("hist");
  c_PU200_prompt_type_v2_EE_vtx->Print("plots/trktype_v2_PU200_prompt_EE_vtx.pdf");
  cout << "Endcap prompt    |" << h_PU200_prompt_type_v2_EE_vtx->GetBinContent(1) << " : " << h_PU200_prompt_type_v2_EE_vtx->GetBinContent(2) << " : " << h_PU200_prompt_type_v2_EE_vtx->GetBinContent(3) << " : " << h_PU200_prompt_type_v2_EE_vtx->GetBinContent(4) << endl;
    // nonprompt
      // Barrel
  TCanvas* c_PU200_nonprompt_type_v2_EB_vtx = new TCanvas("c_PU200_nonprompt_type_v2_EB_vtx", "c_PU200_nonprompt_type_v2_EB_vtx", 1500, 1500);
  c_PU200_nonprompt_type_v2_EB_vtx->cd();
  c_PU200_nonprompt_type_v2_EB_vtx->SetLeftMargin(0.12);
  TH1D* h_PU200_nonprompt_type_v2_EB_incl_vtx = new TH1D("h_PU200_nonprompt_type_v2_EB_incl_vtx","h_PU200_nonprompt_type_v2_EB_incl_vtx",4,0,4);
  h_PU200_nonprompt_type_v2_EB_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_PU200_nonprompt_type_v2_EB_incl_vtx->SetTitle("");
  h_PU200_nonprompt_type_v2_EB_incl_vtx->SetMinimum(0.);
  h_PU200_nonprompt_type_v2_EB_incl_vtx->SetBinContent(1,h_PU200_nonprompt_type_v2_EB_vtx->GetBinContent(1));
  h_PU200_nonprompt_type_v2_EB_incl_vtx->SetBinContent(2,h_PU200_nonprompt_type_v2_EB_vtx->GetBinContent(2));
  h_PU200_nonprompt_type_v2_EB_incl_vtx->SetBinContent(3,h_PU200_nonprompt_type_v2_EB_vtx->GetBinContent(3));
  h_PU200_nonprompt_type_v2_EB_incl_vtx->SetBinContent(4,h_PU200_nonprompt_type_v2_EB_vtx->GetBinContent(4));
  h_PU200_nonprompt_type_v2_EB_incl_vtx->Draw("hist");
  c_PU200_nonprompt_type_v2_EB_vtx->Print("plots/trktype_v2_PU200_nonprompt_EB_vtx.pdf");
  cout << "Barrel nonprompt |" << h_PU200_nonprompt_type_v2_EB_vtx->GetBinContent(1) << " : " << h_PU200_nonprompt_type_v2_EB_vtx->GetBinContent(2) << " : " << h_PU200_nonprompt_type_v2_EB_vtx->GetBinContent(3) << " : " << h_PU200_nonprompt_type_v2_EB_vtx->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_PU200_nonprompt_type_v2_EE_vtx = new TCanvas("c_PU200_nonprompt_type_v2_EE_vtx", "c_PU200_nonprompt_type_v2_EE_vtx", 1500, 1500);
  c_PU200_nonprompt_type_v2_EE_vtx->cd();
  c_PU200_nonprompt_type_v2_EE_vtx->SetLeftMargin(0.12);
  TH1D* h_PU200_nonprompt_type_v2_EE_incl_vtx = new TH1D("h_PU200_nonprompt_type_v2_EE_incl_vtx","h_PU200_nonprompt_type_v2_EE_incl_vtx",4,0,4);
  h_PU200_nonprompt_type_v2_EE_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_PU200_nonprompt_type_v2_EE_incl_vtx->SetTitle("");
  h_PU200_nonprompt_type_v2_EE_incl_vtx->SetMinimum(0.);
  h_PU200_nonprompt_type_v2_EE_incl_vtx->SetBinContent(1,h_PU200_nonprompt_type_v2_EE_vtx->GetBinContent(1));
  h_PU200_nonprompt_type_v2_EE_incl_vtx->SetBinContent(2,h_PU200_nonprompt_type_v2_EE_vtx->GetBinContent(2));
  h_PU200_nonprompt_type_v2_EE_incl_vtx->SetBinContent(3,h_PU200_nonprompt_type_v2_EE_vtx->GetBinContent(3));
  h_PU200_nonprompt_type_v2_EE_incl_vtx->SetBinContent(4,h_PU200_nonprompt_type_v2_EE_vtx->GetBinContent(4));
  h_PU200_nonprompt_type_v2_EE_incl_vtx->Draw("hist");
  c_PU200_nonprompt_type_v2_EE_vtx->Print("plots/trktype_v2_PU200_nonprompt_EE_vtx.pdf");
  cout << "Endcap nonprompt |" << h_PU200_nonprompt_type_v2_EE_vtx->GetBinContent(1) << " : " << h_PU200_nonprompt_type_v2_EE_vtx->GetBinContent(2) << " : " << h_PU200_nonprompt_type_v2_EE_vtx->GetBinContent(3) << " : " << h_PU200_nonprompt_type_v2_EE_vtx->GetBinContent(4) << endl;
  // noPU
    // prompt
      // Barrel
  TCanvas* c_noPU_prompt_type_v2_EB_vtx = new TCanvas("c_noPU_prompt_type_v2_EB_vtx", "c_noPU_prompt_type_v2_EB_vtx", 1500, 1500);
  c_noPU_prompt_type_v2_EB_vtx->cd();
  c_noPU_prompt_type_v2_EB_vtx->SetLeftMargin(0.12);
  TH1D* h_noPU_prompt_type_v2_EB_incl_vtx = new TH1D("h_noPU_prompt_type_v2_EB_incl_vtx","h_noPU_prompt_type_v2_EB_incl_vtx",4,0,4);
  h_noPU_prompt_type_v2_EB_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_noPU_prompt_type_v2_EB_incl_vtx->SetTitle("");
  h_noPU_prompt_type_v2_EB_incl_vtx->SetMinimum(0.);
  h_noPU_prompt_type_v2_EB_incl_vtx->SetBinContent(1,h_noPU_prompt_type_v2_EB_vtx->GetBinContent(1));
  h_noPU_prompt_type_v2_EB_incl_vtx->SetBinContent(2,h_noPU_prompt_type_v2_EB_vtx->GetBinContent(2));
  h_noPU_prompt_type_v2_EB_incl_vtx->SetBinContent(3,h_noPU_prompt_type_v2_EB_vtx->GetBinContent(3));
  h_noPU_prompt_type_v2_EB_incl_vtx->SetBinContent(4,h_noPU_prompt_type_v2_EB_vtx->GetBinContent(4));
  h_noPU_prompt_type_v2_EB_incl_vtx->Draw("hist");
  c_noPU_prompt_type_v2_EB_vtx->Print("plots/trktype_v2_noPU_prompt_EB_vtx.pdf");
  cout << endl;
  cout << "noPU" << endl;
  cout << "Barrel prompt    |" << h_noPU_prompt_type_v2_EB_vtx->GetBinContent(1) << " : " << h_noPU_prompt_type_v2_EB_vtx->GetBinContent(2) << " : " << h_noPU_prompt_type_v2_EB_vtx->GetBinContent(3) << " : " << h_noPU_prompt_type_v2_EB_vtx->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_noPU_prompt_type_v2_EE_vtx = new TCanvas("c_noPU_prompt_type_v2_EE_vtx", "c_noPU_prompt_type_v2_EE_vtx", 1500, 1500);
  c_noPU_prompt_type_v2_EE_vtx->cd();
  c_noPU_prompt_type_v2_EE_vtx->SetLeftMargin(0.12);
  TH1D* h_noPU_prompt_type_v2_EE_incl_vtx = new TH1D("h_noPU_prompt_type_v2_EE_incl_vtx","h_noPU_prompt_type_v2_EE_incl_vtx",4,0,4);
  h_noPU_prompt_type_v2_EE_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_noPU_prompt_type_v2_EE_incl_vtx->SetTitle("");
  h_noPU_prompt_type_v2_EE_incl_vtx->SetMinimum(0.);
  h_noPU_prompt_type_v2_EE_incl_vtx->SetBinContent(1,h_noPU_prompt_type_v2_EE_vtx->GetBinContent(1));
  h_noPU_prompt_type_v2_EE_incl_vtx->SetBinContent(2,h_noPU_prompt_type_v2_EE_vtx->GetBinContent(2));
  h_noPU_prompt_type_v2_EE_incl_vtx->SetBinContent(3,h_noPU_prompt_type_v2_EE_vtx->GetBinContent(3));
  h_noPU_prompt_type_v2_EE_incl_vtx->SetBinContent(4,h_noPU_prompt_type_v2_EE_vtx->GetBinContent(4));
  h_noPU_prompt_type_v2_EE_incl_vtx->Draw("hist");
  c_noPU_prompt_type_v2_EE_vtx->Print("plots/trktype_v2_noPU_prompt_EE_vtx.pdf");
  cout << "Endcap prompt    |" << h_noPU_prompt_type_v2_EE_vtx->GetBinContent(1) << " : " << h_noPU_prompt_type_v2_EE_vtx->GetBinContent(2) << " : " << h_noPU_prompt_type_v2_EE_vtx->GetBinContent(3) << " : " << h_noPU_prompt_type_v2_EE_vtx->GetBinContent(4) << endl;
    // nonprompt
      // Barrel
  TCanvas* c_noPU_nonprompt_type_v2_EB_vtx = new TCanvas("c_noPU_nonprompt_type_v2_EB_vtx", "c_noPU_nonprompt_type_v2_EB_vtx", 1500, 1500);
  c_noPU_nonprompt_type_v2_EB_vtx->cd();
  c_noPU_nonprompt_type_v2_EB_vtx->SetLeftMargin(0.12);
  TH1D* h_noPU_nonprompt_type_v2_EB_incl_vtx = new TH1D("h_noPU_nonprompt_type_v2_EB_incl_vtx","h_noPU_nonprompt_type_v2_EB_incl_vtx",4,0,4);
  h_noPU_nonprompt_type_v2_EB_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_noPU_nonprompt_type_v2_EB_incl_vtx->SetTitle("");
  h_noPU_nonprompt_type_v2_EB_incl_vtx->SetMinimum(0.);
  h_noPU_nonprompt_type_v2_EB_incl_vtx->SetBinContent(1,h_noPU_nonprompt_type_v2_EB_vtx->GetBinContent(1));
  h_noPU_nonprompt_type_v2_EB_incl_vtx->SetBinContent(2,h_noPU_nonprompt_type_v2_EB_vtx->GetBinContent(2));
  h_noPU_nonprompt_type_v2_EB_incl_vtx->SetBinContent(3,h_noPU_nonprompt_type_v2_EB_vtx->GetBinContent(3));
  h_noPU_nonprompt_type_v2_EB_incl_vtx->SetBinContent(4,h_noPU_nonprompt_type_v2_EB_vtx->GetBinContent(4));
  h_noPU_nonprompt_type_v2_EB_incl_vtx->Draw("hist");
  c_noPU_nonprompt_type_v2_EB_vtx->Print("plots/trktype_v2_noPU_nonprompt_EB_vtx.pdf");
  cout << "Barrel nonprompt |" << h_noPU_nonprompt_type_v2_EB_vtx->GetBinContent(1) << " : " << h_noPU_nonprompt_type_v2_EB_vtx->GetBinContent(2) << " : " << h_noPU_nonprompt_type_v2_EB_vtx->GetBinContent(3) << " : " << h_noPU_nonprompt_type_v2_EB_vtx->GetBinContent(4) << endl;
      // Endcap
  TCanvas* c_noPU_nonprompt_type_v2_EE_vtx = new TCanvas("c_noPU_nonprompt_type_v2_EE_vtx", "c_noPU_nonprompt_type_v2_EE_vtx", 1500, 1500);
  c_noPU_nonprompt_type_v2_EE_vtx->cd();
  c_noPU_nonprompt_type_v2_EE_vtx->SetLeftMargin(0.12);
  TH1D* h_noPU_nonprompt_type_v2_EE_incl_vtx = new TH1D("h_noPU_nonprompt_type_v2_EE_incl_vtx","h_noPU_nonprompt_type_v2_EE_incl_vtx",4,0,4);
  h_noPU_nonprompt_type_v2_EE_incl_vtx->GetYaxis()->SetTitle("Counts");
  h_noPU_nonprompt_type_v2_EE_incl_vtx->SetTitle("");
  h_noPU_nonprompt_type_v2_EE_incl_vtx->SetMinimum(0.);
  h_noPU_nonprompt_type_v2_EE_incl_vtx->SetBinContent(1,h_noPU_nonprompt_type_v2_EE_vtx->GetBinContent(1));
  h_noPU_nonprompt_type_v2_EE_incl_vtx->SetBinContent(2,h_noPU_nonprompt_type_v2_EE_vtx->GetBinContent(2));
  h_noPU_nonprompt_type_v2_EE_incl_vtx->SetBinContent(3,h_noPU_nonprompt_type_v2_EE_vtx->GetBinContent(3));
  h_noPU_nonprompt_type_v2_EE_incl_vtx->SetBinContent(4,h_noPU_nonprompt_type_v2_EE_vtx->GetBinContent(4));
  h_noPU_nonprompt_type_v2_EE_incl_vtx->Draw("hist");
  c_noPU_nonprompt_type_v2_EE_vtx->Print("plots/trktype_v2_noPU_nonprompt_EE_vtx.pdf");
  cout << "Endcap nonprompt |" << h_noPU_nonprompt_type_v2_EE_vtx->GetBinContent(1) << " : " << h_noPU_nonprompt_type_v2_EE_vtx->GetBinContent(2) << " : " << h_noPU_nonprompt_type_v2_EE_vtx->GetBinContent(3) << " : " << h_noPU_nonprompt_type_v2_EE_vtx->GetBinContent(4) << endl;
  cout << endl;

}

void track_sigma_type_v2() {
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


  // (muon, track)
    // PU200_prompt
  TH1D* h_PU200_prompt_PV_EB 	  = new TH1D("h_PU200_prompt_PV_EB", 	  "h_PU200_prompt_PV_EB",      10, 0, 10);
  TH1D* h_PU200_prompt_SV_EB 	  = new TH1D("h_PU200_prompt_SV_EB", 	  "h_PU200_prompt_SV_EB",      10, 0, 10);
  TH1D* h_PU200_prompt_PU_EB	  = new TH1D("h_PU200_prompt_PU_EB", 	  "h_PU200_prompt_PU_EB",      10, 0, 10);
  TH1D* h_PU200_prompt_fake_EB 	  = new TH1D("h_PU200_prompt_fake_EB",    "h_PU200_prompt_fake_EB",    10, 0, 10);
  TH1D* h_PU200_prompt_no_tErr_EB = new TH1D("h_PU200_prompt_no_tErr_EB", "h_PU200_prompt_no_tErr_EB", 10, 0, 10);
  TH1D* h_PU200_prompt_tot_reco2sim_EB     = new TH1D("h_PU200_prompt_tot_reco2sim_EB",     "h_PU200_prompt_tot_reco2sim_EB",     10, 0, 10);
  TH1D* h_PU200_prompt_tot_EB     = new TH1D("h_PU200_prompt_tot_EB",     "h_PU200_prompt_tot_EB",     10, 0, 10);
  TH1D* h_PU200_prompt_PV_EE 	  = new TH1D("h_PU200_prompt_PV_EE", 	  "h_PU200_prompt_PV_EE",      10, 0, 10);
  TH1D* h_PU200_prompt_SV_EE 	  = new TH1D("h_PU200_prompt_SV_EE", 	  "h_PU200_prompt_SV_EE",      10, 0, 10);
  TH1D* h_PU200_prompt_PU_EE	  = new TH1D("h_PU200_prompt_PU_EE", 	  "h_PU200_prompt_PU_EE",      10, 0, 10);
  TH1D* h_PU200_prompt_fake_EE 	  = new TH1D("h_PU200_prompt_fake_EE",    "h_PU200_prompt_fake_EE",    10, 0, 10);
  TH1D* h_PU200_prompt_no_tErr_EE = new TH1D("h_PU200_prompt_no_tErr_EE", "h_PU200_prompt_no_tErr_EE", 10, 0, 10);
  TH1D* h_PU200_prompt_tot_reco2sim_EE     = new TH1D("h_PU200_prompt_tot_reco2sim_EE",     "h_PU200_prompt_tot_reco2sim_EE",     10, 0, 10);
  TH1D* h_PU200_prompt_tot_EE     = new TH1D("h_PU200_prompt_tot_EE",     "h_PU200_prompt_tot_EE",     10, 0, 10);
    // PU200_nonprompt
  TH1D* h_PU200_nonprompt_PV_EB      = new TH1D("h_PU200_nonprompt_PV_EB", 	"h_PU200_nonprompt_PV_EB",      10, 0, 10);
  TH1D* h_PU200_nonprompt_SV_EB      = new TH1D("h_PU200_nonprompt_SV_EB", 	"h_PU200_nonprompt_SV_EB",      10, 0, 10);
  TH1D* h_PU200_nonprompt_PU_EB	     = new TH1D("h_PU200_nonprompt_PU_EB", 	"h_PU200_nonprompt_PU_EB",      10, 0, 10);
  TH1D* h_PU200_nonprompt_fake_EB    = new TH1D("h_PU200_nonprompt_fake_EB",    "h_PU200_nonprompt_fake_EB",    10, 0, 10);
  TH1D* h_PU200_nonprompt_no_tErr_EB = new TH1D("h_PU200_nonprompt_no_tErr_EB", "h_PU200_nonprompt_no_tErr_EB", 10, 0, 10);
  TH1D* h_PU200_nonprompt_tot_reco2sim_EB     = new TH1D("h_PU200_nonprompt_tot_reco2sim_EB",     "h_PU200_nonprompt_tot_reco2sim_EB",     10, 0, 10);
  TH1D* h_PU200_nonprompt_tot_EB     = new TH1D("h_PU200_nonprompt_tot_EB",     "h_PU200_nonprompt_tot_EB",     10, 0, 10);
  TH1D* h_PU200_nonprompt_PV_EE      = new TH1D("h_PU200_nonprompt_PV_EE", 	"h_PU200_nonprompt_PV_EE",      10, 0, 10);
  TH1D* h_PU200_nonprompt_SV_EE      = new TH1D("h_PU200_nonprompt_SV_EE", 	"h_PU200_nonprompt_SV_EE",      10, 0, 10);
  TH1D* h_PU200_nonprompt_PU_EE	     = new TH1D("h_PU200_nonprompt_PU_EE", 	"h_PU200_nonprompt_PU_EE",      10, 0, 10);
  TH1D* h_PU200_nonprompt_fake_EE    = new TH1D("h_PU200_nonprompt_fake_EE",    "h_PU200_nonprompt_fake_EE",    10, 0, 10);
  TH1D* h_PU200_nonprompt_no_tErr_EE = new TH1D("h_PU200_nonprompt_no_tErr_EE", "h_PU200_nonprompt_no_tErr_EE", 10, 0, 10);
  TH1D* h_PU200_nonprompt_tot_reco2sim_EE     = new TH1D("h_PU200_nonprompt_tot_reco2sim_EE",     "h_PU200_nonprompt_tot_reco2sim_EE",     10, 0, 10);
  TH1D* h_PU200_nonprompt_tot_EE     = new TH1D("h_PU200_nonprompt_tot_EE",     "h_PU200_nonprompt_tot_EE",     10, 0, 10);
    // noPU_prompt
  TH1D* h_noPU_prompt_PV_EB 	  = new TH1D("h_noPU_prompt_PV_EB", 	  "h_noPU_prompt_PV_EB",      10, 0, 10);
  TH1D* h_noPU_prompt_SV_EB 	  = new TH1D("h_noPU_prompt_SV_EB", 	  "h_noPU_prompt_SV_EB",      10, 0, 10);
  TH1D* h_noPU_prompt_PU_EB	  = new TH1D("h_noPU_prompt_PU_EB", 	  "h_noPU_prompt_PU_EB",      10, 0, 10);
  TH1D* h_noPU_prompt_fake_EB 	  = new TH1D("h_noPU_prompt_fake_EB",     "h_noPU_prompt_fake_EB",    10, 0, 10);
  TH1D* h_noPU_prompt_no_tErr_EB  = new TH1D("h_noPU_prompt_no_tErr_EB",  "h_noPU_prompt_no_tErr_EB", 10, 0, 10);
  TH1D* h_noPU_prompt_tot_reco2sim_EB     = new TH1D("h_noPU_prompt_tot_reco2sim_EB",     "h_noPU_prompt_tot_reco2sim_EB",     10, 0, 10);
  TH1D* h_noPU_prompt_tot_EB     = new TH1D("h_noPU_prompt_tot_EB",     "h_noPU_prompt_tot_EB",     10, 0, 10);
  TH1D* h_noPU_prompt_PV_EE 	  = new TH1D("h_noPU_prompt_PV_EE", 	  "h_noPU_prompt_PV_EE",      10, 0, 10);
  TH1D* h_noPU_prompt_SV_EE 	  = new TH1D("h_noPU_prompt_SV_EE", 	  "h_noPU_prompt_SV_EE",      10, 0, 10);
  TH1D* h_noPU_prompt_PU_EE	  = new TH1D("h_noPU_prompt_PU_EE", 	  "h_noPU_prompt_PU_EE",      10, 0, 10);
  TH1D* h_noPU_prompt_fake_EE 	  = new TH1D("h_noPU_prompt_fake_EE",     "h_noPU_prompt_fake_EE",    10, 0, 10);
  TH1D* h_noPU_prompt_no_tErr_EE  = new TH1D("h_noPU_prompt_no_tErr_EE",  "h_noPU_prompt_no_tErr_EE", 10, 0, 10);
  TH1D* h_noPU_prompt_tot_reco2sim_EE     = new TH1D("h_noPU_prompt_tot_reco2sim_EE",     "h_noPU_prompt_tot_reco2sim_EE",     10, 0, 10);
  TH1D* h_noPU_prompt_tot_EE     = new TH1D("h_noPU_prompt_tot_EE",     "h_noPU_prompt_tot_EE",     10, 0, 10);
    // noPU_nonprompt
  TH1D* h_noPU_nonprompt_PV_EB      = new TH1D("h_noPU_nonprompt_PV_EB", 	"h_noPU_nonprompt_PV_EB",      10, 0, 10);
  TH1D* h_noPU_nonprompt_SV_EB      = new TH1D("h_noPU_nonprompt_SV_EB", 	"h_noPU_nonprompt_SV_EB",      10, 0, 10);
  TH1D* h_noPU_nonprompt_PU_EB	    = new TH1D("h_noPU_nonprompt_PU_EB", 	"h_noPU_nonprompt_PU_EB",      10, 0, 10);
  TH1D* h_noPU_nonprompt_fake_EB    = new TH1D("h_noPU_nonprompt_fake_EB",      "h_noPU_nonprompt_fake_EB",    10, 0, 10);
  TH1D* h_noPU_nonprompt_no_tErr_EB = new TH1D("h_noPU_nonprompt_no_tErr_EB",   "h_noPU_nonprompt_no_tErr_EB", 10, 0, 10);
  TH1D* h_noPU_nonprompt_tot_reco2sim_EB     = new TH1D("h_noPU_nonprompt_tot_reco2sim_EB",     "h_noPU_nonprompt_tot_reco2sim_EB",     10, 0, 10);
  TH1D* h_noPU_nonprompt_tot_EB     = new TH1D("h_noPU_nonprompt_tot_EB",     "h_noPU_nonprompt_tot_EB",     10, 0, 10);
  TH1D* h_noPU_nonprompt_PV_EE      = new TH1D("h_noPU_nonprompt_PV_EE", 	"h_noPU_nonprompt_PV_EE",      10, 0, 10);
  TH1D* h_noPU_nonprompt_SV_EE      = new TH1D("h_noPU_nonprompt_SV_EE", 	"h_noPU_nonprompt_SV_EE",      10, 0, 10);
  TH1D* h_noPU_nonprompt_PU_EE	    = new TH1D("h_noPU_nonprompt_PU_EE", 	"h_noPU_nonprompt_PU_EE",      10, 0, 10);
  TH1D* h_noPU_nonprompt_fake_EE    = new TH1D("h_noPU_nonprompt_fake_EE",      "h_noPU_nonprompt_fake_EE",    10, 0, 10);
  TH1D* h_noPU_nonprompt_no_tErr_EE = new TH1D("h_noPU_nonprompt_no_tErr_EE",   "h_noPU_nonprompt_no_tErr_EE", 10, 0, 10);
  TH1D* h_noPU_nonprompt_tot_reco2sim_EE     = new TH1D("h_noPU_nonprompt_tot_reco2sim_EE",     "h_noPU_nonprompt_tot_reco2sim_EE",     10, 0, 10);
  TH1D* h_noPU_nonprompt_tot_EE     = new TH1D("h_noPU_nonprompt_tot_EE",     "h_noPU_nonprompt_tot_EE",     10, 0, 10);

  // (PV, track)
    // PU200_prompt
  TH1D* h_PU200_prompt_PV_EB_vtx      = new TH1D("h_PU200_prompt_PV_EB_vtx", 	  "h_PU200_prompt_PV_EB_vtx",      10, 0, 10);
  TH1D* h_PU200_prompt_SV_EB_vtx      = new TH1D("h_PU200_prompt_SV_EB_vtx", 	  "h_PU200_prompt_SV_EB_vtx",      10, 0, 10);
  TH1D* h_PU200_prompt_PU_EB_vtx      = new TH1D("h_PU200_prompt_PU_EB_vtx", 	  "h_PU200_prompt_PU_EB_vtx",      10, 0, 10);
  TH1D* h_PU200_prompt_fake_EB_vtx    = new TH1D("h_PU200_prompt_fake_EB_vtx",    "h_PU200_prompt_fake_EB_vtx",    10, 0, 10);
  TH1D* h_PU200_prompt_no_tErr_EB_vtx = new TH1D("h_PU200_prompt_no_tErr_EB_vtx", "h_PU200_prompt_no_tErr_EB_vtx", 10, 0, 10);
  TH1D* h_PU200_prompt_tot_reco2sim_EB_vtx     = new TH1D("h_PU200_prompt_tot_reco2sim_EB_vtx",     "h_PU200_prompt_tot_reco2sim_EB_vtx",     10, 0, 10);
  TH1D* h_PU200_prompt_tot_EB_vtx     = new TH1D("h_PU200_prompt_tot_EB_vtx",     "h_PU200_prompt_tot_EB_vtx",     10, 0, 10);
  TH1D* h_PU200_prompt_PV_EE_vtx      = new TH1D("h_PU200_prompt_PV_EE_vtx", 	  "h_PU200_prompt_PV_EE_vtx",      10, 0, 10);
  TH1D* h_PU200_prompt_SV_EE_vtx      = new TH1D("h_PU200_prompt_SV_EE_vtx", 	  "h_PU200_prompt_SV_EE_vtx",      10, 0, 10);
  TH1D* h_PU200_prompt_PU_EE_vtx      = new TH1D("h_PU200_prompt_PU_EE_vtx", 	  "h_PU200_prompt_PU_EE_vtx",      10, 0, 10);
  TH1D* h_PU200_prompt_fake_EE_vtx    = new TH1D("h_PU200_prompt_fake_EE_vtx",    "h_PU200_prompt_fake_EE_vtx",    10, 0, 10);
  TH1D* h_PU200_prompt_no_tErr_EE_vtx = new TH1D("h_PU200_prompt_no_tErr_EE_vtx", "h_PU200_prompt_no_tErr_EE_vtx", 10, 0, 10);
  TH1D* h_PU200_prompt_tot_reco2sim_EE_vtx     = new TH1D("h_PU200_prompt_tot_reco2sim_EE_vtx",     "h_PU200_prompt_tot_reco2sim_EE_vtx",     10, 0, 10);
  TH1D* h_PU200_prompt_tot_EE_vtx     = new TH1D("h_PU200_prompt_tot_EE_vtx",     "h_PU200_prompt_tot_EE_vtx",     10, 0, 10);
    // PU200_nonprompt
  TH1D* h_PU200_nonprompt_PV_EB_vtx      = new TH1D("h_PU200_nonprompt_PV_EB_vtx", 	"h_PU200_nonprompt_PV_EB_vtx",      10, 0, 10);
  TH1D* h_PU200_nonprompt_SV_EB_vtx      = new TH1D("h_PU200_nonprompt_SV_EB_vtx", 	"h_PU200_nonprompt_SV_EB_vtx",      10, 0, 10);
  TH1D* h_PU200_nonprompt_PU_EB_vtx	 = new TH1D("h_PU200_nonprompt_PU_EB_vtx", 	"h_PU200_nonprompt_PU_EB_vtx",      10, 0, 10);
  TH1D* h_PU200_nonprompt_fake_EB_vtx    = new TH1D("h_PU200_nonprompt_fake_EB_vtx",    "h_PU200_nonprompt_fake_EB_vtx",    10, 0, 10);
  TH1D* h_PU200_nonprompt_no_tErr_EB_vtx = new TH1D("h_PU200_nonprompt_no_tErr_EB_vtx", "h_PU200_nonprompt_no_tErr_EB_vtx", 10, 0, 10);
  TH1D* h_PU200_nonprompt_tot_reco2sim_EB_vtx     = new TH1D("h_PU200_nonprompt_tot_reco2sim_EB_vtx",     "h_PU200_nonprompt_tot_reco2sim_EB_vtx",     10, 0, 10);
  TH1D* h_PU200_nonprompt_tot_EB_vtx     = new TH1D("h_PU200_nonprompt_tot_EB_vtx",     "h_PU200_nonprompt_tot_EB_vtx",     10, 0, 10);
  TH1D* h_PU200_nonprompt_PV_EE_vtx      = new TH1D("h_PU200_nonprompt_PV_EE_vtx", 	"h_PU200_nonprompt_PV_EE_vtx",      10, 0, 10);
  TH1D* h_PU200_nonprompt_SV_EE_vtx      = new TH1D("h_PU200_nonprompt_SV_EE_vtx", 	"h_PU200_nonprompt_SV_EE_vtx",      10, 0, 10);
  TH1D* h_PU200_nonprompt_PU_EE_vtx	 = new TH1D("h_PU200_nonprompt_PU_EE_vtx", 	"h_PU200_nonprompt_PU_EE_vtx",      10, 0, 10);
  TH1D* h_PU200_nonprompt_fake_EE_vtx    = new TH1D("h_PU200_nonprompt_fake_EE_vtx",    "h_PU200_nonprompt_fake_EE_vtx",    10, 0, 10);
  TH1D* h_PU200_nonprompt_no_tErr_EE_vtx = new TH1D("h_PU200_nonprompt_no_tErr_EE_vtx", "h_PU200_nonprompt_no_tErr_EE_vtx", 10, 0, 10);
  TH1D* h_PU200_nonprompt_tot_reco2sim_EE_vtx     = new TH1D("h_PU200_nonprompt_tot_reco2sim_EE_vtx",     "h_PU200_nonprompt_tot_reco2sim_EE_vtx",     10, 0, 10);
  TH1D* h_PU200_nonprompt_tot_EE_vtx     = new TH1D("h_PU200_nonprompt_tot_EE_vtx",     "h_PU200_nonprompt_tot_EE_vtx",     10, 0, 10);
    // noPU_prompt
  TH1D* h_noPU_prompt_PV_EB_vtx      = new TH1D("h_noPU_prompt_PV_EB_vtx", 	 "h_noPU_prompt_PV_EB_vtx",      10, 0, 10);
  TH1D* h_noPU_prompt_SV_EB_vtx      = new TH1D("h_noPU_prompt_SV_EB_vtx", 	 "h_noPU_prompt_SV_EB_vtx",      10, 0, 10);
  TH1D* h_noPU_prompt_PU_EB_vtx	     = new TH1D("h_noPU_prompt_PU_EB_vtx", 	 "h_noPU_prompt_PU_EB_vtx",      10, 0, 10);
  TH1D* h_noPU_prompt_fake_EB_vtx    = new TH1D("h_noPU_prompt_fake_EB_vtx",     "h_noPU_prompt_fake_EB_vtx",    10, 0, 10);
  TH1D* h_noPU_prompt_no_tErr_EB_vtx = new TH1D("h_noPU_prompt_no_tErr_EB_vtx",  "h_noPU_prompt_no_tErr_EB_vtx", 10, 0, 10);
  TH1D* h_noPU_prompt_tot_reco2sim_EB_vtx     = new TH1D("h_noPU_prompt_tot_reco2sim_EB_vtx",     "h_noPU_prompt_tot_reco2sim_EB_vtx",     10, 0, 10);
  TH1D* h_noPU_prompt_tot_EB_vtx     = new TH1D("h_noPU_prompt_tot_EB_vtx",     "h_noPU_prompt_tot_EB_vtx",     10, 0, 10);
  TH1D* h_noPU_prompt_PV_EE_vtx      = new TH1D("h_noPU_prompt_PV_EE_vtx", 	 "h_noPU_prompt_PV_EE_vtx",      10, 0, 10);
  TH1D* h_noPU_prompt_SV_EE_vtx      = new TH1D("h_noPU_prompt_SV_EE_vtx", 	 "h_noPU_prompt_SV_EE_vtx",      10, 0, 10);
  TH1D* h_noPU_prompt_PU_EE_vtx	     = new TH1D("h_noPU_prompt_PU_EE_vtx", 	 "h_noPU_prompt_PU_EE_vtx",      10, 0, 10);
  TH1D* h_noPU_prompt_fake_EE_vtx    = new TH1D("h_noPU_prompt_fake_EE_vtx",     "h_noPU_prompt_fake_EE_vtx",    10, 0, 10);
  TH1D* h_noPU_prompt_no_tErr_EE_vtx = new TH1D("h_noPU_prompt_no_tErr_EE_vtx",  "h_noPU_prompt_no_tErr_EE_vtx", 10, 0, 10);
  TH1D* h_noPU_prompt_tot_reco2sim_EE_vtx     = new TH1D("h_noPU_prompt_tot_reco2sim_EE_vtx",     "h_noPU_prompt_tot_reco2sim_EE_vtx",     10, 0, 10);
  TH1D* h_noPU_prompt_tot_EE_vtx     = new TH1D("h_noPU_prompt_tot_EE_vtx",     "h_noPU_prompt_tot_EE_vtx",     10, 0, 10);
    // noPU_nonprompt
  TH1D* h_noPU_nonprompt_PV_EB_vtx      = new TH1D("h_noPU_nonprompt_PV_EB_vtx",      "h_noPU_nonprompt_PV_EB_vtx",      10, 0, 10);
  TH1D* h_noPU_nonprompt_SV_EB_vtx      = new TH1D("h_noPU_nonprompt_SV_EB_vtx",      "h_noPU_nonprompt_SV_EB_vtx",      10, 0, 10);
  TH1D* h_noPU_nonprompt_PU_EB_vtx	= new TH1D("h_noPU_nonprompt_PU_EB_vtx",      "h_noPU_nonprompt_PU_EB_vtx",      10, 0, 10);
  TH1D* h_noPU_nonprompt_fake_EB_vtx    = new TH1D("h_noPU_nonprompt_fake_EB_vtx",    "h_noPU_nonprompt_fake_EB_vtx",    10, 0, 10);
  TH1D* h_noPU_nonprompt_no_tErr_EB_vtx = new TH1D("h_noPU_nonprompt_no_tErr_EB_vtx", "h_noPU_nonprompt_no_tErr_EB_vtx", 10, 0, 10);
  TH1D* h_noPU_nonprompt_tot_reco2sim_EB_vtx     = new TH1D("h_noPU_nonprompt_tot_reco2sim_EB_vtx",     "h_noPU_nonprompt_tot_reco2sim_EB_vtx",     10, 0, 10);
  TH1D* h_noPU_nonprompt_tot_EB_vtx     = new TH1D("h_noPU_nonprompt_tot_EB_vtx",     "h_noPU_nonprompt_tot_EB_vtx",     10, 0, 10);
  TH1D* h_noPU_nonprompt_PV_EE_vtx      = new TH1D("h_noPU_nonprompt_PV_EE_vtx",      "h_noPU_nonprompt_PV_EE_vtx",      10, 0, 10);
  TH1D* h_noPU_nonprompt_SV_EE_vtx      = new TH1D("h_noPU_nonprompt_SV_EE_vtx",      "h_noPU_nonprompt_SV_EE_vtx",      10, 0, 10);
  TH1D* h_noPU_nonprompt_PU_EE_vtx	= new TH1D("h_noPU_nonprompt_PU_EE_vtx",      "h_noPU_nonprompt_PU_EE_vtx",      10, 0, 10);
  TH1D* h_noPU_nonprompt_fake_EE_vtx    = new TH1D("h_noPU_nonprompt_fake_EE_vtx",    "h_noPU_nonprompt_fake_EE_vtx",    10, 0, 10);
  TH1D* h_noPU_nonprompt_no_tErr_EE_vtx = new TH1D("h_noPU_nonprompt_no_tErr_EE_vtx", "h_noPU_nonprompt_no_tErr_EE_vtx", 10, 0, 10);
  TH1D* h_noPU_nonprompt_tot_reco2sim_EE_vtx     = new TH1D("h_noPU_nonprompt_tot_reco2sim_EE_vtx",     "h_noPU_nonprompt_tot_reco2sim_EE_vtx",     10, 0, 10);
  TH1D* h_noPU_nonprompt_tot_EE_vtx     = new TH1D("h_noPU_nonprompt_tot_EE_vtx",     "h_noPU_nonprompt_tot_EE_vtx",     10, 0, 10);


  ///////////////////
  // (muon, track) //
  ///////////////////
  // PU200
    // Barrel
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_PV_EB",      "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_SV_EB",      "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_PU_EB",      "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_fake_EB",    "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_no_tErr_EB", "muon_prompt_==1 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && (muon_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_tot_reco2sim_EB",     "muon_prompt_==1 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_tot_EB",     "muon_prompt_==1 && muon_isBarrel_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_PV_EB",      "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_SV_EB",      "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_PU_EB",      "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_fake_EB",    "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_no_tErr_EB", "muon_prompt_==0 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && (muon_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_tot_reco2sim_EB",     "muon_prompt_==0 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_tot_EB",     "muon_prompt_==0 && muon_isBarrel_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
    // Endcap
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_PV_EE",      "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_SV_EE",      "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_PU_EE",      "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_fake_EE",    "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_no_tErr_EE", "muon_prompt_==1 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && (muon_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_tot_reco2sim_EE",     "muon_prompt_==1 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_muon_track_>>h_PU200_prompt_tot_EE",     "muon_prompt_==1 && muon_isBarrel_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_PV_EE",      "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_SV_EE",      "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_PU_EE",      "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_fake_EE",    "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_no_tErr_EE", "muon_prompt_==0 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && (muon_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_tot_reco2sim_EE",     "muon_prompt_==0 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_muon_track_>>h_PU200_nonprompt_tot_EE",     "muon_prompt_==0 && muon_isBarrel_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  
  // noPU
    // Barrel
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_PV_EB",      "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_SV_EB",      "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_PU_EB",      "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_fake_EB",    "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_no_tErr_EB", "muon_prompt_==1 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && (muon_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_tot_reco2sim_EB",     "muon_prompt_==1 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_tot_EB",     "muon_prompt_==1 && muon_isBarrel_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_PV_EB",      "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_SV_EB",      "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_PU_EB",      "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_fake_EB",    "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_no_tErr_EB", "muon_prompt_==0 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && (muon_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_tot_reco2sim_EB",     "muon_prompt_==0 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_tot_EB",     "muon_prompt_==0 && muon_isBarrel_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
    // Endcap
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_PV_EE",      "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_SV_EE",      "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_PU_EE",      "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_fake_EE",    "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_no_tErr_EE", "muon_prompt_==1 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && (muon_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_tot_reco2sim_EE",     "muon_prompt_==1 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_muon_track_>>h_noPU_prompt_tot_EE",     "muon_prompt_==1 && muon_isBarrel_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_PV_EE",      "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_SV_EE",      "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_PU_EE",      "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_fake_EE",    "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_no_tErr_EE", "muon_prompt_==0 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && (muon_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_tot_reco2sim_EE",     "muon_prompt_==0 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_muon_track_>>h_noPU_nonprompt_tot_EE",     "muon_prompt_==0 && muon_isBarrel_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  
  
  ///////////////////
  /// (PV, track) ///
  ///////////////////
  // PU200
    // Barrel
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_PV_EB_vtx",      "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_SV_EB_vtx",      "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_PU_EB_vtx",      "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_fake_EB_vtx",    "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_no_tErr_EB_vtx", "muon_prompt_==1 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && (vtx_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_tot_reco2sim_EB_vtx",     "muon_prompt_==1 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_tot_EB_vtx",     "muon_prompt_==1 && muon_isBarrel_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_PV_EB_vtx",      "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_SV_EB_vtx",      "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_PU_EB_vtx",      "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_fake_EB_vtx",    "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_no_tErr_EB_vtx", "muon_prompt_==0 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && (vtx_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_tot_reco2sim_EB_vtx",     "muon_prompt_==0 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_tot_EB_vtx",     "muon_prompt_==0 && muon_isBarrel_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
    // Endcap
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_PV_EE_vtx",      "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_SV_EE_vtx",      "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_PU_EE_vtx",      "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_fake_EE_vtx",    "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_no_tErr_EE_vtx", "muon_prompt_==1 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && (vtx_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_tot_reco2sim_EE_vtx",     "muon_prompt_==1 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_prompt->Draw("dtsig_vtx_track_>>h_PU200_prompt_tot_EE_vtx",     "muon_prompt_==1 && muon_isBarrel_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_PV_EE_vtx",      "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_SV_EE_vtx",      "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_PU_EE_vtx",      "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_fake_EE_vtx",    "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_no_tErr_EE_vtx", "muon_prompt_==0 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && (vtx_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_tot_reco2sim_EE_vtx",     "muon_prompt_==0 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_PU200_nonprompt->Draw("dtsig_vtx_track_>>h_PU200_nonprompt_tot_EE_vtx",     "muon_prompt_==0 && muon_isBarrel_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  
  // noPU
    // Barrel
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_PV_EB_vtx",      "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_SV_EB_vtx",      "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_PU_EB_vtx",      "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_fake_EB_vtx",    "muon_prompt_==1 && muon_isBarrel_==1 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_no_tErr_EB_vtx", "muon_prompt_==1 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && (vtx_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_tot_reco2sim_EB_vtx",     "muon_prompt_==1 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_tot_EB_vtx",     "muon_prompt_==1 && muon_isBarrel_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_PV_EB_vtx",      "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_SV_EB_vtx",      "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_PU_EB_vtx",      "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_fake_EB_vtx",    "muon_prompt_==0 && muon_isBarrel_==1 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_no_tErr_EB_vtx", "muon_prompt_==0 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && (vtx_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_tot_reco2sim_EB_vtx",     "muon_prompt_==0 && muon_isBarrel_==1 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_tot_EB_vtx",     "muon_prompt_==0 && muon_isBarrel_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
    // Endcap
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_PV_EE_vtx",      "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_SV_EE_vtx",      "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_PU_EE_vtx",      "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_fake_EE_vtx",    "muon_prompt_==1 && muon_isBarrel_==0 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_no_tErr_EE_vtx", "muon_prompt_==1 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && (vtx_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_tot_reco2sim_EE_vtx",     "muon_prompt_==1 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_prompt->Draw("dtsig_vtx_track_>>h_noPU_prompt_tot_EE_vtx",     "muon_prompt_==1 && muon_isBarrel_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_PV_EE_vtx",      "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_SV_EE_vtx",      "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==3 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_PU_EE_vtx",      "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_fake_EE_vtx",    "muon_prompt_==0 && muon_isBarrel_==0 && track_type_==2 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_no_tErr_EE_vtx", "muon_prompt_==0 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && (vtx_time_err_<=0 || track_time_err_<=0) && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_tot_reco2sim_EE_vtx",     "muon_prompt_==0 && muon_isBarrel_==0 && match_vtx_reco2sim_==1 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");
  ch_noPU_nonprompt->Draw("dtsig_vtx_track_>>h_noPU_nonprompt_tot_EE_vtx",     "muon_prompt_==0 && muon_isBarrel_==0 && selectedLV_==1 && match_vtx_sim2reco_ && match_vtx_reco2sim_==1", "goff");


  // Re-binning (tracks not having time_Err were stored in underflow)
    // muon track
  h_PU200_prompt_no_tErr_EB->SetBinContent(1,h_PU200_prompt_no_tErr_EB->GetBinContent(0)); h_PU200_prompt_no_tErr_EB->SetBinContent(0,0);
  h_PU200_prompt_no_tErr_EE->SetBinContent(1,h_PU200_prompt_no_tErr_EE->GetBinContent(0)); h_PU200_prompt_no_tErr_EE->SetBinContent(0,0);
  h_PU200_nonprompt_no_tErr_EB->SetBinContent(1,h_PU200_nonprompt_no_tErr_EB->GetBinContent(0)); h_PU200_nonprompt_no_tErr_EB->SetBinContent(0,0);
  h_PU200_nonprompt_no_tErr_EE->SetBinContent(1,h_PU200_nonprompt_no_tErr_EE->GetBinContent(0)); h_PU200_nonprompt_no_tErr_EE->SetBinContent(0,0);
  h_noPU_prompt_no_tErr_EB->SetBinContent(1,h_noPU_prompt_no_tErr_EB->GetBinContent(0)); h_noPU_prompt_no_tErr_EB->SetBinContent(0,0);
  h_noPU_prompt_no_tErr_EE->SetBinContent(1,h_noPU_prompt_no_tErr_EE->GetBinContent(0)); h_noPU_prompt_no_tErr_EE->SetBinContent(0,0);
  h_noPU_nonprompt_no_tErr_EB->SetBinContent(1,h_noPU_nonprompt_no_tErr_EB->GetBinContent(0)); h_noPU_nonprompt_no_tErr_EB->SetBinContent(0,0);
  h_noPU_nonprompt_no_tErr_EE->SetBinContent(1,h_noPU_nonprompt_no_tErr_EE->GetBinContent(0)); h_noPU_nonprompt_no_tErr_EE->SetBinContent(0,0);
    // vtx
  h_PU200_prompt_no_tErr_EB_vtx->SetBinContent(1,h_PU200_prompt_no_tErr_EB_vtx->GetBinContent(0)); h_PU200_prompt_no_tErr_EB_vtx->SetBinContent(0,0);
  h_PU200_prompt_no_tErr_EE_vtx->SetBinContent(1,h_PU200_prompt_no_tErr_EE_vtx->GetBinContent(0)); h_PU200_prompt_no_tErr_EE_vtx->SetBinContent(0,0);
  h_PU200_nonprompt_no_tErr_EB_vtx->SetBinContent(1,h_PU200_nonprompt_no_tErr_EB_vtx->GetBinContent(0)); h_PU200_nonprompt_no_tErr_EB_vtx->SetBinContent(0,0);
  h_PU200_nonprompt_no_tErr_EE_vtx->SetBinContent(1,h_PU200_nonprompt_no_tErr_EE_vtx->GetBinContent(0)); h_PU200_nonprompt_no_tErr_EE_vtx->SetBinContent(0,0);
  h_noPU_prompt_no_tErr_EB_vtx->SetBinContent(1,h_noPU_prompt_no_tErr_EB_vtx->GetBinContent(0)); h_noPU_prompt_no_tErr_EB_vtx->SetBinContent(0,0);
  h_noPU_prompt_no_tErr_EE_vtx->SetBinContent(1,h_noPU_prompt_no_tErr_EE_vtx->GetBinContent(0)); h_noPU_prompt_no_tErr_EE_vtx->SetBinContent(0,0);
  h_noPU_nonprompt_no_tErr_EB_vtx->SetBinContent(1,h_noPU_nonprompt_no_tErr_EB_vtx->GetBinContent(0)); h_noPU_nonprompt_no_tErr_EB_vtx->SetBinContent(0,0);
  h_noPU_nonprompt_no_tErr_EE_vtx->SetBinContent(1,h_noPU_nonprompt_no_tErr_EE_vtx->GetBinContent(0)); h_noPU_nonprompt_no_tErr_EE_vtx->SetBinContent(0,0);


  ///////////////
  // Cosmetics //
  ///////////////
  // (muon, track)
  h_PU200_prompt_PV_EB->SetFillColor(kGray+1); h_PU200_prompt_PU_EB->SetFillColor(kAzure+7); h_PU200_prompt_fake_EB->SetFillColor(kGreen+2); h_PU200_prompt_no_tErr_EB->SetFillColor(kYellow-7); h_PU200_prompt_SV_EB->SetFillColor(kOrange+7);
  h_PU200_prompt_PV_EE->SetFillColor(kGray+1); h_PU200_prompt_PU_EE->SetFillColor(kAzure+7); h_PU200_prompt_fake_EE->SetFillColor(kGreen+2); h_PU200_prompt_no_tErr_EE->SetFillColor(kYellow-7); h_PU200_prompt_SV_EE->SetFillColor(kOrange+7);
  h_PU200_nonprompt_PV_EB->SetFillColor(kGray+1); h_PU200_nonprompt_PU_EB->SetFillColor(kAzure+7); h_PU200_nonprompt_fake_EB->SetFillColor(kGreen+2); h_PU200_nonprompt_no_tErr_EB->SetFillColor(kYellow-7); h_PU200_nonprompt_SV_EB->SetFillColor(kOrange+7);
  h_PU200_nonprompt_PV_EE->SetFillColor(kGray+1); h_PU200_nonprompt_PU_EE->SetFillColor(kAzure+7); h_PU200_nonprompt_fake_EE->SetFillColor(kGreen+2); h_PU200_nonprompt_no_tErr_EE->SetFillColor(kYellow-7); h_PU200_nonprompt_SV_EE->SetFillColor(kOrange+7);
  h_noPU_prompt_PV_EB->SetFillColor(kGray+1); h_noPU_prompt_PU_EB->SetFillColor(kAzure+7); h_noPU_prompt_fake_EB->SetFillColor(kGreen+2); h_noPU_prompt_no_tErr_EB->SetFillColor(kYellow-7); h_noPU_prompt_SV_EB->SetFillColor(kOrange+7);
  h_noPU_prompt_PV_EE->SetFillColor(kGray+1); h_noPU_prompt_PU_EE->SetFillColor(kAzure+7); h_noPU_prompt_fake_EE->SetFillColor(kGreen+2); h_noPU_prompt_no_tErr_EE->SetFillColor(kYellow-7); h_noPU_prompt_SV_EE->SetFillColor(kOrange+7);
  h_noPU_nonprompt_PV_EB->SetFillColor(kGray+1); h_noPU_nonprompt_PU_EB->SetFillColor(kAzure+7); h_noPU_nonprompt_fake_EB->SetFillColor(kGreen+2); h_noPU_nonprompt_no_tErr_EB->SetFillColor(kYellow-7); h_noPU_nonprompt_SV_EB->SetFillColor(kOrange+7);
  h_noPU_nonprompt_PV_EE->SetFillColor(kGray+1); h_noPU_nonprompt_PU_EE->SetFillColor(kAzure+7); h_noPU_nonprompt_fake_EE->SetFillColor(kGreen+2); h_noPU_nonprompt_no_tErr_EE->SetFillColor(kYellow-7); h_noPU_nonprompt_SV_EE->SetFillColor(kOrange+7);
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
  TLegend *leg_PU200_prompt_EB = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_PV_EB, "Track from PV", "F");
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_SV_EB, "Track from SV", "F");
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_PU_EB, "PU Track", "F");
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_fake_EB, "Fake Track", "F");
  leg_PU200_prompt_EB->AddEntry(h_PU200_prompt_no_tErr_EB, "Track without tErr", "F");
  TLegend *leg_PU200_prompt_EE = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_PV_EE, "Track from PV", "F");
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_SV_EE, "Track from SV", "F");
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_PU_EE, "PU Track", "F");
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_fake_EE, "Fake Track", "F");
  leg_PU200_prompt_EE->AddEntry(h_PU200_prompt_no_tErr_EE, "Track without tErr", "F");
  TLegend *leg_PU200_nonprompt_EB = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_PV_EB, "Track from PV", "F");
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_SV_EB, "Track from SV", "F");
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_PU_EB, "PU Track", "F");
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_fake_EB, "Fake Track", "F");
  leg_PU200_nonprompt_EB->AddEntry(h_PU200_nonprompt_no_tErr_EB, "Track without tErr", "F");
  TLegend *leg_PU200_nonprompt_EE = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_PV_EE, "Track from PV", "F");
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_SV_EE, "Track from SV", "F");
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_PU_EE, "PU Track", "F");
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_fake_EE, "Fake Track", "F");
  leg_PU200_nonprompt_EE->AddEntry(h_PU200_nonprompt_no_tErr_EE, "Track without tErr", "F");
  TLegend *leg_noPU_prompt_EB = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_PV_EB, "Track from PV", "F");
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_SV_EB, "Track from SV", "F");
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_PU_EB, "PU Track", "F");
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_fake_EB, "Fake Track", "F");
  leg_noPU_prompt_EB->AddEntry(h_noPU_prompt_no_tErr_EB, "Track without tErr", "F");
  TLegend *leg_noPU_prompt_EE = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_PV_EE, "Track from PV", "F");
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_SV_EE, "Track from SV", "F");
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_PU_EE, "PU Track", "F");
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_fake_EE, "Fake Track", "F");
  leg_noPU_prompt_EE->AddEntry(h_noPU_prompt_no_tErr_EE, "Track without tErr", "F");
  TLegend *leg_noPU_nonprompt_EB = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_PV_EB, "Track from PV", "F");
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_SV_EB, "Track from SV", "F");
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_PU_EB, "PU Track", "F");
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_fake_EB, "Fake Track", "F");
  leg_noPU_nonprompt_EB->AddEntry(h_noPU_nonprompt_no_tErr_EB, "Track without tErr", "F");
  TLegend *leg_noPU_nonprompt_EE = new TLegend(0.65, 0.65, 0.88, 0.85);
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_PV_EE, "Track from PV", "F");
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_SV_EE, "Track from SV", "F");
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_PU_EE, "PU Track", "F");
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_fake_EE, "Fake Track", "F");
  leg_noPU_nonprompt_EE->AddEntry(h_noPU_nonprompt_no_tErr_EE, "Track without tErr", "F");
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

  st_PU200_prompt_EB->Add(h_PU200_prompt_PV_EB); st_PU200_prompt_EB->Add(h_PU200_prompt_PU_EB); st_PU200_prompt_EB->Add(h_PU200_prompt_fake_EB); st_PU200_prompt_EB->Add(h_PU200_prompt_SV_EB); st_PU200_prompt_EB->Add(h_PU200_prompt_no_tErr_EB); st_PU200_prompt_EB->SetTitle(";#sigma_{t};Counts");
  st_PU200_prompt_EE->Add(h_PU200_prompt_PV_EE); st_PU200_prompt_EE->Add(h_PU200_prompt_PU_EE); st_PU200_prompt_EE->Add(h_PU200_prompt_fake_EE); st_PU200_prompt_EE->Add(h_PU200_prompt_SV_EE); st_PU200_prompt_EE->Add(h_PU200_prompt_no_tErr_EE); st_PU200_prompt_EE->SetTitle(";#sigma_{t};Counts");
  st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_PV_EB); st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_PU_EB); st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_fake_EB); st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_SV_EB); st_PU200_nonprompt_EB->Add(h_PU200_nonprompt_no_tErr_EB); st_PU200_nonprompt_EB->SetTitle(";#sigma_{t};Counts");
  st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_PV_EE); st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_PU_EE); st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_fake_EE); st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_SV_EE); st_PU200_nonprompt_EE->Add(h_PU200_nonprompt_no_tErr_EE); st_PU200_nonprompt_EE->SetTitle(";#sigma_{t};Counts");
  st_noPU_prompt_EB->Add(h_noPU_prompt_PV_EB); st_noPU_prompt_EB->Add(h_noPU_prompt_PU_EB); st_noPU_prompt_EB->Add(h_noPU_prompt_fake_EB); st_noPU_prompt_EB->Add(h_noPU_prompt_SV_EB); st_noPU_prompt_EB->Add(h_noPU_prompt_no_tErr_EB); st_noPU_prompt_EB->SetTitle(";#sigma_{t};Counts");
  st_noPU_prompt_EE->Add(h_noPU_prompt_PV_EE); st_noPU_prompt_EE->Add(h_noPU_prompt_PU_EE); st_noPU_prompt_EE->Add(h_noPU_prompt_fake_EE); st_noPU_prompt_EE->Add(h_noPU_prompt_SV_EE); st_noPU_prompt_EE->Add(h_noPU_prompt_no_tErr_EE); st_noPU_prompt_EE->SetTitle(";#sigma_{t};Counts");
  st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_PV_EB); st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_PU_EB); st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_fake_EB); st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_SV_EB); st_noPU_nonprompt_EB->Add(h_noPU_nonprompt_no_tErr_EB); st_noPU_nonprompt_EB->SetTitle(";#sigma_{t};Counts");
  st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_PV_EE); st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_PU_EE); st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_fake_EE); st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_SV_EE); st_noPU_nonprompt_EE->Add(h_noPU_nonprompt_no_tErr_EE); st_noPU_nonprompt_EE->SetTitle(";#sigma_{t};Counts");
  // (PV, track)
  THStack* st_PU200_prompt_EB_vtx    = new THStack("st_PU200_prompt_EB_vtx","st_PU200_prompt_EB_vtx");
  THStack* st_PU200_prompt_EE_vtx    = new THStack("st_PU200_prompt_EE_vtx","st_PU200_prompt_EE_vtx");
  THStack* st_PU200_nonprompt_EB_vtx = new THStack("st_PU200_nonprompt_EB_vtx","st_PU200_nonprompt_EB_vtx");
  THStack* st_PU200_nonprompt_EE_vtx = new THStack("st_PU200_nonprompt_EE_vtx","st_PU200_nonprompt_EE_vtx");
  THStack* st_noPU_prompt_EB_vtx     = new THStack("st_noPU_prompt_EB_vtx","st_noPU_prompt_EB_vtx");
  THStack* st_noPU_prompt_EE_vtx     = new THStack("st_noPU_prompt_EE_vtx","st_noPU_prompt_EE_vtx");
  THStack* st_noPU_nonprompt_EB_vtx  = new THStack("st_noPU_nonprompt_EB_vtx","st_noPU_nonprompt_EB_vtx");
  THStack* st_noPU_nonprompt_EE_vtx  = new THStack("st_noPU_nonprompt_EE_vtx","st_noPU_nonprompt_EE_vtx");

  st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_PV_EB_vtx); st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_PU_EB_vtx); st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_fake_EB_vtx); st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_SV_EB_vtx); st_PU200_prompt_EB_vtx->Add(h_PU200_prompt_no_tErr_EB_vtx); st_PU200_prompt_EB_vtx->SetTitle(";#sigma_{t};Counts");
  st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_PV_EE_vtx); st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_PU_EE_vtx); st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_fake_EE_vtx); st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_SV_EE_vtx); st_PU200_prompt_EE_vtx->Add(h_PU200_prompt_no_tErr_EE_vtx); st_PU200_prompt_EE_vtx->SetTitle(";#sigma_{t};Counts");
  st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_PV_EB_vtx); st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_PU_EB_vtx); st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_fake_EB_vtx); st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_SV_EB_vtx); st_PU200_nonprompt_EB_vtx->Add(h_PU200_nonprompt_no_tErr_EB_vtx); st_PU200_nonprompt_EB_vtx->SetTitle(";#sigma_{t};Counts");
  st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_PV_EE_vtx); st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_PU_EE_vtx); st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_fake_EE_vtx); st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_SV_EE_vtx); st_PU200_nonprompt_EE_vtx->Add(h_PU200_nonprompt_no_tErr_EE_vtx); st_PU200_nonprompt_EE_vtx->SetTitle(";#sigma_{t};Counts");
  st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_PV_EB_vtx); st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_PU_EB_vtx); st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_fake_EB_vtx); st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_SV_EB_vtx); st_noPU_prompt_EB_vtx->Add(h_noPU_prompt_no_tErr_EB_vtx); st_noPU_prompt_EB_vtx->SetTitle(";#sigma_{t};Counts");
  st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_PV_EE_vtx); st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_PU_EE_vtx); st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_fake_EE_vtx); st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_SV_EE_vtx); st_noPU_prompt_EE_vtx->Add(h_noPU_prompt_no_tErr_EE_vtx); st_noPU_prompt_EE_vtx->SetTitle(";#sigma_{t};Counts");
  st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_PV_EB_vtx); st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_PU_EB_vtx); st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_fake_EB_vtx); st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_SV_EB_vtx); st_noPU_nonprompt_EB_vtx->Add(h_noPU_nonprompt_no_tErr_EB_vtx); st_noPU_nonprompt_EB_vtx->SetTitle(";#sigma_{t};Counts");
  st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_PV_EE_vtx); st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_PU_EE_vtx); st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_fake_EE_vtx); st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_SV_EE_vtx); st_noPU_nonprompt_EE_vtx->Add(h_noPU_nonprompt_no_tErr_EE_vtx); st_noPU_nonprompt_EE_vtx->SetTitle(";#sigma_{t};Counts");

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
  c_st_PU200_prompt_EB_v2->Print("plots/track_sigma_PU200_prompt_EB_v2.pdf");
  cout << "///////////////////////" << endl;
  cout << "///// muon, track /////" << endl;
  cout << "///////////////////////" << endl;
  cout << "/////    2sigma   /////" << endl;
  cout << "///////////////////////" << endl;
  cout << "PU200 Barrel prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EB->Integral(1,-1) << " : " << h_PU200_prompt_PV_EB->Integral(1,2) << " : " << h_PU200_prompt_PV_EB->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EB->Integral(3,-1)/h_PU200_prompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EB->Integral(1,-1) << " : " << h_PU200_prompt_PU_EB->Integral(1,2) << " : " << h_PU200_prompt_PU_EB->Integral(3,-1) << "  " << "(" << 100 -h_PU200_prompt_PU_EB->Integral(3,-1)/h_PU200_prompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EB->Integral(1,-1) << " : " << h_PU200_prompt_fake_EB->Integral(1,2) << " : " << h_PU200_prompt_fake_EB->Integral(3,-1) << "  " << "(" << 100 -h_PU200_prompt_fake_EB->Integral(3,-1)/h_PU200_prompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EB->Integral(1,-1) << " : " << h_PU200_prompt_SV_EB->Integral(1,2) << " : " << h_PU200_prompt_SV_EB->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EB->Integral(3,-1)/h_PU200_prompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EB->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_prompt_no_tErr_EB->Integral() << " / " << h_PU200_prompt_tot_reco2sim_EB->Integral(0,-1) << " / " << h_PU200_prompt_tot_EB->Integral(0,-1) << endl;
  cout << h_PU200_prompt_tot_EB->Integral(0,-1)-h_PU200_prompt_tot_reco2sim_EB->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_prompt_no_tErr_EB->Integral()/h_PU200_prompt_tot_reco2sim_EB->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_PU200_prompt_EE_v2 = new TCanvas("c_st_PU200_prompt_EE_v2", "c_st_PU200_prompt_EE_v2", 1500, 1500);
  c_st_PU200_prompt_EE_v2->cd();
  c_st_PU200_prompt_EE_v2->SetLeftMargin(0.12);
  st_PU200_prompt_EE->Draw();
  leg_PU200_prompt_EE->Draw();
  c_st_PU200_prompt_EE_v2->Print("plots/track_sigma_PU200_prompt_EE_v2.pdf");
  cout << "PU200 Endcap prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EE->Integral(1,-1) << " : " << h_PU200_prompt_PV_EE->Integral(1,2) << " : " << h_PU200_prompt_PV_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EE->Integral(3,-1)/h_PU200_prompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EE->Integral(1,-1) << " : " << h_PU200_prompt_PU_EE->Integral(1,2) << " : " << h_PU200_prompt_PU_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EE->Integral(3,-1)/h_PU200_prompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EE->Integral(1,-1) << " : " << h_PU200_prompt_fake_EE->Integral(1,2) << " : " << h_PU200_prompt_fake_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EE->Integral(3,-1)/h_PU200_prompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EE->Integral(1,-1) << " : " << h_PU200_prompt_SV_EE->Integral(1,2) << " : " << h_PU200_prompt_SV_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EE->Integral(3,-1)/h_PU200_prompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EE->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_prompt_no_tErr_EE->Integral() << " / " << h_PU200_prompt_tot_reco2sim_EE->Integral(0,-1) << " / " << h_PU200_prompt_tot_EE->Integral(0,-1) << endl;
  cout << h_PU200_prompt_tot_EE->Integral(0,-1)-h_PU200_prompt_tot_reco2sim_EE->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_prompt_no_tErr_EE->Integral()/h_PU200_prompt_tot_reco2sim_EE->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
    // nonprompt
      // Barrel
  TCanvas* c_st_PU200_nonprompt_EB_v2 = new TCanvas("c_st_PU200_nonprompt_EB_v2", "c_st_PU200_nonprompt_EB_v2", 1500, 1500);
  c_st_PU200_nonprompt_EB_v2->cd();
  c_st_PU200_nonprompt_EB_v2->SetLeftMargin(0.12);
  st_PU200_nonprompt_EB->Draw();
  leg_PU200_nonprompt_EB->Draw();
  c_st_PU200_nonprompt_EB_v2->Print("plots/track_sigma_PU200_nonprompt_EB_v2.pdf");
  cout << "PU200 Barrel nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EB->Integral(1,2) << " : " << h_PU200_nonprompt_PV_EB->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EB->Integral(3,-1)/h_PU200_nonprompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EB->Integral(1,2) << " : " << h_PU200_nonprompt_PU_EB->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EB->Integral(3,-1)/h_PU200_nonprompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EB->Integral(1,2) << " : " << h_PU200_nonprompt_fake_EB->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EB->Integral(3,-1)/h_PU200_nonprompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EB->Integral(1,2) << " : " << h_PU200_nonprompt_SV_EB->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EB->Integral(3,-1)/h_PU200_nonprompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EB->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_nonprompt_no_tErr_EB->Integral() << " / " << h_PU200_nonprompt_tot_reco2sim_EB->Integral(0,-1) << " / " << h_PU200_nonprompt_tot_EB->Integral(0,-1) << endl;
  cout << h_PU200_nonprompt_tot_EB->Integral(0,-1)-h_PU200_nonprompt_tot_reco2sim_EB->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_nonprompt_no_tErr_EB->Integral()/h_PU200_nonprompt_tot_reco2sim_EB->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_PU200_nonprompt_EE_v2 = new TCanvas("c_st_PU200_nonprompt_EE_v2", "c_st_PU200_nonprompt_EE_v2", 1500, 1500);
  c_st_PU200_nonprompt_EE_v2->cd();
  c_st_PU200_nonprompt_EE_v2->SetLeftMargin(0.12);
  st_PU200_nonprompt_EE->Draw();
  leg_PU200_nonprompt_EE->Draw();
  c_st_PU200_nonprompt_EE_v2->Print("plots/track_sigma_PU200_nonprompt_EE_v2.pdf");
  cout << "PU200 Endcap nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EE->Integral(1,2) << " : " << h_PU200_nonprompt_PV_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EE->Integral(3,-1)/h_PU200_nonprompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EE->Integral(1,2) << " : " << h_PU200_nonprompt_PU_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EE->Integral(3,-1)/h_PU200_nonprompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EE->Integral(1,2) << " : " << h_PU200_nonprompt_fake_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EE->Integral(3,-1)/h_PU200_nonprompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EE->Integral(1,2) << " : " << h_PU200_nonprompt_SV_EE->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EE->Integral(3,-1)/h_PU200_nonprompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EE->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_nonprompt_no_tErr_EE->Integral() << " / " << h_PU200_nonprompt_tot_reco2sim_EE->Integral(0,-1) << " / " << h_PU200_nonprompt_tot_EE->Integral(0,-1) << endl;
  cout << h_PU200_nonprompt_tot_EE->Integral(0,-1)-h_PU200_nonprompt_tot_reco2sim_EE->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_nonprompt_no_tErr_EE->Integral()/h_PU200_nonprompt_tot_reco2sim_EE->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;

  // noPU
    // prompt
      // Barrel
  TCanvas* c_st_noPU_prompt_EB_v2 = new TCanvas("c_st_noPU_prompt_EB_v2", "c_st_noPU_prompt_EB_v2", 1500, 1500);
  c_st_noPU_prompt_EB_v2->cd();
  c_st_noPU_prompt_EB_v2->SetLeftMargin(0.12);
  st_noPU_prompt_EB->Draw();
  leg_noPU_prompt_EB->Draw();
  c_st_noPU_prompt_EB_v2->Print("plots/track_sigma_noPU_prompt_EB_v2.pdf");
  cout << "noPU Barrel prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EB->Integral(1,-1) << " : " << h_noPU_prompt_PV_EB->Integral(1,2) << " : " << h_noPU_prompt_PV_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EB->Integral(3,-1)/h_noPU_prompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EB->Integral(1,-1) << " : " << h_noPU_prompt_PU_EB->Integral(1,2) << " : " << h_noPU_prompt_PU_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EB->Integral(3,-1)/h_noPU_prompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EB->Integral(1,-1) << " : " << h_noPU_prompt_fake_EB->Integral(1,2) << " : " << h_noPU_prompt_fake_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EB->Integral(3,-1)/h_noPU_prompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EB->Integral(1,-1) << " : " << h_noPU_prompt_SV_EB->Integral(1,2) << " : " << h_noPU_prompt_SV_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EB->Integral(3,-1)/h_noPU_prompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EB->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_prompt_no_tErr_EB->Integral() << " / " << h_noPU_prompt_tot_reco2sim_EB->Integral(0,-1) << " / " << h_noPU_prompt_tot_EB->Integral(0,-1) << endl;
  cout << h_noPU_prompt_tot_EB->Integral(0,-1)-h_noPU_prompt_tot_reco2sim_EB->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_prompt_no_tErr_EB->Integral()/h_noPU_prompt_tot_reco2sim_EB->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_noPU_prompt_EE_v2 = new TCanvas("c_st_noPU_prompt_EE_v2", "c_st_noPU_prompt_EE_v2", 1500, 1500);
  c_st_noPU_prompt_EE_v2->cd();
  c_st_noPU_prompt_EE_v2->SetLeftMargin(0.12);
  st_noPU_prompt_EE->Draw();
  leg_noPU_prompt_EE->Draw();
  c_st_noPU_prompt_EE_v2->Print("plots/track_sigma_noPU_prompt_EE_v2.pdf");
  cout << "noPU Endcap prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EE->Integral(1,-1) << " : " << h_noPU_prompt_PV_EE->Integral(1,2) << " : " << h_noPU_prompt_PV_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EE->Integral(3,-1)/h_noPU_prompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EE->Integral(1,-1) << " : " << h_noPU_prompt_PU_EE->Integral(1,2) << " : " << h_noPU_prompt_PU_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EE->Integral(3,-1)/h_noPU_prompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EE->Integral(1,-1) << " : " << h_noPU_prompt_fake_EE->Integral(1,2) << " : " << h_noPU_prompt_fake_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EE->Integral(3,-1)/h_noPU_prompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EE->Integral(1,-1) << " : " << h_noPU_prompt_SV_EE->Integral(1,2) << " : " << h_noPU_prompt_SV_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EE->Integral(3,-1)/h_noPU_prompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EE->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_prompt_no_tErr_EE->Integral() << " / " << h_noPU_prompt_tot_reco2sim_EE->Integral(0,-1) << " / " << h_noPU_prompt_tot_EE->Integral(0,-1) << endl;
  cout << h_noPU_prompt_tot_EE->Integral(0,-1)-h_noPU_prompt_tot_reco2sim_EE->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_prompt_no_tErr_EE->Integral()/h_noPU_prompt_tot_reco2sim_EE->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
    // nonprompt
      // Barrel
  TCanvas* c_st_noPU_nonprompt_EB_v2 = new TCanvas("c_st_noPU_nonprompt_EB_v2", "c_st_noPU_nonprompt_EB_v2", 1500, 1500);
  c_st_noPU_nonprompt_EB_v2->cd();
  c_st_noPU_nonprompt_EB_v2->SetLeftMargin(0.12);
  st_noPU_nonprompt_EB->Draw();
  leg_noPU_nonprompt_EB->Draw();
  c_st_noPU_nonprompt_EB_v2->Print("plots/track_sigma_noPU_nonprompt_EB_v2.pdf");
  cout << "noPU Barrel nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EB->Integral(1,2) << " : " << h_noPU_nonprompt_PV_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EB->Integral(3,-1)/h_noPU_nonprompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EB->Integral(1,2) << " : " << h_noPU_nonprompt_PU_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EB->Integral(3,-1)/h_noPU_nonprompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EB->Integral(1,2) << " : " << h_noPU_nonprompt_fake_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EB->Integral(3,-1)/h_noPU_nonprompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EB->Integral(1,2) << " : " << h_noPU_nonprompt_SV_EB->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EB->Integral(3,-1)/h_noPU_nonprompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EB->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_nonprompt_no_tErr_EB->Integral() << " / " << h_noPU_nonprompt_tot_reco2sim_EB->Integral(0,-1) << " / " << h_noPU_nonprompt_tot_EB->Integral(0,-1) << endl;
  cout << h_noPU_nonprompt_tot_EB->Integral(0,-1)-h_noPU_nonprompt_tot_reco2sim_EB->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_nonprompt_no_tErr_EB->Integral()/h_noPU_nonprompt_tot_reco2sim_EB->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_noPU_nonprompt_EE_v2 = new TCanvas("c_st_noPU_nonprompt_EE_v2", "c_st_noPU_nonprompt_EE_v2", 1500, 1500);
  c_st_noPU_nonprompt_EE_v2->cd();
  c_st_noPU_nonprompt_EE_v2->SetLeftMargin(0.12);
  st_noPU_nonprompt_EE->Draw();
  leg_noPU_nonprompt_EE->Draw();
  c_st_noPU_nonprompt_EE_v2->Print("plots/track_sigma_noPU_nonprompt_EE_v2.pdf");
  cout << "noPU Endcap nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EE->Integral(1,2) << " : " << h_noPU_nonprompt_PV_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EE->Integral(3,-1)/h_noPU_nonprompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EE->Integral(1,2) << " : " << h_noPU_nonprompt_PU_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EE->Integral(3,-1)/h_noPU_nonprompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EE->Integral(1,2) << " : " << h_noPU_nonprompt_fake_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EE->Integral(3,-1)/h_noPU_nonprompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EE->Integral(1,2) << " : " << h_noPU_nonprompt_SV_EE->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EE->Integral(3,-1)/h_noPU_nonprompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EE->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_nonprompt_no_tErr_EE->Integral() << " / " << h_noPU_nonprompt_tot_reco2sim_EE->Integral(0,-1) << " / " << h_noPU_nonprompt_tot_EE->Integral(0,-1) << endl;
  cout << h_noPU_nonprompt_tot_EE->Integral(0,-1)-h_noPU_nonprompt_tot_reco2sim_EE->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_nonprompt_no_tErr_EE->Integral()/h_noPU_nonprompt_tot_reco2sim_EE->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;

  cout << "///////////////////////" << endl;
  cout << "/////    3sigma   /////" << endl;
  cout << "///////////////////////" << endl;
  cout << "PU200 Barrel prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EB->Integral(1,-1) << " : " << h_PU200_prompt_PV_EB->Integral(1,3) << " : " << h_PU200_prompt_PV_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EB->Integral(4,-1)/h_PU200_prompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EB->Integral(1,-1) << " : " << h_PU200_prompt_PU_EB->Integral(1,3) << " : " << h_PU200_prompt_PU_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EB->Integral(4,-1)/h_PU200_prompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EB->Integral(1,-1) << " : " << h_PU200_prompt_fake_EB->Integral(1,3) << " : " << h_PU200_prompt_fake_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EB->Integral(4,-1)/h_PU200_prompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EB->Integral(1,-1) << " : " << h_PU200_prompt_SV_EB->Integral(1,3) << " : " << h_PU200_prompt_SV_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EB->Integral(4,-1)/h_PU200_prompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EB->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_prompt_no_tErr_EB->Integral() << " / " << h_PU200_prompt_tot_reco2sim_EB->Integral(0,-1) << " / " << h_PU200_prompt_tot_EB->Integral(0,-1) << endl;
  cout << h_PU200_prompt_tot_EB->Integral(0,-1)-h_PU200_prompt_tot_reco2sim_EB->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_prompt_no_tErr_EB->Integral()/h_PU200_prompt_tot_reco2sim_EB->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  cout << "PU200 Endcap prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EE->Integral(1,-1) << " : " << h_PU200_prompt_PV_EE->Integral(1,3) << " : " << h_PU200_prompt_PV_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EE->Integral(4,-1)/h_PU200_prompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EE->Integral(1,-1) << " : " << h_PU200_prompt_PU_EE->Integral(1,3) << " : " << h_PU200_prompt_PU_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EE->Integral(4,-1)/h_PU200_prompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EE->Integral(1,-1) << " : " << h_PU200_prompt_fake_EE->Integral(1,3) << " : " << h_PU200_prompt_fake_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EE->Integral(4,-1)/h_PU200_prompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EE->Integral(1,-1) << " : " << h_PU200_prompt_SV_EE->Integral(1,3) << " : " << h_PU200_prompt_SV_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EE->Integral(4,-1)/h_PU200_prompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EE->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_prompt_no_tErr_EE->Integral() << " / " << h_PU200_prompt_tot_reco2sim_EE->Integral(0,-1) << " / " << h_PU200_prompt_tot_EE->Integral(0,-1) << endl;
  cout << h_PU200_prompt_tot_EE->Integral(0,-1)-h_PU200_prompt_tot_reco2sim_EE->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_prompt_no_tErr_EE->Integral()/h_PU200_prompt_tot_reco2sim_EE->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
    // nonprompt
      // Barrel
  cout << "PU200 Barrel nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EB->Integral(1,3) << " : " << h_PU200_nonprompt_PV_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EB->Integral(4,-1)/h_PU200_nonprompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EB->Integral(1,3) << " : " << h_PU200_nonprompt_PU_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EB->Integral(4,-1)/h_PU200_nonprompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EB->Integral(1,3) << " : " << h_PU200_nonprompt_fake_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EB->Integral(4,-1)/h_PU200_nonprompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EB->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EB->Integral(1,3) << " : " << h_PU200_nonprompt_SV_EB->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EB->Integral(4,-1)/h_PU200_nonprompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EB->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_nonprompt_no_tErr_EB->Integral() << " / " << h_PU200_nonprompt_tot_reco2sim_EB->Integral(0,-1) << " / " << h_PU200_nonprompt_tot_EB->Integral(0,-1) << endl;
  cout << h_PU200_nonprompt_tot_EB->Integral(0,-1)-h_PU200_nonprompt_tot_reco2sim_EB->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_nonprompt_no_tErr_EB->Integral()/h_PU200_nonprompt_tot_reco2sim_EB->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  cout << "PU200 Endcap nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EE->Integral(1,3) << " : " << h_PU200_nonprompt_PV_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EE->Integral(4,-1)/h_PU200_nonprompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EE->Integral(1,3) << " : " << h_PU200_nonprompt_PU_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EE->Integral(4,-1)/h_PU200_nonprompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EE->Integral(1,3) << " : " << h_PU200_nonprompt_fake_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EE->Integral(4,-1)/h_PU200_nonprompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EE->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EE->Integral(1,3) << " : " << h_PU200_nonprompt_SV_EE->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EE->Integral(4,-1)/h_PU200_nonprompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EE->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_nonprompt_no_tErr_EE->Integral() << " / " << h_PU200_nonprompt_tot_reco2sim_EE->Integral(0,-1) << " / " << h_PU200_nonprompt_tot_EE->Integral(0,-1) << endl;
  cout << h_PU200_nonprompt_tot_EE->Integral(0,-1)-h_PU200_nonprompt_tot_reco2sim_EE->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_nonprompt_no_tErr_EE->Integral()/h_PU200_nonprompt_tot_reco2sim_EE->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;

  // noPU
    // prompt
      // Barrel
  cout << "noPU Barrel prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EB->Integral(1,-1) << " : " << h_noPU_prompt_PV_EB->Integral(1,3) << " : " << h_noPU_prompt_PV_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EB->Integral(4,-1)/h_noPU_prompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EB->Integral(1,-1) << " : " << h_noPU_prompt_PU_EB->Integral(1,3) << " : " << h_noPU_prompt_PU_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EB->Integral(4,-1)/h_noPU_prompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EB->Integral(1,-1) << " : " << h_noPU_prompt_fake_EB->Integral(1,3) << " : " << h_noPU_prompt_fake_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EB->Integral(4,-1)/h_noPU_prompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EB->Integral(1,-1) << " : " << h_noPU_prompt_SV_EB->Integral(1,3) << " : " << h_noPU_prompt_SV_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EB->Integral(4,-1)/h_noPU_prompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EB->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_prompt_no_tErr_EB->Integral() << " / " << h_noPU_prompt_tot_reco2sim_EB->Integral(0,-1) << " / " << h_noPU_prompt_tot_EB->Integral(0,-1) << endl;
  cout << h_noPU_prompt_tot_EB->Integral(0,-1)-h_noPU_prompt_tot_reco2sim_EB->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_prompt_no_tErr_EB->Integral()/h_noPU_prompt_tot_reco2sim_EB->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  cout << "noPU Endcap prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EE->Integral(1,-1) << " : " << h_noPU_prompt_PV_EE->Integral(1,3) << " : " << h_noPU_prompt_PV_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EE->Integral(4,-1)/h_noPU_prompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EE->Integral(1,-1) << " : " << h_noPU_prompt_PU_EE->Integral(1,3) << " : " << h_noPU_prompt_PU_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EE->Integral(4,-1)/h_noPU_prompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EE->Integral(1,-1) << " : " << h_noPU_prompt_fake_EE->Integral(1,3) << " : " << h_noPU_prompt_fake_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EE->Integral(4,-1)/h_noPU_prompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EE->Integral(1,-1) << " : " << h_noPU_prompt_SV_EE->Integral(1,3) << " : " << h_noPU_prompt_SV_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EE->Integral(4,-1)/h_noPU_prompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EE->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_prompt_no_tErr_EE->Integral() << " / " << h_noPU_prompt_tot_reco2sim_EE->Integral(0,-1) << " / " << h_noPU_prompt_tot_EE->Integral(0,-1) << endl;
  cout << h_noPU_prompt_tot_EE->Integral(0,-1)-h_noPU_prompt_tot_reco2sim_EE->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_prompt_no_tErr_EE->Integral()/h_noPU_prompt_tot_reco2sim_EE->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
    // nonprompt
      // Barrel
  cout << "noPU Barrel nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EB->Integral(1,3) << " : " << h_noPU_nonprompt_PV_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EB->Integral(4,-1)/h_noPU_nonprompt_PV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EB->Integral(1,3) << " : " << h_noPU_nonprompt_PU_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EB->Integral(4,-1)/h_noPU_nonprompt_PU_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EB->Integral(1,3) << " : " << h_noPU_nonprompt_fake_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EB->Integral(4,-1)/h_noPU_nonprompt_fake_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EB->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EB->Integral(1,3) << " : " << h_noPU_nonprompt_SV_EB->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EB->Integral(4,-1)/h_noPU_nonprompt_SV_EB->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EB->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_nonprompt_no_tErr_EB->Integral() << " / " << h_noPU_nonprompt_tot_reco2sim_EB->Integral(0,-1) << " / " << h_noPU_nonprompt_tot_EB->Integral(0,-1) << endl;
  cout << h_noPU_nonprompt_tot_EB->Integral(0,-1)-h_noPU_nonprompt_tot_reco2sim_EB->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_nonprompt_no_tErr_EB->Integral()/h_noPU_nonprompt_tot_reco2sim_EB->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  cout << "noPU Endcap nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EE->Integral(1,3) << " : " << h_noPU_nonprompt_PV_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EE->Integral(4,-1)/h_noPU_nonprompt_PV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EE->Integral(1,3) << " : " << h_noPU_nonprompt_PU_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EE->Integral(4,-1)/h_noPU_nonprompt_PU_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EE->Integral(1,3) << " : " << h_noPU_nonprompt_fake_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EE->Integral(4,-1)/h_noPU_nonprompt_fake_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EE->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EE->Integral(1,3) << " : " << h_noPU_nonprompt_SV_EE->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EE->Integral(4,-1)/h_noPU_nonprompt_SV_EE->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EE->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_nonprompt_no_tErr_EE->Integral() << " / " << h_noPU_nonprompt_tot_reco2sim_EE->Integral(0,-1) << " / " << h_noPU_nonprompt_tot_EE->Integral(0,-1) << endl;
  cout << h_noPU_nonprompt_tot_EE->Integral(0,-1)-h_noPU_nonprompt_tot_reco2sim_EE->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_nonprompt_no_tErr_EE->Integral()/h_noPU_nonprompt_tot_reco2sim_EE->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;


  // (PV, track)
  // PU200
    // prompt
      // Barrel
  TCanvas* c_st_PU200_prompt_EB_vtx_v2 = new TCanvas("c_st_PU200_prompt_EB_vtx_v2", "c_st_PU200_prompt_EB_vtx_v2", 1500, 1500);
  c_st_PU200_prompt_EB_vtx_v2->cd();
  c_st_PU200_prompt_EB_vtx_v2->SetLeftMargin(0.12);
  st_PU200_prompt_EB_vtx->Draw();
  leg_PU200_prompt_EB_vtx->Draw();
  c_st_PU200_prompt_EB_vtx_v2->Print("plots/track_sigma_PU200_prompt_EB_vtx_v2.pdf");
  cout << "///////////////////////" << endl;
  cout << "///// vtx, track //////" << endl;
  cout << "///////////////////////" << endl;
  cout << "/////   2sigma   //////" << endl;
  cout << "///////////////////////" << endl;
  cout << "PU200 Barrel prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PV_EB_vtx->Integral(1,2) << " : " << h_PU200_prompt_PV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EB_vtx->Integral(3,-1)/h_PU200_prompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PU_EB_vtx->Integral(1,2) << " : " << h_PU200_prompt_PU_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EB_vtx->Integral(3,-1)/h_PU200_prompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_fake_EB_vtx->Integral(1,2) << " : " << h_PU200_prompt_fake_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EB_vtx->Integral(3,-1)/h_PU200_prompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_SV_EB_vtx->Integral(1,2) << " : " << h_PU200_prompt_SV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EB_vtx->Integral(3,-1)/h_PU200_prompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EB_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_prompt_no_tErr_EB_vtx->Integral() << " / " << h_PU200_prompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " / " << h_PU200_prompt_tot_EB_vtx->Integral(0,-1) << endl;
  cout << h_PU200_prompt_tot_EB_vtx->Integral(0,-1)-h_PU200_prompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_prompt_no_tErr_EB_vtx->Integral()/h_PU200_prompt_tot_reco2sim_EB_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_PU200_prompt_EE_vtx_v2 = new TCanvas("c_st_PU200_prompt_EE_vtx_v2", "c_st_PU200_prompt_EE_vtx_v2", 1500, 1500);
  c_st_PU200_prompt_EE_vtx_v2->cd();
  c_st_PU200_prompt_EE_vtx_v2->SetLeftMargin(0.12);
  st_PU200_prompt_EE_vtx->Draw();
  leg_PU200_prompt_EE_vtx->Draw();
  c_st_PU200_prompt_EE_vtx_v2->Print("plots/track_sigma_PU200_prompt_EE_vtx_v2.pdf");
  cout << "PU200 Endcap prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PV_EE_vtx->Integral(1,2) << " : " << h_PU200_prompt_PV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EE_vtx->Integral(3,-1)/h_PU200_prompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PU_EE_vtx->Integral(1,2) << " : " << h_PU200_prompt_PU_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EE_vtx->Integral(3,-1)/h_PU200_prompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_fake_EE_vtx->Integral(1,2) << " : " << h_PU200_prompt_fake_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EE_vtx->Integral(3,-1)/h_PU200_prompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_SV_EE_vtx->Integral(1,2) << " : " << h_PU200_prompt_SV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EE_vtx->Integral(3,-1)/h_PU200_prompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EE_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_prompt_no_tErr_EE_vtx->Integral() << " / " << h_PU200_prompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " / " << h_PU200_prompt_tot_EE_vtx->Integral(0,-1) << endl;
  cout << h_PU200_prompt_tot_EE_vtx->Integral(0,-1)-h_PU200_prompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_prompt_no_tErr_EE_vtx->Integral()/h_PU200_prompt_tot_reco2sim_EE_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
    // nonprompt
      // Barrel
  TCanvas* c_st_PU200_nonprompt_EB_vtx_v2 = new TCanvas("c_st_PU200_nonprompt_EB_vtx_v2", "c_st_PU200_nonprompt_EB_vtx_v2", 1500, 1500);
  c_st_PU200_nonprompt_EB_vtx_v2->cd();
  c_st_PU200_nonprompt_EB_vtx_v2->SetLeftMargin(0.12);
  st_PU200_nonprompt_EB_vtx->Draw();
  leg_PU200_nonprompt_EB_vtx->Draw();
  c_st_PU200_nonprompt_EB_vtx_v2->Print("plots/track_sigma_PU200_nonprompt_EB_vtx_v2.pdf");
  cout << "PU200 Barrel nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EB_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_PV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EB_vtx->Integral(3,-1)/h_PU200_nonprompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EB_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_PU_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EB_vtx->Integral(3,-1)/h_PU200_nonprompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EB_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_fake_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EB_vtx->Integral(3,-1)/h_PU200_nonprompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EB_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_SV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EB_vtx->Integral(3,-1)/h_PU200_nonprompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EB_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_nonprompt_no_tErr_EB_vtx->Integral() << " / " << h_PU200_nonprompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " / " << h_PU200_nonprompt_tot_EB_vtx->Integral(0,-1) << endl;
  cout << h_PU200_nonprompt_tot_EB_vtx->Integral(0,-1)-h_PU200_nonprompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_nonprompt_no_tErr_EB_vtx->Integral()/h_PU200_nonprompt_tot_reco2sim_EB_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_PU200_nonprompt_EE_vtx_v2 = new TCanvas("c_st_PU200_nonprompt_EE_vtx_v2", "c_st_PU200_nonprompt_EE_vtx_v2", 1500, 1500);
  c_st_PU200_nonprompt_EE_vtx_v2->cd();
  c_st_PU200_nonprompt_EE_vtx_v2->SetLeftMargin(0.12);
  st_PU200_nonprompt_EE_vtx->Draw();
  leg_PU200_nonprompt_EE_vtx->Draw();
  c_st_PU200_nonprompt_EE_vtx_v2->Print("plots/track_sigma_PU200_nonprompt_EE_vtx_v2.pdf");
  cout << "PU200 Endcap nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EE_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_PV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EE_vtx->Integral(3,-1)/h_PU200_nonprompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EE_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_PU_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EE_vtx->Integral(3,-1)/h_PU200_nonprompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EE_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_fake_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EE_vtx->Integral(3,-1)/h_PU200_nonprompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EE_vtx->Integral(1,2) << " : " << h_PU200_nonprompt_SV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EE_vtx->Integral(3,-1)/h_PU200_nonprompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EE_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_nonprompt_no_tErr_EE_vtx->Integral() << " / " << h_PU200_nonprompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " / " << h_PU200_nonprompt_tot_EE_vtx->Integral(0,-1) << endl;
  cout << h_PU200_nonprompt_tot_EE_vtx->Integral(0,-1)-h_PU200_nonprompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_nonprompt_no_tErr_EE_vtx->Integral()/h_PU200_nonprompt_tot_reco2sim_EE_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;

  // noPU
    // prompt
      // Barrel
  TCanvas* c_st_noPU_prompt_EB_vtx_v2 = new TCanvas("c_st_noPU_prompt_EB_vtx_v2", "c_st_noPU_prompt_EB_vtx_v2", 1500, 1500);
  c_st_noPU_prompt_EB_vtx_v2->cd();
  c_st_noPU_prompt_EB_vtx_v2->SetLeftMargin(0.12);
  st_noPU_prompt_EB_vtx->Draw();
  leg_noPU_prompt_EB_vtx->Draw();
  c_st_noPU_prompt_EB_vtx_v2->Print("plots/track_sigma_noPU_prompt_EB_vtx_v2.pdf");
  cout << "noPU Barrel prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PV_EB_vtx->Integral(1,2) << " : " << h_noPU_prompt_PV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EB_vtx->Integral(3,-1)/h_noPU_prompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PU_EB_vtx->Integral(1,2) << " : " << h_noPU_prompt_PU_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EB_vtx->Integral(3,-1)/h_noPU_prompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_fake_EB_vtx->Integral(1,2) << " : " << h_noPU_prompt_fake_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EB_vtx->Integral(3,-1)/h_noPU_prompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_SV_EB_vtx->Integral(1,2) << " : " << h_noPU_prompt_SV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EB_vtx->Integral(3,-1)/h_noPU_prompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EB_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_prompt_no_tErr_EB_vtx->Integral() << " / " << h_noPU_prompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " / " << h_noPU_prompt_tot_EB_vtx->Integral(0,-1) << endl;
  cout << h_noPU_prompt_tot_EB_vtx->Integral(0,-1)-h_noPU_prompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_prompt_no_tErr_EB_vtx->Integral()/h_noPU_prompt_tot_reco2sim_EB_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_noPU_prompt_EE_vtx_v2 = new TCanvas("c_st_noPU_prompt_EE_vtx_v2", "c_st_noPU_prompt_EE_vtx_v2", 1500, 1500);
  c_st_noPU_prompt_EE_vtx_v2->cd();
  c_st_noPU_prompt_EE_vtx_v2->SetLeftMargin(0.12);
  st_noPU_prompt_EE_vtx->Draw();
  leg_noPU_prompt_EE_vtx->Draw();
  c_st_noPU_prompt_EE_vtx_v2->Print("plots/track_sigma_noPU_prompt_EE_vtx_v2.pdf");
  cout << "noPU Endcap prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PV_EE_vtx->Integral(1,2) << " : " << h_noPU_prompt_PV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EE_vtx->Integral(3,-1)/h_noPU_prompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PU_EE_vtx->Integral(1,2) << " : " << h_noPU_prompt_PU_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EE_vtx->Integral(3,-1)/h_noPU_prompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_fake_EE_vtx->Integral(1,2) << " : " << h_noPU_prompt_fake_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EE_vtx->Integral(3,-1)/h_noPU_prompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_SV_EE_vtx->Integral(1,2) << " : " << h_noPU_prompt_SV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EE_vtx->Integral(3,-1)/h_noPU_prompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EE_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_prompt_no_tErr_EE_vtx->Integral() << " / " << h_noPU_prompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " / " << h_noPU_prompt_tot_EE_vtx->Integral(0,-1) << endl;
  cout << h_noPU_prompt_tot_EE_vtx->Integral(0,-1)-h_noPU_prompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_prompt_no_tErr_EE_vtx->Integral()/h_noPU_prompt_tot_reco2sim_EE_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
    // nonprompt
      // Barrel
  TCanvas* c_st_noPU_nonprompt_EB_vtx_v2 = new TCanvas("c_st_noPU_nonprompt_EB_vtx_v2", "c_st_noPU_nonprompt_EB_vtx_v2", 1500, 1500);
  c_st_noPU_nonprompt_EB_vtx_v2->cd();
  c_st_noPU_nonprompt_EB_vtx_v2->SetLeftMargin(0.12);
  st_noPU_nonprompt_EB_vtx->Draw();
  leg_noPU_nonprompt_EB_vtx->Draw();
  c_st_noPU_nonprompt_EB_vtx_v2->Print("plots/track_sigma_noPU_nonprompt_EB_vtx_v2.pdf");
  cout << "noPU Barrel nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EB_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_PV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EB_vtx->Integral(3,-1)/h_noPU_nonprompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EB_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_PU_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EB_vtx->Integral(3,-1)/h_noPU_nonprompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EB_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_fake_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EB_vtx->Integral(3,-1)/h_noPU_nonprompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EB_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_SV_EB_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EB_vtx->Integral(3,-1)/h_noPU_nonprompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EB_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_nonprompt_no_tErr_EB_vtx->Integral() << " / " << h_noPU_nonprompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " / " << h_noPU_nonprompt_tot_EB_vtx->Integral(0,-1) << endl;
  cout << h_noPU_nonprompt_tot_EB_vtx->Integral(0,-1)-h_noPU_nonprompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_nonprompt_no_tErr_EB_vtx->Integral()/h_noPU_nonprompt_tot_reco2sim_EB_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  TCanvas* c_st_noPU_nonprompt_EE_vtx_v2 = new TCanvas("c_st_noPU_nonprompt_EE_vtx_v2", "c_st_noPU_nonprompt_EE_vtx_v2", 1500, 1500);
  c_st_noPU_nonprompt_EE_vtx_v2->cd();
  c_st_noPU_nonprompt_EE_vtx_v2->SetLeftMargin(0.12);
  st_noPU_nonprompt_EE_vtx->Draw();
  leg_noPU_nonprompt_EE_vtx->Draw();
  c_st_noPU_nonprompt_EE_vtx_v2->Print("plots/track_sigma_noPU_nonprompt_EE_vtx_v2.pdf");
  cout << "noPU Endcap nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EE_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_PV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EE_vtx->Integral(3,-1)/h_noPU_nonprompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EE_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_PU_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EE_vtx->Integral(3,-1)/h_noPU_nonprompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EE_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_fake_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EE_vtx->Integral(3,-1)/h_noPU_nonprompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EE_vtx->Integral(1,2) << " : " << h_noPU_nonprompt_SV_EE_vtx->Integral(3,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EE_vtx->Integral(3,-1)/h_noPU_nonprompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EE_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_nonprompt_no_tErr_EE_vtx->Integral() << " / " << h_noPU_nonprompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " / " << h_noPU_nonprompt_tot_EE_vtx->Integral(0,-1) << endl;
  cout << h_noPU_nonprompt_tot_EE_vtx->Integral(0,-1)-h_noPU_nonprompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_nonprompt_no_tErr_EE_vtx->Integral()/h_noPU_nonprompt_tot_reco2sim_EE_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;

  cout << "///////////////////////" << endl;
  cout << "/////   3sigma   //////" << endl;
  cout << "///////////////////////" << endl;
  cout << "PU200 Barrel prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PV_EB_vtx->Integral(1,3) << " : " << h_PU200_prompt_PV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EB_vtx->Integral(4,-1)/h_PU200_prompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PU_EB_vtx->Integral(1,3) << " : " << h_PU200_prompt_PU_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EB_vtx->Integral(4,-1)/h_PU200_prompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_fake_EB_vtx->Integral(1,3) << " : " << h_PU200_prompt_fake_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EB_vtx->Integral(4,-1)/h_PU200_prompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EB_vtx->Integral(1,-1) << " : " << h_PU200_prompt_SV_EB_vtx->Integral(1,3) << " : " << h_PU200_prompt_SV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EB_vtx->Integral(4,-1)/h_PU200_prompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EB_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_prompt_no_tErr_EB_vtx->Integral() << " / " << h_PU200_prompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " / " << h_PU200_prompt_tot_EB_vtx->Integral(0,-1) << endl;
  cout << h_PU200_prompt_tot_EB_vtx->Integral(0,-1)-h_PU200_prompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_prompt_no_tErr_EB_vtx->Integral()/h_PU200_prompt_tot_reco2sim_EB_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  cout << "PU200 Endcap prompt" << endl;
  cout << "[PV]      " << h_PU200_prompt_PV_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PV_EE_vtx->Integral(1,3) << " : " << h_PU200_prompt_PV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PV_EE_vtx->Integral(4,-1)/h_PU200_prompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_prompt_PU_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_PU_EE_vtx->Integral(1,3) << " : " << h_PU200_prompt_PU_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_PU_EE_vtx->Integral(4,-1)/h_PU200_prompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_prompt_fake_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_fake_EE_vtx->Integral(1,3) << " : " << h_PU200_prompt_fake_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_fake_EE_vtx->Integral(4,-1)/h_PU200_prompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_prompt_SV_EE_vtx->Integral(1,-1) << " : " << h_PU200_prompt_SV_EE_vtx->Integral(1,3) << " : " << h_PU200_prompt_SV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_prompt_SV_EE_vtx->Integral(4,-1)/h_PU200_prompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_prompt_no_tErr_EE_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_prompt_no_tErr_EE_vtx->Integral() << " / " << h_PU200_prompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " / " << h_PU200_prompt_tot_EE_vtx->Integral(0,-1) << endl;
  cout << h_PU200_prompt_tot_EE_vtx->Integral(0,-1)-h_PU200_prompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_prompt_no_tErr_EE_vtx->Integral()/h_PU200_prompt_tot_reco2sim_EE_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
    // nonprompt
      // Barrel
  cout << "PU200 Barrel nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EB_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_PV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EB_vtx->Integral(4,-1)/h_PU200_nonprompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EB_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_PU_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EB_vtx->Integral(4,-1)/h_PU200_nonprompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EB_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_fake_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EB_vtx->Integral(4,-1)/h_PU200_nonprompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EB_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EB_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_SV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EB_vtx->Integral(4,-1)/h_PU200_nonprompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EB_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_nonprompt_no_tErr_EB_vtx->Integral() << " / " << h_PU200_nonprompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " / " << h_PU200_nonprompt_tot_EB_vtx->Integral(0,-1) << endl;
  cout << h_PU200_nonprompt_tot_EB_vtx->Integral(0,-1)-h_PU200_nonprompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_nonprompt_no_tErr_EB_vtx->Integral()/h_PU200_nonprompt_tot_reco2sim_EB_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  cout << "PU200 Endcap nonprompt" << endl;
  cout << "[PV]      " << h_PU200_nonprompt_PV_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PV_EE_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_PV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PV_EE_vtx->Integral(4,-1)/h_PU200_nonprompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_PU200_nonprompt_PU_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_PU_EE_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_PU_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_PU_EE_vtx->Integral(4,-1)/h_PU200_nonprompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_PU200_nonprompt_fake_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_fake_EE_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_fake_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_fake_EE_vtx->Integral(4,-1)/h_PU200_nonprompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_PU200_nonprompt_SV_EE_vtx->Integral(1,-1) << " : " << h_PU200_nonprompt_SV_EE_vtx->Integral(1,3) << " : " << h_PU200_nonprompt_SV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_PU200_nonprompt_SV_EE_vtx->Integral(4,-1)/h_PU200_nonprompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_PU200_nonprompt_no_tErr_EE_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_PU200_nonprompt_no_tErr_EE_vtx->Integral() << " / " << h_PU200_nonprompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " / " << h_PU200_nonprompt_tot_EE_vtx->Integral(0,-1) << endl;
  cout << h_PU200_nonprompt_tot_EE_vtx->Integral(0,-1)-h_PU200_nonprompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_PU200_nonprompt_no_tErr_EE_vtx->Integral()/h_PU200_nonprompt_tot_reco2sim_EE_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;

  // noPU
    // prompt
      // Barrel
  cout << "noPU Barrel prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PV_EB_vtx->Integral(1,3) << " : " << h_noPU_prompt_PV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EB_vtx->Integral(4,-1)/h_noPU_prompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PU_EB_vtx->Integral(1,3) << " : " << h_noPU_prompt_PU_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EB_vtx->Integral(4,-1)/h_noPU_prompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_fake_EB_vtx->Integral(1,3) << " : " << h_noPU_prompt_fake_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EB_vtx->Integral(4,-1)/h_noPU_prompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EB_vtx->Integral(1,-1) << " : " << h_noPU_prompt_SV_EB_vtx->Integral(1,3) << " : " << h_noPU_prompt_SV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EB_vtx->Integral(4,-1)/h_noPU_prompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EB_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_prompt_no_tErr_EB_vtx->Integral() << " / " << h_noPU_prompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " / " << h_noPU_prompt_tot_EB_vtx->Integral(0,-1) << endl;
  cout << h_noPU_prompt_tot_EB_vtx->Integral(0,-1)-h_noPU_prompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_prompt_no_tErr_EB_vtx->Integral()/h_noPU_prompt_tot_reco2sim_EB_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  cout << "noPU Endcap prompt" << endl;
  cout << "[PV]      " << h_noPU_prompt_PV_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PV_EE_vtx->Integral(1,3) << " : " << h_noPU_prompt_PV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PV_EE_vtx->Integral(4,-1)/h_noPU_prompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_prompt_PU_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_PU_EE_vtx->Integral(1,3) << " : " << h_noPU_prompt_PU_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_PU_EE_vtx->Integral(4,-1)/h_noPU_prompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_prompt_fake_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_fake_EE_vtx->Integral(1,3) << " : " << h_noPU_prompt_fake_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_fake_EE_vtx->Integral(4,-1)/h_noPU_prompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_prompt_SV_EE_vtx->Integral(1,-1) << " : " << h_noPU_prompt_SV_EE_vtx->Integral(1,3) << " : " << h_noPU_prompt_SV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_prompt_SV_EE_vtx->Integral(4,-1)/h_noPU_prompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_prompt_no_tErr_EE_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_prompt_no_tErr_EE_vtx->Integral() << " / " << h_noPU_prompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " / " << h_noPU_prompt_tot_EE_vtx->Integral(0,-1) << endl;
  cout << h_noPU_prompt_tot_EE_vtx->Integral(0,-1)-h_noPU_prompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_prompt_no_tErr_EE_vtx->Integral()/h_noPU_prompt_tot_reco2sim_EE_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
    // nonprompt
      // Barrel
  cout << "noPU Barrel nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EB_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_PV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EB_vtx->Integral(4,-1)/h_noPU_nonprompt_PV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EB_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_PU_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EB_vtx->Integral(4,-1)/h_noPU_nonprompt_PU_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EB_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_fake_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EB_vtx->Integral(4,-1)/h_noPU_nonprompt_fake_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EB_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EB_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_SV_EB_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EB_vtx->Integral(4,-1)/h_noPU_nonprompt_SV_EB_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EB_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_nonprompt_no_tErr_EB_vtx->Integral() << " / " << h_noPU_nonprompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " / " << h_noPU_nonprompt_tot_EB_vtx->Integral(0,-1) << endl;
  cout << h_noPU_nonprompt_tot_EB_vtx->Integral(0,-1)-h_noPU_nonprompt_tot_reco2sim_EB_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_nonprompt_no_tErr_EB_vtx->Integral()/h_noPU_nonprompt_tot_reco2sim_EB_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;
      // Endcap
  cout << "noPU Endcap nonprompt" << endl;
  cout << "[PV]      " << h_noPU_nonprompt_PV_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PV_EE_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_PV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PV_EE_vtx->Integral(4,-1)/h_noPU_nonprompt_PV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[PU]      " << h_noPU_nonprompt_PU_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_PU_EE_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_PU_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_PU_EE_vtx->Integral(4,-1)/h_noPU_nonprompt_PU_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[fake]    " << h_noPU_nonprompt_fake_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_fake_EE_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_fake_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_fake_EE_vtx->Integral(4,-1)/h_noPU_nonprompt_fake_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[SV]      " << h_noPU_nonprompt_SV_EE_vtx->Integral(1,-1) << " : " << h_noPU_nonprompt_SV_EE_vtx->Integral(1,3) << " : " << h_noPU_nonprompt_SV_EE_vtx->Integral(4,-1) << "  " << "(" << 100 - h_noPU_nonprompt_SV_EE_vtx->Integral(4,-1)/h_noPU_nonprompt_SV_EE_vtx->Integral(1,-1)*100 << " % is considered)" << endl;
  cout << "[no tErr] " << h_noPU_nonprompt_no_tErr_EE_vtx->Integral() << endl;
  cout << "[no tErr / reco2sim / total] : " << h_noPU_nonprompt_no_tErr_EE_vtx->Integral() << " / " << h_noPU_nonprompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " / " << h_noPU_nonprompt_tot_EE_vtx->Integral(0,-1) << endl;
  cout << h_noPU_nonprompt_tot_EE_vtx->Integral(0,-1)-h_noPU_nonprompt_tot_reco2sim_EE_vtx->Integral(0,-1) << " tracks are not reco2sim matched,  " << h_noPU_nonprompt_no_tErr_EE_vtx->Integral()/h_noPU_nonprompt_tot_reco2sim_EE_vtx->Integral(0,-1)*100 << " tracks are not having time_err" << endl;
  cout << endl;


}

void muon_mother_pdgId() {
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


  TH1D* h_PU200_prompt    = new TH1D("h_PU200_prompt",    "h_PU200_prompt",    45, -15, 30);
  //TH1D* h_PU200_nonprompt = new TH1D("h_PU200_nonprompt", "h_PU200_nonprompt", 12000, -6000, 6000);
  TH1D* h_PU200_nonprompt = new TH1D("h_PU200_nonprompt", "h_PU200_nonprompt", 1000, -500, 500);
  TH1D* h_noPU_prompt     = new TH1D("h_noPU_prompt",     "h_noPU_prompt",     45, -15, 30);
  //TH1D* h_noPU_nonprompt  = new TH1D("h_noPU_nonprompt",  "h_noPU_nonprompt",  12000, -6000, 6000);
  TH1D* h_noPU_nonprompt  = new TH1D("h_noPU_nonprompt",  "h_noPU_nonprompt",  1000, -500, 500);
  // PU200
    // Barrel
  ch_PU200_prompt->Draw("muon_mother_pdgId_>>h_PU200_prompt",       "muon_prompt_==1", "goff");
  ch_PU200_nonprompt->Draw("muon_mother_pdgId_>>h_PU200_nonprompt", "muon_prompt_==0", "goff");
  ch_noPU_prompt->Draw("muon_mother_pdgId_>>h_noPU_prompt",         "muon_prompt_==1", "goff");
  ch_noPU_nonprompt->Draw("muon_mother_pdgId_>>h_noPU_nonprompt",   "muon_prompt_==0", "goff");

  cout << endl;
  cout << endl;
  cout << "//////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "////////////////////////////  Check pdgId of muons  //////////////////////////" << endl;
  cout << "//////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "[PU200 prompt]" << endl;
  for(int i=1; i<h_PU200_prompt->GetNbinsX(); i++) {
    if(h_PU200_prompt->GetBinContent(i)!=0) {
	  cout << pdgId(h_PU200_prompt->GetXaxis()->GetBinLowEdge(i)) << ", " << h_PU200_prompt->GetBinContent(i) << endl;
	}
  }
  cout << endl;
  cout << "[PU200 nonprompt]" << endl;
  for(int i=1; i<h_PU200_nonprompt->GetNbinsX(); i++) {
    if(h_PU200_nonprompt->GetBinContent(i)!=0) {
	  cout << pdgId(h_PU200_nonprompt->GetXaxis()->GetBinLowEdge(i)) << ", " << h_PU200_nonprompt->GetBinContent(i) << endl;
	}
  }
  cout << endl;
  cout << "[noPU prompt]" << endl;
  for(int i=1; i<h_noPU_prompt->GetNbinsX(); i++) {
    if(h_noPU_prompt->GetBinContent(i)!=0) {
	  cout << pdgId(h_noPU_prompt->GetXaxis()->GetBinLowEdge(i)) << ", " << h_noPU_prompt->GetBinContent(i) << endl;
	}
  }
  cout << endl;
  cout << "[noPU nonprompt]" << endl;
  for(int i=1; i<h_noPU_nonprompt->GetNbinsX(); i++) {
    if(h_noPU_nonprompt->GetBinContent(i)!=0) {
	  cout << pdgId(h_noPU_nonprompt->GetXaxis()->GetBinLowEdge(i)) << ", " << h_noPU_nonprompt->GetBinContent(i) << endl;
	}
  }

  cout << "//////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "////////////////////////////  Get under/over flow  ///////////////////////////" << endl;
  cout << "//////////////////////////////////////////////////////////////////////////////" << endl;
  // Get underflow and overflow
  cout << "PU200_prompt: " << h_PU200_prompt->GetBinContent(0) + h_PU200_prompt->GetBinContent(31) << endl;
  cout << "PU200_nonprompt: " << h_PU200_nonprompt->GetBinContent(0) + h_PU200_nonprompt->GetBinContent(601) << endl;
  cout << "noPU_prompt: " << h_noPU_prompt->GetBinContent(0) + h_noPU_prompt->GetBinContent(31) << endl;
  cout << "noPU_nonprompt: " << h_noPU_nonprompt->GetBinContent(0) + h_noPU_nonprompt->GetBinContent(601) << endl;

  // Normalization
  h_PU200_prompt->Scale(1./h_PU200_prompt->Integral(0,-1));
  h_PU200_nonprompt->Scale(1./h_PU200_nonprompt->Integral(0,-1));
  h_noPU_prompt->Scale(1./h_noPU_prompt->Integral(0,-1));
  h_noPU_nonprompt->Scale(1./h_noPU_nonprompt->Integral(0,-1));

  // Cosmetics
  h_PU200_prompt->SetLineColor(kBlack); h_PU200_prompt->GetXaxis()->SetTitle("pdgId of mother particle of muon"); h_PU200_prompt->GetYaxis()->SetTitle("% Counts"); h_PU200_prompt->SetTitle("");
  h_PU200_nonprompt->SetLineColor(kBlack); h_PU200_nonprompt->GetXaxis()->SetTitle("pdgId of mother particle of muon"); h_PU200_nonprompt->GetYaxis()->SetTitle("% Counts"); h_PU200_nonprompt->SetTitle("");
  h_noPU_prompt->SetLineColor(kBlack); h_noPU_prompt->GetXaxis()->SetTitle("pdgId of mother particle of muon"); h_noPU_prompt->GetYaxis()->SetTitle("% Counts"); h_noPU_prompt->SetTitle("");
  h_noPU_nonprompt->SetLineColor(kBlack); h_noPU_nonprompt->GetXaxis()->SetTitle("pdgId of mother particle of muon"); h_noPU_nonprompt->GetYaxis()->SetTitle("% Counts"); h_noPU_nonprompt->SetTitle("");


  TCanvas* c_PU200_prompt_pdgId = new TCanvas("c_PU200_prompt_pdgId", "c_PU200_prompt_pdgId", 1500, 1500);
  c_PU200_prompt_pdgId->cd();
  c_PU200_prompt_pdgId->SetLeftMargin(0.12);
  h_PU200_prompt->Draw("hist");
  c_PU200_prompt_pdgId->Print("plots/PU200_prompt_muon_mother_pdgId.pdf");

  TCanvas* c_PU200_nonprompt_pdgId = new TCanvas("c_PU200_nonprompt_pdgId", "c_PU200_nonprompt_pdgId", 1500, 1500);
  c_PU200_nonprompt_pdgId->cd();
  c_PU200_nonprompt_pdgId->SetLeftMargin(0.12);
  h_PU200_nonprompt->Draw("hist");
  c_PU200_nonprompt_pdgId->Print("plots/PU200_nonprompt_muon_mother_pdgId.pdf");

  TCanvas* c_noPU_prompt_pdgId = new TCanvas("c_noPU_prompt_pdgId", "c_noPU_prompt_pdgId", 1500, 1500);
  c_noPU_prompt_pdgId->cd();
  c_noPU_prompt_pdgId->SetLeftMargin(0.12);
  h_noPU_prompt->Draw("hist");
  c_noPU_prompt_pdgId->Print("plots/noPU_prompt_muon_mother_pdgId.pdf");

  TCanvas* c_noPU_nonprompt_pdgId = new TCanvas("c_noPU_nonprompt_pdgId", "c_noPU_nonprompt_pdgId", 1500, 1500);
  c_noPU_nonprompt_pdgId->cd();
  c_noPU_nonprompt_pdgId->SetLeftMargin(0.12);
  h_noPU_nonprompt->Draw("hist");
  c_noPU_nonprompt_pdgId->Print("plots/noPU_nonprompt_muon_mother_pdgId.pdf");

}

void fraction_mva_cut() {
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


  TH1D* h_PU200_prompt_muon             = new TH1D("h_PU200_prompt_muon",            "h_PU200_prompt_muon",            1,0,1);
  TH1D* h_PU200_prompt_muon_mva         = new TH1D("h_PU200_prompt_muon_mva",        "h_PU200_prompt_muon_mva",        1,0,1);
  TH1D* h_PU200_prompt_muon_time        = new TH1D("h_PU200_prompt_muon_time",       "h_PU200_prompt_muon_time",       1,0,1);
  TH1D* h_PU200_prompt_muon_time_err    = new TH1D("h_PU200_prompt_muon_time_err",   "h_PU200_prompt_muon_time_err",   1,0,1);
  TH1D* h_PU200_prompt_muon_mva_cut     = new TH1D("h_PU200_prompt_muon_mva_cut",    "h_PU200_prompt_muon_mva_cut",    1,0,1);
  TH1D* h_PU200_prompt_vtx             = new TH1D("h_PU200_prompt_vtx",            "h_PU200_prompt_vtx",            1,0,1);
  TH1D* h_PU200_prompt_vtx_time        = new TH1D("h_PU200_prompt_vtx_time",       "h_PU200_prompt_vtx_time",       1,0,1);
  TH1D* h_PU200_prompt_vtx_time_err    = new TH1D("h_PU200_prompt_vtx_time_err",   "h_PU200_prompt_vtx_time_err",   1,0,1);
  TH1D* h_PU200_prompt_tot_vtx_track      = new TH1D("h_PU200_prompt_tot_vtx_track",     "h_PU200_prompt_tot_vtx_track",     1,0,1);
  TH1D* h_PU200_prompt_track            = new TH1D("h_PU200_prompt_track",           "h_PU200_prompt_track",           1,0,1);
  TH1D* h_PU200_prompt_track_mva        = new TH1D("h_PU200_prompt_track_mva",       "h_PU200_prompt_track_mva",       1,0,1);
  TH1D* h_PU200_prompt_track_time       = new TH1D("h_PU200_prompt_track_time",      "h_PU200_prompt_track_time",      1,0,1);
  TH1D* h_PU200_prompt_track_time_err   = new TH1D("h_PU200_prompt_track_time_err",  "h_PU200_prompt_track_time_err",  1,0,1);
  TH1D* h_PU200_prompt_track_mva_cut    = new TH1D("h_PU200_prompt_track_mva_cut",   "h_PU200_prompt_track_mva_cut",   1,0,1);
  TH1D* h_PU200_prompt_tot_mva_cut      = new TH1D("h_PU200_prompt_tot_mva_cut",     "h_PU200_prompt_tot_mva_cut",     1,0,1);

  TH1D* h_PU200_nonprompt_muon             = new TH1D("h_PU200_nonprompt_muon",            "h_PU200_nonprompt_muon",            1,0,1);
  TH1D* h_PU200_nonprompt_muon_mva         = new TH1D("h_PU200_nonprompt_muon_mva",        "h_PU200_nonprompt_muon_mva",        1,0,1);
  TH1D* h_PU200_nonprompt_muon_time        = new TH1D("h_PU200_nonprompt_muon_time",       "h_PU200_nonprompt_muon_time",       1,0,1);
  TH1D* h_PU200_nonprompt_muon_time_err    = new TH1D("h_PU200_nonprompt_muon_time_err",   "h_PU200_nonprompt_muon_time_err",   1,0,1);
  TH1D* h_PU200_nonprompt_muon_mva_cut     = new TH1D("h_PU200_nonprompt_muon_mva_cut",    "h_PU200_nonprompt_muon_mva_cut",    1,0,1);
  TH1D* h_PU200_nonprompt_vtx             = new TH1D("h_PU200_nonprompt_vtx",            "h_PU200_nonprompt_vtx",            1,0,1);
  TH1D* h_PU200_nonprompt_vtx_time        = new TH1D("h_PU200_nonprompt_vtx_time",       "h_PU200_nonprompt_vtx_time",       1,0,1);
  TH1D* h_PU200_nonprompt_vtx_time_err    = new TH1D("h_PU200_nonprompt_vtx_time_err",   "h_PU200_nonprompt_vtx_time_err",   1,0,1);
  TH1D* h_PU200_nonprompt_tot_vtx_track      = new TH1D("h_PU200_nonprompt_tot_vtx_track",     "h_PU200_nonprompt_tot_vtx_track",     1,0,1);
  TH1D* h_PU200_nonprompt_track            = new TH1D("h_PU200_nonprompt_track",           "h_PU200_nonprompt_track",           1,0,1);
  TH1D* h_PU200_nonprompt_track_mva        = new TH1D("h_PU200_nonprompt_track_mva",       "h_PU200_nonprompt_track_mva",       1,0,1);
  TH1D* h_PU200_nonprompt_track_time       = new TH1D("h_PU200_nonprompt_track_time",      "h_PU200_nonprompt_track_time",      1,0,1);
  TH1D* h_PU200_nonprompt_track_time_err   = new TH1D("h_PU200_nonprompt_track_time_err",  "h_PU200_nonprompt_track_time_err",  1,0,1);
  TH1D* h_PU200_nonprompt_track_mva_cut    = new TH1D("h_PU200_nonprompt_track_mva_cut",   "h_PU200_nonprompt_track_mva_cut",   1,0,1);
  TH1D* h_PU200_nonprompt_tot_mva_cut      = new TH1D("h_PU200_nonprompt_tot_mva_cut",     "h_PU200_nonprompt_tot_mva_cut",     1,0,1);

  TH1D* h_noPU_prompt_muon             = new TH1D("h_noPU_prompt_muon",            "h_noPU_prompt_muon",            1,0,1);
  TH1D* h_noPU_prompt_muon_mva         = new TH1D("h_noPU_prompt_muon_mva",        "h_noPU_prompt_muon_mva",        1,0,1);
  TH1D* h_noPU_prompt_muon_time        = new TH1D("h_noPU_prompt_muon_time",       "h_noPU_prompt_muon_time",       1,0,1);
  TH1D* h_noPU_prompt_muon_time_err    = new TH1D("h_noPU_prompt_muon_time_err",   "h_noPU_prompt_muon_time_err",   1,0,1);
  TH1D* h_noPU_prompt_muon_mva_cut     = new TH1D("h_noPU_prompt_muon_mva_cut",    "h_noPU_prompt_muon_mva_cut",    1,0,1);
  TH1D* h_noPU_prompt_vtx             = new TH1D("h_noPU_prompt_vtx",            "h_noPU_prompt_vtx",            1,0,1);
  TH1D* h_noPU_prompt_vtx_time        = new TH1D("h_noPU_prompt_vtx_time",       "h_noPU_prompt_vtx_time",       1,0,1);
  TH1D* h_noPU_prompt_vtx_time_err    = new TH1D("h_noPU_prompt_vtx_time_err",   "h_noPU_prompt_vtx_time_err",   1,0,1);
  TH1D* h_noPU_prompt_tot_vtx_track      = new TH1D("h_noPU_prompt_tot_vtx_track",     "h_noPU_prompt_tot_vtx_track",     1,0,1);
  TH1D* h_noPU_prompt_track            = new TH1D("h_noPU_prompt_track",           "h_noPU_prompt_track",           1,0,1);
  TH1D* h_noPU_prompt_track_mva        = new TH1D("h_noPU_prompt_track_mva",       "h_noPU_prompt_track_mva",       1,0,1);
  TH1D* h_noPU_prompt_track_time       = new TH1D("h_noPU_prompt_track_time",      "h_noPU_prompt_track_time",      1,0,1);
  TH1D* h_noPU_prompt_track_time_err   = new TH1D("h_noPU_prompt_track_time_err",  "h_noPU_prompt_track_time_err",  1,0,1);
  TH1D* h_noPU_prompt_track_mva_cut    = new TH1D("h_noPU_prompt_track_mva_cut",   "h_noPU_prompt_track_mva_cut",   1,0,1);
  TH1D* h_noPU_prompt_tot_mva_cut      = new TH1D("h_noPU_prompt_tot_mva_cut",     "h_noPU_prompt_tot_mva_cut",     1,0,1);

  TH1D* h_noPU_nonprompt_muon             = new TH1D("h_noPU_nonprompt_muon",            "h_noPU_nonprompt_muon",            1,0,1);
  TH1D* h_noPU_nonprompt_muon_mva         = new TH1D("h_noPU_nonprompt_muon_mva",        "h_noPU_nonprompt_muon_mva",        1,0,1);
  TH1D* h_noPU_nonprompt_muon_time        = new TH1D("h_noPU_nonprompt_muon_time",       "h_noPU_nonprompt_muon_time",       1,0,1);
  TH1D* h_noPU_nonprompt_muon_time_err    = new TH1D("h_noPU_nonprompt_muon_time_err",   "h_noPU_nonprompt_muon_time_err",   1,0,1);
  TH1D* h_noPU_nonprompt_muon_mva_cut     = new TH1D("h_noPU_nonprompt_muon_mva_cut",    "h_noPU_nonprompt_muon_mva_cut",    1,0,1);
  TH1D* h_noPU_nonprompt_vtx             = new TH1D("h_noPU_nonprompt_vtx",            "h_noPU_nonprompt_vtx",            1,0,1);
  TH1D* h_noPU_nonprompt_vtx_time        = new TH1D("h_noPU_nonprompt_vtx_time",       "h_noPU_nonprompt_vtx_time",       1,0,1);
  TH1D* h_noPU_nonprompt_vtx_time_err    = new TH1D("h_noPU_nonprompt_vtx_time_err",   "h_noPU_nonprompt_vtx_time_err",   1,0,1);
  TH1D* h_noPU_nonprompt_tot_vtx_track      = new TH1D("h_noPU_nonprompt_tot_vtx_track",     "h_noPU_nonprompt_tot_vtx_track",     1,0,1);
  TH1D* h_noPU_nonprompt_track            = new TH1D("h_noPU_nonprompt_track",           "h_noPU_nonprompt_track",           1,0,1);
  TH1D* h_noPU_nonprompt_track_mva        = new TH1D("h_noPU_nonprompt_track_mva",       "h_noPU_nonprompt_track_mva",       1,0,1);
  TH1D* h_noPU_nonprompt_track_time       = new TH1D("h_noPU_nonprompt_track_time",      "h_noPU_nonprompt_track_time",      1,0,1);
  TH1D* h_noPU_nonprompt_track_time_err   = new TH1D("h_noPU_nonprompt_track_time_err",  "h_noPU_nonprompt_track_time_err",  1,0,1);
  TH1D* h_noPU_nonprompt_track_mva_cut    = new TH1D("h_noPU_nonprompt_track_mva_cut",   "h_noPU_nonprompt_track_mva_cut",   1,0,1);
  TH1D* h_noPU_nonprompt_tot_mva_cut      = new TH1D("h_noPU_nonprompt_tot_mva_cut",     "h_noPU_nonprompt_tot_mva_cut",     1,0,1);


  // PU200 prompt
    // muon
  ch_PU200_prompt->Draw("muon_pt_>>h_PU200_prompt_muon",             "muon_prompt_==1", 				      "goff");
  ch_PU200_prompt->Draw("muon_pt_>>h_PU200_prompt_muon_mva",         "muon_prompt_==1 && muon_mva_>0",        "goff");
  ch_PU200_prompt->Draw("muon_pt_>>h_PU200_prompt_muon_time",        "muon_prompt_==1 && muon_time_!=0",      "goff");
  ch_PU200_prompt->Draw("muon_pt_>>h_PU200_prompt_muon_time_err",    "muon_prompt_==1 && muon_time_err_i_>0", "goff");
  ch_PU200_prompt->Draw("muon_pt_>>h_PU200_prompt_muon_mva_cut",     "muon_prompt_==1 && muon_mva_>0.5",      "goff");
    // track
  ch_PU200_prompt->Draw("track_pt_>>h_PU200_prompt_track",             "muon_prompt_==1", 				         "goff");
  ch_PU200_prompt->Draw("track_pt_>>h_PU200_prompt_track_mva",         "muon_prompt_==1 && track_mva_>0",        "goff");
  ch_PU200_prompt->Draw("track_pt_>>h_PU200_prompt_track_time",        "muon_prompt_==1 && track_time_!=0",       "goff");
  ch_PU200_prompt->Draw("track_pt_>>h_PU200_prompt_track_time_err",    "muon_prompt_==1 && track_time_err_i_>0", "goff");
  ch_PU200_prompt->Draw("track_pt_>>h_PU200_prompt_track_mva_cut",     "muon_prompt_==1 && track_mva_>0.5",      "goff");
    // vertex
  ch_PU200_prompt->Draw("muon_pt_>>h_PU200_prompt_vtx",             "muon_prompt_==1", 				      "goff");
  ch_PU200_prompt->Draw("muon_pt_>>h_PU200_prompt_vtx_time",        "muon_prompt_==1 && vtx_time_!=0",      "goff");
  ch_PU200_prompt->Draw("muon_pt_>>h_PU200_prompt_vtx_time_err",    "muon_prompt_==1 && vtx_time_err_>0", "goff");
  ch_PU200_prompt->Draw("track_pt_>>h_PU200_prompt_tot_vtx_track",     "muon_prompt_==1 && vtx_time_err_>0 && track_time_err_>0",      "goff");
    // total
  ch_PU200_prompt->Draw("track_pt_>>h_PU200_prompt_tot_mva_cut",     "muon_prompt_==1 && muon_time_err_>0 && track_time_err_>0",      "goff");
  // PU200 nonprompt
    // muon
  ch_PU200_nonprompt->Draw("muon_pt_>>h_PU200_nonprompt_muon",             "muon_prompt_==0", 				        "goff");
  ch_PU200_nonprompt->Draw("muon_pt_>>h_PU200_nonprompt_muon_mva",         "muon_prompt_==0 && muon_mva_>0",        "goff");
  ch_PU200_nonprompt->Draw("muon_pt_>>h_PU200_nonprompt_muon_time",        "muon_prompt_==0 && muon_time_!=0",       "goff");
  ch_PU200_nonprompt->Draw("muon_pt_>>h_PU200_nonprompt_muon_time_err",    "muon_prompt_==0 && muon_time_err_i_>0", "goff");
  ch_PU200_nonprompt->Draw("muon_pt_>>h_PU200_nonprompt_muon_mva_cut",     "muon_prompt_==0 && muon_mva_>0.5",      "goff");
    // vertex
  ch_PU200_nonprompt->Draw("muon_pt_>>h_PU200_nonprompt_vtx",             "muon_prompt_==0", 				      "goff");
  ch_PU200_nonprompt->Draw("muon_pt_>>h_PU200_nonprompt_vtx_time",        "muon_prompt_==0 && vtx_time_!=0",      "goff");
  ch_PU200_nonprompt->Draw("muon_pt_>>h_PU200_nonprompt_vtx_time_err",    "muon_prompt_==0 && vtx_time_err_>0", "goff");
  ch_PU200_nonprompt->Draw("track_pt_>>h_PU200_nonprompt_tot_vtx_track",     "muon_prompt_==0 && vtx_time_err_>0 && track_time_err_>0",      "goff");
    // track
  ch_PU200_nonprompt->Draw("track_pt_>>h_PU200_nonprompt_track",             "muon_prompt_==0", 				       "goff");
  ch_PU200_nonprompt->Draw("track_pt_>>h_PU200_nonprompt_track_mva",         "muon_prompt_==0 && track_mva_>0",        "goff");
  ch_PU200_nonprompt->Draw("track_pt_>>h_PU200_nonprompt_track_time",        "muon_prompt_==0 && track_time_!=0",       "goff");
  ch_PU200_nonprompt->Draw("track_pt_>>h_PU200_nonprompt_track_time_err",    "muon_prompt_==0 && track_time_err_i_>0", "goff");
  ch_PU200_nonprompt->Draw("track_pt_>>h_PU200_nonprompt_track_mva_cut",     "muon_prompt_==0 && track_mva_>0.5",      "goff");
    // total
  ch_PU200_nonprompt->Draw("track_pt_>>h_PU200_nonprompt_tot_mva_cut",     "muon_prompt_==0 && muon_time_err_>0 && track_time_err_>0",      "goff");

  // noPU prompt
    // muon
  ch_noPU_prompt->Draw("muon_pt_>>h_noPU_prompt_muon",             "muon_prompt_==1", 				        "goff");
  ch_noPU_prompt->Draw("muon_pt_>>h_noPU_prompt_muon_mva",         "muon_prompt_==1 && muon_mva_>0",        "goff");
  ch_noPU_prompt->Draw("muon_pt_>>h_noPU_prompt_muon_time",        "muon_prompt_==1 && muon_time_!=0",       "goff");
  ch_noPU_prompt->Draw("muon_pt_>>h_noPU_prompt_muon_time_err",    "muon_prompt_==1 && muon_time_err_i_>0", "goff");
  ch_noPU_prompt->Draw("muon_pt_>>h_noPU_prompt_muon_mva_cut",     "muon_prompt_==1 && muon_mva_>0.5",      "goff");
    // vertex
  ch_noPU_prompt->Draw("muon_pt_>>h_noPU_prompt_vtx",             "muon_prompt_==1", 				      "goff");
  ch_noPU_prompt->Draw("muon_pt_>>h_noPU_prompt_vtx_time",        "muon_prompt_==1 && vtx_time_!=0",      "goff");
  ch_noPU_prompt->Draw("muon_pt_>>h_noPU_prompt_vtx_time_err",    "muon_prompt_==1 && vtx_time_err_>0", "goff");
  ch_noPU_prompt->Draw("track_pt_>>h_noPU_prompt_tot_vtx_track",     "muon_prompt_==1 && vtx_time_err_>0 && track_time_err_>0",      "goff");
    // track
  ch_noPU_prompt->Draw("track_pt_>>h_noPU_prompt_track",             "muon_prompt_==1", 				       "goff");
  ch_noPU_prompt->Draw("track_pt_>>h_noPU_prompt_track_mva",         "muon_prompt_==1 && track_mva_>0",        "goff");
  ch_noPU_prompt->Draw("track_pt_>>h_noPU_prompt_track_time",        "muon_prompt_==1 && track_time_!=0",       "goff");
  ch_noPU_prompt->Draw("track_pt_>>h_noPU_prompt_track_time_err",    "muon_prompt_==1 && track_time_err_i_>0", "goff");
  ch_noPU_prompt->Draw("track_pt_>>h_noPU_prompt_track_mva_cut",     "muon_prompt_==1 && track_mva_>0.5",      "goff");
    // total
  ch_noPU_prompt->Draw("track_pt_>>h_noPU_prompt_tot_mva_cut",     "muon_prompt_==1 && muon_time_err_>0 && track_time_err_>0",      "goff");
  // noPU nonprompt
    // muon
  ch_noPU_nonprompt->Draw("muon_pt_>>h_noPU_nonprompt_muon",             "muon_prompt_==0", 				      "goff");
  ch_noPU_nonprompt->Draw("muon_pt_>>h_noPU_nonprompt_muon_mva",         "muon_prompt_==0 && muon_mva_>0",        "goff");
  ch_noPU_nonprompt->Draw("muon_pt_>>h_noPU_nonprompt_muon_time",        "muon_prompt_==0 && muon_time_!=0",       "goff");
  ch_noPU_nonprompt->Draw("muon_pt_>>h_noPU_nonprompt_muon_time_err",    "muon_prompt_==0 && muon_time_err_i_>0", "goff");
  ch_noPU_nonprompt->Draw("muon_pt_>>h_noPU_nonprompt_muon_mva_cut",     "muon_prompt_==0 && muon_mva_>0.5",      "goff");
    // vertex
  ch_noPU_nonprompt->Draw("muon_pt_>>h_noPU_nonprompt_vtx",             "muon_prompt_==0", 				      "goff");
  ch_noPU_nonprompt->Draw("muon_pt_>>h_noPU_nonprompt_vtx_time",        "muon_prompt_==0 && vtx_time_!=0",      "goff");
  ch_noPU_nonprompt->Draw("muon_pt_>>h_noPU_nonprompt_vtx_time_err",    "muon_prompt_==0 && vtx_time_err_>0", "goff");
  ch_noPU_nonprompt->Draw("track_pt_>>h_noPU_nonprompt_tot_vtx_track",     "muon_prompt_==0 && vtx_time_err_>0 && track_time_err_>0",      "goff");
    // track
  ch_noPU_nonprompt->Draw("track_pt_>>h_noPU_nonprompt_track",             "muon_prompt_==0", 				         "goff");
  ch_noPU_nonprompt->Draw("track_pt_>>h_noPU_nonprompt_track_mva",         "muon_prompt_==0 && track_mva_>0",        "goff");
  ch_noPU_nonprompt->Draw("track_pt_>>h_noPU_nonprompt_track_time",        "muon_prompt_==0 && track_time_!=0",       "goff");
  ch_noPU_nonprompt->Draw("track_pt_>>h_noPU_nonprompt_track_time_err",    "muon_prompt_==0 && track_time_err_i_>0", "goff");
  ch_noPU_nonprompt->Draw("track_pt_>>h_noPU_nonprompt_track_mva_cut",     "muon_prompt_==0 && track_mva_>0.5",      "goff");
    // total
  ch_noPU_nonprompt->Draw("track_pt_>>h_noPU_nonprompt_tot_mva_cut",     "muon_prompt_==0 && muon_time_err_>0 && track_time_err_>0",      "goff");


  // Get yields in each case
  int num_PU200_prompt_muon_tot,     num_PU200_prompt_muon_mva,     num_PU200_prompt_muon_time,     num_PU200_prompt_muon_time_err,     num_PU200_prompt_muon_mva_cut;
  int num_PU200_prompt_track_tot,    num_PU200_prompt_track_mva,    num_PU200_prompt_track_time,    num_PU200_prompt_track_time_err,    num_PU200_prompt_track_mva_cut,    num_PU200_prompt_tot_mva_cut;
  int num_PU200_nonprompt_muon_tot,  num_PU200_nonprompt_muon_mva,  num_PU200_nonprompt_muon_time,  num_PU200_nonprompt_muon_time_err,  num_PU200_nonprompt_muon_mva_cut;
  int num_PU200_nonprompt_track_tot, num_PU200_nonprompt_track_mva, num_PU200_nonprompt_track_time, num_PU200_nonprompt_track_time_err, num_PU200_nonprompt_track_mva_cut, num_PU200_nonprompt_tot_mva_cut;
  int num_noPU_prompt_muon_tot,      num_noPU_prompt_muon_mva,      num_noPU_prompt_muon_time,      num_noPU_prompt_muon_time_err,      num_noPU_prompt_muon_mva_cut;
  int num_noPU_prompt_track_tot,     num_noPU_prompt_track_mva,     num_noPU_prompt_track_time,     num_noPU_prompt_track_time_err,     num_noPU_prompt_track_mva_cut,     num_noPU_prompt_tot_mva_cut;
  int num_noPU_nonprompt_muon_tot,   num_noPU_nonprompt_muon_mva,   num_noPU_nonprompt_muon_time,   num_noPU_nonprompt_muon_time_err,   num_noPU_nonprompt_muon_mva_cut;
  int num_noPU_nonprompt_track_tot,  num_noPU_nonprompt_track_mva,  num_noPU_nonprompt_track_time,  num_noPU_nonprompt_track_time_err,  num_noPU_nonprompt_track_mva_cut,  num_noPU_nonprompt_tot_mva_cut;
    // vertex
  int num_PU200_prompt_vtx_tot,                                     num_PU200_prompt_vtx_time,      num_PU200_prompt_vtx_time_err,                                         num_PU200_prompt_tot_vtx_track;
  int num_PU200_nonprompt_vtx_tot,                                  num_PU200_nonprompt_vtx_time,   num_PU200_nonprompt_vtx_time_err,                                       num_PU200_nonprompt_tot_vtx_track;
  int num_noPU_prompt_vtx_tot,                                      num_noPU_prompt_vtx_time,       num_noPU_prompt_vtx_time_err,                                           num_noPU_prompt_tot_vtx_track;
  int num_noPU_nonprompt_vtx_tot,                                   num_noPU_nonprompt_vtx_time,    num_noPU_nonprompt_vtx_time_err,                                        num_noPU_nonprompt_tot_vtx_track;

  num_PU200_prompt_muon_tot = h_PU200_prompt_muon->Integral(0,-1); num_PU200_prompt_muon_mva = h_PU200_prompt_muon_mva->Integral(0,-1); num_PU200_prompt_muon_time = h_PU200_prompt_muon_time->Integral(0,-1); num_PU200_prompt_muon_time_err = h_PU200_prompt_muon_time_err->Integral(0,-1); num_PU200_prompt_muon_mva_cut = h_PU200_prompt_muon_mva_cut->Integral(0,-1);
  num_PU200_prompt_track_tot = h_PU200_prompt_track->Integral(0,-1); num_PU200_prompt_track_mva = h_PU200_prompt_track_mva->Integral(0,-1); num_PU200_prompt_track_time = h_PU200_prompt_track_time->Integral(0,-1); num_PU200_prompt_track_time_err = h_PU200_prompt_track_time_err->Integral(0,-1); num_PU200_prompt_track_mva_cut = h_PU200_prompt_track_mva_cut->Integral(0,-1); num_PU200_prompt_tot_mva_cut = h_PU200_prompt_tot_mva_cut->Integral(0,-1);
  num_PU200_nonprompt_muon_tot = h_PU200_nonprompt_muon->Integral(0,-1); num_PU200_nonprompt_muon_mva = h_PU200_nonprompt_muon_mva->Integral(0,-1); num_PU200_nonprompt_muon_time = h_PU200_nonprompt_muon_time->Integral(0,-1); num_PU200_nonprompt_muon_time_err = h_PU200_nonprompt_muon_time_err->Integral(0,-1); num_PU200_nonprompt_muon_mva_cut = h_PU200_nonprompt_muon_mva_cut->Integral(0,-1);
  num_PU200_nonprompt_track_tot = h_PU200_nonprompt_track->Integral(0,-1); num_PU200_nonprompt_track_mva = h_PU200_nonprompt_track_mva->Integral(0,-1); num_PU200_nonprompt_track_time = h_PU200_nonprompt_track_time->Integral(0,-1); num_PU200_nonprompt_track_time_err = h_PU200_nonprompt_track_time_err->Integral(0,-1); num_PU200_nonprompt_track_mva_cut = h_PU200_nonprompt_track_mva_cut->Integral(0,-1); num_PU200_nonprompt_tot_mva_cut = h_PU200_nonprompt_tot_mva_cut->Integral(0,-1);
  num_noPU_prompt_muon_tot = h_noPU_prompt_muon->Integral(0,-1); num_noPU_prompt_muon_mva = h_noPU_prompt_muon_mva->Integral(0,-1); num_noPU_prompt_muon_time = h_noPU_prompt_muon_time->Integral(0,-1); num_noPU_prompt_muon_time_err = h_noPU_prompt_muon_time_err->Integral(0,-1); num_noPU_prompt_muon_mva_cut = h_noPU_prompt_muon_mva_cut->Integral(0,-1);
  num_noPU_prompt_track_tot = h_noPU_prompt_track->Integral(0,-1); num_noPU_prompt_track_mva = h_noPU_prompt_track_mva->Integral(0,-1); num_noPU_prompt_track_time = h_noPU_prompt_track_time->Integral(0,-1); num_noPU_prompt_track_time_err = h_noPU_prompt_track_time_err->Integral(0,-1); num_noPU_prompt_track_mva_cut = h_noPU_prompt_track_mva_cut->Integral(0,-1); num_noPU_prompt_tot_mva_cut = h_noPU_prompt_tot_mva_cut->Integral(0,-1);
  num_noPU_nonprompt_muon_tot = h_noPU_nonprompt_muon->Integral(0,-1); num_noPU_nonprompt_muon_mva = h_noPU_nonprompt_muon_mva->Integral(0,-1); num_noPU_nonprompt_muon_time = h_noPU_nonprompt_muon_time->Integral(0,-1); num_noPU_nonprompt_muon_time_err = h_noPU_nonprompt_muon_time_err->Integral(0,-1); num_noPU_nonprompt_muon_mva_cut = h_noPU_nonprompt_muon_mva_cut->Integral(0,-1);
  num_noPU_nonprompt_track_tot = h_noPU_nonprompt_track->Integral(0,-1); num_noPU_nonprompt_track_mva = h_noPU_nonprompt_track_mva->Integral(0,-1); num_noPU_nonprompt_track_time = h_noPU_nonprompt_track_time->Integral(0,-1); num_noPU_nonprompt_track_time_err = h_noPU_nonprompt_track_time_err->Integral(0,-1); num_noPU_nonprompt_track_mva_cut = h_noPU_nonprompt_track_mva_cut->Integral(0,-1); num_noPU_nonprompt_tot_mva_cut = h_noPU_nonprompt_tot_mva_cut->Integral(0,-1);
    // vtx
  num_PU200_prompt_vtx_tot = h_PU200_prompt_vtx->Integral(0,-1); num_PU200_prompt_vtx_time = h_PU200_prompt_vtx_time->Integral(0,-1); num_PU200_prompt_vtx_time_err = h_PU200_prompt_vtx_time_err->Integral(0,-1); num_PU200_prompt_tot_vtx_track = h_PU200_prompt_tot_vtx_track->Integral(0,-1);
  num_PU200_nonprompt_vtx_tot = h_PU200_nonprompt_vtx->Integral(0,-1); num_PU200_nonprompt_vtx_time = h_PU200_nonprompt_vtx_time->Integral(0,-1); num_PU200_nonprompt_vtx_time_err = h_PU200_nonprompt_vtx_time_err->Integral(0,-1); num_PU200_nonprompt_tot_vtx_track = h_PU200_nonprompt_tot_vtx_track->Integral(0,-1);
  num_noPU_prompt_vtx_tot = h_noPU_prompt_vtx->Integral(0,-1); num_noPU_prompt_vtx_time = h_noPU_prompt_vtx_time->Integral(0,-1); num_noPU_prompt_vtx_time_err = h_noPU_prompt_vtx_time_err->Integral(0,-1); num_noPU_prompt_tot_vtx_track = h_noPU_prompt_tot_vtx_track->Integral(0,-1);
  num_noPU_nonprompt_vtx_tot = h_noPU_nonprompt_vtx->Integral(0,-1); num_noPU_nonprompt_vtx_time = h_noPU_nonprompt_vtx_time->Integral(0,-1); num_noPU_nonprompt_vtx_time_err = h_noPU_nonprompt_vtx_time_err->Integral(0,-1); num_noPU_nonprompt_tot_vtx_track = h_noPU_nonprompt_tot_vtx_track->Integral(0,-1);


  cout << endl;
  cout << endl;
  cout << "//////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "////Check on whether muon/track has MVA/time/tErr or passes MVA cut or not////" << endl;
  cout << "//////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "[PU200 prompt]" << endl;
  cout << "               |           MVA           |           time           |           tErr           |       MVA > 0.5 / total       |        MVA > 0.5 / MVA       |" << endl;
  cout << "      muon     | " 
	  << num_PU200_prompt_muon_mva << "/" << num_PU200_prompt_muon_tot << " = " << 1.*num_PU200_prompt_muon_mva/num_PU200_prompt_muon_tot << " % | " 
	  << num_PU200_prompt_muon_time << "/" << num_PU200_prompt_muon_tot << " = " << 1.*num_PU200_prompt_muon_time/num_PU200_prompt_muon_tot << " % | "
	  << num_PU200_prompt_muon_time_err << "/" << num_PU200_prompt_muon_tot << " = " << 1.*num_PU200_prompt_muon_time_err/num_PU200_prompt_muon_tot << " % | "
	  << num_PU200_prompt_muon_mva_cut << "/" << num_PU200_prompt_muon_tot << " = " << 1.*num_PU200_prompt_muon_mva_cut/num_PU200_prompt_muon_tot << " % | "
	  << num_PU200_prompt_muon_mva_cut << "/" << num_PU200_prompt_muon_mva << " = " << 1.*num_PU200_prompt_muon_mva_cut/num_PU200_prompt_muon_mva << " % | "
	  << endl;
  cout << "      track    | " 
	  << num_PU200_prompt_track_mva << "/" << num_PU200_prompt_track_tot << " = " << 1.*num_PU200_prompt_track_mva/num_PU200_prompt_track_tot << " % | " 
	  << num_PU200_prompt_track_time << "/" << num_PU200_prompt_track_tot << " = " << 1.*num_PU200_prompt_track_time/num_PU200_prompt_track_tot << " % | "
	  << num_PU200_prompt_track_time_err << "/" << num_PU200_prompt_track_tot << " = " << 1.*num_PU200_prompt_track_time_err/num_PU200_prompt_track_tot << " % | "
	  << num_PU200_prompt_track_mva_cut << "/" << num_PU200_prompt_track_tot << " = " << 1.*num_PU200_prompt_track_mva_cut/num_PU200_prompt_track_tot << " % | "
	  << num_PU200_prompt_track_mva_cut << "/" << num_PU200_prompt_track_mva << " = " << 1.*num_PU200_prompt_track_mva_cut/num_PU200_prompt_track_mva << " % | "
	  << endl;
  cout << "  muon&&track  | " 
	  << num_PU200_prompt_tot_mva_cut << "/" << num_PU200_prompt_track_tot << " = " << 1.*num_PU200_prompt_tot_mva_cut/num_PU200_prompt_track_tot << " % | " 
	  << endl;

  cout << endl;
  cout << "[PU200 nonprompt]" << endl;
  cout << "               |           MVA           |           time           |           tErr           |       MVA > 0.5 / total       |        MVA > 0.5 / MVA       |" << endl;
  cout << "      muon     | " 
	  << num_PU200_nonprompt_muon_mva << "/" << num_PU200_nonprompt_muon_tot << " = " << 1.*num_PU200_nonprompt_muon_mva/num_PU200_nonprompt_muon_tot << " % | " 
	  << num_PU200_nonprompt_muon_time << "/" << num_PU200_nonprompt_muon_tot << " = " << 1.*num_PU200_nonprompt_muon_time/num_PU200_nonprompt_muon_tot << " % | "
	  << num_PU200_nonprompt_muon_time_err << "/" << num_PU200_nonprompt_muon_tot << " = " << 1.*num_PU200_nonprompt_muon_time_err/num_PU200_nonprompt_muon_tot << " % | "
	  << num_PU200_nonprompt_muon_mva_cut << "/" << num_PU200_nonprompt_muon_tot << " = " << 1.*num_PU200_nonprompt_muon_mva_cut/num_PU200_nonprompt_muon_tot << " % | "
	  << num_PU200_nonprompt_muon_mva_cut << "/" << num_PU200_nonprompt_muon_mva << " = " << 1.*num_PU200_nonprompt_muon_mva_cut/num_PU200_nonprompt_muon_mva << " % | "
	  << endl;
  cout << "      track    | " 
	  << num_PU200_nonprompt_track_mva << "/" << num_PU200_nonprompt_track_tot << " = " << 1.*num_PU200_nonprompt_track_mva/num_PU200_nonprompt_track_tot << " % | " 
	  << num_PU200_nonprompt_track_time << "/" << num_PU200_nonprompt_track_tot << " = " << 1.*num_PU200_nonprompt_track_time/num_PU200_nonprompt_track_tot << " % | "
	  << num_PU200_nonprompt_track_time_err << "/" << num_PU200_nonprompt_track_tot << " = " << 1.*num_PU200_nonprompt_track_time_err/num_PU200_nonprompt_track_tot << " % | "
	  << num_PU200_nonprompt_track_mva_cut << "/" << num_PU200_nonprompt_track_tot << " = " << 1.*num_PU200_nonprompt_track_mva_cut/num_PU200_nonprompt_track_tot << " % | "
	  << num_PU200_nonprompt_track_mva_cut << "/" << num_PU200_nonprompt_track_mva << " = " << 1.*num_PU200_nonprompt_track_mva_cut/num_PU200_nonprompt_track_mva << " % | "
	  << endl;
  cout << "  muon&&track  | " 
	  << num_PU200_nonprompt_tot_mva_cut << "/" << num_PU200_nonprompt_track_tot << " = " << 1.*num_PU200_nonprompt_tot_mva_cut/num_PU200_nonprompt_track_tot << " % | " 
	  << endl;

  cout << endl;
  cout << "[noPU prompt]" << endl;
  cout << "               |           MVA           |           time           |           tErr           |       MVA > 0.5 / total       |        MVA > 0.5 / MVA       |" << endl;
  cout << "      muon     | " 
	  << num_noPU_prompt_muon_mva << "/" << num_noPU_prompt_muon_tot << " = " << 1.*num_noPU_prompt_muon_mva/num_noPU_prompt_muon_tot << " % | " 
	  << num_noPU_prompt_muon_time << "/" << num_noPU_prompt_muon_tot << " = " << 1.*num_noPU_prompt_muon_time/num_noPU_prompt_muon_tot << " % | "
	  << num_noPU_prompt_muon_time_err << "/" << num_noPU_prompt_muon_tot << " = " << 1.*num_noPU_prompt_muon_time_err/num_noPU_prompt_muon_tot << " % | "
	  << num_noPU_prompt_muon_mva_cut << "/" << num_noPU_prompt_muon_tot << " = " << 1.*num_noPU_prompt_muon_mva_cut/num_noPU_prompt_muon_tot << " % | "
	  << num_noPU_prompt_muon_mva_cut << "/" << num_noPU_prompt_muon_mva << " = " << 1.*num_noPU_prompt_muon_mva_cut/num_noPU_prompt_muon_mva << " % | "
	  << endl;
  cout << "      track    | " 
	  << num_noPU_prompt_track_mva << "/" << num_noPU_prompt_track_tot << " = " << 1.*num_noPU_prompt_track_mva/num_noPU_prompt_track_tot << " % | " 
	  << num_noPU_prompt_track_time << "/" << num_noPU_prompt_track_tot << " = " << 1.*num_noPU_prompt_track_time/num_noPU_prompt_track_tot << " % | "
	  << num_noPU_prompt_track_time_err << "/" << num_noPU_prompt_track_tot << " = " << 1.*num_noPU_prompt_track_time_err/num_noPU_prompt_track_tot << " % | "
	  << num_noPU_prompt_track_mva_cut << "/" << num_noPU_prompt_track_tot << " = " << 1.*num_noPU_prompt_track_mva_cut/num_noPU_prompt_track_tot << " % | "
	  << num_noPU_prompt_track_mva_cut << "/" << num_noPU_prompt_track_mva << " = " << 1.*num_noPU_prompt_track_mva_cut/num_noPU_prompt_track_mva << " % | "
	  << endl;
  cout << "  muon&&track  | " 
	  << num_noPU_prompt_tot_mva_cut << "/" << num_noPU_prompt_track_tot << " = " << 1.*num_noPU_prompt_tot_mva_cut/num_noPU_prompt_track_tot << " % | " 
	  << endl;

  cout << endl;
  cout << "[noPU nonprompt]" << endl;
  cout << "               |           MVA           |           time           |           tErr           |       MVA > 0.5 / total       |        MVA > 0.5 / MVA       |" << endl;
  cout << "      muon     | " 
	  << num_noPU_nonprompt_muon_mva << "/" << num_noPU_nonprompt_muon_tot << " = " << 1.*num_noPU_nonprompt_muon_mva/num_noPU_nonprompt_muon_tot << " % | " 
	  << num_noPU_nonprompt_muon_time << "/" << num_noPU_nonprompt_muon_tot << " = " << 1.*num_noPU_nonprompt_muon_time/num_noPU_nonprompt_muon_tot << " % | "
	  << num_noPU_nonprompt_muon_time_err << "/" << num_noPU_nonprompt_muon_tot << " = " << 1.*num_noPU_nonprompt_muon_time_err/num_noPU_nonprompt_muon_tot << " % | "
	  << num_noPU_nonprompt_muon_mva_cut << "/" << num_noPU_nonprompt_muon_tot << " = " << 1.*num_noPU_nonprompt_muon_mva_cut/num_noPU_nonprompt_muon_tot << " % | "
	  << num_noPU_nonprompt_muon_mva_cut << "/" << num_noPU_nonprompt_muon_mva << " = " << 1.*num_noPU_nonprompt_muon_mva_cut/num_noPU_nonprompt_muon_mva << " % | "
	  << endl;
  cout << "      track    | " 
	  << num_noPU_nonprompt_track_mva << "/" << num_noPU_nonprompt_track_tot << " = " << 1.*num_noPU_nonprompt_track_mva/num_noPU_nonprompt_track_tot << " % | " 
	  << num_noPU_nonprompt_track_time << "/" << num_noPU_nonprompt_track_tot << " = " << 1.*num_noPU_nonprompt_track_time/num_noPU_nonprompt_track_tot << " % | "
	  << num_noPU_nonprompt_track_time_err << "/" << num_noPU_nonprompt_track_tot << " = " << 1.*num_noPU_nonprompt_track_time_err/num_noPU_nonprompt_track_tot << " % | "
	  << num_noPU_nonprompt_track_mva_cut << "/" << num_noPU_nonprompt_track_tot << " = " << 1.*num_noPU_nonprompt_track_mva_cut/num_noPU_nonprompt_track_tot << " % | "
	  << num_noPU_nonprompt_track_mva_cut << "/" << num_noPU_nonprompt_track_mva << " = " << 1.*num_noPU_nonprompt_track_mva_cut/num_noPU_nonprompt_track_mva << " % | "
	  << endl;
  cout << "  muon&&track  | " 
	  << num_noPU_nonprompt_tot_mva_cut << "/" << num_noPU_nonprompt_track_tot << " = " << 1.*num_noPU_nonprompt_tot_mva_cut/num_noPU_nonprompt_track_tot << " % | " 
	  << endl;
    // vtx
  cout << endl;
  cout << endl;
  cout << "//////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "             //// Check on whether vtx has time/tErr ////" << endl;
  cout << "//////////////////////////////////////////////////////////////////////////////" << endl;
  cout << "[PU200 prompt]" << endl;
  cout << "               |           time           |           tErr           |" << endl;
  cout << "      muon     | " 
	  << num_PU200_prompt_vtx_time << "/" << num_PU200_prompt_vtx_tot << " = " << 1.*num_PU200_prompt_vtx_time/num_PU200_prompt_vtx_tot << " % | "
	  << num_PU200_prompt_vtx_time_err << "/" << num_PU200_prompt_vtx_tot << " = " << 1.*num_PU200_prompt_vtx_time_err/num_PU200_prompt_vtx_tot << " % | "
	  << endl;
  cout << "   vtx&&track  | " 
	  << num_PU200_prompt_tot_vtx_track << "/" << num_PU200_prompt_track_tot << " = " << 1.*num_PU200_prompt_tot_vtx_track/num_PU200_prompt_track_tot << " % | " 
	  << endl;
  cout << "[PU200 nonprompt]" << endl;
  cout << "               |           time           |           tErr           |" << endl;
  cout << "      muon     | " 
	  << num_PU200_nonprompt_vtx_time << "/" << num_PU200_nonprompt_vtx_tot << " = " << 1.*num_PU200_nonprompt_vtx_time/num_PU200_nonprompt_vtx_tot << " % | "
	  << num_PU200_nonprompt_vtx_time_err << "/" << num_PU200_nonprompt_vtx_tot << " = " << 1.*num_PU200_nonprompt_vtx_time_err/num_PU200_nonprompt_vtx_tot << " % | "
	  << endl;
  cout << "   vtx&&track  | " 
	  << num_PU200_nonprompt_tot_vtx_track << "/" << num_PU200_nonprompt_track_tot << " = " << 1.*num_PU200_nonprompt_tot_vtx_track/num_PU200_nonprompt_track_tot << " % | " 
	  << endl;
  cout << "[noPU prompt]" << endl;
  cout << "               |           time           |           tErr           |" << endl;
  cout << "      muon     | " 
	  << num_noPU_prompt_vtx_time << "/" << num_noPU_prompt_vtx_tot << " = " << 1.*num_noPU_prompt_vtx_time/num_noPU_prompt_vtx_tot << " % | "
	  << num_noPU_prompt_vtx_time_err << "/" << num_noPU_prompt_vtx_tot << " = " << 1.*num_noPU_prompt_vtx_time_err/num_noPU_prompt_vtx_tot << " % | "
	  << endl;
  cout << "   vtx&&track  | " 
	  << num_noPU_prompt_tot_vtx_track << "/" << num_noPU_prompt_track_tot << " = " << 1.*num_noPU_prompt_tot_vtx_track/num_noPU_prompt_track_tot << " % | " 
	  << endl;
  cout << "[noPU nonprompt]" << endl;
  cout << "               |           time           |           tErr           |" << endl;
  cout << "      muon     | " 
	  << num_noPU_nonprompt_vtx_time << "/" << num_noPU_nonprompt_vtx_tot << " = " << 1.*num_noPU_nonprompt_vtx_time/num_noPU_nonprompt_vtx_tot << " % | "
	  << num_noPU_nonprompt_vtx_time_err << "/" << num_noPU_nonprompt_vtx_tot << " = " << 1.*num_noPU_nonprompt_vtx_time_err/num_noPU_nonprompt_vtx_tot << " % | "
	  << endl;
  cout << "   vtx&&track  | " 
	  << num_noPU_nonprompt_tot_vtx_track << "/" << num_noPU_nonprompt_track_tot << " = " << 1.*num_noPU_nonprompt_tot_vtx_track/num_noPU_nonprompt_track_tot << " % | " 
	  << endl;

}

void evtId_pvtrk() {
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


  TH1D* h_PU200_prompt    = new TH1D("h_PU200_prompt",    "h_PU200_prompt",    300, 0, 300);
  TH1D* h_PU200_nonprompt = new TH1D("h_PU200_nonprompt", "h_PU200_nonprompt", 300, 0, 300);
  TH1D* h_noPU_prompt     = new TH1D("h_noPU_prompt",     "h_noPU_prompt",     300, 0, 300);
  TH1D* h_noPU_nonprompt  = new TH1D("h_noPU_nonprompt",  "h_noPU_nonprompt",  300, 0, 300);
  // PU200
    // Barrel
  ch_PU200_prompt->Draw("track_evtId_>>h_PU200_prompt",       "muon_prompt_==1 && track_type_==0", "goff");
  ch_PU200_nonprompt->Draw("track_evtId_>>h_PU200_nonprompt", "muon_prompt_==0 && track_type_==0", "goff");
  ch_noPU_prompt->Draw("track_evtId_>>h_noPU_prompt",         "muon_prompt_==1 && track_type_==0", "goff");
  ch_noPU_nonprompt->Draw("track_evtId_>>h_noPU_nonprompt",   "muon_prompt_==0 && track_type_==0", "goff");


  TCanvas* c_PU200_prompt_evtId = new TCanvas("c_PU200_prompt_evtId", "c_PU200_prompt_evtId", 1500, 1500);
  c_PU200_prompt_evtId->cd();
  c_PU200_prompt_evtId->SetLogy();
  h_PU200_prompt->GetXaxis()->SetTitle("track_eventId");
  h_PU200_prompt->GetYaxis()->SetTitle("Counts");
  h_PU200_prompt->SetTitle("");
  h_PU200_prompt->Draw("hist");
  c_PU200_prompt_evtId->Print("plots/evtId_PU200_prompt_pvtrk.pdf");

  TCanvas* c_PU200_nonprompt_evtId = new TCanvas("c_PU200_nonprompt_evtId", "c_PU200_nonprompt_evtId", 1500, 1500);
  c_PU200_nonprompt_evtId->cd();
  c_PU200_nonprompt_evtId->SetLogy();
  h_PU200_nonprompt->GetXaxis()->SetTitle("track_eventId");
  h_PU200_nonprompt->GetYaxis()->SetTitle("Counts");
  h_PU200_nonprompt->SetTitle("");
  h_PU200_nonprompt->Draw("hist");
  c_PU200_nonprompt_evtId->Print("plots/evtId_PU200_nonprompt_pvtrk.pdf");

  TCanvas* c_noPU_prompt_evtId = new TCanvas("c_noPU_prompt_evtId", "c_noPU_prompt_evtId", 1500, 1500);
  c_noPU_prompt_evtId->cd();
  c_noPU_prompt_evtId->SetLogy();
  h_noPU_prompt->GetXaxis()->SetTitle("track_eventId");
  h_noPU_prompt->GetYaxis()->SetTitle("Counts");
  h_noPU_prompt->SetTitle("");
  h_noPU_prompt->Draw("hist");
  c_noPU_prompt_evtId->Print("plots/evtId_noPU_prompt_pvtrk.pdf");

  TCanvas* c_noPU_nonprompt_evtId = new TCanvas("c_noPU_nonprompt_evtId", "c_noPU_nonprompt_evtId", 1500, 1500);
  c_noPU_nonprompt_evtId->cd();
  c_noPU_nonprompt_evtId->SetLogy();
  h_noPU_nonprompt->GetXaxis()->SetTitle("track_eventId");
  h_noPU_nonprompt->GetYaxis()->SetTitle("Counts");
  h_noPU_nonprompt->SetTitle("");
  h_noPU_nonprompt->Draw("hist");
  c_noPU_nonprompt_evtId->Print("plots/evtId_noPU_nonprompt_pvtrk.pdf");

 
}

void pt_eta_phi() {
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


  // PU200 prompt
    // muon
  TH1D* h_PU200_prompt_muon_pt     = new TH1D("h_PU200_prompt_muon_pt",     "h_PU200_prompt_muon_pt",    20, 0, 100);
  TH1D* h_PU200_prompt_muon_pt_EB  = new TH1D("h_PU200_prompt_muon_pt_EB",  "h_PU200_prompt_muon_pt_EB", 20, 0, 100);
  TH1D* h_PU200_prompt_muon_pt_EE  = new TH1D("h_PU200_prompt_muon_pt_EE",  "h_PU200_prompt_muon_pt_EE", 20, 0, 100);
  TH1D* h_PU200_prompt_muon_eta    = new TH1D("h_PU200_prompt_muon_eta",    "h_PU200_prompt_muon_eta",    48, -2.4, 2.4);
  TH1D* h_PU200_prompt_muon_eta_EB = new TH1D("h_PU200_prompt_muon_eta_EB", "h_PU200_prompt_muon_eta_EB", 48, -2.4, 2.4);
  TH1D* h_PU200_prompt_muon_eta_EE = new TH1D("h_PU200_prompt_muon_eta_EE", "h_PU200_prompt_muon_eta_EE", 48, -2.4, 2.4);
  TH1D* h_PU200_prompt_muon_phi    = new TH1D("h_PU200_prompt_muon_phi",    "h_PU200_prompt_muon_phi",    64, -3.2, 3.2);
  TH1D* h_PU200_prompt_muon_phi_EB = new TH1D("h_PU200_prompt_muon_phi_EB", "h_PU200_prompt_muon_phi_EB", 64, -3.2, 3.2);
  TH1D* h_PU200_prompt_muon_phi_EE = new TH1D("h_PU200_prompt_muon_phi_EE", "h_PU200_prompt_muon_phi_EE", 64, -3.2, 3.2);
    // track
  TH1D* h_PU200_prompt_track_pt     = new TH1D("h_PU200_prompt_track_pt",     "h_PU200_prompt_track_pt",    20, 0, 100);
  TH1D* h_PU200_prompt_track_pt_EB  = new TH1D("h_PU200_prompt_track_pt_EB",  "h_PU200_prompt_track_pt_EB", 20, 0, 100);
  TH1D* h_PU200_prompt_track_pt_EE  = new TH1D("h_PU200_prompt_track_pt_EE",  "h_PU200_prompt_track_pt_EE", 20, 0, 100);
  TH1D* h_PU200_prompt_track_eta    = new TH1D("h_PU200_prompt_track_eta",    "h_PU200_prompt_track_eta",    48, -2.4, 2.4);
  TH1D* h_PU200_prompt_track_eta_EB = new TH1D("h_PU200_prompt_track_eta_EB", "h_PU200_prompt_track_eta_EB", 48, -2.4, 2.4);
  TH1D* h_PU200_prompt_track_eta_EE = new TH1D("h_PU200_prompt_track_eta_EE", "h_PU200_prompt_track_eta_EE", 48, -2.4, 2.4);
  TH1D* h_PU200_prompt_track_phi    = new TH1D("h_PU200_prompt_track_phi",    "h_PU200_prompt_track_phi",    64, -3.2, 3.2);
  TH1D* h_PU200_prompt_track_phi_EB = new TH1D("h_PU200_prompt_track_phi_EB", "h_PU200_prompt_track_phi_EB", 64, -3.2, 3.2);
  TH1D* h_PU200_prompt_track_phi_EE = new TH1D("h_PU200_prompt_track_phi_EE", "h_PU200_prompt_track_phi_EE", 64, -3.2, 3.2);
  // PU200 nonprompt
    // muon
  TH1D* h_PU200_nonprompt_muon_pt     = new TH1D("h_PU200_nonprompt_muon_pt",     "h_PU200_nonprompt_muon_pt",    20, 0, 100);
  TH1D* h_PU200_nonprompt_muon_pt_EB  = new TH1D("h_PU200_nonprompt_muon_pt_EB",  "h_PU200_nonprompt_muon_pt_EB", 20, 0, 100);
  TH1D* h_PU200_nonprompt_muon_pt_EE  = new TH1D("h_PU200_nonprompt_muon_pt_EE",  "h_PU200_nonprompt_muon_pt_EE", 20, 0, 100);
  TH1D* h_PU200_nonprompt_muon_eta    = new TH1D("h_PU200_nonprompt_muon_eta",    "h_PU200_nonprompt_muon_eta",    48, -2.4, 2.4);
  TH1D* h_PU200_nonprompt_muon_eta_EB = new TH1D("h_PU200_nonprompt_muon_eta_EB", "h_PU200_nonprompt_muon_eta_EB", 48, -2.4, 2.4);
  TH1D* h_PU200_nonprompt_muon_eta_EE = new TH1D("h_PU200_nonprompt_muon_eta_EE", "h_PU200_nonprompt_muon_eta_EE", 48, -2.4, 2.4);
  TH1D* h_PU200_nonprompt_muon_phi    = new TH1D("h_PU200_nonprompt_muon_phi",    "h_PU200_nonprompt_muon_phi",    64, -3.2, 3.2);
  TH1D* h_PU200_nonprompt_muon_phi_EB = new TH1D("h_PU200_nonprompt_muon_phi_EB", "h_PU200_nonprompt_muon_phi_EB", 64, -3.2, 3.2);
  TH1D* h_PU200_nonprompt_muon_phi_EE = new TH1D("h_PU200_nonprompt_muon_phi_EE", "h_PU200_nonprompt_muon_phi_EE", 64, -3.2, 3.2);
    // track
  TH1D* h_PU200_nonprompt_track_pt     = new TH1D("h_PU200_nonprompt_track_pt",     "h_PU200_nonprompt_track_pt",    20, 0, 100);
  TH1D* h_PU200_nonprompt_track_pt_EB  = new TH1D("h_PU200_nonprompt_track_pt_EB",  "h_PU200_nonprompt_track_pt_EB", 20, 0, 100);
  TH1D* h_PU200_nonprompt_track_pt_EE  = new TH1D("h_PU200_nonprompt_track_pt_EE",  "h_PU200_nonprompt_track_pt_EE", 20, 0, 100);
  TH1D* h_PU200_nonprompt_track_eta    = new TH1D("h_PU200_nonprompt_track_eta",    "h_PU200_nonprompt_track_eta",    48, -2.4, 2.4);
  TH1D* h_PU200_nonprompt_track_eta_EB = new TH1D("h_PU200_nonprompt_track_eta_EB", "h_PU200_nonprompt_track_eta_EB", 48, -2.4, 2.4);
  TH1D* h_PU200_nonprompt_track_eta_EE = new TH1D("h_PU200_nonprompt_track_eta_EE", "h_PU200_nonprompt_track_eta_EE", 48, -2.4, 2.4);
  TH1D* h_PU200_nonprompt_track_phi    = new TH1D("h_PU200_nonprompt_track_phi",    "h_PU200_nonprompt_track_phi",    64, -3.2, 3.2);
  TH1D* h_PU200_nonprompt_track_phi_EB = new TH1D("h_PU200_nonprompt_track_phi_EB", "h_PU200_nonprompt_track_phi_EB", 64, -3.2, 3.2);
  TH1D* h_PU200_nonprompt_track_phi_EE = new TH1D("h_PU200_nonprompt_track_phi_EE", "h_PU200_nonprompt_track_phi_EE", 64, -3.2, 3.2);
  // noPU prompt
    // muon
  TH1D* h_noPU_prompt_muon_pt     = new TH1D("h_noPU_prompt_muon_pt",     "h_noPU_prompt_muon_pt",    20, 0, 100);
  TH1D* h_noPU_prompt_muon_pt_EB  = new TH1D("h_noPU_prompt_muon_pt_EB",  "h_noPU_prompt_muon_pt_EB", 20, 0, 100);
  TH1D* h_noPU_prompt_muon_pt_EE  = new TH1D("h_noPU_prompt_muon_pt_EE",  "h_noPU_prompt_muon_pt_EE", 20, 0, 100);
  TH1D* h_noPU_prompt_muon_eta    = new TH1D("h_noPU_prompt_muon_eta",    "h_noPU_prompt_muon_eta",    48, -2.4, 2.4);
  TH1D* h_noPU_prompt_muon_eta_EB = new TH1D("h_noPU_prompt_muon_eta_EB", "h_noPU_prompt_muon_eta_EB", 48, -2.4, 2.4);
  TH1D* h_noPU_prompt_muon_eta_EE = new TH1D("h_noPU_prompt_muon_eta_EE", "h_noPU_prompt_muon_eta_EE", 48, -2.4, 2.4);
  TH1D* h_noPU_prompt_muon_phi    = new TH1D("h_noPU_prompt_muon_phi",    "h_noPU_prompt_muon_phi",    64, -3.2, 3.2);
  TH1D* h_noPU_prompt_muon_phi_EB = new TH1D("h_noPU_prompt_muon_phi_EB", "h_noPU_prompt_muon_phi_EB", 64, -3.2, 3.2);
  TH1D* h_noPU_prompt_muon_phi_EE = new TH1D("h_noPU_prompt_muon_phi_EE", "h_noPU_prompt_muon_phi_EE", 64, -3.2, 3.2);
    // track
  TH1D* h_noPU_prompt_track_pt     = new TH1D("h_noPU_prompt_track_pt",     "h_noPU_prompt_track_pt",    20, 0, 100);
  TH1D* h_noPU_prompt_track_pt_EB  = new TH1D("h_noPU_prompt_track_pt_EB",  "h_noPU_prompt_track_pt_EB", 20, 0, 100);
  TH1D* h_noPU_prompt_track_pt_EE  = new TH1D("h_noPU_prompt_track_pt_EE",  "h_noPU_prompt_track_pt_EE", 20, 0, 100);
  TH1D* h_noPU_prompt_track_eta    = new TH1D("h_noPU_prompt_track_eta",    "h_noPU_prompt_track_eta",    48, -2.4, 2.4);
  TH1D* h_noPU_prompt_track_eta_EB = new TH1D("h_noPU_prompt_track_eta_EB", "h_noPU_prompt_track_eta_EB", 48, -2.4, 2.4);
  TH1D* h_noPU_prompt_track_eta_EE = new TH1D("h_noPU_prompt_track_eta_EE", "h_noPU_prompt_track_eta_EE", 48, -2.4, 2.4);
  TH1D* h_noPU_prompt_track_phi    = new TH1D("h_noPU_prompt_track_phi",    "h_noPU_prompt_track_phi",    64, -3.2, 3.2);
  TH1D* h_noPU_prompt_track_phi_EB = new TH1D("h_noPU_prompt_track_phi_EB", "h_noPU_prompt_track_phi_EB", 64, -3.2, 3.2);
  TH1D* h_noPU_prompt_track_phi_EE = new TH1D("h_noPU_prompt_track_phi_EE", "h_noPU_prompt_track_phi_EE", 64, -3.2, 3.2);
  // noPU nonprompt
    // muon
  TH1D* h_noPU_nonprompt_muon_pt     = new TH1D("h_noPU_nonprompt_muon_pt",     "h_noPU_nonprompt_muon_pt",    20, 0, 100);
  TH1D* h_noPU_nonprompt_muon_pt_EB  = new TH1D("h_noPU_nonprompt_muon_pt_EB",  "h_noPU_nonprompt_muon_pt_EB", 20, 0, 100);
  TH1D* h_noPU_nonprompt_muon_pt_EE  = new TH1D("h_noPU_nonprompt_muon_pt_EE",  "h_noPU_nonprompt_muon_pt_EE", 20, 0, 100);
  TH1D* h_noPU_nonprompt_muon_eta    = new TH1D("h_noPU_nonprompt_muon_eta",    "h_noPU_nonprompt_muon_eta",    48, -2.4, 2.4);
  TH1D* h_noPU_nonprompt_muon_eta_EB = new TH1D("h_noPU_nonprompt_muon_eta_EB", "h_noPU_nonprompt_muon_eta_EB", 48, -2.4, 2.4);
  TH1D* h_noPU_nonprompt_muon_eta_EE = new TH1D("h_noPU_nonprompt_muon_eta_EE", "h_noPU_nonprompt_muon_eta_EE", 48, -2.4, 2.4);
  TH1D* h_noPU_nonprompt_muon_phi    = new TH1D("h_noPU_nonprompt_muon_phi",    "h_noPU_nonprompt_muon_phi",    64, -3.2, 3.2);
  TH1D* h_noPU_nonprompt_muon_phi_EB = new TH1D("h_noPU_nonprompt_muon_phi_EB", "h_noPU_nonprompt_muon_phi_EB", 64, -3.2, 3.2);
  TH1D* h_noPU_nonprompt_muon_phi_EE = new TH1D("h_noPU_nonprompt_muon_phi_EE", "h_noPU_nonprompt_muon_phi_EE", 64, -3.2, 3.2);
    // track
  TH1D* h_noPU_nonprompt_track_pt     = new TH1D("h_noPU_nonprompt_track_pt",     "h_noPU_nonprompt_track_pt",    20, 0, 100);
  TH1D* h_noPU_nonprompt_track_pt_EB  = new TH1D("h_noPU_nonprompt_track_pt_EB",  "h_noPU_nonprompt_track_pt_EB", 20, 0, 100);
  TH1D* h_noPU_nonprompt_track_pt_EE  = new TH1D("h_noPU_nonprompt_track_pt_EE",  "h_noPU_nonprompt_track_pt_EE", 20, 0, 100);
  TH1D* h_noPU_nonprompt_track_eta    = new TH1D("h_noPU_nonprompt_track_eta",    "h_noPU_nonprompt_track_eta",    48, -2.4, 2.4);
  TH1D* h_noPU_nonprompt_track_eta_EB = new TH1D("h_noPU_nonprompt_track_eta_EB", "h_noPU_nonprompt_track_eta_EB", 48, -2.4, 2.4);
  TH1D* h_noPU_nonprompt_track_eta_EE = new TH1D("h_noPU_nonprompt_track_eta_EE", "h_noPU_nonprompt_track_eta_EE", 48, -2.4, 2.4);
  TH1D* h_noPU_nonprompt_track_phi    = new TH1D("h_noPU_nonprompt_track_phi",    "h_noPU_nonprompt_track_phi",    64, -3.2, 3.2);
  TH1D* h_noPU_nonprompt_track_phi_EB = new TH1D("h_noPU_nonprompt_track_phi_EB", "h_noPU_nonprompt_track_phi_EB", 64, -3.2, 3.2);
  TH1D* h_noPU_nonprompt_track_phi_EE = new TH1D("h_noPU_nonprompt_track_phi_EE", "h_noPU_nonprompt_track_phi_EE", 64, -3.2, 3.2);


  // PU200 prompt
    // muon
  ch_PU200_prompt->Draw("muon_pt_>>h_PU200_prompt_muon_pt",                  "muon_prompt_==1",                      "goff");
  ch_PU200_prompt->Draw("muon_pt_>>h_PU200_prompt_muon_pt_EB",               "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_PU200_prompt->Draw("muon_pt_>>h_PU200_prompt_muon_pt_EE",               "muon_prompt_==1 && muon_isBarrel_==0", "goff");
  ch_PU200_prompt->Draw("muon_eta_>>h_PU200_prompt_muon_eta",                "muon_prompt_==1",                      "goff");
  ch_PU200_prompt->Draw("muon_eta_>>h_PU200_prompt_muon_eta_EB",             "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_PU200_prompt->Draw("muon_eta_>>h_PU200_prompt_muon_eta_EE",             "muon_prompt_==1 && muon_isBarrel_==0", "goff");
  ch_PU200_prompt->Draw("muon_phi_>>h_PU200_prompt_muon_phi",                "muon_prompt_==1",                      "goff");
  ch_PU200_prompt->Draw("muon_phi_>>h_PU200_prompt_muon_phi_EB",             "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_PU200_prompt->Draw("muon_phi_>>h_PU200_prompt_muon_phi_EE",             "muon_prompt_==1 && muon_isBarrel_==0", "goff");
    // track
  ch_PU200_prompt->Draw("track_pt_>>h_PU200_prompt_track_pt",                  "muon_prompt_==1",                      "goff");
  ch_PU200_prompt->Draw("track_pt_>>h_PU200_prompt_track_pt_EB",               "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_PU200_prompt->Draw("track_pt_>>h_PU200_prompt_track_pt_EE",               "muon_prompt_==1 && muon_isBarrel_==0", "goff");
  ch_PU200_prompt->Draw("track_eta_>>h_PU200_prompt_track_eta",                "muon_prompt_==1",                      "goff");
  ch_PU200_prompt->Draw("track_eta_>>h_PU200_prompt_track_eta_EB",             "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_PU200_prompt->Draw("track_eta_>>h_PU200_prompt_track_eta_EE",             "muon_prompt_==1 && muon_isBarrel_==0", "goff");
  ch_PU200_prompt->Draw("track_phi_>>h_PU200_prompt_track_phi",                "muon_prompt_==1",                      "goff");
  ch_PU200_prompt->Draw("track_phi_>>h_PU200_prompt_track_phi_EB",             "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_PU200_prompt->Draw("track_phi_>>h_PU200_prompt_track_phi_EE",             "muon_prompt_==1 && muon_isBarrel_==0", "goff");
  // PU200 nonprompt
    // muon
  ch_PU200_nonprompt->Draw("muon_pt_>>h_PU200_nonprompt_muon_pt",            "muon_prompt_==0",                      "goff");
  ch_PU200_nonprompt->Draw("muon_pt_>>h_PU200_nonprompt_muon_pt_EB",         "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_PU200_nonprompt->Draw("muon_pt_>>h_PU200_nonprompt_muon_pt_EE",         "muon_prompt_==0 && muon_isBarrel_==0", "goff");
  ch_PU200_nonprompt->Draw("muon_eta_>>h_PU200_nonprompt_muon_eta",          "muon_prompt_==0",                      "goff");
  ch_PU200_nonprompt->Draw("muon_eta_>>h_PU200_nonprompt_muon_eta_EB",       "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_PU200_nonprompt->Draw("muon_eta_>>h_PU200_nonprompt_muon_eta_EE",       "muon_prompt_==0 && muon_isBarrel_==0", "goff");
  ch_PU200_nonprompt->Draw("muon_phi_>>h_PU200_nonprompt_muon_phi",          "muon_prompt_==0",                      "goff");
  ch_PU200_nonprompt->Draw("muon_phi_>>h_PU200_nonprompt_muon_phi_EB",       "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_PU200_nonprompt->Draw("muon_phi_>>h_PU200_nonprompt_muon_phi_EE",       "muon_prompt_==0 && muon_isBarrel_==0", "goff");
    // track
  ch_PU200_nonprompt->Draw("track_pt_>>h_PU200_nonprompt_track_pt",            "muon_prompt_==0",                      "goff");
  ch_PU200_nonprompt->Draw("track_pt_>>h_PU200_nonprompt_track_pt_EB",         "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_PU200_nonprompt->Draw("track_pt_>>h_PU200_nonprompt_track_pt_EE",         "muon_prompt_==0 && muon_isBarrel_==0", "goff");
  ch_PU200_nonprompt->Draw("track_eta_>>h_PU200_nonprompt_track_eta",          "muon_prompt_==0",                      "goff");
  ch_PU200_nonprompt->Draw("track_eta_>>h_PU200_nonprompt_track_eta_EB",       "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_PU200_nonprompt->Draw("track_eta_>>h_PU200_nonprompt_track_eta_EE",       "muon_prompt_==0 && muon_isBarrel_==0", "goff");
  ch_PU200_nonprompt->Draw("track_phi_>>h_PU200_nonprompt_track_phi",          "muon_prompt_==0",                      "goff");
  ch_PU200_nonprompt->Draw("track_phi_>>h_PU200_nonprompt_track_phi_EB",       "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_PU200_nonprompt->Draw("track_phi_>>h_PU200_nonprompt_track_phi_EE",       "muon_prompt_==0 && muon_isBarrel_==0", "goff");
  // noPU prompt
    // muon
  ch_noPU_prompt->Draw("muon_pt_>>h_noPU_prompt_muon_pt",                    "muon_prompt_==1",                      "goff");
  ch_noPU_prompt->Draw("muon_pt_>>h_noPU_prompt_muon_pt_EB",                 "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_noPU_prompt->Draw("muon_pt_>>h_noPU_prompt_muon_pt_EE",                 "muon_prompt_==1 && muon_isBarrel_==0", "goff");
  ch_noPU_prompt->Draw("muon_eta_>>h_noPU_prompt_muon_eta",                  "muon_prompt_==1",                      "goff");
  ch_noPU_prompt->Draw("muon_eta_>>h_noPU_prompt_muon_eta_EB",               "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_noPU_prompt->Draw("muon_eta_>>h_noPU_prompt_muon_eta_EE",               "muon_prompt_==1 && muon_isBarrel_==0", "goff");
  ch_noPU_prompt->Draw("muon_phi_>>h_noPU_prompt_muon_phi",                  "muon_prompt_==1",                      "goff");
  ch_noPU_prompt->Draw("muon_phi_>>h_noPU_prompt_muon_phi_EB",               "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_noPU_prompt->Draw("muon_phi_>>h_noPU_prompt_muon_phi_EE",               "muon_prompt_==1 && muon_isBarrel_==0", "goff");
    // track
  ch_noPU_prompt->Draw("track_pt_>>h_noPU_prompt_track_pt",                    "muon_prompt_==1",                      "goff");
  ch_noPU_prompt->Draw("track_pt_>>h_noPU_prompt_track_pt_EB",                 "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_noPU_prompt->Draw("track_pt_>>h_noPU_prompt_track_pt_EE",                 "muon_prompt_==1 && muon_isBarrel_==0", "goff");
  ch_noPU_prompt->Draw("track_eta_>>h_noPU_prompt_track_eta",                  "muon_prompt_==1",                      "goff");
  ch_noPU_prompt->Draw("track_eta_>>h_noPU_prompt_track_eta_EB",               "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_noPU_prompt->Draw("track_eta_>>h_noPU_prompt_track_eta_EE",               "muon_prompt_==1 && muon_isBarrel_==0", "goff");
  ch_noPU_prompt->Draw("track_phi_>>h_noPU_prompt_track_phi",                  "muon_prompt_==1",                      "goff");
  ch_noPU_prompt->Draw("track_phi_>>h_noPU_prompt_track_phi_EB",               "muon_prompt_==1 && muon_isBarrel_==1", "goff");
  ch_noPU_prompt->Draw("track_phi_>>h_noPU_prompt_track_phi_EE",               "muon_prompt_==1 && muon_isBarrel_==0", "goff");
  // noPU nonprompt
    // muon
  ch_noPU_nonprompt->Draw("muon_pt_>>h_noPU_nonprompt_muon_pt",              "muon_prompt_==0",                      "goff");
  ch_noPU_nonprompt->Draw("muon_pt_>>h_noPU_nonprompt_muon_pt_EB",           "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_noPU_nonprompt->Draw("muon_pt_>>h_noPU_nonprompt_muon_pt_EE",           "muon_prompt_==0 && muon_isBarrel_==0", "goff");
  ch_noPU_nonprompt->Draw("muon_eta_>>h_noPU_nonprompt_muon_eta",            "muon_prompt_==0",                      "goff");
  ch_noPU_nonprompt->Draw("muon_eta_>>h_noPU_nonprompt_muon_eta_EB",         "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_noPU_nonprompt->Draw("muon_eta_>>h_noPU_nonprompt_muon_eta_EE",         "muon_prompt_==0 && muon_isBarrel_==0", "goff");
  ch_noPU_nonprompt->Draw("muon_phi_>>h_noPU_nonprompt_muon_phi",            "muon_prompt_==0",                      "goff");
  ch_noPU_nonprompt->Draw("muon_phi_>>h_noPU_nonprompt_muon_phi_EB",         "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_noPU_nonprompt->Draw("muon_phi_>>h_noPU_nonprompt_muon_phi_EE",         "muon_prompt_==0 && muon_isBarrel_==0", "goff");
    // track
  ch_noPU_nonprompt->Draw("track_pt_>>h_noPU_nonprompt_track_pt",              "muon_prompt_==0",                      "goff");
  ch_noPU_nonprompt->Draw("track_pt_>>h_noPU_nonprompt_track_pt_EB",           "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_noPU_nonprompt->Draw("track_pt_>>h_noPU_nonprompt_track_pt_EE",           "muon_prompt_==0 && muon_isBarrel_==0", "goff");
  ch_noPU_nonprompt->Draw("track_eta_>>h_noPU_nonprompt_track_eta",            "muon_prompt_==0",                      "goff");
  ch_noPU_nonprompt->Draw("track_eta_>>h_noPU_nonprompt_track_eta_EB",         "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_noPU_nonprompt->Draw("track_eta_>>h_noPU_nonprompt_track_eta_EE",         "muon_prompt_==0 && muon_isBarrel_==0", "goff");
  ch_noPU_nonprompt->Draw("track_phi_>>h_noPU_nonprompt_track_phi",            "muon_prompt_==0",                      "goff");
  ch_noPU_nonprompt->Draw("track_phi_>>h_noPU_nonprompt_track_phi_EB",         "muon_prompt_==0 && muon_isBarrel_==1", "goff");
  ch_noPU_nonprompt->Draw("track_phi_>>h_noPU_nonprompt_track_phi_EE",         "muon_prompt_==0 && muon_isBarrel_==0", "goff");

  // Normalize
    // pt
  h_PU200_prompt_muon_pt->Scale(1./h_PU200_prompt_muon_pt->Integral(0,-1)); h_PU200_prompt_muon_pt_EB->Scale(1./h_PU200_prompt_muon_pt_EB->Integral(0,-1)); h_PU200_prompt_muon_pt_EE->Scale(1./h_PU200_prompt_muon_pt_EE->Integral(0,-1));
  h_PU200_nonprompt_muon_pt->Scale(1./h_PU200_nonprompt_muon_pt->Integral(0,-1)); h_PU200_nonprompt_muon_pt_EB->Scale(1./h_PU200_nonprompt_muon_pt_EB->Integral(0,-1)); h_PU200_nonprompt_muon_pt_EE->Scale(1./h_PU200_nonprompt_muon_pt_EE->Integral(0,-1));
  h_noPU_prompt_muon_pt->Scale(1./h_noPU_prompt_muon_pt->Integral(0,-1)); h_noPU_prompt_muon_pt_EB->Scale(1./h_noPU_prompt_muon_pt_EB->Integral(0,-1)); h_noPU_prompt_muon_pt_EE->Scale(1./h_noPU_prompt_muon_pt_EE->Integral(0,-1));
  h_noPU_nonprompt_muon_pt->Scale(1./h_noPU_nonprompt_muon_pt->Integral(0,-1)); h_noPU_nonprompt_muon_pt_EB->Scale(1./h_noPU_nonprompt_muon_pt_EB->Integral(0,-1)); h_noPU_nonprompt_muon_pt_EE->Scale(1./h_noPU_nonprompt_muon_pt_EE->Integral(0,-1));
  h_PU200_prompt_track_pt->Scale(1./h_PU200_prompt_track_pt->Integral(0,-1)); h_PU200_prompt_track_pt_EB->Scale(1./h_PU200_prompt_track_pt_EB->Integral(0,-1)); h_PU200_prompt_track_pt_EE->Scale(1./h_PU200_prompt_track_pt_EE->Integral(0,-1));
  h_PU200_nonprompt_track_pt->Scale(1./h_PU200_nonprompt_track_pt->Integral(0,-1)); h_PU200_nonprompt_track_pt_EB->Scale(1./h_PU200_nonprompt_track_pt_EB->Integral(0,-1)); h_PU200_nonprompt_track_pt_EE->Scale(1./h_PU200_nonprompt_track_pt_EE->Integral(0,-1));
  h_noPU_prompt_track_pt->Scale(1./h_noPU_prompt_track_pt->Integral(0,-1)); h_noPU_prompt_track_pt_EB->Scale(1./h_noPU_prompt_track_pt_EB->Integral(0,-1)); h_noPU_prompt_track_pt_EE->Scale(1./h_noPU_prompt_track_pt_EE->Integral(0,-1));
  h_noPU_nonprompt_track_pt->Scale(1./h_noPU_nonprompt_track_pt->Integral(0,-1)); h_noPU_nonprompt_track_pt_EB->Scale(1./h_noPU_nonprompt_track_pt_EB->Integral(0,-1)); h_noPU_nonprompt_track_pt_EE->Scale(1./h_noPU_nonprompt_track_pt_EE->Integral(0,-1));
    // eta
  h_PU200_prompt_muon_eta->Scale(1./h_PU200_prompt_muon_eta->Integral(0,-1)); h_PU200_prompt_muon_eta_EB->Scale(1./h_PU200_prompt_muon_eta_EB->Integral(0,-1)); h_PU200_prompt_muon_eta_EE->Scale(1./h_PU200_prompt_muon_eta_EE->Integral(0,-1));
  h_PU200_nonprompt_muon_eta->Scale(1./h_PU200_nonprompt_muon_eta->Integral(0,-1)); h_PU200_nonprompt_muon_eta_EB->Scale(1./h_PU200_nonprompt_muon_eta_EB->Integral(0,-1)); h_PU200_nonprompt_muon_eta_EE->Scale(1./h_PU200_nonprompt_muon_eta_EE->Integral(0,-1));
  h_noPU_prompt_muon_eta->Scale(1./h_noPU_prompt_muon_eta->Integral(0,-1)); h_noPU_prompt_muon_eta_EB->Scale(1./h_noPU_prompt_muon_eta_EB->Integral(0,-1)); h_noPU_prompt_muon_eta_EE->Scale(1./h_noPU_prompt_muon_eta_EE->Integral(0,-1));
  h_noPU_nonprompt_muon_eta->Scale(1./h_noPU_nonprompt_muon_eta->Integral(0,-1)); h_noPU_nonprompt_muon_eta_EB->Scale(1./h_noPU_nonprompt_muon_eta_EB->Integral(0,-1)); h_noPU_nonprompt_muon_eta_EE->Scale(1./h_noPU_nonprompt_muon_eta_EE->Integral(0,-1));
  h_PU200_prompt_track_eta->Scale(1./h_PU200_prompt_track_eta->Integral(0,-1)); h_PU200_prompt_track_eta_EB->Scale(1./h_PU200_prompt_track_eta_EB->Integral(0,-1)); h_PU200_prompt_track_eta_EE->Scale(1./h_PU200_prompt_track_eta_EE->Integral(0,-1));
  h_PU200_nonprompt_track_eta->Scale(1./h_PU200_nonprompt_track_eta->Integral(0,-1)); h_PU200_nonprompt_track_eta_EB->Scale(1./h_PU200_nonprompt_track_eta_EB->Integral(0,-1)); h_PU200_nonprompt_track_eta_EE->Scale(1./h_PU200_nonprompt_track_eta_EE->Integral(0,-1));
  h_noPU_prompt_track_eta->Scale(1./h_noPU_prompt_track_eta->Integral(0,-1)); h_noPU_prompt_track_eta_EB->Scale(1./h_noPU_prompt_track_eta_EB->Integral(0,-1)); h_noPU_prompt_track_eta_EE->Scale(1./h_noPU_prompt_track_eta_EE->Integral(0,-1));
  h_noPU_nonprompt_track_eta->Scale(1./h_noPU_nonprompt_track_eta->Integral(0,-1)); h_noPU_nonprompt_track_eta_EB->Scale(1./h_noPU_nonprompt_track_eta_EB->Integral(0,-1)); h_noPU_nonprompt_track_eta_EE->Scale(1./h_noPU_nonprompt_track_eta_EE->Integral(0,-1));
    // phi
  h_PU200_prompt_muon_phi->Scale(1./h_PU200_prompt_muon_phi->Integral(0,-1)); h_PU200_prompt_muon_phi_EB->Scale(1./h_PU200_prompt_muon_phi_EB->Integral(0,-1)); h_PU200_prompt_muon_phi_EE->Scale(1./h_PU200_prompt_muon_phi_EE->Integral(0,-1));
  h_PU200_nonprompt_muon_phi->Scale(1./h_PU200_nonprompt_muon_phi->Integral(0,-1)); h_PU200_nonprompt_muon_phi_EB->Scale(1./h_PU200_nonprompt_muon_phi_EB->Integral(0,-1)); h_PU200_nonprompt_muon_phi_EE->Scale(1./h_PU200_nonprompt_muon_phi_EE->Integral(0,-1));
  h_noPU_prompt_muon_phi->Scale(1./h_noPU_prompt_muon_phi->Integral(0,-1)); h_noPU_prompt_muon_phi_EB->Scale(1./h_noPU_prompt_muon_phi_EB->Integral(0,-1)); h_noPU_prompt_muon_phi_EE->Scale(1./h_noPU_prompt_muon_phi_EE->Integral(0,-1));
  h_noPU_nonprompt_muon_phi->Scale(1./h_noPU_nonprompt_muon_phi->Integral(0,-1)); h_noPU_nonprompt_muon_phi_EB->Scale(1./h_noPU_nonprompt_muon_phi_EB->Integral(0,-1)); h_noPU_nonprompt_muon_phi_EE->Scale(1./h_noPU_nonprompt_muon_phi_EE->Integral(0,-1));
  h_PU200_prompt_track_phi->Scale(1./h_PU200_prompt_track_phi->Integral(0,-1)); h_PU200_prompt_track_phi_EB->Scale(1./h_PU200_prompt_track_phi_EB->Integral(0,-1)); h_PU200_prompt_track_phi_EE->Scale(1./h_PU200_prompt_track_phi_EE->Integral(0,-1));
  h_PU200_nonprompt_track_phi->Scale(1./h_PU200_nonprompt_track_phi->Integral(0,-1)); h_PU200_nonprompt_track_phi_EB->Scale(1./h_PU200_nonprompt_track_phi_EB->Integral(0,-1)); h_PU200_nonprompt_track_phi_EE->Scale(1./h_PU200_nonprompt_track_phi_EE->Integral(0,-1));
  h_noPU_prompt_track_phi->Scale(1./h_noPU_prompt_track_phi->Integral(0,-1)); h_noPU_prompt_track_phi_EB->Scale(1./h_noPU_prompt_track_phi_EB->Integral(0,-1)); h_noPU_prompt_track_phi_EE->Scale(1./h_noPU_prompt_track_phi_EE->Integral(0,-1));
  h_noPU_nonprompt_track_phi->Scale(1./h_noPU_nonprompt_track_phi->Integral(0,-1)); h_noPU_nonprompt_track_phi_EB->Scale(1./h_noPU_nonprompt_track_phi_EB->Integral(0,-1)); h_noPU_nonprompt_track_phi_EE->Scale(1./h_noPU_nonprompt_track_phi_EE->Integral(0,-1));


  // Cosmetics
    // pt
  h_PU200_prompt_muon_pt->SetLineColor(kBlack); h_PU200_prompt_muon_pt_EB->SetLineColor(kBlack); h_PU200_prompt_muon_pt_EE->SetLineColor(kBlack);
  h_PU200_nonprompt_muon_pt->SetLineColor(kBlack); h_PU200_nonprompt_muon_pt_EB->SetLineColor(kBlack); h_PU200_nonprompt_muon_pt_EE->SetLineColor(kBlack);
  h_noPU_prompt_muon_pt->SetLineColor(kBlack); h_noPU_prompt_muon_pt_EB->SetLineColor(kBlack); h_noPU_prompt_muon_pt_EE->SetLineColor(kBlack);
  h_noPU_nonprompt_muon_pt->SetLineColor(kBlack); h_noPU_nonprompt_muon_pt_EB->SetLineColor(kBlack); h_noPU_nonprompt_muon_pt_EE->SetLineColor(kBlack);
  h_PU200_prompt_track_pt->SetLineColor(kBlack); h_PU200_prompt_track_pt_EB->SetLineColor(kBlack); h_PU200_prompt_track_pt_EE->SetLineColor(kBlack);
  h_PU200_nonprompt_track_pt->SetLineColor(kBlack); h_PU200_nonprompt_track_pt_EB->SetLineColor(kBlack); h_PU200_nonprompt_track_pt_EE->SetLineColor(kBlack);
  h_noPU_prompt_track_pt->SetLineColor(kBlack); h_noPU_prompt_track_pt_EB->SetLineColor(kBlack); h_noPU_prompt_track_pt_EE->SetLineColor(kBlack);
  h_noPU_nonprompt_track_pt->SetLineColor(kBlack); h_noPU_nonprompt_track_pt_EB->SetLineColor(kBlack); h_noPU_nonprompt_track_pt_EE->SetLineColor(kBlack);
    // eta
  h_PU200_prompt_muon_eta->SetLineColor(kBlack); h_PU200_prompt_muon_eta_EB->SetLineColor(kBlack); h_PU200_prompt_muon_eta_EE->SetLineColor(kBlack);
  h_PU200_nonprompt_muon_eta->SetLineColor(kBlack); h_PU200_nonprompt_muon_eta_EB->SetLineColor(kBlack); h_PU200_nonprompt_muon_eta_EE->SetLineColor(kBlack);
  h_noPU_prompt_muon_eta->SetLineColor(kBlack); h_noPU_prompt_muon_eta_EB->SetLineColor(kBlack); h_noPU_prompt_muon_eta_EE->SetLineColor(kBlack);
  h_noPU_nonprompt_muon_eta->SetLineColor(kBlack); h_noPU_nonprompt_muon_eta_EB->SetLineColor(kBlack); h_noPU_nonprompt_muon_eta_EE->SetLineColor(kBlack);
  h_PU200_prompt_track_eta->SetLineColor(kBlack); h_PU200_prompt_track_eta_EB->SetLineColor(kBlack); h_PU200_prompt_track_eta_EE->SetLineColor(kBlack);
  h_PU200_nonprompt_track_eta->SetLineColor(kBlack); h_PU200_nonprompt_track_eta_EB->SetLineColor(kBlack); h_PU200_nonprompt_track_eta_EE->SetLineColor(kBlack);
  h_noPU_prompt_track_eta->SetLineColor(kBlack); h_noPU_prompt_track_eta_EB->SetLineColor(kBlack); h_noPU_prompt_track_eta_EE->SetLineColor(kBlack);
  h_noPU_nonprompt_track_eta->SetLineColor(kBlack); h_noPU_nonprompt_track_eta_EB->SetLineColor(kBlack); h_noPU_nonprompt_track_eta_EE->SetLineColor(kBlack);
    // phi
  h_PU200_prompt_muon_phi->SetLineColor(kBlack); h_PU200_prompt_muon_phi_EB->SetLineColor(kBlack); h_PU200_prompt_muon_phi_EE->SetLineColor(kBlack);
  h_PU200_nonprompt_muon_phi->SetLineColor(kBlack); h_PU200_nonprompt_muon_phi_EB->SetLineColor(kBlack); h_PU200_nonprompt_muon_phi_EE->SetLineColor(kBlack);
  h_noPU_prompt_muon_phi->SetLineColor(kBlack); h_noPU_prompt_muon_phi_EB->SetLineColor(kBlack); h_noPU_prompt_muon_phi_EE->SetLineColor(kBlack);
  h_noPU_nonprompt_muon_phi->SetLineColor(kBlack); h_noPU_nonprompt_muon_phi_EB->SetLineColor(kBlack); h_noPU_nonprompt_muon_phi_EE->SetLineColor(kBlack);
  h_PU200_prompt_track_phi->SetLineColor(kBlack); h_PU200_prompt_track_phi_EB->SetLineColor(kBlack); h_PU200_prompt_track_phi_EE->SetLineColor(kBlack);
  h_PU200_nonprompt_track_phi->SetLineColor(kBlack); h_PU200_nonprompt_track_phi_EB->SetLineColor(kBlack); h_PU200_nonprompt_track_phi_EE->SetLineColor(kBlack);
  h_noPU_prompt_track_phi->SetLineColor(kBlack); h_noPU_prompt_track_phi_EB->SetLineColor(kBlack); h_noPU_prompt_track_phi_EE->SetLineColor(kBlack);
  h_noPU_nonprompt_track_phi->SetLineColor(kBlack); h_noPU_nonprompt_track_phi_EB->SetLineColor(kBlack); h_noPU_nonprompt_track_phi_EE->SetLineColor(kBlack);

  // PU200 prompt pt
    // muon
  TCanvas* c_PU200_prompt_muon_pt = new TCanvas("c_PU200_prompt_muon_pt", "c_PU200_prompt_muon_pt", 1500, 1500);
  c_PU200_prompt_muon_pt->cd();
  h_PU200_prompt_muon_pt->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_prompt_muon_pt->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_muon_pt->SetTitle("");
  h_PU200_prompt_muon_pt->Draw("hist");
  c_PU200_prompt_muon_pt->Print("plots/PU200_prompt_muon_pt.pdf");
  TCanvas* c_PU200_prompt_muon_pt_EB = new TCanvas("c_PU200_prompt_muon_pt_EB", "c_PU200_prompt_muon_pt_EB", 1500, 1500);
  c_PU200_prompt_muon_pt_EB->cd();
  h_PU200_prompt_muon_pt_EB->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_prompt_muon_pt_EB->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_muon_pt_EB->SetTitle("");
  h_PU200_prompt_muon_pt_EB->Draw("hist");
  c_PU200_prompt_muon_pt_EB->Print("plots/PU200_prompt_muon_pt_EB.pdf");
  TCanvas* c_PU200_prompt_muon_pt_EE = new TCanvas("c_PU200_prompt_muon_pt_EE", "c_PU200_prompt_muon_pt_EE", 1500, 1500);
  c_PU200_prompt_muon_pt_EE->cd();
  h_PU200_prompt_muon_pt_EE->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_prompt_muon_pt_EE->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_muon_pt_EE->SetTitle("");
  h_PU200_prompt_muon_pt_EE->Draw("hist");
  c_PU200_prompt_muon_pt_EE->Print("plots/PU200_prompt_muon_pt_EE.pdf");
    // track
  TCanvas* c_PU200_prompt_track_pt = new TCanvas("c_PU200_prompt_track_pt", "c_PU200_prompt_track_pt", 1500, 1500);
  c_PU200_prompt_track_pt->cd();
  h_PU200_prompt_track_pt->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_prompt_track_pt->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_track_pt->SetTitle("");
  h_PU200_prompt_track_pt->Draw("hist");
  c_PU200_prompt_track_pt->Print("plots/PU200_prompt_track_pt.pdf");
  TCanvas* c_PU200_prompt_track_pt_EB = new TCanvas("c_PU200_prompt_track_pt_EB", "c_PU200_prompt_track_pt_EB", 1500, 1500);
  c_PU200_prompt_track_pt_EB->cd();
  h_PU200_prompt_track_pt_EB->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_prompt_track_pt_EB->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_track_pt_EB->SetTitle("");
  h_PU200_prompt_track_pt_EB->Draw("hist");
  c_PU200_prompt_track_pt_EB->Print("plots/PU200_prompt_track_pt_EB.pdf");
  TCanvas* c_PU200_prompt_track_pt_EE = new TCanvas("c_PU200_prompt_track_pt_EE", "c_PU200_prompt_track_pt_EE", 1500, 1500);
  c_PU200_prompt_track_pt_EE->cd();
  h_PU200_prompt_track_pt_EE->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_prompt_track_pt_EE->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_track_pt_EE->SetTitle("");
  h_PU200_prompt_track_pt_EE->Draw("hist");
  c_PU200_prompt_track_pt_EE->Print("plots/PU200_prompt_track_pt_EE.pdf");

  // PU200 prompt eta
    // muon
  TCanvas* c_PU200_prompt_muon_eta = new TCanvas("c_PU200_prompt_muon_eta", "c_PU200_prompt_muon_eta", 1500, 1500);
  c_PU200_prompt_muon_eta->cd();
  h_PU200_prompt_muon_eta->GetXaxis()->SetTitle("#eta (GeV)");
  h_PU200_prompt_muon_eta->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_muon_eta->SetTitle("");
  h_PU200_prompt_muon_eta->Draw("hist");
  c_PU200_prompt_muon_eta->Print("plots/PU200_prompt_muon_eta.pdf");
  TCanvas* c_PU200_prompt_muon_eta_EB = new TCanvas("c_PU200_prompt_muon_eta_EB", "c_PU200_prompt_muon_eta_EB", 1500, 1500);
  c_PU200_prompt_muon_eta_EB->cd();
  h_PU200_prompt_muon_eta_EB->GetXaxis()->SetTitle("#eta (GeV)");
  h_PU200_prompt_muon_eta_EB->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_muon_eta_EB->SetTitle("");
  h_PU200_prompt_muon_eta_EB->Draw("hist");
  c_PU200_prompt_muon_eta_EB->Print("plots/PU200_prompt_muon_eta_EB.pdf");
  TCanvas* c_PU200_prompt_muon_eta_EE = new TCanvas("c_PU200_prompt_muon_eta_EE", "c_PU200_prompt_muon_eta_EE", 1500, 1500);
  c_PU200_prompt_muon_eta_EE->cd();
  h_PU200_prompt_muon_eta_EE->GetXaxis()->SetTitle("#eta (GeV)");
  h_PU200_prompt_muon_eta_EE->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_muon_eta_EE->SetTitle("");
  h_PU200_prompt_muon_eta_EE->Draw("hist");
  c_PU200_prompt_muon_eta_EE->Print("plots/PU200_prompt_muon_eta_EE.pdf");
    // track
  TCanvas* c_PU200_prompt_track_eta = new TCanvas("c_PU200_prompt_track_eta", "c_PU200_prompt_track_eta", 1500, 1500);
  c_PU200_prompt_track_eta->cd();
  h_PU200_prompt_track_eta->GetXaxis()->SetTitle("#eta (GeV)");
  h_PU200_prompt_track_eta->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_track_eta->SetTitle("");
  h_PU200_prompt_track_eta->Draw("hist");
  c_PU200_prompt_track_eta->Print("plots/PU200_prompt_track_eta.pdf");
  TCanvas* c_PU200_prompt_track_eta_EB = new TCanvas("c_PU200_prompt_track_eta_EB", "c_PU200_prompt_track_eta_EB", 1500, 1500);
  c_PU200_prompt_track_eta_EB->cd();
  h_PU200_prompt_track_eta_EB->GetXaxis()->SetTitle("#eta (GeV)");
  h_PU200_prompt_track_eta_EB->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_track_eta_EB->SetTitle("");
  h_PU200_prompt_track_eta_EB->Draw("hist");
  c_PU200_prompt_track_eta_EB->Print("plots/PU200_prompt_track_eta_EB.pdf");
  TCanvas* c_PU200_prompt_track_eta_EE = new TCanvas("c_PU200_prompt_track_eta_EE", "c_PU200_prompt_track_eta_EE", 1500, 1500);
  c_PU200_prompt_track_eta_EE->cd();
  h_PU200_prompt_track_eta_EE->GetXaxis()->SetTitle("#eta (GeV)");
  h_PU200_prompt_track_eta_EE->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_track_eta_EE->SetTitle("");
  h_PU200_prompt_track_eta_EE->Draw("hist");
  c_PU200_prompt_track_eta_EE->Print("plots/PU200_prompt_track_eta_EE.pdf");

  // PU200 prompt phi
    // muon
  TCanvas* c_PU200_prompt_muon_phi = new TCanvas("c_PU200_prompt_muon_phi", "c_PU200_prompt_muon_phi", 1500, 1500);
  c_PU200_prompt_muon_phi->cd();
  h_PU200_prompt_muon_phi->GetXaxis()->SetTitle("#phi (GeV)");
  h_PU200_prompt_muon_phi->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_muon_phi->SetTitle("");
  h_PU200_prompt_muon_phi->Draw("hist");
  c_PU200_prompt_muon_phi->Print("plots/PU200_prompt_muon_phi.pdf");
  TCanvas* c_PU200_prompt_muon_phi_EB = new TCanvas("c_PU200_prompt_muon_phi_EB", "c_PU200_prompt_muon_phi_EB", 1500, 1500);
  c_PU200_prompt_muon_phi_EB->cd();
  h_PU200_prompt_muon_phi_EB->GetXaxis()->SetTitle("#phi (GeV)");
  h_PU200_prompt_muon_phi_EB->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_muon_phi_EB->SetTitle("");
  h_PU200_prompt_muon_phi_EB->Draw("hist");
  c_PU200_prompt_muon_phi_EB->Print("plots/PU200_prompt_muon_phi_EB.pdf");
  TCanvas* c_PU200_prompt_muon_phi_EE = new TCanvas("c_PU200_prompt_muon_phi_EE", "c_PU200_prompt_muon_phi_EE", 1500, 1500);
  c_PU200_prompt_muon_phi_EE->cd();
  h_PU200_prompt_muon_phi_EE->GetXaxis()->SetTitle("#phi (GeV)");
  h_PU200_prompt_muon_phi_EE->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_muon_phi_EE->SetTitle("");
  h_PU200_prompt_muon_phi_EE->Draw("hist");
  c_PU200_prompt_muon_phi_EE->Print("plots/PU200_prompt_muon_phi_EE.pdf");
    // track
  TCanvas* c_PU200_prompt_track_phi = new TCanvas("c_PU200_prompt_track_phi", "c_PU200_prompt_track_phi", 1500, 1500);
  c_PU200_prompt_track_phi->cd();
  h_PU200_prompt_track_phi->GetXaxis()->SetTitle("#phi (GeV)");
  h_PU200_prompt_track_phi->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_track_phi->SetTitle("");
  h_PU200_prompt_track_phi->Draw("hist");
  c_PU200_prompt_track_phi->Print("plots/PU200_prompt_track_phi.pdf");
  TCanvas* c_PU200_prompt_track_phi_EB = new TCanvas("c_PU200_prompt_track_phi_EB", "c_PU200_prompt_track_phi_EB", 1500, 1500);
  c_PU200_prompt_track_phi_EB->cd();
  h_PU200_prompt_track_phi_EB->GetXaxis()->SetTitle("#phi (GeV)");
  h_PU200_prompt_track_phi_EB->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_track_phi_EB->SetTitle("");
  h_PU200_prompt_track_phi_EB->Draw("hist");
  c_PU200_prompt_track_phi_EB->Print("plots/PU200_prompt_track_phi_EB.pdf");
  TCanvas* c_PU200_prompt_track_phi_EE = new TCanvas("c_PU200_prompt_track_phi_EE", "c_PU200_prompt_track_phi_EE", 1500, 1500);
  c_PU200_prompt_track_phi_EE->cd();
  h_PU200_prompt_track_phi_EE->GetXaxis()->SetTitle("#phi (GeV)");
  h_PU200_prompt_track_phi_EE->GetYaxis()->SetTitle("% Counts");
  h_PU200_prompt_track_phi_EE->SetTitle("");
  h_PU200_prompt_track_phi_EE->Draw("hist");
  c_PU200_prompt_track_phi_EE->Print("plots/PU200_prompt_track_phi_EE.pdf");

  // PU200 nonprompt pt
    // muon
  TCanvas* c_PU200_nonprompt_muon_pt = new TCanvas("c_PU200_nonprompt_muon_pt", "c_PU200_nonprompt_muon_pt", 1500, 1500);
  c_PU200_nonprompt_muon_pt->cd();
  h_PU200_nonprompt_muon_pt->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_nonprompt_muon_pt->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_muon_pt->SetTitle("");
  h_PU200_nonprompt_muon_pt->Draw("hist");
  c_PU200_nonprompt_muon_pt->Print("plots/PU200_nonprompt_muon_pt.pdf");
  TCanvas* c_PU200_nonprompt_muon_pt_EB = new TCanvas("c_PU200_nonprompt_muon_pt_EB", "c_PU200_nonprompt_muon_pt_EB", 1500, 1500);
  c_PU200_nonprompt_muon_pt_EB->cd();
  h_PU200_nonprompt_muon_pt_EB->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_nonprompt_muon_pt_EB->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_muon_pt_EB->SetTitle("");
  h_PU200_nonprompt_muon_pt_EB->Draw("hist");
  c_PU200_nonprompt_muon_pt_EB->Print("plots/PU200_nonprompt_muon_pt_EB.pdf");
  TCanvas* c_PU200_nonprompt_muon_pt_EE = new TCanvas("c_PU200_nonprompt_muon_pt_EE", "c_PU200_nonprompt_muon_pt_EE", 1500, 1500);
  c_PU200_nonprompt_muon_pt_EE->cd();
  h_PU200_nonprompt_muon_pt_EE->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_nonprompt_muon_pt_EE->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_muon_pt_EE->SetTitle("");
  h_PU200_nonprompt_muon_pt_EE->Draw("hist");
  c_PU200_nonprompt_muon_pt_EE->Print("plots/PU200_nonprompt_muon_pt_EE.pdf");
    // track
  TCanvas* c_PU200_nonprompt_track_pt = new TCanvas("c_PU200_nonprompt_track_pt", "c_PU200_nonprompt_track_pt", 1500, 1500);
  c_PU200_nonprompt_track_pt->cd();
  h_PU200_nonprompt_track_pt->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_nonprompt_track_pt->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_track_pt->SetTitle("");
  h_PU200_nonprompt_track_pt->Draw("hist");
  c_PU200_nonprompt_track_pt->Print("plots/PU200_nonprompt_track_pt.pdf");
  TCanvas* c_PU200_nonprompt_track_pt_EB = new TCanvas("c_PU200_nonprompt_track_pt_EB", "c_PU200_nonprompt_track_pt_EB", 1500, 1500);
  c_PU200_nonprompt_track_pt_EB->cd();
  h_PU200_nonprompt_track_pt_EB->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_nonprompt_track_pt_EB->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_track_pt_EB->SetTitle("");
  h_PU200_nonprompt_track_pt_EB->Draw("hist");
  c_PU200_nonprompt_track_pt_EB->Print("plots/PU200_nonprompt_track_pt_EB.pdf");
  TCanvas* c_PU200_nonprompt_track_pt_EE = new TCanvas("c_PU200_nonprompt_track_pt_EE", "c_PU200_nonprompt_track_pt_EE", 1500, 1500);
  c_PU200_nonprompt_track_pt_EE->cd();
  h_PU200_nonprompt_track_pt_EE->GetXaxis()->SetTitle("pT (GeV)");
  h_PU200_nonprompt_track_pt_EE->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_track_pt_EE->SetTitle("");
  h_PU200_nonprompt_track_pt_EE->Draw("hist");
  c_PU200_nonprompt_track_pt_EE->Print("plots/PU200_nonprompt_track_pt_EE.pdf");

  // PU200 nonprompt eta
    // muon
  TCanvas* c_PU200_nonprompt_muon_eta = new TCanvas("c_PU200_nonprompt_muon_eta", "c_PU200_nonprompt_muon_eta", 1500, 1500);
  c_PU200_nonprompt_muon_eta->cd();
  h_PU200_nonprompt_muon_eta->GetXaxis()->SetTitle("#eta (GeV)");
  h_PU200_nonprompt_muon_eta->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_muon_eta->SetTitle("");
  h_PU200_nonprompt_muon_eta->Draw("hist");
  c_PU200_nonprompt_muon_eta->Print("plots/PU200_nonprompt_muon_eta.pdf");
  TCanvas* c_PU200_nonprompt_muon_eta_EB = new TCanvas("c_PU200_nonprompt_muon_eta_EB", "c_PU200_nonprompt_muon_eta_EB", 1500, 1500);
  c_PU200_nonprompt_muon_eta_EB->cd();
  h_PU200_nonprompt_muon_eta_EB->GetXaxis()->SetTitle("#eta (GeV)");
  h_PU200_nonprompt_muon_eta_EB->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_muon_eta_EB->SetTitle("");
  h_PU200_nonprompt_muon_eta_EB->Draw("hist");
  c_PU200_nonprompt_muon_eta_EB->Print("plots/PU200_nonprompt_muon_eta_EB.pdf");
  TCanvas* c_PU200_nonprompt_muon_eta_EE = new TCanvas("c_PU200_nonprompt_muon_eta_EE", "c_PU200_nonprompt_muon_eta_EE", 1500, 1500);
  c_PU200_nonprompt_muon_eta_EE->cd();
  h_PU200_nonprompt_muon_eta_EE->GetXaxis()->SetTitle("#eta (GeV)");
  h_PU200_nonprompt_muon_eta_EE->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_muon_eta_EE->SetTitle("");
  h_PU200_nonprompt_muon_eta_EE->Draw("hist");
  c_PU200_nonprompt_muon_eta_EE->Print("plots/PU200_nonprompt_muon_eta_EE.pdf");
    // track
  TCanvas* c_PU200_nonprompt_track_eta = new TCanvas("c_PU200_nonprompt_track_eta", "c_PU200_nonprompt_track_eta", 1500, 1500);
  c_PU200_nonprompt_track_eta->cd();
  h_PU200_nonprompt_track_eta->GetXaxis()->SetTitle("#eta (GeV)");
  h_PU200_nonprompt_track_eta->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_track_eta->SetTitle("");
  h_PU200_nonprompt_track_eta->Draw("hist");
  c_PU200_nonprompt_track_eta->Print("plots/PU200_nonprompt_track_eta.pdf");
  TCanvas* c_PU200_nonprompt_track_eta_EB = new TCanvas("c_PU200_nonprompt_track_eta_EB", "c_PU200_nonprompt_track_eta_EB", 1500, 1500);
  c_PU200_nonprompt_track_eta_EB->cd();
  h_PU200_nonprompt_track_eta_EB->GetXaxis()->SetTitle("#eta (GeV)");
  h_PU200_nonprompt_track_eta_EB->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_track_eta_EB->SetTitle("");
  h_PU200_nonprompt_track_eta_EB->Draw("hist");
  c_PU200_nonprompt_track_eta_EB->Print("plots/PU200_nonprompt_track_eta_EB.pdf");
  TCanvas* c_PU200_nonprompt_track_eta_EE = new TCanvas("c_PU200_nonprompt_track_eta_EE", "c_PU200_nonprompt_track_eta_EE", 1500, 1500);
  c_PU200_nonprompt_track_eta_EE->cd();
  h_PU200_nonprompt_track_eta_EE->GetXaxis()->SetTitle("#eta (GeV)");
  h_PU200_nonprompt_track_eta_EE->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_track_eta_EE->SetTitle("");
  h_PU200_nonprompt_track_eta_EE->Draw("hist");
  c_PU200_nonprompt_track_eta_EE->Print("plots/PU200_nonprompt_track_eta_EE.pdf");

  // PU200 nonprompt phi
    // muon
  TCanvas* c_PU200_nonprompt_muon_phi = new TCanvas("c_PU200_nonprompt_muon_phi", "c_PU200_nonprompt_muon_phi", 1500, 1500);
  c_PU200_nonprompt_muon_phi->cd();
  h_PU200_nonprompt_muon_phi->GetXaxis()->SetTitle("#phi (GeV)");
  h_PU200_nonprompt_muon_phi->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_muon_phi->SetTitle("");
  h_PU200_nonprompt_muon_phi->Draw("hist");
  c_PU200_nonprompt_muon_phi->Print("plots/PU200_nonprompt_muon_phi.pdf");
  TCanvas* c_PU200_nonprompt_muon_phi_EB = new TCanvas("c_PU200_nonprompt_muon_phi_EB", "c_PU200_nonprompt_muon_phi_EB", 1500, 1500);
  c_PU200_nonprompt_muon_phi_EB->cd();
  h_PU200_nonprompt_muon_phi_EB->GetXaxis()->SetTitle("#phi (GeV)");
  h_PU200_nonprompt_muon_phi_EB->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_muon_phi_EB->SetTitle("");
  h_PU200_nonprompt_muon_phi_EB->Draw("hist");
  c_PU200_nonprompt_muon_phi_EB->Print("plots/PU200_nonprompt_muon_phi_EB.pdf");
  TCanvas* c_PU200_nonprompt_muon_phi_EE = new TCanvas("c_PU200_nonprompt_muon_phi_EE", "c_PU200_nonprompt_muon_phi_EE", 1500, 1500);
  c_PU200_nonprompt_muon_phi_EE->cd();
  h_PU200_nonprompt_muon_phi_EE->GetXaxis()->SetTitle("#phi (GeV)");
  h_PU200_nonprompt_muon_phi_EE->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_muon_phi_EE->SetTitle("");
  h_PU200_nonprompt_muon_phi_EE->Draw("hist");
  c_PU200_nonprompt_muon_phi_EE->Print("plots/PU200_nonprompt_muon_phi_EE.pdf");
    // track
  TCanvas* c_PU200_nonprompt_track_phi = new TCanvas("c_PU200_nonprompt_track_phi", "c_PU200_nonprompt_track_phi", 1500, 1500);
  c_PU200_nonprompt_track_phi->cd();
  h_PU200_nonprompt_track_phi->GetXaxis()->SetTitle("#phi (GeV)");
  h_PU200_nonprompt_track_phi->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_track_phi->SetTitle("");
  h_PU200_nonprompt_track_phi->Draw("hist");
  c_PU200_nonprompt_track_phi->Print("plots/PU200_nonprompt_track_phi.pdf");
  TCanvas* c_PU200_nonprompt_track_phi_EB = new TCanvas("c_PU200_nonprompt_track_phi_EB", "c_PU200_nonprompt_track_phi_EB", 1500, 1500);
  c_PU200_nonprompt_track_phi_EB->cd();
  h_PU200_nonprompt_track_phi_EB->GetXaxis()->SetTitle("#phi (GeV)");
  h_PU200_nonprompt_track_phi_EB->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_track_phi_EB->SetTitle("");
  h_PU200_nonprompt_track_phi_EB->Draw("hist");
  c_PU200_nonprompt_track_phi_EB->Print("plots/PU200_nonprompt_track_phi_EB.pdf");
  TCanvas* c_PU200_nonprompt_track_phi_EE = new TCanvas("c_PU200_nonprompt_track_phi_EE", "c_PU200_nonprompt_track_phi_EE", 1500, 1500);
  c_PU200_nonprompt_track_phi_EE->cd();
  h_PU200_nonprompt_track_phi_EE->GetXaxis()->SetTitle("#phi (GeV)");
  h_PU200_nonprompt_track_phi_EE->GetYaxis()->SetTitle("% Counts");
  h_PU200_nonprompt_track_phi_EE->SetTitle("");
  h_PU200_nonprompt_track_phi_EE->Draw("hist");
  c_PU200_nonprompt_track_phi_EE->Print("plots/PU200_nonprompt_track_phi_EE.pdf");

  // noPU prompt pt
    // muon
  TCanvas* c_noPU_prompt_muon_pt = new TCanvas("c_noPU_prompt_muon_pt", "c_noPU_prompt_muon_pt", 1500, 1500);
  c_noPU_prompt_muon_pt->cd();
  h_noPU_prompt_muon_pt->GetXaxis()->SetTitle("pT (GeV)");
  h_noPU_prompt_muon_pt->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_muon_pt->SetTitle("");
  h_noPU_prompt_muon_pt->Draw("hist");
  c_noPU_prompt_muon_pt->Print("plots/noPU_prompt_muon_pt.pdf");
  TCanvas* c_noPU_prompt_muon_pt_EB = new TCanvas("c_noPU_prompt_muon_pt_EB", "c_noPU_prompt_muon_pt_EB", 1500, 1500);
  c_noPU_prompt_muon_pt_EB->cd();
  h_noPU_prompt_muon_pt_EB->GetXaxis()->SetTitle("pT (GeV)");
  h_noPU_prompt_muon_pt_EB->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_muon_pt_EB->SetTitle("");
  h_noPU_prompt_muon_pt_EB->Draw("hist");
  c_noPU_prompt_muon_pt_EB->Print("plots/noPU_prompt_muon_pt_EB.pdf");
  TCanvas* c_noPU_prompt_muon_pt_EE = new TCanvas("c_noPU_prompt_muon_pt_EE", "c_noPU_prompt_muon_pt_EE", 1500, 1500);
  c_noPU_prompt_muon_pt_EE->cd();
  h_noPU_prompt_muon_pt_EE->GetXaxis()->SetTitle("pT (GeV)");
  h_noPU_prompt_muon_pt_EE->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_muon_pt_EE->SetTitle("");
  h_noPU_prompt_muon_pt_EE->Draw("hist");
  c_noPU_prompt_muon_pt_EE->Print("plots/noPU_prompt_muon_pt_EE.pdf");
    // track
  TCanvas* c_noPU_prompt_track_pt = new TCanvas("c_noPU_prompt_track_pt", "c_noPU_prompt_track_pt", 1500, 1500);
  c_noPU_prompt_track_pt->cd();
  h_noPU_prompt_track_pt->GetXaxis()->SetTitle("pT (GeV)");
  h_noPU_prompt_track_pt->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_track_pt->SetTitle("");
  h_noPU_prompt_track_pt->Draw("hist");
  c_noPU_prompt_track_pt->Print("plots/noPU_prompt_track_pt.pdf");
  TCanvas* c_noPU_prompt_track_pt_EB = new TCanvas("c_noPU_prompt_track_pt_EB", "c_noPU_prompt_track_pt_EB", 1500, 1500);
  c_noPU_prompt_track_pt_EB->cd();
  h_noPU_prompt_track_pt_EB->GetXaxis()->SetTitle("pT (GeV)");
  h_noPU_prompt_track_pt_EB->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_track_pt_EB->SetTitle("");
  h_noPU_prompt_track_pt_EB->Draw("hist");
  c_noPU_prompt_track_pt_EB->Print("plots/noPU_prompt_track_pt_EB.pdf");
  TCanvas* c_noPU_prompt_track_pt_EE = new TCanvas("c_noPU_prompt_track_pt_EE", "c_noPU_prompt_track_pt_EE", 1500, 1500);
  c_noPU_prompt_track_pt_EE->cd();
  h_noPU_prompt_track_pt_EE->GetXaxis()->SetTitle("pT (GeV)");
  h_noPU_prompt_track_pt_EE->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_track_pt_EE->SetTitle("");
  h_noPU_prompt_track_pt_EE->Draw("hist");
  c_noPU_prompt_track_pt_EE->Print("plots/noPU_prompt_track_pt_EE.pdf");

  // noPU prompt eta
    // muon
  TCanvas* c_noPU_prompt_muon_eta = new TCanvas("c_noPU_prompt_muon_eta", "c_noPU_prompt_muon_eta", 1500, 1500);
  c_noPU_prompt_muon_eta->cd();
  h_noPU_prompt_muon_eta->GetXaxis()->SetTitle("#eta (GeV)");
  h_noPU_prompt_muon_eta->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_muon_eta->SetTitle("");
  h_noPU_prompt_muon_eta->Draw("hist");
  c_noPU_prompt_muon_eta->Print("plots/noPU_prompt_muon_eta.pdf");
  TCanvas* c_noPU_prompt_muon_eta_EB = new TCanvas("c_noPU_prompt_muon_eta_EB", "c_noPU_prompt_muon_eta_EB", 1500, 1500);
  c_noPU_prompt_muon_eta_EB->cd();
  h_noPU_prompt_muon_eta_EB->GetXaxis()->SetTitle("#eta (GeV)");
  h_noPU_prompt_muon_eta_EB->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_muon_eta_EB->SetTitle("");
  h_noPU_prompt_muon_eta_EB->Draw("hist");
  c_noPU_prompt_muon_eta_EB->Print("plots/noPU_prompt_muon_eta_EB.pdf");
  TCanvas* c_noPU_prompt_muon_eta_EE = new TCanvas("c_noPU_prompt_muon_eta_EE", "c_noPU_prompt_muon_eta_EE", 1500, 1500);
  c_noPU_prompt_muon_eta_EE->cd();
  h_noPU_prompt_muon_eta_EE->GetXaxis()->SetTitle("#eta (GeV)");
  h_noPU_prompt_muon_eta_EE->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_muon_eta_EE->SetTitle("");
  h_noPU_prompt_muon_eta_EE->Draw("hist");
  c_noPU_prompt_muon_eta_EE->Print("plots/noPU_prompt_muon_eta_EE.pdf");
    // track
  TCanvas* c_noPU_prompt_track_eta = new TCanvas("c_noPU_prompt_track_eta", "c_noPU_prompt_track_eta", 1500, 1500);
  c_noPU_prompt_track_eta->cd();
  h_noPU_prompt_track_eta->GetXaxis()->SetTitle("#eta (GeV)");
  h_noPU_prompt_track_eta->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_track_eta->SetTitle("");
  h_noPU_prompt_track_eta->Draw("hist");
  c_noPU_prompt_track_eta->Print("plots/noPU_prompt_track_eta.pdf");
  TCanvas* c_noPU_prompt_track_eta_EB = new TCanvas("c_noPU_prompt_track_eta_EB", "c_noPU_prompt_track_eta_EB", 1500, 1500);
  c_noPU_prompt_track_eta_EB->cd();
  h_noPU_prompt_track_eta_EB->GetXaxis()->SetTitle("#eta (GeV)");
  h_noPU_prompt_track_eta_EB->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_track_eta_EB->SetTitle("");
  h_noPU_prompt_track_eta_EB->Draw("hist");
  c_noPU_prompt_track_eta_EB->Print("plots/noPU_prompt_track_eta_EB.pdf");
  TCanvas* c_noPU_prompt_track_eta_EE = new TCanvas("c_noPU_prompt_track_eta_EE", "c_noPU_prompt_track_eta_EE", 1500, 1500);
  c_noPU_prompt_track_eta_EE->cd();
  h_noPU_prompt_track_eta_EE->GetXaxis()->SetTitle("#eta (GeV)");
  h_noPU_prompt_track_eta_EE->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_track_eta_EE->SetTitle("");
  h_noPU_prompt_track_eta_EE->Draw("hist");
  c_noPU_prompt_track_eta_EE->Print("plots/noPU_prompt_track_eta_EE.pdf");

  // noPU prompt phi
    // muon
  TCanvas* c_noPU_prompt_muon_phi = new TCanvas("c_noPU_prompt_muon_phi", "c_noPU_prompt_muon_phi", 1500, 1500);
  c_noPU_prompt_muon_phi->cd();
  h_noPU_prompt_muon_phi->GetXaxis()->SetTitle("#phi (GeV)");
  h_noPU_prompt_muon_phi->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_muon_phi->SetTitle("");
  h_noPU_prompt_muon_phi->Draw("hist");
  c_noPU_prompt_muon_phi->Print("plots/noPU_prompt_muon_phi.pdf");
  TCanvas* c_noPU_prompt_muon_phi_EB = new TCanvas("c_noPU_prompt_muon_phi_EB", "c_noPU_prompt_muon_phi_EB", 1500, 1500);
  c_noPU_prompt_muon_phi_EB->cd();
  h_noPU_prompt_muon_phi_EB->GetXaxis()->SetTitle("#phi (GeV)");
  h_noPU_prompt_muon_phi_EB->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_muon_phi_EB->SetTitle("");
  h_noPU_prompt_muon_phi_EB->Draw("hist");
  c_noPU_prompt_muon_phi_EB->Print("plots/noPU_prompt_muon_phi_EB.pdf");
  TCanvas* c_noPU_prompt_muon_phi_EE = new TCanvas("c_noPU_prompt_muon_phi_EE", "c_noPU_prompt_muon_phi_EE", 1500, 1500);
  c_noPU_prompt_muon_phi_EE->cd();
  h_noPU_prompt_muon_phi_EE->GetXaxis()->SetTitle("#phi (GeV)");
  h_noPU_prompt_muon_phi_EE->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_muon_phi_EE->SetTitle("");
  h_noPU_prompt_muon_phi_EE->Draw("hist");
  c_noPU_prompt_muon_phi_EE->Print("plots/noPU_prompt_muon_phi_EE.pdf");
    // track
  TCanvas* c_noPU_prompt_track_phi = new TCanvas("c_noPU_prompt_track_phi", "c_noPU_prompt_track_phi", 1500, 1500);
  c_noPU_prompt_track_phi->cd();
  h_noPU_prompt_track_phi->GetXaxis()->SetTitle("#phi (GeV)");
  h_noPU_prompt_track_phi->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_track_phi->SetTitle("");
  h_noPU_prompt_track_phi->Draw("hist");
  c_noPU_prompt_track_phi->Print("plots/noPU_prompt_track_phi.pdf");
  TCanvas* c_noPU_prompt_track_phi_EB = new TCanvas("c_noPU_prompt_track_phi_EB", "c_noPU_prompt_track_phi_EB", 1500, 1500);
  c_noPU_prompt_track_phi_EB->cd();
  h_noPU_prompt_track_phi_EB->GetXaxis()->SetTitle("#phi (GeV)");
  h_noPU_prompt_track_phi_EB->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_track_phi_EB->SetTitle("");
  h_noPU_prompt_track_phi_EB->Draw("hist");
  c_noPU_prompt_track_phi_EB->Print("plots/noPU_prompt_track_phi_EB.pdf");
  TCanvas* c_noPU_prompt_track_phi_EE = new TCanvas("c_noPU_prompt_track_phi_EE", "c_noPU_prompt_track_phi_EE", 1500, 1500);
  c_noPU_prompt_track_phi_EE->cd();
  h_noPU_prompt_track_phi_EE->GetXaxis()->SetTitle("#phi (GeV)");
  h_noPU_prompt_track_phi_EE->GetYaxis()->SetTitle("% Counts");
  h_noPU_prompt_track_phi_EE->SetTitle("");
  h_noPU_prompt_track_phi_EE->Draw("hist");
  c_noPU_prompt_track_phi_EE->Print("plots/noPU_prompt_track_phi_EE.pdf");

  // noPU nonprompt pt
    // muon
  TCanvas* c_noPU_nonprompt_muon_pt = new TCanvas("c_noPU_nonprompt_muon_pt", "c_noPU_nonprompt_muon_pt", 1500, 1500);
  c_noPU_nonprompt_muon_pt->cd();
  h_noPU_nonprompt_muon_pt->GetXaxis()->SetTitle("pT (GeV)");
  h_noPU_nonprompt_muon_pt->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_muon_pt->SetTitle("");
  h_noPU_nonprompt_muon_pt->Draw("hist");
  c_noPU_nonprompt_muon_pt->Print("plots/noPU_nonprompt_muon_pt.pdf");
  TCanvas* c_noPU_nonprompt_muon_pt_EB = new TCanvas("c_noPU_nonprompt_muon_pt_EB", "c_noPU_nonprompt_muon_pt_EB", 1500, 1500);
  c_noPU_nonprompt_muon_pt_EB->cd();
  h_noPU_nonprompt_muon_pt_EB->GetXaxis()->SetTitle("pT (GeV)");
  h_noPU_nonprompt_muon_pt_EB->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_muon_pt_EB->SetTitle("");
  h_noPU_nonprompt_muon_pt_EB->Draw("hist");
  c_noPU_nonprompt_muon_pt_EB->Print("plots/noPU_nonprompt_muon_pt_EB.pdf");
  TCanvas* c_noPU_nonprompt_muon_pt_EE = new TCanvas("c_noPU_nonprompt_muon_pt_EE", "c_noPU_nonprompt_muon_pt_EE", 1500, 1500);
  c_noPU_nonprompt_muon_pt_EE->cd();
  h_noPU_nonprompt_muon_pt_EE->GetXaxis()->SetTitle("pT (GeV)");
  h_noPU_nonprompt_muon_pt_EE->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_muon_pt_EE->SetTitle("");
  h_noPU_nonprompt_muon_pt_EE->Draw("hist");
  c_noPU_nonprompt_muon_pt_EE->Print("plots/noPU_nonprompt_muon_pt_EE.pdf");
    // track
  TCanvas* c_noPU_nonprompt_track_pt = new TCanvas("c_noPU_nonprompt_track_pt", "c_noPU_nonprompt_track_pt", 1500, 1500);
  c_noPU_nonprompt_track_pt->cd();
  h_noPU_nonprompt_track_pt->GetXaxis()->SetTitle("pT (GeV)");
  h_noPU_nonprompt_track_pt->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_track_pt->SetTitle("");
  h_noPU_nonprompt_track_pt->Draw("hist");
  c_noPU_nonprompt_track_pt->Print("plots/noPU_nonprompt_track_pt.pdf");
  TCanvas* c_noPU_nonprompt_track_pt_EB = new TCanvas("c_noPU_nonprompt_track_pt_EB", "c_noPU_nonprompt_track_pt_EB", 1500, 1500);
  c_noPU_nonprompt_track_pt_EB->cd();
  h_noPU_nonprompt_track_pt_EB->GetXaxis()->SetTitle("pT (GeV)");
  h_noPU_nonprompt_track_pt_EB->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_track_pt_EB->SetTitle("");
  h_noPU_nonprompt_track_pt_EB->Draw("hist");
  c_noPU_nonprompt_track_pt_EB->Print("plots/noPU_nonprompt_track_pt_EB.pdf");
  TCanvas* c_noPU_nonprompt_track_pt_EE = new TCanvas("c_noPU_nonprompt_track_pt_EE", "c_noPU_nonprompt_track_pt_EE", 1500, 1500);
  c_noPU_nonprompt_track_pt_EE->cd();
  h_noPU_nonprompt_track_pt_EE->GetXaxis()->SetTitle("pT (GeV)");
  h_noPU_nonprompt_track_pt_EE->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_track_pt_EE->SetTitle("");
  h_noPU_nonprompt_track_pt_EE->Draw("hist");
  c_noPU_nonprompt_track_pt_EE->Print("plots/noPU_nonprompt_track_pt_EE.pdf");

  // noPU nonprompt eta
    // muon
  TCanvas* c_noPU_nonprompt_muon_eta = new TCanvas("c_noPU_nonprompt_muon_eta", "c_noPU_nonprompt_muon_eta", 1500, 1500);
  c_noPU_nonprompt_muon_eta->cd();
  h_noPU_nonprompt_muon_eta->GetXaxis()->SetTitle("#eta (GeV)");
  h_noPU_nonprompt_muon_eta->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_muon_eta->SetTitle("");
  h_noPU_nonprompt_muon_eta->Draw("hist");
  c_noPU_nonprompt_muon_eta->Print("plots/noPU_nonprompt_muon_eta.pdf");
  TCanvas* c_noPU_nonprompt_muon_eta_EB = new TCanvas("c_noPU_nonprompt_muon_eta_EB", "c_noPU_nonprompt_muon_eta_EB", 1500, 1500);
  c_noPU_nonprompt_muon_eta_EB->cd();
  h_noPU_nonprompt_muon_eta_EB->GetXaxis()->SetTitle("#eta (GeV)");
  h_noPU_nonprompt_muon_eta_EB->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_muon_eta_EB->SetTitle("");
  h_noPU_nonprompt_muon_eta_EB->Draw("hist");
  c_noPU_nonprompt_muon_eta_EB->Print("plots/noPU_nonprompt_muon_eta_EB.pdf");
  TCanvas* c_noPU_nonprompt_muon_eta_EE = new TCanvas("c_noPU_nonprompt_muon_eta_EE", "c_noPU_nonprompt_muon_eta_EE", 1500, 1500);
  c_noPU_nonprompt_muon_eta_EE->cd();
  h_noPU_nonprompt_muon_eta_EE->GetXaxis()->SetTitle("#eta (GeV)");
  h_noPU_nonprompt_muon_eta_EE->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_muon_eta_EE->SetTitle("");
  h_noPU_nonprompt_muon_eta_EE->Draw("hist");
  c_noPU_nonprompt_muon_eta_EE->Print("plots/noPU_nonprompt_muon_eta_EE.pdf");
    // track
  TCanvas* c_noPU_nonprompt_track_eta = new TCanvas("c_noPU_nonprompt_track_eta", "c_noPU_nonprompt_track_eta", 1500, 1500);
  c_noPU_nonprompt_track_eta->cd();
  h_noPU_nonprompt_track_eta->GetXaxis()->SetTitle("#eta (GeV)");
  h_noPU_nonprompt_track_eta->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_track_eta->SetTitle("");
  h_noPU_nonprompt_track_eta->Draw("hist");
  c_noPU_nonprompt_track_eta->Print("plots/noPU_nonprompt_track_eta.pdf");
  TCanvas* c_noPU_nonprompt_track_eta_EB = new TCanvas("c_noPU_nonprompt_track_eta_EB", "c_noPU_nonprompt_track_eta_EB", 1500, 1500);
  c_noPU_nonprompt_track_eta_EB->cd();
  h_noPU_nonprompt_track_eta_EB->GetXaxis()->SetTitle("#eta (GeV)");
  h_noPU_nonprompt_track_eta_EB->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_track_eta_EB->SetTitle("");
  h_noPU_nonprompt_track_eta_EB->Draw("hist");
  c_noPU_nonprompt_track_eta_EB->Print("plots/noPU_nonprompt_track_eta_EB.pdf");
  TCanvas* c_noPU_nonprompt_track_eta_EE = new TCanvas("c_noPU_nonprompt_track_eta_EE", "c_noPU_nonprompt_track_eta_EE", 1500, 1500);
  c_noPU_nonprompt_track_eta_EE->cd();
  h_noPU_nonprompt_track_eta_EE->GetXaxis()->SetTitle("#eta (GeV)");
  h_noPU_nonprompt_track_eta_EE->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_track_eta_EE->SetTitle("");
  h_noPU_nonprompt_track_eta_EE->Draw("hist");
  c_noPU_nonprompt_track_eta_EE->Print("plots/noPU_nonprompt_track_eta_EE.pdf");

  // noPU nonprompt phi
    // muon
  TCanvas* c_noPU_nonprompt_muon_phi = new TCanvas("c_noPU_nonprompt_muon_phi", "c_noPU_nonprompt_muon_phi", 1500, 1500);
  c_noPU_nonprompt_muon_phi->cd();
  h_noPU_nonprompt_muon_phi->GetXaxis()->SetTitle("#phi (GeV)");
  h_noPU_nonprompt_muon_phi->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_muon_phi->SetTitle("");
  h_noPU_nonprompt_muon_phi->Draw("hist");
  c_noPU_nonprompt_muon_phi->Print("plots/noPU_nonprompt_muon_phi.pdf");
  TCanvas* c_noPU_nonprompt_muon_phi_EB = new TCanvas("c_noPU_nonprompt_muon_phi_EB", "c_noPU_nonprompt_muon_phi_EB", 1500, 1500);
  c_noPU_nonprompt_muon_phi_EB->cd();
  h_noPU_nonprompt_muon_phi_EB->GetXaxis()->SetTitle("#phi (GeV)");
  h_noPU_nonprompt_muon_phi_EB->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_muon_phi_EB->SetTitle("");
  h_noPU_nonprompt_muon_phi_EB->Draw("hist");
  c_noPU_nonprompt_muon_phi_EB->Print("plots/noPU_nonprompt_muon_phi_EB.pdf");
  TCanvas* c_noPU_nonprompt_muon_phi_EE = new TCanvas("c_noPU_nonprompt_muon_phi_EE", "c_noPU_nonprompt_muon_phi_EE", 1500, 1500);
  c_noPU_nonprompt_muon_phi_EE->cd();
  h_noPU_nonprompt_muon_phi_EE->GetXaxis()->SetTitle("#phi (GeV)");
  h_noPU_nonprompt_muon_phi_EE->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_muon_phi_EE->SetTitle("");
  h_noPU_nonprompt_muon_phi_EE->Draw("hist");
  c_noPU_nonprompt_muon_phi_EE->Print("plots/noPU_nonprompt_muon_phi_EE.pdf");
    // track
  TCanvas* c_noPU_nonprompt_track_phi = new TCanvas("c_noPU_nonprompt_track_phi", "c_noPU_nonprompt_track_phi", 1500, 1500);
  c_noPU_nonprompt_track_phi->cd();
  h_noPU_nonprompt_track_phi->GetXaxis()->SetTitle("#phi (GeV)");
  h_noPU_nonprompt_track_phi->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_track_phi->SetTitle("");
  h_noPU_nonprompt_track_phi->Draw("hist");
  c_noPU_nonprompt_track_phi->Print("plots/noPU_nonprompt_track_phi.pdf");
  TCanvas* c_noPU_nonprompt_track_phi_EB = new TCanvas("c_noPU_nonprompt_track_phi_EB", "c_noPU_nonprompt_track_phi_EB", 1500, 1500);
  c_noPU_nonprompt_track_phi_EB->cd();
  h_noPU_nonprompt_track_phi_EB->GetXaxis()->SetTitle("#phi (GeV)");
  h_noPU_nonprompt_track_phi_EB->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_track_phi_EB->SetTitle("");
  h_noPU_nonprompt_track_phi_EB->Draw("hist");
  c_noPU_nonprompt_track_phi_EB->Print("plots/noPU_nonprompt_track_phi_EB.pdf");
  TCanvas* c_noPU_nonprompt_track_phi_EE = new TCanvas("c_noPU_nonprompt_track_phi_EE", "c_noPU_nonprompt_track_phi_EE", 1500, 1500);
  c_noPU_nonprompt_track_phi_EE->cd();
  h_noPU_nonprompt_track_phi_EE->GetXaxis()->SetTitle("#phi (GeV)");
  h_noPU_nonprompt_track_phi_EE->GetYaxis()->SetTitle("% Counts");
  h_noPU_nonprompt_track_phi_EE->SetTitle("");
  h_noPU_nonprompt_track_phi_EE->Draw("hist");
  c_noPU_nonprompt_track_phi_EE->Print("plots/noPU_nonprompt_track_phi_EE.pdf");


}



int main(int argc, char **argv)
{
//  N_muon();
  draw_iso_efficiency();
  draw_iso_distribution();
  draw_reliso_roc();
  draw_pt_roc();
  pTeff();
//  N_genMatched();
//  track_type_v1();
//  track_sigma_type_v1();
  track_type_v2();
  track_sigma_type_v2();
  muon_mother_pdgId();
  fraction_mva_cut();
  evtId_pvtrk();
  pt_eta_phi();


  return 0;

}
