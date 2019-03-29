// ------------------------------------------------------------------------------------
//  ROOT macro that produces average RecHit energy from PFG ntuples
//
//  Author : Ken H
//  Written on May 24, 2018
// ------------------------------------------------------------------------------------
//  
// Pre-requisite :
//
//   You should have the PFG ntuple for the Run from which you want to do a measurement. 
//   Instruction on how to make PFG ntuples can be found here : FIXME link here 
//
//   You should have "Fig" directory for plots 
//
// Usage : 
//
//   $ root -b  
//   root> .L ana_PFStudy.C+
//   root> ana_PFStudy("/cms/data/store/user/hatake/HCAL/ntuples/10_2_x/pi50_trees_MCfull_CMSSW_10_2_0_pre3_*.root","hcal_timestudy_pi50_histograms.root")
//   or
//   root> ana_PFStudy("list_trees_pi50_MCfull_CMSSW_10_2_0_pre3.txt","hcal_timestudy_pi50_histograms.root")
//   or
//   from command line:
/*
     root.exe -b -q 'ana_PFStudy.C++("trees_relval_ttbar_phase2_age_new2_4500ultimate.root","hcal_noisestudy_histograms_age_new2_4500ultimate.root")'
     root.exe -b -q 'ana_PFStudy.C++("trees_relval_ttbar_phase2_age_org.root","hcal_noisestudy_histograms_age_org.root")'
     root.exe -b -q 'ana_PFStudy.C++("trees_relval_ttbar_phase2_noage.root","hcal_noisestudy_histograms_noage.root")'

 */
//    
// -----------------------------------------------------------------------------------
// 

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip> // for setw()
#include <algorithm> 

#include "TROOT.h"
#include "TF1.h"
#include "TMath.h"
#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TCanvas.h"
#include "TDirectory.h"
#include "TBranch.h"
#include "TString.h"
#include "TStyle.h"
#include "TInterpreter.h"
#include "TStyle.h"
#include "TLorentzVector.h"

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TTreeReaderArray.h>

// In order to use vector of vectors : vector<vector<data type> >
// ACLiC makes dictionary for this
// [ref] http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=10236&p=44117#p44117
#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<int> >+;
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

using namespace std;

bool DRAWPLOTS  = true;  // draw plots or not (make "Fig" directory first before turning this on)
bool VERBOSE    = false;  // print out mean +/- sigma for each channel or not



//
// Book histograms
//
void book1D(TList *v_hist, std::string name, int n, double min, double max);
void book1DProf(TList *v_hist, std::string name, int n, double min, double max, double ymin, double ymax, Option_t *option);
void bookHistograms(TList *v_hist);

//
// Fill histograms
//
void fill1D(TList *v_hist, std::string name, double value);
void fill1DProf(TList *v_hist, std::string name, double value, double valuey);

//
// Aux
//
void relabelProfA(TList *v_hist, std::string name);

//
// Main analyzer
//
void PFCheckRun(TString rootfile, TString outfile, int maxevents=-1, int option=2) 
{ 

   cout << "[PF analyzer] Running option " << option << " for " << endl; 

   // fit pannel display option
   gStyle->SetOptFit(1011);

   //
   // Get the tree from the PFG ntuple 
   //

   TChain *ch = new TChain("hcalTupleTree/tree");

   std::string filename(rootfile);
   std::string::size_type idx;
   idx = filename.rfind('.');
   std::string extension = filename.substr(idx+1);
   std::string line;
   
   if(idx != std::string::npos && extension=="txt")
     {
       std::cout << rootfile << " " << extension << std::endl;
       std::ifstream in(rootfile);
       while (std::getline(in, line)) {     // Process line
	 if (line.size()>0) ch->Add(line.c_str());
       }
     }
   else
     {
       // No extension found
       ch->Add(rootfile);
     }

   printf("%d;\n",ch->GetNtrees());
   printf("%lld;\n",ch->GetEntries());

   TTreeReader     fReader(ch);  //!the tree reader

   //
   // Set up TTreeReader's
   // -- use MakeSelector of root
   //
   // Readers to access the data (delete the ones you do not need).
   TTreeReaderArray<double> GenParEta = {fReader, "GenParEta"};
   TTreeReaderArray<double> GenParM = {fReader, "GenParM"};
   TTreeReaderArray<double> GenParPhi = {fReader, "GenParPhi"};
   TTreeReaderArray<double> GenParPt = {fReader, "GenParPt"};
   TTreeReaderArray<double> PFParEta = {fReader, "PFParEta"};
   TTreeReaderArray<double> PFParM = {fReader, "PFParM"};
   TTreeReaderArray<double> PFParPhi = {fReader, "PFParPhi"};
   TTreeReaderArray<double> PFParPt = {fReader, "PFParPt"};
   TTreeReaderArray<int> PFParPdgId = {fReader, "PFParPdgId"};
   TTreeReaderArray<double> PFSimParEta = {fReader, "PFSimParEta"};
   TTreeReaderArray<double> PFSimParPhi = {fReader, "PFSimParPhi"};
   TTreeReaderArray<double> PFSimParPt = {fReader, "PFSimParPt"};
   // TTreeReaderArray<float> SimTracksPt = {fReader, "SimTracksPt"};

   // Define histograms to fill
   //
   TList *v_hist = new TList();
   
   bookHistograms(v_hist); // most of histograms booked here
   
   //
   // Loop over entries
   //
   unsigned int nentries = (Int_t)ch->GetEntries();
   cout << "[Hcal analyzer] The number of entries is: " << nentries << endl;

   bool debug=true;
   
   //---------------------------------------------------------------------------------------------------------
   // main event loop
   //---------------------------------------------------------------------------------------------------------

   int ievent=0;
   
   while (fReader.Next()) {
  
     // Progress indicator 
     ievent++;
     if(ievent%1000==0) cout << "[HCAL analyzer] Processed " << ievent << " out of " << nentries << " events" << endl; 
     if (maxevents>0 && ievent>maxevents) break;
     
     //--------------------
     // Loop over PF candidates
     //--------------------
     std::string strtmp;
     TLorentzVector v_gen;
     TLorentzVector v_pfsim;
     v_gen.SetPtEtaPhiM(0,0,0,0);
     v_pfsim.SetPtEtaPhiM(0,0,0,0);
     
     std::vector<int> index_pfsim_max;
     index_pfsim_max.clear();
     
     for (int igen = 0, ngen = GenParPt.GetSize(); igen<ngen; ++igen){
      
       v_gen.SetPtEtaPhiM(GenParPt[igen],GenParEta[igen],GenParPhi[igen],GenParM[igen]);

       double PFSimParPtMax =0;
       int index_pfsim_max_tmp = -1;

       // PF sim particle loop -- first
       
       for (int ipfsim = 0, npfsim = PFSimParPt.GetSize(); ipfsim<npfsim; ++ipfsim){
	 
	 v_pfsim.SetPtEtaPhiM(PFSimParPt[ipfsim],PFSimParEta[ipfsim],PFSimParPhi[ipfsim],0);

	 if (v_gen.DeltaR(v_pfsim)<0.4) {
	 if (PFSimParPt[ipfsim] > PFSimParPtMax){
	   PFSimParPtMax = PFSimParPt[ipfsim];
	   index_pfsim_max_tmp = ipfsim;
	 }
	 }
	 
       } // pfsim-loop ends      
      
       if (index_pfsim_max_tmp!=-1)
	 index_pfsim_max.push_back(index_pfsim_max_tmp);
       
     } // gen-loop ends
    
     //printf("index_max is: %8.3f",in);
     // PF sim particle loop - 2nd
     
     for (int ipfsim = 0, npfsim = PFSimParPt.GetSize(); ipfsim<npfsim; ++ipfsim){
       if (fabs(PFSimParEta[ipfsim])>2.4) continue;
       if (PFSimParPt[ipfsim]<20) continue;
      
       double dR = .1;
       double dRmin = 1.; // dummy. big value.
       v_pfsim.SetPtEtaPhiM(PFSimParPt[ipfsim],PFSimParEta[ipfsim],PFSimParPhi[ipfsim],0);
       
       TLorentzVector v_temp;
       TLorentzVector v_pf_all;
       TLorentzVector v_pf_ch_nh;
       TLorentzVector v_pf_dRmin;
       int PFParPdgId_dRmin=-1;
       v_temp.SetPtEtaPhiM(0,0,0,0);
       v_pf_all.SetPtEtaPhiM(0,0,0,0);
       v_pf_ch_nh.SetPtEtaPhiM(0,0,0,0);
       v_pf_dRmin.SetPtEtaPhiM(0,0,0,0);
       bool cHad = false;
       bool cHad_nHad = false;
       
       
       for (int ipf = 0, npf =  PFParPt.GetSize(); ipf < npf ; ++ipf) {
	 v_temp.SetPtEtaPhiM(PFParPt[ipf],PFParEta[ipf],PFParPhi[ipf],PFParM[ipf]);
	 if (v_temp.DeltaR(v_pfsim) < dR ){ 

	   v_pf_all += v_temp;	  
 
	   if (v_temp.DeltaR(v_pfsim) < dRmin){
	     dRmin = v_temp.DeltaR(v_pfsim);
	     v_pf_dRmin = v_temp;
	     PFParPdgId_dRmin = PFParPdgId[ipf];
	   }

	   if ((PFParPdgId[ipf] == 5) || (PFParPdgId[ipf] == 1)) v_pf_ch_nh += v_temp;
	   
	 }

       } //  loop-over PF candidate ends

       if (fabs(PFParPdgId_dRmin) == 1) cHad = true;
       
       if (cHad){
	 fill1D(v_hist, "All PFPar Response", (v_pf_all.Pt()-v_pfsim.Pt())/v_pfsim.Pt());
	 fill1D(v_hist, "cHad_nHad Response", (v_pf_ch_nh.Pt()-v_pfsim.Pt())/v_pfsim.Pt());
	 fill1D(v_hist, "dRmin Response", (v_pf_dRmin.Pt()-v_pfsim.Pt())/v_pfsim.Pt());
	 
	 if ( std::find(index_pfsim_max.begin(), index_pfsim_max.end(), ipfsim) != index_pfsim_max.end() ) {
	   fill1D(v_hist, "All PFPar with Max Sim Par", (v_pf_all.Pt()-v_pfsim.Pt())/v_pfsim.Pt());
	   fill1D(v_hist, "cHad_nHad with Max Sim Par", (v_pf_ch_nh.Pt()-v_pfsim.Pt())/v_pfsim.Pt());
	   fill1D(v_hist, "dRmin Response with Max Sim Par", (v_pf_dRmin.Pt()-v_pfsim.Pt())/v_pfsim.Pt());
	 }
    	   
       } // closing cHad
     
     } // ipfsim loop

   } // Event loop ends
   
   //---------------------------------------------------------------------------------------------------------
   // main event loop ends
   //---------------------------------------------------------------------------------------------------------

   // output file for histograms
   TFile file_out(outfile,"RECREATE");
   
   v_hist->Write();
    
   file_out.ls();
   file_out.Close();

}

//
// Main function
//
void ana_PFStudy_PFgeneral(TString rootfile="SinglePi_PGun_step3_RECO_10_4_0_E2_500_NoPU_v3.root",TString outfile="NoPU_SinglePi_RECO_histograms.root",int maxevents=-1)
{
  PFCheckRun(rootfile, outfile, maxevents, 0);
}



//
// --- Aux ---
//

//
// Book 1D histograms
//
void book1D(TList *v_hist, std::string name, int n, double min, double max)
{
  TH1D *htemp = new TH1D(name.c_str(), name.c_str(), n, min, max);
  v_hist->Add(htemp);
}
//
// Book 1D profile histograms
//
void book1DProf(TList *v_hist, std::string name, int n, double min, double max, double ymin, double ymax, Option_t *option="")
{
  TProfile *htemp = new TProfile(name.c_str(), name.c_str(), n, min, max, ymin, ymax, option);
  v_hist->Add(htemp);
}
//
// Book histograms
//
void bookHistograms(TList *v_hist)
{

  Char_t histo[100];
  std::string strtmp;
  
 
  //
  // Booking histograms
  // 
  //book1D(v_hist, "Response distribution", 150, -1.5, 1.5);
  book1D(v_hist, "All PFPar Response", 150, -1.5, 1.5);
  book1D(v_hist, "cHad_nHad Response", 150, -1.5, 1.5);
  book1D(v_hist, "dRmin Response", 150, -1.5, 1.5);
  book1D(v_hist, "All PFPar with Max Sim Par", 150, -1.5, 1.5);
  book1D(v_hist, "cHad_nHad with Max Sim Par", 150, -1.5, 1.5);
  book1D(v_hist, "dRmin Response with Max Sim Par", 150, -1.5, 1.5);
 

}
//
// Fill 1D histograms
//
void relabelProfA(TList *v_hist, std::string name)
{
  TProfile* htemp = (TProfile*) v_hist->FindObject(name.c_str());
  htemp->GetXaxis()->SetBinLabel(1,"Track");
  htemp->GetXaxis()->SetBinLabel(2,"ECAL");
  htemp->GetXaxis()->SetBinLabel(3,"HCAL");
  htemp->GetXaxis()->SetBinLabel(5,"HCAL d1");
  htemp->GetXaxis()->SetBinLabel(6,"HCAL d2");
  htemp->GetXaxis()->SetBinLabel(7,"HCAL d3");
  htemp->GetXaxis()->SetBinLabel(8,"HCAL d4");
  htemp->GetXaxis()->SetBinLabel(9,"HCAL d5");
  htemp->GetXaxis()->SetBinLabel(10,"HCAL d6");
  htemp->GetXaxis()->SetBinLabel(11,"HCAL d7");
}
//
// Fill 1D histograms
//
void fill1D(TList *v_hist, std::string name, double value)
{
  TH1F* htemp = (TH1F*) v_hist->FindObject(name.c_str());
  htemp->Fill(value);
}
//
// Fill 1D Profile histograms
//
void fill1DProf(TList *v_hist, std::string name, double value, double valuey)
{
  TProfile* htemp = (TProfile*) v_hist->FindObject(name.c_str());
  htemp->Fill(value,valuey);
}


