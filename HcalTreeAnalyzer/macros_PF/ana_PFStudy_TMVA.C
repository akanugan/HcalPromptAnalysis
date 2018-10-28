// ------------------------------------------------------------------------------------
//  ROOT macro example
//
//  Author : Ken H
//  Written on May 24, 2018
// ------------------------------------------------------------------------------------
//  
// Usage : 
//
//   $ root -b  
//   root> .L ana_PFStudy_TMVA.C+
//   root> ana_PFStudy_TMVA("relval_ttbar_2019.root","pfstudy_histograms_2019.root");
//   or
//   root> ana_PFStudy_TMVA();
//   or
//   from command line:
/*
     root.exe -b -q 'ana_PFStudy_TMVA.C+("relval_ttbar_2019.root","pfstudy_histograms_2019.root")'
 */
//    
// -----------------------------------------------------------------------------------
// 

#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <iomanip> // test
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

#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/TMVARegGui.h"

// In order to use vector of vectors : vector<vector<data type> >
// ACLiC makes dictionary for this
// [ref] http://root.cern.ch/phpBB3/viewtopic.php?f=3&t=10236&p=44117#p44117
#ifdef __MAKECINT__
#pragma link C++ class std::vector < std::vector<int> >+;
#pragma link C++ class std::vector < std::vector<float> >+;
#endif

using namespace std;

bool DRAWPLOTS  = false;  // draw plots or not (make "Fig" directory first before turning this on)
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
   /*
   TTreeReaderArray<float> HBHERecHitEnergy = {fReader, "HBHERecHitEnergy"};
   TTreeReaderArray<float> HBHERecHitEta = {fReader, "HBHERecHitEta"};
   TTreeReaderArray<float> HBHERecHitPhi = {fReader, "HBHERecHitPhi"};
   TTreeReaderArray<float> HBHERecHitTime = {fReader, "HBHERecHitTime"};
   */
   TTreeReaderArray<float> PFParEcalEnergyFrac = {fReader, "PFParEcalEnergyFrac"};
   TTreeReaderArray<float> PFParHOEnergyFrac = {fReader, "PFParHOEnergyFrac"};
   TTreeReaderArray<float> PFParHcalEnergyFrac = {fReader, "PFParHcalEnergyFrac"};
   TTreeReaderArray<float> PFParHcalFrac1 = {fReader, "PFParHcalFrac1"};
   TTreeReaderArray<float> PFParHcalFrac2 = {fReader, "PFParHcalFrac2"};
   TTreeReaderArray<float> PFParHcalFrac3 = {fReader, "PFParHcalFrac3"};
   TTreeReaderArray<float> PFParHcalFrac4 = {fReader, "PFParHcalFrac4"};
   TTreeReaderArray<float> PFParHcalFrac5 = {fReader, "PFParHcalFrac5"};
   TTreeReaderArray<float> PFParHcalFrac6 = {fReader, "PFParHcalFrac6"};
   TTreeReaderArray<float> PFParHcalFrac7 = {fReader, "PFParHcalFrac7"};
   TTreeReaderArray<float> PFParTrackPt = {fReader, "PFParTrackPt"};
   TTreeReaderArray<int> GenParPdgId = {fReader, "GenParPdgId"};
   TTreeReaderArray<int> GenParStatus = {fReader, "GenParStatus"};
   /*
   TTreeReaderArray<int> HBHERecHitAux = {fReader, "HBHERecHitAux"};
   TTreeReaderArray<int> HBHERecHitDepth = {fReader, "HBHERecHitDepth"};
   TTreeReaderArray<int> HBHERecHitFlags = {fReader, "HBHERecHitFlags"};
   TTreeReaderArray<int> HBHERecHitHPDid = {fReader, "HBHERecHitHPDid"};
   TTreeReaderArray<int> HBHERecHitIEta = {fReader, "HBHERecHitIEta"};
   TTreeReaderArray<int> HBHERecHitIPhi = {fReader, "HBHERecHitIPhi"};
   TTreeReaderArray<int> HBHERecHitRBXid = {fReader, "HBHERecHitRBXid"};
   */
   TTreeReaderArray<int> PFParPdgId = {fReader, "PFParPdgId"};
   TTreeReaderArray<int> PFParStatus = {fReader, "PFParStatus"};
   TTreeReaderValue<UInt_t> bx = {fReader, "bx"};
   TTreeReaderValue<UInt_t> event = {fReader, "event"};
   TTreeReaderValue<UInt_t> ls = {fReader, "ls"};
   TTreeReaderValue<UInt_t> orbit = {fReader, "orbit"};
   TTreeReaderValue<UInt_t> run = {fReader, "run"};

   //
   // Define histograms to fill
   //
   TList *v_hist = new TList();
   
   bookHistograms(v_hist); // most of histograms booked here
   
   //
   // Loop over entries
   //
   unsigned int nentries = (Int_t)ch->GetEntries();
   cout << "[Hcal analyzer] The number of entries is: " << nentries << endl;

   //

   //TString outfileName2( "SingleK0L.root" );
   //TFile* outputFile2 = TFile::Open( outfileName2, "RECREATE" ); 

   TTree t1("t1","a simple Tree with simple variables");
 
   // Create one branch. If splitlevel is set, event is a superbranch
   // creating a sub branch for each data member of the Event object.

   Float_t gen_pt, gen_e, gen_eta, gen_phi;   
   t1.Branch("gen_pt", &gen_pt, "gen_pt/F");
   t1.Branch("gen_e",  &gen_e,  "gen_e/F");
   t1.Branch("gen_eta",&gen_eta,"gen_eta/F");
   t1.Branch("gen_phi",&gen_phi,"gen_phi/F");

   Float_t pf_pt, pf_e, pf_eta, pf_phi;   
   t1.Branch("pf_pt", &pf_pt, "pf_pt/F");
   t1.Branch("pf_e",  &pf_e,  "pf_e/F");
   t1.Branch("pf_eta",&pf_eta,"pf_eta/F");
   t1.Branch("pf_phi",&pf_phi,"pf_phi/F");

   Float_t npf, npf_nh;
   t1.Branch("npf",   &npf,   "npf/F");
   t1.Branch("npf_nh",&npf_nh,"npf_nh/F");
   
   Float_t pf_ecalFrac, pf_hcalFrac, pf_hcalFrac1, pf_hcalFrac2, pf_hcalFrac3, pf_hcalFrac4, pf_hcalFrac5, pf_hcalFrac6, pf_hcalFrac7;
   t1.Branch("pf_ecalFrac",&pf_ecalFrac,"pf_ecalFrac/F");
   t1.Branch("pf_hcalFrac",&pf_hcalFrac,"pf_hcalFrac/F");
   t1.Branch("pf_hcalFrac1",&pf_hcalFrac1,"pf_hcalFrac1/F");
   t1.Branch("pf_hcalFrac2",&pf_hcalFrac2,"pf_hcalFrac2/F");
   t1.Branch("pf_hcalFrac3",&pf_hcalFrac3,"pf_hcalFrac3/F");
   t1.Branch("pf_hcalFrac4",&pf_hcalFrac4,"pf_hcalFrac4/F");
   t1.Branch("pf_hcalFrac5",&pf_hcalFrac5,"pf_hcalFrac5/F");
   t1.Branch("pf_hcalFrac6",&pf_hcalFrac6,"pf_hcalFrac6/F");
   t1.Branch("pf_hcalFrac7",&pf_hcalFrac7,"pf_hcalFrac7/F");
   
   //---------------------------------------------------------------------------------------------------------
   // main event loop
   //---------------------------------------------------------------------------------------------------------

   int ievent=0;
   while (fReader.Next()) {
  
     // Progress indicator 
     ievent++;
     if(ievent%10==0) cout << "[HCAL analyzer] Processed " << ievent << " out of " << nentries << " events" << endl; 
     if (maxevents>0 && ievent>maxevents) break;

     bool debug=false;
     
     //--------------------
     // Loop over GEN particles
     //--------------------
     for (int ipar = 0, npar =  GenParPt.GetSize(); ipar < npar; ++ipar) {
     if (GenParPdgId[ipar]==130){
	 //std::cout << GenParPdgId[ipar] << std::endl;
       
       if (debug) printf("ipar,pt,eta,phi:      %4d %8.2f %8.2f %8.2f %8.2f %5d\n",ipar,GenParPt[ipar],GenParEta[ipar],GenParPhi[ipar],GenParM[ipar],GenParPdgId[ipar]);
       TLorentzVector tlv_k0l; tlv_k0l.SetPtEtaPhiM(GenParPt[ipar],GenParEta[ipar],GenParPhi[ipar],GenParM[ipar]);

       gen_pt  = GenParPt[ipar];
       gen_e   = tlv_k0l.E();
       gen_eta = GenParEta[ipar];
       gen_phi = GenParPhi[ipar];
       
       //--------------------
       // Loop over PF candidates
       //--------------------
       TLorentzVector tlv_pf_sum(0.,0.,0.,0.);
       TLorentzVector tlv_pf_sum_ecal(0.,0.,0.,0.);
       TLorentzVector tlv_pf_sum_hcal(0.,0.,0.,0.);
       TLorentzVector tlv_pf_sum_hcal1(0.,0.,0.,0.);
       TLorentzVector tlv_pf_sum_hcal2(0.,0.,0.,0.);
       TLorentzVector tlv_pf_sum_hcal3(0.,0.,0.,0.);
       TLorentzVector tlv_pf_sum_hcal4(0.,0.,0.,0.);
       TLorentzVector tlv_pf_sum_hcal5(0.,0.,0.,0.);
       TLorentzVector tlv_pf_sum_hcal6(0.,0.,0.,0.);
       TLorentzVector tlv_pf_sum_hcal7(0.,0.,0.,0.);
       bool has_charge = false;
       bool has_photon = false;
       int Npf=0;
       int Npf_nh=0;
       for (int ipfcand = 0, npfcand =  PFParPt.GetSize(); ipfcand < npfcand; ++ipfcand) {
	 TLorentzVector tlv_pf; tlv_pf.SetPtEtaPhiM(PFParPt[ipfcand],PFParEta[ipfcand],PFParPhi[ipfcand],PFParM[ipfcand]);

	 /*
	 if (fabs(PFParEta[ipfcand])>2.5 && fabs(PFParEta[ipfcand])<2.8 && tlv_k0l.DeltaR(tlv_pf)<0.1){
	   printf("aaa  ipfcand,pt,eta,phi: %10d %4d %8.2f %8.2f %8.2f %8.2f %5d\n",ievent,ipfcand,PFParPt[ipfcand],PFParEta[ipfcand],PFParPhi[ipfcand],PFParM[ipfcand],PFParPdgId[ipfcand]);
	   std::cout << PFParTrackPt[ipfcand] << " " << PFParEcalEnergyFrac[ipfcand] << " "
		     << PFParHcalEnergyFrac[ipfcand] << " " << PFParHOEnergyFrac[ipfcand] << std::endl;
	 }
	 */

	 //if (fabs(PFParPdgId[ipfcand])!=5) continue;
	 if (tlv_k0l.DeltaR(tlv_pf)<0.1){
	   Npf++;	   
	   if (fabs(PFParPdgId[ipfcand])==5) Npf_nh++;
	   tlv_pf_sum+=tlv_pf;
	   tlv_pf_sum_ecal+=tlv_pf*PFParEcalEnergyFrac[ipfcand];
	   tlv_pf_sum_hcal+=tlv_pf*PFParHcalEnergyFrac[ipfcand];
	   tlv_pf_sum_hcal1+=tlv_pf*PFParHcalEnergyFrac[ipfcand]*PFParHcalFrac1[ipfcand];
	   tlv_pf_sum_hcal2+=tlv_pf*PFParHcalEnergyFrac[ipfcand]*PFParHcalFrac2[ipfcand];
	   tlv_pf_sum_hcal3+=tlv_pf*PFParHcalEnergyFrac[ipfcand]*PFParHcalFrac3[ipfcand];
	   tlv_pf_sum_hcal4+=tlv_pf*PFParHcalEnergyFrac[ipfcand]*PFParHcalFrac4[ipfcand];
	   tlv_pf_sum_hcal5+=tlv_pf*PFParHcalEnergyFrac[ipfcand]*PFParHcalFrac5[ipfcand];
	   tlv_pf_sum_hcal6+=tlv_pf*PFParHcalEnergyFrac[ipfcand]*PFParHcalFrac6[ipfcand];
	   tlv_pf_sum_hcal7+=tlv_pf*PFParHcalEnergyFrac[ipfcand]*PFParHcalFrac7[ipfcand];
	   if (fabs(PFParPdgId[ipfcand])<=3) has_charge=true;
	   if (fabs(PFParPdgId[ipfcand])==4) has_photon=true;
	   if (debug) { 
	     printf("  ipfcand,pt,eta,phi: %4d %8.2f %8.2f %8.2f %8.2f %5d\n",ipfcand,PFParPt[ipfcand],PFParEta[ipfcand],PFParPhi[ipfcand],PFParM[ipfcand],PFParPdgId[ipfcand]);
	     std::cout << PFParTrackPt[ipfcand] << " " << PFParEcalEnergyFrac[ipfcand] << " "
		       << PFParHcalEnergyFrac[ipfcand] << " " << PFParHOEnergyFrac[ipfcand] << std::endl;
	   } // debug
	 }   // matched to gen par
       }     // loop over PF candidates   

       if (has_charge) continue; // don't want to use this K0L that has a matched charged PF candidate

       pf_pt  = tlv_pf_sum.Pt();
       pf_e   = tlv_pf_sum.E();
       pf_eta = tlv_pf_sum.Eta();
       pf_phi = tlv_pf_sum.Phi();
       pf_ecalFrac  = tlv_pf_sum_ecal.Pt()/tlv_pf_sum.Pt();
       pf_hcalFrac  = tlv_pf_sum_hcal.Pt()/tlv_pf_sum.Pt();
       pf_hcalFrac1 = tlv_pf_sum_hcal1.Pt()/tlv_pf_sum.Pt();
       pf_hcalFrac2 = tlv_pf_sum_hcal2.Pt()/tlv_pf_sum.Pt();
       pf_hcalFrac3 = tlv_pf_sum_hcal3.Pt()/tlv_pf_sum.Pt();
       pf_hcalFrac4 = tlv_pf_sum_hcal4.Pt()/tlv_pf_sum.Pt();
       pf_hcalFrac5 = tlv_pf_sum_hcal5.Pt()/tlv_pf_sum.Pt();
       pf_hcalFrac6 = tlv_pf_sum_hcal6.Pt()/tlv_pf_sum.Pt();
       pf_hcalFrac7 = tlv_pf_sum_hcal7.Pt()/tlv_pf_sum.Pt();
       if (debug) {
	 if (fabs(tlv_pf_sum.Eta())>2.5) std::cout << npf << std::endl;
       }
	 
       // plotting
       std::string strtmp;
       if (fabs(GenParEta[ipar])<1.0){
	 strtmp = "PFTask_response_vs_pt_eta0to1";
	 fill1DProf(v_hist, strtmp, GenParPt[ipar],tlv_pf_sum.Pt()/GenParPt[ipar]);
	 strtmp = "PFTask_response_vs_pt_eta0to1_PFK0L";
	 if (!has_photon) fill1DProf(v_hist, strtmp, GenParPt[ipar],tlv_pf_sum.Pt()/GenParPt[ipar]);
       }
       if (fabs(GenParEta[ipar])>1.4 && fabs(GenParEta[ipar])<2.2){
	 strtmp = "PFTask_response_vs_pt_eta14to22";	 
	 fill1DProf(v_hist, strtmp, GenParPt[ipar],tlv_pf_sum.Pt()/GenParPt[ipar]);
	 strtmp = "PFTask_response_vs_pt_eta14to22_PFK0L";	 
	 if (!has_photon) fill1DProf(v_hist, strtmp, GenParPt[ipar],tlv_pf_sum.Pt()/GenParPt[ipar]);
       }
       if (fabs(GenParEta[ipar])>2.5 && fabs(GenParEta[ipar])<2.8){
	 strtmp = "PFTask_response_vs_pt_eta25to28";	 
	 fill1DProf(v_hist, strtmp, GenParPt[ipar],tlv_pf_sum.Pt()/GenParPt[ipar]);
	 strtmp = "PFTask_response_vs_pt_eta25to28_PFK0L";	 
	 if (pf_ecalFrac<0.01) fill1DProf(v_hist, strtmp, GenParPt[ipar],tlv_pf_sum.Pt()/GenParPt[ipar]);
       }

       npf=Npf;
       npf_nh=Npf_nh;
       
       if (pf_pt>0.) t1.Fill(); // fill tree
       
     } // GenK0L
     } // gen particle loop

     //--------------------
     // Loop over PF candidates
     //--------------------
     for (int ipfcand = 0, npfcand =  PFParPt.GetSize(); ipfcand < npfcand; ++ipfcand) {
       
       std::string strtmp;

       strtmp = "PFTask_hcalFrac1Zero_vs_pt";
       float zero1=0.;
       if (PFParHcalFrac1[ipfcand]==0.) zero1=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero1);

       strtmp = "PFTask_hcalFrac2Zero_vs_pt";
       float zero2=0.;
       if (PFParHcalFrac2[ipfcand]==0.) zero2=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero2);

       strtmp = "PFTask_hcalFrac3Zero_vs_pt";
       float zero3=0.;
       if (PFParHcalFrac3[ipfcand]==0.) zero3=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero3);

       strtmp = "PFTask_hcalFrac4Zero_vs_pt";
       float zero4=0.;
       if (PFParHcalFrac4[ipfcand]==0.) zero4=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero4);

       strtmp = "PFTask_hcalFrac5Zero_vs_pt";
       float zero5=0.;
       if (PFParHcalFrac5[ipfcand]==0.) zero5=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero5);

       strtmp = "PFTask_hcalFrac6Zero_vs_pt";
       float zero6=0.;
       if (PFParHcalFrac6[ipfcand]==0.) zero6=1.;
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero6);

       strtmp = "PFTask_hcalFrac7Zero_vs_pt";
       float zero7=0.;
       if (PFParHcalFrac7[ipfcand]==0.) zero7=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]),zero7);

       strtmp = "PFTask_hcalFracAllZero_vs_pt";
       float zeroAll=0.;
       if (PFParHcalFrac1[ipfcand]==0.
	   && PFParHcalFrac2[ipfcand]==0. 
	   && PFParHcalFrac3[ipfcand]==0. 
	   && PFParHcalFrac4[ipfcand]==0. 
	   && PFParHcalFrac5[ipfcand]==0. 
	   && PFParHcalFrac6[ipfcand]==0. 
	   && PFParHcalFrac7[ipfcand]==0.) zeroAll=1.; 
       fill1DProf(v_hist, strtmp, log10(PFParPt[ipfcand]), zeroAll);

       //
       // Endcap
       // 
       if ( fabs(PFParEta[ipfcand])<2.9&&fabs(PFParEta[ipfcand])>1.5){

	 //
	 // Charged hadrons
	 //
	 if ( PFParPdgId[ipfcand]==1&&fabs(PFParEta[ipfcand])<2.5 ){ 
	   if ( PFParPt[ipfcand]>5.){
	     /*
	     std::cout << PFParPt[ipfcand] << " " << PFParEta[ipfcand] << " " << PFParPdgId[ipfcand] << std::endl;
	     std::cout << PFParTrackPt[ipfcand] << " " << PFParEcalEnergyFrac[ipfcand] << " "
		       << PFParHcalEnergyFrac[ipfcand] << " " << PFParHOEnergyFrac[ipfcand] << std::endl;
	     std::cout << PFParHcalFrac1[ipfcand] << " "
		       << PFParHcalFrac2[ipfcand] << " " 
		       << PFParHcalFrac3[ipfcand] << " " 
		       << PFParHcalFrac4[ipfcand] << " " 
		       << PFParHcalFrac5[ipfcand] << " " 
		       << PFParHcalFrac6[ipfcand] << " " 
		       << PFParHcalFrac7[ipfcand] << std::endl;
	     std::cout << std::endl;
	     */
	   }

	   //--- 
	   if ( PFParPt[ipfcand]>5. ){
	     strtmp = "PFTask_Profile_ChargedHadron_Endcap_PtAbove5";
	     fill1DProf(v_hist, strtmp, -3, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	     fill1DProf(v_hist, strtmp, -2, PFParEcalEnergyFrac[ipfcand]);
	     fill1DProf(v_hist, strtmp, -1, PFParHcalEnergyFrac[ipfcand]);
	     fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	     fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	     fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	     fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	     fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	     fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	     fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);

	     /*
	     std::cout << PFParPt[ipfcand] << " " << PFParEta[ipfcand] << " " << PFParPdgId[ipfcand] << std::endl;	   
	     std::cout << PFParTrackPt[ipfcand] << " " << PFParEcalEnergyFrac[ipfcand] << " "
		       << PFParHcalEnergyFrac[ipfcand] << " " << PFParHOEnergyFrac[ipfcand] << std::endl;
	     std::cout << PFParHcalFrac1[ipfcand] << " "
		       << PFParHcalFrac2[ipfcand] << " " 
		       << PFParHcalFrac3[ipfcand] << " " 
		       << PFParHcalFrac4[ipfcand] << " " 
		       << PFParHcalFrac5[ipfcand] << " " 
		       << PFParHcalFrac6[ipfcand] << " " 
		     << PFParHcalFrac7[ipfcand] << std::endl;
	     std::cout << std::endl;
	     */

	     strtmp = "PFTask_ChargedHadron_TrackPtRatio_PtAbove5";
	     fill1D(v_hist, strtmp, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);

	     if (PFParTrackPt[ipfcand]!=PFParPt[ipfcand]){
	       strtmp = "PFTask_Profile_ChargedHadron_Endcap_PtAbove5_Special";
	       fill1DProf(v_hist, strtmp, -3, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	       fill1DProf(v_hist, strtmp, -2, PFParEcalEnergyFrac[ipfcand]);
	       fill1DProf(v_hist, strtmp, -1, PFParHcalEnergyFrac[ipfcand]);
	       fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	       fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	       fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	       fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	       fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	       fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	       fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);
	     }
	     
	   } else if ( PFParPt[ipfcand]>1.){
	     
	     strtmp = "PFTask_Profile_ChargedHadron_Endcap_Pt1To5";
	     fill1DProf(v_hist, strtmp, -3, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	     fill1DProf(v_hist, strtmp, -2, PFParEcalEnergyFrac[ipfcand]);
	     fill1DProf(v_hist, strtmp, -1, PFParHcalEnergyFrac[ipfcand]);
	     fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	     fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	     fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	     fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	     fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	     fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	     fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);

	     strtmp = "PFTask_ChargedHadron_TrackPtRatio_Pt1To5";
	     fill1D(v_hist, strtmp, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	     
	     if (PFParTrackPt[ipfcand]!=PFParPt[ipfcand]){
	       strtmp = "PFTask_Profile_ChargedHadron_Endcap_Pt1To5_Special";
	       fill1DProf(v_hist, strtmp, -3, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	       fill1DProf(v_hist, strtmp, -2, PFParEcalEnergyFrac[ipfcand]);
	       fill1DProf(v_hist, strtmp, -1, PFParHcalEnergyFrac[ipfcand]);
	       fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	       fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	       fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	       fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	       fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	       fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	       fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);
	     }

	   } else {
	     
	     strtmp = "PFTask_Profile_ChargedHadron_Endcap_PtBelow1";
	     fill1DProf(v_hist, strtmp, -3, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	     fill1DProf(v_hist, strtmp, -2, PFParEcalEnergyFrac[ipfcand]);
	     fill1DProf(v_hist, strtmp, -1, PFParHcalEnergyFrac[ipfcand]);
	     fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	     fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	     fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	     fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	     fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	     fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	     fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);
	     
	     strtmp = "PFTask_ChargedHadron_TrackPtRatio_PtBelow1";
	     fill1D(v_hist, strtmp, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);

	     if (PFParTrackPt[ipfcand]!=PFParPt[ipfcand]){
	       strtmp = "PFTask_Profile_ChargedHadron_Endcap_PtBelow1_Special";
	       fill1DProf(v_hist, strtmp, -3, PFParTrackPt[ipfcand]/PFParPt[ipfcand]);
	       fill1DProf(v_hist, strtmp, -2, PFParEcalEnergyFrac[ipfcand]);
	       fill1DProf(v_hist, strtmp, -1, PFParHcalEnergyFrac[ipfcand]);
	       fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	       fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	       fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	       fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	       fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	       fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	       fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);
	     }

	   }
	   
	 //
	 // Neutral hadrons
	 //
	 } else if ( PFParPdgId[ipfcand]==5 ){

	   if ( PFParPt[ipfcand]>5.){
	     strtmp = "PFTask_hcalProfile_NeutralHadron_Endcap_PtAbove5";
	     fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	     fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	     fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	     fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	     fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	     fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	     fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);
	     
	     /*
	     std::cout << PFParPt[ipfcand] << " " << PFParEta[ipfcand] << " " << PFParPdgId[ipfcand] << std::endl;	   
	     std::cout << PFParTrackPt[ipfcand] << " " << PFParEcalEnergyFrac[ipfcand] << " "
		       << PFParHcalEnergyFrac[ipfcand] << " " << PFParHOEnergyFrac[ipfcand] << std::endl;
	     std::cout << PFParHcalFrac1[ipfcand] << " "
		       << PFParHcalFrac2[ipfcand] << " " 
		       << PFParHcalFrac3[ipfcand] << " " 
		       << PFParHcalFrac4[ipfcand] << " " 
		       << PFParHcalFrac5[ipfcand] << " " 
		       << PFParHcalFrac6[ipfcand] << " " 
		     << PFParHcalFrac7[ipfcand] << std::endl;
	     std::cout << std::endl;
	     */
	     
	   } else if ( PFParPt[ipfcand]>1.){
	     
	     strtmp = "PFTask_hcalProfile_NeutralHadron_Endcap_Pt1To5";
	     fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	     fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	     fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	     fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	     fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	     fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	     fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);

	   } else {
	     
	     strtmp = "PFTask_hcalProfile_NeutralHadron_Endcap_PtBelow1";
	     fill1DProf(v_hist, strtmp, 1., PFParHcalFrac1[ipfcand]);
	     fill1DProf(v_hist, strtmp, 2., PFParHcalFrac2[ipfcand]);
	     fill1DProf(v_hist, strtmp, 3., PFParHcalFrac3[ipfcand]);
	     fill1DProf(v_hist, strtmp, 4., PFParHcalFrac4[ipfcand]);
	     fill1DProf(v_hist, strtmp, 5., PFParHcalFrac5[ipfcand]);
	     fill1DProf(v_hist, strtmp, 6., PFParHcalFrac6[ipfcand]);
	     fill1DProf(v_hist, strtmp, 7., PFParHcalFrac7[ipfcand]);
	     
	   }
	   
	 } // Neutral Hadron
       } // Endcap
     } // PF candidiate loop
     
   }   // Event loop ends
   //---------------------------------------------------------------------------------------------------------
   // main event loop ends
   //---------------------------------------------------------------------------------------------------------

   // output file for histograms
   TFile file_out(outfile,"RECREATE");
   
   t1.Write();     
   //outputFile2->Close();
   
   v_hist->Write();
   
   file_out.ls();
   file_out.Close();

}

//
// Main function
//
//void ana_PFStudy_TMVA(TString rootfile="relval_ttbar_2018_pmx25ns.root",TString outfile="pfstudy_histograms.root",int maxevents=-1)
//void ana_PFStudy_TMVA(TString rootfile="relval_ttbar_2018_pmx25ns.root",TString outfile="pfstudy_histograms.root",int maxevents=100)
void ana_PFStudy_TMVA(TString rootfile="filelist_K0L.txt",TString outfile="pfstudy_histograms_K0S.root",int maxevents=-1)
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
  sprintf(histo, "PFTask_hcalFrac1Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFrac2Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFrac3Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFrac4Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFrac5Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFrac6Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFrac7Zero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);
  sprintf(histo, "PFTask_hcalFracAllZero_vs_pt");
  book1DProf(v_hist, histo, 50., -1, 4., -1., 2.);

  sprintf(histo, "PFTask_ChargedHadron_TrackPtRatio_PtAbove5");
  book1D(v_hist, histo, 70., -0.2, 1.2);
  sprintf(histo, "PFTask_ChargedHadron_TrackPtRatio_Pt1To5");
  book1D(v_hist, histo, 70., -0.2, 1.2);
  sprintf(histo, "PFTask_ChargedHadron_TrackPtRatio_PtBelow1");
  book1D(v_hist, histo, 70., -0.2, 1.2);
  
  sprintf(histo,"PFTask_Profile_ChargedHadron_Endcap_PtAbove5");
  book1DProf(v_hist, histo, 11, -3.5, 7.5, -1., 2., "S");
  TProfile* htemp = (TProfile*) v_hist->FindObject(histo);
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
  htemp->Print("all");
  sprintf(histo,"PFTask_Profile_ChargedHadron_Endcap_Pt1To5");
  book1DProf(v_hist, histo, 11, -3.5, 7.5, -1., 2., "S");
  htemp = (TProfile*) v_hist->FindObject(histo);
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
  htemp->Print("all");
  sprintf(histo,"PFTask_Profile_ChargedHadron_Endcap_PtBelow1");
  book1DProf(v_hist, histo, 11, -3.5, 7.5, -1., 2., "S");
  htemp = (TProfile*) v_hist->FindObject(histo);
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
  htemp->Print("all");

  // PF candidate Pt not 
  sprintf(histo,"PFTask_Profile_ChargedHadron_Endcap_PtAbove5_Special");
  book1DProf(v_hist, histo, 11, -3.5, 7.5, -1., 2., "S");
  htemp = (TProfile*) v_hist->FindObject(histo);
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
  htemp->Print("all");
  sprintf(histo,"PFTask_Profile_ChargedHadron_Endcap_Pt1To5_Special");
  book1DProf(v_hist, histo, 11, -3.5, 7.5, -1., 2., "S");
  htemp = (TProfile*) v_hist->FindObject(histo);
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
  htemp->Print("all");
  sprintf(histo,"PFTask_Profile_ChargedHadron_Endcap_PtBelow1_Special");
  book1DProf(v_hist, histo, 11, -3.5, 7.5, -1., 2., "S");
  htemp = (TProfile*) v_hist->FindObject(histo);
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
  htemp->Print("all");
  
  sprintf(histo,"PFTask_hcalProfile_NeutralHadron_Endcap_PtAbove5");
  book1DProf(v_hist, histo, 7, 0.5,7.5, -1., 2., "S");
  sprintf(histo,"PFTask_hcalProfile_NeutralHadron_Endcap_Pt1To5");
  book1DProf(v_hist, histo, 7, 0.5,7.5, -1., 2., "S");
  sprintf(histo,"PFTask_hcalProfile_NeutralHadron_Endcap_PtBelow1");
  book1DProf(v_hist, histo, 7, 0.5,7.5, -1., 2., "S");

  sprintf(histo,"PFTask_response_vs_pt_eta0to1");
  book1DProf(v_hist, histo, 20, 0., 200., -1., 2., "S");
  sprintf(histo,"PFTask_response_vs_pt_eta14to22");	 
  book1DProf(v_hist, histo, 20, 0., 200., -1., 2., "S");
  sprintf(histo,"PFTask_response_vs_pt_eta25to28");	 
  book1DProf(v_hist, histo, 20, 0., 200., -1., 2., "S");
  
  sprintf(histo,"PFTask_response_vs_pt_eta0to1_PFK0L");
  book1DProf(v_hist, histo, 20, 0., 200., -1., 2., "S");
  sprintf(histo,"PFTask_response_vs_pt_eta14to22_PFK0L");	 
  book1DProf(v_hist, histo, 20, 0., 200., -1., 2., "S");
  sprintf(histo,"PFTask_response_vs_pt_eta25to28_PFK0L");	 
  book1DProf(v_hist, histo, 20, 0., 200., -1., 2., "S");
  
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
