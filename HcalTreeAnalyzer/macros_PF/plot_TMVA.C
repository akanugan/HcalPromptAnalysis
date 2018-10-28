{

  TFile *_file0 = TFile::Open("TMVAReg_PF.root");

  TTree *TestTree = (TTree*)_file0->Get("dataset/TestTree");

  //-----
  TCanvas *c1 = new TCanvas("c1","c1",800,600);
  gPad->SetGrid();
  TProfile *tprof = new TProfile("tprof","tprof",20,0.,100.,-1.,+2.);
  TestTree->Draw("(pf_pt/gen_pt):gen_pt>>tprof","","prof");
  TProfile *tprof2 = new TProfile("tprof2","tprof2",20,0.,100.,-1.,+2.);
  TestTree->Draw("(pf_pt/gen_pt):gen_pt>>tprof2","pf_ecalFrac<0.01","prof");
  tprof->SetMaximum(1.5);
  tprof->SetMinimum(0.0);
  tprof->Draw();
  tprof2->SetLineColor(2);
  tprof2->Draw("same");
  TLegend *legend = new TLegend(0.5,0.3,0.8,0.50);
  legend->AddEntry(tprof,"All","l");
  legend->AddEntry(tprof2,"ecalFrac<0.01","l");
  legend->Draw();
  c1->SaveAs("K0L_response_genpt.pdf");

  //-----
  TCanvas *c1b = new TCanvas("c1b","c1b",800,600);
  gPad->SetGrid();
  //TProfile *tprof = new TProfile("tprof","tprof",20,0.,100.,-1.,+2.);
  TestTree->Draw("(DNN/gen_pt):gen_pt>>tprof","","prof");
  //TProfile *tprof2 = new TProfile("tprof2","tprof2",20,0.,100.,-1.,+2.);
  TestTree->Draw("(DNN/gen_pt):gen_pt>>tprof2","pf_ecalFrac<0.01","prof");
  tprof->SetMaximum(1.5);
  tprof->SetMinimum(0.0);
  tprof->Draw();
  tprof2->SetLineColor(2);
  tprof2->Draw("same");
  //TLegend *legend = new TLegend(0.5,0.3,0.8,0.50);
  //legend->AddEntry(tprof,"All","l");
  //legend->AddEntry(tprof2,"ecalFrac<0.01","l");
  legend->Draw();
  c1b->SaveAs("K0L_response_genpt_DNN.pdf");

  //-----
  TCanvas *c2 = new TCanvas("c2","c2",800,600);
  gPad->SetGrid();
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1011);
  // Set stat options
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);                
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);
  // Set height of stat-box (fraction of pad size)
  TH1F *th1b = new TH1F("th1b","th1b",100,-0.5,0.5);
  TestTree->Draw("(DNN-gen_pt)/gen_pt>>th1b","gen_pt>50.");
  TH1F *th1 = new TH1F("th1","th1",100,-0.5,0.5);
  TestTree->Draw("(pf_pt-gen_pt)/gen_pt>>th1","gen_pt>50.");
  th1b->SetLineColor(2);
  th1b->Fit("gaus");
  th1b->GetFunction("gaus")->SetLineColor(2);
  th1->Fit("gaus");
  th1->GetFunction("gaus")->SetLineColor(4);
  th1b->Draw();
  th1->Draw("sames");
  c2->Update();
  TPaveStats *sta1  = (TPaveStats*)th1->FindObject("stats");
  TPaveStats *sta1b = (TPaveStats*)th1b->FindObject("stats");
  sta1->SetX1NDC(0.55);
  sta1->SetX2NDC(0.75);
  sta1b->SetX1NDC(0.75);
  sta1b->SetX2NDC(0.95);
  th1b->Draw();
  th1->Draw("sames");
  
  c2->SaveAs("resol_genpt50.pdf");

  //-----
  TCanvas *c2b = new TCanvas("c2b","c2b",800,600);
  gPad->SetGrid();
  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1011);
  // Set stat options
  gStyle->SetStatY(0.9);                
  gStyle->SetStatX(0.9);                
  gStyle->SetStatW(0.2);                
  gStyle->SetStatH(0.2);
  // Set height of stat-box (fraction of pad size)
  //TH1F *th1b = new TH1F("th1b","th1b",100,-0.5,0.5);
  TestTree->Draw("(DNN-gen_pt)/gen_pt>>th1b","gen_pt>30.");
  //TH1F *th1 = new TH1F("th1","th1",100,-0.5,0.5);
  TestTree->Draw("(pf_pt-gen_pt)/gen_pt>>th1","gen_pt>30.");
  th1b->SetLineColor(2);
  th1b->Fit("gaus");
  th1b->GetFunction("gaus")->SetLineColor(2);
  th1->Fit("gaus");
  th1->GetFunction("gaus")->SetLineColor(4);
  th1b->Draw();
  th1->Draw("sames");
  c2->Update();
  //TPaveStats *sta1  = (TPaveStats*)th1->FindObject("stats");
  //TPaveStats *sta1b = (TPaveStats*)th1b->FindObject("stats");
  //sta1->SetX1NDC(0.55);
  //sta1->SetX2NDC(0.75);
  //sta1b->SetX1NDC(0.75);
  //sta1b->SetX2NDC(0.95);
  th1b->Draw();
  th1->Draw("sames");
  
  c2b->SaveAs("resol_genpt30.pdf");

}
