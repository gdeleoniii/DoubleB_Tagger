#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <THnSparse.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "untuplizer.h"
#include "readSample.h"


void QCD_HT() {

  TFile* file[5];


  file[0]  = TFile::Open("HT500to700_3.0.1.root");
  file[1]  = TFile::Open("HT700to1000_3.0.1.root");
  file[2]  = TFile::Open("HT1000to1500_3.0.1.root");
  file[3]  = TFile::Open("HT1500to2000_3.0.1.root");
  file[4]  = TFile::Open("HT2000toInf_3.0.1.root");

  TH1F *h_leadDSV[5] = {NULL};
  TH1F *h_sublDSV[5] = {NULL};
  TH1F *h_leadPt[5]  = {NULL};
  TH1F *h_sublPt[5]  = {NULL};
  TH1F *h_leadPR[5]  = {NULL};
  TH1F *h_sublPR[5]  = {NULL};
  TH1F *h_leadEta[5] = {NULL};
  TH1F *h_sublEta[5] = {NULL};
  TH1F *h_Msubt[5]   = {NULL};
  TH1F *h_Mjj[5]     = {NULL};
  TH1F *h_DelEta[5]  = {NULL};


  for(int i=0;i<5;i++) {

    h_leadDSV[i]  = (TH1F*)(file[i]->Get("leadDSV"));
    h_sublDSV[i]  = (TH1F*)(file[i]->Get("sublDSV"));
    h_leadPt[i]   = (TH1F*)(file[i]->Get("leadPt"));
    h_sublPt[i]   = (TH1F*)(file[i]->Get("sublPt"));
    h_leadPR[i]   = (TH1F*)(file[i]->Get("leadPR"));
    h_sublPR[i]   = (TH1F*)(file[i]->Get("sublPR"));
    h_leadEta[i]  = (TH1F*)(file[i]->Get("leadEta"));
    h_sublEta[i]  = (TH1F*)(file[i]->Get("sublEta"));
    h_Msubt[i]    = (TH1F*)(file[i]->Get("Msubt"));
    h_Mjj[i]      = (TH1F*)(file[i]->Get("Mjj"));
    h_DelEta[i]   = (TH1F*)(file[i]->Get("DelEta"));
  }

  Double_t norm[5];
  norm[0] = 32100/14882885;
  norm[1] = 6831/14050788;
  norm[2] = 1207/14613039;
  norm[3] = 119.9/3097583;
  norm[4] = 25.24/5911831;

  TH1F* hleadDSV = (TH1F*)h_leadDSV[1]->Clone("hleadDSV");
  hleadDSV->Reset();
  TH1F* hsublDSV = (TH1F*)h_sublDSV[1]->Clone("hsublDSV");
  hsublDSV->Reset();
  TH1F* hleadPt = (TH1F*)h_leadPt[1]->Clone("hleadPt");
  hleadPt->Reset();
  TH1F* hsublPt = (TH1F*)h_sublPt[1]->Clone("hsublPt");
  hsublPt->Reset();
  TH1F* hleadPR = (TH1F*)h_leadPR[1]->Clone("hleadPR");
  hleadPR->Reset();
  TH1F* hsublPR = (TH1F*)h_sublPR[1]->Clone("hsublPR");
  hsublPR->Reset();
  TH1F* hleadEta = (TH1F*)h_leadEta[1]->Clone("hleadEta");
  hleadEta->Reset();
  TH1F* hsublEta = (TH1F*)h_sublEta[1]->Clone("hsublEta");
  hsublEta->Reset();
  TH1F* hMsubt = (TH1F*)h_Msubt[1]->Clone("hMsubt");
  hMsubt->Reset();
  TH1F* hMjj = (TH1F*)h_Mjj[1]->Clone("hMjj");
  hMjj->Reset();
  TH1F* hDelEta = (TH1F*)h_DelEta[1]->Clone("hDelEta");
  hDelEta->Reset();

  for(int x=1;x<5;x++) {
    hleadDSV->Add(h_leadDSV[x],norm[x]);
    hsublDSV->Add(h_sublDSV[x],norm[x]);
    hleadPt->Add(h_leadPt[x],norm[x]);
    hsublPt->Add(h_sublPt[x],norm[x]);
    hleadPR->Add(h_leadPR[x],norm[x]);
    hsublPR->Add(h_sublPR[x],norm[x]);
    hleadEta->Add(h_leadEta[x],norm[x]);
    hsublEta->Add(h_sublEta[x],norm[x]);
    hMsubt->Add(h_Msubt[x],norm[x]);
    hMjj->Add(h_Mjj[x],norm[x]);
    hDelEta->Add(h_DelEta[x],norm[x]);
  }

  TFile* outfile = new TFile("QCD_PythiaMC_3.0.2.root","recreate");
  hleadDSV->Write("leadDSV");
  hsublDSV->Write("sublDSV");
  hleadPt->Write("leadPt");
  hsublPt->Write("sublPt");
  hleadPR->Write("leadPR");
  hsublPR->Write("sublPR");
  hleadEta->Write("leadEta");
  hsublEta->Write("sublEta");
  hMsubt->Write("Msubt");
  hMjj->Write("Mjj");
  hDelEta->Write("DelEta");
  outfile->Write();

}
