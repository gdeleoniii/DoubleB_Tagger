#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF2.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TPad.h>
#include <TStyle.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <setNCUStyle.C>

void eff() {

  TFile* file[2];

  file[0]  = TFile::Open("SignalEff.root");
  file[1]  = TFile::Open("BkgEff_all.root");

  TH2D* h0  = (TH2D*)(file[0]->Get("signal"));     
  TH2D* h1  = (TH2D*)(file[1]->Get("bkg"));
 
  float bkg[19], sig[19];
  float densig = h0->Integral(1,20,1,20);
  float denbkg = h1->Integral(1,20,1,20);
  for(int i=0; i<19;i++) {
    float numsig1 = h0->Integral(i+1,20,i+1,20);
    float numbkg1 = h1->Integral(i+1,20,i+1,20);
    bkg[i]=1-(numbkg1/denbkg);
    sig[i]=numsig1/densig;
    std::cout << "Signal Eff of " << i << " to 20 = " << sig[i] << std::endl;
    std::cout << "Bkg Eff of " << i  << " to 20 = " << bkg[i] << std::endl;
  }
  
  setNCUStyle();
  TCanvas* c3 =  new TCanvas("c3","c3",0,0,600,600);                                                           
  c3->cd();
  TGraph *roc_curve = new TGraph(19,sig,bkg);
  roc_curve->SetLineWidth();
  roc_curve->SetMarkerStyle(23);
  roc_curve->SetMarkerSize(1);
  roc_curve->SetMarkerColor(kSpring-2);
  roc_curve->SetLineColor(kSpring-2);
  roc_curve->SetTitle("DoubleSV");
  roc_curve->GetXaxis()->SetLimits(0,1);
  roc_curve->SetMinimum(0);
  roc_curve->SetMaximum(1);
  roc_curve->GetYaxis()->SetTitle("1 - #varepsilon_{Bkg}");
  roc_curve->GetXaxis()->SetTitle("#varepsilon_{Sig}");
  roc_curve->Draw();
  
}
