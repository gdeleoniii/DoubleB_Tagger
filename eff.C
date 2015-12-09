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

  TH2D* h0  = (TH2D*)(file[0]->Get("signal")); //2D histogram of the lead addjet doubleSV vs sublead addjet doubleSV for signal Mx=1.8 TeV     
  TH2D* h1  = (TH2D*)(file[1]->Get("bkg")); // 2D histogram of the lead addjet doubleSV vs sublead addjet doubleSV for QCD bkg


  float bkg[19],sig[19];
  double n = 19;
  float densig = h0->Integral(1,20,1,20);
  float denbkg = h1->Integral(1,20,1,20);
  for(int i=1;i<=19;i++) {
    if(i==1) {
      float numsig1 = h0->Integral(n,20,n,20);
      float numbkg1 = h1->Integral(n,20,n,20);
      bkg[i]= 1-(numbkg1/denbkg);
      sig[i]=numsig2/densig;
      std::cout << "Signal Eff of " << n  << " to 20 = " << sig[i] << std::endl;
      std::cout << "Bkg Eff of " << n  << " to 20 = " << bkg[i] << std::endl;
    }
    else {
      n-=1;
      float numsig2 = h0->Integral(n,20,n,20);
      float numbkg2 = h1->Integral(n,20,n,20);
      bkg[i]= 1-(numbkg2/denbkg);
      sig[i]=numsig2/densig;
      std::cout << "Signal Eff of " << n  << " to 20 = " << sig[i] << std::endl;
      std::cout << "Bkg Eff of " << n  << " to 20 = " << bkg[i] << std::endl;
    }
  }
 

  //Values for signal efficiency = {0.315137, 0.424725, 0.499836, 0.565178, 0.615334, 0.655804, 0.692743, 0.725332, 0.759235, 0.786899, 0.811854, 0.836808, 0.86045, 0.880397, 0.902725, 0.925628, 0.95124, 0.977426, 1};
  //Values for 1 - background efficiency = {0.999511, 0.998887, 0.997856, 0.99696,  0.995386, 0.993785, 0.991424, 0.988629, 0.985345, 0.981247, 0.975955, 0.968573, 0.959725, 0.947323, 0.928542, 0.901322, 0.843624, 0.692024,0};

  //setNCUStyle();
  TCanvas* c3 =  new TCanvas("c3","c3",0,0,600,600);                                                           
  c3->cd();
  TGraph *roc_curve = new TGraph(19,sig,bkg);
  roc_curve->SetLineWidth(3);
  roc_curve->SetMarkerStyle(23);
  roc_curve->SetMarkerSize(1);
  roc_curve->SetMarkerColor(kSpring-2);
  roc_curve->SetLineColor(kSpring-2);
  roc_curve->SetTitle("DoubleSV");
  roc_curve->GetXaxis()->SetLimits(0,1);
  roc_curve->SetMinimum(0);
  roc_curve->SetMaximum(1);
  roc_curve->GetXaxis()->SetTitle("varepsilon_{sig}");
  roc_curve->GetYaxis()->SetTitle("1 - #varepsilon_{bkg}");
  roc_curve->Draw("ac");
}