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
#include <THnSparse.h>
#include <TProfile.h>
#include <TPad.h>
#include <TStyle.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <setNCUStyle.C>

void eff() {

  TFile* file[6];

  file[0]  = TFile::Open("SignalEff.root");
  file[1]  = TFile::Open("BkgEff_all.root");

  TH2D* h0  = (TH2D*)(file[0]->Get("doublebtagging"));  //dbt signal   
  TH2D* h1  = (TH2D*)(file[0]->Get("fatjetcsv")); //csv signal
  TH2D* h2  = (TH2D*)(file[1]->Get("doublebtagging")); //dbt background
  TH2D* h3  = (TH2D*)(file[1]->Get("fatjetcsv")); //csv background
  THnSparseD* h4  = (THnSparseD*)(file[1]->Get("subjetcsv"));
  THnSparseD* h5  = (THnSparseD*)(file[0]->Get("subjetcsv"));
 
  //--------------For DBT and fatjetCSV------------
  float bkg0[19], sig0[19],bkg1[19], sig1[19],bkg2[9], sig2[9];
  float densig0 = h0->Integral(1,20,1,20); //denominator dbt signal
  float denbkg0 = h2->Integral(1,20,1,20); //denominator dbt background
  float densig1 = h1->Integral(1,20,1,20); //denominator csv signal
  float denbkg1 = h3->Integral(1,20,1,20); //denominator csv background
   for(int i=0; i<19;i++) {
    float numsig0 = h0->Integral(i+1,20,i+1,20);
    float numbkg0 = h2->Integral(i+1,20,i+1,20);
    float numsig1 = h1->Integral(i+1,20,i+1,20);
    float numbkg1 = h3->Integral(i+1,20,i+1,20);
    bkg0[i]=1-(numbkg0/denbkg0);
    sig0[i]=numsig0/densig0;
    bkg1[i]=1-(numbkg1/denbkg1);
    sig1[i]=numsig1/densig1;
   }


  //---------For subjetcsv-----------------
  Long64_t ncutsig[9]= {0};
  Long64_t ncutbkg[9]= {0};
  for(int q=0;q<10;q++)
    for(int i=q+1;i<=10; i++)
      for(int j=q+1; j<=10; j++)
        for(int k=q+1; k<=10; k++)
          for(int m=q+1; m<=10; m++)
            {
              int index[4]={i,j,k,m};
              ncutsig[q] += h4->GetBinContent(index);
	      ncutbkg[q] += h5->GetBinContent(index);
            }

  float densig2 = ncutsig[0];
  float denbkg2 = ncutbkg[0];
  for(int q=0;q<9;q++) {
    float numsig2 = ncutsig[q+1];
    float numbkg2 = ncutbkg[q+1];
    sig2[q] = numsig2/densig2;
    bkg2[q] = 1-(numbkg2/denbkg2);
  }

  setNCUStyle();
  TCanvas* c3 =  new TCanvas("c3","c3",0,0,600,600);                                                           
  c3->cd();
  TGraph *roc_curve = new TGraph(19,sig0,bkg0);
  TGraph *roc2 = new TGraph(19,sig1,bkg1);
  TGraph *roc3 = new TGraph(9,sig2,bkg2);
  roc_curve->SetLineWidth(3);
  roc_curve->SetLineStyle(8);
  roc_curve->SetMarkerStyle(23);
  roc_curve->SetMarkerSize(1.2);
  roc_curve->SetMarkerColor(kAzure+2);
  roc_curve->SetLineColor(kAzure-3);
  roc_curve->SetTitle("DoubleSV");
  roc_curve->GetXaxis()->SetLimits(0,1);
  roc_curve->SetMinimum(0);
  roc_curve->SetMaximum(1);
  roc_curve->GetYaxis()->SetTitle("1 - #varepsilon_{Bkg}");
  roc_curve->GetXaxis()->SetTitle("#varepsilon_{Sig}");
  roc_curve->Draw("apc");
  roc2->SetLineWidth(3);
  roc2->SetLineStyle(2);
  roc2->SetMarkerStyle(20);
  roc2->SetMarkerSize(1.2);
  roc2->SetMarkerColor(kBlue+2);
  roc2->SetLineColor(kBlue-3);
  roc2->Draw("pc");
  roc3->SetLineWidth(3);
  roc3->SetLineStyle(6);
  roc3->SetMarkerStyle(21);
  roc3->SetMarkerSize(1.2);
  roc3->SetMarkerColor(kTeal+5);
  roc3->SetLineColor(kTeal-3);
  roc3->Draw("pc");

  TLegend *legend2 = new TLegend(0.16,0.40,.61,0.49);
  legend2->AddEntry(roc_curve,"Double B-tagging","lep");
  legend2->AddEntry(roc2,"FAT Jet CSV", "lep");
  legend2->AddEntry(roc3,"with 4 FAT Sub-jet CSV > 0.605", "lep");
  legend2->SetFillStyle(0);
  legend2->Draw();

}
