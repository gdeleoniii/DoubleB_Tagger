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
  THnSparseD* h4  = (THnSparseD*)(file[0]->Get("subjetcsv"));
  THnSparseD* h5  = (THnSparseD*)(file[1]->Get("subjetcsv"));
  THnSparseD* h6  = (THnSparseD*)(file[0]->Get("subjetcsv1"));
  THnSparseD* h7  = (THnSparseD*)(file[1]->Get("subjetcsv1"));
  THnSparseD* h8  = (THnSparseD*)(file[0]->Get("subjetcsv2"));
  THnSparseD* h9  = (THnSparseD*)(file[1]->Get("subjetcsv2"));
 
  float bkg0[19], sig0[19],bkg1[19], sig1[19],bkg2[10], sig2[10], bkg3[10], sig3[10], bkg4[10], sig4[10];
  
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
  Long64_t ncutsig[10]= {0};
  Long64_t ncutbkg[10]= {0};
  Long64_t ncutsig1[10]= {0};
  Long64_t ncutbkg1[10]= {0};
  Long64_t ncutsig2[10]= {0};
  Long64_t ncutbkg2[10]= {0};
  for(int q=0;q<10;q++) {
    for(int i=q+1;i<=10; i++) {
      for(int j=q+1; j<=10; j++) {

	int index2[2]={i,j};
	ncutsig2[q] += h8->GetBinContent(index2);
	ncutbkg2[q] += h9->GetBinContent(index2);

        for(int k=q+1; k<=10; k++) {
	  
	  int index1[3]={i,j,k};
	  ncutsig1[q] += h6->GetBinContent(index1);
	  ncutbkg1[q] += h7->GetBinContent(index1);
          
	  for(int m=q+1; m<=10; m++) {
              int index[4]={i,j,k,m};
              ncutsig[q] += h4->GetBinContent(index);
	      ncutbkg[q] += h5->GetBinContent(index);
	  }
	}
      }
    }
  }
  float densig2 = ncutsig[0];
  float denbkg2 = ncutbkg[0];
  float densig3 = ncutsig1[0];
  float denbkg3 = ncutbkg1[0];
  float densig4 = ncutsig2[0];
  float denbkg4 = ncutbkg2[0];
  for(int q=0;q<10;q++) {
    float numsig2 = ncutsig[q];
    float numbkg2 = ncutbkg[q];
    float numsig3 = ncutsig1[q];
    float numbkg3 = ncutbkg1[q];
    float numsig4 = ncutsig2[q];
    float numbkg4 = ncutbkg2[q];
    sig2[q] = numsig2/densig2;
    bkg2[q] = 1-(numbkg2/denbkg2);
    sig3[q] = numsig3/densig3;
    bkg3[q] = 1-(numbkg3/denbkg3);
    sig4[q] = numsig4/densig4;
    bkg4[q] = 1-(numbkg4/denbkg4);
  }



  setNCUStyle();
  TCanvas* c3 =  new TCanvas("c3","c3",0,0,600,600);                                                           
  c3->cd();
  TGraph *roc_curve = new TGraph(19,sig0,bkg0);
  TGraph *roc2 = new TGraph(19,sig1,bkg1);
  TGraph *roc3 = new TGraph(10,sig2,bkg2);
  TGraph *roc4 = new TGraph(10,sig3,bkg3);
  TGraph *roc5 = new TGraph(10,sig4,bkg4);
  roc_curve->SetLineWidth(2);
  roc_curve->SetLineStyle(5);
  roc_curve->SetMarkerStyle(23);
  roc_curve->SetMarkerSize(1.1);
  roc_curve->SetMarkerColor(kRed+2);
  roc_curve->SetLineColor(kRed-3);
  roc_curve->SetTitle("DoubleSV");
  roc_curve->GetXaxis()->SetLimits(0,1);
  roc_curve->SetMinimum(0);
  roc_curve->SetMaximum(1);
  roc_curve->GetYaxis()->SetTitle("1 - #varepsilon_{Bkg}");
  roc_curve->GetXaxis()->SetTitle("#varepsilon_{Sig}");
  roc_curve->Draw("apc");
  roc2->SetLineWidth(2);
  roc2->SetLineStyle(3);
  roc2->SetMarkerStyle(20);
  roc2->SetMarkerSize(0.8);
  roc2->SetMarkerColor(kBlue+2);
  roc2->SetLineColor(kBlue-3);
  roc2->Draw("pc");
  roc3->SetLineWidth(2);
  roc3->SetLineStyle(6);
  roc3->SetMarkerStyle(21);
  roc3->SetMarkerSize(0.8);
  roc3->SetMarkerColor(kGreen+2);
  roc3->SetLineColor(kGreen-3);
  roc3->Draw("pc");
  roc4->SetLineWidth(2);
  roc4->SetLineStyle(6);
  roc4->SetMarkerStyle(22);
  roc4->SetMarkerSize(1.1);
  roc4->SetMarkerColor(kYellow);
  roc4->SetLineColor(kYellow-7);
  roc4->Draw("pc");
  roc5->SetLineWidth(2);
  roc5->SetLineStyle(2);
  roc5->SetMarkerStyle(29);
  roc5->SetMarkerSize(1.3);
  roc5->SetMarkerColor(kOrange+8);
  roc5->SetLineColor(kOrange+7);
  roc5->Draw("pc");

  TLegend *legend2 = new TLegend(0.16,0.40,.61,0.49);
  //legend2->SetHeader("ADD Jet Double SV");
  legend2->AddEntry(roc_curve,"Double B-tagging","lep");
  legend2->AddEntry(roc2,"Fatjet CSV", "lep");
  legend2->AddEntry(roc3,"4 Sub-jet CSV", "lep");
  legend2->AddEntry(roc4,"3 Sub-jet CSV", "lep");
  legend2->AddEntry(roc5,"2 Sub-jet CSV", "lep");
  legend2->SetFillStyle(0);
  legend2->Draw();

}
