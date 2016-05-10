#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1.h>
#include <TH2.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TPad.h>
#include <THnSparse.h>
#include <TStyle.h>
#include <TStyle.h>
#include "setNCUStyle.C"
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>

void dataMCratio(std::string variable) {
 
  TFile* file[4];

  file[0]  = TFile::Open("data_1.root");
  file[1]  = TFile::Open("data_2.root");
  file[2]  = TFile::Open("data_3.root");
  file[3]  = TFile::Open("qcd.root");

  TH1F* h_data1     = (TH1F*)(file[0]->Get(Form("%s",variable.data())));
  TH1F* h_data2     = (TH1F*)(file[1]->Get(Form("%s",variable.data())));
  TH1F* h_data3     = (TH1F*)(file[2]->Get(Form("%s",variable.data())));
  TH1F* h_bkg       = (TH1F*)(file[3]->Get(Form("%s",variable.data())));
 
  TH1F* h_data = (TH1F*)h_data1->Clone("h_data");
  h_data->Reset();
  h_data->Add(h_data1);
  h_data->Add(h_data2);
  h_data->Add(h_data3);


  TH1F* h_ratio = (TH1F*)h_data->Clone("h_ratio");

  h_ratio->Reset();

  Int_t nbin = h_ratio->GetNbinsX();
  Double_t ratio[nbin];
  Double_t error[nbin];
  Double_t numer_nbincontent[nbin];
  Double_t denom_nbincontent[nbin];
  Double_t numer_binerror[nbin];
  Double_t denom_binerror[nbin];

  for(Int_t i = 1; i <= nbin; i++){

    numer_nbincontent[i] = h_data->GetBinContent(i);
    denom_nbincontent[i] = h_bkg ->GetBinContent(i);
    numer_binerror[i]    = h_data->GetBinError(i);
    denom_binerror[i]    = h_bkg->GetBinError(i);

    if( denom_nbincontent[i] <= 0 || numer_nbincontent[i] <= 0 ) continue;
    if( denom_binerror[i] <= 0 || numer_binerror[i] <= 0 ) continue;

    ratio[i] = (Double_t)numer_nbincontent[i]/denom_nbincontent[i];
    error[i] = (ratio[i])*sqrt(pow(numer_binerror[i]/numer_nbincontent[i],2)+pow(denom_binerror[i]/denom_nbincontent[i],2));

    h_ratio->SetBinContent(i,ratio[i]);
    h_ratio->SetBinError(i,error[i]);

  }

  setNCUStyle(true);

  Double_t up_height     = 0.8;
  Double_t dw_correction = 1.455;
  Double_t dw_height     = (1-up_height)*dw_correction;

  TCanvas c("c","",0,0,900,900);
  c.Divide(1,2);

  TPad* c_up = (TPad*) c.GetListOfPrimitives()->FindObject("c_1");
  TPad* c_dw = (TPad*) c.GetListOfPrimitives()->FindObject("c_2"); 

  c_up->SetPad(0,1-up_height,1,1);
  c_dw->SetPad(0,0,1,dw_height);
  c_dw->SetBottomMargin(0.25);

  c_up->cd();
  //c_up->SetLogy();

  h_data->SetLineColor(kBlack);
  h_data->SetMarkerStyle(8);
  h_data->SetMarkerSize(0.8);
  h_data->GetYaxis()->SetTitleOffset(1.3);
  h_data->GetYaxis()->SetTitle("Number of Events");
  h_data->GetXaxis()->SetTitle("");
  h_data->GetXaxis()->SetLabelOffset(999);
  h_data->GetXaxis()->SetLabelSize(0);
  h_data->Draw("el");
  h_bkg->SetLineColor(kAzure+2);
  h_bkg->SetFillColor(kAzure+2);
  h_bkg->Draw("same");
  c_up->RedrawAxis();


  TLegend *leg = new TLegend(0.81, 0.53, 0.94, 0.63);
  leg->SetBorderSize(0);
  leg->SetFillColor(0);
  leg->SetFillStyle(0);
  leg->SetTextSize(0.04);
  leg->AddEntry(h_bkg, "QCD", "f");
  leg->AddEntry(h_data, "Data", "lp");
  leg->Draw();

  c_dw->cd();
  h_ratio->SetLineColor(kBlack);
  h_ratio->SetMarkerStyle(8);
  h_ratio->SetMarkerSize(0.8);
  h_ratio->SetTitle("");
  h_ratio->GetYaxis()->SetTitle("data/MC");
  h_ratio->GetYaxis()->SetTitleOffset(0.45);
  h_ratio->GetXaxis()->SetLabelSize(0.1);
  h_ratio->GetXaxis()->SetLabelOffset(0.005);
  h_ratio->GetXaxis()->SetTitleSize(0.125);
  h_ratio->GetXaxis()->SetTitleOffset(0.8);
  h_ratio->GetYaxis()->SetLabelSize(0.1);
  h_ratio->GetYaxis()->SetTitleSize(0.1);
  //h_ratio->GetYaxis()->SetNdivisions(505);
  //h_ratio->GetYaxis()->SetRangeUser(0,2);
  h_ratio->Draw();


  c.Print(Form("%s.pdf",variable.data()));


}
