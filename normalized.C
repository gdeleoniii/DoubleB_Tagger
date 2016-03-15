#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TH2.h>
#include <THnSparse.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "untuplizer.h"
#include "readSample.h"

void bkg() {

  TFile* file[8];

  //file[0]  = TFile::Open("BkgEff_1.root");
  //file[1]  = TFile::Open("BkgEff_2.root");
  //file[2]  = TFile::Open("BkgEff_3.root");
  file[3]  = TFile::Open("BkgEff_4.root");
  file[4]  = TFile::Open("BkgEff_5.root");
  file[5]  = TFile::Open("BkgEff_6.root");
  file[6]  = TFile::Open("BkgEff_7.root");
  file[7]  = TFile::Open("BkgEff_8.root");
 
  //TH2D* h2_dbtn1  = (TH2D*)(file[0]->Get("doublebtagging"));
  //TH2D* h2_fatn1  = (TH2D*)(file[0]->Get("fatjetcsv"));
  //TH2D* h2_dbtn2  = (TH2D*)(file[1]->Get("doublebtagging"));
  //TH2D* h2_fatn2  = (TH2D*)(file[1]->Get("fatjetcsv")); 
  //TH2D* h2_dbtn3  = (TH2D*)(file[2]->Get("doublebtagging"));
  //TH2D* h2_fatn3  = (TH2D*)(file[2]->Get("fatjetcsv"));
  TH2D* h2_dbtn4  = (TH2D*)(file[3]->Get("doublebtagging"));
  TH2D* h2_fatn4  = (TH2D*)(file[3]->Get("fatjetcsv"));
  TH2D* h2_dbtn5  = (TH2D*)(file[4]->Get("doublebtagging"));
  TH2D* h2_fatn5  = (TH2D*)(file[4]->Get("fatjetcsv"));
  TH2D* h2_dbtn6  = (TH2D*)(file[5]->Get("doublebtagging"));
  TH2D* h2_fatn6  = (TH2D*)(file[5]->Get("fatjetcsv"));
  TH2D* h2_dbtn7  = (TH2D*)(file[6]->Get("doublebtagging"));
  TH2D* h2_fatn7  = (TH2D*)(file[6]->Get("fatjetcsv"));
  TH2D* h2_dbtn8  = (TH2D*)(file[7]->Get("doublebtagging"));
  TH2D* h2_fatn8  = (TH2D*)(file[7]->Get("fatjetcsv"));

  //THnSparseD* hs_subn1  = (THnSparseD*)(file[0]->Get("subjetcsv"));
  //THnSparseD* hs_subn2  = (THnSparseD*)(file[1]->Get("subjetcsv"));
  //THnSparseD* hs_subn3  = (THnSparseD*)(file[2]->Get("subjetcsv"));
  THnSparseD* hs_subn4  = (THnSparseD*)(file[3]->Get("subjetcsv"));
  THnSparseD* hs_subn5  = (THnSparseD*)(file[4]->Get("subjetcsv"));
  THnSparseD* hs_subn6  = (THnSparseD*)(file[5]->Get("subjetcsv"));
  THnSparseD* hs_subn7  = (THnSparseD*)(file[6]->Get("subjetcsv"));
  THnSparseD* hs_subn8  = (THnSparseD*)(file[7]->Get("subjetcsv"));

  //----------double b-tagging------

  //Double_t sdbt1=27990000/h2_dbtn1->Integral(1,20,1,20);
  //Double_t sdbt2=1712000/h2_dbtn2->Integral(1,20,1,20);
  //Double_t sdbt3=347700/h2_dbtn3->Integral(1,20,1,20);
  Double_t sdbt4=32100/h2_dbtn4->Integral(1,20,1,20);
  Double_t sdbt5=6831/h2_dbtn5->Integral(1,20,1,20);
  Double_t sdbt6=1207/h2_dbtn6->Integral(1,20,1,20);
  Double_t sdbt7=119.9/h2_dbtn7->Integral(1,20,1,20);
  Double_t sdbt8=25.24/h2_dbtn8->Integral(1,20,1,20);

  TH2D* h2_dbtnadd = (TH2D*)h2_dbtn4->Clone("h2_dbtnadd");
  h2_dbtnadd->Reset();

  //h2_dbtnadd->Add(h2_dbtn1,sdbt1);
  //h2_dbtnadd->Add(h2_dbtn2,sdbt2);
  //h2_dbtnadd->Add(h2_dbtn3,sdbt3);
  h2_dbtnadd->Add(h2_dbtn4,sdbt4);
  h2_dbtnadd->Add(h2_dbtn5,sdbt5);
  h2_dbtnadd->Add(h2_dbtn6,sdbt6);
  h2_dbtnadd->Add(h2_dbtn7,sdbt7);
  h2_dbtnadd->Add(h2_dbtn8,sdbt8);

  //-------fatjet csv-------
  //Double_t sfat1=27990000/h2_fatn1->Integral(1,20,1,20);
  //Double_t sfat2=1712000/h2_fatn2->Integral(1,20,1,20);
  //Double_t sfat3=347700/h2_fatn3->Integral(1,20,1,20);
  Double_t sfat4=32100/h2_fatn4->Integral(1,20,1,20);
  Double_t sfat5=6831/h2_fatn5->Integral(1,20,1,20);
  Double_t sfat6=1207/h2_fatn6->Integral(1,20,1,20);
  Double_t sfat7=119.9/h2_fatn7->Integral(1,20,1,20);
  Double_t sfat8=25.24/h2_fatn8->Integral(1,20,1,20);

  TH2D* h2_fatnadd = (TH2D*)h2_fatn4->Clone("h2_fatnadd");
  h2_fatnadd->Reset();

  //h2_fatnadd->Add(h2_fatn1,sfat1);
  //h2_fatnadd->Add(h2_fatn2,sfat2);
  //h2_fatnadd->Add(h2_fatn3,sfat3);
  h2_fatnadd->Add(h2_fatn4,sfat4);
  h2_fatnadd->Add(h2_fatn5,sfat5);
  h2_fatnadd->Add(h2_fatn6,sfat6);
  h2_fatnadd->Add(h2_fatn7,sfat7);
  h2_fatnadd->Add(h2_fatn8,sfat8);

  //---------subjet csv-------------
  //Double_t ssub1=27990000/hs_subn1->GetEntries();
  //Double_t ssub2=1712000/hs_subn2->GetEntries();
  //Double_t ssub3=347700/hs_subn3->GetEntries();
  Double_t ssub4=32100/hs_subn4->GetEntries();
  Double_t ssub5=6831/hs_subn5->GetEntries();
  Double_t ssub6=1207/hs_subn6->GetEntries();
  Double_t ssub7=119.9/hs_subn7->GetEntries();
  Double_t ssub8=25.24/hs_subn8->GetEntries();

  THnSparseD* hs_subnadd = (THnSparseD*)hs_subn4->Clone("hs_subnadd");
  hs_subnadd->Reset();

  //hs_subnadd->Add(hs_subn1,ssub1);
  //hs_subnadd->Add(hs_subn2,ssub2);
  //hs_subnadd->Add(hs_subn3,ssub3);
  hs_subnadd->Add(hs_subn4,ssub4);
  hs_subnadd->Add(hs_subn5,ssub5);
  hs_subnadd->Add(hs_subn6,ssub6);
  hs_subnadd->Add(hs_subn7,ssub7);
  hs_subnadd->Add(hs_subn8,ssub8);


  TFile* outfile = new TFile("BkgEffnadd.root","recreate");
  h2_dbtnadd->Write("doublebtagging");
  h2_fatnadd->Write("fatjetcsv");
  hs_subnadd->Write("subjetcsv");
  outfile->Write();

}
