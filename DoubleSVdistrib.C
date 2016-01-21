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

void DoubleSV(std::string inputFile,char name) {

  //read the ntuples (in pcncu)

  std::vector<string> infiles;

  readSample(inputFile, infiles);
  
  TreeReader data(infiles);
   
  // Declare the histogram

  //TH1D* h_leaddoubleSV=new TH1D("","",40,-1.2,1.2);
  //TH1D* h_subldoubleSV=new TH1D("","",40,-1.2,1.2);
  TH2F *h2_doubleSV   = new TH2F("","",20,-1,1,20,-1,1);
  TH2F *h2_fatjetCSV = new TH2F("","",20,0,1,20,0,1);

  const int nDim=4;
  int nbins[nDim]={10,10,10,10};
  double xmin[nDim]={0,0,0,0};
  double xmax[nDim]={1,1,1,1};
  THnSparseD* hs = new THnSparseD("hs", "hs", nDim, nbins, xmin, xmax);

  // begin of event loop

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev % 10000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

    int nFATJet         = data.GetInt("FATnJet");
    const int nFJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCSV    = data.GetPtrFloat("FATjetCSV");
    //Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    Int_t*   FATnSubSDJet   = data.GetPtrInt("FATnSubSDJet");
    vector<float>* FATsubjetSDCSV       = data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJet);

    int nADDJet         = data.GetInt("ADDnJet");
    const int nAJets=nADDJet;
    TClonesArray* addjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    Float_t*  addjet_doubleSV = data.GetPtrFloat("ADDjet_DoubleSV");

    vector<int> fatjet;
    vector<pair<int,int>> Mjj;
    for(int ij=0; ij<nFJets; ij++) {
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
      if(thisJet->Pt()<170)continue;
      if(fabs(thisJet->Eta())>2.4)continue;
      if(!passFatJetLooseID[ij])continue;
      if(fatjetPRmassL2L3Corr[ij]<105 || fatjetPRmassL2L3Corr[ij]>135)continue;
      //if(fatjetCSV[ij]<0.605)continue;                                                                                         
      //if( FATnSubSDJet[ij] != 2 ) continue;
      //if( FATsubjetSDCSV[ij][0] < 0.605 || FATsubjetSDCSV[ij][1] < 0.605 ) continue;

      fatjet.push_back(ij);
    }

    if(fatjet.size()<2)continue;

    for(unsigned int i=0; i<fatjet.size(); i++) {
      for(unsigned int j=0; j<i; j++) {
        int index_that = fatjet[i];
        int index_those = fatjet[j];
        TLorentzVector* thatJet  = (TLorentzVector*)fatjetP4->At(index_that);
        TLorentzVector* thoseJet = (TLorentzVector*)fatjetP4->At(index_those);
        float dEta = fabs(thatJet->Eta() - thoseJet->Eta());
        if(dEta>1.3)continue;

        Double_t mjj = (*thatJet+*thoseJet).M();
        if(mjj<1000)continue;
        Mjj.push_back(make_pair(index_that,index_those));
      }
    }

    if(Mjj.size()<1)continue;
    int addJetIndex[2]={-1,-1};
    int aa = Mjj[0].second;
    int ee = Mjj[0].first;
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa);
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);
    h2_fatjetCSV->Fill(fatjetCSV[aa],fatjetCSV[ee]);
    
    for(int ad=0; ad<nAJets; ad++) {
      TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
      if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;}
      if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;}
    }
    
    
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;
    h2_doubleSV->Fill(addjet_doubleSV[ addJetIndex[0]], addjet_doubleSV[addJetIndex[1]]);

    int fatjetIndex[2]={Mjj[0].second, Mjj[0].first};
    double subjetCSV[4]={-1,-1,-1,-1};
    for(int i=0; i<2; i++)
      {
        int ijet = fatjetIndex[i];
        for(int isub=0; isub < FATnSubSDJet[ijet]; isub++)
          {
            int vectorID = i*2 + isub;
            subjetCSV[vectorID]=    FATsubjetSDCSV[ijet][isub];
	  }
      }
    hs->Fill(subjetCSV);
  

  } // end of event loop
  fprintf(stderr, "Processed all events\n");

  Double_t scale;
  Double_t scale1;
  Double_t scale2;

  if(name == 1) {
    scale=2799000/h2_doubleSV->Integral();
    scale1=27990000/h2_fatjetCSV->Integral();
    scale2=27990000/hs->GetEntries();
  }
  else if(name == 2) {
    scale=1712000/h2_doubleSV->Integral();
    scale1=1712000/h2_fatjetCSV->Integral();
    scale2=1712000/hs->GetEntries();
  }
  else if(name == 3) {
    scale=347700/h2_doubleSV->Integral();
    scale1=347700/h2_fatjetCSV->Integral();
    scale2=347700/hs->GetEntries();
  }
  else if(name == 4) {
    scale=32100/h2_doubleSV->Integral();
    scale1=32100/h2_fatjetCSV->Integral();
    scale2=32100/hs->GetEntries();
  }
  else if(name == 5) {
    scale=6831/h2_doubleSV->Integral();
    scale1=6831/h2_fatjetCSV->Integral();
    scale2=6831/hs->GetEntries();
  }
  else if(name == 6) {
    scale=1207/h2_doubleSV->Integral();
    scale1=1207/h2_fatjetCSV->Integral();
    scale2=1207/hs->GetEntries();
  }
  else if(name == 7) {
    scale=119.9/h2_doubleSV->Integral();
    scale1=119.9/h2_fatjetCSV->Integral();
    scale2=119.9/hs->GetEntries();
  }
  else if(name == 8) {
    scale=25.24/h2_doubleSV->Integral();
    scale1=25.24/h2_fatjetCSV->Integral();
    scale2=25.24/hs->GetEntries();
  }

  h2_doubleSV->Scale(scale);
  h2_fatjetCSV->Scale(scale1);
  hs->Scale(scale2);

  TFile* outfile = new TFile(Form("BkgEff_%d.root",name),"recreate");
  h2_doubleSV->Write("doublebtagging");
  h2_fatjetCSV->Write("fatjetcsv");
  hs->Write("subjetcsv");
  outfile->Write();
}
