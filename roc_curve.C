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
#include <TStyle.h>
#include <TStyle.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "setNCUStyle.C"
 
using namespace std;
void bkgsig(std::string inputFile) {
  
  //get TTree from file ...
  TreeReader data(inputFile.data());

  TH2F *h2_doubleSV   = new TH2F("","",20,-1,1,20,-1,1);

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
    
    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);
    
    int nFATJet         = data.GetInt("FATnJet");
    const int nFJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCSV    = data.GetPtrFloat("FATjetCSV");
    Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
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
      if(thisJet->Pt()<200)continue;
      if(fabs(thisJet->Eta())>2.4)continue;
      if(!passFatJetLooseID[ij])continue;
      if(fatjetPRmassL2L3Corr[ij]<105 || fatjetPRmassL2L3Corr[ij]>135)continue;
      //if(fatjetCSV[ij]<0.605)continue; //2nd roc curve
      if( FATnSubSDJet[ij] != 2 ) continue; //3rd roc curve
      if( FATsubjetSDCSV[ij][0] < 0.605 || FATsubjetSDCSV[ij][1] < 0.605 ) continue; //3rd roc curve 
      
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
    
    int addJetIndex[2]={-1,-1};
    for(unsigned int ae = 0;ae<Mjj.size(); ae++) {
      int aa = Mjj[0].second;
      int ee = Mjj[0].first;
      TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa); 
      TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);
      for(int ad=0; ad<nAJets; ad++) {
	TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
	if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;} // first add jet to pass the delta r cut
	if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;} // first add jet to pass the delta r cut
      }
    }
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;
    h_leaddoubleSV->Fill(addjet_doubleSV[addJetIndex[0]]);
    h_subldoubleSV->Fill(addjet_doubleSV[addJetIndex[1]]);
    h2_doubleSV->Fill(addjet_doubleSV[ addJetIndex[0]], addjet_doubleSV[addJetIndex[1]]);
  
  
  }
  
  TFile* outfile = new TFile("SignalEff3.root","recreate");    
  h2_doubleSV->Write("signal3");
  outfile->Write();
  
}
