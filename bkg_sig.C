#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TString.h>
#include <map>
#include <TH1.h>
#include <TFile.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TStyle.h>
#include <TProfile.h>
#include <TPad.h>
#include "TEfficiency.h"
#include <TStyle.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "setNCUStyle.C"
 
using namespace std;
void bkgsig(std::string inputFile) {
  
  //get TTree from file ...
  TreeReader data(inputFile.data());

  Long64_t nPass[20]={0};

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
    
    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);
    
    int nFATJet         = data.GetInt("FATnJet");
    const int nFJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
    Float_t*  fatjetPRmass = data.GetPtrFloat("FATjetPRmass");
    vector<bool>    &passFatJetLooseID = *((vector<bool>*) data.GetPtr("FATjetPassIDLoose"));
    Float_t*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    
    int nADDJet         = data.GetInt("ADDnJet");
    const int nAJets=nADDJet;
    TClonesArray* addjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    Float_t*  addjet_doubleSV = data.GetPtrFloat("ADDjet_DoubleSV");
    
    int ss=0;
    vector<int> fatjet;
    vector<pair<int,int>> Mjj;
    for(int ij=0; ij<nFJets; ij++) {
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
      if(thisJet->Pt()<200)continue;
      if(fabs(thisJet->Eta())>2.4)continue;
      if(!passFatJetLooseID[ij])continue;
      if(fatjetPRmassL2L3Corr[ij]<105 || fatjetPRmassL2L3Corr[ij]>135)continue;
      
      fatjet.push_back(ij);	
      ss++;
    }
    
    if(ss>=2)nPass[0]++; //number of events passing pre-selection
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
    
    for(unsigned int ae = 0;ae<Mjj.size(); ae++) {
      int aa = Mjj[0].second;
      int ee = Mjj[0].first;
      TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa); 
      TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);
      for(int ad=0; ad<nAJets; ad++) {
	TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
	if(Jet1->DeltaR(*Jet3)>0.1 && Jet2->DeltaR(*Jet3)>0.1)continue;

	double n = -0.9;
	for(int i=1;i<=19;i++) {
	  if(i==1) {
	    if(addjet_doubleSV[ad]>n)continue;
	    nPass[i]++;
	  }
	  else {
	    n+=0.1;
	    if(addjet_doubleSV[ad]>n)continue;
	    nPass[i]++;
          }
	}
      }
    }
    
  }

  for(int i=0;i<20;i++)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;

}
