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
 
using namespace std;
void bkgsig(std::string inputFile) {
  
  //get TTree from file ...
  TreeReader data(inputFile.data());

  Int_t before = 0;
  Int_t after  = 0; 

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);
    
    before++;

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

    vector<int> fatjet;
    for(int ij=0; ij<nFJets; ij++)
      {
	TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
	if(thisJet->Pt()<300)continue;
	if(fabs(thisJet->Eta())>2.4)continue;
	if(!passFatJetLooseID[ij])continue;
		
	fatjet.push_back(ij);
	
      } 
    if(fatjet.size()<2)continue;
    
    for(int ad=0; ad<nAJets; ad++)
      {
        TLorentzVector* theseJet = (TLorentzVector*)addjetP4->At(ad);
	//doubleSV

        for(unsigned int ae=0; ae<fatjet.size(); ae++) 
	  {
	    int aa = fatjet[0];
	    int ab = fatjet[1];
	    if(fatjetPRmassL2L3Corr[aa]<105 || fatjetPRmassL2L3Corr[aa]>135)continue;
	    if(fatjetPRmassL2L3Corr[ab]<105 || fatjetPRmassL2L3Corr[ab]>135)continue; //cut on corrected pruned mass for both H
	 
	    TLorentzVector* thatJet  = (TLorentzVector*)fatjetP4->At(aa);
	    TLorentzVector* thoseJet = (TLorentzVector*)fatjetP4->At(ab);
	    float dEta = fabs(thatJet->Eta() - thoseJet->Eta());
	    if(dEta>1.3)continue; //|delta Eta| < 1.3 cut between leading and subleading FAT jet 	      
	    
	    Double_t mjj = (*thatJet+*thoseJet).M(); //Mx = leading FAT jet + subleading FAT jet
	    if(mjj<1000)continue; //Mjj < 1000 GeV Cut
	       
	    if(thatJet->DeltaR(*theseJet)>0.1)continue; //delta R < 0.1 cut between leading FAT jet and ADD jet 
	    after++;
	  }

      }   

  }
  std::cout <<"before = "<<before<<std::endl;
  std::cout <<"after = "<<after<<std::endl;      

}
