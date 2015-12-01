#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "untuplizer.h"
#include "readSample.h"

void DoubleSV(std::string inputFile, char name){

  //read the ntuples (in pcncu)

  std::vector<string> infiles;

  readSample(inputFile, infiles);
  
  TreeReader data(infiles);
  //TreeReader data(inputFile.data());
  
  // Declare the histogram

  TH1D* h_leaddoubleSV=new TH1D("","",40,-1.2,1.2);
  TH1D* h_subldoubleSV=new TH1D("","",40,-1.2,1.2);

  // begin of event loop

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev % 10000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);

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
    vector<pair<int,int>> Mjj;
    for(int ij=0; ij<nFJets; ij++) {
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
      if(thisJet->Pt()<170)continue;
      if(fabs(thisJet->Eta())>2.4)continue;
      if(!passFatJetLooseID[ij])continue;
      if(fatjetPRmassL2L3Corr[ij]<105 || fatjetPRmassL2L3Corr[ij]>135)continue;

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

    for(unsigned int ae = 0;ae<Mjj.size(); ae++) {
      int aa = Mjj[0].second;
      int ee = Mjj[0].first;
      TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa);
      TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);
      for(int ad=0; ad<nAJets; ad++) {
        TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
        if(Jet1->DeltaR(*Jet3)<0.1) {
          h_leaddoubleSV->Fill(addjet_doubleSV[ad]);                                           
	}
        if(Jet2->DeltaR(*Jet3)<0.1) {
          h_subldoubleSV->Fill(addjet_doubleSV[ad]);                                                                    
	}
      }
    }

  } // end of event loop
  fprintf(stderr, "Processed all events\n");

  TFile* outfile = new TFile(Form("doubleSVBKG_%d.root",name),"recreate");
  h_leaddoubleSV->Write(Form("dsvBKGlead_%d",name));
  h_subldoubleSV->Write(Form("dsvBKGsubl_%d",name));
  outfile->Write();
}
