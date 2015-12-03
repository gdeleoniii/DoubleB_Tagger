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

void DoubleSV(std::string inputFile,char name) {

  //read the ntuples (in pcncu)

  std::vector<string> infiles;

  readSample(inputFile, infiles);
  
  TreeReader data(infiles);
  //TreeReader data(inputFile.data());
  
  // Declare the histogram

  Long64_t nTotal=0;
  Long64_t nPass[20]={0};

  TH1D* h_leaddoubleSV=new TH1D("","",40,-1.2,1.2);
  TH1D* h_subldoubleSV=new TH1D("","",40,-1.2,1.2);

  // begin of event loop

  for( Long64_t ev = 0; ev < data.GetEntriesFast(); ev++ ){

    if( ev % 10000 == 0 )
      fprintf(stderr, "Processing event %lli of %lli\n", ev + 1, data.GetEntriesFast());

    data.GetEntry(ev);
    nTotal++;

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
    int fjetselection= 0;
    vector<pair<int,int>> Mjj;
    for(int ij=0; ij<nFJets; ij++) {
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
      if(thisJet->Pt()<170)continue;
      if(fabs(thisJet->Eta())>2.4)continue;
      if(!passFatJetLooseID[ij])continue;
      if(fatjetPRmassL2L3Corr[ij]<105 || fatjetPRmassL2L3Corr[ij]>135)continue;

      fatjet.push_back(ij);
      fjetselection++;
    }

    if(fjetselection>=2)nPass[0]++;
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
	nPass[1]++;
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
	  nPass[2]++;                                           
	}
        if(Jet2->DeltaR(*Jet3)<0.1) {
          h_subldoubleSV->Fill(addjet_doubleSV[ad]);                                                                    
	  nPass[3]++;
	}
      }
    }

  } // end of event loop
  fprintf(stderr, "Processed all events\n");

  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++) 
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;

  TH1D* h_CutFlow = new TH1D(" ", "" ,4,0,4);
     
  char* cut_name[4] = {"FATJet","Mjj","LeadDSV","SublDSV"};  

  for(int i=1;i<=4;i++){ // i is the index of column of cut flow plot
    if(i==1) {h_CutFlow->SetBinContent(i,nPass[0]); }
    else {h_CutFlow->SetBinContent(i,nPass[i-1]); }
    h_CutFlow->GetXaxis()->SetBinLabel( i , cut_name[i-1] );
  }

  TFile* outfile = new TFile(Form("counting_%d.root",name),"recreate");
  h_CutFlow->Write(Form("count_%d",name));
  outfile->Write();

  /*
  TFile* outfile = new TFile(Form("doubleSVBKG_%d.root",name),"recreate");
  h_leaddoubleSV->Write(Form("dsvBKGlead_%d",name));
  h_subldoubleSV->Write(Form("dsvBKGsubl_%d",name));
  outfile->Write();
  */
}
