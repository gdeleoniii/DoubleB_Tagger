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
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "readSample.h"
 
using namespace std;
void dataMC(std::string inputFile, char name) {

    //read the ntuples (in pcncu)

    std::vector<string> infiles;

    readSample(inputFile, infiles);
  
    TreeReader data(infiles);
  

  TH1F* h_leadDSV=new TH1F("","",20,-1,1);
  TH1F* h_sublDSV=new TH1F("","",20,-1,1);
  TH1F* h_leadPR=new TH1F("","",20,50,250);
  TH1F* h_sublPR=new TH1F("","",20,50,250);
  TH1F* h_leadPt=new TH1F("","",40,200,1000);
  TH1F* h_sublPt=new TH1F("","",40,200,1000);
  TH1F* h_leadEta=new TH1F("","",50,-2.5,2.5);
  TH1F* h_sublEta=new TH1F("","",50,-2.5,2.5);
  TH1F* h_Mjj=new TH1F("","",32,900,2500);
  TH1F* h_DelEta=new TH1F("","",28,0,1.4);

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
      if(thisJet->Pt()<300)continue;
      if(fabs(thisJet->Eta())>2.4)continue;
      if(!passFatJetLooseID[ij])continue;
      if(fatjetPRmassL2L3Corr[ij]<70 || fatjetPRmassL2L3Corr[ij]>200)continue;
      
      fatjet.push_back(ij);     
    }
    
    if(fatjet.size()<2)continue;

    int fat0=fatjet[0];
    int fat1=fatjet[1];

    h_leadPR->Fill(fatjetPRmass[fat0]);
    h_sublPR->Fill(fatjetPRmass[fat1]);

    TLorentzVector* fatjet0 = (TLorentzVector*)fatjetP4->At(fat0);
    TLorentzVector* fatjet1 = (TLorentzVector*)fatjetP4->At(fat1);

    h_leadPt->Fill(fatjet0->Pt());
    h_sublPt->Fill(fatjet1->Pt());
    h_leadEta->Fill(fatjet0->Eta());
    h_sublEta->Fill(fatjet1->Eta());
    
    for(unsigned int i=0; i<fatjet.size(); i++) {
      for(unsigned int j=0; j<i; j++) {
        int index_that = fatjet[i];
        int index_those = fatjet[j];
        TLorentzVector* thatJet  = (TLorentzVector*)fatjetP4->At(index_that);
        TLorentzVector* thoseJet = (TLorentzVector*)fatjetP4->At(index_those);
        float dEta = fabs(thatJet->Eta() - thoseJet->Eta());
        if(dEta>1.3)continue;

        h_DelEta->Fill(dEta);
        
        Double_t mjj = (*thatJet+*thoseJet).M();
        if(mjj<1000)continue; 
        Mjj.push_back(make_pair(index_that,index_those));

      }
    }

    
    if(Mjj.size()<1)continue;   

    int aa = Mjj[0].second;
    int ee = Mjj[0].first;
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa); 
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);
      
    Float_t mff=(*Jet1+*Jet2).M();
    h_Mjj->Fill(mff);

    int addJetIndex[2]={-1,-1}; 
    for(int ad=0; ad<nAJets; ad++) {
      TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
      if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;} // first add jet to pass the delta r cut
      if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;} // first add jet to pass the delta r cut
    }
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;
    h_leadDSV->Fill(addjet_doubleSV[addJetIndex[0]]);
    h_sublDSV->Fill(addjet_doubleSV[addJetIndex[1]]);   

  } //end of the event loop

  TFile* outfile = new TFile(Form("qcd_%d.root",name),"recreate");
  h_leadDSV->Write("leadDSV");
  h_sublDSV->Write("sublDSV");
  h_leadPt->Write("leadPt");
  h_sublPt->Write("sublPt");
  h_leadPR->Write("leadPR");
  h_sublPR->Write("sublPR");
  h_leadEta->Write("leadEta");
  h_sublEta->Write("sublEta");
  h_Mjj->Write("Mjj");
  h_DelEta->Write("DelEta");
  outfile->Write();

  std::cout << "Events on lead DSV = " << h_leadDSV->Integral() << std::endl;

}
