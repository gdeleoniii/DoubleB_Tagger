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
#include "setNCUStyle.C"
 
using namespace std;
void roc_curve(std::string inputFile) {
  
  //get TTree from file ...
  TreeReader data(inputFile.data());


  TH2F *h2_doubleSV   = new TH2F("","",20,-1,1,20,-1,1);
  TH2F *h2_fatjetCSV = new TH2F("","",20,0,1,20,0,1);

  const int nDim=4;
  int nbins[nDim]={10,10,10,10};
  double xmin[nDim]={0,0,0,0};
  double xmax[nDim]={1,1,1,1};
  THnSparseD* hs = new THnSparseD("hs", "hs", nDim, nbins, xmin, xmax);

  Long64_t DENOM = 0;
  Long64_t DBTnum[21] = {0};
  Long64_t FATnum[21] = {0};
  Long64_t SUBnum[21] = {0};

  Long64_t nExact4[10]={0};
  Long64_t nExact3[10]={0};
  Long64_t nExact2[10]={0};

  //Event loop
  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){
    
    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());
    
    data.GetEntry(jEntry);
    
    int nFATJet         = data.GetInt("FATnJet");
    const int nFJets=nFATJet;
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    Float_t*  fatjetCSV    = data.GetPtrFloat("FATjetCSV");
    Float_t*  fatjetCISVV2 = data.GetPtrFloat("FATjetCISVV2");
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
      if(thisJet->Pt()<300 || thisJet->Pt()400)continue; //this part will be changed for different pt cut
      if(fabs(thisJet->Eta())>2.4)continue;
      if(!passFatJetLooseID[ij])continue;
      if(fatjetPRmassL2L3Corr[ij]<70 || fatjetPRmassL2L3Corr[ij]>200)continue;
      
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

    int aa = Mjj[0].second;
    int ee = Mjj[0].first;
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa); 
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);
    h2_fatjetCSV->Fill(fatjetCISVV2[aa],fatjetCISVV2[ee]);
   
    int addJetIndex[2]={-1,-1}; 
    for(int ad=0; ad<nAJets; ad++) {
      TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
      if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;} // first add jet to pass the delta r cut
      if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;} // first add jet to pass the delta r cut
    }
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;
    h2_doubleSV->Fill(addjet_doubleSV[ addJetIndex[0]], addjet_doubleSV[addJetIndex[1]]);

    Long64_t nSubJetPass[10]={0};

    int fatjetIndex[2]={Mjj[0].second, Mjj[0].first};
    double subjetCSV[4]={-1,-1,-1,-1};

    //---------------per jet counting-----------------\
    
    Long64_t nDBTPass[21] = {0};
    Long64_t nFATPass[21] = {0};
    DENOM += 2;

    Float_t add = -1;
    Float_t fat = 0;
    for(int i=0;i<21;i++) {
     for(int j=0;j<2;j++) {
       int ajet = addJetIndex[j];
       int fjet = fatjetIndex[j];
        if(addjet_doubleSV[ajet]>add)nDBTPass[i]++;
	if(fatjetCISVV2[fjet]>fat)nFATPass[i]++;
     }
    
    DBTnum[i] += nDBTPass[i];
    FATnum[i] += nFATPass[i];
    add += 0.1;
    fat += 0.05;
    }

    Long64_t nSUB1Pass[21] = {0};
    Long64_t nSUB2Pass[21] = {0};

    Float_t sub = 0;
    int fjet1 = fatjetIndex[0];
    int fjet2 = fatjetIndex[1];
    for(int i=0;i<21;i++) {
      for(int sub1=0;sub1<FATnSubSDJet[fjet1];sub1++) {
	if(FATsubjetSDCSV[fjet1][sub1]>sub)nSUB1Pass[i]++;
      }
      for(int sub2=0;sub2<FATnSubSDJet[fjet2];sub2++) {
	if(FATsubjetSDCSV[fjet2][sub2]>sub)nSUB2Pass[i]++;
      }

      if(nSUB1Pass[i] == 1) nSUB1Pass[i] = 0;
      if(nSUB2Pass[i] == 1) nSUB2Pass[i] = 0;
      SUBnum[i] += (nSUB1Pass[i]/2) + (nSUB2Pass[i]/2);
      sub +=0.05;
    }
 
    //------------------------------------------------



      for(int i=0; i<2; i++)
	{
	  int ijet = fatjetIndex[i];
	  for(int isub=0; isub < FATnSubSDJet[ijet]; isub++)
	    {
	      int vectorID = i*2 + isub;
	      subjetCSV[vectorID]=    FATsubjetSDCSV[ijet][isub];
	      if(FATsubjetSDCSV[ijet][isub]>n)nSubJetPass[q]++;
	      
	    }
	}
      hs->Fill(subjetCSV);


  } //end of the event loop

  std::cout << "Denominator= "<< DENOM << std::endl;
  std::cout << "DBT Generated Events= " << h2_doubleSV->Integral(1,20,1,20) << std::endl;
  for(int i=0;i<21;i++) {
    std::cout <<  DBTnum[i] << std::endl;
  }

  std::cout << "FAT Generated Events= " << h2_fatjetCSV->Integral(1,20,1,20) << std::endl;
  for(int i=0;i<21;i++) {
    std::cout <<  FATnum[i] << std::endl;
  }

  std::cout << "SUB Generated Events= " << hs->GetEntries() << std::endl;
  for(int i=0;i<21;i++) {
    std::cout <<  SUBnum[i] << std::endl;
  }
  

}
