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
void dataMC(std::string inputFile, std::string name) {

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
  TH1F* h_Msubt=new TH1F("","",36,700,2500);
  TH1F* h_Mjj=new TH1F("","",32,900,2500);
  TH1F* h_DelEta=new TH1F("","",28,0,1.4);

  
  h_leadDSV->Sumw2();
  h_sublDSV->Sumw2();
  h_leadPR->Sumw2();
  h_sublPR->Sumw2();
  h_leadPt->Sumw2();
  h_sublPt->Sumw2();
  h_leadEta->Sumw2();
  h_sublEta->Sumw2();
  h_Msubt->Sumw2();
  h_Mjj->Sumw2();
  h_DelEta->Sumw2();
  

  Long64_t nPass[20] = {0};

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
    Float_t*  fatjetTau1 = data.GetPtrFloat("FATjetTau1");
    Float_t*  fatjetTau2 = data.GetPtrFloat("FATjetTau2");
    Float_t mcWeight  = data.GetFloat("mcWeight");
    //Float_t*  fatjet_doubleSV = data.GetPtrFloat("FATjet_DoubleSV");
    
    int nADDJet         = data.GetInt("ADDnJet");
    const int nAJets=nADDJet;
    TClonesArray* addjetP4 = (TClonesArray*) data.GetPtrTObject("ADDjetP4");
    Float_t*  addjet_doubleSV = data.GetPtrFloat("ADDjet_DoubleSV");
    
    //nPass[0] =data.GetEntriesFast();

    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));
    const Int_t nsize = data.GetPtrStringSize();

    bool passTrigger=false;
    for(int it=0; it< nsize; it++)
      {
	std::string thisTrig= trigName[it];
	bool results = trigResult[it];

	// std::cout << thisTrig << " : " << results << std::endl;
	
	if( (thisTrig.find("HLT_PFHT800")!= std::string::npos && results==1)
	    )
	  {
	    passTrigger=true;
	    break;
	  }


      }


    if(!passTrigger)continue;
    nPass[1]++;


    //3. has a good vertex
    Int_t nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;
    nPass[2]++;

    vector<int> fatjet;
    vector<pair<int,int>> Mjj;
    for(int ij=0; ij<nFJets; ij++) {
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
      if(thisJet->Pt()<200)continue;
      if(fabs(thisJet->Eta())>2.4)continue;
      if(!passFatJetLooseID[ij])continue;
      if(fatjetPRmassL2L3Corr[ij]<70 || fatjetPRmassL2L3Corr[ij]>200)continue;

      Double_t tau21 = (fatjetTau2[ij]/fatjetTau1[ij]);
      if(tau21>0.6)continue;
      
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

	h_DelEta->Fill(dEta,mcWeight);
	h_Mjj->Fill(mjj,mcWeight);
        Mjj.push_back(make_pair(index_that,index_those));

      }
    }

    
    if(Mjj.size()<1)continue;   

    int aa = Mjj[0].second;
    int ee = Mjj[0].first;
    TLorentzVector* Jet1 = (TLorentzVector*)fatjetP4->At(aa); 
    TLorentzVector* Jet2 = (TLorentzVector*)fatjetP4->At(ee);
  
    
    Float_t mff=(*Jet1+*Jet2).M();
    Float_t msubt = mff-(Jet1->M()-125)-(Jet2->M()-125);
    if(msubt<800)continue;


    h_Msubt->Fill(msubt,mcWeight);
    h_leadPR->Fill(fatjetPRmass[aa],mcWeight);
    h_sublPR->Fill(fatjetPRmass[ee],mcWeight);
    h_leadPt->Fill(Jet1->Pt(),mcWeight);
    h_sublPt->Fill(Jet2->Pt(),mcWeight);
    h_leadEta->Fill(Jet1->Eta(),mcWeight);
    h_sublEta->Fill(Jet2->Eta(),mcWeight);
    //h_leadDSV->Fill(fatjet_doubleSV[aa]); //,mcWeight);
    //h_sublDSV->Fill(fatjet_doubleSV[ee]); //,mcWeight);

    
    int addJetIndex[2]={-1,-1}; 
    for(int ad=0; ad<nAJets; ad++) {
      TLorentzVector* Jet3 = (TLorentzVector*)addjetP4->At(ad);
      if(Jet1->DeltaR(*Jet3)<0.1 && addJetIndex[0] < 0) { addJetIndex[0]=ad;} // first add jet to pass the delta r cut
      if(Jet2->DeltaR(*Jet3)<0.1 && addJetIndex[1] < 0) { addJetIndex[1]=ad;} // first add jet to pass the delta r cut
    }
    if(addJetIndex[0]<0 || addJetIndex[1]<0)continue;
    h_leadDSV->Fill(addjet_doubleSV[addJetIndex[0]],mcWeight);
    h_sublDSV->Fill(addjet_doubleSV[addJetIndex[1]],mcWeight);  
     

  } //end of the event loop

  TFile* outfile = new TFile(Form("%s_3.0.root",name.data()),"recreate");
  h_leadDSV->Write("leadDSV");
  h_sublDSV->Write("sublDSV");
  h_leadPt->Write("leadPt");
  h_sublPt->Write("sublPt");
  h_leadPR->Write("leadPR");
  h_sublPR->Write("sublPR");
  h_leadEta->Write("leadEta");
  h_sublEta->Write("sublEta");
  h_Msubt->Write("Msubt");
  h_Mjj->Write("Mjj");
  h_DelEta->Write("DelEta");
  outfile->Write();

  std::cout << "Events on lead DSV = " << h_leadDSV->Integral() << std::endl;
  //std::cout << "Total number of Events = " << nPass[0] << std::endl;
  std::cout << nPass[1] << " " << nPass[2] << std::endl;
}
