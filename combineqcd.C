#include <vector>
#include <string>
#include <iostream>
#include <TH1D.h>
#include <TMath.h>
#include <TFile.h>
#include <TH1.h>
#include <THnSparse.h>
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include "untuplizer.h"
#include "readSample.h"

void qcd() {

  TFile* file[12];

  //file[0]  = TFile::Open("data_1.root");
  //file[1]  = TFile::Open("data_2.root");
  //file[2]  = TFile::Open("data_3.root");
  //file[3]  = TFile::Open("data_4.root");
  //file[4]  = TFile::Open("data_5.root");
  //file[5]  = TFile::Open("data_6.root");
  file[6]  = TFile::Open("qcd_7.root");
  file[7]  = TFile::Open("qcd_8.root");
  file[8]  = TFile::Open("qcd_9.root");
  file[9]  = TFile::Open("qcd_10.root");
  file[10]  = TFile::Open("qcd_11.root");
  file[11]  = TFile::Open("qcd_12.root");
  
  TH1F* h_leadDSV7  = (TH1F*)(file[6]->Get("leadDSV"));
  TH1F* h_sublDSV7  = (TH1F*)(file[6]->Get("sublDSV"));
  TH1F* h_leadPt7   = (TH1F*)(file[6]->Get("leadPt"));
  TH1F* h_sublPt7   = (TH1F*)(file[6]->Get("sublPt"));
  TH1F* h_leadPR7   = (TH1F*)(file[6]->Get("leadPR"));
  TH1F* h_sublPR7   = (TH1F*)(file[6]->Get("sublPR"));
  TH1F* h_leadEta7  = (TH1F*)(file[6]->Get("leadEta"));
  TH1F* h_sublEta7  = (TH1F*)(file[6]->Get("sublEta"));
  TH1F* h_Mjj7      = (TH1F*)(file[6]->Get("Mjj"));
  TH1F* h_DelEta7   = (TH1F*)(file[6]->Get("DelEta"));

  TH1F* h_leadDSV8  = (TH1F*)(file[7]->Get("leadDSV"));
  TH1F* h_sublDSV8  = (TH1F*)(file[7]->Get("sublDSV"));
  TH1F* h_leadPt8   = (TH1F*)(file[7]->Get("leadPt"));
  TH1F* h_sublPt8   = (TH1F*)(file[7]->Get("sublPt"));
  TH1F* h_leadPR8   = (TH1F*)(file[7]->Get("leadPR"));
  TH1F* h_sublPR8   = (TH1F*)(file[7]->Get("sublPR"));
  TH1F* h_leadEta8  = (TH1F*)(file[7]->Get("leadEta"));
  TH1F* h_sublEta8  = (TH1F*)(file[7]->Get("sublEta"));
  TH1F* h_Mjj8      = (TH1F*)(file[7]->Get("Mjj"));
  TH1F* h_DelEta8   = (TH1F*)(file[7]->Get("DelEta"));

  TH1F* h_leadDSV9  = (TH1F*)(file[8]->Get("leadDSV"));
  TH1F* h_sublDSV9  = (TH1F*)(file[8]->Get("sublDSV"));
  TH1F* h_leadPt9   = (TH1F*)(file[8]->Get("leadPt"));
  TH1F* h_sublPt9   = (TH1F*)(file[8]->Get("sublPt"));
  TH1F* h_leadPR9   = (TH1F*)(file[8]->Get("leadPR"));
  TH1F* h_sublPR9   = (TH1F*)(file[8]->Get("sublPR"));
  TH1F* h_leadEta9  = (TH1F*)(file[8]->Get("leadEta"));
  TH1F* h_sublEta9  = (TH1F*)(file[8]->Get("sublEta"));
  TH1F* h_Mjj9      = (TH1F*)(file[8]->Get("Mjj"));
  TH1F* h_DelEta9   = (TH1F*)(file[8]->Get("DelEta"));

  TH1F* h_leadDSV10  = (TH1F*)(file[9]->Get("leadDSV"));
  TH1F* h_sublDSV10  = (TH1F*)(file[9]->Get("sublDSV"));
  TH1F* h_leadPt10   = (TH1F*)(file[9]->Get("leadPt"));
  TH1F* h_sublPt10   = (TH1F*)(file[9]->Get("sublPt"));
  TH1F* h_leadPR10   = (TH1F*)(file[9]->Get("leadPR"));
  TH1F* h_sublPR10   = (TH1F*)(file[9]->Get("sublPR"));
  TH1F* h_leadEta10  = (TH1F*)(file[9]->Get("leadEta"));
  TH1F* h_sublEta10  = (TH1F*)(file[9]->Get("sublEta"));
  TH1F* h_Mjj10      = (TH1F*)(file[9]->Get("Mjj"));
  TH1F* h_DelEta10   = (TH1F*)(file[9]->Get("DelEta"));

  TH1F* h_leadDSV11  = (TH1F*)(file[10]->Get("leadDSV"));
  TH1F* h_sublDSV11  = (TH1F*)(file[10]->Get("sublDSV"));
  TH1F* h_leadPt11   = (TH1F*)(file[10]->Get("leadPt"));
  TH1F* h_sublPt11   = (TH1F*)(file[10]->Get("sublPt"));
  TH1F* h_leadPR11   = (TH1F*)(file[10]->Get("leadPR"));
  TH1F* h_sublPR11   = (TH1F*)(file[10]->Get("sublPR"));
  TH1F* h_leadEta11  = (TH1F*)(file[10]->Get("leadEta"));
  TH1F* h_sublEta11  = (TH1F*)(file[10]->Get("sublEta"));
  TH1F* h_Mjj11      = (TH1F*)(file[10]->Get("Mjj"));
  TH1F* h_DelEta11   = (TH1F*)(file[10]->Get("DelEta"));

  TH1F* h_leadDSV12  = (TH1F*)(file[11]->Get("leadDSV"));
  TH1F* h_sublDSV12  = (TH1F*)(file[11]->Get("sublDSV"));
  TH1F* h_leadPt12   = (TH1F*)(file[11]->Get("leadPt"));
  TH1F* h_sublPt12   = (TH1F*)(file[11]->Get("sublPt"));
  TH1F* h_leadPR12   = (TH1F*)(file[11]->Get("leadPR"));
  TH1F* h_sublPR12   = (TH1F*)(file[11]->Get("sublPR"));
  TH1F* h_leadEta12  = (TH1F*)(file[11]->Get("leadEta"));
  TH1F* h_sublEta12  = (TH1F*)(file[11]->Get("sublEta"));
  TH1F* h_Mjj12      = (TH1F*)(file[11]->Get("Mjj"));
  TH1F* h_DelEta12   = (TH1F*)(file[11]->Get("DelEta"));


  //-------Lead DSV-------
  Double_t norm7  = 2684.*8654.49315/7910182;
  Double_t norm8  = 2684.*797.3526900/7205073;
  Double_t norm9  = 2684.*79.02553776/3841262;
  Double_t norm10 = 2684.*25.09505908/3984898;
  Double_t norm11 = 2684.*4.707368272/3566431;
  Double_t norm12 = 2684.*1.621316920/3714871;

  TH1F* h_leadDSV = (TH1F*)h_leadDSV7->Clone("h_leadDSV");
  h_leadDSV->Reset();
  
  h_leadDSV->Add(h_leadDSV7,norm7);
  h_leadDSV->Add(h_leadDSV8,norm8);
  h_leadDSV->Add(h_leadDSV9,norm9);
  h_leadDSV->Add(h_leadDSV10,norm10);
  h_leadDSV->Add(h_leadDSV11,norm11);
  h_leadDSV->Add(h_leadDSV12,norm12);

  TH1F* h_sublDSV = (TH1F*)h_sublDSV7->Clone("h_sublDSV");
  h_sublDSV->Reset();

  h_sublDSV->Add(h_sublDSV7,norm7);
  h_sublDSV->Add(h_sublDSV8,norm8);
  h_sublDSV->Add(h_sublDSV9,norm9);
  h_sublDSV->Add(h_sublDSV10,norm10);
  h_sublDSV->Add(h_sublDSV11,norm11);
  h_sublDSV->Add(h_sublDSV12,norm12);

  TH1F* h_leadPt = (TH1F*)h_leadPt7->Clone("h_leadPt");
  h_leadPt->Reset();

  h_leadPt->Add(h_leadPt7,norm7);
  h_leadPt->Add(h_leadPt8,norm8);
  h_leadPt->Add(h_leadPt9,norm9);
  h_leadPt->Add(h_leadPt10,norm10);
  h_leadPt->Add(h_leadPt11,norm11);
  h_leadPt->Add(h_leadPt12,norm12);

  TH1F* h_sublPt = (TH1F*)h_sublPt7->Clone("h_sublPt");
  h_sublPt->Reset();

  h_sublPt->Add(h_sublPt7,norm7);
  h_sublPt->Add(h_sublPt8,norm8);
  h_sublPt->Add(h_sublPt9,norm9);
  h_sublPt->Add(h_sublPt10,norm10);
  h_sublPt->Add(h_sublPt11,norm11);
  h_sublPt->Add(h_sublPt12,norm12);

  TH1F* h_leadPR = (TH1F*)h_leadPR7->Clone("h_leadPR");
  h_leadPR->Reset();

  h_leadPR->Add(h_leadPR7,norm7);
  h_leadPR->Add(h_leadPR8,norm8);
  h_leadPR->Add(h_leadPR9,norm9);
  h_leadPR->Add(h_leadPR10,norm10);
  h_leadPR->Add(h_leadPR11,norm11);
  h_leadPR->Add(h_leadPR12,norm12);

  TH1F* h_sublPR = (TH1F*)h_sublPR7->Clone("h_sublPR");
  h_sublPR->Reset();

  h_sublPR->Add(h_sublPR7,norm7);
  h_sublPR->Add(h_sublPR8,norm8);
  h_sublPR->Add(h_sublPR9,norm9);
  h_sublPR->Add(h_sublPR10,norm10);
  h_sublPR->Add(h_sublPR11,norm11);
  h_sublPR->Add(h_sublPR12,norm12);

  TH1F* h_leadEta = (TH1F*)h_leadEta7->Clone("h_leadEta");
  h_leadEta->Reset();

  h_leadEta->Add(h_leadEta7,norm7);
  h_leadEta->Add(h_leadEta8,norm8);
  h_leadEta->Add(h_leadEta9,norm9);
  h_leadEta->Add(h_leadEta10,norm10);
  h_leadEta->Add(h_leadEta11,norm11);
  h_leadEta->Add(h_leadEta12,norm12);

  TH1F* h_sublEta = (TH1F*)h_sublEta7->Clone("h_sublEta");
  h_sublEta->Reset();

  h_sublEta->Add(h_sublEta7,norm7);
  h_sublEta->Add(h_sublEta8,norm8);
  h_sublEta->Add(h_sublEta9,norm9);
  h_sublEta->Add(h_sublEta10,norm10);
  h_sublEta->Add(h_sublEta11,norm11);
  h_sublEta->Add(h_sublEta12,norm12);

  TH1F* h_Mjj = (TH1F*)h_Mjj7->Clone("h_Mjj");
  h_Mjj->Reset();

  h_Mjj->Add(h_Mjj7,norm7);
  h_Mjj->Add(h_Mjj8,norm8);
  h_Mjj->Add(h_Mjj9,norm9);
  h_Mjj->Add(h_Mjj10,norm10);
  h_Mjj->Add(h_Mjj11,norm11);
  h_Mjj->Add(h_Mjj12,norm12);

  TH1F* h_DelEta = (TH1F*)h_DelEta7->Clone("h_DelEta");
  h_DelEta->Reset();

  h_DelEta->Add(h_DelEta7,norm7);
  h_DelEta->Add(h_DelEta8,norm8);
  h_DelEta->Add(h_DelEta9,norm9);
  h_DelEta->Add(h_DelEta10,norm10);
  h_DelEta->Add(h_DelEta11,norm11);
  h_DelEta->Add(h_DelEta12,norm12);



  TFile* outfile = new TFile("qcd_luminized.root","recreate");
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
}
