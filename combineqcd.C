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
  Double_t ldsv7  = 8654.49315/h_leadDSV7->Integral();
  Double_t ldsv8  = 797.3526900/h_leadDSV8->Integral();
  Double_t ldsv9  = 79.02553776/h_leadDSV9->Integral();
  Double_t ldsv10 = 25.09505908/h_leadDSV10->Integral();
  Double_t ldsv11 = 25.09505908/h_leadDSV11->Integral();
  Double_t ldsv12 = 1.621316920/h_leadDSV12->Integral();

  TH1F* h_leadDSV = (TH1F*)h_leadDSV7->Clone("h_leadDSV");
  h_leadDSV->Reset();
  
  h_leadDSV->Add(h_leadDSV7,ldsv7);
  h_leadDSV->Add(h_leadDSV8,ldsv8);
  h_leadDSV->Add(h_leadDSV9,ldsv9);
  h_leadDSV->Add(h_leadDSV10,ldsv10);
  h_leadDSV->Add(h_leadDSV11,ldsv11);
  h_leadDSV->Add(h_leadDSV12,ldsv12);

  //-------Sublead DSV-------                                                                                                                             
  Double_t sdsv7  = 8654.49315/h_sublDSV7->Integral();
  Double_t sdsv8  = 797.3526900/h_sublDSV8->Integral();
  Double_t sdsv9  = 79.02553776/h_sublDSV9->Integral();
  Double_t sdsv10 = 25.09505908/h_sublDSV10->Integral();
  Double_t sdsv11 = 25.09505908/h_sublDSV11->Integral();
  Double_t sdsv12 = 1.621316920/h_sublDSV12->Integral();

  TH1F* h_sublDSV = (TH1F*)h_sublDSV7->Clone("h_sublDSV");
  h_sublDSV->Reset();

  h_sublDSV->Add(h_sublDSV7,sdsv7);
  h_sublDSV->Add(h_sublDSV8,sdsv8);
  h_sublDSV->Add(h_sublDSV9,sdsv9);
  h_sublDSV->Add(h_sublDSV10,sdsv10);
  h_sublDSV->Add(h_sublDSV11,sdsv11);
  h_sublDSV->Add(h_sublDSV12,sdsv12);

  //-------lead Pt-------                                                                                                                            
  Double_t lpt7  = 8654.49315/h_leadPt7->Integral();
  Double_t lpt8  = 797.3526900/h_leadPt8->Integral();
  Double_t lpt9  = 79.02553776/h_leadPt9->Integral();
  Double_t lpt10 = 25.09505908/h_leadPt10->Integral();
  Double_t lpt11 = 25.09505908/h_leadPt11->Integral();
  Double_t lpt12 = 1.621316920/h_leadPt12->Integral();

  TH1F* h_leadPt = (TH1F*)h_leadPt7->Clone("h_leadPt");
  h_leadPt->Reset();

  h_leadPt->Add(h_leadPt7,lpt7);
  h_leadPt->Add(h_leadPt8,lpt8);
  h_leadPt->Add(h_leadPt9,lpt9);
  h_leadPt->Add(h_leadPt10,lpt10);
  h_leadPt->Add(h_leadPt11,lpt11);
  h_leadPt->Add(h_leadPt12,lpt12);

  //-------Sublead Pt------                                                                                                                            
  Double_t spt7  = 8654.49315/h_sublPt7->Integral();
  Double_t spt8  = 797.3526900/h_sublPt8->Integral();
  Double_t spt9  = 79.02553776/h_sublPt9->Integral();
  Double_t spt10 = 25.09505908/h_sublPt10->Integral();
  Double_t spt11 = 25.09505908/h_sublPt11->Integral();
  Double_t spt12 = 1.621316920/h_sublPt12->Integral();

  TH1F* h_sublPt = (TH1F*)h_sublPt7->Clone("h_sublPt");
  h_sublPt->Reset();

  h_sublPt->Add(h_sublPt7,spt7);
  h_sublPt->Add(h_sublPt8,spt8);
  h_sublPt->Add(h_sublPt9,spt9);
  h_sublPt->Add(h_sublPt10,spt10);
  h_sublPt->Add(h_sublPt11,spt11);
  h_sublPt->Add(h_sublPt12,spt12);

  //-------lead PR-------                                                                                                                             
  Double_t lpr7  = 8654.49315/h_leadPR7->Integral();
  Double_t lpr8  = 797.3526900/h_leadPR8->Integral();
  Double_t lpr9  = 79.02553776/h_leadPR9->Integral();
  Double_t lpr10 = 25.09505908/h_leadPR10->Integral();
  Double_t lpr11 = 25.09505908/h_leadPR11->Integral();
  Double_t lpr12 = 1.621316920/h_leadPR12->Integral();

  TH1F* h_leadPR = (TH1F*)h_leadPR7->Clone("h_leadPR");
  h_leadPR->Reset();

  h_leadPR->Add(h_leadPR7,lpr7);
  h_leadPR->Add(h_leadPR8,lpr8);
  h_leadPR->Add(h_leadPR9,lpr9);
  h_leadPR->Add(h_leadPR10,lpr10);
  h_leadPR->Add(h_leadPR11,lpr11);
  h_leadPR->Add(h_leadPR12,lpr12);

  //-------Sublead PR------                                                                                                                             
  Double_t spr7  = 8654.49315/h_sublPR7->Integral();
  Double_t spr8  = 797.3526900/h_sublPR8->Integral();
  Double_t spr9  = 79.02553776/h_sublPR9->Integral();
  Double_t spr10 = 25.09505908/h_sublPR10->Integral();
  Double_t spr11 = 25.09505908/h_sublPR11->Integral();
  Double_t spr12 = 1.621316920/h_sublPR12->Integral();

  TH1F* h_sublPR = (TH1F*)h_sublPR7->Clone("h_sublPR");
  h_sublPR->Reset();

  h_sublPR->Add(h_sublPR7,spr7);
  h_sublPR->Add(h_sublPR8,spr8);
  h_sublPR->Add(h_sublPR9,spr9);
  h_sublPR->Add(h_sublPR10,spr10);
  h_sublPR->Add(h_sublPR11,spr11);
  h_sublPR->Add(h_sublPR12,spr12);

  //-------lead eta-------                                                                                                                              
  Double_t leta7  = 8654.49315/h_leadEta7->Integral();
  Double_t leta8  = 797.3526900/h_leadEta8->Integral();
  Double_t leta9  = 79.02553776/h_leadEta9->Integral();
  Double_t leta10 = 25.09505908/h_leadEta10->Integral();
  Double_t leta11 = 25.09505908/h_leadEta11->Integral();
  Double_t leta12 = 1.621316920/h_leadEta12->Integral();

  TH1F* h_leadEta = (TH1F*)h_leadEta7->Clone("h_leadEta");
  h_leadEta->Reset();

  h_leadEta->Add(h_leadEta7,leta7);
  h_leadEta->Add(h_leadEta8,leta8);
  h_leadEta->Add(h_leadEta9,leta9);
  h_leadEta->Add(h_leadEta10,leta10);
  h_leadEta->Add(h_leadEta11,leta11);
  h_leadEta->Add(h_leadEta12,leta12);

  //-------lead eta-------
  Double_t seta7  = 8654.49315/h_sublEta7->Integral();
  Double_t seta8  = 797.3526900/h_sublEta8->Integral();
  Double_t seta9  = 79.02553776/h_sublEta9->Integral();
  Double_t seta10 = 25.09505908/h_sublEta10->Integral();
  Double_t seta11 = 25.09505908/h_sublEta11->Integral();
  Double_t seta12 = 1.621316920/h_sublEta12->Integral();

  TH1F* h_sublEta = (TH1F*)h_sublEta7->Clone("h_sublEta");
  h_sublEta->Reset();

  h_sublEta->Add(h_sublEta7,seta7);
  h_sublEta->Add(h_sublEta8,seta8);
  h_sublEta->Add(h_sublEta9,seta9);
  h_sublEta->Add(h_sublEta10,seta10);
  h_sublEta->Add(h_sublEta11,seta11);
  h_sublEta->Add(h_sublEta12,seta12);

  //-------Mjj-------                                                            
  Double_t mjj7  = 8654.49315/h_Mjj7->Integral();
  Double_t mjj8  = 797.3526900/h_Mjj8->Integral();
  Double_t mjj9  = 79.02553776/h_Mjj9->Integral();
  Double_t mjj10 = 25.09505908/h_Mjj10->Integral();
  Double_t mjj11 = 25.09505908/h_Mjj11->Integral();
  Double_t mjj12 = 1.621316920/h_Mjj12->Integral();

  TH1F* h_Mjj = (TH1F*)h_Mjj7->Clone("h_Mjj");
  h_Mjj->Reset();

  h_Mjj->Add(h_Mjj7,mjj7);
  h_Mjj->Add(h_Mjj8,mjj8);
  h_Mjj->Add(h_Mjj9,mjj9);
  h_Mjj->Add(h_Mjj10,mjj10);
  h_Mjj->Add(h_Mjj11,mjj11);
  h_Mjj->Add(h_Mjj12,mjj12);

  //------Delta eta-------                                                                                                                         
  Double_t deta7  = 8654.49315/h_DelEta7->Integral();
  Double_t deta8  = 797.3526900/h_DelEta8->Integral();
  Double_t deta9  = 79.02553776/h_DelEta9->Integral();
  Double_t deta10 = 25.09505908/h_DelEta10->Integral();
  Double_t deta11 = 25.09505908/h_DelEta11->Integral();
  Double_t deta12 = 1.621316920/h_DelEta12->Integral();

  TH1F* h_DelEta = (TH1F*)h_DelEta7->Clone("h_DelEta");
  h_DelEta->Reset();

  h_DelEta->Add(h_DelEta7,deta7);
  h_DelEta->Add(h_DelEta8,deta8);
  h_DelEta->Add(h_DelEta9,deta9);
  h_DelEta->Add(h_DelEta10,deta10);
  h_DelEta->Add(h_DelEta11,deta11);
  h_DelEta->Add(h_DelEta12,deta12);



  TFile* outfile = new TFile("qcd.root","recreate");
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
