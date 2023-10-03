#if !defined(__CINT__) || defined(__MAKECINT__)
#include <TROOT.h>
#include <TSystem.h>
#include <TFile.h>
#include <TChain.h>
#include <TTree.h>
#include <TMath.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <vector>
#include <iostream>
#include <iomanip>

#endif

  TCanvas *c1 = new TCanvas("MC_Data","",900,600);

void Pulse_shape(TString infile) {

  gDirectory->DeleteAll();

  TChain *chain = new TChain("ExportTree/HcalTree");
  chain->AddFile(infile);
  
  Int_t           PulseCount;
  Double_t        Charge[5184][10];
  Double_t        Pedestal[5184][10];
  Int_t           IEta[5184];

  chain->SetBranchAddress("PulseCount", &PulseCount);
  chain->SetBranchAddress("Charge", Charge);
  chain->SetBranchAddress("Pedestal", Pedestal);
  chain->SetBranchAddress("IEta", IEta);

  double TS[10];
  double sumQ;
  double ped;
  
  TProfile* PulseShape =  new TProfile("PulseShape", "PulseShape;Time Slice;Charge (fC)", 10,0,10);

  long int nentries = std::min(chain->GetEntries(), (long long int) 100000);

  for (UInt_t i=0; i<nentries; i++) {
    chain->GetEntry(i);

    for (int j = 0; j < (int)PulseCount; j++) {
      if (abs(IEta[j])>27 || ( abs(IEta[j])>14 && abs(IEta[j])<19)) continue;
      if(Charge[j][0]==Charge[j][4] && Charge[j][4]==Charge[j][5]) continue;

      sumQ=0; ped=0;

      for (UInt_t k=0; k<3; k++) { ped+=Charge[j][k]+Pedestal[j][k];}
      for (UInt_t k=3; k<7; k++) { sumQ+=Charge[j][k]+Pedestal[j][k];}
      
      ped/=3;
      sumQ-=4*ped;

      for ( int k=0; k < 10; k++)
         TS[k] = Charge[j][k]+Pedestal[j][k]-ped;
      
      if (sumQ<5 || TS[4]<5 || ped<0 || TS[3]<0 || TS[5]<0 || TS[6]<0) continue;

      for ( int k=0; k < 10; k++)
         PulseShape->Fill(k, TS[k]);
    }
  }
 
  PulseShape->Draw("HIST");
  c1->SaveAs("PulseShape.pdf");
}
