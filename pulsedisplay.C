#include <iostream>
#include <vector>

#include "TSystem.h"
#include "TStyle.h"
#include "TH1F.h"
#include "TFile.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TString.h"
#include "TPad.h"
#include "TLine.h"
#include "TLatex.h"
#include "TChain.h"
#include "TError.h"
#include "THStack.h"
#include "TAxis.h"

//#include "utilities.hpp"

using namespace std;

TH1D* InitTH1D(TString Name, TString Title, int Nbins, double XMin, double XMax)
{
  TH1D *h1 = new TH1D(Name, Title, Nbins, XMin, XMax);
  h1->Sumw2();
  return h1;
}

TLegend* drawPulsePlot(TChain *ch, int printIndex, THStack* &st, TH1D* &h1_data, TH1D* &h1_all, TString title, float &max) 
                  //int run, int ls, int evt, int ieta, int iphi, int depth, 
{ 
  float inputTS[10];
  float itPulse[10];
  float p1Pulse[10];
  float n1Pulse[10];
  float p2Pulse[10];
  float n2Pulse[10];
  float p3Pulse[10];
  float n3Pulse[10];
  float n4Pulse[10];
  float totalUCNoise[10];
  float p1Energy=0;
  float n1Energy=0;
  float p2Energy=0;
  float n2Energy=0;
  float p3Energy=0;
  float n3Energy=0;
  float n4Energy=0;
  float mahiEnergy=0;
  float pedEnergy=0;
  float inGain=0;
  float chiSq=0;

  ch->SetBranchAddress("inputTS",    &inputTS);
  ch->SetBranchAddress("itPulse",    &itPulse);
  ch->SetBranchAddress("p1Pulse",    &p1Pulse);
  ch->SetBranchAddress("n1Pulse",    &n1Pulse);
  ch->SetBranchAddress("p2Pulse",    &p2Pulse);
  ch->SetBranchAddress("n2Pulse",    &n2Pulse);
  ch->SetBranchAddress("p3Pulse",    &p3Pulse);
  ch->SetBranchAddress("n3Pulse",    &n3Pulse);
  ch->SetBranchAddress("n4Pulse",    &n4Pulse);
  ch->SetBranchAddress("totalUCNoise",    &totalUCNoise);
  ch->SetBranchAddress("p1Energy",   &p1Energy);
  ch->SetBranchAddress("n1Energy",   &n1Energy);
  ch->SetBranchAddress("p2Energy",   &p2Energy);
  ch->SetBranchAddress("n2Energy",   &n2Energy);
  ch->SetBranchAddress("p3Energy",   &p3Energy);
  ch->SetBranchAddress("n3Energy",   &n3Energy);
  ch->SetBranchAddress("n4Energy",   &n4Energy);
  ch->SetBranchAddress("mahiEnergy", &mahiEnergy);
  ch->SetBranchAddress("pedEnergy",  &pedEnergy);
  ch->SetBranchAddress("inGain",     &inGain);
  ch->SetBranchAddress("chiSq",      &chiSq);
      
  TH1D *h1_it   = InitTH1D("h1_it",   "h1_it",   8, -3.5, 4.5); 
  TH1D *h1_p1   = InitTH1D("h1_p1",   "h1_p1",   8, -3.5, 4.5); 
  TH1D *h1_p2   = InitTH1D("h1_p2",   "h1_p2",   8, -3.5, 4.5); 
  TH1D *h1_p3   = InitTH1D("h1_p3",   "h1_p3",   8, -3.5, 4.5); 
  TH1D *h1_n1   = InitTH1D("h1_n1",   "h1_n1",   8, -3.5, 4.5); 
  TH1D *h1_n2   = InitTH1D("h1_n2",   "h1_n2",   8, -3.5, 4.5); 
  TH1D *h1_n3   = InitTH1D("h1_n3",   "h1_n3",   8, -3.5, 4.5); 
  TH1D *h1_n4   = InitTH1D("h1_n4",   "h1_n4",   8, -3.5, 4.5); 
  TH1D *h1_ped  = InitTH1D("h1_ped",  "h1_ped",  8, -3.5, 4.5); 
  //TH1D *h1_all  = InitTH1D("h1_all",  "h1_all",  8, -3.5, 4.5); 

  ch->GetEntry(printIndex);

  title.ReplaceAll("E= GeV", Form("E= %.1f GeV", mahiEnergy*inGain)); 
  title.ReplaceAll("chi2=", Form("#chi^{2}= %.1f", chiSq)); 

  if(0)
 {
    cout << title << endl; 
    cout << p3Energy << endl; 
    cout << p2Energy << endl; 
    cout << p1Energy << endl; 
    cout << mahiEnergy << endl; 
    cout << n1Energy << endl; 
    cout << n2Energy << endl; 
    cout << n3Energy << endl; 
    cout << n4Energy << endl; 
    cout << pedEnergy << endl; 
  }
  

  for(int its=0; its<8; its++) 
  {
    if(p2Energy<0) p2Energy=0; 
    if(p3Energy<0) p3Energy=0; 
    if(n2Energy<0) n2Energy=0; 
    if(n3Energy<0) n3Energy=0; 
    if(n4Energy<0) n4Energy=0; 

    //h1_data->SetBinContent(its+1, inputTS[its]>0?inputTS[its]*inGain:0);
    h1_data->SetBinContent(its+1, inputTS[its]*inGain);
    h1_it->SetBinContent(its+1, itPulse[its]*mahiEnergy*inGain);
    h1_p1->SetBinContent(its+1, p1Pulse[its]*p1Energy*inGain);
    h1_p2->SetBinContent(its+1, p2Pulse[its]*p2Energy*inGain);
    h1_p3->SetBinContent(its+1, p3Pulse[its]*p3Energy*inGain);
    h1_n1->SetBinContent(its+1, n1Pulse[its]*n1Energy*inGain);
    h1_n2->SetBinContent(its+1, n2Pulse[its]*n2Energy*inGain);
    h1_n3->SetBinContent(its+1, n3Pulse[its]*n3Energy*inGain);
    h1_n4->SetBinContent(its+1, n4Pulse[its]*n4Energy*inGain);
    h1_ped->SetBinContent(its+1, pedEnergy*inGain);
  } 
  
  //  
  if(mahiEnergy>0) h1_it->Scale(mahiEnergy*inGain/h1_it->Integral(1,8));
  if(p1Energy>0) h1_p1->Scale(p1Energy*inGain/h1_p1->Integral(1,8));
  if(p2Energy>0) h1_p2->Scale(p2Energy*inGain/h1_p2->Integral(1,8));
  if(p3Energy>0) h1_p3->Scale(p3Energy*inGain/h1_p3->Integral(1,8));
  if(n1Energy>0) h1_n1->Scale(n1Energy*inGain/h1_n1->Integral(1,8));
  if(n2Energy>0) h1_n2->Scale(n2Energy*inGain/h1_n2->Integral(1,8));
  if(n3Energy>0) h1_n3->Scale(n3Energy*inGain/h1_n3->Integral(1,8));
  if(n4Energy>0) h1_n4->Scale(n4Energy*inGain/h1_n4->Integral(1,8));

  // all 
  for(int its=0; its<8; its++)
  {
    float all = 0;
    all = h1_it->GetBinContent(its+1);  
    all += h1_p1->GetBinContent(its+1);
    all += h1_p2->GetBinContent(its+1);
    all += h1_p3->GetBinContent(its+1);
    all += h1_n1->GetBinContent(its+1);
    all += h1_n2->GetBinContent(its+1);
    all += h1_n3->GetBinContent(its+1);
    all += h1_n4->GetBinContent(its+1);
    all += h1_ped->GetBinContent(its+1);
    h1_all->SetBinContent(its+1, all); 
    h1_all->SetBinError(its+1, TMath::Sqrt(totalUCNoise[its])*inGain); 
    //cout << "noise: " << its <<" " << TMath::Sqrt(totalUCNoise[its])*inGain << endl;
  }

  if(h1_all->GetMaximum()>max) max = h1_all->GetMaximum();

  h1_data->SetLineColor(kBlack); h1_data->SetMarkerColor(kBlack);
  h1_data->SetMarkerStyle(20);   h1_data->SetMarkerSize(1); 
  h1_it->SetLineColor(kBlack);   h1_it->SetFillColor(kYellow-10);
  h1_p1->SetLineColor(kBlack);   h1_p1->SetFillColor(kAzure+1);
  h1_p2->SetLineColor(kBlack);   h1_p2->SetFillColor(kAzure+6);
  h1_p3->SetLineColor(kBlack);   h1_p3->SetFillColor(kAzure);
  h1_n1->SetLineColor(kBlack);   h1_n1->SetFillColor(kPink+1);
  h1_n2->SetLineColor(kBlack);   h1_n2->SetFillColor(kPink+6);
  h1_n3->SetLineColor(kBlack);   h1_n3->SetFillColor(kPink);
  h1_n4->SetLineColor(kBlack);   h1_n4->SetFillColor(kPink-3);
  h1_ped->SetLineColor(kBlack);  h1_ped->SetFillColor(kGreen-10);
  h1_all->SetLineColor(kBlack);  h1_all->SetLineWidth(2);

  // legend 
  TLegend *l1 = new TLegend(0.15, 0.3, 0.25, 0.85);
  l1->SetNColumns(1);
  l1->SetBorderSize(0);
  l1->SetFillColor(0);
  l1->SetFillStyle(0);
  l1->SetTextFont(42);
  l1->SetTextAlign(12);
  l1->SetTextSize(0.06);
  l1->SetFillColor(kWhite);
  l1->SetLineColor(kWhite);
  l1->SetShadowColor(kWhite);
  l1->AddEntry(h1_data, " digi", "p");
  if(h1_ped->Integral(1,8)>0) l1->AddEntry(h1_ped, Form(" base (%.1f GeV)",h1_ped->Integral(1,8)), "f");
  if(h1_p3->Integral(1,8)>0) l1->AddEntry(h1_p3,   Form(" -3BX (%.1f GeV)",h1_p3->Integral(1,8)),  "f");
  if(h1_p2->Integral(1,8)>0) l1->AddEntry(h1_p2,   Form(" -2BX (%.1f GeV)",h1_p2->Integral(1,8)),  "f");
  if(h1_p1->Integral(1,8)>0) l1->AddEntry(h1_p1,   Form(" -1BX (%.1f GeV)",h1_p1->Integral(1,8)),  "f");
  if(h1_it->Integral(1,8)>0) l1->AddEntry(h1_it,   Form("  0BX (%.1f GeV)",h1_it->Integral(1,8)),  "f");
  if(h1_n1->Integral(1,8)>0) l1->AddEntry(h1_n1,   Form(" +1BX (%.1f GeV)",h1_n1->Integral(1,8)),  "f");
  if(h1_n2->Integral(1,8)>0) l1->AddEntry(h1_n2,   Form(" +2BX (%.1f GeV)",h1_n2->Integral(1,8)),  "f");
  if(h1_n3->Integral(1,8)>0) l1->AddEntry(h1_n3,   Form(" +3BX (%.1f GeV)",h1_n3->Integral(1,8)),  "f"); 
  if(h1_n4->Integral(1,8)>0) l1->AddEntry(h1_n4,   Form(" +4BX (%.1f GeV)",h1_n4->Integral(1,8)),  "f");

  if(h1_ped->Integral(1,8)>0) st->Add(h1_ped);
  if(h1_p3->Integral(1,8)>0) st->Add(h1_p3);
  if(h1_p2->Integral(1,8)>0) st->Add(h1_p2);
  if(h1_p1->Integral(1,8)>0) st->Add(h1_p1);
  if(h1_n1->Integral(1,8)>0) st->Add(h1_n1);
  if(h1_n2->Integral(1,8)>0) st->Add(h1_n2);
  if(h1_n3->Integral(1,8)>0) st->Add(h1_n3); 
  if(h1_n4->Integral(1,8)>0) st->Add(h1_n4); 
  if(h1_it->Integral(1,8)>0) st->Add(h1_it);
  st->SetTitle(title); 
  st->SetMaximum(h1_data->GetMaximum()*1.2); 
  st->SetMinimum(h1_data->GetMinimum()>0?0:h1_data->GetMinimum()*1.2); 
  return l1;
}

void pulsedisplay()
{
    gErrorIgnoreLevel=kError+1;
    TChain* ch_8p       = new TChain("HcalTree"); 
    ch_8p->Add("input_pulsedisplay.root");  
   
    int run;
    int ls;
    int evt;
    int ieta;
    int iphi;
    int depth;
    float mahiEnergy;
    float inGain;
    ch_8p->SetBranchAddress("run",        &run);
    ch_8p->SetBranchAddress("ls",         &ls);
    ch_8p->SetBranchAddress("evt",        &evt);
    ch_8p->SetBranchAddress("ieta",       &ieta);
    ch_8p->SetBranchAddress("iphi",       &iphi);
    ch_8p->SetBranchAddress("depth",      &depth);
    ch_8p->SetBranchAddress("mahiEnergy", &mahiEnergy);
    ch_8p->SetBranchAddress("inGain",     &inGain);
 
    vector<int> vec_printIndex;
    vector<TString> vec_title1;
    vector<TString> vec_title2;

    int nentries = ch_8p->GetEntries();
    for(int i=0; i<nentries; i++)
    {
      ch_8p->GetEntry(i);  

      vec_printIndex.push_back(i);
      vec_title1.push_back(Form("Run %i  LS %i Event %i", run, ls, evt));
      vec_title2.push_back(Form("IEta %i, IPhi %i, Depth %i", ieta, iphi, depth));
    }

    for(unsigned int i=0; i<vec_printIndex.size(); i++) 
    {
      THStack *st_8p        = new THStack("st_8p",        Form("8p   %s", vec_title1.at(i).Data()));
      TH1D *h1_8p_data      = InitTH1D("h1_8p_data",      "h1_8p_data",       8, -3.5, 4.5); 
      TH1D *h1_8p_all       = InitTH1D("h1_8p_all",      "h1_8p_all",       8, -3.5, 4.5); 

      float max=-1;
      TLegend *leg_8p   = drawPulsePlot(ch_8p,   vec_printIndex.at(i), st_8p,   h1_8p_data, h1_8p_all,  
                          "", max);

      h1_8p_all->SetMarkerSize(0);
      h1_8p_all->SetFillColor(kBlack);
      h1_8p_all->SetLineColor(kBlack);
      h1_8p_all->SetFillStyle(3354);

      if(h1_8p_data->GetMaximum()>max) max=h1_8p_data->GetMaximum();
 
      //
      TCanvas *c = new TCanvas("c","c",800,400); 
      c->SetBottomMargin(0.15); 
      c->SetRightMargin(0.05); 
      st_8p->Draw("hist"); 
      st_8p->SetMaximum(max*1.2);
      h1_8p_all->Draw("e2 same");
      h1_8p_data->Draw("p same");
      // x 
      st_8p->GetXaxis()->SetLabelSize(0.06);
      st_8p->GetXaxis()->SetLabelOffset(0.01);
      st_8p->GetXaxis()->SetTitle("Time slice");
      st_8p->GetXaxis()->SetTitleSize(0.06);
      // y
      st_8p->GetYaxis()->SetLabelSize(0.06);
      st_8p->GetYaxis()->SetLabelOffset(0.01);
      st_8p->GetYaxis()->SetTitle("E [GeV]");
      st_8p->GetYaxis()->SetTitleSize(0.06); 
      st_8p->GetYaxis()->SetTitleOffset(0.7);
      leg_8p->Draw();
      leg_8p->Draw();

      float textSize = 0.06;
      TLatex *TexEnergyLumi = new TLatex(0.95,0.92,"#font[42]{13 TeV}");
      TexEnergyLumi->SetNDC();
      TexEnergyLumi->SetTextSize(textSize);
      TexEnergyLumi->SetTextAlign (31);
      TexEnergyLumi->SetLineWidth(2);

      TLatex *TexCMS = new TLatex(0.1,0.92,"CMS #font[52]{Preliminary} 2018");
      TexCMS->SetNDC();
      TexCMS->SetTextSize(textSize);
      TexCMS->SetLineWidth(2);
      
      TLatex *TexLabel1 = new TLatex(0.94,0.83,vec_title1.at(i));
      TexLabel1->SetNDC();
      TexLabel1->SetTextSize(textSize-0.01);
      TexLabel1->SetLineWidth(2);
      TexLabel1->SetTextAlign(31);
      TexLabel1->SetTextFont(42);
      
      TLatex *TexLabel2 = new TLatex(0.94,0.83-textSize,vec_title2.at(i));
      TexLabel2->SetNDC();
      TexLabel2->SetTextSize(textSize-0.01);
      TexLabel2->SetLineWidth(2);
      TexLabel2->SetTextAlign(31);
      TexLabel2->SetTextFont(42);
      
      TexEnergyLumi->Draw("same");
      TexCMS->Draw("same");
      TexLabel1->Draw("same");
      TexLabel2->Draw("same");

      c->Print(Form("pulse_run%i_ls%i_evt%i_ieta%i_iphi%i_depth%i.pdf", run, ls, evt, ieta, iphi, depth));
    }

}
