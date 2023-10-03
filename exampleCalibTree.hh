//////////////////////////////////////////////////////////
// Header file
// with CalibTree class 
// for isotrack calibration
//////////////////////////////////////////////////////////

#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TProfile.h>
#include <TLegend.h>
#include <TString.h>
#include <TF1.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>


//**********************************************************
// Class with TTree containing parameters of selected events
//**********************************************************
class CalibTree {
public :
  TChain          *fChain;   //!pointer to the analyzed TTree
  //TChain          *inChain;   //!pointer to the analyzed TChain
  Int_t           fCurrent; //!current Tree number in a TChain

  // Declaration of leaf types
  Int_t           t_Run;
  Int_t           t_Event;
  Int_t           t_nVtx;
  Int_t           t_nTrk;
  Double_t        t_EventWeight;
  Double_t        t_p;
  Double_t        t_pt;
  Int_t           t_ieta;
  Double_t        t_phi;
  Double_t        t_eMipDR;
  Double_t        t_eHcal;
  Double_t        t_eHcal10;
  Double_t        t_eHcal30;
  Double_t        t_hmaxNearP;
  Bool_t          t_selectTk;
  Bool_t          t_qltyMissFlag;
  Bool_t          t_qltyPVFlag;
  std::vector<unsigned int> *t_DetIds;
  std::vector<double>  *t_HitEnergies;
  
  // List of branches
  TBranch        *b_t_Run;   //!
  TBranch        *b_t_Event;   //!
  TBranch        *b_t_nVtx;
  TBranch        *b_t_nTrk;
  TBranch        *b_t_EventWeight;   //!
  TBranch        *b_t_p;   //!
  TBranch        *b_t_pt;   //!
  TBranch        *b_t_ieta;   //!
  TBranch        *b_t_phi;   //!
  TBranch        *b_t_eMipDR;   //!
  TBranch        *b_t_eHcal;   //!
  TBranch        *b_t_eHcal10;   //!
  TBranch        *b_t_eHcal30;   //!
  TBranch        *b_t_hmaxNearP;   //!
  TBranch        *b_t_selectTk;   //!
  TBranch        *b_t_qltyMissFlag;   //!
  TBranch        *b_t_qltyPVFlag;   //!
  TBranch        *b_t_DetIds;   //!
  TBranch        *b_t_HitEnergies;   //!

  //--- constructor & destructor
  //CalibTree(TTree *tree=0);
  CalibTree(TChain *tree,
	    double min_enrHcal, double min_pt,
	    double lim_mipEcal, double lim_charIso,
	    double min_trackMom, double max_trackMom);
  virtual ~CalibTree();
  
  //--- functions
  virtual Int_t      GetEntry(Long64_t entry);
  virtual Long64_t   LoadTree(Long64_t entry);
  virtual void       Init(TChain *tree);
  virtual Int_t      LoopOverEvents(int, bool, double);
  Bool_t             goodTrack();
  Bool_t             openOutputRootFile(std::string);

  //--- variables for statistics
  int maxNumOfTracksForIeta;
  std::map<unsigned int, int> subDetector_trk;
  std::map<unsigned int, int> subDetector_final;
  std::map<unsigned int, int> nTrks;
  std::map<unsigned int, int> nSubdetInEvent;
  std::map<unsigned int, int> nPhiMergedInEvent;

  //--- variables for selection
  double minEnrHcal;
  double minTrackPt;
  double minTrackMom;
  double maxTrackMom;
  double limMipEcal;
  double limCharIso;
  
  //--- variables for plotting
  TFile *foutRootFile;
  TH1D* h_nvtx;
  TProfile* h_cluster; // number of cells in cluster
  TH1F* e2pALL;
  TH1F* e2pHB;
  TH1F* e2pTR;
  TH1F* e2pHE;
  };

//**********************************************************
// CalibTree constructor
//**********************************************************
//CalibTree::CalibTree(TTree *tree) : fChain(0) {
CalibTree::CalibTree(TChain *tree,
		     double min_enrHcal,
		     double min_pt,
		     double lim_mipEcal,
		     double lim_charIso,
		     double min_trackMom,
		     double max_trackMom )
{ //: fChain(0) {
  // if parameter tree is not specified (or zero), connect the file
  // used to generate this class and read the Tree.
  if (tree == 0) {
    TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("output.root");
    if (!f || !f->IsOpen()) {
      f = new TFile("output.root");
    }
    TDirectory * dir = (TDirectory*)f->Get("IsoTrackCalibration");
    dir->GetObject("CalibTree",tree);
  }
  Init(tree);

  maxNumOfTracksForIeta = 0;
  // initialization of maps
  subDetector_trk.clear();
  subDetector_final.clear();
  nTrks.clear();
  nSubdetInEvent.clear();
  nPhiMergedInEvent.clear();
  
  // selection
  minEnrHcal = min_enrHcal;
  minTrackPt = min_pt;
  minTrackMom = min_trackMom;
  maxTrackMom = max_trackMom;
  limMipEcal = lim_mipEcal;
  limCharIso = lim_charIso;
}

//**********************************************************
// CalibTree destructor
//**********************************************************
CalibTree::~CalibTree() {

    foutRootFile->cd();
    foutRootFile->Write();
    foutRootFile->Close();

  if (!fChain) return;
  delete fChain->GetCurrentFile();
}

//**********************************************************
// Get entry function
//**********************************************************
Int_t CalibTree::GetEntry(Long64_t entry) {
  // Read contents of entry.
  if (!fChain) return 0;
  return fChain->GetEntry(entry);
}

//**********************************************************
// Load tree function
//**********************************************************
Long64_t CalibTree::LoadTree(Long64_t entry) {
  // Set the environment to read one entry
  if (!fChain) return -5;
  Long64_t centry = fChain->LoadTree(entry);
  if (centry < 0) return centry;
  if (fChain->GetTreeNumber() != fCurrent) {
    fCurrent = fChain->GetTreeNumber();
  }
  return centry;
}

//**********************************************************
// Initialisation of TTree
//**********************************************************
void CalibTree::Init(TChain *tree) {
  // Set object pointer
  t_DetIds = 0;
  t_HitEnergies = 0;
  // Set branch addresses and branch pointers
  if (!tree) return;
  fChain = tree;
  fCurrent = -1;
  fChain->SetMakeClass(1);
  
  fChain->SetBranchAddress("t_Run", &t_Run, &b_t_Run);
  fChain->SetBranchAddress("t_Event", &t_Event, &b_t_Event);
  fChain->SetBranchAddress("t_nVtx", &t_nVtx, &b_t_nVtx);
  fChain->SetBranchAddress("t_nTrk", &t_nTrk, &b_t_nTrk);
  fChain->SetBranchAddress("t_EventWeight", &t_EventWeight, &b_t_EventWeight);
  fChain->SetBranchAddress("t_p", &t_p, &b_t_p);
  fChain->SetBranchAddress("t_pt", &t_pt, &b_t_pt);
  fChain->SetBranchAddress("t_ieta", &t_ieta, &b_t_ieta);
  fChain->SetBranchAddress("t_phi", &t_phi, &b_t_phi);
  fChain->SetBranchAddress("t_eMipDR", &t_eMipDR, &b_t_eMipDR);
  fChain->SetBranchAddress("t_eHcal", &t_eHcal, &b_t_eHcal);
  fChain->SetBranchAddress("t_eHcal10", &t_eHcal10, &b_t_eHcal10);
  fChain->SetBranchAddress("t_eHcal30", &t_eHcal30, &b_t_eHcal30);
  fChain->SetBranchAddress("t_hmaxNearP", &t_hmaxNearP, &b_t_hmaxNearP);
  fChain->SetBranchAddress("t_selectTk", &t_selectTk, &b_t_selectTk);
  fChain->SetBranchAddress("t_qltyMissFlag", &t_qltyMissFlag, &b_t_qltyMissFlag);
  fChain->SetBranchAddress("t_qltyPVFlag", &t_qltyPVFlag, &b_t_qltyPVFlag);
  fChain->SetBranchAddress("t_DetIds", &t_DetIds, &b_t_DetIds);
  fChain->SetBranchAddress("t_HitEnergies", &t_HitEnergies, &b_t_HitEnergies);
}
//**********************************************************
// Open file and book histograms
//**********************************************************
bool CalibTree::openOutputRootFile(std::string fname)
{
  bool decision = false;
  
  foutRootFile = new TFile(fname.c_str(), "RECREATE");
  if ( foutRootFile != NULL ) decision = true;  
  foutRootFile->cd();

  return decision;
}
//**********************************************************


