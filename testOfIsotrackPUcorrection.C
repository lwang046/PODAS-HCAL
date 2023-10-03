//////////////////////////////////////////////////////////
//---- Example ROOT macro ----------------
// Correction for pileup for isolated charged hadrons
// Requires header file exampleCalibTree.hh
// CalibTree class contains ROOT-tree
// generated with IsoTrackCalibration plugin
////////////////////////////////////////////////////////////

#include <TSystem.h>
#include <TStyle.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TString.h>
#include <TF1.h>
#include <TMath.h>

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <map>
#include <utility>

//**********************************************************
// Constants
//**********************************************************
const int NVERTEX_MIN = 1;   // minimum number of PV in event
const int NVERTEX_MAX = 100; // maximum number of PV in event

const unsigned int MAXNUM_SUBDET = 250;
const int FIRST_IETA_TR = 15;
const int FIRST_IETA_HE = 18;
const unsigned int N_DEPTHS = 8;

// detID format for CMSSW_8_X and +
const unsigned int PHI_MASK       = 0x3FF;
const unsigned int ETA_OFFSET     = 10;
const unsigned int ETA_MASK       = 0x1FF;
const unsigned int ZSIDE_MASK     = 0x80000;
const unsigned int DEPTH_OFFSET   = 20;
const unsigned int DEPTH_MASK     = 0xF;
const unsigned int DEPTH_SET      = 0xF00000;

//--------------------------------------------
// merging depths

const unsigned int MERGE_PHI_AND_DEPTHS = 1;
const unsigned int MASK(0xFFC00); //merge phi and depth
/*
const unsigned int MERGE_PHI_AND_DEPTHS = 0;
const unsigned int MASK(0xFFFC00); //merge phi
*/
//-------------------------------------------

// individual ieta rings
const unsigned int MASK2(0); // no second mask
const int N_ETA_RINGS_PER_BIN = 1;
/*
// twin (even+odd) ieta rings
const unsigned int MASK2(0x80);
const int N_ETA_RINGS_PER_BIN = 2;

// 4-fold ieta rings
const unsigned int MASK2(0x180);
const int N_ETA_RINGS_PER_BIN = 4;
*/
const int MAX_ONESIDE_ETA_RINGS = 30;
const int HALF_NUM_ETA_BINS =
  (MAX_ONESIDE_ETA_RINGS + 1*(N_ETA_RINGS_PER_BIN>1))/N_ETA_RINGS_PER_BIN;
const int NUM_ETA_BINS = 2*HALF_NUM_ETA_BINS + 1;

const int MIN_N_ENTRIES_FOR_FIT = 100;

//--- tuned with MC 2018 using pion samples with and w/o pileup
//--- the description of the method can be found in DN-2016/029
//--- parametrisation is modified compared to the Note 
const double DELTA_CUT = 0.03;
const double CONST_COR_COEF[4]  = { 0.973, 0.998,  0.992,  0.965 };
const double LINEAR_COR_COEF[4] = { 0,    -0.318, -0.261, -0.406 };
const double SQUARE_COR_COEF[4] = { 0,     0,      0.047,  0.089 };
const int PU_IETA_1 = 7;
const int PU_IETA_2 = 16;
const int PU_IETA_3 = 25;

//--- restrictions on response range to avoid outliers
const double MIN_RESPONSE_HIST = 0.0;
const double MAX_RESPONSE_HIST = 3.0;
const int NBIN_RESPONSE_HIST = 1200;
const double RMSRANGE4FIRSTFIT = 1.0;

//**********************************************************
// Polynomial parametrization
// of correction for pileup
//**********************************************************
double corPU(int jeta, double d2p)
{
  // jeta - absolute value of ieta
  unsigned icor = unsigned(jeta >= PU_IETA_1)
    + unsigned(jeta >= PU_IETA_2)
    + unsigned(jeta >= PU_IETA_3);
    
  double cor = CONST_COR_COEF[icor] + LINEAR_COR_COEF[icor]*d2p + SQUARE_COR_COEF[icor]*d2p*d2p;

  return cor;
}
//**********************************************************
// Header with CalibTree class definition
//**********************************************************

#include "exampleCalibTree.hh"

//**********************************************************
// Description of function to run iteration
//**********************************************************
int calculateResponse(const char *inFileDir = "./data",                 // folder with input files
		      const char *inFileNamePrefix = "Commissioning_2023EraD",    // input file prefix
		      const int firstInputFileEnum = 1,                 // numbers of input files (see input folder)
		      const int lastInputFileEnum = 8,
		      const int maxTrackIeta = 25,                     // max abs(ieta), up to which plots will be produced 
		      const bool correctForPU = true,                  // turn on and off correction for pileup
		      const double fitRangeInRMS = 2.0,                // fit range for the second fit (in 2-step fit)
		      const double limitForChargeIsolation = 10.0,     // constraint on charge isolation in GeV
		      const char *outFilePrefix = "example",
		      const double minTrackMomentum = 40.0,            // in GeV
		      const double maxTrackMomentum = 60.0,            // in GeV
		      const double limitForMipInEcal = 1.0,            // in GeV
		      const double minHcalEnergy = 10.0,               // in GeV (to remove mip-like events)
		      const double minPt = 7.0,                        // in GeV
		      const char *treeName = "hcalIsoTrkAnalyzer/CalibTree"
		      )
{
  
  std::cout << "Response to isolated charged hadrons..." << std::endl;

  char corprefix[20] = "_noCor";
  if ( correctForPU ) sprintf(corprefix,"_cor%02d",
					 int(100*DELTA_CUT)
					 );				    
  char fnameInput[120];
  char fnameOutRoot[120];
  char tname[100];
  
  sprintf(tname, "%s", treeName );
  TChain tree(tname);

  //--- combine tree from several enumerated files with the same prefix
  //    or one file w/o number (firstInputFileEnum = lastInputFileEnum < 0 )

  for ( int ik = firstInputFileEnum; ik <= lastInputFileEnum; ik++ ) {
    if ( ik < 0 ) 
      sprintf(fnameInput, "%s/%s.root", inFileDir, inFileNamePrefix);
    else if (ik < 10 )
      sprintf(fnameInput, "%s/%s_%1d.root", inFileDir, inFileNamePrefix, ik);
    else if (ik < 100 )
      sprintf(fnameInput, "%s/%s_%2d.root", inFileDir, inFileNamePrefix, ik);
    else if (ik < 1000 )
      sprintf(fnameInput, "%s/%s_%3d.root", inFileDir, inFileNamePrefix, ik);
    else
      sprintf(fnameInput, "%s/%s_%4d.root", inFileDir, inFileNamePrefix, ik);

    if ( !gSystem->Which("./", fnameInput ) ) { // check file availability
      std::cout << "File " << fnameInput << " doesn't exist." << std::endl;
    }
    else {
      tree.Add(fnameInput);
      std::cout << "Add tree from " << fnameInput 
	        << "   total number of entries (tracks): "
		<< tree.GetEntries() << std::endl;
    }
  }
  if ( tree.GetEntries() == 0 ) {
    std:: cout << "Tree is empty." << std::endl;
    return -2;
  }

  //--- Initialize tree
  CalibTree t(&tree,
	      minHcalEnergy, minPt,
	      limitForMipInEcal, limitForChargeIsolation,
	      minTrackMomentum, maxTrackMomentum);
  
  char isoPrefix[14];
  sprintf(isoPrefix, "const%02d",
	  int(t.limCharIso)
	  );
    
  //--- Define output file
  sprintf(fnameOutRoot,
	  "%s%s_%s_%s_p%02d-%02d_pt%02d_eh%02d_ee%1d_rings%1d_%3.1frms.root",
	  outFilePrefix, corprefix,
	  inFileNamePrefix,
	  isoPrefix,
	  int(minTrackMomentum), int(maxTrackMomentum), int(minPt),
	  int(minHcalEnergy), int(limitForMipInEcal),
	  N_ETA_RINGS_PER_BIN,
	  fitRangeInRMS
	  );

  if ( !t.openOutputRootFile(fnameOutRoot) ) {
    std::cout << "Problems with booking output file " << fnameOutRoot << std::endl;
    return -1;
  }
  std::cout << "Correction for PU: ";
  if ( correctForPU ) {
    std::cout << " applied for delta > " << DELTA_CUT << std::endl;
    std::cout << " [ 1" << " ; " << PU_IETA_1-1 << " ]: "
              << CONST_COR_COEF[0]
              << " ; " << LINEAR_COR_COEF[0]
              << " ; " << SQUARE_COR_COEF[0]
              << std::endl;
    std::cout << " [ " << PU_IETA_1 << " ; " << PU_IETA_2-1 << " ]: "
              << CONST_COR_COEF[1]
              << " ; " << LINEAR_COR_COEF[1]
              << " ; " << SQUARE_COR_COEF[1]
              << std::endl;
    std::cout << " [ " << PU_IETA_2 << " ; " << PU_IETA_3-1 << " ]: "
              << CONST_COR_COEF[2]
              << " ; " << LINEAR_COR_COEF[2]
              << " ; " << SQUARE_COR_COEF[2]
              << std::endl;
    std::cout << " [ " << PU_IETA_3 << " ; 26 ]: "
              << CONST_COR_COEF[3]
              << " ; " << LINEAR_COR_COEF[3]
              << " ; " << SQUARE_COR_COEF[3]
              << std::endl;
  }
  else
    std::cout << " no " << std::endl;

    //--- Prepare initial histograms and count good track
  Int_t nEventsWithGoodTrack = t.LoopOverEvents(maxTrackIeta, correctForPU, fitRangeInRMS);
    std::cout << "Number of events with good track = "
	      << nEventsWithGoodTrack << std::endl;
    std::cout << "Plots saved in " << fnameOutRoot << std::endl;

  return 0;
}

//**********************************************************
// Initial loop over events in the tree
//**********************************************************
Int_t CalibTree::LoopOverEvents(int maxIeta,
				bool applyCorPU,
				double fitRange
				)
{
  char name[100];
  double maxRespForGoodTrack(0);
  double minRespForGoodTrack(1000);
  
  TH1F* e2p[NUM_ETA_BINS];
  int ntrk_ieta[NUM_ETA_BINS];
  for ( int j = 0; j < NUM_ETA_BINS; j++ ) {
    ntrk_ieta[j] = 0;
  }
  
  char scorr[80] = "correction for PU";
  char sxlabel[80] ="E^{cor}_{hcal}/(p_{track} - E_{ecal})"; 
  if ( !applyCorPU ) {
    sprintf(scorr,"no correction for PU");
    sprintf(sxlabel,"E_{hcal}/(p_{track} - E_{ecal})");
  }
  
  TF1* f1 = new TF1("f1","gaus", MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);

  h_nvtx = new TH1D("h_nvtx","Number of vertices in selected events",100,0,100);
  h_nvtx->GetXaxis()->SetTitle("N_{vtx}");
  h_cluster = new TProfile("h_cluster","Number of subdetectors in cluster",
			   2*MAX_ONESIDE_ETA_RINGS,
			   -MAX_ONESIDE_ETA_RINGS, MAX_ONESIDE_ETA_RINGS); 
  h_cluster->GetXaxis()->SetTitle("i#eta of track");
  h_cluster->GetYaxis()->SetTitle("<N_{subdet}>");

  //--------- initialize histograms for response -----------------------------------------
  sprintf(name,"HB+HE: %s", scorr);
  e2pALL = new TH1F("e2pALL", name,
		      NBIN_RESPONSE_HIST, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pALL->Sumw2();
  e2pALL->GetXaxis()->SetTitle(sxlabel);
  sprintf(name,"Events/%5.3f", e2pALL->GetBinWidth(1));
  e2pALL->GetYaxis()->SetTitle(name);
  
  sprintf(name,"Hcal Barrel: %s", scorr);
  e2pHB = new TH1F("e2pHB", name,
			NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pHB->Sumw2();
  e2pHB->GetXaxis()->SetTitle(sxlabel);
  sprintf(name,"Events/%5.3f", e2pHB->GetBinWidth(1));
  e2pHB->GetYaxis()->SetTitle(name);

  sprintf(name,"Transition Region: %s", scorr);
  e2pTR = new TH1F("e2pTR", name,
			NBIN_RESPONSE_HIST/10, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pTR->Sumw2();
  e2pTR->GetXaxis()->SetTitle(sxlabel);
  sprintf(name,"Events/%5.3f", e2pTR->GetBinWidth(1));
  e2pTR->GetYaxis()->SetTitle(name);

  sprintf(name,"Hcal Endcap: %s", scorr);
  e2pHE = new TH1F("e2pHE", name,
			NBIN_RESPONSE_HIST/2, MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
  e2pHE->Sumw2();
  e2pHE->GetXaxis()->SetTitle(sxlabel);
  sprintf(name,"Events/%5.3f", e2pHE->GetBinWidth(1));
  e2pHE->GetYaxis()->SetTitle(name);

//--- initialize chain ----------------------------------------
  if (fChain == 0) return 0;
  Long64_t nentries = fChain->GetEntriesFast();
  Long64_t nb = 0;
  
  int nSelectedEvents(0);

// ----------------------- first loop over events ---------------------------- 
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if ( ientry < 0 ) break;   
    nb = fChain->GetEntry(jentry);   //nbytes += nb;
     
// --------------- selection of good track --------------------    

    if ( t_nVtx < NVERTEX_MIN  ||  t_nVtx > NVERTEX_MAX ) continue;
    if ( t_ieta < -maxIeta || t_ieta > maxIeta ) continue;
    
    if ( !goodTrack() ) continue;
    
    h_nvtx->Fill(t_nVtx,1.0);
          
    // ---- loop over active subdetectors in the event for total energy ---
    unsigned int nDets = (*t_DetIds).size();
    h_cluster->Fill(t_ieta, nDets);
      
    // check for possibility to correct for PU
    double correctionForPU(1.0);
    int abs_t_ieta = abs(t_ieta);
    
    if ( applyCorPU ) { 
      double de2p = (t_eHcal30 - t_eHcal10)/t_p;
      if ( de2p > DELTA_CUT )
	correctionForPU = corPU(abs_t_ieta,de2p);
    } 
   
    //--- check whether the calculated correction for PU makes sense
    if ( correctionForPU <= 0 || correctionForPU > 1 ) continue;
    nSelectedEvents++;

    std::map<unsigned int, bool> sameSubdet;
    sameSubdet.clear();

    for (unsigned int idet = 0; idet < nDets; idet++) { 
      unsigned int detId = ( (*t_DetIds)[idet] & MASK ) | MASK2 ;

      if (nPhiMergedInEvent.find(detId) != nPhiMergedInEvent.end()) 
	nPhiMergedInEvent[detId]++;
      else 
	nPhiMergedInEvent.insert(std::pair<unsigned int,int>(detId, 1));
		
      if (nTrks.find(detId) != nTrks.end()) {
	if ( sameSubdet.find(detId) == sameSubdet.end() ) {
	  nTrks[detId]++;
	  nSubdetInEvent[detId] += nDets;
	  sameSubdet.insert(std::pair<unsigned int,bool>(detId, true));
	}
      }
      else {
	nTrks.insert(std::pair<unsigned int,int>(detId, 1));
	nSubdetInEvent.insert(std::pair<unsigned int,int>(detId, nDets));
	sameSubdet.insert(std::pair<unsigned int,bool>(detId, true));
	subDetector_trk.insert(std::pair<unsigned int,
			       int>( detId,((*t_DetIds)[idet] &0xe000000) / 0x2000000 ));
      }
      
    }

    //--- Number of tracks per ieta ------------------------------
    int jj = HALF_NUM_ETA_BINS + int(t_ieta/N_ETA_RINGS_PER_BIN);
    ntrk_ieta[jj]++;

  } // ------------------- end of first loop over events -------------------------------------

// ------------------- Calculate max number of tracks per ieta ------------
  for ( int j = 0; j < NUM_ETA_BINS; j++ ) {
    if ( maxNumOfTracksForIeta < ntrk_ieta[j] ) maxNumOfTracksForIeta = ntrk_ieta[j];
  }
//---------------------- Book individual e2p[ieta] ----------------------
  int n_ieta_bins = 0;
  n_ieta_bins = 2*2.5*pow(maxNumOfTracksForIeta,1/3.0);
  for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
    sprintf(name,"e2p[%02d]", i);
    e2p[i] = new TH1F(name, "",
		      n_ieta_bins,
		      MIN_RESPONSE_HIST, MAX_RESPONSE_HIST);
    e2p[i]->Sumw2();
    e2p[i]->GetXaxis()->SetTitle(sxlabel);
    sprintf(name,"Events/%5.3f", e2p[i]->GetBinWidth(1));
    e2p[i]->GetYaxis()->SetTitle(name);
  }
  
// ----------------------- second loop over events ---------------------------- 
  for (Long64_t jentry=0; jentry<nentries; jentry++) {
    Long64_t ientry = LoadTree(jentry);
    if ( ientry < 0 ) break;   
    nb = fChain->GetEntry(jentry);   //nbytes += nb;
     
// --------------- selection of good track --------------------    

    if ( t_nVtx < NVERTEX_MIN  ||  t_nVtx > NVERTEX_MAX ) continue;
    if ( t_ieta < -maxIeta || t_ieta > maxIeta ) continue;
    
    if ( !goodTrack() ) continue;
    
// --- Correction for PU  --------
    double eTotalCor(t_eHcal);
    double correctionForPU(1.0);
    int abs_t_ieta = abs(t_ieta);
    
    if ( applyCorPU ) { 
      double de2p = (t_eHcal30 - t_eHcal10)/t_p;
      if ( de2p > DELTA_CUT )
	correctionForPU = corPU(abs_t_ieta,de2p);
    }    

    eTotalCor = t_eHcal*correctionForPU;
    double response = eTotalCor/(t_p - t_eMipDR);

// --- Fill merged histograms ---------------------------      
    e2pALL->Fill(response ,1.0);

    if ( abs_t_ieta < FIRST_IETA_TR )
      e2pHB->Fill(response ,1.0);
    else if ( abs_t_ieta < FIRST_IETA_HE )
      e2pTR->Fill(response ,1.0);
    else
      e2pHE->Fill(response ,1.0);

    if ( response > maxRespForGoodTrack  )
      maxRespForGoodTrack = response;
    if ( response < minRespForGoodTrack )
      minRespForGoodTrack = response;

// --- Fill ieta histograms ---------------------------      
    int jj = HALF_NUM_ETA_BINS + int(t_ieta/N_ETA_RINGS_PER_BIN);
    e2p[jj]->Fill(response,1.0);
    
  } // ------------------- end of second loop over events -------------------------------------

// ------------------- Fill graphs with number of tracks and subdetectors ------------
  double jeta[N_DEPTHS][MAXNUM_SUBDET];
  double nTrk[N_DEPTHS][MAXNUM_SUBDET];
  double nSub[N_DEPTHS][MAXNUM_SUBDET];
  double nPhi[N_DEPTHS][MAXNUM_SUBDET];
  unsigned int kdep[N_DEPTHS];
  for ( unsigned ik = 0; ik < N_DEPTHS; ik++ ) { kdep[ik] = 0; }

  // fill number of tracks
  std::map <unsigned int,int>::iterator nTrksItr = nTrks.begin();
  for (nTrksItr = nTrks.begin(); nTrksItr != nTrks.end(); nTrksItr++ ) {
    unsigned int detId = nTrksItr->first;
    int depth= ((detId>>DEPTH_OFFSET) & DEPTH_MASK) + int(MERGE_PHI_AND_DEPTHS);
    int zside= (detId&ZSIDE_MASK) ? 1 : -1;
    unsigned int kcur = kdep[depth-1];
    
    jeta[depth-1][kcur] = int((detId>>ETA_OFFSET) & ETA_MASK)*zside;
    nTrk[depth-1][kcur] = nTrksItr->second;
    nSub[depth-1][kcur] = double(nSubdetInEvent[detId])/double(nTrksItr->second);
    nPhi[depth-1][kcur] = double(nPhiMergedInEvent[detId])/double(nTrksItr->second);
    kdep[depth-1]++;
  }
  unsigned ngraphs = N_DEPTHS;
  if ( MERGE_PHI_AND_DEPTHS ) ngraphs = 1;
  
  for ( unsigned ik = 0; ik < ngraphs; ik++ ) {
    double x[MAXNUM_SUBDET];
    double ytrk[MAXNUM_SUBDET], ysub[MAXNUM_SUBDET], yphi[MAXNUM_SUBDET];
    for ( unsigned im = 0; im < MAXNUM_SUBDET; im++ ) {
      x[im] = jeta[ik][im];
      ytrk[im] = nTrk[ik][im];
      ysub[im] = nSub[ik][im];
      yphi[im] = nPhi[ik][im];
    }
    TGraph*  g_ntrk = new TGraph(kdep[ik], x, ytrk);
    sprintf(name, "Number of tracks for depth %1d", ik+1);
    g_ntrk->SetTitle(name);
    g_ntrk->Draw("a*");
    sprintf(name, "nTrk_depth%1d", ik+1);
    foutRootFile->WriteTObject(g_ntrk, name);

    TGraph*  g_nsub = new TGraph(kdep[ik], x, ysub);
    sprintf(name, "Mean number of active subdetectors, depth %1d", ik+1);
    g_nsub->SetTitle(name);
    g_nsub->Draw("a*");
    sprintf(name, "nSub_depth%1d", ik+1);
    foutRootFile->WriteTObject(g_nsub, name);

    TGraph*  g_nphi = new TGraph(kdep[ik], x, yphi);
    sprintf(name, "Mean number of phi-merged subdetectors, depth %1d", ik+1);
    g_nphi->SetTitle(name);
    g_nphi->Draw("a*");
    sprintf(name, "nPhi_depth%1d", ik+1);
    foutRootFile->WriteTObject(g_nphi, name);
  }

  //--- estimate ratio mean/MPV
  double xl = e2pALL->GetMean() - RMSRANGE4FIRSTFIT*e2pALL->GetRMS();
  double xr = e2pALL->GetMean() + RMSRANGE4FIRSTFIT*e2pALL->GetRMS();
  f1->SetRange(xl,xr);
  e2pALL->Fit("f1","QR");
  xl = f1->GetParameter(1) - fitRange*f1->GetParameter(2);
  xr = f1->GetParameter(1) + fitRange*f1->GetParameter(2);
  f1->SetRange(xl,xr);
  e2pALL->Fit("f1","QR");
    
  std::cout << "<Mean from sample>/<MPV from fit> = "
	    << e2pALL->GetMean()/f1->GetParameter(1)
	    << "  (chi2ndf = " << f1->GetChisquare()/f1->GetNDF() << ")"
	    << std::endl;

  //---- for HB
  xl = e2pHB->GetMean() - RMSRANGE4FIRSTFIT*e2pHB->GetRMS();
  xr = e2pHB->GetMean() + RMSRANGE4FIRSTFIT*e2pHB->GetRMS();
  f1->SetRange(xl,xr);
  e2pHB->Fit("f1","QR");
  xl = f1->GetParameter(1) - fitRange*f1->GetParameter(2);
  xr = f1->GetParameter(1) + fitRange*f1->GetParameter(2);
  f1->SetRange(xl,xr);
  e2pHB->Fit("f1","QR");

  std::cout << "In HB: <Mean from sample>/<MPV from fit> = "
	    << e2pHB->GetMean()/f1->GetParameter(1)
	    << "  (chi2ndf = " << f1->GetChisquare()/f1->GetNDF() << ")"
	    << std::endl;

  //---- for TR
  xl = e2pTR->GetMean() - RMSRANGE4FIRSTFIT*e2pTR->GetRMS();
  xr = e2pTR->GetMean() + RMSRANGE4FIRSTFIT*e2pTR->GetRMS();
  f1->SetRange(xl,xr);
  e2pTR->Fit("f1","QR");
  xl = f1->GetParameter(1) - fitRange*f1->GetParameter(2);
  xr = f1->GetParameter(1) + fitRange*f1->GetParameter(2);
  f1->SetRange(xl,xr);
  e2pTR->Fit("f1","QR");

  std::cout << "In TR: <Mean from sample>/<MPV from fit> = "
	    << e2pTR->GetMean()/f1->GetParameter(1)
	    << "  (chi2ndf = " << f1->GetChisquare()/f1->GetNDF() << ")"
	    << std::endl;

  //---- for HE
  xl = e2pHE->GetMean() - RMSRANGE4FIRSTFIT*e2pHE->GetRMS();
  xr = e2pHE->GetMean() + RMSRANGE4FIRSTFIT*e2pHE->GetRMS();
  f1->SetRange(xl,xr);
  e2pHE->Fit("f1","QR");
  xl = f1->GetParameter(1) - fitRange*f1->GetParameter(2);
  xr = f1->GetParameter(1) + fitRange*f1->GetParameter(2);
  f1->SetRange(xl,xr);
  e2pHE->Fit("f1","QR");

  std::cout << "In HE: <Mean from sample>/<MPV from fit> = "
	    << e2pHE->GetMean()/f1->GetParameter(1)
	    << "  (chi2ndf = " << f1->GetChisquare()/f1->GetNDF() << ")"
	    << std::endl;

  //----- print additional info
  std::cout << "Minimal response for good tracks = " 
	    << minRespForGoodTrack
	    << std::endl;
  std::cout << "Maximal response for good tracks = " 
	    << maxRespForGoodTrack
	    << std::endl;
  std::cout << "Maximum number of selected tracks per ieta bin = " 
	    << maxNumOfTracksForIeta
	    << std::endl;
  std::cout << "Number of selected tracks in HB = "
	    << e2pHB->GetEntries()
	    << std::endl;
  std::cout << "Number of selected tracks in TR = "
	    << e2pTR->GetEntries()
	    << std::endl;
  std::cout << "Number of selected tracks in HE = "
	    << e2pHE->GetEntries()
	    << std::endl;

  //--- Response versus ieta --------------------------
  TGraph *g_chi = new TGraph(NUM_ETA_BINS);  
  TGraphErrors* g_e2pFit = new TGraphErrors(NUM_ETA_BINS);
  TGraphErrors* g_e2pMean = new TGraphErrors(NUM_ETA_BINS);
  
  int ipointF(0);
  int ipointM(0);
  for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
    int ieta = (i - HALF_NUM_ETA_BINS)*N_ETA_RINGS_PER_BIN;
    if ( N_ETA_RINGS_PER_BIN > 1 ) {
      ieta = (i > HALF_NUM_ETA_BINS) ? ieta+1 : ieta-1;
    }
    if ( abs(ieta) > maxIeta ) continue;

    int nhistentries = e2p[i]->GetEntries();
    if ( nhistentries < 1 ) continue;
    else {
      g_e2pMean->SetPoint(ipointM, ieta, e2p[i]->GetMean());
      g_e2pMean->SetPointError(ipointM, 0, e2p[i]->GetMeanError());
      ipointM++;

      if ( nhistentries >= MIN_N_ENTRIES_FOR_FIT ) {
	
	double xl = e2p[i]->GetMean() - RMSRANGE4FIRSTFIT*e2p[i]->GetRMS();
	double xr = e2p[i]->GetMean() + RMSRANGE4FIRSTFIT*e2p[i]->GetRMS();
	f1->SetRange(xl,xr);
	e2p[i]->Fit("f1","QR");
	xl = f1->GetParameter(1) - fitRange*f1->GetParameter(2);
	xr = f1->GetParameter(1) + fitRange*f1->GetParameter(2);
	f1->SetRange(xl,xr);
	e2p[i]->Fit("f1","QR");
	g_e2pFit->SetPoint(ipointF, ieta, f1->GetParameter(1));
	g_e2pFit->SetPointError(ipointF, 0, f1->GetParError(1));
	g_chi->SetPoint(ipointF, ieta, f1->GetChisquare()/f1->GetNDF());
	ipointF++;
      }
    }
  }
  // remove empty points
  for ( int k = ipointF; k < NUM_ETA_BINS; k++ ) {
    g_e2pFit->RemovePoint(ipointF);
    g_chi->RemovePoint(ipointF);
  }
  for ( int k = ipointM; k < NUM_ETA_BINS; k++ ) {
    g_e2pMean->RemovePoint(ipointM);
    }
  
  sprintf(name, "Response from fit");
  g_e2pFit->SetTitle(name);
  g_e2pFit->GetXaxis()->SetTitle("i#eta");
  sprintf(name, "respFit");
  foutRootFile->WriteTObject(g_e2pFit, name);

  sprintf(name, "Mean response");
  g_e2pMean->SetTitle(name);
  g_e2pMean->GetXaxis()->SetTitle("i#eta");
  sprintf(name, "respMean");
  foutRootFile->WriteTObject(g_e2pMean, name);

  sprintf(name, "Chi2/NDF");
  g_chi->SetTitle(name);
  g_chi->GetXaxis()->SetTitle("i#eta");
  sprintf(name, "chi2ndf");
  foutRootFile->WriteTObject(g_chi, name); 

// ----- Clone ieta hists with non-zero entries for saving ------------
  for ( int i = 0; i < NUM_ETA_BINS; i++ ) {
    int inum = (i - HALF_NUM_ETA_BINS)*N_ETA_RINGS_PER_BIN;
    sprintf(name,"ieta_%d", inum);
    if ( e2p[i]->GetEntries() > 0 ) e2p[i]->Clone(name);
      delete e2p[i];
  }

  return nSelectedEvents;
}

//**********************************************************
// Isolated track selection
//**********************************************************
Bool_t CalibTree::goodTrack()
{

  bool ok = ( (t_selectTk)          // track quality
	      && (t_qltyMissFlag)
	      && (t_hmaxNearP < limCharIso)  // charge isolation
	      && (t_eMipDR < limMipEcal) 
	      && (t_p > minTrackMom) && (t_p < maxTrackMom)  // track momentum range
	      && (t_pt >= minTrackPt)               // constraint on track pt
	      && (t_eHcal >= minEnrHcal)            // constraint on Hcal energy
	      && (t_eHcal/t_p < MAX_RESPONSE_HIST)
		 // reject events with too big cluster energy
	     );
  return ok;
}
//********************************************************************************
