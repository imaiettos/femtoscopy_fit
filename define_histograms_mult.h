#include "call_libraries.h"  // call libraries from ROOT and C++
//#include "random_mixing_5.h" // random mixing
//#include "random_mixing_2.h"
#include "eta_mixing.h"

// define the bins
// qinv
const int nQBins = 1600;   // number of qinv bins
const double minQ = 0.0;  // minimum qinv
const double maxQ = 8.0;  // maximumm qinv

//const int nQBins = 600;   // number of qinv bins
//const double minQ = 0.0;  // minimum qinv 
//const double maxQ = 3.0;  // maximumm qinv

// q3D
const int nQBins3D = 100;   // number of q3D bins
const double minQ3D = 0.0;  // minimum q3D
const double maxQ3D = 2.0;  // maximumm q3D

//kT
//const int nKtBins = 11; // number of average transverse momentum bins
//double KtBins[nKtBins+1] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.5}; 

const int nKtBins = 5;
const double KtBins[nKtBins+1] = {0.1,0.3,0.5,0.7,0.9,1.5};
//const double KtBins[nKtBins+1] = {0.0,0.2,0.3,0.4,0.5,0.6,0.8,1.0};


// multiplicity
//const int nCentBins = 11; // number of multiplicity bins
//double CentBins[nCentBins+1] = {0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0, 70.0, 100.0, 150.0, 200.0, 250.0};  // multiplicity bins to match PbPb
const int nCentBins = 6; // number of multiplicity bins
double CentBins[nCentBins+1] = {0.0, 5.0, 10.0, 15.0,  25.0,  35.0, 100.0};  // multiplicity bins to match PbPb

// Event histograms
TH1I *Nevents = new TH1I("Nevents", "Nevents", 10, 0, 10);
TH1D *vzhist_beforecuts = new TH1D("vzhist_beforecuts", "vzhist_beforecuts", 80, -20., 20.);
TH1D *vzhist = new TH1D("vzhist", "vzhist", 80, -20., 20.);
TH1D *multiplicity_beforecuts = new TH1D("multiplicity_beforecuts", "multiplicity_beforecuts", 400, 0.0, 400.0);
TH1D *multiplicity = new TH1D("multiplicity", "multiplicity", 400, 0.0, 400.0);
TH1D *multposplus = new TH1D("multposplus", "multposplus", 400, 0.0, 400.0);
TH1D *multposminus = new TH1D("multposminus", "multposminus", 400, 0.0, 400.0);
TH1I *NeventsAss = new TH1I("NeventsAss", "NeventsAss", 20, 0, 20);
TH1D *hkT = new TH1D("hkT", "hkT", 100, 0, 1.5);
TH1D* hkT_pos = new TH1D("hkT_pos", "Relative p_{T} (#eta_{#gamma} > 0);p_{T} [GeV/c];Counts", 100, 0, 1.5);
TH1D* hkT_neg = new TH1D("hkT_neg", "Relative p_{T} (#eta_{#gamma} < 0);p_{T} [GeV/c];Counts", 100, 0, 1.5);

TH1D *hpT = new TH1D("hpT", "hpT", 100, 0, 1.5);
TH1D* hpT_pos = new TH1D("hpT_pos", "Relative p_{T} (#eta_{#gamma} > 0);p_{T} [GeV/c];Counts", 100, 0, 1.5);
TH1D* hpT_neg = new TH1D("hpT_neg", "Relative p_{T} (#eta_{#gamma} < 0);p_{T} [GeV/c];Counts", 100, 0, 1.5);

TH1D *heta = new TH1D("heta", "heta", 100, -2.4, 2.4);
TH1D* heta_pos = new TH1D("heta_pos", "Relative #eta (#eta_{#gamma} > 0); #eta;Counts", 100, -2.4, 2.4); 
TH1D* heta_neg = new TH1D("heta_neg", "Relative #eta (#eta_{#gamma} < 0); #eta;Counts", 100, -2.4, 2.4);

// ZDC
TH1D *hZDCPlusZeroBias_beforecuts  = new TH1D( "hZDCPlusZeroBias_beforecuts",  "", 150, 0, 15000 );
TH1D *hZDCMinusZeroBias_beforecuts = new TH1D( "hZDCMinusZeroBias_beforecuts", "", 150, 0, 15000 );
TH1D *hZDCPlusZDCOr_beforecuts  = new TH1D( "hZDCPlusZDCOr_beforecuts", "", 150, 0, 15000 );
TH1D *hZDCMinusZDCOr_beforecuts = new TH1D( "hZDCMinusZDCOr_beforecuts", "", 150, 0, 15000 );
TH1D *hZDCPlusZeroBias  = new TH1D( "hZDCPlusZeroBias",  "Zero Bias Trigger ZDC Plus Distribution;ZDC Offline Energy Sum (GeV);", 150, 0, 15000 );
TH1D *hZDCMinusZeroBias = new TH1D( "hZDCMinusZeroBias", "Zero Bias Trigger ZDC Minus Distribution;ZDC Offline Energy Sum (GeV);", 150, 0, 15000 );
TH1D *hZDCPlusZDCOr  = new TH1D( "hZDCPlusZDCOr", "ZDCOr Trigger ZDC Plus Distribution;ZDC Offline Energy Sum (GeV);", 150, 0, 15000 );
TH1D *hZDCMinusZDCOr = new TH1D( "hZDCMinusZDCOr", "ZDCOr Trigger ZDC Minus Distribution;ZDC Offline Energy Sum (GeV);", 150, 0, 15000 );

// PF Candidates
TH2D *h_h      = new TH2D( "h_h",     "Charged Hadrons; #eta; p_{T}", 1100, -5.5, 5.5, 500, 0, 50 );
TH2D *h_e      = new TH2D( "h_e",     "Electrons; #eta; p_{T}", 1100, -5.5, 5.5, 500, 0, 50 );
TH2D *h_mu     = new TH2D( "h_mu",    "Muons; #eta; p_{T}", 1100, -5.5, 5.5, 500, 0, 50 );
TH2D *h_gamma  = new TH2D( "h_gamma", "Photons; #eta; E_{T}", 1100, -5.5, 5.5, 500, 0, 50 );
TH2D *h_h0     = new TH2D( "h_h0",    "Neutral Hadrons; #eta; E_{T}", 1100, -5.5, 5.5, 500, 0, 50 );
TH2D *h_HFhad  = new TH2D( "h_HFhad", "HF Hadronic; #eta; E", 1100, -5.5, 5.5, 500, 0, 50 );
TH2D *h_HFem   = new TH2D( "h_HFem",  "HF Electromagnetic; #eta; E", 1100, -5.5, 5.5, 500, 0, 50 );

TH2D *trigGammaAMult   = new TH2D( "trigGammaAMult",  "trigGammaAMult", 400, 0, 400, 800, 0, 400 );
TH2D *trigGammaAMult_pos   = new TH2D( "trigGammaAMult_pos",  "trigGammaAMult_pos", 400, 0, 400, 800, 0, 400 );
TH2D *trigGammaAMult_neg   = new TH2D( "trigGammaAMult_neg",  "trigGammaAMult_neg", 400, 0, 400, 800, 0, 400 );

//Nch in function of trigger

// For gap info
// Axis : 0 -> etagapplus, 1 -> etagapminus, 2 -> multiplicity , 3 -> ETSum, 4 -> itrigger, 5 -> photon side
int    bins_gap[6]   =   { 60    ,  60   ,   250    ,  800 , 9    , 3 };
double xmin_gap[6]   =   { 0.0   ,  0.0  ,   0.0    ,  0.0 , -0.5 , -1.5};
double xmax_gap[6]   =   { 6.0  ,  6.0 ,   250.0  ,  400 , 8.5  , 1.5};
// --> Reco
THnSparseD *hist_gaps = new THnSparseD("hist_gaps", "hist_gaps", 6, bins_gap, xmin_gap, xmax_gap);

// Track/Particle histograms
// Axis : 0 -> track pT, 1 -> trk eta, 2 -> trk phi, 3 -> trk charge, 4 -> multiplicity bin, 5 -> photon side
int	bins_trk[6]      =   { 200   ,  24  ,   30				      , 3   , nCentBins  			, 3 };
double xmin_trk[6]   =   { 0.0   , -2.4 ,   -TMath::Pi()  		  , -1.5, CentBins[0]			, -1.5};
double xmax_trk[6]   =   { 50.0  ,  2.4 ,   TMath::Pi()  		  ,  1.5, CentBins[nCentBins]	, 1.5};
// --> Reco
THnSparseD *hist_reco_trk = new THnSparseD("hist_reco_trk", "hist_reco_trk", 6, bins_trk, xmin_trk, xmax_trk);
THnSparseD *hist_reco_trk_corr = new THnSparseD("hist_reco_trk_corr", "hist_reco_trk_corr", 6, bins_trk, xmin_trk, xmax_trk);

// HBT histograms for tests
// Axis : 0 -> qinv, 1 -> kT, 2 -> multiplicity bin, 3 -> photon side
int	bins_qinv[4]  =   { nQBins , nKtBins  		   ,   nCentBins, 3};
double xmin_qinv[4]   =   { minQ   , KtBins[0] 		   ,   CentBins[0], -1.5};
double xmax_qinv[4]   =   { maxQ   , KtBins[nKtBins] ,   CentBins[nCentBins], 1.5};
THnSparseD *hist_qinv_SS = new THnSparseD("hist_qinv_SS", "hist_qinv_SS", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_INV = new THnSparseD("hist_qinv_SS_INV", "hist_qinv_SS_INV", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_ROT = new THnSparseD("hist_qinv_SS_ROT", "hist_qinv_SS_ROT", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_MIX = new THnSparseD("hist_qinv_SS_MIX", "hist_qinv_SS_MIX", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS = new THnSparseD("hist_qinv_OS", "hist_qinv_OS", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_INV = new THnSparseD("hist_qinv_OS_INV", "hist_qinv_OS_INV", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_ROT = new THnSparseD("hist_qinv_OS_ROT", "hist_qinv_OS_ROT", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_MIX = new THnSparseD("hist_qinv_OS_MIX", "hist_qinv_OS_MIX", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_gen = new THnSparseD("hist_qinv_SS_gen", "hist_qinv_SS_gen", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_gen_INV = new THnSparseD("hist_qinv_SS_gen_INV", "hist_qinv_SS_gen_INV", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_gen_ROT = new THnSparseD("hist_qinv_SS_gen_ROT", "hist_qinv_SS_gen_ROT", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_SS_gen_MIX = new THnSparseD("hist_qinv_SS_gen_MIX", "hist_qinv_SS_gen_MIX", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_gen = new THnSparseD("hist_qinv_OS_gen", "hist_qinv_OS_gen", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_gen_INV = new THnSparseD("hist_qinv_OS_gen_INV", "hist_qinv_OS_gen_INV", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_gen_ROT = new THnSparseD("hist_qinv_OS_gen_ROT", "hist_qinv_OS_gen_ROT", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qinv_OS_gen_MIX = new THnSparseD("hist_qinv_OS_gen_MIX", "hist_qinv_OS_gen_MIX", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_SS = new THnSparseD("hist_qlcms_SS", "hist_qlcms_SS", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_SS_INV = new THnSparseD("hist_qlcms_SS_INV", "hist_qlcms_SS_INV", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_SS_ROT = new THnSparseD("hist_qlcms_SS_ROT", "hist_qlcms_SS_ROT", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_SS_MIX = new THnSparseD("hist_qlcms_SS_MIX", "hist_qlcms_SS_MIX", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_OS = new THnSparseD("hist_qlcms_OS", "hist_qlcms_OS", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_OS_INV = new THnSparseD("hist_qlcms_OS_INV", "hist_qlcms_OS_INV", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_OS_ROT = new THnSparseD("hist_qlcms_OS_ROT", "hist_qlcms_OS_ROT", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_OS_MIX = new THnSparseD("hist_qlcms_OS_MIX", "hist_qlcms_OS_MIX", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_SS_gen = new THnSparseD("hist_qlcms_SS_gen", "hist_qlcms_SS_gen", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_SS_gen_INV = new THnSparseD("hist_qlcms_SS_gen_INV", "hist_qlcms_SS_gen_INV", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_SS_gen_ROT = new THnSparseD("hist_qlcms_SS_gen_ROT", "hist_qlcms_SS_gen_ROT", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_SS_gen_MIX = new THnSparseD("hist_qlcms_SS_gen_MIX", "hist_qlcms_SS_gen_MIX", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_OS_gen = new THnSparseD("hist_qlcms_OS_gen", "hist_qlcms_OS_gen", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_OS_gen_INV = new THnSparseD("hist_qlcms_OS_gen_INV", "hist_qlcms_OS_gen_INV", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_OS_gen_ROT = new THnSparseD("hist_qlcms_OS_gen_ROT", "hist_qlcms_OS_gen_ROT", 4, bins_qinv, xmin_qinv, xmax_qinv);
THnSparseD *hist_qlcms_OS_gen_MIX = new THnSparseD("hist_qlcms_OS_gen_MIX", "hist_qlcms_OS_gen_MIX", 4, bins_qinv, xmin_qinv, xmax_qinv);

// HBT histograms for tests
// Axis : 0 -> qlong, 1 -> qout, 2 -> qside, 3 -> kT, 4 -> multiplicity bin, 5 -> photon side
int	bins_q3D[6]      =   { nQBins3D, nQBins3D, nQBins3D,  nKtBins 		      , nCentBins            , 3};
double xmin_q3D[6]   =   { minQ3D  , minQ3D  , minQ3D  ,  KtBins[0] 		  , CentBins[0]          , -1.5};
double xmax_q3D[6]   =   { maxQ3D  , maxQ3D  , maxQ3D  ,  KtBins[nKtBins]     , CentBins[nCentBins]  , 1.5};
THnSparseD *hist_q3D_SS = new THnSparseD("hist_q3D_SS", "hist_q3D_SS", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_SS_INV = new THnSparseD("hist_q3D_SS_INV", "hist_q3D_SS_INV", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_SS_ROT = new THnSparseD("hist_q3D_SS_ROT", "hist_q3D_SS_ROT", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_SS_MIX = new THnSparseD("hist_q3D_SS_MIX", "hist_q3D_SS_MIX", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_OS = new THnSparseD("hist_q3D_OS", "hist_q3D_OS", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_OS_INV = new THnSparseD("hist_q3D_OS_INV", "hist_q3D_OS_INV", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_OS_ROT = new THnSparseD("hist_q3D_OS_ROT", "hist_q3D_OS_ROT", 6, bins_q3D, xmin_q3D, xmax_q3D);
THnSparseD *hist_q3D_OS_MIX = new THnSparseD("hist_q3D_OS_MIX", "hist_q3D_OS_MIX", 6, bins_q3D, xmin_q3D, xmax_q3D);

TH2D *hist_dpt_cos_SS = new TH2D("hist_dpt_cos_SS", "hist_dpt_cos_SS",1000, 0.99910, 1.0001, 100, 0, 0.25);
TH2D *hist_dpt_cos_OS = new TH2D("hist_dpt_cos_OS", "hist_dpt_cos_OS",1000, 0.99910, 1.0001, 100, 0, 0.25);
TH1D *hist_pairSS_Mass = new TH1D("hist_pairSS_Mass", "Invariant mass same-sign tracks", 500, 0, 1.0);
TH1D *hist_pairOS_Mass = new TH1D("hist_pairOS_Mass", "Invariant mass opposite-sign tracks", 500, 0, 1.0);
TH2D *hist_detadphi_SS = new TH2D("hist_detadphi_SS", "hist_detadphi_SS",1000, 0.0, 0.1, 1000, 0, 0.1);
TH2D *hist_detadphi_OS = new TH2D("hist_detadphi_OS", "hist_detadphi_OS",1000, 0.0, 0.1, 1000, 0, 0.1);

void sw2(const double* KtBinning, const double* CentBinning){

hist_reco_trk->GetAxis(4)->Set(bins_trk[4],CentBinning);	
hist_reco_trk_corr->GetAxis(4)->Set(bins_trk[4],CentBinning);
hist_qinv_SS->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qinv_SS_INV->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qinv_SS_ROT->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qinv_SS_MIX->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qinv_OS->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qinv_OS_INV->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qinv_OS_ROT->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qinv_OS_MIX->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qinv_SS->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qinv_SS_INV->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qinv_SS_ROT->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qinv_SS_MIX->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qinv_OS->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qinv_OS_INV->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qinv_OS_ROT->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qinv_OS_MIX->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qlcms_SS->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qlcms_SS_INV->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qlcms_SS_ROT->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qlcms_SS_MIX->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qlcms_OS->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qlcms_OS_INV->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qlcms_OS_ROT->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qlcms_OS_MIX->GetAxis(1)->Set(bins_qinv[1],KtBinning);
hist_qlcms_SS->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qlcms_SS_INV->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qlcms_SS_ROT->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qlcms_SS_MIX->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qlcms_OS->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qlcms_OS_INV->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qlcms_OS_ROT->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_qlcms_OS_MIX->GetAxis(2)->Set(bins_qinv[2],CentBinning);
hist_q3D_SS->GetAxis(3)->Set(bins_q3D[3],KtBinning);
hist_q3D_SS_INV->GetAxis(3)->Set(bins_q3D[3],KtBinning);
hist_q3D_SS_ROT->GetAxis(3)->Set(bins_q3D[3],KtBinning);
hist_q3D_SS_MIX->GetAxis(3)->Set(bins_q3D[3],KtBinning);
hist_q3D_OS->GetAxis(3)->Set(bins_q3D[3],KtBinning);
hist_q3D_OS_INV->GetAxis(3)->Set(bins_q3D[3],KtBinning);
hist_q3D_OS_ROT->GetAxis(3)->Set(bins_q3D[3],KtBinning);
hist_q3D_OS_MIX->GetAxis(3)->Set(bins_q3D[3],KtBinning);
hist_q3D_SS->GetAxis(4)->Set(bins_q3D[4],CentBinning);
hist_q3D_SS_INV->GetAxis(4)->Set(bins_q3D[4],CentBinning);
hist_q3D_SS_ROT->GetAxis(4)->Set(bins_q3D[4],CentBinning);
hist_q3D_SS_MIX->GetAxis(4)->Set(bins_q3D[4],CentBinning);
hist_q3D_OS->GetAxis(4)->Set(bins_q3D[4],CentBinning);
hist_q3D_OS_INV->GetAxis(4)->Set(bins_q3D[4],CentBinning);
hist_q3D_OS_ROT->GetAxis(4)->Set(bins_q3D[4],CentBinning);
hist_q3D_OS_MIX->GetAxis(4)->Set(bins_q3D[4],CentBinning);

/*
hist_qinv_SS->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qinv_SS_INV->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qinv_SS_ROT->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qinv_SS_MIX->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qinv_OS->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qinv_OS_INV->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qinv_OS_ROT->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qinv_OS_MIX->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qinv_SS->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qinv_SS_INV->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qinv_SS_ROT->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qinv_SS_MIX->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qinv_OS->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qinv_OS_INV->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qinv_OS_ROT->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qinv_OS_MIX->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qlcms_SS->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qlcms_SS_INV->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qlcms_SS_ROT->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qlcms_SS_MIX->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qlcms_OS->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qlcms_OS_INV->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qlcms_OS_ROT->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qlcms_OS_MIX->GetAxis(1)->Set(bins_qinv[1],KtBins);
hist_qlcms_SS->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qlcms_SS_INV->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qlcms_SS_ROT->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qlcms_SS_MIX->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qlcms_OS->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qlcms_OS_INV->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qlcms_OS_ROT->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_qlcms_OS_MIX->GetAxis(2)->Set(bins_qinv[2],CentBins);
hist_q3D_SS->GetAxis(3)->Set(bins_q3D[3],KtBins);
hist_q3D_SS_INV->GetAxis(3)->Set(bins_q3D[3],KtBins);
hist_q3D_SS_ROT->GetAxis(3)->Set(bins_q3D[3],KtBins);
hist_q3D_SS_MIX->GetAxis(3)->Set(bins_q3D[3],KtBins);
hist_q3D_OS->GetAxis(3)->Set(bins_q3D[3],KtBins);
hist_q3D_OS_INV->GetAxis(3)->Set(bins_q3D[3],KtBins);
hist_q3D_OS_ROT->GetAxis(3)->Set(bins_q3D[3],KtBins);
hist_q3D_OS_MIX->GetAxis(3)->Set(bins_q3D[3],KtBins);
hist_q3D_SS->GetAxis(4)->Set(bins_q3D[4],CentBins);
hist_q3D_SS_INV->GetAxis(4)->Set(bins_q3D[4],CentBins);
hist_q3D_SS_ROT->GetAxis(4)->Set(bins_q3D[4],CentBins);
hist_q3D_SS_MIX->GetAxis(4)->Set(bins_q3D[4],CentBins);
hist_q3D_OS->GetAxis(4)->Set(bins_q3D[4],CentBins);
hist_q3D_OS_INV->GetAxis(4)->Set(bins_q3D[4],CentBins);
hist_q3D_OS_ROT->GetAxis(4)->Set(bins_q3D[4],CentBins);
hist_q3D_OS_MIX->GetAxis(4)->Set(bins_q3D[4],CentBins);
*/
//void sw2()
	// Event histograms
Nevents->Sumw2();
vzhist_beforecuts->Sumw2();
vzhist->Sumw2();
multiplicity_beforecuts->Sumw2();
multiplicity->Sumw2();
multposplus->Sumw2();
multposminus->Sumw2();
NeventsAss->Sumw2();
hkT->Sumw2();
hkT_pos->Sumw2();
hkT_neg->Sumw2();
hpT->Sumw2();
hpT_pos->Sumw2();
hpT_neg->Sumw2();
heta->Sumw2();
heta_pos->Sumw2();
heta_neg->Sumw2();
// ZDC
hZDCPlusZeroBias_beforecuts->Sumw2();
hZDCMinusZeroBias_beforecuts->Sumw2();
hZDCPlusZDCOr_beforecuts->Sumw2();
hZDCMinusZDCOr_beforecuts->Sumw2();
hZDCPlusZeroBias->Sumw2();
hZDCMinusZeroBias->Sumw2();
hZDCPlusZDCOr->Sumw2();
hZDCMinusZDCOr->Sumw2();
// PF Candidates
h_h->Sumw2();
h_e->Sumw2();
h_mu->Sumw2();
h_gamma->Sumw2();
h_h0->Sumw2();
h_HFhad->Sumw2();
h_HFem->Sumw2();
hist_gaps->Sumw2();
trigGammaAMult->Sumw2();
trigGammaAMult_pos->Sumw2();
trigGammaAMult_neg->Sumw2();

/*hist_reco_trk->Sumw2();
hist_reco_trk_corr->Sumw2();
hist_qinv_SS->Sumw2();
hist_qinv_SS_INV->Sumw2();
hist_qinv_SS_ROT->Sumw2();
hist_qinv_SS_MIX->Sumw2();
hist_qinv_OS->Sumw2();
hist_qinv_OS_INV->Sumw2();
hist_qinv_OS_ROT->Sumw2();
hist_qinv_OS_MIX->Sumw2();
hist_qlcms_SS->Sumw2();
hist_qlcms_SS_INV->Sumw2();
hist_qlcms_SS_ROT->Sumw2();
hist_qlcms_SS_MIX->Sumw2();
hist_qlcms_OS->Sumw2();
hist_qlcms_OS_INV->Sumw2();
hist_qlcms_OS_ROT->Sumw2();
hist_qlcms_OS_MIX->Sumw2();
hist_q3D_SS->Sumw2();
hist_q3D_SS_INV->Sumw2();
hist_q3D_SS_ROT->Sumw2();
hist_q3D_SS_MIX->Sumw2();
hist_q3D_OS->Sumw2();
hist_q3D_OS_INV->Sumw2();
hist_q3D_OS_ROT->Sumw2();
hist_q3D_OS_MIX->Sumw2();
hist_pairSS_Mass->Sumw2();
hist_pairOS_Mass->Sumw2();
hist_detadphi_SS->Sumw2();
hist_detadphi_OS->Sumw2();
*/

hist_reco_trk->Sumw2();
hist_qinv_SS->Sumw2();
hist_qinv_SS_INV->Sumw2();
hist_qinv_SS_ROT->Sumw2();
hist_qinv_SS_MIX->Sumw2();
hist_qinv_OS->Sumw2();
hist_qinv_OS_INV->Sumw2();
hist_qinv_OS_ROT->Sumw2();
hist_qinv_OS_MIX->Sumw2();
hist_qlcms_SS->Sumw2();
hist_qlcms_SS_INV->Sumw2();
hist_qlcms_SS_ROT->Sumw2();
hist_qlcms_SS_MIX->Sumw2();
hist_qlcms_OS->Sumw2();
hist_qlcms_OS_INV->Sumw2();
hist_qlcms_OS_ROT->Sumw2();
hist_qlcms_OS_MIX->Sumw2();
hist_q3D_SS->Sumw2();
hist_q3D_SS_INV->Sumw2();
hist_q3D_SS_ROT->Sumw2();
hist_q3D_SS_MIX->Sumw2();
hist_q3D_OS->Sumw2();
hist_q3D_OS_INV->Sumw2();
hist_q3D_OS_ROT->Sumw2();
hist_q3D_OS_MIX->Sumw2();
hist_pairSS_Mass->Sumw2();
hist_pairOS_Mass->Sumw2();
hist_detadphi_SS->Sumw2();
hist_detadphi_OS->Sumw2();

}

// histogram writting

void write_eventQA(){
	// Event histograms
	Nevents->Write();
	vzhist_beforecuts->Write();
	vzhist->Write();
	multiplicity_beforecuts->Write();
	multiplicity->Write();
	multposplus->Write();
	multposminus->Write();
	NeventsAss->Write();
	hkT->Write();
	hkT_pos->Write();
	hkT_neg->Write();
	
        hpT->Write();
        hpT_pos->Write();
        hpT_neg->Write();

        heta->Write();
        heta_pos->Write();
        heta_neg->Write();

	// ZDC
	hZDCPlusZeroBias_beforecuts->Write();
	hZDCMinusZeroBias_beforecuts->Write();
	hZDCPlusZDCOr_beforecuts->Write();
	hZDCMinusZDCOr_beforecuts->Write();
	hZDCPlusZeroBias->Write();
	hZDCMinusZeroBias->Write();
	hZDCPlusZDCOr->Write();
	hZDCMinusZDCOr->Write();
	// PF Candidates
	h_h->Write();
	h_e->Write();
	h_mu->Write();
	h_gamma->Write();
	h_h0 ->Write();
	h_HFhad->Write();
	h_HFem ->Write();
	trigGammaAMult->Write();
	trigGammaAMult_pos->Write();
	trigGammaAMult_neg->Write();
	// Gaps
	hist_gaps->Write();
	// Tracks
	hist_reco_trk->Write();
        hist_reco_trk_corr->Write();

}

void write_HBT1D(){
	hist_qinv_SS->Write();
	hist_qinv_SS_INV->Write();
	hist_qinv_SS_ROT->Write();
	hist_qinv_SS_MIX->Write();
	hist_qinv_OS->Write();
	hist_qinv_OS_INV->Write();
	hist_qinv_OS_ROT->Write();
	hist_qinv_OS_MIX->Write();
	hist_qlcms_SS->Write();
	hist_qlcms_SS_INV->Write();
	hist_qlcms_SS_ROT->Write();
	hist_qlcms_SS_MIX->Write();
	hist_qlcms_OS->Write();
	hist_qlcms_OS_INV->Write();
	hist_qlcms_OS_ROT->Write();
	hist_qlcms_OS_MIX->Write();
	hist_dpt_cos_SS->Write();
	hist_dpt_cos_OS->Write();
	hist_pairSS_Mass->Write();
	hist_pairOS_Mass->Write();
	hist_detadphi_SS->Write();
	hist_detadphi_OS->Write();
}

void write_HBT3D(){
	hist_q3D_SS->Write();
	hist_q3D_SS_INV->Write();
	hist_q3D_SS_ROT->Write();
	hist_q3D_SS_MIX->Write();
	hist_q3D_OS->Write();
	hist_q3D_OS_INV->Write();
	hist_q3D_OS_ROT->Write();
	hist_q3D_OS_MIX->Write();
}

