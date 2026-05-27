#include "read_tree.h" // read the TChains
#include "tracking_correction.h" // tracking correction
#include "define_histograms_mult.h" // histogram definition
#include "Constants.h"
#include "findBin.h"
#define pi_mass 0.1396

int Encode(int w){
    int Range[] = {3,5,8,10,13,16,20,25,30,200,10000};
    int i(0), j(0);
    while(w >= Range[i++]) j++;
    return j;
}

Double_t ComputeEventWeight(std::vector<ROOT::Math::PtEtaPhiMVector>& tracks){
    int nFwd(0), nBkw(0), nCtr(0);
    for(auto& trk : tracks){
        double eta = trk.Eta();
        if(eta < -0.8)     nBkw++;
        else if(eta > 0.8) nFwd++;
        else               nCtr++;
    }
    return (100*Encode(nFwd) + 10*Encode(nCtr) + Encode(nBkw));
}

void correlation_UPC_new(TString input_file, TString ouputfile, int doquicktest, int domixing, int Nmixevents, int mincentormult, float minvz, int hbt3d, int gamov, int syst){

    clock_t sec_start, sec_end;
    sec_start = clock(); // start timing measurement
    TDatime* date = new TDatime(); // to add date in the output file

    // ------------------------
    // Check input file
    // ------------------------
    if (gSystem->AccessPathName(input_file)) {
        std::cerr << "ERROR: Input file does not exist: " << input_file << std::endl;
        return;
    }

    std::ifstream fin(input_file.Data());
    if (!fin.is_open()) {
        std::cerr << "ERROR: Could not open input file: " << input_file << std::endl;
        return;
    }

    std::vector<std::string> lines;
    std::string line;
    while (std::getline(fin, line)) {
        if (!line.empty()) lines.push_back(line);
    }
    fin.close();

    if (lines.empty()) {
        std::cerr << "WARNING: Input file is empty: " << input_file << std::endl;
        return;
    }

    std::cout << "Successfully read " << lines.size() << " lines from " << input_file << std::endl;

    bool do_quicktest = (doquicktest != 0);
    bool do_mixing    = (domixing == 0);
    bool do_hbt3d     = (hbt3d == 0);
    //bool do_gamov     = (gamov != 0);

    //if(syst == 4 || syst == 5) do_gamov = true;


    bool do_gamov = true;
    
    bool dosplit = true;

    TString systematics = "nonapplied_nominal";
    if(syst == 0) systematics =  "nominal";
    if(syst == 1) systematics =  "vznarrow";
    if(syst == 2) systematics =  "vzwide";
    if(syst == 3) systematics =  "removeduplicatedcut";
    if(syst == 4) systematics =  "gamovplus15";
    if(syst == 5) systematics =  "gamovminus15";
    if(syst == 6){ systematics =  "Nmix20"; Nmixevents += 10; }
    if(syst == 7){ systematics =  "Nmix05"; Nmixevents -= 5; }
    if(syst == 8){ systematics =  "minvz3"; minvz += 1.0; }
    if(syst == 9){ systematics =  "minvz1"; minvz -= 1.0; }
    if(syst == 10){ systematics =  "mincentormult10"; mincentormult += 5; }
    if(syst == 11){ systematics =  "mincentormult3"; mincentormult -= 2; }
    if(syst == 12){ systematics =  "hfthres4"; }
    if(syst == 13){ systematics =  "hfthres6"; }
    if(syst == 14){ systematics =  "hfthres8"; }
    if(syst == 15){ systematics =  "lowpfthres"; }
    if(syst == 16){ systematics =  "highpfthres"; }

    // Read the list of input file(s)
    std::vector<TString> file_name_vector;
    for(const auto& entry : lines) {
        file_name_vector.push_back(entry.c_str());
    }

    if(file_name_vector.empty()) {
        std::cerr << "ERROR: No valid input files found in " << input_file << std::endl;
        return;
    }

    // Create chain and add files
    TChain *hea_tree = new TChain("analyzer/eventTree");
    for (const auto& fileName : file_name_vector) {
        TFile *testfile = TFile::Open(fileName,"READ");
        if(testfile && !testfile->IsZombie() && !testfile->TestBit(TFile::kRecovered)){
            std::cout << "Adding file " << fileName << " to the chain" << std::endl;
            hea_tree->Add(fileName);
        } else {
            std::cout << "File: " << fileName << " failed!" << std::endl;
        }
    }
    file_name_vector.clear();

    // Read branches
    read_tree(hea_tree);

    // Histogram setup
//    const int nKtBinsF = 11;
//    const double KtBinsF[nKtBinsF+1] = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.5};
    const int nKtBinsF = 5;
    const double KtBinsF[nKtBinsF+1] = {0.1,0.3,0.5,0.7,0.9,1.5};
    //const double KtBinsF[nKtBinsF+1] = {0.0,0.2,0.3,0.4,0.5,0.6,0.8,1.0};

    const int nCentBinsF = 6;
    const double CentBinsF[nCentBinsF+1] = {0.0, 5.0, 10.0, 15.0, 25.0, 35.0, 100.0};
    sw2(KtBinsF, CentBinsF);

    int nevents = hea_tree->GetEntries();
    std::cout << "Total number of events: " << nevents << std::endl;

std::vector<std::pair<Double_t, std::vector<ROOT::Math::PtEtaPhiMVector>>> track_4vector;
std::vector<std::pair<Double_t, std::vector<double>>>                       track_weights_vector;
std::vector<std::pair<Double_t, std::vector<int>>>                          track_charge_vector;
std::vector<std::pair<Double_t, int>>                                       multiplicity_vector;
std::vector<double> vz_vector;
std::vector<int>    phpos_vector;
std::vector<int>    trigger_vector;

    // Open efficiency file once
    TFile *effFile = TFile::Open("GeneralTracks_2DEFFTable_2023pbpb_upc.root");
    if (!effFile || effFile->IsZombie()) {
    	std::cerr << "ERROR: Could not open efficiency file!" << std::endl;
    	return;
	}
    TH2D *hTrkEff = (TH2D*)effFile->Get("hEff_2D");
    TH2D *hTrkFak = (TH2D*)effFile->Get("hFak_2D");

    double nev = (double)nevents;

    // Event loop
    for (int i = 0; i < nevents; i++){
        hea_tree->GetEntry(i);

        if(i != 0 && i % 100000 == 0) {
            std::cout << "Running -> percentage: " << std::setprecision(3) << ((double)i/nev*100) << "%" << std::endl;
        }
       	if(i != 0 && (i % 100000) == 0){double alpha = (double)i; cout << " Running -> percentage: " << std::setprecision(3) << ((alpha / nev) * 100) << "%" << endl;} // % processed
 
        if(do_quicktest && i != 0 && i % 1000 == 0) break;

        int Ntroff = mult;
        double vertexz = (double)vtx[2];


        // fill quantities before any cut
        vzhist_beforecuts->Fill(vertexz);
        multiplicity_beforecuts->Fill(Ntroff);
				
	Nevents->Fill(0); // Filled after each event cut -> 0 means no cuts	
	// vertex selection 
        // Vertex and ZDC cuts
        if(syst == 1 && fabs(vertexz) > 3.0) continue;
        else if(syst == 2 && (fabs(vertexz) < 3.0 || fabs(vertexz) > 15.0)) continue;
        else if(syst != 1 && syst != 2 && fabs(vertexz) > 15.0) continue;

        Nevents->Fill(1); // Vertex Z cut	
   
	if(std::isnan(zdcSumPlus) || std::isnan(zdcSumMinus) || Ntroff < 0) continue;

	if( zdcSumPlus <= -50000 || zdcSumPlus > 300000 || zdcSumMinus <= -50000 || zdcSumMinus > 300000 ) continue; // clean up XZDC
	Nevents->Fill(2); // ZDC cut


	if( triggers[4] || triggers[5] ) { hZDCPlusZDCOr_beforecuts->Fill( zdcSumPlus ); hZDCMinusZDCOr_beforecuts->Fill( zdcSumMinus );}
	if( triggers[6] || triggers[7] ) { hZDCPlusZeroBias_beforecuts->Fill( zdcSumPlus ); hZDCMinusZeroBias_beforecuts->Fill( zdcSumMinus ); }

	float ZDC_MINUS_CUT = 800;
	float ZDC_PLUS_CUT = 1000;
	//If needed for systematics we can add it here using the systematic flag , e.g.
	//if( syst == XX ) {ZDC_MINUS_CUT = XXXX; ZDC_PLUS_CUT = XXXX;} 

	bool posPhoton;
	if( zdcSumMinus > ZDC_MINUS_CUT && zdcSumPlus < ZDC_PLUS_CUT ){ //photon traveling in + direction
		posPhoton = true;
		} else if( zdcSumMinus < ZDC_MINUS_CUT && zdcSumPlus > ZDC_PLUS_CUT ){ //photon traveling in - direction
			posPhoton = false;
		} else { continue; }
		Nevents->Fill(3); // Photon travel side cuts

	// HF threshold
	double HFThreshold = 0.0; // have to decide the nominal
	if(syst == 12 ) HFThreshold = 2.0;
	if(syst == 13 ) HFThreshold = 4.0;
	if(syst == 14 ) HFThreshold = 6.0;

	//Find Eta Gaps
	double etaGapPos = 0;
	double etaGapNeg = 0;
	double ETSum = 0;
	int numNegHFClusters=0, numPosHFClusters=0;
	etagaps(syst, HFThreshold, ETSum, numNegHFClusters, numPosHFClusters, pfPt, pfEta, pfE, pfID, isPrimary, etaGapPos, etaGapNeg, h_h, h_e, h_mu, h_gamma, h_h0, h_HFhad, h_HFem, posPhoton);

	if( HFThreshold > 0.0 && (numNegHFClusters < 1 || numPosHFClusters < 1) ) continue; // HF Coincidence
	Nevents->Fill(4);
	double sumGapCut = 3.0;
	//if( Ntroff > 35 ) sumGapCut = 2.8;
	//if( Ntroff < 0 ) continue; // remove events with multiplicity < 0
	//if( Ntroff > 100 ) continue; // remove events with multiplicity > 250
        //if( posPhoton && etaGapPos <= sumGapCut ) continue;
        //else if( !posPhoton && etaGapNeg <= sumGapCut ) continue;

	Nevents->Fill(5); // Multiplicity cut

	const int NUM_TRIGS = 8;
	string TRIG_STRINGS[ NUM_TRIGS ] =
  		        { "HLT_HIUPC_ZDC1nXOR_MBHF1AND_PixelTrackMultiplicity20_v",        //0
    			  "HLT_HIUPC_ZDC1nXOR_MBHF2AND_PixelTrackMultiplicity20_v",        //1
    			  "HLT_HIUPC_ZDC1nAsymXOR_MBHF1AND_PixelTrackMultiplicity20_v",    //2
    			  "HLT_HIUPC_ZDC1nAsymXOR_MBHF2AND_PixelTrackMultiplicity20_v",    //3
    			  "HLT_HIUPC_ZDC1nOR_SinglePixelTrackLowPt_MaxPixelCluster400_v",  //4
    			  "HLT_HIUPC_ZDC1nOR_MinPixelCluster400_MaxPixelCluster10000_v",   //5
    			  "HLT_HIUPC_ZeroBias_SinglePixelTrackLowPt_MaxPixelCluster400_v", //6
    			  "HLT_HIUPC_ZeroBias_MinPixelCluster400_MaxPixelCluster10000_v" };//7


	for( int iTrig = 0; iTrig < NUM_TRIGS - 2; iTrig++ ){
			bool firedTrigger = false;
			if( iTrig < 4 && triggers[ iTrig ] ) { firedTrigger = true; } else if( iTrig == 4 && (triggers[4] || triggers[5]) ){ firedTrigger = true; } else if( iTrig == 5 && (triggers[6] || triggers[7]) ) { firedTrigger = true; }
			if( !firedTrigger ) continue;
			double x_gap[6]={etaGapPos,etaGapNeg,(double) Ntroff, ETSum,(double) iTrig, (double) posPhoton};
			hist_gaps->Fill(x_gap);
			if( posPhoton && etaGapPos > sumGapCut )//SUM_GAP_CUT )
		            trigGammaAMult_pos->Fill( mult, ETSum );
          		else if( !posPhoton && etaGapNeg > sumGapCut )//SUM_GAP_CUT )
            		    trigGammaAMult_neg->Fill( mult, ETSum );
		}

	// can be also modified as systematics (if needed)
	//cout << "posPhoton: "<< (int) posPhoton << endl;
	if( posPhoton && etaGapPos <= sumGapCut ) continue;
	//if( !posPhoton && etaGapNeg <= sumGapCut ) continue;
        else if( !posPhoton && etaGapNeg <= sumGapCut ) continue;
	Nevents->Fill(6); // Multiplicity cut
      	
	if( Ntroff <= 15 && !triggers[4] && !triggers[5] ) continue;
	Nevents->Fill(7); // Multiplicity cut
	if( Ntroff > 15 && !triggers[0] && !triggers[1] && !triggers[2] && !triggers[3] ) continue;
	Nevents->Fill(8); // Multiplicity cut

	if ( triggers[7] ) { trigger_vector.push_back(7); }
	else if ( triggers[6] )  { trigger_vector.push_back(6); }
	else if ( Ntroff <= 15 && triggers[0] ) { trigger_vector.push_back(0); }
	else if ( Ntroff <= 15 && triggers[1] ) {trigger_vector.push_back(1); }
	else if ( Ntroff <= 15 && triggers[2] ) {trigger_vector.push_back(2); }
	else if ( Ntroff <= 15 && triggers[3] ) {trigger_vector.push_back(3); }
	else if ( Ntroff > 15 && triggers[4] ) {trigger_vector.push_back(4); }
	else if ( Ntroff > 15 && triggers[5] ) {trigger_vector.push_back(5); }

	// Fill event histograms after all cuts
	vzhist->Fill(vertexz);
	multiplicity->Fill( Ntroff );
	if(posPhoton) multposplus->Fill( Ntroff );
	if(!posPhoton) multposminus->Fill( Ntroff );
	if( triggers[4] || triggers[5] ) { hZDCPlusZDCOr->Fill( zdcSumPlus ); hZDCMinusZDCOr->Fill( zdcSumMinus );}
	if( triggers[6] || triggers[7] ) { hZDCPlusZeroBias->Fill( zdcSumPlus ); hZDCMinusZeroBias->Fill( zdcSumMinus ); }

        // ... [event cuts remain unchanged]

        std::vector<ROOT::Math::PtEtaPhiMVector> tracks_reco;
        std::vector<double> track_weight_reco;
        std::vector<int> track_charge_reco;

        int numPF = pfID->size();
	if(numPF == 0 || pfPt->size() != numPF || pfEta->size() != numPF || pfPhi->size() != numPF || pfCharge->size() != numPF || isPrimary->size() != numPF) {
    	continue; // skip this event if vectors are inconsistent
	}

        if(pfID->size() > 0){
            for (size_t j=0; j<pfID->size(); j++){
                double trkEta = pfEta->at(j);
                double trkPt  = pfPt->at(j);
                int trkCharge = pfCharge->at(j);
                double trkPhi  = pfPhi->at(j);

                if(trkPt <= 0.3) continue;
                if(pfID->at(j) != 1) continue;
                if(!isPrimary->at(j)) continue;
                if(fabs(trkEta) > 2.4) continue;
                if(trkNPixelHits->at(j) < 1) continue;

		//std::cout << "trkNPixelHits" << trkNPixelHits->at(j) << std::endl;


                int effEtaBin = hTrkEff->GetXaxis()->FindBin(trkEta);
                int effPtBin  = hTrkEff->GetYaxis()->FindBin(trkPt);
                double trgEff = hTrkEff->GetBinContent(effEtaBin, effPtBin);
                double trgFak = hTrkFak->GetBinContent(effEtaBin, effPtBin);
                double trk_weight = (trgEff>0) ? (1-trgFak)/(trgEff) : 1.0;
                double trkPt_corr = trk_weight*trkPt;

                double x_reco_trk[6]={trkPt,trkEta,pfPhi->at(j),(double)trkCharge,(double)Ntroff,(double)posPhoton};
                hist_reco_trk->Fill(x_reco_trk);
                hpT->Fill(trkPt,trk_weight);
                        if (posPhoton == 1) hpT_pos->Fill(trkPt, trk_weight);
                        else                hpT_neg->Fill(trkPt, trk_weight);

                heta->Fill(trkEta,trk_weight);
                        if (posPhoton == 1) heta_pos->Fill(trkEta, trk_weight);
                        else                heta_neg->Fill(trkEta, trk_weight);

/*                hpT->Fill(trkPt);
                        if (posPhoton == 1) hpT_pos->Fill(trkPt);
                        else                hpT_neg->Fill(trkPt);
*/
		for (size_t i = j + 1; i < pfID->size(); i++) {

			double trkEta_i  = pfEta->at(i);
        		double trkPt_i   = pfPt->at(i);
        		double trkPhi_i  = pfPhi->at(i);
        		int    trkCharge_i = pfCharge->at(i);
                        
			if(trkPt_i <= 0.3) continue;
                        if(pfID->at(i) != 1) continue;
                        if(!isPrimary->at(i)) continue;
                        if(fabs(trkEta_i) > 2.4) continue;
                        if(trkNPixelHits->at(i) < 1) continue;
        		if (trkCharge_i == 0) continue;


                int effEtaBin_i = hTrkEff->GetXaxis()->FindBin(trkEta_i);
                int effPtBin_i  = hTrkEff->GetYaxis()->FindBin(trkPt_i);
                double trgEff_i = hTrkEff->GetBinContent(effEtaBin_i, effPtBin_i);
                double trgFak_i = hTrkFak->GetBinContent(effEtaBin_i, effPtBin_i);
                double trk_weight_i = (trgEff_i>0) ? (1-trgFak_i)/(trgEff_i) : 1.0;

                //double trkPt_corr_i = trk_weight_i*trkPt_i;


        		// --- Compute ΔR ---
        		double dEta = trkEta - trkEta_i;
        		double dPhi = TVector2::Phi_mpi_pi(trkPhi - trkPhi_i); // handles wrap-around
        		double dR   = TMath::Sqrt(dEta*dEta + dPhi*dPhi);

        		// --- Compute relative kT ---
        		// Use the softer track's pT (common convention)
        		double minPt = TMath::Min(trkPt, trkPt_i);
        		double kT    = minPt * TMath::Sin(dR);
          		hkT->Fill(kT,trk_weight);

			// Fill the correct histogram based on photon direction
        		if (posPhoton == 1) hkT_pos->Fill(kT,trk_weight);
        		else                hkT_neg->Fill(kT,trk_weight);
                        
			/*hkT->Fill(kT);

                        // Fill the correct histogram based on photon direction
                        if (posPhoton == 1) hkT_pos->Fill(kT);
                        else                hkT_neg->Fill(kT);
*/		}	

		/*int effEtaBin = hTrkEff->GetXaxis()->FindBin(trkEta);
                int effPtBin  = hTrkEff->GetYaxis()->FindBin(trkPt);
                double trgEff = hTrkEff->GetBinContent(effEtaBin, effPtBin);
                double trgFak = hTrkFak->GetBinContent(effEtaBin, effPtBin);
                double trk_weight = (trgEff>0) ? (1-trgFak)/(4.*trgEff) : 1.0;
*/
                hist_reco_trk_corr->Fill(x_reco_trk, trk_weight);

                ROOT::Math::PtEtaPhiMVector TrackFourVector;
                TrackFourVector.SetM(pi_mass);
                TrackFourVector.SetPt(trkPt);
                TrackFourVector.SetEta(trkEta);
                TrackFourVector.SetPhi(pfPhi->at(j));

                tracks_reco.push_back(TrackFourVector);
                track_charge_reco.push_back(trkCharge);
                track_weight_reco.push_back(trk_weight);
            }
        }

        if(tracks_reco.size() > 1){
            twoparticlecorrelation(tracks_reco, track_charge_reco, track_weight_reco,
                                   hist_pairSS_Mass, hist_dpt_cos_SS, hist_detadphi_SS,
                                   hist_qinv_SS, hist_qinv_SS_INV, hist_qinv_SS_ROT,
                                   hist_qlcms_SS, hist_qlcms_SS_INV, hist_qlcms_SS_ROT,
                                   hist_q3D_SS, hist_q3D_SS_INV, hist_q3D_SS_ROT,
                                   hist_pairOS_Mass, hist_dpt_cos_OS, hist_detadphi_OS,
                                   hist_qinv_OS, hist_qinv_OS_INV, hist_qinv_OS_ROT,
                                   hist_qlcms_OS, hist_qlcms_OS_INV, hist_qlcms_OS_ROT,
                                   hist_q3D_OS, hist_q3D_OS_INV, hist_q3D_OS_ROT,
                                   Ntroff, (int)posPhoton, dosplit, do_hbt3d, do_gamov, syst);

            /*track_4vector.push_back(tracks_reco);
            track_charge_vector.push_back(track_charge_reco);
            track_weights_vector.push_back(track_weight_reco);
            multiplicity_vector.push_back(Ntroff);
            vz_vector.push_back(vertexz);
            phpos_vector.push_back((int)posPhoton);
        */
    Double_t etaWeight = ComputeEventWeight(tracks_reco);
    track_4vector.push_back({etaWeight, tracks_reco});
    track_charge_vector.push_back({etaWeight, track_charge_reco});
    track_weights_vector.push_back({etaWeight, track_weight_reco});
    multiplicity_vector.push_back({etaWeight, Ntroff});
    vz_vector.push_back(vertexz);
    phpos_vector.push_back((int)posPhoton);	  
      
	}
    }

    // Safe mixing
/*    if(do_mixing && !track_4vector.empty()){
        std::cout << "Time aor mixing" << std::endl;
        MixEvents(mincentormult, Nmixevents, multiplicity_vector, vz_vector, phpos_vector, trigger_vector,
                  minvz, track_4vector, track_charge_vector, track_weights_vector,
                  hist_qinv_SS_MIX, hist_qlcms_SS_MIX, hist_q3D_SS_MIX,
                  hist_qinv_OS_MIX, hist_qlcms_OS_MIX, hist_q3D_OS_MIX,
                  dosplit, do_hbt3d, do_gamov, syst, NeventsAss);
    }
*/

if(do_mixing) {
    if(track_4vector.empty()) {
        std::cout << "Warning: no events to mix!" << std::endl;
    } else {
        std::cout << "Time for mixing" << std::endl;
        /*MixEvents(mincentormult, Nmixevents, multiplicity_vector, vz_vector, phpos_vector, trigger_vector,
                  minvz, track_4vector, track_charge_vector, track_weights_vector,
                  hist_qinv_SS_MIX, hist_qlcms_SS_MIX, hist_q3D_SS_MIX,
                  hist_qinv_OS_MIX, hist_qlcms_OS_MIX, hist_q3D_OS_MIX,
                  dosplit, do_hbt3d, do_gamov, syst, NeventsAss);
    */
          //MixEvents_eta_random(mincentormult, Nmixevents, multiplicity_vector, vz_vector, phpos_vector, trigger_vector, minvz, track_4vector, track_charge_vector, track_weights_vector, hist_qinv_SS_MIX, hist_qlcms_SS_MIX, hist_q3D_SS_MIX, hist_qinv_OS_MIX, hist_qlcms_OS_MIX, hist_q3D_OS_MIX, dosplit, do_hbt3d, false, syst, NeventsAss);

for(int ich = 0; ich < nCentBinsF; ich++){
    MixEvents_eta_random(
        (int)CentBinsF[ ich ],       // ntrkoff_min
        (int)CentBinsF[ ich + 1 ] - 1,   // ntrkoff_max
        mincentormult,             // centrality_or_ntrkoff_int
        Nmixevents,                // nEvt_to_mix
        track_4vector,             // ev_GoodTrackFourVector_etaMixWeight_vec
        track_charge_vector,       // ev_GoodTrackCharge_etaMixWeight_vec
        track_weights_vector,      // ev_GoodTrackEff_etaMixWeight_vec
        multiplicity_vector,       // ev_ntrkoff_etaMixWeight_vec
        vz_vector,                 // vtx_z_vec
        phpos_vector,              // ev_phpos
        trigger_vector,            // ev_trg
        minvz,                     // vzcut
        hist_qinv_SS_MIX,          // histo_SS
        hist_qlcms_SS_MIX,         // histo_SSLCMS
        hist_q3D_SS_MIX,           // histo_SS3D
        hist_qinv_OS_MIX,          // histo_OS
        hist_qlcms_OS_MIX,         // histo_OSLCMS
        hist_q3D_OS_MIX,           // histo_OS3D
        dosplit,                   // docostdptcut
        do_hbt3d,                  // do_hbt3d
        false,                     // dogamovcorrection
        syst,                      // systematic
        NeventsAss                 // NeventsAss
    );
}


		  }
}


    // Output
    
// Split path and filename
string full_path = ouputfile.Data();
string dir_path = full_path.substr(0, full_path.find_last_of('/') + 1);
string base_name = full_path.substr(full_path.find_last_of('/') + 1);

string file_output = Form("%s_syst_%s_Nmix_%i_Mixint_%i_Vzint_%.f_on_%i_corr_pT0p3",
                          base_name.c_str(), systematics.Data(), Nmixevents,
                          mincentormult, minvz, date->GetDate());
std::replace(file_output.begin(), file_output.end(), '.', 'p');

// Prepend the directory path back
file_output = dir_path + file_output;

TFile *MyFile = new TFile(Form("%s.root", file_output.c_str()), "RECREATE");

//string file_output = Form("%s_syst_%s_Nmix_%i_Mixint_%i_Vzint_%.f_on_%i_corr_pT0p2",
//                              ouputfile.Data(), systematics.Data(), Nmixevents,
//                              mincentormult, minvz, date->GetDate());
//    std::replace(file_output.begin(), file_output.end(), '.', 'p');

//    TFile *MyFile = new TFile(Form("%s.root", file_output.c_str()), "RECREATE");
    if(MyFile->IsOpen()) std::cout << "Output file: " << file_output << ".root" << std::endl;

    // Write histograms
    MyFile->mkdir("EventQA_histograms"); MyFile->cd("EventQA_histograms"); write_eventQA();
    MyFile->mkdir("HBT_1D"); MyFile->cd("HBT_1D"); write_HBT1D();
    if(do_hbt3d){ MyFile->mkdir("HBT_3D"); MyFile->cd("HBT_3D"); write_HBT3D(); }

    MyFile->Close();
    std::cout << "Job finished successfully. Output: " << file_output << std::endl;

    sec_end = clock();
    std::cout << "Total running time: " << (double)(sec_end - sec_start)/CLOCKS_PER_SEC << " [s]" << std::endl;
}

