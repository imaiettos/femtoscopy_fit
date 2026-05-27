#include "call_libraries.h"  // call libraries from ROOT and C++
//#define coscut 0.99996
//#define dptcut 0.04
#define coscut 0.99991
//#define coscut_i 0.999985
#define dptcut 0.045
//#define dptcut_i 0.055598
#define pimass 0.1396

/*
Calculate q invariant
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
cos_cut: cut used in the cosine of angle between the particles
dpt_cut: cut used in the difference in pT of the particles
*/
bool splitcomb(ROOT::Math::PtEtaPhiMVector &vec1,ROOT::Math::PtEtaPhiMVector &vec2, double cos_cut, double dpt_cut){
//bool splitcomb(ROOT::Math::PtEtaPhiMVector &vec1,ROOT::Math::PtEtaPhiMVector &vec2, double cos_cut, double dpt_cut, double cos_cut_i, double dpt_cut_i){
   bool issplit=false;
   double cosa = TMath::Abs(vec1.Px()*vec2.Px() + vec1.Py()*vec2.Py() + vec1.Pz()*vec2.Pz())/(vec1.P()*vec2.P());
   double deltapt = fabs(vec1.Pt() - vec2.Pt());
//   if ( ((cosa > cos_cut) && (deltapt < dpt_cut)) || ((cosa > cos_cut_i) && (deltapt < dpt_cut_i)) ) { issplit = true;}
   if  ((cosa > cos_cut) && (deltapt < dpt_cut)) { issplit = true;}  
  
   return issplit;
}

/*
Calculate q long in the LCMS
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
float GetQlongLCMS(ROOT::Math::PtEtaPhiMVector &p1, ROOT::Math::PtEtaPhiMVector &p2){
   Double_t num = 2.0*( (p1.Pz())*(p2.E()) - (p2.Pz())*(p1.E()) );
   Double_t den = TMath::Sqrt( (p1.E()+p2.E())*(p1.E()+p2.E()) - (p1.Pz()+p2.Pz())*(p1.Pz()+p2.Pz()) );
   Double_t qlongLCMS = 0.0;
   if(den != 0) qlongLCMS = fabs(num/den);
   return qlongLCMS;
}

/*
Calculate q invariant
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
float GetQ(ROOT::Math::PtEtaPhiMVector &p1, ROOT::Math::PtEtaPhiMVector &p2){
   	ROOT::Math::PtEtaPhiMVector Sum4V = p1+p2;
   	Double_t q = Sum4V.M2() - 4.0*p1.mass()*p2.mass();
   	return (  q > 0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}

/*
Calculate q invariant
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
float GetQLCMS(ROOT::Math::PtEtaPhiMVector &p1, ROOT::Math::PtEtaPhiMVector &p2){
	Double_t qlong = GetQlongLCMS(p1,p2);
	Double_t q = (p1.Px()-p2.Px())*(p1.Px()-p2.Px()) + (p1.Py()-p2.Py())*(p1.Py()-p2.Py()) + qlong*qlong;
   	return (  q > 0 ?  TMath::Sqrt(q) : -TMath::Sqrt(-q)  );
}

/*
Calculate q out in the LCMS
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
float GetQout(ROOT::Math::PtEtaPhiMVector &p1, ROOT::Math::PtEtaPhiMVector &p2){
   TVector3 qT;
   qT.SetXYZ(p1.Px()-p2.Px(),p1.Py()-p2.Py(),0.0);
   TVector3 kT;
   kT.SetXYZ( (p1.Px()+p2.Px())/2.0 , (p1.Py()+p2.Py())/2.0 ,0.0);
   TVector3 qout;
   qout = qT.Dot(kT.Unit())*kT.Unit();
   Double_t absValue = qout.Mag();
   return absValue; 
}

/*
Calculate q side in the LCMS
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
float GetQside(ROOT::Math::PtEtaPhiMVector &p1, ROOT::Math::PtEtaPhiMVector &p2){
   TVector3 qT;
   qT.SetXYZ(p1.Px()-p2.Px(),p1.Py()-p2.Py(),0.0);
   TVector3 kT;
   kT.SetXYZ( (p1.Px()+p2.Px())/2.0 , (p1.Py()+p2.Py())/2.0 ,0.0);
   TVector3 qout;
   qout = qT.Dot(kT.Unit())*kT.Unit();
   TVector3 qsid;
   qsid = qT - qout;
   Double_t absValue = qsid.Mag();
   return absValue;
}

/*
Invert Px, Py and Pz for reference sample
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
ROOT::Math::PtEtaPhiMVector InvertPVector( ROOT::Math::PtEtaPhiMVector &vec){
   ROOT::Math::PtEtaPhiMVector ovec = vec;
   ovec.SetPxPyPzE(-vec.Px(),-vec.Py(),-vec.Pz(),vec.E());
   return ovec;
}

/*
Rotate X and Y for reference sample
--> Arguments
p1: particle 1 4-vector
p2: particle 2 4-vector
*/
ROOT::Math::PtEtaPhiMVector InvertXYVector( ROOT::Math::PtEtaPhiMVector &vec){
  ROOT::Math::PtEtaPhiMVector ovec = vec;
  ovec.SetXYZT(-vec.X(),-vec.Y(),vec.Z(),vec.T());
  return ovec;
}

/*
Return the weight factor due to Coloumb repulsion [Gamow] same charge
--> Arguments
q: q invariant
systematic: systematic number
*/
const double CoulombSS(const double &q, int systematic){
   if(q == 0) return 0.0;
   const double alpha = 1./137.;
   double x = 1.0*TMath::Pi()*(alpha*pimass/q);
   double weight = 1.0;
   if(systematic == 9) weight = 1.15; //for syst. +15%
   if(systematic == 10) weight = 0.85; //for syst. -15%
   return weight*( (TMath::Exp(x) - 1.0)/ x - 1.0 ) + 1.0;
}

/*
Return the weight factor due to Coloumb attraction [Gamow] opposite charge
--> Arguments
q: q invariant
systematic: systematic number
*/
const double CoulombOS(const double &q, int systematic){
   if(q == 0) return 0.0;
   const double alpha = 1./137.;
   double x = 1.0*TMath::Pi()*(alpha*pimass/q);
   double weight = 1.0;
   if(systematic == 9) weight = 1.15; //for syst. +15%
   if(systematic == 10) weight = 0.85; //for syst. -15%
   return weight*( (1.-TMath::Exp(-x)) / x - 1.0 ) + 1.0;
}

/*
Find the eta gap and total energy
--> Arguments
syst: for systematics
PF info: pT, Eta, E and ID
isPrimary?
etaGapPos: eta gap positive
etaGapNeg: eta gap negative
2D histograms for multiple PF particles
*/
void etagaps(int syst, double hfthreshold, double &ETSum, int &numPosHFClusters, int &numNegHFClusters, std::vector<float>* pfPt, std::vector<float>* pfEta, std::vector<float>* pfE, std::vector<int>* pfID, std::vector<bool>* isPrimary, double &etaGapPos, double &etaGapNeg, TH2D* h_h,  TH2D* h_e, TH2D* h_mu, TH2D* h_gamma, TH2D* h_h0, TH2D* h_HFhad, TH2D* h_HFem, bool posPhoton){

      int numPF = pfID->size();	  
      vector<float> pfEtaPosVector;
      vector<float> pfEtaNegVector;
      //double ETSum = 0;
      float multCorr = 0;      
      for( int iPF = 0; iPF < numPF; iPF++ ) {
        double Et = pfE->at( iPF )/TMath::CosH( pfEta->at( iPF ) );

	// Find et sum on Pb going side
	if( pfID->at( iPF ) == 6 || pfID->at( iPF ) == 7 ){
        if( ((!posPhoton && pfEta->at( iPF ) > 0) || (posPhoton && pfEta->at( iPF ) < 0)) && pfE->at( iPF ) > 3.0 ) ETSum += Et;
		if( pfE->at( iPF ) > hfthreshold && pfEta->at( iPF ) > 0 ) numPosHFClusters++;
		if( pfE->at( iPF ) > hfthreshold && pfEta->at( iPF ) < 0 ) numNegHFClusters++;
		}
		
        if( pfID->at( iPF ) == 1 ) {      //h+/-
          h_h->Fill( pfEta->at( iPF ), pfPt->at( iPF ) );
        }else if( pfID->at( iPF ) == 2 ) { //e
          h_e->Fill( pfEta->at( iPF ), pfPt->at( iPF ) );
        }else if( pfID->at( iPF ) == 3 ) { //mu
          h_mu->Fill( pfEta->at( iPF ), pfPt->at( iPF ) );
        }else if( pfID->at( iPF ) == 4 ) { //gamma
          h_gamma->Fill( pfEta->at( iPF ), Et );
        }else if( pfID->at( iPF ) == 5 ) { //h0
          h_h0->Fill( pfEta->at( iPF ), Et );
        }else if( pfID->at( iPF ) == 6 ) { //Had HF
          h_HFhad->Fill( pfEta->at( iPF ), pfE->at( iPF ) ); 
        }else if( pfID->at( iPF ) == 7 ) { //EM HF
          h_HFem->Fill( pfEta->at( iPF ), pfE->at( iPF ) ); }
        if(syst == 15){ 
	        //Lower threshold cuts
    	    //if( ( pfID->at( iPF ) == 1 && pfPt->at( iPF ) < 0.3 && isPrimary->at( iPF ) ) || //h+/-
        if( ( pfID->at( iPF ) == 1 && (pfPt->at( iPF ) < 0.3 || !isPrimary->at( iPF )) ) || //h+/-	
		( pfID->at( iPF ) == 2 && pfPt->at( iPF ) < 1.0 ) || //e
     	        ( pfID->at( iPF ) == 3 && pfPt->at( iPF ) < 0.5 ) || //mu
     	        ( pfID->at( iPF ) == 4 && Et < 0.5 ) || //gamma
     	        ( pfID->at( iPF ) == 5 && Et < 2.0 && fabs( pfEta->at( iPF ) ) < 1.4 ) || //h0 barrel
     	        ( pfID->at( iPF ) == 5 && Et < 0.5 && fabs( pfEta->at( iPF ) ) > 1.4 ) || //h0 endcaps
      	        ( pfID->at( iPF ) == 6 && pfE->at( iPF ) < 2.0 ) || //Had HF
     	        ( pfID->at( iPF ) == 7 && pfE->at( iPF ) < 2.0 ) ) {  //EM HF
     	     		continue; }
         }else if (syst == 16){
	        //High threshold cuts
            	        if( ( pfID->at( iPF ) == 1 && (pfPt->at( iPF ) < 0.5 || !isPrimary->at( iPF )) ) || //h+/-
		( pfID->at( iPF ) == 2 && pfPt->at( iPF ) < 3.0 ) || //e
            	( pfID->at( iPF ) == 3 && pfPt->at( iPF ) < 2.0 ) || //mu
            	( pfID->at( iPF ) == 4 && Et < 1.5 ) || //gamma
            	( pfID->at( iPF ) == 5 && Et < 4.0 && fabs( pfEta->at( iPF ) ) < 1.4 ) || //h0 barrel
            	( pfID->at( iPF ) == 5 && Et < 1.5 && fabs( pfEta->at( iPF ) ) > 1.4 ) || //h0 endcaps
            	( pfID->at( iPF ) == 6 && pfE->at( iPF ) < 4.0 ) || //Had HF
           		( pfID->at( iPF ) == 7 && pfE->at( iPF ) < 4.0 ) ) {  //EM HF
          			continue; }
         }else{
	        //Default cuts
    	    //if( ( pfID->at( iPF ) == 1 && pfPt->at( iPF ) < 0.4 ) || //h+/-
        if( ( pfID->at( iPF ) == 1 && (pfPt->at( iPF ) < 0.4 || !isPrimary->at( iPF )) ) || //h+/- 	
	        ( pfID->at( iPF ) == 2 && pfPt->at( iPF ) < 2.0 ) || //e
            	( pfID->at( iPF ) == 3 && pfPt->at( iPF ) < 1.0 ) || //mu
            	( pfID->at( iPF ) == 4 && Et < 1.0 ) || //gamma
            	( pfID->at( iPF ) == 5 && Et < 3.0 && fabs( pfEta->at( iPF ) ) < 1.4 ) || //h0 barrel
            	( pfID->at( iPF ) == 5 && Et < 1.0 && fabs( pfEta->at( iPF ) ) > 1.4 ) || //h0 endcaps
            	( pfID->at( iPF ) == 6 && pfE->at( iPF ) < 3.0 ) || //Had HF
            	( pfID->at( iPF ) == 7 && pfE->at( iPF ) < 3.0 ) ) {  //EM HF
		          	continue; }
		}
        if( pfEta->at( iPF ) > 0 ){ pfEtaPosVector.push_back( pfEta->at( iPF ) );
        }else { pfEtaNegVector.push_back( pfEta->at( iPF ) ); }

        //Find ETSum of HF on photon-going side and see if HF cluster passes threshold
        if( pfID->at( iPF ) == 6 || pfID->at( iPF ) == 7 )
        { 
          if( ((!posPhoton && pfEta->at( iPF ) > 0) || (posPhoton && pfEta->at( iPF ) < 0))   
              && pfE->at( iPF ) > 3.0 )
            ETSum += Et;
        }



      }


      pfEtaPosVector.push_back( 0.0 );//Add middle of detector as a particle
      pfEtaPosVector.push_back( 5.3 );//Add edge of detector as a particle
      sort( pfEtaPosVector.begin(), pfEtaPosVector.end() );
      for( int i = 1; i < pfEtaPosVector.size(); i++ ) {
        double deltaEta = fabs( pfEtaPosVector[i] - pfEtaPosVector[i - 1] );
        if( deltaEta > 0.5 ) etaGapPos += deltaEta;
      }
      pfEtaPosVector.clear();

      pfEtaNegVector.push_back( 0.0 );//Add middle of detector as a particle
      pfEtaNegVector.push_back( -5.3 );//Add edge of detector as a particle
      sort( pfEtaNegVector.begin(), pfEtaNegVector.end() );
      for( int i = 1; i < pfEtaNegVector.size(); i++ ){
        double deltaEta = fabs( pfEtaNegVector[i] - pfEtaNegVector[i - 1] );
        if( deltaEta > 0.5 ) etaGapNeg += deltaEta;
      }
      pfEtaNegVector.clear();

} 


/*
Measure the 2 particle correlation
--> Arguments
tracks: vector with track informations
tracks_charge: vector with track charge informations
tracks_weight: vector with track efficiency informations
*/
void twoparticlecorrelation(std::vector<ROOT::Math::PtEtaPhiMVector> tracks, std::vector<int> tracks_charge, std::vector<double> tracks_weight, TH1D* pairmass_samesign, TH2D* costhetadpt_samesign, TH2D* dphideta_samesign, THnSparse* histo_2pcorr_samesign,  THnSparse* histo_2pcorr_samesign_inverted,  THnSparse* histo_2pcorr_samesign_rotated, THnSparse* histo_2pcorrlcms_samesign,  THnSparse* histo_2pcorrlcms_samesign_inverted,  THnSparse* histo_2pcorrlcms_samesign_rotated, THnSparse* histo_2pcorr_samesign3D,  THnSparse* histo_2pcorr_samesign3D_inverted,  THnSparse* histo_2pcorr_samesign3D_rotated, TH1D* pairmass_oppsign, TH2D* costhetadpt_oppsign, TH2D* dphideta_oppsign, THnSparse* histo_2pcorr_oppsign, THnSparse* histo_2pcorr_oppsign_inverted, THnSparse* histo_2pcorr_oppsign_rotated, THnSparse* histo_2pcorrlcms_oppsign,  THnSparse* histo_2pcorrlcms_oppsign_inverted,  THnSparse* histo_2pcorrlcms_oppsign_rotated, THnSparse* histo_2pcorr_oppsign3D, THnSparse* histo_2pcorr_oppsign3D_inverted, THnSparse* histo_2pcorr_oppsign3D_rotated, int cent, int posPhoton, bool docostdptcut, bool do_hbt3d, bool dogamovcorrection, int systematic){
	// get correlation histograms
	for (int ipair = 0; ipair < (tracks.size()*tracks.size()); ipair++){ // start loop over tracks
	
        int a = ipair / tracks.size();
        int b = ipair % tracks.size();
        if( b <= a ) continue; // This makes it correlate like: b = a + 1;  

		double eff_trk_a = tracks_weight[a];
		double eff_trk_b = tracks_weight[b];
		double tot_eff = eff_trk_a*eff_trk_b;

		if(docostdptcut){if(splitcomb(tracks[a],tracks[b],coscut,dptcut)) continue;}

        ROOT::Math::PtEtaPhiMVector psum2 = tracks[a] + tracks[b];
        double kt = (psum2.Pt())/2.;
        ROOT::Math::PtEtaPhiMVector psum2_inverted = tracks[a] + InvertPVector(tracks[b]);
        double kt_inverted = (psum2_inverted.Pt())/2.;
        ROOT::Math::PtEtaPhiMVector psum2_rotated = tracks[a] + InvertXYVector(tracks[b]);
        double kt_rotated = (psum2_rotated.Pt())/2.;

		ROOT::Math::PtEtaPhiMVector trackb_inverted = InvertPVector(tracks[b]);
		ROOT::Math::PtEtaPhiMVector trackb_rotated = InvertXYVector(tracks[b]);

		double qinv = GetQ(tracks[a],tracks[b]);
		double qinv_inverted = GetQ(tracks[a], trackb_inverted);
		double qinv_rotated = GetQ(tracks[a], trackb_rotated);
		double x_2pc_hbt[4]={qinv, kt, (double)cent, (double) posPhoton}; 
		double x_2pc_hbt_inv[4]={qinv_inverted, kt_inverted, (double)cent, (double) posPhoton}; 
		double x_2pc_hbt_rot[4]={qinv_rotated, kt_rotated, (double)cent, (double) posPhoton}; 

		double qlcms = GetQLCMS(tracks[a],tracks[b]);
		double qlcms_inverted = GetQLCMS(tracks[a], trackb_inverted);
		double qlcms_rotated = GetQLCMS(tracks[a], trackb_rotated);
		double x_2pc_hbtlcms[4]={qlcms, kt, (double)cent, (double) posPhoton};
		double x_2pc_hbtlcms_inv[4]={qlcms_inverted, kt_inverted, (double)cent, (double) posPhoton};
		double x_2pc_hbtlcms_rot[4]={qlcms_rotated, kt_rotated, (double)cent, (double) posPhoton};

		double qlong = GetQlongLCMS(tracks[a],tracks[b]);
		double qlong_inverted = GetQlongLCMS(tracks[a], trackb_inverted);
		double qlong_rotated = GetQlongLCMS(tracks[a], trackb_rotated);
		double qout = GetQout(tracks[a],tracks[b]);
		double qout_inverted = GetQout(tracks[a], trackb_inverted);
		double qout_rotated = GetQout(tracks[a], trackb_rotated);
		double qside = GetQside(tracks[a],tracks[b]);
		double qside_inverted = GetQside(tracks[a], trackb_inverted);
		double qside_rotated = GetQside(tracks[a], trackb_rotated);
		double x_2pc_hbt_3D[6]={qlong, qout, qside, kt, (double)cent, (double) posPhoton}; 
		double x_2pc_hbt_3D_inv[6]={qlong_inverted, qout_inverted, qside_inverted, kt_inverted, (double)cent, (double) posPhoton}; 
		double x_2pc_hbt_3D_rot[6]={qlong_rotated, qout_rotated, qside_rotated, kt_rotated, (double)cent, (double) posPhoton}; 
			
		double coulomb_ss = 1.0;
		double coulomb_os = 1.0;
		if(dogamovcorrection){
			coulomb_ss = CoulombSS(qinv,systematic);
			coulomb_os = CoulombOS(qinv,systematic);		
		}
		
		//if(fabs(tracks[a].Eta() - tracks[b].Eta()) == 0 && fabs(tracks[a].Phi() - tracks[b].Phi()) == 0) continue;
		double deltaetaabs = fabs(tracks[a].Eta() - tracks[b].Eta());
		double deltaphiabs = fabs(TVector2::Phi_mpi_pi(tracks[a].Phi() - tracks[b].Phi()));
		double costheta = TMath::Abs(tracks[a].Px()*tracks[b].Px() + tracks[a].Py()*tracks[b].Py() + tracks[a].Pz()*tracks[b].Pz())/(tracks[a].P()*tracks[b].P());
   		double dpt = fabs(tracks[a].Pt() - tracks[b].Pt());
   		double pairmass = psum2.M();

		if(tracks_charge[a]*tracks_charge[b] > 0){
			pairmass_samesign->Fill(pairmass,tot_eff);
			costhetadpt_samesign->Fill(costheta,dpt,tot_eff);
			dphideta_samesign->Fill(deltaetaabs, deltaphiabs);
			histo_2pcorr_samesign->Fill(x_2pc_hbt,coulomb_ss*tot_eff);
			histo_2pcorr_samesign_inverted->Fill(x_2pc_hbt_inv,coulomb_ss*tot_eff);
			histo_2pcorr_samesign_rotated->Fill(x_2pc_hbt_rot,coulomb_ss*tot_eff);
			histo_2pcorrlcms_samesign->Fill(x_2pc_hbtlcms,coulomb_ss*tot_eff);
			histo_2pcorrlcms_samesign_inverted->Fill(x_2pc_hbtlcms_inv,coulomb_ss*tot_eff);
			histo_2pcorrlcms_samesign_rotated->Fill(x_2pc_hbtlcms_rot,coulomb_ss*tot_eff);
			if(do_hbt3d){
				histo_2pcorr_samesign3D->Fill(x_2pc_hbt_3D,coulomb_ss*tot_eff);
				histo_2pcorr_samesign3D_inverted->Fill(x_2pc_hbt_3D_inv,coulomb_ss*tot_eff);
				histo_2pcorr_samesign3D_rotated->Fill(x_2pc_hbt_3D_rot,coulomb_ss*tot_eff);
			}
		}else{
			pairmass_oppsign->Fill(pairmass,tot_eff);
			costhetadpt_oppsign->Fill(costheta,dpt,tot_eff);
			dphideta_oppsign->Fill(deltaetaabs, deltaphiabs);
			histo_2pcorr_oppsign->Fill(x_2pc_hbt,coulomb_os*tot_eff);			
			histo_2pcorr_oppsign_inverted->Fill(x_2pc_hbt_inv,coulomb_os*tot_eff);			
			histo_2pcorr_oppsign_rotated->Fill(x_2pc_hbt_rot,coulomb_os*tot_eff);			
			histo_2pcorrlcms_oppsign->Fill(x_2pc_hbtlcms,coulomb_os*tot_eff);
			histo_2pcorrlcms_oppsign_inverted->Fill(x_2pc_hbtlcms_inv,coulomb_os*tot_eff);
			histo_2pcorrlcms_oppsign_rotated->Fill(x_2pc_hbtlcms_rot,coulomb_os*tot_eff);
			if(do_hbt3d){
				histo_2pcorr_oppsign3D->Fill(x_2pc_hbt_3D,coulomb_os*tot_eff);			
				histo_2pcorr_oppsign3D_inverted->Fill(x_2pc_hbt_3D_inv,coulomb_os*tot_eff);			
				histo_2pcorr_oppsign3D_rotated->Fill(x_2pc_hbt_3D_rot,coulomb_os*tot_eff);			
			}
		}
	} // 2PC loop
}


