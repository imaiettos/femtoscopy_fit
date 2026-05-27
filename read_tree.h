#include "call_libraries.h"  // call libraries from ROOT and C++

// declare variables
float zdcSumPlus, zdcSumMinus;
int mult, npixcl;
float vtx[3];
bool triggers[8];
vector<int> *pfID = 0;
vector<int> *pfCharge = 0;
vector<bool> *isPrimary = 0;
vector<float> *pfPhi = 0;
vector<float> *pfEta = 0;
vector<float> *pfPt = 0;
vector<float> *pfE = 0;
vector<char> *trkNPixelHits = 0;


//All variables listed above are readed in the function bellow
/*
Function to read the Forest/Skim tree
Arguments ->  transfer quantities from trees to our variables
tree: input TChain from jet_analyzer.C file
is_MC: true -> MC; false -> Data
*/
void read_tree(TChain *tree){

    tree->SetBranchStatus("*", 0); // disable all branches - this is important while reading big files

//    tree->SetBranchStatus( "nPixelClusters"  , 1);
    tree->SetBranchStatus( "zdcSumPlus"	     , 1);
    tree->SetBranchStatus( "zdcSumMinus"	 , 1);
    tree->SetBranchStatus( "multiplicity"    , 1);
    tree->SetBranchStatus( "vtx"			 , 1);
    tree->SetBranchStatus( "triggersFired"   , 1);
    tree->SetBranchStatus( "pfCharge"                         , 1);
    tree->SetBranchStatus( "pfID"			 , 1);
    tree->SetBranchStatus( "isPrimary"	  	 , 1);
    tree->SetBranchStatus( "trkNPixelHits"           , 1);
    tree->SetBranchStatus( "pfEta"			 , 1);
    tree->SetBranchStatus( "pfPhi"		 	 , 1);
    tree->SetBranchStatus( "pfPt"			 , 1);
    tree->SetBranchStatus( "pfE"			 , 1);

//    tree->SetBranchAddress( "nPixelClusters"  ,   &npixcl );
    tree->SetBranchAddress( "zdcSumPlus"	  ,   &zdcSumPlus );
    tree->SetBranchAddress( "zdcSumMinus"	  ,   &zdcSumMinus );
    tree->SetBranchAddress( "multiplicity"    ,   &mult );
    tree->SetBranchAddress( "vtx"			  ,   vtx );
    tree->SetBranchAddress( "triggersFired"   ,   triggers );
    tree->SetBranchAddress( "pfCharge"                         , &pfCharge);
    tree->SetBranchAddress( "pfID"			  ,   &pfID );
    tree->SetBranchAddress( "isPrimary"	  	  ,   &isPrimary );
    tree->SetBranchAddress( "trkNPixelHits"           ,   &trkNPixelHits );
    tree->SetBranchAddress( "pfEta"			  ,   &pfEta );
    tree->SetBranchAddress( "pfPhi"		 	  ,   &pfPhi );
    tree->SetBranchAddress( "pfPt"			  ,   &pfPt );
    tree->SetBranchAddress( "pfE"			  ,   &pfE );

}
