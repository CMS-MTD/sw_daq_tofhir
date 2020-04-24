//================================================================================================================
// Done by abdollah.mohammadi@cern.ch and Si.Xie@cern.ch
//================================================================================================================
// include std libraries
#include <iostream>
#include <fstream>
#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <cstring>
#include <string.h>
#include <sstream>
// include ROOT libraries 
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TTree.h"
#include "TChain.h"
#include "TGraph.h"
#include "TStyle.h"
#include "TFolder.h"
#include "TCanvas.h"
#include "TRandom.h"
#include "TMath.h"
#include "TFile.h"
#include "TSystem.h"
#include "TProfile.h"

#include "Tofhir.h"
//#include "Tracker.h"
#include "TrackerSlow.h"


#define MAX_TOFPET_CHANNEL 400
#define SEARCHWINDOWSIZE 15000 //in ps units  (5000 was changed to 15000 after an optimization study)
#define MATCH_THRESHOLD 0.000005 // in seconds
#define DEBUG 0
#define ClockScaleFactor 1.04617

using namespace std;

int main(int argc, char* argv[]){
//    if(argc < 3) {
//        cerr << "Please give 3 arguments " << "inputTOFHIRFileName " << " inputTrackerFileName " << "outputFileName" <<endl;
//        return -1;
//    }
    if(argc < 1) {
        cerr << "Please give the run number " <<endl;
        return -1;
    }

    const char *Run = argv[1];
    std::string runNumber(Run);

    
    //=====================Open Tracker input files=========================
//    const char *inputTrackerFileName = argv[2];
    string trackingFile = "/home/daq/2019_04_April_CMSTiming/Tracks/Run"+runNumber+"_CMSTiming_SlowTriggerStream_converted.root";
    TFile *trackerFile = new TFile(trackingFile.c_str(),"READ");
    TTree *trackerTreeSlow = (TTree*)trackerFile->Get("CMSTiming");
    
    if( trackerTreeSlow != NULL ) cout << ">>> got track tree from file " << trackingFile << endl;
    else exit(-1);

//    =====================Open Slow Tracker input files=========================

    SlowTrigTracker TRK_Slow(trackerTreeSlow);


        for(int k=0; k<trackerTreeSlow->GetEntries(); k++) {

            trackerTreeSlow->GetEntry( 0 );
            Long64_t TrackerBCOThisEventZero = TRK_Slow.trackerEventSlow->bco;

            trackerTreeSlow->GetEntry( k );
            Long64_t TrackerBCOThisEvent = TRK_Slow.trackerEventSlow->bco;

            if (TRK_Slow.trackerEventSlow->bco== 4294967295) continue;

            if (1) cout<<"k = "<<k<<"  time is= "<<  TrackerBCOThisEvent  * 144e-9 * ClockScaleFactor <<" - "<< TrackerBCOThisEventZero * 144e-9 * ClockScaleFactor <<" ---> " <<(TrackerBCOThisEvent - TrackerBCOThisEventZero ) * 144e-9 * ClockScaleFactor<< "\n";

        }

    }
