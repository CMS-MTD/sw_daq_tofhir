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
#include "Tracker.h"
//#include "TrackerSlow.h"


#define MAX_TOFPET_CHANNEL 400
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


    //=====================Open TOFHIR input files=========================
//    const char *inputTofhirFileName = argv[1];
    string singleFile="/home/daq/2020_02_cmstiming_BTL/TOFHIR/RecoData/v1/RecoWithTracks/run"+runNumber+"_singles.root";
  // "/home/daq/2019_04_April_CMSTiming/TOFHIR/RecoData/v1/RecoWithTracks/run"+runNumber+"_singles.root";
    TFile *tofhirFile=new TFile(singleFile.c_str(),"read");
    TTree *tofhirTree = (TTree*)tofhirFile->Get("data");
    if( tofhirTree != NULL ) cout << "\n>>> got TOFHIR tree from file " << singleFile << endl;
    else exit(-1);
    
    
    TOFHIR TOF_(tofhirTree);
    
    TOFHIR::SpillInfo  Ev1stSpill1st = TOF_.get1stEvent1stSplill(0);
    TOFHIR::SpillInfo  EvLastSpill1st =  TOF_.getLastEvent1stSplill(Ev1stSpill1st.Index );
    TOFHIR::SpillInfo  Ev1stSpill2nd =  TOF_.get1stEvent2ndSplill(EvLastSpill1st.Index);
    TOFHIR::SpillInfo  EvLastSpill2nd =  TOF_.getLastEvent2ndSplill(Ev1stSpill2nd.Index);
    TOFHIR::SpillInfo  Ev1stSpill3rd =  TOF_.get1stEvent3rdSplill(EvLastSpill2nd.Index);
    
    
    cout << "TOFHIR: first event's index and time in spill 1: " << Ev1stSpill1st.Index << " : " << Ev1stSpill1st.Time << "\n";
    cout << "TOFHIR: last event's index and time in spill 1:  " << EvLastSpill1st.Index << " : " << EvLastSpill1st.Time << "\n";
    cout << "TOFHIR: first event's index and time in spill 2: " << Ev1stSpill2nd.Index << " : " << Ev1stSpill2nd.Time << "\n";
    cout << "TOFHIR: last event's index and time in spill 2:  " << EvLastSpill2nd.Index << " : " << EvLastSpill2nd.Time << "\n";
    cout << "TOFHIR: first event's index and time in spill 3: " << Ev1stSpill3rd.Index << " : " << Ev1stSpill3rd.Time << "\n";

    
    
    //=====================Open Tracker input files=========================
//    const char *inputTrackerFileName = argv[2];
    string trackingFile = "/home/daq/2019_04_April_CMSTiming/Tracks/Run"+runNumber+"_CMSTiming_FastTriggerStream_converted.root";
    TFile *trackerFile = new TFile(trackingFile.c_str(),"READ");
    TTree *trackerTree = (TTree*)trackerFile->Get("CMSTiming");
    
    if( trackerTree != NULL ) cout << ">>> got track tree from file " << trackingFile << endl;
    else exit(-1);


    TRACKER TRK_(trackerTree);
    
    TRACKER::SpillInfo  Ev1stSpill1st_trk = TRK_.get1stEvent1stSplill(trackerTree, 0);
    TRACKER::SpillInfo  EvLastSpill1st_trk =  TRK_.getLastEvent1stSplill(trackerTree, Ev1stSpill1st_trk.Index );
    TRACKER::SpillInfo  Ev1stSpill2nd_trk =  TRK_.get1stEvent2ndSplill(trackerTree, EvLastSpill1st_trk.Index);
    TRACKER::SpillInfo  EvLastSpill2nd_trk =  TRK_.getLastEvent2ndSplill(trackerTree, Ev1stSpill2nd_trk.Index);
    TRACKER::SpillInfo  Ev1stSpill3rd_trk =  TRK_.get1stEvent3rdSplill(trackerTree, EvLastSpill2nd_trk.Index);


    cout << "TRACKER: first event's index and time in spill 1: " << Ev1stSpill1st_trk.Index << " : " << Ev1stSpill1st_trk.Time << "\n";
    cout << "TRACKER: last event's index and time in spill 1:  " << EvLastSpill1st_trk.Index << " : " << EvLastSpill1st_trk.Time << "\n";
    cout << "TRACKER: first event's index and time in spill 2: " << Ev1stSpill2nd_trk.Index << " : " << Ev1stSpill2nd_trk.Time << "\n";
    cout << "TRACKER: last event's index and time in spill 2:  " << EvLastSpill2nd_trk.Index << " : " << EvLastSpill2nd_trk.Time << "\n";
    cout << "TRACKER: first event's index and time in spill 3: " << Ev1stSpill3rd_trk.Index << " : " << Ev1stSpill3rd_trk.Time << "\n";


    //=====================Open Slow Tracker input files=========================

//    SlowTrigTracker TRK_(trackerTree);
//
//    SlowTrigTracker::SpillInfo  Ev1stSpill1st_trk = TRK_.get1stEvent1stSplill(trackerTree, 0);
//    SlowTrigTracker::SpillInfo  EvLastSpill1st_trk =  TRK_.getLastEvent1stSplill(trackerTree, Ev1stSpill1st_trk.Index );
//    SlowTrigTracker::SpillInfo  Ev1stSpill2nd_trk =  TRK_.get1stEvent2ndSplill(trackerTree, EvLastSpill1st_trk.Index);
//    SlowTrigTracker::SpillInfo  EvLastSpill2nd_trk =  TRK_.getLastEvent2ndSplill(trackerTree, Ev1stSpill2nd_trk.Index);
//    SlowTrigTracker::SpillInfo  Ev1stSpill3rd_trk =  TRK_.get1stEvent3rdSplill(trackerTree, EvLastSpill2nd_trk.Index);
//
//
//    cout << "SlowTrigTracker: first event's index and time in spill 1: " << Ev1stSpill1st_trk.Index << " : " << Ev1stSpill1st_trk.Time << "\n";
//    cout << "SlowTrigTracker: last event's index and time in spill 1:  " << EvLastSpill1st_trk.Index << " : " << EvLastSpill1st_trk.Time << "\n";
//    cout << "SlowTrigTracker: first event's index and time in spill 2: " << Ev1stSpill2nd_trk.Index << " : " << Ev1stSpill2nd_trk.Time << "\n";
//    cout << "SlowTrigTracker: last event's index and time in spill 2:  " << EvLastSpill2nd_trk.Index << " : " << EvLastSpill2nd_trk.Time << "\n";
//    cout << "SlowTrigTracker: first event's index and time in spill 3: " << Ev1stSpill3rd_trk.Index << " : " << Ev1stSpill3rd_trk.Time << "\n";



    //sanity check
//    Int_t totalNumEve=TOFHIR_TriggerIndices.size();
//    if (trackerTree->GetEntries()/TOFHIR_TriggerIndices.size() < 1.0 or trackerTree->GetEntries()/TOFHIR_TriggerIndices.size() > 1.5 ){
//        cout <<"Tracker has much less or much more events than TOFHIR !!!\n";
//        exit(-1);
//    }

    //=====================RECREATE outputFile=========================
//    const char *outFileName   = argv[3];
    std::string OutName="outFile_"+runNumber+".root";
    TFile *outFile = new TFile(OutName.c_str(),"recreate");
    TTree *outTree = new TTree("data","data");
    outTree->SetAutoSave();
//
//    //=====================Config Variables=========================
    const UInt_t triggerChannelID = 384;
//    const double  ClockScaleFactor=1.04617;

//    //=====================Initialize input tree variables=========================
    UInt_t channelID;
    Long64_t time;
    float tot;
    float xIntercept;
    float yIntercept;
    float xSlope;
    float ySlope;
    float x_dut;
    float y_dut;
    float chi2;
    int ntracks;
    int nplanes;
    float step1_;
    float step2_;
    int matchEff;

    channelID=-9999;
    time=-9999;
    tot=-9999;
    xIntercept=-9999;
    yIntercept=-9999;
    xSlope=-9999;
    ySlope=-9999;
    x_dut=-9999;
    y_dut=-9999;
    chi2 = -9999;
    ntracks = -1;
    nplanes = -1;
    step1_=-99;
    step2_=-99;
    matchEff=0;

    int totalNumEve=TOF_.triggeredTofhirEv.size();
    int totalNumEveMatched=1;

    int numSpill=0;
//    //==============set Branch addresses for all the input variables================
    tofhirTree->SetBranchAddress("channelID",&channelID);
    tofhirTree->SetBranchAddress("time",&time);
    tofhirTree->SetBranchAddress("tot",&tot);
    tofhirTree->SetBranchAddress("step1",&step1_);
    tofhirTree->SetBranchAddress("step2",&step2_);
//    //=====================Initialize output tree variables=========================
    Long64_t chTime[400];
    float chtot[400];
    Int_t event;
    float step1;
    float step2;
    //initialize for event 1
    event=1;
    for(int k=0;k<400;k++){
        chTime[k]=-9999;
        chtot[k]=-9999;
    }
    Long64_t TofhirTimeZero=0;
    Long64_t TrackerTimeZero=0;

    //==============set Branch addresses for all the output variables================
    outTree->Branch("event",&event,"event/I");
    outTree->Branch("step1", &step1, "step1/F");
    outTree->Branch("step2", &step2, "step2/F");
    outTree->Branch("chTime",&chTime,"chTime[400]/L");
    outTree->Branch("chtot",&chtot,"chtot[400]/F");
    outTree->Branch("xIntercept", &xIntercept, "xIntercept/F");
    outTree->Branch("yIntercept", &yIntercept, "yIntercept/F");
    outTree->Branch("xSlope", &xSlope, "xSlope/F");
    outTree->Branch("ySlope", &ySlope, "ySlope/F");
    outTree->Branch("x_dut", &x_dut, "x_dut/F");
    outTree->Branch("y_dut", &y_dut, "y_dut/F");
    outTree->Branch("chi2", &chi2, "chi2/F");
    outTree->Branch("ntracks", &ntracks, "ntracks/I");
    outTree->Branch("nplanes", &nplanes, "nplanes/I");
    outTree->Branch("matchEff", &matchEff, "matchEff/I");

//    //================================================================================================================
//    //================================================================================================================
//    //======================== Look for coincidences =================================================================
//    //================================================================================================================
//    //================================================================================================================
//
    int PreviousMatchIndex = Ev1stSpill1st_trk.Index;
    cout << std::setprecision(10);

    for (Int_t q=0;q<TOF_.triggeredTofhirEv.size(); q++) {
    

        tofhirTree->GetEntry(TOF_.triggeredTofhirEv[q].Index);
        
        //================================================================================================================
        //======================== reset the clock for 2nd and 3rd spill =================================================
        //================================================================================================================
        TofhirTimeZero=Ev1stSpill1st.Time;
        TrackerTimeZero=Ev1stSpill1st_trk.Time;
            if ( q >= Ev1stSpill2nd.Index){
                TofhirTimeZero=Ev1stSpill2nd.Time;
                TrackerTimeZero=Ev1stSpill2nd_trk.Time;
            }
            if ( q >= Ev1stSpill3rd.Index){
                TofhirTimeZero=Ev1stSpill3rd.Time;
                TrackerTimeZero=Ev1stSpill3rd_trk.Time;
            }
        //================================================================================================================
        if (1) cout<<"q= "<<q<<"   time is= "<<TOF_.triggeredTofhirEv[q].Time*1e-12 <<"  dif=" <<(TOF_.triggeredTofhirEv[q].Time-TofhirTimeZero)*1e-12<<"\n";
//
//
        //initialize at beginning of each event
        for(int pp=0;pp<400;pp++){
            chTime[pp]=-9999;
            chtot[pp]=-9999;
        }

        //populate trigger data
        chTime[triggerChannelID] = time;
        chtot[triggerChannelID] = tot;
//
//
        //find signal hits that correspond to the given trigger
        // We find that sipms are typically around 170ns EARLIER than the trigger timestamp.
        // So look within a 10ns window around that.
        Long64_t tTrigger = time;  // This is the time for channelID 384

        // Just look at the hits that are 100 before the trigger or 100 after the trigger. If their time matches with trigger time then write them donw in the chTime/chtot vectors, otherwsie they will be filled with default -999; Now it is insured that the time and tot of all hits from the the single events are the same
        int j_start = TOF_.triggeredTofhirEv[q].Index-100;
        if (j_start < 0) j_start = 0;
        int j_end = TOF_.triggeredTofhirEv[q].Index+100;

        for (Int_t j=j_start;j<j_end; j++) {
            tofhirTree->GetEntry(j);

            Long64_t tdiff = tTrigger - time;

            //If channel falls within our search window, then populate the data for that channel
            if (tdiff >= 170000 - 20000 &&  tdiff <= 170000 + 5000) {
                chTime[channelID] = time;
                chtot[channelID] = tot;
            }
        }

        double TOFHIRTriggerTimestamp = (tTrigger - TofhirTimeZero)*1e-12;
//
        int NMatchedTracks = 0;
        matchEff=0;


//        for(int k=PreviousMatchIndex; k<trackerTree->GetEntries(); k++) {
//
//            trackerTree->GetEntry( k );
//            Long64_t TrackerBCOThisEvent = TRK_.trackerEvent->bco;
//
//            if (TRK_.trackerEvent->bco== 4294967295) continue;
//
//            double trackerTime = (TrackerBCOThisEvent - TrackerTimeZero) * 144e-9 * ClockScaleFactor;

//            if (1) cout<<"k = "<<k<<"  time is= "<<TrackerBCOThisEvent<< " - "<<TrackerTimeZero<<" =>   dif= "<<trackerTime<<"\n";

////            //================================================================================================================
////            //now find a match with TOFHIR trigger time
//            if (fabs (TOFHIRTriggerTimestamp - trackerTime ) < MATCH_THRESHOLD  ){
//                cout <<"Event with index of Tofhir= "<<q <<"  MATCHED w/ event w/ index of Tracker= " << k << " : " << TOFHIRTriggerTimestamp << "   | "<< trackerTime<<"\n";
//                totalNumEveMatched++;
//                NMatchedTracks++;
//                PreviousMatchIndex = k;
//
//                //populate tracking data
//                xIntercept=TRK_.trackerEvent->xIntercept * 1e-3;
//                yIntercept=TRK_.trackerEvent->yIntercept * 1e-3;
//                xSlope=TRK_.trackerEvent->xSlope * 1e-3;
//                ySlope=TRK_.trackerEvent->ySlope * 1e-3;
//                x_dut= (TRK_.trackerEvent->xIntercept + TRK_.trackerEvent->xSlope * 2.0e5) * 1e-3;
//                y_dut= (TRK_.trackerEvent->yIntercept + TRK_.trackerEvent->ySlope * 2.0e5) * 1e-3;
//                chi2 = TRK_.trackerEvent->chi2;
//                ntracks = NMatchedTracks;
//                nplanes = TRK_.trackerEvent->nPlanes;
//                step1=step1_;
//                step2=step2_;
//                matchEff=1;
//
//            }
//            else if (TOFHIRTriggerTimestamp - trackerTime < -0.001 || (TOFHIRTriggerTimestamp==0 && trackerTime==0)) {
//            
//                break;
//            }
//        }

        outTree->Fill();
        event++;
    }

    outFile->Write();
    outFile->Close();
    tofhirFile->Close();
//    cout<<"\n============>  Matching efficiency is "<<totalNumEveMatched <<"/" <<totalNumEve <<" = "<< 1.0*totalNumEveMatched/totalNumEve<<"\n";
//    cout<<"====================================================================================================================================\n\n";
}
