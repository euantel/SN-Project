#include "TStyle.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TPaveText.h"
#include <iostream>
#include <memory>
#define _USE_MATH_DEFINES
#include <fstream>
#include <sstream>
#include <vector>
#include <TH1D.h>
#include <TH2D.h>
#include <TH3D.h>
#include <TF1.h>
#include <TStyle.h>
#include "TKey.h"
#include "TFile.h"
#include "TTree.h"
#include "TLine.h"
#include "TROOT.h"
#include <TText.h>
#include <TLatex.h>
#include <TRandom3.h>
#include <TLegend.h>
#include <TSystem.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <TParameter.h>
#include <TPaletteAxis.h>
using namespace std;

/* checks for events with >= 1 active calorimeter and and adjacent 
layer 8 or 9 tracker hit

om id = 260*side + column*13 + row
tracker id = side + column*9 + layer           (maybe wrong)
*/

bool within_two(int a, int b) {
    //check if: b-2 <= a <= b+2
    if (a >= (b-2) && a <= (b+2)) {return 1;}
    return 0;
}

void tracker_OM_adjacent() {

    //get tree and setup relevant branches
    TFile *f = new TFile("snemo_run-1166_udd.root", "READ");
    TTree *tree = (TTree*)f->Get("SimData");

    int event = 0;
    int calohits = 0;
    int trackerhits = 0;
    vector<int> *caloid = new vector<int>;
    vector<int> *caloside = new vector<int>;
    vector<int> *calocolumn = new vector<int>;
    vector<int> *calorow = new vector<int>;
    vector<long> *timestamp = new vector<long>;
    vector<int> *charge = new vector<int>;
    vector<int> *trackerid = new vector<int>;
    vector<int> *trackerside = new vector<int>;
    vector<int> *trackercolumn = new vector<int>;
    vector<int> *trackerlayer = new vector<int>;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("header.eventnumber", 1);
    tree->SetBranchAddress("header.eventnumber", &event);
    tree->SetBranchStatus("digicalo.id", 1);                    //calorimeters
    tree->SetBranchAddress("digicalo.id", &caloid);
    tree->SetBranchStatus("digicalo.nohits", 1);                   
    tree->SetBranchAddress("digicalo.nohits", &calohits);
    tree->SetBranchStatus("digicalo.side", 1);
    tree->SetBranchAddress("digicalo.side", &caloside);
    tree->SetBranchStatus("digicalo.column", 1);
    tree->SetBranchAddress("digicalo.column", &calocolumn);
    tree->SetBranchStatus("digicalo.row", 1);
    tree->SetBranchAddress("digicalo.row", &calorow);
    tree->SetBranchStatus("digicalo.timestamp", 1);
    tree->SetBranchAddress("digicalo.timestamp", &timestamp);
    tree->SetBranchStatus("digicalo.charge", 1);
    tree->SetBranchAddress("digicalo.charge", &charge);
    tree->SetBranchStatus("digitracker.id", 1);                    //trackers
    tree->SetBranchAddress("digitracker.id", &trackerid);
    tree->SetBranchStatus("digitracker.nohits", 1);                 
    tree->SetBranchAddress("digitracker.nohits", &trackerhits);
    tree->SetBranchStatus("digitracker.side", 1);
    tree->SetBranchAddress("digitracker.side", &trackerside);
    tree->SetBranchStatus("digitracker.column", 1);
    tree->SetBranchAddress("digitracker.column", &trackercolumn);
    tree->SetBranchStatus("digitracker.layer", 1);
    tree->SetBranchAddress("digitracker.layer", &trackerlayer);

    //get cut of events with a calorimeter hit, maybe add tracker cells > 3
    TCut cut_calohit = "digicalo.nohits > 3";    //pairing

    tree->Draw(">>elist", cut_calohit, "entrylist");
    TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
    int totalentries = elist->GetN();
    cout << totalentries << " entries with 4+ calorimeter hits\n";

    //create outfile and tree for particle tracks
    TFile *out = new TFile("tracks.root", "RECREATE");
    TTree *outtree = new TTree("tracks", "SimData Tracks");
    int hitentry = 0;
    int hit_eventid = 0;
    int hit_caloid = 0;
    double hit_energy = 0;
    long hit_time = 0;
    vector<int> *hittrackSide = new vector<int>;
    vector<int> *hittrackColumn = new vector<int>;
    vector<int> *hittrackLayer = new vector<int>;
    outtree->Branch("entry", &hitentry, "hitentry/I");
    outtree->Branch("hit_event_id", &hit_eventid, "hit_eventid/I");
    outtree->Branch("hit_timestamp", &hit_time);
    outtree->Branch("hit_calo_id", &hit_caloid, "hit_caloid/I");
    outtree->Branch("hit_energy", &hit_energy, "hit_energy/D");
    outtree->Branch("hit_track_side", &hittrackSide);
    outtree->Branch("hit_track_column", &hittrackColumn);
    outtree->Branch("hit_track_layer", &hittrackLayer);

    //comparison vector, +- 4 inclusive gives tracker columns 
    vector<double> tab_column = {0.8, 6.8, 12.5, 18.2, 24.2, 29.9, 35.9, 41.8, 46.7, 53.6, 59.4, 65.3, 71.1, 77, 83.1, 88.8, 94.8, 100.6, 106.4, 112.2};

    //input charge-energy calibration from text
    vector<double> calib;
    ifstream calib_file("run-1351_fee-charge-to-energy.txt");
    int n1;
    double n2;
    while (calib_file >> n1 >> n2) {
        calib.push_back(n2);
    }

    //check calorimeters for nearby active tracker cells
    int good_events = 0;

    for (int i=0; i < totalentries; i++) {
        int entryno = elist->GetEntry(i);
        tree->GetEntry(entryno);

        for (int j=0; j<calohits; j++) {          //for each hit calorimeter j
            int col = calocolumn->at(j); 
            float tcol_min = tab_column.at(col) - 4; 
            float tcol_max = tab_column.at(col) + 4;

            hit_caloid = (13*col) + (260*caloside->at(j)) + calorow->at(j);          
            hit_energy = (charge->at(j))*calib[hit_caloid]*(-1./1000.);             //changed to MeV and flipped sign
            if (hit_energy < 0.3) {continue;}          

            //cout << "Event: " << event << "\tCaloID: " << hit_caloid << "\tCharge: " << charge->at(j) << "\tE: " << hit_energy << "\n";

            for (int k=0; k<trackercolumn->size(); k++) {
                if (trackerside->at(k) != caloside->at(j)) {continue;}          //check same side
                if (trackerlayer->at(k) < 7) {continue;}                        //check layer 8 or 9

                int tcol = trackercolumn->at(k);
                int tlayer = trackerlayer->at(k);
                int tside = trackerside->at(k);
                if (tcol <= tcol_max && tcol >= tcol_min) {      
                    
                    hitentry = good_events;
                    good_events += 1;

                    hit_eventid = event;
                    hit_time = timestamp->at(j);

                    hittrackSide->clear();
                    hittrackColumn->clear();
                    hittrackLayer->clear();
                    hittrackSide->push_back(tside);
                    hittrackColumn->push_back(tcol);
                    hittrackLayer->push_back(tlayer);

                    /*         FIND BETTER WAY TO DO THIS, maybe cuts? --------------------------------------

                    while (1) {
                        int start = 0;
                        int tracksize = hittrackColumn->size();
                        for (int l=start;l<tracksize; l++) {           //loop through current track

                            for (int m=0; m<trackercolumn->size();m++) {

                                bool duplicate = 0;                 //check not already in track vector
                                for (int n=0;n<tracksize;n++){
                                    if (trackerside->at(m) == hittrackSide->at(n) &&
                                        trackercolumn->at(m) == hittrackColumn->at(n) &&
                                        trackerlayer->at(m) == hittrackLayer->at(n)) {duplicate = 1;}
                                }

                                if (duplicate) {continue;}
                                if (trackerside->at(m) != hittrackSide->at(l)) {continue;}
                                if (within_two(trackercolumn->at(m), hittrackColumn->at(l)) == 0) {continue;}
                                if (within_two(trackerlayer->at(m), hittrackLayer->at(l)) == 0) {continue;}

                                //append track vector
                                hittrackSide->push_back(trackerside->at(m));
                                hittrackColumn->push_back(trackercolumn->at(m));
                                hittrackLayer->push_back(trackerlayer->at(m));
                                
                            }

                            start += 1;

                        }
                        
                        if (tracksize == hittrackColumn->size()) {break;}       //break if nothing changes
                    }      
                    */

                    outtree->Fill();

                    break;                                                 
                }
            }
        }
    }

    out->cd();
    outtree->Write();
    cout << "Adjacent events: " << good_events << "\n";

}