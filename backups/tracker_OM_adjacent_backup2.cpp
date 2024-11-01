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

    //get cut of events with a calorimeter hit, maybe add more cuts?
    TCut cut_calohit = "digicalo.nohits > 0";

    tree->Draw(">>elist", cut_calohit, "entrylist");
    TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
    int totalentries = elist->GetN();
    cout << totalentries << " entries with 1+ calorimeter hit\n";

    //create outfile and tree for particle tracks
    TFile *out = new TFile("tracks.root", "RECREATE");
    TTree *outtree = new TTree("tracks", "SimData Tracks");
    int hit_eventid = 0;
    int hit_caloid = 0;
    vector<int> *hittrackSide = new vector<int>;
    vector<int> *hittrackColumn = new vector<int>;
    vector<int> *hittrackLayer = new vector<int>;
    outtree->Branch("hit_event_id", &hit_eventid, "hit_eventid/I");
    outtree->Branch("hit_calo_id", &hit_caloid, "hit_caloid/I");
    outtree->Branch("hit_track_side", &hittrackSide);
    outtree->Branch("hit_track_column", &hittrackColumn);
    outtree->Branch("hit_track_layer", &hittrackLayer);

    //comparison vector, +- 4 inclusive gives tracker columns 
    vector<double> tab_column = {0.8, 6.8, 12.5, 18.2, 24.2, 29.9, 35.9, 41.8, 46.7, 53.6, 59.4, 65.3, 71.1, 77, 83.1, 88.8, 94.8, 100.6, 106.4, 112.2};

    //check calorimeters for nearby active tracker cells
    int good_events = 0;

    for (int i=0; i < 1000; i++) {
        int entryno = elist->GetEntry(i);
        tree->GetEntry(entryno);

        for (int j=0; j<calohits; j++) {          //for each hit calorimeter j
            int col = calocolumn->at(j); 
            float tcol_min = tab_column.at(col) - 4; 
            float tcol_max = tab_column.at(col) + 4;

            for (int k=0; k<trackercolumn->size(); k++) {
                if (trackerside->at(k) != caloside->at(j)) {continue;}          //check same side
                if (trackerlayer->at(k) < 7) {continue;}                        //check layer 8 or 9

                int tcol = trackercolumn->at(k);
                if (tcol <= tcol_max && tcol >= tcol_min) {                     //check within OM range
                    //do something or output something
                    good_events += 1;
                    hit_eventid = event;
                    hit_caloid = (13*col) + (260*caloside->at(j)) + calorow->at(j);
                    hittrackSide->push_back(trackerside->at(k));
                    hittrackColumn->push_back(tcol);
                    hittrackLayer->push_back(trackerlayer->at(k));

                    int l = 0;
                    while (0) {
                        //TODO track reconstruction, append vectors then save

                        l += 1;
                    }

                    //write branch here 
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