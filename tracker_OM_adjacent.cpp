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


*/

bool within_two(vector<int> a, vector<int> b) {
    //check if 3D vectors are within a 2x2 box on the same side. takes {side, col, layer}. 
    if (a[0] != b[0]) {return 0;}
    if (a[1] < (b[1]-2) || a[1] > (b[1]+2)) {return 0;}
    if (a[2] < (b[2]-2) || a[1] > (b[2]+2)) {return 0;}
    return 1;
}

void tracker_OM_adjacent() {

    //get tree and setup relevant branches
    TFile *f = new TFile("snemo_run-1166_udd.root", "READ");
    TTree *tree = (TTree*)f->Get("SimData");

    gInterpreter->GenerateDictionary("vector<vector<int> >","vector");          //trying this
    gInterpreter->GenerateDictionary("vector<vector<long> >","vector");

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
    vector<int> *wall = new vector<int>;
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
    tree->SetBranchStatus("digicalo.wall", 1);
    tree->SetBranchAddress("digicalo.wall", &wall);
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

    std::vector<vector<long>> *anode_R0 = new std::vector<vector<long>>;
    tree->SetBranchStatus("digitracker.anodetimestampR0", 1);
    tree->SetBranchAddress("digitracker.anodetimestampR0", &anode_R0);

    
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

    vector<vector<int>> *hit_track = new vector<vector<int>>;
    outtree->Branch("entry", &hitentry, "hitentry/I");
    outtree->Branch("hit_event_id", &hit_eventid, "hit_eventid/I");
    outtree->Branch("hit_timestamp", &hit_time);
    outtree->Branch("hit_calo_id", &hit_caloid, "hit_caloid/I");
    outtree->Branch("hit_energy", &hit_energy, "hit_energy/D");
    outtree->Branch("hit_track", &hit_track);

    //comparison vector, +- 4 inclusive gives tracker columns 
    vector<double> tab_column = {0.8, 6.8, 12.5, 18.2, 24.2, 29.9, 35.9, 41.8, 46.7, 53.6, 59.4, 65.3, 71.1, 77, 83.1, 88.8, 94.8, 100.6, 106.4, 112.2};

    //input charge-energy calibration from text
    double energy_conv = -1./4194.304;
    vector<double> calib;
    ifstream calib_file("run-1351_fee-charge-to-energy.txt");
    int n1;
    double n2;
    while (calib_file >> n1 >> n2) {
        calib.push_back(n2);
    }

    //check specific calorimeter for energy spectrum
    TH1D *spectrum = new TH1D("spectrum", "Energies for given OM", 100, 0, 10);
    TH1D *timehist = new TH1D("time", "OM to tracker delta_t", 120, -20, 100);

    //check calorimeters for nearby active tracker cells
    int good_events = 0;

    for (int i=0; i < totalentries; i++) {
        int entryno = elist->GetEntry(i);
        tree->GetEntry(entryno);

        if (calohits < 4) {continue;}
        if (trackercolumn->size() < 4) {continue;}

        for (int j=0; j<calohits; j++) {            //for each hit calorimeter j
            if (wall->at(j) != -1) {continue;}      //main wall only?

            int col = calocolumn->at(j); 
            float tcol_min = tab_column.at(col) - 4.; 
            float tcol_max = tab_column.at(col) + 4.;

            vector<vector<int>> *track = new vector<vector<int>>;

            hit_caloid = (13*col) + (260*caloside->at(j)) + calorow->at(j);          
            hit_energy = (charge->at(j))*calib[hit_caloid]*energy_conv;             //changed to MeV and flipped sign
            if (hit_energy < 0.3) {continue;}          

            if (hit_caloid == 123) {spectrum->Fill(hit_energy);}                //record energies at specific OM 

            for (int k=0; k<trackercolumn->size(); k++) {
                if (trackerside->at(k) != caloside->at(j)) {continue;}          //check same side
                if (trackerlayer->at(k) < 7) {continue;}                        //check layer 8 or 9

                long delta_t = (2*anode_R0->at(k).at(0) - timestamp->at(j))*6.25/1000.;         //microseconds
                timehist->Fill(delta_t);

                if (delta_t < 0.2 || delta_t > 50) {continue;}                                  //time cut 

                int tcol = trackercolumn->at(k);
                int tlayer = trackerlayer->at(k);
                int tside = trackerside->at(k);

                if (tcol <= tcol_max && tcol >= tcol_min) {              
                    track->push_back({tside, tcol, tlayer});
                    break;
                }
            }
            //loop again to reconstruct then save 
            if (track->size() == 0) {continue;}
            while (true == true) {
                int pre_size = track->size();
                for (int l=0; l<track->size();l++) {
                    for (int m=0;m<trackercolumn->size();m++) {
                        vector<int> next_tracker = {trackerside->at(m), trackercolumn->at(m), trackerlayer->at(m)};
                        
                        bool dupe = 0;                                      //check not already in track vector
                        for (int n=0;n<track->size();n++) {
                            if (next_tracker == track->at(n)) {dupe = 1;}
                        }
                        if (dupe == 1) {continue;}

                        if (within_two(track->at(l), next_tracker)) {       //append to track vector if within 2x2 box
                            track->push_back(next_tracker);               
                        }
                    }
                }

                if (pre_size == track->size()) {
                    if (track->size() > 3) {                        //save as "good" event if track length > 3 

                        hitentry = good_events;
                        good_events += 1;

                        hit_eventid = event;
                        hit_time = timestamp->at(j);

                        hit_track = track;                          //add this to outtree
                        outtree->Fill();
                    }

                    break;
                }
            }
        }
    }

    out->cd();
    outtree->Write();
    cout << "Adjacent events: " << good_events << "\n";

    spectrum->SetTitle("Energy spectrum for specific OM;Energy (MeV);Count");
    timehist->SetTitle("Time difference between OM and adj tracker;delta_t (us);Count");
    spectrum->Write();
    timehist->Write();

    //spectrum->Draw();
    timehist->Draw();

}