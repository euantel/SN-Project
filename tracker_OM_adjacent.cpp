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

bool within_x(vector<int> a, vector<int> b, int x) {
    //check if 3D vectors are within a x by x box on the same side. takes {side, col, layer}. 
    if (a[0] != b[0]) {return 0;}
    if (a[1] < (b[1]-x) || a[1] > (b[1]+x)) {return 0;}
    if (a[2] < (b[2]-x) || a[1] > (b[2]+x)) {return 0;}
    return 1;
}

void tracker_OM_adjacent() {

    //get tree and setup relevant branches
    TFile *f = new TFile("snemo_run-1166_udd.root", "READ");
    TTree *tree = (TTree*)f->Get("SimData");

    gInterpreter->GenerateDictionary("vector<vector<int>>","vector");          //trying this
    gInterpreter->GenerateDictionary("vector<vector<long>>","vector");

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

    /* alternate way of cutting, leaving out for now 

    //get cut of events with a calorimeter hit, maybe add tracker cells > 3
    TCut cut_calohit = "digicalo.nohits > 3";    //pairing

    tree->Draw(">>elist", cut_calohit, "entrylist");
    TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
    int totalentries = elist->GetN();
    cout << totalentries << " entries with 4+ calorimeter hits\n";

    */

    //create outfile and tree for particle tracks
    TFile *out = new TFile("tracks.root", "RECREATE");
    TTree *outtree = new TTree("tracks", "SimData Tracks");
    int hitentry = 0;
    int hit_eventid = 0;
    int hit_caloid, hit_caloside, hit_calorow, hit_calocol = 0;
    double hit_energy = 0;
    long hit_time = 0;

    vector<vector<int>> *hit_track = new vector<vector<int>>;
    outtree->Branch("entry", &hitentry, "hitentry/I");
    outtree->Branch("hit_event_id", &hit_eventid, "hit_eventid/I");
    outtree->Branch("hit_timestamp", &hit_time);
    outtree->Branch("hit_energy", &hit_energy, "hit_energy/D");
    outtree->Branch("hit_calo_id", &hit_caloid, "hit_caloid/I");
    outtree->Branch("hit_calo_side", &hit_caloside, "hit_caloside/I");
    outtree->Branch("hit_calo_row", &hit_calorow, "hit_calorow/I");
    outtree->Branch("hit_calo_column", &hit_calocol, "hit_calocol/I");
    outtree->Branch("hit_track", &hit_track);

    //out TTree for final gamma entries, may integrate into main code at some point
    TTree *gammas = new TTree("tracks_gammas", "Electron-like tracks with correlated gamma");
    int g_event = 0;
    int e_caloid, gamma_caloid = -1;
    double e_energy, gamma_energy = 0;
    vector<int> e_calo, gamma_calo = {};
    vector<vector<int>> *e_track = new vector<vector<int>>; 
    long gamma_timestamp, e_timestamp;
    gammas->Branch("event_id", &g_event, "g_event/I");
    gammas->Branch("gamma_energy", &gamma_energy, "gamma_energy/D");
    gammas->Branch("gamma_caloid", &gamma_caloid, "gamma_caloid/I");
    gammas->Branch("gamma_calo", &gamma_calo);
    gammas->Branch("gamma_timestamp", &gamma_timestamp, "gamma_timestamp/L");
    gammas->Branch("e_energy", &e_energy, "e_energy/D");
    gammas->Branch("e_caloid", &e_caloid, "e_caloid/I");
    gammas->Branch("e_calo", &e_calo);
    gammas->Branch("e_timestamp", &e_timestamp, "e_timestamp/L");
    gammas->Branch("e_track", &e_track);

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

    int good_events = 0;
    int totalentries = tree->GetEntries();

    //FLAGS                                
    bool flag_cut_calohits = 1;             //enforce 2 calo hits (4 due to pairing)
    bool flag_cut_tracklength = 1;          //enforce electron track length 4+ cells
    bool flag_cut_e_energy = 1;             //enforce electron energy > 0.3 MeV
    bool flag_cut_OM_delta_t = 1;           //enforce OM-tracker time difference of -0.2 to 50 us
    bool flag_track_recon = 1;              //reconstruct tracks and save 
    bool flag_gammas = 1;                   //run extra loop to seek correlated gammas (within 50ns)

    for (int i=0; i < totalentries; i++) {
        tree->GetEntry(i);

        if (calohits < 4 && flag_cut_calohits == 1) {continue;}
        if (trackercolumn->size() < 4 && flag_cut_tracklength == 1) {continue;}

        vector<vector<int>> *gamma_candidate = new vector<vector<int>>;
        vector<double> *gamma_candidate_e = new vector<double>;
        vector<long> *gamma_candidate_t = new vector<long>;

        for (int j=0; j<calohits; j++) {            //for each hit calorimeter j
            if (wall->at(j) != -1) {continue;}      //main wall only?

            int col = calocolumn->at(j); 
            float tcol_min = tab_column.at(col) - 4.; 
            float tcol_max = tab_column.at(col) + 4.;

            vector<vector<int>> *track = new vector<vector<int>>;

            hit_calocol = col;
            hit_caloside = caloside->at(j);
            hit_calorow = calorow->at(j);
            hit_caloid = (13*col) + (260*caloside->at(j)) + calorow->at(j);          
            hit_energy = (charge->at(j))*calib[hit_caloid]*energy_conv;             //changed to MeV and flipped sign
            if (hit_energy < 0.3 && flag_cut_e_energy == 1) {continue;}          

            if (hit_caloid == 123) {spectrum->Fill(hit_energy);}                //record energies at specific OM 

            bool adj_tracker = 0;
            for (int k=0; k<trackercolumn->size(); k++) {
                if (trackerside->at(k) != caloside->at(j)) {continue;}          //check same side
                if (trackerlayer->at(k) < 7) {continue;}                        //check layer 8 or 9

                long delta_t = (2*anode_R0->at(k).at(0) - timestamp->at(j))*6.25/1000.;         //microseconds
                timehist->Fill(delta_t);

                if ((delta_t < 0.2 || delta_t > 50) && flag_cut_OM_delta_t == 1) {continue;}    //time cut 

                int tcol = trackercolumn->at(k);
                int tlayer = trackerlayer->at(k);
                int tside = trackerside->at(k);

                if (tcol <= tcol_max && tcol >= tcol_min) {              
                    track->push_back({tside, tcol, tlayer});
                    adj_tracker = 1; 
                    break;
                }
            }

            if (adj_tracker == 0) {
                gamma_candidate->push_back({hit_caloid, hit_caloside, hit_calorow, hit_calocol});
                gamma_candidate_e->push_back(hit_energy);
                gamma_candidate_t->push_back(timestamp->at(j));
                }

            //loop again to reconstruct then save 
            if (track->size() == 0) {continue;}
            while (true == true) {
                int pre_size = track->size();
                for (int l=0; l<track->size();l++) {
                    if (flag_track_recon == 0) {break;}
                    for (int m=0;m<trackercolumn->size();m++) {
                        vector<int> next_tracker = {trackerside->at(m), trackercolumn->at(m), trackerlayer->at(m)};
                        
                        bool dupe = 0;                                      //check not already in track vector
                        for (int n=0;n<track->size();n++) {
                            if (next_tracker == track->at(n)) {dupe = 1;}
                        }
                        if (dupe == 1) {continue;}

                        if (within_x(track->at(l), next_tracker, 2)) {       //append to track vector if within 2x2 box
                            track->push_back(next_tracker);               
                        }
                    }
                }

                if (pre_size == track->size()) {
                    if (track->size() > 3 || flag_track_recon == 0) {       //save as "good" event if track length > 3 

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

        //per event code, try putting gamma stuff here?            TODO: move into it's own function or even new macro. 
        //get entries at this event no in outtree
        for (int j = 0; j<outtree->GetEntries(); j++) {
            if (flag_gammas == 0) {break;}

            outtree->GetEntry(outtree->GetEntries() - j);           //not even sure you can do this 
            if (hit_eventid != i) {break;}

            for (int k=0; k<gamma_candidate->size(); k++) {

                if (within_x({hit_caloside, hit_calorow, hit_calocol}, 
                    {gamma_candidate->at(k).at(1), gamma_candidate->at(k).at(2), gamma_candidate->at(k).at(3)}, 1)) {continue;}
                
                if (abs(6.25*(hit_time - gamma_candidate_t->at(k))) < 50
                && hit_calocol != gamma_candidate->at(k).at(3)) {                //exclude same column for now

                    //record event
                    g_event = hit_eventid;
                    e_caloid = hit_caloid;
                    gamma_caloid = gamma_candidate->at(k).at(0);
                    e_energy = hit_energy; 
                    gamma_energy = gamma_candidate_e->at(k);
                    e_calo = {hit_caloside, hit_calorow, hit_calocol};
                    gamma_calo = {gamma_candidate->at(k).at(1), gamma_candidate->at(k).at(2), gamma_candidate->at(k).at(3)};
                    e_track = hit_track;
                    e_timestamp = hit_time;
                    gamma_timestamp = gamma_candidate_t->at(k);
                    gammas->Fill();
                }

            }
            
        }
    }

    out->cd();
    outtree->Write();
    gammas->Write();
    cout << "Adjacent events: " << good_events << "\n";

    spectrum->SetTitle("Energy spectrum for specific OM;Energy (MeV);Count");
    timehist->SetTitle("Time difference between OM and adj tracker;delta_t (us);Count");
    spectrum->Write();
    timehist->Write();

    //spectrum->Draw();
    timehist->Draw();

}