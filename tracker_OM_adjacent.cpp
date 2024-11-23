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
    int e_hit_caloid = 0, gamma_caloid;
    double e_hit_energy, gamma_energy = 0;
    long e_hit_time = 0, gamma_timestamp = 0;
    vector<int> *hit_track = new vector<int>;                         //vector of IDs

    //FLAGS                                
    bool flag_cut_calohits = 0;             //enforce 2 calo hits (4 due to pairing)
    bool flag_cut_tracklength = 0;          //enforce electron track length 4+ cells
    bool flag_cut_e_energy = 0;             //enforce electron energy > 0.3 MeV
    bool flag_cut_OM_delta_t = 0;           //enforce OM-tracker time difference of -0.2 to 50 us

    //re-saving data from original events 
    outtree->Branch("header.eventnumber", &event);
    outtree->Branch("digicalo.nohits", &calohits);
    outtree->Branch("digicalo.side", &caloside);
    outtree->Branch("digicalo.column", &calocolumn);
    outtree->Branch("digicalo.row", &calorow);
    outtree->Branch("digicalo.timestamp", &timestamp);
    outtree->Branch("digicalo.charge", &charge);
    outtree->Branch("digicalo.wall", &wall);
    outtree->Branch("digitracker.nohits", &trackerhits);
    outtree->Branch("digitracker.side", &trackerside);
    outtree->Branch("digitracker.column", &trackercolumn);
    outtree->Branch("digitracker.layer", &trackerlayer);
    outtree->Branch("digitracker.anodetimestampR0", &anode_R0);

    //cut flags, 0 = failed cut 
    outtree->Branch("pass_cut_calohits", &flag_cut_calohits, "flag_cut_calohits/O");
    outtree->Branch("pass_cut_energy", &flag_cut_e_energy, "flag_cut_e_energy/O");
    outtree->Branch("pass_cut_tracklength", &flag_cut_tracklength, "flag_cut_tracklength/O");
    outtree->Branch("pass_cut_OM_delta_t", &flag_cut_OM_delta_t, "flag_cut_OM_delta_t/O");

    //additional data 
    outtree->Branch("e_hit_timestamp", &e_hit_time);
    outtree->Branch("e_hit_energy", &e_hit_energy, "hit_energy/D");
    outtree->Branch("e_hit_calo_id", &e_hit_caloid, "hit_caloid/I");
    outtree->Branch("e_hit_track", &hit_track);
    outtree->Branch("gamma_timestamp", &gamma_timestamp, "gamma_timestamp/L");
    outtree->Branch("gamma_energy", &gamma_energy, "gamma_energy/D");
    outtree->Branch("gamma_caloid", &gamma_caloid, "gamma_caloid/I");

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
    int cut_calohits = 0, cut_e_energy = 0, cut_OM_deltat = 0;

    for (int i=0; i < totalentries/10; i++) {
        tree->GetEntry(i);

        //set default values
        e_hit_time = -1;
        e_hit_energy = -1;
        e_hit_caloid = -1;
        hit_track->clear();
        gamma_timestamp = -1;
        gamma_energy = -1;
        gamma_caloid = -1;
        flag_cut_calohits = 0;
        flag_cut_e_energy = 0;
        flag_cut_OM_delta_t = 0;
        flag_cut_tracklength = 0;

        if (calohits == 4) {
            flag_cut_calohits = 1;   //passed first cut
        } else {
            outtree->Fill();
            continue;
        }

        for (int j=0; j<calohits; j++) {            //for each hit calorimeter j
            if (wall->at(j) != -1) {continue;}      //main wall only?

            int col = calocolumn->at(j); 
            float tcol_min = tab_column.at(col) - 4.; 
            float tcol_max = tab_column.at(col) + 4.;

            vector<vector<int>> *track = new vector<vector<int>>;
            track->clear();

            bool adj_tracker = 0;

            int hit_caloid = (13*col) + (260*caloside->at(j)) + calorow->at(j);          
            double hit_energy = (charge->at(j))*calib[hit_caloid]*energy_conv;             //changed to MeV and flipped sign

            if (hit_energy > 0.3) {         //energy cut
                flag_cut_e_energy = 1;
            } else {
                continue;
            }

            if (hit_caloid == 123) {spectrum->Fill(hit_energy);}                //record energies at specific OM 

            for (int k=0; k<trackercolumn->size(); k++) {
                if (trackerside->at(k) != caloside->at(j)) {continue;}          //check same side
                if (trackerlayer->at(k) < 7) {continue;}                        //check layer 8 or 9

                int tcol = trackercolumn->at(k);
                int tlayer = trackerlayer->at(k);
                int tside = trackerside->at(k);

                if (tcol <= tcol_max && tcol >= tcol_min) {                    //within +- 4 range, start track 

                    long delta_t = (2.*anode_R0->at(k).at(0) - timestamp->at(j))*6.25/1000.;         //microseconds
                    timehist->Fill(delta_t);

                    if ((delta_t > -0.2 && delta_t < 50)) {     //time cut 
                        flag_cut_OM_delta_t = 1;
                    }

                    track->push_back({tside, tcol, tlayer});
                    adj_tracker = 1; 
                    break;
                }
            }

            if (adj_tracker == 0) {
                //record possible gamma hits (for now just any hit with no time correlated adjacent track
                gamma_energy = hit_energy;
                gamma_caloid = hit_caloid;
                gamma_timestamp = timestamp->at(j);
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

                        if (within_x(track->at(l), next_tracker, 2)) {       //append to track vector if within 2x2 box
                            track->push_back(next_tracker);               
                        }
                    }
                }

                if (pre_size == track->size()) {
                    e_hit_time = timestamp->at(j);
                    e_hit_energy = hit_energy;
                    e_hit_caloid = hit_caloid;

                    for (int p = 0; p<track->size(); p++) {
                        hit_track->push_back(track->at(p).at(0)*1017 + track->at(p).at(1)*9 + track->at(p).at(2));
                    }

                    break; 
                }
            }
        }

        if (hit_track->size() > 3) {
            flag_cut_tracklength = 1;
            good_events += 1;
            }
        if (flag_cut_calohits == 1) {cut_calohits += 1;}
        if (flag_cut_e_energy == 1) {cut_e_energy += 1;}
        if (flag_cut_OM_delta_t == 1) {cut_OM_deltat += 1;}
        outtree->Fill();
    }

    out->cd();
    outtree->Write();

    //output number of events cut, some events may have multiple recorded tracks 
    cout << "Initial events: \t\t" << totalentries/10 << "\n";
    cout << "Events with 4 OM hits: \t\t" << cut_calohits << "\n";
    cout << "Events with > 0.3MeV hits: \t" << cut_e_energy << "\n";
    cout << "Events with -0.2 < dt < 50us: \t" << cut_OM_deltat << "\n";
    cout << "Events with track length > 3: \t" << good_events << "\n";

    //quick text output to file 
    ofstream outtxt;
    outtxt.open("cuts.txt");
    outtxt << "Initial events: \t\t" << totalentries << "\n";
    outtxt << "Events with 4+ OM hits: \t" << cut_calohits << "\n";
    outtxt << "Events with > 0.3MeV hits: \t" << cut_e_energy << "\n";
    outtxt << "Events with -0.2 < dt < 50us: \t" << cut_OM_deltat << "\n";
    outtxt << "Events with track length > 3: \t" << good_events << "\n";
    outtxt.close();

    spectrum->SetTitle("Energy spectrum for specific OM;Energy (MeV);Count");
    timehist->SetTitle("Time difference between OM and adj tracker;delta_t (us);Count");
    spectrum->Write();
    timehist->Write();

    //spectrum->Draw();
    //timehist->Draw();

}

void find_gammas() {
    /*
    
    trying gamma stuff here, to be ran after the main macro
    look for events with an electron-like track and OM hit within 50ns
    match-up then save to a new file
    
    */

    gInterpreter->GenerateDictionary("vector<vector<int>>","vector"); 

    //get trees from tracks file 
    TFile *g = new TFile("tracks.root", "READ");
    TTree *gammas = (TTree*)g->Get("gammas");
    int g_eventid = 0, g_caloside = 0, g_calorow = 0, g_calocol = 0;
    Long64_t g_time = 0;
    double g_energy = 0;
    gammas->SetBranchAddress("event_id", &g_eventid);
    gammas->SetBranchAddress("gamma_timestamp", &g_time);
    gammas->SetBranchAddress("gamma_energy", &g_energy);
    gammas->SetBranchAddress("gamma_caloside", &g_caloside);
    gammas->SetBranchAddress("gamma_calorow", &g_calorow);
    gammas->SetBranchAddress("gamma_calocol", &g_calocol);

    TTree *electrons = (TTree*)g->Get("tracks");
    int e_eventid = 0, e_caloside = 0, e_calorow = 0, e_calocol = 0;
    long e_time = 0;
    double e_energy = 0;
    vector<vector<int>> *e_track = new vector<vector<int>>;
    electrons->SetBranchAddress("hit_event_id", &e_eventid);
    electrons->SetBranchAddress("hit_timestamp", &e_time);
    electrons->SetBranchAddress("hit_energy", &e_energy);
    electrons->SetBranchAddress("hit_calo_side", &e_caloside);
    electrons->SetBranchAddress("hit_calo_row", &e_calorow);
    electrons->SetBranchAddress("hit_calo_column", &e_calocol);
    electrons->SetBranchAddress("hit_track", &e_track);

    //output seperate file for now
    TFile *out = new TFile("e_gamma_events.root", "RECREATE");
    TTree *out_events = new TTree("e_gammas", "Correlated electron + gammas");

    int out_eventid = 0 ;
    vector<int> *e_OM = new vector<int>;
    vector<int> *g_OM = new vector<int>;
    double out_e_energy = 0, out_g_energy = 0;
    long out_e_time = 0, out_g_time = 0;
    vector<vector<int>> *out_e_track = new vector<vector<int>>;

    out_events->Branch("eventID", &out_eventid, "out_eventid/I");
    out_events->Branch("electron_OM", &e_OM);
    out_events->Branch("gamma_OM", &g_OM);
    out_events->Branch("electron_energy", &out_e_energy, "out_e_energy/D");
    out_events->Branch("electron_timestamp", &out_e_time);
    out_events->Branch("gamma_energy", &out_g_energy, "out_g_energy/D");
    out_events->Branch("gamma_timestamp", &out_g_time);
    out_events->Branch("electron_track", &out_e_track);

    int events_found = 0;

    TH1D *gamma_spectrum = new TH1D("gamma_energies", "Gamma Energies", 140, 0, 7);

    //loop through both and match up timestamps 
    for (int i=0; i<electrons->GetEntries(); i++) {
        electrons->GetEntry(i);

        if (i % 250 == 0) {cout << i << "\n";}
        if (i == 25000) {break;}                                         //break for testing
        
        //using entrylist seems faster than looping for comparison
        TString cutstring = TString::Format("event_id == %d", e_eventid);

        gammas->Draw(">>elist", cutstring, "entrylist");
        TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
        int matches = elist->GetN();
        

        for (int j=0; j<matches; j++) {
            gammas->GetEntry(elist->GetEntry(j));
            
            e_OM->clear();
            g_OM->clear();
            
            if (e_calocol == g_calocol) {continue;}                     //avoid same column for now 
            if (abs(6.25*(e_time - g_time)) > 50) {continue;}           //enforce time correlation 50ns

            //discard adjacent OM hits
            if (within_x({e_caloside, e_calorow, e_calocol}, {g_caloside, g_calorow, g_calocol}, 1)) {continue;}

            //record events that pass the cuts
            out_eventid = e_eventid;
            e_OM->push_back(e_caloside);
            e_OM->push_back(e_calorow);
            e_OM->push_back(e_calocol);
            g_OM->push_back(g_caloside);
            g_OM->push_back(g_calorow);
            g_OM->push_back(g_calocol);
            out_e_energy = e_energy;
            out_e_time = e_time;
            out_g_energy = g_energy;
            out_g_time = g_time;
            out_e_track = e_track;
            out_events->Fill();

            events_found += 1;
            gamma_spectrum->Fill(g_energy);

            break;                     //will only take one from each for now
        }
    }

    out->cd();
    out_events->Write();
    gamma_spectrum->SetTitle("Energies of Correlated Gammas;E (MeV);Count");
    gamma_spectrum->Write();
    cout << "Electron + Gamma-like events: " << events_found << "\n";
}