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

/* 
for application to a simulated dataset
*/

bool within_x(vector<int> a, vector<int> b, int x) {
    //check if 3D vectors are within a x by x box on the same side. takes {side, col, layer}. 
    if (a[0] != b[0]) {return 0;}
    if (a[1] < (b[1]-x) || a[1] > (b[1]+x)) {return 0;}
    if (a[2] < (b[2]-x) || a[2] > (b[2]+x)) {return 0;}
    return 1;
}

void e_gamma_sim() {

    //get tree and setup relevant branches
    TFile *f = new TFile("Falaise_UDD_5_1_2.root", "READ");
    TTree *tree = (TTree*)f->Get("SimData");

    gInterpreter->GenerateDictionary("vector<vector<int>>","vector");          //seems to fix 2D vectors
    gInterpreter->GenerateDictionary("vector<vector<long>>","vector");

    int event = 0;
    int calohits = 0;
    int trackerhits = 0;
    vector<int> *caloid = new vector<int>;
    vector<int> *caloside = new vector<int>;
    vector<int> *calocolumn = new vector<int>;
    vector<int> *calorow = new vector<int>;
    vector<double> *timestamp = new vector<double>;
    vector<double> *caloenergy = new vector<double>;
    vector<int> *calocharge = new vector<int>; //temporary check
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
    tree->SetBranchStatus("calo.time", 1);                          //changed this to calo.time 
    tree->SetBranchAddress("calo.time", &timestamp);
    tree->SetBranchStatus("calo.energy", 1);                            //Simuation uses calo.energy rather than digicalo.charge
    tree->SetBranchAddress("calo.energy", &caloenergy);
    tree->SetBranchStatus("digicalo.charge", 1);                          
    tree->SetBranchAddress("digicalo.charge", &calocharge);
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
    std::vector<vector<long>> *R5 = new std::vector<vector<long>>;
    std::vector<vector<long>> *R6 = new std::vector<vector<long>>;
    tree->SetBranchStatus("digitracker.anodetimestampR0", 1);
    tree->SetBranchAddress("digitracker.anodetimestampR0", &anode_R0);
    tree->SetBranchStatus("digitracker.bottomcathodetimestamp", 1);
    tree->SetBranchAddress("digitracker.bottomcathodetimestamp", &R5);
    tree->SetBranchStatus("digitracker.topcathodetimestamp", 1);
    tree->SetBranchAddress("digitracker.topcathodetimestamp", &R6);

    std::vector<double> *trackertime = new std::vector<double>;
    tree->SetBranchStatus("truetracker.time", 1);
    tree->SetBranchAddress("truetracker.time", &trackertime);

    /* alternate way of cutting, leaving out for now 

    //get cut of events with a calorimeter hit, maybe add tracker cells > 3
    TCut cut_calohit = "digicalo.nohits > 3";    //pairing

    tree->Draw(">>elist", cut_calohit, "entrylist");
    TEntryList *elist = (TEntryList*)gDirectory->Get("elist");
    int totalentries = elist->GetN();
    cout << totalentries << " entries with 4+ calorimeter hits\n";

    */

    //create outfile and tree for particle tracks
    TFile *out = new TFile("tracks_SIM.root", "RECREATE");
    TTree *outtree = new TTree("tracks", "SimData Tracks");
    int e_hit_caloid = 0, gamma_caloid;
    double e_hit_energy, gamma_energy = 0;
    double e_hit_time = 0, gamma_timestamp = 0;
    vector<int> *hit_track = new vector<int>;                         //vector of IDs

    //FLAGS                                
    bool flag_cut_calohits = 0;             //enforce 2 calo hits (4 due to pairing)
    bool flag_cut_tracklength = 0;          //enforce electron track length 4+ cells
    bool flag_cut_e_energy = 0;             //enforce electron energy > 0.3 MeV
    bool flag_cut_OM_delta_t = 0;           //enforce OM-tracker time difference of -0.2 to 50 us
    bool flag_e_g_correlated = 0;           //time correlation between electron and gamma-like hits

    //re-saving data from original events 
    outtree->Branch("header.eventnumber", &event);
    outtree->Branch("digicalo.nohits", &calohits);
    outtree->Branch("digicalo.side", &caloside);
    outtree->Branch("digicalo.column", &calocolumn);
    outtree->Branch("digicalo.row", &calorow);
    outtree->Branch("calo.time", &timestamp);
    outtree->Branch("calo.energy", &caloenergy);
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
    outtree->Branch("electron_gamma_correlated", &flag_e_g_correlated, "flag_e_g_correlated/O");

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

    //input tracker activity to replicate dead cells
    vector<long> tracker_activity = {};
    ifstream tracker_file("tracker_activity.txt");
    long n3;
    while (tracker_file >> n3) {
        tracker_activity.push_back(n3);
    }

    //output correlated events to text for plotting/testing
    ofstream eventtxt;
    eventtxt.open("e_gamma_events_SIM.txt");

    //check specific calorimeter for energy spectrum
    TH1D *spectrum = new TH1D("spectrum", "Energies for given OM", 100, 0, 5);
    TH1D *timehist = new TH1D("time", "OM to tracker delta_t", 120, -20, 100);
    TH1D *gamma_spectrum = new TH1D("gamma_energies", "Gamma Energies", 140, 0, 3.5);
    TH1D *zposhist = new TH1D("z_positions", "z-positions around OM centre", 100, -5, 5);
    TH1D *corr_hist = new TH1D("time_correlation", "delta t between electron and gamma", 200, -100, 100);

    //recording energy at only OM ID 123
    TH1D *spectrum_123 = new TH1D("spectrum_123", "Electron and gamma energies for OM 123", 100, 0, 5);

    int good_events = 0;
    int totalentries = tree->GetEntries();
    int cut_calohits = 0, cut_e_energy = 0, cut_OM_deltat = 0, cut_correlated = 0;

    //for z-pos test
    bool flag_cut_zpos;
    int cut_zpos = 0;
    double z_k = 0.01;
    double z_H = 2.95;

    for (int i=0; i < totalentries; i++) {
        tree->GetEntry(i);

        if (i % 10000 == 0) {cout << i << " out of " << totalentries << "\n";}

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
        flag_e_g_correlated = 0;
        flag_cut_zpos = 0;

        //temporary locations for time correlation
        int e_side = 0, e_row = 0, e_col = 0;
        int g_side = 0, g_row = 0, g_col = 0;

        if (calohits >= 2) {
            flag_cut_calohits = 1;   //passed first cut
        } else {
            outtree->Fill();
            continue;
        }

        int no_high_energy = 0;

        for (int j=0; j<calohits; j++) {            //for each hit calorimeter j

            int col = calocolumn->at(j); 
            float tcol_min = tab_column.at(col) - 4.; 
            float tcol_max = tab_column.at(col) + 4.;

            vector<vector<int>> *track = new vector<vector<int>>;
            track->clear();

            bool adj_tracker = 0;

            int hit_caloid = (13*col) + (260*caloside->at(j)) + calorow->at(j);  
            
            double hit_energy = 0.;
            if (calib[hit_caloid] != 0) {
                hit_energy = caloenergy->at(j);//(calocharge->at(j))/(-1000.);          //should omit dead OMs 
            }                       
        

            if (hit_caloid > 519) {continue;}    //main wall only

            if (hit_energy > 0.3) {         //energy cut
                no_high_energy += 1;
            } else {
                continue;
            }

            //calc OM z-range
            double z_min = -1.1165 + 0.18714*calorow->at(j) - 1;
            double z_max = -1.1165 + 0.18714*calorow->at(j) + 1;

            for (int k=0; k<trackercolumn->size(); k++) {
                if (trackerside->at(k) != caloside->at(j)) {continue;}          //check same side
                if (trackerlayer->at(k) < 7) {continue;}                        //check layer 8 or 9

                int tcol = trackercolumn->at(k);
                int tlayer = trackerlayer->at(k);
                int tside = trackerside->at(k);

                //skip trackers with low activity in the real run
                if (tracker_activity[tside*1017 + tcol*9 + tlayer] < 25000) {continue;}

                //enforce only some known good cells (z-pos check only)
                //if (tside != 0 || tcol < 9 || tcol > 37 || tlayer != 8) {continue;} revisit this later ----------------------

                if (tcol <= tcol_max && tcol >= tcol_min) {                    //within +- 4 range, start track 

                    double delta_t = trackertime->at(k) - timestamp->at(j);    //need to clarify this 
                    timehist->Fill(delta_t);

                    delta_t = 0; //TEMPORARY TO GET SOME DATA
                    
                    if ((delta_t > -0.2 && delta_t < 50)) {     //time cut 
                        flag_cut_OM_delta_t = 1;
                    }

                    //z-pos calculation
                    double t_top = (R6->at(k).at(0) - anode_R0->at(k).at(0))*12.5E-3;          // Put it in µs
                    double t_bottom = (R5->at(k).at(0) - anode_R0->at(k).at(0))*12.5E-3;
                    double z_gg = -99999; 

                    if (t_top > 0 && t_bottom > 0 && flag_cut_OM_delta_t == 1) {
                        double t_ratio = ((t_bottom - t_top)/(t_top + t_bottom));
                        //z_gg = t_ratio;

                        //newer calulation
                        z_gg = (z_H/2)*t_ratio - (z_k*(z_H*z_H)/4)*t_ratio*(1-abs(t_ratio));

                        zposhist->Fill(z_gg - (z_min + z_max)/2.);
                    }

                    if (z_gg >= z_min && z_gg <= z_max) {flag_cut_zpos = 1;}     //z-pos cut 

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

                g_side = caloside->at(j);
                g_row = calorow->at(j);
                g_col = calocolumn->at(j);
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

                        //dont include inactive trackers (need to think more on if this belongs here)
                        if (tracker_activity[trackerside->at(m)*1017 + trackercolumn->at(m)*9 + trackerlayer->at(m)] < 25000) {continue;}

                        if (within_x(track->at(l), next_tracker, 2)) {       //append to track vector if within 2x2 box
                            track->push_back(next_tracker);               
                        }
                    }
                }

                if (pre_size == track->size()) {        //stop reconstruction if not changing anymore
                    e_hit_time = timestamp->at(j);
                    e_hit_energy = hit_energy;
                    e_hit_caloid = hit_caloid;

                    e_side = caloside->at(j);
                    e_row = calorow->at(j);
                    e_col = calocolumn->at(j);

                    for (int p = 0; p<track->size(); p++) {
                        hit_track->push_back(track->at(p).at(0)*1017 + track->at(p).at(1)*9 + track->at(p).at(2));
                    }

                    break;  
                }
            }
        }

        if (no_high_energy == 2) {flag_cut_e_energy = 1;} //enforce only two particles above 0.3MeV
        if (hit_track->size() > 3) {
            flag_cut_tracklength = 1;
        }

        //add electron-gamma correlation here?
        if (flag_cut_calohits && flag_cut_e_energy && flag_cut_e_energy && flag_cut_tracklength) {
            //check for adjacency or same column
            if (within_x({e_side, e_row, e_col}, {g_side, g_row, g_col}, 1) == 0 && e_col != g_col) {
                //enforce time correlation 
                if (abs((e_hit_time - gamma_timestamp)) < 20) {
                    flag_e_g_correlated = 1; 

                    corr_hist->Fill((e_hit_time-gamma_timestamp));

                    //write to histograms
                    gamma_spectrum->Fill(gamma_energy);
                    spectrum->Fill(e_hit_energy);
                    
                    if (gamma_caloid == 123) {spectrum_123->Fill(gamma_energy);}
                    if (e_hit_caloid == 123) {spectrum_123->Fill(e_hit_energy);}

                    eventtxt << i << "\n";
                }
            }
        }

        //trying out a nested if to count passes of ALL cuts
        if (flag_cut_calohits == 1) {
            cut_calohits += 1;
            if (flag_cut_e_energy == 1) {
                cut_e_energy += 1;
                if (flag_cut_OM_delta_t == 1) {
                    cut_OM_deltat += 1;
                    if (flag_cut_tracklength == 1) {
                        good_events += 1;
                        if (flag_cut_zpos == 1) {
                            cut_zpos += 1;
                            if (flag_e_g_correlated == 1) {
                                cut_correlated += 1;
                            }
                        }
                    }
                }
            }
        }

        outtree->Fill();
    }

    out->cd();
    outtree->Write();

    //output number of events cut, some events may have multiple recorded tracks 
    cout << "Initial events:                 " << totalentries << "\n";
    cout << "Events with >= 2 OM hits:       " << cut_calohits << "\n";
    cout << "Events with two > 0.3MeV hits:  " << cut_e_energy << "\n";
    cout << "Events with -0.2 < dt < 50us:   " << cut_OM_deltat << "\n";
    cout << "Events with track length > 3:   " << good_events << "\n";
    cout << "Correlated z-position:          " << cut_zpos << "\n";
    cout << "Correlated electron and gamma:  " << cut_correlated << "\n\n";

    //quick text output to file 
    ofstream outtxt;
    outtxt.open("cuts_SIM.txt");
    outtxt << "Initial events:                " << totalentries << "\n";
    outtxt << "Events with >= 2 OM hits:      " << cut_calohits << "\n";
    outtxt << "Events with two > 0.3MeV hits: " << cut_e_energy << "\n";
    outtxt << "Events with -0.2 < dt < 50us:  " << cut_OM_deltat << "\n";
    outtxt << "Events with track length > 3:  " << good_events << "\n";
    outtxt << "Correlated z-position:         " << cut_zpos << "\n";
    outtxt << "Correlated electron and gamma: " << cut_correlated << "\n\n";
    outtxt.close();

    eventtxt.close();

    spectrum->SetTitle("Energy of Correlated Electrons;Energy (MeV);Count");
    timehist->SetTitle("Time difference between OM and adjacent tracker;delta_t (us);Count");
    gamma_spectrum->SetTitle("Energy of Correlated Photons;Energy (MeV);Count");
    spectrum_123->SetTitle("Energy of correlated electrons and photons on OM 123;Energy (MeV);Count");
    corr_hist->SetTitle("Time difference between electron and gamma OM;Time difference (ns);Count");

    spectrum->Write();
    timehist->Write();
    gamma_spectrum->Write();
    corr_hist->Write();
    spectrum_123->Write();

    corr_hist->Draw();
    //gamma_spectrum->Draw();
    //spectrum->Draw();
    //timehist->Draw();
    //zposhist->Draw();

}