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
checks a given file for electron+gamma like events
electron and gamma are time correlated by < 50ns 

records OM location/energy/timestamp for both at a given event
and records basic track reconstruction of electron

E. Telfer 2024
*/

bool within_x(vector<int> a, vector<int> b, int x) {
    //check if 3D vectors are within a x by x box on the same side. takes {side, col, layer}. 
    if (a[0] != b[0]) {return 0;}
    if (a[1] < (b[1]-x) || a[1] > (b[1]+x)) {return 0;}
    if (a[2] < (b[2]-x) || a[2] > (b[2]+x)) {return 0;}
    return 1;
}

void e_gamma() {

    //get tree and setup relevant branches
    TFile *f = new TFile("snemo_run-1101_udd.root", "READ");                    //change to whichever run required 
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
    std::vector<vector<long>> *anode_R1 = new std::vector<vector<long>>;
    std::vector<vector<long>> *anode_R2 = new std::vector<vector<long>>;
    std::vector<vector<long>> *anode_R3 = new std::vector<vector<long>>;
    std::vector<vector<long>> *anode_R4 = new std::vector<vector<long>>;
    std::vector<vector<long>> *R5 = new std::vector<vector<long>>;
    std::vector<vector<long>> *R6 = new std::vector<vector<long>>;
    tree->SetBranchStatus("digitracker.anodetimestampR0", 1);
    tree->SetBranchAddress("digitracker.anodetimestampR0", &anode_R0);
    tree->SetBranchStatus("digitracker.anodetimestampR1", 1);
    tree->SetBranchAddress("digitracker.anodetimestampR1", &anode_R1);
    tree->SetBranchStatus("digitracker.anodetimestampR2", 1);
    tree->SetBranchAddress("digitracker.anodetimestampR2", &anode_R2);
    tree->SetBranchStatus("digitracker.anodetimestampR3", 1);
    tree->SetBranchAddress("digitracker.anodetimestampR3", &anode_R3);
    tree->SetBranchStatus("digitracker.anodetimestampR4", 1);
    tree->SetBranchAddress("digitracker.anodetimestampR4", &anode_R4);
    tree->SetBranchStatus("digitracker.bottomcathodetimestamp", 1);
    tree->SetBranchAddress("digitracker.bottomcathodetimestamp", &R5);
    tree->SetBranchStatus("digitracker.topcathodetimestamp", 1);
    tree->SetBranchAddress("digitracker.topcathodetimestamp", &R6);
    
    //falling cell time for time calibration
    vector<int> *falling_cell = new vector<int>;
    tree->SetBranchStatus("digicalo.falling_cell", 1);
    tree->SetBranchAddress("digicalo.falling_cell", &falling_cell);

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
    bool flag_e_g_correlated = 0;           //time correlation between electron and gamma-like hits

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
    vector<double> calib = {};
    ifstream calib_file("run-1351_fee-charge-to-energy.txt");
    int n1;
    double n2;
    while (calib_file >> n1 >> n2) {
        calib.push_back(n2);
    }

    //import time calibration
    vector<double> time_calibration = {};
    ifstream time_file("output_data_merged_2.txt");
    int n3;
    double n4;
    double n5;    //std error on the calibration, unused for now
    while(time_file >> n3 >> n4 >> n5) {
        time_calibration.push_back(n4);
    }

    //output correlated events to text for plotting/testing
    ofstream eventtxt;
    eventtxt.open("e_gamma_events.txt");

    //check specific calorimeter for energy spectrum
    TH1D *spectrum = new TH1D("spectrum", "Energies for correlated electrons", 100, 0, 5);
    TH1D *timehist = new TH1D("time", "OM to tracker delta_t", 120, -20, 100);
    TH1D *gamma_spectrum = new TH1D("gamma_energies", "Gamma Energies", 140, 0, 3.5);
    TH1D *zposhist = new TH1D("z_positions", "z-positions around OM centre", 100, -5, 5);
    TH1D *trackerhist = new TH1D("trackers", "tracker activity distribution", 2034, 0, 2034);
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

    //record first and last timestamp for total runtime
    long time0 = 0, time1 = 0; 

    //record active trackers
    long tracker_array[2034] = {}; 
    int affected = 0, unaffected = 0;

    for (int i=0; i < 50000; i++) {
        tree->GetEntry(i);

        //monitor tracker activity before any cuts, might slow code quite a bit
        for (int J=0; J<trackercolumn->size(); J++) {
            tracker_array[trackerside->at(J)*1017 + trackercolumn->at(J)*9 + trackerlayer->at(J)] += 1;
            trackerhist->Fill(trackerside->at(J)*1017 + trackercolumn->at(J)*9 + trackerlayer->at(J));
        }

        if (i % 10000 == 0) {cout << i << " out of " << totalentries << "\n";}

        //record time of run 
        if (i == 0) {time0 = timestamp->at(0);}
        if (i == (totalentries-1)) {time1 = timestamp->at(timestamp->size()-1);}

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

        //for z=pos investigation
        vector<int> tracker_zpos = {0, 0};

        //trying new energy cut strategy
        int no_high_energy = 0;

        int e_calo_index = 9999;
        int gamma_calo_index = 9999;

        for (int j=0; j<calohits; j++) {            //for each hit calorimeter j

            int col = calocolumn->at(j); 
            float tcol_min = tab_column.at(col) - 4.; 
            float tcol_max = tab_column.at(col) + 4.;

            vector<vector<int>> *track = new vector<vector<int>>;
            track->clear();

            bool adj_tracker = 0;

            int hit_caloid = (13*col) + (260*caloside->at(j)) + calorow->at(j);          
            double hit_energy = (charge->at(j))*calib[hit_caloid]*energy_conv;             //changed to MeV and flipped sign

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
                int tcol = trackercolumn->at(k);
                int tlayer = trackerlayer->at(k);
                int tside = trackerside->at(k);

                if (trackerside->at(k) != caloside->at(j)) {continue;}          //check same side
                if (trackerlayer->at(k) < 7) {continue;}                        //check layer 8 or 9

                //enforce only some known good cells (z-pos check only)
                //if (tside != 0 || tcol < 9 || tcol > 37 || tlayer != 8) {continue;} revisit this later ----------------------

                if (tcol <= tcol_max && tcol >= tcol_min) {                    //within +- 4 range, start track 

                    long delta_t = (2.*anode_R0->at(k).at(0) - timestamp->at(j))*6.25/1000.;         //microseconds
                    timehist->Fill(delta_t);

                    if ((delta_t > -0.2 && delta_t < 50)) {     //time cut 
                        flag_cut_OM_delta_t = 1;
                    }

                    //check for good timestamps
                    long top_timestamp = 0, bottom_timestamp = 0;
                    if (R6->at(k).at(0) > 0) {top_timestamp = R6->at(k).at(0);} else {
                        //substitute with anode times within ~5us for now 
                        if ((abs(R5->at(k).at(0)-anode_R1->at(k).at(0)) < 500) || (abs(R5->at(k).at(0)-anode_R3->at(k).at(0)) < 500)) {
                            if (R5->at(k).at(0) > 0) {
                                if (anode_R2->at(k).at(0) > 0) {top_timestamp = anode_R2->at(k).at(0);}
                                else if (anode_R4->at(k).at(0) > 0) {top_timestamp = anode_R4->at(k).at(0);}
                            }
                        }
                        if ((abs(R5->at(k).at(0)-anode_R2->at(k).at(0)) < 500) || (abs(R5->at(k).at(0)-anode_R4->at(k).at(0)) < 500)) {
                            if (R5->at(k).at(0) > 0) {
                                if (anode_R1->at(k).at(0) > 0) {top_timestamp = anode_R1->at(k).at(0);}
                                else if (anode_R3->at(k).at(0) > 0) {top_timestamp = anode_R3->at(k).at(0);}
                            }
                        }
                    } 
                    if (R5->at(k).at(0) > 0) {bottom_timestamp = R5->at(k).at(0);} else {
                        //substitute with anode times within ~5us for now 
                        if ((abs(R6->at(k).at(0)-anode_R1->at(k).at(0)) < 500) || (abs(R6->at(k).at(0)-anode_R3->at(k).at(0)) < 500)) {
                            if (R6->at(k).at(0) > 0) {
                                if (anode_R2->at(k).at(0) > 0) {bottom_timestamp = anode_R2->at(k).at(0);}
                                else if (anode_R4->at(k).at(0) > 0) {bottom_timestamp = anode_R4->at(k).at(0);}
                            }
                        }
                        if ((abs(R6->at(k).at(0)-anode_R2->at(k).at(0)) < 500) || (abs(R6->at(k).at(0)-anode_R4->at(k).at(0)) < 500)) {
                            if (R6->at(k).at(0) > 0) {
                                if (anode_R1->at(k).at(0) > 0) {bottom_timestamp = anode_R1->at(k).at(0);}
                                else if (anode_R3->at(k).at(0) > 0) {bottom_timestamp = anode_R3->at(k).at(0);}
                            }
                        }  
                    }

                    //z-pos calculation
                    double t_top = (top_timestamp - anode_R0->at(k).at(0))*12.5E-3;          // Put it in µs
                    double t_bottom = (bottom_timestamp - anode_R0->at(k).at(0))*12.5E-3;
                    double z_gg = -99999; 

                    if (t_top > 0 && t_bottom > 0 && flag_cut_OM_delta_t == 1) {
                        double t_ratio = ((t_bottom - t_top)/(t_top + t_bottom));
                        //z_gg = t_ratio;

                        //newer calulation
                        z_gg = (z_H/2)*t_ratio - (z_k*(z_H*z_H)/4)*t_ratio*(1-abs(t_ratio));

                        zposhist->Fill(z_gg - (z_min + z_max)/2.);
                    }

                    if (z_gg >= z_min && z_gg <= z_max) {
                        flag_cut_zpos = 1;                      //z-pos cut 
                        tracker_zpos[1] = 1;                    //for z-pos investigation
                    }    

                    //for z-pos investigation
                    if (tside == 0 && tcol == 14 && tlayer == 8) {
                        tracker_zpos[0] = 1;
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
                gamma_calo_index = j;

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

                        if (within_x(track->at(l), next_tracker, 2)) {       //append to track vector if within 2x2 box
                            track->push_back(next_tracker);               
                        }
                    }
                }

                if (pre_size == track->size()) {        //stop reconstruction if not changing anymore
                    e_hit_time = timestamp->at(j);
                    e_hit_energy = hit_energy;
                    e_hit_caloid = hit_caloid;
                    e_calo_index = j;

                    e_side = caloside->at(j);
                    e_row = calorow->at(j);
                    e_col = calocolumn->at(j);

                    for (int p = 0; p<track->size(); p++) {
                        //record track
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

                double delta_T = 999999;

                if (e_hit_time > 0 && gamma_timestamp > 0) {
                //apply time calibration correction
                    double calib_e_time = 6.25*e_hit_time - time_calibration[e_hit_caloid] + ((falling_cell->at(e_calo_index)/256)*(0.390625));                 
                    double calib_gamma_time = 6.25*gamma_timestamp - time_calibration[gamma_caloid] + ((falling_cell->at(gamma_calo_index)/256)*(0.390625));

                    delta_T = (calib_e_time - calib_gamma_time);            //calibrated time difference
                    corr_hist->Fill(delta_T);
                }

                //enforce time correlation
                if (abs(delta_T) < 50) {
                    flag_e_g_correlated = 1; 
                    cut_correlated += 1;

                    //fill histograms
                    gamma_spectrum->Fill(gamma_energy);
                    spectrum->Fill(e_hit_energy);

                    if (gamma_caloid == 123) {spectrum_123->Fill(gamma_energy);}
                    if (e_hit_caloid == 123) {spectrum_123->Fill(e_hit_energy);}

                    eventtxt << i << "\n";

                    //z-pos investigation
                    if (tracker_zpos[0] == 1) {
                        if (flag_cut_zpos == 1) {unaffected += 1;} else {affected += 1;}
                    }
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
    cout << "time of run: " << (time1-time0)*6.25/1000000000. << "s\n\n";
    cout << "z-pos test: " << affected << " of " << unaffected << "\n";

    //quick text output to file 
    ofstream outtxt;
    outtxt.open("cuts.txt");
    outtxt << "Initial events:                " << totalentries << "\n";
    outtxt << "Events with >= 2 OM hits:      " << cut_calohits << "\n";
    outtxt << "Events with two > 0.3MeV hits: " << cut_e_energy << "\n";
    outtxt << "Events with -0.2 < dt < 50us:  " << cut_OM_deltat << "\n";
    outtxt << "Events with track length > 3:  " << good_events << "\n";
    outtxt << "Correlated z-position:         " << cut_zpos << "\n";
    outtxt << "Correlated electron and gamma: " << cut_correlated << "\n\n";
    outtxt << "time of run: " << (time1-time0)*6.25/1000000000. << "s\n\n";
    outtxt << "z-pos: " << affected << " events cut, " << unaffected << "\n";
    outtxt.close();

    eventtxt.close();

    //write tracker activity
    ofstream trackertxt;
    trackertxt.open("tracker_activity.txt");
    for (int i=0; i<2034; i++) {
        trackertxt << tracker_array[i] << "\n";
    }
    trackertxt.close();

    spectrum->SetTitle("Energy of Correlated Electrons;Energy (MeV);Count");
    timehist->SetTitle("Time difference between OM and adjacent tracker;delta_t (us);Count");
    gamma_spectrum->SetTitle("Energy of Correlated Photons;Energy (MeV);Count");
    spectrum_123->SetTitle("Energy of correlated electrons and photons on OM 123;Energy (MeV);Count");
    trackerhist->SetTitle("Tracker activity over entire run;Tracker ID;Count");
    corr_hist->SetTitle("Time difference between electron and gamma OM;Time difference (ns);Count");
    spectrum->Write();
    timehist->Write();
    gamma_spectrum->Write();
    corr_hist->Write();
    spectrum_123->Write();
    trackerhist->Write();

    corr_hist->Draw();
    //gamma_spectrum->Draw();
    //spectrum->Draw();
    //timehist->Draw();
    //zposhist->Draw();

}