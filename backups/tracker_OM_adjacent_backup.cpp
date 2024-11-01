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
    TFile *f = new TFile("snemo_run-1166_udd.root", "read");
    TTree *tree = (TTree*)f->Get("SimData");

    int event = 0;
    int calohits = 0;
    int trackerhits = 0;
    vector<int> *caloside = new vector<int>;
    vector<int> *calocolumn = new vector<int>;
    vector<int> *calorow = new vector<int>;
    vector<int> *trackerside = new vector<int>;
    vector<int> *trackercolumn = new vector<int>;
    vector<int> *trackerlayer = new vector<int>;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("header.eventnumber", 1);
    tree->SetBranchAddress("header.eventnumber", &event);
    tree->SetBranchStatus("digicalo.nohits", 1);                    //calorimeters
    tree->SetBranchAddress("digicalo.nohits", &calohits);
    tree->SetBranchStatus("digicalo.side", 1);
    tree->SetBranchAddress("digicalo.side", &caloside);
    tree->SetBranchStatus("digicalo.column", 1);
    tree->SetBranchAddress("digicalo.column", &calocolumn);
    tree->SetBranchStatus("digicalo.row", 1);
    tree->SetBranchAddress("digicalo.row", &calorow);
    tree->SetBranchStatus("digitracker.nohits", 1);                    //trackers
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

    //comparison vector, +- 4 inclusive gives tracker columns 
    vector<double> tab_column = {0.8, 6.8, 12.5, 18.2, 24.2, 29.9, 35.9, 41.8, 46.7, 53.6, 59.4, 65.3, 71.1, 77, 83.1, 88.8, 94.8, 100.6, 106.4, 112.2};

    //check calorimeters for nearby active tracker cells
    int good_events = 0;
    int almost_events = 0;

    for (int i=0; i < totalentries; i++) {
        int entryno = elist->GetEntry(i);
        tree->GetEntry(entryno);

        for (int j=0; j<calohits; j++) {          //for each hit calorimeter
            int col = calocolumn->at(j); 
            float tcol_min = tab_column.at(col) - 4;            //adjacent tracker column range
            float tcol_max = tab_column.at(col) + 4;

            for (int k=0; k<trackercolumn->size(); k++) {
                if (trackerside->at(k) != caloside->at(j)) {continue;}          //check same side
                if (trackerlayer->at(k) < 7) {continue;}                        //check layer 8 or 9
                almost_events += 1;

                int tcol = trackercolumn->at(k);
                if (tcol <= tcol_max && tcol >= tcol_min) {
                    //do something or output something...
                    good_events += 1;
                    break;                                                   
                }
            }
        }
    }

    cout << "Almost events: " << almost_events << "\n";
    cout << "Adjacent events: " << good_events << "\n";

}