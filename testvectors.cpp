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

void testvectors() {
    TFile *f = new TFile("1166_out.root", "read");

    TTree *T;
    gDirectory->GetObject("arbre_alpha", T);

    //try using a vector 
    vector<double> *vec = new vector<double>; //remember this!
    T->SetBranchStatus("*", 0);
    T->SetBranchStatus("energy_calo_gamma", 1);
    T->SetBranchAddress("energy_calo_gamma", &vec);

    int nentries = T->GetEntries();
    TH1D *hist = new TH1D("hist", "mag hist", 100, -1, 2);

    for (int j=0;j<10;j++){
        T->GetEntry(j);

        for (int i=0;i<vec->size();i++) {
            cout << vec->at(i) << "\n";
            hist->Fill(vec->at(i));
            //hist->Fill(i);

        }
    }

    cout << "finished\n";
    hist->Draw();                    //crashing
}