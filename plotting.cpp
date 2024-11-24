void plotting() {
    TFile *infile = new TFile("tracks.root", "READ");
    TTree *tree = (TTree*)infile->Get("tracks");

    int g_OM = 0;
    bool corr;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("gamma_caloid", 1);
    tree->SetBranchAddress("gamma_caloid", &g_OM);
    tree->SetBranchStatus("electron_gamma_correlated", 1);
    tree->SetBranchAddress("electron_gamma_correlated", &corr);

    TH1D *calocolumns = new TH1D("columns", "Hits per OM column", 20, 0, 20);

    int col;
    for (int i=0; i<tree->GetEntries(); i++) {
        tree->GetEntry(i);

        if (corr == 0) {continue;}

        if (g_OM/260 == 0) {
            col = g_OM/13;
        } else {
            col = (g_OM - 260)/13;
        }

        calocolumns->Fill(col);
    }

    calocolumns->Draw(); 
}