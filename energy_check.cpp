void energy_check(){
    TFile *f = new TFile("snemo_run-1166_udd.root", "READ");
    TTree *tree = (TTree*)f->Get("SimData");

    vector<int> *charge = new vector<int>;
    vector<int> *col = new vector<int>;
    vector<int> *side = new vector<int>;
    vector<int> *row = new vector<int>;
    int calohits = 0;
    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("digicalo.charge", 1);
    tree->SetBranchAddress("digicalo.charge", &charge);
    tree->SetBranchStatus("digicalo.nohits", 1);                   
    tree->SetBranchAddress("digicalo.nohits", &calohits);
    tree->SetBranchStatus("digicalo.side", 1);
    tree->SetBranchAddress("digicalo.side", &side);
    tree->SetBranchStatus("digicalo.column", 1);
    tree->SetBranchAddress("digicalo.column", &col);
    tree->SetBranchStatus("digicalo.row", 1);
    tree->SetBranchAddress("digicalo.row", &row);

    //input charge-energy calibration from text
    vector<double> calib;
    ifstream calib_file("run-1351_fee-charge-to-energy.txt");
    int n1;
    double n2;
    while (calib_file >> n1 >> n2) {
        calib.push_back(n2);
    }

    int entries = tree->GetEntries();
    int caloid;
    double energy;

    TH1D *totalspectrum = new TH1D("totalspectrum", "all OM energies", 100, 0, 10);

    for (int i=0;i<entries;i++) {
        tree->GetEntry(i); 
        for (int j=0;j<calohits;j++) {
            caloid = (13*col->at(j)) + (260*side->at(j)) + row->at(j);
            energy = charge->at(j)*calib[caloid]*(-1./1000.);
            if (energy > 0.3) {totalspectrum->Fill(energy);}
        }
        
    }

    totalspectrum->Draw();
}