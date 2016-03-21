void PMTPosition() {
    TString filename = "25-18345-PMT-Location";
    TString inputfile = filename + ".txt";
    ifstream sfile;
    sfile.open(inputfile,ios::in);
    int line = 0;
    vector<int> copyno;
    vector<double> theta, phi;
    while (!sfile.eof()) {
        int tmp_lineno, tmp_copyno;
        double tmp_theta, tmp_phi;
        sfile >> tmp_lineno >> tmp_copyno >> tmp_theta >> tmp_phi;
        if (!sfile.eof()) {
            copyno.push_back(tmp_copyno);
            theta.push_back(tmp_theta);
            phi.push_back(tmp_phi);
            line++;
        }
    }
    sfile.close();

    int m_copyno;
    int m_pmtID;
    double m_theta, m_phi;

    TFile* PMTPos = new TFile(filename+".root","RECREATE");
    TTree* PMTPosition = new TTree("PMTPosition","PMTPosition");
    PMTPosition->Branch("pmtID",&m_pmtID,"m_pmtID/I");
    PMTPosition->Branch("theta",&m_theta,"m_theta/D");
    PMTPosition->Branch("phi",&m_phi,"m_phi/D");

    m_pmtID = 0;
    for (int i=line-1;i>=0;i--) {
        m_copyno = copyno[i]; 
        for (int j=0;j<m_copyno;j++) {
            m_theta = theta[i]*180/TMath::Pi();
            m_phi = (phi[i] + j*2*TMath::Pi()/m_copyno)*180/TMath::Pi();
            PMTPosition->Fill();
            m_pmtID++;
        }
    }
    for (int i=line-1;i>0;i--) {
        m_copyno = copyno[i];
        for (int j=0;j<m_copyno;j++) {
            m_theta = 180 - theta[i]*180/TMath::Pi();
            m_phi = (phi[i] + j*2*TMath::Pi()/m_copyno)*180/TMath::Pi();
            PMTPosition->Fill();
            m_pmtID++;
        }
    }
    PMTPosition->Write();
}
