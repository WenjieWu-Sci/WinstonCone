void PMTconfigure() {
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
    cout << line << endl;
    sfile.close();

    int m_copyno;
    double m_theta, m_phi;
    TString outfile = filename + ".csv";
    ofstream rfile;
    rfile.open(outfile,ios::out);
    rfile.setf(ios::fixed);
    for (int i=line-1;i>=0;i--) {
        m_copyno = copyno[i]; 
        for (int j=0;j<m_copyno;j++) {
            m_theta = theta[i]*180/TMath::Pi();
            m_phi = (phi[i] + j*2*TMath::Pi()/m_copyno)*180/TMath::Pi();
            rfile.precision(10);
            rfile << j << "\t" << m_theta << "\t";
            rfile.precision(10);
            rfile << m_phi << endl;
            rfile << endl;
        }
    }
    for (int i=line-1;i>0;i--) {
        m_copyno = copyno[i];
        for (int j=0;j<m_copyno;j++) {
            m_theta = 180 - theta[i]*180/TMath::Pi();
            m_phi = (phi[i] + j*2*TMath::Pi()/m_copyno)*180/TMath::Pi();
            rfile.precision(10);
            rfile << j << "\t" << m_theta << "\t";
            rfile.precision(10);
            rfile << m_phi << endl;
            rfile << endl;
        }
    }
    rfile.close();
}
