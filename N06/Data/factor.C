#include <fstream>
#include <vector>
#include <iostream>
#include "TString.h"
#include "TCanvas.h"
#include "TGraph.h"
#include "TAxis.h"

using namespace std;

int main (int argc, char* argv[]) {
    ifstream sfile;
    TString filename(argv[1]);
    sfile.open(filename,ios::in);
    if (!sfile.good()) {
        cout << "Can not open this file!" << endl;
        return -1;
    }
    double angle, enterNum, exitNum;
    vector<double> theta, Factor;
    int line = 0;
    while (!sfile.eof()) {
        sfile >> angle >> enterNum >> exitNum;
        if (!sfile.eof()) {
            line++;
            theta.push_back(angle);
            Factor.push_back(exitNum/enterNum);
        }
    }
    TCanvas* cc = new TCanvas("cc","cc",800,600);
    TGraph* tg = new TGraph(line,&theta[0],&Factor[0]);
    tg->SetMarkerStyle(8);
    tg->SetMarkerSize(1);
    tg->GetXaxis()->SetTitle("angle (deg)");
    tg->GetYaxis()->SetTitle("Concentration Factor");
    TString output = filename.Remove(filename.Length()-4);
    tg->Draw("APL");
    cc->SaveAs(output+".png");
}
