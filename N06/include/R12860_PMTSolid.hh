#ifndef R12860_PMTSolid_hh
#define R12860_PMTSolid_hh

#include "globals.hh"

class G4VSolid;

class R12860_PMTSolid {

public:
    R12860_PMTSolid(double R1, double R1z, double R2, double R3, double H, double H3);
    ~R12860_PMTSolid();

    G4VSolid* GetSolid(G4String solidname, double thickness=0.0);

    G4VSolid* GetCoverTop(G4String solidname, 
                            G4double c_h1,
                            G4double c_r1,
                            G4double c_h2,
                            G4double c_r2,
                            G4double thickness=1*mm);
    G4VSolid* GetCoverTopCurve(G4String solidname,
                            G4double c_h1,
                            G4double c_r1,
                            G4double c_h2,
                            G4double c_r2,
                            G4double center2equator,
                            G4double ball_r);
    G4VSolid* GetCoverBottom(G4String solidname, 
                            G4double c_h2,
                            G4double c_r2,
                            G4double c_h3,
                            G4double c_r3,
                            G4double thick_out=11*mm,
                            G4double thick_in=1*mm);

    G4VSolid* GetCoverSolid(G4String solidname, 
                            G4double c_h1,
                            G4double c_r1,
                            G4double c_h2,
                            G4double c_r2,
                            G4double c_h3,
                            G4double c_r3,
                            G4double thickness=11.01*mm);
    G4VSolid* GetCoverSolidWithTopCurve(G4String solidname, 
                            G4double c_h0,
                            G4double c_r0,
                            G4double c_h1,
                            G4double c_r1,
                            G4double c_h2,
                            G4double c_r2,
                            G4double c_h3,
                            G4double c_r3,
                            G4double center2equator,
                            G4double ball_r,
                            G4double thickness=11.01*mm);

    void SetCoverEdge(G4int n) {
        numSide = n;
    }

private:
    double m_R1;
    double m_R1z;
    double m_R1p;
    double m_R2;
    double m_R3;
    double m_theta;
    double m_H;
    double m_H_1_2;

    double m_H3;

    static double delta_torlerance;
    G4int numSide;

};

#endif
