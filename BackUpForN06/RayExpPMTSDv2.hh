#ifndef RAYEXPPMTSDv2_h
#define RAYEXPPMTSDv2_h 1

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "RayExpPMTHitv2.hh"

class RayExpPMTSDv2 : public G4VSensitiveDetector {
    public:
        RayExpPMTSDv2(G4String name);
        ~RayExpPMTSDv2();

        void Initialize(G4HCofThisEvent* HCE);
        G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist);
        void EndOfEvent(G4HCofThisEvent* HCE);

    private:
        RayExpPMTHitsCollectionv2* PMTCollectionv2;
        bool iftrig;
        double energy;

    private:
        G4int N;
        G4double TID;
        G4double PID;
        G4double StepLength;
        G4double velocity;
};

#endif
