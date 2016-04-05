#include "RayExpRun.hh"
#include "G4Event.hh"
#include "G4HCofThisEvent.hh"
#include "G4SDManager.hh"
//#include "RayExpPMTSD.hh"
//#include "RayExpPMTHit.hh"
#include "dywSD_PMT_v2.hh"
#include "dywHit_PMT.hh"
#include <iostream>
#include <fstream>
//#include "RayExpPMTSDv2.hh"
//#include "RayExpPMTHitv2.hh"

RayExpRun::RayExpRun() {
//    EnterHCID = -1;
    ExitHCID = -1;
//    HitEnterNum = 0;
    HitExitNum = 0;
//    EnterNum = 0;
    ExitNum = 0;
}

RayExpRun::~RayExpRun() {
}

void RayExpRun::RecordEvent(const G4Event* evt) {
//    if (EnterHCID<0) {
//        EnterHCID = G4SDManager::GetSDMpointer()->GetCollectionID("PMTCollection");
//        G4cout << "entranSD : " << EnterHCID << G4endl;
//    }
    if (ExitHCID<0) {
        ExitHCID = G4SDManager::GetSDMpointer()->GetCollectionID("hitCollection");
        G4cout << "ExitSD : " << ExitHCID << G4endl;
    }

//    RayExpPMTHitsCollection* EnterCHC = 0;
    dywHit_PMT_Collection* ExitCHC = 0;
//    RayExpPMTHitsCollectionv2* ExitCHC = 0;
    G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
    if (HCE) {
        G4cout << "Get hit collection" << G4endl;
//        EnterCHC = static_cast<RayExpPMTHitsCollection*>(HCE->GetHC(EnterHCID));
        ExitCHC = static_cast<dywHit_PMT_Collection*>(HCE->GetHC(ExitHCID));
    }
//    if (EnterCHC) {
//        G4cout << "Get enter hit collection" << G4endl;
//        for (G4int i=0;i<EnterCHC->entries();i++) {
//            G4String thePrePVname = (*EnterCHC)[i]->GetPrePVname();
//            G4String thePostPVname = (*EnterCHC)[i]->GetPostPVname();
//            G4double thePreZ = (*EnterCHC)[i]->GetPreZcoordinate();
//            G4double thePostZ = (*EnterCHC)[i]->GetPostZcoordinate();
//            G4cout << "Enter aperture : " << G4endl;
//            G4cout << "thePrePVname : " << thePrePVname << G4endl;
//            G4cout << "thePostPVname : " << thePostPVname << G4endl;
//            HitEnterNum++;
//            if (thePrePVname=="entra_phys" && thePostPVname=="World" && thePreZ > thePostZ) {
//                EnterNum++;
//            }
//        }
//    }
    if (ExitCHC) {
        G4cout << "Get PMT hit collection" << G4endl;
        G4cout << ExitCHC->entries() << G4endl;
        for (G4int i=0;i<ExitCHC->entries();i++) {
            G4String thePrePVname = (*ExitCHC)[i]->GetPrePVname();
            G4String thePostPVname = (*ExitCHC)[i]->GetPostPVname();
            G4double thePreZ = (*ExitCHC)[i]->GetPreZcoordinate();
            G4double thePostZ = (*ExitCHC)[i]->GetPostZcoordinate();
            G4cout << "pmt : " << G4endl;
            G4cout << "thePrePVname : " << thePrePVname << G4endl;
            G4cout << "thePostPVname : " << thePostPVname << G4endl;
            HitExitNum++;
            if (thePrePVname=="body_phys" && thePostPVname=="inner1_phys" ) {
                ExitNum++;
            }
        }
    }

    G4cout << "##################################################" << G4endl;
//    G4cout << HitEnterNum << G4endl;
//    G4cout << "Number of entering photon : " << EnterNum << G4endl;
    G4cout << HitExitNum << G4endl;
    G4cout << "Number of exiting photon : " << ExitNum << G4endl;
    G4cout << "##################################################" << G4endl;

//    std::ofstream sfile;
//    sfile.open("WithoutLC.txt",std::ios::app);
//    sfile << ExitNum << G4endl;
//    sfile.close();

    G4Run::RecordEvent(evt);
}
