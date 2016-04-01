#include "RayExpPMTSD.hh"

#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4HCofThisEvent.hh"
#include "G4Track.hh"
#include "G4ParticleDefinition.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4LossTableManager.hh"
#include "G4Gamma.hh"
#include "G4Electron.hh"

RayExpPMTSD::RayExpPMTSD( G4String name ) : G4VSensitiveDetector(name) {
    G4String HCname = "PMTCollection";
    collectionName.insert(HCname);
}

RayExpPMTSD::~RayExpPMTSD() {
}

void RayExpPMTSD::Initialize(G4HCofThisEvent* HCE) {
    G4cout << "RayExpPMTSD::Initialize" << G4endl;
    PMTCollection = new RayExpPMTHitsCollection(SensitiveDetectorName,collectionName[0]);
    static G4int HCID = -1;
    if (HCID < 0) {
        HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
    }
    HCE->AddHitsCollection(HCID, PMTCollection);

    iftrig = false;
}

G4bool RayExpPMTSD::ProcessHits(G4Step* aStep, G4TouchableHistory*) {
    G4Track* fTrack = aStep->GetTrack();
    const G4DynamicParticle* fParticle = fTrack->GetDynamicParticle();
    RayExpPMTHit* hit = new RayExpPMTHit();

    energy = fTrack->GetDynamicParticle()->GetTotalEnergy();
  //  G4cout << "This is the energy of this track : " << energy << G4endl;
    hit->SetPDG(fTrack->GetDefinition()->GetPDGEncoding());
    G4int m_ihit = PMTCollection->insert(hit);
    G4cout << "This is the " << m_ihit << "th Hit" << G4endl;

    G4double edep = aStep->GetTotalEnergyDeposit()/keV;
    hit->SetEdep(edep);
    //if ( edep==0 ) {
    //    G4cout << "edep = 0" << G4endl;
    //    return false;
    //}
    G4String Name = fTrack->GetDefinition()->GetParticleName();
    N = fTrack->GetCurrentStepNumber(); 
    TID = fTrack->GetTrackID();
    hit->SetTrackID(TID);
    PID = fTrack->GetParentID();
    StepLength = fTrack->GetStepLength()/mm;
    velocity = fTrack->GetVelocity();
    
    // check whether it is entering to the entrance aperture
    G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
    G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();
    G4ThreeVector pre = thePrePoint->GetPosition()/mm;
    G4double PreZ = pre.getZ();
    G4String thePrePVname = thePrePV->GetName();

    G4StepPoint* thePostPoint = aStep->GetPostStepPoint();
    G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();
    G4ThreeVector post = thePostPoint->GetPosition()/mm;
    G4double PostZ = post.getZ();
    G4String thePostPVname = thePostPV->GetName();

    hit->SetPreZcoordinate(PreZ);
    hit->SetPostZcoordinate(PostZ);
    hit->SetPrePVname(thePrePVname);
    hit->SetPostPVname(thePostPVname);
    //if (thePrePVname == "entra_phys" && thePostPVname == "World") {
    //    EnterPhoton++;
    //}
    //if (thePrePVname == "collec_phys" && thePostPVname == "World") {
    //    ExitPhoton++;
    //}

//    G4cout << "Name                : " << Name          << G4endl;
//    G4cout << "Current Step Number : " << N             << G4endl;
//    G4cout << "Track ID            : " << TID           << G4endl;
//    G4cout << "Parent ID           : " << PID           << G4endl;
//    G4cout << "edep                : " << edep          << G4endl;
//    G4cout << "PreStep             : " << pre           << G4endl;
//    G4cout << "PostStep            : " << post          << G4endl;
//    G4cout << "StepLength          : " << StepLength    << G4endl;
//    G4cout << "Velocity            : " << velocity      << G4endl;
//    G4cout << "thePrePV            : " << thePrePVname  << G4endl;
//    G4cout << "thePostPV           : " << thePostPVname << G4endl;

    return true;
}

void RayExpPMTSD::EndOfEvent(G4HCofThisEvent* HCE) {
    G4cout << "RayExpPMTSD::EndofEvent" << G4endl;
}
