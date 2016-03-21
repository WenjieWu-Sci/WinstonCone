#include "RayExpPMTHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<RayExpPMTHit> RayExpPMTHitAllocator;

RayExpPMTHit::RayExpPMTHit() {}

RayExpPMTHit::~RayExpPMTHit() {}

RayExpPMTHit::RayExpPMTHit(const RayExpPMTHit& right) : G4VHit() {
  m_detType   = right.m_detType;
  m_detID     = right.m_detID;
  m_trackID   = right.m_trackID;
  m_pdg       = right.m_pdg;
  m_KE        = right.m_KE;
  m_edep      = right.m_edep;
  m_quenchedE = right.m_quenchedE;
}

const RayExpPMTHit& RayExpPMTHit::operator = (const RayExpPMTHit& right) {
  m_detType   = right.m_detType;
  m_detID     = right.m_detID;
  m_trackID   = right.m_trackID;
  m_pdg       = right.m_pdg;
  m_KE        = right.m_KE;
  m_edep      = right.m_edep;
  m_quenchedE = right.m_quenchedE;
    
  return *this;
}

int RayExpPMTHit::operator == (const RayExpPMTHit& right) const {
    return (this==&right) ? 1 : 0;
}
