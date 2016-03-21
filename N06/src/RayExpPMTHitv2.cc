#include "RayExpPMTHitv2.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

G4Allocator<RayExpPMTHitv2> RayExpPMTHitv2Allocator;

RayExpPMTHitv2::RayExpPMTHitv2() {}

RayExpPMTHitv2::~RayExpPMTHitv2() {}

RayExpPMTHitv2::RayExpPMTHitv2(const RayExpPMTHitv2& right) : G4VHit() {
  m_detType   = right.m_detType;
  m_detID     = right.m_detID;
  m_trackID   = right.m_trackID;
  m_pdg       = right.m_pdg;
  m_KE        = right.m_KE;
  m_edep      = right.m_edep;
  m_quenchedE = right.m_quenchedE;
}

const RayExpPMTHitv2& RayExpPMTHitv2::operator = (const RayExpPMTHitv2& right) {
  m_detType   = right.m_detType;
  m_detID     = right.m_detID;
  m_trackID   = right.m_trackID;
  m_pdg       = right.m_pdg;
  m_KE        = right.m_KE;
  m_edep      = right.m_edep;
  m_quenchedE = right.m_quenchedE;
    
  return *this;
}

int RayExpPMTHitv2::operator == (const RayExpPMTHitv2& right) const {
    return (this==&right) ? 1 : 0;
}
