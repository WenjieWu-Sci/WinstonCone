//------------------------------------------------------------------
//                       dywHit_PMT
//                       
// The data members of a hit on sensitive detector is defined.
// Their values are obtained in dywSD_PMT using the information of
// steps which hit the sensitive detector.
// -----------------------------------------------------------------
//  Author: Liang Zhan,  2006/01/27
//  Modified by: Weili Zhong, 2006/03/01
// -----------------------------------------------------------------

#ifndef dywHit_PMT_h
#define dywHit_PMT_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

////////////////////////////////////////////////////////////////////////

class dywHit_PMT : public G4VHit
{
  public:
//    dywHit_PMT(G4int z);  
//    dywHit_PMT(G4int i,G4double t);
    dywHit_PMT();
    ~dywHit_PMT();
    dywHit_PMT(const dywHit_PMT &);
    const dywHit_PMT& operator=(const dywHit_PMT &);
    int operator==(const dywHit_PMT &/*right*/) const; //{return 0;};
//    void Init();

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);

//    // for visulization of hits
//    void Draw();
//    void Print();
//    void Dump(const char* prefix);
    

  private:
    
    G4String      m_thePrePVname;
    G4String      m_thePostPVname;
    G4double      m_PreZ;
    G4double      m_PostZ;
    G4double edep;
    G4double energy;     // the energy of the photon hitting PMT
    G4ThreeVector pos;   // the photon hitting position
    G4ThreeVector mom;   // the photon momentum when hitting PMT
    G4ThreeVector pol;   // the photon polarization when hitting PMT
//    G4int pmtID;         // the ID of the PMT the photon hits
//    G4double time;       // the time when photon hitting PMT
    G4int iHitCount;     // the Number of p.e. from a Hit
    G4int weight;        // if the number of photons is scaled to speed up
                         // processing, this is the weight of the photon. 
                         // (ie for a muon, might track fewer photons and
                         // compensate with increased weight.)
    G4double wavelength; // wavelength of photon hitting PMT

//    G4int producerID;    // track ID of parent of optical photon.  See
//                         // dywTrajectory.
//    G4bool isFromCerenkov;
//    G4bool isReemission;
//
//    G4bool isOriginalOP;
//    G4double OriginalOPStartT;
//
//    G4ThreeVector boundary_pos;

  public:
    void SetEdep(G4double de) { edep = de; };
    G4double GetEdep() { return edep; };
      
  public: 
        void SetPrePVname (G4String thePrePVname) {m_thePrePVname = thePrePVname; };
        void SetPostPVname (G4String thePostPVname) {m_thePostPVname = thePostPVname; };
        void SetPreZcoordinate(G4double PreZ) { m_PreZ = PreZ; };
        void SetPostZcoordinate(G4double PostZ) {m_PostZ = PostZ;};
        G4String GetPrePVname() {return m_thePrePVname;};
        G4String GetPostPVname() {return m_thePostPVname;};
        G4double GetPreZcoordinate() { return m_PreZ; };
        G4double GetPostZcoordinate() { return m_PostZ; };
//    inline G4double GetTime() const { return time; }
//    inline void SetTime(G4double val) { time = val; }

    G4double GetWavelength() const { return wavelength; }
    void SetWavelength(G4double val) { wavelength = val; }

    G4int GetWeight() const { return weight; }
    void SetWeight(G4int val) { weight = val; }
    
    G4double GetKineticEnergy() const { return energy; }
    void SetKineticEnergy(G4double val) { energy = val; }
    
    void SetPosition(G4ThreeVector xyz) { pos = xyz; }
    G4ThreeVector GetPosition() const { return pos; }
    
    void SetMomentum(G4ThreeVector mxyz) { mom = mxyz; }
    G4ThreeVector GetMomentum() const { return mom; }
    
    void SetPolarization(G4ThreeVector pxyz) { pol=pxyz; }
    G4ThreeVector GetPolarization() const {return pol;}   
    
//    inline void SetPMTID(G4int z) { pmtID = z; }
//    inline G4int GetPMTID() const { return pmtID; }

    void SetCount(G4int n) { iHitCount = n; }
    G4int GetCount() const { return iHitCount; }

//    G4int GetProducerID() const { return producerID; }
//    void SetProducerID(G4int pid) { producerID = pid; }
//
//    G4bool IsFromCerenkov() const { return isFromCerenkov;}
//    void SetFromCerenkov(G4bool flag) { isFromCerenkov = flag; }
//
//    G4bool IsReemission() const { return isReemission;}
//    void SetReemission(G4bool flag) { isReemission = flag; }
//
//    G4bool IsOriginalOP() const { return isOriginalOP; }
//    void SetOriginalOP(G4bool flag) { isOriginalOP = flag; }
//
//    G4double GetOriginalOPStartT() { return OriginalOPStartT; }
//    void SetOriginalOPStartT(double t) { OriginalOPStartT = t; }
//
//    inline void SetBoundaryPosition(G4ThreeVector xyz) { boundary_pos = xyz; }
//    inline G4ThreeVector GetBoundaryPosition() const { return boundary_pos; }

};

// dywHit_PMT_Collection is a vector of hits
typedef G4THitsCollection<dywHit_PMT> dywHit_PMT_Collection;

extern G4Allocator<dywHit_PMT> dywHit_PMT_Allocator;

inline void* dywHit_PMT::operator new(size_t)
{
  void *aHit;
  aHit = (void *) dywHit_PMT_Allocator.MallocSingle();
  return aHit;
}

inline void dywHit_PMT::operator delete(void *aHit)
{
  dywHit_PMT_Allocator.FreeSingle((dywHit_PMT*) aHit);
}

#endif
