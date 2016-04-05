#ifndef RayExpRun_h
#define RayExpRun_h 1

#include "globals.hh"
#include "G4Run.hh"
#include "G4THitsCollection.hh"

class G4Event;

class RayExpRun : public G4Run {
    public:
        RayExpRun();
        virtual ~RayExpRun();
        G4int GetExitNum(); 

    public:
        virtual void RecordEvent(const G4Event*);

    private:
//        G4int EnterHCID;
        G4int ExitHCID;
//        G4int EnterNum;
//        G4int HitEnterNum;
        G4int HitExitNum;
        G4int ExitNum;

};

#endif
