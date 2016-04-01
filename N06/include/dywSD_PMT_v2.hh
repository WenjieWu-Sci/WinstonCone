//--------------------------------------------------------------------------
//                            dywSD_PMT_v2
//
// PMTs are difined as sensitive detector. They collect hits on them.
// The data members of hits are set up here using the information of G4Step.
// -------------------------------------------------------------------------
// Author: Liang Zhan, 2006/01/27
// Modified by: Weili Zhong, 2006/03/01
// -------------------------------------------------------------------------

#ifndef dywSD_PMT_v2_h
#define dywSD_PMT_v2_h 1

#include "G4VSensitiveDetector.hh"
#include "G4ThreeVector.hh"
#include "dywHit_PMT.hh"
//#include "dywHit_PMT_muon.hh"
//#include "IToolForSD_PMT.h"
//#include "PMTHitMerger.hh"
#include <map>
#include <vector>
#include "G4Step.hh"
//#include <TF1.h>

//////////////////////////////////////////////////////////////////////////

class G4Step;
class G4Track;
class G4HCofThisEvent;

class dywSD_PMT_v2 : public G4VSensitiveDetector //, public IToolForSD_PMT
{
    public:
        //dywSD_PMT_v2(const std::string& name);
        dywSD_PMT_v2(G4String name);
        ~dywSD_PMT_v2();

        void Initialize(G4HCofThisEvent*HCE);
        G4bool ProcessHits(G4Step*aStep,G4TouchableHistory*ROhist);
        void EndOfEvent(G4HCofThisEvent*HCE);

//        void clear();
//        void DrawAll();
//        void PrintAll();
//        void SimpleHit( const ParamsForSD_PMT& );

        void setCEMode(const std::string& mode);
        void setCEFlatValue(double v) {m_ce_flat_value = v;}
//        void setMergeFlag(bool f) { m_merge_flag = f; }
//        void setMergeWindows(double t) { m_time_window = t; }
//        void setMerger(PMTHitMerger* phm) { m_pmthitmerger=phm; }

        void setHitType(int i) { m_hit_type = i; }
        int getHitType() { return m_hit_type; }

        void disableSD() { m_disable = true; }
        void enableSD() { m_disable = false; }

        //void setCEFunc(const std::string& func, const std::vector<double>& param);
    private:
//        int get_pmtid(G4Track*);
        double get_ce(const std::string& volname, const G4ThreeVector& localpos);
    private:
        dywHit_PMT_Collection* hitCollection;
    private:
        bool m_debug;
        std::string m_ce_mode;

        // if flat mode enabled, this is used to set the fixed number
        double m_ce_flat_value;

        // 20inchfunc mode, function mode
        std::string m_ce_func_str;
        std::vector<double> m_ce_func_params;
       // TF1* m_ce_func;

        // flag to enable/disable the sensitive detector
        // disable is true, means disable the SD.
        bool m_disable;

    private:
        // ========================================================================
        // merge related
        // ========================================================================
        // = merge flag 
 //       bool m_merge_flag;
 //       // = the time window is used when merge is enabled.
 //       double m_time_window;

        typedef std::multimap<int, int> PMTID2COLIDS;
        typedef std::pair< PMTID2COLIDS::iterator, PMTID2COLIDS::iterator > PMTITER;
        PMTID2COLIDS m_pmtid2idincol;

        // new merger
//        PMTHitMerger* m_pmthitmerger;

        // ========================================================================
        // Hit Collection switcher
        // * normal, dywHit_PMT
        // * for muon, dywHit_PMT_muon
        // ========================================================================
        G4int m_hit_type;
};

#endif

