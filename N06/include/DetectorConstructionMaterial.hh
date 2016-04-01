#ifndef DetectorConstructionMaterial_h
#define DetectorConstructionMaterial_h 1

#include "globals.hh"

class G4Material;
class G4Element;

class DetectorConstructionMaterial {
    public:
        DetectorConstructionMaterial();
        ~DetectorConstructionMaterial();

    public:
        G4Material* GetAl(); 
        G4Material* GetBlacksheet(); 
        G4Material* GetAir(); 
        G4Material* GetPyrex(); 
        G4Material* GetVaccum();
        
    private:
        G4Element* Si;
        G4Element* O;
        G4Element* B;
        G4Element* Na;
        G4Element* N;

    private:
        G4Material* Al;
        G4Material* Blacksheet;
        G4Material* Air;
        G4Material* SiO2;
        G4Material* B2O2;
        G4Material* Na2O;
        G4Material* Pyrex;
        G4Material* Vacuum;
};

#endif
