#include "DetectorConstructionMaterial.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "OpticalProperty.icc"
#include "G4NistManager.hh"

DetectorConstructionMaterial::DetectorConstructionMaterial(){
    Si = new G4Element("Silicon", "Si", 14., 28.09*g/mole);
    O = new G4Element("Oxygen", "O", 8., 16.00*g/mole); 
    B = new  G4Element("Boron",  "B", 5, 10.811*g/mole);
    Na = new G4Element("Sodium", "Na", 11., 22.98977*g/mole);
    N = new G4Element("Nitrogen", "N", 7 , 14.01*g/mole);
}

DetectorConstructionMaterial::~DetectorConstructionMaterial() {;}

    // ---------------- Mateiral definition -------------------

G4Material*
DetectorConstructionMaterial::GetAl() {
    // Al
    G4double z, a, density;
    Al = new G4Material("Aluminum", z= 13., a= 26.98*g/mole, density= 2.7*g/cm3);

    return Al;
}

G4Material*
DetectorConstructionMaterial::GetBlacksheet() {
    // Blacksheet
    G4double density;
    G4NistManager* nistManager = G4NistManager::Instance();
    Blacksheet = new G4Material("Blacksheet",density=0.95*g/cm3,2);
    Blacksheet->AddElement(nistManager->FindOrBuildElement(6), 1);
    Blacksheet->AddElement(nistManager->FindOrBuildElement(1), 2);

    return Blacksheet;
}

G4Material*
DetectorConstructionMaterial::GetAir() {
    // Air
    G4double density;
    G4int nelements;
    Air = new G4Material("Air", density = 1.29*mg/cm3, nelements=2);
    Air->AddElement(N, 70.*perCent);
    Air->AddElement(O, 30.*perCent);
    G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
    myMPT2->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex2, 32);
    Air->SetMaterialPropertiesTable(myMPT2);

    return Air;
}

G4Material* 
DetectorConstructionMaterial::GetPyrex() {
    // Borosilicate glass (Pyrex) for the pmt faces
    // from PDG: 
    // 80% SiO2 + 13% B2O2 + 7% Na2O
    // by fractional mass?
    SiO2 = new G4Material("SiO2", 2.23*g/cm3, 2);
    SiO2->AddElement(Si, 1);
    SiO2->AddElement(O , 2);

    B2O2 = new G4Material("B2O2", 2.23*g/cm3, 2);
    B2O2->AddElement(B,  2);
    B2O2->AddElement(O,  2);

    G4Material* Na2O = new G4Material("Na2O", 2.23*g/cm3, 2);
    Na2O->AddElement(Na, 2);
    Na2O->AddElement(O,  1);
    
    Pyrex = new G4Material("Pyrex", 2.23*g/cm3, 3);
    Pyrex->AddMaterial(SiO2, .80);
    Pyrex->AddMaterial(B2O2, .13);
    Pyrex->AddMaterial(Na2O, .07);
    G4MaterialPropertiesTable* PyrexMPT = new G4MaterialPropertiesTable();
    PyrexMPT->AddProperty("RINDEX", fPP_Pyrex, fPyrexRINDEX, 6);
    PyrexMPT->AddProperty("ABSLENGTH", fPP_PyrexABS, fPyrexABSORPTION, 9);
    Pyrex->SetMaterialPropertiesTable(PyrexMPT);

    return Pyrex;
}

G4Material*
DetectorConstructionMaterial::GetVaccum() {
    // Vacuum
    // --- PMT vacuum is very dilute air -------
    G4double density     =  1e-3 * kGasThreshold;         //from PhysicalConstants.h
    G4double temperature = STP_Temperature;         //from PhysicalConstants.h
    G4double pressure    = STP_Pressure * density / (1.29e-3*g/cm3);
    Vacuum = new G4Material("Vacuum", density, 1, kStateGas,temperature,pressure);
    Vacuum->AddMaterial(Air, 1.);

    G4double VacPP[4] =
    {
        1.55*eV, 6.20*eV, 10.33*eV, 15.5*eV
    };
    G4double VacRINDEX[4] =
    {
        1.000001, 1.000001, 1.000001, 1.000001
    };
    G4double VacABSLENGTH[4] =
    {
        1.0e6*m, 1.0e6*m, 1.0e6*m, 1.0e6*m
    };
    G4MaterialPropertiesTable* VacMPT = new G4MaterialPropertiesTable();
    VacMPT->AddProperty("RINDEX",    VacPP, VacRINDEX,     4);
    VacMPT->AddProperty("ABSLENGTH", VacPP, VacABSLENGTH,  4);
    Vacuum->SetMaterialPropertiesTable(VacMPT);

    return Vacuum;
}







