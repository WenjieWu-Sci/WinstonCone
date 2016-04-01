//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: ExN06DetectorConstruction.cc,v 1.18 2010-10-23 19:27:38 gum Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN06DetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"
#include "G4OpticalSurface.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4Polycone.hh"
#include "LC_JUNO_coordinates_Design2_480mm.icc"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4NistManager.hh"
#include "G4SDManager.hh"
#include "RayExpPMTSD.hh"
#include "G4Ellipsoid.hh"
#include "G4SubtractionSolid.hh"
#include "R12860_PMTSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4VSolid.hh"
#include "DetectorConstructionMaterial.hh"
#include "dywSD_PMT_v2.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06DetectorConstruction::ExN06DetectorConstruction()
{
    expHall_x = expHall_y = 1.0*m;
    expHall_z = 5.0*m;
    m_pmt_r = 254.*mm;
    m_pmt_h = 523.*mm;
    m_z_equator = 184.*mm;  // From top to equator
    m_pmtsolid_maker = new R12860_PMTSolid(
                            m_pmt_r,
                            m_z_equator,
                            47*mm,
                            50*mm,
                            m_pmt_h,
                            120.*mm);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06DetectorConstruction::~ExN06DetectorConstruction(){
    if (m_pmtsolid_maker) {
        delete m_pmtsolid_maker;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExN06DetectorConstruction::Construct()
{

    DetectorConstructionMaterial* MatTable = new DetectorConstructionMaterial();
//
//	------------- Volumes --------------

    G4double DeployPos = 0.5*m - expHall_z;
    G4ThreeVector Translation(0.,0.,DeployPos);
// The experimental Hall
//
    G4Box* expHall_box = new G4Box("World",expHall_x,expHall_y,expHall_z);
  
    G4LogicalVolume* expHall_log
                    = new G4LogicalVolume(expHall_box,MatTable->GetAir(),"World",0,0,0);
  
    G4VPhysicalVolume* expHall_phys
                    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);

// light concentrator
//
    G4bool fCheckOverlaps = true;
  
    NumOfZ = 200;
    for (int i=0;i<NumOfZ;i++) {
        Zconc[i] = height[i]*cm;
        RconcMin[i] = radius[i]*cm;
        RconcMax[i] = (radius[i] + 0.1)*cm;
    //  RExtconcMin[i] = (radius[i] + 0.11)*cm;
    //  RExtconcMax[i] = (radius[i] + 0.16)*cm;
    //  G4cout << "RconcMin =  " << RconcMin[i] << G4endl;
    //  G4cout << "RconcMax =  " << RconcMax[i] << G4endl;
    //  G4cout << "RExtconcMin =  " << RExtconcMin[i] << G4endl;
    //  G4cout << "RExtconcMax =  " << RExtconcMax[i] << G4endl;
    }
    G4cout << "NumOfZ = " << NumOfZ << G4endl;
    G4Polycone* Concentrator_solid = new G4Polycone("Conc_solid",0.0*deg,360.0*deg,NumOfZ,Zconc,RconcMin,RconcMax);
    G4LogicalVolume* Concentrator_log = new G4LogicalVolume(Concentrator_solid,MatTable->GetAl(),"Conc_log",0,0,0);
    G4double ConcentratorPos_x = 0.*cm;
    G4double ConcentratorPos_y = 0.*cm;
    G4double ConcentratorPos_z = 0.*cm;
    G4VPhysicalVolume* Concentrator_phys = new G4PVPlacement(0,
            Translation, //G4ThreeVector(ConcentratorPos_x,ConcentratorPos_y,ConcentratorPos_z),
            Concentrator_log,     // its logical volume
            "Conc_phys",          // its name
            expHall_log,          // its mother volume
            false,                // no boolean os
            0,                    // no particular field
            fCheckOverlaps);

// Entrance aperture

//    G4double innerRadius = 0.*cm;
//    G4double outerRadius = (radius[NumOfZ-1] - 0.005)*cm;
//    G4double hight = (height[NumOfZ-1] - height[NumOfZ-2])/2*cm;
//    G4double startAngle = 0.*deg;
//    G4double spanningAngle = 360.*deg;
//    G4Tubs* entrance_tube = new G4Tubs("entra_tube",innerRadius,
//                                      outerRadius,hight,
//                                      startAngle,spanningAngle);
    
//    G4double space = 0.005*cm;
//    G4int NumOfZentrance = 2;
//    G4double Zentrance[2] = {height[NumOfZ-2]*cm, height[NumOfZ-1]*cm};
//    G4double RMinentrance[2] = {0, 0};
//    G4double RMaxentrance[2] = {radius[NumOfZ-2]*cm-space, radius[NumOfZ-1]*cm-space};
//    G4Polycone* entrance_tube = new G4Polycone("entra_tube",0.0*deg,360.0*deg,NumOfZentrance,Zentrance,RMinentrance,RMaxentrance);
//    
//    G4LogicalVolume* entrance_log = new G4LogicalVolume(entrance_tube,MatTable->GetAir(),"entra_log",0,0,0);
//    G4double Pos_x = 0.*cm;
//    G4double Pos_y = 0.*cm;
//    G4double Pos_z = 0.*cm;
//    G4VPhysicalVolume* entrance_phys = new G4PVPlacement(0,
//               Translation, //G4ThreeVector(Pos_x,Pos_y,Pos_z),
//               entrance_log,
//               "entra_phys",
//               expHall_log,
//               false,
//               0,
//               fCheckOverlaps);

    // pmt

    pmt_solid = m_pmtsolid_maker->GetSolid("pmt_solid",1E-3*mm);
    body_solid = m_pmtsolid_maker->GetSolid("body_solid");
    inner_solid = m_pmtsolid_maker->GetSolid("inner_solid",-5*mm);
    G4double helper_sep_tube_r = m_pmt_r;
    G4double helper_sep_tube_h = m_z_equator;
    G4double helper_sep_tube_hh = helper_sep_tube_h/2;

    G4VSolid* pInnerSep = new G4Tubs("Inner_Separator",
            0.,
            helper_sep_tube_r+1E-9*mm,
            helper_sep_tube_hh+1E-9*mm,
            0.,360.*degree);
    G4ThreeVector innerSepDispl(0.,0.,helper_sep_tube_hh-1E-9*mm);
    inner1_solid = new G4IntersectionSolid("inner1_solid",
            inner_solid, pInnerSep, NULL, innerSepDispl);
    inner2_solid = new G4SubtractionSolid("inner2_solid",
            inner_solid, pInnerSep, NULL, innerSepDispl);

    m_logical_pmt = new G4LogicalVolume(pmt_solid, MatTable->GetPyrex(), "pmt_log");
    body_log = new G4LogicalVolume(body_solid, MatTable->GetPyrex(), "body_log");
    inner1_log = new G4LogicalVolume(inner1_solid, MatTable->GetVaccum(), "inner1_log");
    inner2_log = new G4LogicalVolume(inner2_solid, MatTable->GetVaccum(), "inner2_log");

    G4ThreeVector noTranslation(0.,0.,0.);
    
    // place outer solids in envelope
    m_phys_pmt = new G4PVPlacement(0,
                                   Translation, //noTranslation,
                                   m_logical_pmt,
                                   "pmt_phys",
                                   expHall_log,
                                   false,
                                   0,
                                   fCheckOverlaps);
    body_phys = new G4PVPlacement(0,
                                  noTranslation,
                                  body_log,
                                  "body_phys",
                                  m_logical_pmt,
                                  false,
                                  0,
                                  fCheckOverlaps);
    // place inner solids in outer solid
    inner1_phys = new G4PVPlacement(0,
                                    noTranslation,
                                    inner1_log,
                                    "inner1_phys",
                                    body_log,
                                    false,
                                    0,
                                    fCheckOverlaps);
    inner2_phys = new G4PVPlacement(0,
                                    noTranslation,
                                    inner2_log,
                                    "inner2_phys",
                                    body_log,
                                    false,
                                    0,
                                    fCheckOverlaps);
                                    
//
//	------------- Surfaces --------------
//
// surface of light concentrator
//

    const G4int NUM = 2;
    G4double PP[NUM] = { 1.4E-9*GeV,6.2E-9*GeV};     
    G4double REFLECTIVITY_conc[NUM] = {1,1}; 
    G4cout<<"ref conc = "<<REFLECTIVITY_conc[0]<<G4endl; 
    G4double EFFICIENCY_conc[NUM] = { 0.,0.};
    G4OpticalSurface* OpSurfConcentrator = new G4OpticalSurface("OpSurfConcentrator");
    OpSurfConcentrator->SetType(dielectric_metal);
    OpSurfConcentrator->SetModel(glisur);
    OpSurfConcentrator->SetFinish(polished);
    G4MaterialPropertiesTable *myST3 = new G4MaterialPropertiesTable();
    myST3->AddProperty("REFLECTIVITY", PP, REFLECTIVITY_conc, NUM);
    myST3->AddProperty("EFFICIENCY", PP, EFFICIENCY_conc, NUM);
    OpSurfConcentrator->SetMaterialPropertiesTable(myST3);

    G4LogicalSkinSurface* surfConcentrator;
    surfConcentrator  = new G4LogicalSkinSurface("surfConcentrator",Concentrator_log,OpSurfConcentrator);

// surface of absorber

    //const G4int NUMENTRIES_water=60;
    //G4double ENERGY_water[NUMENTRIES_water] =
    //{ 1.56962e-09*GeV, 1.58974e-09*GeV, 1.61039e-09*GeV, 1.63157e-09*GeV, 
    //    1.65333e-09*GeV, 1.67567e-09*GeV, 1.69863e-09*GeV, 1.72222e-09*GeV, 
    //    1.74647e-09*GeV, 1.77142e-09*GeV, 1.7971e-09*GeV,  1.82352e-09*GeV, 
    //    1.85074e-09*GeV, 1.87878e-09*GeV, 1.90769e-09*GeV, 1.93749e-09*GeV, 
    //    1.96825e-09*GeV, 1.99999e-09*GeV, 2.03278e-09*GeV, 2.06666e-09*GeV,
    //    2.10169e-09*GeV, 2.13793e-09*GeV, 2.17543e-09*GeV, 2.21428e-09*GeV, 
    //    2.25454e-09*GeV, 2.29629e-09*GeV, 2.33962e-09*GeV, 2.38461e-09*GeV, 
    //    2.43137e-09*GeV, 2.47999e-09*GeV, 2.53061e-09*GeV, 2.58333e-09*GeV, 
    //    2.63829e-09*GeV, 2.69565e-09*GeV, 2.75555e-09*GeV, 2.81817e-09*GeV, 
    //    2.88371e-09*GeV, 2.95237e-09*GeV, 3.02438e-09*GeV, 3.09999e-09*GeV,
    //    3.17948e-09*GeV, 3.26315e-09*GeV, 3.35134e-09*GeV, 3.44444e-09*GeV, 
    //    3.54285e-09*GeV, 3.64705e-09*GeV, 3.75757e-09*GeV, 3.87499e-09*GeV, 
    //    3.99999e-09*GeV, 4.13332e-09*GeV, 4.27585e-09*GeV, 4.42856e-09*GeV, 
    //    4.59258e-09*GeV, 4.76922e-09*GeV, 4.95999e-09*GeV, 5.16665e-09*GeV, 
    //    5.39129e-09*GeV, 5.63635e-09*GeV, 5.90475e-09*GeV, 6.19998e-09*GeV };
    //G4double RINDEX_blacksheet[NUM] = { 1.6, 1.6 };
    //G4double SPECULARLOBECONSTANT[NUM] = { 0.3, 0.3 };
    //G4double SPECULARSPIKECONSTANT[NUM] ={ 0.2, 0.2 };
    //G4double BACKSCATTERCONSTANT[NUM] =  { 0.2, 0.2 };
    //G4double EFFICIENCY_blacksheet[NUMENTRIES_water] = { 0.0 };

    //G4double SK1SK2FF = 1.0;
    //G4bool BlackSheetFudgeFactor=true;
    //if (BlackSheetFudgeFactor) SK1SK2FF=SK1SK2FF*1.55;

    //G4double REFLECTIVITY_blacksheet[NUMENTRIES_water] = {0.0};
    ////{ 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    ////    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    ////    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    ////    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    ////    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    ////    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    ////    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    ////    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    ////    0.055*SK1SK2FF, 0.057*SK1SK2FF, 0.059*SK1SK2FF, 0.060*SK1SK2FF, 
    ////    0.059*SK1SK2FF, 0.058*SK1SK2FF, 0.057*SK1SK2FF, 0.055*SK1SK2FF, 
    ////    0.050*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 
    ////    0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF,
    ////    0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF,
    ////    0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF,
    ////    0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF ,
    ////    0.045*SK1SK2FF, 0.045*SK1SK2FF };

    //G4OpticalSurface *OpSurfBlinder = new G4OpticalSurface("OpSurfBlinder");		
    //OpSurfBlinder->SetType(dielectric_dielectric);
    //OpSurfBlinder->SetModel(unified);
    //OpSurfBlinder->SetFinish(groundfrontpainted);
    //OpSurfBlinder->SetSigmaAlpha(0.1);
    //G4MaterialPropertiesTable *myST1 = new G4MaterialPropertiesTable();  			   
    //myST1->AddProperty("RINDEX", ENERGY_water, RINDEX_blacksheet, NUMENTRIES_water);
    //myST1->AddProperty("SPECULARLOBECONSTANT", PP, SPECULARLOBECONSTANT, NUM);
    //myST1->AddProperty("SPECULARSPIKECONSTANT", PP, SPECULARSPIKECONSTANT, NUM);
    //myST1->AddProperty("BACKSCATTERCONSTANT", PP, BACKSCATTERCONSTANT, NUM);
    //myST1->AddProperty("REFLECTIVITY", ENERGY_water, REFLECTIVITY_blacksheet, NUMENTRIES_water);
    //myST1->AddProperty("EFFICIENCY", ENERGY_water, EFFICIENCY_blacksheet, NUMENTRIES_water);
    //OpSurfBlinder->SetMaterialPropertiesTable(myST1);

    //G4LogicalSkinSurface* surfBlinder;					
    //surfBlinder  = new G4LogicalSkinSurface("surfBlinder",absorber_log,OpSurfBlinder);

//  OpWaterSurface->SetType(dielectric_dielectric);
//  OpWaterSurface->SetFinish(ground);
//  OpWaterSurface->SetModel(unified);
//
//  new G4LogicalBorderSurface("WaterSurface",
//                                 waterTank_phys,expHall_phys,OpWaterSurface);
//
////
//// Generate & Add Material Properties Table attached to the optical surfaces
////
//  const G4int num = 2;
//  G4double Ephoton[num] = {2.034*eV, 4.136*eV};
//
//  //OpticalWaterSurface 
//  G4double RefractiveIndex[num] = {1.35, 1.40};
//  G4double SpecularLobe[num]    = {0.3, 0.3};
//  G4double SpecularSpike[num]   = {0.2, 0.2};
//  G4double Backscatter[num]     = {0.2, 0.2};
//
//  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();
//  
//  myST1->AddProperty("RINDEX",                Ephoton, RefractiveIndex, num);
//  myST1->AddProperty("SPECULARLOBECONSTANT",  Ephoton, SpecularLobe,    num);
//  myST1->AddProperty("SPECULARSPIKECONSTANT", Ephoton, SpecularSpike,   num);
//  myST1->AddProperty("BACKSCATTERCONSTANT",   Ephoton, Backscatter,     num);
//
//  OpWaterSurface->SetMaterialPropertiesTable(myST1);
//
//  //OpticalAirSurface
//  G4double Reflectivity[num] = {0.3, 0.5};
//  G4double Efficiency[num]   = {0.8, 1.0};
//
//  G4MaterialPropertiesTable *myST2 = new G4MaterialPropertiesTable();
//
//  myST2->AddProperty("REFLECTIVITY", Ephoton, Reflectivity, num);
//  myST2->AddProperty("EFFICIENCY",   Ephoton, Efficiency,   num);
//
//  OpAirSurface->SetMaterialPropertiesTable(myST2);

//	------------- Sensitive Detectors --------------
    
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

//    RayExpPMTSD* entranceSD = new RayExpPMTSD("entranceSD");
//    SDman->AddNewDetector( entranceSD );
//    entrance_log->SetSensitiveDetector( entranceSD );

    dywSD_PMT_v2* PMTSD = new dywSD_PMT_v2("PMTSD");
    SDman->AddNewDetector( PMTSD );
    body_log->SetSensitiveDetector( PMTSD );

//--------- Visualization attributes -------------------------------

    G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    expHall_log->SetVisAttributes(BoxVisAtt);  
  
    G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(0.5,0,0));
    VisAtt->SetForceWireframe(true);
    VisAtt->SetForceAuxEdgeVisible(true);
    entrance_log->SetVisAttributes(VisAtt);
  
    G4VisAttributes* ConeVisAtt = new G4VisAttributes(G4Colour(0.75,0.75,0.75));
    ConeVisAtt->SetForceWireframe(true);
    ConeVisAtt->SetForceAuxEdgeVisible(true);
    Concentrator_log->SetVisAttributes(ConeVisAtt);
  
    G4VisAttributes* visAtt;
    visAtt = new G4VisAttributes(G4Color(0.7,0.5,0.3));
    visAtt->SetForceSolid(true);
    inner1_log->SetVisAttributes(visAtt);
    visAtt = new G4VisAttributes(G4Color(0.6,0.7,0.8));
    visAtt->SetForceSolid(true);
    inner2_log->SetVisAttributes(visAtt);
    //visAtt = new G4VisAttributes(G4Color(0.3,0.4,0.7));
    //visAtt->SetForceSolid(true);
    //body_log->SetVisAttributes(visAtt);
    //m_logical_pmt->SetVisAttributes(visAtt);
  
//always return the physical World
    return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
