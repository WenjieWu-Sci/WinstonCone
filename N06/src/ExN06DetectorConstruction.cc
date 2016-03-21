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
#include "RayExpPMTSDv2.hh"
#include "G4Ellipsoid.hh"
#include "G4SubtractionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06DetectorConstruction::ExN06DetectorConstruction()
{
  expHall_x = expHall_y = expHall_z = 3.0*m;
//  tank_x    = tank_y    = tank_z    =  5.0*m;
//  bubble_x  = bubble_y  = bubble_z  =  0.5*m;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06DetectorConstruction::~ExN06DetectorConstruction(){;}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* ExN06DetectorConstruction::Construct()
{

//	------------- Materials -------------

  G4double a, z, density;
  G4int nelements;

// Aluminum
  G4Material* Al = new G4Material("Aluminum", z= 13., a= 26.98*g/mole, density= 2.7*g/cm3);

// Blacksheet
    G4NistManager* nistManager = G4NistManager::Instance();
    G4Material* Blacksheet = new G4Material("Blacksheet",density=0.95*g/cm3,2);
    Blacksheet->AddElement(nistManager->FindOrBuildElement(6), 1);
    Blacksheet->AddElement(nistManager->FindOrBuildElement(1), 2);

// Air
// 
  G4Element* N = new G4Element("Nitrogen", "N", z=7 , a=14.01*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=16.00*g/mole);
  G4Material* Air = new G4Material("Air", density=1.29*mg/cm3, nelements=2);
  Air->AddElement(N, 70.*perCent);
  Air->AddElement(O, 30.*perCent);

  const G4int nEntries = 32;
  G4double PhotonEnergy[nEntries] =
            { 2.034*eV, 2.068*eV, 2.103*eV, 2.139*eV,
              2.177*eV, 2.216*eV, 2.256*eV, 2.298*eV,
              2.341*eV, 2.386*eV, 2.433*eV, 2.481*eV,
              2.532*eV, 2.585*eV, 2.640*eV, 2.697*eV,
              2.757*eV, 2.820*eV, 2.885*eV, 2.954*eV,
              3.026*eV, 3.102*eV, 3.181*eV, 3.265*eV,
              3.353*eV, 3.446*eV, 3.545*eV, 3.649*eV,
              3.760*eV, 3.877*eV, 4.002*eV, 4.136*eV };
  G4double RefractiveIndex2[nEntries] =
            { 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
              1.00, 1.00, 1.00, 1.00 };
  G4MaterialPropertiesTable* myMPT2 = new G4MaterialPropertiesTable();
  myMPT2->AddProperty("RINDEX", PhotonEnergy, RefractiveIndex2, nEntries);
  Air->SetMaterialPropertiesTable(myMPT2);

//
//	------------- Volumes --------------

// The experimental Hall
//
  G4Box* expHall_box = new G4Box("World",expHall_x,expHall_y,expHall_z);

  G4LogicalVolume* expHall_log
    = new G4LogicalVolume(expHall_box,Air,"World",0,0,0);

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
  G4LogicalVolume* Concentrator_log = new G4LogicalVolume(Concentrator_solid,Al,"Conc_log",0,0,0);
  G4double ConcentratorPos_x = 0.*cm;
  G4double ConcentratorPos_y = 0.*cm;
  G4double ConcentratorPos_z = 0.*cm;
  G4VPhysicalVolume* Concentrator_phys = new G4PVPlacement(0,
          G4ThreeVector(ConcentratorPos_x,ConcentratorPos_y,ConcentratorPos_z),
          Concentrator_log,     // its logical volume
          "Conc_phys",          // its name
          expHall_log, // its mother volume
          false,                // no boolean os
          0,                    // no particular field
          fCheckOverlaps);

// Entrance aperture

  G4double innerRadius = 0.*cm;
  G4double outerRadius = (radius[NumOfZ-1] - 0.05)*cm;
  G4double hight = (height[NumOfZ-1] - height[NumOfZ-2])/2*cm;
  G4double startAngle = 0.*deg;
  G4double spanningAngle = 360.*deg;
  G4Tubs* entrance_tube = new G4Tubs("entra_tube",innerRadius,
                                    outerRadius,hight,
                                    startAngle,spanningAngle);
  G4LogicalVolume* entrance_log = new G4LogicalVolume(entrance_tube,Air,"entra_log",0,0,0);
  G4double Pos_x = 0.*cm;
  G4double Pos_y = 0.*cm;
  G4double Pos_z = height[NumOfZ-1]*cm;
  G4VPhysicalVolume* entrance_phys = new G4PVPlacement(0,
             G4ThreeVector(Pos_x,Pos_y,Pos_z),
             entrance_log,
             "entra_phys",
             expHall_log,
             false,
             0,
             fCheckOverlaps);

// photon collector: exit aperture

//  G4double innerRadiusOfTheTube = 0.*cm;
//  G4double outerRadiusOfTheTube = (radius[0] - 0.05)*cm;
//  G4double hightOfTheTube = (height[1] - height[0])/2*cm;
//  G4double startAngleOfTheTube = 0.*deg;
//  G4double spanningAngleOfTheTube = 360.*deg;
//  G4Tubs* collector_tube = new G4Tubs("collec_tube",innerRadiusOfTheTube,
//                                    outerRadiusOfTheTube,hightOfTheTube,
//                                    startAngleOfTheTube,spanningAngleOfTheTube);

    double rLPmt = 25.4*cm;
    double gap = 0.5*mm;
    double rLPmt_2 = rLPmt - gap; 
    double hhPmt = 18.4*cm;
    double hhPmt_2 = hhPmt - gap;
    double BottomCut = Zconc[0] - 0.1*cm; 

    double gap2 = 1.01*mm;
    double gap3 = 2.*mm;
    double rLPmt_3 = rLPmt - gap2;
    double rLPmt_4 = rLPmt - gap3;
    double hhPmt_3 = hhPmt - gap2;
    double hhPmt_4 = hhPmt - gap3;
    G4Ellipsoid* collector_part1 = new G4Ellipsoid("collector_part1",
                                                        rLPmt_2,
                                                        rLPmt_2,
                                                        hhPmt_2,
                                                        BottomCut,
                                                        hhPmt_2);
    G4Ellipsoid* collector_part2 = new G4Ellipsoid("collector_part2",
                                                        rLPmt,
                                                        rLPmt,
                                                        hhPmt,
                                                        BottomCut,
                                                        hhPmt);
    G4SubtractionSolid* collector_tube = new G4SubtractionSolid("collec_tube",
                                                        collector_part2,
                                                        collector_part1,
                                                        0,
                                                        G4ThreeVector());
    G4LogicalVolume* collector_log = new G4LogicalVolume(collector_tube,Air,"collec_log",0,0,0);
    G4double collectorPos_x = 0.*cm;
    G4double collectorPos_y = 0.*cm;
    G4double collectorPos_z = 0.*cm;
    G4VPhysicalVolume* collector_phys = new G4PVPlacement(0,
               G4ThreeVector(collectorPos_x,collectorPos_y,collectorPos_z),
               collector_log,
               "collec_phys",
               expHall_log,
               false,
               0,
               fCheckOverlaps);

    G4Ellipsoid* absorber_part1 = new G4Ellipsoid("absorb_part1",
                                                        rLPmt_3,
                                                        rLPmt_3,
                                                        hhPmt_3,
                                                        BottomCut,
                                                        hhPmt_3);
    G4Ellipsoid* absorber_part2 = new G4Ellipsoid("absorb_part2",
                                                        rLPmt_4,
                                                        rLPmt_4,
                                                        hhPmt_4,
                                                        BottomCut,
                                                        hhPmt_4);
    G4SubtractionSolid* absorber_ellip = new G4SubtractionSolid("absorb_ellip",
                                                        absorber_part1,
                                                        absorber_part2,
                                                        0,
                                                        G4ThreeVector());
    G4LogicalVolume* absorber_log = new G4LogicalVolume(absorber_ellip,Blacksheet,"absorb_log",0,0,0);
    G4double absorberPos_x = 0.*cm;
    G4double absorberPos_y = 0.*cm;
    G4double absorberPos_z = 0.*cm;
    G4VPhysicalVolume* absorber_phys = new G4PVPlacement(0,
               G4ThreeVector(absorberPos_x,absorberPos_y,absorberPos_z),
               absorber_log,
               "absorb_phys",
               expHall_log,
               false,
               0,
               fCheckOverlaps);

//// The Water Tank
////	
//  G4Box* waterTank_box = new G4Box("Tank",tank_x,tank_y,tank_z);
//
//  G4LogicalVolume* waterTank_log
//    = new G4LogicalVolume(waterTank_box,Water,"Tank",0,0,0);
//
//  G4VPhysicalVolume* waterTank_phys
//    = new G4PVPlacement(0,G4ThreeVector(),waterTank_log,"Tank",
//                        expHall_log,false,0);
//
//// The Air Bubble
////   
//  G4Box* bubbleAir_box = new G4Box("Bubble",bubble_x,bubble_y,bubble_z);
//
//  G4LogicalVolume* bubbleAir_log
//    = new G4LogicalVolume(bubbleAir_box,Air,"Bubble",0,0,0);
//
////G4VPhysicalVolume* bubbleAir_phys =
//      new G4PVPlacement(0,G4ThreeVector(0,2.5*m,0),bubbleAir_log,"Bubble",
//                        waterTank_log,false,0);
//
//	------------- Surfaces --------------
//
// surface of light concentrator
//
    const G4int NUM = 2;
    G4double PP[NUM] = { 1.4E-9*GeV,6.2E-9*GeV};     
    G4double REFLECTIVITY_conc[NUM] = {0.9,0.9}; 
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

    const G4int NUMENTRIES_water=60;
    G4double ENERGY_water[NUMENTRIES_water] =
    { 1.56962e-09*GeV, 1.58974e-09*GeV, 1.61039e-09*GeV, 1.63157e-09*GeV, 
        1.65333e-09*GeV, 1.67567e-09*GeV, 1.69863e-09*GeV, 1.72222e-09*GeV, 
        1.74647e-09*GeV, 1.77142e-09*GeV, 1.7971e-09*GeV,  1.82352e-09*GeV, 
        1.85074e-09*GeV, 1.87878e-09*GeV, 1.90769e-09*GeV, 1.93749e-09*GeV, 
        1.96825e-09*GeV, 1.99999e-09*GeV, 2.03278e-09*GeV, 2.06666e-09*GeV,
        2.10169e-09*GeV, 2.13793e-09*GeV, 2.17543e-09*GeV, 2.21428e-09*GeV, 
        2.25454e-09*GeV, 2.29629e-09*GeV, 2.33962e-09*GeV, 2.38461e-09*GeV, 
        2.43137e-09*GeV, 2.47999e-09*GeV, 2.53061e-09*GeV, 2.58333e-09*GeV, 
        2.63829e-09*GeV, 2.69565e-09*GeV, 2.75555e-09*GeV, 2.81817e-09*GeV, 
        2.88371e-09*GeV, 2.95237e-09*GeV, 3.02438e-09*GeV, 3.09999e-09*GeV,
        3.17948e-09*GeV, 3.26315e-09*GeV, 3.35134e-09*GeV, 3.44444e-09*GeV, 
        3.54285e-09*GeV, 3.64705e-09*GeV, 3.75757e-09*GeV, 3.87499e-09*GeV, 
        3.99999e-09*GeV, 4.13332e-09*GeV, 4.27585e-09*GeV, 4.42856e-09*GeV, 
        4.59258e-09*GeV, 4.76922e-09*GeV, 4.95999e-09*GeV, 5.16665e-09*GeV, 
        5.39129e-09*GeV, 5.63635e-09*GeV, 5.90475e-09*GeV, 6.19998e-09*GeV };
    G4double RINDEX_blacksheet[NUM] = { 1.6, 1.6 };
    G4double SPECULARLOBECONSTANT[NUM] = { 0.3, 0.3 };
    G4double SPECULARSPIKECONSTANT[NUM] ={ 0.2, 0.2 };
    G4double BACKSCATTERCONSTANT[NUM] =  { 0.2, 0.2 };
    G4double EFFICIENCY_blacksheet[NUMENTRIES_water] = { 0.0 };

    G4double SK1SK2FF = 1.0;
    G4bool BlackSheetFudgeFactor=true;
    if (BlackSheetFudgeFactor) SK1SK2FF=SK1SK2FF*1.55;

    G4double REFLECTIVITY_blacksheet[NUMENTRIES_water] = {0.0};
    //{ 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    //    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    //    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    //    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    //    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    //    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    //    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    //    0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 0.055*SK1SK2FF, 
    //    0.055*SK1SK2FF, 0.057*SK1SK2FF, 0.059*SK1SK2FF, 0.060*SK1SK2FF, 
    //    0.059*SK1SK2FF, 0.058*SK1SK2FF, 0.057*SK1SK2FF, 0.055*SK1SK2FF, 
    //    0.050*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 
    //    0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF,
    //    0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF,
    //    0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF,
    //    0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF, 0.045*SK1SK2FF ,
    //    0.045*SK1SK2FF, 0.045*SK1SK2FF };

    G4OpticalSurface *OpSurfBlinder = new G4OpticalSurface("OpSurfBlinder");		
    OpSurfBlinder->SetType(dielectric_dielectric);
    OpSurfBlinder->SetModel(unified);
    OpSurfBlinder->SetFinish(groundfrontpainted);
    OpSurfBlinder->SetSigmaAlpha(0.1);
    G4MaterialPropertiesTable *myST1 = new G4MaterialPropertiesTable();  			   
    myST1->AddProperty("RINDEX", ENERGY_water, RINDEX_blacksheet, NUMENTRIES_water);
    myST1->AddProperty("SPECULARLOBECONSTANT", PP, SPECULARLOBECONSTANT, NUM);
    myST1->AddProperty("SPECULARSPIKECONSTANT", PP, SPECULARSPIKECONSTANT, NUM);
    myST1->AddProperty("BACKSCATTERCONSTANT", PP, BACKSCATTERCONSTANT, NUM);
    myST1->AddProperty("REFLECTIVITY", ENERGY_water, REFLECTIVITY_blacksheet, NUMENTRIES_water);
    myST1->AddProperty("EFFICIENCY", ENERGY_water, EFFICIENCY_blacksheet, NUMENTRIES_water);
    OpSurfBlinder->SetMaterialPropertiesTable(myST1);

    G4LogicalSkinSurface* surfBlinder;					
    surfBlinder  = new G4LogicalSkinSurface("surfBlinder",absorber_log,OpSurfBlinder);

//  OpWaterSurface->SetType(dielectric_dielectric);
//  OpWaterSurface->SetFinish(ground);
//  OpWaterSurface->SetModel(unified);
//
//  new G4LogicalBorderSurface("WaterSurface",
//                                 waterTank_phys,expHall_phys,OpWaterSurface);
//
//// Air Bubble
////
//  G4OpticalSurface* OpAirSurface = new G4OpticalSurface("AirSurface");
//  OpAirSurface->SetType(dielectric_dielectric);
//  OpAirSurface->SetFinish(polished);
//  OpAirSurface->SetModel(glisur);
//
//  G4LogicalSkinSurface* AirSurface = 
//	  new G4LogicalSkinSurface("AirSurface", bubbleAir_log, OpAirSurface);
//
//  G4OpticalSurface* opticalSurface = dynamic_cast <G4OpticalSurface*>
//        (AirSurface->GetSurface(bubbleAir_log)->GetSurfaceProperty());
//
//  if (opticalSurface) opticalSurface->DumpInfo();
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

    RayExpPMTSD* entranceSD = new RayExpPMTSD("entranceSD");
    SDman->AddNewDetector( entranceSD );
    entrance_log->SetSensitiveDetector( entranceSD );

    RayExpPMTSDv2* collectorSD = new RayExpPMTSDv2("collectorSD");
    SDman->AddNewDetector( collectorSD );
    collector_log->SetSensitiveDetector( collectorSD );

//--------- Visualization attributes -------------------------------

  G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
  expHall_log->SetVisAttributes(BoxVisAtt);  

  G4VisAttributes* VisAtt = new G4VisAttributes(G4Colour(0.5,0,0));
//  VisAtt->SetForceWireframe(true);
//  VisAtt->SetForceAuxEdgeVisible(true);
  entrance_log->SetVisAttributes(VisAtt);

  G4VisAttributes* ConeVisAtt = new G4VisAttributes(G4Colour(1.0,0.5,0));
  ConeVisAtt->SetForceWireframe(true);
  ConeVisAtt->SetForceAuxEdgeVisible(true);
  Concentrator_log->SetVisAttributes(ConeVisAtt);

  G4VisAttributes* TubeVisAtt = new G4VisAttributes(G4Colour(1.0,0,0));
//  TubeVisAtt->SetForceWireframe(true);
//  TubeVisAtt->SetForceAuxEdgeVisible(true);
  collector_log->SetVisAttributes(TubeVisAtt);
  
  G4VisAttributes* EllipVisAtt = new G4VisAttributes(G4Colour(0,0.5,0.5));
//  TubeVisAtt->SetForceWireframe(true);
//  TubeVisAtt->SetForceAuxEdgeVisible(true);
  absorber_log->SetVisAttributes(EllipVisAtt);

//always return the physical World
  return expHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
