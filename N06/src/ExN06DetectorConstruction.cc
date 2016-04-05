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
#include "G4Ellipsoid.hh"
#include "G4SubtractionSolid.hh"
#include "R12860_PMTSolid.hh"
#include "G4IntersectionSolid.hh"
#include "G4VSolid.hh"
#include "DetectorConstructionMaterial.hh"
#include "dywSD_PMT_v2.hh"
#include "G4MaterialPropertyVector.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"
//#include "G4GeometryManager.hh"
//#include "G4PhysicalVolumeStore.hh"
//#include "G4LogicalVolumeStore.hh"
//#include "G4SolidStore.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06DetectorConstruction::ExN06DetectorConstruction() 
: m_angle(0.*deg), m_WithLC(true)
{
    messenger = new ExN06DetectorConstructionMessenger(this);

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

    RotatePMT = new G4RotationMatrix();
    RotatePMT->rotateX(m_angle);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06DetectorConstruction::~ExN06DetectorConstruction(){
    delete RotatePMT;
    if (m_pmtsolid_maker) {
        delete m_pmtsolid_maker;
    }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//G4VPhysicalVolume* ExN06DetectorConstruction::Construct() {
//    // Get all Material through MatTable
//    MatTable = new DetectorConstructionMaterial();
////    ConstructDetector();
////    SetSensitiveDet();
//    return ConstructDetector();
//}

G4VPhysicalVolume* ExN06DetectorConstruction::Construct() {
    // Get all Material through MatTable
    MatTable = new DetectorConstructionMaterial();

    G4double DeployPos = 0.5*m - expHall_z;
    G4ThreeVector Translation(0.,0.,DeployPos);  // deployed position of PMT
    G4bool fCheckOverlaps = true;                // whether check overlap of physical volumes

// there may exist little overlap between pmt/body/inner1/inner2, that's ok, cause they are in one log volume

//	------------- Volumes --------------
// The experimental Hall
//
    G4Box* expHall_box = new G4Box("World",expHall_x,expHall_y,expHall_z);
  
    G4LogicalVolume* expHall_log
                    = new G4LogicalVolume(expHall_box,MatTable->GetAir(),"World",0,0,0);
  
    G4VPhysicalVolume* expHall_phys
                    = new G4PVPlacement(0,G4ThreeVector(),expHall_log,"World",0,false,0);

// light concentrator
//
    if (m_WithLC) {
        double gap = 0.1*mm;
        double rLPmt_2 = m_pmt_r + gap;
        double hhPmt_2 = m_z_equator + gap;
        double dR_LPmtBlinder = 0.05*cm;
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
        //G4double ConcentratorPos_x = 0.*cm;
        //G4double ConcentratorPos_y = 0.*cm;
        //G4double ConcentratorPos_z = 0.*cm;
        G4VPhysicalVolume* Concentrator_phys = new G4PVPlacement(
                RotatePMT,
                Translation, //G4ThreeVector(ConcentratorPos_x,ConcentratorPos_y,ConcentratorPos_z),
                Concentrator_log,     // its logical volume
                "Conc_phys",          // its name
                expHall_log,          // its mother volume
                false,                // no boolean os
                0,                    // no particular field
                fCheckOverlaps);

        // visualization of LC
        G4VisAttributes* ConeVisAtt = new G4VisAttributes(G4Colour(0.75,0.75,0.75));
        ConeVisAtt->SetForceWireframe(true);
        ConeVisAtt->SetForceAuxEdgeVisible(true);
        Concentrator_log->SetVisAttributes(ConeVisAtt);

        // surface of light concentrator
        const G4int NUM = 2;
        G4double PP[NUM] = { 1.4E-9*GeV,6.2E-9*GeV};     
        G4double REFLECTIVITY_conc[NUM] = {1.,1. }; 
        G4cout<<"ref conc = "<<REFLECTIVITY_conc[0]<<G4endl; 
        G4double EFFICIENCY_conc[NUM] = { 0., 0.};
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

        // pmt blinder
        G4Ellipsoid* solidLPmtBlinder_part1 = new G4Ellipsoid("LPmtBlinder_part1",
                                                            rLPmt_2,
                                                            rLPmt_2,
                                                            hhPmt_2,
                                                            0,
                                                            Zconc[0]-gap);

        G4Ellipsoid* solidLPmtBlinder_part2 = new G4Ellipsoid("LPmtBlinder_part2",
                                                            rLPmt_2 + dR_LPmtBlinder,
                                                            rLPmt_2 + dR_LPmtBlinder,
                                                            hhPmt_2 + dR_LPmtBlinder,
                                                            0,
                                                            Zconc[0]-gap);
        G4SubtractionSolid* solidLPmtBlinder = new G4SubtractionSolid("LPmtBlinder",
                                                            solidLPmtBlinder_part2,
                                                            solidLPmtBlinder_part1,
                                                            0,
                                                            G4ThreeVector());
        G4cout << "Exit aperture is not at the equator" << G4endl;

        G4LogicalVolume* logicLPmtBlinder = new G4LogicalVolume(solidLPmtBlinder,MatTable->GetBlacksheet(),"LPmtBlinder");
        G4VPhysicalVolume* physicLPmtBlinder = new G4PVPlacement(
                RotatePMT,
                Translation, //G4ThreeVector(0,0,0),
                logicLPmtBlinder,     // its logical volume
                "LPmtBlinder",	  // its name
                expHall_log,	  // its mother volume
                false,			  // no boolean os
                0,		  // no particular field
                fCheckOverlaps);

        // visualization of pmt blinder
        G4VisAttributes* bottom_visatt = new G4VisAttributes(G4Colour(0, 0, 0));
        bottom_visatt -> SetForceWireframe(true);  
        bottom_visatt -> SetForceAuxEdgeVisible(true);
        logicLPmtBlinder -> SetVisAttributes(bottom_visatt);

        // surface of pmt blinder
        G4OpticalSurface *OpSurfPMTBlinder;		
        OpSurfPMTBlinder = new G4OpticalSurface("OpSurfPMTBlinder");
        OpSurfPMTBlinder->SetType(dielectric_dielectric);
        OpSurfPMTBlinder->SetModel(unified);
        OpSurfPMTBlinder->SetFinish(groundfrontpainted);
        OpSurfPMTBlinder->SetSigmaAlpha(0.1);
        G4MaterialPropertiesTable* PMTBlinderMPT = new G4MaterialPropertiesTable();
        PMTBlinderMPT = MatTable->GetBlacksheet()->GetMaterialPropertiesTable();
        OpSurfPMTBlinder->SetMaterialPropertiesTable(PMTBlinderMPT);
    }

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
    m_phys_pmt = new G4PVPlacement(RotatePMT,
                                   Translation, //noTranslation,
                                   m_logical_pmt,
                                   "pmt_phys",
                                   expHall_log,
                                   false,
                                   0,
                                   fCheckOverlaps);
    body_phys = new G4PVPlacement(0, //RotatePMT,
                                  noTranslation,
                                  body_log,
                                  "body_phys",
                                  m_logical_pmt,
                                  false,
                                  0,
                                  fCheckOverlaps);
    // place inner solids in outer solid
    inner1_phys = new G4PVPlacement(0, //RotatePMT,
                                    noTranslation,
                                    inner1_log,
                                    "inner1_phys",
                                    body_log,
                                    false,
                                    0,
                                    fCheckOverlaps);
    inner2_phys = new G4PVPlacement(0, //RotatePMT,
                                    noTranslation,
                                    inner2_log,
                                    "inner2_phys",
                                    body_log,
                                    false,
                                    0,
                                    fCheckOverlaps);
                                    
    // visualization of PMT
    G4VisAttributes* visAtt;
    visAtt = new G4VisAttributes(G4Color(0.7,0.5,0.3));
    visAtt->SetForceSolid(true);
    inner1_log->SetVisAttributes(visAtt);
    visAtt = new G4VisAttributes(G4Color(0.6,0.7,0.8));
    visAtt->SetForceSolid(true);
    inner2_log->SetVisAttributes(visAtt);
    
    // surface of Photocathode
    G4OpticalSurface* OpSurfPhotocathode = new G4OpticalSurface("OpSurfPhotocathode");
    OpSurfPhotocathode->SetType(dielectric_metal); // ignored if RINDEX defined
    OpSurfPhotocathode->SetModel(glisur);
    OpSurfPhotocathode->SetFinish(polished);
    G4MaterialPropertiesTable* PhotocathodeMPT = new G4MaterialPropertiesTable();
    PhotocathodeMPT = MatTable->GetPhotocathode()->GetMaterialPropertiesTable();
    OpSurfPhotocathode->SetMaterialPropertiesTable(PhotocathodeMPT);

    G4LogicalBorderSurface* photocathodeSurf1 = 
        new G4LogicalBorderSurface("photocathode_logsurf1",inner1_phys,body_phys,OpSurfPhotocathode);
    G4LogicalBorderSurface* photocathodeSurf2 = 
        new G4LogicalBorderSurface("photocathode_logsurf2",body_phys,inner1_phys,OpSurfPhotocathode);

    // surface of mirror
    // construct a static mirror surface with idealized properties
    G4OpticalSurface* OpSurfMirror = new G4OpticalSurface("OpSurfMirror");
    OpSurfMirror->SetType(dielectric_metal);
    OpSurfMirror->SetFinish(polishedfrontpainted); // needed for mirror
    OpSurfMirror->SetModel(glisur);
    OpSurfMirror->SetPolish(0.999); 
    G4cout << "Warning: setting PMT mirror reflectivity to 0.9999 "
           << "because no PMT_Mirror material properties defined" << G4endl;
    G4MaterialPropertiesTable* propMirror = new G4MaterialPropertiesTable();
    propMirror->AddProperty("REFLECTIVITY", new G4MaterialPropertyVector());
    propMirror->AddEntry("REFLECTIVITY", 1.55*eV, 0.9999);
    propMirror->AddEntry("REFLECTIVITY", 15.5*eV, 0.9999);
    OpSurfMirror->SetMaterialPropertiesTable( propMirror );

    G4LogicalBorderSurface* mirrorSurf1 = 
        new G4LogicalBorderSurface("mirror_logsurf1",inner2_phys,body_phys,OpSurfMirror);
    G4LogicalBorderSurface* mirrorSurf2 = 
        new G4LogicalBorderSurface("mirror_logsurf2",body_phys,inner2_phys,OpSurfMirror);

    //	------------- Sensitive Detectors --------------
    G4SDManager* SDman = G4SDManager::GetSDMpointer();

    dywSD_PMT_v2* PMTSD = new dywSD_PMT_v2("PMTSD");
    SDman->AddNewDetector( PMTSD );
    body_log->SetSensitiveDetector( PMTSD );

    //--------- Visualization attributes -------------------------------
    G4VisAttributes* BoxVisAtt= new G4VisAttributes(G4Colour(1.0,1.0,1.0));
    expHall_log->SetVisAttributes(BoxVisAtt);  
  
    //always return the physical World
    return expHall_phys;
}

//void ExN06DetectorConstruction::SetSensitiveDet() {
//    //	------------- Sensitive Detectors --------------
//    G4SDManager* SDman = G4SDManager::GetSDMpointer();
//
//    dywSD_PMT_v2* PMTSD = new dywSD_PMT_v2("PMTSD");
//    SDman->AddNewDetector( PMTSD );
//    body_log->SetSensitiveDetector( PMTSD );
//}

void ExN06DetectorConstruction::setAngle(G4double angle) {
    if (!m_phys_pmt) {
        G4cerr << "Detector has not yet been constructed" << G4endl;
        return;
    }

    m_angle = angle;
    *RotatePMT = G4RotationMatrix();
    RotatePMT->rotateX(m_angle);
    G4RunManager::GetRunManager()->GeometryHasBeenModified();
}

//void ExN06DetectorConstruction::setLC(G4bool WithLC) {
//    m_WithLC = WithLC;
//}

//void ExN06DetectorConstruction::UpdateGeometry() {
//    // clean-up previous geometry
//    G4GeometryManager::GetInstance()->OpenGeometry();
//
//    G4PhysicalVolumeStore::GetInstance()->Clean();
//    G4LogicalVolumeStore::GetInstance()->Clean();
//    G4SolidStore::GetInstance()->Clean();
//    G4LogicalSkinSurface::CleanSurfaceTable();
//    G4LogicalBorderSurface::CleanSurfaceTable();
//    G4SurfaceProperty::CleanSurfacePropertyTable();
//
//    // define new one
//    G4RunManager::GetRunManager()->DefineWorldVolume(Construct());
//    G4RunManager::GetRunManager()->GeometryHasBeenModified();
//}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
