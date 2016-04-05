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
// $Id: ExN06DetectorConstruction.hh,v 1.5 2006-06-29 17:53:55 gunter Exp $
// GEANT4 tag $Name: not supported by cvs2svn $
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#ifndef ExN06DetectorConstruction_h
#define ExN06DetectorConstruction_h 1

#include "globals.hh"
#include "G4VUserDetectorConstruction.hh"

#include "DetectorConstructionMaterial.hh"
#include "ExN06DetectorConstructionMessenger.hh"
#include <iostream>
#include <string>
#include "G4RotationMatrix.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class R12860_PMTSolid;
class G4VSolid;
class G4LogicalVolume;
class G4PVPlacement;

class ExN06DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
    ExN06DetectorConstruction();
   ~ExN06DetectorConstruction();

  public:
    G4VPhysicalVolume* Construct();
//    G4VPhysicalVolume* ConstructDetector();

  public:
//    void SetSensitiveDet();
    void setAngle(G4double angle);
//    void setLC(G4bool WithLC);
//    void UpdateGeometry();

  private:
    /* solid maker*/
    R12860_PMTSolid* m_pmtsolid_maker;
    // * pmt solid (a little big than body solid)
    //  * body solid
    //   + inner1
    //   + inner2
    G4VSolid* pmt_solid;
    G4VSolid* body_solid;
    G4VSolid* inner_solid;
    G4VSolid* inner1_solid;
    G4VSolid* inner2_solid;
    /* logical volumes */
    G4LogicalVolume* m_logical_pmt;
    G4LogicalVolume* body_log;
    G4LogicalVolume* inner1_log;
    G4LogicalVolume* inner2_log;
    /* physical volumes */
    G4PVPlacement* m_phys_pmt;
    G4PVPlacement* body_phys;
    G4PVPlacement* inner1_phys;
    G4PVPlacement* inner2_phys;

  private:
    DetectorConstructionMaterial* MatTable;

    G4double expHall_x;
    G4double expHall_y;
    G4double expHall_z;

    G4double m_pmt_r;
    G4double m_pmt_h;
    G4double m_z_equator;

    G4int NumOfZ;
    G4double Zconc[200];
    G4double RconcMin[200];
    G4double RconcMax[200];
    G4double RExtconcMin[200];
    G4double RExtconcMax[200];

    G4RotationMatrix* RotatePMT;
    G4double m_angle;
    G4bool m_WithLC;

    ExN06DetectorConstructionMessenger* messenger;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif /*ExN06DetectorConstruction_h*/
