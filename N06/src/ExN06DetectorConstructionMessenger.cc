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
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "ExN06DetectorConstructionMessenger.hh"

#include <sstream>

#include "ExN06DetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcommand.hh"
#include "G4UIparameter.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithABool.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06DetectorConstructionMessenger::ExN06DetectorConstructionMessenger(ExN06DetectorConstruction* manager)
:det(manager)
{
  detDir = new G4UIdirectory("/det/");
  detDir->SetGuidance("detector control");

    armCmd = new G4UIcmdWithADoubleAndUnit("/det/angle",this);
    armCmd->SetGuidance("Set rotation angle of PMT");
    armCmd->SetParameterName("angle",true);
    armCmd->SetRange("angle>=0. && angle<=90");
    armCmd->SetDefaultValue(0.);
    armCmd->SetDefaultUnit("deg");
    
//    LCCmd = new G4UIcmdWithABool("/det/WithLC",this);
//    LCCmd->SetGuidance("Simulation W/ or W/o LC");
//    LCCmd->SetParameterName("LC",true);
//    LCCmd->SetDefaultValue(true);
//
//    updateCmd = new G4UIcommand("/det/update",this);
//    updateCmd->SetGuidance("Update the detector geometry with changed value");
//    updateCmd->SetGuidance("Must be run before beamOn if detector has been changed");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

ExN06DetectorConstructionMessenger::~ExN06DetectorConstructionMessenger()
{
 // delete setheight;
 // delete setseperate;
 // delete setradius;
 // delete setref;
 // delete setpolished;
 // delete setmetal;
  delete armCmd;
//  delete LCCmd;
//  delete updateCmd;
  delete detDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void ExN06DetectorConstructionMessenger::SetNewValue(G4UIcommand* command,G4String newValues)
{
//  if(command == setheight) det->setHeight(newValues);
//  if(command == setseperate) det->setSeperate(newValues);
//  if(command == setradius) det->setRadius(newValues);
//  if(command == setref) det->setRef(newValues);
//  if(command == setpolished) det->setPolished(setpolished->GetNewBoolValue(newValues));
//  if(command == setmetal) det->setMetal(setmetal->GetNewBoolValue(newValues));
    if(command == armCmd) det->setAngle(armCmd->GetNewDoubleValue(newValues));
//    if(command == LCCmd) det->setLC(LCCmd->GetNewBoolValue(newValues));
//    if(command == updateCmd) det->UpdateGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
