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
// $Id: TntEventMessenger.cc 70256 2013-05-28 07:29:30Z gcosmo $
//
/// \file optical/Tnt/src/TntEventMessenger.cc
/// \brief Implementation of the TntEventMessenger class
//
//
#include "TntEventMessenger.hh"
#include "TntEventAction.hh"

#include "G4UIcmdWithABool.hh"
#include "G4UIcmdWithAnInteger.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntEventMessenger::TntEventMessenger(TntEventAction* event)
 : fTntEvent(event)
{
  fSaveThresholdCmd = new G4UIcmdWithAnInteger("/Tnt/saveThreshold",this);
  fSaveThresholdCmd->
SetGuidance("Set the photon count threshold for saving the random number seed");
  fSaveThresholdCmd->SetParameterName("photons",true);
  fSaveThresholdCmd->SetDefaultValue(4500);
  fSaveThresholdCmd->AvailableForStates(G4State_PreInit,G4State_Idle);

  fVerboseCmd = new G4UIcmdWithAnInteger("/Tnt/eventVerbose",this);
  fVerboseCmd->SetGuidance("Set the verbosity of event data.");
  fVerboseCmd->SetParameterName("verbose",true);
  fVerboseCmd->SetDefaultValue(1);

  fPmtThresholdCmd = new G4UIcmdWithAnInteger("/Tnt/pmtThreshold",this);
  fPmtThresholdCmd->SetGuidance("Set the pmtThreshold (in # of photons)");

  fForceDrawPhotonsCmd=new G4UIcmdWithABool("/Tnt/forceDrawPhotons",this);
  fForceDrawPhotonsCmd->SetGuidance("Force drawing of photons.");
  fForceDrawPhotonsCmd
    ->SetGuidance("(Higher priority than /Tnt/forceDrawNoPhotons)");

  fForceDrawNoPhotonsCmd=new G4UIcmdWithABool("/Tnt/forceDrawNoPhotons",this);
  fForceDrawNoPhotonsCmd->SetGuidance("Force no drawing of photons.");
  fForceDrawNoPhotonsCmd
    ->SetGuidance("(Lower priority than /Tnt/forceDrawPhotons)");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntEventMessenger::~TntEventMessenger(){
  delete fSaveThresholdCmd;
  delete fVerboseCmd;
  delete fPmtThresholdCmd;
  delete fForceDrawPhotonsCmd;
  delete fForceDrawNoPhotonsCmd;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntEventMessenger::SetNewValue(G4UIcommand* command, G4String newValue){
  if( command == fSaveThresholdCmd ){
    fTntEvent->SetSaveThreshold(fSaveThresholdCmd->GetNewIntValue(newValue));
  }
  else if( command == fVerboseCmd ){
    fTntEvent->SetEventVerbose(fVerboseCmd->GetNewIntValue(newValue));
  }
  else if( command == fPmtThresholdCmd ){
    fTntEvent->SetPMTThreshold(fPmtThresholdCmd->GetNewIntValue(newValue));
  }
  else if(command == fForceDrawPhotonsCmd){
    fTntEvent->SetForceDrawPhotons(fForceDrawPhotonsCmd
                                  ->GetNewBoolValue(newValue));
  }
  else if(command == fForceDrawNoPhotonsCmd){
    fTntEvent->SetForceDrawNoPhotons(fForceDrawNoPhotonsCmd
                                  ->GetNewBoolValue(newValue));
    G4cout<<"TEST"<<G4endl;
  }
}
