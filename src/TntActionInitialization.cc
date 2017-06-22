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
// $Id: TntActionInitialization.cc 68058 2013-03-13 14:47:43Z gcosmo $
//
/// \file TntActionInitialization.cc
/// \brief Implementation of the TntActionInitialization class

#include "TntActionInitialization.hh"

#include "TntPrimaryGeneratorAction.hh"

#include "TntRunAction.hh"
#include "TntEventAction.hh"
#include "TntTrackingAction.hh"
#include "TntSteppingAction.hh"
#include "TntStackingAction.hh"
#include "TntSteppingVerbose.hh"

#include "TntRecorderBase.hh"
#include "TntGlobalParams.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntActionInitialization::TntActionInitialization(TntRecorderBase* recorder)
 : G4VUserActionInitialization(), fRecorder(recorder)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntActionInitialization::~TntActionInitialization()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntActionInitialization::BuildForMaster() const
{
  SetUserAction(new TntRunAction(fRecorder));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntActionInitialization::Build() const
{
	if(TntGlobalParams::Instance()->GetReacFile() == "0") {
		SetUserAction(new TntPrimaryGeneratorAction());
	} else {
		SetUserAction(new TntPGAReaction());
	}

  SetUserAction(new TntStackingAction());

  SetUserAction(new TntRunAction(fRecorder));
  SetUserAction(new TntEventAction(fRecorder));
  SetUserAction(new TntTrackingAction(fRecorder));
  SetUserAction(new TntSteppingAction(fRecorder));
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VSteppingVerbose* TntActionInitialization::InitializeSteppingVerbose() const
{
  return new TntSteppingVerbose();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
