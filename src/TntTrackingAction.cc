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
// $Id: TntTrackingAction.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/Tnt/src/TntTrackingAction.cc
/// \brief Implementation of the TntTrackingAction class
//
//
#include "TntTrajectory.hh"
#include "TntTrackingAction.hh"
#include "TntUserTrackInformation.hh"
#include "TntDetectorConstruction.hh"
#include "TntRecorderBase.hh"

#include "G4TrackingManager.hh"
#include "G4Track.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntTrackingAction::TntTrackingAction(TntRecorderBase* r)
  : fRecorder(r) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntTrackingAction::PreUserTrackingAction(const G4Track* aTrack)
{
  //Let this be up to the user via vis.mac
  //  fpTrackingManager->SetStoreTrajectory(true);

  //Use custom trajectory class
  fpTrackingManager->SetTrajectory(new TntTrajectory(aTrack));

  //This user track information is only relevant to the photons
  fpTrackingManager->SetUserTrackInformation(new TntUserTrackInformation);

  /*  const G4VProcess* creator = aTrack->GetCreatorProcess();
  if(creator)
    G4cout<<creator->GetProcessName()<<G4endl;
  */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntTrackingAction::PostUserTrackingAction(const G4Track* aTrack){
  TntTrajectory* trajectory=(TntTrajectory*)fpTrackingManager->GimmeTrajectory();
  TntUserTrackInformation*
    trackInformation=(TntUserTrackInformation*)aTrack->GetUserInformation();

//by Shuya 160502. Initially declared in Tnt.cc (main file).
    extern G4int NumOfCreatedPhotons;

  //Lets choose to draw only the photons that hit the sphere and a pmt
  if(aTrack->GetDefinition()==G4OpticalPhoton::OpticalPhotonDefinition()){

	//by Shuya 160502. To count the number of created photons whether detected or not by PMTs.
	//G4cout << "TEST OF CREATED PHOTON NUMBERS " << G4endl;
	NumOfCreatedPhotons++;

    const G4VProcess* creator=aTrack->GetCreatorProcess();
    if(creator && creator->GetProcessName()=="OpWLS"){
      trajectory->WLS();
      trajectory->SetDrawTrajectory(true);
    }

    if(TntDetectorConstruction::GetSphereOn()){
      if((trackInformation->GetTrackStatus()&hitPMT)&&
         (trackInformation->GetTrackStatus()&hitSphere)){
        trajectory->SetDrawTrajectory(true);
      }
    }
    else{
      if(trackInformation->GetTrackStatus()&hitPMT)
        trajectory->SetDrawTrajectory(true);
    }
  }
  else //draw all other trajectories
    trajectory->SetDrawTrajectory(true);

  if(trackInformation->GetForceDrawTrajectory())
    trajectory->SetDrawTrajectory(true);

  if(fRecorder)fRecorder->RecordTrack(aTrack);
}
