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
// $Id: TntSteppingAction.cc 73915 2013-09-17 07:32:26Z gcosmo $
//
/// \file optical/Tnt/src/TntSteppingAction.cc
/// \brief Implementation of the TntSteppingAction class
//
//
#include "TntSteppingAction.hh"
#include "TntEventAction.hh"
#include "TntTrackingAction.hh"
#include "TntTrajectory.hh"
#include "TntPMTSD.hh"
#include "TntUserTrackInformation.hh"
#include "TntUserEventInformation.hh"
#include "TntSteppingMessenger.hh"
#include "TntRecorderBase.hh"

#include "G4SteppingManager.hh"
#include "G4SDManager.hh"
#include "G4EventManager.hh"
#include "G4ProcessManager.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4Event.hh"
#include "G4StepPoint.hh"
#include "G4TrackStatus.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntSteppingAction::TntSteppingAction(TntRecorderBase* r)
  : fRecorder(r),fOneStepPrimaries(false)
{
  fSteppingMessenger = new TntSteppingMessenger(this);

  fExpectedNextStatus = Undefined;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntSteppingAction::~TntSteppingAction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntSteppingAction::UserSteppingAction(const G4Step * theStep){

  G4Track* theTrack = theStep->GetTrack();

  if ( theTrack->GetCurrentStepNumber() == 1 ) fExpectedNextStatus = Undefined;
 
  TntUserTrackInformation* trackInformation
    =(TntUserTrackInformation*)theTrack->GetUserInformation();
  TntUserEventInformation* eventInformation
    =(TntUserEventInformation*)G4EventManager::GetEventManager()
    ->GetConstCurrentEvent()->GetUserInformation();

  G4StepPoint* thePrePoint = theStep->GetPreStepPoint();
  G4VPhysicalVolume* thePrePV = thePrePoint->GetPhysicalVolume();

  G4StepPoint* thePostPoint = theStep->GetPostStepPoint();
  G4VPhysicalVolume* thePostPV = thePostPoint->GetPhysicalVolume();

  G4OpBoundaryProcessStatus boundaryStatus=Undefined;
  static G4ThreadLocal G4OpBoundaryProcess* boundary=NULL;

  //find the boundary process only once
  if(!boundary){
    G4ProcessManager* pm
      = theStep->GetTrack()->GetDefinition()->GetProcessManager();
    G4int nprocesses = pm->GetProcessListLength();
    G4ProcessVector* pv = pm->GetProcessList();
    G4int i;
    for( i=0;i<nprocesses;i++){
      if((*pv)[i]->GetProcessName()=="OpBoundary"){
        boundary = (G4OpBoundaryProcess*)(*pv)[i];
        break;
      }
    }
  }

//by Shuya 160420
//G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!! " << theStep->GetTrack()->GetTrackID() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;
	//if(theStep->GetTrack()->GetTrackID()==0)
	//if(theStep->GetTrack()->GetParentID()==0)
	if(theStep->GetTrack()->GetDefinition()->GetParticleName()!="opticalphoton")
	{
# if 0
		G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;
		//G4cout << "Hit Number: " << i << G4endl;
		G4cout << "the Particle Name: " << theStep->GetTrack()->GetDefinition()->GetParticleName() << G4endl;
		G4cout << "Prestep kinetic energy: " << theStep->GetTrack()->GetKineticEnergy() + theStep->GetTotalEnergyDeposit() << G4endl;
		G4cout << "kinetic energy (MeV, Poststep): " << theStep->GetTrack()->GetKineticEnergy() << G4endl;
		G4cout << "edep in total: " << theStep->GetDeltaEnergy() << G4endl;
		G4cout << "edep in medium: " << theStep->GetTotalEnergyDeposit() << G4endl;
		G4cout << "the Parent TrackID: " << theStep->GetTrack()->GetParentID()  << G4endl;
		G4cout << "the TrackID: " << theStep->GetTrack()->GetTrackID() << G4endl;
		//Comment by Shuya 160502. I couldn't make it (Got Segmentation fault). If you want to make it, see below (inside if(theTrack->GetParentID()..., creator, creatorName will be a hint)
		//However, without this, you can still guess what made this reaction from other outputs.
		//G4cout << "the CreatorProcess: " << theStep->GetTrack()->GetCreatorProcess()->GetProcessName() << G4endl;
		G4cout << "Global Time: " << theStep->GetPreStepPoint()->GetGlobalTime()  << G4endl;
		G4cout << "Prestep Position (mm): " << theStep->GetPreStepPoint()->GetPosition() << G4endl;
		G4cout << "Poststep Position (mm): " << theStep->GetPostStepPoint()->GetPosition() << G4endl;
#endif
	}


  if(theTrack->GetParentID()==0){
    //This is a primary track
 
    G4TrackVector* fSecondary=fpSteppingManager->GetfSecondary();
    G4int tN2ndariesTot = fpSteppingManager->GetfN2ndariesAtRestDoIt()
      + fpSteppingManager->GetfN2ndariesAlongStepDoIt()
      + fpSteppingManager->GetfN2ndariesPostStepDoIt();

    //If we havent already found the conversion position and there were
    //secondaries generated, then search for it
    if(!eventInformation->IsConvPosSet() && tN2ndariesTot>0 ){
      for(size_t lp1=(*fSecondary).size()-tN2ndariesTot;
          lp1<(*fSecondary).size(); lp1++){
        const G4VProcess* creator=(*fSecondary)[lp1]->GetCreatorProcess();
        if(creator){
          G4String creatorName=creator->GetProcessName();
          if(creatorName=="phot"||creatorName=="compt"||creatorName=="conv"){
            //since this is happening before the secondary is being tracked
            //the Vertex position has not been set yet(set in initial step)
            eventInformation->SetConvPos((*fSecondary)[lp1]->GetPosition());
	//by Shuya 160419
	//G4cout << "!! AAAAAAA " << G4endl;
          }
        }
      }
    }

//Comments by Shuya 160418. I have no idea how this statement is working so far. 
//See SteppingManager class to see how set fOneStepPrimaries.
    if(fOneStepPrimaries&&thePrePV->GetName()=="scintillator")
	{
	//Comments by Shuya 160418. I checked with this G4cout, and this if statement is not working.
	//G4cout << "!! OOOOOO " << G4endl;
      theTrack->SetTrackStatus(fStopAndKill);
	}

//by Shuya 160418. To cut off all the secondary reactions from the primary particle (neutron).
//but it didn't work....
/*
	//if(tN2ndariesTot > 0)
	if((*fSecondary).size()>0)ost{
	G4cout << "!! AAAAAAA " << G4endl;
      theTrack->SetTrackStatus(fStopAndKill);
      //theTrack->SetTrackStatus(fStopButAlive);
      //theTrack->SetTrackStatus(fKillTrackAndSecondaries);
	}
*/

  }

  if(!thePostPV){//out of world
    fExpectedNextStatus=Undefined;
    return;
  }

  G4ParticleDefinition* particleType = theTrack->GetDefinition();
  if(particleType==G4OpticalPhoton::OpticalPhotonDefinition()){
    //Optical photon only

    if(thePrePV->GetName()=="Slab")
      //force drawing of photons in WLS slab
      trackInformation->SetForceDrawTrajectory(true);
    else if(thePostPV->GetName()=="expHall")
      //Kill photons entering expHall from something other than Slab
      theTrack->SetTrackStatus(fStopAndKill);

    //Was the photon absorbed by the absorption process
    if(thePostPoint->GetProcessDefinedStep()->GetProcessName()
       =="OpAbsorption"){
      eventInformation->IncAbsorption();
      trackInformation->AddTrackStatusFlag(absorbed);
    }

    boundaryStatus=boundary->GetStatus();

    //Check to see if the partcile was actually at a boundary
    //Otherwise the boundary status may not be valid
    //Prior to Geant4.6.0-p1 this would not have been enough to check
    if(thePostPoint->GetStepStatus()==fGeomBoundary){
      if(fExpectedNextStatus==StepTooSmall){
        if(boundaryStatus!=StepTooSmall){
          G4ExceptionDescription ed;
          ed << "TntSteppingAction::UserSteppingAction(): "
                << "No reallocation step after reflection!"
                << G4endl;
          G4Exception("TntSteppingAction::UserSteppingAction()", "TntExpl01",
          FatalException,ed,
          "Something is wrong with the surface normal or geometry");
        }
      }

//COMMENT by Shuya 160504. All the statements below are included if statement above (if(thePostpoint->.....==fGeomBoundary){}), so opticalphotons under consideration here are in boundary.)
      fExpectedNextStatus=Undefined;
      switch(boundaryStatus){
      case Absorption:
        trackInformation->AddTrackStatusFlag(boundaryAbsorbed);
        eventInformation->IncBoundaryAbsorption();
        break;
      case Detection: //Note, this assumes that the volume causing detection
                      //is the photocathode because it is the only one with
                      //non-zero efficiency
                      //<---Suppled by Shuya (note explanation above is by the original author) 160502. non-zero efficiency means if(efficiency=true).
        {
        //Triger sensitive detector manually since photon is
        //absorbed but status was Detection
        G4SDManager* SDman = G4SDManager::GetSDMpointer();
        G4String sdName="/TntDet/pmtSD";
        TntPMTSD* pmtSD = (TntPMTSD*)SDman->FindSensitiveDetector(sdName);
        if(pmtSD)pmtSD->ProcessHits_constStep(theStep,NULL);
        trackInformation->AddTrackStatusFlag(hitPMT);
        break;
        }
      case FresnelReflection:
      case TotalInternalReflection:
      case LambertianReflection:
      case LobeReflection:
      case SpikeReflection:
      case BackScattering:
        trackInformation->IncReflections();
        fExpectedNextStatus=StepTooSmall;
        break;
      default:
        break;
      }
      if(thePostPV->GetName()=="sphere")
        trackInformation->AddTrackStatusFlag(hitSphere);
    }
  }

  if(fRecorder)fRecorder->RecordStep(theStep);
}
