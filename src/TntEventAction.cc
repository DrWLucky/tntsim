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
// $Id: TntEventAction.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/Tnt/src/TntEventAction.cc
/// \brief Implementation of the TntEventAction class
//
//
#include "TntEventAction.hh"
#include "TntScintHit.hh"
#include "TntPMTHit.hh"
#include "TntUserEventInformation.hh"
#include "TntTrajectory.hh"
#include "TntRecorderBase.hh"
#include "TntGlobalParams.hh"

#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "G4UImanager.hh"
#include "G4SystemOfUnits.hh"

//by Shuya 160408
#include "TntDataRecordTree.hh"




//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntEventAction::TntEventAction(TntRecorderBase* r)
  : fRecorder(r),fSaveThreshold(0),fScintCollID(-1),fPMTCollID(-1),fVerbose(1),
		fPMTThreshold(1),fForcedrawphotons(false),fForcenophotons(false),
//by Shuya 160407
		numberOfEvent(-1)
{
  fEventMessenger = new TntEventMessenger(this);

  TntGetDataEV = TntDataRecordTree::TntPointer;
//by Shuya 160408
  TntDataOutEV = TntDataRecordTree::TntPointer;
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntEventAction::~TntEventAction(){}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntEventAction::BeginOfEventAction(const G4Event* anEvent){
 

  //New event, add the user information object
  G4EventManager::
    GetEventManager()->SetUserInformation(new TntUserEventInformation);

  G4SDManager* SDman = G4SDManager::GetSDMpointer();
  if(fScintCollID<0)
    fScintCollID=SDman->GetCollectionID("scintCollection");
  if(fPMTCollID<0)
    fPMTCollID=SDman->GetCollectionID("pmtHitCollection");

  if(fRecorder)fRecorder->RecordBeginOfEvent(anEvent);
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntEventAction::EndOfEventAction(const G4Event* anEvent){


//by Shuya 160421
	extern G4int Counter;
	//G4cout << "COUNTER " << Counter << G4endl;
//by Shuya 160502
	extern G4int NumOfCreatedPhotons;
//by Shuya 160509
	G4int npmtX = TntGlobalParams::Instance()->GetNumPmtX();
	G4int npmtY = TntGlobalParams::Instance()->GetNumPmtY();
	

//by Shuya 160421
  numberOfEvent++;
	//G4cout << "NUMOFEVE " << numberOfEvent << G4endl;

  TntUserEventInformation* eventInformation
    =(TntUserEventInformation*)anEvent->GetUserInformation();
 
  G4TrajectoryContainer* trajectoryContainer=anEvent->GetTrajectoryContainer();
 
  G4int n_trajectories = 0;
  if (trajectoryContainer) n_trajectories = trajectoryContainer->entries();
//G4cout << "!!!  " << n_trajectories << G4endl;

  // extract the trajectories and draw them
  if (G4VVisManager::GetConcreteInstance()){
    for (G4int i=0; i<n_trajectories; i++){
      TntTrajectory* trj = (TntTrajectory*)
        ((*(anEvent->GetTrajectoryContainer()))[i]);
      if(trj->GetParticleName()=="opticalphoton"){
        trj->SetForceDrawTrajectory(fForcedrawphotons);
        trj->SetForceNoDrawTrajectory(fForcenophotons);
      }
      trj->DrawTrajectory();
    }
  }
 
  TntScintHitsCollection* scintHC = 0;
  TntPMTHitsCollection* pmtHC = 0;
  G4HCofThisEvent* hitsCE = anEvent->GetHCofThisEvent();
 
  //Get the hit collections
  if(hitsCE){
    if(fScintCollID>=0)scintHC = (TntScintHitsCollection*)(hitsCE->GetHC(fScintCollID));
    if(fPMTCollID>=0)pmtHC = (TntPMTHitsCollection*)(hitsCE->GetHC(fPMTCollID));
//by Shuya 160322.
    //if(fPMTCollID>=0)G4cout << "PMT No!!! = " << fPMTCollID << G4endl;
  }

  //Hits in scintillator
  if(scintHC){
    int n_hit = scintHC->entries();
//by Shuya 160406.
    //G4cout << "n_hit SCINT!!! = " << n_hit << G4endl;
    G4ThreeVector  eWeightPos(0.);
    G4double edep;
    G4double edepMax=0;

    for(int i=0;i<n_hit;i++){ //gather info on hits in scintillator
      edep=(*scintHC)[i]->GetEdep();
      eventInformation->IncEDep(edep); //sum up the edep
      eWeightPos += (*scintHC)[i]->GetPos()*edep;//calculate energy weighted pos
      if(edep>edepMax){
        edepMax=edep;//store max energy deposit
        G4ThreeVector posMax=(*scintHC)[i]->GetPos();
        eventInformation->SetPosMax(posMax,edep);
      }
    }

//by Shuya 160502
		if(fVerbose >0){
			G4cout << "================OUTPUT SENT FROM TntEventAction Part1====================================" << G4endl;
		}

    if(eventInformation->GetEDep()==0.){
      if(fVerbose>0) G4cout<<"No hits in the scintillator this event."<<G4endl;
    }
    else{
      //Finish calculation of energy weighted position
      eWeightPos/=eventInformation->GetEDep();
      eventInformation->SetEWeightPos(eWeightPos);
      if(fVerbose>0){
        G4cout << "\tEnergy weighted position of hits in Tnt : "
               << eWeightPos/mm << G4endl;
      }
    }
    if(fVerbose>0){
//Comment by Shuya 160502. This should be equal to Tot_sumDep in DataRecord.cc (and I confirmed it is).
			G4cout << "\tTotal energy deposition in scintillator : "
						 << eventInformation->GetEDep() / keV << " (keV)" << G4endl;
    }
  }
 
//by Shuya 160502. I put this here in case PMT count is 0.
  TntDataOutEV->createdataPMT(Counter);

//by Shuya 160502. I moved these from inside if statement of (pmtHC).
	G4int pmtphotonfrontsum = 0;
	G4int pmtphotonbacksum = 0;

  if(pmtHC){
    G4ThreeVector reconPos(0.,0.,0.);
    G4int pmts=pmtHC->entries();
//by Shuya 160404.
//G4cout << pmts << "!" << G4endl;
//by Shuya 160322
    //G4int pmtphotonfrontsum = 0;
    //G4int pmtphotonbacksum = 0;

    //Gather info from all PMTs
    for(G4int i=0;i<pmts;i++){
//by Shuya 160408
			G4int pmtnumber=(*pmtHC)[i]->GetPMTNumber();
//by Shuya 160421
			//if(i==0) TntDataOutEV->createdataPMT(pmtnumber, anEvent->GetEventID());
			//if(i==0) TntDataOutEV->createdataPMT(pmtnumber, numberOfEvent);
			//by Shuya 160502
			//if(i==0) TntDataOutEV->createdataPMT(Counter);

      eventInformation->IncHitCount((*pmtHC)[i]->GetPhotonCount());
//by Shuya 160322
			//G4cout << i << " " << (*pmtHC)[i]->GetPhotonCount() << G4endl;
			//pmtphotonsum += (*pmtHC)[i]->GetPhotonCount();
			//by Shuya 160427
			//if(pmtnumber<100)	pmtphotonfrontsum += (*pmtHC)[i]->GetPhotonCount();
			//else if(pmtnumber>=100 && pmtnumber<200)	pmtphotonbacksum += (*pmtHC)[i]->GetPhotonCount();
			//by Shuya 160509
			if(pmtnumber<(npmtX*npmtY))	pmtphotonfrontsum += (*pmtHC)[i]->GetPhotonCount();
			else if(pmtnumber>=(npmtX*npmtY) && pmtnumber<(2*npmtX*npmtY))	pmtphotonbacksum += (*pmtHC)[i]->GetPhotonCount();
			//Comment by Shuya 160428. If you want to check the pmt count, remove the comment out below.
			//G4cout << pmtnumber << " " << (*pmtHC)[i]->GetPhotonCount() << G4endl;

//by Shuya 160421
			//TntDataOutEV->senddataPMT(pmtnumber,(*pmtHC)[i]->GetPhotonCount());
			//TntDataOutEV->senddataPMT(pmtnumber,(*pmtHC)[i]->GetPhotonCount(),anEvent->GetEventID());
			//TntDataOutEV->senddataPMT(pmtnumber,(*pmtHC)[i]->GetPhotonCount(),numberOfEvent);
			TntDataOutEV->senddataPMT(pmtnumber,(*pmtHC)[i]->GetPhotonCount(),Counter);

			for(std::vector<G4double>::const_iterator iPhot = (*pmtHC)[i]->GetPhotonTime().begin();
					iPhot != (*pmtHC)[i]->GetPhotonTime().end(); ++iPhot) {
				TntDataOutEV->senddataPMT_Time(pmtnumber, *iPhot);
			}
			
      reconPos+=(*pmtHC)[i]->GetPMTPos()*(*pmtHC)[i]->GetPhotonCount();
			//G4cout << "CHECKING RECONPOS " << reconPos << " CHECKING PMTPOS " << (*pmtHC)[i]->GetPMTPos() << " " << "CHECKING PHOTONCOUNTS " << (*pmtHC)[i]->GetPhotonCount() << G4endl;
      if((*pmtHC)[i]->GetPhotonCount()>=fPMTThreshold){
        eventInformation->IncPMTSAboveThreshold();
      }
      else{//wasnt above the threshold, turn it back off
        (*pmtHC)[i]->SetDrawit(false);
      }
    }
		//G4cout << "PMT Front Sum " << pmtphotonfrontsum << G4endl;
		//G4cout << "PMT Back Sum " << pmtphotonbacksum << G4endl;
		//by Shuya 160426
		//TntDataOutEV->senddataEV(7,(double)pmtphotonfrontsum);
		//TntDataOutEV->senddataEV(8,(double)pmtphotonbacksum);
 
    if(eventInformation->GetHitCount()>0){//dont bother unless there were hits
      reconPos/=eventInformation->GetHitCount();
      if(fVerbose>0){
        G4cout << "\tReconstructed position of hits in Tnt : "
               << reconPos/mm << G4endl;
      }
      eventInformation->SetReconPos(reconPos);
    }
    pmtHC->DrawAllHits();
  }

	//by Shuya 160502. I moved these from inside if statement of (pmtHC).
 	TntDataOutEV->senddataEV(7,(double)pmtphotonfrontsum);
 	TntDataOutEV->senddataEV(8,(double)pmtphotonbacksum);

	//by Shuya 160502. NumOfCreatedPhotons are counted in TrackingAction and now sending data to DataRecord.cc, and then initialization
 	TntDataOutEV->senddataEV(9,(double)NumOfCreatedPhotons);

 	TntDataOutEV->senddataEV(16,(double)eventInformation->GetAbsorptionCount());
 	TntDataOutEV->senddataEV(17,(double)eventInformation->GetBoundaryAbsorptionCount());

//by Shuya 160427. I moved the function here to fill the pmtphotonsum data which is stored above.
	TntDataOutEV->FillTree();
//by Shuya 160422
 	TntDataOutEV->FillTree2(Counter);

  if(fVerbose>0){
//by Shuya 160407
		//G4int event_show = 10000;
		G4int event_show = 1;   
		if( ((numberOfEvent) % event_show == 0)  )
		{
			// G4cout << "====================================================" << G4endl;
      G4cout << "================OUTPUT SENT FROM TntEventAction Part2====================================" << G4endl;
      G4cout << ">>> Event (from GetEventID())" << anEvent->GetEventID() << G4endl;
      G4cout << ">>> Event (By Shuya 160426)" << Counter << G4endl;
      G4cout << ">>> Note Event from GetEventID() would't be correct because of multithreading" << G4endl;

      //G4int n_hits = myCollection->entries();
      G4int n_hits = scintHC->entries();
      G4cout << n_hits << " Hits were registered in this event! " << G4endl; 
      TntGetDataEV->ShowDataFromEvent();
      //  G4cout << "The Event Action Tnt Data Out pointer is: " << TntDataOutEV << G4endl;
		}

//by Shuya 160502
		G4cout << "================OUTPUT SENT FROM TntEventAction Part3====================================" << G4endl;
    //End of event output. later to be controlled by a verbose level
    G4cout << "\tNumber of photons that hit PMTs in this event : "
           << eventInformation->GetHitCount() << G4endl;
    G4cout << "\tNumber of PMTs above threshold("<<fPMTThreshold<<") : "
           << eventInformation->GetPMTSAboveThreshold() << G4endl;
//by Shuya 160502
    G4cout << "\tNumber of photons produced in this event (counting in TrackAction, shoud be equal to one below): "
           << NumOfCreatedPhotons << G4endl;
    G4cout << "\tNumber of photons produced by scintillation in this event : "
           << eventInformation->GetPhotonCount_Scint() << G4endl;
    G4cout << "\tNumber of photons produced by cerenkov in this event : "
           << eventInformation->GetPhotonCount_Ceren() << G4endl;
    G4cout << "\tNumber of photons absorbed (OpAbsorption) in this event : "
           << eventInformation->GetAbsorptionCount() << G4endl;
    G4cout << "\tNumber of photons absorbed at boundaries (OpBoundary) in "
           << "this event : " << eventInformation->GetBoundaryAbsorptionCount()
           << G4endl;
    G4cout << "Unacounted for photons in this event : "
           << (eventInformation->GetPhotonCount_Scint() +
               eventInformation->GetPhotonCount_Ceren() -
               eventInformation->GetAbsorptionCount() -
               eventInformation->GetHitCount() -
               eventInformation->GetBoundaryAbsorptionCount())
           << G4endl;

//by Shuya 160502
		G4cout << "================END OF EVENT " << Counter << " ====================================" << G4endl;
		G4cout << G4endl;
		G4cout << G4endl;
		G4cout << G4endl;

		if(Counter%1000 == 0) {
			G4cerr << "Event " << Counter << G4endl;
		}
  }

//by Shuya 160502
	NumOfCreatedPhotons = 0;

  //If we have set the flag to save 'special' events, save here
  if(fSaveThreshold&&eventInformation->GetPhotonCount() <= fSaveThreshold)
    G4RunManager::GetRunManager()->rndmSaveThisEvent();

  if(fRecorder)fRecorder->RecordEndOfEvent(anEvent);

//by Shuya 160421
	Counter++;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntEventAction::SetSaveThreshold(G4int save){
/*Sets the save threshold for the random number seed. If the number of photons
	generated in an event is lower than this, then save the seed for this event
	in a file called run###evt###.rndm
*/
  fSaveThreshold=save;
  G4RunManager::GetRunManager()->SetRandomNumberStore(true);
  G4RunManager::GetRunManager()->SetRandomNumberStoreDir("random/");
  //  G4UImanager::GetUIpointer()->ApplyCommand("/random/setSavingFlag 1");
}
