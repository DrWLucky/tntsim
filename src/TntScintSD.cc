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
// $Id: TntScintSD.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/Tnt/src/TntScintSD.cc
/// \brief Implementation of the TntScintSD class
//
//
#include <vector>
#include "TntScintSD.hh"
#include "TntScintHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4LogicalVolume.hh"
#include "G4Track.hh"
#include "G4Step.hh"
#include "G4ParticleDefinition.hh"
#include "G4VTouchable.hh"
#include "G4TouchableHistory.hh"
#include "G4ios.hh"
#include "G4VProcess.hh"
//by Shuya 160407
#include "G4SDManager.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntScintSD::TntScintSD(G4String name, G4String Light)
  : G4VSensitiveDetector(name), Light_Conv(Light)
{
  fScintCollection = NULL;
  collectionName.insert("scintCollection");

//by Shuya 160407
  TntDataOutEV = TntDataRecordTree::TntPointer;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntScintSD::~TntScintSD() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntScintSD::Initialize(G4HCofThisEvent* hitsCE){
  fScintCollection = new TntScintHitsCollection
		(SensitiveDetectorName,collectionName[0]);
  //A way to keep all the hits of this event in one place if needed
  static G4int hitsCID = -1;
  if(hitsCID<0){
    hitsCID = GetCollectionID(0);
  }
  hitsCE->AddHitsCollection( hitsCID, fScintCollection );
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool TntScintSD::ProcessHits(G4Step* aStep,G4TouchableHistory* ){
  G4double edep = aStep->GetTotalEnergyDeposit();
  if(edep==0.) return false; //No edep so dont count as hit

  G4StepPoint* thePrePoint = aStep->GetPreStepPoint();
  G4TouchableHistory* theTouchable =
    (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* thePrePV = theTouchable->GetVolume();
 
  G4StepPoint* thePostPoint = aStep->GetPostStepPoint();

  //Get the average position of the hit
  //G4ThreeVector pos = thePrePoint->GetPosition() + thePostPoint->GetPosition();
  //pos/=2.;
	//by Shuya 160428. To check the pencil beam source has no spreadout. 
  G4ThreeVector pos = thePrePoint->GetPosition();

  TntScintHit* scintHit = new TntScintHit(thePrePV);

  scintHit->SetEdep(edep);
  scintHit->SetPos(pos);

//by Shuya 160408
  //newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() ); // Note - collects reaction point! 
  //newHit->SetPos( aStep->GetTrack()->GetVertexPosition());
  scintHit->SetTrackID( aStep->GetTrack()->GetTrackID() );
  scintHit->SetParentTrackID( aStep->GetTrack()->GetParentID() );
  scintHit->SetTOF( aStep->GetPreStepPoint()->GetGlobalTime() );

  scintHit->SetParticleName( aStep->GetTrack()->GetDefinition()->GetParticleName() );
  scintHit->SetParticleCharge( aStep->GetTrack()->GetDefinition()->GetPDGCharge() );
  scintHit->SetParticleA( aStep->GetTrack()->GetDefinition()->GetBaryonNumber() );

  const G4VProcess* theProcess = aStep->GetTrack()->GetCreatorProcess();
  if(theProcess != 0)
	{
		const G4String theProcessName = theProcess->GetProcessName();
		scintHit->SetParticleProcess(theProcessName);
	}
  else
	{ scintHit->SetParticleProcess("NoReaction"); }
  /* 
		 G4cout << "The Particle Type was : " << newHit->GetParticleName() << G4endl;
		 G4cout << "The Hit occured at " << newHit->GetPos() << G4endl;
		 G4cout << "The Energy Deposited at the Hit was " << edep/MeV << " MeV" << G4endl;
  */

  fScintCollection->insert(scintHit);

//by Shuya 160420
/*
//G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!! " << aStep->GetTrack()->GetTrackID() << "!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;
//if(aStep->GetTrack()->GetTrackID()==1)
//if(aStep->GetTrack()->GetDefinition()->GetParticleName()=="neutron")
//{
G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;
//G4cout << "Hit Number: " << i << G4endl;
G4cout << "the Particle Name: " << aStep->GetTrack()->GetDefinition()->GetParticleName() << G4endl;
G4cout << "edep: " << edep << G4endl;
G4cout << "the Parent TrackID: " << aStep->GetTrack()->GetParentID()  << G4endl;
G4cout << "the TrackID: " << aStep->GetTrack()->GetTrackID() << G4endl;
G4cout << "the CreatorProcess Name: " << aStep->GetTrack()->GetCreatorProcess()->GetProcessName() << G4endl;
G4cout << "Global Time: " << aStep->GetPreStepPoint()->GetGlobalTime()  << G4endl;
G4cout << "First Hit Position: " << pos << G4endl;
//}
*/
	
  return true;

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//by Shuya 160407
//void TntScintSD::EndOfEvent(G4HCofThisEvent* ) 
void TntScintSD::EndOfEvent(G4HCofThisEvent* hitsCE) 
{
//by Shuya 160407
  // Note "EndOfEvent()" is name, NOT EndofEvent !

	G4SDManager* SDMan = G4SDManager::GetSDMpointer();
	G4String detName = "TntScintSD";
	//G4int collectionID = SDMan->GetCollectionID("TntScintSD/TntScintHitsCollection");
	//G4int collectionID = SDMan->GetCollectionID("/TntDet/scintSD/scintCollection");
	G4int collectionID = SDMan->GetCollectionID("scintCollection");
//G4cout << collectionID << "!!!!!!!!!" << G4endl;

	//TntScintHitsCollection* myCollection = (TntScintHitsCollection*)(HCE->GetHC(collectionID));
	//Comment by Shuya. myCollection corresponds to scintHC in EventAction.
	TntScintHitsCollection* myCollection = (TntScintHitsCollection*)(hitsCE->GetHC(collectionID));

	/*
		G4cout << "The EndOfEvent HCE is : " << HCE << G4endl;
		G4cout << "The Event Collection ID is : " << myCollection << G4endl;
		G4cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>End of Event Method is going!" << G4endl;
		G4cout << "The Tnt Pointer Address is" << TntDataOutEV << G4endl;
	*/

	G4double totE = 0;
	G4double EsumProton = 0;
	G4double EsumDeuteron = 0;
	G4double EsumTriton = 0;
	G4double EsumHe3 = 0;
	G4double EsumAlpha = 0;
	G4double EsumC12 = 0;
	G4double EsumC13 = 0;
	G4double EsumEorGamma = 0;
	G4double EsumExotic = 0;

//by Shuya 160502
	G4double totEDep = 0;
	G4double EsumDepProton = 0;
	G4double EsumDepDeuteron = 0;
	G4double EsumDepTriton = 0;
	G4double EsumDepHe3 = 0;
	G4double EsumDepAlpha = 0;
	G4double EsumDepC12 = 0;
	G4double EsumDepC13 = 0;
	G4double EsumDepGamma = 0;
	G4double EsumDepExotic = 0;

//by Shuya 160407 to count the number of photon.
	G4double totPh = 0;
/*
	G4double PhsumProton = 0;
	G4double PhsumDeuteron = 0;
	G4double PhsumTriton = 0;
	G4double PhsumHe3 = 0;
	G4double PhsumAlpha = 0;
	G4double PhsumC12 = 0;
	G4double PhsumC13 = 0;
	G4double PhsumEorGamma = 0;
	G4double PhsumExotic = 0;
*/

	G4double GlobalTime = 0.;
	G4double PulseTime = 0.;  // Time after first hit
	G4ThreeVector FirstHitPos(0.,0.,0.);

	std::vector<TntDataRecordTree::Hit_t> AllHitPosTime(0);
 
 
	G4int number_of_gammahits = 0;
	G4int NumberOfTracks = 0;
	G4int theLastTrackID = 0;

//G4cout << "!!" << G4endl;

	if(myCollection)
	{
    G4int n_hit = myCollection->entries();
//G4cout << n_hit << "!!!!!!" << G4endl;

    // Event readout loop!

    for(int i=0;i<n_hit;i++)
    {
       
      // (*myCollection)[i] is the pointer to the i-th hit of the event.
      // All data analysis output is read out here and then sent to the DataTree!

      TntScintHit* theCurrentHit = (*myCollection)[i];
     
      G4String theParticleName = theCurrentHit->GetParticleName();
      G4double edep = theCurrentHit->GetEdep();
      G4int theParentTrackID = theCurrentHit->GetParentTrackID();
      G4int theTrackID = theCurrentHit->GetTrackID();
      G4int theParticleA = theCurrentHit->GetParticleA();
			bool isNewTrack = false;

			
      if(theTrackID != theLastTrackID && theParticleA >=1)
			{
				NumberOfTracks++;
				theLastTrackID = theTrackID;
				isNewTrack = true;
			}

      G4String theCurrentProcess = theCurrentHit->GetParticleProcess();
      
//COMMENT BY SHUYA 160502. Note this may not correctly give the first hit's information. This is information of one of secondary products created by a first particle (see explanation below). 
      if(i == 0)
			{
				/*   Record Time of "first hit"!  */
				GlobalTime = (theCurrentHit->GetTOF())/ns;
				FirstHitPos = theCurrentHit->GetPos();

				//by Shuya 160428... This should be something about neutron, and I confrimed it was (menate_neutron).
				//G4cout << "CHECK FIRST HIT PROCESS!! " << theCurrentProcess << G4endl;
				//G4cout << "CHECK FIRST HIT POSITION!! " << FirstHitPos << G4endl;
	
//COMMENT BY SHUYA 160502. Note this doesn't necessarily give the first hit energy deposit. This is just an energy deposit by a secondary particles created by a first particle. 
//For example, n+12C->alpha(1st)+alpha(2nd)+alpha(3rd)+n', then this just gives energy deposit of alpha(1st) in its first hit.
//Run and check (by comment out) if you want to understand about this...

				//G4cout << "CHECK FIRST HIT ENERGY DEPOSIT!! " << theCurrentHit->GetEdep() << G4endl;
	
			}  

      PulseTime = (theCurrentHit->GetTOF())/ns;

			if(theParticleName != "opticalphoton") {

				TntDataRecordTree::Hit_t thisHit = {
					theCurrentHit->GetPos()[0],
					theCurrentHit->GetPos()[1],
					theCurrentHit->GetPos()[2],
					PulseTime,
					edep,
					theTrackID
				};

				if(isNewTrack) {
					AllHitPosTime.push_back( thisHit );
				}
			}

			
//by Shuya 160415.

			//if(theParentTrackID==1)
			if(theParticleName != "opticalphoton")
				//if(theTrackID==1)
			{
				G4cout << "!!!!!!!!!!!!!!!!!!!!!!!!! TEST !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << G4endl;
				G4cout << "Hit Number: " << i << G4endl;
				G4cout << "the Particle Name: " << theParticleName << G4endl;
				G4cout << "the Particle A: " << theParticleA << G4endl;
				G4cout << "edep: " << edep << G4endl;
				G4cout << "the Parent TrackID: " << theParentTrackID << G4endl;
				G4cout << "the TrackID: " << theTrackID << G4endl;
				G4cout << "the CurrentProcess: " << theCurrentProcess << G4endl;
				G4cout << "Global Time: " << PulseTime << G4endl;
				G4cout << "Hit Position: " <<  theCurrentHit->GetPos() << G4endl;
			}


//G4cout << theParticleName << G4endl;
			// 3/27/2009 - BTR
			// Need to Convert to Light step by step to get proper response function
			// 
			// Should I add a condition like Pozzi et al. here to say only take
			// first 3 hadron tracks?
			//if(NumberOfTracks <=3)
			//	
			if(theParticleName=="proton")
			{EsumProton += ConvertToLight("proton",1,edep,Light_Conv);
				//G4cout << "!!!!!!!!!!!!" << G4endl
				//by Shuya 160502
				EsumDepProton += edep;
			}
			else if(theParticleName=="deuteron")
			{EsumDeuteron += ConvertToLight("deuteron",1,edep,Light_Conv);
				//by Shuya 160502
				EsumDepDeuteron += edep;
			}
			else if(theParticleName=="triton")
			{EsumTriton += ConvertToLight("triton",1,edep,Light_Conv);
				//by Shuya 160502
				EsumDepTriton += edep;
			}
			else if(theParticleName=="He3")
			{EsumHe3 += ConvertToLight("He3",2,edep,Light_Conv);
				//by Shuya 160502
				EsumDepHe3 += edep;
			}
			else if(theParticleName=="alpha")
			{EsumAlpha += ConvertToLight("alpha",2,edep,Light_Conv);
				//by Shuya 160502
				EsumDepAlpha += edep;
			}
			//Comment by Shuya 160504. NOTE that this is only for ground state and energy loss by excited particles is counted as by Exotic Particles...
			//by Shuya 160504. I changed this to include excited states.
			//else if(theParticleName=="C12[0.0]")
			else if(theParticleName=="C12")
			{EsumC12 += ConvertToLight("C12",6,edep,Light_Conv);
				//by Shuya 160502
				EsumDepC12 += edep;
			}
			//Comment by Shuya 160504. NOTE that this is only for ground state and energy loss by excited particles is counted as by Exotic Particles...
			//by Shuya 160504. I changed this to include excited states.
			//else if(theParticleName=="C13[0.0]")
			else if(theParticleName=="C13")
			{EsumC13 += ConvertToLight("C13",6,edep,Light_Conv);
				//by Shuya 160502
				EsumDepC13 += edep;
			}
			else if(theParticleName=="e-" || 
							theParticleName=="e+" ||
							theParticleName=="gamma")
			{
				EsumEorGamma += ConvertToLight("e-",-1,edep,Light_Conv);
				number_of_gammahits++;
				//by Shuya 160502
				EsumDepGamma += edep;
			}
//by Shuya 160407
			else if(theParticleName=="opticalphoton")
			{totPh++;}
			else
			{
//Comment by Shuya 160502. I assume these particles are Be9, etc...

				EsumExotic += edep;  
				//G4cout << "Exotic Particle Created!!!!! >>> ID = "<< theParticleName << G4endl;
				//G4cout << "The Energy Deposited by the Exotic Particle : " << edep << G4endl; 
	      
				//by Shuya 160502
				EsumDepExotic += edep;
			}
	        
		}   
	}
	else
	{
		/* No hit collection! (although, always a HitCollection even if no Hits ! */
		G4cout << "Warning! No Hits Collection! " << G4endl;
	}

	//G4double BaryonEngEE = EsumProton+EsumDeuteron+EsumTriton+EsumHe3+EsumAlpha+EsumC12+EsumC13+EsumExotic;
	G4double BaryonEngEE = EsumProton+EsumDeuteron+EsumTriton+EsumHe3+EsumAlpha+EsumC12+EsumC13;
	G4double GammaEngEE = EsumEorGamma;

//by Shuya 160502
	G4double BaryonEngEEDep = EsumDepProton+EsumDepDeuteron+EsumDepTriton+EsumDepHe3+EsumDepAlpha+EsumDepC12+EsumDepC13+EsumDepExotic;
	G4double GammaEngEEDep = EsumDepGamma;
  
	//EsumExotic = EsumProton+EsumC12; // Test for 4.4 MeV escape peak - comes from added C12 residual
	// and proton recoil after 12C(n,n')12C* 4.4 gamma decay.

	G4bool gamma_flag = false;

	if(GammaEngEE >= BaryonEngEE)
	{gamma_flag = true;}

	//Add up light energy from hit collection
	//protons, alphas, C12 only
	//for n+gamma discrim.,veto event if there is a gamma hit (for NE213)
	if(number_of_gammahits == 0 || gamma_flag == false)
	{
		//totE = BaryonEngEE;
		totE = BaryonEngEE+GammaEngEE;

//by Shuya 160502
//Comment by Shuya 160502. I have no idea if this should be in the if statement.
		totEDep = BaryonEngEEDep + GammaEngEEDep;
	}
	else
	{
		//G4cout << "totE = 0" << endl;
		totE = 0.; // event gated out

		totEDep = 0.; // event gated out
	}


		
//G4cout << "!!!!" << G4endl;
	// Send data from hit collection to DataRecordTree and FillTree and text files!
 
	TntDataOutEV->senddataEV(1,totE);
	TntDataOutEV->senddataEV(2,EsumProton);
	TntDataOutEV->senddataEV(3,EsumAlpha);
	TntDataOutEV->senddataEV(4,EsumC12);
	TntDataOutEV->senddataEV(5,EsumEorGamma);
	TntDataOutEV->senddataEV(6,EsumExotic);
//by Shuya 160407
//by Shuya 160426. I moved this to the end of EventAction.
// TntDataOutEV->senddataEV(7,totPh);

//by Shuya 160502. These are ENERGY (not converted to light output) deposited in Scintillation Detector.
	TntDataOutEV->senddataEV(10,totEDep);
	TntDataOutEV->senddataEV(11,EsumDepProton);
	TntDataOutEV->senddataEV(12,EsumDepAlpha);
	TntDataOutEV->senddataEV(13,EsumDepC12);
	TntDataOutEV->senddataEV(14,EsumDepGamma);
	TntDataOutEV->senddataEV(15,EsumDepExotic);

	TntDataOutEV->senddataTOF(GlobalTime);
	TntDataOutEV->senddataPosition(FirstHitPos);
	TntDataOutEV->senddataHits(AllHitPosTime, true);
 
//by Shuya 160427. I restored this function.
//by Shuya 160426. I removed this and Set FillTree at the end of EventAction, where PMT photons data are collected and sent TntDataRecord. (because in an event, ScintSD->EventAction is the order)
// TntDataOutEV->FillTree();


	// G4cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>End of Event Analysis!" << G4endl;

}

G4double TntScintSD::ConvertToLight(G4String theName, G4double theCharge, G4double edep,G4String Light_Conv)
{
	//Converts Energy deposited in hit using Cecil's formula.

//G4cout << "!!!! " << Light_Conv << G4endl;
	G4double theLightEng = 0.;

	if(Light_Conv == "none")
	{
		theLightEng = edep;
	}
	else if(Light_Conv == "light-conv" || Light_Conv == "light+resol")
	{     
//G4cout << "!!!!!!!!!!!!!!  " << edep << G4endl;
		if(theCharge == 1. && theName != "e+")
		{ theLightEng = 0.83*edep-2.82*(1-exp(-0.25*pow(edep,0.93))); }
		else if(theCharge == 2.)
		{ theLightEng = 0.41*edep-5.9*(1-exp(-0.065*pow(edep,1.01)));}
		else if(theCharge == 3.)
		{ theLightEng = 0.1795*( edep );}  // Obtained from EXP fit of measured leading coeffs
		else if(theCharge == 4.)
		{ theLightEng = 0.0821*( edep );}  // Obtained from EXP fit of measured leading coeffs
		else if(theCharge == 5.)
		{ theLightEng =  0.0375*( edep );} // Obtained from EXP fit of measured leading coeffs
		else if(theCharge == 6.)
		{ theLightEng = 0.017*( edep );}
		else if(theName =="e-" || theName =="e+")
		{ theLightEng = edep;}
		else if (theName =="gamma")
		{ theLightEng = edep; }

		// Include Detector light collection resolution
		// Follows method of Daniel Cano student thesis
		// and Dekempener et al. NIM A 256 (1987) 489-498

		if(Light_Conv == "light+resol")
		{
			G4double AlphaRes = 0.045;
			G4double BetaRes = 0.075;
			G4double GammaRes = 0.002;
			G4double Argument = pow(AlphaRes,2)*pow(theLightEng,2)+pow(BetaRes,2)*theLightEng+pow(GammaRes,2);    
			if(Argument > 0.)
			{
				G4double FWHMFact = 2.35482;
				// 3/27/2009 - Should this be the resolution?
				// Need to look in Dekempenner et al., NIM A 256 (1987) 489-498
				// or S. Pozzi NIM A 582 (2007)
				G4double FWHMRes = (sqrt(Argument)/FWHMFact);
				theLightEng += CLHEP::RandGauss::shoot(0.,FWHMRes);
			}
		}
		// Remove events where you get Light<0 due to resolution
		if(theLightEng <= 0.*MeV)
		{theLightEng = 0.*MeV;}
	}
	else
	{
		G4cout << "*********************************************" << G4endl;
		G4cerr << "Bad choice for light conversion method in TntScintSD::ConvertToLight() method." << G4endl;
		G4cerr << "Proper Values include \"none\", \"light-conv\", or \"light+resol\"." << G4endl;
		G4ExceptionSeverity severity=FatalException;
		G4Exception("Exception thrown -- Program stopped in TntScintSD::ConvertToLight() method.", "Error",severity,"Error!");
		//Note: the enum G4ExceptionSeverity has been defined since Geant4 release 5.0 and its values are: FatalException, FatalErrorInArgument, EventMustBeAborted, FatalException, and JustWarning. 
		//G4Exception("Exception thrown -- Program stopped in TntSD::ConvertToLight() method.");
	} 
    
	return theLightEng;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntScintSD::clear() {} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntScintSD::DrawAll() {} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntScintSD::PrintAll() {} 
