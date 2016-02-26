//==========================================================================
// TntSD.cc
// Based on DemonScintSD.cc and ExN04TrackerSD.cc
//
// Written/Modified by: Brian Roeder, LPC Caen 02/14/07
//                      email - roeder@lpccaen.in2p3.fr
//
// Purpose: Defines TntSD member functions for simulation data readout
//          Sets and accesses Data in TntDetHit.hh methods
//
//==========================================================================
//
// - See UM Hits Presentation and JLab Scoring 2 Talk for more info.
//
//


#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"

#include "TntDetHit.hh"
#include "G4VPhysicalVolume.hh"
#include "G4TouchableHistory.hh"

#include "G4RunManager.hh"
#include "Randomize.hh"
#include "G4VProcess.hh"

#include "TntSD.hh"


TntSD::TntSD(G4String name, G4String Light) : G4VSensitiveDetector(name), Light_Conv(Light)
{
  /* Constructor of Detector Method */
  // Note: Constructor runs before "run" begins
  // G4cout << "You are Using a Tnt SD! " << G4endl;
  G4String HCname;
  collectionName.insert(HCname="TntHitsCollection"); // Set name of Event Hits Collection
  TntDataOutEV = TntDataRecordTree::TntPointer;
}

TntSD::~TntSD()
{/* Destructor */ }

void TntSD::Initialize(G4HCofThisEvent* HCE)    // Member Function of TntSD
{
  //- This method is run at the beginning of each event
  //- Sets up HitsCollection, etc.
  //- Note that "GetCollectionID()" is a slow op. - Following Demon Setup!

  TntHitsCollection = new TntDetHitsCollection(SensitiveDetectorName,collectionName[0]);
  static int DetCollectionID = -1;

  if(DetCollectionID<0)
    { DetCollectionID = GetCollectionID(0); }  // Sets Hit CollectionID

  HCE->AddHitsCollection( DetCollectionID, TntHitsCollection ); // Sets Events Hits Collection

  //The below is for debugging purposes
  /*
   G4cout << "DetCollection ID = " << DetCollectionID << G4endl;
   G4cout << "The Detector was initialized!" << G4endl;
   G4cout << "This is the address of the Tnt Hits Collection : " << TntHitsCollection << G4endl;
   G4cout << "This is the address of the  ->HitsCollection (HCE) : " << HCE << G4endl;
  */
}

G4bool TntSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
  // Not too much different from ExN04, see also JLab ScoringII talk
  // Should go automatically once detector is defined as sensitive in DetectorConstruction!

  G4double edep = aStep->GetTotalEnergyDeposit();

  if(edep == 0.)
    {return false;}

  TntDetHit* newHit = new TntDetHit();

  newHit->SetEdep( edep );
  newHit->SetPos( aStep->GetPreStepPoint()->GetPosition() ); // Note - collects reaction point! 
  //newHit->SetPos( aStep->GetTrack()->GetVertexPosition());
  newHit->SetTrackID( aStep->GetTrack()->GetTrackID() );
  newHit->SetParentTrackID( aStep->GetTrack()->GetParentID() );
  newHit->SetTOF( aStep->GetPreStepPoint()->GetGlobalTime() );

  newHit->SetParticleName( aStep->GetTrack()->GetDefinition()->GetParticleName() );
  newHit->SetParticleCharge( aStep->GetTrack()->GetDefinition()->GetPDGCharge() );
  newHit->SetParticleA( aStep->GetTrack()->GetDefinition()->GetBaryonNumber() );

  const G4VProcess* theProcess = aStep->GetTrack()->GetCreatorProcess();
  if(theProcess != 0)
    {
      const G4String theProcessName = theProcess->GetProcessName();
      newHit->SetParticleProcess(theProcessName);
    }
  else
    { newHit->SetParticleProcess("NoReaction"); }
  /* 
  G4cout << "The Particle Type was : " << newHit->GetParticleName() << G4endl;
  G4cout << "The Hit occured at " << newHit->GetPos() << G4endl;
  G4cout << "The Energy Deposited at the Hit was " << edep/MeV << " MeV" << G4endl;
  */
  
  TntHitsCollection->insert( newHit );

  return true;
}

void TntSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  // Note "EndOfEvent()" is name, NOT EndofEvent !

 G4SDManager* SDMan = G4SDManager::GetSDMpointer();
 G4String detName = "TntSD";
 G4int collectionID = SDMan->GetCollectionID("TntSD/TntHitsCollection");

 TntDetHitsCollection* myCollection = (TntDetHitsCollection*)(HCE->GetHC(collectionID));

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

 G4double GlobalTime = 0.;
 G4double PulseTime = 0.;  // Time after first hit
 G4ThreeVector FirstHitPos(0.,0.,0.);

 G4int number_of_gammahits = 0;
 G4int NumberOfTracks = 0;
 G4int theLastTrackID = 0;

 if(myCollection)
  {
    G4int n_hit = myCollection->entries();

    // Event readout loop!

    for(int i=0;i<n_hit;i++)
    {
       
      // (*myCollection)[i] is the pointer to the i-th hit of the event.
      // All data analysis output is read out here and then sent to the DataTree!

      TntDetHit* theCurrentHit = (*myCollection)[i];
     
      G4String theParticleName = theCurrentHit->GetParticleName();
      G4double edep = theCurrentHit->GetEdep();
      G4int theParentTrackID = theCurrentHit->GetParentTrackID();
      G4int theTrackID = theCurrentHit->GetTrackID();
      G4int theParticleA = theCurrentHit->GetParticleA();

      if(theTrackID != theLastTrackID && theParticleA >=1)
	{
	  NumberOfTracks++;
	  theLastTrackID = theTrackID;
	}

      G4String theCurrentProcess = theCurrentHit->GetParticleProcess();
      
      if(i == 0)
	{
	  /*   Record Time of "first hit"!  */
	  GlobalTime = (theCurrentHit->GetTOF())/ns;
          FirstHitPos = theCurrentHit->GetPos();
	}    

      PulseTime = (theCurrentHit->GetTOF())/ns;

      // 3/27/2009 - BTR
      // Need to Convert to Light step by step to get proper response function
      // 
      // Should I add a condition like Pozzi et al. here to say only take
      // first 3 hadron tracks?
      //if(NumberOfTracks <=3)
      //	
	  if(theParticleName=="proton")
	    {EsumProton += ConvertToLight("proton",1,edep,Light_Conv);}
	  else if(theParticleName=="deuteron")
	    {EsumDeuteron += ConvertToLight("deuteron",1,edep,Light_Conv);}
	  else if(theParticleName=="triton")
	    {EsumTriton += ConvertToLight("triton",1,edep,Light_Conv);}
	  else if(theParticleName=="He3")
	    {EsumHe3 += ConvertToLight("He3",2,edep,Light_Conv);}
	  else if(theParticleName=="alpha")
	    {EsumAlpha += ConvertToLight("alpha",2,edep,Light_Conv);}
	  else if(theParticleName=="C12[0.0]")
	    {EsumC12 += ConvertToLight("C12",6,edep,Light_Conv);}
	  else if(theParticleName=="C13[0.0]")
	    {EsumC13 += ConvertToLight("C13",6,edep,Light_Conv);}
	  else if(theParticleName=="e-" || 
		  theParticleName=="e+" ||
		  theParticleName=="gamma")
	    {
	      EsumEorGamma += ConvertToLight("e-",-1,edep,Light_Conv);
	      number_of_gammahits++;
	    }
	  else
	    {
	      EsumExotic += edep;  
	      //G4cout << "Exotic Particle Created!!!!! >>> ID = "<< theParticleName << G4endl;
	      //G4cout << "The Energy Deposited by the Exotic Particle : " << edep << G4endl; 
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
   }
 else
   {
     //G4cout << "totE = 0" << endl;
     totE = 0.; // event gated out
   }

 // Send data from hit collection to DataRecordTree and FillTree and text files!
 TntDataOutEV->senddataEV(1,totE);
 TntDataOutEV->senddataEV(2,EsumProton);
 TntDataOutEV->senddataEV(3,EsumAlpha);
 TntDataOutEV->senddataEV(4,EsumC12);
 TntDataOutEV->senddataEV(5,EsumEorGamma);
 TntDataOutEV->senddataEV(6,EsumExotic);

 TntDataOutEV->senddataTOF(GlobalTime);
 TntDataOutEV->senddataPosition(FirstHitPos);

 TntDataOutEV->FillTree();

 // G4cout << ">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>End of Event Analysis!" << G4endl;

} // End of EndOfEvent() Method


G4double TntSD::ConvertToLight(G4String theName, G4double theCharge, G4double edep,G4String Light_Conv)
{
  //Converts Energy deposited in hit using Cecil's formula.

  G4double theLightEng = 0.;

  if(Light_Conv == "none")
    {
      theLightEng = edep;
    }
  else if(Light_Conv == "light-conv" || Light_Conv == "light+resol")
    {     
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
      G4cerr << "Bad choice for light conversion method in TntSD::ConvertToLight() method." << G4endl;
      G4cerr << "Proper Values include \"none\", \"light-conv\", or \"light+resol\"." << G4endl;
	G4ExceptionSeverity severity=FatalException;
	G4Exception("Exception thrown -- Program stopped in TntSD::ConvertToLight() method.", "Error",severity,"Error!");
	//Note: the enum G4ExceptionSeverity has been defined since Geant4 release 5.0 and its values are: FatalException, FatalErrorInArgument, EventMustBeAborted, FatalException, and JustWarning. 
      //G4Exception("Exception thrown -- Program stopped in TntSD::ConvertToLight() method.");
    } 

 
    
  return theLightEng;
}
