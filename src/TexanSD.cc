/// \file TexanSD.cc
/// \brief Implementation of the TexanSD class
///
#include <cassert>

#include "texansim/TexanSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4VProcess.hh"
#include "G4Step.hh"
#include "G4ios.hh"




namespace txs = texansim;






txs::TexanSD::TexanSD(const G4String& name, const G4String& hitsCollectionName)
  : G4VSensitiveDetector(name),
		fHitsCollection(NULL)
{
	collectionName.insert(hitsCollectionName);
}




txs::TexanSD::~TexanSD()
{
}




void txs::TexanSD::Initialize(G4HCofThisEvent* hce)
{
  /// - Create hits collection
  fHitsCollection 
    = new TexanHitsCollection(SensitiveDetectorName, collectionName[0]); 

  /// - Add this collection in hce
  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}





G4bool txs::TexanSD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
#if 0
	if ( aStep->GetTrack()->GetTrackID() != 1 )
		return false;

	/// - Read deposited energy
	G4double edep = aStep->GetTotalEnergyDeposit();

	/// - Create new TexanHit object to fill with data and add to hits collection
	TexanHit* newHit = new TexanHit();

  // newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  // newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
	// 										 ->GetCopyNumber());

	newHit->SetEdep(edep);
	newHit->SetTime(aStep->GetTrack()->GetGlobalTime());
	newHit->SetPosition(aStep->GetPostStepPoint()->GetPosition());

	G4String pname = aStep->GetPostStepPoint()->GetProcessDefinedStep()->GetProcessName();
	newHit->SetProcessName(pname);
	newHit->SetParticleName(aStep->GetTrack()->GetDefinition()->GetParticleName());


	newHit->SetStep(aStep);

  fHitsCollection->insert(newHit);
#endif


	TexanHit* newHit = new TexanHit(aStep);
  fHitsCollection->insert(newHit);

  if(verboseLevel > 1)
		newHit->Print();
  return true;
}




void txs::TexanSD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
		G4int nofHits = fHitsCollection->entries();
		G4cout << G4endl
					 << "-------->Hits Collection: in this event they are " << nofHits 
					 << " hits in the tracker chambers: " << G4endl;
		for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}
