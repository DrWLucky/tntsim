/// \file ArraySD.cc
/// \brief Implementation of the ArraySD class
///
#include <cassert>

#include "texansim/ArraySD.hh"

#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"



namespace txs = texansim;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::ArraySD::ArraySD(const G4String& name, const G4String& hitsCollectionName)
  : G4VSensitiveDetector(name),
		fHitsCollection(NULL)
{
	collectionName.insert(hitsCollectionName);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::ArraySD::~ArraySD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::ArraySD::Initialize(G4HCofThisEvent* hce)
{
  /// - Create hits collection
  fHitsCollection 
    = new ArrayHitsCollection(SensitiveDetectorName, collectionName[0]); 

  /// - Add this collection in hce
  G4int hcID 
    = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hce->AddHitsCollection( hcID, fHitsCollection ); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool txs::ArraySD::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	/// - Read deposited energy, abort if 0
	G4double edep = aStep->GetTotalEnergyDeposit();
	if(edep == 0)
		return false;

	/// - Create new hit object, fill with data from detector, add to hits collection
	ArrayHit* newHit = new ArrayHit();
  newHit->SetTrackID  (aStep->GetTrack()->GetTrackID());
  newHit->SetChamberNb(aStep->GetPreStepPoint()->GetTouchableHandle()
											 ->GetCopyNumber());
  newHit->SetEdep(edep);
  newHit->SetPos (aStep->GetPostStepPoint()->GetPosition());

	newHit->fMass   = aStep->GetTrack()->GetDynamicParticle()->GetMass();
	newHit->fCharge = aStep->GetTrack()->GetDynamicParticle()->GetCharge();

  fHitsCollection->insert( newHit );
  newHit->Print();

  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::ArraySD::EndOfEvent(G4HCofThisEvent*)
{
  if ( verboseLevel>1 ) { 
		G4int nofHits = fHitsCollection->entries();
		G4cout << G4endl
					 << "-------->Hits Collection: in this event they are " << nofHits 
					 << " hits in the tracker chambers: " << G4endl;
		for ( G4int i=0; i<nofHits; i++ ) (*fHitsCollection)[i]->Print();
  }
}
