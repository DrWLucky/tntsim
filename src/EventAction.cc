/// \file EventAction.cc
/// \brief Implementation of the EventAction class
///
#include "texansim/EventAction.hh"
#include "texansim/Analysis.hh"
#include "texansim/Utils.hh"
#include "texansim/ArrayHit.hh"
// #include "texansim/Run.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
// #include "Randomize.hh"


namespace txs = texansim;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::EventAction::EventAction()
: G4UserEventAction()
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::EventAction::BeginOfEventAction(const G4Event*)
{    
	;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::EventAction::EndOfEventAction(const G4Event* event)
{   
	G4VHitsCollection* hc = event->GetHCofThisEvent()->GetHC(0);
	for(G4int i=0; i< TXS_MAX_HITS; ++i) {

		G4String colname = FormatStr1("fEdep", i);	
		G4double edep = (i < (G4int)hc->GetSize()) ?
			dynamic_cast<ArrayHit*>(hc->GetHit(i))->GetEdep() : 0;

		Analysis::FillNtupleColumn(colname, edep);
	}

	Analysis::FillNtupleColumn("fNumHits", (G4int)hc->GetSize());
	Analysis::AddNtupleRow();
}
