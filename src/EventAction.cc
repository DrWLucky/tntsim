/// \file EventAction.cc
/// \brief Implementation of the EventAction class
///
#include "texansim/EventAction.hh"
#include "G4Event.hh"


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

void txs::EventAction::EndOfEventAction(const G4Event*)
{   
	/// Presently empty - analysis stuff done in G4Run::RecordEvent
	/// \todo I think we can get rid of this class entirely...
	;
}
