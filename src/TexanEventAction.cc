/// \file TexanEventAction.cc
/// \brief Implementation of the TexanEventAction class
///
#include "TexanEventAction.hh"
#include "TexanAnalysis.hh"
// #include "TexanRun.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "Randomize.hh"


namespace txs = texansim;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::EventAction::EventAction()
: G4UserEventAction(),
  fEdep(0.)
{} 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::EventAction::BeginOfEventAction(const G4Event*)
{    
  fEdep = G4MTRandFlat::shoot(0., 10.);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::EventAction::EndOfEventAction(const G4Event*)
{   
	Analysis::FillH1("hval1", fEdep);
	Analysis::FillNtupleColumn("val1", fEdep);
	Analysis::AddNtupleRow();

  // accumulate statistics in B1Run
  // B1Run* run 
  //   = static_cast<B1Run*>(
  //       G4RunManager::GetRunManager()->GetNonConstCurrentRun());
  // run->AddEdep(fEdep);

	
}
