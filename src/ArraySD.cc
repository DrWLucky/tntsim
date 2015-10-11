/// \file ArraySD.cc
/// \brief Implementation of the ArraySD class
///
#include <cassert>
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4SystemOfUnits.hh"

#include "texansim/ArraySD.hh"



namespace txs = texansim;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::ArraySD::ArraySD(const G4String& name)
  : G4VSensitiveDetector(name)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::ArraySD::~ArraySD()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::ArraySD::Initialize(G4HCofThisEvent*)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool txs::ArraySD::ProcessHits(G4Step* step, G4TouchableHistory*)
{
	G4double edep = step->GetTotalEnergyDeposit()/MeV;
  G4cout << "Processing hits ....  Energy deposited: " << edep << " MeV" << G4endl;
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::ArraySD::EndOfEvent(G4HCofThisEvent*)
{
}
