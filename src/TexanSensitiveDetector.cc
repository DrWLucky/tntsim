/// \file TexanSensitiveDetector.cc
/// \brief Implementation of the SensitiveDetector class
///
#include <cassert>
#include "TexanSensitiveDetector.hh"
#include "G4HCofThisEvent.hh"
#include "G4Step.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"


namespace txs = texansim;



//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::SensitiveDetector::SensitiveDetector(const G4String& name)
  : G4VSensitiveDetector(name)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::SensitiveDetector::~SensitiveDetector()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::SensitiveDetector::Initialize(G4HCofThisEvent*)
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4bool txs::SensitiveDetector::ProcessHits(G4Step*, G4TouchableHistory*)
{
	assert(0 && "SensitiveDetector");
  G4cout << "Processing hits ...." << G4endl; 
  return true;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::SensitiveDetector::EndOfEvent(G4HCofThisEvent*)
{
}
