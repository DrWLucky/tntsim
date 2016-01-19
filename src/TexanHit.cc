/// \file TexanHit.cc
/// \brief Implementation of the TexanHit class

#include "tntsim/TexanHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4Step.hh"
#include <iomanip>

namespace tnt = tntsim;

namespace tntsim {
G4ThreadLocal G4Allocator<tntsim::TexanHit>* TexanHitAllocator = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

tnt::TexanHit::TexanHit()
	: G4VHit(), fStep(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

tnt::TexanHit::TexanHit(G4Step* step)
	: G4VHit()
{
	SetStep(step);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

tnt::TexanHit::~TexanHit()
{ if(fStep) delete fStep; }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int tnt::TexanHit::operator==(const tnt::TexanHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void tnt::TexanHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
		const G4ThreeVector &pos = fStep->GetPostStepPoint()->GetPosition();
    G4Circle circle(pos);
    circle.SetScreenSize(4.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void tnt::TexanHit::Print()
{
	const G4ThreeVector &pos = fStep->GetPostStepPoint()->GetPosition();
  G4cout
		<< "  trackID: " << fStep->GetTrack()->GetTrackID() // << " chamberNb: " << fChamberNb
		<< "Edep: "
		<< std::setw(7) << G4BestUnit(fStep->GetTotalEnergyDeposit(), "Energy") //GetData().fEdep,"Energy")
		<< " Position: "
		<< std::setw(7) << G4BestUnit( pos,"Length")
		<< G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
