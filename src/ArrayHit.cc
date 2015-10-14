/// \file ArrayHit.cc
/// \brief Implementation of the ArrayHit class

#include "texansim/ArrayHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

namespace txs = texansim;

namespace texansim {
G4ThreadLocal G4Allocator<texansim::ArrayHit>* ArrayHitAllocator = 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::ArrayHit::ArrayHit()
 : G4VHit()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::ArrayHit::~ArrayHit() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// txs::ArrayHit::ArrayHit(const txs::ArrayHit& right)
//   : G4VHit()
// {
//   fTrackID   = right.fTrackID;
//   fChamberNb = right.fChamberNb;
//   fEdep      = right.fEdep;
//   fPos       = right.fPos;
// }

// //....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

// const txs::ArrayHit& txs::ArrayHit::operator=(const txs::ArrayHit& right)
// {
//   fTrackID   = right.fTrackID;
//   fChamberNb = right.fChamberNb;
//   fEdep      = right.fEdep;
//   fPos       = right.fPos;

//   return *this;
// }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4int txs::ArrayHit::operator==(const txs::ArrayHit& right) const
{
  return ( this == &right ) ? 1 : 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::ArrayHit::Draw()
{
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
		G4ThreeVector pos(fData.fX, fData.fY, fData.fZ);
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

void txs::ArrayHit::Print()
{
	G4ThreeVector pos(fData.fX, fData.fY, fData.fZ);
  G4cout
     // << "  trackID: " << fTrackID << " chamberNb: " << fChamberNb
     << "Edep: "
     << std::setw(7) << G4BestUnit(fData.fEdep,"Energy")
     << " Position: "
     << std::setw(7) << G4BestUnit( pos,"Length")
     << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
