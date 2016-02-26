//==========================================================================
// TntDetHit.cc
// Based on DemonScintHit.cc and ExN04TrackerHit.cc
//
// Written/Modified by: Brian Roeder, LPC Caen 02/14/07
//                      email - roeder@lpccaen.in2p3.fr
//
// Purpose: Defines TntDetHit class for simulation data readout
//          Sets and accesses Data in TntDetHit.hh methods
//
//==========================================================================
//
// - See UM Hits Presentation and JLab Scoring 2 Talk for more info.
//

#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include "TntDetHit.hh"

G4Allocator<TntDetHit> TntDetHitAllocator;


TntDetHit::TntDetHit() 
{ /* Constructor */ }


TntDetHit::~TntDetHit() 
{ /* Destructor */ }


// The below taken from ExNO4TrackerHit.cc, but is not needed!
/*
TntDetHit::TntDetHit(const TntDetHit &right) : G4VHit()
{
  G4cout << ">>>>>>>>>>>>>>>>>>>>>A Hit was seen" << G4endl;
  edep = right.edep;
  pos = right.pos;
}

const TntDetHit& TntDetHit::operator=(const TntDetHit &right)
{
  edep = right.edep;
  pos = right.pos;
  return *this;
}

G4int TntDetHit::operator==(const TntDetHit &right) const
{
  return (this==&right) ? 1 : 0;
}
*/

// These Methods were taken directly from ExN04 and Demon

void TntDetHit::Draw()
{
  
  G4VVisManager* pVVisManager = G4VVisManager::GetConcreteInstance();
  if(pVVisManager)
  {
    G4Circle circle(pos);
    circle.SetScreenSize(10.);
    circle.SetFillStyle(G4Circle::filled);
    G4Colour colour(1.,0.,0.);
    G4VisAttributes attribs(colour);
    circle.SetVisAttributes(attribs);
    pVVisManager->Draw(circle);
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetHit::Print()
{ /* Empty for now */ }
