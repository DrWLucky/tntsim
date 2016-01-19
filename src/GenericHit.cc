/// \file GenericHit.cc
/// \brief Implements class to implement generic hit data-storage functions.
#include "tntsim/GenericHit.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"


namespace tnt = tntsim;


G4ThreeVector tnt::GenericHit::GetPosition() const
{
	return G4ThreeVector(fData.fX * cm, fData.fY * cm, fData.fZ *cm);
}

void tnt::GenericHit::SetPosition(const G4ThreeVector& pos)
{
	fData.fX = pos.getX() / cm;
	fData.fY = pos.getY() / cm;
	fData.fZ = pos.getZ() / cm;

	fData.fRho   = pos.getRho()   / cm;
	fData.fPhi   = pos.getPhi()   / deg;
	fData.fTheta = pos.getTheta() / deg;

	fData.fR = pos.getR() / cm;
}
