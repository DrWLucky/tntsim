/// \file GenericHit.cc
/// \brief Implements class to implement generic hit data-storage functions.
#include "texansim/GenericHit.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"


namespace txs = texansim;


G4ThreeVector txs::GenericHit::GetPosition() const
{
	return G4ThreeVector(fData.fX * cm, fData.fY * cm, fData.fZ *cm);
}

void txs::GenericHit::SetPosition(const G4ThreeVector& pos)
{
	fData.fX = pos.getX() / cm;
	fData.fY = pos.getY() / cm;
	fData.fZ = pos.getZ() / cm;

	fData.fRho   = pos.getRho()   / cm;
	fData.fPhi   = pos.getPhi()   / deg;
	fData.fTheta = pos.getTheta() / deg;

	fData.fR = pos.getR() / cm;
}
