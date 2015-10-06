/// \file TexanDetectorConstruction.cc
/// \brief Implements detector construction (GDML) class.

#include "TexanDetectorConstruction.hh"


namespace txs = texansim;

txs::DetectorConstruction::DetectorConstruction(G4VPhysicalVolume *setWorld):
	fWorld(setWorld)
{ }


G4VPhysicalVolume* txs::DetectorConstruction::Construct()
{
	return fWorld;
}
