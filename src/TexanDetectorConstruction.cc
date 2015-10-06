/// \file TexanDetectorConstruction.cc
/// \brief Implements detector construction (GDML) class.

#include "TexanDetectorConstruction.hh"


TexanDetectorConstruction::TexanDetectorConstruction(G4VPhysicalVolume *setWorld):
	fWorld(setWorld)
{ }


G4VPhysicalVolume* TexanDetectorConstruction::Construct()
{
	return fWorld;
}
