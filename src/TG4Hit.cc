#include <TClonesArray.h>
#include "G4ThreeVector.hh"
#include "tntsim/TG4Hit.h"






TG4Hit::TG4Hit():
	TNamed("g4hit", "GEANT4 hit")
{
	Reset();
}


TG4Hit::TG4Hit(const char* name, const char* title):
	TNamed(name, title)
{
	Reset();
}

void TG4Hit::Reset()
{
	fValid = false;
	SetXYZ(0, 0, 0);
	SetEdep(0);
	SetTime(0);
}


const CLHEP::Hep3Vector& TG4Hit::GetCLHEPPosition() const
{
	static CLHEP::Hep3Vector v;
	v.set(fPosition.x(), fPosition.y(), fPosition.z());
	return v;
}


void TG4Hit::SetCLHEPPosition(const CLHEP::Hep3Vector& pos)
{
	fPosition.SetXYZ(pos.x(), pos.y(), pos.z());
}


TG4Hit* TG4Hit::ReadFromArray(TClonesArray& array, Int_t i)
{
	return i < array.GetEntries() ?
		static_cast<TG4Hit*>(array[i]) : 0;
}


TG4Hit* TG4Hit::ReadFromArray(TClonesArray* array, Int_t i)
{
	return ReadFromArray(*array, i);
}
