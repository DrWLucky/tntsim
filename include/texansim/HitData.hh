#ifndef HIT_DATA_HEADER_FILE_INCUDE_GUARD
#define HIT_DATA_HEADER_FILE_INCUDE_GUARD
#include "G4Types.hh"
#include "TObject.h"




namespace texansim {

/// Simple container for hit data
struct HitData
	: public TObject
{
	G4double fEdep; /// Energy deposited in the hit (MeV)
	G4double fX;    /// X position (cm)
	G4double fY;    /// Y position (cm)
	G4double fZ;    /// Z position (cm)

	ClassDef(HitData, 1);
};

}


#endif
#if 0
  fTrackID   = right.fTrackID;
  fChamberNb = right.fChamberNb;
  fEdep      = right.fEdep;
  fPos       = right.fPos;
#endif
