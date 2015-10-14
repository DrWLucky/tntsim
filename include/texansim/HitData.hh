#ifndef HIT_DATA_HEADER_FILE_INCUDE_GUARD
#define HIT_DATA_HEADER_FILE_INCUDE_GUARD
#include "G4Types.hh"
#include "TObject.h"




namespace texansim {

/// Simple container for hit data
/*! Units are:
 *  - Energy: MeV
 *  - Position: cm
 *  - Angle: degrees
 *  - Time: ns
 */
struct HitData
	: public TObject
{
	G4double fEdep;  /// Energy deposited in the hit (MeV)
	G4double fTime;  /// Global time (since event creation)

	G4double fX;     /// X position
	G4double fY;     /// Y position
	G4double fZ;     /// Z position

	G4double fR;     /// Spherical coordinates radius
	G4double fTheta; /// Spherical coordinates radius polar angle
	G4double fPhi;   /// Spherical coordinates radius azimuthal angle

	G4double fRho;   /// Cylindrical coordinates radius


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
