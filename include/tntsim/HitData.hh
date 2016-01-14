#ifndef HIT_DATA_HEADER_FILE_INCUDE_GUARD
#define HIT_DATA_HEADER_FILE_INCUDE_GUARD
#include "G4Types.hh"
#include "TClonesArray.h"
#include <string>

#include "TVector3.h"
namespace CLHEP { class Hep3Vector; }

namespace tntsim {

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
	G4double fEdep;  ///< Energy deposited in the hit (MeV)
	G4double fTime;  ///< Global time (since event creation)

	G4double fX;     ///< X position
	G4double fY;     ///< Y position
	G4double fZ;     ///< Z position

	G4double fR;     ///< Spherical coordinates radius
	G4double fTheta; ///< Spherical coordinates radius polar angle
	G4double fPhi;   ///< Spherical coordinates radius azimuthal angle
	G4double fRho;   ///< Cylindrical coordinates radius

	std::string fProcessName; ///< Name of the process following the hit
	std::string fParticleName; ///< Name of the track particle

	TVector3 fPosition;

	ClassDef(HitData, 1);
};

/// Read hit from a TClonesArray reference
inline HitData* ReadHit(TClonesArray& hitArray, G4int i)
{ return i < hitArray.GetEntries() ? static_cast<HitData*>(hitArray[i]) : 0; }

/// Read hit from a TClonesArray pointer
inline HitData* ReadHit(TClonesArray* hitArray, G4int i)
{ return ReadHit(*hitArray, i); }

class ThreeVector : public TVector3
{
public:
	ThreeVector(double x_=0, double y_=0, double z_=0):
		TVector3(x_,y_,z_) { }
	ThreeVector& operator=(const CLHEP::Hep3Vector& rhs);
	CLHEP::Hep3Vector CopyAsCLHEP();
	ClassDef(ThreeVector,1);
};



}


#ifndef __MAKECINT__
#include "G4ThreeVector.hh"

inline tntsim::ThreeVector&
tntsim::ThreeVector::operator=(const CLHEP::Hep3Vector& rhs)
{
	SetXYZ(rhs.x(), rhs.y(), rhs.z());
	return *this;
}
inline CLHEP::Hep3Vector tntsim::ThreeVector::CopyAsCLHEP()
{
	return CLHEP::Hep3Vector(X(), Y(), Z());
}

#endif



#endif
