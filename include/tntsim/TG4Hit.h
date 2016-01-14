/// \file TG4Hit.hh
/// \brief Defines ROOT class to save Geant4 hit information
///
#ifndef TG4HIT_HEADER_FILE_INCUDE_GUARD
#define TG4HIT_HEADER_FILE_INCUDE_GUARD
#include <TNamed.h>
#include <TVector3.h>


class TClonesArray;
namespace CLHEP { class Hep3Vector; }


/// Container for hit data
/*! Units are:
 *  - Energy: MeV
 *  - Position: cm
 *  - Angle: degrees
 *  - Time: ns
 */
class TG4Hit
	: public TNamed
{
public:
	TG4Hit();
	TG4Hit(const char* name, const char* title);

	void Reset();
	
	Double_t GetEdep() const { return fEdep; }
	Int_t GetTrackId() const { return fTrackId; }
	Double_t GetTime() const { return fTime; }
	const TVector3& GetPosition() const { return fPosition; }
	Bool_t GetValid() const { return fValid; }

	void SetEdep(Double_t edep) { fEdep = edep; }
	void SetPosition(const TVector3& pos) { fPosition = pos; }
	void SetTime(Double_t time) { fTime = time; }
	void SetTrackId(Int_t id) { fTrackId = id; }
	void SetValid(Bool_t valid = true) { fValid = valid; }
	void SetXYZ(Double_t x, Double_t y, Double_t z)
		{ fPosition.SetXYZ(x, y, z); }

	const CLHEP::Hep3Vector& GetCLHEPPosition() const;
	void SetCLHEPPosition(const CLHEP::Hep3Vector& pos);	

	static TG4Hit* ReadFromArray(TClonesArray& array, Int_t i);
	static TG4Hit* ReadFromArray(TClonesArray* array, Int_t i);

public:
  /// Process by which this hit interacts
	TString fProcess;
	TString fParticle;

private:
	/// Valid flag
	Bool_t fValid;
	/// Track ID
	Int_t fTrackId;
  /// Energy deposited in the hit (MeV)
	Double_t fEdep;
	/// Global time (since event creation)
	Double_t fTime;  
	/// Position of the hit
	TVector3 fPosition;

	ClassDef(TG4Hit, 1);
};





#endif

/* Local Variables:  */
/* mode: c++         */
/* End:              */
