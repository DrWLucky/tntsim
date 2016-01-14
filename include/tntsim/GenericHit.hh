/// \file GenericHit.hh
/// \brief Defines class to implement generic hit data-storage functions.
#ifndef TXS_GENERIC_HIT_HEADER
#define TXS_GENERIC_HIT_HEADER
#include "tntsim/HitData.hh"
#include "tntsim/Units.hh"
#include "G4ThreeVector.hh"



namespace tntsim {

/// Implements data storage functions common to any type of hit.
class GenericHit {
public:
	/// Empty
	GenericHit() { }
	/// Empty
	virtual ~GenericHit() { }

	// 'Setter methods'
	/// Get read-only access to all hit data
	const HitData& GetData() const { return fData; }
	/// Get vector of hit position
	G4ThreeVector GetPosition() const;

	// 'Setter methods'
	/// Set position vector
	void SetPosition(const G4ThreeVector& pos);
	/// Set deposited energy
	void SetEdep(G4double edep) { fData.fEdep = edep / MeV; }
	/// Set time
	void SetTime(G4double time) { fData.fTime = time / ns; }
	/// Set process name
	void SetProcessName(const char* name) { fData.fProcessName = name; }
	/// Set particle name
	void SetParticleName(const char* name) { fData.fParticleName = name; }

private:
	/// Contains hit data
	HitData fData;
};

}



#endif
