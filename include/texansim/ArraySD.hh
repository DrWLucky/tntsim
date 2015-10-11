/// \file ArraySD.hh
/// \brief Definition of the TexanArraySD class
///
#ifndef TEXAN_ARRAY_SD_HH
#define TEXAN_ARRAY_SD_HH 1

#include "G4VSensitiveDetector.hh"

#include "texansim/ArrayHit.hh"

class G4Step;

namespace texansim {

/// Sensitive detector to be attached to the GDML geometry
class ArraySD : public G4VSensitiveDetector
{
public:
	/// Ctor
	ArraySD(const G4String& name, const G4String& hitsCollectionName);
	/// Dtor
	~ArraySD();

	/// Initialization actions
	virtual void Initialize(G4HCofThisEvent*);
	/// Process a hit
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	/// Executed at end of event
	virtual void EndOfEvent(G4HCofThisEvent*);

private:
	texansim::ArrayHitsCollection* fHitsCollection;
};

}

#endif

