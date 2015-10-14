/// \file TexanSD.hh
/// \brief Definition of the TexanTexanSD class
///
#ifndef TEXAN_ARRAY_SD_HH
#define TEXAN_ARRAY_SD_HH 1

#include "G4VSensitiveDetector.hh"
#include "texansim/TexanHit.hh"



class G4Step;

namespace texansim {

/// Sensitive detector for TEXAN array.
/** Sensitive detector to be attached to the GDML geometry
 * describing the neutron array TEXAN.
 */
class TexanSD
	: public G4VSensitiveDetector
{
public:
	/// Ctor
	TexanSD(const G4String& name, const G4String& hitsCollectionName);
	/// Dtor
	~TexanSD();

	/// Initialization actions
	virtual void Initialize(G4HCofThisEvent*);
	/// Process a hit
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	/// Executed at end of event
	virtual void EndOfEvent(G4HCofThisEvent*);

private:
	/// Collection of hits
	texansim::TexanHitsCollection* fHitsCollection;
};

}

#endif

