/// \file SensitiveDetector.hh
/// \brief Definition of the TexanSensitiveDetector class
///
#ifndef TexanSensitiveDetector_h
#define TexanSensitiveDetector_h 1

#include "G4VSensitiveDetector.hh"

class G4Step;

namespace texansim {

/// Sensitive detector to be attached to the GDML geometry
class SensitiveDetector : public G4VSensitiveDetector
{
public:
	SensitiveDetector(const G4String&);
	~SensitiveDetector();

	virtual void Initialize(G4HCofThisEvent*);
	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);
	virtual void EndOfEvent(G4HCofThisEvent*);
};

}

#endif

