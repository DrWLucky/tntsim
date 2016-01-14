/// \file Run.hh
/// \brief defines custom G4run class
#ifndef TEXAN_RUN_HEADER_FILE_12345678
#define TEXAN_RUN_HEADER_FILE_12345678
#include "G4Run.hh"


class TTexan;
class TClonesArray;

namespace texansim {

class VPersistenceManager;
class RunAction;

/// Custom Run class, we want to override the RecordEvent() function for anaysis
class Run : public G4Run
{
public:
	/// Ctor
	Run();
	/// Dtor
	~Run();
	/// Override RecordEvent() to perform analysis functions
	virtual void RecordEvent(const G4Event* event);

private:
	Run(const Run&):G4Run() { }
	Run& operator=(const Run&) { return *this; }

private:
	/// For output
	VPersistenceManager* fPersistence;
	/// Total number of hits
	G4int fNumHits;
	/// Sum of deposited energy
	G4double fEdep;
	/// Array of all hits in the event
	TClonesArray* fHitArray;
	///
	G4int fNumPrimaryHits;
	///
	TClonesArray* fPrimaryHits;
	/// "Experimental" analysis class
	TTexan* fTexan;
};

}



#endif
