/// \file Run.hh
/// \brief defines custom G4run class
#ifndef TEXAN_RUN_HEADER_FILE_12345678
#define TEXAN_RUN_HEADER_FILE_12345678
#include "G4Run.hh"


class TTexan;

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
	/// For output
	VPersistenceManager* fPersistence;
	G4double fEdep;

	TTexan* fTexan;
};

}



#endif
