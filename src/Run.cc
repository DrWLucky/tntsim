/// \file Run.cc
/// \brief Implements texansim::Run class
#include "texansim/Run.hh"
#include "G4Event.hh"


texansim::Run::Run():
	G4Run()
{
	/// Empty
	;
}

texansim::Run::~Run()
{
	/// Empty
	;
}



void texansim::Run::RecordEvent(const G4Event* event)
{
	G4cerr << "analyzing event " << event << G4endl;
	G4Run::RecordEvent(event);
}
