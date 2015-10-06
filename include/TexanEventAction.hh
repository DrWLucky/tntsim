/// \file TexanEventAction.hh
/// \brief Definition of the EventAction class

#ifndef TexanEventAction_h
#define TexanEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"


namespace texansim {

/// Event action class
///
class EventAction : public G4UserEventAction
{
public:
	EventAction();
	virtual ~EventAction();
    
	virtual void BeginOfEventAction(const G4Event* event);
	virtual void EndOfEventAction(const G4Event* event);

	void AddEdep(G4double edep) { fEdep += edep; }

private:
	G4double  fEdep;
};

}


#endif

    
