/// \file TexanEventAction.hh
/// \brief Definition of the TexanEventAction class

#ifndef TexanEventAction_h
#define TexanEventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

/// Event action class
///
class TexanEventAction : public G4UserEventAction
{
public:
	TexanEventAction();
	virtual ~TexanEventAction();
    
	virtual void BeginOfEventAction(const G4Event* event);
	virtual void EndOfEventAction(const G4Event* event);

	void AddEdep(G4double edep) { fEdep += edep; }

private:
	G4double  fEdep;
};



#endif

    
