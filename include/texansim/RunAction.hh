/// \file TexanRunAction.hh
/// \brief Definition of the TexanRunAction class

#ifndef TexanRunAction_h
#define TexanRunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class G4LogicalVolume;

// class TTree;
// class TFile;

namespace texansim {

/// Run action class
///
/// In EndOfRunAction(), it calculates the dose in the selected volume 
/// from the energy deposit accumulated via stepping and event actions.
/// The computed dose is then printed on the screen.
class RunAction : public G4UserRunAction
{
public:
	RunAction();
	virtual ~RunAction();

	virtual G4Run* GenerateRun();
	virtual void BeginOfRunAction(const G4Run*);
	virtual void EndOfRunAction(const G4Run*);

// private:
// 	TFile* fFile;
// 	TTree* fTree;
};

}

#endif

