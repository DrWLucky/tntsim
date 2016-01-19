/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "tntsim/ActionInitialization.hh"
#include "tntsim/PrimaryGeneratorAction.hh"
#include "tntsim/RunAction.hh"


namespace tnt = tntsim;

tnt::ActionInitialization::ActionInitialization()
	: G4VUserActionInitialization()
{}



tnt::ActionInitialization::~ActionInitialization()
{}



void tnt::ActionInitialization::BuildForMaster() const
{
  SetUserAction(new tnt::RunAction);
}


void tnt::ActionInitialization::Build() const
{
  SetUserAction(
		new tnt::PrimaryGeneratorAction("neutron", 10*MeV,
																		G4ThreeVector(0,0,0),
																		G4ThreeVector(0,0,1))
		);
  SetUserAction(new tnt::RunAction);
}
