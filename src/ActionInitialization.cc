/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "tntsim/ActionInitialization.hh"
#include "tntsim/PrimaryGeneratorAction.hh"
#include "tntsim/RunAction.hh"


namespace txs = tntsim;

txs::ActionInitialization::ActionInitialization()
	: G4VUserActionInitialization()
{}



txs::ActionInitialization::~ActionInitialization()
{}



void txs::ActionInitialization::BuildForMaster() const
{
  SetUserAction(new txs::RunAction);
}


void txs::ActionInitialization::Build() const
{
  SetUserAction(
		new txs::PrimaryGeneratorAction("neutron", 10*MeV,
																		G4ThreeVector(0,0,0),
																		G4ThreeVector(0,0,1))
		);
  SetUserAction(new txs::RunAction);
}
