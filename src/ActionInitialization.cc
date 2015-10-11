/// \file ActionInitialization.cc
/// \brief Implementation of the ActionInitialization class

#include "texansim/ActionInitialization.hh"
#include "texansim/PrimaryGeneratorAction.hh"
#include "texansim/RunAction.hh"
#include "texansim/EventAction.hh"


namespace txs = texansim;

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
  
  txs::EventAction* eventAction = new txs::EventAction;
  SetUserAction(eventAction);
}
