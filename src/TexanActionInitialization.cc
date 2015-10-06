/// \file TexanActionInitialization.cc
/// \brief Implementation of the TexanActionInitialization class

#include "TexanActionInitialization.hh"
#include "TexanPrimaryGeneratorAction.hh"
#include "TexanRunAction.hh"
#include "TexanEventAction.hh"


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
  SetUserAction(new txs::PrimaryGeneratorAction);
  SetUserAction(new txs::RunAction);
  
  txs::EventAction* eventAction = new txs::EventAction;
  SetUserAction(eventAction);
}  
