/// \file TexanActionInitialization.cc
/// \brief Implementation of the TexanActionInitialization class

#include "TexanActionInitialization.hh"
#include "TexanPrimaryGeneratorAction.hh"
#include "TexanRunAction.hh"
#include "TexanEventAction.hh"



TexanActionInitialization::TexanActionInitialization()
 : G4VUserActionInitialization()
{}



TexanActionInitialization::~TexanActionInitialization()
{}



void TexanActionInitialization::BuildForMaster() const
{
  SetUserAction(new TexanRunAction);
}


void TexanActionInitialization::Build() const
{
  SetUserAction(new TexanPrimaryGeneratorAction);
  SetUserAction(new TexanRunAction);
  
  TexanEventAction* eventAction = new TexanEventAction;
  SetUserAction(eventAction);
}  
