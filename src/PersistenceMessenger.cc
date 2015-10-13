/// \file PersistenceMessenger.cc
/// \brief Implements PersistenceManager class
///
#include "texansim/PersistenceMessenger.hh"
#include "texansim/VPersistenceManager.hh"
#include "texansim/Utils.hh"

#include "G4UIcmdWithAString.hh"


namespace txs = texansim;


txs::PersistenceMessenger::PersistenceMessenger(VPersistenceManager* manager):
	G4UImessenger(),
	fManager(manager),
	fDirectory(0),
	fSetFilenameCmd(0)
{
	fDirectory = new G4UIdirectory("/persistence/");
	fDirectory->SetGuidance("Set persistent event output parameters");
      
	fSetFilenameCmd 
		= new G4UIcmdWithAString("/placement/setFilename", this);
	fSetFilenameCmd->SetGuidance("Change persistent event output filename");
	fSetFilenameCmd->SetParameterName("Filename", false);
}


txs::PersistenceMessenger::~PersistenceMessenger()
{
	Zap(fDirectory);
	Zap(fSetFilenameCmd);
}


void txs::PersistenceMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
	if(command == fSetFilenameCmd) { 
		G4cout << "Setting persistence output file name to \"" << newValue << "\"" << G4endl;
		fManager->SetFilename(newValue);
	} else {
		G4cerr << "Ignoring invalid UI command \"" << command->GetCommandPath() << "\"" << G4endl;
	}
}
