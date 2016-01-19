/// \file PersistenceMessenger.cc
/// \brief Implements PersistenceManager class
///
#include "tntsim/PersistenceMessenger.hh"
#include "tntsim/VPersistenceManager.hh"
#include "tntsim/Utils.hh"

#include "G4UIcmdWithAString.hh"

#include <algorithm>


namespace tnt = tntsim;




tnt::PersistenceMessenger::PersistenceMessenger():
	DelayedMessenger(),
	fDirectory(0),
	fSetFilenameCmd(0)
{
	fDirectory = new G4UIdirectory("/persistence/");
	fDirectory->SetGuidance("Set persistent event output parameters");
      
	fSetFilenameCmd 
		= new G4UIcmdWithAString("/persistence/setFilename", this);
	fSetFilenameCmd->SetGuidance("Change persistent event output filename");
	fSetFilenameCmd->SetParameterName("Filename", false);
}


tnt::PersistenceMessenger::~PersistenceMessenger()
{
	Zap(fDirectory);
	Zap(fSetFilenameCmd);
}


void tnt::PersistenceMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
	if ( command == fSetFilenameCmd ) {
		fCommands.insert(std::make_pair(command, newValue));
	} else {
		G4cerr << "Ignoring invalid UI command \"" << command->GetCommandPath() << "\"" << G4endl;
	}
}



namespace { struct DoApplyCommand {
	tnt::VPersistenceManager* fManager;
	DoApplyCommand(tnt::VPersistenceManager* manager):
		fManager(manager) { }
	void operator() (const std::pair<G4UIcommand*, G4String>& element)
		{
			const G4String& newValue = element.second;
			G4cout << "Setting new persistence output file name to \"" << newValue << "\"" << G4endl;
			fManager->SetFilename(newValue);		
		}
}; }

void tnt::PersistenceMessenger::ApplyCommands(void* managerAddr)
{
	VPersistenceManager* manager =
		reinterpret_cast<VPersistenceManager*>(managerAddr);
	DoApplyCommand applyCommand(manager);
	std::for_each(fCommands.begin(), fCommands.end(), applyCommand);
}
