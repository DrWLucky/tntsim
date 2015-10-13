/// \file PersistenceMessenger.hh
/// \brief Defines class to pass messages to VPersistenceManager instance
///
#ifndef TXS_PERSISTENCE_MESSENGER_HH_8675309
#define TXS_PERSISTENCE_MESSENGER_HH_8675309
#include "G4UImessenger.hh"


class G4UIcmdWithAString;

namespace texansim {

class VPersistenceManager;


/// Persistence messenger class
class PersistenceMessenger : public G4UImessenger
{
public:
	/// Ctor, associate to specific persistence manager instance
	PersistenceMessenger(VPersistenceManager* manager);
	/// Dtor, empty
	~PersistenceMessenger();
	/// Set new value
	virtual void SetNewValue(G4UIcommand *command, G4String newValue);

private:
	VPersistenceManager* fManager;
	G4UIdirectory* fDirectory;
	G4UIcmdWithAString* fSetFilenameCmd;
};

}


#endif
