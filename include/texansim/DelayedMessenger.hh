/// \file DelayedMessenger.hh
/// \brief Defines delayed messenger class
///
#ifndef TXS_DELAYED_MESSENGER_HEADER_FILE
#define TXS_DELAYED_MESSENGER_HEADER_FILE
#include "G4UImessenger.hh"


namespace texansim {

/// Delayed messenger class
/*! Abstract class for prcessing UI messenges in a delayed manner.
 *  The problem to solve is this: macro files are all read in before
 *  /run/beamOn is called. So if we want to use the macro to do something
 *  with a class constructed _after_ beamOn the messenger still needs to be
 *  initialized _before_ beamOn. What this class is intended for is
 *  initializing the messenger class before beamOn, then saving
 *  a list of actions to later be performed when the time is right.
 */
class DelayedMessenger : public G4UImessenger
{
public:
	/// Ctor, call base only
	DelayedMessenger(): G4UImessenger() { }
	/// Dtor, empty
	virtual ~DelayedMessenger() { }
	/// Set new value
	/*! Defined in derived class, no different from normal
	 * SetNewValue
	 */
	virtual void SetNewValue(G4UIcommand* command, G4String newValue) = 0;
	/// Apply commands
	/*! Tells how to apply commands once it's time to run them.
	 *  \param p Pointer to the object used to apply the commands.
	 */
	virtual void ApplyCommands(void* p) = 0;
};

}

#endif
