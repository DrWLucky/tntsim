/// \file VPersistenceManager.cc
/// \brief Implements VPersistence manager class
///
#include "tntsim/VPersistenceManager.hh"
#include "tntsim/PersistenceMessenger.hh"
#include "tntsim/Utils.hh"

#include <cassert>
#include "G4Threading.hh"



tntsim::VPersistenceManager::VPersistenceManager()
{
	;
}



tntsim::VPersistenceManager::~VPersistenceManager()
{
	;
}


void tntsim::VPersistenceManager::InitializeBase()
{
	/// Takes care of filename stuff:
	/// - Set default filename ("g4output")
	SetFilename("g4output");

	/// - Check for UI message to change file name
	if (gMessenger != NULL) {
		GetMessenger()->ApplyCommands(this);
	}
}


tntsim::DelayedMessenger* tntsim::VPersistenceManager::GetMessenger()
{
	return gMessenger;
}



void tntsim::VPersistenceManager::SetMessenger(DelayedMessenger* messenger)
{
#ifdef G4MULTITHREADED
	assert(G4Threading::IsWorkerThread() == false);
#endif

	gMessenger = messenger;
}


tntsim::DelayedMessenger* tntsim::VPersistenceManager::gMessenger = NULL;
