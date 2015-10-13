/// \file VPersistenceManager.cc
/// \brief Implements VPersistence manager class
///
#include "texansim/VPersistenceManager.hh"
#include "texansim/PersistenceMessenger.hh"
#include "texansim/Utils.hh"

#include <cassert>
#include "G4Threading.hh"



texansim::VPersistenceManager::VPersistenceManager()
{
	;
}



texansim::VPersistenceManager::~VPersistenceManager()
{
	;
}


void texansim::VPersistenceManager::InitializeBase()
{
	/// Takes care of filename stuff:
	/// - Set default filename ("g4output")
	SetFilename("g4output");

	/// - Check for UI message to change file name
	if (gMessenger != NULL) {
		GetMessenger()->ApplyCommands(this);
	}
}


texansim::DelayedMessenger* texansim::VPersistenceManager::GetMessenger()
{
	return gMessenger;
}



void texansim::VPersistenceManager::SetMessenger(DelayedMessenger* messenger)
{
#ifdef G4MULTITHREADED
	assert(G4Threading::IsWorkerThread() == false);
#endif

	gMessenger = messenger;
}


texansim::DelayedMessenger* texansim::VPersistenceManager::gMessenger = NULL;
