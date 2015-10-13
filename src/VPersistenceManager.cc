/// \file VPersistenceManager.cc
/// \brief Implements VPersistence manager class
///
#include "texansim/VPersistenceManager.hh"
#include "texansim/PersistenceMessenger.hh"
#include "texansim/Utils.hh"




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
	/// Set default filename ("g4output")
	SetFilename("g4output");

	/// Check for UI messenges and if found, do them
	if (gMessenger != NULL) {
		static_cast<PersistenceMessenger*>(gMessenger)->ApplyCommands(this);
	}
}


G4UImessenger* texansim::VPersistenceManager::GetMessenger()
{
	return gMessenger;
}



void texansim::VPersistenceManager::SetMessenger(G4UImessenger* messenger)
{
	/// \todo assert that not in thread
	gMessenger = messenger;
}


G4UImessenger* texansim::VPersistenceManager::gMessenger = NULL;
