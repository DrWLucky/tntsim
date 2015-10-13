/// \file VPersistenceManager.cc
/// \brief Implements VPersistence manager class
///
#include "texansim/VPersistenceManager.hh"
#include "texansim/PersistenceMessenger.hh"
#include "texansim/Utils.hh"




texansim::VPersistenceManager::VPersistenceManager()
{
	fMessenger = new PersistenceMessenger(this);
}



texansim::VPersistenceManager::~VPersistenceManager()
{
	Zap(fMessenger);
}



G4UImessenger* texansim::VPersistenceManager::GetMessenger() const
{
	return fMessenger;
}



void texansim::VPersistenceManager::SetMessenger(G4UImessenger* messenger)
{
	ResetPointer(fMessenger, messenger);
}



// G4String& texansim::VPersistenceManager::GetFilenameBase()
// {
// 	static G4String* fname = 0;
// 	if(fname == 0)
// 		fname = new G4String();
// 	return *fname;
// }
