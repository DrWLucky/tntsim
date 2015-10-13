/// \file Run.cc
/// \brief Implements texansim::Run class
#include "texansim/Run.hh"
#include "texansim/ArrayHit.hh"
#include "texansim/RootPersistenceManager.hh"
#include "G4Event.hh"
#include "G4Threading.hh"
#include "G4RunManager.hh"


texansim::Run::Run():
	G4Run()
{
	/// Set up persistence output
	fPersistence = new RootPersistenceManager();
	fPersistence->OpenFile();

	fPersistence->AddHistogram1d("hEdep", "Deposited energy", 20, 0, 20, &fEdep, VPersistenceManager::kDouble);
	fPersistence->AddPrimitive("edep", &fEdep, VPersistenceManager::kDouble);
}

texansim::Run::~Run()
{
	/// Take care of persistence closing and output
	fPersistence->Write();
	fPersistence->Close();
	delete fPersistence;
}

void texansim::Run::RecordEvent(const G4Event* event)
{
	/// Get event data and save, also call G4Run version.
	ArrayHitsCollection& hc =
		static_cast<ArrayHitsCollection&>(*(event->GetHCofThisEvent()->GetHC(0)));

	fEdep = ((G4int)hc.GetSize() > 0) ? hc[0]->GetEdep() : 0;

	fPersistence->SaveEvent();
	G4Run::RecordEvent(event);
}
