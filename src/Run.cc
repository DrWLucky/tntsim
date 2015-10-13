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

	/// - Add histograms
	fPersistence->AddHistogram1d("hEdep", "Deposited energy", 20, 0, 20, &fEdep, VPersistenceManager::kDouble);


	/// - Add primitives
	fPersistence->AddPrimitive("edep", &fEdep, VPersistenceManager::kDouble);

	/// - Add classes
}

texansim::Run::~Run()
{
	/// Take care of persistence closing and output: write/close/delete output
	fPersistence->Write();
	fPersistence->Close();
	delete fPersistence;
}

void texansim::Run::RecordEvent(const G4Event* event)
{
	/// Get hit data and save from event's hit collection, save event to disk.
	ArrayHitsCollection& hc =
		static_cast<ArrayHitsCollection&>(*(event->GetHCofThisEvent()->GetHC(0)));

	fEdep = ((G4int)hc.GetSize() > 0) ? hc[0]->GetEdep() : 0;

	fPersistence->SaveEvent();

	/// Also need to call base class version
	G4Run::RecordEvent(event);
}
