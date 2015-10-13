/// \file Run.cc
/// \brief Implements texansim::Run class
///
#include "texansim/Run.hh"
#include "texansim/ArrayHit.hh"
#include "texansim/RootPersistenceManager.hh"

#include "G4Event.hh"
#include "G4Threading.hh"
#include "G4RunManager.hh"
#include "TTree.h"
#include "texana/TTexan.h"



texansim::Run::Run():
	G4Run()
{
	fTexan = new TTexan("texan", "");

	/// Set up persistence output
	fPersistence = new RootPersistenceManager();
	fPersistence->OpenFile();

	/// - Add histograms
	fPersistence->AddHistogram1d("hEdep", "Deposited energy", 20, 0, 20, &fEdep, VPersistenceManager::kDouble);


	/// - Add primitives
	fPersistence->AddPrimitive("edep", &fEdep, VPersistenceManager::kDouble);

	/// - Add classes
	fPersistence->AddObject(fTexan->GetName(), fTexan->ClassName(), &fTexan);
}

texansim::Run::~Run()
{
	/// Take care of persistence closing and output: write/close/delete output
	fPersistence->Write();
	fPersistence->Close();
	delete fPersistence;
	delete fTexan;
}

void texansim::Run::RecordEvent(const G4Event* event)
{
	/// Get hit data and save from event's hit collection, save event to disk.
	fTexan->Reset();

	ArrayHitsCollection& hc =
		static_cast<ArrayHitsCollection&>(*(event->GetHCofThisEvent()->GetHC(0)));

	for(int i=0; i< 32; ++i) {
		if(i >= (G4int)hc.GetSize())
			continue;
		fTexan->SetEcal(i, hc[i]->GetEdep());
	}
		
	fEdep = ((G4int)hc.GetSize() > 0) ? hc[0]->GetEdep() : 0;

	fPersistence->SaveEvent();

	/// Also need to call base class version
	G4Run::RecordEvent(event);
}
