/// \file Run.cc
/// \brief Implements texansim::Run class
///
#include "texansim/Run.hh"
#include "texansim/HitData.hh"
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
	fHitData = new HitData();

	/// Set up persistence output
	fPersistence = new RootPersistenceManager();
	fPersistence->OpenFile();

	/// - Add histograms
	fPersistence->AddHistogram1d("hEdep", "Deposited energy", 20, 0, 20, &fEdep, VPersistenceManager::kDouble);


	/// - Add primitives
	fPersistence->AddPrimitive("numhits", &fNumHits,        VPersistenceManager::kInt);

	/// - Add classes
	fPersistence->AddObject("hit", fHitData->ClassName(), &fHitData);
	fPersistence->AddObject(fTexan->GetName(), fTexan->ClassName(), &fTexan);
}

texansim::Run::~Run()
{
	/// Take care of persistence closing and output: write/close/delete output
	fPersistence->Write();
	fPersistence->Close();
	delete fPersistence;
	delete fTexan;
	delete fHitData;
}

void texansim::Run::RecordEvent(const G4Event* event)
{
	/// Get hit data and save from event's hit collection, save event to disk.
	fTexan->Reset();

	ArrayHitsCollection& hc =
		static_cast<ArrayHitsCollection&>(*(event->GetHCofThisEvent()->GetHC(0)));

	fHitData->Initialize( hc.GetSize() );

	for(int i=0; i< (G4int)hc.GetSize(); ++i) {
		if(i < 32) {
			fTexan->fParticleMass[i]   = hc[i]->fMass;
			fTexan->fParticleCharge[i] = hc[i]->fCharge;
			fTexan->SetEcal(i, hc[i]->GetEdep());
		}

		fHitData->fEdep[i] = hc[i]->GetEdep();
		fHitData->fSum += fHitData->fEdep[i];
	}
		
	fEdep = ((G4int)hc.GetSize() > 0) ? hc[0]->GetEdep() : 0;
	fNumHits = hc.GetSize();

	fPersistence->SaveEvent();

	/// Also need to call base class version
	G4Run::RecordEvent(event);
}
