/// \file Run.cc
/// \brief Implements texansim::Run class
///
#include "texansim/Run.hh"
#include "texansim/HitData.hh"
#include "texansim/TexanHit.hh"
#include "texansim/RootPersistenceManager.hh"

#include "G4Event.hh"
#include "G4Threading.hh"
#include "G4RunManager.hh"

#include "TClass.h"
#include "TClonesArray.h"

#include "texana/TTexan.h"





texansim::Run::Run():
	G4Run()
{
	fTexan = new TTexan("texan", "");
	fHitArray = new TClonesArray( HitData::Class()->GetName() );

	/// Set up persistence output
	fPersistence = new RootPersistenceManager();
	fPersistence->OpenFile();

	/// - Add histograms
	fPersistence->AddHistogram1d("hEdep", "Total Deposited energy", 200, 0, 20,
															 &fEdep, VPersistenceManager::kDouble);

	/// - Add primitives
	fPersistence->AddPrimitive("fNumHits", &fNumHits,  VPersistenceManager::kInt);
	fPersistence->AddPrimitive("fEdep", &fEdep, VPersistenceManager::kDouble);

	/// - Add classes
	fPersistence->AddObject(fTexan->GetName(), fTexan->ClassName(), &fTexan);
	fPersistence->AddObject("hit", "TClonesArray", &fHitArray);
	fHitArray->BypassStreamer();
}

texansim::Run::~Run()
{
	/// Take care of persistence closing and output: write/close/delete output
	fPersistence->Write();
	fPersistence->Close();
	delete fPersistence;
	delete fTexan;
	delete fHitArray;
}

void texansim::Run::RecordEvent(const G4Event* event)
{
	/// Get hit data and save from event's hit collection, save event to disk.
	fTexan->Reset();

	TexanHitsCollection& hc =
		static_cast<TexanHitsCollection&>(*(event->GetHCofThisEvent()->GetHC(0)));

	fEdep = 0;
	fNumHits = hc.GetSize();
	fHitArray->Clear();

	for(G4int i=0; i< fNumHits; ++i) {
		// Experimental class
		if(i < 32) {
			// fTexan->fParticleMass[i]   = hc[i]->GetHitData().fMass;
			// fTexan->fParticleCharge[i] = hc[i]->fCharge;
			fTexan->SetEcal(i, hc[i]->GetData().fEdep);
		}

		// Hit array
		TObject*& hitMemory = (*fHitArray)[i];
		HitData* hitdata = new(hitMemory) HitData(hc[i]->GetData());
		fEdep += hitdata->fEdep;
	}

	fPersistence->SaveEvent();

	/// Also need to call base class version
	G4Run::RecordEvent(event);
}
