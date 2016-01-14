/// \file Run.cc
/// \brief Implements tntsim::Run class
///
#include "tntsim/Run.hh"
#include "tntsim/TexanHit.hh"
#include "tntsim/RootPersistenceManager.hh"
#include "tntsim/TG4Hit.h"

#include <TClonesArray.h>

#include "G4Event.hh"
#include "G4Threading.hh"
#include "G4RunManager.hh"

#include "texana/TTexan.h"





tntsim::Run::Run():
	G4Run()
{
	fTexan = new TTexan("texan", "");
	fHitArray = new TClonesArray("TG4Hit");
	fPrimaryHits = new TClonesArray("TG4Hit");
	fNumPrimaryHits = 0;

	/// Set up persistence output
	fPersistence = new RootPersistenceManager();
	fPersistence->OpenFile();

	/// - Add histograms
	fPersistence->AddHistogram1d("hEdep", "Total Deposited energy", 200, 0, 20,
															 &fEdep, VPersistenceManager::kDouble);

	/// - Add primitives
	fPersistence->AddPrimitive("fNumHits", &fNumHits,  VPersistenceManager::kInt);
	fPersistence->AddPrimitive("fNumPrimaryHits", &fNumPrimaryHits,  VPersistenceManager::kInt);
	fPersistence->AddPrimitive("fEdep", &fEdep, VPersistenceManager::kDouble);
	
	/// - Add classes
	fPersistence->AddObject(fTexan->GetName(), fTexan->ClassName(), &fTexan);
	fPersistence->AddObject("hit", "TClonesArray", &fHitArray);
	fHitArray->BypassStreamer();
	fPersistence->AddObject("primaryhit", "TClonesArray", &fPrimaryHits);
	fPrimaryHits->BypassStreamer();
}

tntsim::Run::~Run()
{
	/// Take care of persistence closing and output: write/close/delete output
	fPersistence->Write();
	fPersistence->Close();
	delete fPersistence;
	delete fTexan;
	delete fHitArray;
	delete fPrimaryHits;
}




void tntsim::Run::RecordEvent(const G4Event* event)
{
	/// Get hit data and save from event's hit collection, save event to disk.
	fTexan->Reset();

	TexanHitsCollection& hc =
		static_cast<TexanHitsCollection&>(*(event->GetHCofThisEvent()->GetHC(0)));

	fEdep = 0;
	fNumHits = hc.GetSize();
	fHitArray->Clear();

	fNumPrimaryHits = 0;
	for(G4int i=0; i< fNumHits; ++i) {
		const G4Step* aStep = hc[i]->GetStep();
		// if ( aStep->GetTrack()->GetTrackID() != 1 )
		// 	continue;

		/// - Read deposited energy
		G4double edep = aStep->GetTotalEnergyDeposit();
		// Hit array
		TG4Hit* g4hit = new( (*fHitArray)[i] ) TG4Hit(Form("Hit%i", i), "");

		g4hit->SetTrackId(aStep->GetTrack()->GetTrackID());
		g4hit->SetEdep(edep);
		g4hit->SetTime(aStep->GetTrack()->GetGlobalTime());
		g4hit->SetCLHEPPosition(aStep->GetPostStepPoint()->GetPosition());
		g4hit->fParticle = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
		g4hit->fProcess  = aStep->GetTrack()->GetParentID() == 0 ?
			"primary" : aStep->GetTrack()->GetCreatorProcess()->GetProcessName();

		// increment total edep
		fEdep += g4hit->GetEdep();

		// Primary hits //
		if(aStep->GetTrack()->GetParentID() == 0) {

			TG4Hit* primaryHit = new( (*fPrimaryHits)[fNumPrimaryHits] ) TG4Hit();
//			primaryHit->SetTitle(primaryHit->GetName());
			primaryHit->SetName(Form("PrimaryHit_%i", fNumPrimaryHits));

			// Get total energy deposited by primaries in this step
			G4double edepTot = 0;
			const std::vector<const G4Track*>* secondaries = 
				aStep->GetSecondaryInCurrentStep();

			for(size_t j=0; j< secondaries->size(); ++j) {
				if(!j) {
					primaryHit->fProcess = secondaries->at(j)->GetCreatorProcess()->GetProcessName();
					primaryHit->fParticle = secondaries->at(j)->GetParticleDefinition()->GetParticleName();
				}
				// else
				// 	if( std::string(secondaries->at(j)->GetCreatorProcess()->GetProcessName().data()) != primaryHit->fProcess )
				// 	{
				// 		G4cerr << "DIFFERENT Process: Current \"" << secondaries->at(j)->GetCreatorProcess()->GetProcessName()
				// 					 << "\", set: \"" << primaryHit->fProcess << "\"" << G4endl;
				// 	}
				edepTot += secondaries->at(j)->GetStep()->GetTotalEnergyDeposit();
			}
			
			primaryHit->SetEdep(edepTot);
			++fNumPrimaryHits;
		}

#if 0
		// Experimental class
		if(i < 32) {
			// fTexan->fParticleMass[i]   = hc[i]->GetHitData().fMass;
			// fTexan->fParticleCharge[i] = hc[i]->fCharge;
			fTexan->SetEcal(i, g4hit->GetEdep());
		}
#endif

	}

	fPersistence->SaveEvent();

	/// Also need to call base class version
	G4Run::RecordEvent(event);
}
