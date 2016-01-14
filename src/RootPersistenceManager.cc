/// \file RootPersistenceManager.cc
/// \brief Implements RootPersistenceManager class
#include "tntsim/RootPersistenceManager.hh"
#include "tntsim/Utils.hh"

#include "TH1D.h"
#include "TH1F.h"
#include "TH1I.h"
#include "TH2D.h"
#include "TH2F.h"
#include "TH2I.h"
#include "TFile.h"
#include "TTree.h"
#include "TROOT.h"
#ifdef G4MULTITHREADED
#include "TChain.h"
#endif
#include "TSystem.h"
#include "TObjArray.h"


#include "G4ios.hh"
#include "G4Threading.hh"
#include "G4AutoLock.hh"
#include "G4RunManager.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#endif

namespace txs =  tntsim;


#ifdef G4MULTITHREADED
//Create a global auto lock/unlock device
namespace{ G4Mutex mutex = G4MUTEX_INITIALIZER; }
#define TXS_THREADLOCK G4AutoLock LockLockLock(&mutex)
#else
#define TXS_THREADLOCK do {} while(0)
#endif


namespace {
G4int gNumWorkerFiles = 0;
G4int gNumClosedWorkerFiles = 0;
}


txs::RootPersistenceManager::RootPersistenceManager():
	fFile(0)
{
	TXS_THREADLOCK;
	
	/// Initialize base (chooses the correct file name)
	InitializeBase();

	/// Initialize arrays
	fHists1d = new TObjArray();
	fHists2d = new TObjArray();
}


void txs::RootPersistenceManager::SetFilename(const G4String& name)
{
	/// - Append extension and thread markers
#ifdef G4MULTITHREADED
	fFilename = name;
	if(G4Threading::IsWorkerThread())
		fFilename += G4String(Form("_t%d", G4Threading::G4GetThreadId()));
#else
	fFilename = name;
#endif

	/// - If file is already open, re-name it
	if(fFile) {
		TXS_THREADLOCK;
		G4String newname = GetFilename() + ".root";
		gSystem->Rename(fFile->GetName(), newname.data());
	}
}

G4bool txs::RootPersistenceManager::OpenFile()
{
	TXS_THREADLOCK;

	Zap(fFile);
	G4String name = GetFilename() + ".root";
	fFile = new TFile(name.data(), "recreate");
	if(fFile == 0 || fFile->IsZombie()) {
		Zap(fFile);
		return false;
	}

	fTree = new TTree("t", "Geant4 saved objects");

	for(int i=0; i< fHists1d->GetEntries(); ++i) {
		static_cast<TH1*>(fHists1d->At(i))->SetDirectory(fFile);
	}
	for(int i=0; i< fHists2d->GetEntries(); ++i) {
		static_cast<TH1*>(fHists1d->At(i))->SetDirectory(fFile);
	}

	if(G4Threading::IsWorkerThread())
		++gNumWorkerFiles;

	return true;
}


G4bool txs::RootPersistenceManager::AddObject(const G4String& name, const G4String& classname, void* obj)
{
	TXS_THREADLOCK;
	return fTree ? fTree->Branch(name.data(), classname.data(), obj) : false;
}

G4bool txs::RootPersistenceManager::AddPrimitive(const G4String& name, void* p, Type_t type)
{
	TXS_THREADLOCK;
	if(!fTree)
		return false;
	G4String leaflist;
	switch(type) {
	case kDouble:
		leaflist = name + G4String("/D");
		break;
	case kFloat:
		leaflist = name + G4String("/F");
		break;
	case kInt:
		leaflist = name + G4String("/I");
		break;
	default:
		G4cerr << "Invalid primitive type: " << type << G4endl;
		return false;
		break;
	}

	return fTree->Branch(name.data(), p, leaflist.data());
}

namespace { struct HistDir0 {
	HistDir0()  { TH1::AddDirectory(false); }
	~HistDir0() { TH1::AddDirectory(true ); }
}; }

G4bool txs::RootPersistenceManager::AddHistogram1d(
	const G4String& name, const G4String& title,
	G4int bins, G4double xlow, G4double xhigh,
	void* valuePointer, Type_t type)
{
	TXS_THREADLOCK;
	::HistDir0 hd0;

	TH1* hst = 0;
	switch(type) {
	case kDouble:
		hst = new TH1D(name.data(), title.data(), bins, xlow, xhigh);
		break;
	case kFloat:
		hst = new TH1F(name.data(), title.data(), bins, xlow, xhigh);
		break;
	case kInt:
		hst = new TH1I(name.data(), title.data(), bins, xlow, xhigh);
		break;
	default:
		G4cerr << "Invalid 1d histogram type: " << type << G4endl;
		return false;
		break;
	}
	if(fFile) hst->SetDirectory(fFile);
	fHists1d->Add(hst);
	fHistData1d.push_back(valuePointer);

	return true;
}

G4bool txs::RootPersistenceManager::AddHistogram2d(
	const G4String& name, const G4String& title,
	G4int xbins, G4double xlow, G4double xhigh,
	G4int ybins, G4double ylow, G4double yhigh,
	void* valuePointerx, void* valuePointery,
	Type_t type)
{
	TXS_THREADLOCK;
	::HistDir0 hd0;

	TH2* hst = 0;
	switch(type) {
	case kDouble:
		hst = new TH2D(name.data(), title.data(), xbins, xlow, xhigh, ybins, ylow, yhigh);
		break;																											                    
	case kFloat:																									                    
		hst = new TH2F(name.data(), title.data(), xbins, xlow, xhigh, ybins, ylow, yhigh);
		break;																											                    
	case kInt:
		hst = new TH2I(name.data(), title.data(), xbins, xlow, xhigh, ybins, ylow, yhigh);
		break;
	default:
		G4cerr << "Invalid 2d histogram type: " << type << G4endl;
		return false;
		break;
	}
	if(fFile) hst->SetDirectory(fFile);
	fHists2d->Add(hst);
	fHistData2d.push_back(std::make_pair(valuePointerx, valuePointery));

	return true;
}



void txs::RootPersistenceManager::SaveEvent()
{
	fTree->Fill();

	for(int i=0; i< fHists1d->GetEntries(); ++i) {
		TH1* hst = static_cast<TH1*>(fHists1d->At(i));

		if (0) { }
		else if(hst->IsA() == TH1D::Class()) {
			hst->Fill( *(reinterpret_cast<G4double*>(fHistData1d.at(i))) );
		}
		else if(hst->IsA() == TH1F::Class()) {
			hst->Fill( *(reinterpret_cast<G4float*>(fHistData1d.at(i))) );
		}
		else if(hst->IsA() == TH1I::Class()) {
			hst->Fill( *(reinterpret_cast<G4int*>(fHistData1d.at(i))) );
		}
		else {
			G4cerr << "Wrong histogram type!" << G4endl;
		}
	}

	for(int i=0; i< fHists2d->GetEntries(); ++i) {
		TH2* hst = static_cast<TH2*>(fHists2d->At(i));
		void* px = fHistData2d.at(i).first;
		void* py = fHistData2d.at(i).second;

		if (0) { }
		else if(hst->IsA() == TH2D::Class()) {
			hst->Fill( *(reinterpret_cast<G4double*>(px)),
								 *(reinterpret_cast<G4double*>(py)) );
		}
		else if(hst->IsA() == TH2F::Class()) {
			hst->Fill( *(reinterpret_cast<G4float*>(px)),
								 *(reinterpret_cast<G4float*>(py)) );
		}
		else if(hst->IsA() == TH2I::Class()) {
			hst->Fill( *(reinterpret_cast<G4int*>(px)),
								 *(reinterpret_cast<G4int*>(py)) );
		}
		else {
			G4cerr << "Wrong histogram type!" << G4endl;
		}
	}	
}


void txs::RootPersistenceManager::Write()
{
	TXS_THREADLOCK;
	if(G4Threading::IsMultithreadedApplication()) {
		if(G4Threading::IsWorkerThread() == false) {
			///
			/// Take care of in Merge() if on main thread
			return; 
		}
	}
	fFile->cd();
	if(fTree)
		fTree->AutoSave();
	for(int i=0; i< fHists1d->GetEntries(); ++i) {
		fHists1d->At(i)->Write();
	}
	for(int i=0; i< fHists2d->GetEntries(); ++i) {
		fHists2d->At(i)->Write();
	}
}


txs::RootPersistenceManager::~RootPersistenceManager()
{
	Close();

	TXS_THREADLOCK;
	Zap(fHists1d);
	Zap(fHists2d);
}

#include <cassert>
void txs::RootPersistenceManager::Close()
{
	TXS_THREADLOCK;
	if(fFile) {
#ifdef G4MULTITHREADED
		if(G4Threading::IsMultithreadedApplication()) {
			if(G4Threading::IsWorkerThread() == false) {
				Merge();
			}
		}
#endif
		fFile->Close();
		Zap(fFile);
		++gNumClosedWorkerFiles;
	}
}

void txs::RootPersistenceManager::Merge()
{
#ifdef G4MULTITHREADED
	assert(gNumClosedWorkerFiles == gNumWorkerFiles);

	gROOT->cd();
	TChain* chain = new TChain("t", "Geant4 saved objects (all threads)");

	for(G4int i=0; i< gNumWorkerFiles; ++i) {
		G4String fname = GetFilename() + Form("_t%d.root", i);
		TFile* file = TFile::Open(fname.data());
		if(file == 0)
			continue;

		/// - TChain for event-by-event data
		chain->AddFile(fname.data());

		/// - Then histograms
		for(int j=0; j< fHists1d->GetEntries(); ++j) {
			TH1* hthis = dynamic_cast<TH1*>(fHists1d->At(j));
			TH1* hthat = dynamic_cast<TH1*>(file->Get(hthis->GetName()));
			hthis->Add(hthat);
		}
		for(int j=0; j< fHists2d->GetEntries(); ++j) {
			TH2* hthis = dynamic_cast<TH2*>(fHists2d->At(j));
			TH2* hthat = dynamic_cast<TH2*>(file->Get(hthis->GetName()));
			hthis->Add(hthat);
		}

		file->Delete();
	}

	fFile->cd();
	for(int j=0; j< fHists1d->GetEntries(); ++j) {
		fHists1d->At(j)->Write();
	}
	for(int j=0; j< fHists2d->GetEntries(); ++j) {
		fHists2d->At(j)->Write();
	}
	
	chain->Write();
	chain->Delete();
		
#else
	return;
#endif
}
