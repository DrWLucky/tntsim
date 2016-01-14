/// \file RootPersistenceManager.hh
/// \brief ROOT implementation of VPersistenceManager
#ifndef TXS_ROOT_PERSISTENCE_HEADER_FILE_8675309
#define TXS_ROOT_PERSISTENCE_HEADER_FILE_8675309
#include "tntsim/VPersistenceManager.hh"
#include <vector>
#include <utility>

class TH1;
class TFile;
class TTree;
class TObjArray;

namespace tntsim {

class RootPersistenceManager : public VPersistenceManager
{
public:
	/// Ctor
	RootPersistenceManager();
	/// Dtor
	virtual ~RootPersistenceManager();
	/// Set output file name (Base - no extension)
	virtual void SetFilename(const G4String& name);
	/// Get output file name
	virtual const G4String& GetFilename() const { return fFilename; }
	/// Open output file
	virtual G4bool OpenFile();
	/// Add object to be saved
	virtual G4bool AddObject(const G4String& name, const G4String& classname, void*);
	/// Add primitive type (double, float, int) to be saved
	virtual G4bool AddPrimitive(const G4String& name, void* p, Type_t type);
	/// Add 1d histogram
	virtual G4bool AddHistogram1d(
		const G4String& name, const G4String& title,
		G4int bins, G4double xlow, G4double xhigh,
		void* valuePointer, Type_t type);
	/// Add 2d histogram
	virtual G4bool AddHistogram2d(
		const G4String& name, const G4String& title,
		G4int xbins, G4double xlow, G4double xhigh,
		G4int ybins, G4double ylow, G4double yhigh,
		void* valuePointerx, void* valuePointery,
		Type_t type);
	/// Save an event of histograms and persistent objects
	virtual void SaveEvent();
	/// Write to disk
	virtual void Write();
	/// Close everything out
	virtual void Close();

protected:
	/// Merge files if multithreaded
	void Merge();

private:
	RootPersistenceManager(const RootPersistenceManager&) { }
	RootPersistenceManager& operator=(const RootPersistenceManager&) { return *this; }

private:
	G4String fFilename;
	TFile* fFile;
	TTree* fTree;
	TObjArray* fHists1d;
	TObjArray* fHists2d;
	std::vector<void*> fHistData1d;
	std::vector<std::pair<void*, void*> > fHistData2d;
};

}




#endif
