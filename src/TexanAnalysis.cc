#include "G4UImanager.hh"
#include "G4RootAnalysisManager.hh"
#include "TexanAnalysis.hh"


namespace txs = texansim;


txs::AnalysisManager*& txs::AnalysisManager::Instance()
{
	static txs::AnalysisManager* ana = 0;
	if(!ana) {
		ana = new txs::AnalysisManager();
	}
	return ana;
}

G4VAnalysisManager* txs::AnalysisManager::G4()
{
	static G4VAnalysisManager* ana = 0;
	if(!ana) {
		// Read from .mac file to decide analyzer type //
		ana = new G4RootAnalysisManager();
	}
	return ana;
}

txs::AnalysisManager::AnalysisManager()
{ 
}

txs::AnalysisManager::~AnalysisManager()
{
	/// \warning Do not call delete G4VRootAnalysisManager anywhere else!
	/// (same for other G4VAnalysisManager classes)
	delete G4();
}


G4int txs::AnalysisManager::BookH1(const G4String &name, const G4String &title, G4int nbins, G4double xmin, G4double xmax, const G4String &unitName, const G4String &fcnName)
{
	/// \warning All names must be unique
	G4int id = G4()->CreateH1(name, title, nbins, xmin, xmax, unitName, fcnName);
	fMap[name] = id;
	return id;
}

G4int txs::AnalysisManager::BookH2(const G4String &name, const G4String &title, G4int nxbins, G4double xmin, G4double xmax, G4int nybins, G4double ymin, G4double ymax,  const G4String &unitName, const G4String &fcnName)
{
	/// \warning All names must be unique
	G4int id = G4()->CreateH2(name, title, nxbins, xmin, xmax, nybins, ymin, ymax, unitName, fcnName);
	fMap[name] = id;
	return id;
}
