#include "G4UImanager.hh"
#include "G4RootAnalysisManager.hh"
#include "TexanAnalysis.hh"


namespace txs = texansim;



G4int txs::Analysis::BookH1(const G4String &name, const G4String &title, G4int nbins, G4double xmin, G4double xmax, const G4String &unitName, const G4String &fcnName)
{
	/// \warning All names must be unique
	G4int id = gAna->CreateH1(name, title, nbins, xmin, xmax, unitName, fcnName);
	gBookedMap()[name] = id;
	return id;
}

G4int txs::Analysis::BookH2(const G4String &name, const G4String &title, G4int nxbins, G4double xmin, G4double xmax, G4int nybins, G4double ymin, G4double ymax,  const G4String &unitName, const G4String &fcnName)
{
	/// \warning All names must be unique
	G4int id = gAna->CreateH2(name, title, nxbins, xmin, xmax, nybins, ymin, ymax, unitName, fcnName);
	gBookedMap()[name] = id;
	return id;
}

std::map<G4String, G4int>& txs::Analysis::gBookedMap()
{
	std::map<G4String, G4int>* m = 0;
	if(!m) {
		m = new std::map<G4String, G4int>();
	}
	return *m;
}
