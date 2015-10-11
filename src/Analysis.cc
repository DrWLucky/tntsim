/// \file Analysis.cc
/// \brief Implements analysis classes
#include <cassert>
#include "G4UImanager.hh"
#include "G4RootAnalysisManager.hh"
#include "texansim/Analysis.hh"


namespace txs = texansim;


void txs::Analysis::Cleanup()
{
	/// \attention Call from RunAction destructor
	delete G4AnalysisManager::Instance(); 
}

G4int txs::Analysis::BookH1(const G4String &name, const G4String &title, G4int nbins, G4double xmin, G4double xmax, const G4String &unitName, const G4String &fcnName)
{
	/// \warning All names must be unique
	G4int id = gAna()->CreateH1(name, title, nbins, xmin, xmax, unitName, fcnName);
	gHistMap()[name] = id;
	return id;
}

G4int txs::Analysis::BookH2(const G4String &name, const G4String &title, G4int nxbins, G4double xmin, G4double xmax, G4int nybins, G4double ymin, G4double ymax,  const G4String &unitName, const G4String &fcnName)
{
	/// \warning All names must be unique
	G4int id = gAna()->CreateH2(name, title, nxbins, xmin, xmax, nybins, ymin, ymax, unitName, fcnName);
	gHistMap()[name] = id;
	return id;
}

txs::Analysis::IdLookup_t& txs::Analysis::gHistMap()
{
	static IdLookup_t* m = 0;
	if(!m) {
		m = new IdLookup_t();
	}
	return *m;
}

txs::Analysis::IdLookup_t& txs::Analysis::gColumnMap()
{
	static IdLookup_t* m = 0;
	if(!m) {
		m = new IdLookup_t();
	}
	return *m;
}


G4int txs::Analysis::GetHistId(const G4String &name)
{
	IdLookup_t::const_iterator it = gHistMap().find(name);
	if(it != gHistMap().end())
		return it->second;

	G4cerr << "Error getting hostogram id: \"" << name << "\"\n";
	return -1;
}

G4int txs::Analysis::GetColumnId(const G4String &name)
{
	IdLookup_t::const_iterator it = gColumnMap().find(name);
	if(it != gColumnMap().end())
		return it->second;

	G4cerr << "Error getting column id: \"" << name << "\"\n";
	return -1;
}
