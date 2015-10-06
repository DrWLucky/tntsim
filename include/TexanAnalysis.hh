#ifndef HAVE_TEXAN_ANALYSIS_HH
#define HAVE_TEXAN_ANALYSIS_HH 1
/// \file TexanAnalysis.hh
/// \brief Class for handling analysis stuff
#include <map>
#include <typeinfo>

#include "g4root.hh"




namespace txprivate
{
/// Private map 'singleton' instance
/** \attention Do not use this in your code */
std::map<G4String, G4int>& gBookedMap();
}



namespace texansim {

/// Book 1d histogram
G4int BookH1(const G4String &name, const G4String &title, G4int nbins, G4double xmin, G4double xmax, const G4String &unitName="none", const G4String &fcnName="none");

/// Book 2d histogram
G4int BookH2(const G4String &name, const G4String &title, G4int nxbins, G4double xmin, G4double xmax, G4int nybins, G4double ymin, G4double ymax,  const G4String &unitName="none", const G4String &fcnName="none");

/// Book NTuple column
/** Valid template types G4double, G4float, G4int, defaults to G4double if other type specified */
template <class T> 
G4int BookNtupleColumn(const G4String& name);

/// Look up object ID by name
inline G4int GetId(const G4String &name) { return txprivate::gBookedMap()[name];	}

/// Fill NTuple column by name
template <class T>
G4bool FillNtupleColumn (const G4String& name, const T& value);

/// Alias the global analysis manager instance
static G4AnalysisManager* gAna = G4AnalysisManager::Instance();

}



namespace {
template <class T> G4bool do_fill(G4int id, const T& value)
{
	G4cerr << "Warning: Defaulting fill to double\n";
	return texansim::gAna->FillNtupleDColumn(id, value);
}

template <> G4bool do_fill<G4double>(G4int id, const G4double& value)
{
	return texansim::gAna->FillNtupleDColumn(id, value);
}

template <> G4bool do_fill<G4float>(G4int id, const G4float& value)
{
	return texansim::gAna->FillNtupleFColumn(id, value);
}

template <> G4bool do_fill<G4int>(G4int id, const G4int& value)
{
	return texansim::gAna->FillNtupleIColumn(id, value);
}

}


template <class T>
G4bool texansim::FillNtupleColumn(const G4String& name, const T& value)
{
	/// \warning All names must be unique
	G4int id = GetId(name);
	G4cerr << id << "  <<<<<< ID\n";
	return do_fill<T> (id, value);
}



template <class T>
inline G4int texansim::BookNtupleColumn(const G4String& name)
{
	/// \warning All names must be unique
	G4int retval;
	
	if (0) { }
	else if (typeid(T) == typeid(G4double)) {
		retval = gAna->CreateNtupleDColumn(name);
	}
	else if (typeid(T) == typeid(G4float)) {
		retval = gAna->CreateNtupleFColumn(name);
	}
	else if (typeid(T) == typeid(G4int)) {
		retval = gAna->CreateNtupleIColumn(name);
	}
	else {
		G4cerr << "Warning in BookNtupleColumn: Defaulting column type to double for ntuple column \"" << name << "\"\n";
		retval = gAna->CreateNtupleDColumn(name);
	}

	txprivate::gBookedMap()[name] = retval;
	return retval;
}


#endif
