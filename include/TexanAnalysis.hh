#if 0
// #ifndef HAVE_TEXAN_ANALYSIS_HH
// #define HAVE_TEXAN_ANALYSIS_HH 1
/// \file TexanAnalysis.hh
/// \brief Class for handling analysis stuff
#include <map>
#include <typeinfo>

#include "G4VAnalysisManager.hh"



namespace texansim {

/// Class to manage data analysis.
/*! Hooks into the G4VAnalysisManagerFramework, also 
	providing some ways to simplify access. Allows lookup
	of object (histogram, ntuple column) IDs by name.
	\warning All given names using the Book*() functions must be unique
*/
class AnalysisManager {
public:
	/// Global instance
	static AnalysisManager*& Instance();
	/// Global G4Vanalysis manager instance
	G4VAnalysisManager*& G4();
	/// Deletes global G4vanalysismanager
	virtual ~AnalysisManager();

	/// Book 1d histogram
	G4int BookH1(const G4String &name, const G4String &title, G4int nbins, G4double xmin, G4double xmax, const G4String &unitName="none", const G4String &fcnName="none");

	/// Book 2d histogram
	G4int BookH2(const G4String &name, const G4String &title, G4int nxbins, G4double xmin, G4double xmax, G4int nybins, G4double ymin, G4double ymax,  const G4String &unitName="none", const G4String &fcnName="none");

	/// Create Ntuple
	void CreateNtuple(const G4String& name, const G4String& title)
		{ G4()->CreateNtuple(name, title); }

	/// Book NTuple column
	/* Valid template types G4double, G4float, G4int, defaults to G4double if other type specified */
	template <class T> 
	G4int BookNtupleColumn(const G4String& name);
	
	/// Look up object ID by name
	G4int GetId(const G4String &name) { return fMap[name];	}

	/// Fill NTuple column by name
	template <class T>
	G4bool FillNtupleColumn (const G4String& name, const T& value);
	
protected:
	/// Empty
	AnalysisManager();

private:
	std::map<G4String, G4int> fMap;
};

/// Shorthand access to global instance
inline AnalysisManager*& Ana()
{
	return AnalysisManager::Instance();
}

}

namespace {
template <class T> G4bool do_fill(G4int id, const T& value)
{
	G4cerr << "Warning: Defaulting fill to double\n";
	return texansim::Ana()->G4()->FillNtupleDColumn(id, value);
}

template <> G4bool do_fill<G4double>(G4int id, const G4double& value)
{
	return texansim::Ana()->G4()->FillNtupleDColumn(id, value);
}

template <> G4bool do_fill<G4float>(G4int id, const G4float& value)
{
	return texansim::Ana()->G4()->FillNtupleFColumn(id, value);
}

template <> G4bool do_fill<G4int>(G4int id, const G4int& value)
{
	return texansim::Ana()->G4()->FillNtupleIColumn(id, value);
}

}


template <class T>
G4bool texansim::AnalysisManager::FillNtupleColumn(const G4String& name, const T& value)
{
	/// \warning All names must be unique
	G4int id = GetId(name);
	return do_fill<T> (id, value);
}



template <class T>
inline G4int texansim::AnalysisManager::BookNtupleColumn(const G4String& name)
{
	/// \warning All names must be unique
	G4int retval;
	
	if (0) { }
	// else if (typeid(T) == typeid(G4double)) {
	// 	retval = G4()->CreateNtupleDColumn(name);
	// }
	// else if (typeid(T) == typeid(G4float)) {
	// 	retval = G4()->CreateNtupleFColumn(name);
	// }
	// else if (typeid(T) == typeid(G4int)) {
	// 	retval = G4()->CreateNtupleIColumn(name);
	// }
	else {
		G4cerr << "Warning: Defaulting column type to double for ntuple column \"" << name << "\"\n";
		retval = G4()->CreateNtupleDColumn(name);
	}

	fMap[name] = retval;
	return retval;
}


#endif
