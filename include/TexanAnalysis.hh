#ifndef HAVE_TEXAN_ANALYSIS_HH
#define HAVE_TEXAN_ANALYSIS_HH 1
/// \file TexanAnalysis.hh
/// \brief Class for handling analysis stuff
#include <map>


#include "g4root.hh"




namespace txprivate
{
/// Private map 'singleton' instance
/** \attention Do not use this in your code */
std::map<G4String, G4int>& gBookedMap();
}



namespace texansim {

/// Alias the global analysis manager instance
static G4AnalysisManager* gAna = G4AnalysisManager::Instance();



class Analysis {
public:
  /// Book 1d histogram
	static G4int BookH1(const G4String &name, const G4String &title, G4int nbins, G4double xmin, G4double xmax, const G4String &unitName="none", const G4String &fcnName="none");

  /// Book 2d histogram
	static G4int BookH2(const G4String &name, const G4String &title, G4int nxbins, G4double xmin, G4double xmax, G4int nybins, G4double ymin, G4double ymax,  const G4String &unitName="none", const G4String &fcnName="none");

  /// Look up object ID by name
	static G4int GetId(const G4String &name) { return txprivate::gBookedMap()[name];	}

	/// Add ntuple row
	static G4bool AddNtupleRow() { return gAna->AddNtupleRow(); }

  /// Fill NTuple column by name
	template <class T>
	static G4bool FillNtupleColumn(const G4String& name, const T& value)
		{
			G4cerr << "Warning in FillNtupleColumn: Defaulting fill to double\n";
			return texansim::gAna->FillNtupleDColumn(GetId(name), value);
		}

  /// Book NTuple column
  /** Valid template types G4double, G4float, G4int, defaults to G4double
			(with warning) if other type specified */
	template <class T>
	static G4int BookNtupleColumn(const G4String& name)
		{
			/// \warning All names must be unique
			G4cerr << "Warning in BookNtupleColumn: Defaulting column type to double for ntuple column \""
						 << name << "\"\n";
			G4int retval = gAna->CreateNtupleDColumn(name);

			txprivate::gBookedMap()[name] = retval;
			return retval;
		}
};


/// Specialization for double
template <>
inline G4bool Analysis::FillNtupleColumn<G4double>(const G4String& name, const G4double& value)
{
	return texansim::gAna->FillNtupleDColumn(Analysis::GetId(name), value);
}

/// Specialization for float
template <>
inline G4bool Analysis::FillNtupleColumn<G4float>(const G4String& name, const G4float& value)
{
	return texansim::gAna->FillNtupleFColumn(Analysis::GetId(name), value);
}

/// Specialization for int
template <>
inline G4bool Analysis::FillNtupleColumn<G4int>(const G4String& name, const G4int& value)
{
	return texansim::gAna->FillNtupleIColumn(Analysis::GetId(name), value);
}



/// Specialization for double
template <>
inline G4int Analysis::BookNtupleColumn<G4double>(const G4String& name)
{
	G4int retval = gAna->CreateNtupleDColumn(name);
	txprivate::gBookedMap()[name] = retval;
	return retval;
}

/// Specialization for float
template <>
inline G4int Analysis::BookNtupleColumn<G4float>(const G4String& name)
{
	G4int retval = gAna->CreateNtupleFColumn(name);
	txprivate::gBookedMap()[name] = retval;
	return retval;
}

/// Specialization for int
template <>
inline G4int Analysis::BookNtupleColumn<G4int>(const G4String& name)
{
	G4int retval = gAna->CreateNtupleIColumn(name);
	txprivate::gBookedMap()[name] = retval;
	return retval;
}

} // namespace txsim

#endif
