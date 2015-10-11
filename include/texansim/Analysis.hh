#ifndef HAVE_TEXAN_ANALYSIS_HH
#define HAVE_TEXAN_ANALYSIS_HH 1
/// \file TexanAnalysis.hh
/// \brief Class for handling analysis stuff
#include <map>

#if   defined (TEXAN_ANALYZE_ROOT)
#include "g4root.hh"
#elif defined (TEXAN_ANALYZE_XML)
#include "g4xml.hh"
#elif defined (TEXAN_ANALYZE_CSV)
#include "g4csv.hh"
#else 
#error("Analyzer not defined")
#endif




namespace texansim {

/// Analysis class ("static")
/*! Provides function interface mimicing that of G4VAnalysisManager.
	  Except this is in a static class, so there's no singleton instance.
		Also provides additional functionality: lookup of histogram and NTuple
		IDs from their name strings. To use this, call the Book*() functions to
		create histograms and NTuples, and Fill*() functions that take 'string name'
		as an argument to add data. For ntuple columns, make use of templates to
		deterine column type rather than having separate function names.
		
		\attention All names given to objects must be unique (though histograms and Ntuple
		columns can have the same name).

		\note Undocumented functions are simply copies of the corresponding G4VAnalysisManager
		functions.
 */
class Analysis {
public:
	/// Clean up - delete G4AnalysisManager
	static void Cleanup();

  /// Book 1d histogram
	static G4int BookH1(const G4String &name, const G4String &title, G4int nbins, G4double xmin, G4double xmax, const G4String &unitName="none", const G4String &fcnName="none");

  /// Book 2d histogram
	static G4int BookH2(const G4String &name, const G4String &title, G4int nxbins, G4double xmin, G4double xmax, G4int nybins, G4double ymin, G4double ymax,  const G4String &unitName="none", const G4String &fcnName="none");

  /// Look up Histogram ID by name
	static G4int GetHistId(const G4String &name);
	
  /// Look up Ntuple column ID by name
	static G4int GetColumnId(const G4String &name);

	/// Fill 1d histogram by name
	static G4bool FillH1(const G4String& name, G4double value, G4double weight=1.0)
		{ return gAna()->FillH1(GetHistId(name), value, weight); }

	/// Fill 2d histogram by name
	static G4bool FillH2(const G4String& name, G4double xvalue, G4double yvalue, G4double weight=1.0)
		{ return gAna()->FillH2(GetHistId(name), xvalue, yvalue, weight); }

  /// Fill NTuple column by name
	template <class T>
	static G4bool FillNtupleColumn(const G4String& name, const T& value)
		{
			G4cerr << "Warning in FillNtupleColumn: Defaulting fill to double\n";
			return gAna()->FillNtupleDColumn(GetColumnId(name), value);
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
			G4int retval = gAna()->CreateNtupleDColumn(name);

			gColumnMap()[name] = retval;
			return retval;
		}



// COPIES OF G4VAnalysisManager functions //

	static G4bool   	OpenFile () 
		{ return gAna()->OpenFile(); }
	static G4bool   	OpenFile (const G4String &fileName) 
		{ return gAna()->OpenFile(fileName); }
	static G4bool   	Write () 
		{ return gAna()->Write(); }
	static G4bool   	CloseFile () 
		{ return gAna()->CloseFile(); }
	static G4bool   	SetFileName (const G4String &fileName) 
		{ return gAna()->SetFileName(fileName); }
	static G4bool   	SetHistoDirectoryName (const G4String &dirName) 
		{ return gAna()->SetHistoDirectoryName(dirName); }
	static G4bool   	SetNtupleDirectoryName (const G4String &dirName) 
		{ return gAna()->SetNtupleDirectoryName(dirName); }
	static G4String 	GetFileName ()   
		{ return gAna()->GetFileName(); }
	static G4String 	GetHistoDirectoryName ()  
		{ return gAna()->GetHistoDirectoryName(); }
	static G4String 	GetNtupleDirectoryName ()  
		{ return gAna()->GetNtupleDirectoryName(); }
	static G4int 	    CreateH1 (const G4String &name, const G4String &title,
															G4int nbins, G4double xmin, G4double xmax,
															const G4String &unitName="none", const G4String &fcnName="none")
		{ return gAna()->CreateH1(name, title, nbins, xmin, xmax, unitName, fcnName); }
	static G4int 	    CreateH2 (const G4String &name, const G4String &title,
															G4int nxbins, G4double xmin, G4double xmax,
															G4int nybins, G4double ymin, G4double ymax,
															const G4String &xunitName="none", const G4String &yunitName="none",
															const G4String &xfcnName="none", const G4String &yfcnName="none")
		{ return gAna()->CreateH2(name, title, nxbins, xmin, xmax, nybins, ymin, ymax, xunitName, xfcnName, yunitName, yfcnName ); }
	static G4bool   	SetH1 (G4int id, G4int nbins, G4double xmin, G4double xmax,
													 const G4String &unitName="none", const G4String &fcnName="none")
		{ return gAna()->SetH1(id, nbins, xmin, xmax, unitName, fcnName); }
	static G4bool   	SetH2 (G4int id, G4int nxbins, G4double xmin, G4double xmax,
													 G4int nybins, G4double ymin, G4double ymax,
													 const G4String &xunitName="none", const G4String &yunitName="none",
													 const G4String &xfcnName="none", const G4String &yfcnName="none")
		{ return gAna()->SetH2(id, nxbins, xmin, xmax, nybins, ymin, ymax, xunitName, yunitName, xfcnName, yfcnName); }
	static G4bool   	ScaleH1 (G4int id, G4double factor) 
		{ return gAna()->ScaleH1(id, factor); }
	static G4bool   	ScaleH2 (G4int id, G4double factor) 
		{ return gAna()->ScaleH2(id, factor); }
	static void 	    CreateNtuple (const G4String &name, const G4String &title) 
		{ gAna()->CreateNtuple(name, title); }
	static G4int 	    CreateNtupleIColumn (const G4String &name) 
		{ return gAna()->CreateNtupleIColumn(name); }
	static G4int      CreateNtupleFColumn (const G4String &name) 
		{ return gAna()->CreateNtupleFColumn(name); }
	static G4int     	CreateNtupleDColumn (const G4String &name) 
		{ return gAna()->CreateNtupleDColumn(name); }
	static void     	FinishNtuple () 
		{ gAna()->FinishNtuple(); }
	static G4bool   	SetFirstHistoId (G4int firstId) 
		{ return gAna()->SetFirstHistoId(firstId); }
	static G4bool   	SetFirstNtupleColumnId (G4int firstId) 
		{ return gAna()->SetFirstNtupleColumnId(firstId); }
	static G4bool   	FillH1 (G4int id, G4double value, G4double weight=1.0) 
		{ return gAna()->FillH1(id, value, weight); }
	static G4bool   	FillH2 (G4int id, G4double xvalue, G4double yvalue, G4double weight=1.0) 
		{ return gAna()->FillH2(id, xvalue, yvalue, weight); }
	static G4bool   	FillNtupleIColumn (G4int id, G4int value)    
		{ return gAna()->FillNtupleIColumn(id, value); }
	static G4bool   	FillNtupleFColumn (G4int id, G4float value)  
		{ return gAna()->FillNtupleFColumn(id, value); }
	static G4bool   	FillNtupleDColumn (G4int id, G4double value) 
		{ return gAna()->FillNtupleDColumn(id, value); }
	static G4bool   	AddNtupleRow () 
		{ return gAna()->AddNtupleRow (); }
	static void     	SetActivation (G4bool activation) 
		{ gAna()->SetActivation(activation); }
	static G4bool   	GetActivation () 
		{ return gAna()->GetActivation(); }
	static G4bool   	IsActive () 
		{ return gAna()->IsActive(); }
	static G4bool   	IsAscii () 
		{ return gAna()->IsAscii(); }
	static G4int 	    GetNofH1s () 
		{ return gAna()->GetNofH1s (); }
	static G4int 	    GetNofH2s ()  
		{ return gAna()->GetNofH2s (); }
	static G4int 	    GetH1Nbins (G4int id)      
		{ return gAna()->GetH1Nbins(id); }
	static G4double 	GetH1Xmin (G4int id)       
		{ return gAna()->GetH1Xmin(id); }
	static G4double 	GetH1Xmax (G4int id)       
		{ return gAna()->GetH1Xmax (id); }      
	static G4double 	GetH1Width (G4int id)      
		{ return gAna()->GetH1Width (id); }     
	static G4int 	    GetH2Nxbins (G4int id)     
		{ return gAna()->GetH2Nxbins (id); }    
	static G4double 	GetH2Xmin (G4int id)       
		{ return gAna()->GetH2Xmin (id); }      
	static G4double 	GetH2Xmax (G4int id)       
		{ return gAna()->GetH2Xmax (id); }      
	static G4double 	GetH2XWidth (G4int id)     
		{ return gAna()->GetH2XWidth (id); }    
	static G4int     	GetH2Nybins (G4int id)     
		{ return gAna()->GetH2Nybins (id); }    
	static G4double 	GetH2Ymin (G4int id)       
		{ return gAna()->GetH2Ymin (id); }      
	static G4double 	GetH2Ymax (G4int id)       
		{ return gAna()->GetH2Ymax (id); }      
	static G4double 	GetH2YWidth (G4int id)     
		{ return gAna()->GetH2YWidth (id); }    
	static G4String 	GetH1Name (G4int id)       
		{ return gAna()->GetH1Name (id); }      
	static G4double 	GetH1Unit (G4int id)       
		{ return gAna()->GetH1Unit (id); }      
	static G4bool 	  GetH1Activation (G4int id) 
		{ return gAna()->GetH1Activation (id); }
	static G4bool 	  GetH1Ascii (G4int id)      
		{ return gAna()->GetH1Ascii (id); }     
	static G4String 	GetH2Name (G4int id)       
		{ return gAna()->GetH2Name (id); }      
	static G4double 	GetH2XUnit (G4int id)      
		{ return gAna()->GetH2XUnit (id); }     
	static G4double 	GetH2YUnit (G4int id)      
		{ return gAna()->GetH2YUnit (id); }     
	static G4bool 	  GetH2Activation (G4int id) 
		{ return gAna()->GetH2Activation (id); }
	static G4bool 	  GetH2Ascii (G4int id)      
		{ return gAna()->GetH2Ascii (id); }     
	static G4bool 	  SetH1Title (G4int id, const G4String &title)
		{ return gAna()->SetH1Title      (id, title); }
	static G4bool 	  SetH1XAxisTitle (G4int id, const G4String &title)
		{ return gAna()->SetH1XAxisTitle (id, title); }
	static G4bool 	  SetH1YAxisTitle (G4int id, const G4String &title)
		{ return gAna()->SetH1YAxisTitle (id, title); }
	static G4bool 	  SetH2Title (G4int id, const G4String &title)
		{ return gAna()->SetH2Title      (id, title); }
	static G4bool 	  SetH2XAxisTitle (G4int id, const G4String &title)
		{ return gAna()->SetH2XAxisTitle (id, title); }
	static G4bool 	  SetH2YAxisTitle (G4int id, const G4String &title)
		{ return gAna()->SetH2YAxisTitle (id, title); }
	static G4bool 	  SetH2ZAxisTitle (G4int id, const G4String &title)
		{ return gAna()->SetH2ZAxisTitle (id, title); }
	static G4String 	GetH1Title (G4int id) 
		{ return gAna()->GetH1Title(id); } 
	static G4String 	GetH1XAxisTitle (G4int id) 
		{ return gAna()->GetH1XAxisTitle(id); } 
	static G4String 	GetH1YAxisTitle (G4int id) 
		{ return gAna()->GetH1YAxisTitle(id); } 
	static G4String 	GetH2Title (G4int id) 
		{ return gAna()->GetH2Title(id); } 
	static G4String 	GetH2XAxisTitle (G4int id) 
		{ return gAna()->GetH2XAxisTitle(id); } 
	static G4String 	GetH2YAxisTitle (G4int id) 
		{ return gAna()->GetH2YAxisTitle(id); } 
	static G4String 	GetH2ZAxisTitle (G4int id) 
		{ return gAna()->GetH2ZAxisTitle(id); } 
	static G4int 	    GetVerboseLevel ()
		{ return gAna()->GetVerboseLevel(); }
	static void 	    SetVerboseLevel (G4int verboseLevel)
		{ gAna()->SetVerboseLevel(verboseLevel); }
	static G4String 	GetType ()
		{ return gAna()->GetType(); }
	static G4String 	GetFileType ()
		{ return gAna()->GetFileType(); }

// END COPIES //

private:
	typedef std::map<G4String, G4int> IdLookup_t;

/// Singleton map for histograms (1d and 2d)
	static IdLookup_t& gHistMap();
/// Singleton map for NTuple columns (all types)
	static IdLookup_t& gColumnMap();
/// Singleton G4AnalysManager instance (just an alias)
	static G4AnalysisManager* gAna() { return G4AnalysisManager::Instance(); }
};






// TEMPLATE SPECIALIZATIONS //


/// Specialization for double
template <>
inline G4bool Analysis::FillNtupleColumn<G4double>(const G4String& name, const G4double& value)
{
	return gAna()->FillNtupleDColumn(GetColumnId(name), value);
}

/// Specialization for float
template <>
inline G4bool Analysis::FillNtupleColumn<G4float>(const G4String& name, const G4float& value)
{
	return gAna()->FillNtupleFColumn(GetColumnId(name), value);
}

/// Specialization for int
template <>
inline G4bool Analysis::FillNtupleColumn<G4int>(const G4String& name, const G4int& value)
{
	return gAna()->FillNtupleIColumn(GetColumnId(name), value);
}



/// Specialization for double
template <>
inline G4int Analysis::BookNtupleColumn<G4double>(const G4String& name)
{
	G4int retval = gAna()->CreateNtupleDColumn(name);
	gColumnMap()[name] = retval;
	return retval;
}

/// Specialization for float
template <>
inline G4int Analysis::BookNtupleColumn<G4float>(const G4String& name)
{
	G4int retval = gAna()->CreateNtupleFColumn(name);
	gColumnMap()[name] = retval;
	return retval;
}

/// Specialization for int
template <>
inline G4int Analysis::BookNtupleColumn<G4int>(const G4String& name)
{
	G4int retval = gAna()->CreateNtupleIColumn(name);
	gColumnMap()[name] = retval;
	return retval;
}

} // namespace txsim

#endif
