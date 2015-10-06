#ifndef HAVE_TEXAN_ANALYSIS_HH
#define HAVE_TEXAN_ANALYSIS_HH 1
/// \file TexanAnalysis.hh
/// \brief Class for handling analysis stuff
#include "G4VAnalysisManager.hh"


namespace texansim {

/// Class to manage data analysis.
/*! Hooks into the G4VAnalysisManagerFramework, also 
	providing some ways to simplify access
*/
class AnalysisManager {
public:
	/// Global instance
	static AnalysisManager*& Instance();
	/// Global G4Vanalysis manager instance
	G4VAnalysisManager*& G4();
	/// Deletes global G4vanalysismanager
	~AnalysisManager();
protected:
	/// Empty
	AnalysisManager();
};

/// Shorthand access to global instance
inline AnalysisManager*& Ana()
{
	return AnalysisManager::Instance();
}

}


#endif
