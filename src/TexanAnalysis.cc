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

G4VAnalysisManager*& txs::AnalysisManager::G4()
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
