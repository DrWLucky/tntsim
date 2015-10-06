/// \file texansim.cc
/// \brief Defines main() program
///
#include <memory>
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4GDMLParser.hh"

// User files //
#include "TexanDetectorConstruction.hh"  // Detector construction
#include "QGSP_BIC_HP.hh"                // Physics list
#include "TexanActionInitialization.hh"  // Action initialization


/// Namespace for all project classes
namespace texansim {    }
namespace txs = texansim;


int main(int argc, char** argv)
{
	// GDML parser
	G4GDMLParser parser;

	if(argc > 1)
		parser.Read(argv[1]);
	else
		parser.Read("detectors.gdml");

	// construct the default run manager
#ifdef G4MULTITHREADED
	std::auto_ptr<G4RunManager> runManager (new G4MTRunManager);
#else
	std::auto_ptr<G4RunManager> runManager (new G4RunManager);
#endif

	// set mandatory initialization classes
	 
	/// - Detector construction from GDML file (XML)
	runManager->SetUserInitialization(new txs::DetectorConstruction(parser.GetWorldVolume()));
	/// - TEMPORARY physics list
	/// \todo implement permanent (or changable) physics list
	runManager->SetUserInitialization(new QGSP_BIC_HP);
	/// - Action initializatioj
	runManager->SetUserInitialization(new txs::ActionInitialization);

	// initialize G4 kernel
	runManager->Initialize();

	// get the pointer to the UI manager and set verbosities
	G4UImanager* UI = G4UImanager::GetUIpointer();
	UI->ApplyCommand("/run/verbose 1");
	UI->ApplyCommand("/event/verbose 1");
	UI->ApplyCommand("/tracking/verbose 1");

	// start a run
	int numberOfEvent = 3;
	runManager->BeamOn(numberOfEvent);

	// job termination
	return 0;
}
