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
	// Check for macro file argument
	if(argc < 2) {
		G4cerr << "usage: texansim <run.mac>\n";
		return 1;
	}

	// GDML parser
	G4GDMLParser parser;
	parser.Read("/home/gacgroup/gchristian/packages/simulation/texansim/build/detectors.gdml");

	/// Construct the default run manager
#ifdef G4MULTITHREADED
 	std::auto_ptr<G4RunManager> runManager (new G4MTRunManager);
#else
	std::auto_ptr<G4RunManager> runManager (new G4RunManager);
#endif

	/// - Set mandatory initialization classes
	///
	/// - Detector construction from GDML file (XML)
	runManager->SetUserInitialization(new txs::DetectorConstruction(parser.GetWorldVolume()));
	/// - TEMPORARY physics list
	/// \todo implement permanent (or changable) physics list
	runManager->SetUserInitialization(new QGSP_BIC_HP);
	/// - Action initialization
	runManager->SetUserInitialization(new txs::ActionInitialization);

	/// - Initialize G4 kernel
	runManager->Initialize();

  /// - Read a macro file of commands
	G4UImanager* UI = G4UImanager::GetUIpointer();
  G4String command  = "/control/execute ";
  G4String fileName = argv[1];
	UI->ApplyCommand(command+fileName); 

	return 0;
}
