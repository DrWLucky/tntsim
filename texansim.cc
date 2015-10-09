/// \file texansim.cc
/// \brief Defines main() program
///
#include <memory>
#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"


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


	/// Construct the default run manager
#ifdef G4MULTITHREADED
 	std::auto_ptr<G4RunManager> runManager (new G4MTRunManager);
#else
	std::auto_ptr<G4RunManager> runManager (new G4RunManager);
#endif

	G4UImanager* UI = G4UImanager::GetUIpointer();
	// UI->ApplyCommand("/persistency/gdml/read detector1.gdml");
	// UI->ApplyCommand("/texansim/detectorFile detector1.gdml");


	/// - Set mandatory initialization classes
	///
	/// - Detector construction from GDML file (XML)
	txs::DetectorConstruction* det = new txs::DetectorConstruction();
	runManager->SetUserInitialization(det);

	/// - TEMPORARY physics list
	/// \todo implement permanent (or changable) physics list
	runManager->SetUserInitialization(new QGSP_BIC_HP);
	/// - Action initialization
	runManager->SetUserInitialization(new txs::ActionInitialization);

	/// - Initialize G4 kernel
	runManager->Initialize();

  /// - Read a macro file of commands
  G4String command  = "/control/execute ";
  G4String fileName = argv[1];
	UI->ApplyCommand(command+fileName); 

	return 0;
}
