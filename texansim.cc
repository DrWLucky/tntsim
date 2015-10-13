/// \file texansim.cc
/// \brief Defines main() program
///
#include <memory>
#include <string>

#include "G4RunManager.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#endif
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

// User files //
#include "QGSP_BIC_HP.hh"                 // Physics list
#include "texansim/DetectorConstruction.hh"  // Detector construction
#include "texansim/ActionInitialization.hh"  // Action initialization
#include "texansim/Analysis.hh"
#include "texansim/VPersistenceManager.hh"
#include "texansim/PersistenceMessenger.hh"



namespace txs = texansim;


namespace {
int usage()
{
	G4cerr << "usage: texansim <run*.mac> [--geo[metry]=*.gdml] [--threads=N] [--vis[ualize]]\n\n";
	return 1;
}
int novis()
{
	G4cerr << "visualization was not enabled at compile time\n"
				 << "re-compile with the proper flags to make use of visualization\n\n";
	return 1;
}

}


/// Namespace for all project classes
namespace texansim {

/// The main program
int main(int argc, char** argv)
{
	/// I'm just putting this in the texansim namespace so that it
	/// shows up under the namespace section of doxygen. The real
	/// main() is a one-liner simply calling this.

	/// - Check arguments
	if(argc < 2) // no macro file specified
		return usage();

	std::string arg1(argv[1]);
	if(arg1 == "--help" || arg1 == "-h")
		return usage();
	if(arg1.substr(0,10) == "--geometry"  || arg1.substr(0,5) == "--geo")
		return usage();
	if(arg1.substr(0,9) == "--threads")
		return usage();
	

	/// - Needed for thread safety: call Analysis::CallSingletons();
	Analysis::ConstructSingletons();

	G4int nthreads = 0;
	bool visualize = false;
	G4String geofile = TEXAN_BUILD_DIR + G4String("/test.gdml");

	for(int i = 1; i< argc; ++i) {
		std::string arg = argv[i];
		if(0) { }
		else if(arg == "--visualize" || arg == "--vis") {
			visualize = true;
		}
		else if(arg.substr(0,6) == "--geo=") {
			geofile = arg.substr(6);
		}
		else if(arg.substr(0,11) == "--geometry=") {
			geofile = arg.substr(11);
		}
		else if(arg.substr(0, 10) == "--threads=") {
			nthreads = atoi(arg.substr(10).c_str());
			if(nthreads == 1) nthreads = 0;
		}
	}
	
	if(visualize) {
#if defined(G4UI_USE) && defined(G4VIS_USE)
		;
#else
		return novis();
#endif
	}


	/// - Construct the default run manager
	std::auto_ptr<G4RunManager> runManager(0);

	if(nthreads != 0) { // Multi-threaded
#ifdef G4MULTITHREADED
		runManager.reset(new G4MTRunManager());
		static_cast<G4MTRunManager*>(runManager.get())->SetNumberOfThreads(nthreads);
		G4cout << "Enabling multi-threaded support, number of threads = " << nthreads << G4endl;
#else
		G4cerr << "Multi-threading support not enabled!\n";
		return 1;
#endif
	} else { // Single thread
		runManager.reset(new G4RunManager());
	}

	txs::VPersistenceManager::SetMessenger(new PersistenceMessenger());


	/// - Set mandatory initialization classes
	///
	/// - Detector construction from GDML file (XML)
	txs::DetectorConstruction* det = new txs::DetectorConstruction(geofile);
	runManager->SetUserInitialization(det);

	/// - TEMPORARY physics list
	/// \todo implement permanent (or changable) physics list
	runManager->SetUserInitialization(new QGSP_BIC_HP);
	/// - Action initialization
	runManager->SetUserInitialization(new txs::ActionInitialization);

	/// - Initialize G4 kernel
	runManager->Initialize();


	/// - UI stuff
	G4UImanager* UI = G4UImanager::GetUIpointer();

	if(!visualize) { // run simulation
		G4String command  = "/control/execute ";
		G4String fileName = argv[1];
		UI->ApplyCommand(command + fileName); 
	} else { // visualize
		std::auto_ptr<G4UIExecutive> ui(new G4UIExecutive(argc, argv));
		{
			std::auto_ptr<G4VisManager> visManager(new G4VisExecutive);
			visManager->Initialize();
			UI->ApplyCommand("/control/execute vis.mac");
			ui->SessionStart();
		}
	}

	return 0;
}

}


int main(int argc, char** argv)
{ return txs::main(argc, argv); }
