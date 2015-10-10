/// \file texansim.cc
/// \brief Defines main() program
///
#include <memory>

#include "G4RunManager.hh"
#include "G4MTRunManager.hh"
#include "G4UImanager.hh"
#include "G4VisManager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"



// User files //
#include "TexanDetectorConstruction.hh"  // Detector construction
#include "QGSP_BIC_HP.hh"                // Physics list
#include "TexanActionInitialization.hh"  // Action initialization
#include "TexanSensitiveDetector.hh"

/// Namespace for all project classes
namespace texansim {    }
namespace txs = texansim;


namespace {
int usage()
{
	G4cerr << "usage: texansim <run*.mac> [--geo[metry] <*.gdml>] [--vis[ualize]]\n\n";
	return 1;
}
int novis()
{
	G4cerr << "visualization was not enabled at compile time\n"
				 << "re-compile with the proper flags to make use of visualization\n\n";
	return 1;
}

}


int main(int argc, char** argv)
{
	/// - Check arguments
	if(argc < 2) // no macro file specified
		return usage();

	if(G4String(argv[1]) == "--help" || G4String(argv[1]) == "-h")
		return usage();
	if(G4String(argv[1]) == "--visualize" || G4String(argv[1]) == "--vis")
		return usage();
	if(G4String(argv[1]) == "--geometry"  || G4String(argv[1]) == "--geo")
		return usage();
	

	bool visualize = false;
	G4String geofile = TEXAN_BUILD_DIR + G4String("/empty.gdml");

	for(int i = 1; i< argc; ++i) {
		if(0) { }
		else if(G4String(argv[i]) == "--visualize") {
			visualize = true;
		}
		else if(G4String(argv[i]) == "--geo") {
			if(i+1 < argc)
				geofile = argv[++i];
			else
				return usage();
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
#ifdef G4MULTITHREADED
 	std::auto_ptr<G4RunManager> runManager (new G4MTRunManager);
#else
	std::auto_ptr<G4RunManager> runManager (new G4RunManager);
#endif


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


# if 0

	/// - Sensitive detectors
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	txs::SensitiveDetector* scintSD = new txs::SensitiveDetector("Scint");
	SDman->AddNewDetector( scintSD );



	///////////////////////////////////////////////////////////////////////
   //
   // Example how to retrieve Auxiliary Information for sensitive detector
   //
	const G4GDMLAuxMapType* auxmap = parser.GetAuxMap();
	G4cout << "Found " << auxmap->size()
				 << " volume(s) with auxiliary information."
				 << G4endl << G4endl;
	for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin();
			iter!=auxmap->end(); iter++) 
	{
		G4cout << "Volume " << ((*iter).first)->GetName()
					 << " has the following list of auxiliary information: "
					 << G4endl << G4endl;
		for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin();
				 vit!=(*iter).second.end(); vit++)
		{
			G4cout << "--> Type: " << (*vit).type
						 << " Value: " << (*vit).value << G4endl;
		}
	}
	G4cout << G4endl;

	// The same as above, but now we are looking for
	// sensitive detectors setting them for the volumes

	for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin();
			iter!=auxmap->end(); iter++) 
	{
		G4cout << "Volume " << ((*iter).first)->GetName()
					 << " has the following list of auxiliary information: "
					 << G4endl << G4endl;
		for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin();
				 vit!=(*iter).second.end();vit++)
		{
			if ((*vit).type=="SensDet")
			{
				G4cout << "Attaching sensitive detector " << (*vit).value
							 << " to volume " << ((*iter).first)->GetName()
							 <<  G4endl << G4endl;

				G4VSensitiveDetector* mydet = 
					SDman->FindSensitiveDetector((*vit).value);
				if(mydet) 
				{
					G4LogicalVolume* myvol = (*iter).first;
					myvol->SetSensitiveDetector(mydet);
				}
				else
				{
					G4cout << (*vit).value << " detector not found" << G4endl;
				}
			}
		}
	}
	//
	// End of Auxiliary Information block
	//
	////////////////////////////////////////////////////////////////////////

#endif


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
