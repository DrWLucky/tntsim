//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: Tnt_EU.cc,v 2.0 2008/10/01 17:47:10 roeder Exp $
// GEANT4 tag $Name: geant4-09-02.p01 $
//
// Modified by Brian Roeder, TAMU on 7/7/2009
// email - broeder@comp.tamu.edu
// 
//--------------------------------------------------------------------
// Latest version of Tnt Simulation for a DEMON module
//--------------------------------------------------------------------
//


#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

#include "G4ThreeVector.hh"
#include "TntDetectorConstruction.hh"
#include "TntPrimaryGeneratorAction.hh"
#include "TntEventAction.hh"
#include "TntSD.hh"
#include "TntPhysicsList.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#include "TntDataRecordTree.hh"

int main()
{
  // Set beamType for Primary Generator action
  G4String beamType;
  beamType = "pencil";
  //beamType = "diffuse";
  //beamType = "conic";
  //beamType = "iso";

  G4double Neutron_Eng = 15.0*MeV;    // Starting Energy
  G4double Det_Threshold = 0.5*MeV;  // Det. Threshold in MeVee

  // Detector-Types Defined (In TntDetectorConstruction.cc):
  // cylinder - DetDim[0] = rad. of cyl., 
  //            DetDim[1] = DetDepth/2., 
  //            DetDim[2]->NotUsed
  // bar - DetDim[0] = X Dim, DetDim[1] = Y Dim, DetDim[2] = ZDim

# if 0
  // Cylinder defined is DEMON detector module (NE213)
	G4String DetMaterial = "NE213";
  G4String DetType = "cylinder";
  G4ThreeVector DetDim;
  DetDim[0] = 8.0*cm;     // rad. of cyl.
  DetDim[1] = 10.0*cm;    //DetDepth/2.
  DetDim[2] = 0.0*cm;     // z-dim
#else
  // TNT detector
	G4String DetMaterial = "BC505";
  G4String DetType = "bar";
  G4ThreeVector DetDim;
  DetDim[0] = 200.0*cm; // x
  DetDim[1] = 200.0*cm; // y
  DetDim[2] = 30.00*cm; // z
#endif

  // Possible LightConv Settings:
  // "none" - no light conversion - get only energy deposited
  // "light-conv" - convert energy dep. to light following method of Cecil et al.
  // "light+resol" - convert to light with Cecil's formula and Dekempner's resolution 

  G4String LightConv = "light+resol";

  // Set VerboseFlag to "1" to run in verbose mode
  G4int VerboseFlag = 0;
  // Set VisFlag to "1" to run in visualization mode
  G4int VisFlag = 1;
  //**** Note! You should run fewer events (such as < 1000) when using the viewer!
  G4String VisType = "OPENGL";   // Set vis type. OPENGL works in UBUNTU 12.04 and later
   //G4String VisType = "VRMLVIEW";
  G4int SaveVisFile = 0;         // Set to 1 to save OPENGL picture to G4OpenGL.eps file

  //Set number of events per energy in sim.
  G4int ch_eng = 0;
  //G4int ch_eng = 20000;
  G4int numberOfEvent;
  if (ch_eng == 0)
    {numberOfEvent = 10;}
  else
    { numberOfEvent = ch_eng*200; }  // 200 energies from 1 to 200 MeV

//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
// Construct the default run manager
//-----------------------------------------------------------------------
//-----------------------------------------------------------------------
  G4RunManager* runManager = new G4RunManager;

  // Start Time and Random Number Engine
  G4int start_time = time(NULL);  // stores initial time of program start
  CLHEP::RanecuEngine* theEngine = new CLHEP::RanecuEngine();
  CLHEP::HepRandom::setTheEngine(theEngine);
  CLHEP::HepRandom::setTheSeed(start_time);
  CLHEP::HepRandom::showEngineStatus();

  // set mandatory initialization classes
  G4VUserDetectorConstruction* detector = new TntDetectorConstruction(DetMaterial,DetType,DetDim, LightConv);
  runManager->SetUserInitialization(detector);

  G4VUserPhysicsList* physics = new TntPhysicsList;
  runManager->SetUserInitialization(physics);

  // Visualization, if you choose to have it!
#ifdef G4VIS_USE
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // Generate Analysis Pointer Class
  TntDataRecordTree* TntPointer = new TntDataRecordTree(Det_Threshold);
  cout << TntPointer << endl;

  G4VUserPrimaryGeneratorAction* gen_action = new TntPrimaryGeneratorAction(ch_eng,beamType,Neutron_Eng);  
  runManager->SetUserAction(gen_action);

  //set optional user action classes

  //Include Event Action Sequence to invoke event by event analysis
  G4UserEventAction* event_action = new TntEventAction();
  runManager->SetUserAction(event_action);
  
  // Initialize G4 kernel
  //
  runManager->Initialize();

  // Get the pointer to the UI manager and set verbosities
  G4UImanager* UI = G4UImanager::GetUIpointer();

  if (VerboseFlag == 1)
    {
  UI->ApplyCommand("/run/verbose 0");
  UI->ApplyCommand("/event/verbose 1");  
  UI->ApplyCommand("/tracking/verbose 2");
    }
  else
    {/* silent running */}
 
  // Open visualization window to see stuff made; particles.
 
  if (VisFlag == 1)
    {
      if(VisType == "VRMLVIEW")
	{
	  // commands to open vrmlview file
	  UI->ApplyCommand("/vis/scene/create");
	  UI->ApplyCommand("/vis/open VRML2FILE");
	  UI->ApplyCommand("/tracking/storeTrajectory 1");
	  UI->ApplyCommand("/vis/scene/endOfEventAction accumulate");
	  UI->ApplyCommand("/vis/scene/add/trajectories");
	}
 else if(VisType == "OPENGL")
	{
	  UI->ApplyCommand("/vis/open OGL 900x600-0+0");
	  //UI->ApplyCommand("/vis/open OGLI 900x600-0+0");
	  //UI->ApplyCommand("/vis/open DAWNFILE");
	  UI->ApplyCommand("/vis/viewer/set/autoRefresh true");
	  //UI->ApplyCommand("/vis/viewer/set/background red ! ! 0.2");
	  UI->ApplyCommand("/vis/viewer/set/background 0.5 0.5 0.5 0.1");
	  UI->ApplyCommand("/vis/verbose errors");
	  UI->ApplyCommand("/vis/drawVolume");
	  UI->ApplyCommand("/vis/viewer/set/viewpointVector -1 0 0");
	  UI->ApplyCommand("/vis/viewer/set/lightsVector -1 0 0");
	  UI->ApplyCommand("/vis/viewer/set/style wireframe");
	  //UI->ApplyCommand("/vis/viewer/set/style surface");
	  // UI->ApplyCommand("/vis/viewer/set/auxiliaryEdge true");
	  UI->ApplyCommand("/vis/viewer/set/lineSegmentsPerCircle 100");
	  UI->ApplyCommand("/tracking/storeTrajectory 1");
	  UI->ApplyCommand("/vis/scene/endOfEventAction accumulate");
	  UI->ApplyCommand("/vis/scent/endOfEventAction accumulate 2000"); // view more events!
	  UI->ApplyCommand("/vis/scene/add/trajectories");
	  UI->ApplyCommand("/vis/modeling/trajectories/create/drawByCharge");
	  UI->ApplyCommand("/vis/modeling/trajectories/drawByCharge-0/default/setDrawStepPts true");
	  UI->ApplyCommand("/vis/modeling/trajectories/drawByCharge-0/default/setStepPtsSize 2");
	  //UI->ApplyCommand("/vis/scene/add/hits");
	  UI->ApplyCommand("/vis/set/textColour green");
	  UI->ApplyCommand("/vis/set/textLayout right");
	  UI->ApplyCommand("/vis/scene/add/text2D 0.9 -.9 24 ! ! Demon Simulation"); 
	  UI->ApplyCommand("/vis/scene/add/text 0 10 35 cm 18 10 10 Det Module");
	  UI->ApplyCommand("/vis/scene/add/scale 25 cm");   //Simple scale line
	  //UI->ApplyCommand("/vis/scene/add/axes");    //Simple axes: x=red, y=green, z=blue.
	  UI->ApplyCommand("/vis/scene/add/logo2D");  //Simple logo
	  //UI->ApplyCommand("/vis/geometry/set/visibility Envelope 0 true");	
	  //UI->ApplyCommand("/vis/viewer/set/hiddenMarker true");

	  // Set viewing angle!
	  //UI->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 135 150");
	  UI->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 90 180");
	  //UI->ApplyCommand("/vis/viewer/set/viewpointThetaPhi 135 135");

          UI->ApplyCommand("/vis/viewer/pan -10 0 cm");

	  UI->ApplyCommand("/vis/viewer/zoom 0.5");  // Zoom in, > 1, Zoom out < 1	
	}
    }
  else
    {/*No Visualization! */}

  runManager->BeamOn(numberOfEvent);

 if(VisType == "OPENGL" && SaveVisFile == 1)
    {
    // print vis output
      UI->ApplyCommand("/vis/ogl/set/printMode pixmap"); // note: vectored mode crashes! 
      UI->ApplyCommand("/vis/ogl/set/printSize 900 600"); 
      UI->ApplyCommand("/vis/ogl/printEPS"); 
      G4cout << "OPENGL Visual Scene saved to G4OpenGL_x.eps. " << G4endl;
    }


  // Job termination
  //
  // Free the store: user actions, physics_list and detector_description are
  //                 owned and deleted by the run manager, so they should not
  //                 be deleted in the main() program !
  //
  G4cout
    << "############################################################" << G4endl;
  G4cout 
    << " Run Summary" << G4endl;
  G4cout
    << "############################################################" << G4endl;
  TntPointer->GetParticleTotals();
  delete TntPointer;
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete theEngine;
  delete runManager;

 cout << "End of Simulation! " << endl;

G4int end_time = time(NULL);
G4double run_time = static_cast<double>(end_time) - static_cast<double>(start_time);
  if(run_time < 60.)
    {G4cout << "The total run time was " << run_time << " seconds." << G4endl;}
  else
    {
      run_time /= 60.;
      G4cout << "The total run time was " << run_time << " minutes." << G4endl;
    }

 // Ascii and Root data files managed by TntDataRecordTree.hh class
 
  return 0;
}


