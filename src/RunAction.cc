/// \file RunAction.cc
/// \brief Implementation of the RunAction class
#include <cassert>

#include "tntsim/RunAction.hh"
#include "tntsim/PrimaryGeneratorAction.hh"
#include "tntsim/DetectorConstruction.hh"
#include "tntsim/Utils.hh"
#include "tntsim/Run.hh"

#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


namespace txs = tntsim;




txs::RunAction::RunAction():
	G4UserRunAction()
{ }


txs::RunAction::~RunAction()
{ }


G4Run* txs::RunAction::GenerateRun()
{
	// This method is invoked at the beginning of BeamOn. Because the user can inherit
  // the class G4Run and create his/her own concrete class to store some information about 
  // the run, the GenerateRun() method is the place to instantiate such an object. It is
  // also the ideal place to set variables which affect the physics table (such as 
  // production thresholds) for a particular run, because GenerateRun() is invoked before 
  // the calculation of the physics table.
	//
	/// Return an instance of our custom class tntsim::Run
	return new txs::Run();
}


void txs::RunAction::BeginOfRunAction(const G4Run*)
{
	// This method is invoked before entering the event loop. A typical use of this
	// method would be to initialize and/or book histograms for a particular run. 
	// This method is invoked after the calculation of the physics tables.
}



void txs::RunAction::EndOfRunAction(const G4Run*)
{
	// This method is invoked at the very end of the run processing. It is typically
	// used for a simple analysis of the processed run.
}

