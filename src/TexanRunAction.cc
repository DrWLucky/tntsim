/// \file TexanRunAction.cc
/// \brief Implementation of the TexanRunAction class

#include "TexanRunAction.hh"
#include "TexanPrimaryGeneratorAction.hh"
#include "TexanDetectorConstruction.hh"
#include "TexanRootAnalyzer.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include <TTree.h>


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TexanRunAction::TexanRunAction()
: G4UserRunAction()
{
	TexanRootAnalyzer::Instance()->OpenFile("output.root", TexanRootAnalyzer::kAlways);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TexanRunAction::~TexanRunAction()
{
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* TexanRunAction::GenerateRun()
{
	/// This method is invoked at the beginning of BeamOn. Because the user can inherit
  /// the class G4Run and create his/her own concrete class to store some information about 
  /// the run, the GenerateRun() method is the place to instantiate such an object. It is
  /// also the ideal place to set variables which affect the physics table (such as 
  /// production thresholds) for a particular run, because GenerateRun() is invoked before 
  /// the calculation of the physics table.
	return new G4Run;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace { double val; }
void TexanRunAction::BeginOfRunAction(const G4Run*)
{
	/// This method is invoked before entering the event loop. A typical use of this
	/// method would be to initialize and/or book histograms for a particular run. 
	/// This method is invoked after the calculation of the physics tables.

	TTree* t1 = TexanRootAnalyzer::Instance()->CreateTree("t1", "t1");
	t1->Branch("val", &val, "val/D");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TexanRunAction::EndOfRunAction(const G4Run*)
{
	/// This method is invoked at the very end of the run processing. It is typically
	/// used for a simple analysis of the processed run.
	
	TexanRootAnalyzer::Instance()->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
