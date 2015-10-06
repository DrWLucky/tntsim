/// \file TexanRunAction.cc
/// \brief Implementation of the TexanRunAction class

#include "TexanRunAction.hh"
#include "TexanPrimaryGeneratorAction.hh"
#include "TexanDetectorConstruction.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"
#include "TexanAnalysis.hh"


namespace txs = texansim;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::RunAction::RunAction():
	G4UserRunAction()
{
	txs::Ana()->G4()->SetVerboseLevel(1);
  txs::Ana()->G4()->SetFirstHistoId(1);

  // Creating histograms
  txs::Ana()->G4()->CreateH1("1","Edep in absorber", 100, 0., 800*MeV);
  txs::Ana()->G4()->CreateH1("2","Edep in gap", 100, 0., 100*MeV);
	txs::Ana()->G4()->CreateNtuple("t1", "Test ntuple");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::RunAction::~RunAction()
{
	delete txs::Ana();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4Run* txs::RunAction::GenerateRun()
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


void txs::RunAction::BeginOfRunAction(const G4Run*)
{
	/// This method is invoked before entering the event loop. A typical use of this
	/// method would be to initialize and/or book histograms for a particular run. 
	/// This method is invoked after the calculation of the physics tables.
  
  // Open an output file
	G4String fname =
		G4UImanager::GetUIpointer()->GetCurrentStringValue("/analysis/setFileName");
  txs::Ana()->G4()->OpenFile(fname);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::RunAction::EndOfRunAction(const G4Run*)
{
	/// This method is invoked at the very end of the run processing. It is typically
	/// used for a simple analysis of the processed run.

  // Save histograms
  txs::Ana()->G4()->Write();
  txs::Ana()->G4()->CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
