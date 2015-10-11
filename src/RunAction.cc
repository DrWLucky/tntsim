/// \file RunAction.cc
/// \brief Implementation of the RunAction class
#include <cassert>

#include "texansim/Analysis.hh"
#include "texansim/RunAction.hh"
#include "texansim/PrimaryGeneratorAction.hh"
#include "texansim/DetectorConstruction.hh"
#include "texansim/Utils.hh"

#include "G4Run.hh"
#include "G4UImanager.hh"
#include "G4RunManager.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4LogicalVolume.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"


namespace txs = texansim;


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::RunAction::RunAction():
	G4UserRunAction()
{
	/// Set up all of the analysis stuff
	
	/// Default file name (overwrite with macro /analysis/setFileName)
	Analysis::SetFileName("texansim_output");

  /// Create ntuple and rows
	Analysis::CreateNtuple("t", "TEXAN Geant4 simulation event data");

	Analysis::BookNtupleColumn<G4int>("fNumHits");
	for(G4int i=0; i< TXS_MAX_HITS; ++i) {
		Analysis::BookNtupleColumn<G4double>(FormatStr1<G4int>("fEdep", i));
	}
	Analysis::FinishNtuple();

	/// Create histograms
	Analysis::BookH1("hst1", "", 100, 0, 10, "MeV");
	Analysis::BookH2("hst2", "", 100, 0, 10, 100, 0, 10, "MeV");
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::RunAction::~RunAction()
{
	/// Clean up analysis stuff
	Analysis::Cleanup();
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

	/// Open analysis file
	/// Default file name is set in the constructor.
	/// Can be overwritten with the macro /analysis/setFileName
	Analysis::OpenFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::RunAction::EndOfRunAction(const G4Run*)
{
	/// This method is invoked at the very end of the run processing. It is typically
	/// used for a simple analysis of the processed run.

  /// Write analysis stuff
	Analysis::Write();
	Analysis::CloseFile();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
