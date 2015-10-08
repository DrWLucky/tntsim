/// \file TexanDetectorMessenger.cc
/// \brief Implements messenger class to read detector files

#include "globals.hh"

#include "TexanDetectorMessenger.hh"
#include "TexanDetectorConstruction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithAnInteger.hh"

namespace txs = texansim;


txs::DetectorMessenger::DetectorMessenger( txs::DetectorConstruction* myDet )
  : G4UImessenger(),
    fTheDetector( myDet ),
    fTheDetectorDir(0),
    fTheReadCommand(0)
{ 
  fTheDetectorDir = new G4UIdirectory( "/texan/detector/" );
  fTheDetectorDir->SetGuidance("Detector control.");

  fTheReadCommand = new G4UIcmdWithAString("/texan/detector/readFile", this);
  fTheReadCommand ->SetGuidance("READ GDML file with given name");
  fTheReadCommand ->SetParameterName("FileRead", false);
  fTheReadCommand ->SetDefaultValue("detectors.gdml");
  fTheReadCommand ->AvailableForStates(G4State_PreInit);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

txs::DetectorMessenger::~DetectorMessenger()
{
  delete fTheReadCommand;
  delete fTheDetectorDir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void txs::DetectorMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{ 
  if ( command == fTheReadCommand )
  { 
    fTheDetector->SetReadFile(newValue );
  }
}
