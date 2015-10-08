/// \file TexanDetectorMessenger.hh
/// \brief Defines messenger class to read detector files

#ifndef TexanDetectorMessenger_h
#define TexanDetectorMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithAnInteger;

// ----------------------------------------------------------------------------


namespace texansim {

class DetectorConstruction;


/// Detector messenger class used to set GDML file
class DetectorMessenger: public G4UImessenger
{
public:
	DetectorMessenger( texansim::DetectorConstruction* );
	~DetectorMessenger();
    
	virtual void SetNewValue( G4UIcommand*, G4String );

private:
	DetectorConstruction*      fTheDetector;
	G4UIdirectory*             fTheDetectorDir;
	G4UIcmdWithAString*        fTheReadCommand;
};

}

// ----------------------------------------------------------------------------

#endif
