/// \file DetectorConstruction.hh
/// \brief Definition of detector construction class using GDML
//
#ifndef HAVE_TexanDETECTORCONSTRUCTION_H_
#define HAVE_TexanDETECTORCONSTRUCTION_H_

#include "G4VUserDetectorConstruction.hh"


class G4GDMLParser;
class G4UImessenger;

namespace tntsim {


/// Detector construction using the geometry read from a GDML (XML) file
class DetectorConstruction : public G4VUserDetectorConstruction
{
public: 
	/// Ctor
	DetectorConstruction(const G4String& file);
	/// Delete allocated stuff
	virtual ~DetectorConstruction();
	/// Construct the geometry
	virtual G4VPhysicalVolume *Construct();
	/// Construct sensitive volumes
	virtual void ConstructSDandField();
	/// Set GDML file
	void SetReadFile(const G4String& File);

private:
	G4VPhysicalVolume *fWorld;
	G4GDMLParser *fParser;
	G4String fReadFile;
	G4UImessenger* fMessenger;
};

}


#endif
