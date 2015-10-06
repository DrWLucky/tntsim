/// \file TexanDetectorConstruction.hh
/// \brief Definition of detector construction class using GDML
//
#ifndef HAVE_TexanDETECTORCONSTRUCTION_H_
#define HAVE_TexanDETECTORCONSTRUCTION_H_


#include "G4VUserDetectorConstruction.hh"


namespace texansim {

/// Detector construction using the geometry read from a GDML (XML) file
class DetectorConstruction : public G4VUserDetectorConstruction
{
public: 
	DetectorConstruction(G4VPhysicalVolume *setWorld = 0);
	virtual G4VPhysicalVolume *Construct();

private:
	G4VPhysicalVolume *fWorld;
};

}


#endif
