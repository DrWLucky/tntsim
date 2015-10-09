/// \file TexanDetectorConstruction.cc
/// \brief Implements detector construction (GDML) class.
#include "G4GDMLParser.hh"
#include "TexanDetectorConstruction.hh"
#include "TexanDetectorMessenger.hh"
#include "TexanUtils.hh"

namespace txs = texansim;

txs::DetectorConstruction::DetectorConstruction()
{ 
	fReadFile  = TEXAN_BUILD_DIR + G4String("/detectors.gdml");
	fParser    = new G4GDMLParser();		
	fMessenger = new G4GDMLMessenger(fParser);
}

txs::DetectorConstruction::~DetectorConstruction()
{
	Zap(fParser);
	Zap(fMessenger);
}

G4VPhysicalVolume* txs::DetectorConstruction::Construct()
{
	fParser->Read(fReadFile);	
	fWorld = fParser->GetWorldVolume();
	return fWorld;
}

void txs::DetectorConstruction::SetReadFile( const G4String& File )
{
  fReadFile=File;
}
