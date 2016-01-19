/// \file DetectorConstruction.cc
/// \brief Implements detector construction (GDML) class.
#include "G4GDMLParser.hh"
#include "G4SDManager.hh"

#include "tntsim/DetectorConstruction.hh"
#include "tntsim/TexanSD.hh"
#include "tntsim/Utils.hh"




namespace tnt = tntsim;

tnt::DetectorConstruction::DetectorConstruction(const G4String& file)
{ 
	fReadFile  = file;
	fParser    = new G4GDMLParser();		
	fMessenger = new G4GDMLMessenger(fParser);
}

tnt::DetectorConstruction::~DetectorConstruction()
{
	Zap(fParser);
	Zap(fMessenger);
}

void tnt::DetectorConstruction::SetReadFile( const G4String& File )
{
  fReadFile=File;
}


G4VPhysicalVolume* tnt::DetectorConstruction::Construct()
{
	fParser->Read(fReadFile);	
	fWorld = fParser->GetWorldVolume();
	return fWorld;
}

void tnt::DetectorConstruction::ConstructSDandField()
{
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	tnt::TexanSD* texanSD = new tnt::TexanSD("TEXAN", "TexanHitsCollection");
	SDman->AddNewDetector( texanSD );

	const G4GDMLAuxMapType* auxmap = fParser->GetAuxMap();
	G4cout << "Found " << auxmap->size()
				 << " volume(s) with auxiliary information."
				 << G4endl << G4endl;
	for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin();
			iter!=auxmap->end(); iter++) 
	{
		G4cout << "Volume " << ((*iter).first)->GetName()
					 << " has the following list of auxiliary information: "
					 << G4endl << G4endl;
		for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin();
				 vit!=(*iter).second.end(); vit++)
		{
			G4cout << "--> Type: " << (*vit).type
						 << " Value: " << (*vit).value << G4endl;
		}
	}
	G4cout << G4endl;


	// The same as above, but now we are looking for
	// sensitive detectors setting them for the volumes

	for(G4GDMLAuxMapType::const_iterator iter=auxmap->begin();
			iter!=auxmap->end(); iter++) 
	{
		G4cout << "Volume " << ((*iter).first)->GetName()
					 << " has the following list of auxiliary information: "
					 << G4endl << G4endl;
		for (G4GDMLAuxListType::const_iterator vit=(*iter).second.begin();
				 vit!=(*iter).second.end();vit++)
		{
			if ((*vit).type=="SensDet")
			{
				G4cout << "Attaching sensitive detector " << (*vit).value
							 << " to volume " << ((*iter).first)->GetName()
							 <<  G4endl << G4endl;

				G4VSensitiveDetector* mydet = 
					SDman->FindSensitiveDetector((*vit).value);
				if(mydet) 
				{
					G4LogicalVolume* myvol = (*iter).first;
					myvol->SetSensitiveDetector(mydet);
				}
				else
				{
					G4cout << (*vit).value << " detector not found" << G4endl;
				}
			}
		}
	}
}
