/// \file DetectorConstruction.cc
/// \brief Implements detector construction (GDML) class.
#include "G4GDMLParser.hh"
#include "G4SDManager.hh"

#include "texansim/DetectorConstruction.hh"
#include "texansim/TexanSD.hh"
#include "texansim/Utils.hh"




namespace txs = texansim;

txs::DetectorConstruction::DetectorConstruction(const G4String& file)
{ 
	fReadFile  = file;
	fParser    = new G4GDMLParser();		
	fMessenger = new G4GDMLMessenger(fParser);
}

txs::DetectorConstruction::~DetectorConstruction()
{
	Zap(fParser);
	Zap(fMessenger);
}

void txs::DetectorConstruction::SetReadFile( const G4String& File )
{
  fReadFile=File;
}


G4VPhysicalVolume* txs::DetectorConstruction::Construct()
{
	fParser->Read(fReadFile);	
	fWorld = fParser->GetWorldVolume();
	return fWorld;
}

void txs::DetectorConstruction::ConstructSDandField()
{
	G4SDManager* SDman = G4SDManager::GetSDMpointer();
	txs::TexanSD* texanSD = new txs::TexanSD("TEXAN", "TexanHitsCollection");
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
