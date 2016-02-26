//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
//
// $Id: TntDetectorConstruction.cc,v 1.9 2007/02/14 17:47:19 roeder Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// Modified by Brian Roeder LPC Caen on 02/14/2007
// email - roeder@lpccaen.in2p3.fr
// 
// -------------------------------------------------------------------
// Based on     GEANT 4 - exampleN01, adding elements from other exps.
// -------------------------------------------------------------------
//
// 12/13/2006 - Created Demon Scintillator Tube - Measures Energy Loss only
// 01/16/2007 - Adding Flux and Filter classes to differentiate between particles
// 02/14/2007 - Copied from Serra Sim. - switching to "Sensitive Detector" - 
//            - Tracking Detector definitions. 

#include "TntDetectorConstruction.hh"

#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "G4RotationMatrix.hh"
#include "globals.hh"
#include "G4ios.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4SDManager.hh"

// The methods below implement a user defined sensitive detector

#include "G4String.hh"

#include "TntSD.hh"            // User defined Sensitive Detector.
#include "G4VReadOutGeometry.hh"




namespace {

// calculate mass fraction from atomic ratios
inline void calcMassFraction(G4double atomicRatio, // ratio N1/N2 
														 G4double mass1, G4double mass2,
														 G4double* fraction1, G4double* fraction2)
{
	G4double af1 = atomicRatio / (1 + atomicRatio); // atomic fraction
	G4double af2 = 1. - af1; //atomic fraction

	G4double p1 = mass1*af1;
	G4double p2 = mass2*af2;

	*fraction1 = p1 / (p1 + p2);
	*fraction2 = p2 / (p1 + p2);
}

inline G4Material* createHydrocarbon(
	const char* name,
	G4double density,
	G4double HtoC_ratio,
	G4Element* pH,
	G4Element* pC)
{
	G4Material* hc = new G4Material(name, density, 2);
		
	G4double mf1, mf2;
	calcMassFraction(HtoC_ratio, pH->GetN(), pC->GetN(), &mf1, &mf2);
	
	hc->AddElement(pH, mf1);
	hc->AddElement(pC, mf2);

	return hc;
}

class MaterialNameIs {
	G4String name_;
public:
	MaterialNameIs(const G4String& name):
		name_(name) { }
	bool operator() (G4Material* m)
		{
			return m->GetName() == name_;
		}
};

G4Material* findMaterial(const G4String& name)
{
	G4MaterialTable* mt = G4Material::GetMaterialTable();
	G4MaterialTable::iterator it =
		std::find_if(mt->begin(), mt->end(), MaterialNameIs(name));

	return it != mt->end() ? *it : 0;
}
}

TntDetectorConstruction::TntDetectorConstruction(
	G4String Material, G4String Type,G4ThreeVector Dims,G4String Light)
	:solidWorld(0),  logicWorld(0),  physiWorld(0),
	 solidTarget(0), logicTarget(0), physiTarget(0),
	 solidDetector(0),logicDetector(0), physiDetector(0),
	 DetMaterial(Material),DetType(Type),DetDim(Dims),Light_Conv_Method(Light)
{;}

TntDetectorConstruction::~TntDetectorConstruction()
{ /* Destructor */ }

G4VPhysicalVolume* TntDetectorConstruction::Construct()
{

  // Elements and  materials used ---------------------------------------

  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density, temperature, pressure;
  G4int nel; // number of elements in compound

  G4Element* N  = new G4Element("Nitrogen", "N",  z=7.,  a= 14.01*g/mole);
  G4Element* O  = new G4Element("Oxygen"  , "O",  z=8.,  a= 16.00*g/mole);
  G4Element* H  = new G4Element("Hydrogen", "H" , z= 1., a= 1.01*g/mole);
  G4Element* C  = new G4Element("Carbon"  , "C" , z= 6., a= 12.011*g/mole); 

  // G4Material* Ar = 
  // new G4Material("ArgonGas", z= 18., a= 39.95*g/mole, density= 1.782*mg/cm3);

  // G4Material* Al = 
  //new G4Material("Aluminum", z= 13., a= 26.98*g/mole, density= 2.7*g/cm3);

  // G4Material* Pb = 
  // new G4Material("Lead", z= 82., a= 207.19*g/mole, density= 11.35*g/cm3);

  // Gold (for target material)
  // density = 19.3*g/cm3;
  // G4Material* Au = new G4Material("Au", z=81., a=196.97*g/mole, density);

  // Iron (for target material)
  // density = 7.87*g/cm3;
  // G4Material* Fe = new G4Material("Fe", z=26., a=55.845*g/mole, density);

	//  Vacuum

	density = 0.00000001*mg/cm3;
	G4Material* Vacuum = new G4Material("Vacuum", density, nel=2, kStateGas, 
																			temperature= 273.15*kelvin, pressure=0.0001*pascal);
	Vacuum->AddElement(N, .7);
	Vacuum->AddElement(O, .3);

// NE213
	density = 0.874*g/cm3 ;
	G4Material* NE213 = new G4Material("NE213", density, nel=2);
	NE213->AddElement(C, 90.779*perCent);    // New Calc w/ numbers from BC-501A
	NE213->AddElement(H, 9.221*perCent);
	//NE213->AddElement(C, 90.82*perCent);   // Marc's Calculation
	//NE213->AddElement(H, 9.18*perCent);

	{
		const G4double* numAtomsVect = NE213->GetVecNbOfAtomsPerVolume();
		G4double numC = numAtomsVect[0];
		G4double numH = numAtomsVect[1];
		G4cout << "NE-213: The Num. of C atoms is : " << numC << G4endl;
		G4cout << "NE-213: The Num. of H atoms is : " << numH << G4endl;
	}

	// Add all of the plastic & liquid scintillators
		
	G4Material* BC505 = 
		createHydrocarbon("BC505", 0.877*g/cm3, 1.331, H, C);

	// Print all the materials defined.
  //
  G4cout << G4endl << "The materials defined are : " << G4endl << G4endl;
  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  //World and Detector Volumes-----------------------------------------

  //World Volume - A cube of 10 m^3 
  //
  //Creates cubic world, note pointer "solidWorld" predefined in header
  //G4double WorldSize = 1.0*m;
  G4double WorldSize = 11.0*m;
  solidWorld = new G4Box("World",WorldSize/2.,WorldSize/2.,WorldSize/2.); 

  //World logic volume - filled with vacuum 
  logicWorld = new G4LogicalVolume(solidWorld,         //solid of logical World
																	 Vacuum,             //material filling logic volume
																	 "World",0,0,0);     //Name, other properties 0.

  //World Physical volume (placement)
  physiWorld = new G4PVPlacement(0,                // no rotation
																 G4ThreeVector(),  // position at 0,0,0
																 logicWorld,       // its logical volume
																 "World",          // its name
																 0,                // its mother volume;(0) for world
																 false,            // no boolean op.
																 0);               // not a copy
  //World Visualization (see Demon)
	G4VisAttributes* WorldVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.0,0.0));

	WorldVisAtt->SetForceWireframe(true);
	//logicWorld->SetVisAttributes(WorldVisAtt);
	logicWorld->SetVisAttributes (G4VisAttributes::Invisible);
    
  //
  // Detector Tube of NE213, placed along z-axis (Demon dimensions)
  //

	G4double innerRadius = 0.*cm;
	G4double startAngle = 0.*deg;
	G4double spanAngle = 360.*deg;
    
	G4double outerRadius;
	G4double halfHeight;
 
	if(DetType == "cylinder")
	{
		outerRadius = DetDim[0];
		halfHeight = DetDim[1];
		solidDetector = new G4Tubs("detector",innerRadius,outerRadius,halfHeight,startAngle,spanAngle);
	}
	else if(DetType == "bar")
	{
		G4double XDim = DetDim[0];
		G4double YDim = DetDim[1];
		G4double ZDim = DetDim[2];
		G4cout << "The Det. Dims were : " << XDim/cm << " cm, " << YDim/cm << " cm, " << ZDim/cm << " cm. " << G4endl;
		solidDetector = new G4Box("detector",XDim/2.,YDim/2.,ZDim/2.);
	}
	else 
	{ 
		G4cerr << "DetType unknown !!" << G4endl;
		G4cerr << "Please add type to Detector Construction!" << G4endl;
		G4ExceptionSeverity severity=FatalException;
		G4Exception("Program aborted in TntDetectorConstruction::Construct() method.", "Error",severity,"Error!");
		//Note: the enum G4ExceptionSeverity has been defined since Geant4 release 5.0 and its values are: FatalException, FatalErrorInArgument, EventMustBeAborted, FatalException, and JustWarning. 
		//G4Exception("Program aborted in TntDetectorConstruction::Construct() method.");
	}
 
 

  //Detector Logic Volume (tube of NE213 w/ no shell)
  //Eventually need to add shell as a separate logical volume for multi-detector array

	G4Material* detectorMaterial = findMaterial(DetMaterial);
	if(!detectorMaterial)
	{
		G4cerr << "Couldn't find detector material: \"" << DetMaterial << "\"" << G4endl;

		G4ExceptionSeverity severity=FatalException;
		G4Exception("Program aborted in TntDetectorConstruction::Construct() method.", "Error",severity,"Error!");
	}

	logicDetector = new G4LogicalVolume(solidDetector,         // Detector's solid
																			NE213,                 // It's Material
																			"detector",0,0,0);     // It's name, other props. 0

  //Detector Placement (at 0 degrees)
	G4ThreeVector positionDetector = G4ThreeVector(0.,0.,25.*cm);
	//G4RotationMatrix *rot1=new G4RotationMatrix();
	//rot1->rotateX(90.*deg);

	physiDetector = new G4PVPlacement(0,                  // rotated 0
																		positionDetector,  // position defined above
																		logicDetector,     // its logical volume				  
																		"detector",        // its name
																		logicWorld,      // its mother  volume
																		false,           // no boolean operations
																		0);              // copy number 
  
  // Detector Visualization (is it there?)
  // colour is r,g,b between 0 and 1.
  G4VisAttributes* DetectorVisAtt
    = new G4VisAttributes(G4Colour(0.0,0.5,0.0));  //Green Det.
  DetectorVisAtt->SetForceSolid(true);
  logicDetector->SetVisAttributes(DetectorVisAtt);

  // Let's try to make the detector sensitive

  // 2/14/07
  // User Defined Sensitive Detector Defs.
  // See Example ExN04, Demon

  // Create Sensitive Detector Manager

  G4SDManager* SDMan = G4SDManager::GetSDMpointer();

  // This follows from Scoring School Talks and Demon....

  G4String TntSDname = "TntSD";   // To access use "Tnt/TntSD"
  TntSD* theTntSD = new TntSD(TntSDname, Light_Conv_Method);
  SDMan->AddNewDetector( theTntSD );

  logicDetector->SetSensitiveDetector(theTntSD);      // Sets user defined detector as sensitive

  return physiWorld;    // Returns Detector Construction to RunManager

  // End of Detector Construction
 
}

