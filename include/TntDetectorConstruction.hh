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
// $Id: TntDetectorConstruction.hh,v 1.6 2007/02/14 17:47:13 roeder Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// Modified by Brian Roeder LPC Caen on 02/14/2007
// email - roeder@lpccaen.in2p3.fr
// 
// --------------------------------------------------------------
// Based on     GEANT 4 - exampleN01, adding elements from other examples.
// --------------------------------------------------------------
//
//

#ifndef TntDetectorConstruction_H
#define TntDetectorConstruction_H 1

class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VSolid;
class G4Box;
class G4Tubs;
class G4Material;

#include "globals.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"

class TntDetectorConstruction : public G4VUserDetectorConstruction
{
  public:

  TntDetectorConstruction(G4String Material, G4String Type, G4ThreeVector Dims, G4String Light);
  ~TntDetectorConstruction();

    G4VPhysicalVolume* Construct();

  private:
  TntDetectorConstruction(){;}

  //World
    G4Box*             solidWorld;    // pointer to the solid envelope 
    G4LogicalVolume*   logicWorld;    // pointer to the logical envelope
    G4VPhysicalVolume* physiWorld;    // pointer to the physical envelope
  //
  // Elastic Scattering Target
  //  
    G4Box*             solidTarget;   // pointer to the solid Target
    G4LogicalVolume*   logicTarget;   // pointer to the logical Target
    G4VPhysicalVolume* physiTarget;   // pointer to the physical Target

  //
  // Silicon Detector
  //
    G4VSolid*          solidDetector;   // pointer to the solid Detector
    G4LogicalVolume*   logicDetector;   // pointer to the logical Detector
    G4VPhysicalVolume* physiDetector;   // pointer to the physical Detector

  //
  // Pointers - (See ExN02)
  //
  G4Material*         TargetMatter;  // pointer to the target  material

  // Class Variables
	G4String DetMaterial;
  G4String DetType;
  G4ThreeVector DetDim;
  G4String Light_Conv_Method;

};

#endif

