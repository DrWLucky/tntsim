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
// $Id: TntPrimaryGeneratorAction.cc,v 1.6 2006/12/13 17:47:21 roeder Exp $
// GEANT4 tag $Name: geant4-08-01-patch-02 $
//
// Modified by Brian Roeder LPC Caen on 12/13/2006
// email - roeder@lpccaen.in2p3.fr
// 
// -------------------------------------------------------------------
// Based on     GEANT 4 - exampleN01, adding elements from other exps.
// -------------------------------------------------------------------
// 
// 12/13/06 - Added code in "Generate Primaries" to generate a random 
//            energy particle and save the initial energy in a file
// 12/13/06 - Pencil Beam of neutrons at random energy in x-direction
//
// 09/05/07 - Random number generator now defined in function randnum() 
//            in header "randCLHEP.h"
//  

#include "TntPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "globals.hh"

// need the below for random theta angle source (from Demon)
#include "Randomize.hh"
#include "G4UnitsTable.hh"

// 19 Jan 07 - C++ data io goes to TntDataRecordTree class 
// through "TntDataOut" pointer (assigned to "TntPointer")


TntPrimaryGeneratorAction::TntPrimaryGeneratorAction(G4int eng_step, G4String beamType, G4double Eng) : 
  ch_eng(eng_step), BeamType(beamType)
{
  TntDataOutPG = TntDataRecordTree::TntPointer;
  G4int n_particle = 1;
  particleGun = new G4ParticleGun(n_particle);

  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  //particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="geantino"));
  particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="neutron"));
  //particleGun->SetParticleDefinition(particleTable->FindParticle(particleName="proton"));

  particle_energy = Eng;
  dataevent = 0;
}

TntPrimaryGeneratorAction::~TntPrimaryGeneratorAction()
{
  delete particleGun;
}

void TntPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{
  dataevent++;
  
 // Energy Setup, ch_eng > 0->monoenergetic for "ch_eng" number of events, then increases by 1 MeV
 //               else (ch_eng <= 0) generates random energy from 1 to 100 MeV
 //               "ch_eng" variable defined in Tnt.cc "main" routine

  if (ch_eng > 0) 
    {
      if ((dataevent) % ch_eng == 0)
	{
	  particle_energy++;
	  //particle_energy += 0.1*MeV;
          // Calculate Detection efficiency at this energy.
	  cout << endl;
	  TntDataOutPG->CalculateEff(ch_eng); // Calculate Detection efficiency at this energy.
          cout << "****************************************************************" << endl;
          G4cout << "Energy is now: "<< particle_energy << "!" << endl;
	  cout << "****************************************************************" << endl;
	  dataevent = 0;
        }
    particleGun->SetParticleEnergy(particle_energy);  // particle_energy starts at 1*MeV; (in constructor)
    }
  else
    {
     // Use to set to a constant energy for the whole sim.
     // Value set in main() program 
     particleGun->SetParticleEnergy(particle_energy);
      /*
      // Generates random energy between 1 and 100 MeV
     G4double ran1;
     ran1=G4UniformRand(); //  random number between 0 and 1
     particle_energy=1.+99.*ran1;   
     particleGun->SetParticleEnergy(particle_energy);
      */
    }
   
  TntDataOutPG->senddataPG(particle_energy);

  //****************************************************
  // Now select type of beam
  //****************************************************

  G4double momentum_x = 0.;   // Default momentum is z-direction for
  G4double momentum_y = 0.;   // Pencil beams
  G4double momentum_z = 1.;
  
  const G4double Pi = CLHEP::pi;
  
  particleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, -5.0*cm));

  if(BeamType == "pencil")
    {
      //position of source of particles
      //Straight Pencil Beam
      particleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, -25.*cm));
      //particleGun->SetParticlePosition(G4ThreeVector(0.0, -20.0*cm, 25.*cm));
      momentum_x = 0.;
      momentum_y = 0.;
      momentum_z = 1.;  
    }
  else if(BeamType == "diffuse")
    {
      // Edited to give a "true" diffuse beam spot - 11/01/08 BTR
      // Follows method of "BeamPosDet" simulation      
      //Straight Pencil Beam with a diffuse beam spot. 
      //Particle starts at a random position
      //in x-y plane.
      momentum_x = 0.;
      momentum_y = 0.;
      momentum_z = 1.;

      G4double theRadius = 8.0*cm;
      theRadius *= 22.5;          // correction factor 
    
      // Since diffuse beam is forward focused in small theta,
      // Need Pi-small number shaped like cosine to get circular beam

     G4double ran = G4UniformRand();
     G4double theta = acos(1.-0.001*ran);

     ran = G4UniformRand();
     G4double phi = twopi*ran; // Flat in phi

     G4double x_pos = cos(phi)*sin(theta);
     G4double y_pos = sin(phi)*sin(theta);

     G4ThreeVector theStartPos(0.,0.,0.);
      theStartPos[0] = theRadius*x_pos;
      theStartPos[1] = theRadius*y_pos;
      theStartPos[2] = -25.*cm;

      particleGun->SetParticlePosition(theStartPos);
    }
  else if(BeamType == "conic")   // Conic source (from Sega1, DEMON sims)
   {

    // 3D-conic source set x,y,z directions (see Demon for more info)
    // Convert to "Conic Source" (cone with bottom radius same as detector at front of detector)
    G4double OpenAngle=180.*deg; 
    G4double theta;
    G4double phi;
    
    //Calculate aperature based on source distance from detector ...
    G4double dist = 5.0*m;
    G4double z_pos_det = 25.*cm-10.*cm;  // distance on z-axis-halfheight of module
    //G4double z_pos_det = 25.*cm-2.5*cm; // Eden test source pos.
    G4double z_pos = z_pos_det - dist;

    particleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, z_pos)); //Default 
    OpenAngle = atan(8.*cm/dist);

    theta = acos(1.+(cos(OpenAngle)-1.)*G4UniformRand());
  
    phi = 2.*Pi*G4UniformRand(); // flat

    momentum_x = sin(theta)*cos(phi);  // source in z-direction
    momentum_y = sin(theta)*sin(phi);
    momentum_z = cos(theta);
   }

  // Set particle direction
  G4ThreeVector v(momentum_x,momentum_y,momentum_z);
  
  particleGun->SetParticleMomentumDirection(v);
  particleGun->GeneratePrimaryVertex(anEvent);
}


