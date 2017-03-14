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
// $Id: TntPrimaryGeneratorAction.cc 68752 2013-04-05 10:23:47Z gcosmo $
//
/// \file optical/Tnt/src/TntPrimaryGeneratorAction.cc
/// \brief Implementation of the TntPrimaryGeneratorAction class
//
//
#include <cassert>
#include "TntPrimaryGeneratorAction.hh"

#include "G4Event.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4SystemOfUnits.hh"
#include "globals.hh"

#include "TntDataRecordTree.hh"
#include "TntGlobalParams.hh"

//by Shuya 160510. Just copied from Tntsim.
// need the below for random theta angle source (from Demon)
#include "Randomize.hh"
#include "G4UnitsTable.hh"
extern G4int Counter;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntPrimaryGeneratorAction::TntPrimaryGeneratorAction(){
//by Shuya 160407.
  TntDataOutPG = TntDataRecordTree::TntPointer;

//by Shuya 160510
  // Set beamType for Primary Generator action
  //BeamType = "pencil";
  //BeamType = "diffuse";
  //BeamType = "conic";
	BeamType = TntGlobalParams::Instance()->GetBeamType();

  G4int n_particle = 1;
  fParticleGun = new G4ParticleGun(n_particle);
 
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
 
  G4String particleName;
  fParticleGun->SetParticleDefinition(particleTable->
//by Shuya 160404
                                     //FindParticle(particleName="gamma"));
																			FindParticle(particleName="neutron"));
	// 																		FindParticle(particleName="e-"));
//	 FindParticle(particleName="proton"));
                                     //FindParticle(particleName="deuteron"));
                                     //FindParticle(particleName="alpha"));
  //Default energy,position,momentum
//by Shuya 160404
  //fParticleGun->SetParticleEnergy(511.*keV);
//  fParticleGun->SetParticleEnergy(15.*MeV);
//by Shuya 160510
  //fParticleGun->SetParticleEnergy(1.*MeV);
	fParticleGun->SetParticleEnergy(TntGlobalParams::Instance()->GetNeutronEnergy());

	
//by Shuya 160510. I incorporated these in below GeneratePrimaries().
//Comments by Shuya 160427. This is a pencil beam. 
  //fParticleGun->SetParticlePosition(G4ThreeVector(0.0 , 0.0, -100.0*cm));
  //fParticleGun->SetParticleMomentumDirection(G4ThreeVector(0.,0.,1.));

  TntDataOutPG->senddataPG(fParticleGun->GetParticleEnergy());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntPrimaryGeneratorAction::~TntPrimaryGeneratorAction(){
    delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

//by Shuya 160427. I upgraded this as below so I can use different beam distributions.
//by Shuya 160509. I continued on upgrading it.
/*
void TntPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
*/

void TntPrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent){

//  dataevent++;
  
 // Energy Setup, ch_eng > 0->monoenergetic for "ch_eng" number of events, then increases by 1 MeV
 //               else (ch_eng <= 0) generates random energy from 1 to 100 MeV
 //               "ch_eng" variable defined in Tnt.cc "main" routine

////////////Comments by Shuya 160427. In the original Tnt code, below is not commented out.
/*

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

      // Generates random energy between 1 and 100 MeV
     //G4double ran1;
     //ran1=G4UniformRand(); //  random number between 0 and 1
     //particle_energy=1.+99.*ran1;   
     //particleGun->SetParticleEnergy(particle_energy);
    }
   
  TntDataOutPG->senddataPG(particle_energy);
 */

  //****************************************************
  // Now select type of beam
  //****************************************************
  G4double momentum_x = 0.;   // Default momentum is z-direction for
  G4double momentum_y = 0.;   // Pencil beams
  G4double momentum_z = 1.;
  
  const G4double Pi = CLHEP::pi;

//by Shuya 160510.
  G4double beam_z = -100.0*cm;
  
//by Shuya 160510
  //fparticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, -5.0*cm));
  fParticleGun->SetParticlePosition(G4ThreeVector(0.0 , 0.0, beam_z));

  if(BeamType == "pencil")
    {
      //position of source of particles
      //Straight Pencil Beam
      //particleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, -25.*cm));
      //by Shuya 160509
      fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, beam_z));
      //particleGun->SetParticlePosition(G4ThreeVector(0.0, -20.0*cm, 25.*cm));
      momentum_x = 0.;
      momentum_y = 0.;
      momentum_z = 1.;  
    }
	else if(BeamType == "rectangle")
	{
		momentum_x = 0.;
		momentum_y = 0.;
		momentum_z = 1.;
		G4double posx = G4UniformRand() - 0.5; // -0.5 -> 0.5
		G4double posy = G4UniformRand() - 0.5; // -0.5 -> 0.5
		posx *= 28; // -14 -> 14
		posy *= 28; // -14 -> 14
		fParticleGun->SetParticlePosition(G4ThreeVector(posx*cm, posy*cm, beam_z));									
	}
	else if(BeamType == "scan")
	{
		momentum_x = 0.;
		momentum_y = 0.;
		momentum_z = 1.;
		G4int ix = Counter/100;
		G4int iz = Counter%100;
		assert(ix<280);
		assert(iz<100);
		G4double posx = ix*0.1 - 14.; // cm
		G4double posz = iz*0.1 - 5.; // cm
		fParticleGun->SetParticlePosition(G4ThreeVector(posx*cm, 0., posz*cm));
		G4cout << "GENREATED EVENT:: (x,y) position = " << posx << " mm, " << posz << " mm" << G4endl;
	}
//COMMENT by Shuya 160510. Diffuse beam is parallel beam to z-axis but have random x-y (on theRadius circle) as a source position.
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

      G4double theRadius = 15*cm;
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
//by Shuya 160510
      //theStartPos[2] = -25.*cm;
      theStartPos[2] = beam_z;

	//by Shuya 160510
      //particleGun->SetParticlePosition(theStartPos);
      fParticleGun->SetParticlePosition(theStartPos);
    }
  else if(BeamType == "conic")   // Conic source (from Sega1, DEMON sims)
   {

    // 3D-conic source set x,y,z directions (see Demon for more info)
    // Convert to "Conic Source" (cone with bottom radius same as detector at front of detector)
    G4double OpenAngle=180.*deg; 
    G4double theta;
    G4double phi;
    
//Comment by Shuya 160510... This is inexplicable to me....
/*
    //Calculate aperature based on source distance from detector ...
    G4double dist = 5.0*m;
    G4double z_pos_det = 25.*cm-10.*cm;  // distance on z-axis-halfheight of module
    //G4double z_pos_det = 25.*cm-2.5*cm; // Eden test source pos.
    G4double z_pos = z_pos_det - dist;
*/

    //by Shuya 160510
    //fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, z_pos)); //Default 
    fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, beam_z)); //Default 
    //OpenAngle = atan(8.*cm/dist);
    //by Shuya 160510. Detector height is 100 cm, so height from z-axis is half of it (50cm).
    OpenAngle = atan(50.*cm/beam_z);

    theta = acos(1.+(cos(OpenAngle)-1.)*G4UniformRand());
  
    phi = 2.*Pi*G4UniformRand(); // flat

    momentum_x = sin(theta)*cos(phi);  // source in z-direction
    momentum_y = sin(theta)*sin(phi);
    momentum_z = cos(theta);
   }

  // Set particle direction
  G4ThreeVector v(momentum_x,momentum_y,momentum_z);
  
//by Shuya 160510
  //fparticleGun->SetParticleMomentumDirection(v);
  fParticleGun->SetParticleMomentumDirection(v);

  fParticleGun->GeneratePrimaryVertex(anEvent);
}
