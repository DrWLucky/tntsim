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

#include "G4GenPhaseSpace.hh"
#include "TntDataRecordTree.hh"
#include "TntGlobalParams.hh"
#include "TntNeutronDecay.hh"
#include "TntNuclearMasses.hh"
#include "TntInputFileParser.hh"
#include "TntReaction.hh"

//by Shuya 160510. Just copied from Tntsim.
// need the below for random theta angle source (from Demon)
#include "Randomize.hh"
#include "G4UnitsTable.hh"
extern G4int Counter;


#if 0
namespace { 
void generate_from_7he(const G4double& ebeam,
											 G4double& px,
											 G4double& py,
											 G4double& pz, 
											 G4double& eneut,
											 G4LorentzVector& p_6he)
{
	// Generate neutron from breakup of 7he g.s.
	// Unbound state at 410 keV, width 150 keV

	// 7He @beam energy, 100% along beam axis
	const G4double he6mass = 5.60553446318; // GeV/c^2
	const G4double neutronmass = 9.39565378e-01; // GeV/c^2
	const G4double he7mass = he6mass + neutronmass + 410e-6; // S_b = -410 keV
	// \todo Add in width

	const G4double etot = he7mass + 7*ebeam/1e3; // total energy, GeV/c^2
	const G4double ptot = sqrt(etot*etot - he7mass*he7mass); // GeV/c

	// Decay products
	const G4double decayProductMasses[2] = { he6mass, neutronmass };

	static G4GenPhaseSpace* gen = 0;
	if(!gen) {
		gen = new G4GenPhaseSpace();
		G4LorentzVector p_7he(0, 0, ptot, etot);
		gen->SetDecay(p_7he, 2, decayProductMasses);
	}

	gen->Generate();

	G4LorentzVector* p_f = gen->GetDecay(0); // 6he fragment from 7he decay
	p_6he.set(p_f->px()*1e3, p_f->py()*1e3, p_f->pz()*1e3, p_f->e()*1e3);

	G4LorentzVector* p_neutron = gen->GetDecay(1); // neutron from 7he decay
	px = p_neutron->px() * 1e3;
	py = p_neutron->py() * 1e3;
	pz = p_neutron->pz() * 1e3;
	eneut = 1e3*(p_neutron->e() - p_neutron->m());

}


void generate_from_6he(const G4double& ebeam,
											 G4double& px1,
											 G4double& py1,
											 G4double& pz1, 
											 G4double& eneut1,
											 G4double& px2,
											 G4double& py2,
											 G4double& pz2, 
											 G4double& eneut2,
											 G4LorentzVector& p_6he)
{
	// Generate neutron from breakup of 6he first excited state @ 1797 keV

	// 7He @beam energy, 100% along beam axis
	const G4double he6mass = 5.60553446318 + 1797e-6; // GeV/c^2
	const G4double he4mass = 3.72737916179; // GeV/c^2
	const G4double neutronmass = 0.939565378; // GeV/c^2
	// \todo Add in width

	const G4double etot = he6mass + 6*ebeam/1e3; // total energy, GeV/c^2
	const G4double ptot = sqrt(etot*etot - he6mass*he6mass); // GeV/c

	// Decay products
	const G4double decayProductMasses[3] = { he4mass, neutronmass, neutronmass };

	static G4GenPhaseSpace* gen = 0;
	if(!gen) {
		gen = new G4GenPhaseSpace();
		G4LorentzVector p_6he1(0, 0, ptot, etot);
		bool okay = gen->SetDecay(p_6he1, 3, decayProductMasses);
		assert(okay);
	}

	while(1) {
		G4double w = gen->Generate();
		G4double r = G4UniformRand();
		if(r < w) { break; }
	}

	G4LorentzVector* p_f = gen->GetDecay(0); // 6he fragment from 7he decay
	p_6he.set(p_f->px()*1e3, p_f->py()*1e3, p_f->pz()*1e3, p_f->e()*1e3);

	G4LorentzVector* p_neutron1 = gen->GetDecay(1); // neutron from 6he decay
	px1 = p_neutron1->px() * 1e3;
	py1 = p_neutron1->py() * 1e3;
	pz1 = p_neutron1->pz() * 1e3;
	eneut1 = 1e3*(p_neutron1->e() - p_neutron1->m());

	G4LorentzVector* p_neutron2 = gen->GetDecay(2); // neutron from 6he decay
	px2 = p_neutron2->px() * 1e3;
	py2 = p_neutron2->py() * 1e3;
	pz2 = p_neutron2->pz() * 1e3;
	eneut2 = 1e3*(p_neutron2->e() - p_neutron2->m());
}

} // namespace
#endif



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
  fParticleGun->SetParticleDefinition(particleTable->FindParticle(particleName="neutron"));

	if(BeamType != "he7")
	{
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
		G4cout << "TntPrimaryGeneratorAction:: neutron energy " << fParticleGun->GetParticleEnergy() << G4endl;
	}
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
  //****************************************************
  // Now select type of beam
  //****************************************************
  G4double momentum_x = 0.;   // Default momentum is z-direction for
  G4double momentum_y = 0.;   // Pencil beams
  G4double momentum_z = 1.;
  
  const G4double Pi = CLHEP::pi;

//by Shuya 160510.
	G4double detector_thickness = TntGlobalParams::Instance()->GetDetectorZ();
  G4double beam_z = (-1*TntGlobalParams::Instance()->GetSourceZ() - (detector_thickness/2))*cm;
	
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
		posx *= (TntGlobalParams::Instance()->GetDetectorX()*cm);
		posy *= (TntGlobalParams::Instance()->GetDetectorX()*cm);
		fParticleGun->SetParticlePosition(G4ThreeVector(posx, posy, beam_z));
		// momentum_x = atan(posx/beam_z);
		// momentum_y = atan(posy/beam_z);
		// momentum_z = sqrt(1 - pow(momentum_x, 2) - pow(momentum_y, 2));
		G4ThreeVector v(momentum_x, momentum_y, momentum_z);
		fParticleGun->SetParticleMomentumDirection(v);
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
#if 0
	else if(BeamType == "he7")
	{
		G4double eneut = 0;
		G4LorentzVector p_6he;
    fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, beam_z)); //Default target posn (0,0)
		::generate_from_7he(TntGlobalParams::Instance()->GetNeutronEnergy(), //< beam energy
												momentum_x, momentum_y, momentum_z, eneut, p_6he);
		fParticleGun->SetParticleEnergy(eneut);
		TntDataOutPG->senddataPG(fParticleGun->GetParticleEnergy());
		TntDataOutPG->senddataSecondary(G4ThreeVector(0,0,beam_z), p_6he);
		
		G4cout << "TntPrimaryGeneratorAction:: neutron energy " <<
			fParticleGun->GetParticleEnergy() << G4endl;
	}
	else if(BeamType == "he6")
	{
		static TntNeutronDecay* decay = 0;
		if(!decay) {
			G4double m_6He = TntNuclearMasses::GetNuclearMass(2, 6)*MeV + 1.797*MeV; // EXCITED MASS
			decay = new TntTwoNeutronDecayPhaseSpace();
			G4double beamEnergy = TntGlobalParams::Instance()->GetNeutronEnergy() * 6; // MeV
			G4double beamMomentum = sqrt(pow(beamEnergy + m_6He, 2) - pow(m_6He, 2));
			decay->SetInitial(2, 6, G4LorentzVector(0, 0, beamMomentum, beamEnergy+m_6He));
			// decay->SetDecayParameter("energy", 1797*keV);
			// decay->SetDecayParameter("width", 0*keV); // actually 113 keV
		}
	
		static int whichNeutron = 0;

		fParticleGun->SetParticlePosition(G4ThreeVector(0.0, 0.0, beam_z)); //Default target posn (0,0)
		if(whichNeutron == 0) {
			decay->Generate();

			G4double en1 = decay->GetFinal(2).e() - decay->GetFinal(2).m();
			G4LorentzVector p_4he = decay->GetFinal(1);
			
			fParticleGun->SetParticleEnergy(en1);
			TntDataOutPG->senddataPG(fParticleGun->GetParticleEnergy());
			TntDataOutPG->senddataSecondary(G4ThreeVector(0,0,beam_z), p_4he);
			
			momentum_x = decay->GetFinal(2).px();
			momentum_y = decay->GetFinal(2).py();
			momentum_z = decay->GetFinal(2).pz();
 
			whichNeutron = 1;
		} else {
			G4double en2 = decay->GetFinal(3).e() - decay->GetFinal(3).m();
			G4LorentzVector p_4he = decay->GetFinal(1);
			
			fParticleGun->SetParticleEnergy(en2);
			TntDataOutPG->senddataPG(fParticleGun->GetParticleEnergy());
			TntDataOutPG->senddataSecondary(G4ThreeVector(0,0,beam_z), p_4he);
			
			momentum_x = decay->GetFinal(3).px();
			momentum_y = decay->GetFinal(3).py();
			momentum_z = decay->GetFinal(3).pz();
 
			whichNeutron = 0;
		}
		
		
		// if(whichNeutron == 0) {
		// 	generate_from_6he(TntGlobalParams::Instance()->GetNeutronEnergy(),
		// 										px1, py1, pz1, e1, px2,	py2, pz2, e2,
		// 										p_6he);
		// 	fParticleGun->SetParticleEnergy(e1);
		// 	TntDataOutPG->senddataPG(fParticleGun->GetParticleEnergy());
		// 	TntDataOutPG->senddataSecondary(G4ThreeVector(0,0,beam_z), p_6he);
		// 	momentum_x = px1;
		// 	momentum_y = py1;
		// 	momentum_z = pz1;
 
		// 	whichNeutron = 1;
		// } else {
		// 	fParticleGun->SetParticleEnergy(e2);
		// 	TntDataOutPG->senddataPG(fParticleGun->GetParticleEnergy());
		// 	TntDataOutPG->senddataSecondary(G4ThreeVector(0,0,beam_z), p_6he);
		// 	momentum_x = px2;
		// 	momentum_y = py2;
		// 	momentum_z = pz2;
 
		// 	whichNeutron = 0;
		// }
	}
#endif
	else 
	{
		G4cerr << "ERROR<TntPrimaryGeneratorAction.cc>:: Invalid BeamType: " << BeamType << G4endl;
		assert( 0 && "Invalid Beam Type" );

#ifdef NOT_COMPILED_
		
		// Attempt to reac from reaction file //
		G4String reacFile = BeamType;
		static TntReaction* reac = 0;
		static TntNeutronDecay* decay = 0;
		if(!reac) {
			// Set reaction parameters
			//
			TntReactionFactory factory;
			TntInputFileParser<TntReactionFactory> reacParser(&factory);
			reacParser.AddInput("beam",     &TntReactionFactory::SetBeam);
			reacParser.AddInput("target",   &TntReactionFactory::SetTarget);
			reacParser.AddInput("ejectile", &TntReactionFactory::SetEjectile);
			reacParser.AddInput("ebeam",    &TntReactionFactory::SetEbeamPerA);
			reacParser.AddInput("ex",       &TntReactionFactory::SetEx);
			reacParser.AddInput("width",    &TntReactionFactory::SetWidth);
			reacParser.AddInput("angdist",  &TntReactionFactory::SetAngDistFile);
			try 
			{ 
				reacParser.Parse(reacFile); 
			}
			catch (std::string s)  // NOT A VALID FILE
			{
				G4cerr << "ERROR<TntPrimaryGeneratorAction.cc>:: Invalid BeamType: " << BeamType << G4endl;
				assert( 0 && "Invalid Beam Type" );
			}

			reac = factory.CreateReaction();


			// Neutron decay parameters
			//
			TntNeutronDecayFactory decayFactory;
			TntInputFileParser<TntNeutronDecayFactory> decayParser(&decayFactory);
			decayParser.AddInput("decaytype", &TntNeutronDecayFactory::SetDecayType);
			decayParser.AddInput("decayoption", &TntNeutronDecayFactory::SetDecayOption);
			decayParser.Parse(reacFile);

			decay = decayFactory.Create();

			// Beam Parameters
			//
			{ // X //
				TntBeamEmittance *emX = new TntBeamEmittance();
				TntInputFileParser<TntBeamEmittance> parseEmX(emX);
				parseEmX.AddInput("embeamx", &TntBeamEmittance::SetTwist);
				parseEmX.AddInput("beamx", &TntBeamEmittance::SetX0);
				parseEmX.Parse(reacFile);
				if(emX->GetEpsilon()) { reac->SetEmittanceX(emX); } // takes ownership!
				else { delete emX; }
			}
			{ // Y //
				TntBeamEmittance *emY = new TntBeamEmittance();
				TntInputFileParser<TntBeamEmittance> parseEmY(emY);
				parseEmY.AddInput("embeamy", &TntBeamEmittance::SetTwist);
				parseEmY.AddInput("beamy", &TntBeamEmittance::SetX0);
				parseEmY.Parse(reacFile);
				if(emY->GetEpsilon()) { reac->SetEmittanceY(emY); } // takes ownership!
				else { delete emY; }
			}	
		}

		// Set up done, now generate event-by-event reaction
		//
		int nNeut = decay->GetNumberOfNeutrons();
		static G4int whichNeutron = 0;
		if(whichNeutron == 0) {
			// Now Generate reaction
			reac->Generate();

			// Neutron Decay
			decay->SetInitial(reac->GetZ4(), reac->GetA4(), reac->GetRecoil());
			decay->Generate();
		}
		
		// Save beam position
		G4ThreeVector beamPos(reac->GetBeamPosition().x(), reac->GetBeamPosition().y(), beam_z);
		fParticleGun->SetParticlePosition(beamPos);	

		// Set neutron energy+momentum
		G4double eNeut = decay->GetFinal(whichNeutron+2).e() - decay->GetFinal(whichNeutron+2).m();
		momentum_x = decay->GetFinal(whichNeutron+2).px();
		momentum_y = decay->GetFinal(whichNeutron+2).py();
		momentum_z = decay->GetFinal(whichNeutron+2).pz();

		// Recoil momentum (4-vector)
		G4LorentzVector p_recoil = decay->GetFinal(whichNeutron+1);

		// Send to data record class
		fParticleGun->SetParticleEnergy(eNeut);
		TntDataOutPG->senddataPG(fParticleGun->GetParticleEnergy());
		TntDataOutPG->senddataSecondary(beamPos, p_recoil);
		TntDataOutPG->senddataEjectile(beamPos,
																	 reac->GetEjectile(), 
																	 reac->GetThetaCM());
		TntDataOutPG->senddataReaction(beamPos,
																	 reac->GetBeam(),
																	 reac->GetTarget().m());
																	 

		if(whichNeutron == nNeut-1) { whichNeutron = 0; }
		else                        { ++whichNeutron;		}

#endif
	}
	
  // Set particle direction
  G4ThreeVector v(momentum_x,momentum_y,momentum_z);
  
//by Shuya 160510
  //fparticleGun->SetParticleMomentumDirection(v);
  fParticleGun->SetParticleMomentumDirection(v);

  fParticleGun->GeneratePrimaryVertex(anEvent);

	TntDataRecordTree::TntPointer->senddataPrimary(fParticleGun->GetParticlePosition(), 
																								 fParticleGun->GetParticleMomentumDirection());
}

// Utility class to parse reaction files
//
namespace { struct reac_file_params {
	G4String beam, target, ejectile;
	G4double ebeam, debeam;
	G4double ex, width;
	G4String angdist;
	G4double epsx, alphax, sigmax; // twist
	G4double epsy, alphay, sigmay;
	G4double x0,y0; // beam
	
	void set_beam(G4String b) { beam=b; }
	void set_target(G4String t) { target=t; }
	void set_ejectile(G4String e) { ejectile=e; }
	void set_ebeam(G4double e, G4double de) { ebeam=e; debeam=de; }
	void set_ex(G4double e) {ex=e;}
	void set_width(G4double w) {width=w;}
	void set_angdist(G4String a) {angdist=a;}
	void set_twist_x(G4double e, G4double a, G4double s)
		{ epsx=e; alphax=a; sigmax=s; }
	void set_twist_y(G4double e, G4double a, G4double s)
		{ epsy=e; alphay=a; sigmay=s; }
	void set_beam_x0(G4double x) { x0=x; }
	void set_beam_y0(G4double y) { y0=y; }
}; }
	
	

TntPGAReaction::TntPGAReaction():
	TntPrimaryGeneratorAction(),
	fReac(0), fDecay(0),
	fRngEbeam(0), fRngEx3(0), fRngEx4(0), fRngTheta(0), fRngPhi(0),
	fEmX(0), fEmY(0)
{
	// Parse reaction file
	//
	fReacFile = TntGlobalParams::Instance()->GetReacFile();

	reac_file_params rfp;
	TntInputFileParser<reac_file_params> parser(&rfp);
	parser.AddInput("beam",        &reac_file_params::set_beam);
	parser.AddInput("target",      &reac_file_params::set_target);
	parser.AddInput("ejectile",    &reac_file_params::set_ejectile);
	parser.AddInput("ebeam",       &reac_file_params::set_ebeam);
	parser.AddInput("ex",          &reac_file_params::set_ex);
	parser.AddInput("width",       &reac_file_params::set_width);
	parser.AddInput("angdist",     &reac_file_params::set_angdist);
	parser.AddInput("embeamx",     &reac_file_params::set_twist_x);
	parser.AddInput("embeamy",     &reac_file_params::set_twist_y);
	parser.AddInput("beamx",       &reac_file_params::set_beam_x0);
	parser.AddInput("beamy",       &reac_file_params::set_beam_y0);

	try { parser.Parse(fReacFile); }
	catch (std::string s) {
		G4cerr << "ERROR:: Invalid reaction file:: " << s << G4endl;
		exit(1);
	}

	// (+handle decay separately)
	TntNeutronDecayFactory decayFactory;
	TntInputFileParser<TntNeutronDecayFactory> decayParser(&decayFactory);
	decayParser.AddInput("decaytype", &TntNeutronDecayFactory::SetDecayType);
	decayParser.AddInput("decayoption", &TntNeutronDecayFactory::SetDecayOption);
	decayParser.Parse(fReacFile);


	// Initialize Generators
	//
	// Beam & Energy
	TntNuclearMasses::GetZAFromSymbol(rfp.beam, fZ, fA);
	TntNuclearMasses::GetZAFromSymbol(rfp.target, fZ+1, fA+1);
	TntNuclearMasses::GetZAFromSymbol(rfp.ejectile, fZ+2, fA+2);
	fZ[3] = fZ[0]+fZ[1] - fZ[2];
	fA[3] = fA[0]+fA[1] - fA[2];
	
	fRngEbeam.reset(new TntRngGaus(rfp.ebeam*fA[0], rfp.debeam*fA[0])); // beam energy

	// Reaction
	fReac.reset(new TntTwoBodyReactionGenerator());
	fReac->SetBeamTargetEjectile(rfp.beam, rfp.target, rfp.ejectile);
	// (+Emittance)
	fEmX.reset(new TntBeamEmittance(rfp.epsx, rfp.alphax, rfp.sigmax, rfp.x0));
	fEmY.reset(new TntBeamEmittance(rfp.epsy, rfp.alphay, rfp.sigmay, rfp.y0));						 
	fReac->SetEmittance(fEmX.get(), fEmY.get());
	// (+Angle Gen.)
	fRngTheta.reset(new TntRngCustomAngDist(rfp.angdist));
	fRngPhi.reset(new TntRngUniform(0, 2*CLHEP::pi));
	// fReac->SetRNGs():: DO LATER ::

	// Decay
	fDecay.reset(decayFactory.Create());
	if(typeid(&fDecay) == typeid(TntTwoNeutronDecayDiNeutron)) {
		fRngEx4.reset(new TntRngVolyaDiNeutronEx(rfp.ex, rfp.width, -18.7, fA[0]));
	}
	else {
		fRngEx4.reset(new TntRngBreitWigner(rfp.ex, rfp.width));
	}
	fReac->SetRNGs(fRngEbeam.get(), 0, fRngEx4.get(), fRngTheta.get(), fRngPhi.get());
	
}

TntPGAReaction::~TntPGAReaction()
{ }

void TntPGAReaction::GeneratePrimaries(G4Event* anEvent)
{
  //****************************************************
  // Now select type of beam
  //****************************************************
  G4double momentum_x = 0.;   // Default momentum is z-direction for
  G4double momentum_y = 0.;   // Pencil beams
  G4double momentum_z = 1.;
  
//by Shuya 160510.
	G4double detector_thickness = TntGlobalParams::Instance()->GetDetectorZ();
  G4double beam_z = (-1*TntGlobalParams::Instance()->GetSourceZ() - (detector_thickness/2))*cm;

	// Generate event-by-event reaction & neutron decay
	// Treat n>1 decays as separate 'events' (saved w/ same frag. data)
	// but neutron data from the corresponding neutron
	//
	int nNeut = fDecay->GetNumberOfNeutrons();
	static G4int whichNeutron = 0;
	if(whichNeutron == 0) {
		// We are on the first neutron, 
		// so generate a new reaction + decay
		//
		G4bool success = fReac->Generate();
		assert(success);
		
		// Neutron Decay
		fDecay->SetInputReaction(fReac.get());
		success = fDecay->Generate();
		assert(success);		
	}
		
	// Save beam position
	G4ThreeVector beamPos(fReac->GetPos().x(), fReac->GetPos().y(), beam_z);

	// Set neutron energy+momentum (for 'whichNeutron')
	// Offset in fDecay->GetFinal() is +2 (initial beam, fragment)
	//
	G4double eNeut = fDecay->GetFinal(whichNeutron+2).e() - fDecay->GetFinal(whichNeutron+2).m();
	momentum_x = fDecay->GetFinal(whichNeutron+2).px();
	momentum_y = fDecay->GetFinal(whichNeutron+2).py();
	momentum_z = fDecay->GetFinal(whichNeutron+2).pz();

	// Recoil ('fragment') momentum (4-vector)
	G4LorentzVector p_recoil = fDecay->GetFinal(whichNeutron+1);
	
	// Send to data record class
	TntDataOutPG->senddataPG(fParticleGun->GetParticleEnergy());
	TntDataOutPG->senddataSecondary(beamPos, p_recoil);
	TntDataOutPG->senddataEjectile(beamPos,
																 fReac->GetReactant(3),
																 fReac->GetThetaCM());
	TntDataOutPG->senddataReaction(beamPos,
																 fReac->GetReactant(1),
																 fReac->GetReactant(2).m());

	
	// Iterate through successive neutrons
	//
	if(whichNeutron == nNeut-1) { whichNeutron = 0; }
	else                        { ++whichNeutron;		}


	// Set particle gun paramters
	//
	// Energy and position
	fParticleGun->SetParticlePosition(beamPos);	
	fParticleGun->SetParticleEnergy(eNeut);

	// Direction
  G4ThreeVector v(momentum_x,momentum_y,momentum_z);
  fParticleGun->SetParticleMomentumDirection(v);

	// Generate event (neutron...)
  fParticleGun->GeneratePrimaryVertex(anEvent);
	TntDataRecordTree::TntPointer->senddataPrimary(fParticleGun->GetParticlePosition(), 
																								 fParticleGun->GetParticleMomentumDirection());
}
