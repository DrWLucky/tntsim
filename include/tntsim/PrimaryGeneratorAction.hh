/// \file PrimaryGeneratorAction.hh
/// \brief Primary generator action
///
#ifndef TexanPrimaryGeneratorAction_h
#define TexanPrimaryGeneratorAction_h 1

#include "G4SystemOfUnits.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"

class G4ParticleGun;
class G4Event;


namespace tntsim {

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
public:
	PrimaryGeneratorAction(
		const G4String& particleName,
		G4double energy,
		G4ThreeVector position,
		G4ThreeVector momentumDirection);
	~PrimaryGeneratorAction();

	// methods
	virtual void GeneratePrimaries(G4Event*);

private:
	// data members
	G4ParticleGun*  fParticleGun; //pointer a to G4 service class
};

}


#endif
