#include <vector>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "TntNuclearMasses.hh"
#include "TntReactionKinematics.hh"
#include "TntReactionGenerator.hh"


TntTwoBodyReactionGenerator::TntTwoBodyReactionGenerator():
	fBeam(0,0,0,0), fTarget(0,0,0,0), fEjectile(0,0,0,0), fRecoil(0,0,0,0),
	fPos(0,0,0),
	fRngEbeam(0), fRngEx3(0), fRngEx4(0), fRngTheta(0), fRngPhi(0),
	fEmX(0), fEmY(0)
{ fM1=fM2=fM3=fM4=fTheta=fPhi=0; }

TntTwoBodyReactionGenerator::~TntTwoBodyReactionGenerator()
{ }

void TntTwoBodyReactionGenerator::SetBeamTargetEjectile(const G4String& beam, 
																												const G4String& target, 
																												const G4String& ejectile)
{
	int z1,a1; TntNuclearMasses::GetZAFromSymbol(beam,&z1,&a1);
	int z2,a2; TntNuclearMasses::GetZAFromSymbol(target,&z2,&a2);
	int z3,a3; TntNuclearMasses::GetZAFromSymbol(ejectile,&z3,&a3);

	SetBeamTargetEjectile(z1,a1,z2,a2,z3,a3);
}

void TntTwoBodyReactionGenerator::SetBeamTargetEjectile(G4int Zbeam, G4int Abeam, 
																												G4int Ztrgt, G4int Atrgt,
																												G4int Zejectile, G4int Aejectile)
{
	G4int Zrecoil = Zbeam + Ztrgt - Zejectile;
	G4int Arecoil = Abeam + Atrgt - Aejectile;
	fM1 = TntNuclearMasses::GetNuclearMass(Zbeam,     Abeam)     * MeV;
	fM2 = TntNuclearMasses::GetNuclearMass(Ztrgt,     Atrgt)     * MeV;
	fM3 = TntNuclearMasses::GetNuclearMass(Zejectile, Aejectile) * MeV;
	fM4 = TntNuclearMasses::GetNuclearMass(Zrecoil,   Arecoil)   * MeV;
}

void TntTwoBodyReactionGenerator::SetEmittance(TntBeamEmittance* emX, TntBeamEmittance* emY)
{
	fEmX = emX;
	fEmY = emY;
}

void TntTwoBodyReactionGenerator::SetRNGs(TntRng* rngEbeam, 
																					TntRng* rngEx3, 
																					TntRng* rngEx4, 
																					TntRng* rngTheta,
																					TntRng* rngPhi)
{
	fRngEbeam = rngEbeam;
	fRngEx3 = rngEx3;
	fRngEx4 = rngEx4;
	fRngTheta = rngTheta;
	fRngPhi = rngPhi;
}


G4double TntTwoBodyReactionGenerator::GetReactantMass(G4int i) const
{
	switch(i) {
	case 1: return fM1;
	case 2: return fM2;
	case 3: return fM3;
	case 4: return fM4;
	default:
		{
			G4cerr << "ERROR:: TntTwoBodyReactionGenerator::GetReactantMass:: Invalid indx: " << i
					 << ". Valid arguments are 1,2,3,4. Returning DUMMY value: 0" << G4endl;
			return 0;
		}
	}
}
const G4LorentzVector& TntTwoBodyReactionGenerator::GetReactant(G4int i) const
{
	switch(i) {
	case 1: return fBeam;
	case 2: return fTarget;
	case 3: return fEjectile;
	case 4: return fRecoil;
	default:
		{
			G4cerr << "ERROR:: TntTwoBodyReactionGenerator::GetReactant:: Invalid indx: " << i
					 << ". Valid arguments are 1,2,3,4. Returning DUMMY vector: (0,0,0,0)" << G4endl;
			static G4LorentzVector dummy(0,0,0,0);
			return dummy;
		}
	}
}

G4bool TntTwoBodyReactionGenerator::Generate()
{
	// Generate incoming beam particle enegy
	//
	G4double ebeam = fRngEbeam->GenerateAbove(0); // beam kinetic energy
	G4double pbeam = sqrt(pow(fM1+ebeam, 2) - fM1*fM1); // beam momentum;
	
	// Generate Beam Angle & Position //
	//
	fPos.set(0,0,0);
	G4double pbeam_xyz[3] = {0,0,pbeam};
	G4int iloop = 0; 
	for (auto* em : { fEmX, fEmY }) {
		if(em) {
			TntRngGaus2d rng(em->GetSigmaX(), em->GetSigmaTX(), em->GetRho());
			auto xtx = rng.Generate();
			G4double x  = xtx.first * mm  +  em->GetX0() * mm; // position
			fPos[iloop] = x;

			G4double tx = xtx.second * mrad; // angle
			pbeam_xyz[iloop] = pbeam*sin(tx/rad);
		}
		++iloop;
	}	
	pbeam_xyz[2] = sqrt(pow(pbeam, 2) - pow(pbeam_xyz[0], 2) - pow(pbeam_xyz[1], 2));
	fBeam.set(pbeam_xyz[0], pbeam_xyz[1], pbeam_xyz[2], fM1+ebeam);

	// Set target
	//
	fTarget.set(0,0,0,fM2);
	
	// Generate Excitation Energies
	//
	G4double ex3 = fRngEx3 ? fRngEx3->GenerateAbove(0) : 0;
	G4double ex4 = fRngEx4 ? fRngEx4->GenerateAbove(0) : 0;

	// Generate CM angles
	//
	fTheta = fRngTheta ? fRngTheta->Generate() : 0;
	fPhi = fRngPhi ? fRngPhi->Generate() : 0;

	// Now calculate reaction
	G4double masses[2] = { fM3+ex3, fM4+ex4 };
	TntTwoBodyReactionKinematics reacKin(fBeam, fTarget, std::vector<G4double>(masses, masses+2));
	G4bool isGood = reacKin.Calculate(fTheta, fPhi);
	if(isGood) {
		fEjectile = reacKin.GetProduct(3);
		fRecoil = reacKin.GetProduct(4);
	} else {
		fEjectile.set(0,0,0,0);
		fRecoil.set(0,0,0,0);
	}
	
	return TRUE;
}
