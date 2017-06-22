#include <cassert>
#include <vector>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "TntReactionKinematics.hh"
#include "TntReactionGenerator.hh"


namespace {
template<typename T> T POW2(const T& t) { return t*t; }
}


TntTwoBodyReactionGenerator::TntTwoBodyReactionGenerator():
	fRngEbeam(0), fRngEx3(0), fRngEx4(0), fRngTheta(0), fRngPhi(0),
	fEmX(0), fEmY(0)
{ fTheta=fPhi=0; }

TntTwoBodyReactionGenerator::~TntTwoBodyReactionGenerator()
{ }

void TntTwoBodyReactionGenerator::SetBeamTargetEjectile(const G4String& beam, 
																												const G4String& target, 
																												const G4String& ejectile)
{
	fP1.SetNucleus(beam);
	fP2.SetNucleus(target);
	fP3.SetNucleus(ejectile);
	G4int Z4 = fP1.Z() + fP2.Z() - fP3.Z();
	G4int A4 = fP1.A() + fP2.A() - fP3.A();
	fP4.SetNucleus(Z4, A4);
}

void TntTwoBodyReactionGenerator::SetBeamTargetEjectile(G4int Zbeam, G4int Abeam, 
																												G4int Ztrgt, G4int Atrgt,
																												G4int Zejectile, G4int Aejectile)
{
	fP1.SetNucleus(Zbeam, Abeam);
	fP2.SetNucleus(Ztrgt, Atrgt);
	fP3.SetNucleus(Zejectile, Aejectile);
	G4int Z4 = fP1.Z() + fP2.Z() - fP3.Z();
	G4int A4 = fP1.A() + fP2.A() - fP3.A();
	fP4.SetNucleus(Z4, A4);
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
	fRngEx3   = rngEx3;
	fRngEx4   = rngEx4;
	fRngTheta = rngTheta;
	fRngPhi   = rngPhi;
}

const TntParticle& TntTwoBodyReactionGenerator::GetReactant(G4int i) const
{
	switch(i) {
	case 1: return fP1;
	case 2: return fP2;
	case 3: return fP3;
	case 4: return fP4;
	default:
		{
			G4cerr << "ERROR:: TntTwoBodyReactionGenerator::GetReactant:: Invalid indx: " << i
					 << ". Valid arguments are 1,2,3,4. Returning DUMMY vector: (0,0,0,0)" << G4endl;
			static TntParticle dummy;
			return dummy;
		}
	}
}

G4bool TntTwoBodyReactionGenerator::Generate()
{
	// Generate incoming beam particle enegy
	//
	G4double ebeam = fRngEbeam->GenerateAbove(0); // beam kinetic energy
	G4double pbeam = sqrt( POW2(fP1.M()+ebeam) - fP1.M2() ); // beam momentum;
	
	// Generate Beam Angle & Position //
	//
	G4ThreeVector pos(0,0,0);
	G4ThreeVector pbeam3(0,0,pbeam);
	G4int iloop = 0; 
	for (auto* em : { fEmX, fEmY }) {
		if(em) {
			TntRngGaus2d rng(em->GetSigmaX(), em->GetSigmaTX(), em->GetRho());
			auto xtx = rng.Generate();
			G4double x  = xtx.first * mm  +  em->GetX0() * mm; // position
			pos[iloop] = x;

			G4double tx = xtx.second * mrad; // angle
			pbeam3[iloop] = pbeam*sin(tx/rad);
		}
		++iloop;
	}	
	pbeam3[2] = sqrt( POW2(pbeam) - POW2(pbeam3[0]) - POW2(pbeam3[1]) );
	fP1.SetPosition(pos);
	fP1.SetP3(pbeam3);
	assert( fabs(fP1.E() - fP1.M()) - ebeam < 1e-8 );
	assert( fabs(fP1.P() - pbeam) < 1e-8 );
	assert( fabs(fP1.Momentum().m() - fP1.M()) < 1e-8 );
	
	// Set target
	//
	fP2.SetP3XYZ(0,0,0);
	fP2.SetPosition(fP1.Position());
	
	// Generate Excitation Energies
	//
	G4double ex3 = fRngEx3 ? fRngEx3->GenerateAbove(0) : 0;
	G4double ex4 = fRngEx4 ? fRngEx4->GenerateAbove(0) : 0;

	// Generate CM angles
	//
	fTheta = fRngTheta ? fRngTheta->Generate() : 0;
	fPhi = fRngPhi ? fRngPhi->Generate() : 0;

	// Now calculate reaction
	G4double masses[2] = { fP3.M()+ex3, fP4.M()+ex4 };
	TntTwoBodyReactionKinematics 
		reacKin(fP1.Momentum(), fP2.Momentum(), std::vector<G4double>(masses, masses+2));
	G4bool isGood = reacKin.Calculate(fTheta, fPhi);
	if(isGood) {
		fP3.SetPosition(fP1.Position());
		fP4.SetPosition(fP1.Position());
		fP3.SetP3(reacKin.GetProduct(3).vect());
		fP4.SetP3(reacKin.GetProduct(4).vect());
	} else {
		fP3.SetPosition(0,0,0);
		fP4.SetPosition(0,0,0);
		fP3.SetP3XYZ(0,0,0);
		fP4.SetP3XYZ(0,0,0);
	}
	
	return isGood;
}
