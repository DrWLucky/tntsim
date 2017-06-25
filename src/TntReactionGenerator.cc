#include <cassert>
#include <vector>
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "G4GenPhaseSpace.hh"
#include "TntError.hh"
#include "TntReactionKinematics.hh"
#include "TntReactionGenerator.hh"


namespace {
template<typename T> T POW2(const T& t) { return t*t; }


// Helper function to generate beam + target distributions
// shared implementation common to different classes
void generate_beam_target(TntRng* fRngEbeam,
													TntParticle& fP1,
													TntParticle& fP2,
													TntBeamEmittance *fEmX,
													TntBeamEmittance *fEmY)
{
	// Generate incoming beam particle energy
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
}

} // namespace


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

void TntTwoBodyReactionGenerator::SetRNGs ( std::initializer_list<TntRng*> rngList )
{
	/** List items are:	rngEbeam, rngEx3, rngEx4, rngTheta, rngPhi
	 */
	int i=0;
	TntRng** fRng[5] =
		{ &fRngEbeam, &fRngEx3, &fRngEx4, &fRngTheta, &fRngPhi };
	for(const auto& rng : rngList) { *(fRng[i++]) = rng; }
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
	// Generate beam and target
	//
	generate_beam_target(fRngEbeam,fP1,fP2,fEmX,fEmY);

	// Generate Excitation Energies
	//
	G4double ex3 = fRngEx3 ? fRngEx3->GenerateAbove(0) : 0;
	fP3.SetEx(ex3);
	G4double ex4 = fRngEx4 ? fRngEx4->GenerateAbove(0) : 0;
	fP4.SetEx(ex4);
	
	// Generate CM angles
	//
	fTheta = fRngTheta ? fRngTheta->Generate() : 0;
	fPhi = fRngPhi ? fRngPhi->Generate() : 0;

	// Now calculate reaction
	G4double masses[2] = { fP3.MplusEx() , fP4.MplusEx() };
	TntTwoBodyReactionKinematics 
		reacKin(fP1.Momentum(), fP2.Momentum(), std::vector<G4double>(masses, masses+2));
	G4bool isGood = reacKin.Calculate(fTheta, fPhi);
	if(isGood) {
		fP3.SetPosition(fP1.Position());
		fP4.SetPosition(fP1.Position());
		fP3.SetP3(reacKin.GetProduct(0).vect());
		fP4.SetP3(reacKin.GetProduct(1).vect());
	} else {
		fP3.SetPosition(0,0,0);
		fP4.SetPosition(0,0,0);
		fP3.SetP3XYZ(0,0,0);
		fP4.SetP3XYZ(0,0,0);
	}
	
	return isGood;
}



//===============================================
//== TntNeutronPhaseSpaceReactionGenerator
//==

TntNeutronPhaseSpaceReactionGenerator::TntNeutronPhaseSpaceReactionGenerator(G4int n_neut):
	fNeutrons(n_neut), 
	fRngEbeam(0),
	fEmX(0), 
	fEmY(0)
{ fTheta=fPhi=0; }
	
TntNeutronPhaseSpaceReactionGenerator::~TntNeutronPhaseSpaceReactionGenerator()
{ }

void TntNeutronPhaseSpaceReactionGenerator::SetBeamTargetEjectile(const G4String& beam, 
																																	const G4String& target, 
																																	const G4String& ejectile)
{
	fP1.SetNucleus(beam);
	fP2.SetNucleus(target);
	fP3.SetNucleus(ejectile);
	G4int Z4 = fP1.Z() + fP2.Z() - fP3.Z();
	G4int A4 = fP1.A() + fP2.A() - fP3.A();
	fP4.SetNucleus(Z4, A4 - fNeutrons.size());
	for(TntParticle& P : fNeutrons) { P.SetNucleus(0, 1); }
}

void TntNeutronPhaseSpaceReactionGenerator::SetBeamTargetEjectile(G4int Zbeam, G4int Abeam, 
																																	G4int Ztrgt, G4int Atrgt,
																																	G4int Zejectile, G4int Aejectile)
{
	fP1.SetNucleus(Zbeam, Abeam);
	fP2.SetNucleus(Ztrgt, Atrgt);
	fP3.SetNucleus(Zejectile, Aejectile);
	G4int Z4 = fP1.Z() + fP2.Z() - fP3.Z();
	G4int A4 = fP1.A() + fP2.A() - fP3.A();
	fP4.SetNucleus(Z4, A4 - fNeutrons.size());
	for(TntParticle& P : fNeutrons) { P.SetNucleus(0, 1); }
}

void TntNeutronPhaseSpaceReactionGenerator::SetEmittance(TntBeamEmittance* emX, TntBeamEmittance* emY)
{
	fEmX = emX;
	fEmY = emY;
}

void TntNeutronPhaseSpaceReactionGenerator::SetRNGs( std::initializer_list<TntRng*> rngList )
{
	fRngEbeam = *rngList.begin();
}

const TntParticle& TntNeutronPhaseSpaceReactionGenerator::GetReactant(G4int i) const
{
	switch(i) {
	case 1: return fP1;
	case 2: return fP2;
	case 3: return fP3;
	case 4: return fP4;
	default:
		try { return fNeutrons.at(i-5); }
		catch (std::exception&)
		{
			G4cerr << "ERROR:: TntTwoBodyReactionGenerator::GetReactant:: Invalid indx: " << i
					 << ". Valid arguments are 1,2,3,4. Returning DUMMY vector: (0,0,0,0)" << G4endl;
			static TntParticle dummy;
			return dummy;
		}
	}
}


G4bool TntNeutronPhaseSpaceReactionGenerator::Generate()
{
	generate_beam_target(fRngEbeam, fP1, fP2, fEmX, fEmY);

	// Generate phase space
	//
	G4GenPhaseSpace gen;
	std::vector<G4double> finals = { fP3.M() , fP4.M() };
	for(const TntParticle& p : fNeutrons) { finals.push_back(p.M()); }
	G4bool isGood = gen.SetDecay(fP1.Momentum() + fP2.Momentum(), finals.size(), &finals[0]);
	if(!isGood) 
	{ // not enough energy for decay!
		fP3.SetPosition(0,0,0);
		fP4.SetPosition(0,0,0);
		fP3.SetP3XYZ(0,0,0);
		fP4.SetP3XYZ(0,0,0);
		for(TntParticle& P : fNeutrons) {
			P.SetPosition(0,0,0);
			P.SetP3XYZ(0,0,0);
		}
		TNTERR << "TntNeutronPhaseSpaceReactionGenerator::Generate() :: "
					 << "Not enough energy for decay!\nEnergies are (in | {out} | out sum):\n\t"
			;
		G4double sum = 0; for(G4double mf : finals) { sum += mf; }
		G4cerr << fP1.Momentum().e() << " | { ";
		for(G4double mf : finals) { G4cerr << mf << " "; }
		G4cerr << "} | " << sum << G4endl;
	}
	else
	{	
		G4double relwt, ran;
		do {
			relwt = gen.Generate() / gen.GetWtMax();
			ran = TntRngUniform().Generate();
		} while(relwt < ran);
	
		// Calculate CM angles
		//
		G4LorentzVector v3CM = *gen.GetDecay(0);
		G4ThreeVector boost3 = v3CM.boostVector();
		v3CM.boost(-boost3);
		fTheta = v3CM.theta();
		fPhi = v3CM.phi();

		// Set outputs
		//
		fP3.SetPosition(fP1.Position());
		fP4.SetPosition(fP1.Position());
		fP3.SetP3(gen.GetDecay(0)->vect());
		fP4.SetP3(gen.GetDecay(1)->vect());
		int iP = 2;
		for(TntParticle& P : fNeutrons) {
			P.SetPosition(fP1.Position());
			P.SetP3(gen.GetDecay(iP++)->vect());
		}
	}
	
	return isGood;
}
