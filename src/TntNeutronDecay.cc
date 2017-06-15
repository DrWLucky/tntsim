#include <map>
#include <CLHEP/Random/RandBreitWigner.h>
#include <G4SystemOfUnits.hh>
#include <Randomize.hh>
#include "G4GenPhaseSpace.hh"
#include "TntNuclearMasses.hh"
#include "TntNeutronDecay.hh"

namespace {

const G4double kNeutronMass =
	TntNuclearMasses::GetNuclearMass(0, 1)*MeV;

} //namespace


/////////////////////////////////////////////////////////////////
// TntNeutronDecayIntermediate
//

TntNeutronDecayIntermediate::TntNeutronDecayIntermediate(G4int number_of_neutrons_emitted):
	mNumberOfNeutronsEmitted(number_of_neutrons_emitted),
	mFinal(number_of_neutrons_emitted + 2),
	mFragMass(0),
	mInitial(0,0,0,0)
{ }

TntNeutronDecayIntermediate::~TntNeutronDecayIntermediate()
{ }

void TntNeutronDecayIntermediate::SetInitial(G4int Z, G4int A, const G4LorentzVector& momentum)
{
	mBeamMass = TntNuclearMasses::GetNuclearMass(Z, A)*MeV;
	mFragMass = TntNuclearMasses::GetNuclearMass(Z, A - mNumberOfNeutronsEmitted)*MeV;
	mInitial = momentum;
}

void TntNeutronDecayIntermediate::SetFinal(G4int indx, const G4LorentzVector& v)
{
	if(size_t(indx) < mFinal.size())	{
		mFinal[indx].set(v.x(), v.y(), v.z(), v.t());
	}
	else { 		
		G4cerr << "ERROR:: Invalid index to TntNeutronDecayIntermediate::SetFinal:: " << indx << G4endl;
		mFinal.at(indx).set(0,0,0,0); // throws 
	}
}

const G4LorentzVector& TntNeutronDecayIntermediate::GetFinal(G4int indx) const
{
	if(size_t(indx) < mFinal.size())	{
		return mFinal[indx]; 
	}
	else {
		G4cerr << "ERROR:: Invalid index to TntNeutronDecayIntermediate::GetFinal:: " << indx << G4endl;
		return mFinal.at(indx); // throws
	}
}

namespace {
G4double do_weighted_generation(bool uniform, G4GenPhaseSpace& gen)
{
	G4double weight = gen.Generate();
	if(uniform) {
		while(weight < G4UniformRand()) {
			weight = gen.Generate();
		}
		weight = 1;
	}
	return weight;
} }


/////////////////////////////////////////////////////////////////
// TntOneNeutronDecay
//

TntOneNeutronDecay::TntOneNeutronDecay():
	TntNeutronDecayIntermediate(1)
{ }

TntOneNeutronDecay::~TntOneNeutronDecay()
{ }


G4double TntOneNeutronDecay::Generate(G4bool uniformWeight)
{
	G4double mOut[2] = { mFragMass, kNeutronMass };

	if(mInitial.m() < mFragMass + kNeutronMass) {
		G4cerr << "ERROR:: TntTwoNeutronDecayPhaseSpace:: Not enough energy for decay!" << G4endl;
		G4cerr << "MASSES (MBeam+ex, MF, MN):: " << mInitial.m() << ", "
					 << mOut[0] << ", " << mOut[1] << G4endl;
		throw false;
	}
	SetFinal(0, mInitial); // Initial excited nucleus

	G4GenPhaseSpace gen;
	G4bool possible = gen.SetDecay(mInitial, 2, mOut);
	if(!possible) {
		G4cerr << "ERROR:: TntOneNeutronDecay:: Not enough energy for decay!" << G4endl;
		G4cerr << "MASSES (MBeam+ex, MF, MN):: " << mInitial.m() << ", "
					 << mOut[0] << ", " << mOut[1] << G4endl;
		throw possible;
	}

	G4double weight = do_weighted_generation(uniformWeight, gen);
	SetFinal(1, *(gen.GetDecay(0))); // fragment
	SetFinal(2, *(gen.GetDecay(1))); // neutron
	
	return weight;
}



/////////////////////////////////////////////////////////////////
// TntTwoNeutronDecayPhaseSpace
//

TntTwoNeutronDecayPhaseSpace::TntTwoNeutronDecayPhaseSpace():
	TntNeutronDecayIntermediate(2)
{ }

TntTwoNeutronDecayPhaseSpace::~TntTwoNeutronDecayPhaseSpace()
{ }


G4double TntTwoNeutronDecayPhaseSpace::Generate(G4bool uniformWeight)
{
	G4double mOut[3] = { mFragMass, kNeutronMass, kNeutronMass };
	if(mInitial.m() < mFragMass + 2*kNeutronMass) {
		G4cerr << "ERROR:: TntTwoNeutronDecayPhaseSpace:: Not enough energy for decay!" << G4endl;
		G4cerr << "MASSES (MBeam+ex, MF, 2*MN):: " << mInitial.m() << ", "
					 << mOut[0] << ", " << 2*mOut[1] << G4endl;
		throw false;
	}
	SetFinal(0, mInitial);

	G4GenPhaseSpace gen;
	G4bool possible = gen.SetDecay(mInitial, 3, mOut);
	if(!possible) {
		G4cerr << "ERROR:: TntTwoNeutronDecayPhaseSpace:: Not enough energy for decay!" << G4endl;
		G4cerr << "MASSES (MBeam+ex, MF, 2*MN):: " << mInitial.m() << ", "
					 << mOut[0] << ", " << 2*mOut[1] << G4endl;
		throw possible;
	}

	G4double weight = do_weighted_generation(uniformWeight, gen);
	SetFinal(1, *(gen.GetDecay(0))); // fragment
	SetFinal(2, *(gen.GetDecay(1))); // neutron 1
	SetFinal(3, *(gen.GetDecay(2))); // neutron 1

	return weight;
}

















TntNeutronDecay* TntNeutronDecayFactory::Create()
{
	if(false) { }
	else if(GetDecayType() == "1n") { return new TntOneNeutronDecay(); }
	else if(GetDecayType() == "2nPhaseSpace") { return new TntTwoNeutronDecayPhaseSpace(); }
//	else if(GetDecayType() == "2nDiNeutron")  { return new TntTwoNeutronDecayDiNeutron();  }
	else {
		G4cerr << "TntNeutronDecayFactory::Create:: Invalid GetDecayType():: " << GetDecayType()
					 << G4endl;
		throw GetDecayType();
	}
}
