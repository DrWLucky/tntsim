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
	mBeamMass(0), 
	mFragMass(0),
	mInitialMomentum(0,0,0)
{ }

TntNeutronDecayIntermediate::~TntNeutronDecayIntermediate()
{ }

void TntNeutronDecayIntermediate::SetBeam(G4int Z, G4int A, const G4ThreeVector& momentum)
{
	mBeamMass = TntNuclearMasses::GetNuclearMass(Z, A)*MeV;
	mFragMass = TntNuclearMasses::GetNuclearMass(Z, A - mNumberOfNeutronsEmitted)*MeV;
	mInitialMomentum.set(momentum.x(), momentum.y(), momentum.z());
}

void TntNeutronDecayIntermediate::SetDecayParameter(const G4String& parname, G4double value)
{
	mParams[parname] = value;
}

G4double TntNeutronDecayIntermediate::GetDecayParameter(const G4String& parname)
{
	std::map<G4String, G4double>::iterator it = mParams.find(parname);
	if(it == mParams.end()) {
		G4cerr << "ERROR:: Invalid parameter name: " << parname << G4endl;
		throw parname;
	}
	return it->second;
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
	G4double e0 = GetDecayParameter("energy");
	G4double ex = e0;
	G4double mOut[2] = { mFragMass, kNeutronMass };

	if(ex+mBeamMass < mFragMass + kNeutronMass) {
		G4cerr << "ERROR:: TntTwoNeutronDecayPhaseSpace:: Not enough energy for decay!" << G4endl;
		G4cerr << "MASSES (MBeam, ex, MF, MN):: " << mBeamMass << ", " << ex  << ", "
					 << mOut[0] << ", " << mOut[1] << G4endl;
		throw false;
	}
if(GetDecayParameter("width") > 1e-6) {
		do {
			ex = CLHEP::RandBreitWigner::shoot(e0, GetDecayParameter("width"));
		} while(ex+mBeamMass < mFragMass+kNeutronMass);
	}
	
	G4double M0 = mBeamMass + ex;
	G4double Etot = sqrt(mInitialMomentum.mag2() + M0*M0);
	G4LorentzVector vBeam(mInitialMomentum, Etot);
	SetFinal(0, vBeam); // Beam

	G4GenPhaseSpace gen;
	G4bool possible = gen.SetDecay(vBeam, 2, mOut);
	if(!possible) {
		G4cerr << "ERROR:: TntOneNeutronDecay:: Not enough energy for decay!" << G4endl;
		G4cerr << "MASSES (MBeam, ex, MF, MN):: " << mBeamMass << ", " << ex  << ", "
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
	G4double e0 = GetDecayParameter("energy");
	G4double ex = e0;
	G4double mOut[3] = { mFragMass, kNeutronMass, kNeutronMass };
	if(ex+mBeamMass < mFragMass + 2*kNeutronMass) {
		G4cerr << "ERROR:: TntTwoNeutronDecayPhaseSpace:: Not enough energy for decay!" << G4endl;
		G4cerr << "MASSES (MBeam, ex, MF, 2*MN):: " << mBeamMass << ", " << ex  << ", "
					 << mOut[0] << ", " << 2*mOut[1] << G4endl;
		throw false;
	}
	if(GetDecayParameter("width") > 1e-6) {
		do {
			ex = CLHEP::RandBreitWigner::shoot(e0, GetDecayParameter("width"));
		} while(ex+mBeamMass < mFragMass + 2*kNeutronMass);
	}
	
	G4double M0 = mBeamMass + ex;
	G4double Etot = sqrt(mInitialMomentum.mag2() + M0*M0);
	G4LorentzVector vBeam(mInitialMomentum, Etot);
	SetFinal(0, vBeam); // Beam

	G4GenPhaseSpace gen;
	G4bool possible = gen.SetDecay(vBeam, 3, mOut);
	if(!possible) {
		G4cerr << "ERROR:: TntTwoNeutronDecayPhaseSpace:: Not enough energy for decay!" << G4endl;
		G4cerr << "MASSES (MBeam, ex, MF, 2*MN):: " << mBeamMass << ", " << ex  << ", "
					 << mOut[0] << ", " << 2*mOut[1] << G4endl;
		throw possible;
	}

	G4double weight = do_weighted_generation(uniformWeight, gen);
	SetFinal(1, *(gen.GetDecay(0))); // fragment
	SetFinal(2, *(gen.GetDecay(1))); // neutron 1
	SetFinal(3, *(gen.GetDecay(2))); // neutron 1
	return weight;
}
