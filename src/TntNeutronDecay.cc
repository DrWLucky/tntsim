#include <map>
#include <CLHEP/Random/RandBreitWigner.h>
#include <G4SystemOfUnits.hh>
#include "G4GenPhaseSpace.hh"
#include "TntNuclearMasses.hh"
#include "TntNeutronDecay.hh"
#include "TntRng.hh"

namespace {

const G4double kNeutronMass =
	TntNuclearMasses::GetNuclearMass(0, 1)*MeV;

TntRngUniform kRngUniform(0,1); // RNG uniform 0->1
TntRngUniform kRngUniform_neg1_pos1(-1,1); // RNG uniform -1->1
TntRngUniform kRngUniform_0_2pi(0, 2*CLHEP::pi); // RNG uniform 0->2pi
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
		while(weight < kRngUniform.Generate()) {
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





/////////////////////////////////////////////////////////////////
// TntTwoNeutronDecayDiNeutron
//

TntTwoNeutronDecayDiNeutron::TntTwoNeutronDecayDiNeutron():
	TntNeutronDecayIntermediate(2)
{ }

TntTwoNeutronDecayDiNeutron::~TntTwoNeutronDecayDiNeutron()
{ }


G4double TntTwoNeutronDecayDiNeutron::Generate(G4bool)
{
	////////////////////////////////////////////////
  // initial unbound state decay: A -> (A-2) + (2n)
	
	if(mInitial.m() < mFragMass + 2*kNeutronMass) {
		G4cerr << "ERROR:: TntTwoNeutronDecayPhaseSpace:: Not enough energy for decay!" << G4endl;
		G4cerr << "MASSES (MBeam+ex, MF, 2*MN):: " << mInitial.m() << ", "
					 << mFragMass << ", " << 2*kNeutronMass << G4endl;
		throw false;
	}

	// heavy fragment
	G4LorentzVector lvFrag(0,0,0,mFragMass);

	// di-neutron
	G4LorentzVector lv2N(0,0,0,2*kNeutronMass);

	///TODO:: Need Generators for exenTotal, exenDiNeutron, incl. Volya stuff.
	/// (see st_reaction.cc, line 967)
	double exenTotal = mInitial.m() - lvFrag.m() - lv2N.m(); // TOTAL decay energy
	double exenDiNeutron = 118.5*keV;  // NOMINAL dineutron breakup energy, and POSSIBLY WRONG!!
	exenDiNeutron = TntRngBreitWigner(118.5*keV, 100*keV).GenerateAbove(0);
	double exen12_left = exenTotal - exenDiNeutron;

	double e2N, eF;   // total neutron and fragment energy
  double eCM;   // total CM energy

	eCM = lv2N.m() + lvFrag.m() + exen12_left;  // total CM energy
	e2N = pow(eCM, 2) + pow(lv2N.m(), 2) - pow(lvFrag.m(), 2);
	e2N = e2N/(2*eCM); // total energy neutron
  eF  = pow(eCM, 2) - pow(lv2N.m(), 2) + pow(lvFrag.m(), 2);
  eF = eF/(2.*eCM);  // total energy of fragment
	
  double p2N,pF; // fragment and neutron momentum
  p2N = e2N*e2N - pow(lv2N.m(), 2);
  p2N = sqrt(p2N);
  pF = eF*eF - pow(lvFrag.m(), 2);
  if (pF < 0) { pF = 0; }
  else        { pF = sqrt(pF); }
  pF = -pF;  // fragment goes in opposite direction

  lvFrag.set(0,0,pF,eF);
  lv2N.set(0,0,p2N,e2N);

	// set theta and phi for the first decay to some random value
  double cosTheta = kRngUniform_neg1_pos1.Generate();  // cos(theta)
  double theta    = acos(cosTheta);
  double phi      = kRngUniform_0_2pi.Generate();

  lvFrag.setTheta(theta); // Set fragment angle
  lvFrag.setPhi(phi);     // Set fragment phi angle
  lv2N.setTheta(CLHEP::pi - theta); // Set di-neutron angle 180-theta
  lv2N.setPhi(phi + CLHEP::pi);     // set di-neutron phi angle

  G4ThreeVector t4Boost = lv2N.boostVector();  // boost vector for di-neutron



	////////////////////////////////////////////////
  // start of the di-neutron decay: 2n -> n+n

  G4LorentzVector lvN1(0,0,0,kNeutronMass);  // particle 3 (neutron 1)
  G4LorentzVector lvN2(0,0,0,kNeutronMass);  // particle 4 (neutron 2)

  double eN_1, eN_2; // total neutron and fragment energy
  double eCM2;       // total CM energy

  eCM2 = 2*kNeutronMass + exenDiNeutron;  // total CM energy
	eN_1 = eN_2 = pow(eCM2, 2);
	eN_1 = eN_1/(2.*eCM2); // total energy of neutron 1
	eN_2 = eN_2/(2.*eCM2); // total energy of neutron 2
	
  double pN_2, pN_1; // fragment and neutron momentum
  pN_2 = +sqrt(eN_2*eN_2 - pow(kNeutronMass, 2));
  pN_1 = -sqrt(eN_1*eN_1 - pow(kNeutronMass, 2)); // goes in opposite direction

  lvN1.set(0,0,pN_1,eN_1);
  lvN2.set(0,0,pN_2,eN_2);

  // set theta and phi to some random value
  cosTheta = kRngUniform_neg1_pos1.Generate(); // cos(theta)
  theta    = acos(cosTheta);
  phi      = kRngUniform_0_2pi.Generate();     // phi
  lvN1.setTheta(theta);       // Set neutron1 angle
  lvN1.setPhi(phi);           // Set neutron1 phi angle
  lvN2.setTheta(CLHEP::pi - theta); // Set neutron2 angle 180-theta
  lvN2.setPhi(phi + CLHEP::pi);     // Set neutron2 phi angle

  //Boosting into the frame where the fragment is at rest after the first decay
  lvN1.boost(t4Boost);
  lvN2.boost(t4Boost);

//		G4cerr << eCM2 << "\t" << lvN1.e()-kNeutronMass << "\t" << lvN1.m() << "<<<<\n";

	

	////////////////////////////////////////////////
  // final boost into lab frame + set return values
	
  G4ThreeVector t3Boost = mInitial.boostVector();
  lvFrag.boost(t3Boost);
  lvN1.boost(t3Boost);
  lvN2.boost(t3Boost);

	SetFinal(0, mInitial); // initial nucleus, before decay
	SetFinal(1, lvFrag);   // fragment
	SetFinal(2, lvN1);     // neutron 1
	SetFinal(3, lvN2);     // neutron 2

//	G4cerr << lvFrag.e() - lvFrag.m() << "\t" << lvN1.e() - lvN1.m() << "\t" << lvN2.e() - lvN2.m() << G4endl;
	
	return 1.; // Always constant weight value = 1
}















/////////////////////////////////////////////////////////////////
// TntNeutronDecayFactory
//

TntNeutronDecay* TntNeutronDecayFactory::Create()
{
	if(false) { }
	else if(GetDecayType() == "1n") { return new TntOneNeutronDecay(); }
	else if(GetDecayType() == "2nPhaseSpace") { return new TntTwoNeutronDecayPhaseSpace(); }
	else if(GetDecayType() == "2nDiNeutron")  { return new TntTwoNeutronDecayDiNeutron();  }
	else {
		G4cerr << "TntNeutronDecayFactory::Create:: Invalid GetDecayType():: " << GetDecayType()
					 << G4endl;
		throw GetDecayType();
	}
}
