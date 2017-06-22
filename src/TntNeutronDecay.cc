#include <cassert>
#include <map>
#include <CLHEP/Random/RandBreitWigner.h>
#include <G4SystemOfUnits.hh>
#include "G4GenPhaseSpace.hh"
#include "TntNuclearMasses.hh"
#include "TntNeutronDecay.hh"
#include "TntReactionGenerator.hh"
#include "TntRng.hh"

namespace {

const G4double kNeutronMass =
	TntNuclearMasses::GetNuclearMass(0, 1)*MeV;

TntRngUniform kRngUniform(0,1); // RNG uniform 0->1
TntRngUniform kRngUniform_neg1_pos1(-1,1); // RNG uniform -1->1
TntRngUniform kRngUniform_0_2pi(0, 2*CLHEP::pi); // RNG uniform 0->2pi

template<typename T> T POW2(const T& t) { return t*t; }
} //namespace


/////////////////////////////////////////////////////////////////
// TntNeutronDecayIntermediate
//

TntNeutronDecayIntermediate::TntNeutronDecayIntermediate(G4int number_of_neutrons_emitted):
	mNumberOfNeutronsEmitted(number_of_neutrons_emitted),
	mFinal(number_of_neutrons_emitted + 2),
	mInitialMass(0),
	mFinalFragMass(0),
	mInitialA(0),
	mInitialZ(0), 
	mInitial(0,0,0,0),
	fReaction(0)
{ }

TntNeutronDecayIntermediate::~TntNeutronDecayIntermediate()
{ }

void TntNeutronDecayIntermediate::SetInputReaction(const TntReactionGenerator* r)
{
	assert(0 && "NEED TO IMPLEMENT TntNeutronDecayIntermediate::SetInputReaction!!!!");
	// mInitialA = A;
	// mInitialZ = Z;
	// mInitialMass = TntNuclearMasses::GetNuclearMass(Z, A)*MeV;
	// mFinalFragMass = TntNuclearMasses::GetNuclearMass(Z, A - mNumberOfNeutronsEmitted)*MeV;
	// mInitial = momentum;
	fReaction = r;
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

void TntNeutronDecayIntermediate::SetParam(const G4String& par, G4double val)
{
	mParams[par] = val;
}

G4double TntNeutronDecayIntermediate::GetParam(const G4String& par)
{
	std::map<G4String, G4double>::iterator it = mParams.find(par);
	if(it != mParams.end()) { return it->second; }
	G4cerr << "ERROR:: TntNeutronDecayIntermediate::GetParam:: Invalid Parameter "
				 << par << G4endl;
	throw par;
}



/////////////////////////////////////////////////////////////////
// TntOneNeutronDecay
//

TntOneNeutronDecay::TntOneNeutronDecay():
	TntNeutronDecayIntermediate(1)
{ }

TntOneNeutronDecay::~TntOneNeutronDecay()
{ }


G4bool TntOneNeutronDecay::Generate()
{
	G4LorentzVector lvF, lvN;
	TntNeutronEvaporation evap(mInitial.m(), mFinalFragMass, kNeutronMass);
	evap(&lvF, &lvN);

	G4ThreeVector t3boost = mInitial.boostVector();
	lvN.boost(t3boost); // neutron 1
	lvF.boost(t3boost); // final fragment

	SetFinal(0, mInitial);
	SetFinal(1, lvF);
	SetFinal(2, lvN);

	return true;
}


/////////////////////////////////////////////////////////////////
// TntTwoNeutronDecayPhaseSpace
//

TntTwoNeutronDecayPhaseSpace::TntTwoNeutronDecayPhaseSpace(G4bool fsi):
	TntNeutronDecayIntermediate(2),
	fFSI(fsi)
{
	SetParam("r0", 6.6);
}

TntTwoNeutronDecayPhaseSpace::~TntTwoNeutronDecayPhaseSpace()
{ }


namespace {
//Function from Marques for FSI three body
//email from MT to JKS on 3/18/2014 (later forwarded by
//JKS to GAC June 2017)
double Cnn0(double X, double r0)
{
	double hc,f0,d0;
	hc = 197.3;   // hc [MeVfm]
	f0 =  18.5;   // nn scatt. length [fm]
	d0 =   2.8;   // effective range  [fm]
	double    a,b,Ref,Imf,f_2,B0,Bi,pi,F1,F2,CNN0;

	if( r0 == 0. )
	{
		CNN0 = 1.;
		return CNN0;
	}

	pi = acos(-1.);
	a = 1./f0+d0*pow(X/hc,2)/2.; //
	b = -X/hc;                 // f = 1/(a+ib)
	Ref = a/(a*a+b*b);       //
	Imf = -b/(a*a+b*b);      //
	f_2 = Ref*Ref+Imf*Imf;
	B0 = -.5*exp(-4.*pow(r0*X/hc,2));    // t0 = 0
	if( X == 0. )
	{
		F1 = 1.;
		F2 = 0.;
	}
	else
	{
		F1 = (erf(2.*r0*X/hc)*sqrt(pi)/2.)*exp(-4.*pow(r0*X/hc,2))/(2.*r0*X/hc);
		F2 = (1.-exp(-4.*pow(r0*X/hc,2)))/(2.*r0*X/hc);
	}
	Bi = .25*f_2/pow(r0,2)*(1.-d0/(2.*sqrt(pi)*r0))+Ref*F1/(sqrt(pi)*r0)-Imf*F2/(2.*r0);

	CNN0 = 1.+B0+Bi;

	return CNN0;
} }

G4bool TntTwoNeutronDecayPhaseSpace::Generate()
{
	G4double mOut[3] = { mFinalFragMass, kNeutronMass, kNeutronMass };
	if(mInitial.m() < mFinalFragMass + 2*kNeutronMass) {
		G4cerr << "ERROR:: TntTwoNeutronDecayPhaseSpace:: Not enough energy for decay!" << G4endl;
		G4cerr << "MASSES (MBeam+ex, MF, 2*MN):: " << mInitial.m() << ", "
					 << mOut[0] << ", " << 2*mOut[1] << G4endl;
		return false;
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

	if(fFSI) { // Use Final State Interaction

		// Not 100% sure what this parameter is, need to ask JKS.
		// Likely it's the "source size" (avg. distance between the
		// two neutrons in the original nucleus). This was 6.6 fm in
		// the JKS NIM paper (and also the original paper on 11Li FSI
		// by Marquis).
		//
		// If this is correct, then the value needs to be dynamic...
		// and I should look up the correct value for 6He (if it exists)
		G4double r0 = GetParam("r0");
		G4double CnnM = Cnn0(0, r0);
		while(1) {
			G4double rel_weight = gen.Generate() / gen.GetWtMax();
			G4double ran = kRngUniform.Generate();

			G4LorentzVector ptemp1 = *gen.GetDecay(0);
			G4LorentzVector ptemp2 = *gen.GetDecay(1);
			G4LorentzVector ptemp3 = *gen.GetDecay(2);
			G4LorentzVector n2sys = ptemp2 + ptemp3;
			G4ThreeVector n2vel = n2sys.boostVector();
			ptemp2.boost(-n2vel);
			ptemp3.boost(-n2vel);

			// n.b. already in MeV here, but JKS code was ptemp2.P()*1000
			G4double Cnn = Cnn0(ptemp2.v().mag(), r0);
			G4double ran2 = TntRngUniform(0, CnnM).Generate();

			if( (ran < rel_weight) && (ran2 < Cnn) ) { break; }
		}
	}	else { // NO FSI
		G4double relwt, ran;
		do {
			relwt = gen.Generate() / gen.GetWtMax();
			ran = kRngUniform.Generate();
		} while(relwt < ran);
	}
	
	SetFinal(1, *(gen.GetDecay(0))); // fragment
	SetFinal(2, *(gen.GetDecay(1))); // neutron 1
	SetFinal(3, *(gen.GetDecay(2))); // neutron 1

	return true;
}



/////////////////////////////////////////////////////////////////
// TntTwoNeutronDecaySequential
//

TntTwoNeutronDecaySequential::TntTwoNeutronDecaySequential():
	TntNeutronDecayIntermediate(2),
	mIntermediateFragMass(0)
{ }

TntTwoNeutronDecaySequential::~TntTwoNeutronDecaySequential()
{ }

void TntTwoNeutronDecaySequential::SetInputReaction(const TntReactionGenerator* r)
{
	TntNeutronDecayIntermediate::SetInputReaction(r);
	mIntermediateFragMass = TntNuclearMasses::GetNuclearMass(mInitialZ, mInitialA - 1)*MeV;
}


/////////////////
// Code taken from st_reaction.cc out of st_mona simulation in
// use by the MoNA collaboration.
G4bool TntTwoNeutronDecaySequential::Generate()
{
	// Choose decay energies
	// NOTE: This is my own interpretation of double BW.
	// It is not necessarily correct. Once I have the RNGs from
	// Jenna, update to do a 'Volya' sequential decay.

	G4double edecay1, edecay2;
	const G4double ex1 = GetParam("ex1");  // Intermediate state EXCITATION energy
	const G4double  w1 = GetParam("width1"); // Intermediate state width
	const G4double edecay = mInitial.m() - mFinalFragMass - 2*kNeutronMass; // TOTAL decay energy
	do {
		G4double exIntermediate = TntRngBreitWigner(ex1, w1).Generate();
		edecay1 = mInitial.m() - (mIntermediateFragMass+exIntermediate) - kNeutronMass;
		edecay2 = edecay - edecay1;
	} while(edecay1 < 0 || edecay2 < 0);


	////////////////////////////////////////
	// Evaportation of neutrons

	G4LorentzVector lvN1, lvF1, lvN2, lvF2;
	G4double intermediateMass = mInitial.m() - kNeutronMass - edecay1;

	// Do first neutron evaporation in COM frame
	{
		TntNeutronEvaporation evap1(mInitial.m(), intermediateMass, kNeutronMass);
		evap1(&lvF1, &lvN1);
	}

	// Second neutron evaporation (COM frame)
	{
		TntNeutronEvaporation evap2(intermediateMass, mFinalFragMass, kNeutronMass);
		evap2(&lvF2, &lvN2);
	}

	
  ////////////////////////////////////////
	// Boost into lab frame

	// First boost to frame where frag is at rest after first decay
	// n.b. neutron 1 (lvN1) is already in this frame
	G4ThreeVector t4boost = lvF1.boostVector();
	lvF2.boost(t4boost); // FINAL fragment
	lvN2.boost(t4boost); // neutron from SECOND decay

	// Now boost all three into lab frame
	G4ThreeVector t3boost = mInitial.boostVector();
	lvN1.boost(t3boost); // neutron 1
	lvN2.boost(t3boost); // neutron 2
	lvF2.boost(t3boost); // final fragment

	
	SetFinal(0, mInitial); // initial nucleus, before decay
	SetFinal(1, lvF2);     // fragment
	SetFinal(2, lvN1);     // neutron 1
	SetFinal(3, lvN2);     // neutron 2

	return true;
}



/////////////////////////////////////////////////////////////////
// TntTwoNeutronDecayDiNeutron
//

TntTwoNeutronDecayDiNeutron::TntTwoNeutronDecayDiNeutron():
	TntNeutronDecayIntermediate(2)
{
	//fRngVolya = 
}

TntTwoNeutronDecayDiNeutron::~TntTwoNeutronDecayDiNeutron()
{
//	if(fRngVolya) { delete fRngVolya; fRngVolya = 0; }
}

/////////////////
// Code taken from st_reaction.cc out of st_mona simulation in
// use by the MoNA collaboration.
G4bool TntTwoNeutronDecayDiNeutron::Generate()
{
	return true;

#if 0
	////////////////////////////////////////////////
  // initial unbound state decay: A -> (A-2) + (2n)
	
	if(mInitial.m() < mFinalFragMass + 2*kNeutronMass) {
		G4cerr << "ERROR:: TntTwoNeutronDecayPhaseSpace:: Not enough energy for decay!" << G4endl;
		G4cerr << "MASSES (MBeam+ex, MF, 2*MN):: " << mInitial.m() << ", "
					 << mFinalFragMass << ", " << 2*kNeutronMass << G4endl;
		throw false;
	}

	// heavy fragment
	G4LorentzVector lvFrag(0,0,0,mFinalFragMass);

	// di-neutron
	G4LorentzVector lv2N(0,0,0,2*kNeutronMass);

	///TODO:: Need Generators for exenTotal, exenDiNeutron, incl. Volya stuff.
	/// (see st_reaction.cc, line 967)
	double exenTotal = mInitial.m() - lvFrag.m() - lv2N.m(); // TOTAL decay energy
	double exenDiNeutron = 118.5*keV;  // NOMINAL dineutron breakup energy, and POSSIBLY WRONG!!
	exenDiNeutron = TntRngBreitWigner(118.5*keV, 100*keV).GenerateAbove(0);
	double exen12_left = exenTotal - exenDiNeutron;

	if(!fRngVolya) {
//		fRngVolya = new TntRngVolyaDiNeutron(exenTotal, 
																				 // (G4double Ei, G4double Gi, G4double a_s, G4int Ai):
	}

	
	double e2N, eF;   // total neutron and fragment energy
  double eCM;   // total CM energy

	eCM = lv2N.m() + lvFrag.m() + exen12_left;  // total CM energy
	e2N = POW2(eCM) + POW2(lv2N.m()) - POW2(lvFrag.m());
	e2N = e2N/(2*eCM); // total energy neutron
  eF  = POW2(eCM) - POW2(lv2N.m()) + POW2(lvFrag.m());
  eF = eF/(2.*eCM);  // total energy of fragment
	
  double p2N,pF; // fragment and neutron momentum
  p2N = e2N*e2N - POW2(lv2N.m());
  p2N = sqrt(p2N);
  pF = eF*eF - POW2(lvFrag.m());
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
	eN_1 = eN_2 = POW2(eCM2);
	eN_1 = eN_1/(2.*eCM2); // total energy of neutron 1
	eN_2 = eN_2/(2.*eCM2); // total energy of neutron 2
	
  double pN_2, pN_1; // fragment and neutron momentum
  pN_2 = +sqrt(POW2(eN_2) - POW2(kNeutronMass));
  pN_1 = -sqrt(POW2(eN_1) - POW2(kNeutronMass)); // goes in opposite direction

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

	return true;
#endif
}




///////////////////////////////////////////////////////////////
// TntNeutronEvaporation
// 
void TntNeutronEvaporation::operator() (G4LorentzVector* Frag, G4LorentzVector* Neut)
{
	// Method for simulating neutron evaporation in center of mass frame
	// Copied from st_mona code. Results still need to be boosted into the
	// lab frame.

	// total CM energy
	G4double eCM = mM0;

	// total neutron, frag energies
	G4double eN = (POW2(eCM) + POW2(mMn) - POW2(mMf)) / (2*eCM);
	G4double eF = (POW2(eCM) - POW2(mMn) + POW2(mMf)) / (2*eCM);

	// neutron, frag momenta
	G4double pN = POW2(eN) - POW2(mMn);
	pN = pN < 0 ? 0 : sqrt(pN);
	
	G4double pF = POW2(eF) - POW2(mMf);
	pF = pF < 0 ? 0 : -1*sqrt(pF); // n.b.: opposite direction

	// Momentum laong z-direction
	Frag->set(0,0,pF,eF);
	Neut->set(0,0,pN,eN);

	// Generate & set angles
	G4double cosTheta = kRngUniform_neg1_pos1.Generate();
	G4double theta = acos(cosTheta);
	G4double phi = kRngUniform_0_2pi.Generate();

	Frag->setTheta(theta);
	Frag->setPhi(phi);
	Neut->setTheta(CLHEP::pi - theta);
	Neut->setPhi(phi + CLHEP::pi);
}


///////////////////////////////////////////////////////////////
// TntNeutronDecayFactory
//

G4double TntNeutronDecayFactory::GetDecayOption(G4String option) const
{
	std::map<G4String, G4double>::const_iterator it = mOptions.find(option);
	return it == mOptions.end() ? 0 : it->second;
}

namespace {
struct ToLower { void operator() (char& c) {
	const std::string up = "ABCDEFGHIJKLMNOPQRSTUVWXYZ";
	const std::string lo = "abcdefghijklmnopqrstuvwxyz";
	size_t i = up.find(c);
	if(i  < up.size() ) c = lo[i];
} };
	
void to_lower(std::string& str) {
	std::for_each(str.begin(), str.end(), ::ToLower());
} }

TntNeutronDecay* TntNeutronDecayFactory::Create()
{
	TntNeutronDecay* decay = 0;
	G4String decayType = this->GetDecayType();
	::to_lower(decayType);
	
	if(decayType == "1n") {
		decay = new TntOneNeutronDecay();      
	}
	else if(decayType == "2nphasespace") {
		decay = new TntTwoNeutronDecayPhaseSpace(FALSE); 
	}
	else if(decayType == "2nphasespacefsi") {
		decay = new TntTwoNeutronDecayPhaseSpace(TRUE); 
	}
	else if(decayType == "2ndineutron") {
		decay = new TntTwoNeutronDecayDiNeutron();
	}
	else if(decayType == "2nsequential") {
		decay = new TntTwoNeutronDecaySequential();
	}
	else {
		G4cerr << "TntNeutronDecayFactory::Create:: Invalid GetDecayType():: " << GetDecayType()
					 << G4endl;
		throw GetDecayType();
	}

	TntNeutronDecayIntermediate *di = dynamic_cast<TntNeutronDecayIntermediate*>(decay);
	if(di) {
		for(const auto& it : mOptions) {
			di->SetParam(it.first, it.second);
		}
	}
	return decay;
}
