#include <fstream>
#include <cassert>
#include <map>
#include <CLHEP/Random/RandBreitWigner.h>
#include <G4SystemOfUnits.hh>
#include "G4GenPhaseSpace.hh"
#include "TntNuclearMasses.hh"
#include "TntNeutronDecay.hh"
#include "TntReactionGenerator.hh"
#include "TntRng.hh"
#include "TntError.hh"


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
	fVerb(2),
	mFinalFragMass(0),
	mInitial(),
	fReaction(0)
{ }

TntNeutronDecayIntermediate::~TntNeutronDecayIntermediate()
{ }

void TntNeutronDecayIntermediate::SetInputReaction(const TntReactionGenerator* r)
{
	fReaction = r;
	mInitial  = &(r->GetReactant(4));
	mFinalFragMass = 
		TntNuclearMasses::GetNuclearMass(mInitial->Z(), mInitial->A() - mNumberOfNeutronsEmitted)*MeV;
}

void TntNeutronDecayIntermediate::SetFinal(G4int indx, const G4LorentzVector& v)
{
	if(size_t(indx) < mFinal.size())	{
		mFinal[indx].set(v.x(), v.y(), v.z(), v.t());
	}
	else {
		if(GetVerboseLevel() > 0) {
			TNTERR << "Invalid index to TntNeutronDecayIntermediate::SetFinal:: " << indx << G4endl;
		}
		exit(1);
	}
}

const G4LorentzVector& TntNeutronDecayIntermediate::GetFinal(G4int indx) const
{
	if(size_t(indx) < mFinal.size())	{
		return mFinal[indx]; 
	}
	else {
		if(GetVerboseLevel() > 0) {
			TNTERR << "Invalid index to TntNeutronDecayIntermediate::GetFinal:: " << indx << G4endl;
		}
		exit(1);
	}
}

void TntNeutronDecayIntermediate::SetDecayParam(const G4String& par, G4double val)
{
	mParams[par] = val;
}

G4double TntNeutronDecayIntermediate::GetDecayParam(const G4String& par)
{
	std::map<G4String, G4double>::iterator it = mParams.find(par);
	if(it != mParams.end()) { return it->second; }

	if(GetVerboseLevel() > 0) {
		TNTERR << "TntNeutronDecayIntermediate::GetDecayParam:: Invalid Parameter "
					 << par << G4endl;
	}
	exit(1);
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
	TntNeutronEvaporation evap(mInitial->MplusEx(), mFinalFragMass, kNeutronMass);
	evap(&lvF, &lvN);

	G4ThreeVector t3boost = mInitial->Momentum().boostVector();
	lvN.boost(t3boost); // neutron 1
	lvF.boost(t3boost); // final fragment

	SetFinal(0, mInitial->Momentum());
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
	SetDecayParam("r0", 2.4);
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

	G4GenPhaseSpace gen;
	G4bool possible = gen.SetDecay(mInitial->Momentum(), 3, mOut);
	if(!possible) {
		if(GetVerboseLevel() > 1) {
			TNTERR << "TntTwoNeutronDecayPhaseSpace:: Not enough energy for decay!" << G4endl;
			TNTERR << "MASSES (MBeam+ex, MF, 2*MN):: " << mInitial->MplusEx() << ", "
						 << mOut[0] << ", " << 2*mOut[1] << G4endl;
		}
		return false;
	}

	if(fFSI) { // Use Final State Interaction
		/**
		* \note The 'r0' parameter is the "source size" from
		* PLB 476, 219 (2000). It is set in the 'reaction.dat' 
		* input file. The PLB paper gives respective source sizes for 
		* 6He, 11Li, 14Be of r0 = 2.4, 2.7, and 2.2 fm.
		*/
		G4double r0 = GetDecayParam("r0");
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

	SetFinal(0, mInitial->Momentum()); // recoil
	SetFinal(1, *(gen.GetDecay(0)));  // fragment
	SetFinal(2, *(gen.GetDecay(1)));  // neutron 1
	SetFinal(3, *(gen.GetDecay(2)));  // neutron 1

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
	mIntermediateFragMass =
		TntNuclearMasses::GetNuclearMass(mInitial->Z(), mInitial->A() - 1)*MeV;
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
	const G4double ex1 = GetDecayParam("ex1");  // Intermediate state EXCITATION energy
	const G4double  w1 = GetDecayParam("width1"); // Intermediate state width
	const G4double edecay = mInitial->MplusEx() - mFinalFragMass - 2*kNeutronMass; // TOTAL decay energy
	do {
		G4double exIntermediate = TntRngBreitWigner(ex1, w1).Generate();
		edecay1 = mInitial->MplusEx() - (mIntermediateFragMass+exIntermediate) - kNeutronMass;
		edecay2 = edecay - edecay1;
	} while(edecay1 < 0 || edecay2 < 0);


	////////////////////////////////////////
	// Evaportation of neutrons

	G4LorentzVector lvN1, lvF1, lvN2, lvF2;
	G4double intermediateMass = mInitial->MplusEx() - kNeutronMass - edecay1;

	// Do first neutron evaporation in COM frame
	{
		TntNeutronEvaporation evap1(mInitial->MplusEx(), intermediateMass, kNeutronMass);
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
	G4ThreeVector t3boost = mInitial->Momentum().boostVector();
	lvN1.boost(t3boost); // neutron 1
	lvN2.boost(t3boost); // neutron 2
	lvF2.boost(t3boost); // final fragment

	
	SetFinal(0, mInitial->Momentum()); // initial nucleus, before decay
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
{ }

TntTwoNeutronDecayDiNeutron::~TntTwoNeutronDecayDiNeutron()
{ }

/////////////////
// Code taken from st_reaction.cc out of st_mona simulation in
// use by the MoNA collaboration.
G4bool TntTwoNeutronDecayDiNeutron::Generate()
{
	////////////////////////////////////////////////
  // figure out decay energy, etc
	//
	const G4double s2n = // 2n separation energy
		mInitial->M() - mFinalFragMass - GetNumberOfNeutrons()*kNeutronMass;
	G4double ed2n; // total DECAY energy:  A -> [(A-2) + (2n)]
	G4double ei2n; // dineutron INTRINSIC energy
	try
	{
		/** \note Using 'Volya' dineutron decay.
		 *  Generates intrinsic and kinetic energies of dineutron.
		 *  The generation happens in TntRngTwoBodyGenerator, so we
		 *  need to read the last generated values from the
		 *  RNG used to pick the recoil excitation energy, when
		 *  the reaction generator was run.
		 */
		const TntRngVolyaDiNeutronEx& rngVolya =
			dynamic_cast<const TntRngVolyaDiNeutronEx&>(*fReaction->GetRngEx4());

		ei2n = rngVolya.GetRng2d()->GetLast().first;
		ed2n = s2n + rngVolya.GetRng2d()->GetLast().second;
	}
	catch (std::exception& e) 
	{
		// Wrong RNG class set in TntReactionGenerator
		if(GetVerboseLevel() > 0) {
			TNTERR << "TntTwoNeutronDecayDiNeutron:: TntReactionGenerator has the wrong RNG class"
						 << " for dineuton decay!" << G4endl;
		}
		exit(1);
	}
	
	if(ed2n < 0) {
		if(GetVerboseLevel() > 1) {
			TNTERR << "TntTwoNeutronDecayDiNeutron:: Not enough energy for decay!" << G4endl
						 << "MASSES (MBeam, ex, MF, 2*MN):: " << mInitial->M() << ", " 
						 << mInitial->Ex() << ", "
						 << mFinalFragMass << ", " << 2*kNeutronMass << G4endl
						 << "Edecay (2n), s2n, Eintrinsic (2n):" 
						 << ed2n << ", " << s2n << ", " << ei2n << G4endl;
		}
		return false;
	}

	////////////////////////////////////////////////
  // initial unbound state decay: A -> (A-2) + (2n)
	//
	G4LorentzVector lvFrag, lv2N; // fragment, dineutron
	{
		G4double eCM = 2*kNeutronMass + mFinalFragMass + ed2n; // total initial energy
		TntNeutronEvaporation evap(eCM, mFinalFragMass, 2*kNeutronMass);
		evap(&lvFrag, &lv2N);
	}

	////////////////////////////////////////////////
  // di-neutron decay: 2n -> n+n
	//
	G4LorentzVector lvN1, lvN2; // neutron 1, neutron 2
	{
		G4double eCM = 2*kNeutronMass + ei2n; // total CM energy
		TntNeutronEvaporation evap(eCM, kNeutronMass, kNeutronMass);
		evap(&lvN1, &lvN2);
	}
	
	////////////////////////////////////////////////
  // Boosting into the frame where the fragment is
	// at rest after the first decay
	//
  lvN1.boost(lv2N.boostVector());
  lvN2.boost(lv2N.boostVector());

	////////////////////////////////////////////////
  // final boost into lab frame + set return values
	//
  const G4ThreeVector& labBoost = mInitial->Momentum().boostVector();
  lvFrag.boost(labBoost);
  lvN1.boost(labBoost);
  lvN2.boost(labBoost);

	SetFinal(0, mInitial->Momentum()); // initial nucleus, before decay
	SetFinal(1, lvFrag);   // fragment
	SetFinal(2, lvN1);     // neutron 1
	SetFinal(3, lvN2);     // neutron 2

	return true;
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

	// Momentum along z-direction
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
		TNTERR << "TntNeutronDecayFactory::Create:: Invalid GetDecayType():: " << GetDecayType()
					 << G4endl;
		exit(1);
	}

	TntNeutronDecayIntermediate *di = dynamic_cast<TntNeutronDecayIntermediate*>(decay);
	if(di) {
		for(const auto& it : mOptions) {
			di->SetDecayParam(it.first, it.second);
		}
	}
	return decay;
}
