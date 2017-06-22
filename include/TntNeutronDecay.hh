#ifndef TNT_NEUTRON_DECAY_HH
#define TNT_NEUTRON_DECAY_HH
#include <map>
#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4LorentzVector.hh>


class TntReactionGenerator;

// Abstract base class for generic neutron decay
//
class TntNeutronDecay {
public:
	TntNeutronDecay() { }
	virtual ~TntNeutronDecay() { }
	virtual void SetInputReaction(const TntReactionGenerator* r) = 0;
	virtual G4int GetNumberOfNeutrons() const = 0;
	virtual G4bool Generate() = 0;
	virtual const G4LorentzVector& GetFinal(G4int indx) const = 0;
};

// Factory class
class TntNeutronDecayFactory {
public:
	void SetDecayType(G4String type) { fDecayType = type; }
	const G4String& GetDecayType() const { return fDecayType; }
	void SetDecayOption(G4String option, G4double value) { mOptions[option] = value; }
	G4double GetDecayOption(G4String option) const;
	TntNeutronDecay* Create();
private:
	G4String fDecayType;
	std::map<G4String, G4double> mOptions;
};


// Intermediate class, handling lots of the data stuff
// Still abstract, need to implement Generate()
//
class TntNeutronDecayIntermediate : public TntNeutronDecay {
public:
	TntNeutronDecayIntermediate(G4int number_of_neutrons_emitted);
	virtual ~TntNeutronDecayIntermediate();
	virtual void SetInputReaction(const TntReactionGenerator* r);
	virtual G4int GetNumberOfNeutrons() const { return mNumberOfNeutronsEmitted; }
	virtual G4bool Generate() = 0;
	void SetParam(const G4String& par, G4double val);
	G4double GetParam(const G4String& par);
protected:
	void SetFinal(G4int indx, const G4LorentzVector& v);
	virtual const G4LorentzVector& GetFinal(G4int indx) const;
private:
	G4int mNumberOfNeutronsEmitted;
	std::map<G4String, G4double> mParams;
	std::vector<G4LorentzVector> mFinal;
protected:
	G4double mInitialMass;   // Ground State Mass of initial nucleus
	G4double mFinalFragMass; // Rest mass of final decay fragment
	G4int mInitialA, mInitialZ;
	G4LorentzVector mInitial;
	const TntReactionGenerator* fReaction;
};

// Helper class to calculate neutron evaporation results
//
class TntNeutronEvaporation {
public:
	TntNeutronEvaporation(G4double m0, G4double mf, G4double mn):
		mM0(m0), mMf(mf), mMn(mn) { }
	~TntNeutronEvaporation() {}
	void operator() (G4LorentzVector* Frag, G4LorentzVector* Neut);
private:
	G4double mM0, mMf, mMn;
};


// Concrete class for single neutron decay, Breit Wigner
// Parameters to set are "energy" and "width"
// Setting "width" to zero returns a uniform (spike) decay energy
//
class TntOneNeutronDecay : public TntNeutronDecayIntermediate {
public:
	TntOneNeutronDecay();
	virtual ~TntOneNeutronDecay();
	virtual G4bool Generate();
};

// Concrete class for two neutron
// phase space decay. Includes optional
// final state interaction (FSI).
class TntTwoNeutronDecayPhaseSpace : public TntNeutronDecayIntermediate {
public:
	TntTwoNeutronDecayPhaseSpace(G4bool fsi = FALSE);
	virtual ~TntTwoNeutronDecayPhaseSpace();
	virtual G4bool Generate();
private:
	G4bool fFSI;
};

// Concrete class for two neutron
// 'dineutron' decay
//
class TntTwoNeutronDecayDiNeutron : public TntNeutronDecayIntermediate {
public:
	TntTwoNeutronDecayDiNeutron();
	virtual ~TntTwoNeutronDecayDiNeutron();
	virtual G4bool Generate();
private:
	//TntRngVolyaDiNeutronEx* fRngVolya;
};

// Concrete class for two neutron
// sequential decay
//
class TntTwoNeutronDecaySequential : public TntNeutronDecayIntermediate {
public:
	TntTwoNeutronDecaySequential();
	virtual ~TntTwoNeutronDecaySequential();
	virtual void SetInputReaction(const TntReactionGenerator* r);
	virtual G4bool Generate();
public:
	G4double mIntermediateFragMass;
};


#endif
