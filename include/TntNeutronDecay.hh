#ifndef TNT_NEUTRON_DECAY_HH
#define TNT_NEUTRON_DECAY_HH
#include <map>
#include <G4String.hh>
#include <G4ThreeVector.hh>
#include <G4LorentzVector.hh>

// Abstract base class for generic neutron decay
//
class TntNeutronDecay {
public:
	TntNeutronDecay() { }
	virtual ~TntNeutronDecay() { }
	virtual void SetBeam(G4int Z, G4int A, const G4ThreeVector& momentum) = 0;
	virtual void SetDecayParameter(const G4String& parname, G4double value) = 0;	
	virtual G4double GetDecayParameter(const G4String& parname) = 0;
	virtual G4double Generate(G4bool uniformWeight) = 0;
	virtual const G4LorentzVector& GetFinal(G4int indx) const = 0;
};


// Intermediate class, handling lots of the data stuff
// Still abstract, need to implement Generate()
//
class TntNeutronDecayIntermediate : public TntNeutronDecay {
public:
	TntNeutronDecayIntermediate(G4int number_of_neutrons_emitted);
	virtual ~TntNeutronDecayIntermediate();
	virtual void SetBeam(G4int Z, G4int A, const G4ThreeVector& momentum);
	virtual void SetDecayParameter(const G4String& parname, G4double value);
	virtual G4double GetDecayParameter(const G4String& parname);
	virtual G4double Generate(G4bool uniformWeight) = 0;
	virtual const G4LorentzVector& GetFinal(G4int indx) const;
protected:
	void SetFinal(G4int indx, const G4LorentzVector& v);
private:
	G4int mNumberOfNeutronsEmitted;
	std::map<G4String, G4double> mParams;
	std::vector<G4LorentzVector> mFinal;
protected:
	G4double mBeamMass; // Ground State Mass of initial nucleus
	G4double mFragMass; // Rest mass of decay fragment
	G4ThreeVector mInitialMomentum;
};

// Concrete class for single neutron decay, Breit Wigner
// Parameters to set are "energy" and "width"
// Setting "width" to zero returns a uniform (spike) decay energy
//
class TntOneNeutronDecay : public TntNeutronDecayIntermediate {
public:
	TntOneNeutronDecay();
	virtual ~TntOneNeutronDecay();
	virtual G4double Generate(G4bool uniformWeight);
};

// Concrete class for two neutron decay, Breit Wigner resonance energy,
// phase space decay.
// Parameters to set are "energy" and "width"
// Setting "width" to zero returns a uniform (spike) decay energy
//
class TntTwoNeutronDecayPhaseSpace : public TntNeutronDecayIntermediate {
public:
	TntTwoNeutronDecayPhaseSpace();
	virtual ~TntTwoNeutronDecayPhaseSpace();
	virtual G4double Generate(G4bool uniformWeight);
};


#endif
