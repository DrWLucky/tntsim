#ifndef TNT_RNG_HH
#define TNT_RNG_HH
#include <vector>
#include <globals.hh>

class TntRng {
public:
	TntRng() { }
	virtual ~TntRng() { }
	virtual G4double Generate() = 0;
	G4double GenerateAbove(G4double low);
	G4double GenerateBelow(G4double high);
	G4double GenerateBetween(G4double low, G4double high);
};

class TntRngGaus : public TntRng {
public:
	TntRngGaus(G4double mean, G4double sigma);
	~TntRngGaus();
	G4double Generate();
private:
	G4double mMean, mSigma;
};

class TntRngBreitWigner : public TntRng {
public:
	TntRngBreitWigner(G4double mean, G4double sigma);
	~TntRngBreitWigner();
	G4double Generate();
private:
	G4double mMean, mSigma;
};

class TntRngUniform : public TntRng {
public:
	TntRngUniform(G4double low = 0, G4double high = 1);
	~TntRngUniform();
	G4double Generate();
private:
	G4double mLow, mHigh;
};

class TntRngCustom : public TntRng {
public:
	TntRngCustom(const G4String& filename);
	~TntRngCustom();
	G4double Generate();
protected:
	std::vector<G4double> mXlow, mCdf;
};

class TntRngCustomAngDist : public TntRngCustom {
public:
	TntRngCustomAngDist(const G4String& filename);
	~TntRngCustomAngDist();
};

#endif
