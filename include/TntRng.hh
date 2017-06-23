#ifndef TNT_RNG_HH
#define TNT_RNG_HH
#include <vector>
#include <utility>
#include <gsl/gsl_histogram2d.h>
#include <globals.hh>


class TntRng {
public:
	TntRng() { }
	virtual ~TntRng() { }
	G4double Generate() { return fR = DoGenerate(); }
	G4double GenerateAbove(G4double low);
	G4double GenerateBelow(G4double high);
	G4double GenerateBetween(G4double low, G4double high);
	G4double GetLast() const { return fR; }
private:
	virtual G4double DoGenerate() = 0;
	G4double fR;
};

class TntRng2d {
public:
	TntRng2d() { }
	virtual ~TntRng2d() { }
	std::pair<G4double, G4double> Generate() { return fR2 = DoGenerate(); }
	std::pair<G4double, G4double> GetLast() const { return fR2; }
private:
	virtual std::pair<G4double, G4double> DoGenerate() = 0;
	std::pair<G4double, G4double> fR2;
};

class TntRngConstant : public TntRng {
public:
	TntRngConstant(G4double val): fVal(val) { }
private:
	G4double DoGenerate() { return fVal; }
	G4double fVal;
};	

class TntRngGaus : public TntRng {
public:
	TntRngGaus(G4double mean, G4double sigma);
	~TntRngGaus();
private:
	G4double DoGenerate();
private:
	G4double mMean, mSigma;
};

class TntRngBreitWigner : public TntRng {
public:
	TntRngBreitWigner(G4double mean, G4double sigma);
	~TntRngBreitWigner();
private:
	G4double DoGenerate();
private:
	G4double mMean, mSigma;
};

class TntRngUniform : public TntRng {
public:
	TntRngUniform(G4double low = 0, G4double high = 1);
	~TntRngUniform();
private:
	G4double DoGenerate();
private:
	G4double mLow, mHigh;
};

class TntRngCustom : public TntRng {
public:
	TntRngCustom(const G4String& filename);
	~TntRngCustom();
private:
	G4double DoGenerate();
protected:
	std::vector<G4double> mXlow, mCdf;
};

class TntRngCustomAngDist : public TntRngCustom {
public:
	TntRngCustomAngDist(const G4String& filename);
	~TntRngCustomAngDist();
};

class TntRngGaus2d : public TntRng2d {
public:
	TntRngGaus2d(G4double sigma_x, G4double sigma_y, G4double rho);
private:
	virtual std::pair<G4double, G4double> DoGenerate();
private:
	G4double mSigmaX, mSigmaY, mRho;
};

class TntRngVolyaDiNeutron : public TntRng2d {
public:
	TntRngVolyaDiNeutron(G4double Ei, G4double Gi, G4double a_s, G4int Ai);
	virtual ~TntRngVolyaDiNeutron();
	const gsl_histogram2d_pdf* GetPdfHist() const { return fP2d; }
	const gsl_histogram2d* GetHist() const { return fH2d; }
private:
	virtual std::pair<G4double, G4double> DoGenerate();
	G4double GetDineutronDecayRate(G4double Ebw, G4double di_ei, G4double r0, 
																 G4double RR, G4double as);
private:
	G4double fEi; // initial state energy
	G4double fGi; // initial state width
	G4double fAs; // n-n scattering length
	G4int fA;     // initial state mass number

	gsl_histogram2d_pdf* fP2d;
	gsl_histogram2d* fH2d;	
};

class TntRngVolyaDiNeutronEx : public TntRng {
public:
	TntRngVolyaDiNeutronEx(G4double e, G4double w, G4double a_s, G4double Ai);
	virtual ~TntRngVolyaDiNeutronEx();
	const TntRngVolyaDiNeutron* GetRng2d() const { return fR2d.get(); }
private:	
	G4double DoGenerate();
private:
	std::auto_ptr<TntRngVolyaDiNeutron> fR2d;
};


/// Class to prevent infinite loops when trying to
/// generate RNGs within a range
class TntCheckMaxTries {
public:
	TntCheckMaxTries(G4int max=10000): kMaxTries(max) { }
	void operator() (G4int& n, const char* fct, G4double* low = 0, G4double* high = 0);
private:
	const G4int kMaxTries;
};


#endif
