/// \brief Defines random number generator classes
#ifndef TNT_RNG_HH
#define TNT_RNG_HH
#include <vector>
#include <utility>
#include <gsl/gsl_histogram2d.h>
#include <globals.hh>



/// Set Global RNG seed
extern void TntSetRngSeed(G4int seed);

/// Get Global RNG seed
extern G4int TntGetRngSeed();

/// Abstract random number generator class
/** Derived classes implement details of generating random
 *  numbers from set distributions.
 *  \attention All derived classes should use the CLHEP random number
 *   generators ONLY. This ensures consistency throughout the simulation,
 *   in particular it means only one seed value needs to be changed to
 *   ensire an independent simulation.
 */
class TntRng {
public:
	/// Ctor
	TntRng() { }
	/// Dtor
	virtual ~TntRng() { }
	/// Generate a random number, and return the value
	G4double Generate() { return fR = DoGenerate(); }
	/// Generate a random number above some value, and return it
	G4double GenerateAbove(G4double low);
	/// Generate a random number below some value, and return it
	G4double GenerateBelow(G4double high);
	/// Generate a random number between two values, and return it
	G4double GenerateBetween(G4double low, G4double high);
	/// Get the most recently generated random number
	G4double GetLast() const { return fR; }
private:
	/// Implements the actual random number generation
	/** \returns The generated random number
	 */
	virtual G4double DoGenerate() = 0;
	/// Most recently generated random number
	G4double fR;
};

/// Abstract random number generator class, generating two
/// correlated random numbers from a distribution
/** Derived classes implement details of generating random
 *  numbers from set distributions.
 *  \attention All derived classes should use the CLHEP random number
 *   generators ONLY. This ensures consistency throughout the simulation,
 *   in particular it means only one seed value needs to be changed to
 *   ensire an independent simulation.
 */
class TntRng2d {
public:
	/// Ctor
	TntRng2d() { }
	/// Dtor
	virtual ~TntRng2d() { }
	/// Generate two correlated random numbers and return then
	/// as a `std::pair`
	std::pair<G4double, G4double> Generate() { return fR2 = DoGenerate(); }
	/// Get the most recently generated pair of random numbers
	std::pair<G4double, G4double> GetLast() const { return fR2; }
private:
	/// Implements the actual random number generation
	/** \returns The generated random number pair
	 */
	virtual std::pair<G4double, G4double> DoGenerate() = 0;
	/// Most recently generated pair
	std::pair<G4double, G4double> fR2;
};


///  Always returns a constant value
/** Not really a random number generator, but made to look like one
 */
class TntRngConstant : public TntRng {
public:
	///  Ctor
	/** \param [in] val The constant number to return */
	TntRngConstant(G4double val): fVal(val) { }
private:
	G4double DoGenerate() { return fVal; }
	G4double fVal;
};	

/// Gaussian [normal] distribution
class TntRngGaus : public TntRng {
public:
	/// Ctor
	/** \param [in] mean Mean value of gaussian
	 *  \param [in] sigma Standard deviation of gaussian
	 */
	TntRngGaus(G4double mean, G4double sigma);
	~TntRngGaus();
private:
	G4double DoGenerate();
private:
	G4double mMean, mSigma;
};

/// Breit Wigner distribution, as in resonance decay
class TntRngBreitWigner : public TntRng {
public:
	/// Ctor
	/** \param [in] mean Mean value
	 *  \param [in] width Width
	 */
	TntRngBreitWigner(G4double mean, G4double width);
	~TntRngBreitWigner();
private:
	G4double DoGenerate();
private:
	G4double mMean, mSigma;
};

/// Uniform distribution
class TntRngUniform : public TntRng {
public:
	/// Ctor
	/** \param [in] low Low value of uniform range (inclusive)
	 *  \param [in] high High value of uniform range (exclusive)
	 */
	TntRngUniform(G4double low = 0, G4double high = 1);
	~TntRngUniform();
private:
	G4double DoGenerate();
private:
	G4double mLow, mHigh;
};

/// Uses custom distrubition from a text file
class TntRngCustom : public TntRng {
public:
	/// Ctor
	/** \param [in] filename Name of the file containing the custom
	 *  distribution.
	 *
	 *  Should be structured like a histogram with the format:
	 *  \code
	 *  bin_low <whitespace> bin_value
	 *  \endcode
	 */
	TntRngCustom(const G4String& filename);
	~TntRngCustom();
private:
	G4double DoGenerate();
protected:
	std::vector<G4double> mXlow, mCdf;
};

/// Custom angular distribution for nuclear reactions
class TntRngCustomAngDist : public TntRngCustom {
public:
	/// Ctor
	/** \param [in] filename Name of the file containing the custom
	 *  distribution.
	 *  
	 *  Should be structured like a histogram with the format:
	 *  \code
	 *  bin_low <whitespace> bin_value
	 *  \endcode
	 *  The bin values should be dSigma/dOmega. Multiplication by
	 *  sin(theta) is taken care of automatically.
	 */
	TntRngCustomAngDist(const G4String& filename);
	~TntRngCustomAngDist();
};

/// Two-dimensional Gaussian [normal] with correlated x,y
/** No option to shift the means, so just add an offset value
 *  manually.
 */
class TntRngGaus2d : public TntRng2d {
public:
	/// Ctor
	/** \param [in] sigma_x standard deviation in x
	 *  \param [in] sigma_y standard deviation in y
	 *  \param [in] rho  Correlation coefficient
	 */
	TntRngGaus2d(G4double sigma_x, G4double sigma_y, G4double rho);
private:
	virtual std::pair<G4double, G4double> DoGenerate();
private:
	G4double mSigmaX, mSigmaY, mRho;
};

/// Custom generator for 'Volya' dineutron decay model
/** Generates pairs of <dineutron intrinsic energy, dineutron kinetic energy>
 *  The total decay energy is the sum of the two.
 *  \attention The constructor is highly non-trivial and very computationally
 *  expensive. It is recommended to create instances as few times as possible!
 */
class TntRngVolyaDiNeutron : public TntRng2d {
public:
	/// Ctor
	/** \param [in] Ei Central decay energy
	 *  \param [in] Gi Decay width
	 *  \param [in] a_s n-n scattering length (-18.7 fm is canonical)
	 *  \param [in] Ai Mass number of initial nucleus
	 */
	TntRngVolyaDiNeutron(G4double Ei, G4double Gi, G4double a_s, G4int Ai);
	virtual ~TntRngVolyaDiNeutron();
	/// Get internal histogram of CDF for picking random numbers
	const gsl_histogram2d_pdf* GetPdfHist() const { return fP2d; }
	/// Get internal histogram of PDF for picking random numbers
	const gsl_histogram2d* GetHist() const { return fH2d; }
private:
	virtual std::pair<G4double, G4double> DoGenerate();
	G4double GetDineutronDecayRate(G4double Ebw, G4double di_ei, G4double r0, 
																 G4double RR, G4double as);
	gsl_histogram2d_pdf* fP2d;
	gsl_histogram2d* fH2d;	
};

/// Custom excitation energy generator for 'Volya' dineutron decay model
/** Generates the total decay energy.
 *  Uses TntRngVolyaDiNeutron internally; pair of <EI, EK>
 *  can be accessed through internal 2d generator, using GetRng2d()
 *  \attention The constructor is highly non-trivial and very computationally
 *  expensive. It is recommended to create instances as few times as possible!
 */
class TntRngVolyaDiNeutronEx : public TntRng {
public:
	/// 
	TntRngVolyaDiNeutronEx(G4double e, G4double w, G4double a_s, G4double Ai);
	virtual ~TntRngVolyaDiNeutronEx();
	/// Returns internal 2d generator
	const TntRngVolyaDiNeutron* GetRng2d() const { return fR2d.get(); }
private:	
	G4double DoGenerate();
private:
	std::auto_ptr<TntRngVolyaDiNeutron> fR2d;
};

/// Custom generator for 'Volya' sequentual 2n decay model
/** Generates pairs of <dineutron intrinsic energy, dineutron kinetic energy>
 *  The total decay energy is the sum of the two.
 *  \attention The constructor is highly non-trivial and very computationally
 *  expensive. It is recommended to create instances as few times as possible!
 */
class TntRngVolyaSequential : public TntRng2d {
public:
	/// Ctor
	/** \param [in] Ei    Initial state energy
	 *  \param [in] Ev    Intermediate state energy, if Ev < 0 then it is assumed to be scattering length for L=0 decay
	 *  \param [in] sI    Spectroscopic factor for the I->V single particle decay
	 *  \param [in] sV    Spectroscopic factor for the V->F single particle decay
	 *  \param [in] L     Orbital angular momentum
	 *  \param [in] Gamma_in  Initial width
	 *  \param [in] fragA Initial fragment mass number
	 */
	TntRngVolyaSequential(double Ei, double Ev, double sI, double sV,
												int L, double Gamma_in, int fragA);
	virtual ~TntRngVolyaSequential();
	/// Get internal histogram of CDF for picking random numbers
	const gsl_histogram2d_pdf* GetPdfHist() const { return fP2d; }
	/// Get internal histogram of PDF for picking random numbers
	const gsl_histogram2d* GetHist() const { return fH2d; }
private:
	/// Returns pair: <EBW , EREL_N>
	virtual std::pair<G4double, G4double> DoGenerate();
	/// Get single particle decay width
	G4double SPDW(double ee, int L, double fragA, int add);
	/// Get gamma value
	G4double Gamma(double ee, int L, double fragA, int add, double S);
	/// Get decay rate
	G4double GetDecayRate(int L, double Ei, double S, double Ev, double Er, double e, double ep, double sI, double sV, double fragA);

private:
	gsl_histogram2d_pdf* fP2d;
	gsl_histogram2d* fH2d;	
};

/// Custom excitation energy generator for 'Volya' dineutron decay model
/** Generates the total decay energy.
 *  Uses TntRngVolyaSequential internally; pair of <EBW , EREL_N>
 *  can be accessed through internal 2d generator, using GetRng2d()
 *  \attention The constructor is highly non-trivial and very computationally
 *  expensive. It is recommended to create instances as few times as possible!
 */
class TntRngVolyaSequentialEx : public TntRng {
public:
	/// Ctor, same parameters as TntRngVolyaSequential
	TntRngVolyaSequentialEx(double Ei, double Ev, double sI, double sV,
													int L, double Gamma_in, int fragA);
	virtual ~TntRngVolyaSequentialEx();
	/// Returns internal 2d generator
	const TntRngVolyaSequential* GetRng2d() const { return fR2d.get(); }
private:	
	G4double DoGenerate();
private:
	std::auto_ptr<TntRngVolyaSequential> fR2d;
};




/// Class to prevent infinite loops when trying to
/// generate RNGs within a range
class TntCheckMaxTries {
public:
	/// Ctor
	/** \param [in] max Max number of tries before complaining.
	 */
	TntCheckMaxTries(G4int max=10000): kMaxTries(max) { }
	/// Check if a loop has reached the maximum number of tries.
	/** If it has, give user the option to exit the program, or to continue
	 */
	void operator() (G4int& n, const char* fct, G4double* low = 0, G4double* high = 0);
private:
	const G4int kMaxTries;
};



#endif
