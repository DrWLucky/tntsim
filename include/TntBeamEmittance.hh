#ifndef EMITTANCE_HEADER_INCLUDE_GUARD
#define EMITTANCE_HEADER_INCLUDE_GUARD
#include <vector>
#include <utility>

/// Class to calculate beam emittance using Courant-Snyder equation
/// 
class TntBeamEmittance {
public:
	TntBeamEmittance();
	TntBeamEmittance(double epsilon, double alpha,  double sigma_x);
	TntBeamEmittance(double epsilon, double alpha,  double sigma_x, double x0);
	void SetTwist(double epsilon, double alpha, double sigma_x);
	double GetEpsilon() const { return fEpsilon; }
	double GetAlpha()   const { return fAlpha;   }
	double GetBeta()    const { return fBeta;    }
	double GetGamma()   const { return fGamma;   }
	double GetSigmaX()  const;
	double GetSigmaTX() const;
	double GetRho()     const;
	void   SetX0(double x) { fX0 = x; }
	double GetX0()   const { return fX0; }
private:
	double fEpsilon, fBeta, fAlpha, fGamma, fX0;
};



#endif
