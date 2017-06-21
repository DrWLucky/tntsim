#include <cmath>
#include "TntBeamEmittance.hh"

TntBeamEmittance::TntBeamEmittance():
	fEpsilon(0), fBeta(0), fAlpha(0), fGamma(0), fX0(0)
{}

TntBeamEmittance::TntBeamEmittance(double epsilon, double alpha, double sigmaX)
{
	SetTwist(epsilon,alpha,sigmaX);
}

void TntBeamEmittance::SetTwist(double epsilon, double alpha, double sigmaX)
{
	fEpsilon = epsilon;
	fAlpha   = alpha;
	fBeta    = pow(sigmaX, 2) / fEpsilon;
	fGamma   = (1 + pow(fAlpha, 2)) / fBeta;
}

double TntBeamEmittance::GetSigmaX()const 
{
	return sqrt(fBeta*fEpsilon);  
}

double TntBeamEmittance::GetSigmaTX() const
{
	return sqrt(fGamma*fEpsilon);
}

double TntBeamEmittance::GetRho() const
{
	return -1*fAlpha*fEpsilon / (GetSigmaX()*GetSigmaTX()); 
}

