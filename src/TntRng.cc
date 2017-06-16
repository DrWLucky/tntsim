#include <cassert>
#include <fstream>
#include <sstream>
#include <algorithm>
#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "TntRng.hh"

#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

namespace { 

int kMaxTries = 10000;


void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
}

void check_max_tries(G4int& n, const char* fct, G4double* low, G4double* high)
{
	if(n++ == kMaxTries) {
		G4cerr << "WARNING in TntRng::" << fct << ":: Failed to generate RNG ";
		if(low)  { G4cerr << "above " << *low; }
		if(high) { if(low) { G4cerr << " and"; } G4cerr << " below " << *high; }
		G4cerr << " after " << kMaxTries << " attempts. Enter 1 to abort program "
					 << "(and print stack trace), or 0 to continue..." << G4endl;
		G4int val;
		G4cin >> val;
		if(val != 0) { handler(0); exit(1); }
	}
}

}

// TNT RNG (BASE CLASS) //

G4double TntRng::GenerateAbove(G4double low)
{
	G4int n = 0;
	G4double value;
	do { 
		value = Generate();
		check_max_tries(n, "GenerateAbove", &low, 0);
	} while (value < low);
	return value;
}

G4double TntRng::GenerateBelow(G4double high)
{
	G4int n = 0;
	G4double value;
	do {
		value = Generate();
		check_max_tries(n, "GenerateBelow", 0, &high);
	} while (value >= high);
	return value;
}

G4double TntRng::GenerateBetween(G4double low, G4double high)
{
	G4int n = 0;
	G4double value;
	do { 
		value = Generate();
		check_max_tries(n, "GenerateBetween", &low, &high);
	} while (value < low || value >= high);
	return value;
}


// TNT RNG UNIFORM //

TntRngUniform::TntRngUniform(G4double low, G4double high):
	mLow(low), mHigh(high)
{ }

TntRngUniform::~TntRngUniform()
{ }

G4double TntRngUniform::Generate()
{
	G4double rndm = G4UniformRand();
	return mLow + (mHigh-mLow)*rndm;
}


// TNT RNG GAUS //

TntRngGaus::TntRngGaus(G4double mean, G4double sigma):
	mMean(mean), mSigma(sigma)
{ }

TntRngGaus::~TntRngGaus()
{ }

G4double TntRngGaus::Generate()
{
	return mSigma > 1e-9 ? CLHEP::RandGauss::shoot(mMean, mSigma) : mMean;
}


// TNT RNG BREIT WIGNER //

TntRngBreitWigner::TntRngBreitWigner(G4double mean, G4double sigma):
	mMean(mean), mSigma(sigma)
{ }

TntRngBreitWigner::~TntRngBreitWigner()
{ }

G4double TntRngBreitWigner::Generate()
{
	return mSigma > 1e-9 ? CLHEP::RandGauss::shoot(mMean, mSigma) : mMean;
}


namespace { void do_file_init(const G4String& filename,
															std::vector<G4double>& mXlow,
															std::vector<G4double>& mCdf,
															std::vector<G4double>* mPdf)
{
	if(filename == "NULL") return;
	
	std::ifstream ifs(filename.c_str());
  if(!ifs.good()) {
		G4cerr<< "ERROR:: Invalid filename passed to TntRngCustom ctor:: "
					<< filename << G4endl;
		throw ifs;
	}

  std::string line;
  std::vector<double> pdf;
  double integral = 0;

  double x, w;
  while(std::getline(ifs, line)) {
    std::stringstream iss(line);
    if (!(iss >> x >> w)) { continue; } // skip comments
		mXlow.push_back(x);
		pdf.push_back(w);
    integral += w;
  }

  double cdfi = 0;
  mCdf.resize(pdf.size());
  for(size_t i=0; i< pdf.size(); ++i) {
    pdf[i] /= integral;
    cdfi += pdf[i];
    mCdf[i] = cdfi;
  }	

	if(mPdf) {
		mPdf->resize(pdf.size());
		std::copy(pdf.begin(), pdf.end(), mPdf->begin());
	}
} }

// TNT RNG CUSTOM //

TntRngCustom::TntRngCustom(const G4String& filename)
{
	do_file_init(filename, mXlow, mCdf, 0);
}

TntRngCustom::~TntRngCustom()
{ }

G4double TntRngCustom::Generate()
{
	G4double r = TntRngUniform(0,1).Generate();
	std::vector<double>::const_iterator it =
		std::lower_bound(mCdf.begin(), mCdf.end(), r);
	size_t indx = it - mCdf.begin();

	if(indx+1 < mXlow.size()) {
		return TntRngUniform(mXlow[indx], mXlow[indx+1]).Generate();
	} else {
		G4cerr << "ERROR:: TntRngCustom::Generate:: "
					 << "Found out-of-range value in CDF. (indx, r): "
					 << indx << ", " << r << G4endl;
		throw indx;
		return 0;
	}
}



// TNT RNG CUSTOM ANG DIST //

TntRngCustomAngDist::TntRngCustomAngDist(const G4String& filename):
	TntRngCustom("NULL")
{
	std::vector<G4double> pdf;
	do_file_init(filename, mXlow, mCdf, &pdf);

	G4double integral = 0;
	for(size_t i=0; i< mXlow.size() - 1; ++i) {
		G4double diff = mXlow.at(i+1) - mXlow.at(i);
		G4double Angle = mXlow.at(i) + diff/2;
		pdf.at(i) = pdf[i]*sin(Angle*CLHEP::pi/180);
		integral += pdf.at(i);
	}

  double cdfi = 0;
  mCdf.resize(pdf.size());
  for(size_t i=0; i< pdf.size(); ++i) {
    pdf[i] /= integral;
    cdfi += pdf[i];
    mCdf[i] = cdfi;
  }
}

//      // multiply CStemp by sin(theta)
//       TF1* fsin = new TF1("sin",Form("1/(sin(x*%f/180.))",M_PI),0,180);
//       CStemp->Divide(fsin,1);
//       SetCrossSectionHist(CStemp);
//       delete fsin;
//     }

TntRngCustomAngDist::~TntRngCustomAngDist()
{ }
