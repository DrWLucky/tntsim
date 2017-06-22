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

G4double TntRngUniform::DoGenerate()
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

G4double TntRngGaus::DoGenerate()
{
	return mSigma > 1e-9 ? CLHEP::RandGauss::shoot(mMean, mSigma) : mMean;
}


// TNT RNG BREIT WIGNER //

TntRngBreitWigner::TntRngBreitWigner(G4double mean, G4double sigma):
	mMean(mean), mSigma(sigma)
{ }

TntRngBreitWigner::~TntRngBreitWigner()
{ }

G4double TntRngBreitWigner::DoGenerate()
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

G4double TntRngCustom::DoGenerate()
{
	G4double r = TntRngUniform(0,1).Generate();
	std::vector<double>::const_iterator it =
		std::lower_bound(mCdf.begin(), mCdf.end(), r);
	size_t indx = it - mCdf.begin();

	if(indx+1 < mXlow.size()) {
		return TntRngUniform(mXlow[indx], mXlow[indx+1]).Generate();
	} else {
		G4cerr << "ERROR:: TntRngCustom::DoGenerate:: "
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

TntRngCustomAngDist::~TntRngCustomAngDist()
{ }



// TNT RNG VOLYA DINEUTRON EXCITATION ENERGY //

TntRngVolyaDiNeutronEx::TntRngVolyaDiNeutronEx(G4double e, G4double w, G4double a_s, G4double Ai):
	fR2d(new TntRngVolyaDiNeutron(e,w,a_s,Ai))
{ }

TntRngVolyaDiNeutronEx::~TntRngVolyaDiNeutronEx()
{ }

G4double TntRngVolyaDiNeutronEx::DoGenerate()
{
	auto gen2d = fR2d->Generate(); // dineutron ei, ek
	return gen2d.first + gen2d.second; // sum of energies is total decay energy
}



// TNT RNG GAUS 2D //

TntRngGaus2d::TntRngGaus2d(G4double sigma_x, G4double sigma_y, G4double rho):
	mSigmaX(sigma_x),
	mSigmaY(sigma_y),
	mRho(rho)
{  }

std::pair<G4double, G4double> TntRngGaus2d::DoGenerate()
{
	// Routine taken from GSL source code (GPL licensed)
	// Available online at:
	// https://github.com/LuaDist/gsl/blob/master/randist/bigauss.c
	//
	double u, v, r2, scale;
	TntRngUniform rng;
	
  do
    {
      /* choose x,y in uniform square (-1,-1) to (+1,+1) */

      u = -1 + 2 * rng.Generate();
      v = -1 + 2 * rng.Generate();

      /* see if it is in the unit circle */
      r2 = u * u + v * v;
    }
  while (r2 > 1.0 || r2 == 0);

  scale = sqrt (-2.0 * log (r2) / r2);

  double x = mSigmaX * u * scale;
  double y = mSigmaY * (mRho * u + sqrt(1 - mRho*mRho) * v) * scale;

	return std::make_pair(x, y);
}


// TNT RNG VOLYA DINEUTRON //

TntRngVolyaDiNeutron::TntRngVolyaDiNeutron(G4double E_i, G4double G_i, G4double a_s, G4int A_i):
	fEi(E_i), fGi(G_i), fAs(a_s), fA(A_i)
{
	// Set Paramaters
  G4double FragA = fA - 2; //Residual fragments after dineutron decay
	G4double Ei = fEi; //MeV Initial State
  G4double as = fAs; //n-n scattering length (fm)

  //NewParameter for addition of BW
  G4double Gamma_in = fGi;

  //Constants
  G4double r0 = 1.4; //fm
  //fixing the line below - JKS, 3 October 2013
  //G4double RR = r0 * pow(FragA + 2, 0.333333); //fm
  G4double RR = r0 * ( pow(FragA,0.333333) + pow(2,0.333333) ); //fm // channel radius for epsilon_K, first step
  // next line changed by JKS to account for channel radius in dineutron break-up
  G4double r = r0 * ( pow(1,0.333333) + pow(1,0.333333) ); //fm // channel radius for second step - breakup of dineutron

  const G4int nSteps = 6000;
  G4double min=0, max = 10.;
  G4double stepsize = (max-min) / nSteps;

	{
		char buf[4096*2], buf1[4096*2];
		sprintf(buf,"Inside StRNGVolya_Dineutron: Ei:%f GammaIn:%f as:%f FragA:%f\n",Ei,Gamma_in,as,FragA);
		sprintf(buf1,"  nSteps:%d  stepsize:%f  Range:0-%f MeV\n\n",nSteps,stepsize,stepsize*nSteps);
		G4cout << buf << G4endl;
		G4cout << buf1 << G4endl;
	}

  if((1.5*Ei) > max) { 
		G4cerr << "\n\nNeed to increase maxbin in StRNGVolya_2nseqGSL\n" << G4endl;
	}

  //Initialize the hists and pdfs.
  fH2d = gsl_histogram2d_alloc (nSteps,nSteps); // define histo gram and alloc
  gsl_histogram2d_set_ranges_uniform (fH2d, min, max, min, max); // set range and bins
  fP2d = gsl_histogram2d_pdf_alloc (nSteps,nSteps); // init PDF


  //Loop of BW
  for(int j=0; j<nSteps; j++){
    if((j%500)==0) {
			G4cout << "Creating 2-D PDF for Volya+BW distribution.  Step: " 
						 <<  j << " of " << nSteps <<  G4endl;
		}
    G4double Ebw0 = j*stepsize;
    G4double Ebw1 = (j+1)*stepsize;
    //Evaluate at midpoint
    G4double Ebw = (Ebw0+Ebw1)/2.;
    
    G4double DecayRateArray[nSteps]={0}; 
    G4double Total2nDecayWidth = 0;
    //Loop over Erel
    //JKS - does this have to be re-done just like the 2n sequential?
    for(int i=0; i<nSteps; i++){
      G4double di_ei0 = i*stepsize;  //intrinsic energy of the dineutron
      G4double di_ei1 = (i+1)*stepsize; 
      //Evaluate at midpoint
      G4double di_ei = (di_ei1+di_ei0)/2.;
      
      G4double DecayRate = -1;
      //why this if statement?
      if(di_ei<=Ebw) DecayRate = this->GetDineutronDecayRate(FragA,Ebw, di_ei, r, RR, as);
      else DecayRate = 0;
      DecayRate = DecayRate * stepsize;

      DecayRateArray[i] = DecayRate;
      Total2nDecayWidth = Total2nDecayWidth + DecayRate;
    }//i Erel   

    // This is different than the sequential decay. - JKS, 3 October 2013
    // Gamma_in is the input total decay width
    G4double Gamma_Total = Total2nDecayWidth * Gamma_in;
     //G4double Gamma_Total = Total2nDecayWidth + Gamma_in;

    G4double BW = Gamma_in / ( pow(Ebw - Ei, 2.0) + 0.25*pow(Gamma_Total, 2.0));
    BW = BW * stepsize;

    //loop again to calc prob for each erel of each Ebw
    for(int i=0; i<nSteps; i++){
      G4double di_ei0 = i*stepsize;  //intrinsic energy of the dineutron
      G4double di_ei1 = (i+1)*stepsize; 
      //Evaluate at midpoint
      G4double di_ei = (di_ei1+di_ei0)/2.;

      G4double Prob = BW * DecayRateArray[i];
      
      G4double di_ek = Ebw - di_ei; //kinetic energy of dineutron   

      //FILL HISTOGRAM FOR THE 2D PDF
      gsl_histogram2d_accumulate(fH2d, di_ei, di_ek, Prob);

    }//i second loop

  }//j Ebw
	
  
  //CREATE 2-D PDF
  gsl_histogram2d_pdf_init(fP2d, fH2d);
}


TntRngVolyaDiNeutron::~TntRngVolyaDiNeutron()
{
	gsl_histogram2d_free(fH2d);
	gsl_histogram2d_pdf_free(fP2d);
}

G4double TntRngVolyaDiNeutron::GetDineutronDecayRate(G4double FragA, 
																										 G4double Ebw, 
																										 G4double di_ei, 
																										 G4double r0,
																										 G4double RR, 
																										 G4double as) 
{
  G4double redMass = 931.494 * ((1.008*1.008) / (1.008 + 1.008)); //Reduced dineutron mass
  G4double Epsilon0 = (197*197) / (2*redMass*as*as); //MeV
  G4double term1 = sqrt( (Ebw - di_ei)*di_ei );
  G4double term2 = 1 + ( (r0*di_ei) / (2*as*Epsilon0) );
  term2 = term2*term2*Epsilon0 + di_ei;
  G4double DR = (1./CLHEP::pi) * (term1 / term2) * (r0 / RR);
  return DR;
}


std::pair<G4double, G4double> TntRngVolyaDiNeutron::DoGenerate()
{
  double RanNum =  TntRngUniform().Generate();
  double RanNum2 =  TntRngUniform().Generate();
   
  double EI, EK;
  gsl_histogram2d_pdf_sample(fP2d, RanNum, RanNum2, &EI, &EK);

	return std::make_pair(EI, EK);
}
