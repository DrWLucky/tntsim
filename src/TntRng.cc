#include <stdio.h>
#include <execinfo.h>
#include <signal.h>
#include <stdlib.h>
#include <unistd.h>

#include <cassert>
#include <fstream>
#include <sstream>
#include <algorithm>

#include "G4PhysicalConstants.hh"
#include "Randomize.hh"
#include "TntError.hh"
#include "TntRng.hh"


namespace { 

void handler(int sig) {
  void *array[10];
  size_t size;

  // get void*'s for all entries on the stack
  size = backtrace(array, 10);

  // print out all the frames to stderr
  fprintf(stderr, "Error: signal %d:\n", sig);
  backtrace_symbols_fd(array, size, STDERR_FILENO);
  exit(1);
} }

void TntCheckMaxTries::operator() (G4int& n, const char* fct, G4double* low, G4double* high)
{
	if(n++ == kMaxTries) {
		TNTWAR << "TntRng::" << fct << ":: Failed to generate RNG fulfilling requirements";
		if(low)  { 
			G4cerr << "(above " << *low; 
		}
		if(high) { 
			if(low) { G4cerr << " and"; } else { G4cerr << "("; }
			G4cerr << " below " << *high << ")"; 
		}
		G4cerr << " after " << kMaxTries << " attempts. Enter 1 to abort program "
					 << "(and print stack trace), or 0 to continue..." << G4endl;
		G4int val;
		G4cin >> val;
		if(val != 0) { handler(0); exit(1); }
	}
}

// TNT RNG (BASE CLASS) //

G4double TntRng::GenerateAbove(G4double low)
{
	G4int n = 0;
	G4double value;
	do { 
		value = Generate();
		TntCheckMaxTries() (n, "GenerateAbove", &low, 0);
	} while (value < low);
	return value;
}

G4double TntRng::GenerateBelow(G4double high)
{
	G4int n = 0;
	G4double value;
	do {
		value = Generate();
		TntCheckMaxTries() (n, "GenerateBelow", 0, &high);
	} while (value >= high);
	return value;
}

G4double TntRng::GenerateBetween(G4double low, G4double high)
{
	G4int n = 0;
	G4double value;
	do { 
		value = Generate();
		TntCheckMaxTries() (n, "GenerateBetween", &low, &high);
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

TntRngVolyaDiNeutron::TntRngVolyaDiNeutron(G4double E_i, G4double G_i, G4double as, G4int A_i)
{
	// Set Paramaters
  G4double FragA = A_i - 2; //Residual fragments after dineutron decay
	G4double Ei = E_i;   //MeV Initial State

  //NewParameter for addition of BW
  G4double Gamma_in = G_i;

  //Constants
  const G4double r0 = 1.4; //fm
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
      if(di_ei<=Ebw) DecayRate = this->GetDineutronDecayRate(Ebw, di_ei, r, RR, as);
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

G4double TntRngVolyaDiNeutron::GetDineutronDecayRate(G4double Ebw, 
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
	if(fP2d == 0) {
		TNTERR << "TntRngVolyaDineutron::DoGenerate :: PDF histogram not set!" << G4endl;
		return std::make_pair(0,0);
	}
	
  double RanNum  = TntRngUniform().Generate();
  double RanNum2 = TntRngUniform().Generate();
   
  double EI, EK;
  gsl_histogram2d_pdf_sample(fP2d, RanNum, RanNum2, &EI, &EK);

	return std::make_pair(EI, EK);
}



///////////////////////////////////////////////////////
// TNT RNG VOLYA SEQUENTIAL
//

TntRngVolyaSequential::TntRngVolyaSequential(double Ei, double Ev, double sI,
																						 double sV, int L, double Gamma_in,
																						 int fragA)
{
#if 0	
  //public variable
  fInitialStateEnergy = par1;
#endif
  double FragA = fragA - 2; //remove 2 neutrons
	double S = Ev - Ei;  //MeV Threshold
	{
		char buf[4096*2];
		sprintf(buf,"Inside StRNGVolya_2nseqGSL: Ei:%f L:%d S:%f Ev:%f sI:%f sV:%f GammaIn:%f FragA:%f\n",Ei,L,S,Ev,sI,sV,Gamma_in,FragA);
		G4cout << buf << G4endl;
	}

 //-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*
  //   For L=0 seq. decays need to re-calcualte Ev such that virtual state 
  //   (based on s-wave scattering length is correcly calculated with no R dependence
  //   So, we need to calculated a new Ev that corresponds to the correct virtual state 
  //   of the input Ev above.  First lets set the current Ev to the desired Evirt;
  if(Ev<0 && L!=0) {
		TNTERR << "TntRngVolyaSequential:: EV CANNOT BE NEGATIVE WITH L!=0 !!!!!!!!" << G4endl;
		exit(1);
	}
  if(L==0){
    G4double Evirt = Ev;
    //Now calculate the scattering length corresponding to the virtual state
    G4double redMass = 931.494 * ((1.008 * (FragA)) / (1.008 + (FragA)));  //MeV/c2    
    G4double as = sqrt( (197.3269*197.3269) / (2*redMass * Evirt) );

    bool Input_as = false;
    if(Ev < 0){
      as = -Ev;
      Evirt = (197.3269*197.3269) / (2*redMass * as*as);
      Input_as = true;
    }
    
    //Next we can determine the energy of the intermiedate state (Ev) from the
    //scattering length based on the scattering length of a Breit-Wigner
    G4double R = 1.4 * (pow(1.008, 1./3.)+pow(FragA,1./3.)); //fm
    Ev = (sV*197.3269*197.3269) / (redMass*R*as);

    G4cout << "  L=0 requires re-calculation of intermiediate state" << G4endl;
    if(Input_as){
			char buf1[4096], buf2[4096];
      sprintf(buf1,"  Input Scattering Length of Intermediate State: -%f fm\n",as);
      sprintf(buf2,"  Corresponding Virtual State: %f MeV \n",Evirt);
			G4cout << buf1 << buf2 << G4endl;
    }else{
			char buf1[4096], buf2[4096];
      sprintf(buf1, "  Input Virtual State: %f MeV\n", Evirt);
      sprintf(buf2, "  Corresponding Scattering Length: -%f fm\n",as);
			G4cout << buf1 << buf2 << G4endl;
    }
		{
			char buf[4096];
			sprintf(buf, "  Ev used for calcualtion is: %f MeV\n",Ev);
			G4cout << buf << G4endl;
		}
  }

  //-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*-*


  const G4int nSteps = 6000;
  G4double min=0, max = 10.;
  G4double stepsize = (max-min) / nSteps;
	{
		char buf[4096];
		sprintf(buf, "nSteps:%d  stepsize:%f  Range:0-%f MeV\n\n",nSteps,stepsize,stepsize*nSteps);
		G4cout << buf << G4endl;
	}

  if((1.5*Ei) > max) {
		TNTERR << "TntRngVolyaSequential :: Need to increase maxbin in StRNGVolya_2nseqGSL\n" << G4endl;
		exit(1);
	}

  //Initialize the hists and pdfs.
//  fH2d = gsl_histogram2d_alloc (nSteps,nSteps); // define histo gram and alloc
//  second index (Erel) should be 2*nSteps
  fH2d = gsl_histogram2d_alloc (nSteps,2*nSteps); // define histo gram and alloc
//  gsl_histogram2d_set_ranges_uniform (fH2d, min, max, min, max); // set range and bins
//  second range (Erel) should be -max to max
  gsl_histogram2d_set_ranges_uniform (fH2d, min, max, -max, max); // set range and bins
//  fP2d = gsl_histogram2d_pdf_alloc (nSteps,nSteps); // init PDF
  fP2d = gsl_histogram2d_pdf_alloc (nSteps,2*nSteps); // init PDF


  //Loop of BW
  for(int j=0; j<nSteps; j++){
    if((j%500)==0) G4cout << "Creating 2-D PDF for Volya+BW distribution.  Step: " <<  j << " of " << nSteps <<  G4endl;
    G4double Ebw0 = j*stepsize;
    G4double Ebw1 = (j+1)*stepsize;
    //Evaluate at midpoint
    G4double Ebw = (Ebw0+Ebw1)/2.;
    
    G4double DecayRateArray[2*nSteps]={0}; 
    G4double Total2nDecayWidth = 0;
    //Loop over Erel
    for(int i=-nSteps; i<nSteps; i++){
      G4double Er0 = i*stepsize;  //Relative energy of e and e' (2 neutrons)
      G4double Er1 = (i+1)*stepsize;  //Relative energy of e and e' (2 neutrons)
      //Evaluate at midpoint
      G4double Er = (Er1+Er0)/2.;

      G4double e = 0.5 * (Ebw +Er);
      G4double ep = 0.5 * (Ebw -Er);

      //Reset Threshold to keep the intermediate state at a constant energy.
      S = Ev - Ebw;

      G4double DecayRate = -1;
      if(fabs(Er)<=Ebw) DecayRate = GetDecayRate(L,Ebw, S, Ev, Er, e, ep, sI,sV, FragA);
      else DecayRate = 0;
      DecayRate = DecayRate * stepsize;

      DecayRateArray[i+nSteps] = DecayRate;
      Total2nDecayWidth = Total2nDecayWidth + DecayRate;
    }//i Erel   

    G4double Gamma_Total = Total2nDecayWidth + Gamma_in;

    G4double BW = Gamma_in / ( pow(Ebw - Ei,2.0) + 0.25*pow(Gamma_Total,2.0));
    BW = BW * stepsize;

    //loop again to calc prob for each erel of each Ebw
    for(int i=-nSteps; i<nSteps; i++){
      G4double Er0 = i*stepsize;  //Relative energy of e and e' (2 neutrons)
      G4double Er1 = (i+1)*stepsize;  //Relative energy of e and e' (2 neutrons)
      //Evaluate at midpoint
      G4double Er = (Er1+Er0)/2.;
      
      G4double Prob = BW * DecayRateArray[i+nSteps];

      //FILL HISTOGRAM FOR THE 2D PDF
      gsl_histogram2d_accumulate(fH2d, Ebw, Er, Prob);

    }//i second loop

  }//j Ebw

  
  //CREATE 2-D PDF
  gsl_histogram2d_pdf_init(fP2d, fH2d);
}

TntRngVolyaSequential::~TntRngVolyaSequential()
{
	gsl_histogram2d_free(fH2d);
	gsl_histogram2d_pdf_free(fP2d);
}

std::pair<G4double, G4double> TntRngVolyaSequential::DoGenerate()
{
	if(fP2d == 0) {
		TNTERR << "TntRngVolyaSequential::DoGenerate :: PDF histogram not set!" << G4endl;
		return std::make_pair(0,0);
	}
	
  double RanNum  = TntRngUniform().Generate();
  double RanNum2 = TntRngUniform().Generate();
   
  double ebw, ereln;
  gsl_histogram2d_pdf_sample(fP2d, RanNum, RanNum2, &ebw, &ereln);

	return std::make_pair(ebw, ereln);
}

G4double TntRngVolyaSequential::SPDW(double ee, int L, double fragA, int add)
{
  //TODO: should have 2 reduced masses? for fragA and fragA+1
  double redMass = 931.494 * ((1.008 * (fragA+add)) / (1.008 + (fragA+add)));  //MeV/c2
  double hbarC = 197.3269; //MeV fm
  double R = 1.4 * (pow(1.008, 1./3.)+pow(fragA+add,1./3.)); //fm

  double k = sqrt(2*redMass*ee) / hbarC;

  double spdw = 2*hbarC*hbarC / (redMass*R*R);

  spdw = spdw * k * R * fabs( (2.*L-1.)/(2.*L+1.) );

  double transProb = 0;
  if(L==0) transProb = 1;
  else if(L==1) transProb = pow(k*R, 2.0) / (1 + pow(k*R, 2.0));
  else if(L==2) transProb = pow(k*R, 4.0) / (9 + 3*pow(k*R,2.0) + pow(k*R, 4.0));
  else {
		TNTERR << "TntRngVolyaSequential::SPDW :: Values of L=0,1,2 are only accepted at this time." << G4endl;
		exit(1);
	}
  
  spdw = spdw * transProb;
  return spdw;
}

G4double TntRngVolyaSequential::Gamma(double ee, int L, double fragA, int add, double S)
{
  double gamma = SPDW(ee,L,fragA,add) * S;
  //Add small const. to gamma for L>0 to account for a mathematical singularity that can occur,
  //which then leads to some unrealistic distributions.
  //if(L>0) gamma = gamma + 0.001;
  //Take it away because it leads to calculations that don't make sense
  //when energy of the intermediate state gets close to zero.
  //if(L>0) gamma = gamma;
  return gamma;
}

G4double TntRngVolyaSequential::GetDecayRate(int L, double Ei, double S, double /*Ev*/, double /*Er*/,
																						 double e, double ep, double sI, double sV, double fragA)
{
  double DR = 1./(8*CLHEP::pi);
  DR = DR * SPDW(e,L,fragA,1) * SPDW(ep,L,fragA,0);

  double A = Ei+2*S;
  double B = 0.5 * ( Gamma(e,L,fragA,1,sV) + Gamma(ep,L,fragA,0,sV) );

  //Numerator is now A - iB

  double C = S + e;
  double D = 0.5 * Gamma(ep,L,fragA,1,sV);

  double X = S + ep;
  double Y = 0.5 * Gamma(e,L,fragA,0,sV);

  //Denominator is now [C-iD] [X-iY]
  //Complex Modulus gives:
  //            sqrt(A^2 + B^2)
  //       -------------------------
  //       sqrt(C^2+D^2)*sqrt(X^2+Y^2)
  //
  //and then square everything and remove the sqrts

  double modsq = (A*A + B*B) / ( (C*C + D*D) * (X*X + Y*Y) );

  DR = DR * modsq;

  DR = DR * sI * sV;

  return DR;
}
