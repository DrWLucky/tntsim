// Porting of ROOT's TGenPhaseSpace code to the GEANT4 architecture
// Motivations are twofold: 1) Allow phase space calcs. w/o ROOT
//                          2) Use GEANT4 RNGs for phase space calcs, for consistency
// Note that Momentum, Energy units are Gev/C, GeV
//

#include <cstdlib>
#include "G4GenPhaseSpace.hh"
#include "Randomize.hh"

namespace {
const G4int kMAXP = 18;
const G4double kPI = 3.14159265358979;
}

////////////////////////////////////////////////////////////////////////////////
/// The PDK function.

G4double G4GenPhaseSpace::PDK(G4double a, G4double b, G4double c)
{
	G4double x = (a-b-c)*(a+b+c)*(a-b+c)*(a+b-c);
	x = sqrt(x)/(2*a);
	return x;
}

////////////////////////////////////////////////////////////////////////////////
/// Special max function

G4int DoubleMax(const void *a, const void *b)
{
	G4double aa = * ((G4double *) a);
	G4double bb = * ((G4double *) b);
	if (aa > bb) return  1;
	if (aa < bb) return -1;
	return 0;

}

////////////////////////////////////////////////////////////////////////////////
/// Copy constructor

G4GenPhaseSpace::G4GenPhaseSpace(const G4GenPhaseSpace &gen)
{
	fNt      = gen.fNt;
	fWtMax   = gen.fWtMax;
	fTeCmTm  = gen.fTeCmTm;
	fBeta[0] = gen.fBeta[0];
	fBeta[1] = gen.fBeta[1];
	fBeta[2] = gen.fBeta[2];
	for (G4int i=0;i<fNt;i++) {
		fMass[i]   = gen.fMass[i];
		fDecPro[i] = gen.fDecPro[i];
	}
}


////////////////////////////////////////////////////////////////////////////////
/// Assignment operator

G4GenPhaseSpace& G4GenPhaseSpace::operator=(const G4GenPhaseSpace &gen)
{
	fNt      = gen.fNt;
	fWtMax   = gen.fWtMax;
	fTeCmTm  = gen.fTeCmTm;
	fBeta[0] = gen.fBeta[0];
	fBeta[1] = gen.fBeta[1];
	fBeta[2] = gen.fBeta[2];
	for (G4int i=0;i<fNt;i++) {
		fMass[i]   = gen.fMass[i];
		fDecPro[i] = gen.fDecPro[i];
	}
	return *this;
}

////////////////////////////////////////////////////////////////////////////////
///  Generate a random final state.
///  The function returns the weight of the current event.
///  The G4LorentzVector of each decay product can be obtained using GetDecay(n).
///
/// Note that Momentum, Energy units are Gev/C, GeV

G4double G4GenPhaseSpace::Generate()
{
	G4double rno[kMAXP];
	rno[0] = 0;
	G4int n;
	if (fNt>2) {
		for (n=1; n<fNt-1; n++)  rno[n]=G4UniformRand();   // fNt-2 random numbers
		qsort(rno+1 ,fNt-2 ,sizeof(G4double) ,DoubleMax);  // sort them
	}
	rno[fNt-1] = 1;

	G4double invMas[kMAXP], sum=0;
	for (n=0; n<fNt; n++) {
		sum      += fMass[n];
		invMas[n] = rno[n]*fTeCmTm + sum;
	}

	//
	//-----> compute the weight of the current event
	//
	G4double wt=fWtMax;
	G4double pd[kMAXP];
	for (n=0; n<fNt-1; n++) {
		pd[n] = PDK(invMas[n+1],invMas[n],fMass[n+1]);
		wt *= pd[n];
	}

	//
	//-----> complete specification of event (Raubold-Lynch method)
	//
	fDecPro[0].set(0, pd[0], 0 , sqrt(pd[0]*pd[0]+fMass[0]*fMass[0]) );

	G4int i=1;
	G4int j;
	while (1) {
		fDecPro[i].set(0, -pd[i-1], 0 , sqrt(pd[i-1]*pd[i-1]+fMass[i]*fMass[i]) );

		G4double cZ   = 2*G4UniformRand() - 1;
		G4double sZ   = sqrt(1-cZ*cZ);
		G4double angY = 2*kPI * G4UniformRand();
		G4double cY   = cos(angY);
		G4double sY   = sin(angY);
		for (j=0; j<=i; j++) {
			G4LorentzVector *v = fDecPro+j;
			G4double x = v->px();
			G4double y = v->py();
			v->setPx( cZ*x - sZ*y );
			v->setPy( sZ*x + cZ*y );   // rotation around Z
			x = v->px();
			G4double z = v->pz();
			v->setPx( cY*x - sY*z );
			v->setPz( sY*x + cY*z );   // rotation around Y
		}

		if (i == (fNt-1)) break;

		G4double beta = pd[i] / sqrt(pd[i]*pd[i] + invMas[i]*invMas[i]);
		for (j=0; j<=i; j++) fDecPro[j].boost(0,beta,0);
		i++;
	}

	//
	//---> final boost of all particles
	//
	for (n=0;n<fNt;n++) fDecPro[n].boost(fBeta[0],fBeta[1],fBeta[2]);

	//
	//---> return the weight of event
	//
	return wt;
}

////////////////////////////////////////////////////////////////////////////////
/// Return Lorentz vector corresponding to decay n

G4LorentzVector *G4GenPhaseSpace::GetDecay(G4int n)
{
	if (n>fNt) return 0;
	return fDecPro+n;
}


////////////////////////////////////////////////////////////////////////////////
/// Input:
///  - G4LorentzVector &P:    decay particle (Momentum, Energy units are Gev/C, GeV)
///  - G4int nt:             number of decay products
///  - G4double *mass:       array of decay product masses
///  - const char *opt:        default -> constant cross section
///                       "Fermi" -> Fermi energy dependence
/// Return value:
///  - true:      the decay is permitted by kinematics
///  - false:     the decay is forbidden by kinematics
///

bool G4GenPhaseSpace::SetDecay(G4LorentzVector &P, G4int nt,
															 const G4double *mass, const char *opt)
{
	G4int n;
	fNt = nt;
	if (fNt<2 || fNt>18) return false;  // no more then 18 particle

	//
	//
	//
	fTeCmTm = P.mag();           // total energy in C.M. minus the sum of the masses
	for (n=0;n<fNt;n++) {
		fMass[n]  = mass[n];
		fTeCmTm  -= mass[n];
	}

	if (fTeCmTm<=0) return false;    // not enough energy for this decay

	//
	//------> the max weight depends on opt:
	//   opt == "Fermi"  --> fermi energy dependence for cross section
	//   else            --> constant cross section as function of TECM (default)
	//
	if (strcasecmp(opt,"fermi")==0) {
		// ffq[] = pi * (2*pi)**(FNt-2) / (FNt-2)!
		G4double ffq[] = {0
											,3.141592, 19.73921, 62.01255, 129.8788, 204.0131
											,256.3704, 268.4705, 240.9780, 189.2637
											,132.1308,  83.0202,  47.4210,  24.8295
											,12.0006,   5.3858,   2.2560,   0.8859 };
		fWtMax = pow(fTeCmTm,fNt-2) * ffq[fNt-1] / P.mag();

	} else {
		G4double emmax = fTeCmTm + fMass[0];
		G4double emmin = 0;
		G4double wtmax = 1;
		for (n=1; n<fNt; n++) {
			emmin += fMass[n-1];
			emmax += fMass[n];
			wtmax *= PDK(emmax, emmin, fMass[n]);
		}
		fWtMax = 1/wtmax;
	}

	//
	//---->  save the betas of the decaying particle
	//
	if (P.beta()) {
		G4double w = P.beta()/P.rho();
		fBeta[0] = P(0)*w;
		fBeta[1] = P(1)*w;
		fBeta[2] = P(2)*w;
	}
	else fBeta[0]=fBeta[1]=fBeta[2]=0;

	return true;
}
