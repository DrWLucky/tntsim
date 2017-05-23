// Porting of ROOT's TGenPhaseSpace code to the GEANT4 architecture
// Motivations are twofold: 1) Allow phase space calcs. w/o ROOT
//                          2) Use GEANT4 RNGs for phase space calcs, for consistency
// Note that Momentum, Energy units are Gev/C, GeV
//
#ifndef G4_GenPhaseSpace_GAC
#define G4_GenPhaseSpace_GAC

#include "G4LorentzVector.hh"

class G4GenPhaseSpace {
private:
	G4int        fNt;             // number of decay particles
	G4double     fMass[18];       // masses of particles
	G4double     fBeta[3];        // betas of decaying particle
	G4double     fTeCmTm;         // total energy in the C.M. minus the total mass
	G4double     fWtMax;          // maximum weigth
	G4LorentzVector  fDecPro[18];  //kinematics of the generated particles

	G4double PDK(G4double a, G4double b, G4double c);

public:
	G4GenPhaseSpace(): fNt(0), fMass(), fBeta(), fTeCmTm(0.), fWtMax(0.) {}
	G4GenPhaseSpace(const G4GenPhaseSpace &gen);
	virtual ~G4GenPhaseSpace() {}
	G4GenPhaseSpace& operator=(const G4GenPhaseSpace &gen);

	bool          SetDecay(G4LorentzVector &P, G4int nt, const G4double *mass, const char *opt="");
	G4double        Generate();
	G4LorentzVector *GetDecay(G4int n);

	G4int    GetNt()      const { return fNt;}
	G4double GetWtMax()   const { return fWtMax;}
};

#endif

