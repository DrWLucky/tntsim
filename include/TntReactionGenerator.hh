#ifndef TNT_REACTION_HH
#define TNT_REACTION_HH
#include <memory>
#include <vector>
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "TntRng.hh"
#include "TntBeamEmittance.hh"



/// Class to GENERATE Two Body nuclear reactions
/// Includes beam emittance and excitation energy of recoil
class TntTwoBodyReactionGenerator {
public:
	TntTwoBodyReactionGenerator();
	virtual ~TntTwoBodyReactionGenerator();
	void SetBeamTargetEjectile(const G4String& beam, 
														 const G4String& target, 
														 const G4String& ejectile);
	void SetBeamTargetEjectile(G4int Zbeam, G4int Abeam, 
														 G4int Ztrgt, G4int Atrgt,
														 G4int Zejectile, G4int Aejectile);
	void SetEmittance(TntBeamEmittance* emX, TntBeamEmittance* emY);
	void SetRNGs(TntRng* rngEbeam, TntRng* rngEx3, TntRng* rngEx4, TntRng* rngTheta, TntRng* rngPhi);
	G4double GetReactantMass(G4int i) const;
	const G4LorentzVector& GetReactant(G4int i) const;
	const G4ThreeVector& GetPos() const { return fPos; }
	const TntBeamEmittance* GetEmittanceX() const { return fEmX.get(); }
	const TntBeamEmittance* GetEmittanceY() const { return fEmY.get(); }
	const TntRng* GetRngEbeam()   const { return fRngEbeam.get(); }
	const TntRng* GetRngEx3()     const { return fRngEx3.get(); } /// EJECTILE
	const TntRng* GetRngEx4()     const { return fRngEx4.get(); } /// RECOIL
	const TntRng* GetRngTheta()   const { return fRngTheta.get();}
	const TntRng* GetRngPhi()     const { return fRngPhi.get();  }

	G4bool Generate();
private:
	/// Beam, Target, Ejectile, Recoil
	G4double fM1,fM2,fM3,fM4;
	/// Beam, Target, Ejectile, Recoil
	G4LorentzVector fBeam, fTarget, fEjectile, fRecoil;
	/// Position
	G4ThreeVector fPos;
	/// beam energy, ex particle 3 (ejectile), ex particle 4 (recoil), dSigma/dOmega
	std::auto_ptr<TntRng> fRngEbeam, fRngEx3, fRngEx4, fRngTheta, fRngPhi;
	/// Beam emittance x, y
	std::auto_ptr<TntBeamEmittance> fEmX, fEmY;
};

#endif
