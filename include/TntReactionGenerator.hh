#ifndef TNT_REACTION_HH
#define TNT_REACTION_HH
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "TntRng.hh"
#include "TntBeamEmittance.hh"



class TntReactionGenerator {
public:
	TntReactionGenerator() { }
	virtual ~TntReactionGenerator() { }
	virtual void SetBeamTargetEjectile(const G4String& beam, 
														 const G4String& target, 
														 const G4String& ejectile) = 0;
	virtual void SetBeamTargetEjectile(G4int Zbeam, G4int Abeam, 
														 G4int Ztrgt, G4int Atrgt,
														 G4int Zejectile, G4int Aejectile) = 0;
	virtual void SetEmittance(TntBeamEmittance* emX, TntBeamEmittance* emY) = 0;
	virtual void SetRNGs(TntRng* rngEbeam, TntRng* rngEx3, TntRng* rngEx4, TntRng* rngTheta, TntRng* rngPhi) = 0;
	virtual G4double GetReactantMass(G4int i) const = 0;
	virtual G4double GetThetaCM() const = 0;
	virtual G4double GetPhiCM() const = 0;
	virtual const G4LorentzVector& GetReactant(G4int i) const = 0;
	virtual const G4ThreeVector& GetPos() const = 0;
	virtual const TntBeamEmittance* GetEmittanceX() const = 0;
	virtual const TntBeamEmittance* GetEmittanceY() const = 0;
	virtual const TntRng* GetRngEbeam() const  = 0;
	virtual const TntRng* GetRngEx3()   const  = 0;
	virtual const TntRng* GetRngEx4()   const  = 0;
	virtual const TntRng* GetRngTheta() const  = 0;
	virtual const TntRng* GetRngPhi()   const  = 0;
	virtual G4bool Generate() = 0;
};


/// Class to GENERATE Two Body nuclear reactions
/// Includes beam emittance and excitation energy of recoil
class TntTwoBodyReactionGenerator : public TntReactionGenerator {
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
	G4double GetThetaCM() const { return fTheta; }
	G4double GetPhiCM() const { return fPhi; }
	const G4LorentzVector& GetReactant(G4int i) const;
	const G4ThreeVector& GetPos() const { return fPos; }
	const TntBeamEmittance* GetEmittanceX() const { return fEmX; }
	const TntBeamEmittance* GetEmittanceY() const { return fEmY; }
	const TntRng* GetRngEbeam() const { return fRngEbeam; }
	const TntRng* GetRngEx3()   const { return fRngEx3;   } /// EJECTILE
	const TntRng* GetRngEx4()   const { return fRngEx4;   } /// RECOIL
	const TntRng* GetRngTheta() const { return fRngTheta; }
	const TntRng* GetRngPhi()   const { return fRngPhi;   }

	G4bool Generate();
private:
	/// Beam, Target, Ejectile, Recoil
	G4double fM1,fM2,fM3,fM4;
	/// Beam, Target, Ejectile, Recoil
	G4LorentzVector fBeam, fTarget, fEjectile, fRecoil;
	/// Position
	G4ThreeVector fPos;
	/// beam energy, ex particle 3 (ejectile), ex particle 4 (recoil), dSigma/dOmega
	TntRng *fRngEbeam, *fRngEx3, *fRngEx4, *fRngTheta, *fRngPhi;
	/// Beam emittance x, y
	TntBeamEmittance *fEmX, *fEmY;
	/// Theta, Phi from last generated Event
	G4double fTheta, fPhi;
};

#endif
