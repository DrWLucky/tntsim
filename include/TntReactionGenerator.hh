#ifndef TNT_REACTION_HH
#define TNT_REACTION_HH
#include <initializer_list>
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "TntRng.hh"
#include "TntParticle.hh"
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
	virtual void SetRNGs( std::initializer_list<TntRng*> rngList ) = 0;
	virtual G4double GetThetaCM() const = 0;
	virtual G4double GetPhiCM() const = 0;
	virtual const TntParticle& GetReactant(G4int i) const = 0;
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
	/// rngEbeam, rngEx3, rngEx4, rngTheta, rngPhi
	void SetRNGs( std::initializer_list<TntRng*> rngList );
	G4double GetThetaCM() const { return fTheta; }
	G4double GetPhiCM() const { return fPhi; }
	const TntParticle& GetReactant(G4int i) const;
	const TntBeamEmittance* GetEmittanceX() const { return fEmX; }
	const TntBeamEmittance* GetEmittanceY() const { return fEmY; }
	const TntRng* GetRngEbeam() const { return fRngEbeam; }
	const TntRng* GetRngEx3()   const { return fRngEx3;   } /// EJECTILE
	const TntRng* GetRngEx4()   const { return fRngEx4;   } /// RECOIL
	const TntRng* GetRngTheta() const { return fRngTheta; }
	const TntRng* GetRngPhi()   const { return fRngPhi;   }

	G4bool Generate();
private:
	TntParticle fP1, fP2, fP3, fP4; /// Beam, target, ejectile, recoil
	/// Beam energy, ex particle 3 (ejectile), ex particle 4 (recoil), dSigma/dOmega
	TntRng *fRngEbeam, *fRngEx3, *fRngEx4, *fRngTheta, *fRngPhi;
	/// Beam emittance x, y
	TntBeamEmittance *fEmX, *fEmY;
	/// Theta, Phi from last generated Event
	G4double fTheta, fPhi;
};


/// N-body "neutron" phase space reaction generator
/** Generates final states using a n-body phase space
 *  calculation (e.g. TGenPhaseSpace from ROOT).
 *  The size of N is set based on the arguments to
 *  SetBeamTargetEjectile()
 */
class TntNeutronPhaseSpaceReactionGenerator : public TntReactionGenerator {
public:
	TntNeutronPhaseSpaceReactionGenerator(G4int n_neut);
	virtual ~TntNeutronPhaseSpaceReactionGenerator();
	virtual void SetBeamTargetEjectile(const G4String& beam, 
																		 const G4String& target, 
																		 const G4String& ejectile);
	virtual void SetBeamTargetEjectile(G4int Zbeam, G4int Abeam, 
																		 G4int Ztrgt, G4int Atrgt,
																		 G4int Zejectile, G4int Aejectile);
	virtual void SetEmittance(TntBeamEmittance* emX, TntBeamEmittance* emY);
	/// rngEbeam
	void SetRNGs( std::initializer_list<TntRng*> rngList );
	virtual G4double GetThetaCM() const { return fTheta; }
	virtual G4double GetPhiCM() const { return fPhi;   }
	virtual const TntParticle& GetReactant(G4int i) const;
	const TntBeamEmittance* GetEmittanceX() const { return fEmX; }
	const TntBeamEmittance* GetEmittanceY() const { return fEmY; }
	virtual const TntRng* GetRngEbeam() const  { return fRngEbeam; }
	virtual const TntRng* GetRngEx3()   const  { return 0; }
	virtual const TntRng* GetRngEx4()   const  { return 0; }
	virtual const TntRng* GetRngTheta() const  { return 0; }
	virtual const TntRng* GetRngPhi()   const  { return 0; }

	virtual G4bool Generate();
private:
  /// Beam, target, ejectile, recoil
	TntParticle fP1, fP2, fP3, fP4;
	/// Neutrons
	std::vector<TntParticle> fNeutrons;
	/// Beam energy, ex particle 3 (ejectile), ex particle 4 (recoil), dSigma/dOmega
	TntRng *fRngEbeam;
	/// Beam emittance x, y
	TntBeamEmittance *fEmX, *fEmY;
	/// Theta, Phi from last generated Event
	G4double fTheta, fPhi;
};


#endif
