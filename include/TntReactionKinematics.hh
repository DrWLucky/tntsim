#ifndef TNT_REACTION_KINEMATICS_HH
#define TNT_REACTION_KINEMATICS_HH
#include <memory>
#include <vector>
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "TntRng.hh"
#include "TntBeamEmittance.hh"



/// Class to calculate nuclear reaction kinematics
class TntReactionKinematics {
public:
	TntReactionKinematics();
	
	/// Set Initial Beam + Excited State Four-Vectors (LAB frame), and output masses [MeV]
	/// Outpus masses should include excitation energy
	TntReactionKinematics(const G4LorentzVector& beam, 
												const G4LorentzVector& target,
												const std::vector<G4double>& finalProductMasses);
	/// Empty
	virtual ~TntReactionKinematics();
	
	/// Set Initial Beam + Excited State Four-Vectors (LAB frame), and output masses [MeV]
	/// Outpus masses should include excitation energy
	void SetInputs(const G4LorentzVector& beam,
								 const G4LorentzVector& target, 
								 const std::vector<G4double>& finalProductMasses);

	/// Get Output Reaction product With index i, corresponding to the index in the input mass vector
	const G4LorentzVector& GetProduct(size_t i);
	/// Get Beam
	const G4LorentzVector& GetBeam() const { return fBeam; }
	/// Get Target
	const G4LorentzVector& GetTarget() const { return fTarget; }
	// Have Inputs?
	G4bool HaveInputs() const { return fHaveInputs; }
	/// Get Output Product Masses
	G4double GetProductMass(size_t i);
	/// Get Number of output products
	G4int GetNumProducts() const { return fOutputs.size(); }

	/// Calculate Output Lorentz Vectors, for a given theta, phi (both in center of mass)
	/// Pure virtual, must be implemented in child class
	virtual G4bool Calculate(G4double thetaCM, G4double phiCM) = 0;

protected:
	G4bool SetProduct(size_t i, const G4LorentzVector& v);
	
private:
	G4bool fHaveInputs;
	G4LorentzVector fBeam, fTarget;
	std::vector<G4double> fOutputMasses;
	std::vector<G4LorentzVector> fOutputs;
};


/// Class for Two Body Reaction Kinematics
class TntTwoBodyReactionKinematics : public TntReactionKinematics {
public:
	TntTwoBodyReactionKinematics() :
		TntReactionKinematics() { }
	TntTwoBodyReactionKinematics(const G4LorentzVector& beam,
															 const G4LorentzVector& target, 
															 const std::vector<G4double>& finalProductMasses) 
		:	TntReactionKinematics(beam, target, finalProductMasses) {  }
	virtual ~TntTwoBodyReactionKinematics() { }
	virtual G4bool Calculate(G4double thetaCM, G4double phiCM);
};



#endif
