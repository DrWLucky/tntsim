#ifndef TNT_PARTICLE_HEADER
#define TNT_PARTICLE_HEADER
#include <G4PhysicalConstants.hh>
#include <G4ThreeVector.hh>
#include <G4LorentzVector.hh>


class TntParticle {
public:
	/// Ctor
	TntParticle();
	/// Dtor
	~TntParticle();
	/// Sets A,Z, and GROUND STATE Mass
	/** Doesn't touch kinematic properties */
	void SetNucleus(G4int Z, G4int A);  /// Sets A,Z, and Mass
	/// Sets A,Z, and GROUND STATE Mass
	/** Doesn't touch kinematic properties */
	void SetNucleus(const G4String& symbol);
	/// Set GROUND state mass
	/** Doesn't touch any kinematic properties, or A, Z */
	void SetMass(G4double mass);
	void SetEx(G4double ex);
	/// Set mass number
	/** Doesn't touch any kinematic properties, or the mass */
	void SetA(G4int a) { fA = a; }
	/// Set charge
	/** Doesn't touch any kinematic properties, or the mass */
	void SetZ(G4int z) { fZ = z; }
	/// Mass number
	G4int    A()         const { return fA; }
	/// Charge
	G4int    Z()         const { return fZ; }
	/// GROUND STATE mass
	G4double M()         const { return fM;      }
	/// GROUND STATE mass^2
	G4double M2()        const { return fM*fM;   }
	/// Excitation energy
	G4double Ex()        const { return fEx;     }
	/// TOTAL intrinsic energy (or rest energy),
	/// e.g. g.s. mass + excitation
	G4double MplusEx()   const { return fM+fEx;}
	/// GROUND STATE Mass in AMU
	G4double M_amu()     const { return fM*CLHEP::amu_c2; }
	/// TOTAL rest mass (g.s.+ex) in AMU
	G4double MplusEx_amu() const { return MplusEx()*CLHEP::amu_c2; }
	/// TOTAL ENERGY, g.s. + ex + kinetic
	G4double E()         const { return fP.e();  }
	/// Magnitude of momentum
	G4double P()         const { return fP.vect().mag();  }
	/// X-dir momentum
	G4double Px()        const { return fP.px(); }
	/// Y-dir momentum
	G4double Py()        const { return fP.py(); }
	/// Z-dir momentum
	G4double Pz()        const { return fP.pz(); }
	/// Kinetic energy
	G4double Ekin()      const { return fP.e() - fP.m(); }
	/// Theta angle of momentum, relative to z-axis
	G4double Theta()     const { return fP.theta(); }
	/// Phi angle of momentum, relative to x-axis
	G4double Phi()       const { return fP.phi();   }
	/// Get the momentum 4-vector
	/** \attention The mass stored in the four-vector is the
	 *  TOTAL rest mass (g.s. + ex). It is not the same as this->M() !!
	 */
	const G4LorentzVector& Momentum() const { return fP; }
	/// Set 3-momentum vector
	/** Changes kinematics only, rest energy stays the same */
	void SetP3(const G4ThreeVector& p);
	/// Set 3-momentum vector, from components
	/** Changes kinematics only, rest energy stays the same */
	void SetP3XYZ(G4double px, G4double py, G4double pz);
	/// Set 3-momentum vector, from theta and phi angles + momentum magnitude
	/** Changes kinematics only, rest energy stays the same */
	void SetP3ThetaPhi(G4double p, G4double Theta, G4double phi);
	/// Set 3-momentum vector, from theta and phi angles + kinetic energy
	/** Changes kinematics only, rest energy stays the same */
	void SetEkinThetaPhi(G4double ekin, G4double Theta, G4double phi);
	/// Boost into a different frame
	void Boost(const G4ThreeVector& bv) { fP.boost(bv); }
	/// Get the position of the particle
	const G4ThreeVector& Position() const { return fPos; }
	/// X-position
	G4double PosX()   const { return fP.px(); }
	/// Y-position
	G4double PosY()   const { return fP.py(); }
	/// Z-position
	G4double PosZ()   const { return fP.pz(); }
	/// Set position
	void SetPosition(const G4ThreeVector& pos) { fPos = pos; }
	/// Set position, by components
	void SetPosition(G4double x, G4double y, G4double z) { fPos.set(x,y,z); }
	/// Set x-position
	void SetPosX(G4double x) { fP[0] = x; }
	/// Set y-position
	void SetPosY(G4double y) { fP[1] = y; }
	/// Set z-position
	void SetPosZ(G4double z) { fP[2] = z; }
	
private:
	G4int fA, fZ;
	G4double fM, fEx;
	G4LorentzVector fP;
	G4ThreeVector fPos;
};


#endif
