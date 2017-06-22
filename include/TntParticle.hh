#ifndef TNT_PARTICLE_HEADER
#define TNT_PARTICLE_HEADER
#include <G4PhysicalConstants.hh>
#include <G4ThreeVector.hh>
#include <G4LorentzVector.hh>


class TntParticle {
public:
	TntParticle();
	~TntParticle();
	void SetNucleus(G4int Z, G4int A);  /// Sets A,Z, and Mass
	void SetNucleus(const G4String& symbol);  /// Sets A,Z, and Mass
	void SetMass(G4double mass) { fM = mass; }
	void SetA(G4int a) { fA = a; }
	void SetZ(G4int z) { fZ = z; }
	
	G4int    A()     const { return fA; }
	G4int    Z()     const { return fZ; }
	G4double M()     const { return fM;      }
	G4double M2()    const { return fM*fM;   }
	G4double Amu()   const { return fM*CLHEP::amu_c2; }
	G4double E()     const { return fP.e();  }
	G4double P()     const { return fP.vect().mag();  }
	G4double Px()    const { return fP.px(); }
	G4double Py()    const { return fP.py(); }
	G4double Pz()    const { return fP.pz(); }
	G4double Ekin()  const { return fP.e() - fP.m(); }
	G4double Theta() const { return fP.theta(); }
	G4double Phi()   const { return fP.phi();   }
	const G4LorentzVector& Momentum() const { return fP; }
	void SetP3(const G4ThreeVector& p);
	void SetP3XYZ(G4double px, G4double py, G4double pz);
	void SetP3ThetaPhi(G4double p, G4double Theta, G4double phi);
	void SetEkinThetaPhi(G4double ekin, G4double Theta, G4double phi);
	void Boost(const G4ThreeVector& bv) { fP.boost(bv); }
	
	const G4ThreeVector& Position() const { return fPos; }
	G4double PosX()   const { return fP.px(); }
	G4double PosY()   const { return fP.py(); }
	G4double PosZ()   const { return fP.pz(); }
	void SetPosition(const G4ThreeVector& pos) { fPos = pos; }
	void SetPosition(G4double x, G4double y, G4double z) { fPos.set(x,y,z); }
	
private:
	G4int fA, fZ;
	G4double fM;
	G4LorentzVector fP;
	G4ThreeVector fPos;
};


#endif
