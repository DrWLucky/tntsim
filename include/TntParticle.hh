#ifndef TNT_PARTICLE_HEADER
#define TNT_PARTICLE_HEADER
#include <G4PhysicalConstants.hh>
#include <G4ThreeVector.hh>
#include <G4LorentzVector.hh>


class TntParticle {
public:
	TntParticle();
	~TntParticle();
	G4int GetA() const { return fA; }
	G4int GetZ() const { return fZ; }
	void SetNucleus(G4int Z, G4int A);
	void SetNucleus(const G4String& symbol);
	void SetMass(G4double mass) { fM = mass; }
	
	const G4LorentzVector& GetMomentum() const { return fP; }
	G4double GetM()     const { return fM;      }
	G4double GetM2()    const { return fM*fM;   }
	G4double GetAmu()   const { return GetM()*CLHEP::amu_c2; }
	G4double GetE()     const { return fP.e();  }
	G4double GetP()     const { return fP.vect().mag();  }
	G4double GetPx()    const { return fP.px(); }
	G4double GetPy()    const { return fP.py(); }
	G4double GetPz()    const { return fP.pz(); }
	G4double GetEkin()  const { return fP.e() - fP.m(); }
	G4double GetTheta() const { return fP.theta(); }
	G4double GetPhi()   const { return fP.phi();   }
	void SetP3(const G4ThreeVector& p);
	void SetP3XYZ(G4double px, G4double py, G4double pz);
	void SetP3ThetaPhi(G4double p, G4double Theta, G4double phi);
	void SetEkinThetaPhi(G4double ekin, G4double Theta, G4double phi);
	void Boost(const G4ThreeVector& bv) { fP.boost(bv); }
	
	const G4ThreeVector& GetPosition() const { return fPos; }
	G4double GetPosX()   const { return fP.px(); }
	G4double GetPosY()   const { return fP.py(); }
	G4double GetPosZ()   const { return fP.pz(); }
	void SetPosition(const G4ThreeVector& pos) { fPos = pos; }

private:
	G4int fA, fZ;
	G4double fM;
	G4LorentzVector fP;
	G4ThreeVector fPos;
};


#endif
