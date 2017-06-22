#include <G4SystemOfUnits.hh>
#include "TntNuclearMasses.hh"
#include "TntParticle.hh"


TntParticle::TntParticle():
	fA(0), fZ(0), fP(0,0,0,0), fPos(0,0,0)
{ }

TntParticle::~TntParticle()
{ }

void TntParticle::SetNucleus(G4int Z, G4int A)
{
	fZ=Z;
	fA=A;
	SetMass(TntNuclearMasses::GetNuclearMass(Z, A)*MeV);
}

void TntParticle::SetNucleus(const G4String& symbol)
{
	TntNuclearMasses::GetZAFromSymbol(symbol, &fZ, &fA);
	SetMass(TntNuclearMasses::GetNuclearMass(fZ, fA)*MeV);
}

void TntParticle::SetP3(const G4ThreeVector& p)
{
	G4double E = sqrt(p.mag2() + fM*fM);
	fP.set(p, E);
}

void TntParticle::SetP3XYZ(G4double px, G4double py, G4double pz)
{
	SetP3(G4ThreeVector(px,py,pz));
}

void TntParticle::SetP3ThetaPhi(G4double p, G4double theta, G4double phi)
{
	SetP3(G4ThreeVector(p*sin(theta)*cos(phi), 
											p*sin(theta)*sin(phi), 
											p*cos(theta)));
}

void TntParticle::SetEkinThetaPhi(G4double ekin, G4double theta, G4double phi)
{
	G4double E = ekin+GetM();
	G4double p = sqrt(E*E - GetM2());
	SetP3ThetaPhi(p,theta,phi);
}
