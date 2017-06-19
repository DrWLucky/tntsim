#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4ThreeVector.hh"
#include "TntNuclearMasses.hh"
#include "TntReaction.hh"


TntReaction::TntReaction():
	fAngdist(0), fEx(0)
{
	fThetaCM = 0;
	fZ1 = fZ2 = fZ3 = fZ4 = 0;
	fA1 = fA2 = fA3 = fA4 = 0;
	fM1 = fM2 = fM3 = fM4 = fEbeam = 0;
	fEjectile.set(0,0,0,0);
	fRecoil.set(0,0,0,0);
}

TntReaction::TntReaction(G4int beamZ,   G4int beamA,
												 G4int targetZ, G4int targetA,
												 G4int recoilZ, G4int recoilA,
												 G4double ebeamPerA,
												 G4double exciteRecoil,
												 G4double widthRecoil,
												 const G4String& angDistFile):
	fAngdist(0), fEx(0)
{
	fThetaCM = 0;
	SetInputs(beamZ, beamA, targetZ, targetA, recoilZ, recoilA,
						ebeamPerA, exciteRecoil, widthRecoil, angDistFile);
}

TntReaction::~TntReaction()
{ }

void TntReaction::SetInputs(G4int beamZ,   G4int beamA,
														G4int targetZ, G4int targetA,
														G4int recoilZ, G4int recoilA,
														G4double ebeamPerA,
														G4double exciteRecoil,
														G4double widthRecoil,
														const G4String& angDistFile)
{
	fZ1 = beamZ;   fA1 = beamA;
	fZ2 = targetZ; fA2 = targetA;
	fZ4 = recoilZ; fA4 = recoilA;
	fZ3 = beamZ+targetZ-recoilZ;
	fA3 = beamA+targetA-recoilA;
	
	fM1 = TntNuclearMasses::GetNuclearMass(beamZ, beamA) * MeV;
	fM2 = TntNuclearMasses::GetNuclearMass(targetZ, targetA) * MeV;
	fM3 = TntNuclearMasses::GetNuclearMass(beamZ+targetZ-recoilZ, beamA+targetA-recoilA) * MeV;
	fM4 = TntNuclearMasses::GetNuclearMass(recoilZ, recoilA) * MeV;
	fEbeam = ebeamPerA*beamA;
	fEjectile.set(0,0,0,0);
	fRecoil.set(0,0,0,0);

	fEx.reset(new TntRngBreitWigner(exciteRecoil, widthRecoil));
	fAngdist.reset(new TntRngCustomAngDist(angDistFile));
}

namespace {
inline G4LorentzVector createLorentzVector(G4double mass, G4double ekin, 
																					 G4double theta, G4double phi)
{
	G4double etot = mass + ekin;
	G4double ptot = sqrt(etot*etot - mass*mass);
	return G4LorentzVector(ptot*sin(theta)*cos(phi),
												 ptot*sin(theta)*sin(phi),
												 ptot*cos(theta),
												 etot);
} }


G4bool TntReaction::Generate() 
{
	G4LorentzVector v1 = ::createLorentzVector(fM1, fEbeam, 0, 0);
	G4LorentzVector v2 = ::createLorentzVector(fM2, 0, 0, 0);
	G4LorentzVector vinv = v1+v2;
	G4ThreeVector bv = vinv.boostVector();

	G4LorentzVector v1cm = v1.boost(-bv);
	G4LorentzVector v2cm = v1.boost(-bv);

	G4double ex = fEx->GenerateAbove(0);
	
	if(fM1+fM2+fEbeam < fM3+fM4+ex) {
		G4cerr << "ERROR:: TntReaction::SetInputs:: Not enough energy for reaction!" << G4endl;
		fEjectile.set(0,0,0,0);
		fRecoil.set(0,0,0,0);
		return false;
	}

	G4double S = vinv.m2(); // invariant mass squared
	G4double E3 = fM3, E4 = fM4+ex;
	G4double pcm = sqrt((pow(S - E3*E3 - E4*E4, 2) - 4*E3*E3*E4*E4) / (4*S));

	fThetaCM = fAngdist->Generate() * CLHEP::pi / 180;
	G4double phiCM = TntRngUniform(0, 2*CLHEP::pi).Generate();
	fEjectile = ::createLorentzVector(fM3, sqrt(E3*E3 + pcm*pcm) - E3, 
																		CLHEP::pi - fThetaCM, phiCM);
	fRecoil = ::createLorentzVector(E4, sqrt(E4*E4 + pcm*pcm) - E4,
																	fThetaCM, CLHEP::pi - phiCM);
	fEjectile.boost(bv);
	fRecoil.boost(bv);

	fBeam = v1;
	fTarget = v2;
	
	return true;
}																	


TntReaction* TntReactionFactory::CreateReaction()
{
	if(IsComplete() == false) {
		G4cerr << "ERROR:: TntReactionFactory::CreateReaction():: Incomplete"
					 << " reaction specification" << G4endl;
		G4cerr << "What is Set:: \n\t";
		for(int i=0; i< 9; ++i) { G4cerr<<fIsSet[i]<<", "; }
		G4cerr<<G4endl;
		throw false;
		return 0;
	}
	TntReaction* r = 	new TntReaction(fZ1, fA1, fZ2, fA2, fZ4, fA4, 
																		fEbeamPerA, fEx, fWidth,
																		fAngDistFile);	
	return r;
}

G4bool TntReactionFactory::IsComplete()
{
	for(G4bool* p = fIsSet; p< fIsSet+10; ++p) {
		if(*p == false) { return false; }
	}
	return true;
}


namespace {

G4int get_z_from_symbol(G4String symbol)
{
	const G4String elements[] = {
		"n", "H","He","Li","Be","B","C","N","O","F","Ne",
		"Na","Mg","Al","Si","P","S","Cl","Ar","K","Ca",
		"Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn",
		"Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr",
		"Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn",
		"Sb","Te","I","Xe","Cs","Ba","La","Ce","Pr","Nd",
		"Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb",
		"Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg",
		"Tl","Pb","Bi","Po","At","Rn","Fr","Ra","Ac","Th",
		"Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm",
		"Md","No","Lr","Rf","Db","Sg","Bh","Hs","Mt","Ds",
		"Rg","Cn","Nh","Fl","Mc","Lv","Ts","Og"
	};	
	G4int N = sizeof(elements)/sizeof(elements[0]);
	const G4String* p = std::find(elements, elements+N, symbol);
	if(p< elements + N) return p-elements; return -1;
}

void convert_symbol(G4String* symbol)
{
	if(0) { }
	else if(*symbol == "n") *symbol = "1n";
	else if(*symbol == "p") *symbol = "1H";
	else if(*symbol == "d") *symbol = "2H";
	else if(*symbol == "t") *symbol = "3H";
	else if(*symbol == "a") *symbol = "4He";
	else  { }
}

void parse_symbol(G4String symbol, G4int* A, G4int* Z)
{
	convert_symbol(&symbol);
	const char integers[] = {'0','1','2','3','4','5','6','7','8','9'};
	G4int ni = sizeof(integers) / sizeof(integers[0]);
	G4String A_str = "";
	G4int indx = 0;
	while(1) {
		const char *p = std::find(integers,integers+ni,symbol[indx++]);
		if(p<integers+ni) { A_str.push_back(*p); }
		else break;
		if(indx == ni) break;
	}
	*A = atoi(A_str.c_str());
	*Z = get_z_from_symbol(symbol.substr(indx-1));	
}
}

void TntReactionFactory::SetBeam(G4String beam)
{
	G4int A, Z;
	parse_symbol(beam, &A, &Z);
	SetZ1(Z);
	SetA1(A);
}

void TntReactionFactory::SetTarget(G4String target)
{
	G4int A, Z;
	parse_symbol(target, &A, &Z);
	SetZ2(Z);
	SetA2(A);
}

void TntReactionFactory::SetRecoil(G4String recoil)
{
	G4int A, Z;
	parse_symbol(recoil, &A, &Z);
	SetZ4(Z);
	SetA4(A);
}

void TntReactionFactory::SetEjectile(G4String ejectile)
{
	G4int A, Z;
	parse_symbol(ejectile, &A, &Z);
	SetZ3(Z);
	SetA3(A);
}
