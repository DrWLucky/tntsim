#ifndef TNT_REACTION_HH
#define TNT_REACTION_HH
#include <memory>
#include <vector>
#include "globals.hh"
#include "G4LorentzVector.hh"
#include "TntRng.hh"

// class to calculate nuclear reaction kinematics
class TntReaction {
public:
	TntReaction();
	TntReaction(G4int beamZ,   G4int beamA,
							G4int targetZ, G4int targetA,
							G4int recoilZ, G4int recoilA,
							G4double ebeamPerA,
							G4double exciteRecoil,
							G4double widthRecoil,
							const G4String& angDistFile);
	~TntReaction();

	void SetInputs(G4int beamZ,   G4int beamA,
								 G4int targetZ, G4int targetA,
								 G4int recoilZ, G4int recoilA,
								 G4double ebeamPerA,
								 G4double exciteRecoil, 
								 G4double widthRecoil,
								 const G4String& angDistFile);
	
	G4bool Generate();
	const G4LorentzVector& GetEjectile() const { return fEjectile; }
	const G4LorentzVector& GetRecoil() const { return fRecoil; }
	G4double GetThetaCM() const { return fThetaCM; }
	
	G4int GetZ1() const {return fZ1;}
	G4int GetZ2() const {return fZ2;}
	G4int GetZ3() const {return fZ3;}
	G4int GetZ4() const {return fZ4;}

	G4int GetA1() const {return fA1;}
	G4int GetA2() const {return fA2;}
	G4int GetA3() const {return fA3;}
	G4int GetA4() const {return fA4;}
	
private:
	G4int fZ1, fZ2, fA1, fA2, fZ3, fA3, fZ4, fA4;
	G4double fM1, fM2, fM3, fM4, fEbeam, fThetaCM;
	G4LorentzVector fEjectile, fRecoil;							
	std::auto_ptr<TntRng> fAngdist, fEx;
};


class TntReactionFactory {
public:
	TntReactionFactory() { for(int i=0;i<9;++i){fIsSet[i] = false;} }
	~TntReactionFactory() { }
	void SetZ1(G4int z) { fZ1 = z; fIsSet[0] = true; }
	void SetZ2(G4int z) { fZ2 = z; fIsSet[1] = true; }
	void SetZ3(G4int z) { fZ4 = fZ1+fZ2-z; fIsSet[2] = true; }
	void SetZ4(G4int z) { fZ4 = z; fIsSet[2] = true; }

	void SetA1(G4int a) { fA1 = a; fIsSet[3] = true; }
	void SetA2(G4int a) { fA2 = a; fIsSet[4] = true; }
	void SetA3(G4int a) { fA4 = fA1+fA2-a; fIsSet[5] = true; }
	void SetA4(G4int a) { fA4 = a; fIsSet[5] = true; }

	void SetBeam(G4String beam);
	void SetTarget(G4String target);
	void SetRecoil(G4String recoil);
	void SetEjectile(G4String ejectile);
	
	void SetEbeamPerA(G4double e) { fEbeamPerA = e; fIsSet[6] = true; }

	void SetEx(G4double ex) 
		{ fEx = ex; 
			fIsSet[7] = true; }

	void SetWidth(G4double width)
		{ fWidth = width;
			fIsSet[8] = true; }

	void SetAngDistFile(G4String filename) 
		{ fAngDistFile = filename;
			fIsSet[9] = true;	}

	G4bool IsComplete();
	TntReaction* CreateReaction();

private:
	G4bool fIsSet[10];
	G4int fZ1, fZ2, fA1, fA2, fZ4, fA4;
	G4double fEbeamPerA, fEx, fWidth;
	G4String fAngDistFile;
};


#endif
