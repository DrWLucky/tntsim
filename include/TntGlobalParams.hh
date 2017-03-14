/// \file TntGlobalParams.h
/// \brief Define singleton class containing some global values.
///
#ifndef TNT_GLOBAL_PARAMS_
#define TNT_GLOBAL_PARAMS_
#include "globals.hh"

class TntGlobalParams {
public:
	static TntGlobalParams* Instance();

	G4double GetNeutronEnergy() { return fNeutronEnergy; }
	void SetNeutronEnergy(G4double energy) { fNeutronEnergy = energy; }

	G4int GetNumPmtX() { return fNumPmtX; }
	void SetNumPmtX(G4int n) { fNumPmtX = n; }

	G4int GetNumPmtY() { return fNumPmtY; }
	void SetNumPmtY(G4int n) { fNumPmtY = n; }

	G4String GetBeamType() { return fBeamType; }
	void SetBeamType(G4String type) { fBeamType = type; }

	G4String GetRootFileName() { return fRootFileName; }
	void SetRootFileName(G4String name) { fRootFileName = name; }

	G4double GetPhotonResolutionScale() { return fPhotonResolutionScale; }
	void SetPhotonResolutionScale(G4double scale) { fPhotonResolutionScale = scale; }

	
private:
	TntGlobalParams();
	
private:
	G4double fNeutronEnergy;
	G4int fNumPmtX;
	G4int fNumPmtY;
	G4String fBeamType;
	G4String fRootFileName;
	G4double fPhotonResolutionScale;
};


#endif
