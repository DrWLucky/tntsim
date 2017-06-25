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

	G4String GetReacFile() { return fReacFile; }
	void SetReacFile(G4String type) { fReacFile = type; }

	G4String GetInputFile() { return fInputFile; }
	void SetInputFile(G4String type) { fInputFile = type; }

	G4String GetRootFileName() { return fRootFileName; }
	void SetRootFileName(G4String name) { fRootFileName = name; }

	G4double GetPhotonResolutionScale() { return fPhotonResolutionScale; }
	void SetPhotonResolutionScale(G4double scale) { fPhotonResolutionScale = scale; }

	G4int GetMenateR_Tracking();
	void SetMenateR_Tracking(G4int n);

	G4String GetScintMaterial() { return fScintMaterial; }
	void SetScintMaterial(G4String materialName) { fScintMaterial = materialName; }

	G4double GetDetectorX() { return fDetectorX; }
	void SetDetectorX(G4double x) { fDetectorX = x; }

	G4double GetDetectorY() { return fDetectorY; }
	void SetDetectorY(G4double y) { fDetectorY = y; }

	G4double GetDetectorZ() { return fDetectorZ; }
	void SetDetectorZ(G4double z) { fDetectorZ = z; }

	G4double GetSourceZ() { return fSourceZ; }
	void SetSourceZ(G4double z) { fSourceZ = z; }

	G4int GetLightOutput() { return fLightOutput; }
	void SetLightOutput(G4int n) { fLightOutput = n; }
	
	G4double GetQuantumEfficiency() { return fQuantumEfficiency; }
	void SetQuantumEfficiency(G4double e) { fQuantumEfficiency = e; }

	G4String GetAngerAnalysis() { return fAngerAnalysis; }
	void SetAngerAnalysis(G4String name) { fAngerAnalysis = name; }

	
private:
	TntGlobalParams();
	
private:
	G4double fNeutronEnergy;
	G4int fNumPmtX;
	G4int fNumPmtY;
	G4String fBeamType;
	G4String fReacFile;
	G4String fInputFile;
	G4String fRootFileName;
	G4double fPhotonResolutionScale;
	G4int fMenateR_Tracking;
	G4String fScintMaterial;
	G4double fDetectorX, fDetectorY, fDetectorZ, fSourceZ;
	G4int fLightOutput;
	G4double fQuantumEfficiency;
	G4String fAngerAnalysis;
};


#endif
