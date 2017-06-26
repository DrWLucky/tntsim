/// \file TntGlobalParams.h
/// \brief Define singleton class containing some global values.
///
#ifndef TNT_GLOBAL_PARAMS_
#define TNT_GLOBAL_PARAMS_
#include "globals.hh"

class TntGlobalParams {
public:
	static TntGlobalParams* Instance();

	G4double GetNeutronEnergy() const { return fNeutronEnergy; }
	void SetNeutronEnergy(G4double energy) { fNeutronEnergy = energy; }

	G4int GetNumPmtX() const { return fNumPmtX; }
	void SetNumPmtX(G4int n) { fNumPmtX = n; }

	G4int GetNumPmtY() const { return fNumPmtY; }
	void SetNumPmtY(G4int n) { fNumPmtY = n; }

	G4String GetBeamType() const { return fBeamType; }
	void SetBeamType(G4String type) { fBeamType = type; }

	G4String GetReacFile() const { return fReacFile; }
	void SetReacFile(G4String type) { fReacFile = type; }

	G4String GetInputFile() const { return fInputFile; }
	void SetInputFile(G4String type) { fInputFile = type; }

	G4String GetRootFileName() const { return fRootFileName; }
	void SetRootFileName(G4String name) { fRootFileName = name; }

	G4double GetPhotonResolutionScale() const { return fPhotonResolutionScale; }
	void SetPhotonResolutionScale(G4double scale) { fPhotonResolutionScale = scale; }

	G4int GetMenateR_Tracking();
	void SetMenateR_Tracking(G4int n);

	G4String GetScintMaterial() const { return fScintMaterial; }
	void SetScintMaterial(G4String materialName) { fScintMaterial = materialName; }

	G4double GetDetectorX() const { return fDetectorX; }
	void SetDetectorX(G4double x) { fDetectorX = x; }

	G4double GetDetectorY() const { return fDetectorY; }
	void SetDetectorY(G4double y) { fDetectorY = y; }

	G4double GetDetectorZ() const { return fDetectorZ; }
	void SetDetectorZ(G4double z) { fDetectorZ = z; }

	G4double GetSourceZ() const { return fSourceZ; }
	void SetSourceZ(G4double z) { fSourceZ = z; }

	G4int GetLightOutput() const { return fLightOutput; }
	void SetLightOutput(G4int n) { fLightOutput = n; }
	
	G4double GetQuantumEfficiency() const { return fQuantumEfficiency; }
	void SetQuantumEfficiency(G4double e) { fQuantumEfficiency = e; }

	G4String GetAngerAnalysis() const { return fAngerAnalysis; }
	void SetAngerAnalysis(G4String name) { fAngerAnalysis = name; }

	G4double GetNumDetX() const { return fNdetX; }
	void SetNumDetX(G4int i)  { fNdetX = i; }

	G4double GetNumDetY() const { return fNdetY; }
	void SetNumDetY(G4int i)  { fNdetY = i; }

	void SetNumDetXY(G4int nx, G4int ny)
		{ fNdetX = nx; fNdetY = ny; }
	void GetNumDetXY(G4int& nx, G4int& ny)
		{ nx=fNdetX; ny=fNdetY; }
	
private:
	TntGlobalParams();
	
private:
	G4double fNeutronEnergy;
	G4int fNumPmtX;
	G4int fNumPmtY;
	G4int fNdetX;
	G4int fNdetY;
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
