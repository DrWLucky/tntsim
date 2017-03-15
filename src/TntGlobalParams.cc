#include <cassert>
#include "TntGlobalParams.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"



TntGlobalParams::TntGlobalParams(): fNeutronEnergy(1.*MeV),
																		fNumPmtX(4),
																		fNumPmtY(4),
																		fBeamType("pencil"),
																		fRootFileName("TntDataTree.root"),
																		fPhotonResolutionScale(1),
																		fMenateR_Tracking(0)
{ }

TntGlobalParams* TntGlobalParams::Instance()
{
	static TntGlobalParams* instance = 0;
	if(!instance) { instance = new TntGlobalParams(); }
	return instance;
}

G4int TntGlobalParams::GetMenateR_Tracking() 
{
	return fMenateR_Tracking; 
}

void TntGlobalParams::SetMenateR_Tracking(G4int n) 
{
	fMenateR_Tracking = n; 
	assert(fMenateR_Tracking == 0 || fMenateR_Tracking == 1 || fMenateR_Tracking == 2);
}
