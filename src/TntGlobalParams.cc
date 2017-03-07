#include "TntGlobalParams.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"



TntGlobalParams::TntGlobalParams(): fNeutronEnergy(1.*MeV),
																		fNumPmtX(4),
																		fNumPmtY(4),
																		fBeamType("pencil"),
																		fRootFileName("TntDataTree.root")
{ }

TntGlobalParams* TntGlobalParams::Instance()
{
	static TntGlobalParams* instance = 0;
	if(!instance) { instance = new TntGlobalParams(); }
	return instance;
}
