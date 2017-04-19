//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
// $Id: TntDetectorConstruction.cc 82853 2014-07-14 09:07:11Z gcosmo $
//
/// \file optical/Tnt/src/TntDetectorConstruction.cc
/// \brief Implementation of the TntDetectorConstruction class
//
//
//



//////////////////////////////////////////////////////////////////////
//Comments by Shuya
//Major change on 2016/4/21
//BC505 -> BC519
//////////////////////////////////////////////////////////////////////



#include "TntDetectorConstruction.hh"
#include "TntPMTSD.hh"
#include "TntScintSD.hh"
#include "TntDetectorMessenger.hh"
#include "TntMainVolume.hh"
#include "TntWLSSlab.hh"
#include "TntGlobalParams.hh"

#include "G4SDManager.hh"
#include "G4RunManager.hh"

#include "G4GeometryManager.hh"
#include "G4SolidStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4LogicalSkinSurface.hh"

#include "G4OpticalSurface.hh"
#include "G4MaterialTable.hh"
#include "G4VisAttributes.hh"
#include "G4Material.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4PVPlacement.hh"
#include "globals.hh"
#include "G4UImanager.hh"
#include "G4PhysicalConstants.hh"
#include "G4SystemOfUnits.hh"

G4bool TntDetectorConstruction::fSphereOn = true;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntDetectorConstruction::TntDetectorConstruction(G4String Light)
: fTnt_mt(NULL), fMPTPStyrene(NULL), Light_Conv_Method(Light)
{
  fExperimentalHall_box = NULL;
  fExperimentalHall_log = NULL;
  fExperimentalHall_phys = NULL;

  fTnt = fAl = fAir = fVacuum = fGlass = NULL;
  fPstyrene = fPMMA = fPethylene1 = fPethylene2 = NULL;

  fN = fO = fC = fH = NULL;

  SetDefaults();

  fDetectorMessenger = new TntDetectorMessenger(this);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntDetectorConstruction::~TntDetectorConstruction() {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//
//by Shuya 160406
// calculate mass fraction from atomic ratios
inline void calcMassFraction(G4double atomicRatio, // ratio N1/N2 
														 G4double mass1, G4double mass2,
														 G4double* fraction1, G4double* fraction2)
{
	G4double af1 = atomicRatio / (1 + atomicRatio); // atomic fraction
	G4double af2 = 1. - af1; //atomic fraction

	G4double p1 = mass1*af1;
	G4double p2 = mass2*af2;

	*fraction1 = p1 / (p1 + p2);
	*fraction2 = p2 / (p1 + p2);
}

inline G4Material* createHydrocarbon(
	const char* name,
	G4double density,
	G4double HtoC_ratio,
	G4Element* pH,
	G4Element* pC)
{
	G4Material* hc = new G4Material(name, density, 2);
		
	G4double mf1, mf2;
	calcMassFraction(HtoC_ratio, pH->GetN(), pC->GetN(), &mf1, &mf2);
	
//Comment by Shuya 160523. AddElement(A,B). If B is int, number of atoms, if it's double, mass fraction. (A is pointer). So, this case indicated by mass fraction. 
	hc->AddElement(pH, mf1);
	hc->AddElement(pC, mf2);

	return hc;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::DefineMaterials(){
  G4double a;  // atomic mass
  G4double z;  // atomic number
  G4double density;

  G4int polyPMMA = 1;
  G4int nC_PMMA = 3+2*polyPMMA;
  G4int nH_PMMA = 6+2*polyPMMA;

  G4int polyeth = 1;
  G4int nC_eth = 2*polyeth;
  G4int nH_eth = 4*polyeth;

  //***Elements
  fH = new G4Element("H", "H", z=1., a=1.01*g/mole);
  fC = new G4Element("C", "C", z=6., a=12.01*g/mole);
  fN = new G4Element("N", "N", z=7., a= 14.01*g/mole);
  fO = new G4Element("O"  , "O", z=8., a= 16.00*g/mole);

  //***Materials
  //Liquid Xenon
//by Shuya 160406. Just set the name of BC505 to Tnt because it is easier (can use all parameters (reflectivity, etc).
  //fTnt = new G4Material("Tnt",z=54.,a=131.29*g/mole,density=3.020*g/cm3);
  //G4Material* BC505 = 
//by Shuya 160421. 
  //G4Material* BC519 = 
  G4Material* fTnt = 
//////////////////////////// Comment By Shuya 160525. THIS IS TO CHANGE FOR SCINTILLATION MATERIALS (1/8) /////////////////////////////////
		//createHydrocarbon("BC505", 0.877*g/cm3, 1.331, fH, fC);
		//createHydrocarbon("Tnt", 0.893*g/cm3, 1.331, fH, fC);
		//Saint GoBain value
		//createHydrocarbon("Tnt", 0.877*g/cm3, 1.331, fH, fC);
		//by Shuya 160421. BC519.
		//createHydrocarbon("Tnt", 0.875*g/cm3, 1.728, fH, fC);
		//by Shuya 160512. BC404. (# of C:H = 4.68:5.15 -> H/C = 1.100
		createHydrocarbon("Tnt", 1.023*g/cm3, 1.100, fH, fC);
		//by Shuya 160523. EJ309.
		//createHydrocarbon("Tnt", 0.959*g/cm3, 1.25, fH, fC);

  //Aluminum
  fAl = new G4Material("Al",z=13.,a=26.98*g/mole,density=2.7*g/cm3);
  //Vacuum
  fVacuum = new G4Material("Vacuum",z=1.,a=1.01*g/mole,
                          density=universe_mean_density,kStateGas,0.1*kelvin,
                          1.e-19*pascal);

  //Air
  fAir = new G4Material("Air", density= 1.29*mg/cm3, 2);
  fAir->AddElement(fN, 70*perCent);
  fAir->AddElement(fO, 30*perCent);
  //Glass
  fGlass = new G4Material("Glass", density=1.032*g/cm3,2);
  fGlass->AddElement(fC,91.533*perCent);
  fGlass->AddElement(fH,8.467*perCent);
  //Polystyrene
  fPstyrene = new G4Material("Polystyrene", density= 1.03*g/cm3, 2);
  fPstyrene->AddElement(fC, 8);
  fPstyrene->AddElement(fH, 8);
  //Fiber(PMMA)
  fPMMA = new G4Material("PMMA", density=1190*kg/m3,3);
  fPMMA->AddElement(fH,nH_PMMA);
  fPMMA->AddElement(fC,nC_PMMA);
  fPMMA->AddElement(fO,2);
  //Cladding(polyethylene)
  fPethylene1 = new G4Material("Pethylene1", density=1200*kg/m3,2);
  fPethylene1->AddElement(fH,nH_eth);
  fPethylene1->AddElement(fC,nC_eth);
  //Double cladding(flourinated polyethylene)
  fPethylene2 = new G4Material("Pethylene2", density=1400*kg/m3,2);
  fPethylene2->AddElement(fH,nH_eth);
  fPethylene2->AddElement(fC,nC_eth);
//by Shuya 160413.
  //Plexiglass
  /***** REMOVE THIS BECAUSE PMMA = plexiglass.... 
  fPlexiglass = new G4Material("Plexiglass", density=1.19*g/cm3,3);
  fPlexiglass->AddElement(fC,59.9848*perCent);
  fPlexiglass->AddElement(fH,8.0538*perCent);
  fPlexiglass->AddElement(fO,31.9614*perCent);
*//////
 
  //***Material properties tables

//////////////////////////// Comment By Shuya 160525. THIS IS TO CHANGE FOR SCINTILLATION MATERIALS (2/8) /////////////////////////////////
  //G4double lxe_Energy[]    = { 7.0*eV , 7.07*eV, 7.14*eV };
//by Shuya 160414 for BC505, BC519
  //G4double lxe_Energy[]    = { 2.4797*eV , 2.9173*eV, 3.0613*eV };	//500nm (10% of peak), 425nm (peak), 405nm (10% of peak)
//by Shuya 160512 for BC404
  G4double lxe_Energy[]    = { 2.5830*eV , 3.0388*eV, 3.2627*eV };	//500nm (10% of peak), 425nm (peak), 405nm (10% of peak)
//by Shuya 160523 for EJ309
  //G4double lxe_Energy[]    = { 2.4075*eV , 2.9173*eV, 3.2204*eV };	//515nm (10% of peak), 425nm (peak), 385nm (10% of peak)
  const G4int lxenum = sizeof(lxe_Energy)/sizeof(G4double);

//////////////////////////// Comment By Shuya 160525. THIS IS TO CHANGE FOR SCINTILLATION MATERIALS (3/8) /////////////////////////////////
//Comment by Shuya. I didn't change these values because it almost reflects correct values of BC505 and BC519 spectra.
//  G4double lxe_SCINT[] = { 0.1, 1.0, 0.1 };
  //by Shuya 160512 for BC404
  G4double lxe_SCINT[] = { 0.05, 1.0, 0.05 };
  assert(sizeof(lxe_SCINT) == sizeof(lxe_Energy));


//////////////////////////// Comment By Shuya 160525. THIS IS TO CHANGE FOR SCINTILLATION MATERIALS (4/8) /////////////////////////////////
  //G4double lxe_RIND[]  = { 1.59 , 1.57, 1.54 };
//by Shuya 160414 for BC505
  //G4double lxe_RIND[]  = { 1.505 , 1.505, 1.505 };
//by Shuya 160414 for BC519
  //G4double lxe_RIND[]  = { 1.50 , 1.50, 1.50 };
//by Shuya 160512 for BC404
  G4double lxe_RIND[]  = { 1.58 , 1.58, 1.58 };
//by Shuya 160523 for EJ309
//  G4double lxe_RIND[]  = { 1.57 , 1.57, 1.57 };
//by Shuya 160526 for TEST of dependence on photon (angular) distribution (should be independent because of isotropic emission...)
//  G4double lxe_RIND[]  = { 2.0 , 2.0, 2.0 };
  assert(sizeof(lxe_RIND) == sizeof(lxe_Energy));


//////////////////////////// Comment By Shuya 160525. THIS IS TO CHANGE FOR SCINTILLATION MATERIALS (5/8) /////////////////////////////////
//by Shuya 160414 for BC505
  //G4double lxe_ABSL[]  = { 35.*cm, 35.*cm, 35.*cm};
  //G4double lxe_ABSL[]  = { 300.*cm, 300.*cm, 300.*cm};
//by Shuya 160421 for BC519
  //G4double lxe_ABSL[]  = { 100.*cm, 100.*cm, 100.*cm};
//by Shuya 160512 for BC404
	G4double lxe_ABSL[]  = { 160.*cm, 160.*cm, 160.*cm};
//by Shuya 160523 for EJ309
  //G4double lxe_ABSL[]  = { 100.*cm, 100.*cm, 100.*cm};
  assert(sizeof(lxe_ABSL) == sizeof(lxe_Energy));

  fTnt_mt = new G4MaterialPropertiesTable();
  fTnt_mt->AddProperty("FASTCOMPONENT", lxe_Energy, lxe_SCINT, lxenum);
  fTnt_mt->AddProperty("SLOWCOMPONENT", lxe_Energy, lxe_SCINT, lxenum);
  fTnt_mt->AddProperty("RINDEX",        lxe_Energy, lxe_RIND,  lxenum);
  fTnt_mt->AddProperty("ABSLENGTH",     lxe_Energy, lxe_ABSL,  lxenum);


//////////////////////////// Comment By Shuya 160525. THIS IS TO CHANGE FOR SCINTILLATION MATERIALS (6/8) /////////////////////////////////
//by Shuya 160506. To take into account the PMT efficiency (20%) with errors. 
  //fTnt_mt->AddConstProperty("SCINTILLATIONYIELD",(12000.*0.2)/MeV);
//Comment By Shuya 160512. Scintillation Yield: BC505=12000, BC519:9500, BC404=10400, EJ309=11500 (From Ejen catalogue). Anthracene~15000.
  //fTnt_mt->AddConstProperty("SCINTILLATIONYIELD",(9500.*0.2)/MeV);
	G4double nphot = 10400.*0.2;
	fTnt_mt->AddConstProperty("SCINTILLATIONYIELD", nphot/MeV);

  // fTnt_mt->AddConstProperty("SCINTILLATIONYIELD",(10400.*0.2)/MeV);

  //fTnt_mt->AddConstProperty("SCINTILLATIONYIELD",(11500.*0.2)/MeV);
  //fTnt_mt->AddConstProperty("SCINTILLATIONYIELD",(11500.)/MeV);

	G4double resScale = TntGlobalParams::Instance()->GetPhotonResolutionScale();

  fTnt_mt->AddConstProperty("RESOLUTIONSCALE", resScale);


//////////////////////////// Comment By Shuya 160525. THIS IS TO CHANGE FOR SCINTILLATION MATERIALS (7/8) /////////////////////////////////
  //fTnt_mt->AddConstProperty("FASTTIMECONSTANT",20.*ns);
//by Shuya 160414 for BC505
  //fTnt_mt->AddConstProperty("FASTTIMECONSTANT",2.7*ns);
//by Shuya 160421 for BC519
  //fTnt_mt->AddConstProperty("FASTTIMECONSTANT",4.*ns);
//by Shuya 160512 for BC404 (not sure about slow one)
	fTnt_mt->AddConstProperty("FASTSCINTILLATIONRISETIME",0.9*ns); // RISE TIME
  fTnt_mt->AddConstProperty("FASTTIMECONSTANT",2.1*ns); // DECAY TIME

//by Shuya 160523 for EJ309 (not sure about slow one)
  //fTnt_mt->AddConstProperty("FASTTIMECONSTANT",3.5*ns);

//Comment by Shuya 160414. All is Fast component! (Fast/(Fast+Slow) = 1.0),
//So SLOWTIMECONSTANT is actually not needed.
//Except it must be defined sice we added the SLOWCOMPONENT property a few lines
//above. Otherwise the code crashes. But whatever values we set get ignored
  fTnt_mt->AddConstProperty("SLOWTIMECONSTANT",0.9*ns);
	fTnt_mt->AddConstProperty("SLOWSCINTILLATIONRISETIME",2.1*ns);
  fTnt_mt->AddConstProperty("YIELDRATIO",1.0); // Ratio of fast / (fast+slow) ==> 1.0 means all fast
  fTnt->SetMaterialPropertiesTable(fTnt_mt);


  // Set the Birks Constant for the Tnt scintillator

  fTnt->GetIonisation()->SetBirksConstant(0.126*mm/MeV);
 
//by Shuya 160606. To test if the refractive index affects the photon detection at PMTs.
  G4double glass_RIND[]={1.49,1.49,1.49};
  //G4double glass_RIND[]={1.57,1.57,1.57};
  assert(sizeof(glass_RIND) == sizeof(lxe_Energy));
  G4double glass_AbsLength[]={420.*cm,420.*cm,420.*cm};
  assert(sizeof(glass_AbsLength) == sizeof(lxe_Energy));
  G4MaterialPropertiesTable *glass_mt = new G4MaterialPropertiesTable();
  glass_mt->AddProperty("ABSLENGTH",lxe_Energy,glass_AbsLength,lxenum);
  glass_mt->AddProperty("RINDEX",lxe_Energy,glass_RIND,lxenum);
  fGlass->SetMaterialPropertiesTable(glass_mt);

  G4double vacuum_Energy[]={2.0*eV,7.0*eV,7.14*eV};
  const G4int vacnum = sizeof(vacuum_Energy)/sizeof(G4double);
  G4double vacuum_RIND[]={1.,1.,1.};
  assert(sizeof(vacuum_RIND) == sizeof(vacuum_Energy));
  G4MaterialPropertiesTable *vacuum_mt = new G4MaterialPropertiesTable();
  vacuum_mt->AddProperty("RINDEX", vacuum_Energy, vacuum_RIND,vacnum);
  fVacuum->SetMaterialPropertiesTable(vacuum_mt);
  fAir->SetMaterialPropertiesTable(vacuum_mt);//Give air the same rindex

  G4double wls_Energy[] = {2.00*eV,2.87*eV,2.90*eV,3.47*eV};
  const G4int wlsnum = sizeof(wls_Energy)/sizeof(G4double);
 
  G4double rIndexPstyrene[]={ 1.5, 1.5, 1.5, 1.5};
  assert(sizeof(rIndexPstyrene) == sizeof(wls_Energy));
  G4double absorption1[]={2.*cm, 2.*cm, 2.*cm, 2.*cm};
  assert(sizeof(absorption1) == sizeof(wls_Energy));
  G4double scintilFast[]={0.00, 0.00, 1.00, 1.00};
  assert(sizeof(scintilFast) == sizeof(wls_Energy));
  fMPTPStyrene = new G4MaterialPropertiesTable();
  fMPTPStyrene->AddProperty("RINDEX",wls_Energy,rIndexPstyrene,wlsnum);
  fMPTPStyrene->AddProperty("ABSLENGTH",wls_Energy,absorption1,wlsnum);
  fMPTPStyrene->AddProperty("FASTCOMPONENT",wls_Energy, scintilFast,wlsnum);
  fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);
  fMPTPStyrene->AddConstProperty("RESOLUTIONSCALE",1.0);
  fMPTPStyrene->AddConstProperty("FASTTIMECONSTANT", 10.*ns);
  fPstyrene->SetMaterialPropertiesTable(fMPTPStyrene);

  // Set the Birks Constant for the Polystyrene scintillator

  fPstyrene->GetIonisation()->SetBirksConstant(0.126*mm/MeV);

  G4double RefractiveIndexFiber[]={ 1.60, 1.60, 1.60, 1.60};
  assert(sizeof(RefractiveIndexFiber) == sizeof(wls_Energy));
  G4double AbsFiber[]={9.00*m,9.00*m,0.1*mm,0.1*mm};
  assert(sizeof(AbsFiber) == sizeof(wls_Energy));
  G4double EmissionFib[]={1.0, 1.0, 0.0, 0.0};
  assert(sizeof(EmissionFib) == sizeof(wls_Energy));
  G4MaterialPropertiesTable* fiberProperty = new G4MaterialPropertiesTable();
  fiberProperty->AddProperty("RINDEX",wls_Energy,RefractiveIndexFiber,wlsnum);
  fiberProperty->AddProperty("WLSABSLENGTH",wls_Energy,AbsFiber,wlsnum);
  fiberProperty->AddProperty("WLSCOMPONENT",wls_Energy,EmissionFib,wlsnum);
  fiberProperty->AddConstProperty("WLSTIMECONSTANT", 0.5*ns);
  fPMMA->SetMaterialPropertiesTable(fiberProperty);

  G4double RefractiveIndexClad1[]={ 1.49, 1.49, 1.49, 1.49};
  assert(sizeof(RefractiveIndexClad1) == sizeof(wls_Energy));
  G4MaterialPropertiesTable* clad1Property = new G4MaterialPropertiesTable();
  clad1Property->AddProperty("RINDEX",wls_Energy,RefractiveIndexClad1,wlsnum);
  clad1Property->AddProperty("ABSLENGTH",wls_Energy,AbsFiber,wlsnum);
  fPethylene1->SetMaterialPropertiesTable(clad1Property);

  G4double RefractiveIndexClad2[]={ 1.42, 1.42, 1.42, 1.42};
  assert(sizeof(RefractiveIndexClad2) == sizeof(wls_Energy));
  G4MaterialPropertiesTable* clad2Property = new G4MaterialPropertiesTable();
  clad2Property->AddProperty("RINDEX",wls_Energy,RefractiveIndexClad2,wlsnum);
  clad2Property->AddProperty("ABSLENGTH",wls_Energy,AbsFiber,wlsnum);
  fPethylene2->SetMaterialPropertiesTable(clad2Property);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* TntDetectorConstruction::Construct(){

  if (fExperimentalHall_phys) {
     G4GeometryManager::GetInstance()->OpenGeometry();
     G4PhysicalVolumeStore::GetInstance()->Clean();
     G4LogicalVolumeStore::GetInstance()->Clean();
     G4SolidStore::GetInstance()->Clean();
     G4LogicalSkinSurface::CleanSurfaceTable();
     G4LogicalBorderSurface::CleanSurfaceTable();
  }

  DefineMaterials();
  return ConstructDetector();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* TntDetectorConstruction::ConstructDetector()
{
  //The experimental hall walls are all 1m away from housing walls
  G4double expHall_x = fScint_x+fD_mtl+1.*m;
  G4double expHall_y = fScint_y+fD_mtl+1.*m;
  G4double expHall_z = fScint_z+fD_mtl+1.*m;

//by Shuya 160404
/*
G4cout << fScint_x << G4endl;
G4cout << fScint_y << G4endl;
G4cout << fScint_z << G4endl;
G4cout << fNx << G4endl;
G4cout << fNy << G4endl;
G4cout << fNz << G4endl;
*/

  //Create experimental hall
  fExperimentalHall_box
    = new G4Box("expHall_box",expHall_x,expHall_y,expHall_z);
  fExperimentalHall_log = new G4LogicalVolume(fExperimentalHall_box,
                                             fVacuum,"expHall_log",0,0,0);
  fExperimentalHall_phys = new G4PVPlacement(0,G4ThreeVector(),
                              fExperimentalHall_log,"expHall",0,false,0);

  fExperimentalHall_log->SetVisAttributes(G4VisAttributes::Invisible);

  //Place the main volume
  if(fMainVolumeOn){
    fMainVolume
      = new TntMainVolume(0,G4ThreeVector(),fExperimentalHall_log,false,0,this);
  }

  //Place the WLS slab
  if(fWLSslab){
    G4VPhysicalVolume* slab = new TntWLSSlab(0,G4ThreeVector(0.,0.,
                                             -fScint_z/2.-fSlab_z-1.*cm),
                                             fExperimentalHall_log,false,0,
                                             this);

    //Surface properties for the WLS slab
    G4OpticalSurface* scintWrap = new G4OpticalSurface("ScintWrap");
 
    new G4LogicalBorderSurface("ScintWrap", slab,
                               fExperimentalHall_phys,
                               scintWrap);
 
    scintWrap->SetType(dielectric_metal);
    scintWrap->SetFinish(polished);
    scintWrap->SetModel(glisur);

    G4double pp[] = {2.0*eV, 3.5*eV};
    const G4int num = sizeof(pp)/sizeof(G4double);
    G4double reflectivity[] = {1., 1.};
    assert(sizeof(reflectivity) == sizeof(pp));
    G4double efficiency[] = {0.0, 0.0};
    assert(sizeof(efficiency) == sizeof(pp));
    
    G4MaterialPropertiesTable* scintWrapProperty 
      = new G4MaterialPropertiesTable();

    scintWrapProperty->AddProperty("REFLECTIVITY",pp,reflectivity,num);
    scintWrapProperty->AddProperty("EFFICIENCY",pp,efficiency,num);
    scintWrap->SetMaterialPropertiesTable(scintWrapProperty);
  }

  return fExperimentalHall_phys;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::ConstructSDandField() {

  if (!fMainVolume) return;

  // PMT SD

  if (!fPmt_SD.Get()) {
    //Created here so it exists as pmts are being placed
    G4cout << "Construction /TntDet/pmtSD" << G4endl;
    TntPMTSD* pmt_SD = new TntPMTSD("/TntDet/pmtSD");
    fPmt_SD.Put(pmt_SD);

    pmt_SD->InitPMTs((fNx*fNy+fNx*fNz+fNy*fNz)*2); //let pmtSD know # of pmts
    pmt_SD->SetPmtPositions(fMainVolume->GetPmtPositions());
  }

  //sensitive detector is not actually on the photocathode.
  //processHits gets done manually by the stepping action.
  //It is used to detect when photons hit and get absorbed&detected at the
  //boundary to the photocathode (which doesnt get done by attaching it to a
  //logical volume.
  //It does however need to be attached to something or else it doesnt get
  //reset at the begining of events

  SetSensitiveDetector(fMainVolume->GetLogPhotoCath(), fPmt_SD.Get());

  // Scint SD

  if (!fScint_SD.Get()) {
    G4cout << "Construction /TntDet/scintSD" << G4endl;
//by Shuya 160407
    //TntScintSD* scint_SD = new TntScintSD("/TntDet/scintSD");
    TntScintSD* scint_SD = new TntScintSD("/TntDet/scintSD", Light_Conv_Method);
    fScint_SD.Put(scint_SD);
  }
  SetSensitiveDetector(fMainVolume->GetLogScint(), fScint_SD.Get());
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetDimensions(G4ThreeVector dims) {
  this->fScint_x=dims[0];
  this->fScint_y=dims[1];
  this->fScint_z=dims[2];
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}
 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetHousingThickness(G4double d_mtl) {
  this->fD_mtl=d_mtl;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetNX(G4int nx) {
  this->fNx=nx;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetNY(G4int ny) {
  this->fNy=ny;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetNZ(G4int nz) {
  this->fNz=nz;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetPMTRadius(G4double outerRadius_pmt) {
  this->fOuterRadius_pmt=outerRadius_pmt;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//by Shuya 160509
void TntDetectorConstruction::SetPMTSizeX(G4double pmt_x) {
  this->fPmt_x=pmt_x;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//by Shuya 160509
void TntDetectorConstruction::SetPMTSizeY(G4double pmt_y) {
  this->fPmt_y=pmt_y;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetDefaults() {

  //Resets to default values
  fD_mtl=0.0635*cm;

//by Shuya 160414
  //Container_mtl=2.*cm;

//by Shuya 160404
//  fScint_x = 17.8*cm;
//  fScint_y = 17.8*cm;
//  fScint_z = 22.6*cm;

  //fScint_x = 200.0*cm;
  //fScint_y = 200.0*cm;
//by Shuya 160510
  fScint_x = 28.0*cm;
  fScint_y = 28.0*cm;
//by Shuya 160510
//	fScint_z = 30.0*cm;
  fScint_z = 10.0*cm;


//by Shuya 160404
//  fNx = 2;
//  fNy = 2;
//  fNz = 3;
//  by Shuya 160509
  //fNx = 10;
  //fNy = 10;
  extern G4int NX;
  extern G4int NY;

  fNx = NX;
  fNy = NY;
  fNz = 1;


//by Shuya 160404
  //fOuterRadius_pmt = 2.3*cm;
  //fOuterRadius_pmt = 10.*cm;
//by Shuya 160509. Change the PMT shape from Tube to Box.
  //fPmt_x = 20.*cm;
  //fPmt_y = 20.*cm;
  fPmt_x = fScint_x / fNx;
  fPmt_y = fScint_y / fNx;

  //fSphereOn = true;
  fSphereOn = false;
  //fRefl=1.0;
  //by Shuya 160414
  //Comment by Shuya 160502. fRefl might need to be set 1.0 because the (transmission) efficiency must be 0 for photocathod to detect photons (see my explanation of efficiency in MainVolume.cc)
  //Comment by Shuya 160503. Note comment above is wrong. Because efficiency is quantum efficiency and not transmission..., so nothing to do between reflectivity and efficiency.
  fRefl=0.0;

  fNfibers=15;
  fWLSslab=false;
//by Shuya 160414
  //fWLSslab=true;
  fMainVolumeOn=true;
  fMainVolume=NULL;
  fSlab_z=2.5*mm;

  G4UImanager::GetUIpointer()
    ->ApplyCommand("/Tnt/detector/scintYieldFactor 1.");


//////////////////////////// Comment By Shuya 160525. THIS IS TO CHANGE FOR SCINTILLATION MATERIALS (8/8) /////////////////////////////////

//Comment By Shuya 160512. Scintillation Yield: BC505=12000, BC519:9500, BC404=10400, EJ309=11500 (From Ejen catalogue). Anthracene~15000.
  //if(fTnt_mt)fTnt_mt->AddConstProperty("SCINTILLATIONYIELD",(12000.*0.2)/MeV);
  //if(fTnt_mt)fTnt_mt->AddConstProperty("SCINTILLATIONYIELD",(9500.*0.2)/MeV);
  if(fTnt_mt) {
		G4double nphot = 10400.*0.2;
		fTnt_mt->AddConstProperty("SCINTILLATIONYIELD", nphot/MeV);
	}
//  if(fTnt_mt)fTnt_mt->AddConstProperty("SCINTILLATIONYIELD",(11500.*0.2)/MeV);
  //if(fTnt_mt)fTnt_mt->AddConstProperty("SCINTILLATIONYIELD",(11500.)/MeV);

  if(fMPTPStyrene)fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",10./keV);

  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetSphereOn(G4bool b) {
  fSphereOn=b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetHousingReflectivity(G4double r) {
  fRefl=r;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetWLSSlabOn(G4bool b) {
  fWLSslab=b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetMainVolumeOn(G4bool b) {
  fMainVolumeOn=b;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetNFibers(G4int n) {
  fNfibers=n;
  G4RunManager::GetRunManager()->ReinitializeGeometry();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntDetectorConstruction::SetMainScintYield(G4double y) {
  fTnt_mt->AddConstProperty("SCINTILLATIONYIELD",y/MeV);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
 
void TntDetectorConstruction::SetWLSScintYield(G4double y) {
  fMPTPStyrene->AddConstProperty("SCINTILLATIONYIELD",y/MeV);
}

