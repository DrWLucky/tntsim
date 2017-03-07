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
//
// $Id: TntPhysicsList.cc,v 1.6 2007/02/14 17:47:21 roeder Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// Modified by Brian Roeder LPC Caen on 02/14/2007
// email - roeder@lpccaen.in2p3.fr
// 
// -------------------------------------------------------------------
// Based on     GEANT 4 - exampleN01, adding elements from other exps.
// -------------------------------------------------------------------
//
// 11/16/2006 - Added em processes only, for protons, etc.
//            - Added "ConstructEM()" Process to "ConstructProcess()"
//            - Added "G4Proton::ProtonDefinition()" to add proton 
// 12/12/2006 - Adding hadronic processes for proton and neutron
//

//by Shuya 160404
//#include "TntPhysicsList.hh"
#include "TntNuclearPhysics.hh"

#include "G4ParticleTypes.hh"
#include "G4ProcessManager.hh"
#include "G4ParticleDefinition.hh"

// Particle Constructors

#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4Decay.hh"

// Hadron Processes (Add models above to these processes)
#include "G4HadronElasticProcess.hh"
#include "G4HadronFissionProcess.hh"
#include "G4HadronCaptureProcess.hh"

// Scintillation Process by Shuya 160314.
#include "G4Scintillation.hh"
#include "G4Threading.hh"
#include "G4OpAbsorption.hh"
#include "G4OpRayleigh.hh"
#include "G4OpMieHG.hh"
#include "G4OpBoundaryProcess.hh"
#include "G4LossTableManager.hh"
#include "G4EmSaturation.hh"

//by Shuya 160316.
#include "G4ios.hh"
#include <iostream>
using namespace std;
G4int fVerboseLebel;

//by Shuya 160404
#include "G4hMultipleScattering.hh" // For Hardons/Ions only
#include "G4eMultipleScattering.hh"
#include "G4hIonisation.hh"

// Proton only processes


// Neutron only processes
/*
#include "G4NeutronInelasticProcess.hh"

#include "G4ExcitationHandler.hh"
#include "G4PreCompoundModel.hh"
#include "G4CascadeInterface.hh"
#include "G4BinaryCascade.hh"

#include "G4NeutronInelasticCrossSection.hh"   // Data Sets
#include "G4NeutronHPJENDLHEInelasticData.hh"
#include "G4NeutronHPJENDLHEElasticData.hh"
#include "G4HadronElasticDataSet.hh"
#include "G4NeutronInelasticCrossSection.hh"

#include "G4NeutronHPElastic.hh"
#include "G4NeutronHPElasticData.hh"

#include "G4NeutronHPInelastic.hh"
#include "G4NeutronHPInelasticData.hh"

// Other Neutron Elastic Scattering Models - 20/06/2007

#include "G4QElastic.hh"
//#include "G4UHadronElasticProcess.hh"

// A hadronic physics list.... see example/extended/hadronic/Hadr01
#include "HadronPhysicsQGSP_BIC.hh"
*/

#include "menate_R.hh"

// Most important part of program is "ConstructProcess()"
// Add EM physics processes in Construct EM.
// Refer to ExN02PhysicsList to see how to add EM physics in this file.
//
// 12/13/06 - Have only added interaction processes for protons and neutrons.
//

//by Shuya 160404
/*
TntPhysicsList::TntPhysicsList()
{;}

TntPhysicsList::~TntPhysicsList()
{;}
*/

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntNuclearPhysics::TntNuclearPhysics(const G4String& name)
                     :  G4VPhysicsConstructor(name) {}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

TntNuclearPhysics::~TntNuclearPhysics() {
  //fDecayProcess = NULL;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void TntNuclearPhysics::ConstructParticle()
{

  //baryons - Construct all baryons
  G4BaryonConstructor pBaryons;
  pBaryons.ConstructParticle();

  //Construct Generic Ions
  
  // Construct all light nuclei!

  G4Deuteron::DeuteronDefinition();
  G4Triton::TritonDefinition();
  G4Alpha::AlphaDefinition();
  G4He3::He3Definition();

 G4IonConstructor pIons;
 pIons.ConstructParticle();
 
//by Shuya 160404
/*
  //mesons - construct all mesons (strong interactions)
  G4MesonConstructor pMesons;
  pMesons.ConstructParticle();
*/
 
}

void TntNuclearPhysics::ConstructProcess()
{
  // Define physics processes

  ConstructNuclear();
}



// Add Nuclear Processes
// See "TntPhysicsList.hh" for processes included 

void TntNuclearPhysics::ConstructNuclear()
{
  // Proton nuclear scattering processes

  G4ProcessManager *pManager = 0;
  
   pManager = G4Proton::Proton()->GetProcessManager();

// Add proton EM processes
   pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
   pManager->AddProcess(new G4hIonisation,       -1,2,2);
	

// Create Pointers for Proton nuclear processes
/*
   G4HadronElasticProcess* theLE_ElasticProcess = new G4HadronElasticProcess();
   G4HadronFissionProcess* theLE_FissionProcess = new G4HadronFissionProcess();
   G4ProtonInelasticProcess* theLE_Proton_InelasticProcess = new G4ProtonInelasticProcess();


   G4LElastic* theLE_Elastic = new G4LElastic();
   G4LFission* theLE_Fission = new G4LFission();
   G4LEProtonInelastic* theLEProton_Inelastic = new G4LEProtonInelastic();

 // Add proton processes
       	theLE_ElasticProcess->RegisterMe(theLE_Elastic);
	theLE_FissionProcess->RegisterMe(theLE_Fission);
	theLE_Proton_InelasticProcess->RegisterMe(theLEProton_Inelastic);

	pManager->AddDiscreteProcess(theLE_ElasticProcess);
	pManager->AddDiscreteProcess(theLE_Proton_InelasticProcess);
	pManager->AddProcess(theLE_FissionProcess);
*/
//			  
// Neutron nuclear scattering processes
//

 pManager = G4Neutron::Neutron()->GetProcessManager();
	
// Create Pointers for Neutron nuclear processes
// Add neutron processes - Register Model then AddDiscreteProcess to use.
	
 //G4HadronElasticProcess* N_LE_ElasticProcess = new G4HadronElasticProcess();
 //G4HadronCaptureProcess* N_LE_CaptureProcess = new G4HadronCaptureProcess();
 //G4HadronFissionProcess* N_LE_FissionProcess = new G4HadronFissionProcess();
 //G4NeutronInelasticProcess* Neutron_InelasticProcess = new G4NeutronInelasticProcess();
 
 /* 
 G4LCapture* N_LE_Capture = new G4LCapture();
 G4LFission* N_LE_Fission = new G4LFission();
 N_LE_CaptureProcess->RegisterMe(N_LE_Capture);
 N_LE_FissionProcess->RegisterMe(N_LE_Fission);
 

// Neutron Nuclear Scattering Models - G4 - Ghiesha 

 G4LElastic* N_LE_Elastic = new G4LElastic();    
 N_LE_ElasticProcess->RegisterMe(N_LE_Elastic);
 
 G4LENeutronInelastic* N_LE_Neutron_Inelastic = new G4LENeutronInelastic();
 Neutron_InelasticProcess->RegisterMe(N_LE_Neutron_Inelastic);
 */
   // G4JENDL Data Sets
  
   //G4NeutronHPJENDLHEElasticData* theElasticData = new G4NeutronHPJENDLHEElasticData(); 
   //N_LE_ElasticProcess->AddDataSet(theElasticData);
   
 //  G4NeutronHPJENDLHEInelasticData* theInelasticData = new G4NeutronHPJENDLHEInelasticData();
 //  Neutron_InelasticProcess->AddDataSet(theInelasticData);
  
   // NeutronHP Physics models
   /*     
   G4NeutronHPElastic* theHPElastic = new G4NeutronHPElastic;
   G4NeutronHPElasticData* theHPElasticData = new G4NeutronHPElasticData;
   N_LE_ElasticProcess->RegisterMe(theHPElastic);
   N_LE_ElasticProcess->AddDataSet(theHPElasticData);
   
   G4NeutronHPInelastic* theHP_N_Inelastic = new G4NeutronHPInelastic;
   G4NeutronHPInelasticData* theHP_N_InelasticData = new G4NeutronHPInelasticData;
   Neutron_InelasticProcess->RegisterMe(theHP_N_Inelastic);
   Neutron_InelasticProcess->AddDataSet(theHP_N_InelasticData);
   */
   
   // QElastic and Precompound Physics models (above 10 MeV)
   // To use QElastic, should comment out G4HadronElasticProcess line above.
 /*
   const G4String& processName ="G4QElastic";
   G4QElastic* N_LE_ElasticProcess = new G4QElastic(processName);
 
   
   G4ExcitationHandler *theHandler = new G4ExcitationHandler();
   G4PreCompoundModel *thePreCompModel = new G4PreCompoundModel(theHandler);

   Neutron_InelasticProcess->RegisterMe(thePreCompModel);
 */

   // Add discrete process to the neutron process manager (otherwise process not used)   
 
   //pManager->AddDiscreteProcess(N_LE_ElasticProcess);
   //pManager->AddDiscreteProcess(Neutron_InelasticProcess);
   //pManager->AddDiscreteProcess(N_LE_CaptureProcess);
   //pManager->AddDiscreteProcess(N_LE_FissionProcess);

  
 G4String theProcessName = "menate_neutron";
 menate_R* theMENATE = new menate_R(theProcessName);
 pManager->AddDiscreteProcess(theMENATE);

				
// ////// Nuclei ///////////////////////////////////////////////// 
	
// Add Deuteron Physics (EM Stopping Power only)
	pManager = G4Deuteron::Deuteron()->GetProcessManager();
        
        pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
        pManager->AddProcess(new G4ionIonisation,       -1,2,2);

// Add Triton Physics (EM Stopping Power only)
	pManager = G4Triton::Triton()->GetProcessManager();
        
        pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
        pManager->AddProcess(new G4ionIonisation,       -1,2,2);

// Add Helium 3 Physics (EM Stopping Power only)
	pManager = G4He3::He3()->GetProcessManager();
        
        pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
        pManager->AddProcess(new G4ionIonisation,       -1,2,2);
	
		
// Add Alpha Physics (EM Stopping Power only)
	pManager = G4Alpha::Alpha()->GetProcessManager();

        pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
        pManager->AddProcess(new G4ionIonisation,       -1,2,2);
       	 
 // Add Generic Ion Physics (EM and Nuclear Elastic, Stopping Power, etc.)
     pManager = G4GenericIon::GenericIon()->GetProcessManager();
	
 // Nuclear Scattering Process for Generic Ion (see above)
     /*	 
 G4HadronElasticProcess* Ion_LE_ElasticProcess = new G4HadronElasticProcess();
 G4LElastic* Ion_LE_Elastic = new G4LElastic();
 Ion_LE_ElasticProcess->RegisterMe(Ion_LE_Elastic);
 pManager->AddDiscreteProcess(Ion_LE_ElasticProcess);
     */	
 // Add EM Processes to Generic Ions
     	
     pManager->AddProcess(new G4hMultipleScattering,-1,1,1);
     pManager->AddProcess(new G4ionIonisation,       -1,2,2);
     
}


/*
void TntPhysicsList::SetCuts()
{
  // suppress error messages even in case e/gamma/proton do not exist            
  G4int temp = GetVerboseLevel();                                                
  SetVerboseLevel(0);
  //
  //  Sets the cut on the physics interaction calculations                                                       
  //  " G4VUserPhysicsList::SetCutsWithDefault" method sets 
  //   the default cut value for all particle types 
  SetCutsWithDefault(); 

  G4double em_cuts = 1.0*mm;
  SetCutValue(em_cuts,"gamma");
  SetCutValue(em_cuts,"e-");
  SetCutValue(em_cuts,"e+");

  //cuts for ions set to default (1.0mm)
  // G4double ion_cuts = 1.0*mm;
  // SetCutValue(ion_cuts,"proton");

  // Retrieve verbose level
  SetVerboseLevel(temp);  
}

//by Shuya 160316
void TntPhysicsList::SetVerbose(G4int verbose)
{
	fVerboseLebel = verbose;
}
*/
