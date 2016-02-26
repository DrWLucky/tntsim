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
// $Id: TntPhysicsList.hh,v 1.6 2006/12/12 17:47:15 roeder Exp $
// GEANT4 tag $Name: geant4-08-01-patch-02 $
//
// 
//
// TntPhysicsList
//  Construct/define particles and physics processes
//
//  Particles defined in Tnt
//    geantino, protons(em)
//
//  Process defined in Tnt
//    transportation 
//    12/12/2006 -  added EM processes (protons, electrons, positrons, gammas)
//    12/12/2006 -  added Hadronic processes for protons, neutrons
//               -  followed example in Demon
//
//    02/14/2007 - Converted from Serra Simulation to Tnt.
//               - Currently uses GEANT4 Neutron Elastic Parameterisation
//               - and Inelastic Process with unedited JENDL data added. 
//

#ifndef TntPhysicsList_h
#define TntPhysicsList_h 1

#include "G4VUserPhysicsList.hh"
#include "globals.hh"

// Included Physics Lists

// EM processes
// Hadrons
#include "G4eMultipleScattering.hh"
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4ionIonisation.hh"

//Low Energy EM leptons and Boson Processes (following ExNO2)
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4eplusAnnihilation.hh"

 
// Strong force processes ("Low" Energy only)
#include "G4ProtonInelasticProcess.hh"
#include "G4NeutronInelasticProcess.hh"
#include "G4DeuteronInelasticProcess.hh"
#include "G4TritonInelasticProcess.hh"
#include "G4AlphaInelasticProcess.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


class TntPhysicsList: public G4VUserPhysicsList
{
  public:
    TntPhysicsList();
    ~TntPhysicsList();

  protected:
    // Construct particle and physics process
    // 11/16/2006 - Constructed Geantinos and Protons
    void ConstructParticle();
    void ConstructProcess();
    void SetCuts();

  protected:
  // these methods Construct physics processes and register them
  // construct proton em scattering and energy-loss processes
  // construct Low Energy Scattering processes for Neutrons
    void ConstructEM();
    void ConstructNuclear();

  protected:
  //EM process pointers
    // Gamma physics
    G4PhotoElectricEffect thePhotoEffect;
    G4ComptonScattering theComptonEffect;
    G4GammaConversion thePairProduction;
  
    // Electron physics
    G4eMultipleScattering theElectronMultipleScattering;
    G4eIonisation theElectronIonisation;
    G4eBremsstrahlung theElectronBremsStrahlung;
  
    //Positron physics
    G4eMultipleScattering thePositronMultipleScattering;
    G4eIonisation thePositronIonisation; 
    G4eBremsstrahlung thePositronBremsStrahlung;  
    G4eplusAnnihilation theAnnihilation;


};

#endif







