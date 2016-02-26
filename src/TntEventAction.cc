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
// $Id: TntEventAction.cc,v 1.2 2007/02/14 17:45:14 roeder Exp $
// GEANT4 tag $Name: geant4-08-02 $
//
// Modified by Brian Roeder, LPC Caen, monitor of events for Tnt
// email - roeder@lpccaen.in2p3.fr
//
// 2/15/07 - Now reads events generated in TntSD from TntDetHit
 
#include "TntEventAction.hh"

#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4RunManager.hh"
#include "G4SDManager.hh"

#include "G4TrajectoryContainer.hh"
#include "G4Trajectory.hh"

// Added following classes for Hit Map Printing

#include "TntDetHit.hh"
#include "G4HCofThisEvent.hh"
#include "G4THitsMap.hh"
#include "G4UnitsTable.hh"
#include "G4ThreeVector.hh"

// These classes are for scorers (being phased out)
/*
#include "G4MultiFunctionalDetector.hh"
#include "G4VPrimitiveScorer.hh"
*/

// 19 Feb 07 - BTR - Root Files now made in TntDataRecordTree 
//                   Data Analysis/Readout done in TntSD.cc - EndOfEvent()
//

TntEventAction::TntEventAction() :
  numberOfEvent(-1)
{
  /* Constructor - Initializes Counters and Analysis Pointer */
  TntGetDataEV = TntDataRecordTree::TntPointer;
}

//
TntEventAction::~TntEventAction()
{
/* Destructor */
}

//
void TntEventAction::BeginOfEventAction(const G4Event*)
{}

//
void TntEventAction::EndOfEventAction(const G4Event* evt)
{
  
  // Periodically Print out Event Results!
  // Analysis of Event done in TntSD.cc EndOfEvent() Method.

 numberOfEvent++;
 
      //  Periodic Printout of event  //

 G4RunManager* runManager = G4RunManager::GetRunManager();
 G4SDManager* SDMan = G4SDManager::GetSDMpointer();

 G4String detName = "TntSD";
 G4int collectionID = SDMan->GetCollectionID("TntSD/TntHitsCollection");
 const G4Event* currentEvent = runManager->GetCurrentEvent();

 G4HCofThisEvent* HCE = currentEvent->GetHCofThisEvent();
 TntDetHitsCollection* myCollection = (TntDetHitsCollection*)(HCE->GetHC(collectionID));

       G4int event_show = 10000;
      // G4int event_show = 1;   
      if((numberOfEvent) % event_show == 0) 
      {
      G4cout << "====================================================" << G4endl;
      G4cout << ">>> Event " << evt->GetEventID() << G4endl;

      G4int n_hits = myCollection->entries();
      G4cout << n_hits << " Hits were registered in this event! " << G4endl; 
      TntGetDataEV->ShowDataFromEvent();
      //  G4cout << "The Event Action Tnt Data Out pointer is: " << TntDataOutEV << G4endl;
      }
     
}  // End of "End of Event Action"
