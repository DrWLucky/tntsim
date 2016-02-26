//==========================================================================
// TntSD.hh
// Based on DemonScintSD.hh and ExN04TrackerSD.hh
//
// Written/Modified by: Brian Roeder, LPC Caen 02/14/07
//                      email - roeder@lpccaen.in2p3.fr
//
// Purpose: Defines TntSD class for simulation data readout
//          Sets and accesses Data in TntDetHit.hh methods
//
//==========================================================================
//
// - See UM Hits Presentation and JLab Scoring 2 Talk for more info.
//

#ifndef TntSD_h
#define TntSD_h 1

#include "G4VSensitiveDetector.hh"
#include "TntDetHit.hh"
#include "TntDataRecordTree.hh"

#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

class G4Step;

class G4HCofThisEvent;

class G4TouchableHistory;

class TntSD : public G4VSensitiveDetector
{
public:

  TntSD(G4String DetName, G4String Light);
  ~TntSD();
  void Initialize(G4HCofThisEvent *HCE);
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory*);  // Will add TH Later
  void EndOfEvent(G4HCofThisEvent* );
  
  G4double ConvertToLight(G4String theName, G4double theCharge, G4double edep, G4String Light);

private:

  TntDetHitsCollection* TntHitsCollection;
  TntDataRecordTree* TntDataOutEV;

  G4String Light_Conv;
 
       
};

#endif
