//-----------------------------------------------------------------
// TntDataRecordTree.hh
//
// - written by: Brian Roeder, LPC Caen, 18 Jan 07 (in Serra)
// - email - roeder@lpccaen.in2p3.fr
//
// - modified: 14/02/07 for use with Tnt Tracking
//
// - Usage: A C++ class for GEANT4 for creating a ROOT tree for 
// -        the Tnt neutron detector simulation event by event.
//
///////////////////////////////////////////////////////////////////
//
//

#ifndef DATARECORD_H
#define DATARECORD_H

#include <iostream>
using namespace std;

// C++ formatting and file io headers

#include <iomanip>
#include <fstream>

// Root Analysis header files for C++

#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include "TTree.h"
#include "globals.hh"

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

class TntDataRecordTree
{

private:

  // Initialized in class constructor in TntDataRecordTree.cc 
  
  TFile* DataFile;
  TTree* TntEventTree;
  G4double eng_int;
  G4double eng_Tnt;
  G4double eng_Tnt_proton;
  G4double eng_Tnt_alpha;
  G4double eng_Tnt_C12;
  G4double eng_Tnt_EG;
  G4double eng_Tnt_Exotic;

  G4double FirstHitTime;
  G4double FirstHitMag;

  G4double Xpos;
  G4double Ypos;
  G4double Zpos;

  G4double Det_Threshold; // Threshold for Detector in MeVee

  // Particle Counters

  int event_counter;
  int number_total;
  int number_protons;
  int number_alphas;
  int number_C12;
  int number_EG;
  int number_Exotic;

  // Efficiency Calculators
  int number_at_this_energy;
  double efficiency;
  
 public:
    TntDataRecordTree(G4double Threshold);
   ~TntDataRecordTree();
  
 // Creates Pointer to Data Analysis in main() program (Tnt.cc).    
  static TntDataRecordTree* TntPointer;
 
  // Definition of member functions - access and distribute data
  // See "TntDataRecordTree.cc for the function listings.

  void senddataPG(double value1);
  void senddataEV(int type, double value1);
  void senddataPosition(G4ThreeVector pos);
  void senddataTOF(G4double time);
  void ShowDataFromEvent();
  void FillTree();
  void GetParticleTotals();
  void CalculateEff(int ch_eng);

private:
  TntDataRecordTree() {;}   // Hide Default Constructor
}; 
#endif
