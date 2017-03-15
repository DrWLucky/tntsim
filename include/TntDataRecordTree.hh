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
#include <iomanip>
#include <fstream>
#include <vector>

// Root Analysis header files for C++

#include "TROOT.h"
#include "Riostream.h"
#include "TFile.h"
#include "TH1.h"
#include "TNtuple.h"
#include "TTree.h"
#include "globals.hh"

#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"

//by Shuya 160407
#include "TntDataRecordTree.hh"
#include "TH2.h"


class TntDataRecordTree
{
public:
	struct Hit_t { 
		G4double X, Y, Z, T, E;
		G4int TrackID;
		bool operator== (const Hit_t& rhs) {
			if(rhs.X == X && rhs.Y == Y && rhs.Z == Z && 
				 rhs.T == T && rhs.E == E && rhs.TrackID == TrackID) 
			{ 	return true;   }
			else { return false; }
		}
		bool operator!= (const Hit_t& rhs) {
			return !(this->operator==(rhs));
		}
};
	
private:

  // Initialized in class constructor in TntDataRecordTree.cc 
  
  TFile* DataFile;
  TTree* TntEventTree;
//by Shuya 160422.
  TTree* TntEventTree2;
  G4double eng_int;
  G4double eng_Tnt;
  G4double eng_Tnt_proton;
  G4double eng_Tnt_alpha;
  G4double eng_Tnt_C12;
  G4double eng_Tnt_EG;
  G4double eng_Tnt_Exotic;
//by Shuya 160407
  G4int eng_Tnt_PhotonFront;
  G4int eng_Tnt_PhotonBack;
//by Shuya 160502
  G4int eng_Tnt_PhotonTotal;
//by Shuya 160502
  G4double edep_Tnt;
  G4double edep_Tnt_proton;
  G4double edep_Tnt_alpha;
  G4double edep_Tnt_C12;
  G4double edep_Tnt_EG;
  G4double edep_Tnt_Exotic;
//by Shuya 160504
  G4double num_Tnt_NonPMT;
  G4double num_Tnt_Abs;
//by Shuya 160422
  G4int PmtFrontHit[4][4];
  G4int PmtBackHit[4][4];
  // std::vector<std::vector<int> > PmtFront;
  // std::vector<std::vector<int> > PmtBack;
//by Shuya 160509. This is needed since array size is not passed until you get extern G4int in TntDataRecord(), which is originally fixed in Tnt.cc file.
  //G4int** PmtFrontHit;
  //G4int** PmtBackHit;

	std::vector<G4int> PhotonSum;
	std::vector<G4double> HitX;
	std::vector<G4double> HitY;
	std::vector<G4double> HitZ;
	std::vector<G4double> HitT;
	std::vector<G4double> HitE;
	std::vector<G4int>    HitTrackID;
	G4int NumHits;
	

  G4double FirstHitTime;
  G4double FirstHitMag;
//by Shuya 160502
  G4ThreeVector FirstHitPosition;

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
//by Shuya 160407
  int number_Photon;

  // Efficiency Calculators
  int number_at_this_energy;
  double efficiency;


  
 public:
    TntDataRecordTree(G4double Threshold);
   ~TntDataRecordTree();

   TH1F *h_Energy_Initial;  
   TH1D *h_Energy_Tnt;  
   TH1F *h_Energy_Proton;  
   TH1F *h_Energy_Alpha;  
   TH1F *h_Energy_C12;  
   TH1F *h_Energy_EG;  
   TH1F *h_Energy_Exotic;  
   TH1F *h_Energy_Photon;  
   //TH2I *h_Count_PMT_Front; 
   //TH2I *h_Count_PMT_Back; 
	std::vector<TH2I*> h_Count_PMT_Front; 
	std::vector<TH2I*> h_Count_PMT_Back; 


 // Creates Pointer to Data Analysis in main() program (Tnt.cc).    
  static TntDataRecordTree* TntPointer;
 
  // Definition of member functions - access and distribute data
  // See "TntDataRecordTree.cc for the function listings.

//by Shuya 160408
  //void senddataPMT(int id, double value1);
//by Shuya 160421
  void senddataPMT(int id, int value1, int evid);
  void createdataPMT(int evid);

  void senddataPG(double value1);
  void senddataEV(int type, double value1);
  void senddataPosition(G4ThreeVector pos);
	void senddataHits(const std::vector<Hit_t>& hit);
  void senddataTOF(G4double time);
  void ShowDataFromEvent();
  void FillTree();
//by Shuya 160422.
  void FillTree2(int evid);
  void GetParticleTotals();
  void CalculateEff(int ch_eng);

private:
  TntDataRecordTree() {;}   // Hide Default Constructor
}; 
#endif
