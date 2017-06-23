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
#include "TLorentzVector.h"
#include "TClonesArray.h"
#include "TLorentzVector.h"

#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"
#include "G4LorentzVector.hh"

//by Shuya 160407
#include "TntDataRecordTree.hh"
#include "TH2I.h"
class TVector3;


class TntDataRecordTree
{
public:
	struct Hit_t { 
		G4double X, Y, Z, T, E;
		G4int TrackID, ParentTrackID, Type;
		bool operator== (const Hit_t& rhs) {
			if(rhs.X == X && rhs.Y == Y && rhs.Z == Z && 
				 rhs.T == T && rhs.E == E && 
				 rhs.TrackID == TrackID && rhs.ParentTrackID &&
				 rhs.Type == Type) 
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
	TTree* TntInputTree;
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
  G4int PmtFrontHit[64][64];
  G4int PmtBackHit[64][64];
  // std::vector<std::vector<int> > PmtFront;
  // std::vector<std::vector<int> > PmtBack;
//by Shuya 160509. This is needed since array size is not passed until you get extern G4int in TntDataRecord(), which is originally fixed in Tnt.cc file.
  //G4int** PmtFrontHit;
  //G4int** PmtBackHit;
	///
	/// GAC
	std::vector<G4int> PhotonSum; // Sum of all photons incident on a PMT
	std::vector<G4int> PhotonSumFront;
	TH2I* hDigi; // Histogram of digitized time signals for each PMT
	///
	/// Vectors of all hit information
	G4int NumHits;
	std::vector<G4double> HitX;
	std::vector<G4double> HitY;
	std::vector<G4double> HitZ;
	std::vector<G4double> HitT;
	std::vector<G4double> HitE;
	std::vector<G4int>    HitTrackID;
	std::vector<G4int>    HitType;
	TClonesArray* fHits;
	TClonesArray* fHit01;
	Int_t iHit0, iHit1;

	///
	/// MENATA_R hits
	TClonesArray* fMenateHitsPos;
	std::vector<G4double> fMenateHitsE;	
	std::vector<G4int> fMenateHitsType;	
	
	///
	/// Positions of original fired neutron
	G4double PrimaryX;
	G4double PrimaryY;
	G4double PrimaryZ;
	TLorentzVector* PrimaryMomentum;
	TLorentzVector* SecondaryMomentum; // recoil momentum if neutron decay
	TLorentzVector* EjectileMomentum;  // ejectile from population reaction [e.g. (d,3He)]
	TLorentzVector* BeamMomentum;      // beam from population reaction
	G4double mTrgt;                    // mass of target in population reaction
	G4double ReacThetaCM; // COM angle of reaction
	
	TVector3* SecondaryPosition; // other particles involved in reaction
	TVector3* EjectilePosition;  // other particles involved in reaction
	TVector3* BeamPosition;      //
	

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

	// INPUT PARAMETERS //
	G4int npmtX, npmtY;
	G4double eNeut;
	G4double detector_x, detector_y, detector_z;

  
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
	void senddataPMT_Time(int id, G4double time);
  void senddataPMT(int id, int value1, int evid);
  void createdataPMT(int evid);

  void senddataPG(double value1);
	void senddataPrimary(const G4ThreeVector& posn, const G4ThreeVector& momentum);
	void senddataSecondary(const G4ThreeVector& posn, const G4LorentzVector& momentum);
	void senddataEjectile(const G4ThreeVector& posn, const G4LorentzVector& momentum,
												const G4double& ThetaCM);
	void senddataReaction(const G4ThreeVector& beamPosn, const G4LorentzVector& beamMomentum,
												const G4double& targetMass);
  void senddataEV(int type, double value1);
  void senddataPosition(const G4ThreeVector& pos);
	void senddataHits(const std::vector<Hit_t>& hit, bool sortTime);
  void senddataTOF(G4double time);
	void senddataMenateR(G4double ekin, const G4ThreeVector& posn, G4double t, G4int type);
  void ShowDataFromEvent();
  void FillTree();
//by Shuya 160422.
  void FillTree2(int evid);
  void GetParticleTotals();
  void CalculateEff(int ch_eng);

	G4int GetParticleCode(const G4String& name);
	G4int GetReactionCode(const G4String& name);
	G4String GetParticleName(G4int  code);
	G4String GetReactionName(G4int  code);	
	
private:
  TntDataRecordTree() {;}   // Hide Default Constructor
}; 
#endif
