//----------------------------------------------------------------
// TntDataRecordTree.cc
//
// - written by: Brian Roeder, LPC Caen, 18 Jan 07
// - email - roeder@lpccaen.in2p3.fr
//
// - Usage: A C++ class for GEANT4 for creating a ROOT tree for 
// -        the Tnt neutron detector simulation event by event.
//
///////////////////////////////////////////////////////////////////
//
// 19 Jan 07 - For methods and functions associated with TntDataRecordTree
//           - pointers, see TntDataRecordTree.hh 
//
// 25 Apr 08 - Detector Threshold is now set in main() of Tnt_EMIS.cc
// 
// 
//
#include <cassert>
#include <algorithm>
#include "TntGlobalParams.hh"
#include "TntDataRecordTree.hh"

#include <TH2I.h>
#include <TVector3.h>
using namespace std;

//by Shuya 160426.
#ifdef __MAKECINT__
#pragma link C++ class vector<vector<double> >+;
#pragma link C++ class vector<int>+;
#endif

//by Shuya 160509
extern G4int NX;
extern G4int NY;

// Access to Analysis pointer! (see TntSD.cc EndOfEvent() for Example)
TntDataRecordTree* TntDataRecordTree::TntPointer;

TntDataRecordTree::TntDataRecordTree(G4double Threshold) : 
  // Initialized Values
  eng_int(0), eng_Tnt(0), eng_Tnt_alpha(0), eng_Tnt_C12(0), 
  eng_Tnt_EG(0), eng_Tnt_Exotic(0),  FirstHitTime(0), FirstHitMag(0),
//by Shuya 160407
  eng_Tnt_PhotonFront(0), eng_Tnt_PhotonBack(0), eng_Tnt_PhotonTotal(0), number_Photon(0),
  Xpos(0), Ypos(0), Zpos(0), Det_Threshold(Threshold),
  event_counter(0), number_total(0), 
  number_protons(0), number_alphas(0), number_C12(0), number_EG(0), 
  number_Exotic(0), number_at_this_energy(0), efficiency(0),
//by Shuya 160502
  eng_Tnt_proton(0), edep_Tnt(0), edep_Tnt_proton(0), edep_Tnt_alpha(0), edep_Tnt_C12(0), edep_Tnt_EG(0), edep_Tnt_Exotic(0),
//by Shuya 160504
  num_Tnt_NonPMT(0), num_Tnt_Abs(0)
{ /* Constructor */ 
  TntPointer = this;  // When Pointer is constructed, assigns address of this class to it.
  //
  // Create new data storage text file
  // Create new text file for data storage - (Data Recorded by TntDataRecordTree class)
 
  char EffFile[] = "eff_results_file.dat";
  ofstream outfile2(EffFile,ios::trunc);
  outfile2.close();

//by Shuya 160422.
/*
	PmtFrontHit.resize(NX);
	for(int i=0;i<NX;i++)	PmtFrontHit[i].resize(NY);
	PmtBackHit.resize(NX);
	for(int i=0;i<NX;i++)	PmtBackHit[i].resize(NY);
*/


//by Shuya 160418? Initialization just in case...
//by Shuya 160509
//Comment by Shuya 160510. I tried to make a 2D variable vector as below, but this didn't give right values when you filled in TBranch.. (don't know why)
//So I have to set the size of array in TntDataRecord.hh so far...
  //PmtFrontHit = new G4int* [NX];
  for(int i=0;i<NX;i++)
  {
  	//PmtFrontHit[i] = new G4int [NY];
		for(int j=0;j<NY;j++)
		{
			PmtFrontHit[i][j] = 0;
			//G4cout << "TEST PMTFRONTHIT ADDRESS " << &PmtFrontHit[i][j] << G4endl;
		}
  }

  //PmtBackHit = new G4int* [NX];
  for(int i=0;i<NX;i++)
  {
  	//PmtBackHit[i] = new G4int [NY];
		for(int j=0;j<NY;j++)
		{
			PmtBackHit[i][j] = 0;
			//G4cout << "TEST PMTFRONTHIT ADDRESS " << &DAMMY[i][j] << G4endl;
		}
  }

  //
  // Create pointers to ROOT Analysis Tree and Branches
  //
  cout <<"\n Starting Data Tree Constructor" << endl;
  
//  const Char_t* evt_file = "TntDataTree.root";
	
  DataFile = new TFile(TntGlobalParams::Instance()->GetRootFileName().c_str(), "RECREATE");
//G4cout << DataFile << "!!" << G4endl;

//////////////////////////// BY Shuya 160407 TO CHANGE TREE TO HISTOGRAMS IN ROOT ////////////////////////////////////
//TFile *foutR = new TFile(rootFileName,"RECREATE");
	if ( DataFile->IsOpen() )
	{
		printf("root file opened for writing.\n");
	}
	else
	{
		printf("root file didn't open.\n");
		exit (1);
	}

	printf("Output the data files\n");

//By Shuya 160426. I removed all these histograms to replace TTree. 
/*
	TH1F *h_Energy_Initial = new TH1F("Energy_Initial","Energy_Initial", 500, 0, 5000);
	h_Energy_Initial->GetXaxis()->SetTitle("Energy (keV)");
	h_Energy_Initial->GetYaxis()->SetTitle("Counts / bin");

	TH1D *h_Energy_Tnt = new TH1D("Energy_Tnt","Energy_Tnt", 500, 0, 1.);
	h_Energy_Tnt->GetXaxis()->SetTitle("Energy (keV)");
	h_Energy_Tnt->GetYaxis()->SetTitle("Counts / bin");

	TH1F *h_Energy_Proton = new TH1F("Energy_Proton","Energy_Proton", 500, 0, 1.);
	h_Energy_Proton->GetXaxis()->SetTitle("Energy (keV)");
	h_Energy_Proton->GetYaxis()->SetTitle("Counts / bin");
	
	TH1F *h_Energy_Alpha = new TH1F("Energy_Alpha","Energy_Alpha", 500, 0, 1.);
	h_Energy_Alpha->GetXaxis()->SetTitle("Energy (keV)");
	h_Energy_Alpha->GetYaxis()->SetTitle("Counts / bin");

	TH1F *h_Energy_C12 = new TH1F("Energy_C12","Energy_C12", 500, 0, 1.);
	h_Energy_C12->GetXaxis()->SetTitle("Energy (keV)");
	h_Energy_C12->GetYaxis()->SetTitle("Counts / bin");

	TH1F *h_Energy_EG = new TH1F("Energy_EG","Energy_EG", 500, 0, 1.);
	h_Energy_EG->GetXaxis()->SetTitle("Energy (keV)");
	h_Energy_EG->GetYaxis()->SetTitle("Counts / bin");

	TH1F *h_Energy_Exotic = new TH1F("Energy_Exotic","Energy_Exotic", 500, 0, 1.);
	h_Energy_Exotic->GetXaxis()->SetTitle("Energy (keV)");
	h_Energy_Exotic->GetYaxis()->SetTitle("Counts / bin");

	TH1F *h_Energy_Photon = new TH1F("Energy_Photon","Energy_Photon", 500, 0, 1000);
	h_Energy_Photon->GetXaxis()->SetTitle("Energy (keV)");
	h_Energy_Photon->GetYaxis()->SetTitle("Counts / bin");
*/

/*
	TH1F *h_First_Hit_Pos = new TH1F("First_Hit_Pos","First_Hit_Pos", 500, 0, 5000);
	h_First_Hit_Pos->GetXaxis()->SetTitle("Pos (keV)");
	h_First_Hit_Pos->GetYaxis()->SetTitle("Counts / bin");

	TH1F *h_First_Hit_Time = new TH1F("First_Hit_Time","First_Hit_Time", 500, 0, 5000);
	h_First_Hit_Time->GetXaxis()->SetTitle("Pos (keV)");
	h_First_Hit_Time->GetYaxis()->SetTitle("Counts / bin");

	TH1F *h_Xpos = new TH1F("Xpos","Xpos", 500, 0, 5000);
	h_Xpos->GetXaxis()->SetTitle("Pos (keV)");
	h_Xpos->GetYaxis()->SetTitle("Counts / bin");

	TH1F *h_Ypos = new TH1F("Ypos","Ypos", 500, 0, 5000);
	h_Ypos->GetXaxis()->SetTitle("Pos (keV)");
	h_Ypos->GetYaxis()->SetTitle("Counts / bin");

	TH1F *h_Zpos = new TH1F("Zpos","Zpos", 500, 0, 5000);
	h_Zpos->GetXaxis()->SetTitle("Pos (keV)");
	h_Zpos->GetYaxis()->SetTitle("Counts / bin");
*/
//by Shuya 160421
/*
	TH2I *h_Count_PMT_Front[100];
	for(int i=0;i<100;i++)
	{
	char brN[300];
	sprintf(brN, "Count_PMT_Front_run%d",i);
	h_Count_PMT_Front[i] = new TH2I(brN,brN, 10, 0, 10, 10, 0, 10);
	h_Count_PMT_Front[i]->GetXaxis()->SetTitle("X (no.)");
	h_Count_PMT_Front[i]->GetYaxis()->SetTitle("Y (no.)");
	}

	TH2I *h_Count_PMT_Back[100];
	for(int i=0;i<100;i++)
	{
	char brN[300];
	sprintf(brN, "Count_PMT_Back_run%d",i);
	h_Count_PMT_Back[i] = new TH2I(brN,brN, 10, 0, 10, 10, 0, 10);
	h_Count_PMT_Back[i]->GetXaxis()->SetTitle("X (no.)");
	h_Count_PMT_Back[i]->GetYaxis()->SetTitle("Y (no.)");
	}
*/

//by Shuya 160426. Removed the histogram to replace with TTree.
//vector<TH2I*> h_Count_PMT_Front;
	vector<TH2I*> h_Count_PMT_Back;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
  

  TntEventTree = new TTree("t","Tnt Scintillator Simulation Data");
  TntEventTree->Branch("Energy_Initial",&eng_int,"eng_int/D");
  TntEventTree->Branch("LightOutput_Tnt",&eng_Tnt,"eng_Tnt/D");
  TntEventTree->Branch("LightOutput_Proton",&eng_Tnt_proton,"eng_Tnt_proton/D");
  TntEventTree->Branch("LightOutput_Alpha",&eng_Tnt_alpha,"eng_Tnt_alpha/D");
  TntEventTree->Branch("LightOutput_C12", &eng_Tnt_C12,"eng_Tnt_C12/D");
  TntEventTree->Branch("LightOutput_EG",&eng_Tnt_EG,"eng_Tnt_EG/D");
  TntEventTree->Branch("LightOutput_Exotic",&eng_Tnt_Exotic,"eng_Tnt_Exotic/D");
//by Shuya 160502
  TntEventTree->Branch("EnergyDeposit_Total",&edep_Tnt,"edep_Tnt/D");
  TntEventTree->Branch("EnergyDeposit_Proton",&edep_Tnt_proton,"edep_Tnt_proton/D");
  TntEventTree->Branch("EnergyDeposit_Alpha",&edep_Tnt_alpha,"edep_Tnt_alpha/D");
  TntEventTree->Branch("EnergyDeposit_C12", &edep_Tnt_C12,"edep_Tnt_C12/D");
  TntEventTree->Branch("EnergyDeposit_EG",&edep_Tnt_EG,"edep_Tnt_EG/D");
  TntEventTree->Branch("EnergyDeposit_Exotic",&edep_Tnt_Exotic,"edep_Tnt_Exotic/D");
//by Shuya 160407
  TntEventTree->Branch("Num_DetectedPhotonFront",&eng_Tnt_PhotonFront,"eng_Tnt_PhotonFront/I");
//by Shuya 160427
  TntEventTree->Branch("Num_DetectedPhotonBack",&eng_Tnt_PhotonBack,"eng_Tnt_PhotonBack/I");
//by Shuya 160502
  TntEventTree->Branch("Num_CreatedPhotonTotal",&eng_Tnt_PhotonTotal,"eng_Tnt_PhotonTotal/I");
//by Shuya 160504
  TntEventTree->Branch("Num_NonPMTCountTotal",&num_Tnt_NonPMT,"num_Tnt_NonPMT/D");
  TntEventTree->Branch("Num_AbsorptionInDetectorTotal",&num_Tnt_Abs,"num_Tnt_Abs/D");

  TntEventTree->Branch("First_Hit_Pos",&FirstHitMag,"FirstHitMag/D");
  TntEventTree->Branch("First_Hit_Time",&FirstHitTime,"FirstHitTime/D");

  TntEventTree->Branch("Xpos",&Xpos,"Xpos/D");
  TntEventTree->Branch("Ypos",&Ypos,"Ypos/D");
  TntEventTree->Branch("Zpos",&Zpos,"Zpos/D");
	

	// GAC - Array of PMT intensities (photon counts)
	// 
	TntEventTree->Branch("PhotonSum", &PhotonSum);

	// Digitizer histogram (see createdataPMT for more info)
	hDigi = 0;
	TntEventTree->Branch("digi", "TH2I", &hDigi);
	//
	// Array of hit positions and times
	TntEventTree->Branch("HitX", &HitX);
	TntEventTree->Branch("HitY", &HitY);
	TntEventTree->Branch("HitZ", &HitZ);
	TntEventTree->Branch("HitT", &HitT);
	TntEventTree->Branch("HitE", &HitE);
	TntEventTree->Branch("HitTrackID", &HitTrackID);
	TntEventTree->Branch("NumHits", &NumHits);
	//
	//
	fHits = new TClonesArray("TLorentzVector");
	TntEventTree->Branch("HitArray", &fHits, 256000, 0); // splitlevel 0 for custom streamer
	fHits->BypassStreamer();
	fHit01 = new TClonesArray("TLorentzVector");
	TntEventTree->Branch("Hit01", &fHit01, 256000, 0);
	fHit01->BypassStreamer();
	//
	// Original x,y,z positions of the fired neutron
	PrimaryMomentum = 0;
	TntEventTree->Branch("PrimaryX",&PrimaryX,"PrimaryX/D");
  TntEventTree->Branch("PrimaryY",&PrimaryY,"PrimaryY/D");
  TntEventTree->Branch("PrimaryZ",&PrimaryZ,"PrimaryZ/D");
	TntEventTree->Branch("PrimaryMomentum", &PrimaryMomentum);
	//
	// Secondary particles involved in the reaction
	SecondaryMomentum = 0;
	SecondaryPosition = 0;
		

//by Shuya 160422. Making tree for photon hits on each pmt.
  TntEventTree2 = new TTree("t2","Tnt Scintillator Simulation Data");
//by Shuya 160426. NOTE THIS IS REQUIRED otherwise branches are not drawn in the tree. # of entries is 1 this case because each branch just filled one time with [10][10](i.e.,100) data.  
  TntEventTree2->SetEntries(1);

/*
  for(int i=0;i<10;i++)
  {
	for(int j=0;j<10;j++)
	{
	char brN[300];
	char brN2[300];
	sprintf(brN, "PMT_Front_Hit_%d_%d",i,j);
	sprintf(brN2, "PMT_Front_Hit_%d_%d/I",i,j);
	TntEventTree2->Branch(brN,&PmtFrontHit[i][j],brN2);
	sprintf(brN, "PMT_Back_Hit_%d_%d",i,j);
	sprintf(brN2, "PMT_Back_Hit_%d_%d/I",i,j);
	TntEventTree2->Branch(brN,&PmtBackHit[i][j],brN2);
	}
  }
*/

  cout << TntPointer << endl;
}

TntDataRecordTree::~TntDataRecordTree()
{/* Destructor, Close root file */

	hDigi->Delete();
	hDigi = 0;
  DataFile->Write(); 
  DataFile->Close();
  cout << "############################################################" << endl;
  cout << "Created Root Tree File = \"TntDataTree.root\"" << endl;
  cout << "Got the ROOT Tree File from Data Record Routine!" << endl; 
  delete DataFile;

/*
//by Shuya 160509
for(int i=0;i<NX;i++)
{
delete [] PmtFrontHit[i];
delete [] PmtBackHit[i];
}
delete [] PmtFrontHit;
delete [] PmtBackHit;
*/
}

// Member functions of TntDataRecordTree class

void TntDataRecordTree::senddataPG(double value1=0.)
{
	eng_int = value1;
	event_counter++;
	//  cout << "eng_int = " << eng_int << endl;
}

void TntDataRecordTree::senddataPrimary(const G4ThreeVector& pos, const G4ThreeVector& mom)
{
	PrimaryX = pos.x();
	PrimaryY = pos.y();
	PrimaryZ = pos.z();

	TVector3 v(mom.x(), mom.y(), mom.z());
	G4double theta = v.Theta(), phi = v.Phi();

	const G4double MNEUT = 939.565378;
	G4double etot = eng_int + MNEUT;
	G4double ptot = sqrt(etot*etot - MNEUT*MNEUT);
	PrimaryMomentum->SetPxPyPzE(ptot*sin(theta)*cos(phi), 
															ptot*sin(theta)*sin(phi),
															ptot*cos(theta), 
															etot);
}

void TntDataRecordTree::senddataSecondary(const G4ThreeVector& pos, const TLorentzVector& mom)
{
	if(SecondaryMomentum == 0) {
		TntEventTree->Branch("SecondaryMomentum", &SecondaryMomentum);
	}
	if(SecondaryPosition == 0) {
		TntEventTree->Branch("SecondaryPosition", &SecondaryPosition);
	}

	SecondaryPosition->SetXYZ(pos.x(), pos.y(), pos.z());
	*SecondaryMomentum = mom;
}

//by Shuya 160408
//void TntDataRecordTree::senddataPMT(int id, double value1)
//by Shuya 160421
void TntDataRecordTree::senddataPMT(int id, int value1, int evid)
{
	int x, y;
	char brN[300];
	PhotonSum.resize(NX*NY);
	
	//if(id < 100)	//Front side
	//by Shuya 160509
	if(id < (NX*NY))	//Front side
	{
		sprintf(brN, "Count_PMT_Front_Event%d",evid);
		x = id / NX;
		y = id % NY;
		//((TH2I*)DataFile->Get("Count_PMT_Front"))->Fill(x,y,value1);
		//for(int i = 0;i<value1;i++)	((TH2I*)DataFile->Get("Count_PMT_Front"))->Fill(x,y);
		//by Shuya 160426. I removed this to replace with TTree branch.
		//for(int i = 0;i<value1;i++)	((TH2I*)DataFile->Get(brN))->Fill(x,y);

		PmtFrontHit[x][y] = value1;
	}
	//else if(id >= 100 && id < 200)	//Back side
	else if(id >= (NX*NY) && id < (2*NX*NY))	//Back side
	{
		sprintf(brN, "Count_PMT_Back_Event%d",evid);
		//by Shuya 160509
		//x = (id-100) / 10;
		//y = (id-100) % 10;
		x = (id-NX*NY) / NX;
		y = (id-NX*NY) % NY;
		PmtBackHit[x][y] = value1;
		PhotonSum.at(id - NX*NY) = value1;
	}
}

void TntDataRecordTree::senddataPMT_Time(int id, G4double time)
{
	if(id >= (NX*NY) && id < (2*NX*NY)) {	//Back side
		hDigi->Fill(time, id - NX*NY);
	}	
}


//by Shuya 160421
void TntDataRecordTree::createdataPMT(int evid)
{
	char brN[300];

	sprintf(brN, "PMT_Front_Hit_Event_%d", evid);
	//G4cout << "FORMATTEST " << TString::Format("PmtFrontHit[%d][%d]/I", NX,NY) << G4endl;
	//sprintf(brN2, "PMT_Front_Hit_Event_%d/I", evid);
	//by Shuya 160509
	//TntEventTree2->Branch(brN,PmtFrontHit,"PmtFrontHit[8][8]/I");
	TntEventTree2->Branch(brN,PmtFrontHit,TString::Format("PmtFrontHit[%d][%d]/I", NX,NY));
	//TntEventTree->Branch(brN,&PmtFrontHit);
	//G4cout << "TEST OF GETBRANCH FILL " << evid << " " << brN << G4endl;

	sprintf(brN, "PMT_Back_Hit_Event_%d", evid);
	//sprintf(brN2, "PMT_Back_Hit_Event_%d[10][10]/I", evid);
	//by Shuya 160509
	//TntEventTree2->Branch(brN,PmtBackHit,"PmtBackHit[8][8]/I");
	TntEventTree2->Branch(brN,PmtBackHit,TString::Format("PmtBackHit[%d][%d]/I", NX,NY));
	//TntEventTree->Branch(brN,&PmtBackHit);
	//G4cout << "TEST OF GETBRANCH FILL " << evid << " " << brN << G4endl;

//by Shuya 160426. Removed this because of replacement with TTree branch.
/*
//G4cout << "!!!!CREATE!!! " << evid << G4endl;
sprintf(brN, "Count_PMT_Front_run%d",evid);
h_Count_PMT_Front.push_back (new TH2I(brN,brN, 10, 0, 10, 10, 0, 10));
h_Count_PMT_Front[evid]->GetXaxis()->SetTitle("X (no.)");
h_Count_PMT_Front[evid]->GetYaxis()->SetTitle("Y (no.)");

*/

	PhotonSum.clear();

	// re-create digitizer histogram
	// x-axis: signals as recorded by a CAEN V1730 digitizer (bins of 2 ns)
	// y-axis: PMT ID (0->16)
	if(hDigi) { hDigi->Delete(); }
	hDigi = new TH2I("hDigi", "Digitizer Signals", 300, 0, 600, NX*NY, 0, NX*NY);
	hDigi->SetDirectory(0);
}


void TntDataRecordTree::senddataEV(int type, double value1)
{
	// Writes data values from TntEndofEventAction to EventTree
	// Changes from Serra >>> Light Conversion in TntSD.cc, EndOfEvent()
	// Hard coded for Tnt
	switch(type)
	{
	case 1:
		eng_Tnt = value1;  // Sum of all energy in event!
		if (eng_Tnt > Det_Threshold)
		{number_total++;
			//G4cout << "!!! " << eng_Tnt << G4endl;
//G4cout << (TH1D*)DataFile->Get("Energy_Tnt") << "!! !!" << G4endl;
//by Shuya 160426. Removed the histogram to replace with TTree.
			//((TH1D*)DataFile->Get("Energy_Tnt"))->Fill(eng_Tnt);
			//h_Energy_Tnt->Fill(eng_Tnt);
			//h_Energy_Tnt->Fill(10, 20);
		}
		break;
	case 2:
		eng_Tnt_proton = value1; 
		if (eng_Tnt_proton > Det_Threshold)
	  {number_protons++;
//by Shuya 160426. Removed the histogram to replace with TTree.
			//((TH1D*)DataFile->Get("Energy_Proton"))->Fill(eng_Tnt_proton);
			//h_Energy_Proton->Fill(eng_Tnt_proton);
		}
		break;
	case 3:
		eng_Tnt_alpha = value1;
		if (eng_Tnt_alpha > Det_Threshold) 
	  {number_alphas++;
//by Shuya 160426. Removed the histogram to replace with TTree.
			//((TH1D*)DataFile->Get("Energy_Alpha"))->Fill(eng_Tnt_alpha);
			//h_Energy_Alpha->Fill(eng_Tnt_alpha);
		}
		break;
	case 4:
		eng_Tnt_C12 = value1; 
		if (eng_Tnt_C12 > Det_Threshold)  
	  {number_C12++;
//by Shuya 160426. Removed the histogram to replace with TTree.
			//((TH1D*)DataFile->Get("Energy_C12"))->Fill(eng_Tnt_C12);
			//h_Energy_C12->Fill(eng_Tnt_C12);
		}
		break;
	case 5:
		eng_Tnt_EG = value1; 
		if (eng_Tnt_EG > Det_Threshold)  
	  {number_EG++;
//by Shuya 160426. Removed the histogram to replace with TTree.
			//((TH1D*)DataFile->Get("Energy_EG"))->Fill(eng_Tnt_EG);
			//h_Energy_EG->Fill(eng_Tnt_EG);
		}
		break;
	case 6:
		eng_Tnt_Exotic = value1;  
		if (eng_Tnt_Exotic > Det_Threshold)
	  {number_Exotic++;
//by Shuya 160426. Removed the histogram to replace with TTree.
			//((TH1D*)DataFile->Get("Energy_Exotic"))->Fill(eng_Tnt_Exotic);
			//h_Energy_Exotic->Fill(eng_Tnt_Exotic);
		}
		break;
//by Shuya 160407
	case 7:
		eng_Tnt_PhotonFront = (int)value1;  
		if (eng_Tnt_PhotonFront > Det_Threshold)
	  {number_Photon++;
//by Shuya 160426. Removed the histogram to replace with TTree.
			//((TH1D*)DataFile->Get("Energy_Photon"))->Fill(eng_Tnt_Photon);
			//h_Energy_Photon->Fill(eng_Tnt_Photon);
		}
		break;
//by Shuya 160427
	case 8:
		eng_Tnt_PhotonBack = (int)value1;  
		if (eng_Tnt_PhotonBack > Det_Threshold)
	  {number_Photon++;
		}
		break;
//by Shuya 160502
	case 9:
		eng_Tnt_PhotonTotal = (int)value1;  
		if (eng_Tnt_PhotonTotal > Det_Threshold)
	  {number_Photon++;
		}
		break;
//by Shuya 160502
	case 10:
		edep_Tnt = value1;  
		break;
	case 11:
		edep_Tnt_proton = value1;  
		break;
	case 12:
		edep_Tnt_alpha = value1;  
		break;
	case 13:
		edep_Tnt_C12 = value1;  
		break;
	case 14:
		edep_Tnt_EG = value1;  
		break;
	case 15:
		edep_Tnt_Exotic = value1;  
	case 16:
		num_Tnt_NonPMT = value1;  
	case 17:
		num_Tnt_Abs = value1;  
		break;
	default:
		G4cout << "Data Transfer Error!" << G4endl;
		G4ExceptionSeverity severity=FatalException;
		G4Exception("Program Aborted in TntDataRecordTree.cc::senddataEV()", "DataError",severity,"Data Transfer Error!");
		//Note: the enum G4ExceptionSeverity has been defined since Geant4 release 5.0 and its values are: FatalException, FatalErrorInArgument, EventMustBeAborted, FatalException, and JustWarning. 
		//G4Exception("Program Aborted in TntDataRecordTree.cc::senddataEV()");
		break;
	}
	// cout << value1 << " <- This was put in the data file!" << endl;
}

void TntDataRecordTree::senddataPosition(const G4ThreeVector& pos)
{
	// const G4double pi = 3.14159265;
	FirstHitMag = 0;
	// HitAngle = 360;       // Some value that one would never get! 
	Xpos = pos(0);
	Ypos = pos(1);
	Zpos = pos(2);

	FirstHitMag = sqrt(pow(Xpos,2)+pow(Ypos,2)+pow(Zpos,2))/cm;

//by Shuya 160502
	FirstHitPosition = pos;

	// HitAngle = acos(x/HitRad)*180/pi;

	//   G4cout << "Hit Distance s = " << FirstHitMag << G4endl;
}


void TntDataRecordTree::senddataHits(const std::vector<TntDataRecordTree::Hit_t>& hits, bool sortTime)
{
	// Reset hit vectors
	HitX.resize(0);
	HitY.resize(0);
	HitZ.resize(0);
	HitT.resize(0);
	HitE.resize(0);
	HitTrackID.resize(0);
	NumHits = 0;
	fHits->Clear();
	fHit01->Clear();

	if(hits.empty()) { return; } // Nothing to do

	//
	// Non-trivial hit vector
	HitX.reserve(hits.size());
	HitY.reserve(hits.size());
	HitZ.reserve(hits.size());
	HitT.reserve(hits.size());
	HitE.reserve(hits.size());
	HitTrackID.reserve(hits.size());
	NumHits = hits.size();
	
	for(std::vector<Hit_t>::const_iterator it = hits.begin();
			it != hits.end(); ++it)
	{
		if(sortTime == false) { // insert at end of vector
			HitX.push_back(it->X);
			HitY.push_back(it->Y);
			HitZ.push_back(it->Z);
			HitT.push_back(it->T);
			HitE.push_back(it->E);
			HitTrackID.push_back(it->TrackID);
		}	else { // insert, sorted by time vector
			std::vector<G4double>::iterator iT = 
				std::lower_bound(HitT.begin(), HitT.end(), it->T);
			std::vector<G4double>::iterator iX = 
				(iT - HitT.begin()) + HitX.begin();
			std::vector<G4double>::iterator iY = 
				(iT - HitT.begin()) + HitY.begin();
			std::vector<G4double>::iterator iZ = 
				(iT - HitT.begin()) + HitZ.begin();
			std::vector<G4double>::iterator iE = 
				(iT - HitT.begin()) + HitE.begin();
			std::vector<G4int>::iterator iTrackID = 
				(iT - HitT.begin()) + HitTrackID.begin();

			HitX.insert(iX, it->X);
			HitY.insert(iY, it->Y);
			HitZ.insert(iZ, it->Z);
			HitT.insert(iT, it->T);
			HitE.insert(iE, it->E);
			HitTrackID.insert(iTrackID, it->TrackID);
		}
	}

	if(NumHits != HitT.size()) { 
		G4cerr << "ERROR:: NUM HITS, HitT.size():: " << NumHits << ", " << HitT.size() << G4endl; 
	}

	// 
	// Fill TLorrentzVector arrays
	for(size_t i=0; i< HitT.size(); ++i) {
		TLorentzVector* hit4Vector = new( (*fHits)[i] ) TLorentzVector();
		hit4Vector->SetXYZT(HitX.at(i), HitY.at(i), HitZ.at(i), HitT.at(i));
	}
	if(fHits->GetEntries() > 1) {
		TLorentzVector* hit1 = (TLorentzVector*)fHits->At(1);
		TLorentzVector* hit0 = (TLorentzVector*)fHits->At(0);
		new( (*fHit01)[0] ) TLorentzVector( *hit1 - *hit0 );
	}
}


void TntDataRecordTree::senddataTOF(G4double time)
{
	FirstHitTime = time/ns;
	// G4cout << "The first hit time was : " << FirstHitTime << G4endl;
}

 
void TntDataRecordTree::ShowDataFromEvent()
{
#if 0
	G4cout << "================OUTPUT SENT FROM TntDataRecordTree====================================" << G4endl;
	G4cout << "The Energy of the Initial Particle was:    " << eng_int << G4endl;
	G4cout << "The Position of the First Hit was (cm <-note!) (Distance from (x=0,y=0,z=0)):    " << FirstHitMag << G4endl;
	G4cout << "The Position of the First Hit was (mm <-note!) (x,y,z):    " << FirstHitPosition << G4endl;
	G4cout << "The Time of Flight of First Hit was (ns):  " << FirstHitTime << G4endl;
	G4cout << "Measured Total Light Output (unless otherwise light_Conv=NULL) in this Event:  " << eng_Tnt  << G4endl;
	G4cout << "Measured Proton Light Output (unless otherwise light_Conv=NULL) in this Event: " << eng_Tnt_proton << G4endl;
	G4cout << "Measured Alpha Light Output (unless otherwise light_Conv=NULL) in this Event:  " << eng_Tnt_alpha << G4endl;
	G4cout << "Measured C12 Light Output (unless otherwise light_Conv=NULL) in this Event:    " << eng_Tnt_C12 << G4endl;
	G4cout << "Measured Electron (/Gamma) Eng. Loss:               " << eng_Tnt_EG << G4endl;
	G4cout << "Measured Energy Loss from Exotic Particles (Be9,etc):" << eng_Tnt_Exotic << G4endl;
	G4cout << "Measured Energy Loss from Photons (Front):" << eng_Tnt_PhotonFront << G4endl;
	G4cout << "Measured Energy Loss from Photons (Back):" << eng_Tnt_PhotonBack << G4endl;
//by Shuya 160502
	G4cout << "Measured Total Energy Loss in this Event (MeV):  " << edep_Tnt  << G4endl;
	G4cout << "Measured Proton Energy Loss in this Event (MeV): " << edep_Tnt_proton << G4endl;
	G4cout << "Measured Alpha Energy Loss in this Event (MeV):  " << edep_Tnt_alpha << G4endl;
	G4cout << "Measured C12 Energy Loss in this Event (MeV):    " << edep_Tnt_C12 << G4endl;
	G4cout << "Measured Electron/Gamma Eng. Loss (MeV):               " << edep_Tnt_EG << G4endl;
	G4cout << "Measured Energy Loss from Exotic Particles (Be9,etc) (MeV):" << edep_Tnt_Exotic << G4endl;
#endif
}

void TntDataRecordTree::FillTree()
{
	if (eng_Tnt > Det_Threshold)  // Threshold set in main()
	{number_at_this_energy++;}
	TntEventTree->Fill();  

	//G4cout << "FillTree1!" << G4endl;
}

//by Shuya 160422
void TntDataRecordTree::FillTree2(int evid)
{
	char brN[300];

	sprintf(brN, "PMT_Front_Hit_Event_%d", evid);
	//sprintf(brN2, "PMT_Front_Hit_Event_%d/I", evid);
	//G4cout << "TEST OF GETBRANCH FILL " << evid << " " << brN << G4endl;

	///for(int i=0;i<NX;i++)
	//{
	//	for(int j=0;j<NY;j++)	G4cout << PmtFrontHit[i][j] << G4endl;
	//for(int j=0;j<NY;j++)	G4cout << "TEST PMTFRONTHIT ADDRESS " << &PmtFrontHit[i][j] << G4endl;
	//}
	
	TntEventTree2->GetBranch(brN)->Fill();
	//TntEventTree->Branch(brN,&PmtFrontHit);

	sprintf(brN, "PMT_Back_Hit_Event_%d", evid);
	//sprintf(brN2, "PMT_Back_Hit_Event_%d[10][10]/I", evid);
	//G4cout << "TEST OF GETBRANCH FILL " << evid << " " << brN << G4endl;

	//for(int i=0;i<NX;i++)
	//{
	//	for(int j=0;j<NY;j++)	G4cout << PmtBackHit[i][j] << G4endl;
	//}
	
	TntEventTree2->GetBranch(brN)->Fill();
	//TntEventTree->Branch(brN,&PmtBackHit);
  	
	//TntEventTree2->Fill();  

//Initialization.
	//for(int i=0;i<10;i++)
	//by Shuya 160509
	for(int i=0;i<NX;i++)
	{
		//for(int j=0;j<10;j++)
		for(int j=0;j<NY;j++)
		{
			//G4cout << "TEST NO 422 " << PmtFrontHit[i][j] << G4endl;
			PmtFrontHit[i][j] = 0;
			PmtBackHit[i][j] = 0;
		}
	}
}

void TntDataRecordTree::GetParticleTotals()
{
	cout << "The Initial Number of Beam Particles was: " << event_counter << endl;
	cout << "The Detection Threshold is set at : " << Det_Threshold << " MeVee." << endl;
	cout << "The Total Number of Detected Events was:  " << number_total << endl;
	cout << "The Total Number of Protons Detected was: " << number_protons << endl;
	cout << "The Total Number of Alphas Detected was:  " << number_alphas << endl;
	cout << "The Total Number of C12 Detected was:     " << number_C12 << endl;
	cout << "The Total Number of e- or e+    was:      " << number_EG << endl;
	cout << "The Total Number of Exotic Particles was: " << number_Exotic << endl;
	cout << "The Total Number of Photons was: " << number_Photon << endl;
}

void TntDataRecordTree::CalculateEff(int ch_eng)
{
	char EffFile[] = "eff_results_file.dat";
	ofstream outfile2(EffFile,ios::app);

	cout << number_at_this_energy << endl;
	efficiency = 100*(static_cast<double>(number_at_this_energy)/static_cast<double>(ch_eng));
	cout << "Efficiency was: " << number_at_this_energy << "/" << ch_eng << " = " << efficiency << " %" << endl;

	outfile2 << setiosflags(ios::fixed)
					 << setprecision(4)
					 << eng_int << "  "
					 << efficiency 
					 << endl;
	outfile2.close();

	number_at_this_energy = 0;
}
 
