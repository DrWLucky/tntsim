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
 
#include "TntDataRecordTree.hh"


// Access to Analysis pointer! (see TntSD.cc EndOfEvent() for Example)
TntDataRecordTree* TntDataRecordTree::TntPointer;

TntDataRecordTree::TntDataRecordTree(G4double Threshold) : 
  // Initialized Values
  eng_int(0), eng_Tnt(0), eng_Tnt_alpha(0), eng_Tnt_C12(0), 
  eng_Tnt_EG(0), eng_Tnt_Exotic(0),  FirstHitTime(0), FirstHitMag(0),
  Xpos(0), Ypos(0), Zpos(0), Det_Threshold(Threshold),
  event_counter(0), number_total(0), 
  number_protons(0), number_alphas(0), number_C12(0), number_EG(0), 
  number_Exotic(0), number_at_this_energy(0), efficiency(0)

{ /* Constructor */ 
  TntPointer = this;  // When Pointer is constructed, assigns address of this class to it.
  //
  // Create new data storage text file
  // Create new text file for data storage - (Data Recorded by TntDataRecordTree class)
 
  char EffFile[] = "eff_results_file.dat";
  ofstream outfile2(EffFile,ios::trunc);
  outfile2.close();

  //
  // Create pointers to ROOT Analysis Tree and Branches
  //
  cout <<"\n Starting Data Tree Constructor" << endl;
  
  const Char_t* evt_file = "TntDataTree.root";

  DataFile = new TFile(evt_file, "RECREATE");
  
  TntEventTree = new TTree("t","Tnt Scintillator Simulation Data");
  TntEventTree->Branch("Energy_Initial",&eng_int,"eng_int/D");
  TntEventTree->Branch("Energy_Tnt",&eng_Tnt,"eng_Tnt/D");
  TntEventTree->Branch("Energy_Proton",&eng_Tnt_proton,"eng_Tnt_proton/D");
  TntEventTree->Branch("Energy_Alpha",&eng_Tnt_alpha,"eng_Tnt_alpha/D");
  TntEventTree->Branch("Energy_C12", &eng_Tnt_C12,"eng_Tnt_C12/D");
  TntEventTree->Branch("Energy_EG",&eng_Tnt_EG,"eng_Tnt_EG/D");
  TntEventTree->Branch("Energy_Exotic",&eng_Tnt_Exotic,"eng_Tnt_Exotic/D");

  TntEventTree->Branch("First_Hit_Pos",&FirstHitMag,"FirstHitMag/D");
  TntEventTree->Branch("First_Hit_Time",&FirstHitTime,"FirstHitTime/D");

  TntEventTree->Branch("Xpos",&Xpos,"Xpos/D");
  TntEventTree->Branch("Ypos",&Ypos,"Ypos/D");
  TntEventTree->Branch("Zpos",&Zpos,"Zpos/D");

  cout << TntPointer << endl;
}

TntDataRecordTree::~TntDataRecordTree()
{/* Destructor, Close root file */
 
  DataFile->Write(); 
  DataFile->Close();
  cout << "############################################################" << endl;
  cout << "Created Root Tree File = \"TntDataTree.root\"" << endl;
  cout << "Got the ROOT Tree File from Data Record Routine!" << endl; 
  delete DataFile;
}

// Member functions of TntDataRecordTree class

void TntDataRecordTree::senddataPG(double value1=0.)
{
    eng_int = value1;
    event_counter++;
    //  cout << "eng_int = " << eng_int << endl;
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
	{number_total++;}
	break;
      case 2:
	eng_Tnt_proton = value1; 
	if (eng_Tnt_proton > Det_Threshold)
	  {number_protons++;}
	break;
      case 3:
	eng_Tnt_alpha = value1;
	if (eng_Tnt_alpha > Det_Threshold) 
	  {number_alphas++;}
	break;
      case 4:
	eng_Tnt_C12 = value1; 
	if (eng_Tnt_C12 > Det_Threshold)  
	  {number_C12++;}
	break;
      case 5:
	eng_Tnt_EG = value1; 
	if (eng_Tnt_EG > Det_Threshold)  
	  {number_EG++;}
	break;
      case 6:
	eng_Tnt_Exotic = value1;  
	if (eng_Tnt_Exotic > Det_Threshold)
	  {number_Exotic++;}
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

void TntDataRecordTree::senddataPosition(G4ThreeVector pos)
{
  // const G4double pi = 3.14159265;
  FirstHitMag = 0;
  // HitAngle = 360;       // Some value that one would never get! 
  Xpos = pos(0);
  Ypos = pos(1);
  Zpos = pos(2);

    FirstHitMag = sqrt(pow(Xpos,2)+pow(Ypos,2)+pow(Zpos,2))/cm;

    // HitAngle = acos(x/HitRad)*180/pi;

    //   G4cout << "Hit Distance s = " << FirstHitMag << G4endl;
}

void TntDataRecordTree::senddataTOF(G4double time)
{
  FirstHitTime = time/ns;
  // G4cout << "The first hit time was : " << FirstHitTime << G4endl;
}

 
void TntDataRecordTree::ShowDataFromEvent()
{
      G4cout << "====================================================" << G4endl;
      G4cout << "The Energy of the Initial Particle was:    " << eng_int << G4endl;
      G4cout << "The Position of the First Hit was (cm):    " << FirstHitMag << G4endl;
      G4cout << "The Time of Flight of First Hit was (ns):  " << FirstHitTime << G4endl;
      G4cout << "Measured Total Energy Loss in this Event:  " << eng_Tnt  << G4endl;
      G4cout << "Measured Proton Energy Loss in this Event: " << eng_Tnt_proton << G4endl;
      G4cout << "Measured Alpha Energy Loss in this Event:  " << eng_Tnt_alpha << G4endl;
      G4cout << "Measured C12 Energy Loss in this Event:    " << eng_Tnt_C12 << G4endl;
      G4cout << "Measured Electron Eng. Loss:               " << eng_Tnt_EG << G4endl;
      G4cout << "Measured Energy Loss from Exotic Particles:" << eng_Tnt_Exotic << G4endl;
}

void TntDataRecordTree::FillTree()
{
  if (eng_Tnt > Det_Threshold)  // Threshold set in main()
     {number_at_this_energy++;}
    TntEventTree->Fill();  
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
 
