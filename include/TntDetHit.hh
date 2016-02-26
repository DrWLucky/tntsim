//==========================================================================
// TntDetHit.hh
// Based on DemonScintHit.hh and ExN04TrackerHit.hh
//
// Written/Modified by: Brian Roeder, LPC Caen 02/14/07
//                      email - roeder@lpccaen.in2p3.fr
//
// Purpose: Defines hit collection and values to be stored by TntSD.cc
//
//==========================================================================
//
// - See UM Hits Presentation and JLab Scoring 2 Talk for more info.


#ifndef TntDetHit_h
#define TntDetHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "G4PhysicalConstants.hh"


class TntDetHit : public G4VHit
{
  public:

  TntDetHit();        // Constructor
  ~TntDetHit();       // Destructor

  /* //From ExN04TrackerHit.cc - Not needed... 
  TntDetHit(const TntDetHit &right);
  const TntDetHit& operator=(const TntDetHit &right);
  G4int operator==(const TntDetHit &right) const;
  */

  inline void * operator new(size_t);
  inline void operator delete(void *aHit);
  

  void Draw();       // needed to Draw Hit in Visualization
  void Print();      // Print out Event (if not invoked in EndofEventAction)

  private:
  // These are the Data Storage objects of the class.
  // Set/Stored and Accessed through Public Methods below.

  G4double edep;           // Records Energy Deposited in Hit
  G4ThreeVector pos;       // Records Hit Position in Detector
  G4int TrackID;           // Records TrackID of Particle in Hit
  G4int ParentTrackID;     // Records TrackID of Parent Particle in Event
  G4double EvtTOF;         // Records Global Time of "Hit"
  G4String ParticleName;   // Records particle ID in reactions 
  G4double ParticleCharge; // Records Charge of Particle (PDG Charge!)
  G4double ParticleA;      // Records A of Particle (where "A" = Baryon Num)
  G4String CreatorProcess; // Records the "Process" creating the Track in Hit

 public:

  // Set Data Methods
      void SetEdep(G4double de)
      { edep = de; }
      void SetPos(G4ThreeVector xyz)
      { pos = xyz; }
      void SetTrackID(G4int track)
      {TrackID = track;}
      void SetParentTrackID(G4int parentTrack)
      {ParentTrackID = parentTrack;}
      void SetTOF(G4double gTOF)
      {EvtTOF = gTOF;}
      void SetParticleName(G4String theParticleName)
      {ParticleName = theParticleName;}
      void SetParticleCharge(G4double theParticleCharge)
      {ParticleCharge = theParticleCharge;}
      void SetParticleA(G4double theParticleMass)
      {ParticleA = theParticleMass;}
      void SetParticleProcess(G4String theProcess)
      {CreatorProcess = theProcess;}
     
  // Get Data Methods (For EndofEventAction and DataRecordTree)
      G4double GetEdep()
      { return edep; }
      G4ThreeVector GetPos()
      { return pos; }
      G4int GetTrackID()
      { return TrackID; }
      G4int GetParentTrackID()
      {return ParentTrackID; }
      G4double GetTOF()
      {return EvtTOF; }
      G4String GetParticleName()
      { return ParticleName; }
      G4double GetParticleCharge()
      {return ParticleCharge; }
      G4double GetParticleA()
      {return ParticleA; }
      G4String GetParticleProcess()
      {return CreatorProcess;}
};

typedef G4THitsCollection<TntDetHit> TntDetHitsCollection;

// These functions setup the memory functions for the Hit

extern G4Allocator<TntDetHit> TntDetHitAllocator;

// Create the Hit

inline void* TntDetHit::operator new(size_t)
{
  void *aHit;
  aHit = (void *)TntDetHitAllocator.MallocSingle();
  return aHit;
}

// Delete the Hit (at the end of the event, for example)

inline void TntDetHit::operator delete(void *aHit)
{
  TntDetHitAllocator.FreeSingle((TntDetHit*) aHit);
}

#endif
