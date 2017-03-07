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
// $Id: TntScintHit.hh 72250 2013-07-12 08:59:26Z gcosmo $
//
/// \file optical/Tnt/include/TntScintHit.hh
/// \brief Definition of the TntScintHit class
//
//
#ifndef TntScintHit_h
#define TntScintHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "G4LogicalVolume.hh"
#include "G4Transform3D.hh"
#include "G4RotationMatrix.hh"
#include "G4VPhysicalVolume.hh"

#include "tls.hh"

class TntScintHit : public G4VHit
{
  public:
 
    TntScintHit();
    TntScintHit(G4VPhysicalVolume* pVol);
    virtual ~TntScintHit();
    TntScintHit(const TntScintHit &right);
    const TntScintHit& operator=(const TntScintHit &right);
    G4int operator==(const TntScintHit &right) const;

    inline void *operator new(size_t);
    inline void operator delete(void *aHit);
 
    virtual void Draw();
    virtual void Print();

    inline void SetEdep(G4double de) { fEdep = de; }
    inline void AddEdep(G4double de) { fEdep += de; }
    inline G4double GetEdep() { return fEdep; }

    inline void SetPos(G4ThreeVector xyz) { fPos = xyz; }
    inline G4ThreeVector GetPos() { return fPos; }

    inline const G4VPhysicalVolume * GetPhysV() { return fPhysVol; }

//by Shuya 160408
  // Set Data Methods
      //void SetEdep(G4double de)
      //{ edep = de; }
      //void SetPos(G4ThreeVector xyz)
      //{ pos = xyz; }
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

//by Shuya 160407
  // Get Data Methods (For EndofEventAction and DataRecordTree)
//      G4double GetEdep()
//      { return edep; }
//      G4ThreeVector GetPos()
//      { return pos; }
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

  private:
    G4double fEdep;
    G4ThreeVector fPos;
    const G4VPhysicalVolume* fPhysVol;

  //G4double edep;           // Records Energy Deposited in Hit
  //G4ThreeVector pos;       // Records Hit Position in Detector
  G4int TrackID;           // Records TrackID of Particle in Hit
  G4int ParentTrackID;     // Records TrackID of Parent Particle in Event
  G4double EvtTOF;         // Records Global Time of "Hit"
  G4String ParticleName;   // Records particle ID in reactions 
  G4double ParticleCharge; // Records Charge of Particle (PDG Charge!)
  G4double ParticleA;      // Records A of Particle (where "A" = Baryon Num)
  G4String CreatorProcess; // Records the "Process" creating the Track in Hit

};

typedef G4THitsCollection<TntScintHit> TntScintHitsCollection;

extern G4ThreadLocal G4Allocator<TntScintHit>* TntScintHitAllocator;

inline void* TntScintHit::operator new(size_t)
{
  if(!TntScintHitAllocator)
      TntScintHitAllocator = new G4Allocator<TntScintHit>;
  return (void *) TntScintHitAllocator->MallocSingle();
}

inline void TntScintHit::operator delete(void *aHit)
{
  TntScintHitAllocator->FreeSingle((TntScintHit*) aHit);
}

#endif
