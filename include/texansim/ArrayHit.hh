/// \file ArrayHit.hh
/// \brief Definition of the B2TrackerHit class
///
#ifndef TEXAN_ARRAY_HIT_HH
#define TEXAN_ARRAY_HIT_HH 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"
#include "tls.hh"



namespace texansim {

/// Array hit class
/**
	 It defines data members to store the trackID, chamberNb, energy deposit,
	 and position of charged particles in a selected volume:
	 - fTrackID, fChamberNB, fEdep, fPos
*/
class ArrayHit : public G4VHit
{
public:
	ArrayHit();
	ArrayHit(const ArrayHit&);
	virtual ~ArrayHit();

	// operators
	const ArrayHit& operator=(const ArrayHit&);
	G4int operator==(const ArrayHit&) const;

	inline void* operator new(size_t);
	inline void  operator delete(void*);

	// methods from base class
	virtual void Draw();
	virtual void Print();

	// Set methods
	void SetTrackID  (G4int track)      { fTrackID = track; };
	void SetChamberNb(G4int chamb)      { fChamberNb = chamb; };
	void SetEdep     (G4double de)      { fEdep = de; };
	void SetPos      (G4ThreeVector xyz){ fPos = xyz; };

	// Get methods
	G4int GetTrackID() const     { return fTrackID; };
	G4int GetChamberNb() const   { return fChamberNb; };
	G4double GetEdep() const     { return fEdep; };
	G4ThreeVector GetPos() const { return fPos; };

private:

	G4int         fTrackID;
	G4int         fChamberNb;
	G4double      fEdep;
	G4ThreeVector fPos;
};


// typedefs 

typedef G4THitsCollection<texansim::ArrayHit> ArrayHitsCollection;

extern G4ThreadLocal G4Allocator<texansim::ArrayHit>* ArrayHitAllocator;

}


inline void* texansim::ArrayHit::operator new(size_t)
{
  if(!ArrayHitAllocator)
		ArrayHitAllocator = new G4Allocator<ArrayHit>;
  return (void *) ArrayHitAllocator->MallocSingle();
}

inline void texansim::ArrayHit::operator delete(void *hit)
{
  ArrayHitAllocator->FreeSingle((ArrayHit*) hit);
}


#endif
