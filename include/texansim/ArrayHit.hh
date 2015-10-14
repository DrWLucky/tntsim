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

#include "texansim/HitData.hh"


namespace texansim {


class ArraySD;

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
	// ArrayHit(const ArrayHit&);
	virtual ~ArrayHit();

	// operators
	// const ArrayHit& operator=(const ArrayHit&);
	G4int operator==(const ArrayHit&) const;

	inline void* operator new(size_t);
	inline void  operator delete(void*);

	// methods from base class
	virtual void Draw();
	virtual void Print();

	// Get access to hit data
	const HitData& GetData() const { return fData; }

private:
	HitData fData;

	friend class ArraySD;
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
