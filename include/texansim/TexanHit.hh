/// \file TexanHit.hh
/// \brief Definition of the TexanHit class
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


class TexanSD;

/// TEXAN neutron array hit class
class TexanHit : public G4VHit
{
public:
	/// Ctor
	TexanHit();
	/// Dtor
	virtual ~TexanHit();

	/// Special equivalence operator as specified by GEANT
	G4int operator==(const TexanHit&) const;
	/// Special new operator as specified by GEANT
	inline void* operator new(size_t);
	/// Special delete operator as specified by GEANT
	inline void  operator delete(void*);

	/// Draw an event
	virtual void Draw();
	/// Print an event 
	virtual void Print();

	/// Get read-only access to hit data
	const HitData& GetData() const { return fData; }

private:
	/// Contains hit data
	HitData fData;

	/// Let associated sensitive detector class write to this one
	friend class TexanSD;
};


/// Templated hits collection
typedef G4THitsCollection<texansim::TexanHit> TexanHitsCollection;
/// Hit collection allocator
extern G4ThreadLocal G4Allocator<texansim::TexanHit>* TexanHitAllocator;

}


inline void* texansim::TexanHit::operator new(size_t)
{
  if(!TexanHitAllocator)
		TexanHitAllocator = new G4Allocator<TexanHit>;
  return (void *) TexanHitAllocator->MallocSingle();
}

inline void texansim::TexanHit::operator delete(void *hit)
{
  TexanHitAllocator->FreeSingle((TexanHit*) hit);
}


#endif
