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


#include "G4Step.hh"
// class G4Step;

namespace tntsim {


class TexanSD;

/// TEXAN neutron array hit class
class TexanHit : public G4VHit
{
public:
	/// Ctor
	TexanHit();
	/// Sets step
	TexanHit(G4Step*);
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

	/// Get step
	const G4Step* GetStep() const { return fStep; }
	/// Set step
	void SetStep(G4Step* step) { fStep = new G4Step(*step); }

private:
	/// Store pointer to G4Step for later processing
	G4Step* fStep;
};


/// Templated hits collection
typedef G4THitsCollection<tntsim::TexanHit> TexanHitsCollection;
/// Hit collection allocator
extern G4ThreadLocal G4Allocator<tntsim::TexanHit>* TexanHitAllocator;

}


inline void* tntsim::TexanHit::operator new(size_t)
{
  if(!TexanHitAllocator)
		TexanHitAllocator = new G4Allocator<TexanHit>;
  return (void *) TexanHitAllocator->MallocSingle();
}

inline void tntsim::TexanHit::operator delete(void *hit)
{
  TexanHitAllocator->FreeSingle((TexanHit*) hit);
}


#endif
