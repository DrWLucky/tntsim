/// \file TexanActionInitialization.hh
/// \brief Definition of the TexanActionInitialization class

#ifndef TexanActionInitialization_h
#define TexanActionInitialization_h 1

#include "G4VUserActionInitialization.hh"

/// Action initialization class.

class TexanActionInitialization : public G4VUserActionInitialization
{
public:
	TexanActionInitialization();
	virtual ~TexanActionInitialization();

	/// For main thread (ignored if serial)
	virtual void BuildForMaster() const;
	/// For all threads
	virtual void Build() const;
};


#endif
    
