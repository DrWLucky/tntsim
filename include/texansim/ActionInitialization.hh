/// \file ActionInitialization.hh
/// \brief Definition of the TexanActionInitialization class

#ifndef TexanActionInitialization_h
#define TexanActionInitialization_h 1

#include "G4VUserActionInitialization.hh"


namespace texansim {

/// Action initialization class.
class ActionInitialization : public G4VUserActionInitialization
{
public:
	ActionInitialization();
	virtual ~ActionInitialization();

	/// For main thread (ignored if serial)
	virtual void BuildForMaster() const;
	/// For all threads
	virtual void Build() const;
};

}


#endif
    
