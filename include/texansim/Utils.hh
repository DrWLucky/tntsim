/// \file Utils.hh
/// \brief Short inline utility functions
#include <sstream>
#include "G4String.hh"


namespace texansim {

/// Delete a pointer and set to NULL
template <class T>
void Zap(T*& t)
{
	if(t) {
		delete t;
		t = NULL;
	}
}

/// Reset a pointer to a new value, deleting the old one
template<class T>
void ResetPointer(T*& t, T* newValue)
{
	if(t) {
		delete t;
	}

	t = newValue;
}


/// Format string + integer 
template <class T>
G4String FormatStr1(const G4String& str, const T& other)
{
	std::stringstream sstr;
	sstr << str.data() << other;
	return G4String(sstr.str().c_str());
}

}
