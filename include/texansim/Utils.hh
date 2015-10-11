/// \file TexanUtils.hh
/// \brief Short inline utility functions
#include <sstream>
#include "G4String.hh"


namespace {

/// Delete a pointer and set to NULL
template <class T>
void Zap(T*& t)
{
	if(t) {
		delete t;
		t = NULL;
	}
}

template <class T>
G4String FormatStr1(const G4String& str, const T& other)
{
	std::stringstream sstr;
	sstr << str.data() << other;
	return G4String(sstr.str().c_str());
}

}
