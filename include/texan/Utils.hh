/// \file TexanUtils.hh
/// \brief Short inline utility functions


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

}
