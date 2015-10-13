/// \file VPersistenceManager.hh
/// \brief Defines abstract class to implement saving of class objects

/**  Goal is to do this using ROOT/TTree... but for generality's sake,
 *   I'll make an abstract class. Maybe someday someone will make another
 *   implementation.
 */
#ifndef TXS_PERSISTENCE_HEADER_FILE_8675309
#define TXS_PERSISTENCE_HEADER_FILE_8675309
#include "G4String.hh"
#include "G4Types.hh"

namespace texansim {

/// \brief Pure abstract persistence manager class
/// \detailed Also includes some histogram functionality
class VPersistenceManager {
public:
	enum Type_t { kDouble = 0, kFloat = 1, kInt = 2 };

public:
	/// Ctor (empty)
	VPersistenceManager() { }
	/// Dtor (empty)
	virtual ~VPersistenceManager() { }
	/// Set output file name (Base - no extension)
	virtual void SetFilename(const G4String& name) = 0;
	/// Get output file name
	virtual const G4String& GetFilename() const = 0;
	/// Open output file
	virtual G4bool OpenFile() = 0;
	/// Add class object to be saved
	virtual G4bool AddObject(const G4String& name, const G4String& classname, void*) = 0;
	/// Add primitive type (double, float, int) to be saved
	virtual G4bool AddPrimitive(const G4String& name, void* p, Type_t type) = 0;
	/// Add 1d histogram
	virtual G4bool AddHistogram1d(
		const G4String& name, const G4String& title,
		G4int bins, G4double xlow, G4double xhigh,
		void* valuePointer, Type_t type) = 0;
	/// Add 2d histogram
	virtual G4bool AddHistogram2d(
		const G4String& name, const G4String& title,
		G4int xbins, G4double xlow, G4double xhigh,
		G4int ybins, G4double ylow, G4double yhigh,
		void* valuePointerx, void* valuePointery,
		Type_t type) = 0;
	/// Save an event of histograms and persistent objects
	virtual void SaveEvent() = 0;
	/// Write to disk
	virtual void Write() = 0;
	/// Close everything out
	virtual void Close() = 0;
};

}



#endif
