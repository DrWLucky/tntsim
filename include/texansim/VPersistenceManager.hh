/// \file VPersistenceManager.hh
/// \brief Defines abstract class to implement saving of class objects
/**  \details Goal is to do this using ROOT/TTree... but for generality's sake,
 *   I'll make an abstract class. Maybe someday someone will make another
 *   implementation.
 */
#ifndef TXS_PERSISTENCE_HEADER_FILE_8675309
#define TXS_PERSISTENCE_HEADER_FILE_8675309
#include "G4String.hh"
#include "G4Types.hh"



class G4UImessenger;

namespace texansim {


/// Pure abstract persistence manager class
/*! Interface for saving class objects to disk.
 *  Also includes functionality to save histograms
 *  (1d and 2d) as well as "Primitive" numbers (int,
 *  float, double). Implementations shoud be thread-safe
 *  and able to be seamlessly called either from serial
 *  or multi-threaded mode without any visible difference
 *  to the user.
 *
 *  Presently existing implementations are RootPersistenceManager class.
 *
 *  \attention Thread safety assumes that the class is used in a thread-local
 *  part of the code. An ideal place is a custom G4Run class (or perhaps an event
 *  or a hit). For an example, see the texansim::Run class.
 */
class VPersistenceManager {
public:
  /// ID for type of persistent object (e.g. histogram of intergers, floats, etc)
	enum Type_t { kDouble = 0, kFloat = 1, kInt = 2 };

public:
	/// Ctor
	/*! Initialize default messenger */
	VPersistenceManager();

	/// Dtor
	/*! Delete default messenger */
	virtual ~VPersistenceManager();

	/// Set output file name
	/*! Argument should be the desired file name _without_ an extensions.
	 *  Implementation should tack on the correct extension (e.g. .root).
	 *  For worker threads, append _t0 _t1, etc to the name, just before 
	 *  the extension.
	 */
	virtual void SetFilename(const G4String& name) = 0;

	/// Get output file name
	/*! Returns the output filename _after_ thread-local IDs and
	 *  extensions are appended.
	 */
	virtual const G4String& GetFilename() const = 0;


	/// Open the output file
	/*! Open the output file (thread-local) and prepare it for writing.
	 */
	virtual G4bool OpenFile() = 0;


	/// Open the output file with a specific name
	/*! Open the output file (thread-local) and prepare it for writing.
	 */
	G4bool OpenFile(const G4String& name)
		{ SetFilename(name) ; return OpenFile(); }


	/// Add class object to be saved
	/*! Tells the class to serialize a given instance of an object.
	 * Every time SaveEvent() gets called, the current state of the
	 * object should be written to disk.
	 *
	 * \param name Desired name of the serialized object
	 * \param classname Text specifying the class name of the serialized object
	 * \param addr Pointer to the object instance you wish to serialize.
	 * \attention The object at _addr_ must live at least as long as this class
	 */
	virtual G4bool AddObject(const G4String& name, const G4String& classname, void* addr) = 0;


	/// Add primitive type (double, float, int) to be saved
	/*! Allows writing to disk of a basic type.
	 * \param name Name to be associated with the data in the output file
	 * (e.g. like a Ntuple row name).
	 * \param addr Address of the primitive to store. Must point to something
	 * that will last as long as this class. The value pointed to gets saved
	 * to disk upon each call of SaveEvent().
	 * \param type Code telling what type this is. Currently G4double,
	 * G4float, G4int ar supported.
	 * \warning The actual type at the address _addr_ should match that specified
	 *  by _type_. No checking is done to make sure this is the case.
	 */
	virtual G4bool AddPrimitive(const G4String& name, void* addr, Type_t type) = 0;


	/// Add 1d histogram
	/*! Add a 1-d histogram to be saved object. At each call of SaveEvent(),
	 *  the histogram gets filled.
	 *  \param name Histogram name
	 *  \param title Histogram title
	 *  \param bins Number of bins
	 *  \param xlow Low edge of the axis
	 *  \param xhigh High edge of the axis
	 *  \param valuePointer Address of the value to store in the histogram.
	 *   The current value of the object at this address gets filled into the
	 *   histogram at each call to SaveEvent(). The pointed-to data
	 *   must stick around as long as this class does.
	 *  \param type Code indexing the type of histogram (G4int, G4float, G4double).
	 *  \warning The actual type at the address _valuePointer_ should match that specified
	 *  by _type_. No checking is done to make sure this is the case.
	 */
	virtual G4bool AddHistogram1d(
		const G4String& name, const G4String& title,
		G4int bins, G4double xlow, G4double xhigh,
		void* valuePointer, Type_t type) = 0;


	/// Add 2d histogram
	/*! Add a 2-d histogram to be saved object. At each call of SaveEvent(),
	 *  the histogram gets filled.
	 *  \param name Histogram name
	 *  \param title Histogram title
	 *  \param xbins Number of x-axis bins
	 *  \param xlow Low edge of the axis
	 *  \param xhigh High edge of the axis
	 *  \param ybins Number of y-axis bins
	 *  \param ylow Low edge of the y-axis
	 *  \param yhigh High edge of the y-axis
	 *  \param valuePointerx Address of the value to store in the histogram's x-axis.
	 *   The current value of the object at this address gets filled into the
	 *   histogram at each call to SaveEvent(). The pointed-to data
	 *   must stick around as long as this class does.
	 *  \param valuePointery Address of the value to store in the histogram's y-axis.
	 *  \param type Code indexing the type of histogram (G4int, G4float, G4double).
	 *  \warning The actual type at the address _valuePointer_ should match that specified
	 *  by _type_. No checking is done to make sure this is the case.
	 */
	virtual G4bool AddHistogram2d(
		const G4String& name, const G4String& title,
		G4int xbins, G4double xlow, G4double xhigh,
		G4int ybins, G4double ylow, G4double yhigh,
		void* valuePointerx, void* valuePointery,
		Type_t type) = 0;

	/// Save an event of histograms and persistent objects
	/*! Should do the following:
		- Write the current state of all objects registered with AddObject()
		- Write the current value of all basic data registered with AddPrimitive()
		- Fill registered histograms with the present value(s) at their specified
		valuePointer addresses
	*/
	virtual void SaveEvent() = 0;


	/// Write to disk
	/*! Performs final flush of all registered classes, primitive data,
		  and histograms to disk.
	 */
	virtual void Write() = 0;


	/// Close everything out
	/*! Take care of all resource management, freeing files, etc. here.
	 *  After calling this, the destructor should be redundant.
	 */ 
	virtual void Close() = 0;

	
	/// Get pointer to UI messenger
	G4UImessenger* GetMessenger() const;


	/// Set custom UI messenger
	/*! Call this function to use your own derived G4UIMessenger class.
	 *  Default is to use texansim::PersistenceMessenger
	 */
	void SetMessenger(G4UImessenger* messenger);


// protected:
// 	static G4String& GetFilenameBase();
	

private:
	/// UI messenger instance
	G4UImessenger* fMessenger;


	// friend class G4UIMessenger;
};

}



#endif
