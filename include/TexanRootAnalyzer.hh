/// \file RootAnalyzer.hh
/// \brief Defines ROOT analysis class.


class TFile;
class TTree;

namespace texansim {

/// Singleton ROOT analysis class
class RootAnalyzer
{
public:
	/// dtor
	virtual ~RootAnalyzer();
	/// Use to obtain singleton instance
	static RootAnalyzer*& Instance();

	/// Open ROOT file
	bool OpenFile(const char* filename, const char* mode);
	/// Close ROOT file
	void CloseFile();
	/// Get pointer to opened file
	TFile* GetFile() const { return fFile; }

	/// Create a TTree
	TTree* CreateTree(const char* name, const char* title);
	/// Get pointer to TTree
	TTree* GetTree(const char* name);
	
protected:
	RootAnalyzer(): fFile(NULL) { }

private:
	TFile* fFile;
};

}
