/// \file TexanRootAnalyzer.hh
/// \brief Defines ROOT analysis class.

class TFile;
class TTree;

/// Singleton ROOT analysis class
class TexanRootAnalyzer
{
public:
	enum { kNever = 0, kPrompt = 1, kAlways = 2 };

	virtual ~TexanRootAnalyzer();
	static TexanRootAnalyzer*& Instance();

	bool OpenFile(const char* filename, int overwrite = kPrompt);
	void CloseFile();
	TFile* GetFile() const { return fFile; }

	TTree* CreateTree(const char* name, const char* title);
	TTree* GetTree(const char* name);
	
protected:
	TexanRootAnalyzer(): fFile(0) { }

private:
	TFile* fFile;
};
