#include <string>
#include <fstream>
#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include "TexanRootAnalyzer.hh"

template <class T>
void Zap(T*& t)
{
	if(t) {
		delete t;
		t = NULL;
	}
}


TexanRootAnalyzer::~TexanRootAnalyzer()
{
	Zap(fFile);
}

TexanRootAnalyzer*& TexanRootAnalyzer::Instance()
{
	static TexanRootAnalyzer* singleton = 0;
	if(!singleton) {
		singleton = new TexanRootAnalyzer();
	}
	return singleton;
}

bool TexanRootAnalyzer::OpenFile(const char* filename, int overwrite)
{
	{
		std::ifstream ifs(filename);
		if(ifs.good()) {
			switch(overwrite) {
			case kNever:
				return false;
			case kPrompt:
				{
					std::string input;
					std::cout << "Overwrite existing root file? (y/[n])\n";
					std::cin  >> input;
					if( !(input == "y" || input == "Y") )
						return false;
					break;
				}
			case kAlways:
				break;
			default:
				break;
			}
		}
	}

	Zap(fFile);
	fFile = new TFile(filename, "RECREATE");
	return (fFile && !(fFile->IsZombie()));
}

void TexanRootAnalyzer::CloseFile()
{
	if(fFile)
		fFile->Close();
}

TTree* TexanRootAnalyzer::CreateTree(const char* name, const char* title)
{
	if(fFile)
		fFile->cd();
	else
		return NULL;

	TTree* t = new TTree(name, title);
	return t;
}

TTree* TexanRootAnalyzer::GetTree(const char* name)
{
	return fFile ? dynamic_cast<TTree*>(fFile->FindObject(name)) : NULL;
}
