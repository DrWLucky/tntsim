#include <string>
#include <fstream>
#include <iostream>

#include <TFile.h>
#include <TTree.h>

#include "TexanRootAnalyzer.hh"
#include "TexanUtils.hh"



texansim::RootAnalyzer::~RootAnalyzer()
{
	Zap(fFile);
}

texansim::RootAnalyzer*& texansim::RootAnalyzer::Instance()
{
	static texansim::RootAnalyzer* singleton = 0;
	if(!singleton) {
		singleton = new texansim::RootAnalyzer();
	}
	return singleton;
}

bool texansim::RootAnalyzer::OpenFile(const char* filename, const char* mode)
{
	Zap(fFile);
	fFile = new TFile(filename, mode);
	return (fFile && !(fFile->IsZombie()));
}

void texansim::RootAnalyzer::CloseFile()
{
	if(fFile) fFile->Close();
}

TTree* texansim::RootAnalyzer::CreateTree(const char* name, const char* title)
{
	if(fFile)
		fFile->cd();
	else
		return NULL;

	TTree* t = new TTree(name, title);
	return t;
}

TTree* texansim::RootAnalyzer::GetTree(const char* name)
{
	return fFile ? dynamic_cast<TTree*>(fFile->FindObject(name)) : NULL;
}
