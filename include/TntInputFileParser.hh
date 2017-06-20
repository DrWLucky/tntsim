#ifndef TNT_INPUT_FILE_PARSER_HEADER_FILE_12345
#define TNT_INPUT_FILE_PARSER_HEADER_FILE_12345
#include <map>
#include <string>
#include <vector>
#include <cstring>
#include <sstream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdexcept>
#include "globals.hh"


class TntKeyConverterBase {
public:
	virtual void Convert(const std::vector<std::string>&) = 0;
	virtual ~TntKeyConverterBase() { }
};


namespace {
struct cout_par { 
	void operator()(const std::string& par) { G4cerr << "\t" << par << "\n"; } 
}; }
	

template<class T, class ParameterSetter_t> 
class TntKeyConverter_1 : public TntKeyConverterBase {
public:
	TntKeyConverter_1(ParameterSetter_t* SetterClassInstance,
										void (ParameterSetter_t::*setter)(T) ):
		mInstance(SetterClassInstance),
		mSetter(setter) { }

	virtual ~TntKeyConverter_1() 
		{  }

	void Convert(const std::vector<std::string>& args)
		{
			try {
				std::stringstream sstr(args.at(0));
				T t; sstr >> t;
				mSetter(mInstance, t);
			} catch (const std::out_of_range& oor) {
				G4cerr << "ERROR:: TntKeyConverter_1:: " 
							 << "Less than one parameters supplied.\n";
				G4cerr << "Skipping setting these parameters!" << G4endl;
			}		
		}
	
private:
	ParameterSetter_t *mInstance;
	std::mem_fun1_t<void, ParameterSetter_t, T> mSetter;
};


template<class T, class T1, class ParameterSetter_t> 
class TntKeyConverter_2 : public TntKeyConverterBase {
public:
	TntKeyConverter_2(ParameterSetter_t* SetterClassInstance,
										void (ParameterSetter_t::*setter)(T, T1) ):
		mInstance(SetterClassInstance),
		mSetter(setter) { }

	virtual ~TntKeyConverter_2() 
		{  }

	void Convert(const std::vector<std::string>& args)
		{
			try {
				std::stringstream sstr(args.at(0));
				T t; sstr >> t;

				std::stringstream sstr1(args.at(1));
				T1 t1; sstr1 >> t1;
				
				(mInstance->*mSetter)(t, t1);
			} catch (const std::out_of_range& oor) {
				G4cerr << "ERROR:: TntKeyConverter_2:: " 
							 << "Less than two parameters supplied.\n"
							 << "Valid parameters are:\n";
				std::for_each(args.begin(), args.end(), cout_par());
				G4cerr << "Skipping setting these parameters!" << G4endl;
			}
		}
	
private:
	ParameterSetter_t *mInstance;
	void (ParameterSetter_t::*mSetter)(T, T1);
};


template<class T, class T1, class T2, class ParameterSetter_t> 
class TntKeyConverter_3 : public TntKeyConverterBase {
public:
	TntKeyConverter_3(ParameterSetter_t* SetterClassInstance,
										void (ParameterSetter_t::*setter)(T, T1, T2) ):
		mInstance(SetterClassInstance),
		mSetter(setter) { }

	virtual ~TntKeyConverter_3() 
		{  }

	void Convert(const std::vector<std::string>& args)
		{
			try {
				std::stringstream sstr(args.at(0));
				T t; sstr >> t;

				std::stringstream sstr1(args.at(1));
				T1 t1; sstr1 >> t1;

				std::stringstream sstr2(args.at(2));
				T2 t2; sstr2 >> t2;

				(mInstance->*mSetter)(t, t1, t2);
			} catch (const std::out_of_range& oor) {
				G4cerr << "ERROR:: TntKeyConverter_3:: " 
							 << "Less than three parameters supplied.\n"
							 << "Valid parameters are:\n";
				std::for_each(args.begin(), args.end(), cout_par());
				G4cerr << "Skipping setting these parameters!" << G4endl;
			}
		}
	
private:
	ParameterSetter_t *mInstance;
	void (ParameterSetter_t::*mSetter)(T, T1, T2);
};


#if 0
template<class ParameterSetter_t, class MemFun_t>
class TntKeyConverter<std::string, ParameterSetter_t> : public TntKeyConverterBase {
public:
	TntKeyConverter( ParameterSetter_t* SetterClassInstance,
									 void (ParameterSetter_t::*setter)(G4String) ):
		mInstance(SetterClassInstance),
		mSetter(setter) { }
	virtual ~TntKeyConverter() { }
	void Convert(const std::string& str)
	{
		mSetter(mInstance, str.c_str()); 
	}
	
private:
	ParameterSetter_t *mInstance;
	MemFun_t mSetter;
};
#endif

template<class ParameterSetter_t> class TntInputFileParser {
public:
	TntInputFileParser(ParameterSetter_t* setter):
		mSetter(setter) { }

	~TntInputFileParser()
		{
			for(imap_t::iterator it = mInputs.begin(); it != mInputs.end(); ++it) {
				delete it->second;
			}
		}

	void Parse(const std::string& filename)
		{
			std::ifstream ifs(filename.c_str());
			if(!ifs.good()) {
				G4cerr << "ERROR:: TntInputFileParser:: Bad File Name:: " << filename << G4endl;
				throw filename;
			}
			std::string line;
	
			while(std::getline(ifs, line)) {
				line = line.substr(0, line.find('#'));
				tab_to_space(line);
				std::vector<std::string> tokens = tokenize(line);
				std::vector<std::string>::iterator itTokens = tokens.begin();
				if(tokens.size() > 1) {
					imap_t::iterator it = mInputs.find(*itTokens++);
					if(it != mInputs.end()) {
						std::vector<std::string> parameters;
						while(itTokens != tokens.end()) { parameters.push_back(*itTokens++); }
						it->second->Convert(parameters);
					}
				}
			}
		}

	template<class T>
	void AddInput(const std::string& key, void (ParameterSetter_t::*setter) (T))
		{
			TntKeyConverterBase *keyConverter = 
				new TntKeyConverter_1<T, ParameterSetter_t>	(mSetter, setter);
			mInputs.insert(std::make_pair(key, keyConverter));
		}

	template<class T, class T1>
	void AddInput(const std::string& key, void (ParameterSetter_t::*setter) (T, T1))
		{
			TntKeyConverterBase *keyConverter = 
				new TntKeyConverter_2<T, T1, ParameterSetter_t> (mSetter, setter);
			mInputs.insert(std::make_pair(key, keyConverter));
		}

	template<class T, class T1, class T2>
	void AddInput(const std::string& key, void (ParameterSetter_t::*setter) (T, T1, T2))
		{
			TntKeyConverterBase *keyConverter = 
				new TntKeyConverter_3<T, T1, T2, ParameterSetter_t> (mSetter, setter);
			mInputs.insert(std::make_pair(key, keyConverter));
		}

	
private:
	void tab_to_space(std::string& str)
		{
			for(std::string::iterator it = str.begin(); it != str.end(); ++it) {
				if(*it == '\t') *it = ' ';
			}
		}

	std::vector<std::string> tokenize(const std::string& str,
																		char delim = ' ')
		{
			std::vector<std::string> output;
			std::string token = "";
			for(std::string::const_iterator it = str.begin(); 
					it != str.end(); ++it) {
				if(*it != delim) { token.push_back(*it); }
				else { if(!token.empty()) {
						output.push_back(token);
						token.clear();
					}	}
			}
			if(!token.empty()) { output.push_back(token); }
			return output;
		}	

private:
	typedef std::multimap<std::string, TntKeyConverterBase*> imap_t;
	typedef std::pair<imap_t::iterator, imap_t::iterator> irange_t;
	imap_t mInputs;
	ParameterSetter_t *mSetter;
};


#endif
