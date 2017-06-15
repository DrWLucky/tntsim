#ifndef TNT_INPUT_FILE_PARSER_HEADER_FILE_12345
#define TNT_INPUT_FILE_PARSER_HEADER_FILE_12345
#include <map>
#include <string>
#include <cstring>
#include <sstream>
#include <cstring>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include "globals.hh"


class TntKeyConverterBase {
public:
	virtual void Convert(const std::string&) = 0;
	virtual ~TntKeyConverterBase() { }
};

template<class T, class ParameterSetter_t> 
class TntKeyConverter : public TntKeyConverterBase {
public:
	TntKeyConverter(ParameterSetter_t* SetterClassInstance,
									void (ParameterSetter_t::*setter)(T) ):
		mInstance(SetterClassInstance),
		mGlobalSetter(setter) { }

	virtual ~TntKeyConverter() 
		{  }

	void Convert(const std::string& str)
		{
			std::stringstream sstr(str);
			T t; sstr >> t;
			mGlobalSetter(mInstance, t);
		}
	
private:
	ParameterSetter_t *mInstance;
	std::mem_fun1_t<void, ParameterSetter_t, T> mGlobalSetter;
};


template<class ParameterSetter_t>
class TntKeyConverter<std::string, ParameterSetter_t> : public TntKeyConverterBase {
public:
	TntKeyConverter( ParameterSetter_t* SetterClassInstance,
									 void (ParameterSetter_t::*setter)(G4String) ):
		mInstance(SetterClassInstance),
		mGlobalSetter(setter) { }
	virtual ~TntKeyConverter() { }
	void Convert(const std::string& str)
	{
		mGlobalSetter(mInstance, str.c_str());
	}
	
private:
	ParameterSetter_t *mInstance;
	std::mem_fun1_t<void, ParameterSetter_t, G4String> mGlobalSetter;
};

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
					irange_t range = mInputs.equal_range(*itTokens++);
					for(imap_t::iterator it = range.first; it != range.second; ++it) {
						it->second->Convert(*itTokens++);
					}
				}
			}
		}

	template<class T>
	void AddInput(const std::string& key, void (ParameterSetter_t::*setter) (T))
		{
			TntKeyConverter<T, ParameterSetter_t> *keyConverter = 
				new TntKeyConverter<T, ParameterSetter_t> (mSetter, setter);
			mInputs.insert(std::make_pair(key, keyConverter));
		}
	
	
private:
	void tab_to_space(std::string& s)
		{
			for(std::string::iterator it = s.begin(); it != s.end(); ++it) {
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
