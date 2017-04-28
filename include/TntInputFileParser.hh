#ifndef TNT_INPUT_FILE_PARSER_HEADER_FILE_12345
#define TNT_INPUT_FILE_PARSER_HEADER_FILE_12345
#include <map>
#include <string>
#include <cstring>
#include <sstream>

#include "TntGlobalParams.hh"

class TntKeyConverterBase {
public:
	virtual void Convert(const std::string&) = 0;
	virtual ~TntKeyConverterBase() { }
};

template<class T> 
class TntKeyConverter : public TntKeyConverterBase {
public:
	TntKeyConverter( void (TntGlobalParams::*setter)(T) ): 
		mGlobalSetter(setter) { }

	virtual ~TntKeyConverter() 
		{  }

	void Convert(const std::string& str)
		{
			std::stringstream sstr(str);
			T t; sstr >> t;
			mGlobalSetter(TntGlobalParams::Instance(), t);
		}
	
private:
	std::mem_fun1_t<void, TntGlobalParams, T> mGlobalSetter;
};


template<>
class TntKeyConverter<std::string> : public TntKeyConverterBase {
public:
	TntKeyConverter( void (TntGlobalParams::*setter)(G4String) ):
		mGlobalSetter(setter) { }
	virtual ~TntKeyConverter() { }
	void Convert(const std::string& str)
	{
		mGlobalSetter(TntGlobalParams::Instance(), str.c_str());
	}
	
private:
	std::mem_fun1_t<void, TntGlobalParams, G4String> mGlobalSetter;
};

class TntInputFileParser {
public:
	TntInputFileParser();
	~TntInputFileParser();
	void Parse(const std::string& filename);
	template<class T>
	void AddInput(const std::string& key, void (TntGlobalParams::*setter) (T))
		{
			mInputs.insert(std::make_pair(key, new TntKeyConverter<T>(setter)));
		}
private:
	typedef std::multimap<std::string, TntKeyConverterBase*> imap_t;
	typedef std::pair<imap_t::iterator, imap_t::iterator> irange_t;
	imap_t mInputs;
};


#endif
