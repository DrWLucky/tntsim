#include <cstring>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <vector>
#include "TntInputFileParser.hh"


namespace {

inline void tab_to_space(std::string& s)
{
	for(std::string::iterator it = s.begin(); it != s.end(); ++it) {
		if(*it == '\t') *it = ' ';
	}
}

inline std::vector<std::string> tokenize(const std::string& str,
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

}


TntInputFileParser::TntInputFileParser() {  }

TntInputFileParser::~TntInputFileParser()
{
	for(imap_t::iterator it = mInputs.begin(); it != mInputs.end(); ++it) {
		delete it->second;
	}
}
	
void TntInputFileParser::Parse(const std::string& filename)
{
	std::ifstream ifs(filename.c_str());
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
