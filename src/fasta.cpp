
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <boost/regex.hpp>

#include "fasta.h"

////////////////////////////////////////////////////////////////////////////
FastaIO::FastaIO(){
	
}

////////////////////////////////////////////////////////////////////////////
FastaIO::~FastaIO(){
	if(m_input.is_open()) m_input.close();
}

////////////////////////////////////////////////////////////////////////////
int FastaIO::OpenFile(const char *filename){
	m_inputfile = filename;
    m_input.open(m_inputfile.c_str(), std::ios_base::in);

	return(!m_input.good());
}
////////////////////////////////////////////////////////////////////////////
int FastaIO::CloseFile(){
	if(m_input.is_open()) m_input.close();
	return 0;
}

////////////////////////////////////////////////////////////////////////////

int FastaIO::SkipFirstNSeqs(int nseqs){
	for(int i=0;i<nseqs; i++){
		std::vector<std::string>  readNames, readContents;
		ReadNSeqs(readNames,readContents,1);
		if(readNames.size() == 0) break;
	}
	return 0;
}
////////////////////////////////////////////////////////////////////////////
int FastaIO::ReadNSeqs(
	std::vector<std::string> & readNames,
	std::vector<std::string> & readContents,
	int getN){
	
	readNames.clear();
	readContents.clear();

	if(!m_input.is_open()){ return -1;}

	std::string line, name, content;

	boost::regex whitespace("[[:space:]]");
	boost::regex tabsuffix("\\t.*");
	
	for(int i=0;i<getN; i++){
		if(m_input.peek() != '>'){
			return -1;
		}
		//Get the header line
		std::getline( m_input, line);
		name = line.substr(1); //Trim the leading ">"
		// trim any \t suffix
		name = boost::regex_replace(name, tabsuffix, "");
		readNames.push_back(name);

		//Now the sequence
		while(m_input.peek() != '>'){
			line ="";
			std::getline( m_input, line);
			if(! line.empty() ){ 
				line = boost::regex_replace(line,whitespace,"",
					boost::match_default | boost::format_all);
				content += line;
			}else{
				break;
			}
		}
		readContents.push_back(content);

		if((i+1) >= getN){ return 0;}
	}
 
	return 0;
}
//////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////// 
