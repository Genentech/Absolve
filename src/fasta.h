
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//	The FastaIO object provides functions for reading fasta formatted files
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#pragma once


#include <string>
#include <vector>
#include <iostream>
#include <fstream>
 



class FastaIO{
public:
	FastaIO();
	~FastaIO();

	int OpenFile(const char *filename);
	int CloseFile();
	int SkipFirstNSeqs(int nseqs);
	int ReadNSeqs(
		std::vector<std::string> & readNames,
		std::vector<std::string> & readContents,
		int getN);
	
private:
	
	std::string m_inputfile;
    std::ifstream m_input;

};

 
