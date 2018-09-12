
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//	The result class encapsulate the result set for a single query sequence
//	and provides functions for their IO to delimited output files
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <vector>
#include <map>

#include "ighmm.h"
/////////////////////////////////////////////////////////
extern int gClassifierTypesLen;
extern const char* gClassifierTypes[];
	
/////////////////////////////////////////////////////////
class Result {
public:
	Result();
	~Result();

	void Init();
	void GetHMMResults(IgHMM & hmm);
	void AsHeader();
	void Print(FILE *out, char delim);
	void Dump();

	std::vector<std::string> m_resNames;
	std::map<std::string, std::string> m_res;

	int m_HDNAstart;
	int m_LDNAstart;

	std::vector<std::string> m_KabatPos;
	std::vector<std::string> m_KabatSHMPos;

    Domain m_HFW1;
    Domain m_HCDR1;
    Domain m_HFW2;
    Domain m_HCDR2;
    Domain m_HFW3;
    Domain m_HCDR3;
    Domain m_HFW4;

    Domain m_LFW1;
    Domain m_LCDR1;
    Domain m_LFW2;
    Domain m_LCDR2;
    Domain m_LFW3;
    Domain m_LCDR3;
    Domain m_LFW4;

};

/////////////////////////////////////////////////////////
 
