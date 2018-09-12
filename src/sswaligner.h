
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//	The SSWAligner objects encapsulates the use of the SSW algorithm
//	to perform DNA or AA alignment. 
//
//	See http://journals.plos.org/plosone/article?id=10.1371/journal.pone.0082138
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include "ssw_cpp.h"
#include "bwareader.h"

class SSWAligner {
public:
	SSWAligner();
	~SSWAligner();

	// perform a single alignment
	int Align(const char *query, int qlen,
		const char *ref, int rlen,
		std::string & cigar,
		int & rpos,
		double & score,
		int & mismatches,
		int byAA);
 
  
  	// find the best alignment in the library
	int FindBestAlignment( const char *seq, int seqlen,
		std::string & call, std::string & score, 
		std::string & pos, std::string & cigar, int byAA,
		BWAReader & library, int firstpos=0, int lastpos=-1);

	int SetMultiHit(int multihit);

private:	
	int m_multihit;
	StripedSmithWaterman::Aligner *m_NTaligner;
	StripedSmithWaterman::Aligner *m_AAaligner;
	
	// Declares a default filter
	StripedSmithWaterman::Filter m_filter;
	// Declares an alignment that stores the result
	StripedSmithWaterman::Alignment m_alignment;



};



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// 
