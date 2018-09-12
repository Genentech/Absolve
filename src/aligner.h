
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//  The aligner object handles the alignment of query sequences to the 
//  germline databases.  It abstracts the interface to the underlying
//  BWA and SSW aligners.  It expects fasta formatted files as initializers.
//	There is one instance of this class for each classifier, typically 
//	Constant, V and J for both heavy and light.
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include <zlib.h>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <ctype.h>
#include <math.h>
#include <string>
#include <map>


#include "bwaaligner.h"
#include "sswaligner.h"

class Aligner {
public:

	Aligner();
	~Aligner();

	// idxfiles are the fasta files for a given database
	Aligner(const char *idxFile, const char* name);

	int Initialize(const char *idxFile, const char * name,
		BWAalignConfig & config);

	int Initialize(const char *idxFile, const char * name,
		int a = 1, int b = 4, 
		int o_ins = 6, int e_ins = 1,
		int o_del = 6, int e_del = 1,
		int pen_clip3 = 5, int pen_clip5 = 5,
		int min_seed_len = 19,
		int T = 30
	);


	// for a given germline call and query sequence,
	// return the number of mutations and their positions
	int ComputeSHM(std::string call, 
		const char *seq, int seqlen, 
		std::string & cigar,
		std::string & shm,
		std::vector<int> & mutpos
		);

	// Get the best <multihit> hits for a given query sequence
	// return calls, scores, pos, and cigars.
	int RunClassifier(const char *seq, int seqlen,
			std::string & call, std::string & score,
			std::string & pos, std::string & cigar,
			int firstpos = 0, // constrain alignment to firstpos 
			int lastpos = -1 // and lastpos, -1 -> use whole seq
		);

	// return top multihit results, not just the top hit.
	int SetMultiHit(int multihit);

	int m_useSSW;
	int m_byAA;


private:
	int Ready();

	// both types of aligners, each used depending on configuration
	BWAAligner m_bwaAligner;
	SSWAligner m_sswAligner;

	// The germline sequences
	BWAReader m_library;
	std::map<std::string, int> m_libraryEntries;
	

};


 
