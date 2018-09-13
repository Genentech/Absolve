
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: MIT, see included COPYING file
//
//
//	The BWAAligner objects encapsulates the use of the BWA MEM algorithm 
//	to perform only DNA alignment. 
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


#include "bwareader.h"

typedef struct {
	int A;
	int B;
	int OINS;
	int EINS;
	int ODEL;
	int EDEL;
	int CLIP3;
	int CLIP5;
	int MIN_SEED_LEN;
	int MIN_SCORE;
} BWAalignConfig;


class BWAAligner {
public:

	BWAAligner();
	BWAAligner(const char *idxFile, const char* name);
	~BWAAligner();

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
	
	int FindBestAlignment(const char *seq, int seqlen,
		std::string & rname, std::string & score, 
		std::string & pos, std::string & cigar,
		int firstpos = 0,
		int lastpost = -1
		);

	int SetMultiHit(int multihit);


public:
	int m_useSSW;
	int m_byAA;

private:
	int m_multihit;

	int Classify(bseq1_t *tseq,
		std::string & call, std::string & pos, std::string & cigar,
		std::string & score);
	int GetSAMDetails(char *sam,
		std::string & rname, std::string & pos, std::string & cigar,
		std::string & score);

	bseq1_t * CreateBseq(const char *seq, int len);

	int FreeBWASeq(bseq1_t **seq);

	std::string m_name;
	mem_opt_t * m_opt;
	bwaidx_t * m_idx;

	
};


 
