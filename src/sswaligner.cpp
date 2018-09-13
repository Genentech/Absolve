
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: MIT, see included COPYING file
//
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <string>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <boost/algorithm/string.hpp>


#include "sswaligner.h"



////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
static const int8_t ssw_aa_mat50[] = {
	//  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V   B   Z   X   *
     	5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -2, -1, -1, -3, -1,  1,  0, -3, -2,  0, -2, -1, -1, -5,	// A
       -2,  7, -1, -2, -4,  1,  0, -3,  0, -4, -3,  3, -2, -3, -3, -1, -1, -3, -1, -3, -1,  0, -1, -5,	// R
       -1, -1,  7,  2, -2,  0,  0,  0,  1, -3, -4,  0, -2, -4, -2,  1,  0, -4, -2, -3,  5,  0, -1, -5,	// N
       -2, -2,  2,  8, -4,  0,  2, -1, -1, -4, -4, -1, -4, -5, -1,  0, -1, -5, -3, -4,  6,  1, -1, -5,	// D
       -1, -4, -2, -4, 13, -3, -3, -3, -3, -2, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -3, -3, -1, -5,	// C
       -1,  1,  0,  0, -3,  7,  2, -2,  1, -3, -2,  2,  0, -4, -1,  0, -1, -1, -1, -3,  0,  4, -1, -5,	// Q
       -1,  0,  0,  2, -3,  2,  6, -3,  0, -4, -3,  1, -2, -3, -1, -1, -1, -3, -2, -3,  1,  5, -1, -5,	// E
     	0, -3,  0, -1, -3, -2, -3,  8, -2, -4, -4, -2, -3, -4, -2,  0, -2, -3, -3, -4, -1, -2, -1, -5,	// G
       -2,  0,  1, -1, -3,  1,  0, -2, 10, -4, -3,  0, -1, -1, -2, -1, -2, -3,  2, -4,  0,  0, -1, -5,	// H
       -1, -4, -3, -4, -2, -3, -4, -4, -4,  5,  2, -3,  2,  0, -3, -3, -1, -3, -1,  4, -4, -3, -1, -5,	// I
       -2, -3, -4, -4, -2, -2, -3, -4, -3,  2,  5, -3,  3,  1, -4, -3, -1, -2, -1,  1, -4, -3, -1, -5,	// L
       -1,  3,  0, -1, -3,  2,  1, -2,  0, -3, -3,  6, -2, -4, -1,  0, -1, -3, -2, -3,  0,  1, -1, -5,	// K
       -1, -2, -2, -4, -2,  0, -2, -3, -1,  2,  3, -2,  7,  0, -3, -2, -1, -1,  0,  1, -3, -1, -1, -5,	// M
       -3, -3, -4, -5, -2, -4, -3, -4, -1,  0,  1, -4,  0,  8, -4, -3, -2,  1,  4, -1, -4, -4, -1, -5,	// F
       -1, -3, -2, -1, -4, -1, -1, -2, -2, -3, -4, -1, -3, -4, 10, -1, -1, -4, -3, -3, -2, -1, -1, -5,	// P
     	1, -1,  1,  0, -1,  0, -1,  0, -1, -3, -3,  0, -2, -3, -1,  5,  2, -4, -2, -2,  0,  0, -1, -5,	// S
    	0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  2,  5, -3, -2,  0,  0, -1, -1, -5, 	// T
       -3, -3, -4, -5, -5, -1, -3, -3, -3, -3, -2, -3, -1,  1, -4, -4, -3, 15,  2, -3, -5, -2, -1, -5, 	// W
       -2, -1, -2, -3, -3, -1, -2, -3,  2, -1, -1, -2,  0,  4, -3, -2, -2,  2,  8, -1, -3, -2, -1, -5, 	// Y
     	0, -3, -3, -4, -1, -3, -3, -4, -4,  4,  1, -3,  1, -1, -3, -2,  0, -3, -1,  5, -3, -3, -1, -5, 	// V
       -2, -1,  5,  6, -3,  0,  1, -1,  0, -4, -4,  0, -3, -4, -2,  0,  0, -5, -3, -3,  6,  1, -1, -5, 	// B
       -1,  0,  0,  1, -3,  4,  5, -2,  0, -3, -3,  1, -1, -4, -1,  0, -1, -2, -2, -3,  1,  5, -1, -5, 	// Z
       -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -5, 	// X
       -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1 	// *
	};

static int8_t ssw_aa_table[128] = {
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23, 23,
	23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
	14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23,
	23, 0,  20, 4,  3,  6,  13, 7,  8,  9,  23, 11, 10, 12, 2,  23,
	14, 5,  1,  15, 16, 23, 19, 17, 22, 18, 21, 23, 23, 23, 23, 23
};

////////////////////////////////////////////////////////////////////////////

SSWAligner::SSWAligner(){

	m_multihit = 0;

	m_NTaligner = new StripedSmithWaterman::Aligner;
	m_AAaligner = new StripedSmithWaterman::Aligner(
		ssw_aa_mat50, 24, ssw_aa_table, 128);
}

SSWAligner::~SSWAligner(){

	delete m_NTaligner;
	delete m_AAaligner;
}

int SSWAligner::Align(const char *query, int qlen,
	const char *ref, int rlen,
	std::string & cigar,
	int & rpos,
	double & score,
	int & mismatches,
	int byAA){

	if(byAA){
		m_AAaligner->Align(query, ref, rlen, m_filter, &m_alignment);
	}else{
		m_NTaligner->Align(query, ref, rlen, m_filter, &m_alignment);
	}

	score = m_alignment.sw_score;
	mismatches = m_alignment.mismatches;
	//use BAM style cigars with M instead of X and =, more compact:
	cigar = m_alignment.dense_cigar_string;
	
	rpos = m_alignment.ref_begin + 1; // ssw is 0 based

	return 0;
}
////////////////////////////////////////////////////////////////////////////
int SSWAligner::SetMultiHit(int multihit){
	m_multihit = multihit;
	return 0;
}
////////////////////////////////////////////////////////////////////////////
//
typedef struct {
	double score;
	std::string call;
	std::string cigar;
	int pos;
	int mismatches;
} OneHit;


double InsertHit(std::list<OneHit> & allhits,
	int getHits, OneHit & ahit){

	//Where to insert 
	std::list<OneHit>::iterator it = allhits.begin();
	
	if(allhits.size() == 0){
		allhits.insert(it, ahit);
		return ahit.score;
	}
	for(; it != allhits.end(); ++it){
		if(ahit.score > it->score){
			allhits.insert(it, ahit);
			if(allhits.size() > (unsigned int)getHits){ 
				allhits.pop_back(); 
			}
			break;
		}
	}

	return allhits.back().score;
}


int SSWAligner::FindBestAlignment( const char *seq, int seqlen,
	std::string & call, std::string & score, std::string & pos,
	std::string & cigar, int byAA, 
	BWAReader & library, int firstpos, int lastpos ){


	OneHit ahit;
	unsigned int getHits = 1+m_multihit;
	std::list<OneHit> allhits;
	double minScore = 0.0;

	int sl = seqlen;
	if(lastpos != -1){
		sl = lastpos - firstpos + 1;
	}
	for(int i=0; i<library.m_Nbseqs; i++){
		Align(
			seq + firstpos, 
			sl,
			library.m_bseqs[i].seq, 
			library.m_bseqs[i].l_seq,
			ahit.cigar, ahit.pos, ahit.score, ahit.mismatches, byAA);
		if(allhits.size() < getHits || ahit.score > minScore){
			ahit.call = library.m_bseqs[i].name;
			minScore = InsertHit(allhits, getHits, ahit);
		}
	}
	//Assign top hit and concatenate others
	std::list<OneHit>::iterator it = allhits.begin();
	call = it->call;
	cigar = it->cigar;
	score = boost::lexical_cast<std::string>(it->score);
	pos = boost::lexical_cast<std::string>(it->pos);
	//append the rest
	++it;
	for(; it != allhits.end(); ++it){
		call = call + "," + it->call;
		pos = pos + "," + boost::lexical_cast<std::string>(it->pos);
		cigar = cigar + "," + it->cigar;
		score = score + "," + boost::lexical_cast<std::string>(it->score);
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// 
