
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: MIT, see included COPYING file
//
//
//	the IgHMM object encapsulates the use of an HMM to perform:
//		* alignment of a single DNA or AA sequence to the HMM
//		* alignment of a paired-end reads to the HMM
//	
//	The primary outputs of the "alignment" are:
//		* a log odds score
//		* orientation and reading frame
//		* FW and CDR annotations
//		* Kabat numbering
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <string>
#include <vector>

/////////////////////////////////////////////////////////
extern "C" {

	#define new my_new
	#include<p7_config.h>
	#include<hmmer.h>
	#include<easel.h>
	#include<esl_alphabet.h>
	#include<esl_getopts.h>
	#include<esl_msa.h>
	#include<esl_msafile.h>
	#include<esl_sq.h>
	#include<esl_sqio.h>
	#include<esl_vectorops.h>
	#include<esl_translate.h>
	#undef new
}

/////////////////////////////////////////////////////////
#include"util.h"
/////////////////////////////////////////////////////////

enum HMMType { Heavy = 0, Light = 1, ScFV = 2 };
extern const char *HMMTypeS[3];

/////////////////////////////////////////////////////////

class Domain {
public:
	std::string str;
	int AAstart;
	int AAstop;
	int DNAstart;
	int DNAstop;
	void Clear(){ str.clear(); AAstop = AAstart = DNAstart = DNAstop = -1; }
};




/////////////////////////////////////////////////////////
class IgHMM {
public:
	IgHMM();
	int Initialize(const char *hmmfile, int NT_PAIR_MIN);
	~IgHMM();


	// Align a single amino acid sequence
	int AlignSingleAA(std::string & aaseq);
	// Align a single DNA sequence
	int AlignSingle(char *r1, int r1len);
	// Align a paired end read
	int AlignPaired(char *r1, int r1len, char *r2, int r2len,
		int & mergedLen, int & mergedFrame);

////////////////////////////////////////
private:
	int AlignPairedByNT(char *r1, int r1len, char *r2, int r2len,
		char * mergedReads, int & mergedLen);
	
	int AlignPairedByHMM(char *r1, int r1len, char *r2, int r2len,
		char * mergedReads, int & mergedLen, int & mergedFrame);

	int ClearStrings();

	int GetDNAforAA(const char *dna, int dnalen,
		const char *aa, int aalen,
		int frame,
		std::string & subdna);

	int TraceSeq(char *seq, int seqlen,
		P7_TRACE *trace, float *tscore);

	int PairTraces(
			ESL_SQ * r1sq, P7_TRACE *tr1, ESL_SQ *s1, int frame1,
	        ESL_SQ * r2sq, P7_TRACE *tr2, ESL_SQ *s2, int frame2,
			char * mergedReads, int & mergedLen, int & mergedFrame);

	int MergeByMatchPos(int r1NTPos, int r2NTPos,
		char * r1Seq, int r1N,
		char * r2Seq, int r2N,
		char *mergedReads,
		int & mergedLen,
		int & startRead);



	int FindBestORF(ESL_SQ & rsq, 
		int & bestFrame1, int & bestRC1,  float & bestScore1,
		P7_TRACE *trace[6],
		ESL_SQ *raa
		);
	
	int GenHeavyKabat();
	int GenLightKabat();
	int TraceDomains(P7_TRACE *tr, ESL_SQ *sq);
	int TraceHMM(int hoffset, int loffset, int doH, int doL);

	int MapResiduePositions2Kabat();

	int GetMSA(ESL_SQ *seq, P7_TRACE *tr, std::string & msastr);
	int GetMSA2(ESL_SQ *s1, P7_TRACE *tr1, ESL_SQ *s2, P7_TRACE *tr2,
		std::string & msa1, std::string & msa2);

	int STDoutMSA();

	int STDoutMSA(ESL_MSA *somemsa,float somescore);
	int STDoutMSA2(ESL_MSA *somemsa,float score1, float score2);

	int OneLineMSA2(ESL_MSA *r1r2msa, float score1, float score2);


private:


	int m_NT_PAIR_MIN;

	P7_HMM * m_hmm;
	P7_TRACE * m_tr;
	ESL_SQ * m_sq;
	P7_PROFILE * m_gm;
	ESL_ALPHABET * m_abc;

	ESL_ALPHABET * m_dnaabc;

	ESL_MSA * m_msa;

	P7_TRACE * m_r1Traces[6];
	P7_TRACE * m_r2Traces[6];

	ESL_SQ * m_orfSearchAA;

	ESL_SQ * m_r1aa;
	ESL_SQ * m_r2aa;

	// some work buffers
	char *m_aabuf;
	char *m_dnabuf;

	// data for aligning traces
	std::string m_R1prof;
	std::string m_R2prof;
	std::vector<int> m_R1profInd;
	std::vector<int> m_R2profInd;
	int m_R1matchLeft;
	int m_R1matchRight;
	int m_R2matchLeft;
	int m_R2matchRight;


	std::string m_MSA;
    std::string m_Hkabatstr;
    std::string m_Lkabatstr;

public:

	///////////////////////////////	
	HMMType m_hmmType;
	float  m_score;

	char * m_mergedPairs;
	char * m_mergedAA;
	std::string m_alignMethod;
	int m_bestFrame1;
	int m_bestFrame2;

	// domain sequences
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

	std::string m_HAA;
	std::string m_LAA;
	std::string m_HDNA;
	int m_HDNAstart;
	std::string m_LDNA;
	int m_LDNAstart;

	//Map of aa residue position to its assign kabat position
	// (useful for SHM tally)
	std::vector<std::string> m_Hres2kabat;
	std::vector<std::string> m_Lres2kabat;

	// kabat strings
    std::string m_kabatstr;


};

/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// 
