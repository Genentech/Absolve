
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: MIT, see included COPYING file
//
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include<fstream>
#include<utility>
#include<algorithm>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/split.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/regex.hpp>
#include <boost/lexical_cast.hpp>


//////////////////////////////////////////////////////////////////////////////
#include "util.h"
#include "ighmm.h"
#include "translate.h"

const char *HMMTypeS[3] = {
	"Heavy",
	"Light",
	"ScFV"
};


#define MAX(a,b) ((a>b)?a:b)
#define MIN(a,b) ((a<b)?a:b)

////////////////////////////////////////////////////////////////////////////

const char *LETTERS = {
	" ABCDEFGHIJKLMNOPQRSTUVWXYZ"
};

AATranslator gTranslator;

//////////////////////////////////////////////////////////////////////////////
// frameMap, convenience mapping for rc and frameshift to index into r1traces
int frameMap[2][3] = {
	// 0,1,2 (frameshift
	{0,1,2}, // rc0
	{3,4,5}  //rc1
};
//////////////////////////////////////////////////////////////////////////////
// init a HMMER sequence
void SetupSQ(char *s, int l, ESL_SQ *sq){
	memset(sq,0,sizeof(ESL_SQ));
	sq->seq = s; sq->n = l; sq->W = l; sq->L = l;
}

int DumpTrace(ESL_DSQ *dsq, P7_TRACE *tr, P7_PROFILE * gm);
void CountMatches(P7_TRACE *tr, int & numMatches );

int MapAAToNT(int ntlen, int aalen, int aapos, int frame, int rc){
	
	// rc	frame	aapos	ntpos	ntlen	aalen
	// 0	0			0		0
	// 0	0			1		3
	// 0	1			0		1
	// 0	1			1		4
	// 1	0			0		0
	// 1	0			1		3
	// 1	1			0		1
	// 1	1			1		4
	
	return rc ? (ntlen - 1 - (frame + aapos*3)) : (frame + aapos*3);
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

IgHMM::IgHMM(){
	m_hmm = NULL;
	m_tr = NULL;
	m_sq = NULL;
	m_gm = NULL;
	m_abc = NULL;
	m_dnaabc = NULL;
	m_msa = NULL;

	for(int i=0;i<6; i++){ m_r1Traces[i] = NULL; }
	for(int i=0;i<6; i++){ m_r2Traces[i] = NULL; }
	m_r1aa = NULL;
	m_r2aa = NULL;
	m_orfSearchAA = NULL;

	m_mergedPairs = NULL;
	m_mergedAA = NULL;

	m_aabuf = NULL;
	m_dnabuf = NULL;
	
	m_R1prof.clear();
	m_R2prof.clear();
	m_R1profInd.clear();
	m_R2profInd.clear();

}

//////////////////////////////////////////////////////////////////////////////
IgHMM::~IgHMM(){
	if(m_hmm != NULL) p7_hmm_Destroy(m_hmm);
	if(m_tr != NULL) p7_trace_Destroy(m_tr);
	if(m_sq != NULL) esl_sq_Destroy(m_sq);
	if(m_gm != NULL) p7_profile_Destroy(m_gm);
	if(m_abc != NULL) esl_alphabet_Destroy(m_abc);
	if(m_dnaabc != NULL) esl_alphabet_Destroy(m_dnaabc);
	if(m_msa != NULL) esl_msa_Destroy(m_msa);
	
	for(int i=0;i<6; i++){ 
		if(m_r1Traces[i] != NULL) p7_trace_Destroy(m_r1Traces[i]); 
	}
	for(int i=0;i<6; i++){
		if(m_r2Traces[i] != NULL) p7_trace_Destroy(m_r2Traces[i]); 
	}
	if(m_r1aa != NULL) esl_sq_Destroy(m_r1aa);
	if(m_r2aa != NULL) esl_sq_Destroy(m_r2aa);
	if(m_orfSearchAA != NULL) esl_sq_Destroy(m_orfSearchAA);

	if(m_mergedPairs != NULL) free(m_mergedPairs);
	if(m_mergedAA != NULL) free(m_mergedAA);

	if(m_aabuf != NULL) free(m_aabuf);
	if(m_dnabuf != NULL) free(m_dnabuf);
	
	m_R1prof.clear();
	m_R2prof.clear();
	m_R1profInd.clear();
	m_R2profInd.clear();

}
//////////////////////////////////////////////////////////////////////////////
// Create a digital HMMER sequence
// By default CreateDigital only gives us dsq of size 256, make it bigger
ESL_SQ * SQCreateBigDigital(ESL_ALPHABET *abc, int n){
	ESL_SQ *sq = esl_sq_CreateDigital(abc);
	free(sq->dsq);
	sq->salloc = n;
	sq->dsq = (ESL_DSQ *) malloc(sq->salloc * sizeof(ESL_DSQ));
	return sq;
}
//////////////////////////////////////////////////////////////////////////////
int IgHMM::Initialize(const char *hmmfile, int NT_PAIR_MIN){

	m_NT_PAIR_MIN = NT_PAIR_MIN;

	char errbuf[eslERRBUFSIZE];
	P7_HMMFILE *hfp = NULL;

	char buf[8192];
	sprintf(buf,"%s",hmmfile);

	int status;
	status = p7_hmmfile_OpenE(buf, NULL, &hfp, errbuf);
	if      (status == eslENOTFOUND) p7_Fail((char*)"File existence/permissions problem in trying to open HMM file %s.\n%s\n", hmmfile, errbuf);
	else if (status == eslEFORMAT)   p7_Fail((char*)"File format problem in trying to open HMM file %s.\n%s\n",                hmmfile, errbuf);
	else if (status != eslOK)        p7_Fail((char*)"Unexpected error %d in opening HMM file %s.\n%s\n",                       status, hmmfile, errbuf);

	status = p7_hmmfile_Read(hfp, &m_abc, &m_hmm);
	if (status == eslEFORMAT)   p7_Fail((char*)"Bad file format in HMM file %s:\n%s\n",          hfp->fname, hfp->errbuf);
	else if (status == eslEINCOMPAT) p7_Fail((char*)"HMM in %s is not in the expected %s alphabet\n", hfp->fname, esl_abc_DecodeType(m_abc->type));
	else if (status == eslEOF)       p7_Fail((char*)"Empty HMM file %s? No HMM data found.\n",        hfp->fname);
	else if (status != eslOK)        p7_Fail((char*)"Unexpected error in reading HMMs from %s\n",     hfp->fname);

	status = p7_hmmfile_Read(hfp, &m_abc, NULL);
	if      (status != eslEOF)       p7_Fail((char*)"HMM file %s does not contain just one HMM\n",    hfp->fname);
	p7_hmmfile_Close(hfp);

	m_tr = p7_trace_CreateWithPP();
	m_sq = SQCreateBigDigital(m_abc, ABS_MAX_AA_LEN);

	// impl_Init forces faster handling of denormal numbers
	// leads to some over/under flow problems for weird cases,
	// but much faster so keep it.
	impl_Init();
	p7_FLogsumInit();

	P7_BG *bg;
	bg = p7_bg_Create(m_abc);

	m_gm = NULL;
	m_gm = p7_profile_Create (m_hmm->M, m_abc);

	p7_ProfileConfig(m_hmm, bg, m_gm, 100, p7_LOCAL);
	p7_bg_Destroy(bg);

	m_dnaabc = esl_alphabet_Create(eslRNA);
	m_msa = NULL;

	for(int i=0;i<6; i++){ m_r1Traces[i] = p7_trace_CreateWithPP(); }
	for(int i=0;i<6; i++){ m_r2Traces[i] = p7_trace_CreateWithPP(); }
	m_r1aa = SQCreateBigDigital(m_abc, ABS_MAX_AA_LEN);
	m_r2aa = SQCreateBigDigital(m_abc, ABS_MAX_AA_LEN);

	m_orfSearchAA = SQCreateBigDigital(m_abc, ABS_MAX_AA_LEN);//
	if(m_orfSearchAA->seq != NULL) free(m_orfSearchAA->seq);
	m_orfSearchAA->seq = (char *)malloc(ABS_MAX_AA_LEN * sizeof(char));
	
	m_R1prof.resize(m_hmm->M, '-');
	m_R2prof.resize(m_hmm->M, '-');
	m_R1profInd.resize(m_hmm->M,-1);
	m_R2profInd.resize(m_hmm->M,-1);

	m_mergedPairs = (char *)malloc(ABS_MAX_NT_LEN * sizeof(char));

	m_mergedAA = (char *)malloc(ABS_MAX_AA_LEN * sizeof(char));

	m_aabuf = (char *)malloc(ABS_MAX_AA_LEN * sizeof(char));
	m_dnabuf = (char *)malloc(ABS_MAX_NT_LEN * sizeof(char));

	return 0;
}
//////////////////////////////////////////////////////////////////////////////
//  reset everything after each sequence is processed
int IgHMM::ClearStrings(){
	
	m_R1prof = ""; m_R2prof = "";
	m_alignMethod = "";

	m_MSA = "";

    m_HFW1.Clear(); m_HCDR1.Clear(); m_HFW2.Clear(); m_HCDR2.Clear();
    m_HFW3.Clear(); m_HCDR3.Clear(); m_HFW4.Clear(); 
    m_LFW1.Clear(); m_LCDR1.Clear(); m_LFW2.Clear(); m_LCDR2.Clear();
    m_LFW3.Clear(); m_LCDR3.Clear(); m_LFW4.Clear();

    m_Hkabatstr = "";
    m_Lkabatstr = "";
    m_kabatstr = "";

	m_Hres2kabat.clear();
	m_Lres2kabat.clear();

	m_HAA = "";
	m_LAA = "";
	m_HDNA = "";
	m_LDNA = "";
	m_HDNAstart = -1;
	m_LDNAstart = -1;

	return 0;
}
//////////////////////////////////////////////////////////////////////////////
// align a single AA seq to the HMM
int IgHMM::AlignSingleAA(std::string & aaseq){
	m_mergedPairs[0] = '\0';
	m_mergedAA[0] = '\0';
	ClearStrings();
	m_score = -1;
	m_bestFrame1 = -1;
	m_HDNA="";
	m_LDNA="";

	m_alignMethod = "S";

	//Just copy it to mergedAA
	memcpy(m_mergedAA, aaseq.data(), aaseq.size()*sizeof(char));
	m_mergedAA[aaseq.size()] = '\0';

	int mergedAAlen = aaseq.size();
	
	//null out the DNA string so no attempts are made to process it
	m_mergedPairs[0] = '\0';

	TraceSeq(m_mergedAA, mergedAAlen, m_tr, &m_score); // no dependence on DNA
	TraceDomains(m_tr, m_sq);

	return 0;
}
//////////////////////////////////////////////////////////////////////////////
// Align a single DNA sequence to the HMM
int IgHMM::AlignSingle(char *r1, int r1len){

	m_mergedPairs[0] = '\0';
	m_mergedAA[0] = '\0';
	ClearStrings();
	m_score = -1;
	m_bestFrame1 = -1;

	m_alignMethod = "S";
	
	ESL_SQ msq;
	SetupSQ(r1,r1len,&msq);

	int bestFrame;
	int bestRC;
	float bestScore;

	FindBestORF(msq, 
		bestFrame, 
		bestRC,
		bestScore,
		m_r1Traces,
		m_r1aa);

	if(bestRC == 1){ gTranslator.Revcomp(msq.seq,msq.n); }

	//copy the final dna seq to its expected variable
	strncpy(m_mergedPairs, msq.seq, msq.n); m_mergedPairs[msq.n] = '\0';

	if(bestFrame != -1){
		int mergedAAlen;
		gTranslator.Translate(msq.seq, msq.n,
			bestFrame, 0, m_mergedAA, mergedAAlen);
		
		m_bestFrame1 = bestFrame;

		//Run the HMM
		TraceSeq(m_mergedAA, mergedAAlen, m_tr, &m_score);
		TraceDomains(m_tr, m_sq);

	}

	
	return 0;
}
//////////////////////////////////////////////////////////////////////////////
// align a paired end read to the HMM
int IgHMM::AlignPaired(char *r1, int r1len, char *r2, int r2len,
	int & mergedLen, int & mergedFrame){

	ClearStrings();
	m_bestFrame1 = -1;
	m_score = -1;

	if( ((r1len + r2len) >= ABS_MAX_NT_LEN) ||
	  ((r1len + r2len)/3 >= ABS_MAX_AA_LEN)){
		fprintf(stderr, "R1,R2 Sequence too large %d,%d\n",r1len,r2len);
		exit(-1);
	}

	m_alignMethod = "";
	
	mergedLen = 0;
	mergedFrame = -1;
	m_mergedPairs[0] = '\0';
	m_mergedAA[0] = '\0';

	int hit = 0;
	//First try by DNA...
	hit = AlignPairedByNT(r1, r1len, r2, r2len, m_mergedPairs, mergedLen);
	if(hit == 0){
		// if DNA failed, try HMM...
		AlignPairedByHMM(r1, r1len, r2, r2len,
			m_mergedPairs, mergedLen, mergedFrame);
	}else{
		//Find the best orf of the NT merged reads
		ESL_SQ msq;
		SetupSQ(m_mergedPairs,mergedLen,&msq);


		int bestRC;
		float bestScore;
		FindBestORF(msq, 
			mergedFrame, 
			bestRC,
			bestScore,
			m_r1Traces,
			m_r1aa);

		if(bestRC == 1){ gTranslator.Revcomp(m_mergedPairs,mergedLen); }
	}

	////////////
	//
	if(mergedFrame != -1){
		int mergedAAlen;
		gTranslator.Translate(m_mergedPairs, mergedLen,
			mergedFrame, 0, m_mergedAA, mergedAAlen);
		
		m_bestFrame1 = mergedFrame;
		//Run the HMM
		TraceSeq(m_mergedAA, mergedAAlen, m_tr, &m_score);
		TraceDomains(m_tr, m_sq);
	}

	return 0;
}


//////////////////////////////////////////////////////////////////////////////

void Nton(char *s, int n){ for(int i=0;i<n;i++){ if(s[i] == 'N') s[i] = 'n';} }
void ntoN(char *s, int n){ for(int i=0;i<n;i++){ if(s[i] == 'n') s[i] = 'N';} }

//////////////////////////////////////////////////////////////////////////////
// when aligned paired end reads, first look for a substantial overlap
// of DNA sequence
int IgHMM::AlignPairedByNT(char *r1, int r1len, char *r2, int r2len,
		char * mergedReads, int & mergedLen){

	mergedReads[0] = '\0';
	mergedLen = 0;

	//set all 'N' in s1 to 'n' so they don't match
	Nton(r1,r1len);
	
	int hit = 0;
	int bestS1 = -1;
	int bestS2 = -1;
	for(int rc2=1; rc2>=0 && hit<1 ; rc2--){
		gTranslator.Revcomp(r2, r2len);
		for(int s1=0; s1<=(r1len-m_NT_PAIR_MIN) && hit<1 ; s1++){
			for(int s2=0; s2<=(r2len-m_NT_PAIR_MIN) &&  hit<1  ; s2++){
				if(strncmp(r1+s1, r2+s2, m_NT_PAIR_MIN) == 0){
					hit = 1;
					bestS1 = s1;
					bestS2 = s2;
					break;
				}
			}
		}
	}
	ntoN(r1,r1len);
	if(hit != 0){
		int startRead;
		MergeByMatchPos(bestS1, bestS2,
			r1, r1len,
			r2, r2len,
			mergedReads,
			mergedLen,
			startRead);

		m_alignMethod = "NT_" + 
			boost::lexical_cast<std::string>( m_NT_PAIR_MIN );

	}
	return hit;
}

//////////////////////////////////////////////////////////////////////////////
void PrintDigitalSeq(ESL_SQ *sq){
	for(int i=1; i <= sq->n; i++){
		printf("%c", sq->abc->sym[sq->dsq[i]]);
	}
}

//////////////////////////////////////////////////////////////////////////////
int DumpTrace(ESL_DSQ *dsq, P7_TRACE *tr, P7_PROFILE * gm){
	int z;
	FILE *fp = stdout;
	float sc_ = 0.0;
	float accuracy_ = 0.0;
	
	fprintf(fp, "st   k     i      transit emission postprob - traceback len %d\n", tr->N);
	fprintf(fp, "--  ---- ------  -------- -------- --------\n");
	for (z = 0; z < tr->N; z++) {
		float tprob = -1.0;
		float eprob = -1.0;
		float pprob = -1.0;
		char emitted = ' ';

		char state = tr->st[z];
		char nextstate = tr->st[z+1];
		int nextkstate = tr->k[z+1];
		char prevstate = tr->st[z-1];
		char stateType = *p7_hmm_DecodeStatetype(state);
		int kstate = tr->k[z];
		int istate = tr->i[z];

		if(z < (tr->N - 1)){
			p7_profile_GetT(gm, state, kstate, nextstate, nextkstate, &tprob);
		}else tprob = 0.0f;
			
		fprintf(fp, "%c  %4d %6d  %8.4f", stateType,  kstate,istate, tprob );

		std::string dom = "";

		if (tr->pp != NULL) pprob = tr->pp[z];
		if (dsq != NULL) {
			int xi = 0;
			xi = dsq[istate];
			if (state == p7T_M) {
				eprob = p7P_MSC(gm, kstate, xi);
				emitted = gm->abc->sym[xi];
			} else if (state == p7T_I) {
				eprob = p7P_ISC(gm, kstate, xi);
				emitted = (char)tolower((int) gm->abc->sym[xi]);
			} else if ((state == p7T_N && prevstate == p7T_N) ||
				(state == p7T_C && prevstate == p7T_C) ||
				(state == p7T_J && prevstate == p7T_J)) {

				eprob = 0;
				emitted = (char)tolower((int) gm->abc->sym[xi]);
			}
		}
		fprintf(fp, " %8.4f %8.4f %c %s\n", eprob,pprob,emitted, dom.c_str());
	}
	fprintf(fp, "                -------- -------- --------\n");
	fprintf(fp, "                  total: %8.4f %8.4f\n\n", sc_, accuracy_);

	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// generate a vector of states for each AA for a given alignment
void FillProf(P7_TRACE *tr, ESL_SQ *sq, 
	std::string & prof, std::vector<int> & profInd, 
	int & leftMatch,
	int & rightMatch ){

	int size = prof.size();
	char *profsq = (char *) prof.data();
	memset(profsq, '-', (size)*sizeof(char));
	
	leftMatch = -1;
	rightMatch = -1;
	for(int z=0; z < tr->N ; z++){
		char state = tr->st[z];
		int istate = tr->i[z];
		int kstate = tr->k[z];
		int xi = sq->dsq[istate];
		if (state == p7T_M) {
			if(leftMatch == -1){ leftMatch = kstate; }
			rightMatch = kstate;
			char emitted = sq->abc->sym[xi];
			profsq[kstate] = emitted;
			profInd[kstate] = istate - 1; //istates are 1-based, conv to 0-base
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// print an MSA from an alignment
int IgHMM::GetMSA2(ESL_SQ *s1, P7_TRACE *tr1, 
	ESL_SQ *s2, P7_TRACE *tr2,
	std::string & msa1, std::string & msa2
	){

	ESL_SQ *bothseqs[2];
	bothseqs[0] = s1;
	bothseqs[1] = s2;
	P7_TRACE *bothtraces[2];
	bothtraces[0] = tr1;
	bothtraces[1] = tr2;
	int msaopts = 0; msaopts |= p7_ALL_CONSENSUS_COLS;
	ESL_MSA *debugMSA = NULL;
	p7_tracealign_Seqs(bothseqs, bothtraces, 2, m_hmm->M, 
		msaopts, m_hmm, &debugMSA);

	msa1 = debugMSA->aseq[0];
	msa2 = debugMSA->aseq[1];

	esl_msa_Destroy(debugMSA);
	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// given two MSA, find the first overlapping match states shared 
// between them
int	FindFirstOverlappingMatchStates( 
		std::string & msa1,
		std::string & msa2,
		int & r1AApos, int & r2AApos, int & score){

	score = 0;

	r1AApos = 0;
	r2AApos = 0;
	for(unsigned int i=0; i<msa1.size();i++){
		int isaa1 = isalpha(msa1[i]);//a-zA-Z
		int isaa2 = isalpha(msa2[i]);
		int match1 = isupper(msa1[i]); //A-Z
		int match2 = isupper(msa2[i]); //A-Z
		if(match1 && match2){
			if(msa1[i] == msa2[i]){
				score = 1;
				if(i<(msa1.size()-1)){
					if(msa1[i+1] == msa2[i+1]){
						score = 2;
					}
				}
				return 0;
			}
		}
		if(isaa1) r1AApos++;
		if(isaa2) r2AApos++;
	}

	return 0;
}	

////////////////////////////////////////////////////////////////////////////
// merge paired end reads by a shared match position
int IgHMM::MergeByMatchPos(int r1NTPos, int r2NTPos,
	char * r1Seq, int r1N,
	char * r2Seq, int r2N,
	char *mergedReads,
	int & mergedLen,
	int & startRead){

	//Who has more sequence to the left?
	if(r1NTPos > r2NTPos){ //R1 does, start with R1
		//Who has more sequence to right?
		if((r2N - r2NTPos) > (r1N - r1NTPos) ){ //R2 does, end with R2
			mergedLen = r1NTPos + r2N - r2NTPos;
			memcpy(mergedReads, r1Seq, r1NTPos + 1);
			memcpy(mergedReads+r1NTPos, r2Seq + r2NTPos,
				r2N - r2NTPos);
			
		}else{ //juse use all of r1
			mergedLen = r1N;
			memcpy(mergedReads, r1Seq, r1N);
		}
		startRead = 1;
	}else{ // R2 does, start with R2
		if((r1N - r1NTPos) > (r2N - r2NTPos) ){ //R1 does,end with R1
			mergedLen = r2NTPos + r1N - r1NTPos;
			memcpy(mergedReads, r2Seq, r2NTPos + 1);
			memcpy(mergedReads+r2NTPos, r1Seq + r1NTPos,
				r1N - r1NTPos);
		}else{ //juse use all of r2
			mergedLen = r2N;
			memcpy(mergedReads, r2Seq, r2N);
		}
		startRead = 2;
	}
	mergedReads[mergedLen] = '\0';

	return 0;
}

////////////////////////////////////////////////////////////////////////////
// if not enough DNA to mate-pairs, try to rescue with the HMM
int IgHMM::PairTraces(
	ESL_SQ * r1sq, P7_TRACE *tr1, ESL_SQ *s1, int frame1,
	ESL_SQ * r2sq, P7_TRACE *tr2, ESL_SQ *s2, int frame2,
	char * mergedReads, int & mergedLen, int & mergedFrame){

	memset(mergedReads, 'N', ABS_MAX_NT_LEN * sizeof(char));

	std::string msa1,msa2;
	GetMSA2(m_r1aa, tr1, m_r2aa, tr2, msa1, msa2);

	//By this point the sequences are revcomped already
	// so no need to consider rc
	//Find first overlapping match states.
	//
	int r1AApos, r2AApos;
	int score;
	FindFirstOverlappingMatchStates( 
		msa1,msa2,
		r1AApos, r2AApos, score);


	if(score>0){  //matching M states
		m_alignMethod = "HMM_MATCH" + 
			boost::lexical_cast<std::string>( score  );

		int r1NTPos = MapAAToNT(r1sq->n, s1->n,  r1AApos, frame1, 0);
		int r2NTPos = MapAAToNT(r2sq->n, s2->n,  r2AApos, frame2, 0);

		int startRead;
		MergeByMatchPos(r1NTPos, r2NTPos,
			r1sq->seq, r1sq->n,
			r2sq->seq, r2sq->n,
			mergedReads,
			mergedLen,
			startRead);

		if(startRead == 1) mergedFrame = frame1;
		if(startRead == 2) mergedFrame = frame2;
	}else{
		// Otherwise just stitch them, in frame, a run of X amino acids
		// these are basically incomplete sequences and 
		// will likely be discarded anyway
		// so treatment is somewhat arbitrary
		int f1adjust = (3 - (r1sq->n - frame1) % 3) % 3;
		f1adjust = f1adjust + (3-frame2) % 3 ;
		f1adjust = f1adjust % 3;
		int nToAdd = 21 + f1adjust;  // add 7 to 8 N codons depending on frame
		mergedLen = r1sq->n + r2sq->n + 21;
		memcpy(mergedReads, r1sq->seq, r1sq->n*sizeof(char));
		memset(mergedReads + r1sq->n, 'N', nToAdd * sizeof(char)); 
		memcpy(mergedReads + r1sq->n + nToAdd, r2sq->seq, r2sq->n);
		mergedReads[mergedLen] = '\0';
		m_alignMethod = "HMM_FRAME";
		mergedFrame = frame1;
	}


	return 0;
}
//////////////////////////////////////////////////////////////////////////////
void PrintMSA(FILE *f, ESL_MSA *msa){
	printf("pp_cons:%s\n", msa->pp_cons);
	printf("rf     :%s\n", msa->rf);
	printf("aseq   :%s\n", msa->aseq[0]);
	printf("pp     :%s\n", msa->pp[0]);
}
int IgHMM::GetMSA(ESL_SQ *seq, P7_TRACE *tr, std::string & msastr){
	int msaopts = 0; msaopts |= p7_ALL_CONSENSUS_COLS;
	ESL_MSA *msa = NULL;
	p7_tracealign_Seqs(&seq, &tr, 1, m_hmm->M, msaopts, m_hmm, &msa);

	msastr = msa->aseq[0];

	esl_msa_Destroy(msa);
	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// pair end read mating using the HMM
int IgHMM::AlignPairedByHMM(char *r1, int r1len, char *r2, int r2len,
	char * mergedReads, int & mergedLen, int & mergedFrame){

	ESL_SQ r1sq;
	SetupSQ(r1,r1len,&r1sq);
	
	int bestFrame1;
	int bestRC1;
	float bestScore1;

	ESL_SQ r2sq; 
	SetupSQ(r2,r2len,&r2sq);

	int bestFrame2;
	int bestRC2;
	float bestScore2;

	FindBestORF(r1sq, bestFrame1, bestRC1,  bestScore1, m_r1Traces, m_r1aa);
	FindBestORF(r2sq, bestFrame2, bestRC2,  bestScore2, m_r2Traces, m_r2aa);

	P7_TRACE *r1r2BestTraces[2];
	r1r2BestTraces[0] =  m_r1Traces[frameMap[bestRC1][bestFrame1]];
	r1r2BestTraces[1] =  m_r2Traces[frameMap[bestRC2][bestFrame2]];

	//Fix the sequences here so rc doesn't matter anymore
	if(bestRC1 == 1){
		gTranslator.Revcomp(r1sq.seq, r1sq.n);
		bestRC1 = 0;
	}
	if(bestRC2 == 1){
		gTranslator.Revcomp(r2sq.seq, r2sq.n);
		bestRC2 = 0;
	}

	PairTraces(
		&r1sq, r1r2BestTraces[0], m_r1aa, bestFrame1,
		&r2sq, r1r2BestTraces[1], m_r2aa, bestFrame2,
		mergedReads, mergedLen, mergedFrame);

	return 0;
}
//////////////////////////////////////////////////////////////////////////////
// align every ORF and find the ORF with the highest score
int IgHMM::FindBestORF(ESL_SQ & rsq, 
	int & bestFrame, int & bestRc,  float & bestScore,
	P7_TRACE *trace[6],
	ESL_SQ *raa
	){

	if(rsq.seq == NULL || rsq.n < 1 ){ return -1;}

	float curScore;
	int firstLoop = 1;
	for(int rcFlag=0; rcFlag<2; rcFlag++){
		for(int frameshift=0; frameshift < 3; frameshift++){

			int aalen;

			gTranslator.Translate(rsq.seq,rsq.n, frameshift, rcFlag,
				m_orfSearchAA->seq, aalen);

			m_orfSearchAA->n = aalen;
			m_orfSearchAA->end = aalen;
			m_orfSearchAA->W = aalen;
			m_orfSearchAA->L = aalen;

			TraceSeq(m_orfSearchAA->seq,(int) m_orfSearchAA->n, 
				trace[frameMap[rcFlag][frameshift]], &curScore);

			if(firstLoop){
				firstLoop = 0;
				bestScore = curScore;
				bestRc = rcFlag;
				bestFrame = frameshift;
				esl_sq_Copy(m_sq, raa);//Use the builtin digital sq
			}else{
				if(curScore > bestScore){
					bestScore = curScore;
					bestRc = rcFlag;
					bestFrame = frameshift;
					esl_sq_Copy(m_sq, raa);//Use the builtin digital sq
				}
			}

		}
	}


	return 0;
}


//////////////////////////////////////////////////////////////////////////////
void CountMatches(P7_TRACE *tr, int & numMatches ){

	numMatches = 0;
	for(int z=0; z < tr->N ; z++){
		char state = tr->st[z];
		if (state == p7T_M) {
			numMatches++;
		}
	}
}

//////////////////////////////////////////////////////////////////////////////
// align an AA sequence to the HMM
int IgHMM::TraceSeq(char *seq, int seqlen,
	P7_TRACE *trace, float *tscore){

	int retval;

	m_score = 0;

	m_sq->dsq[0] = eslDSQ_SENTINEL;
	int seqn = 1;
	for(int i=0;i<seqlen; i++){
		ESL_DSQ x;
		x = m_sq->abc->inmap[(int) seq[i]];
		if(x<=127){
			m_sq->dsq[seqn] = x;
			seqn++;
		}
	}
	m_sq->dsq[seqn] = eslDSQ_SENTINEL;

	m_sq->n = seqn-1;
	m_sq->start = 1;
	m_sq->end = m_sq->n;
	m_sq->C = 0;
	m_sq->W = m_sq->n;
	m_sq->L = m_sq->n;

	if(trace == NULL) {
		fprintf(stderr, "NULL trace object in IgHMM::TraceSeq\n");
		exit(-1);
	}
	p7_trace_Reuse(trace);
	retval = p7_tracealign_computeTraces(m_hmm, &m_sq, 0, 1, &trace);
	if(retval != eslOK){
		fprintf(stderr,"floating point overflow\n");
	}

	retval = p7_trace_Score(trace, m_sq->dsq, m_gm, tscore);
	if(retval != eslOK){
		fprintf(stderr,"floating point overflow\n");
	}

	return 0;
}
//////////////////////////////////////////////////////////////////////////////
// pretty print functions to examine HMM alignments,
// used for debugging
int IgHMM::STDoutMSA(){
	return STDoutMSA(m_msa,m_score);
}
int IgHMM::STDoutMSA(ESL_MSA *somemsa, float somescore){
	FILE *ofp = stdout;
	int outfmt  = eslMSAFILE_STOCKHOLM;
	eslx_msafile_Write(ofp, somemsa, outfmt);
	fprintf(ofp, "Score:%f\n",somescore);

	return 0;
}

int IgHMM::STDoutMSA2(ESL_MSA *somemsa, float score1, float score2){
	FILE *ofp = stdout;
	int outfmt  = eslMSAFILE_STOCKHOLM;
	eslx_msafile_Write(ofp, somemsa, outfmt);
	fprintf(ofp, "Score1:%f\tScore2:%f\n",score1, score2);

	return 0;
}

int IgHMM::OneLineMSA2(ESL_MSA *r1r2msa, float score1, float score2){
	printf("R1:  %s\nR2:  %s\n",r1r2msa->aseq[0],r1r2msa->aseq[1]);
	printf("Score1:%f\nScore2:%f\n",score1, score2);
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Sanity check a domain start and stop make sure it matches the 
// query sequence,
int SanCheckDomInterval(Domain & dom, char *dna){
	char buf[4096];
	buf[0] = '\0';
	int aalen;
	if(dom.DNAstart != -1 || dom.DNAstop != -1){
		gTranslator.Translate(dna + dom.DNAstart, 
			dom.DNAstop - dom.DNAstart + 1,
			0, 0, buf, aalen);
	}
	std::string tmp = dom.str;
	boost::replace_all(tmp,"-",""); boost::to_upper(tmp);
	std::string tr = buf;
	if(tmp.compare(tr)!=0){
		printf("dna bounds not correct!\n");
		exit(-1);
	}
	return 0;
}
//////////////////////////////////////////////////////////////////////////////
// deteremine DNA start and stop from dom size and location
int GetDomainInterval(Domain & dom, std::string::size_type & loc, int frame){
	int size = 0;
	// only non gap chars count towards size
	for(unsigned int i=0;i<dom.str.size();i++){
		if(dom.str[i] != '-'){ size++;} 
	}
	if(size>0){
		dom.AAstart = loc;
		dom.AAstop = loc + size - 1;
		//convert it to DNA, its always rc=0 so no need for anything
		// but frame.
		dom.DNAstart = MapAAToNT(0, 0, dom.AAstart, frame, 0);
		dom.DNAstop = MapAAToNT(0, 0, dom.AAstop, frame, 0) + 2; //include the whole codon
		// update loc so next domain calc starts from correct position
		loc = loc + size;
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// Get the DNA corresponding to a subset of the AA seq
// e.g. the HFW1 
int GetDomainDNA(std::string dom, int & domOffset, std::string & domNT) {
	return 0;

	domNT = "";
	if(dom.size()>0){
		int loc = 3*domOffset;
		domNT = 
			boost::lexical_cast<std::string>( loc ) + "-" +
			boost::lexical_cast<std::string>( loc+3*dom.size() - 1 ) ;
			
		domOffset += dom.size();	
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Get a domain sequence using the MSA and start and stop match state
//	
//     ||||||||
//   aaF-FFFFFFbbb
//
//   aa-FFFFFFFbbb
//    aAFFFFFFFbbb
//
//   aa--FFFFFFbbb
//     AAFFFFFFbbb
int GetMatchStates(std::string &msa, int start, int stop,
	int ins5p, //include 5p insertions
	int ins3p, //include 3p insertions
	std::string & dom,
	int toUpper){ 

	int firstMatchState = -1;
	
	int curMatchState = 1;
	for(unsigned int i=0;i<msa.size(); i++){
		if( (curMatchState >= start) && (curMatchState <= stop)){
			if(curMatchState == start && !ins5p && islower(msa[i])){
			}else{
				if(firstMatchState == -1) firstMatchState = i;
				if(toUpper){
					dom.append(1,toupper(msa[i]));
				}else{
					dom.append(1,msa[i]);
				}
			}
		}
		if(ins3p && curMatchState == (stop+1) && islower(msa[i]) ) {
			if(toUpper){
				dom.append(1,toupper(msa[i]));
			}else{
				dom.append(1,msa[i]);
			}
		}
		if(isupper(msa[i]) || msa[i] == '-'){
			curMatchState++;
		}
	}

	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// generatlly Kabat insertion go from A-Z but if by chance one is over 26
// use AA, AB , etc.
std::string KabatInsert(int insertpos){
	if(insertpos <= 0){ return "";};
	if(insertpos <= 26){
		return boost::lexical_cast<std::string>(LETTERS[insertpos]);
	}
	return boost::lexical_cast<std::string>(LETTERS[insertpos/26]) + 
		boost::lexical_cast<std::string>(LETTERS[insertpos%26]) ;
}



//////////////////////////////////////////////////////////////////////////////
//  generate the kabat string for a framework region
void KabatFW(std::string & chain, std::string & kabat, 
		int kpos, std::string & fw){
	const char *delim = " ";
	std::string kloc;
	int insertsize = 0;;
	for(unsigned int i=0; i<fw.size(); i++){
		if(islower(fw[i])){
			kloc = chain + boost::lexical_cast<std::string>(kpos-1);
			insertsize++;
			kloc = kloc + KabatInsert(insertsize); 
		}else{
			kloc = chain + boost::lexical_cast<std::string>(kpos);
			kpos++;
			insertsize = 0;
		}

		kabat = kabat + kloc + "." + 
			(char)(toupper(fw[i])) + delim;
	}
}


//////////////////////////////////////////////////////////////////////////////
// generate a kabat string for a CDR region
void KabatCDRInfill(std::string & chain, int startKpos, int endKpos, int wigglepos, int delpos, std::string & cdr, std::string & kabat){

	const char *delim = " ";
	//nominal kabat size, e.g CDRL1 goes from 24 to 34 so it is 11
	int cdrsize = endKpos - startKpos + 1;

	if((int)cdr.size() > cdrsize){ //add wiggles
		int wigglesToAdd = cdr.size() - cdrsize;

		int cdrpos = 0;
		for(int pos=startKpos;pos<=wigglepos;pos++){
			kabat = kabat + chain + boost::lexical_cast<std::string>(pos) +
				"." + cdr[cdrpos] + delim;
			cdrpos++;
		}
		//Add wiggles
		for(int pos=1;pos<=wigglesToAdd;pos++){
			kabat = kabat + 
				chain + boost::lexical_cast<std::string>(wigglepos) +
				KabatInsert(pos) + "." + cdr[cdrpos] + delim;
			cdrpos++;
		}
		for(int pos=wigglepos+1;pos<=endKpos; pos++){
			kabat = kabat + chain + boost::lexical_cast<std::string>(pos) +
				"." + cdr[cdrpos] + delim;
			cdrpos++;
		}
		
	}else{
		//How many positions do we need to delete?
		int npos2del = cdrsize - (int)cdr.size();

		int skippos;
		skippos = std::max(startKpos-1, delpos - npos2del);

		std::string kloc;
		int tpos = startKpos - 1;
		for(unsigned int i=0;i<cdr.size();i++){
			tpos = startKpos+i;
			if(tpos > skippos){
				if(tpos == skippos+1){//fill in deletions
					for(int fill=tpos;fill<(tpos+npos2del);fill++){
						kabat = kabat + chain + 
							boost::lexical_cast<std::string>(fill) +
							".-" + delim;
					}
					
				}
				tpos = tpos + npos2del;
			}
			kloc = chain + boost::lexical_cast<std::string>(tpos);
			kabat = kabat + kloc + "." + (char)(toupper(cdr[i])) + delim;
		}
		if(tpos<endKpos){
			for(int fill=tpos+1;fill<=endKpos;fill++){
				kabat = kabat + chain + 
					boost::lexical_cast<std::string>(fill) + ".-" + delim;
			}
		}
	}

}

//////////////////////////////////////////////////////////////////////////////
void KabatCDRInfill_Test(){
	std::string chain = "L";
	int startKpos = 93;
	int endKpos = 102;
	int wigglepos = 100;
	int delpos = 100;

	std::string cdrlong = "123456789012345678901234567890123456789012345678901234567890";

	for(int i=(endKpos - startKpos + 3);i>=0; i--){
		std::string cdr = cdrlong.substr(0,i);
		std::string kabat = "";
		KabatCDRInfill(chain, startKpos,endKpos,wigglepos,delpos, cdr,kabat);
		printf("%s\n%s\n\n\n",cdr.c_str(),kabat.c_str());
	}

	exit(0);
}

//////////////////////////////////////////////////////////////////////////////
// sanity check the kabat seq against the whole variable domain
void SanCheckKabat(std::string & kstr, std::string vardom){
	int ok = 1;
	boost::replace_all(vardom,"-","");
	for(unsigned int i=0;i<vardom.size();i++){ vardom[i] = toupper(vardom[i]); }

	std::vector<std::string> kstrs;

	boost::split(kstrs, kstr, boost::is_any_of(" "), boost::token_compress_on);

	std::string fk = "";
	boost::regex b4(".*\\.");
	boost::regex after("\\..*");
	boost::regex letters("[A-Z]+");
	boost::regex numbers("[0-9]+");
	int lpos = 0;
	std::string lins="";
	for(unsigned int i=0;i<kstrs.size(); i++){
		if(kstrs[i].size() == 0){continue;}
		std::string prefix = boost::regex_replace(kstrs[i], after, "");
		prefix = prefix.substr(1); //drop the L/H
		std::string num = boost::regex_replace(prefix, letters, "");
		std::string ins = boost::regex_replace(prefix, numbers, "");

		int tpos = boost::lexical_cast<int>(num);
		ok = ok && (
			( tpos == (lpos+1)  && ins.size()==0 ) || //52 to 53
			( tpos == lpos && ins.compare(lins)) //52X to 52Y
			);
		if(!ok){
			fprintf(stderr,"kabatnum problem\n");
			exit(-1);
		}


		lpos = tpos;
		lins = ins;

		std::string aa = boost::regex_replace(kstrs[i], b4, "");
		fk = fk + aa;
	}

	boost::replace_all(fk,"-","");
	ok = ok && (fk.compare(vardom)==0);
	if(!ok){
		fprintf(stderr,"kabat != vardom\n");
		exit(-1);
	}

}

//////////////////////////////////////////////////////////////////////////////
// generate heavy chain kabat numbering
// 
//  The numbers used in the FW and CDR calls here are those
// that happen to be the match states in the HMMs that correspond
// to the domain boundaries.   We rely on these conserved states to 
// extract the domains
int IgHMM::GenHeavyKabat(){
	std::string tmp;
	std::string chain = "H";

	/////////////////////////////////
	//FW1
	KabatFW(chain, m_Hkabatstr, 1, m_HFW1.str);
	
	/////////////////////////////////
	//CDR1 
	KabatCDRInfill(chain, 31, 35, 35, 35, m_HCDR1.str, m_Hkabatstr);
	

	/////////////////////////////////
	//FW2
	KabatFW(chain, m_Hkabatstr, 36, m_HFW2.str);

	/////////////////////////////////
	//CDR2
	KabatCDRInfill(chain, 50, 65, 52, 52, m_HCDR2.str, m_Hkabatstr);

	/////////////////////////////////
	//
	//find the first 17 match pos which = H82 after which we force 82ABC
	int match17pos = 0;
	int matchpos = 0;
	for(unsigned int i=0;i<m_HFW3.str.size(); i++){
		if(isupper(m_HFW3.str[i]) || m_HFW3.str[i] == '-'){ matchpos++; }
		if(matchpos == 17){ match17pos = i; };
	}
	tmp = m_HFW3.str.substr(0,match17pos+1);
	KabatFW(chain, m_Hkabatstr, 66, tmp);
	m_Hkabatstr = m_Hkabatstr + chain + "82A." + m_HFW3.str[match17pos+1] + " ";
	m_Hkabatstr = m_Hkabatstr + chain + "82B." + m_HFW3.str[match17pos+2] + " ";
	m_Hkabatstr = m_Hkabatstr + chain + "82C." + m_HFW3.str[match17pos+3] + " ";
	tmp = m_HFW3.str.substr(match17pos+4);
	KabatFW(chain, m_Hkabatstr, 83, tmp);

	/////////////////////////////////
	//CDR3
	KabatCDRInfill(chain, 93, 102, 100, 100, m_HCDR3.str, m_Hkabatstr);

	/////////////////////////////////
	//FW4
	KabatFW(chain, m_Hkabatstr, 103, m_HFW4.str);

	// for testing
	if(0){
		SanCheckKabat(m_Hkabatstr, 
			m_HFW1.str + m_HCDR1.str +
			m_HFW2.str + m_HCDR2.str +
			m_HFW3.str + m_HCDR3.str +
			m_HFW4.str);
	}
	return 0;
}


//////////////////////////////////////////////////////////////////////////////
// generate heavy chain kabat numbering
//
//  The numbers used in the FW and CDR calls here are those
// that happen to be the match states in the HMMs that correspond
// to the domain boundaries.   We rely on these conserved states to 
// extract the domains

int IgHMM::GenLightKabat(){

	std::string chain = "L";
	/////////////////////////////////
	//FW1
	KabatFW(chain, m_Lkabatstr, 1, m_LFW1.str);
	
	/////////////////////////////////
	//CDR1
	KabatCDRInfill(chain, 24, 34, 27, 27, m_LCDR1.str, m_Lkabatstr);

	/////////////////////////////////
	//FW2
	KabatFW(chain, m_Lkabatstr, 35, m_LFW2.str);

	/////////////////////////////////
	//CDR2
	KabatCDRInfill(chain, 50, 56, 54, 54, m_LCDR2.str, m_Lkabatstr);

	/////////////////////////////////
	//FW3 
	KabatFW(chain, m_Lkabatstr, 57, m_LFW3.str);

	/////////////////////////////////
	//CDR3
	KabatCDRInfill(chain, 89, 97, 95, 95, m_LCDR3.str, m_Lkabatstr);

	/////////////////////////////////
	//FW4
	KabatFW(chain, m_Lkabatstr, 98, m_LFW4.str);

	
	// for testing
	if(0){
		SanCheckKabat(m_Lkabatstr, 
			m_LFW1.str + m_LCDR1.str +
			m_LFW2.str + m_LCDR2.str +
			m_LFW3.str + m_LCDR3.str +
			m_LFW4.str);
	}
	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Keep track of what each residue in HAA and LAA map to in kabat positions
// for SHM tallies.
int IgHMM::MapResiduePositions2Kabat(){
	m_Hres2kabat.clear();
	m_Lres2kabat.clear();

	unsigned int i = 0;
	unsigned int j;
	while(i<m_kabatstr.size()){
		int docheck = 0;

		docheck = (i==0 && (m_kabatstr[i] == 'H' || m_kabatstr[i] == 'L')) ||
			(i>0 && (m_kabatstr[i] == 'H' || m_kabatstr[i] == 'L') &&
				(m_kabatstr[i-1] == ' ' || m_kabatstr[i-1] == ' ') );
		if(docheck){
			// seek forward to .
			for(j=i;j<m_kabatstr.size();j++){ if(m_kabatstr[j] =='.') break; }
			if(j<(m_kabatstr.size()-1)){
				if(m_kabatstr[j+1] != '-'){
					if(m_kabatstr[i] == 'H'){
						m_Hres2kabat.push_back( m_kabatstr.substr(i,j-i) );
					}else{
						m_Lres2kabat.push_back( m_kabatstr.substr(i,j-i) );
					}
				}
			}
			i = j;
		}
		i++;
	}

	return 0;
}

//////////////////////////////////////////////////////////////////////////////
//  Use the HMM to generate the domain (FW and CDR for heavy and light)
//  and generate the Kabat strings
int IgHMM::TraceHMM(int hoffset, int loffset, int doH, int doL){
	

	if(doH){
		GetMatchStates(m_MSA, hoffset +1,hoffset +26, 0,0, m_HFW1.str, 0);
		GetMatchStates(m_MSA, hoffset +27,hoffset +35, 1,1, m_HCDR1.str, 1);
		//transfer 4 aa 
		m_HFW1.str = m_HFW1.str + m_HCDR1.str.substr(0,4);
		m_HCDR1.str = m_HCDR1.str.substr(4);
		boost::replace_all(m_HCDR1.str,"-","");

		GetMatchStates(m_MSA, hoffset +36,hoffset +49, 0,0, m_HFW2.str,0);
		GetMatchStates(m_MSA, hoffset +50,hoffset +66, 1,1, m_HCDR2.str,1);
		boost::replace_all(m_HCDR2.str,"-","");

		GetMatchStates(m_MSA, hoffset +67,hoffset +94, 0,0, m_HFW3.str,0);
		GetMatchStates(m_MSA, hoffset +95,hoffset +111, 1,1, m_HCDR3.str,1);
		//transfer 2 aa
		m_HFW3.str = m_HFW3.str + m_HCDR3.str.substr(0,2);
		// first 5 and last 8 residues of FW3 are pretty conserved

		m_HCDR3.str = m_HCDR3.str.substr(2);
		boost::replace_all(m_HCDR3.str,"-","");

		GetMatchStates(m_MSA, hoffset +112,hoffset +122, 0,0, m_HFW4.str,0);
		GenHeavyKabat();
		m_kabatstr = m_kabatstr + m_Hkabatstr;
	}

	if(doL){
		GetMatchStates(m_MSA, loffset +1,loffset +23, 0,0, m_LFW1.str,0);
		GetMatchStates(m_MSA, loffset +24,loffset +35, 1,1, m_LCDR1.str,1);
		boost::replace_all(m_LCDR1.str,"-","");
		GetMatchStates(m_MSA, loffset +36,loffset +50, 0,0, m_LFW2.str,0);
		GetMatchStates(m_MSA, loffset +51,loffset +57, 1,1, m_LCDR2.str,1);
		boost::replace_all(m_LCDR2.str,"-","");
		GetMatchStates(m_MSA, loffset +58,loffset +89, 0,0, m_LFW3.str,0);
		GetMatchStates(m_MSA, loffset +90,loffset +98, 1,1, m_LCDR3.str,1);
		boost::replace_all(m_LCDR3.str,"-","");
		GetMatchStates(m_MSA, loffset +99,loffset +108, 0,0, m_LFW4.str,0);

		GenLightKabat();
		m_kabatstr = m_kabatstr + m_Lkabatstr;
	}

	// 1,000,000  in 19 seconds for full H and L
	MapResiduePositions2Kabat();

	return 0;
}

//////////////////////////////////////////////////////////////////////////////
// Get match state for each position in an MSA string.
// matchpos = if uppercase or '-' ++
// -matchpos if lowercase (insertion state);
void GetMatchPosVector(std::string & msa, std::vector<int> & matchPos){
	matchPos.resize(msa.size());
	int pos = 1;
	for(unsigned int i=0;i<msa.size(); i++){
		if(!islower(msa[i])){
			matchPos[i] = pos;
			pos++;
		}else{
			matchPos[i] = -pos;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// the hmm will often put initial AA as insert states 
// instead of populating the first match states. we will force it
// to use initial match states
void P3Pad(std::string & msa, int startMatch, int endMatch){

	std::vector<int> matchpos;
	GetMatchPosVector(msa, matchpos);

	int p3ins = 0;
	int p3roiDelstates = 0;
	int p3roiWithinDel = 1;
	int lastp3insPos = -1;
	int in3pinsState = 0;
	int lastp3roiDelPos = -1;
	for(int i=(matchpos.size()-1);i>=0;i--){
		// insertions at 3prime end

		if( i == ((int)matchpos.size()-1)){
			if(matchpos[i] < 0) in3pinsState = matchpos[i];
		}
		if(matchpos[i]  == in3pinsState){
			lastp3insPos = i;
			p3ins++;
		//A match state before our region of interest,
		// they should all be -'s, if not, no padding is done
		}else if(abs(matchpos[i])>endMatch){ // && matchpos[i] > 0){
			if(msa[i] != '-'){ 
				// 
				return;
			}
		// within region of interest
		}else if( (abs(matchpos[i]) >= startMatch) && 
			(abs(matchpos[i]) <= endMatch) )  {
			
			//are there 3p del (-) states? count them
			if(p3roiWithinDel){
				if(msa[i] != '-'){
					p3roiWithinDel = 0;
				}else{
					lastp3roiDelPos = i;
					p3roiDelstates++;
				}
			}
		}else if(matchpos[i] == (startMatch-1)){
			//If the roi is all dels "------"
			// and the very next match state is also a "-"
			// this sequence is likely hosed so don't try to pad/fix it
			if(p3roiDelstates == (endMatch-startMatch+1) && msa[i] == '-'){
				return;
			}
		}

	}

	if(p3roiDelstates == 0) {return;}

	int nchars;
	if(p3roiDelstates >= p3ins){
		nchars = p3ins;
	}else{
		nchars = p3roiDelstates;
	}

	for(int d=0;d<nchars;d++){
		msa[lastp3roiDelPos + d] = toupper( msa[lastp3insPos + d]);
		msa[lastp3insPos + d] = '.';
	}

	boost::erase_all(msa,".");

}
//////////////////////////////////////////////////////////////////////////////
// see P3Pad, same for 5 prime
void P5Pad(std::string & msa, int startMatch, int endMatch){

	std::vector<int> matchpos;
	GetMatchPosVector(msa, matchpos);

	int p5ins = 0;
	int p5roiDelstates = 0;
	int p5roiWithinDel = 1;
	int lastp5insPos = -1;
	int lastp5roiDelPos = -1;
	for(int i=0;i<(int)matchpos.size();i++){
		// insertions at 5prime end
		if(matchpos[i] == -1){
			lastp5insPos = i;
			p5ins++;
		//A match state before our region of interest,
		// they should all be -'s, if not, no padding is done
		}else if(abs(matchpos[i])<startMatch){ 
			if(msa[i] != '-'){ 
				return;
			}
		// within region of interest
		}else if( (abs(matchpos[i]) >= startMatch) && 
			(abs(matchpos[i]) <= endMatch) )  {
			
			//are there 5p del (-) states? count them
			if(p5roiWithinDel){
				if(msa[i] != '-'){
					p5roiWithinDel = 0;
				}else{
					lastp5roiDelPos = i;
					p5roiDelstates++;
				}
			}
		}else if(matchpos[i] == (endMatch+1)){
			//If the roi is all dels "------"
			// and the very next match state is also a "-"
			// this sequence is likely hosed so don't try to pad/fix it
			if(p5roiDelstates == (endMatch-startMatch+1) && msa[i] == '-'){
				return;
			}
		}
	}

	if(p5roiDelstates == 0) {return;}

	int nchars;
	if(p5roiDelstates >= p5ins){
		nchars = p5ins;
	}else{
		nchars = p5roiDelstates;
	}

	for(int d=0;d<nchars;d++){
		msa[lastp5roiDelPos - d] = toupper( msa[lastp5insPos - d]);
		msa[lastp5insPos - d] = '.';
	}

	boost::erase_all(msa,".");
}

//////////////////////////////////////////////////////////////////////////////
void TestP3Pad(){
	const char* msas[] = {
		"F-----cba",
		"F-----cb",
		"F-----a",
		"F-----",
		"F---P-cba",
		"F---p-cba",
		"F----cba",
		"F--cba",
		"F--cb",
		"F--c",
		"F--",
		""
	};

	for(int i=0;i<12;i++){
		std::string msa = msas[i];	
		printf("################\n");
		printf("%s\n",msa.c_str());
		std::string tmp = msa;

		msa = tmp;
		P3Pad(msa, 4,6);
		printf("%s\n",msa.c_str());

		msa = tmp;
		P3Pad(msa, 1,3);
		printf("%s\n",msa.c_str());
	};

	exit(0);

}
void TestP5Pad(){
	const char* msas[] = {
	"abc-----F",
	"bc-----F",
	"a-----F",
	"-----F",
	"abc-P---F",
	"abc-p---F",
	"abc----F",
	"abc--F",
	"bc--F",
	"c--F",
	"--F",
	""
	};

	for(int i=0;i<12;i++){
		std::string msa = msas[i];	
		printf("################\n");
		printf("%s\n",msa.c_str());
		std::string tmp = msa;

		msa = tmp;
		P5Pad(msa, 4,6);
		printf("%s\n",msa.c_str());

		if(i==5){
			printf("break\n");
		}
		msa = tmp;
		P5Pad(msa, 1,3);
		printf("%s\n",msa.c_str());
	};


	exit(0);

}

//////////////////////////////////////////////////////////////////////////////
// find all domains and kabat strings for both heavy and light
int IgHMM::TraceDomains(P7_TRACE *tr, ESL_SQ *sq){

	m_Lkabatstr = ""; m_Hkabatstr = ""; m_kabatstr = "";
	m_HFW1.Clear(); m_HCDR1.Clear(); m_HFW2.Clear(); m_HCDR2.Clear(); 
	m_HFW3.Clear(); m_HCDR3.Clear(); m_HFW4.Clear();
	m_LFW1.Clear(); m_LCDR1.Clear(); m_LFW2.Clear(); m_LCDR2.Clear(); 
	m_LFW3.Clear(); m_LCDR3.Clear(); m_LFW4.Clear();
	m_MSA = "";

	GetMSA(sq, tr, m_MSA);

	int Loffset = 137;
	switch(m_hmmType){
		case Heavy:
			P5Pad(m_MSA, 1, 8); //First 8 of HFW1
			P3Pad(m_MSA, 115, 122); // Last 8 of HFW4
			TraceHMM(0,0,1,0);
		break;
		case Light:
			P5Pad(m_MSA, 1, 8); //First 8 of LFW1
			P3Pad(m_MSA, 101, 108); // Last 8 of LFW4
			TraceHMM(0,0,0,1);
		break;
		case ScFV:
			P5Pad(m_MSA, 1, 8); //First 8 of HFW1
			P3Pad(m_MSA, 115, 122); // Last 8 of HFW4
			P5Pad(m_MSA, Loffset+1, Loffset+8); //First 8 of LFW1
			P3Pad(m_MSA, Loffset+101, Loffset+108); // Last 8 of LFW4
			TraceHMM(0,Loffset,1,1);
	
		break;
		default:
		break;
	}

	m_HAA = m_HFW1.str + m_HCDR1.str + m_HFW2.str + m_HCDR2.str + 
		m_HFW3.str + m_HCDR3.str + m_HFW4.str;
	boost::replace_all(m_HAA,".","");
	boost::replace_all(m_HAA,"-","");
	StringToUpper(m_HAA);

	m_LAA = m_LFW1.str + m_LCDR1.str + m_LFW2.str + m_LCDR2.str + 
		m_LFW3.str + m_LCDR3.str + m_LFW4.str;
	boost::replace_all(m_LAA,".","");
	boost::replace_all(m_LAA,"-","");
	StringToUpper(m_LAA);


	std::string wholeORF = m_mergedAA;
	std::string::size_type loc;

	m_HDNA = "";
	m_HDNAstart = -1;
	if(m_HAA.size() > 0){
		loc = wholeORF.find(m_HAA);
		if(loc != std::string::npos){
			m_HDNA = m_mergedPairs;
			if(m_HDNA.size()>0){
				m_HDNAstart = MapAAToNT(m_HDNA.size(),m_HAA.size(), 
					loc, m_bestFrame1, 0);
				m_HDNA = m_HDNA.substr( m_HDNAstart,
					3*m_HAA.size());  //enough NT for the LAA

				GetDomainInterval(m_HFW1, loc, m_bestFrame1);
				GetDomainInterval(m_HCDR1, loc, m_bestFrame1);
				GetDomainInterval(m_HFW2, loc, m_bestFrame1);
				GetDomainInterval(m_HCDR2, loc, m_bestFrame1);
				GetDomainInterval(m_HFW3, loc, m_bestFrame1);
				GetDomainInterval(m_HCDR3, loc, m_bestFrame1);
				GetDomainInterval(m_HFW4, loc, m_bestFrame1);

				SanCheckDomInterval(m_HFW1, m_mergedPairs);
				SanCheckDomInterval(m_HCDR1, m_mergedPairs);
				SanCheckDomInterval(m_HFW2, m_mergedPairs);
				SanCheckDomInterval(m_HCDR2, m_mergedPairs);
				SanCheckDomInterval(m_HFW3, m_mergedPairs);
				SanCheckDomInterval(m_HCDR3, m_mergedPairs);
				SanCheckDomInterval(m_HFW4, m_mergedPairs);

			}
		}else{
			fprintf(stderr, "ERROR: cannot find HAA in ORF\n");
			exit(-1);
		}
	}

	m_LDNA = "";
	m_LDNAstart = -1;
	if(m_LAA.size() > 0){
		loc = wholeORF.find(m_LAA);
		if(loc != std::string::npos){
			m_LDNA = m_mergedPairs;
			if(m_LDNA.size()){
				m_LDNAstart = MapAAToNT(m_LDNA.size(),m_LAA.size(), 
						loc, m_bestFrame1, 0); 
				m_LDNA = m_LDNA.substr(m_LDNAstart,
					3*m_LAA.size()); //enough NT for the LAA

				GetDomainInterval(m_LFW1, loc, m_bestFrame1);
				GetDomainInterval(m_LCDR1, loc, m_bestFrame1);
				GetDomainInterval(m_LFW2, loc, m_bestFrame1);
				GetDomainInterval(m_LCDR2, loc, m_bestFrame1);
				GetDomainInterval(m_LFW3, loc, m_bestFrame1);
				GetDomainInterval(m_LCDR3, loc, m_bestFrame1);
				GetDomainInterval(m_LFW4, loc, m_bestFrame1);

				SanCheckDomInterval(m_LFW1, m_mergedPairs);
				SanCheckDomInterval(m_LCDR1, m_mergedPairs);
				SanCheckDomInterval(m_LFW2, m_mergedPairs);
				SanCheckDomInterval(m_LCDR2, m_mergedPairs);
				SanCheckDomInterval(m_LFW3, m_mergedPairs);
				SanCheckDomInterval(m_LCDR3, m_mergedPairs);
				SanCheckDomInterval(m_LFW4, m_mergedPairs);
			}
		}else{
			fprintf(stderr, "ERROR: cannot find LAA in ORF\n");
			exit(-1);
		}
	}

	return 0;
}


//////////////////////////////////////////////////////////////////////////////

void TestGetMatchStates(){
	
	//                 12345678901234567
	//                    123456 78901
	std::string msa = "abcDEFHIJkLMONPq";
	std::string dom = "";

	dom = ""; GetMatchStates(msa,1,3,0,0,dom,0);
	printf("GetMatchStates(msa,1,3,0,0,dom); \n %s\n\n",dom.c_str()); 

	dom = ""; GetMatchStates(msa,1,3,0,1,dom,0);
	printf("GetMatchStates(msa,1,3,0,1,dom); \n %s\n\n",dom.c_str()); 

	dom = ""; GetMatchStates(msa,1,3,1,0,dom,0);
	printf("GetMatchStates(msa,1,3,1,0,dom); \n %s\n\n",dom.c_str()); 

	dom = ""; GetMatchStates(msa,1,3,1,1,dom,0);
	printf("GetMatchStates(msa,1,3,1,1,dom); \n %s\n\n",dom.c_str()); 

	dom = ""; GetMatchStates(msa,7,11,0,0,dom,0);
	printf("GetMatchStates(msa,7,11,0,0,dom); \n %s\n\n",dom.c_str()); 

	dom = ""; GetMatchStates(msa,7,11,0,1,dom,0);
	printf("GetMatchStates(msa,7,11,0,1,dom); \n %s\n\n",dom.c_str()); 

	dom = ""; GetMatchStates(msa,7,11,1,0,dom,0);
	printf("GetMatchStates(msa,7,11,1,0,dom); \n %s\n\n",dom.c_str()); 

	dom = ""; GetMatchStates(msa,7,11,1,1,dom,0);
	printf("GetMatchStates(msa,7,11,1,1,dom); \n %s\n\n",dom.c_str()); 


	exit(0);
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////// 
