
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: MIT, see included COPYING file
//
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <fcntl.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <unistd.h>
#include <getopt.h>

#include <string>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include <boost/lexical_cast.hpp>

#include "util.h"
#include "ighmm.h"


#include "bwaaligner.h"

#include "ssw_cpp.h"

////////////////////////////////////////////////////////////////////////////
extern "C" {

extern void mem_process_seqs_nonthreaded(const mem_opt_t *opt, const bwt_t *bwt, const bntseq_t *bns, const uint8_t *pac, int64_t n_processed, int n, bseq1_t *seqs, const mem_pestat_t *pes0);
	
}

////////////////////////////////////////////////////////////////////////////

BWAAligner::BWAAligner(){
	m_opt = mem_opt_init();
	m_idx = NULL;
	m_useSSW = 0; //don't use ssw for classification by default
	m_byAA = 0; // don't classify by AA be default
	m_multihit = 0;
}
////////////////////////////////////////////////////////////////////////////
BWAAligner::~BWAAligner(){
	if(m_opt != NULL) { free(m_opt);}
	if(m_idx != NULL) { bwa_idx_destroy(m_idx);}
}

////////////////////////////////////////////////////////////////////////////
BWAAligner::BWAAligner(const char *idxFile, const char *name){
	m_opt = mem_opt_init();
	m_useSSW = 0; //don't use ssw for classification
	m_byAA = 0; // don't classify by AA be defaul
	m_multihit = 0;

	Initialize(idxFile, name);
}

////////////////////////////////////////////////////////////////////////////
int BWAAligner::SetMultiHit(int multihit){
	m_multihit = multihit;
	if(multihit == 0) return 0;

	m_opt->max_XA_hits = 1000; 
	m_opt->max_XA_hits_alt = 1000;
	m_opt->XA_drop_ratio = 0.2;

	return 0;
}

////////////////////////////////////////////////////////////////////////////
int BWAAligner::Initialize(const char *idxFile, const char * name,
	BWAalignConfig & config){

	return Initialize(idxFile, name,
		config.A, config.B,
		config.OINS, config.EINS,
		config.ODEL,config.EDEL,
		config.CLIP3, config.CLIP5,
		config.MIN_SEED_LEN, config.MIN_SCORE);

}



////////////////////////////////////////////////////////////////////////////
int BWAAligner::Initialize(const char *idxFile, const char * name,
	int a, int b, 
	int o_ins, int e_ins,
	int o_del, int e_del,
	int pen_clip3, int pen_clip5,
	int min_seed_len,
	int T
	){

	m_name = name;
	if(idxFile == NULL) {return -1;}
	std::string idxFileStr = idxFile;
	if(idxFileStr.compare("") == 0){ return -1;}

	m_opt->a = a;
	m_opt->b = b;
	m_opt->o_ins = o_ins;
	m_opt->e_ins = e_ins;
	m_opt->o_del = o_del;
	m_opt->e_del = e_del;
	m_opt->pen_clip3 = pen_clip3;
	m_opt->pen_clip5 = pen_clip5;
	m_opt->min_seed_len = min_seed_len;
	m_opt->T = T;
	m_opt->n_threads = 1;


	bwa_fill_scmat(m_opt->a, m_opt->b, m_opt->mat);

	m_idx = bwa_idx_load_from_disk(idxFile,BWA_IDX_ALL);

	return 0;
}
////////////////////////////////////////////////////////////////////////////
int BWAAligner::Classify(bseq1_t *tseq,
	std::string & call, std::string & pos, std::string & cigar,
	std::string & score
){

	mem_process_seqs_nonthreaded(m_opt, m_idx->bwt, m_idx->bns, m_idx->pac, 0, 1, tseq, NULL);
	GetSAMDetails(tseq->sam, call, pos, cigar, score);

	int curMinSeed = m_opt->min_seed_len;
	// failed to get a hit, keep decreasing the seed size to
	// attempt to find some hit
	// 
	// NB this has a significant time cost, 2 seconds vs 12 seconds 
	// for 344 seqs,
	int init = 1;
	mem_opt_t *relax_opt = NULL;
	while(strncmp(call.c_str(), "*",1) == 0){
		if(curMinSeed < 6){ break;}

		if(init){
			relax_opt = mem_opt_init();
			bwa_fill_scmat(relax_opt->a, relax_opt->b, relax_opt->mat);
			init = 0;
		}
		relax_opt->a = m_opt->a;
		relax_opt->b = m_opt->b;
		relax_opt->o_ins = m_opt->o_ins;
		relax_opt->e_ins = m_opt->e_ins;
		relax_opt->o_del = m_opt->o_del;
		relax_opt->e_del = m_opt->e_del;
		relax_opt->pen_clip3 = m_opt->pen_clip3;
		relax_opt->pen_clip5 = m_opt->pen_clip5;
		relax_opt->T = m_opt->T;
		m_opt->n_threads = 1;

		curMinSeed = 5; //force only two tries
		relax_opt->min_seed_len = curMinSeed;


		mem_process_seqs_nonthreaded(relax_opt, m_idx->bwt, 
			m_idx->bns, m_idx->pac, 0, 1, tseq, NULL);

		GetSAMDetails(tseq->sam, call, pos, cigar, score);
	}

	if(relax_opt != NULL) { free(relax_opt);}

	return 0;
}

////////////////////////////////////////////////////////////////////////////
int BWAAligner::GetSAMDetails(char *sam,
	std::string & rname, std::string & pos, std::string & cigar, 
	std::string & score){
	//QNAME   FLAG    RNAME   POS MAPQ    CIGAR
	
	const char delim[2] = "\t";
	char *token;
	token = strtok(sam, delim);

	std::string xa;
	int field = 0;
	while( token != NULL ) 
	{
		field++;
		if(field == 3){ rname = token; }

		if(field == 4){ pos = token;}
		if(field == 6){ cigar = token; }
		if(field == 14){ 
			score = token;

			size_t prefixIndex = score.find("AS:i:");
			if(prefixIndex != std::string::npos){
				score = score.substr(prefixIndex+5);
			}else{
				score = "";
			}
		}
		if(field == 16) { xa = token; }
		token = strtok(NULL, delim);
		if(field > 16) { break;}
	}

	if(m_multihit){
		//parse the XA field
		std::istringstream xas(xa);
		std::string s;    
		std::string xarn;
		std::string xapos;
		std::string xacigar;
		std::string xascore;

		int tindex = 0;
		while (std::getline(xas, s, ';')) {
			if(tindex == 0){ //trim the "XA:Z:" leader
				size_t xa_z_index = s.find("XA:Z:");
				if(xa_z_index == std::string::npos){
					break;
				}else{
					s = s.substr(xa_z_index+5);
				}
			}
			xarn = s.substr(0,s.find_first_of(','));
			rname = rname + "," + xarn;

			s = s.substr(1 + s.find_first_of(','));

			xapos = s.substr(0,s.find_first_of(','));
			pos = pos + "," + xapos;

			s = s.substr(1 + s.find_first_of(','));
			xacigar = s.substr(0,s.find_first_of(','));

			cigar = cigar + "," + xacigar;

			s = s.substr(1 + s.find_first_of(','));
			xascore = s;

			score = score + "," + xascore;

			tindex++;
			if(tindex >= m_multihit) break;
		}

	}
	return 0;
}
////////////////////////////////////////////////////////////////////////////

typedef struct {
	int sw_score;
} AlignScoreParam;

////////////////////////////////////////////////////////////////////////////
int IsBetterScore(StripedSmithWaterman::Alignment & alignment,
	AlignScoreParam & someScore){

	if(alignment.sw_score > someScore.sw_score){
		return 1;
	}else{
		return 0;
	}
}

////////////////////////////////////////////////////////////////////////////
int BWAAligner::FindBestAlignment(const char *seq, int seqlen,
	std::string & call, std::string & score, 
	std::string & pos, std::string & cigar,
	int firstpos, 
	int lastpos
	){
	call = "";
	pos = "";
	cigar = "";

	bseq1_t *tseq = NULL;
	if(firstpos==0 && lastpos == -1){
		tseq = CreateBseq(seq, seqlen);
	}else if(firstpos>0 && lastpos == -1){
		tseq = CreateBseq(seq+firstpos, seqlen-firstpos);
	}else{
		tseq = CreateBseq(seq+firstpos, lastpos-firstpos+1);
	}
	Classify(tseq,call, pos, cigar, score);
	
	FreeBWASeq(&tseq);

	return 0;
}

////////////////////////////////////////////////////////////////////////////
bseq1_t * BWAAligner::CreateBseq( const char *seq, int len){
	bseq1_t *tseq = (bseq1_t*) malloc(sizeof(bseq1_t));
	
	memset(tseq, 0, sizeof(bseq1_t));
	tseq->name = (char *) malloc(ABS_MAX_NT_LEN);
	tseq->seq = (char *) malloc(ABS_MAX_NT_LEN);
	tseq->comment = (char *) malloc(ABS_MAX_NT_LEN);
	tseq->qual = NULL; 
	tseq->sam = NULL;

	strcpy(tseq->name, "samp");
	
	strncpy(tseq->seq, seq, len);

	tseq->l_seq = len;
	tseq->seq[tseq->l_seq] = '\0';

	return tseq;
}

////////////////////////////////////////////////////////////////////////////
int BWAAligner::FreeBWASeq(bseq1_t **seq){

    free((*seq)->name);
    free((*seq)->seq);
    free((*seq)->comment);
    free((*seq)->qual);
    free((*seq)->sam);
	free((*seq));
	*seq = NULL;
	return 0;
}




////////////////////////////////////////////////////////////////////////////

int AllocBWASeq(bseq1_t **tseq){
	*tseq = (bseq1_t*) malloc(sizeof(bseq1_t));
	memset(*tseq, 0, sizeof(bseq1_t));
	(*tseq)->name = (char *) malloc(ABS_MAX_NT_LEN);
	(*tseq)->seq = (char *) malloc(ABS_MAX_NT_LEN);
	(*tseq)->comment = (char *) malloc(ABS_MAX_NT_LEN);
	(*tseq)->qual = (char *) malloc(ABS_MAX_NT_LEN);
	(*tseq)->sam = NULL;
	return 0;
}

////////////////////////////////////////////////////////////////////////////
int FreeBWASeq(bseq1_t **seq){

    free((*seq)->name);
    free((*seq)->seq);
    free((*seq)->comment);
    free((*seq)->qual);
	free((*seq));
	return 0;
}
////////////////////////////////////////////////////////////////////////////

int InitBWASeq(const char *seqname, int namelen, const char *seq, int seqlen, bseq1_t *tseq){

	return 0;
}

int SetSeq(bseq1_t *tseq,
	const char *sseq,
	int sseqlen,
	const char *seqname,
	const char *squal
	){

	if(sseq != NULL){
		strncpy(tseq->seq, sseq, ABS_MAX_NT_LEN-1);
		tseq->seq[ABS_MAX_NT_LEN-1] = '\0';
		tseq->l_seq = (sseqlen > ABS_MAX_NT_LEN) ? ABS_MAX_NT_LEN :  sseqlen;
	}

	if(seqname != NULL){
		strncpy(tseq->name, seqname, ABS_MAX_NT_LEN-1);
		tseq->name[ABS_MAX_NT_LEN-1] = '\0';
	}

	if(squal != NULL){
		strncpy(tseq->qual, squal, ABS_MAX_NT_LEN-1);
		tseq->qual[ABS_MAX_NT_LEN-1] = '\0';
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////////
int TestBWA(){

	const char *idxFile = "./dbs/species/zoo/Vh.fa";

	bwaidx_t *idx = NULL;
	
	mem_opt_t *opt = mem_opt_init();
	bwa_fill_scmat(opt->a, opt->b, opt->mat);

	idx = bwa_idx_load_from_disk(idxFile,BWA_IDX_ALL);

	bseq1_t *seq = NULL;
	AllocBWASeq(&seq);

	const char *aseq = "GACTACCGTAAAAAGGGTAACTAGAGGTTGAGGTGAT";
	
	int numiters = 3;
	for(int iter=0;iter<numiters;iter++){
		SetSeq(seq, 
			aseq, strlen(aseq),
			"myseq",NULL);

		mem_process_seqs_nonthreaded(opt, idx->bwt, idx->bns, idx->pac, 0, 1, seq, NULL);
		free(seq->sam);
	}

	printf("numiters:%d\n",numiters);

    free(opt);
	bwa_idx_destroy(idx);

	FreeBWASeq(&seq);

	return 0;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

 
