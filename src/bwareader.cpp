
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
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

#include "bwareader.h"


////////////////////////////////////////////////////////////////////////////
//


BWAReader::BWAReader(){
	m_fd = -1;
	m_fd2 = -1;
	m_bseqs = NULL;
	m_Nbseqs = 0;
	m_k1 = NULL;
	m_k2 = NULL;
	m_fp = NULL;
	m_fp2 = NULL;

}
BWAReader::~BWAReader(){
	if(m_k1 != NULL) kseq_destroy(m_k1);
	if(m_k2 != NULL) kseq_destroy(m_k2);
	if(m_fp != NULL) gzclose(m_fp);
	if(m_fp2 != NULL) gzclose(m_fp2);
	if(m_bseqs != NULL) FreeSeqs();

}
////////////////////////////////////////////////////////////////////////////


int BWAReader::OpenSeqFiles(std::string r1file, std::string r2file){
	if(r1file.empty()){
		fprintf(stderr,"R1 file must be specified\n");

		return -1;
	}
	m_fp = m_fp2 = NULL;
	m_k1 = NULL;
	m_k2 = NULL;

	m_fd = open(r1file.c_str(), O_RDONLY);
	m_fp = gzdopen(m_fd, "r");
	m_k1 = kseq_init(m_fp);

    if(!r2file.empty()){
        m_fd2 = open(r2file.c_str(), O_RDONLY);
        m_fp2 = gzdopen(m_fd2, "r");
		m_k2 = kseq_init(m_fp2);
		return 2;
    }
	return 1;
}
////////////////////////////////////////////////////////////////////////////
int BWAReader::SkipFirstNSeqs(int n){

	int nRead  = 0;
	fprintf(stderr,"Skipping %d reads...\n",n);
	if(n>0){
		while(kseq_read(m_k1) >= 0){
			if(m_k2 != NULL) kseq_read(m_k2);
			nRead++;

			if(nRead>=n) break;
		}
	}
	fprintf(stderr,"Skipping %d reads...DONE\n",n);
	return nRead;

}

////////////////////////////////////////////////////////////////////////////

static inline void trim_readno(kstring_t *s)
{
    if (s->l > 2 && s->s[s->l-2] == '/' && isdigit(s->s[s->l-1]))
        s->l -= 2, s->s[s->l] = 0;
}
void CopyKseq2Bseq(const kseq_t *ks, bseq1_t *s){
	strncpy(s->name, ks->name.s, ABS_MAX_NT_LEN);
	s->name[ABS_MAX_NT_LEN-1] = '\0';
	s->comment = 0;
	strncpy(s->seq, ks->seq.s, ABS_MAX_NT_LEN);
	if(ks->qual.l){
		strncpy(s->qual, ks->qual.s, ABS_MAX_NT_LEN);
		s->qual[ABS_MAX_NT_LEN-1] = '\0';
	}
	int rlen = strlen(ks->seq.s);
	if(rlen >= ABS_MAX_NT_LEN){
		rlen = ABS_MAX_NT_LEN - 1;
	}
	s->l_seq = rlen;
}

int BWAReader::ReadNSeqR1(int n, bseq1_t * seqs, int & nRead){
	nRead = 0;
	while(kseq_read(m_k1) >=0 ){
		trim_readno(&m_k1->name);
		nRead++;
		if(nRead >= n) break;
	}
	return nRead;
}
////////////////////////////////////////////////////////////////////////////
int BWAReader::UpperCase(){
	for(int i=0;i<m_Nbseqs;i++){
        if(m_bseqs[i].seq != NULL) ToUpper(m_bseqs[i].seq, m_bseqs[i].l_seq);
	} 
	return 0;
}
////////////////////////////////////////////////////////////////////////////
int BWAReader::ReadNSeqs(int n){

	m_bseqs = bseq_read(n, &m_Nbseqs, m_k1, m_k2);
	UpperCase();

	return m_Nbseqs;
}

////////////////////////////////////////////////////////////////////////////
int BWAReader::FreeSeqs(){
	for(int i=0;i<m_Nbseqs;i++){
        if(m_bseqs[i].name != NULL) free(m_bseqs[i].name);
        if(m_bseqs[i].comment != NULL) free(m_bseqs[i].comment);
        if(m_bseqs[i].seq != NULL) free(m_bseqs[i].seq);
        if(m_bseqs[i].qual != NULL) free(m_bseqs[i].qual);

		m_bseqs[i].name = NULL;
		m_bseqs[i].comment = NULL;
		m_bseqs[i].seq = NULL;
		m_bseqs[i].qual = NULL;
    }
    free(m_bseqs);
	m_bseqs = NULL;
	return 0;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// 
