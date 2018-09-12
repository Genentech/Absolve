
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//	The BWAreader object is based on BWA's tools for reading 
//	fastq files, both compressed and uncompressed and single and paired-end
//	reads.  
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

extern "C" {

#include<bwa.h>
#include<bwamem.h>
#include<kvec.h>
#include<utils.h>
#include<bntseq.h>
#include<kseq.h>
KSEQ_DECLARE(gzFile)

}



class BWAReader {
public:
	BWAReader();
	~BWAReader();

	//r2file is optional
	int OpenSeqFiles(std::string r1file, std::string r2file);
	int ReadNSeqs(int n);
	int FreeSeqs(); // after done with ReadNSeqs, free them

	// force sequences to be upper case.
	int UpperCase();

	int ReadNSeqR1(int n, bseq1_t * seqs, int & nRead);

	int SkipFirstNSeqs(int n);

	bseq1_t * m_bseqs;
	int m_Nbseqs;

private:

	kseq_t * m_k1;
	kseq_t * m_k2;

	gzFile m_fp;
	int m_fd;
	gzFile m_fp2;
	int m_fd2;

};

 
