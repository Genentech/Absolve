
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <string.h>
#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>


#include "aligner.h"

const int MaxLibrarySize = 1000000;

////////////////////////////////////////////////////////////////////////////
Aligner::Aligner(){
	m_useSSW = 0;
	m_byAA = 0;
}

////////////////////////////////////////////////////////////////////////////
Aligner::Aligner(const char *idxFile, const char* name){

	m_useSSW = 0;
	m_byAA = 0;
	m_bwaAligner.Initialize(idxFile, name);
}

////////////////////////////////////////////////////////////////////////////
int Aligner::Initialize(const char *idxFile, const char * name,
	BWAalignConfig & config){

	return Initialize(idxFile, name,
		config.A, config.B,
		config.OINS, config.EINS,
		config.ODEL,config.EDEL,
		config.CLIP3, config.CLIP5,
		config.MIN_SEED_LEN, config.MIN_SCORE);

}



////////////////////////////////////////////////////////////////////////////
int Aligner::Initialize(const char *idxFile, const char * name,
	int a, int b, 
	int o_ins, int e_ins,
	int o_del, int e_del,
	int pen_clip3, int pen_clip5,
	int min_seed_len,
	int T
	){

	if(idxFile == NULL) {return -1;}
	std::string idxFileStr = idxFile;
	if(idxFileStr.compare("") == 0){ return -1;}

	m_bwaAligner.Initialize(idxFile, name,
		a, b, 
		o_ins, e_ins,
		o_del, e_del,
		pen_clip3, pen_clip5,
		min_seed_len,
		T
	);

	// read in all the seqs
	m_library.OpenSeqFiles(idxFile,"");
	int libsize = m_library.ReadNSeqs(MaxLibrarySize); 
	for(int i=0;i<libsize; i++){
		std::string libname = m_library.m_bseqs[i].name;
		m_libraryEntries[libname] = i;
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////
Aligner::~Aligner(){
}
////////////////////////////////////////////////////////////////////////////
int Aligner::SetMultiHit(int multihit){
	m_bwaAligner.SetMultiHit(multihit);
	m_sswAligner.SetMultiHit(multihit);
	return 0;
}
////////////////////////////////////////////////////////////////////////////
void PrettyPrintAlignment(std::string outref, std::string outq){
	std::string matches = std::string(outref.size(),' ');
	int match = 0;
	int mismatch = 0;
	for(unsigned int i=0; i<outref.size(); i++){
		if(outref[i] == outq[i]){
			matches[i] = '|';	
			match++;
		}else{
			mismatch++;
		}
	}

	printf("\tR/Q\n\t%s\n\t%s\n\t%s\n\tMatch:%d\n\tMisMa:%d\n",
		outref.c_str(), matches.c_str(), outq.c_str(), match, mismatch);
}
////////////////////////////////////////////////////////////////////////////
// generate an alignment for ref and query from the cigar string
void GetAlignment(std::string ref, std::string cigar, 
	int pos, std::string q,  std::string & outref, std::string & outq){

	pos = pos - 1; //convert to zero-based
	if(pos<0){fprintf(stderr,"ERROR:Bad pos:%d\n",pos); exit(-1);}

	std::vector<char> cigM;
	std::vector<int> cigI;

	char numbuf[32];
	int numbufi = 0;
	for(unsigned int c=0;c<cigar.size();c++){
		if(isdigit(cigar[c])){
			numbuf[numbufi++] = cigar[c];
		}else{
			numbuf[numbufi] = '\0';
			numbufi = 0;
			cigM.push_back(cigar[c]);
			cigI.push_back(atoi(numbuf));
		}
	}

	outref = "";
	outq = "";
	unsigned int refcursor = 0;
	unsigned int qcursor = 0;
	std::string refadd;
	std::string qadd;
	for(unsigned int i = 0; i< cigM.size(); i++){
		switch(cigM[i]){
			case 'S':
				if(i==0){ // 5' softclip
					if(pos>cigI[i]){ // ref offset bigger
						refadd = ref.substr(0,pos);
						boost::algorithm::to_lower(refadd);
						outref += refadd;

						qadd = std::string(pos-cigI[i],'-') +
							q.substr(0,cigI[i]);
						boost::algorithm::to_lower(qadd);
						outq += qadd;
						refcursor += pos;
						qcursor += cigI[i];
					}else{
						refadd = std::string(cigI[i]-pos,'-') +
							ref.substr(0,pos);
						boost::algorithm::to_lower(refadd);
						outref += refadd;
						qadd = q.substr(0,cigI[i]);
						boost::algorithm::to_lower(qadd);
						outq += qadd;
						refcursor += pos;
						qcursor += cigI[i];
					}
				}else{ // 3' softclip
					if(ref.substr(refcursor).size() > //ref is longer
						q.substr(qcursor).size()){
						
						refadd = ref.substr(refcursor);
						boost::algorithm::to_lower(refadd);
						outref += refadd;
						
						qadd = q.substr(qcursor) + 
							std::string(ref.substr(refcursor).size() -
								q.substr(qcursor).size(),'-');
	
						boost::algorithm::to_lower(qadd);
						outq += qadd;
					}else{
						refadd = ref.substr(refcursor) +
							std::string(q.substr(qcursor).size() -
								ref.substr(refcursor).size(),'-' );
						boost::algorithm::to_lower(refadd);
						outref += refadd;
	
						qadd = q.substr(qcursor);
						boost::algorithm::to_lower(qadd);
						outq += qadd;
					}

					refcursor = ref.size();
					qcursor = q.size();
				}
				assert(outref.size() == outq.size());
			break;
			case 'M':
			case 'X':
			case '=':
				if(i==0){
					outref += ref.substr(refcursor, pos+cigI[i]);
					outq += std::string(pos,'-') + q.substr(qcursor, cigI[i]);
					refcursor += pos + cigI[i];
				}else{
					outref += ref.substr(refcursor, cigI[i]);
					outq += q.substr(qcursor, cigI[i]);
					refcursor += cigI[i];
				}
				qcursor += cigI[i];
				assert(outref.size() == outq.size());
			break;
			case 'I':
				if(i==0){
					outref += ref.substr(0,pos) + std::string(cigI[i],'-');
					refcursor = pos;
					outq += std::string(pos,'-') + q.substr(qcursor, cigI[i]);
					qcursor += cigI[i];
				}else{
					outref += std::string(cigI[i],'-');
					outq += q.substr(qcursor, cigI[i]);
					qcursor += cigI[i];
				}
				assert(outref.size() == outq.size());
			break;
			case 'D':
				if(i==0){
					outq += std::string(pos,'-') + std::string(cigI[i],'-');
					qcursor += 0;
					outref += ref.substr(0,pos) + ref.substr(pos+1, cigI[i]);
					refcursor += pos + cigI[i];
					
				}else{
					outref += ref.substr(refcursor, cigI[i]);
					outq += std::string(cigI[i],'-');
					refcursor += cigI[i];
				}

				assert(outref.size() == outq.size());
			break;
			default:
				fprintf(stderr, "ERROR: unhandled cigar type:%c\n",cigM[i]);
				exit(-1);
			break;
		}
	}

	//fill in the last gaps
	if(refcursor < ref.size() || qcursor < q.size()){
		// C- -- CC CC
		// GG GG G- --
		outq += outq.substr(qcursor);
		outref += ref.substr(refcursor);
		if(outq.size() > outref.size()){
			outref += std::string(outq.size()-outref.size(), '-');
		}else if(outref.size() > outq.size()){
			outq += std::string(outref.size()-outq.size(), '-');
		}
	}

	if(outref.size() != outq.size()){
		fprintf(stderr, "ERROR:alignment string sizes do not match\n");
		exit(-1);
	}

}

////////////////////////////////////////////////////////////////////////////
// compute the SHM and their positions for query against a call
int Aligner::ComputeSHM(std::string call, 
	const char *seq, int seqlen, 
	std::string & cigar,
	std::string & shm,
	std::vector<int> & mutpos){

	if(!Ready()){ return -1; }

	std::vector<std::string> allcalls;
	std::istringstream callstream(call);
	std::string onecall;
	int callindex = 0;
	shm = "";
	while(getline(callstream, onecall, ',')){

		if(m_libraryEntries.find(onecall) == m_libraryEntries.end()){
			return -1;
		}

		int index = m_libraryEntries[onecall];

		char *rseq = m_library.m_bseqs[index].seq;
		
		std::string query = seq;
		std::string ref = rseq;

		std::string sswcigar;
		int rpos, mismatches;
		double score;

		m_sswAligner.Align(
			query.c_str(), query.size(), ref.c_str(), ref.size(), 
			sswcigar, rpos, score, mismatches, m_byAA);

		std::string qalign;
		std::string ralign;

		GetAlignment(rseq, sswcigar, rpos, seq, ralign, qalign);

		//Find the first non "-" 5p char in ref
		//       |
		//   ----aagtTTTTTTT  # ref
		//   gjksacccTTTTTTT  # seq
		//   |
		//   aagtgggTTTTTTT  # ref
		//   ---acccTTTTTTT  # seq
		unsigned int startref = 0;
		int inFlank = 1;
		for(unsigned int i=0;i<ralign.size();i++){
			if(ralign[i] != '-') inFlank = 0;
			if(ralign[i] == '-' && inFlank) startref = i+1;
			if(!inFlank) break;
		}

		//Find the the first non '-' 3p char in ref
		unsigned int endref = ralign.size()-1;
		inFlank = 1;
		for(int i=endref;i>=0; i--){
			if(ralign[i] != '-') inFlank = 0;
			if(ralign[i] == '-' && inFlank) endref = i-1;
			if(!inFlank) break;
		}

		//map each align pos to its position in the query
		//and then mask it
		std::vector<int> querypos;
		int curQpos = 0;
		for(unsigned int i=0;i<qalign.size();i++){
			querypos.push_back(curQpos);
			if(qalign[i] != '-'){
				curQpos++;
				if(curQpos >= seqlen){ curQpos = seqlen-1;}
			}
		}

		int shmscore = 0;
		for(unsigned int i=startref;i<=endref;i++){
			if(ralign[i] != qalign[i]){
				if(callindex == 0){ //only tally mutpos for first call
					mutpos.push_back(querypos[i]);
				}
				shmscore++; 
			}
		}

		if(callindex>0){
			shm = shm + ",";
		}
		shm = shm + boost::lexical_cast<std::string>( shmscore  );
		callindex++;

	}
	return 0;
}
////////////////////////////////////////////////////////////////////////////
int Aligner::Ready(){
	return(m_library.m_Nbseqs > 0);
}

////////////////////////////////////////////////////////////////////////////
//
int Aligner::RunClassifier(const char *seq, int seqlen,
		std::string & call, std::string & score,
		std::string & pos, std::string & cigar,
		int firstpos, // = 0,
		int lastpos // = -1,
	){
	call = "";
	pos = "";
	cigar = "";
	
	if(Ready()){
		if(m_useSSW || m_byAA){
			m_sswAligner.FindBestAlignment(seq, seqlen,
				call, score, pos, cigar, m_byAA, m_library, firstpos, lastpos);
		}else{
			m_bwaAligner.FindBestAlignment(seq, seqlen, 
				call, score, pos, cigar, firstpos, lastpos);
		}
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////// 
