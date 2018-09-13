
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: MIT, see included COPYING file
//
//
//	The pipeline object is the master object that orchestrates the absolve
//	pipeline.  The pipeline object handles:
//		* Program option processing
//		* Parameter initialization
//		* Parameter initialization
//		* Pipeline execution:
//			* Reading sequences from input files
//			* Running the HMM to annotate domains
//			* Running the aligners to perform germline classification
//			* Tracking statistics
//			* Writing results to output files
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#pragma once


#include <string>
#include <vector>
#include <map>

#include "fasta.h"
#include "ighmm.h"
#include "aligner.h"

#include "result.h"
#include "runstats.h"
/////////////////////////////////////////////////////////
//
enum BreakoutLineage {
	H,
	L,
	NONE
};

/////////////////////////////////////////////////////////
typedef struct {
	double ExhaustiveHMMThresh;
	double FilterHMMThresh;
	std::string HeavyHMMFile;
	std::string LightHMMFile;
	std::string ScFVHMMFile;
	std::string TempFolderForDistrib;

} AbsolveRC;

/////////////////////////////////////////////////////////

class AbsolvePipeline {
public:
	AbsolvePipeline();
	~AbsolvePipeline();

	int ProcessOptions(int argc, char *argv[], char **env, const char *exepath);
	int InitPipeline();
	int Run();

private:

	int ReadRCFile(const char *inifilename);
	void LoadInitHMM(std::string & hmmfile, IgHMM & model, HMMType type);
	void ShowHelpExit();

	int SetAligner(std::string & alignerType);
	int SetMultiHit(int multihit);

	// pipeline entry points
	int RunFastQPipeline(); // fastq DNA, paired or single
	int RunFastaPipeline(); // DNA fasta file
	int RunAAPipeline(); // AA fasta file

	// per query sequence routines
	int RunFastaSeq(char *name, char *seq, int len);
	int RunFastaSeqExhaustive(char *name, char *seq, int len);
	int RunAASeq(std::string & name, std::string & seq);
	int RunAASeqExhaustive(std::string & name, std::string & seq);
	int SeqResult(const char *seq, int len, Result & rslt);
	
	// classifiers
	int RunAllClassifiers(const char *seq, int len, Result & rslt);
	int RunAllHeavyClassifiers(const char *seq, int len, Result & rslt);
	int RunAllLightClassifiers(const char *seq, int len, Result & rslt);

	// keep a running tally of all mutations by kabat position
	int TallyMutKabat(Result & rslt, char chain, std::vector<int> & mutpos);

	int WriteResult(Result & result);
	int AppendResults();

	//Determine if a sequence has large deletions
	int IncompleteSequence(Result & rslt);

	std::string GetLineage(std::string & VH, std::string & JH,
		std::string & CDR3);

	int GetOffsetFromCigar(std::string & cigar);

private:
	IgHMM m_ighmmSCFV;
	IgHMM m_ighmmHeavy;
	IgHMM m_ighmmLight;

	IgHMM * m_ighmmDefault;

	std::string m_hmmfileSCFV;
	std::string m_hmmfileHeavy;
	std::string m_hmmfileLight;


	BWAReader m_DNAreader; //
	FastaIO m_AAFastaReader;//

	Aligner CH;
	std::string m_CHdb;
	Aligner VH;
	std::string m_VHdb;
	Aligner JH;
	std::string m_JHdb;

	Aligner CL;
	std::string m_CLdb;
	Aligner VL;
	std::string m_VLdb;
	Aligner JL;
	std::string m_JLdb;

	//map of bwa aligner params read from rcfile
	std::map<std::string, BWAalignConfig> m_bwaconfigs;

	std::string m_R1file;
	std::string m_R2file;
	int m_inputIsAA;

	int m_collectStats;
	int m_profileTiming;


	HMMType m_hmmtype;

	std::string m_startseqArg;
	std::string m_stopseqArg;
	int m_startseq;
	int m_stopseq;

	int m_exhaustive; //Search exhaustively for domains
	double m_exhaustiveHMMThresh;

	// filtering flags
	double m_filterHMMThresh;
	int m_filterHMM;
	int m_filterNoCall;
	int m_filterIncomplete;

	// is this a worker sub process?
	int m_isDistrib;

	std::string m_outdir; // result folder
	std::string m_outdirTmp; // temporary result folder, used by isDistrib.

	std::map<std::string, FILE *> m_outmap;

	BreakoutLineage m_breakout;

	RunStats m_stats;

	std::string m_hostPID;
	std::string m_exepath;
	std::string m_tempfolder;

	int m_NT_PAIR_MIN;

	//use a copy of the seq so classfiers can mask after alignment
	char m_classifierBuf[ABS_MAX_NT_LEN];
};

 
