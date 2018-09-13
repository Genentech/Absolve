
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: MIT, see included COPYING file
//
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <sys/types.h>
#include <unistd.h>

#include <vector>
#include <map>
#include <algorithm>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <fstream>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/ini_parser.hpp>


#include "util.h"

#include "pipeline.h"


///////////////////////////////////////////////////////////
enum SeqFileType {
	sqFASTQ,
	sqFASTA,
	sqUNK
};


///////////////////////////////////////////////////////////

void TestGzipIo(){
	std::ifstream file("hello.gz", std::ios_base::in | std::ios_base::binary);
	boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
	in.push(boost::iostreams::gzip_decompressor());
	in.push(file);

}


///////////////////////////////////////////////////////////

SeqFileType GetSeqFileType(const char *fname){
	boost::regex fqgz(".*fastq.gz");
	boost::regex fq(".*fastq");
	boost::regex fa(".*fa");
	boost::regex fasta(".*fasta");
	boost::regex fna(".*fna");

	if(boost::regex_match(fname, fqgz)) return sqFASTQ;
	if(boost::regex_match(fname, fq)) return sqFASTQ;
	if(boost::regex_match(fname, fa)) return sqFASTA;
	if(boost::regex_match(fname, fasta)) return sqFASTA;
	if(boost::regex_match(fname, fna)) return sqFASTA;

	return sqUNK;

}

///////////////////////////////////////////////////////////////////////////
AbsolvePipeline::AbsolvePipeline(){
	m_R1file = "";
	m_R2file = "";
	m_outdir = "";
	m_outdirTmp = "";
	//m_outfileP = NULL;
	m_startseq = -1;
	m_stopseq = -1;
	m_inputIsAA = 0;
	m_collectStats = 0;
	m_profileTiming = 0;
	m_isDistrib = 0;
	m_exhaustive = 0;
	m_exhaustiveHMMThresh = 70.0; //  55.0 picks ups partials

	m_hmmtype = ScFV; // scfv by default
	m_ighmmDefault = & m_ighmmSCFV;

	m_filterHMM = 0; // filter out by m_filterHMMThresh?
	m_filterHMMThresh = 70.0; //  55.0 picks ups partials

	m_filterNoCall = 0;

	m_filterIncomplete = 0; //flag to filter by incomplete

	char buf[4096];
	gethostname(buf,4096);

	m_hostPID = buf;
	m_hostPID += ".";
	m_hostPID += boost::lexical_cast<std::string>(getpid());
}
AbsolvePipeline::~AbsolvePipeline(){

	std::map<std::string, FILE *>::iterator it;
	for(it = m_outmap.begin(); it!=m_outmap.end(); ++it){
		FILE *ofile = (*it).second;
		if(ofile != NULL) fclose(ofile);
	}

	//clean the tmp dir
	if(m_isDistrib && boost::filesystem::exists(m_outdirTmp)){
		boost::filesystem::remove_all(m_outdirTmp);
	}
}


void AbsolvePipeline::ShowHelpExit(){
	const char *opts[] = {	
" --R1 seqs.(fasta|fastq|fastq.gz), sequences to process, R1 set if paired",
" [--R2 seqs.(fasta|fastq|fastq.gz)], sequences to process, optional R2 if paired",
" [--AA] flag indicating that sequences are amino acid, not DNA",
" [--stats] flag indicating whether statistics should be reported",
"    this will generate TSV files showing counts of germline usage,",
"    CDR3 frequencies, lineages, etc",
" [--profileTiming] flag indicating whether statistics on run times ",
"    should be reported in the stats.  Useful for performance profiling.",
"    will add 'RunTime' and 'AppendTime' to 'statRunStats.tsv'.",
" [--CH dbfile] heavy constant domain database, optional",
" [--VH dbfile] heavy v-segment database, optional",
" [--JH dbfile] heavy j-segment database, optional",
" [--CL dbfile] light constant domain database, optional",
" [--VL dbfile] light v-segment database, optional",
" [--JL dbfile] light j-segment database, optional",
"   all dbfiles should be BWA indices, usually created from a fasta file",
"     and 'bwa index fastafile' (except amino acid dbfiles)",
" [--hmmtype (Heavy| Light|ScFV) what HMM should be run ",
"     Light = just the light chain ",
"     Heavy = just the heavy chain ",
"     ScFV = run the scFv model (heavy followed by light with Gs linker) ", 
"     ScFV is the default ",
"     run times will be faster if you pick the chain you know you want ",
"     i.e. don't just use the default if you know better ",
" [--useAligner SWA|BWA] force the type of aligner to use",
"     By default, SWA will be used for AA and BWA for NT. ",
"     SWA is a full smithwaterman aligner, it works to NT and AA, gives ",
"       full verbose cigar strings and positions for all alignments.",
"     BWA is the bwa mem aligner,  it works only with NT and gives ",
"       non-verbose cigar strings and no cigar strings for alternative hits.",
" [--multihit N] return additional top N hits from database classifiers",
"     for N>0, additional hits are reported as a commma separated list e.g:",
"     'IGHV3-30*01,IGHV3-30*02,IGHV1-1*01' ",
" [--exhaustive] search exhaustively for all heavy and light chain instances ",
"    within a single query sequence.  This is useful for finding variable ",
"    domains within a polycistronic construct for example.  It will ",
"    potentially generate multiple records per accession and is not ",
"    generally used for NGS, scFv or simple antibody sequences.",
" [--exhaustiveThreshold THRESHOLD] minimum HMM score for exhaustive ",
"    searches.  Default is 70.0, 55.0 will give partial Ig like sequences, ",
"    10 to 30 is probably the lowest useful sensitivity.",
" [--filterThresh THRESHOLD] minimum HMM score for regular pipeline ",
"    If this is set, only sequences exceeding this threshold will be ",
"    reported.  If ommitted, no filtering is done.",
" [--filterIncompletes] flag indicating if 'incomplete' sequences",
"    should be filtered out.  Incomplete means more than 3 AA gaps or Zs in ",
"    FW1,CDR1,FW2,CDR2,FW3,FW4.  CDR3 is forgiven.",
" [--startseq] skip (startseq-1) sequences before processing, default 1",
" [--nseqs] process only nseqs, all seqs processed if omitted",
"    Notice that nseqs applies to 'pairs of read' when R2 is specified.",
" [--distrib] flag indicating job is part of distributed batch.",
"    If flag set, /tmp will be used and the results will be finally ",
"    spooled/appended to the --outdir folder with synchronization by filelock.",
"    It is harmless to use this flag on non-distributed jobs, but ",
"    it does add extra steps.",
" [--breakout H|L|NONE] break output into separate files by lineage.",
"    NONE is default, in which case all results go to --outdir/absolve.tsv",
" --outdir folder, folder should be unique to sequencing run",
"    this is where the results will be stored ",
"    as common file names will be used",
"    if folder does not exist, it will be created automatically.",
"    results are generally server .tsv files (tab delimited). ",
" ",
"Example:",
" absolve ",
"  --R1  /tmp/myNGS.R1.fastq.gz",
"  --R2  /tmp/myNGS.R2.fastq.gz",
"  --CH ./dbs/species/human/Ch.fa",
"  --VH ./dbs/species/human/VHv.fa",
"  --JH ./dbs/species/human/VHj.fa",
"  --VL ./dbs/species/human/VLv.fa",
"  --JL ./dbs/species/human/VLj.fa",
"  --startseq 100",
"  --nseqs  10",
"  --outdir /home/user1/myNGSrun77folder",
""
};
	std::vector<std::string> op;
	int i =0;
	while(1){
		op.push_back(opts[i]);
		if(op.back().empty()) break;
		i++;
	}

	fprintf(stderr, "\n\nUsage: absolve \n");
	for(unsigned int i=0;i<op.size();i++){
		fprintf(stderr,"%s\n",op[i].c_str());
	}
	exit(-1);
}

//////////////////////////////////////////////////////////////////////////////
// Read in all environment variables
void ParseEnv(char **env, std::map<std::string, std::string > & envs){

	int i = 0;
	while(env[i] != NULL){
		char *varp = strtok(env[i],"=");
		if(varp != NULL){
			std::string var = varp;
			char * valp = strtok(NULL,"");
			if(valp != NULL){
				std::string val = valp;
				envs[var] = val;
			}
		}
		i++;
	}
}
//////////////////////////////////////////////////////////////////////////////
int AbsolvePipeline::SetAligner(std::string & alignerType){
	if(alignerType.compare("SWA")==0){
		CH.m_useSSW = VH.m_useSSW = JH.m_useSSW = 1;
		CL.m_useSSW = VL.m_useSSW = JL.m_useSSW = 1;
		return 0;
	}else if(alignerType.compare("BWA") == 0){
		CH.m_useSSW = VH.m_useSSW = JH.m_useSSW = 0;
		CL.m_useSSW = VL.m_useSSW = JL.m_useSSW = 0;
		return 0;
	}else{
		fprintf(stderr, "Unrecognized aligner type:%s\n",alignerType.c_str());
		return -1;
	}
}
//////////////////////////////////////////////////////////////////////////////
int AbsolvePipeline::SetMultiHit(int multihit){
	CH.SetMultiHit(multihit);
	VH.SetMultiHit(multihit);
	JH.SetMultiHit(multihit);

	CL.SetMultiHit(multihit);
	VL.SetMultiHit(multihit);
	JL.SetMultiHit(multihit);

	return 0;
}

//////////////////////////////////////////////////////////////////////////////
int AbsolvePipeline::ReadRCFile(const char *inifilename){
	boost::property_tree::ptree pt;
	boost::property_tree::ini_parser::read_ini(inifilename, pt);

	const char *aconfig[] = {"CH", "VH", "JH", "CL", "VL", "JL" };
	for(int i=0; i<6; i++){
		std::string type = aconfig[i];
		std::string pname;
		BWAalignConfig tconfig;

		pname = type + ".A"; tconfig.A = pt.get<int>(pname);
		pname = type + ".B"; tconfig.B = pt.get<int>(pname);
		pname = type + ".OINS"; tconfig.OINS = pt.get<int>(pname);
		pname = type + ".EINS"; tconfig.EINS = pt.get<int>(pname);
		pname = type + ".ODEL"; tconfig.ODEL = pt.get<int>(pname);
		pname = type + ".EDEL"; tconfig.EDEL = pt.get<int>(pname);
		pname = type + ".CLIP3"; tconfig.CLIP3 = pt.get<int>(pname);
		pname = type + ".CLIP5"; tconfig.CLIP5 = pt.get<int>(pname);
		pname = type + ".MIN_SEED_LEN"; tconfig.MIN_SEED_LEN = pt.get<int>(pname);
		pname = type + ".MIN_SCORE"; tconfig.MIN_SCORE = pt.get<int>(pname);

		m_bwaconfigs[type] = tconfig;
	}

	std::string hmmp = m_exepath + "/hmm/";
	m_hmmfileSCFV = hmmp + pt.get<std::string>("HMM.ScFVHMMFile");
	m_hmmfileHeavy = hmmp + pt.get<std::string>("HMM.HeavyHMMFile");
	m_hmmfileLight = hmmp + pt.get<std::string>("HMM.LightHMMFile");
	m_exhaustiveHMMThresh = pt.get<float>("HMM.ExhaustiveHMMThresh");
	m_filterHMMThresh = pt.get<float>("HMM.FilterHMMThresh");

	m_tempfolder = pt.get<std::string>("TEMP.TempFolderForDistrib");

	m_NT_PAIR_MIN = pt.get<int>("pairedend.NT_PAIR_MIN");
	return 0;
}
//////////////////////////////////////////////////////////////////////////////
int AbsolvePipeline::ProcessOptions(int argc, char *argv[], char **env,
	const char *exepath){

	int retval;
	// Get hmm file from the environment
	std::map<std::string, std::string > envs;
	ParseEnv(env, envs);

	m_exepath = exepath;

	m_breakout = NONE;

	std::vector<std::string> args;
	args.clear();
	unsigned int i;
	for(i=1; i<(unsigned int)argc; i++){
		args.push_back(argv[i]);
	}

	std::string rcFile = m_exepath;
	rcFile = rcFile + "/.absolverc";
	//Read the RC file first so --args can override it later
	for(i=0;i<args.size();i++){
		if(args[i].compare("--rcfile")==0 && i<(args.size()-1)){
			rcFile = args[i+1];
			ReadRCFile(rcFile.c_str());
		}
	}
	ReadRCFile(rcFile.c_str());
	
	for(i=0;i<args.size();i++){
		if      (args[i].compare("--R1")==0 && i<(args.size()-1)){
			m_R1file = args[i+1];

		}else if(args[i].compare("--R2")==0 && i<(args.size()-1)){
			m_R2file = args[i+1];

		}else if(args[i].compare("--AA")==0 && i<(args.size()-0)){
			m_inputIsAA = 1;

		}else if(args[i].compare("--stats")==0 && i<(args.size()-0)){
			m_collectStats = 1;
		}else if(args[i].compare("--profileTiming")==0 && i<(args.size()-0)){
			m_profileTiming = 1;
		}else if(args[i].compare("--CH")==0 && i<(args.size()-1)){
			m_CHdb = args[i+1];

		}else if(args[i].compare("--VH")==0 && i<(args.size()-1)){
			m_VHdb = args[i+1];

		}else if(args[i].compare("--JH")==0 && i<(args.size()-1)){
			m_JHdb = args[i+1];

		}else if(args[i].compare("--VL")==0 && i<(args.size()-1)){
			m_VLdb = args[i+1];

		}else if(args[i].compare("--CL")==0 && i<(args.size()-1)){
			m_CLdb = args[i+1];

		}else if(args[i].compare("--JL")==0 && i<(args.size()-1)){
			m_JLdb = args[i+1];

		}else if(args[i].compare("--JL")==0 && i<(args.size()-1)){
			m_JLdb = args[i+1];

		}else if(args[i].compare("--hmmtype")==0 && i<(args.size()-1)){
			if(args[i+1].compare(HMMTypeS[Heavy]) == 0){
				m_hmmtype = Heavy;
				m_ighmmDefault = & m_ighmmHeavy;
			}else if(args[i+1].compare(HMMTypeS[Light]) == 0){
				m_hmmtype = Light;
				m_ighmmDefault = & m_ighmmLight;
			}else if(args[i+1].compare(HMMTypeS[ScFV]) == 0) {
				m_hmmtype = ScFV;
				m_ighmmDefault = & m_ighmmSCFV;
			}else{
				fprintf(stderr,
					"Unrecognized hmmtype %s\nShould be Heavy|Light|ScFV\n",
					args[i+1].c_str());
				ShowHelpExit();
			}
		}else if(args[i].compare("--useAligner")==0 && i<(args.size()-1)){
			retval = SetAligner( args[i+1]);
			if(retval != 0) ShowHelpExit();

		}else if(args[i].compare("--multihit")==0 && i<(args.size()-1)){
			SetMultiHit( atoi(args[i+1].c_str()) );

		}else if(args[i].compare("--startseq")==0 && i<(args.size()-1)){
			m_startseqArg = args[i+1];

		}else if(args[i].compare("--nseqs")==0 && i<(args.size()-1)){
			m_stopseqArg = args[i+1];
		
		}else if(args[i].compare("--distrib")==0 && i<(args.size()-0)){
			m_isDistrib = 1;
		}else if(args[i].compare("--exhaustive")==0 && i<(args.size()-0)){
			m_exhaustive = 1;
		}else if(args[i].compare("--exhaustiveThreshold")==0 && 
			i<(args.size()-1)){
			m_exhaustiveHMMThresh = atof(args[i+1].c_str());
		}else if(args[i].compare("--filterThresh")==0 && 
			i<(args.size()-1)){
			m_filterHMMThresh = atof(args[i+1].c_str());
			m_filterHMM = 1;
		}else if(args[i].compare("--filterIncompletes")==0 && 
			i<(args.size()-0)){
			m_filterIncomplete = 1;
		}else if(args[i].compare("--filterNoCall")==0 && 
				i<(args.size()-0)){
			m_filterNoCall = 1;
		}else if(args[i].compare("--breakout")==0 && i<(args.size()-1)){
			if(args[i+1].compare("H") == 0) m_breakout = H;
			if(args[i+1].compare("L") == 0) m_breakout = L;
			if(args[i+1].compare("NONE") == 0) m_breakout = NONE;

		}else if(args[i].compare("--outdir")==0 && i<(args.size()-1)){
			m_outdir = args[i+1];

		}else if(args[i][0] == '-'){
			fprintf(stderr, "unrecognized option %s\n",args[i].c_str());
			ShowHelpExit();
		}
	}

	return 0;
}


//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

std::string GetTmpFile(std::string tempfolder, 
	std::string outfile, std::string hostPID){

	boost::filesystem::path tmpfile = 
		outfile + ".pid." +  hostPID;
	tmpfile = tmpfile.filename();
	tmpfile = "/tmp" / tmpfile ;
	return(tmpfile.string());
}
//////////////////////////////////////////////////////////////////////////////
int GapCount(std::string & dom){
	return( std::count(dom.begin(), dom.end(), '-') );
}


//////////////////////////////////////////////////////////////////////////////
int HasGapsOrZs(std::string & dom, boost::regex & gaps, 
	boost::regex & zs){
	if(dom.size()==0){return(1);}//empty domains are automatically incomplete

	int hasgaps = boost::regex_match(dom,gaps);
	int hasZs = boost::regex_match(dom,zs);
	return(hasgaps || hasZs);
}
//////////////////////////////////////////////////////////////////////////////
// we define "incomplete" sequences as anything with 3 ore Zs
// or gaps in any domain.  This is used for optional filtering
int AbsolvePipeline::IncompleteSequence(Result & rslt){

	int incH = 0;
	int incL = 0;

	boost::regex gaps3(".*-{3,}.*");
	boost::regex zs3(".*Z{3,}.*");

	if(HasGapsOrZs(rslt.m_res.at("HFW1"),gaps3,zs3) ){ incH++; }
	if(HasGapsOrZs(rslt.m_res.at("HCDR1"),gaps3,zs3) ){ incH++; }
	if(HasGapsOrZs(rslt.m_res.at("HFW2") ,gaps3,zs3) ){ incH++; }
	if(HasGapsOrZs(rslt.m_res.at("HCDR2"),gaps3,zs3) ){ incH++; }
	if(HasGapsOrZs(rslt.m_res.at("HFW3") ,gaps3,zs3) ){ incH++; }
	if(HasGapsOrZs(rslt.m_res.at("HFW4") ,gaps3,zs3) ){ incH++; }

	if(HasGapsOrZs(rslt.m_res.at("LFW1"),gaps3,zs3) ){ incL++; }
	if(HasGapsOrZs(rslt.m_res.at("LCDR1"),gaps3,zs3) ){ incL++; }
	if(HasGapsOrZs(rslt.m_res.at("LFW2") ,gaps3,zs3) ){ incL++; }
	if(HasGapsOrZs(rslt.m_res.at("LCDR2"),gaps3,zs3) ){ incL++; }
	if(HasGapsOrZs(rslt.m_res.at("LFW3") ,gaps3,zs3) ){ incL++; }
	if(HasGapsOrZs(rslt.m_res.at("LFW4") ,gaps3,zs3) ){ incL++; }

	return (incH>0 && incL>0);
}
//////////////////////////////////////////////////////////////////////////////
//
void AbsolvePipeline::LoadInitHMM(
	std::string & hmmfile, IgHMM & model, HMMType type){

	if(hmmfile.empty()){
		fprintf(stderr,"No %s HMM file provided\n",HMMTypeS[type]);
		ShowHelpExit();
	}
	fprintf(stderr,"using %s HMM file: %s\n",HMMTypeS[type], hmmfile.c_str());

	model.Initialize(hmmfile.c_str(), m_NT_PAIR_MIN);
	model.m_hmmType = type;
}
//////////////////////////////////////////////////////////////////////////////

int AbsolvePipeline::InitPipeline(){

	m_stats.m_collectStats = m_collectStats;

	LoadInitHMM(m_hmmfileSCFV, m_ighmmSCFV, ScFV);
	LoadInitHMM(m_hmmfileHeavy, m_ighmmHeavy, Heavy);
	LoadInitHMM(m_hmmfileLight, m_ighmmLight, Light);

	printf("Running only HMM:%s\n", HMMTypeS[m_ighmmDefault->m_hmmType]);

	fprintf(stderr,"using CH: %s\n",m_CHdb.c_str());
	CH.Initialize(m_CHdb.c_str(), "CH",
		m_bwaconfigs.at("CH")
	);

	fprintf(stderr,"using VH: %s\n",m_VHdb.c_str());
	VH.Initialize(m_VHdb.c_str(), "VH",
		m_bwaconfigs.at("VH")
	);
	
	fprintf(stderr,"using JH: %s\n",m_JHdb.c_str());
	JH.Initialize(m_JHdb.c_str(), "JH",
		m_bwaconfigs.at("JH")
	);

	fprintf(stderr,"using CL: %s\n",m_CLdb.c_str());
	CL.Initialize(m_CLdb.c_str(), "CL",
		m_bwaconfigs.at("CL")
	);

	fprintf(stderr,"using VL: %s\n",m_VLdb.c_str());
	VL.Initialize(m_VLdb.c_str(), "VL",
		m_bwaconfigs.at("VL")
	);

	fprintf(stderr,"using JL: %s\n",m_JLdb.c_str());
	JL.Initialize(m_JLdb.c_str(), "JL",
		m_bwaconfigs.at("JL")
	);

	if(m_R1file.empty()){
		fprintf(stderr, "No R1 sequence file provided\n");	
		ShowHelpExit();
	}
	fprintf(stderr,"using R1: %s\n",m_R1file.c_str());

	if(!m_R2file.empty()){
		fprintf(stderr,"using R2: %s\n",m_R2file.c_str());
	}

	if(m_startseqArg.empty()){
		m_startseq = 1;
	}else{
		m_startseq = atoi(m_startseqArg.c_str());
	}

	if(m_stopseqArg.empty()){
		m_stopseq = -1;
	}else{
		m_stopseq = atoi(m_stopseqArg.c_str());
	}


	if(m_outdir.empty()){
		fprintf(stderr,
			"\nERROR: with --out folder must be specified\n");
		ShowHelpExit();
	}else{
		boost::filesystem::create_directories(m_outdir);
		fprintf(stderr,"using out folder:%s\n",m_outdir.c_str());

		if(m_isDistrib){
			//Touch the lock file
			std::string lfilename = m_outdir + "/absolve.lockfile";
			FILE *lfile = fopen(lfilename.c_str(), "a+");
			fclose(lfile);

			m_outdirTmp = GetTmpFile(m_tempfolder, m_outdir, m_hostPID);
			boost::filesystem::create_directories(m_outdirTmp);
			fprintf(stderr,"using tmp out dir:%s\n",m_outdirTmp.c_str());
		}
	}

	return 0;
}

////////////////////////////////////////////////////////////////////////
//

int AppendFile(const char * from, const char *to){
	
	int noheader = boost::filesystem::exists(to);
	FILE *ofile = fopen(to, "a+");
	FILE *tfile = fopen(from, "r");

	if(noheader){ //skip the header line
		char *line = NULL;
		size_t line_n = 0;
		size_t lsize = getline(&line, &line_n, tfile);
		if(lsize>0){ lsize = lsize;}
		if(line!=NULL) free(line);
	}

	int nread;
	do {
		char buf[262144];
		nread = fread(buf, 1, 262144,  tfile);
		if(nread>0){
			fwrite(buf, 1, nread, ofile);
		}
	}while(nread > 0);

	fclose(ofile);
	fclose(tfile);

	return 0;
}
////////////////////////////////////////////////////////////////////////
// distributed processed use a lockfile to synchronize the concatenation
// of results
int AbsolvePipeline::AppendResults(){

	std::string lfilename = m_outdir + "/absolve.lockfile";
	boost::interprocess::file_lock outlock(lfilename.c_str());
	outlock.lock();

	TimePoint startTime = GetCurTime();

	//Close all output file
	std::map<std::string, FILE *>::iterator it;
	for(it = m_outmap.begin(); it!=m_outmap.end(); ++it){
		FILE *ofile = (*it).second;
		if(ofile != NULL) fclose(ofile);
		(*it).second = NULL;

		boost::filesystem::path fromfile = (*it).first;

		boost::filesystem::path jfile = fromfile.filename();
		
		std::string fname = jfile.string();

		std::string outfile = m_outdir + "/" + fname;

		AppendFile(fromfile.string().c_str(), outfile.c_str());
		boost::filesystem::remove(fromfile);
	}
	m_stats.MergeStats(m_outdir);

	int runTime = GetSecondsDiffToNow(startTime);
	if(m_profileTiming){
		m_stats.AddTo(m_hostPID + ".AppendTime", 
				m_stats.m_allStats.at("statRunStats"),
				runTime);
	}
	m_stats.WriteStats(m_outdir);

	outlock.unlock();

	return 0;
}

////////////////////////////////////////////////////////////////////////

int AbsolvePipeline::Run(){

	TimePoint startTime = GetCurTime();

	if(m_inputIsAA == 1){
		RunAAPipeline();
	}else{
		if(GetSeqFileType(m_R1file.c_str()) == sqFASTQ){
			fprintf(stderr,"Running fastq ..\n");
			RunFastQPipeline();
		}else if(GetSeqFileType(m_R1file.c_str()) == sqFASTA){
			fprintf(stderr,"Running FASTA ..\n");
			RunFastaPipeline();
		}else{
			fprintf(stderr,"Unrecognized file type:%s\n",m_R1file.c_str());
		}
	}

	int runTime = GetSecondsDiffToNow(startTime);
	if(m_profileTiming){
		m_stats.AddTo(m_hostPID + ".RunTime", 
				m_stats.m_allStats.at("statRunStats"),
				runTime);
	}
	//Write out stats
	//If distrib, lock the real output file and concat our tmp results
	if(m_isDistrib){
		AppendResults();
	}else{
		m_stats.WriteStats(m_outdir);
	}

	return 0;
}
////////////////////////////////////////////////////////////////////////

int AbsolvePipeline::WriteResult(Result & result){

	std::string fname = "absolve.tsv";
	switch(m_breakout){
		case H:
			fname = result.m_res.at("Hlineage") + ".tsv";

			break;
		case L:
			fname = result.m_res.at("Llineage") + ".tsv";
			break;
		default:
			break;
	}

	// If germline has /'s replace with _
	boost::regex slash("/");
	fname = boost::regex_replace(fname, slash, "_", boost::match_default | boost::format_all);

	std::string outfilename;
	if(m_isDistrib){
		outfilename = m_outdirTmp + "/" + fname;
	}else{
		outfilename = m_outdir + "/" + fname;
	}

	FILE *ofile;
	if(m_outmap.find(outfilename) == m_outmap.end()){
		//Open the file 
		ofile = fopen(outfilename.c_str(),"w+");
		if(ofile == NULL){
			fprintf(stderr,"Failed to open output file:\n%s\n",
				outfilename.c_str());
			exit(-1);
		}

		m_outmap[outfilename] = ofile;

		Result tmp = result;
		tmp.AsHeader();
		tmp.Print(ofile, '\t');
	}else{
		ofile = m_outmap.at(outfilename);
	}

	result.Print(ofile,'\t');
	
	return 0;
}
////////////////////////////////////////////////////////////////////////
// process a single AA sequence
int AbsolvePipeline::RunAASeq(std::string & name, std::string & seq){

	Result tresult;
	tresult.Init();
	tresult.m_res["Accession"] = name; 

	m_ighmmDefault->AlignSingleAA(seq);

	tresult.m_res["DNA"] = "";
	tresult.m_res["ORF"] = m_ighmmDefault->m_mergedAA;
	tresult.m_res["Frame"] = "";

	SeqResult(seq.c_str(), seq.size(), tresult);

	return 0;

}
////////////////////////////////////////////////////////////////////////
// for exhaustive searches, mask out each hit as it is found
void MaskFirstInstance(std::string & str, std::string & toreplace,
	char maskchar){
	std::string mask(toreplace.size(),maskchar);
	str.replace(str.find(toreplace), toreplace.size(), mask);
}

////////////////////////////////////////////////////////////////////////
// for a given query sequence, find ALL heavy and light domains,
// useful for complex constructs with multiple heavy or light chains
int AbsolvePipeline::RunAASeqExhaustive(std::string & name, std::string & seq){

	std::vector<Result> hits;

	std::string maskseq = seq;

	IgHMM *hmms[2];
	hmms[0] = &m_ighmmHeavy;
	hmms[1] = &m_ighmmLight;
	
	for(int i=1;i>=0; i--){
		IgHMM & thisHmm = *hmms[i];
		do{
			Result tresult;
			tresult.Init();
			tresult.m_res["Accession"] = name;

			thisHmm.AlignSingleAA(maskseq);
			tresult.GetHMMResults(thisHmm);
			tresult.m_res["DNA"] = "";
			tresult.m_res["ORF"] = maskseq;

			std::string aa = "";
			if(i==0){
				aa = tresult.m_res.at("HAA");
			}else{
				aa = tresult.m_res.at("LAA");
			}
			//Rerun to get pure score
			thisHmm.AlignSingleAA(aa);
			tresult.m_res["HMMScore"] = 
				boost::lexical_cast<std::string>(thisHmm.m_score);

			if(thisHmm.m_score > m_exhaustiveHMMThresh) {
				hits.push_back(tresult);
				MaskFirstInstance(maskseq,aa,'X');
			}
		}while(thisHmm.m_score > m_exhaustiveHMMThresh);
	}

	for(unsigned int i=0;i<hits.size();i++){
		WriteResult(hits[i]);
	}
	
	return 0;	
}
////////////////////////////////////////////////////////////////////////
// amino acid sequences as input
int AbsolvePipeline::RunAAPipeline(){

	//Force classifiers to use AA
	CH.m_useSSW = CH.m_byAA = 1;
	VH.m_useSSW = VH.m_byAA = 1;
	JH.m_useSSW = JH.m_byAA = 1;
	CL.m_useSSW = CL.m_byAA = 1;
	VL.m_useSSW = VL.m_byAA = 1;
	JL.m_useSSW = JL.m_byAA = 1;
	
	m_AAFastaReader.OpenFile(m_R1file.c_str());	

	if(m_startseq > 1){
		m_AAFastaReader.SkipFirstNSeqs(m_startseq-1);
	}

	int totSeqsProcessed = 0;
	int continueReading = 1;
	int readCount = 0;
	
	while(continueReading){
		std::vector<std::string> readNames, readContents;
		m_AAFastaReader.ReadNSeqs(readNames, readContents, 1);
		if(readNames.size()>0){

			if(m_exhaustive){
				RunAASeqExhaustive(readNames[0], readContents[0]);
			}else{
				RunAASeq(readNames[0], readContents[0]);
			}
	

			totSeqsProcessed++;
			if( (totSeqsProcessed % 1000) == 0){
				fprintf(stderr, "%d seqs so far...\n", totSeqsProcessed);
			}

			readCount++;
			if(m_stopseq > 0){
				if(readCount >= m_stopseq){
					continueReading = 0;
					break;
				}
			}
		}

		if(readNames.size() == 0){
			continueReading = 0;
			break;
		}
	}

	m_AAFastaReader.CloseFile();

	return 0;
}
////////////////////////////////////////////////////////////////////////
//  Keep a running tally of SHM counts by kabat position.
//  gives a measure of SHM location bias.
int AbsolvePipeline::TallyMutKabat(Result & rslt, char chain, 
	std::vector<int> & mutpos){

	//Mutated positions
	int aapos;
	for(unsigned int i=0; i<mutpos.size(); i++){
		if(chain == 'H'){
			if(rslt.m_HDNAstart < 0) {continue;}
			if(rslt.m_res["HAA"].size() == 0 ){continue;}
			aapos = (mutpos[i] - rslt.m_HDNAstart)/3;

			if(aapos >= 0 && aapos < (int)m_ighmmDefault->m_Hres2kabat.size()){
				rslt.m_KabatSHMPos.push_back(
					m_ighmmDefault->m_Hres2kabat[aapos] );
			}
		}else{
			if(rslt.m_LDNAstart < 0) {continue;}
			if(rslt.m_res["LAA"].size() == 0 ){continue;}
			aapos = (mutpos[i] - rslt.m_LDNAstart)/3;

			if(aapos >= 0 && aapos < (int)m_ighmmDefault->m_Lres2kabat.size()){
				rslt.m_KabatSHMPos.push_back(
					m_ighmmDefault->m_Lres2kabat[aapos] );
			}
		}
	}

	return 0;
	
}
////////////////////////////////////////////////////////////////////////
int AbsolvePipeline::RunAllHeavyClassifiers(const char *seq, int len, 
	Result & rslt){
		
	std::vector<int> mutpos;

	int vstart = 0;
	int vend = -1;
	
	VH.RunClassifier( seq,len, rslt.m_res.at("VH"), rslt.m_res.at("VHscore"), 
		rslt.m_res.at("VHpos"), rslt.m_res.at("VHcigar"), vstart, vend);

	VH.ComputeSHM(rslt.m_res.at("VH"), seq, len,
		rslt.m_res.at("VHcigar"), rslt.m_res.at("VHshm"), mutpos);

	int jstart = 0;
	int jend = -1;

	JH.RunClassifier( seq,len, rslt.m_res.at("JH"), rslt.m_res.at("JHscore"), 
		rslt.m_res.at("JHpos"), rslt.m_res.at("JHcigar"), jstart, jend);

	JH.ComputeSHM(rslt.m_res.at("JH"), seq,len,
		rslt.m_res.at("JHcigar"), rslt.m_res.at("JHshm"), mutpos);

	CH.RunClassifier( seq,len, rslt.m_res.at("CH"), rslt.m_res.at("CHscore"),
		rslt.m_res.at("CHpos"), rslt.m_res.at("CHcigar"));

	TallyMutKabat(rslt, 'H', mutpos);

	rslt.m_res.at("Hlineage") = GetLineage(rslt.m_res.at("VH"),
		rslt.m_res.at("JH"), rslt.m_res.at("HCDR3"));

	return 0;
}
////////////////////////////////////////////////////////////////////////
int AbsolvePipeline::RunAllLightClassifiers(const char *seq, int len, 
	Result & rslt){

	std::vector<int> mutpos;

	VL.RunClassifier( seq,len, rslt.m_res.at("VL"), rslt.m_res.at("VLscore"), 
		rslt.m_res.at("VLpos"), rslt.m_res.at("VLcigar"));
	VL.ComputeSHM(rslt.m_res.at("VL"), seq,len,
		rslt.m_res.at("VLcigar"), rslt.m_res.at("VLshm"), mutpos);

	JL.RunClassifier( seq,len, rslt.m_res.at("JL"),  rslt.m_res.at("JLscore"),
		rslt.m_res.at("JLpos"), rslt.m_res.at("JLcigar"));
	JL.ComputeSHM(rslt.m_res.at("JL"), seq, len,
		rslt.m_res.at("JLcigar"),
		rslt.m_res.at("JLshm"), mutpos);

	CL.RunClassifier( seq,len, rslt.m_res.at("CL"),  rslt.m_res.at("CLscore"),
		rslt.m_res.at("CLpos"), rslt.m_res.at("CLcigar"));

	TallyMutKabat(rslt, 'L', mutpos);

	rslt.m_res.at("Llineage") = GetLineage(rslt.m_res.at("VL"), rslt.m_res.at("JL"), 
		rslt.m_res.at("LCDR3"));

	return 0;

}
////////////////////////////////////////////////////////////////////////
int AbsolvePipeline::RunAllClassifiers(const char *seq, int len, 
	Result & rslt){

	AbsolvePipeline::RunAllHeavyClassifiers(
		seq,
		len,
		rslt);
	AbsolvePipeline::RunAllLightClassifiers(
		seq,
		len, 
		rslt);
	
	return 0;
}
////////////////////////////////////////////////////////////////////////
int AbsolvePipeline::RunFastaSeqExhaustive(char *name, char *seq, int len){
	

	std::vector<Result> hits;
	char *maskseq = new char[len+1]; 
	strncpy(maskseq, seq, len);
	maskseq[len] = '\0';

	IgHMM *hmms[2];
	hmms[0] = &m_ighmmHeavy;
	hmms[1] = &m_ighmmLight;
	
	for(int i=1;i>=0; i--){
		IgHMM & thisHmm = *hmms[i];
		do{
			Result tresult;
			tresult.Init();
			tresult.m_res["Accession"] = name;

			thisHmm.AlignSingle(maskseq, len);
			tresult.GetHMMResults(thisHmm);
			tresult.m_res["DNA"] = maskseq;
			tresult.m_res["ORF"] = thisHmm.m_mergedAA;

			std::string dna = "";
			if(i==0){
				dna = tresult.m_res.at("HDNA");
				// realign to get pure score
				thisHmm.AlignSingleAA(tresult.m_res.at("HAA"));
				tresult.m_res["HMMScore"] = 
					boost::lexical_cast<std::string>(thisHmm.m_score);
			}else{
				dna = tresult.m_res.at("LDNA");
				// realign to get pure score
				thisHmm.AlignSingleAA(tresult.m_res.at("LAA"));
				tresult.m_res["HMMScore"] = 
					boost::lexical_cast<std::string>(thisHmm.m_score);
			}
			if(thisHmm.m_score > m_exhaustiveHMMThresh) {
				RunAllClassifiers(dna.c_str(), dna.size(), tresult);
				hits.push_back(tresult);

				std::string newmaskseq = maskseq;
				MaskFirstInstance(newmaskseq, dna, 'N');
				strncpy(maskseq, newmaskseq.c_str(), newmaskseq.size());
			}
		}while(thisHmm.m_score > m_exhaustiveHMMThresh);
	}

	for(unsigned int i=0;i<hits.size();i++){
		WriteResult(hits[i]);
	}
	

	delete [] maskseq;
	return 0;
}
////////////////////////////////////////////////////////////////////////
int AbsolvePipeline::RunFastaSeq(char *name, char *seq, int len){
	Result tresult;
	tresult.Init();
	tresult.m_res["Accession"] = name; 

	m_ighmmDefault->AlignSingle(seq, len);

	tresult.m_res["DNA"] = seq; 
	tresult.m_res["ORF"] = m_ighmmDefault->m_mergedAA;

	SeqResult(seq, len, tresult);

	return 0;
}
////////////////////////////////////////////////////////////////////////
int AbsolvePipeline::RunFastaPipeline(){

	m_DNAreader.OpenSeqFiles(m_R1file.c_str(), "");
	
	if(m_startseq > 1){
		m_DNAreader.SkipFirstNSeqs(m_startseq-1);
	}

	int totSeqsProcessed = 0;
	int continueReading = 1;
	int readCount = 0;
	while(continueReading){
		int nseqs = m_DNAreader.ReadNSeqs(1);
		for(int ssi=0;ssi<nseqs; ssi++){
			if(m_exhaustive){
				RunFastaSeqExhaustive(
					m_DNAreader.m_bseqs[ssi].name,
					m_DNAreader.m_bseqs[ssi].seq,
					m_DNAreader.m_bseqs[ssi].l_seq);
			}else{
				RunFastaSeq(m_DNAreader.m_bseqs[ssi].name, 
					m_DNAreader.m_bseqs[ssi].seq,
					m_DNAreader.m_bseqs[ssi].l_seq);
			}

			totSeqsProcessed++;
			if( (totSeqsProcessed % 1000) == 0){
				fprintf(stderr, "%d seqs so far...\n", totSeqsProcessed);
			}

			readCount++;
			if(m_stopseq > 0){
				if(readCount >= m_stopseq){
					continueReading = 0;
					break;
				}
			}
		}
		m_DNAreader.FreeSeqs();
		if(nseqs == 0) break; //EOF
	}
	fprintf(stderr,"Processed Seqs:%d\n",totSeqsProcessed);


	return 0;
}
////////////////////////////////////////////////////////////////////////
int AbsolvePipeline::RunFastQPipeline(){

	m_DNAreader.OpenSeqFiles(m_R1file.c_str(), m_R2file.c_str());

	Result tresult;

	if(m_startseq > 1){
		m_DNAreader.SkipFirstNSeqs(m_startseq-1);
	}
	int totSeqsProcessed = 0;
	int continueReading = 1;
	int readCount = 0;
	while(continueReading){
		tresult.Init();
		int nseqs = m_DNAreader.ReadNSeqs(1);
		readCount++;
		if(nseqs == 2){
			tresult.m_res["Accession"] = m_DNAreader.m_bseqs[0].name;

			int mergedLen=0; 
			int mergedFrame;
			m_ighmmDefault->AlignPaired(
				m_DNAreader.m_bseqs[0].seq,
				m_DNAreader.m_bseqs[0].l_seq,
				m_DNAreader.m_bseqs[1].seq,
				m_DNAreader.m_bseqs[1].l_seq,
				mergedLen,
				mergedFrame);


			if(mergedLen >  0){ 
				tresult.m_res["DNA"] = m_ighmmDefault->m_mergedPairs;
				tresult.m_res["ORF"] = m_ighmmDefault->m_mergedAA;

				SeqResult(tresult.m_res.at("DNA").c_str(),
					tresult.m_res.at("DNA").size(), tresult);
			}
			totSeqsProcessed++;
			if( (totSeqsProcessed % 1000) == 0){
				fprintf(stderr, "%d seqs so far...\n", totSeqsProcessed);
			}
		}

		m_DNAreader.FreeSeqs();
		if(nseqs == 0) break; //EOF

		if((m_stopseq > 0) && 
			(totSeqsProcessed >= m_stopseq)){
				continueReading = 0;
		}

	}
	fprintf(stderr,"Processed Seqs:%d\n",totSeqsProcessed);

	return 0;
}

////////////////////////////////////////////////////////////////////////
// generate a lineage string from metadata
std::string AbsolvePipeline::GetLineage(std::string & V, std::string & J,
		std::string & CDR3){
	std::string lin;

	// use the allele too.
	std::string v = V;
	v = v.substr(0,v.find_first_of(','));
	std::string j = J;
	j = j.substr(0,j.find_first_of(','));

	boost::regex dash("\\-");
	std::string cdr3 = boost::regex_replace(CDR3, dash, "", boost::match_default | boost::format_all);

	lin = v+"#"+j+"#"+ boost::lexical_cast<std::string>(cdr3.size());

	return lin;

}

/////////////////////////////////////////////////////////
	
int AbsolvePipeline::GetOffsetFromCigar(std::string & cigar){

	// count M and D and first S, but skip the last S.
	char numbuf[256];
	int numbufIndex = 0;
	int offset = 0;
	char lastCigar = ' ';
	for(unsigned int i=0;i<cigar.size();i++){
		char c = cigar[i];
		if(isdigit(c)){
			numbuf[numbufIndex++] = c;
		}else{
			if(c == 'M' || c == 'I' || c == 'S'){
				if(!(lastCigar != ' ' && c == 'S')){ //Don't count last S
					numbuf[numbufIndex++] = '\0';
					offset = offset + atoi(numbuf);
				}
			}
			lastCigar = c;
			numbuf[0] = '\0';
			numbufIndex = 0;
		}
	}
	return offset;
}
////////////////////////////////////////////////////////////////////////
// check if result has no call
int AllCallsEmptyStar(Result & rslt){

	for(int i =0; i<gClassifierTypesLen; i++){
		std::string cname = gClassifierTypes[i];
		std::string call = rslt.m_res.at(cname);
		if(call.length()>0){ //has a call
			if(call.compare("*") != 0){ //its not *
				return 0;
			}
		}
	}

	return 1;
}
////////////////////////////////////////////////////////////////////////
int AbsolvePipeline::SeqResult(const char *seq, int len, Result & rslt){
	
    rslt.GetHMMResults(*m_ighmmDefault);

	//check filter
	int isInc = 0;
	if(m_filterIncomplete){
		isInc = IncompleteSequence(rslt);
	}
	int badHMM = 0;
	if(m_filterHMM && m_ighmmDefault->m_score < m_filterHMMThresh){
		badHMM = 1;
	}
	if(!isInc && !badHMM){
		RunAllClassifiers(seq, len, rslt);

		int noCall = 0;
		if(m_filterNoCall){
			noCall = AllCallsEmptyStar(rslt);		
		}

		if(!noCall){
			WriteResult(rslt);
			m_stats.AddAll(rslt);
			m_stats.AddTo("GoodSeqs", 
				m_stats.m_allStats.at("statRunStats"));
		}else{
			m_stats.AddTo("FilteredSeqs", 
				m_stats.m_allStats.at("statRunStats"));
		}
	}else{
		m_stats.AddTo("FilteredSeqs", 
			m_stats.m_allStats.at("statRunStats"));
	}

	m_stats.AddTo("TotalSeqsProcessed", 
		m_stats.m_allStats.at("statRunStats"));

	return 0;
}


/////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// 
