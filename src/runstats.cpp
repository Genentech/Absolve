
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <sys/types.h>
#include <unistd.h>

#include <vector>
#include <map>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include <iostream>
#include <fstream>
#include <boost/iostreams/filter/gzip.hpp>
#include <boost/iostreams/filtering_streambuf.hpp>

#include <boost/interprocess/sync/file_lock.hpp>
#include <boost/lexical_cast.hpp>


#include "runstats.h"


///////////////////////////////////////////////////////////
RunStats::RunStats(){
	m_collectStats = 0;

	m_statNames.push_back("statRunStats");
	m_statNames.push_back("statPairing");
	m_statNames.push_back("statVH");
	m_statNames.push_back("statVHshm");
	m_statNames.push_back("statJH");
	m_statNames.push_back("statJHshm");
	m_statNames.push_back("statVL");
	m_statNames.push_back("statVLshm");
	m_statNames.push_back("statJL");
	m_statNames.push_back("statJLshm");

	m_statNames.push_back("statKabatPos"); //Wt kabat pos
	m_statNames.push_back("statKabatSHMPos"); //mut kabat pos

	m_statNames.push_back("statHlineage");
	m_statNames.push_back("statLlineage");
	m_statNames.push_back("statHCDR3");
	m_statNames.push_back("statLCDR3");
	m_statNames.push_back("statLCDR3");


	m_statNames.push_back("statHMMScore");


	m_statNames.push_back("statHCDR3Lengths");
	m_statNames.push_back("statLCDR3Lengths");

	IntStat anEmptyMap;
	for(unsigned int i=0; i<m_statNames.size(); i++){
		m_allStats[m_statNames[i]] = anEmptyMap;
	}
}


RunStats::~RunStats(){
}
///////////////////////////////////////////////////////////
//
int RunStats::AddTo(std::string val, IntStat & stat,
	int countToAdd){ //countToAdd defaults to 1
	
	IntStat::const_iterator it = stat.find(val);
	if ( it == stat.end() ) {
		stat[val] = countToAdd;
	}else{
		stat[val] += countToAdd;
	}

	return 0;
}

///////////////////////////////////////////////////////////
int RunStats::AddToMapFromFile(IntStat & dest,
	std::string srcdir, std::string fname){
	
	std::string infilename = srcdir + "/" + fname;
	ReadMap( infilename, dest);
	
	return 0;
}	
int RunStats::AddMapToMap(IntStat & src,
	IntStat & dest){
	return 0;
}
///////////////////////////////////////////////////////////
int RunStats::ReadMap( std::string infilename,
	IntStat & stat){
	std::ifstream in;
	std::string line;
	std::string item;
	std::string counts;

	in.open(infilename.c_str());
	while(std::getline( in, line).good()){
		if(line.empty()){continue;}
		std::size_t tabpos = line.find("\t");
		if( tabpos != std::string::npos ){
			item.clear(); counts.clear();
			item = line.substr(0,tabpos);
			counts = line.substr(tabpos+1, line.size());
			std::size_t count = boost::lexical_cast<size_t>(counts); 
	
			IntStat::const_iterator it = 
				stat.find(item);
			if(it == stat.end() ){
				stat[item] = count;
			}else{
				stat[item] += count;
			}
		}
	}
	in.close();
	return 0;
}

///////////////////////////////////////////////////////////
int RunStats::WriteMap(IntStat & stat,
	std::string outfolder, std::string outfilename){
	
	std::string outf = outfolder + "/" + outfilename;
	std::ofstream out;
	out.open(outf.c_str());
	for (IntStat::iterator it = stat.begin();
		it != stat.end(); ++it){
		out << it->first << "\t" << it->second << "\n";

	}
	out.close();

	return 0;
}

int RunStats::AddToAStat(std::string statName, int toAdd){
	//Is this a valid stat name;
	std::map<std::string, IntStat>::const_iterator it = 
		m_allStats.find(statName);
	if ( it == m_allStats.end() ) {
		return -1;
	}

	return 0;
}
///////////////////////////////////////////////////////////
int RunStats::AddAll(Result & tresult){
	

	if(!m_collectStats) return 0;

	AddTo(tresult.m_res.at("PairStrategy"), m_allStats.at("statPairing"));
	AddTo(tresult.m_res.at("VH"), m_allStats.at("statVH"));
	AddTo(tresult.m_res.at("VHshm"), m_allStats.at("statVHshm"));
	AddTo(tresult.m_res.at("JH"), m_allStats.at("statJH"));
	AddTo(tresult.m_res.at("JHshm"), m_allStats.at("statJHshm"));
	AddTo(tresult.m_res.at("VL"), m_allStats.at("statVL"));
	AddTo(tresult.m_res.at("VLshm"), m_allStats.at("statVLshm"));
	AddTo(tresult.m_res.at("JL"), m_allStats.at("statJL"));
	AddTo(tresult.m_res.at("JLshm"), m_allStats.at("statJLshm"));

	AddTo(tresult.m_res.at("Hlineage"), m_allStats.at("statHlineage"));
	AddTo(tresult.m_res.at("Llineage"), m_allStats.at("statLlineage"));

	AddTo(tresult.m_res.at("HCDR3"), m_allStats.at("statHCDR3"));
	AddTo(tresult.m_res.at("LCDR3"), m_allStats.at("statLCDR3"));

	for(unsigned int i=0; i < tresult.m_KabatSHMPos.size();i++){
		AddTo(tresult.m_KabatSHMPos[i], m_allStats.at("statKabatSHMPos"));
	}

	for(unsigned int i=0; i < tresult.m_KabatPos.size();i++){
		AddTo(tresult.m_KabatPos[i], m_allStats.at("statKabatPos"));
	}

	double hmmScore = boost::lexical_cast<double>(tresult.m_res.at("HMMScore"));
	AddTo(boost::lexical_cast<std::string>((int) hmmScore),
		m_allStats.at("statHMMScore"));

	AddTo(boost::lexical_cast<std::string>(tresult.m_res.at("HCDR3").size()),
		m_allStats.at("statHCDR3Lengths"));
	AddTo(boost::lexical_cast<std::string>(tresult.m_res.at("LCDR3").size()),
		m_allStats.at("statLCDR3Lengths"));


	return 0;
}
///////////////////////////////////////////////////////////
int RunStats::WriteStats(std::string outdir){
	for(unsigned int i=0; i<m_statNames.size(); i++){
		if( (m_statNames[i].compare("statRunStats")==0) || 
			(m_collectStats != 0)) {
			WriteMap(m_allStats.at(m_statNames[i]), outdir, 
				m_statNames[i] + ".tsv");
		}
	}

	return 0;
}

///////////////////////////////////////////////////////////
int RunStats::MergeStats(std::string srcdir){

	for(unsigned int i=0; i<m_statNames.size(); i++){
		if( (m_statNames[i].compare("statRunStats")==0) || 
			(m_collectStats != 0)) {

			AddToMapFromFile(m_allStats.at(m_statNames[i]), 
				srcdir, m_statNames[i] + ".tsv");
		}
	}

	return 0;
}
///////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////// 
