
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//	RunStats is used to accumulate statistics for a whole run and provides
//	file io functions for reading and writing these statistics.
//	See runstats.cpp for full list of statistics collected.
//	
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#pragma once


#include <string>
#include <vector>
#include <map>

#include "fasta.h"
#include "ighmm.h"
#include "bwaaligner.h"

#include "result.h"


typedef std::map<std::string, int> IntStat;
/////////////////////////////////////////////////////////
class RunStats {
public:
	RunStats();
	~RunStats();

	int AddAll(Result & res);
	int AddTo(std::string val, IntStat & stat,
		int countToAdd = 1);
	
	int WriteStats(std::string outdir);
	int MergeStats(std::string srcdir);

public:
	std::vector<std::string> m_statNames;
	std::map<std::string, IntStat > m_allStats;
	int m_collectStats;

private:
	int AddMapToMap(IntStat & src,
		IntStat & dest);
	
	int AddToAStat(std::string statName, int toAdd);
	///////////////////////////////////////////////////////////
	int ReadMap( std::string infilename,
		IntStat & stat);
	///////////////////////////////////////////////////////////
	int WriteMap(IntStat & stat,
		std::string outfolder, std::string outfilename);
		
	int AddToMapFromFile(IntStat & dest,
		std::string srcdir, std::string fname);
	
};


/////////////////////////////////////////////////////////
 
