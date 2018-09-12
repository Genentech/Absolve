
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
#include <boost/lexical_cast.hpp>

#include "result.h"


///////////////////////////////////////////////////////////

int gClassifierTypesLen = 6;
const char* gClassifierTypes[] = {
	"CH", "VH", "JH",
	"CL", "VL", "JL"
};



///////////////////////////////////////////////////////////
Result::Result(){
	m_resNames.push_back("Accession");
	m_resNames.push_back("DNA");
	m_resNames.push_back("ORF");
	m_resNames.push_back("Frame");
	for(int i=0;i<gClassifierTypesLen;i++){
		std::string cname = gClassifierTypes[i];
		m_resNames.push_back(cname);
		m_resNames.push_back(cname + "score");
		m_resNames.push_back(cname + "pos");
		m_resNames.push_back(cname + "cigar");
		if(cname.at(0) != 'C'){ //constant regions don't have SHM
			m_resNames.push_back(cname + "shm");
		}
	}
	m_resNames.push_back("HAA");
	m_resNames.push_back("LAA");
	m_resNames.push_back("HDNA");
	m_resNames.push_back("LDNA");
	m_resNames.push_back("HFW1");
	m_resNames.push_back("HCDR1");
	m_resNames.push_back("HFW2");
	m_resNames.push_back("HCDR2");
	m_resNames.push_back("HFW3");
	m_resNames.push_back("HCDR3");
	m_resNames.push_back("HFW4");
	m_resNames.push_back("LFW1");
	m_resNames.push_back("LCDR1");
	m_resNames.push_back("LFW2");
	m_resNames.push_back("LCDR2");
	m_resNames.push_back("LFW3");
	m_resNames.push_back("LCDR3");
	m_resNames.push_back("LFW4");
	m_resNames.push_back("Llineage");
	m_resNames.push_back("Hlineage");
	m_resNames.push_back("Kabat");
	m_resNames.push_back("HMMScore");
	m_resNames.push_back("PairStrategy");
}


///////////////////////////////////////////////////////////
Result::~Result(){
}
	
///////////////////////////////////////////////////////////
void Result::Init(){
	for(unsigned int i=0; i<m_resNames.size(); i++){
		m_res[m_resNames[i]] = "";
	}
	m_KabatSHMPos.clear();
	m_KabatPos.clear();
}
/////////////////////////////////////////////////////////

void Result::AsHeader(){
	for(unsigned int i=0; i<m_resNames.size(); i++){
		m_res[m_resNames[i]] = m_resNames[i];
	}
}
/////////////////////////////////////////////////////////
// populate relevant result fields from the HMM alignment
// process
void Result::GetHMMResults(IgHMM & hmm){

	m_res["Frame"] = boost::lexical_cast<std::string>( hmm.m_bestFrame1 );
	m_res["HAA"] = hmm.m_HAA;
	
	m_res["LAA"] = hmm.m_LAA;

	m_res["HDNA"] = hmm.m_HDNA;
	m_res["LDNA"] = hmm.m_LDNA;
	m_res["HFW1"] = hmm.m_HFW1.str;
	m_res["HCDR1"] = hmm.m_HCDR1.str;
	m_res["HFW2"] = hmm.m_HFW2.str;
	m_res["HCDR2"] = hmm.m_HCDR2.str;
	m_res["HFW3"] = hmm.m_HFW3.str;
	m_res["HCDR3"] = hmm.m_HCDR3.str;
	m_res["HFW4"] = hmm.m_HFW4.str;
	m_res["LFW1"] = hmm.m_LFW1.str;
	m_res["LCDR1"] = hmm.m_LCDR1.str;
	m_res["LFW2"] = hmm.m_LFW2.str;
	m_res["LCDR2"] = hmm.m_LCDR2.str;
	m_res["LFW3"] = hmm.m_LFW3.str;
	m_res["LCDR3"] = hmm.m_LCDR3.str;
	m_res["LFW4"] = hmm.m_LFW4.str;
	m_res["PairStrategy"] = hmm.m_alignMethod;
	m_res["Kabat"] = hmm.m_kabatstr;
	m_res["HMMScore"] = boost::lexical_cast<std::string>( hmm.m_score );
	m_HDNAstart = hmm.m_HDNAstart;
	m_LDNAstart = hmm.m_LDNAstart;

	for(unsigned int i=0; i<hmm.m_Hres2kabat.size();i++){
		m_KabatPos.push_back(hmm.m_Hres2kabat[i]);
	}
	for(unsigned int i=0; i<hmm.m_Lres2kabat.size();i++){
		m_KabatPos.push_back(hmm.m_Lres2kabat[i]);
	}

	// Get domain start stop details
	m_HFW1 = hmm.m_HFW1;
    m_HCDR1 = hmm.m_HCDR1;
    m_HFW2 = hmm.m_HFW2;
    m_HCDR2 = hmm.m_HCDR2;
    m_HFW3 = hmm.m_HFW3;
    m_HCDR3 = hmm.m_HCDR3;
    m_HFW4 = hmm.m_HFW4;

    m_LFW1 = hmm.m_LFW1;
    m_LCDR1 = hmm.m_LCDR1;
    m_LFW2 = hmm.m_LFW2;
    m_LCDR2 = hmm.m_LCDR2;
    m_LFW3 = hmm.m_LFW3;
    m_LCDR3 = hmm.m_LCDR3;
    m_LFW4 = hmm.m_LFW4;


}
/////////////////////////////////////////////////////////

void Result::Print(FILE *out, char delim){
	for(unsigned int i=0; i<m_resNames.size(); i++){
		if(i < (m_resNames.size()-1)){
			fprintf(out, "%s%c", m_res.at(m_resNames[i]).c_str(), delim);
		}else{
			fprintf(out, "%s\n", m_res.at(m_resNames[i]).c_str());
		}
	}
}

void Result::Dump(){
	for(unsigned int i=0; i<m_resNames.size(); i++){
		printf("%s\t%s\n",m_resNames[i].c_str(),
			m_res.at(m_resNames[i]).c_str());
	}
}



/////////////////////////////////////////////////////////
///////////////////////////////////////////////////////// 
