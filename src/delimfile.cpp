
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include<boost/tokenizer.hpp>
#include "delimfile.h"


////////////////////////////////////////////////////////////////////////////
// separator must be prepended to this.
const char *fieldRegxSuffix = "(?=(?:[^\"]*\"[^\"]*\")*(?![^\"]*\"))";

////////////////////////////////////////////////////////////////////////////

int DelimFileReader::OpenFile(const char *fname, 
		const char *fieldsep, int hasHdr){

	m_fname = fname;
	m_hasHdr = hasHdr;
	m_sep = fieldsep;

	std::string regx = m_sep;
	regx += fieldRegxSuffix;
	m_fieldRegx = regx;

	m_infile.open(m_fname.c_str());
	std::string line;
	if(hasHdr){	
		std::getline(m_infile, line);
		boost::sregex_token_iterator ti(line.begin(), 
			line.end(), m_fieldRegx, -1);
		boost::sregex_token_iterator end2;

		while (ti != end2) {
			std::string token = ti->str();
			++ti;
			m_cols.push_back(token);
		}

		if (line.compare(m_sep) == 0) {
			// last character was a separator
			m_cols.push_back("");
		}
		for(unsigned int i=0;i<m_cols.size();i++){
			m_hdrMap[m_cols[i]] = i;
		}
	}
	return 0;
}
////////////////////////////////////////////////////////////////////////////
int DelimFileReader::SkipRows(int ntoskip){
	int skipped = 0;
	std::string line;
	while(m_infile.good() && !m_infile.eof()){
		std::getline(m_infile, line);
		skipped++;
		if(skipped >= ntoskip) return(skipped);
	}
	return skipped;
}

////////////////////////////////////////////////////////////////////////////
DelimFileRow & DelimFileReader::GetColNames(){
	return m_cols;
}

////////////////////////////////////////////////////////////////////////////
int DelimFileReader::ReadRow(DelimFileRow & row){
	std::string line;
	std::getline(m_infile, line);
   	
	boost::tokenizer<boost::escaped_list_separator<char> > tk(line, 
		boost::escaped_list_separator<char>('\\', m_sep.at(0), '\"'));
	for (boost::tokenizer<boost::escaped_list_separator<char> >::iterator i(tk.begin());
			i!=tk.end();++i) {
		row.push_back(*i);
	}

	return 0;

}

////////////////////////////////////////////////////////////////////////////
int DelimFileReader::ReadRowRegex(DelimFileRow & row){
	std::string line;
	std::getline(m_infile, line);
	boost::sregex_token_iterator ti(line.begin(), 
		line.end(), m_fieldRegx, -1);
	boost::sregex_token_iterator end2;

	while (ti != end2) {
		std::string token = ti->str();
		++ti;
		row.push_back(token);
	}

	if (line.compare(m_sep) == 0) {
		// last character was a separator
		row.push_back("");
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////
int DelimFileReader::ReadRow(DelimFileRow & row,
		std::vector<int> & colsOfInterest){
	DelimFileRow trow;
	ReadRow(trow);
	for(unsigned int i=0;i<colsOfInterest.size();i++){
		if(colsOfInterest[i] != -1){
			if(colsOfInterest[i] < (int)trow.size()){
				row.push_back(trow[colsOfInterest[i]]);
			}
		}else{
			row.push_back("");
		}
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////
int DelimFileReader::ReadRows(std::vector<DelimFileRow> & rows,
		int nRowsToRead){
	
	for(int i=0;i<nRowsToRead;i++){
		if(!m_infile.good() || m_infile.eof()){break;}
		DelimFileRow arow;
		ReadRow(arow);
		if(arow.size()==0 && m_infile.eof()){break;} //empty end line
		rows.push_back(arow);
	}

	return 0;
}
////////////////////////////////////////////////////////////////////////////
int DelimFileReader::ReadRows(std::vector<DelimFileRow> & rows,
		int nRowsToRead, std::vector<int> & colsOfInterest){

	for(int i=0;i<nRowsToRead;i++){
		if(!m_infile.good() || m_infile.eof()){break;}
		DelimFileRow arow;
		ReadRow(arow);
		if(arow.size()==0 && m_infile.eof()){break;} //empty end line
		DelimFileRow trow;
		for(unsigned int i=0;i<colsOfInterest.size();i++){
			if(colsOfInterest[i] != -1){
				if(colsOfInterest[i] < (int)arow.size()){
					trow.push_back(arow[colsOfInterest[i]]);
				}
			}else{
				trow.push_back("");
			}
		}
		rows.push_back(trow);
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////
int DelimFileReader::AtEOF(){
	return(m_infile.eof() || !m_infile.good());
}
////////////////////////////////////////////////////////////////////////////
int DelimFileReader::GetColIndices(std::vector<std::string> & colsOfInterest, 
		std::vector<int> & colIndices){

	for(unsigned int i=0; i<colsOfInterest.size(); i++){
		std::map<std::string,int>::iterator it = 
			m_hdrMap.find(colsOfInterest[i]);
		if(it == m_hdrMap.end()){
			colIndices.push_back(-1);
		}else{
			colIndices.push_back(it->second);
		}
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////////

int DelimFileReader::PrintRow(DelimFileRow & arow, FILE *of){
	for(unsigned int i=0;i<arow.size(); i++){
		if(i>0){fprintf(of, "%s", m_sep.c_str()); }
		fprintf(of, "%s", arow[i].c_str()); 
	}	
	fprintf(of, "\n");

	return 0;
}
 
////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////
