
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//	DelimFileReader provides functions for reading delimited files
//	with column headers.
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include <boost/regex.hpp>
#include<iostream>
#include<fstream>
#include<vector>
#include<string>
#include<map>


typedef std::vector<std::string> DelimFileRow;

class DelimFileReader {
public:
	int OpenFile(const char *fname, const char *fieldsep, int hasHdr);
	
	//slow but robust
	int ReadRowRegex(DelimFileRow & row);
	//fast but less robust
	int ReadRow(DelimFileRow & row);

	int ReadRow(DelimFileRow & row, std::vector<int> & colsOfInterest);
	int ReadRows(std::vector<DelimFileRow> & rows, int nRowsToRead);
	int ReadRows(std::vector<DelimFileRow> & rows, int nRowsToRead, 
		std::vector<int> & colsOfInterest);

	int SkipRows(int ntoskip);

	DelimFileRow & GetColNames();

	int AtEOF();

	int GetColIndices(std::vector<std::string> & colsOfInterest, 
		std::vector<int> & colIndices);

	int PrintRow(DelimFileRow & arow, FILE *of = stdout);

private:
	std::string m_fname;
	std::ifstream m_infile;
	std::string m_sep;
	DelimFileRow m_cols;
	std::map<std::string, int> m_hdrMap;
	int m_hasHdr;

	boost::regex m_fieldRegx;
}; 
