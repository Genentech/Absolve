
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//	A utility class for reverse complementing DNA and AA translation 
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include<stdint.h>

class AATranslator {

public:

	AATranslator();
	~AATranslator();

	int Translate(const char *dna, int dnalen, 
		int frameshift, int rc, char * aa, int & aalen);

	int Revcomp(char *nt, int len);

private:
	char GetAA(const char *codon);
	int MapThisCodon(const char *codon, char AA);
	int DumpTable();

	char * m_Codon2AA;
	int32_t * m_bit;
	char *m_chars;

	char *m_comp;

};


 
