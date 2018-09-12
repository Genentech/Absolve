
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
//
//  License: GPL3, see included COPYING file
//
//
//
//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

#include<stdio.h>
#include <stdlib.h>
#include <string.h>
#include <utility>

#include"translate.h"



//  A 0000 0
//  C 0001 1
//  G 0010 2
//  T 0011 3
//  U 0100 4
//  W 0101 5
//  S 0110 6
//  M 0111 7
//  K 1000 8
//  R 1001 9
//  Y 1010 10
//  B 1011 11
//  D 1100 12
//  H 1101 13
//  V 1110 14
//  N 1111 15


//////////////////////////////////////////////////////////////////////////


char AATranslator::GetAA(const char *codon){
	int index = ( (m_bit[(int)codon[0]] << 8) & ((int32_t) 0xF00) ) |
				( (m_bit[(int)codon[1]] << 4) & ((int32_t) 0x0F0) ) |
				( (m_bit[(int)codon[2]]) & ((int32_t) 0x00F) ) ;
	return(m_Codon2AA[index] );
}

//////////////////////////////////////////////////////////////////////////
int AATranslator::MapThisCodon(const char *codon, char AA){

	int index = ( (m_bit[(int)codon[0]] << 8) & ((int32_t) 0xF00) ) |
				( (m_bit[(int)codon[1]] << 4) & ((int32_t) 0x0F0) ) |
				( (m_bit[(int)codon[2]]) & ((int32_t) 0x00F) ) ;
	m_Codon2AA[index] = AA;
	return 0;
}


//////////////////////////////////////////////////////////////////////////

AATranslator::~AATranslator(){
	delete [] m_bit;
	delete [] m_Codon2AA;
	delete [] m_chars;
	delete [] m_comp;
}


//////////////////////////////////////////////////////////////////////////


AATranslator::AATranslator(){

	// map nucleotide to sub bitpattern
	m_bit = new int32_t[256];
	for(int i=0;i<256; i++){ m_bit[i] = 15; } //001, def map to 'N'
	m_bit['A'] =0;// 0000 
	m_bit['C'] =1;// 0001 
	m_bit['G'] =2;// 0010 
	m_bit['T'] =3;// 0011 
	m_bit['U'] =4;// 0100 
	m_bit['W'] =5;// 0101 
	m_bit['S'] =6;// 0110 
	m_bit['M'] =7;// 0111 
	m_bit['K'] =8;// 1000 
	m_bit['R'] =9;// 1001 
	m_bit['Y'] =10;// 1010 
	m_bit['B'] =11;// 1011 
	m_bit['D'] =12;// 1100 
	m_bit['H'] =13;// 1101 
	m_bit['V'] =14;// 1110 
	m_bit['N'] =15;// 1111 

	m_chars = new char[16];
	m_chars[0] ='A';// 0000 
	m_chars[1] ='C';// 0001 
	m_chars[2] ='G';// 0010 
	m_chars[3] ='T';// 0011 
	m_chars[4] ='U';// 0100 
	m_chars[5] ='W';// 0101 
	m_chars[6] ='S';// 0110 
	m_chars[7] ='M';// 0111 
	m_chars[8] ='K';// 1000 
	m_chars[9] ='R';// 1001 
	m_chars[10] ='Y';// 1010 
	m_chars[11] ='B';// 1011 
	m_chars[12] ='D';// 1100 
	m_chars[13] ='H';// 1101 
	m_chars[14] ='V';// 1110 
	m_chars[15] ='N';// 1111 


	m_comp = new char[256];
	for(int i=0;i<256; i++) {m_comp[i] = 'N';}
	m_comp['A'] = 'T';
	m_comp['G'] = 'C';
	m_comp['C'] = 'G';
	m_comp['T'] = 'A';
	m_comp['U'] = 'A';
	m_comp['Y'] = 'R';
	m_comp['R'] = 'Y';
	m_comp['W'] = 'W';
	m_comp['S'] = 'S';
	m_comp['K'] = 'M';
	m_comp['M'] = 'K';
	m_comp['D'] = 'H';
	m_comp['V'] = 'B';
	m_comp['H'] = 'D';
	m_comp['B'] = 'V';
	m_comp['N'] = 'N';


	// map quad bit pattern to AA
	m_Codon2AA = new char[4096];
	for(int i=0;i<4096; i++){ m_Codon2AA[i] = 'Z'; }

	MapThisCodon("GCU",'A');
	MapThisCodon("GCC",'A');
	MapThisCodon("GCA",'A');
	MapThisCodon("GCG",'A');
	MapThisCodon("GCN",'A');
		MapThisCodon("GCT",'A');


	MapThisCodon("CGU",'R');
	MapThisCodon("CGC",'R');
	MapThisCodon("CGA",'R');
	MapThisCodon("CGG",'R');
	MapThisCodon("AGA",'R');
	MapThisCodon("AGG",'R');
	MapThisCodon("CGN",'R');
	MapThisCodon("MGR",'R');
		MapThisCodon("CGT",'R');

	MapThisCodon("AAU",'N');
	MapThisCodon("AAC",'N');
	MapThisCodon("AAY",'N');
		MapThisCodon("AAT",'N');

	MapThisCodon("GAU",'D');
	MapThisCodon("GAC",'D');
	MapThisCodon("GAY",'D');
		MapThisCodon("GAT",'D');

	MapThisCodon("UGU",'C');
	MapThisCodon("UGC",'C');
	MapThisCodon("UGY",'C');
		MapThisCodon("TGT",'C');
		MapThisCodon("TGC",'C');
		MapThisCodon("TGY",'C');
	
	MapThisCodon("CAA",'Q');
	MapThisCodon("CAG",'Q');
	MapThisCodon("CAR",'Q');
	
	MapThisCodon("GAA",'E');
	MapThisCodon("GAG",'E');
	MapThisCodon("GAR",'E');

	MapThisCodon("GGU",'G');
	MapThisCodon("GGC",'G');
	MapThisCodon("GGA",'G');
	MapThisCodon("GGG",'G');
	MapThisCodon("GGN",'G');
		MapThisCodon("GGT",'G');

	MapThisCodon("CAU",'H');
	MapThisCodon("CAC",'H');
	MapThisCodon("CAY",'H');
		MapThisCodon("CAT",'H');

	MapThisCodon("AUU",'I');
	MapThisCodon("AUC",'I');
	MapThisCodon("AUA",'I');
	MapThisCodon("AUH",'I');
		MapThisCodon("ATT",'I');
		MapThisCodon("ATC",'I');
		MapThisCodon("ATA",'I');
		MapThisCodon("ATH",'I');


	MapThisCodon("UUA",'L');
	MapThisCodon("UUG",'L');
	MapThisCodon("CUU",'L');
	MapThisCodon("CUC",'L');
	MapThisCodon("CUA",'L');
	MapThisCodon("CUG",'L');
	MapThisCodon("YUR",'L');
	MapThisCodon("CUN",'L');
		MapThisCodon("TTA",'L');
		MapThisCodon("TTG",'L');
		MapThisCodon("CTT",'L');
		MapThisCodon("CTC",'L');
		MapThisCodon("CTA",'L');
		MapThisCodon("CTG",'L');
		MapThisCodon("YTR",'L');
		MapThisCodon("CTN",'L');



	MapThisCodon("AAA",'K');
	MapThisCodon("AAG",'K');
	MapThisCodon("AAR",'K');

	MapThisCodon("AUG",'M');
		MapThisCodon("ATG",'M');

	MapThisCodon("UUU",'F');
	MapThisCodon("UUC",'F');
	MapThisCodon("UUY",'F');
		MapThisCodon("TTT",'F');
		MapThisCodon("TTC",'F');
		MapThisCodon("TTY",'F');



	MapThisCodon("CCU",'P');
		MapThisCodon("CCT",'P');
	MapThisCodon("CCC",'P');
	MapThisCodon("CCA",'P');
	MapThisCodon("CCG",'P');
	MapThisCodon("CCN",'P');

	MapThisCodon("UCU",'S');
	MapThisCodon("UCC",'S');
	MapThisCodon("UCA",'S');
	MapThisCodon("UCG",'S');
	MapThisCodon("AGU",'S');
	MapThisCodon("AGC",'S');
	MapThisCodon("UCN",'S');
	MapThisCodon("AGY",'S');
		MapThisCodon("TCT",'S');
		MapThisCodon("TCC",'S');
		MapThisCodon("TCA",'S');
		MapThisCodon("TCG",'S');
		MapThisCodon("AGT",'S');
		MapThisCodon("TCN",'S');


	MapThisCodon("ACU",'T');
		MapThisCodon("ACT",'T');
	MapThisCodon("ACC",'T');
	MapThisCodon("ACA",'T');
	MapThisCodon("ACG",'T');
	MapThisCodon("ACN",'T');

	MapThisCodon("UGG",'W');
		MapThisCodon("TGG",'W');

	MapThisCodon("UAU",'Y');
	MapThisCodon("UAC",'Y');
	MapThisCodon("UAY",'Y');
		MapThisCodon("TAT",'Y');
		MapThisCodon("TAC",'Y');
		MapThisCodon("TAY",'Y');


	MapThisCodon("GUU",'V');
	MapThisCodon("GUC",'V');
	MapThisCodon("GUA",'V');
	MapThisCodon("GUG",'V');
	MapThisCodon("GUN",'V');
		MapThisCodon("GTT",'V');
		MapThisCodon("GTC",'V');
		MapThisCodon("GTA",'V');
		MapThisCodon("GTG",'V');
		MapThisCodon("GTN",'V');

//////////////////////////////////////////////////////
	MapThisCodon("UAA",'X');
	MapThisCodon("UGA",'X');
	MapThisCodon("UAG",'X');
	MapThisCodon("UAR",'X');
	MapThisCodon("URA",'X');

		MapThisCodon("TAA",'X');
		MapThisCodon("TGA",'X');
		MapThisCodon("TAG",'X');
		MapThisCodon("TAR",'X');
		MapThisCodon("TRA",'X');

}


//////////////////////////////////////////////////////////////////////////
int AATranslator::DumpTable(){

	if(1){
		for(int c1 = 0; c1<16; c1++){
		for(int c2 = 0; c2<16; c2++){
		for(int c3 = 0; c3<16; c3++){
			char tcodon[4];
			tcodon[0] = m_chars[c1];
			tcodon[1] = m_chars[c2];
			tcodon[2] = m_chars[c3];
			tcodon[3] = '\0';
			printf("%s : %c\n",tcodon, GetAA(tcodon));
		}
		}
		}
	}

	if(1){
		char *aa = NULL;
		int aalen = 0;

		char *seq = (char*)"GAGGTGCAGCTGGTGGAGTCTGGGGGCGGCTTGGTCCAGCCAGGGGGG";
		int slen = strlen(seq);

		for(int rc=0;rc<2; rc++){
			for(int frameshift=0;frameshift<3; frameshift++){
				printf("rc: %d\tfs: %d\n",rc, frameshift);
				Translate(seq, slen, frameshift, rc, aa, aalen);
				for(int aai=0;aai<aalen; aai++) printf("%c",aa[aai]);
				printf("\n\n");
				free(aa);
				aa = NULL;
				
			}
		}

	}
	return 0;
}



//////////////////////////////////////////////////////////////////////////
// fast ambig translation
int AATranslator::Translate(const char *dna, int dnalen, 
	int frameshift, int rc, char * aa, int & aalen){

	aalen = (dnalen-frameshift)/3;
	if(aa == NULL){
		fprintf(stderr,"NULL aa in AATranslator::Translate\n");
		exit(0);
	}

	int nti;
	int index;
	for(int ai = 0; ai < aalen; ai++){
		if(rc){
			nti = dnalen - (ai*3 + frameshift) - 1;
			index = 
				((m_bit[(int)m_comp[(int)dna[nti]]] << 8) & ((int32_t) 0xF00))|
				((m_bit[(int)m_comp[(int)dna[nti-1]]]<<4) & ((int32_t) 0x0F0))|
				( (m_bit[(int)m_comp[(int)dna[nti-2]]]) & ((int32_t) 0x00F) );

			aa[ai] = m_Codon2AA[index];
	
		}else{
			nti = ai*3 + frameshift;
			index = ( (m_bit[(int)dna[nti]] << 8) & ((int32_t) 0xF00) ) |
				( (m_bit[(int)dna[nti+1]] << 4) & ((int32_t) 0x0F0) ) |
				( (m_bit[(int)dna[nti+2]]) & ((int32_t) 0x00F) );

			aa[ai] = m_Codon2AA[index];
	
		}
	}

	aa[aalen] = '\0';

	return 0;
}

//////////////////////////////////////////////////////////////////////////
int AATranslator::Revcomp(char *nt, int len){
	//reverse
	// 0123   4
	// AGTC
	// TCAT
	// AxxT    i=0 4  4/2  i<2
	// ACAT    i=1 4  4/2  i<2
	// AxT     3  3/2  i<1

	//Complement
	for(int i=0;i<len;i++){
		nt[i] = m_comp[(int)nt[i]];
	}
	//Reverse
	for(int i=0; i<(len/2);i++){
		std::swap(nt[i],nt[len-i-1]);
	}
	

	return 0;
}

//////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////// 
