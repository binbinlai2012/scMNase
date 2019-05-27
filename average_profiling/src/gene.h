#ifndef GENE_H
#define GENE_H

#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <fstream>
#include "operation.h"
using namespace std;



class Transcript
{
public:
	int start;
	int end;
	string chr;
	char strand;
	string geneName;   // Symbol
	string name;       // Refseq id or Ensmbl name
	int cdsStart;
	int cdsEnd;
	int exonCount;
	vector<int > exonStarts;
	vector<int > exonEnds;
	map< Sample_id, double > expression_map;
	double epus_idx;
	bool active;
	
	Transcript()
	{

	}

	Transcript(string inname, string inchr, char instrand, int instart, int inend, int incdsStart, int incdsEnd, int inexonCount, vector<int> &inexonStarts,
		vector<int> &inexonEnds, string ingeneName )
	{
		start = instart;
		end = inend;
		chr = inchr;
		strand = instrand;
		name = inname;
		geneName = ingeneName;
		cdsStart = incdsStart;
		cdsEnd = incdsEnd;
		exonCount = inexonCount;
		exonStarts = inexonStarts;
		exonEnds = inexonEnds;
		
	}
	
	Transcript(string inname, string inchr, char instrand, int instart, int inend, int incdsStart, int incdsEnd, int inexonCount, vector<int> &inexonStarts,
		vector<int> &inexonEnds )
	{
		start = instart;
		end = inend;
		chr = inchr;
		strand = instrand;
		name = inname;
		cdsStart = incdsStart;
		cdsEnd = incdsEnd;
		exonCount = inexonCount;
		exonStarts = inexonStarts;
		exonEnds = inexonEnds;
		
	}

	Transcript(string inname, string inchr, char instrand, int instart, int inend, string ingeneName )
	{
		start = instart;
		end = inend;
		chr = inchr;
		strand = instrand;
		name = inname;
		geneName = ingeneName;
		
	} 
	
	void addTranscript(string inname, string inchr, char instrand, int instart, int inend, int incdsStart, int incdsEnd, int inexonCount, vector<int> &inexonStarts,
		vector<int> &inexonEnds, string ingeneName);
	void addTranscript(string inname, string inchr, char instrand, int instart, int inend, string ingeneName);

	int getTSS();
	int getTES();
	pair<int, int > getPromoter(int extention );
	pair<int, int > getPromoter();
	pair<int, int > getPromoter(int up, int down );
	pair<int, int > getgenebody(int up, int down );
	string getgenetype();
};



class Gene
{
public:
	int start;
	int end;
	string chr;
	char strand;
	set<string> TranscriptName;
	map< Sample_id, double > expression_map;
	string type;
	double epus_idx;
	bool active;
	Gene()
	{
	}

	Gene(int instart, int inend, string inchr)
	{
		start = instart;
		end = inend;
		chr = inchr;
	}

	void addTranscript(string name);
	
};

class LincRNAs
{
public:
	int start;
	int end;
	string chr;
	int win;
	int id;
	
	LincRNAs()
	{
	}
	
	LincRNAs(int inid, string inchr, int instart, int inend, int inwin )
	{
		id = inid;
		chr = inchr;
		start = instart;
		end = inend;
		win = inwin;
	}
	
};

class Activepromoter
{
public:
	string chr;
	int TSS;
	double epus_idx;
	string name;
	map< Sample_id, double > k4me3_map;
};

class Genome
{
public:
	map<string, Transcript > name_Transcript_map;
	map<string, vector<pair<int, int > > > intergenicregion;
	
	map<string, map<pair<int, int >, string > > chr_pos_TranscriptName_map;
	map<string, map<pair<int, int >, vector<string> > > chr_pos_TranscriptNameS_map;   // same pos has multiple distinct trancripts
	map<string, Gene > name_Gene_map;
	map<string, map<pair<int, int >, string > > chr_pos_GeneName_map;
	map<string, map<pair<int, int >, vector<string> > > chr_pos_GeneNameS_map;   
	map<int, LincRNAs > id_LincRNAs_map;
	map<string, map<pair<int, int >, int > > chr_pos_LincRNAsId_map;
	
	vector< Activepromoter > actpro_Vec;
	map<string, map< int, Region_id > > chr_pos_actpro_map;
	
	map<string, string > genomeseq;
	
	map<string, map<pair<int, int>, string > > get_chr_promoter_TranscriptName_map(int up, int down);
	map<string, set<pair<int, int > > > get_chr_promoter( int up, int down );
	map<string, set<pair<int, int > > > get_genebody( int up, int down );
	map<string, set<int > > get_chr_TSS();
	vector<int > get_TSS_from_Gene( string gene);
	void transcripttogene();
	void addgenomeseq(string infile);
	string getsubseq(string chr, int pos, int len );
	
	
	void getintergenicregion( );
	void readtranscriptfromucsc(string &infile );
	
};


class GeneExpression_GW
{
public:
	map<string, double > name_exp_map;

};

class TranscriptExpression_GW
{
public:
	map<string, double > name_exp_map;
};



class lincRNAexpression_GW
{
public:
	map<int, double> id_exp_map;
};



pair<string, int> getnearestgene( string chr, int start, int end, Genome &gn);

#endif

