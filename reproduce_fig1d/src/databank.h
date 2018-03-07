#ifndef DATABANK_H
#define DATABANK_H

#include "gene.h"
#include "enhancer.h"
#include "Histonmark.h"
#include "operation.h"
#include "TF.h"

class PIQcall
{
public:
	int id;
	string chr;
	int coord;
	double purity;
	PIQcall()
	{
	}
	PIQcall(int inid, string inchr, int incoord, double inpurity)
	{
		id = inid;
		chr = inchr;
		coord = incoord;
		purity = inpurity;
	}
};

class PIQcall_GW
{
public:
	vector<PIQcall > callve;
	map<string, map<int, Region_id > > chr_coord_region_map;
};

class Sample_Encode
{
public:
	int day;   // D-3, D0, D2, D7
	string type;    // WAT, BAT
	map<string, TFBinding_GW > name_TFBinding_GW_map;
	map<string, Histonmark_GW > name_Histonmark_GW_map;    //K4me1,2,3, K27ac, FAIRE
	map<string, TFBindingGPS_GW > name_TFBindingGPS_GW_map;
	map<string, TFBindingMACS_GW > name_TFBindingMACS_GW_map;

	GeneExpression_GW geneexpression;
	TranscriptExpression_GW transcriptexpression;
	lincRNAexpression_GW lincrnaexpression;
};

class Bank
{
public:
//	map<string, Transcript > name_Transcript_map;
//	map<string, map<pair<int, int >, string > > chr_pos_TranscriptName_map;
//	map<string, Gene > name_Gene_map;
//	map<string, map<pair<int, int >, string > > chr_pos_GeneName_map;
	Genome genome;
	PIQcall_GW pool_piqcall;
	
	vector<Sample_Encode > encode_Vec;
	map< string, map<int, Sample_id > > cell_time_sample_map;
};


#endif

