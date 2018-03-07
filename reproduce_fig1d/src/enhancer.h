#ifndef ENHANCER_H
#define ENHANCER_H

#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include "operation.h"

using namespace std;

class Enhancer
{
public:
	Position pos;
	int id;
	Enhancer()
	{
	}

	Enhancer(int instart, int inend, string inchr)
	{
		pos.addpos(instart, inend, inchr);
	}

	Enhancer(string instr)
	{
		pos.addpos(instr);
	}

	int start();
	int end();
	string chr();

};

class Summit
{
public:
	string chr;
	int summit;
	Sample_id sid;
	string type;
};

class Narrowpeak   // 1000 rg
{
public:
	string chr;
	int start;
	int end;
	double intensity;
	set<int > summits;
	double epu_idx;  // e.g. 3.5: between epu3 and epu4;  4: epu4.
	int id;
	Narrowpeak()
	{
	}
	
	Narrowpeak( string inchr, int instart, int inend, double intens )
	{
		chr = inchr;
		start = instart;
		end = inend;
		intensity = intens;
	}
	
	map< Sample_id, map< string, vector<Summit > > > merged_summits;
	map< Sample_id, map< string, double > > merged_intensity;
	map< Sample_id, int > type_map;  // 0 none; 1 k4me1; 2 k27ac; 3 both.
};

class Broadpeak
{
public:
	string chr;
	int start;
	int end;
	
	set<int > summits;
	
	Broadpeak()
	{
	}
	Broadpeak( string inchr, int instart, int inend )
	{
		chr = inchr;
		start = instart;
		end = inend;
	}
};

class Enhancer_GW
{
public:
	vector< Narrowpeak > narrowpeak_Vec;
	map< string, map< pair<int, int >, Region_id > > chr_pos_narrowpeak_map;
	vector< Broadpeak > broadpeak_Vec;
	map< string, map< pair<int, int >, Region_id > > chr_pos_broadpeak_map;
	
	
};

class M_peak_mappings
{
public:
	map< Sample_id, vector<Region_id > > sample_peaks_map;
};

class Peak_Bank
{
public:	
	vector<Enhancer_GW > k4me1_GW_Vec;
	map<string, map<int, Sample_id > > cell_time_me1_idmap;
	vector<Enhancer_GW > k27ac_GW_Vec;
	map<string, map<int, Sample_id > > cell_time_ac_idmap;
	
	vector< Narrowpeak > M_narrowpeak_Vec;
	vector< M_peak_mappings > M_narrowpeak_me1_mappings;
	vector< M_peak_mappings > M_narrowpeak_ac_mappings;
	map< string, map< pair<int, int >, Region_id > > chr_pos_M_narrowpeak_map;
	
	vector< Broadpeak > M_broadpeak_Vec;
	vector< M_peak_mappings > M_broadpeak_me1_mappings;
	vector< M_peak_mappings > M_broadpeak_ac_mappings;
	map< string, map< pair<int, int >, Region_id > > chr_pos_M_broadpeak_map;

	void merge_narrowpeak();	
	void merge_broadpeak();
};


#endif

