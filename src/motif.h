#ifndef MOTIF_H
#define MOTIF_H

#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <fstream>

#include "operation.h"

class MotifHit
{
public:
	string name;
	string chr;
	int center;
	char strand;
	
	MotifHit()
	{
	}

	MotifHit(string inchr, int incenter, char instrand )
	{
		chr = inchr;
		center = incenter;
		strand = instrand;
	}
	
};

class MotifHit_GW
{
public:
	vector<MotifHit > motifhit_Vec;
	map<string, map< int, Region_id > > chr_pos_motifhit_map;
	
	vector<Region_id > filter_motifhit( map<string, vector< pair<int, int > > > &regions );
	
};


#endif

