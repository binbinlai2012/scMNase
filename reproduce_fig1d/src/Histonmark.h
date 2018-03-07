#ifndef HISTONMARK_H
#define HISTONMARK_H
#include "operation.h"
#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>

using namespace std;

class Histonmark
{
public:
	Position pos;
	double signal;
	double rank;
	double rcount;
	

	Histonmark()
	{
	}

	Histonmark(int instart, int inend, string inchr)
	{
		pos.addpos(instart, inend, inchr);
	}

	Histonmark(string instr)
	{
		pos.addpos(instr);
	}

	int start();
	int end();
	int win();
	string chr();
	void addsignal(double insignal);
	void addrcount(double inr);
	int getlength();
	
	
};

class Histonmark_GW
{
public:
	vector<Histonmark > histonmark_Vec;
	map<string, map< pair<int, int >, Region_id > > chr_pos_histonmark_map;
	int libsize;
	multimap<double, Region_id > sig_region_map;
	
	void ranknormalize_signal();
	void normalize_rc();
};



#endif

