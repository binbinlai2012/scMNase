#ifndef TF_H
#define TF_H

#include "operation.h"

using namespace std;

class TFBinding
{
public:
	Position pos;
	int signal;
	int rcount;
	double foldchange;
	double rank;
	TFBinding()
	{
	}

	TFBinding(int instart, int inend, string inchr)
	{
		pos.addpos(instart, inend, inchr);
	}

	TFBinding(string instr)
	{
		pos.addpos(instr);
	}

	int start();
	int end();
	int win();
	string chr();
	void addsignal(int insignal);
	void addrcount(int inrc );
	void addfoldchange(double foldchg);
};

class TFBinding_GW
{
public:
	vector<TFBinding > binding_Vec;
	map<string, map< pair<int, int >, Region_id > > chr_pos_binding_map;
	int libsize;
	multimap<int, Region_id > rc_region_map;
	
	void normalize_rc();
	void ranknormalize_rc();
};

class TFBindingGPS
{
public:
	string chr;
	int site;
	double qvalue;
	double pvalue;
	
	TFBindingGPS()
	{
	}

	TFBindingGPS(string inchr, int insite, double inq, double inp )
	{
		chr = inchr;
		site = insite;
		qvalue = inq;
		pvalue = inp;
	}
};

class TFBindingGPS_GW
{
public:
	vector<TFBindingGPS > binding_Vec;
	map<string, map<int, Region_id > > chr_pos_binding_map;
	
};

class TFBindingMACS
{
public:
	string chr;
	string name;
	int summit;
	double qvalue;
	int start;
	int end;
	bool peak;
	
	TFBindingMACS()
	{
	}

	TFBindingMACS(string inchr, int insite, double inq, string inname )
	{
		chr = inchr;
		summit = insite;
		peak = false;
		qvalue = inq;
		name = inname;
	}
};

class TFBindingMACS_GW
{
public:
	vector<TFBindingMACS > binding_Vec;
	map<string, map<int, Region_id > > chr_pos_binding_map;
	map<string, Region_id > peakname_binding_map;
};



#endif

