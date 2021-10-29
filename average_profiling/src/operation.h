#ifndef OPERATION_H
#define OPERATION_H

#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <algorithm>
#include <cstdlib>

using namespace std;

typedef size_t Region_id;
typedef size_t Sample_id;

vector<string > parse_string( string & instr, char spl );
vector<string > parse_string( string & instr );

class Position
{
public:
	int start;
	int end;
	string chr;
	Position()
	{
	}
	Position(int inst, int inend, string inchr)
	{
		start = inst;
		end = inend;
		chr = inchr;
	}
	Position(string str)
	{
		vector<string > splstr1 = parse_string(str, ':');
		if ( splstr1.size() != 2 )
		{
			cout<<"error parse position: str "<< str<<endl;
			exit(1);
		}
		chr = splstr1[0];
		vector<string > splstr2 = parse_string(splstr1[1], '+');
		if ( splstr2.size() != 2 )
		{
			cout<<"error parse position: str "<< str<<endl;
			exit(1);
		}
		start = atoi( splstr2[0].c_str() );
		end = atoi( splstr2[1].c_str() );
	}

	void addpos(int inst, int inend, string inchr);
	void addpos(string str);


};

class Stitched_Region
{
public:
	vector<pair<int, int > > regions;
	string chr;

	int getstart();
	int getend();

	
};

class Dynamics
{
public:
	Sample_id preid;
	Sample_id afid;
	double foldchange;
	int onoff;         // 0: no; 1: on; 2: off
	int updown;        // 0: no; 1: up; 3: down
	int widenarrow;
	int moreless;      
};

class Count_pool
{
public:
	
	vector<string > id_ve;
	vector< vector<int > > counts_table;
	vector< vector<double > > rpkm_table;
	map<string, map<int, char> > tss_strand_map;
	map<string, map<int, size_t> > tss_id_map;
	map<string, map<int, pair<int, int > > > tss_protect_map;
	map<string, map<pair<int, int>, char > > body_strand_map;
	map<string, map<pair<int, int>, size_t > > body_id_map;
	
	map<string, map<pair<int, int >, size_t > > region_id_map;
	map<string, map<pair<int, int >, char > > region_strand_map;
	vector< vector<pair<int, int > > > cut_reg_ve;
	void getregionidmap_tss( int ups, int downs, int win );
	void getregionidmap_body( int ups, int downs, int step );
	void getregionidmap_body_sp( int ups, int downs, int win, int step );
	void getregionidmap_body_sp2( int ups, int downs, int win, int step );
	void smooth_ave_rpkm( int stretch );
	void normrpkm( int total );
	void reverseminusstrand();
	void gettssprotectmap();
	
};

double log_2(double r );

bool overlaptest(Position &p1, Position &p2);

bool overlaptest(Position &p1, Position &p2, int extension);

bool overlaptest(pair<int, int> p1, pair<int, int> p2);

bool overlaptest(pair<int, int> p1, pair<int, int> p2, int extension);

bool overlaptest(pair<int, int> p1, map< pair<int, int >, Region_id > &posmap );

bool overlaptest(pair<int, int> p1, map< pair<int, int >, Region_id > &posmap, vector<Region_id > &idve );

bool overlaptest(pair<int, int> p1, map< pair<int, int >, Region_id > &posmap, vector<Region_id > &idve, int extension );

bool overlaptest( int site1, map<int, Region_id > &posmap, int extension, vector<int > &dis );

bool signalon( int rc1, int rc2 );

bool signaloff( int rc1, int rc2 );

bool singalup( int rc1, int rc2, double chg1, double chg2 );

bool singaldown( int rc1, int rc2, double chg1, double chg2 );

bool expressionup( double rc1, double rc2 );

bool signalup(int rc1, int rc2 );
bool signaldown(int rc1, int rc2 );

bool consregionmore( int rc1, int rc2 );
bool consregionless( int rc1, int rc2 );

bool markon( int length1, int length2 );
bool markoff( int length1, int length2 );
bool markon( double sig1, double sig2 );
bool markoff( double sig1, double sig2 );

bool markwide( int length1, int length2 );

bool marknarrow( int length1, int length2 );

bool markup( double signal1, double signal2 );
bool markdown( double signal1, double signal2 );

int matrixtrans( int stage, int elem, vector<int > &matrix );
vector<int > revtransmatrix( int stage, int elem, int digit );

void getposfrom( map<string, vector< pair<int, int> > > &region, 
				map<string, map<pair<int, int>, string > > &region_map );

void getposfrom( map<string, vector< pair<int, int> > > &region, 
				map<string, map<pair<int, int>, Region_id > > &region_map );

void filter_region_remain( map<string, vector< pair<int, int> > > &filtered, 
				   map<string, vector< pair<int, int > > > &ori_region,
				   map<string, vector< pair<int, int > > > &fr );

void filter_region_minus( map<string, vector< pair<int, int> > > &filtered, 
				   map<string, vector< pair<int, int > > > &ori_region,
				   map<string, vector< pair<int, int > > > &fr );

void filter_region_minus_remain( map<string, vector< pair<int, int> > > &overlapped, 
				   map<string, vector< pair<int, int> > > &nonoverlapped,
				   map<string, vector< pair<int, int > > > &ori_region,
				   map<string, vector< pair<int, int > > > &fr );

void stitch_region( map<string, vector<vector<pair<int, int> > > > &stitched,
				   map<string, vector<pair<int, int > > > &region,
				   int distance );
				   
void find_overlap_bycenterinregionplus(map<string, vector< pair<int, int> > > &overlapped, 
	map<string, vector< pair<int, int> > > &ori_region,
	map<string, vector< pair<int, int > > > &base_region,
	int exten );
	
void region_merge_naive( vector< map<string, vector< pair<int, int> > > > &reg_ve,
	map<string, vector< pair<int, int> > > &merged );

void stitch_region(map<string, vector<Stitched_Region > > &stitched,
				   map<string, vector<pair<int, int > > > &region,
				   int distance );

void cluster_region_CrSample( map<string, vector< vector<pair<int, int> > > > &clustered,
							 vector< map<string, vector<pair<int, int> > > > &posVec,
							 int overlap_extension );
							 
string inttostr(int i );

void outputtable( ofstream &outf, vector<vector<int> > &mat, vector<string > &rowname, vector<string> &colname );

void outputtable( ofstream &outf, vector<vector<double> > &mat, vector<string > &rowname, vector<string> &colname );

vector<vector<double> > transposmat( vector<vector<double> > &mat );

void normalizemat( vector<vector<double > > &mat, vector<double > &thr );

int findnearTSSinchr( set<int > &tssset, pair<int, int > region );

int findnearTSSinchr( set<int > &tssset, pair<int, int > region, vector< int > &containedTSS );

string onofftagmeaning( int tag, int stage );
string updowntagmeaning( int tag, int stage );

vector<double > getqualter( vector< vector<vector<double > > > &mat_vec );
vector<vector<double > > mergecol(vector< vector<vector<double > > > &mat_vec );

int intersection_num( vector<set<Region_id > > &sets );
void venn_quintuple(ofstream &outf, vector<set<Region_id > > &sets, vector<string > &names );
void venn_triple(ofstream &outf, vector<set<Region_id > > &sets, vector<string > &names, string outname );

void venn_threeset( vector<Region_id > &set1, vector<Region_id >& set2, 
	vector<Region_id > &set1_exc, vector<Region_id > &set2_exc, vector<Region_id >& overlap );
	
double match_rank( int rc, multimap<int, Region_id > &sig_r_map );

double match_rank( double sig, multimap<double, Region_id > &sig_r_map );

set<Region_id > vetoset(vector<Region_id > &ve );

vector<pair<int, int > > mergeragions( set<pair<int, int > > &r );

int shift_tagpos( int start, int end, char strand, int seg_len );

int gettotaltagnumber( map<string, vector<int > > &tagpos );

void assigntagpos( map<string, vector<int > > &tagpos, string chr, int start, int end, int win, int total, 
	vector<int > &counts, vector<double > &rpkm );
void assigntagpos( int tagpos, vector<pair<int, int > > &cut_reg, vector<int > &counts);
	
void averagetablecolumn( vector< double > &ave_ve, vector<vector<double > > & rpkm_table );

void assigncountpool( map<string, vector<int > > &tagpos, Count_pool &cpool );

#endif

