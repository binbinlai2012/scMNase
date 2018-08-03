#ifndef NUCL_SPACE_H
#define NUCL_SPACE_H

#include "operation.h"

using namespace std;

class Cell_Nucl_Space_Bank
{
public:
	
	map<string, map<string, set<int > > > cell_nucl_pos;
	map<string, map<string, set<int > > > cell_les78_pos;
	
	void readinnucl( string infile );
	void readinnucl_bunch( string infile );
	void readinles78( string infile );
	void readinles78_bunch( string infile );
	void readinnucl_cut_olp( string infile, 
		vector< map<string, set<int> > > &nucl_group, 
		vector< map<string, vector<pair<int, int > > > > &nucl_cuts );
	void cal_space( map<string, vector<pair<int, int > > > &region, map<string, vector<int > > &cell_distances, int lower, int upper );
	// To reduce ambiguity, select anchor with high PS. In addition, no deviated nucl within single cell.
	void cal_space( map<string, vector<pair<int, int > > > &region, map<string, vector<int > > &cell_distances, int lower, int upper,
		map<string, set<int > > &anchors );
	void cal_space_wholechr( map<string, vector<int > > &cell_distances, int lower, int upper );
	
	// cal pair of nucl and les78 around a center in one cell. If nucl is absent, use -10000 value. les78 is the same.
	void cal_nucl_les78_in_cell( map<string, vector<int > > &centers,
		 int up, int down, map<string, vector<pair<int, int > > > &cell_nucl_les78 );
		 
	// cal nucl space between the -1 highest peak and downstream nucls around dhs. and also +1 and upstream 
	void cal_space_1peak_othersidedownstream( map<string, vector<int > > &centers,
		int peakstart, int peakend, int down, map<string, vector<int > > &cell_space,
		vector<vector<int > > &dist_cell_num );
	void cal_distcellnum_1peak_othersidedownstream( map<string, vector<int > > &centers,
		int peakstart, int peakend, int down,   // 80 255
		map<string, map<int, vector<int > > > &dist_cell_num,
		map<string, map<int, vector<vector<string > > > > &dist_cell_list );
	void cal_nuclvar_1peak_othersidedownstream( map<string, vector<int > > &centers,
		int peakstart, int peakend, int down,   // 80 255
		vector<int > &allelevar_onlyregion1,
		vector<int > &allelevar_onlyregion2,
		vector<int > &allelevar_both );
	void cal_distnuclnum_1peak_othersidedownstream( map<string, vector<int > > &centers,
		int peakstart, int peakend, int down,   // 80 255
		map<string, map<int, vector<int > > > &dist_nucl_num );
	void cal_space_1peak_othersidedownstream_les78filter( map<string, vector<int > > &centers,
		int peakstart, int peakend, int down, map<string, vector<int > > &cell_space );   // require les78 occupy
	void get_space_leftpeak_down( vector<int> &pos, int peakstart, int peakend, int upordown, vector<int> &spaces );
	void get_nuclpos_leftpeak_down( vector<int> &pos, int peakstart, int peakend, int upordown, map<int, int> &relpos_count );
	void get_nuclpos_leftpeak_down_control( vector<int> &pos, int peakstart, int peakend, int upordown, 
	map<int, int> &relpos_count, map<int, int > &relpos_count_2 );
	
	void cal_space_1peak_samesidedownstream( map<string, vector<int > > &centers,
		int peakstart, int peakend, int down, map<string, vector<int > > &cell_space );
		
	void cal_space_1N_downstream( map<string, map<pair<int, int >,  char > > &N1region,
		int down, map<string, vector<int > > &cell_space );
		
	void cal_relpos_count_1N_downstream( map<string, map<pair<int, int >, pair<char, string>  > > &N1region,
		int down, map<string, map<int, int > > &gene_relpos_count );
	void cal_truepos_count_1N_downstream( map<string, map<int, pair<char, string>  > > &tssregion, 
		int up, int down, map<string, map<int, int > > &gene_truepos_count );
	
	void cal_relpos_count_1N_downstream_control( map<string, map<pair<int, int >, pair<char, string>  > > &N1region,
		int down, map<string, map<int, int > > &gene_relpos_count,     // with cell nucl as ref
		map<string, map<int, int > > &gene_relpos_count2 );            // with ref nucl as ref
	
	// cal space between the -1 highest peak nucl and middle les78 around dhs.
	void cal_space_1peak_nucl_middle_les78( map<string, vector<int > > &centers, 
		int peakstart, int peakend, int middle_board, map<string, vector<int > > &cell_space,
		vector<vector<int > > &dist_cell_num );
	void get_space_peaknucl_middleles78( vector<int > &posnucl, vector<int > &posles78,
		int lpeakstart, int lpeakend, int rpeakstart, int rpeakend, int middle_lbound, int middle_rbound, vector<int > &spaces );
};

void get_valid_space( set<int >& nuclpos, int lower, int upper, vector<int > &space );
void get_valid_space( set<int >& anchors, set<int >& nuclpos, int lower, int upper, vector<int > &space );
void pool_space( map<string, vector<int > > &cell_distances, vector<int > &pooled_distances );
void cal_ratio( map<string, vector<int > > &cell_distances, map<string, double > &cell_ratio, int lower, int upper );

void getshift_amongcells( map<string, map<string, set<int > > > &cell_nucl_pos, 
	map<string, map<int, vector< int> > > &chr_allele_shift );
void getshift_withincell( map<string, map<string, set<int > > > &cell_nucl_pos, 
	map<string, map<int, vector<int > > > &chr_allele_shift );

void pool_nuclpos(Cell_Nucl_Space_Bank &bank, map<string, set<int > > &nucl_pos );
void getlendis( map<string, set<int > > &nucl_pos, string prefix, map<string, set<int> > &phased_nucl);
void getlendis( map<string, set<int > > &nucl_pos, string prefix, map<string, set<pair< int, int > > > &region);

void getnuclvar( vector<int > &relpos, vector<int > &var  );
	
#endif

