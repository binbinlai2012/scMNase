#ifndef HM_DYNAM_H
#define HM_DYNAM_H

#include "databank.h"
#include "operation.h"

using namespace std;

class HMark_Cluster
{
public:
	map<Sample_id, vector< Region_id > > sample_region_map;

	map<Sample_id, int > sample_length_map;
	map<Sample_id, double > sample_signal_map;
	
	Position pos;

	void merge_sample_pos( vector<Sample_Encode > &encode_Vec, string hmname );

	vector<Dynamics > dynamics;
	int onoff_tag;
	int widenarrow_tag;
	int updown_tag;

	void generatedynamics( vector<pair<Sample_id, Sample_id > > &stages );

	int getlength( Sample_id id );
	double getsignal( Sample_id id );

};

class Clustered_HMark_GW
{
public:
	vector< HMark_Cluster > clusters;
	map<string, map<pair<int, int>, Region_id > > chr_pos_cluster_map;
	
	map< int, Region_id > onofftag_cluster_map;
	map< int, Region_id > widenarrowtag_cluster_map;
};


void get_hm_cluster_gw( vector<Sample_Encode > &encode_Vec, string &hmname, Clustered_HMark_GW &cluster_gw );

void get_hm_dynamics_gw( Clustered_HMark_GW & cluster_gw, 
						 vector<pair<Sample_id, Sample_id > > &stages,
						map< int, set< Region_id > >& onofftag_cluster_map,
						map< int, set< Region_id > >& widenarrowtag_cluster_map,
						map< int, set< Region_id > >& updowntag_cluster_map );

string hmonofftagmeaning( int tag, int stage );
string hmwidenarrowtagmeaning( int tag, int stage );
string hmupdowntagmeaning( int tag, int stage );

#endif

