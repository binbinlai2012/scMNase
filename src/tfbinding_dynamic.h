#ifndef TF_DYNAM_H
#define TF_DYNAM_H

#include "databank.h"
#include "operation.h"

using namespace std;



class TFBinding_Cluster
{
public:
	map<Sample_id, vector< Region_id > > sample_region_map;

	map<Sample_id, int > sample_rcount_map;
	map<Sample_id, double > sample_foldchange_map;

	Position pos;

	void merge_sample_pos( vector<Sample_Encode > &encode_Vec, string tfname );

	vector<Dynamics > dynamics;
	int onoff_tag;
	int updown_tag;

	void generatedynamics( vector<pair<Sample_id, Sample_id > > &stages );

	int getrcount( Sample_id id );
	double getfoldchange( Sample_id id );

};

class Clustered_TFBinding_GW
{
public:
	vector< TFBinding_Cluster > clusters;
	map<string, map<pair<int, int>, Region_id > > chr_pos_cluster_map;
	
//	map< int, Region_id > onofftag_cluster_map;
//	map< int, Region_id > updowntag_cluster_map;
};

void get_tfb_cluster_gw( vector<Sample_Encode > &encode_Vec, string &tfname, int overlap_extension, Clustered_TFBinding_GW &cluster_gw );

void get_tfb_dynamics_gw( Clustered_TFBinding_GW & cluster_gw, 
						 vector<pair<Sample_id, Sample_id > > &stages,
						map< int, set< Region_id > >& onofftag_cluster_map,
						map< int, set< Region_id > >& updowntag_cluster_map );



#endif

