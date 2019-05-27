#include "histonmark_dynamic.h"


void HMark_Cluster::merge_sample_pos( vector<Sample_Encode > &encode_Vec, string hmname )
{
	int start = -1;
	int end = -1;
	string chr = "";
	map<Sample_id, vector< Region_id > >::iterator ite = sample_region_map.begin();
	for ( ; ite != sample_region_map.end(); ++ite )
	{
		Sample_id sample = ite->first;
		vector< Region_id >::iterator seci = ite->second.begin();

		int addlength = 0;
		double addsignal = 0;
		for ( ; seci != ite->second.end(); ++seci )
		{
			Region_id region = *seci;
			chr = encode_Vec[sample].name_Histonmark_GW_map[hmname].histonmark_Vec[region].pos.chr;
			if ( start == -1 )
			{
				start = encode_Vec[sample].name_Histonmark_GW_map[hmname].histonmark_Vec[region].pos.start;
				end = encode_Vec[sample].name_Histonmark_GW_map[hmname].histonmark_Vec[region].pos.end;
			} else
			{
				int ts = encode_Vec[sample].name_Histonmark_GW_map[hmname].histonmark_Vec[region].pos.start;
				int te = encode_Vec[sample].name_Histonmark_GW_map[hmname].histonmark_Vec[region].pos.end;
				if ( ts < start )
					start = ts;
				if ( te > end )
					end = te;
			}

			addlength += encode_Vec[sample].name_Histonmark_GW_map[hmname].histonmark_Vec[region].getlength();
			addsignal += encode_Vec[sample].name_Histonmark_GW_map[hmname].histonmark_Vec[region].signal;
		}
		
		sample_signal_map[sample] = addsignal/(int)ite->second.size();
		sample_length_map[sample] = addlength;
	}

	pos.addpos( start, end, chr );
}

void HMark_Cluster::generatedynamics(std::vector<pair<Sample_id,Sample_id> > &stages)
{
	dynamics.clear();
	vector< int > onoffmatrix;
	vector< int > updownmatrix;
	vector< int > widenarrowmatrix;
	for ( size_t i = 0; i < stages.size(); ++i )
	{
		Sample_id preid = stages[i].first;
		Sample_id afid = stages[i].second;
		int pre_length = 0;
		int af_length = 0;
		double pre_signal = 0;
		double af_signal = 0;
		if ( sample_length_map.find( preid ) != sample_length_map.end() )
		{
			pre_length = sample_length_map[preid];
		}
		if ( sample_length_map.find( afid ) != sample_length_map.end() )
		{
			af_length = sample_length_map[afid];
		}
		if ( sample_signal_map.find( preid ) != sample_signal_map.end() )
		{
			pre_signal = sample_signal_map[preid];
		}
		if ( sample_signal_map.find( afid ) != sample_signal_map.end() )
		{
			af_signal = sample_signal_map[afid];
		}
		

		Dynamics dob;
		dob.afid = afid;
		dob.preid = preid;
		if ( markon( pre_length, af_length ) )
			dob.onoff = 1;
		else if ( markoff( pre_length, af_length ) )
			dob.onoff = 2;
		else
			dob.onoff = 0;

		if ( markwide( pre_length, af_length ) )
			dob.widenarrow = 1;
		else if ( marknarrow( pre_length, af_length ) )
			dob.widenarrow = 2;
		else
			dob.widenarrow = 0;
			
		if ( markup( pre_length, af_length ) )
			dob.updown = 1;
		else if ( markdown( pre_length, af_length ) )
			dob.updown = 2;
		else
			dob.updown = 0;
			
		
		dynamics.push_back(dob);

		onoffmatrix.push_back(dob.onoff);
		updownmatrix.push_back(dob.updown);
		widenarrowmatrix.push_back(dob.widenarrow);
	}

	onoff_tag = matrixtrans( (int)onoffmatrix.size(), 3, onoffmatrix );
	updown_tag = matrixtrans( (int)updownmatrix.size(), 3, updownmatrix );
	widenarrow_tag = matrixtrans( (int)widenarrowmatrix.size(), 3, widenarrowmatrix );
}

int HMark_Cluster::getlength( Sample_id id )
{
	int tlength = 0;
	if ( sample_length_map.find(id ) != sample_length_map.end() )
	{
		tlength = sample_length_map[id];
	}
	return tlength;
}

double HMark_Cluster::getsignal( Sample_id id )
{
	double tsignal = 0;
	if ( sample_signal_map.find(id ) != sample_signal_map.end() )
	{
		tsignal = sample_signal_map[id];
	}
	return tsignal;
}

void get_hm_cluster_gw( vector<Sample_Encode > &encode_Vec, string &hmname, Clustered_HMark_GW &cluster_gw )
{
	map<string, multimap<pair<int,int>, pair<Sample_id, Region_id > > > chr_pos_region_map;
	for ( size_t i = 0; i < encode_Vec.size(); ++i )
	{
		Sample_id sample = i;
		if ( encode_Vec[i].name_Histonmark_GW_map.find(hmname ) == encode_Vec[i].name_Histonmark_GW_map.end() )
		{
			cout<<"Warning: sample "<<i<<" not have hm data "<<hmname<<endl;
			continue;
		}
		map<string, map< pair<int, int >, Region_id > >::iterator ite = encode_Vec[i].name_Histonmark_GW_map[hmname].chr_pos_histonmark_map.begin();
		for ( ; ite != encode_Vec[i].name_Histonmark_GW_map[hmname].chr_pos_histonmark_map.end(); ++ite )
		{
			string chr = ite->first;
			for ( map< pair<int, int >, Region_id >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
			{
				chr_pos_region_map[chr].insert( make_pair( seci->first, make_pair( sample, seci->second ) ) );
			}
		}
	}


	map<string, multimap<pair<int,int>, pair<Sample_id, Region_id > > >::iterator ite = chr_pos_region_map.begin();
	for ( ; ite != chr_pos_region_map.end(); ++ite )
	{
		string chr = ite->first;
		HMark_Cluster tcluster;
		cluster_gw.clusters.push_back( tcluster );
		multimap<pair<int,int>, pair<Sample_id, Region_id > >::iterator posite = ite->second.begin();
		Sample_id tsample = posite->second.first;
		Region_id tregion = posite->second.second;
		cluster_gw.clusters.back().sample_region_map[tsample].push_back(tregion);
		Position tpos( posite->first.first, posite->first.second, chr );
		++posite;
		for ( ; posite != ite->second.end(); ++posite )
		{
			tsample = posite->second.first;
			tregion = posite->second.second;
			Position nextpos( posite->first.first, posite->first.second, chr );
			if ( overlaptest( tpos, nextpos, 0 ) )
			{
				cluster_gw.clusters.back().sample_region_map[tsample].push_back(tregion);
				if ( nextpos.end > tpos.end )
					tpos = nextpos;
			} else
			{
				cluster_gw.clusters.back().merge_sample_pos(encode_Vec, hmname );
				int pstart = cluster_gw.clusters.back().pos.start;
				int pend = cluster_gw.clusters.back().pos.end;
				cluster_gw.chr_pos_cluster_map[chr][make_pair(pstart, pend)] = cluster_gw.clusters.size()-1;

				// start next cluster
				cluster_gw.clusters.push_back(tcluster);
				cluster_gw.clusters.back().sample_region_map[tsample].push_back(tregion);
				tpos = nextpos;
			}
		}
		cluster_gw.clusters.back().merge_sample_pos(encode_Vec, hmname );
		int pstart = cluster_gw.clusters.back().pos.start;
		int pend = cluster_gw.clusters.back().pos.end;
		cluster_gw.chr_pos_cluster_map[chr][make_pair(pstart, pend)] = cluster_gw.clusters.size()-1;
	}
}

void get_hm_dynamics_gw( Clustered_HMark_GW & cluster_gw, 
						 vector<pair<Sample_id, Sample_id > > &stages,
						map< int, set< Region_id > >& onofftag_cluster_map,
						map< int, set< Region_id > >& widenarrowtag_cluster_map,
						map< int, set< Region_id > >& updowntag_cluster_map )
{
	for ( size_t i = 0; i < cluster_gw.clusters.size(); ++i )
	{
		Region_id id = i;
		cluster_gw.clusters[i].generatedynamics( stages );
		int onofftag = cluster_gw.clusters[i].onoff_tag;
		int widenarrowtag = cluster_gw.clusters[i].widenarrow_tag;
		int updowntag = cluster_gw.clusters[i].updown_tag;
		onofftag_cluster_map[onofftag].insert( id );
		widenarrowtag_cluster_map[widenarrowtag].insert( id );
		updowntag_cluster_map[updowntag].insert(id);
	}
}

string hmonofftagmeaning( int tag, int stage )
{
	vector< int > matrix1;
	matrix1 = revtransmatrix( stage, 3, tag );
	string m = "";
	for ( size_t i = 0; i < matrix1.size(); ++i )
	{
		if ( matrix1[i] == 0 )
		{
			m += "No";
		} else if ( matrix1[i] == 1 )
		{
			m += "On";
		} else if ( matrix1[i] == 2 )
		{
			m += "Off";
		}
		if ( i != matrix1.size()-1)
		{
			m += '_';
		}
	}
	return m;
}

string hmwidenarrowtagmeaning( int tag, int stage )
{
	vector< int > matrix1;
	matrix1 = revtransmatrix( stage, 3, tag );
	string m = "";
	for ( size_t i = 0; i < matrix1.size(); ++i )
	{
		if ( matrix1[i] == 0 )
		{
			m += "No";
		} else if ( matrix1[i] == 1 )
		{
			m += "Wide";
		} else if ( matrix1[i] == 2 )
		{
			m += "Narrow";
		}
		if ( i != matrix1.size()-1)
		{
			m += '_';
		}
	}
	return m;
}


string hmupdowntagmeaning( int tag, int stage )
{
	vector< int > matrix1;
	matrix1 = revtransmatrix( stage, 3, tag );
	string m = "";
	for ( size_t i = 0; i < matrix1.size(); ++i )
	{
		if ( matrix1[i] == 0 )
		{
			m += "No";
		} else if ( matrix1[i] == 1 )
		{
			m += "up";
		} else if ( matrix1[i] == 2 )
		{
			m += "down";
		}
		if ( i != matrix1.size()-1)
		{
			m += '_';
		}
	}
	return m;
}



