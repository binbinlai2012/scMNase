#include "tfbinding_dynamic.h"

void TFBinding_Cluster::merge_sample_pos( vector<Sample_Encode > &encode_Vec, string tfname )
{
	int start = -1;
	int end = -1;
	string chr = "";
	map<Sample_id, vector< Region_id > >::iterator ite = sample_region_map.begin();
	for ( ; ite != sample_region_map.end(); ++ite )
	{
		Sample_id sample = ite->first;
		vector< Region_id >::iterator seci = ite->second.begin();

		int addrcount = 0;
		double addfoldchange = 0;
		for ( ; seci != ite->second.end(); ++seci )
		{
			Region_id region = *seci;
			chr = encode_Vec[sample].name_TFBinding_GW_map[tfname].binding_Vec[region].pos.chr;
			if ( start == -1 )
			{
				start = encode_Vec[sample].name_TFBinding_GW_map[tfname].binding_Vec[region].pos.start;
				end = encode_Vec[sample].name_TFBinding_GW_map[tfname].binding_Vec[region].pos.end;
			} else
			{
				int ts = encode_Vec[sample].name_TFBinding_GW_map[tfname].binding_Vec[region].pos.start;
				int te = encode_Vec[sample].name_TFBinding_GW_map[tfname].binding_Vec[region].pos.end;
				if ( ts < start )
					start = ts;
				if ( te > end )
					end = te;
			}

			addrcount += encode_Vec[sample].name_TFBinding_GW_map[tfname].binding_Vec[region].rcount;
			addfoldchange += encode_Vec[sample].name_TFBinding_GW_map[tfname].binding_Vec[region].foldchange;
		}

		double avgfoldchange = 0;
		if ( !ite->second.empty() )
			avgfoldchange = addfoldchange / (int)ite->second.size();

		sample_rcount_map[sample] = addrcount;
		sample_foldchange_map[sample] = avgfoldchange;

	}

	pos.addpos( start, end, chr );
}

int TFBinding_Cluster::getrcount( Sample_id id )
{
	int rcount = 0;
	if ( sample_rcount_map.find( id ) != sample_rcount_map.end() )
	{
		rcount = sample_rcount_map[id];
	}
	return rcount;
}

double TFBinding_Cluster::getfoldchange( Sample_id id )
{
	double fdch = 0;
	if ( sample_foldchange_map.find( id ) != sample_foldchange_map.end() )
	{
		fdch = sample_foldchange_map[id];
	}
	return fdch;
}

void TFBinding_Cluster::generatedynamics( vector<pair<Sample_id, Sample_id > > &stages )
{
	dynamics.clear();
	vector< int > onoffmatrix;
	vector< int > updownmatrix;
	for ( size_t i = 0; i < stages.size(); ++i )
	{
		Sample_id preid = stages[i].first;
		Sample_id afid = stages[i].second;
		int pre_rcount = 0;
		int af_rcount = 0;
		double pre_foldchange = 0;
		double af_foldchange = 0;
		if ( sample_rcount_map.find( preid ) != sample_rcount_map.end() )
		{
			pre_rcount = sample_rcount_map[preid];
		}
		if ( sample_rcount_map.find( afid ) != sample_rcount_map.end() )
		{
			af_rcount = sample_rcount_map[afid];
		}
		if ( sample_foldchange_map.find( preid ) != sample_foldchange_map.end() )
		{
			pre_foldchange = sample_foldchange_map[preid];
		}
		if ( sample_foldchange_map.find( afid ) == sample_foldchange_map.end() )
		{
			af_foldchange = sample_foldchange_map[afid];
		}

		Dynamics dob;
		dob.afid = afid;
		dob.preid = preid;
		if ( signalon( pre_rcount, af_rcount ) )
			dob.onoff = 1;
		else if ( signaloff( pre_rcount, af_rcount ) )
			dob.onoff = 2;
		else
			dob.onoff = 0;
			
		if ( signalup( pre_rcount, af_rcount ) )
			dob.updown = 1;
		else if ( signaldown( pre_rcount, af_rcount ) )
			dob.updown = 2;
		else
			dob.updown = 0;

	/*	if ( singalup( pre_rcount, af_rcount, pre_foldchange, af_foldchange ) )
			dob.updown = 1;
		else if ( singaldown( pre_rcount, af_rcount, pre_foldchange, af_foldchange ) )
			dob.updown = 2;
		else
			dob.updown = 0; */
		dynamics.push_back(dob);

		onoffmatrix.push_back(dob.onoff);
		updownmatrix.push_back(dob.updown);
	}

	onoff_tag = matrixtrans( (int)onoffmatrix.size(), 3, onoffmatrix );
	updown_tag = matrixtrans( (int)updownmatrix.size(), 3, updownmatrix );

}

void get_tfb_cluster_gw( vector<Sample_Encode > &encode_Vec, string &tfname, int overlap_extension, Clustered_TFBinding_GW &cluster_gw )
{
	map<string, multimap<pair<int,int>, pair<Sample_id, Region_id > > > chr_pos_region_map;
	for ( size_t i = 0; i < encode_Vec.size(); ++i )
	{
		Sample_id sample = i;
		if ( encode_Vec[i].name_TFBinding_GW_map.find(tfname ) == encode_Vec[i].name_TFBinding_GW_map.end() )
		{
			cout<<"Warning: sample "<<i<<" not have tf data "<<tfname<<endl;
			continue;
		}
		map<string, map< pair<int, int >, Region_id > >::iterator ite = encode_Vec[i].name_TFBinding_GW_map[tfname].chr_pos_binding_map.begin();
		for ( ; ite != encode_Vec[i].name_TFBinding_GW_map[tfname].chr_pos_binding_map.end(); ++ite )
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
		TFBinding_Cluster tcluster;
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
			if ( overlaptest( tpos, nextpos, overlap_extension ) )
			{
				cluster_gw.clusters.back().sample_region_map[tsample].push_back(tregion);
				if ( nextpos.end > tpos.end )
					tpos = nextpos;
			} else
			{
				cluster_gw.clusters.back().merge_sample_pos(encode_Vec, tfname );
				int pstart = cluster_gw.clusters.back().pos.start;
				int pend = cluster_gw.clusters.back().pos.end;
				cluster_gw.chr_pos_cluster_map[chr][make_pair(pstart, pend)] = cluster_gw.clusters.size()-1;

				// start next cluster
				cluster_gw.clusters.push_back(tcluster);
				cluster_gw.clusters.back().sample_region_map[tsample].push_back(tregion);
				tpos = nextpos;
			}
		}
		cluster_gw.clusters.back().merge_sample_pos(encode_Vec, tfname );
		int pstart = cluster_gw.clusters.back().pos.start;
		int pend = cluster_gw.clusters.back().pos.end;
		cluster_gw.chr_pos_cluster_map[chr][make_pair(pstart, pend)] = cluster_gw.clusters.size()-1;
	}

}

void get_tfb_dynamics_gw( Clustered_TFBinding_GW & cluster_gw, 
						 vector<pair<Sample_id, Sample_id > > &stages,
						map< int, set< Region_id > > &onofftag_cluster_map,
						map< int, set< Region_id > > &updowntag_cluster_map )
{

	for ( size_t i = 0; i < cluster_gw.clusters.size(); ++i )
	{
		Region_id id = i;
		cluster_gw.clusters[i].generatedynamics( stages );
		int onofftag = cluster_gw.clusters[i].onoff_tag;
		int updowntag = cluster_gw.clusters[i].updown_tag;
		onofftag_cluster_map[onofftag].insert( id );
		updowntag_cluster_map[updowntag].insert( id );
	}
}





