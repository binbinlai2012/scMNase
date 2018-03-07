#include "enhancer.h"

int Enhancer::start()
{
	return pos.start;
}

int Enhancer::end()
{
	return pos.end;
}

string Enhancer::chr()
{
	return pos.chr;
}

void Peak_Bank::merge_narrowpeak()
{
	// first determine the broadest range
	map< string, set<pair<int, int > > > range_pool;
	for ( size_t i = 0; i < k4me1_GW_Vec.size(); ++i )
	{
		for ( map< string, map< pair<int, int >, Region_id > >::iterator ite = k4me1_GW_Vec[i].chr_pos_narrowpeak_map.begin(); 
			ite != k4me1_GW_Vec[i].chr_pos_narrowpeak_map.end(); ++ite )
		{
			for ( map< pair<int, int >, Region_id >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
			{
				range_pool[ite->first].insert( si->first );
			}
		}
	}
	for ( size_t i = 0; i < k27ac_GW_Vec.size(); ++i )
	{
		for ( map< string, map< pair<int, int >, Region_id > >::iterator ite = k27ac_GW_Vec[i].chr_pos_narrowpeak_map.begin(); 
			ite != k27ac_GW_Vec[i].chr_pos_narrowpeak_map.end(); ++ite )
		{
			for ( map< pair<int, int >, Region_id >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
			{
				range_pool[ite->first].insert( si->first );
			}
		}
	}
	map< string, vector<pair<int, int > > > range_sep;
	for ( map< string, set<pair<int, int > > >::iterator ite = range_pool.begin(); ite != range_pool.end(); ++ite )
	{
		vector<pair<int, int > > sve = mergeragions( ite->second );
		range_sep.insert( make_pair(ite->first, sve ) );
	}
	
	//get merged range in each separate range
/*	class Sunit
	{
	public:
		string type;
		Sample_id sid;
		Region_id rid;
		int summit;
		int start;
		int end;
		
	};*/
	
	for ( map< string, vector<pair<int, int > > >::iterator ite = range_sep.begin(); ite != range_sep.end(); ++ite )
	{
		string chr = ite->first;
		vector< Summit > sn_ve;
		vector<multimap<int, size_t > > sep_summit_unit_map;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			multimap<int, size_t > nm;
			sep_summit_unit_map.push_back( nm );
		}
		
		// collect sunit for each separate ragion 
		for ( size_t i = 0; i < k4me1_GW_Vec.size(); ++i )
		{
			for ( map< pair<int, int >, Region_id >::iterator mi = k4me1_GW_Vec[i].chr_pos_narrowpeak_map[chr].begin(); 
				mi != k4me1_GW_Vec[i].chr_pos_narrowpeak_map[chr].end(); ++mi )
			{
				bool fd = false;
				size_t idx = 0;
				for ( size_t j = 0; j < ite->second.size(); ++j )
				{
					if ( mi->first.first >= ite->second[j].first && mi->first.second <= ite->second[j].second )
					{
						fd = true;
						idx = j;
						break;
					}
				}
				if ( !fd )
				{
					cout<<"error peak not find in merge peak"<<endl; exit(1);
				}
				for ( set<int>::iterator smi = k4me1_GW_Vec[i].narrowpeak_Vec[mi->second].summits.begin(); 
					smi != k4me1_GW_Vec[i].narrowpeak_Vec[mi->second].summits.end(); ++smi )
				{
					Summit sn;
					sn.type ="k4me1";  // k4me1
					sn.sid = i;
				//	sn.rid = mi->second;
					sn.summit = *smi;
					sn.chr = chr;
				//	sn.start = mi->first.first;
				//	sn.end = mi->first.second;
					sn_ve.push_back( sn );
					sep_summit_unit_map[idx].insert( make_pair( sn.summit, sn_ve.size()-1 ) );
				}
			}
		} 
		for ( size_t i = 0; i < k27ac_GW_Vec.size(); ++i )
		{
			for ( map< pair<int, int >, Region_id >::iterator mi = k27ac_GW_Vec[i].chr_pos_narrowpeak_map[chr].begin(); 
				mi != k27ac_GW_Vec[i].chr_pos_narrowpeak_map[chr].end(); ++mi )
			{
				bool fd = false;
				size_t idx = 0;
				for ( size_t j = 0; j < ite->second.size(); ++j )
				{
					if ( mi->first.first >= ite->second[j].first && mi->first.second <= ite->second[j].second )
					{
						fd = true;
						idx = j;
						break;
					}
				}
				if ( !fd )
				{
					cout<<"error peak not find in merge peak"<<endl; exit(1);
				}
				for ( set<int>::iterator smi = k27ac_GW_Vec[i].narrowpeak_Vec[mi->second].summits.begin(); 
					smi != k27ac_GW_Vec[i].narrowpeak_Vec[mi->second].summits.end(); ++smi )
				{
					Summit sn;
					sn.type ="k27ac";  // k4me1
					sn.sid = i;
				//	sn.rid = mi->second;
					sn.summit = *smi;
				//	sn.start = mi->first.first;
				//	sn.end = mi->first.second;
					sn.chr = chr;
					sn_ve.push_back( sn );
					sep_summit_unit_map[idx].insert( make_pair( sn.summit, sn_ve.size()-1 ) );
				}
			}
		} 
		
		// get merged narrowpeak
		// 1000bp adjacent bin as unit 
		for ( size_t i = 0; i < sep_summit_unit_map.size(); ++i )
		{
			vector<set< int > > summit_cl;
			set<int > subc;
			for ( multimap<int, size_t >::iterator si = sep_summit_unit_map[i].begin(); si != sep_summit_unit_map[i].end(); ++si )
			{
				if ( subc.empty() )
					subc.insert( si->first );
				else
				{
					if ( *subc.rbegin() + 1000 < si->first )
					{
						summit_cl.push_back( subc );
						subc.clear();
						subc.insert( si->first );
					} else
						subc.insert( si->first );
				}
			}
			if ( !subc.empty() )
				summit_cl.push_back( subc );
			
			for ( size_t j = 0; j < summit_cl.size(); ++j )
			{
				int lowers = *summit_cl[j].begin();
				int uppers = *summit_cl[j].rbegin();
				int len = uppers - lowers + 1;
				int wings = 1000 - (len % 1000);
				int start = lowers - wings/2;
				int steps = len / 1000 + 1;
				int end = start + (steps * 1000) - 1;
				map< size_t, Region_id > idx_rid_map;
				for ( size_t s = 0; (int)s < steps; ++s )
				{
					Narrowpeak np;
					np.chr = chr;
					np.start = start + s*1000;
					np.end = np.start + 999;
					M_narrowpeak_Vec.push_back( np );
					idx_rid_map.insert( make_pair(s, M_narrowpeak_Vec.size()-1 ) );
					chr_pos_M_narrowpeak_map[chr][make_pair(start, end)] = M_narrowpeak_Vec.size()-1;
				}
				for ( multimap<int, size_t >::iterator si = sep_summit_unit_map[i].begin(); si != sep_summit_unit_map[i].end(); ++si )
				{
					if ( si->first < lowers )
						continue;
					if ( si->first > uppers )
						break;
					int s = ( si->first - start ) / 1000;
					if ( s < 0 || s >= steps )
					{
						cout<<"error summit "<<endl; 
						cout<<uppers<<","<<lowers<<","<<si->first<<","<<s<<","<<steps<<","<<wings<<","<<start<<endl;
						exit(1);
					}
					Region_id rid = idx_rid_map[s];
					M_narrowpeak_Vec[rid].summits.insert( si->first );
					string type = sn_ve[si->second].type;
					Sample_id sid = sn_ve[si->second].sid;
					M_narrowpeak_Vec[rid].merged_summits[sid][type].push_back( sn_ve[si->second]);
					
				}
			}
		} 
	}
	
}

void Peak_Bank::merge_broadpeak()
{
	map< string, set<pair<int, int > > > range_pool;
	
//	cout<<"merge broadpeak"<<endl;
	for ( size_t i = 0; i < k4me1_GW_Vec.size(); ++i )
	{
		for ( map< string, map< pair<int, int >, Region_id > >::iterator ite = k4me1_GW_Vec[i].chr_pos_broadpeak_map.begin(); 
			ite != k4me1_GW_Vec[i].chr_pos_broadpeak_map.end(); ++ite )
		{
		//	cout<<ite->first<<","<<ite->second.size()<<endl;
			for ( map< pair<int, int >, Region_id >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
			{
				range_pool[ite->first].insert( si->first );
				
			}
		}
	}
	for ( size_t i = 0; i < k27ac_GW_Vec.size(); ++i )
	{
		for ( map< string, map< pair<int, int >, Region_id > >::iterator ite = k27ac_GW_Vec[i].chr_pos_broadpeak_map.begin(); 
			ite != k27ac_GW_Vec[i].chr_pos_broadpeak_map.end(); ++ite )
		{
			for ( map< pair<int, int >, Region_id >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
			{
				range_pool[ite->first].insert( si->first );
			}
		}
	}
	
	map< string, vector<pair<int, int > > > range_sep;
	for ( map< string, set<pair<int, int > > >::iterator ite = range_pool.begin(); ite != range_pool.end(); ++ite )
	{
	//	cout<<ite->first<<","<<ite->second.size()<<endl;
		vector<pair<int, int > > sve = mergeragions( ite->second );
		range_sep.insert( make_pair(ite->first, sve ) );
	}
	
	for ( map<string, vector<pair<int, int > > >::iterator ite = range_sep.begin(); ite != range_sep.end(); ++ite )
	{
	//	cout<<ite->first<<"\t"<<ite->second.size()<<endl;
		string chr = ite->first;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			Broadpeak bp( chr, ite->second[i].first, ite->second[i].second );
			M_broadpeak_Vec.push_back( bp );
			M_peak_mappings m1;
			M_broadpeak_me1_mappings.push_back( m1 );
			M_peak_mappings m2;
			M_broadpeak_ac_mappings.push_back( m2 );
			chr_pos_M_broadpeak_map[chr][make_pair(bp.start, bp.end) ] = M_broadpeak_Vec.size()-1;
			
		}
		
		// mapping
		for (  size_t i = 0; i < k4me1_GW_Vec.size(); ++i )
		{
			for ( map< pair<int, int >, Region_id >::iterator mi = k4me1_GW_Vec[i].chr_pos_broadpeak_map[chr].begin(); 
				mi != k4me1_GW_Vec[i].chr_pos_broadpeak_map[chr].end(); ++mi )
			{
				bool fd = false;
				size_t idx = 0;
				for ( map< pair<int, int >, Region_id >::iterator si = chr_pos_M_broadpeak_map[chr].begin(); 
					si != chr_pos_M_broadpeak_map[chr].end(); ++si )
				{
					if ( si->first.first <= mi->first.first && si->first.second >= mi->first.second )
					{
						fd = true;
						idx = si->second;
						break;
					}
					
				}
				if ( !fd )
				{
					cout<<"error mapping broad peak"<<endl; exit(1);
				}
				if ( fd )
				{
					M_broadpeak_me1_mappings[idx].sample_peaks_map[i].push_back( mi->second );
				}
			}
		}
		for (  size_t i = 0; i < k27ac_GW_Vec.size(); ++i )
		{
			for ( map< pair<int, int >, Region_id >::iterator mi = k27ac_GW_Vec[i].chr_pos_broadpeak_map[chr].begin(); 
				mi != k27ac_GW_Vec[i].chr_pos_broadpeak_map[chr].end(); ++mi )
			{
				bool fd = false;
				size_t idx = 0;
				for ( map< pair<int, int >, Region_id >::iterator si = chr_pos_M_broadpeak_map[chr].begin(); 
					si != chr_pos_M_broadpeak_map[chr].end(); ++si )
				{
					if ( si->first.first <= mi->first.first && si->first.second >= mi->first.second )
					{
						fd = true;
						idx = si->second;
						break;
					}
					
				}
				if ( !fd )
				{
					cout<<"error mapping broad peak"<<endl; exit(1);
				}
				if ( fd )
				{
					M_broadpeak_ac_mappings[idx].sample_peaks_map[i].push_back( mi->second );
				}
			}
		}
	}
}







