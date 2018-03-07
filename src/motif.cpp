#include "motif.h"

vector<Region_id > MotifHit_GW::filter_motifhit( map<string, vector< pair<int, int > > > &regions )
{
	vector<Region_id > rv;
	for ( map<string, vector< pair<int, int > > >::iterator ite = regions.begin(); ite != regions.end(); ++ite )
	{
		string chr = ite->first;
		if ( chr_pos_motifhit_map.find(chr) == chr_pos_motifhit_map.end() )
		{
			continue;
			
		}
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			int st = ite->second[i].first;
			int ed = ite->second[i].second;
			for ( map< int, Region_id >::iterator seci = chr_pos_motifhit_map[chr].begin(); seci != chr_pos_motifhit_map[chr].end(); ++seci )
			{
				int pos = seci->first;
				if ( pos < st )
					continue;
				if ( pos > ed )
					break;
				rv.push_back( seci->second );
			}
		}
	}
	
	return rv;
	
}



