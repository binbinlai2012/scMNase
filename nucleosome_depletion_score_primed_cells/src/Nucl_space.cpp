#include "Nucl_space.h"

void Cell_Nucl_Space_Bank::readinnucl( string infile )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	
	
	vector<string> filenames = parse_string(infile, '/');
	string name_s1 = filenames.back();
	vector<string> filenames2 = parse_string( name_s1, '.');
	string cell = filenames2[0];
	cout<<"cell "<<cell<<endl;
	string line;
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
			
		vector<string > parsed_items = parse_string(line);
	/*	if ( parsed_items.size() < 3 )
		{
			cout<<"error1 unexpected bed line "<<line<<endl;
			exit(1);
		} */
		
	
		string chr = parsed_items[0];
		/*	if ( chr.substr(0,3) != "chr" )
		{
			cout<<"error2 unexpected bed line "<<line<<endl;
			exit(1);
		} */
		int start = atoi(parsed_items[1].c_str());
		int end = atoi(parsed_items[2].c_str());
		int mid = start + (end-start)/2;
		cell_nucl_pos[cell][chr].insert( mid );
		
	}
	inf.close();

}

void Cell_Nucl_Space_Bank::readinnucl_cut_olp( string infile, 
		vector< map<string, set<int> > > &nucl_group, 
		vector< map<string, vector<pair<int, int > > > > &nucl_cuts )
{
	
	
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	
	for ( size_t i = 0; i < nucl_group.size(); ++i )
	{
		map<string, vector<pair<int, int > > > nul;
		nucl_cuts.push_back(nul);
	}
	
	string line;
	int k = 0;
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
			
		vector<string > parsed_items = parse_string(line);
	/*	if ( parsed_items.size() < 3 )
		{
			cout<<"error1 unexpected bed line "<<line<<endl;
			exit(1);
		} */
		string chr = parsed_items[0];
	/*	if ( chr.substr(0,3) != "chr" )
		{
			cout<<"error2 unexpected bed line "<<line<<endl;
			exit(1);
		} */
		
		int start = atoi(parsed_items[1].c_str());
		int end = atoi(parsed_items[2].c_str());
		int mid = start + (end-start)/2;
		
		for ( size_t i = 0; i < nucl_group.size(); ++i )
		{
			if ( nucl_group[i][chr].find(mid) != nucl_group[i][chr].end() )
			{
				nucl_cuts[i][chr].push_back( make_pair(start, end ) );
				break;
			}
		}
	}
	inf.close();
	
}

void Cell_Nucl_Space_Bank::readinnucl_bunch( string infile )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	
	string line;
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		readinnucl( line );
		
	}
	inf.close();
	
}

void Cell_Nucl_Space_Bank::readinles78( string infile )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	
	vector<string> filenames = parse_string(infile, '/');
	string name_s1 = filenames.back();
	vector<string> filenames2 = parse_string( name_s1, '.');
	string cell = filenames2[0];
	cout<<"cell "<<cell<<endl;
	
	string line;
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
			
		vector<string > parsed_items = parse_string(line);
		if ( parsed_items.size() < 3 )
		{
			cout<<"error1 unexpected bed line "<<line<<endl;
			exit(1);
		}
		string chr = parsed_items[0];
		if ( chr.substr(0,3) != "chr" )
		{
			cout<<"error2 unexpected bed line "<<line<<endl;
			exit(1);
		}
		
		int start = atoi(parsed_items[1].c_str());
		int end = atoi(parsed_items[2].c_str());
		int mid = start + (end-start)/2;
		cell_les78_pos[cell][chr].insert( mid );
		
	}
	inf.close();

}

void Cell_Nucl_Space_Bank::readinles78_bunch( string infile )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	
	string line;
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		readinles78( line );
		
	}
	inf.close();
	
}

void Cell_Nucl_Space_Bank::cal_space_wholechr( map<string, vector<int > > &cell_distances, int lower, int upper )
{

	for ( map<string, map<string, set<int > > >::iterator ci = cell_nucl_pos.begin(); 
				ci != cell_nucl_pos.end(); ++ci )
	{
		string cell = ci->first;
		vector<int > cspace;
	//	cout<<region.size()<<endl;
		cout<<"deal with cell "<<cell<<endl;
		for ( map<string, set<int > >::iterator sci = ci->second.begin(); sci != ci->second.end(); ++sci )
		{
			string chr = sci->first;
			vector<int > vspace;
			get_valid_space( sci->second, lower, upper, vspace );
			
			cspace.insert( cspace.end(), vspace.begin(), vspace.end() );
		}
		
		cell_distances.insert( make_pair(cell, cspace ) );
	}
}

void Cell_Nucl_Space_Bank::cal_space( map<string, vector<pair<int, int > > > &region, 
	map<string, vector<int > > &cell_distances, int lower, int upper )
{

	for ( map<string, map<string, set<int > > >::iterator ci = cell_nucl_pos.begin(); 
				ci != cell_nucl_pos.end(); ++ci )
	{
		string cell = ci->first;
		vector<int > cspace;
	//	cout<<region.size()<<endl;
		cout<<"deal with cell "<<cell<<endl;
		int nn = 0;
		for ( map<string, vector<pair<int, int > > >::iterator ite = region.begin(); ite != region.end(); ++ite )
		{	
		//	cout<<ite->second.size()<<endl;
			
			string chr = ite->first;
			set<pair<int, int > > sorted_region;
			for ( vector<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
			{
				sorted_region.insert( *si );
			}
			set<int > filtered_n;
			
			for ( set<int >::iterator cni = ci->second[chr].begin(); cni != ci->second[chr].end(); ++cni )
			{
				int n = *cni;
				for ( set<pair<int, int > >::iterator si = sorted_region.begin(); si != sorted_region.end(); ++si )
				{
					if ( si->second < n )
						continue;
					if ( si->first > n )
						break;
					filtered_n.insert(n);
					break;
				}
			}
		//	cout<<chr<<" "<<ite->second.size()<<" "<<filtered_n.size()<<" "<<ci->second[chr].size()<<endl;
			for ( vector<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
			{
				int left = si->first;
				int right = si->second;
				set<int > nu_in_r;
				for ( set<int >::iterator cni = filtered_n.begin(); cni != filtered_n.end(); ++cni )
				{
					if ( *cni < left )
						continue;
					if ( *cni > right )
						break;
					nu_in_r.insert( *cni );
				} 
				
				vector<int > vspace;
				get_valid_space( nu_in_r, lower, upper, vspace );
				
				cspace.insert( cspace.end(), vspace.begin(), vspace.end() );
				
				nn += (int)vspace.size();
				if ( nn >= 10000000 )
					break;
			}
			if ( nn >= 10000000 )
				break;
		}
		cell_distances.insert( make_pair(cell, cspace ) );
	}
}

void Cell_Nucl_Space_Bank::cal_space( map<string, vector<pair<int, int > > > &region, map<string, vector<int > > &cell_distances, 
	int lower, int upper, map<string, set<int > > &anchors )
{
	for ( map<string, map<string, set<int > > >::iterator ci = cell_nucl_pos.begin(); 
				ci != cell_nucl_pos.end(); ++ci )
	{
		string cell = ci->first;
		vector<int > cspace;
	//	cout<<region.size()<<endl;
		cout<<"deal with cell "<<cell<<endl;
		int nn = 0;
		for ( map<string, vector<pair<int, int > > >::iterator ite = region.begin(); ite != region.end(); ++ite )
		{	
		//	cout<<ite->second.size()<<endl;
			
			string chr = ite->first;
			set<pair<int, int > > sorted_region;
			for ( vector<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
			{
				sorted_region.insert( *si );
			}
			
			set<int > filtered_n;
			
			for ( set<int >::iterator cni = ci->second[chr].begin(); cni != ci->second[chr].end(); ++cni )
			{
				int n = *cni;
				for ( set<pair<int, int > >::iterator si = sorted_region.begin(); si != sorted_region.end(); ++si )
				{
					if ( si->second < n )
						continue;
					if ( si->first > n )
						break;
					filtered_n.insert(n);
					break;
				}
			}
			
			// Filter anchors Condition 1: within anchors
			set<int > filtered_anchor_n;
			for ( set< int >::iterator ni = filtered_n.begin(); ni != filtered_n.end(); ++ni )
			{
				bool wi = false;
				for ( set<int >::iterator ai = anchors[chr].begin(); ai != anchors[chr].end(); ++ai )
				{
					if ( *ai + 15 < *ni )
						continue;
					if ( *ai - 15 > *ni )
						break;
					wi = true;
					
					break;
				}
				if ( wi )
					filtered_anchor_n.insert( *ni );
			}
			// Filter anchors Condition 2: no deviated nucl
		/*	set<int > filtered_anchor_n_2;
			for ( set<int >::iterator ni = filtered_anchor_n.begin(); ni != filtered_anchor_n.end(); ++ni )
			{
				set<int >::iterator sni = ni;
				bool dev = false;
				while ( sni != filtered_anchor_n.begin() )
				{
					--sni;
					if ( *ni - *sni > 100 )
						break;
					if ( *ni - *sni > 15 )
					{
						dev = true;
						break;
					}
				}
				if ( dev )
					continue;
				sni = ni;
				++sni;
				while ( sni != filtered_anchor_n.end() )
				{
					if ( *sni - *ni > 100 )
						break;
					if ( *sni - *ni > 15 )
					{
						dev = true;
						break;
					}
					++sni;
				}
				if ( !dev )
				{
					filtered_anchor_n_2.insert( *ni );
				}
			}
			*/
			
			// cal distance between nucl and anchors within regions.
			for ( vector<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
			{
				int left = si->first;
				int right = si->second;
				set<int > nu_in_r;
				for ( set<int >::iterator cni = filtered_n.begin(); cni != filtered_n.end(); ++cni )
				{
					if ( *cni < left )
						continue;
					if ( *cni > right )
						break;
					nu_in_r.insert( *cni );
				} 
				set<int > anchor_in_r;
				for ( set<int >::iterator ni = filtered_anchor_n.begin(); ni != filtered_anchor_n.end(); ++ni )
				{
					if ( *ni < left )
						continue;
					if ( *ni > right )
						break;
					anchor_in_r.insert( *ni );
				}
				
				vector<int > vspace;
				get_valid_space( anchor_in_r, nu_in_r, lower, upper, vspace );
				
				cspace.insert( cspace.end(), vspace.begin(), vspace.end() );
				
				nn += (int)vspace.size();
				if ( nn >= 10000000 )
					break;
			}
			if ( nn >= 10000000 )
				break;
			
		}
		cell_distances.insert( make_pair(cell, cspace ) );
	}
}	


void Cell_Nucl_Space_Bank::cal_nucl_les78_in_cell( map<string, vector<int > > &centers,
		 int up, int down, map<string, vector<pair<int, int > > > &cell_nucl_les78 )
{
	map<string, map<string, map<int, set< int> > > > cell_nucl_index_pos;
	map<string, map<string, map<int, set< int> > > > cell_les78_index_pos;
	int win = 100000;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin(); ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_nucl_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	for ( map<string, map<string, set<int > > >::iterator ite = cell_les78_pos.begin(); ite != cell_les78_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_les78_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	
	for ( map<string, vector<int > >::iterator ite = centers.begin(); ite != centers.end(); ++ite )
	{	
		string chr = ite->first;
		cout<<chr<<endl;
		for ( vector<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int c = *si;
			int start = c - up;
			int end = c + down;
			int cindex = c / win;
			for ( map<string, map<string, map<int, set< int> > > >::iterator ci = cell_nucl_index_pos.begin();
				ci != cell_nucl_index_pos.end(); ++ci )
			{
				string cell = ci->first;
				vector<int > nucl_relpos;
				if ( ci->second[chr].find(cindex) != ci->second[chr].end() )
				{
					for ( set<int>::iterator tci = ci->second[chr][cindex].begin(); tci != ci->second[chr][cindex].end(); ++tci )
					{
						if ( *tci < start )
							continue;
						if ( *tci > end )
							break;
						int rlp = *tci - c;
						nucl_relpos.push_back(rlp);
					}
				}
				vector<int > les78_relpos;
				if ( cell_les78_index_pos[cell][chr].find(cindex) != cell_les78_index_pos[cell][chr].end() )
				{
					for ( set<int >::iterator tci = cell_les78_index_pos[cell][chr][cindex].begin(); 
						tci != cell_les78_index_pos[cell][chr][cindex].end(); ++tci )
					{
						if ( *tci < start )
							continue;
						if ( *tci > end )
							break;
						int rlp = *tci - c;
						les78_relpos.push_back(rlp);
					}
				}
				
				if ( !nucl_relpos.empty() && les78_relpos.empty() )
				{
					for ( size_t i = 0; i < nucl_relpos.size(); ++i )
					{
						for ( size_t j = 0; j < les78_relpos.size(); ++j )
						{
							cell_nucl_les78[cell].push_back( make_pair(nucl_relpos[i], les78_relpos[j]) );
						}
					}
				} else if ( !nucl_relpos.empty() )
				{
					for ( size_t i = 0; i < nucl_relpos.size(); ++i )
					{
						cell_nucl_les78[cell].push_back( make_pair(nucl_relpos[i], -10000) );
					}
				} else if ( !les78_relpos.empty() )
				{
					for ( size_t j = 0; j < les78_relpos.size(); ++j )
					{
						cell_nucl_les78[cell].push_back( make_pair(-10000, les78_relpos[j]) );
					}
				}
				
			}
			
			
		}
	}

}

void Cell_Nucl_Space_Bank::get_space_leftpeak_down( vector<int> &pos, int peakstart, int peakend, int upordown, vector<int> &spaces ) // 1: down; 0: up
{
	
	vector<int> nu_inpeak;
	for ( size_t i = 0; i < pos.size(); ++i )
	{
		if ( pos[i] > peakstart && pos[i] < peakend )
			nu_inpeak.push_back( pos[i]); 
	} 
	if ( nu_inpeak.empty() )
		return;
	
	for ( size_t i = 0; i < nu_inpeak.size(); ++i )
	{
		int ni = nu_inpeak[i];
		for ( size_t j = 0; j < pos.size(); ++j )
		{
			int nd = pos[j];
			int s = nd - ni;
			if ( abs(s) < 100 )
				continue;
			if ( upordown == 1 && s > 0 )
				spaces.push_back( s );
			else if ( upordown == 0 && s < 0 )
				spaces.push_back( s*(-1) );
		}
	}
	
}

void Cell_Nucl_Space_Bank::get_nuclpos_leftpeak_down( vector<int> &pos, int peakstart, int peakend, int upordown, 
	map<int, int> &relpos_count ) // 1: down; 0: up
{
	
	vector<int> nu_inpeak;
	for ( size_t i = 0; i < pos.size(); ++i )
	{
		if ( pos[i] > peakstart && pos[i] < peakend )
			nu_inpeak.push_back( pos[i]); 
	} 
	if ( nu_inpeak.empty() )
		return;
	
	for ( size_t i = 0; i < nu_inpeak.size(); ++i )
	{
		int ni = nu_inpeak[i];
		
	/*	if ( relpos_count.find(0) == 0 )
		{
			relpos_count[0] = 1;
		} else
			relpos_count[0] += 1;  */
			
		for ( size_t j = 0; j < pos.size(); ++j )
		{
			int nd = pos[j];
			int s = nd - ni;
			if ( abs(s) < 100 )
			{
				int ts = s;
				if ( upordown == 0 )
				{
					ts = s * (-1);
				}
				if ( relpos_count.find(ts) == relpos_count.end() )
					relpos_count[ts] = 1;
				else
					relpos_count[ts] += 1;
				
				continue;
			}
			
			if ( upordown == 1 && s > 0 )
			{	
			
				if ( relpos_count.find(s) == relpos_count.end() )
					relpos_count[s] = 1;
				else
					relpos_count[s] += 1;
			} else if ( upordown == 0 && s < 0 )
			{
				int ts = s * (-1);
				if ( relpos_count.find(ts) == relpos_count.end() )
					relpos_count[ts] = 1;
				else
					relpos_count[ts] += 1;
				
			}
		}
	}
	
}

void Cell_Nucl_Space_Bank::get_nuclpos_leftpeak_down_control( vector<int> &pos, int peakstart, int peakend, int upordown, 
	map<int, int> &relpos_count, map<int, int > &relpos_count_2 ) // 1: down; 0: up
{
	
	vector<int> nu_inpeak;
	for ( size_t i = 0; i < pos.size(); ++i )
	{
		if ( pos[i] > peakstart && pos[i] < peakend )
			nu_inpeak.push_back( pos[i]); 
	} 
	if ( nu_inpeak.empty() )
		return;
	
	set<size_t > colpos;
	for ( size_t i = 0; i < nu_inpeak.size(); ++i )
	{
		int ni = nu_inpeak[i];
		
	/*	if ( relpos_count.find(0) == 0 )
		{
			relpos_count[0] = 1;
		} else
			relpos_count[0] += 1;  */
			
		for ( size_t j = 0; j < pos.size(); ++j )
		{
			int nd = pos[j];
			int s = nd - ni;
			if ( abs(s) < 100 )
			{
				continue;
			}
			
			if ( upordown == 1 && s > 0 )
			{	
			
				if ( relpos_count.find(s) == relpos_count.end() )
					relpos_count[s] = 1;
				else
					relpos_count[s] += 1;
				
				colpos.insert( j );
			} else if ( upordown == 0 && s < 0 )
			{
				int ts = s * (-1);
				if ( relpos_count.find(ts) == relpos_count.end() )
					relpos_count[ts] = 1;
				else
					relpos_count[ts] += 1;
				colpos.insert( j );
			}
		}
	}
	
	// cal rel pos 2 for control: relative to center i.e. ref.
	int c = peakstart + (peakend-peakstart)/2;
	for ( set<size_t >::iterator ite = colpos.begin(); ite != colpos.end(); ++ite )
	{
		int s = pos[*ite] - c;
		if ( upordown == 1 )
		{
			if ( relpos_count_2.find( s ) == relpos_count_2.end() )
				relpos_count_2[s] = 1;
			else
				relpos_count_2[s] += 1;
		} else if ( upordown == 0 )
		{
			int ts = s * (-1);
			if ( relpos_count_2.find( ts ) == relpos_count_2.end() )
				relpos_count_2[ts] = 1;
			else
				relpos_count_2[ts] += 1;
		}
	}
	
}


void Cell_Nucl_Space_Bank::cal_distcellnum_1peak_othersidedownstream( map<string, vector<int > > &centers,
		int peakstart, int peakend, int down,   // 80 255
		map<string, map<int, vector<int > > > &dist_cell_num,
		map<string, map<int, vector<vector<string > > > > &dist_cell_list )
{
	map<string, map<string, map<int, set< int> > > > cell_nucl_index_pos;
	int win = 100000;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin(); ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_nucl_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	
	for ( map<string, vector<int > >::iterator ite = centers.begin(); ite != centers.end(); ++ite )
	{	
		string chr = ite->first;
		cout<<chr<<endl;
		for ( vector<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int c = *si;

			int cindex = c / win;
			
			int start = c - down;
			int end = c + down;
			int leftpeak_start = c - peakend;
			int leftpeak_end = c - peakstart;
			int rightpeak_start = c + peakstart;
			int rightpeak_end = c + peakend;
			
			int region1_cell = 0;
			int region2_cell = 0;
			int bothregion_cell = 0;
			vector<string > region1_cell_ve;
			vector<string > region2_cell_ve;
			vector<string > bothregion_cell_ve;
			for ( map<string, map<string, map<int, set< int> > > >::iterator ci = cell_nucl_index_pos.begin();
				ci != cell_nucl_index_pos.end(); ++ci )
			{
				string cell = ci->first;
				vector<int > nucl_relpos;
				if ( ci->second[chr].find(cindex) != ci->second[chr].end() )
				{
					for ( set<int>::iterator tci = ci->second[chr][cindex].begin(); tci != ci->second[chr][cindex].end(); ++tci )
					{
						if ( *tci < start )
							continue;
						if ( *tci > end )
							break;
						
						nucl_relpos.push_back(*tci);
					}
				}
				if ( nucl_relpos.empty() )
					continue;
				
				vector< int > t_spaces;
				get_space_leftpeak_down( nucl_relpos, leftpeak_start, leftpeak_end, 1, t_spaces );
				get_space_leftpeak_down( nucl_relpos, rightpeak_start, rightpeak_end, 0, t_spaces );
				
				bool region1 = false;
				bool region2 = false;
				for ( size_t i = 0; i < t_spaces.size(); ++i )
				{
					if ( t_spaces[i] >= 150 && t_spaces[i] <= 230 )
					{
						region1 = true;
					} else if ( t_spaces[i] >= 240 && t_spaces[i] <= 360 )
					{
						region2 = true;
					}
				}
				if ( region1 && region2 )
				{
					bothregion_cell += 1;
					bothregion_cell_ve.push_back(cell);
				} else if ( region1 )
				{
					region1_cell += 1;
					region1_cell_ve.push_back(cell);
				} else if ( region2 )
				{
					region2_cell += 1;
					region2_cell_ve.push_back(cell);
				}
			}
			
			vector<int > cell_n;
			cell_n.push_back( region1_cell );
			cell_n.push_back( region2_cell );
			cell_n.push_back( bothregion_cell );
			dist_cell_num[chr][c] = cell_n;
			vector<vector<string > > cell_list_ve;
			cell_list_ve.push_back(region1_cell_ve);
			cell_list_ve.push_back(region2_cell_ve);
			cell_list_ve.push_back(bothregion_cell_ve );
			dist_cell_list[chr][c] = cell_list_ve;
			
		}
		
	}
}

void Cell_Nucl_Space_Bank::cal_nuclvar_1peak_othersidedownstream( map<string, vector<int > > &centers,
		int peakstart, int peakend, int down,   // 80 255
		vector<int > &allelevar_onlyregion1,
		vector<int > &allelevar_onlyregion2,
		vector<int > &allelevar_both )
{
	map<string, map<string, map<int, set< int> > > > cell_nucl_index_pos;
	int win = 100000;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin(); ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_nucl_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	
	for ( map<string, vector<int > >::iterator ite = centers.begin(); ite != centers.end(); ++ite )
	{	
		string chr = ite->first;
		cout<<chr<<endl;
		for ( vector<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int c = *si;

			int cindex = c / win;
			
			int start = c - down;
			int end = c + down;
			int leftpeak_start = c - peakend;
			int leftpeak_end = c - peakstart;
			int rightpeak_start = c + peakstart;
			int rightpeak_end = c + peakend;
			
			int region1_cell = 0;
			int region2_cell = 0;
			int bothregion_cell = 0;
			vector<string > region1_cell_ve;
			vector<string > region2_cell_ve;
			vector<string > bothregion_cell_ve;
			for ( map<string, map<string, map<int, set< int> > > >::iterator ci = cell_nucl_index_pos.begin();
				ci != cell_nucl_index_pos.end(); ++ci )
			{
				string cell = ci->first;
				vector<int > nucl_relpos;
				if ( ci->second[chr].find(cindex) != ci->second[chr].end() )
				{
					for ( set<int>::iterator tci = ci->second[chr][cindex].begin(); tci != ci->second[chr][cindex].end(); ++tci )
					{
						if ( *tci < start )
							continue;
						if ( *tci > end )
							break;
						
						nucl_relpos.push_back(*tci);
					}
				}
				if ( nucl_relpos.empty() )
					continue;
				
				vector<int > var;
				getnuclvar( nucl_relpos, var );
				
				if ( var.empty() )
					continue;
					
				vector< int > t_spaces;
				get_space_leftpeak_down( nucl_relpos, leftpeak_start, leftpeak_end, 1, t_spaces );
				get_space_leftpeak_down( nucl_relpos, rightpeak_start, rightpeak_end, 0, t_spaces );
				
				bool region1 = false;
				bool region2 = false;
				for ( size_t i = 0; i < t_spaces.size(); ++i )
				{
					if ( t_spaces[i] >= 150 && t_spaces[i] <= 230 )
					{
						region1 = true;
					} else if ( t_spaces[i] >= 240 && t_spaces[i] <= 360 )
					{
						region2 = true;
					}
				}
				
				
				if ( region1 && region2 )
				{
					if ( allelevar_both.size() < 1000000 )
						allelevar_both.insert( allelevar_both.end(), var.begin(), var.end() );
				} else if ( region1 )
				{
					if ( allelevar_onlyregion1.size() < 1000000 )
						allelevar_onlyregion1.insert( allelevar_onlyregion1.end(), var.begin(), var.end() );
				} else if ( region2 )
				{
					if ( allelevar_onlyregion2.size() < 1000000 )
						allelevar_onlyregion2.insert( allelevar_onlyregion2.end(), var.begin(), var.end() );
				}
			}
			
			
		}
		
	}
}

void Cell_Nucl_Space_Bank::cal_distnuclnum_1peak_othersidedownstream( map<string, vector<int > > &centers,
		int peakstart, int peakend, int down,   // 80 255
		map<string, map<int, vector<int > > > &dist_nucl_num )
{
	map<string, map<string, map<int, set< int> > > > cell_nucl_index_pos;
	int win = 100000;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin(); ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_nucl_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	
	for ( map<string, vector<int > >::iterator ite = centers.begin(); ite != centers.end(); ++ite )
	{	
		string chr = ite->first;
		cout<<chr<<endl;
		for ( vector<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int c = *si;

			int cindex = c / win;
			
			int start = c - down;
			int end = c + down;
			int leftpeak_start = c - peakend;
			int leftpeak_end = c - peakstart;
			int rightpeak_start = c + peakstart;
			int rightpeak_end = c + peakend;
			
			int nuclregion1 = 0;
			int nuclregion2 = 0;
			
			for ( map<string, map<string, map<int, set< int> > > >::iterator ci = cell_nucl_index_pos.begin();
				ci != cell_nucl_index_pos.end(); ++ci )
			{
				string cell = ci->first;
				vector<int > nucl_relpos;
				if ( ci->second[chr].find(cindex) != ci->second[chr].end() )
				{
					for ( set<int>::iterator tci = ci->second[chr][cindex].begin(); tci != ci->second[chr][cindex].end(); ++tci )
					{
						if ( *tci < start )
							continue;
						if ( *tci > end )
							break;
						
						nucl_relpos.push_back(*tci);
					}
				}
				if ( nucl_relpos.empty() )
					continue;
				
				vector< int > t_spaces;
				get_space_leftpeak_down( nucl_relpos, leftpeak_start, leftpeak_end, 1, t_spaces );
				get_space_leftpeak_down( nucl_relpos, rightpeak_start, rightpeak_end, 0, t_spaces );
				
				
				for ( size_t i = 0; i < t_spaces.size(); ++i )
				{
					if ( t_spaces[i] >= 150 && t_spaces[i] <= 230 )
					{
						nuclregion1 += 1;
					} else if ( t_spaces[i] >= 240 && t_spaces[i] <= 360 )
					{
						nuclregion2 += 1;
					}
				}
				
			}
			
			vector<int > nucl_n;
			nucl_n.push_back( nuclregion1 );
			nucl_n.push_back( nuclregion2 );
			
			dist_nucl_num[chr][c] = nucl_n;
			
		}
		
	}
}


void Cell_Nucl_Space_Bank::cal_space_1peak_othersidedownstream( map<string, vector<int > > &centers,
		int peakstart, int peakend, int down, map<string, vector<int > > &cell_space,  // 80 255
		vector<vector<int > > &dist_cell_num )   // 150-230,   260-340, 3-ve: region1_cell; region2_cell; bothregion_cell
{
	map<string, map<string, map<int, set< int> > > > cell_nucl_index_pos;
	int win = 100000;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin(); ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_nucl_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	
	for ( map<string, vector<int > >::iterator ite = centers.begin(); ite != centers.end(); ++ite )
	{	
		string chr = ite->first;
		cout<<chr<<endl;
		for ( vector<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int c = *si;

			int cindex = c / win;
			
			int start = c - down;
			int end = c + down;
			int leftpeak_start = c - peakend;
			int leftpeak_end = c - peakstart;
			int rightpeak_start = c + peakstart;
			int rightpeak_end = c + peakend;
			
			int region1_cell = 0;
			int region2_cell = 0;
			int bothregion_cell = 0;
			for ( map<string, map<string, map<int, set< int> > > >::iterator ci = cell_nucl_index_pos.begin();
				ci != cell_nucl_index_pos.end(); ++ci )
			{
				string cell = ci->first;
				vector<int > nucl_relpos;
				if ( ci->second[chr].find(cindex) != ci->second[chr].end() )
				{
					for ( set<int>::iterator tci = ci->second[chr][cindex].begin(); tci != ci->second[chr][cindex].end(); ++tci )
					{
						if ( *tci < start )
							continue;
						if ( *tci > end )
							break;
						
						nucl_relpos.push_back(*tci);
					}
				}
				if ( nucl_relpos.empty() )
					continue;
				
				vector< int > t_spaces;
				get_space_leftpeak_down( nucl_relpos, leftpeak_start, leftpeak_end, 1, t_spaces );
				get_space_leftpeak_down( nucl_relpos, rightpeak_start, rightpeak_end, 0, t_spaces );
				
				cell_space[cell].insert( cell_space[cell].end(), t_spaces.begin(), t_spaces.end() );
				
				bool region1 = false;
				bool region2 = false;
				for ( size_t i = 0; i < t_spaces.size(); ++i )
				{
					if ( t_spaces[i] >= 150 && t_spaces[i] <= 230 )
					{
						region1 = true;
					} else if ( t_spaces[i] >= 240 && t_spaces[i] <= 360 )
					{
						region2 = true;
					}
				}
				if ( region1 && region2 )
				{
					bothregion_cell += 1;
				} else if ( region1 )
				{
					region1_cell += 1;
				} else if ( region2 )
				{
					region2_cell += 1;
				}
				
			}
			
			vector<int > cell_n;
			cell_n.push_back( region1_cell );
			cell_n.push_back( region2_cell );
			cell_n.push_back( bothregion_cell );
			dist_cell_num.push_back(cell_n);
			
		}
		
	}
}

void Cell_Nucl_Space_Bank::cal_space_1peak_samesidedownstream( map<string, vector<int > > &centers,
		int peakstart, int peakend, int down, map<string, vector<int > > &cell_space)  // 80 255 
{
	map<string, map<string, map<int, set< int> > > > cell_nucl_index_pos;
	int win = 100000;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin(); ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_nucl_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	
	for ( map<string, vector<int > >::iterator ite = centers.begin(); ite != centers.end(); ++ite )
	{	
		string chr = ite->first;
		cout<<chr<<endl;
		for ( vector<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int c = *si;

			int cindex = c / win;
			
			int start = c - down;
			int end = c + down;
			int leftpeak_start = c - peakend;
			int leftpeak_end = c - peakstart;
			int rightpeak_start = c + peakstart;
			int rightpeak_end = c + peakend;
			
			for ( map<string, map<string, map<int, set< int> > > >::iterator ci = cell_nucl_index_pos.begin();
				ci != cell_nucl_index_pos.end(); ++ci )
			{
				string cell = ci->first;
				vector<int > nucl_relpos;
				if ( ci->second[chr].find(cindex) != ci->second[chr].end() )
				{
					for ( set<int>::iterator tci = ci->second[chr][cindex].begin(); tci != ci->second[chr][cindex].end(); ++tci )
					{
						if ( *tci < start )
							continue;
						if ( *tci > end )
							break;
						
						nucl_relpos.push_back(*tci);
					}
				}
				if ( nucl_relpos.empty() )
					continue;
				
				vector< int > t_spaces;
				
				get_space_leftpeak_down( nucl_relpos, rightpeak_start, rightpeak_end, 1, t_spaces );
				get_space_leftpeak_down( nucl_relpos, leftpeak_start, leftpeak_end, 0, t_spaces );
				
				cell_space[cell].insert( cell_space[cell].end(), t_spaces.begin(), t_spaces.end() );
			}
		}
	}
				
}

void Cell_Nucl_Space_Bank::cal_space_1N_downstream( map<string, map<pair<int, int >, char  > > &N1region,
		int down, map<string, vector<int > > &cell_space )
{
	map<string, map<string, map<int, set< int> > > > cell_nucl_index_pos;
	int win = 100000;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin(); ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_nucl_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	
	for ( map<string, map<pair<int, int >,  char > >::iterator ite = N1region.begin(); ite != N1region.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, char >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			
			int start = si->first.first;
			int end = si->first.second;
			int c = start + (end-start)/2;
			
			char strand = si->second;
			int cindex = c /win;
			
			int end_ext = end;
			int start_ext = start;
			if ( strand == '+' )
			{
				end_ext = end + down;
			} else
			{
				start_ext = start - down;
			}
			
			for ( map<string, map<string, map<int, set< int> > > >::iterator ci = cell_nucl_index_pos.begin();
				ci != cell_nucl_index_pos.end(); ++ci )
			{
				string cell = ci->first;
				vector<int > nucl_relpos;
				if ( ci->second[chr].find(cindex) != ci->second[chr].end() )
				{
					for ( set<int>::iterator tci = ci->second[chr][cindex].begin(); tci != ci->second[chr][cindex].end(); ++tci )
					{
						if ( *tci < start_ext )
							continue;
						if ( *tci > end_ext )
							break;
						
						nucl_relpos.push_back(*tci);
					}
				}
				if ( nucl_relpos.empty() )
					continue;
					
				vector< int > t_spaces;
				
				if ( strand == '+')
					get_space_leftpeak_down( nucl_relpos, start, end, 1, t_spaces );
				else
					get_space_leftpeak_down( nucl_relpos, start, end, 0, t_spaces );
					
				cell_space[cell].insert( cell_space[cell].end(), t_spaces.begin(), t_spaces.end() );
				
			}
			
		}
	}
	
}

void Cell_Nucl_Space_Bank::cal_relpos_count_1N_downstream( map<string, map<pair<int, int >, pair<char, string>  > > &N1region,
		int down, map<string, map<int, int > > &gene_relpos_count )
{
	map<string, map<string, map<int, set< int> > > > cell_nucl_index_pos;
	int win = 100000;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin(); ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_nucl_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	
	for ( map<string, map<pair<int, int >, pair< char, string > > >::iterator ite = N1region.begin(); ite != N1region.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, pair<char, string> >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			
			int start = si->first.first;
			int end = si->first.second;
			int c = start + (end-start)/2;
			
			char strand = si->second.first;
			string gene = si->second.second;
			int cindex = c /win;
			
			int end_ext = end;
			int start_ext = start;
			if ( strand == '+' )
			{
				end_ext = end + down;
			} else
			{
				start_ext = start - down;
			}
			
			map<int, int > relpos_count;
			
			
			
			for ( map<string, map<string, map<int, set< int> > > >::iterator ci = cell_nucl_index_pos.begin();
				ci != cell_nucl_index_pos.end(); ++ci )
			{
				string cell = ci->first;
				vector<int > nucl_relpos;
				if ( ci->second[chr].find(cindex) != ci->second[chr].end() )
				{
					for ( set<int>::iterator tci = ci->second[chr][cindex].begin(); tci != ci->second[chr][cindex].end(); ++tci )
					{
						if ( *tci < start_ext )
							continue;
						if ( *tci > end_ext )
							break;
						
						nucl_relpos.push_back(*tci);
					}
				}
				if ( nucl_relpos.empty() )
					continue;
					
				
				
				if ( strand == '+')
					get_nuclpos_leftpeak_down( nucl_relpos, start, end, 1, relpos_count );
				else
					get_nuclpos_leftpeak_down( nucl_relpos, start, end, 0, relpos_count );
					
				
			}
			
			gene_relpos_count[gene] = relpos_count;
			
		}
	}
	
}

void Cell_Nucl_Space_Bank::cal_relpos_count_1N_downstream_control( map<string, map<pair<int, int >, pair<char, string>  > > &N1region,
		int down, map<string, map<int, int > > &gene_relpos_count, map<string, map<int, int > > &gene_relpos_count_2 )
{
	map<string, map<string, map<int, set< int> > > > cell_nucl_index_pos;
	int win = 100000;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin(); ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_nucl_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	
	for ( map<string, map<pair<int, int >, pair< char, string > > >::iterator ite = N1region.begin(); ite != N1region.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, pair<char, string> >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			
			int start = si->first.first;
			int end = si->first.second;
			int c = start + (end-start)/2;
			
			char strand = si->second.first;
			string gene = si->second.second;
			int cindex = c /win;
			
			int end_ext = end;
			int start_ext = start;
			if ( strand == '+' )
			{
				end_ext = end + down;
			} else
			{
				start_ext = start - down;
			}
			
			map<int, int > relpos_count;
			
			map<int, int > relpos_count_2;
			
			for ( map<string, map<string, map<int, set< int> > > >::iterator ci = cell_nucl_index_pos.begin();
				ci != cell_nucl_index_pos.end(); ++ci )
			{
				string cell = ci->first;
				vector<int > nucl_relpos;
				if ( ci->second[chr].find(cindex) != ci->second[chr].end() )
				{
					for ( set<int>::iterator tci = ci->second[chr][cindex].begin(); tci != ci->second[chr][cindex].end(); ++tci )
					{
						if ( *tci < start_ext )
							continue;
						if ( *tci > end_ext )
							break;
						
						nucl_relpos.push_back(*tci);
					}
				}
				if ( nucl_relpos.empty() )
					continue;
					
				
				
				if ( strand == '+')
					get_nuclpos_leftpeak_down_control( nucl_relpos, start, end, 1, relpos_count, relpos_count_2 );
				else
					get_nuclpos_leftpeak_down_control( nucl_relpos, start, end, 0, relpos_count, relpos_count_2 );
					
				
			}
			
			gene_relpos_count[gene] = relpos_count;
			gene_relpos_count_2[gene] = relpos_count_2;
			
		}
	}
	
}

void Cell_Nucl_Space_Bank::cal_truepos_count_1N_downstream( map<string, map<int, pair<char, string>  > > &tssregion, 
		int up, int down, map<string, map<int, int > > &gene_truepos_count )
{
	map<string, map<string, map<int, set< int> > > > cell_nucl_index_pos;
	int win = 100000;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin(); ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index1 = (*ti-2000) / win;
				cell_nucl_index_pos[cell][chr][index1].insert(*ti);
				int index2 = (*ti+2000) / win;
				if ( index2 > index1 )
					cell_nucl_index_pos[cell][chr][index2].insert(*ti);
			}
		}
	} 
	
	for ( map<string, map<int, pair< char, string > > >::iterator ite = tssregion.begin(); ite != tssregion.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<int, pair<char, string> >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int tss = si->first;
			
			char strand = si->second.first;
			string gene = si->second.second;
			int cindex = tss /win;
			
			int start = tss - up;
			if ( strand == '-' )
				start = tss - down;
			int end = tss + down;
			if ( strand == '-')
				end = tss + up; 
			
			map<int, int > truepos_count;
			
			for ( map<string, map<string, map<int, set< int> > > >::iterator ci = cell_nucl_index_pos.begin();
				ci != cell_nucl_index_pos.end(); ++ci )
			{
				string cell = ci->first;
				
				if ( ci->second[chr].find(cindex) != ci->second[chr].end() )
				{
					for ( set<int>::iterator tci = ci->second[chr][cindex].begin(); tci != ci->second[chr][cindex].end(); ++tci )
					{
						if ( *tci < start )
							continue;
						if ( *tci > end )
							break;
						
						int pos = *tci - tss;
						if ( strand == '-' )
						{
							pos = tss - *tci;
						}
						if ( truepos_count.find( pos ) == truepos_count.end() )
							truepos_count[pos] = 1;
						else
							truepos_count[pos] += 1;
					//	if ( gene == "NM_001001130" )
					//		cout<<"tmp "<<pos<<endl;
					}
				}
			}
			
			gene_truepos_count[gene] = truepos_count;
		}
	}
				
}

void Cell_Nucl_Space_Bank::cal_space_1peak_othersidedownstream_les78filter( map<string, vector<int > > &centers,
		int peakstart, int peakend, int down, map<string, vector<int > > &cell_space )   // 80 255
{
	map<string, map<string, map<int, set< int> > > > cell_nucl_index_pos;
	map<string, map<string, map<int, set< int> > > > cell_les78_index_pos;
	int win = 100000;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin(); ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_nucl_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	for ( map<string, map<string, set<int > > >::iterator ite = cell_les78_pos.begin(); ite != cell_les78_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_les78_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	
	for ( map<string, vector<int > >::iterator ite = centers.begin(); ite != centers.end(); ++ite )
	{	
		string chr = ite->first;
		cout<<chr<<endl;
		for ( vector<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int c = *si;

			int cindex = c / win;
			
			int start = c - down;
			int end = c + down;
			int leftpeak_start = c - peakend;
			int leftpeak_end = c - peakstart;
			int rightpeak_start = c + peakstart;
			int rightpeak_end = c + peakend;
			
			
			for ( map<string, map<string, map<int, set< int> > > >::iterator ci = cell_nucl_index_pos.begin();
				ci != cell_nucl_index_pos.end(); ++ci )
			{
				string cell = ci->first;
				
				// les78filter 
				bool les78occ = false;
				if ( cell_les78_index_pos[cell][chr].find(cindex ) != cell_les78_index_pos[cell][chr].end() )
				{
					for ( set<int>::iterator tci = cell_les78_index_pos[cell][chr][cindex].begin();
						tci != cell_les78_index_pos[cell][chr][cindex].end(); ++tci )
					{
						if ( *tci < c - 80 )
							continue;
						if ( *tci > c + 80 )
							break;
						les78occ = true;
					}
				}
				if ( !les78occ )
					continue;
				// les78filter end
				
				vector<int > nucl_relpos;
				if ( ci->second[chr].find(cindex) != ci->second[chr].end() )
				{
					for ( set<int>::iterator tci = ci->second[chr][cindex].begin(); tci != ci->second[chr][cindex].end(); ++tci )
					{
						if ( *tci < start )
							continue;
						if ( *tci > end )
							break;
						
						nucl_relpos.push_back(*tci);
					}
				}
				if ( nucl_relpos.empty() )
					continue;
				
				vector< int > t_spaces;
				get_space_leftpeak_down( nucl_relpos, leftpeak_start, leftpeak_end, 1, t_spaces );
				get_space_leftpeak_down( nucl_relpos, rightpeak_start, rightpeak_end, 0, t_spaces );
				
				cell_space[cell].insert( cell_space[cell].end(), t_spaces.begin(), t_spaces.end() );
			}
			
		}
		
	}
}

void Cell_Nucl_Space_Bank::get_space_peaknucl_middleles78( vector<int > &posnucl, vector<int > &posles78,
		int lpeakstart, int lpeakend, int rpeakstart, int rpeakend, int middle_lbound, int middle_rbound, vector<int > &spaces )
{
	vector<int > nu_in_peak;
	for ( size_t i = 0; i < posnucl.size(); ++i )
	{
		if ( (posnucl[i] > lpeakstart && posnucl[i] < lpeakend) || (posnucl[i] > rpeakstart && posnucl[i] < rpeakend) )
			nu_in_peak.push_back( posnucl[i]); 
	} 
	
	vector<int > les78_inmiddle;
	for ( size_t i = 0; i < posles78.size(); ++i )
	{
		if ( posles78[i] > middle_lbound && posles78[i] < middle_rbound )
			les78_inmiddle.push_back( posles78[i] );
	}
//	cout<<nu_in_peak.size()<<" "<<les78_inmiddle.size()<<" "<<lpeakstart<<" "<<lpeakend<<" "<<rpeakstart<<" "<<rpeakend<<" "<<middle_lbound<<" "<<middle_rbound<<endl;
	for ( size_t i = 0; i < nu_in_peak.size(); ++i )
	{
		int ni = nu_in_peak[i];
		for ( size_t j = 0; j < les78_inmiddle.size(); ++j )
		{
			int les78i = les78_inmiddle[j];
			int s = abs(les78i - ni);
			spaces.push_back( s );
		}
	}
}

void Cell_Nucl_Space_Bank::cal_space_1peak_nucl_middle_les78( map<string, vector<int > > &centers, 
		int peakstart, int peakend, int middle_bound, map<string, vector<int > > &cell_space,
		vector<vector<int > > &dist_cell_num )
{
	map<string, map<string, map<int, set< int> > > > cell_nucl_index_pos;
	map<string, map<string, map<int, set< int> > > > cell_les78_index_pos;
	int win = 100000;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin(); ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_nucl_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	for ( map<string, map<string, set<int > > >::iterator ite = cell_les78_pos.begin(); ite != cell_les78_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si  )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / win;
				cell_les78_index_pos[cell][chr][index].insert(*ti);
			}
		}
	} 
	
	for ( map<string, vector<int > >::iterator ite = centers.begin(); ite != centers.end(); ++ite )
	{	
		string chr = ite->first;
		cout<<chr<<endl;
		for ( vector<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int c = *si;

			int cindex = c / win;
			
			int start = c - peakend;
			int end = c + peakend;
			int leftpeak_start = c - peakend;
			int leftpeak_end = c - peakstart;
			int rightpeak_start = c + peakstart;
			int rightpeak_end = c + peakend;
			int leftmiddlebound = c - middle_bound;
			int rightmiddlebound = c + middle_bound;
			
			int pos_cell = 0;
			int neg_cell = 0;
			
			
			for ( map<string, map<string, map<int, set< int> > > >::iterator ci = cell_nucl_index_pos.begin();
				ci != cell_nucl_index_pos.end(); ++ci )
			{
				string cell = ci->first;
				vector<int > nucl_relpos;
				if ( ci->second[chr].find(cindex) != ci->second[chr].end() )
				{
					for ( set<int>::iterator tci = ci->second[chr][cindex].begin(); tci != ci->second[chr][cindex].end(); ++tci )
					{
						if ( *tci < start )
							continue;
						if ( *tci > end )
							break;
						
						nucl_relpos.push_back(*tci);
					}
				}
				if ( nucl_relpos.empty() )
					continue;
				
				vector<int > les78_relpos;
				if ( cell_les78_index_pos[cell][chr].find(cindex) != cell_les78_index_pos[cell][chr].end() )
				{
					for ( set<int >::iterator tci = cell_les78_index_pos[cell][chr][cindex].begin(); 
						tci != cell_les78_index_pos[cell][chr][cindex].end(); ++tci )
					{
						if ( *tci < leftmiddlebound )
							continue;
						if ( *tci > rightmiddlebound )
							break;
						
						les78_relpos.push_back(*tci);
					}
				}
				
				if ( les78_relpos.empty() )
				{
					neg_cell += 1;
					continue;
				}
				
				vector<int > t_spaces;
				get_space_peaknucl_middleles78( nucl_relpos, les78_relpos, leftpeak_start, leftpeak_end, 
					rightpeak_start, rightpeak_end, leftmiddlebound, rightmiddlebound, t_spaces );
				cell_space[cell].insert( cell_space[cell].end(), t_spaces.begin(), t_spaces.end() );
			//	cout<<"A "<<t_spaces.size()<<" "<<nucl_relpos.size()<<" "<<les78_relpos.size()<<endl;
			
				if ( !t_spaces.empty() )
					pos_cell += 1;
				else 
					neg_cell += 1;
				
			}
			vector<int > cell_n;
			cell_n.push_back( pos_cell );
			cell_n.push_back( neg_cell );
			
			dist_cell_num.push_back(cell_n);
		}
	}
	
}

void get_valid_space( set<int >& nuclpos, int lower, int upper, vector<int > &space )
{
	for ( set< int >::iterator ite = nuclpos.begin();  ite != nuclpos.end(); ++ite )
	{
		int n1 = *ite;
		
		set< int >::iterator si = ite;
		
		++si;
		for ( ; si != nuclpos.end(); ++si )
		{
			int n2 = *si;
			int d = n2 - n1;
			if ( d >= lower && d <= upper )
			{
				space.push_back( d );
			} 
			if ( d > upper )
				break;
		}
	}
}

void get_valid_space( set<int >& anchors, set<int >& nuclpos, int lower, int upper, vector<int > &space )
{
	for ( set<int >::iterator ite = anchors.begin(); ite != anchors.end(); ++ite )
	{
		int n1 = *ite;
		
		for ( set<int >::iterator si = nuclpos.begin(); si != nuclpos.end(); ++si )
		{
			int n2 = *si;
			int d = abs(n2 - n1);
			if ( d >= lower && d <= upper )
			{
				space.push_back( d );
			} 
			if ( n2 - n1 > upper )
				break;
		} 
	}
}

void pool_space(map<string, vector<int > > &cell_distances, vector<int > &pooled_distances )
{
	for ( map<string, vector<int > >::iterator ite = cell_distances.begin(); ite != cell_distances.end(); ++ite )
	{
		pooled_distances.insert( pooled_distances.end(), ite->second.begin(), ite->second.end() );
		
	}
}

void cal_ratio( map<string, vector<int > > &cell_distances, map<string, double > &cell_ratio, int lower, int upper )
{
	for ( map<string, vector<int > >::iterator ite = cell_distances.begin(); ite != cell_distances.end(); ++ite )
	{	
		int total = (int)ite->second.size();
		int in = 0;
		for ( vector<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			if ( *si >= lower && *si <= upper )
				in += 1;
		}
		double r = in*1.0/total;
	//	cout<<ite->first<<" "<<total<<" "<<in<<" "<<r<<endl;
		cell_ratio.insert( make_pair(ite->first, r ) );
	}
}

void getshift_amongcells( map<string, map<string, set<int > > > &cell_nucl_pos, 
	map<string, map<int, vector< int> > > &chr_cell_shift )
{
	map<string, map<int, set<string > > > chr_pos_cells;
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin();
		ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			string chr = si->first;
			for ( set<int >::iterator ci = si->second.begin(); ci != si->second.end(); ++ci )
			{
				int pos = *ci;
				chr_pos_cells[chr][pos].insert( cell );
			}
		}
	}
	
	for ( map<string, map<int, set<string > > >::iterator ite = chr_pos_cells.begin(); 
		ite != chr_pos_cells.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<int, set<string > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int pos = si->first;
			map<int, set<string > >::iterator nsi = si;
			++nsi;
			for ( ; nsi != ite->second.end(); ++nsi )
			{
				int l = nsi->first - pos;
				if ( l <= 83 && l > 2 )
				{
					for ( set<string >::iterator fi = si->second.begin(); fi != si->second.end(); ++fi )
					{
						for ( set<string >::iterator sfi = nsi->second.begin(); sfi != nsi->second.end(); ++sfi )
						{	
							if ( *sfi != *fi )
							{
								chr_cell_shift[chr][pos].push_back( l );
								chr_cell_shift[chr][nsi->first].push_back( l );
							}
						}
					}
				} else
					break;
			}
		}
	}
}

void getshift_withincell( map<string, map<string, set<int > > > &cell_nucl_pos, 
	map<string, map<int, vector<int > > > &chr_allele_shift )
{
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin();
		ite != cell_nucl_pos.end(); ++ite )
	{
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				set<int >::iterator nti = ti;
				++nti;
				for ( ; nti != si->second.end(); ++nti )
				{
					int d = *nti - *ti;
					if ( d > 83 )
						break;
					if ( d <= 2 )
						continue;
					chr_allele_shift[chr][*ti].push_back( d );
					chr_allele_shift[chr][*nti].push_back(d );
				}  
			}
		}
	}
}

void pool_nuclpos(Cell_Nucl_Space_Bank &bank, map<string, set<int > > &nucl_pos )
{
	for ( map<string, map<string, set<int > > >::iterator ite = bank.cell_nucl_pos.begin(); 
		ite != bank.cell_nucl_pos.end(); ++ite )
	{
		for ( map<string, set<int > >::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			string chr = si->first;
			nucl_pos[chr].insert( si->second.begin(), si->second.end() );
		}
	}
}

void getlendis( map<string, set<int > > &nucl_pos, string prefix, map<string, set<int> > &phased_nucl)
{
	int L = 100000;
	map<string, map<int, set<int > > > chr_index_site;
	for ( map<string, set<int> >::iterator ite = phased_nucl.begin(); ite != phased_nucl.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int site = *si;
			int index = site / L;
			chr_index_site[chr][index].insert( site );
		}
	}
	
	string outfile = prefix+".link_len.txt";
	ofstream outf( outfile.data() );
	int num = 0;
	for ( map<string, set<int > >::iterator ite = nucl_pos.begin(); ite != nucl_pos.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int n1 = *si;
			
			int tindex = n1/ L;
			bool fd = false;
			if ( chr_index_site[chr].find(tindex) != chr_index_site[chr].end() )
			{
				for ( set<int>::iterator ci = chr_index_site[chr][tindex].begin();
					ci != chr_index_site[chr][tindex].end(); ++ci )
				{
					if ( *ci < n1 - 15 )
						continue;
					if ( *ci > n1 + 15 )
						break;
					fd = true;
					break;
				}
			}
			if ( !fd )
				continue;
				
			num += 1;
			if ( num % 10000==0)
				cout<<num<<endl;
			if ( num >= 1000000 )
				break;
			set<int>::iterator nsi = si;
			++nsi;
			for ( ; nsi != ite->second.end(); ++nsi )
			{
				int n2 = *nsi;
				if ( n2 - n1 <= 30 )
					continue;
				if ( n2 - n1 > 1000 )
					break;
				outf<<n2-n1<<endl;
				
			}
		}
		
	}
	outf.close();
}

void getlendis( map<string, set<int > > &nucl_pos, string prefix, map<string, set<pair< int, int > > > &region)
{
	int L = 100000;
	map<string, map<int, set<pair<int,int> > > > chr_index_site;
	for ( map<string, set<pair<int, int > > >::iterator ite = region.begin(); ite != region.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<pair<int, int> >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			
			int index1 = si->first / L;
			int index2 = si->second / L;
			chr_index_site[chr][index1].insert( *si );
			chr_index_site[chr][index2].insert( *si );
		}
	}
	
	string outfile = prefix+".link_len.txt";
	ofstream outf( outfile.data() );
	int num = 0;
	int nucl_num = 0;
	for ( map<string, set<int > >::iterator ite = nucl_pos.begin(); ite != nucl_pos.end(); ++ite )
	{
		string chr = ite->first;
		nucl_num += ite->second.size();
		for ( set<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int n1 = *si;
			int index = n1/L;
			
			bool fd = false;
			if ( chr_index_site[chr].find( index) != chr_index_site[chr].end() )
			{
				for ( set<pair< int, int > >::iterator ci = chr_index_site[chr][index].begin(); 
					ci != chr_index_site[chr][index].end(); ++ci )
				{
					if ( ci->second < n1 )
						continue;
					if ( ci->first > n1 )
						break;
					fd = true;
					break;
				}
			}
			if ( !fd )
				continue;
			
			num += 1;
			if ( num % 10000 == 0)
				cout<<num<<endl;
			if ( num >= 1000000 )
				break;
			
			set<int>::iterator nsi = si;
			++nsi;
			for ( ; nsi != ite->second.end(); ++nsi )
			{
				int n2 = *nsi;
				if ( n2 - n1 <= 30 )
					continue;
				if ( n2 - n1 > 1000 )
					break;
				outf<<n2-n1<<endl;
			}
		}
		
	}
	outf.close();
	
	cout<<"nucl size " <<nucl_num<<endl;
}

void getnuclvar( vector<int > &relpos, vector<int > &var )
{
	set<int> relposset;
	for ( size_t i = 0; i < relpos.size(); ++i )
	{
		relposset.insert( relpos[i] );
	}
	for ( set<int>::iterator ite = relposset.begin(); ite != relposset.end(); ++ite )
	{
		int n1 = *ite;
		set<int >::iterator si = ite;
		++si;
		for ( ; si != relposset.end(); ++si )
		{
			int n2 = *si;
			int v = n2-n1;
			if ( v >= 83 )
				break;
			if ( v > 2 )
				var.push_back( v );
		}
	}
}









