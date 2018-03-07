#include "operation.h"

void Position::addpos(int inst, int inend, string inchr)
{
	start = inst;
	end = inend;
	chr = inchr;
}

void Position::addpos(string str)
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

int Stitched_Region::getstart()
{
	if ( regions.empty() )
	{
		cout<<"error in Stitched_Region::getstart(): regions empty! not able to get start."<<endl; exit(1);
	}
	
	return regions.front().first;
}

int Stitched_Region::getend()
{
	if ( regions.empty() )
	{
		cout<<"error in Stitched_Region::getend(): regions empty! not able to get end."<<endl; exit(1);
	}
	
	return regions.back().second;
}

void Count_pool::getregionidmap_tss( int ups, int downs, int win )
{
	cut_reg_ve.clear();
	for ( size_t i = 0; i < id_ve.size(); ++i )
	{
		vector<pair<int, int > > nullve;
		cut_reg_ve.push_back(nullve);
	}
	ofstream outf("tmp2.txt");
	for ( map<string, map<int, char> >::iterator ite = tss_strand_map.begin(); ite != tss_strand_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<int, char >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			
			size_t id = tss_id_map[ite->first][si->first];
			int tss = si->first;
			int start = 0;
			int end = 0;
			if ( si->second == '+' )
			{
				start = tss - ups;
				end = tss + downs -1; 
				
			} else if ( si->second == '-' )
			{
				start = tss - downs + 1;
				end = tss + ups;
				
			} else 
			{
				cout<<"error strand in Count_pool::getreagionidmap_tss "<<si->second<<endl; exit(1);
			}
			vector<pair<int, int > > cut_reg; 
			vector<int > ct;
			int sp = start;
			int ep = end;
			if ( tss_protect_map[chr].find( tss ) != tss_protect_map[chr].end() )
			{
				if ( tss_protect_map[chr][tss].first != 0 )
				{
					if ( sp < tss_protect_map[chr][tss].first )
						sp = tss_protect_map[chr][tss].first;
				} 
				if ( tss_protect_map[chr][tss].second != 0 )
				{
					if ( ep > tss_protect_map[chr][tss].second )
						ep = tss_protect_map[chr][tss].second;
				} 
			}
			if ( sp != start || ep != end )
			{
				outf<<chr<<"\t"<<tss<<"\t"<<sp<<"\t"<<ep<<endl;
			}
			region_id_map[chr][make_pair(sp, ep)] = id;
			region_strand_map[chr][make_pair(sp, ep)] = si->second;
			while ( start < end )
			{
				int tend = start + win - 1;
				if ( tend > end )
					tend = end;
				cut_reg.push_back( make_pair(start, tend ) );
				
				ct.push_back(0);
				start = tend + 1;
			}
			cut_reg_ve[id] = cut_reg;
			counts_table[id] = ct;
		}
	}
		
}

void Count_pool::getregionidmap_body_sp2( int ups, int downs, int win, int step )
{
	cut_reg_ve.clear();
	for ( size_t i = 0; i < id_ve.size(); ++i )
	{
		vector<pair<int, int > > nullve;
		cut_reg_ve.push_back(nullve);
	}
	for ( map<string, map<pair<int, int >, char> >::iterator ite = body_strand_map.begin(); ite != body_strand_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, char >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			size_t id = body_id_map[ite->first][si->first];
			if ( si->first.second - si->first.first + 1 < step )
				continue;
			int tss = si->first.first;
			if ( si->second == '-' )
				tss = si->first.second;
			int start = 0;
			int end = 0;
			vector<pair<int, int > > cut_reg; 
			vector<int > ct;
			
			if ( si->second == '+' )
			{	
				start = tss - ups;
				end = tss-1;
				int vstart = start;
				while ( start < end )
				{
					int tend = start + win - 1;
					if ( tend > end )
						tend = end;
					cut_reg.push_back( make_pair(start, tend ) );
				
					ct.push_back(0);
					start = tend + 1;
				}
				cut_reg_ve[id] = cut_reg;
				counts_table[id] = ct;
				
				int gstart = end+1;
				int gend = si->first.second;
				int len = gend-gstart+1;
				double twin = len*1.0 / step;
				if ( twin < 1 )
				{
					cout<<win<<" "<<len<<" "<<step<<endl;
					exit(1);
				}
				
				for ( int i = 0; i < step; ++i )
				{
					double tstart = gstart+ i*twin;
					double tend = gstart + (i+1)*twin;
				
				
					cut_reg.push_back( make_pair((int)tstart, (int)tend-1 ) );
					if ( cut_reg.back().second <= cut_reg.back().first )
					{
						cout<<"error reg "<<cut_reg.back().first<<" "<<cut_reg.back().second<<endl;
					}
					ct.push_back(0);
				
				}
				
				start = si->first.second+1;
				end = start + downs - 1;
				while ( start < end )
				{
					int tend = start + win - 1;
					if ( tend > end )
						tend = end;
					cut_reg.push_back( make_pair(start, tend ) );
				
					ct.push_back(0);
					start = tend + 1;
				}
				cut_reg_ve[id] = cut_reg;
				counts_table[id] = ct;
				region_id_map[chr][make_pair(vstart, end)] = id;
				region_strand_map[chr][make_pair(vstart, end)] = si->second;
			} else
			{
				start = si->first.first - downs;
				end = si->first.first-1;
				int vstart = start;
				while ( start < end )
				{
					int tend = start + win - 1;
					if ( tend > end )
						tend = end;
					cut_reg.push_back( make_pair(start, tend ) );
				
					ct.push_back(0);
					start = tend + 1;
				}
				cut_reg_ve[id] = cut_reg;
				counts_table[id] = ct;
				
				int gstart = end+1;
				int gend = si->first.second;
				int len = gend-gstart+1;
				double twin = len*1.0 / step;
				if ( twin < 1 )
				{
					cout<<win<<" "<<len<<" "<<step<<endl;
					exit(1);
				}
				
				for ( int i = 0; i < step; ++i )
				{
					double tstart = gstart+ i*twin;
					double tend = gstart + (i+1)*twin;
				
				
					cut_reg.push_back( make_pair((int)tstart, (int)tend-1 ) );
					if ( cut_reg.back().second <= cut_reg.back().first )
					{
						cout<<"error reg "<<cut_reg.back().first<<" "<<cut_reg.back().second<<endl;
					}
					ct.push_back(0);
				
				}
				
				start = si->first.second+1;
				end = start + ups - 1;
				while ( start < end )
				{
					int tend = start + win - 1;
					if ( tend > end )
						tend = end;
					cut_reg.push_back( make_pair(start, tend ) );
				
					ct.push_back(0);
					start = tend + 1;
				}
				cut_reg_ve[id] = cut_reg;
				counts_table[id] = ct;
				region_id_map[chr][make_pair(vstart, end)] = id;
				region_strand_map[chr][make_pair(vstart, end)] = si->second;
				
			}
			
			// temp
	/*		ofstream outf("tmp.log");
			outf<<chr<<"\t"<<si->first.first<<"\t"<<si->first.second<<"\t"<<si->second<<"\t"<<tss<<endl;
			for ( size_t t = 0; t != cut_reg.size(); ++t )
			{
				outf<<cut_reg[t].first<<"\t"<<cut_reg[t].second<<"\t"<<ct[t]<<endl;
			}
			outf.close();
	*/	
		}
	}
}

void Count_pool::getregionidmap_body_sp( int ups, int downs, int win, int step )
{
	cut_reg_ve.clear();
	for ( size_t i = 0; i < id_ve.size(); ++i )
	{
		vector<pair<int, int > > nullve;
		cut_reg_ve.push_back(nullve);
	}
	for ( map<string, map<pair<int, int >, char> >::iterator ite = body_strand_map.begin(); ite != body_strand_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, char >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			size_t id = body_id_map[ite->first][si->first];
			
			int tss = si->first.first;
			if ( si->second == '-' )
				tss = si->first.second;
			int start = 0;
			int end = 0;
			vector<pair<int, int > > cut_reg; 
			vector<int > ct;
			int vstart = 0;
			if ( si->second == '+' )
			{
			
				start = tss - ups;
				end = tss + downs -1; 
				vstart = start;
			
				while ( start < end )
				{
					int tend = start + win - 1;
					if ( tend > end )
						tend = end;
					cut_reg.push_back( make_pair(start, tend ) );
				
					ct.push_back(0);
					start = tend + 1;
				}
				cut_reg_ve[id] = cut_reg;
				counts_table[id] = ct;
			
			
				int gstart = end+1;
				int gend = si->first.second;
				int len = gend-gstart+1;
				double twin = len*1.0 / step;
				if ( twin < 1 )
				{
					cout<<twin<<" "<<len<<" "<<step<<endl;
					exit(1);
				}
			/*	start = 0;
				end = 0;
				if ( si->second == '+' )
				{
					start = gstart - (ups*win);
					end = gend + (downs*win); 
				
				} else if ( si->second == '-' )
				{
					start = gstart - downs*win;
					end = gend + ups*win;
				
				} else 
				{
					cout<<"error strand in Count_pool::getreagionidmap_body "<<si->second<<endl; exit(1);
				} */
			//	double ewin = (end - start) / (step+ups+downs);
			//	cout<<gstart<<","<<gend<<","<<win<<"\t"<<start<<","<<end<<","<<ewin<<endl;
				
				for ( int i = 0; i < step; ++i )
				{
					double tstart = gstart+ i*twin;
					double tend = gstart + (i+1)*twin;
				
				
					cut_reg.push_back( make_pair((int)tstart, (int)tend-1 ) );
					if ( cut_reg.back().second <= cut_reg.back().first )
					{
						cout<<"error reg "<<cut_reg.back().first<<" "<<cut_reg.back().second<<endl;
					}
					ct.push_back(0);
				
				
				}
				cut_reg_ve[id] = cut_reg;
				counts_table[id] = ct;
				region_id_map[chr][make_pair(vstart, gend)] = id;
				region_strand_map[chr][make_pair(vstart, gend)] = si->second;
			} else
			{
				start = tss - downs + 1;
				end = tss + ups;
				int gstart = si->first.first;
				int gend = start-1;
				int len = gend-gstart+1;
				double twin = len*1.0 / step;
				if ( twin < 1 )
				{
					cout<<twin<<" "<<len<<" "<<step<<endl;
					
					exit(1);
				}
				for ( int i = 0; i < step; ++i )
				{
					double tstart = gstart+ i*twin;
					double tend = gstart + (i+1)*twin;
				
				
					cut_reg.push_back( make_pair((int)tstart, (int)tend-1 ) );
					if ( cut_reg.back().second <= cut_reg.back().first )
					{
						cout<<"error reg "<<cut_reg.back().first<<" "<<cut_reg.back().second<<endl;
					}
					ct.push_back(0);
				
				
				}
				
				while ( start < end )
				{
					int tend = start + win - 1;
					if ( tend > end )
						tend = end;
					cut_reg.push_back( make_pair(start, tend ) );
				
					ct.push_back(0);
					start = tend + 1;
				}
				cut_reg_ve[id] = cut_reg;
				counts_table[id] = ct;
				region_id_map[chr][make_pair(gstart, end)] = id;
				region_strand_map[chr][make_pair(gstart, end)] = si->second;
				
			}	
			
		/*	if ( id_ve[id] == "ENSMUST00000096232" )
			{
				for ( size_t i = 0; i < cut_reg_ve[id].size(); ++i )
				{
					cout<<cut_reg_ve[id][i].first<<"\t"<<cut_reg_ve[id][i].second<<"\t"<<counts_table[id][i]<<"\t"<<id<<"\t"<<si->second<<"\t"<<chr<<endl;
				}
			} */
		}
	}
//	exit(1);
/*	for ( size_t i = 0; i < counts_table.size(); ++i )
	{
		if ( counts_table[i].size() != step+ups+downs )
		{
			cout<<"error col "<<counts_table[i].size()<< ' '<<i<<endl; 
		}
	}*/
}

void Count_pool::getregionidmap_body( int ups, int downs, int step )
{
	cut_reg_ve.clear();
	for ( size_t i = 0; i < id_ve.size(); ++i )
	{
		vector<pair<int, int > > nullve;
		cut_reg_ve.push_back(nullve);
	}
	for ( map<string, map<pair<int, int >, char> >::iterator ite = body_strand_map.begin(); ite != body_strand_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, char >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			size_t id = body_id_map[ite->first][si->first];
			int gstart = si->first.first;
			int gend = si->first.second;
			int len = gend-gstart+1;
			double win = len*1.0 / step;
			if ( win < 1 )
			{
				cout<<win<<" "<<len<<" "<<step<<" "<<chr<<" "<<gstart<<" "<<gend<<endl;
				exit(1);
			}
			double start = 0;
			double end = 0;
			if ( si->second == '+' )
			{
				start = gstart - (ups*win);
				end = gend + (downs*win); 
				
			} else if ( si->second == '-' )
			{
				start = gstart - downs*win;
				end = gend + ups*win;
				
			} else 
			{
				cout<<"error strand in Count_pool::getreagionidmap_body "<<si->second<<endl; exit(1);
			}
			double ewin = (end - start) / (step+ups+downs);
		//	cout<<gstart<<","<<gend<<","<<win<<"\t"<<start<<","<<end<<","<<ewin<<endl;
			vector<pair<int, int > > cut_reg; 
			vector<int > ct;
			region_id_map[chr][make_pair((int)start, (int)end)] = id;
			region_strand_map[chr][make_pair((int)start, (int)end)] = si->second;
			for ( int i = 0; i < step+ups+downs; ++i )
			{
				double tstart = start+ i*ewin;
				double tend = start + (i+1)*ewin;
				
				
				cut_reg.push_back( make_pair((int)tstart, (int)tend-1 ) );
				if ( cut_reg.back().second <= cut_reg.back().first )
				{
					cout<<"error reg "<<cut_reg.back().first<<" "<<cut_reg.back().second<<endl;
				}
				ct.push_back(0);
				
				
			}
			cut_reg_ve[id] = cut_reg;
			counts_table[id] = ct;
			if ( (int)ct.size() != step+ups+downs )
			{
				cout<<"error ctsize "<<ct.size()<<endl; exit(1);
			}
			
		/*	if ( id == 0 )
			{
				for ( size_t k = 0; k < cut_reg.size(); ++k )
				{
					cout<<cut_reg[k].first<<"\t"<<cut_reg[k].second<<endl; 
				}
				exit(1);
			}*/
			
			
		//	cout<<ct.size()<<endl;
		}
	}
}



void Count_pool::normrpkm( int total )
{
	rpkm_table.clear();
	for ( size_t i = 0; i < cut_reg_ve.size(); ++i )
	{
		vector<double > sv;
		for ( size_t j = 0; j < cut_reg_ve[i].size(); ++j )
		{
			int w = cut_reg_ve[i][j].second - cut_reg_ve[i][j].first + 1;
		//	double nc = ( counts_table[i][j] * 1.0 * 1000 * 1000000 ) / ( w * total );
			double nc = counts_table[i][j] * 1.0 * ( 1000.0 / w ) * ( 1000000.0 / total );
			sv.push_back( nc );
			
		}
		rpkm_table.push_back( sv );
	}
}

void Count_pool::smooth_ave_rpkm( int stretch )
{
	if ( stretch > 0 )
	{
		vector<vector<double > > smoothed_rpkm;
		for ( size_t i = 0; i < rpkm_table.size(); ++i )
		{
			vector<double > sv;
			for ( size_t j = 0; j < rpkm_table[i].size(); ++j )
			{
				int start = max( ((int)j - stretch ), 0 );
				int end = min( ((int)j + stretch ), (int)(rpkm_table[i].size())-1 );
				int n = end-start+1;
				double sm_v = 0;
				for ( int k = start; k <= end; ++k )
				{
					sm_v += rpkm_table[i][k];
				}
				sm_v = sm_v / n;
				sv.push_back( sm_v );
			}
			smoothed_rpkm.push_back( sv );
		}
		rpkm_table = smoothed_rpkm;
	}
}

void Count_pool::reverseminusstrand()
{
	for ( map<string, map<pair<int, int >, char > >::iterator ite = region_strand_map.begin(); ite != region_strand_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, char >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			if ( si->second == '-')
			{
				size_t id = region_id_map[chr][si->first];
				reverse( counts_table[id].begin(), counts_table[id].end() );
				reverse( rpkm_table[id].begin(), rpkm_table[id].end() );
			}
		}
	}
}

void Count_pool::gettssprotectmap()
{
	
	for ( map<string, map<int, size_t> >::iterator ite = tss_id_map.begin(); ite != tss_id_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<int, size_t >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int tss = si->first;
			map<int, size_t >::iterator ls = si;
			int lp = 0;
			while ( ls != ite->second.begin() )
			{
				--ls;
				int l = ls->first+500;
				if ( l < tss - 500 )
				{
					lp = l;
					break;
				}
			}
			map<int, size_t >::iterator rs = si;
			int rp = 0;
			++rs;
			while ( rs != ite->second.end() )
			{
				int r = rs->first - 500;
				if ( r > tss + 500 )
				{
					rp = r;
					break;
				}
				++rs;
			}
			tss_protect_map[chr][si->first] = make_pair( lp, rp );
			
		}
	}
	
}

vector<string > parse_string( string & instr, char spl )
{
	vector<string > strve;
	string s = "";
	for ( size_t i = 0; i < instr.size(); ++i )
	{
		if ( instr[i] == spl )
		{
			if ( !s.empty() )
			{
				strve.push_back(s);
				s = "";
			}
		} else
		{
			s += instr[i];
		}
	}
	if ( !s.empty() )
	{
		strve.push_back(s);
	}
	return strve;
}

vector<string > parse_string( string & instr)
{
	vector<string > strve;
	string s = "";
	for ( size_t i = 0; i < instr.size(); ++i )
	{
		if ( instr[i] == '\t' || instr[i] == ' ' )
		{
			if ( !s.empty() )
			{
				strve.push_back(s);
				s = "";
			}
		} else
		{
			s += instr[i];
		}
	}
	if ( !s.empty() )
	{
		strve.push_back(s);
	}
	return strve;
}

double log_2(double r )
{
	double res = log(r)/log(2.0);
	return res;
}

bool overlaptest(Position &p1, Position &p2)
{
	if ( p1.chr != p2.chr )
		return false;

	if ( p1.start > p2.end || p1.end < p2.start )
		return false;
	else
		return true;
}

bool overlaptest(pair<int, int> p1, pair<int, int> p2)
{
	if ( p1.first > p2.second || p1.second < p2.first )
		return false;
	else
		return true;
}

bool overlaptest(Position &p1, Position &p2, int extension)
{
	if ( p1.chr != p2.chr )
		return false;

	if ( p1.start > (p2.end+extension) || p1.end < (p2.start-extension) )
		return false;
	else
		return true;
}

bool overlaptest(pair<int, int> p1, pair<int, int> p2, int extension)
{
	if ( p1.first > (p2.second+extension) || p1.second < (p2.first-extension) )
		return false;
	else
		return true;
}

bool overlaptest(pair<int, int> p1, map< pair<int, int >, Region_id > &posmap )
{
	bool ovl = false;
	for ( map< pair<int, int >, Region_id >::iterator ite = posmap.begin(); ite != posmap.end(); ++ite )
	{
		if ( ite->first.second < p1.first )
			continue;
		else if ( ite->first.first > p1.second )
			break;
		else
			ovl = true;
			break;
		
	}
	return ovl;
}

bool overlaptest(pair<int, int> p1, map< pair<int, int >, Region_id > &posmap, vector<Region_id > &idve )
{
	bool ovl = false;
	for ( map< pair<int, int >, Region_id >::iterator ite = posmap.begin(); ite != posmap.end(); ++ite )
	{
		if ( ite->first.second < p1.first )
			continue;
		else if ( ite->first.first > p1.second )
			break;
		else
			ovl = true;
			idve.push_back( ite->second );
			
	}
	return ovl;
}

bool overlaptest(pair<int, int> p1, map< pair<int, int >, Region_id > &posmap, vector<Region_id > &idve, int extension )
{
	bool ovl = false;
	for ( map< pair<int, int >, Region_id >::iterator ite = posmap.begin(); ite != posmap.end(); ++ite )
	{
		if ( ite->first.second < p1.first - extension )
			continue;
		else if ( ite->first.first > p1.second + extension )
			break;
		else
			ovl = true;
			idve.push_back( ite->second );
			
	}
	return ovl;
}

bool overlaptest( int site1, map<int, Region_id > &posmap, int extension, vector<int > &dis )
{
	bool ovl = false;
	for ( map<int, Region_id >::iterator ite = posmap.begin(); ite != posmap.end(); ++ite )
	{
		if ( ite->first < site1 - extension )
			continue;
		else if ( ite->first > site1 + extension )
			break;
		else
		{
			ovl = true;
			dis.push_back( ite->first - site1);
		}
	}
	return ovl;
}

bool expressionup( double rc1, double rc2 )
{
	double fd = (rc2 + 0.01) / (rc1+0.01);
	if ( fd > 2 )
		return true;
	else
		return false;
}

bool signalon( int rc1, int rc2 )
{
	if ( rc1 == 0 && rc2 > 0 )
		return true;
	else
		return false;
}

bool signaloff( int rc1, int rc2 )
{
	if ( rc1 > 0 && rc2 == 0 )
		return true;
	else
		return false;
}

bool singalup( int rc1, int rc2, double chg1, double chg2 )
{
	if ( signalon(rc1, rc2 ) )
		return true;
	else
	{
		if ( double(rc2 + 0.01) / double(rc1 + 0.01) > 2 && (chg2+0.01) / (chg1+0.01) > 2 )
			return true;
		else
			return false;
	}
}

bool signalup( int rc1, int rc2 )
{
	if ( signalon(rc1, rc2 ) )
		return true;
	else
	{
		if ( double(rc2 + 0.01) / double(rc1 + 0.01) > 2 )
			return true;
		else
			return false;
	}
}

bool signaldown( int rc1, int rc2 )
{
	return signalup(rc2, rc1);
}

bool singaldown( int rc1, int rc2, double chg1, double chg2 )
{
	if ( singalup( rc2, rc1, chg2, chg1 ) )
		return true;
	else
		return false;
}

bool markon( double sig1, double sig2 )
{
	if ( sig1 == 0 && sig2 > 0 )
		return true;
	else
		return false;
}

bool markoff( double sig1, double sig2 )
{
	return markon( sig2, sig1 );
}

bool markon( int length1, int length2 )
{
	return signalon( length1, length2 );
}

bool markoff( int length1, int length2 )
{
	return signaloff( length1, length2 );
}

bool markwide( int length1, int length2 )
{
	if ( markon(length1, length2 ) )
		return true;
	else
	{
		if ( length2 > length1 )
			return true;
		else
			return false;
	}
}

bool marknarrow( int length1, int length2 )
{
	return markwide( length2, length1 );
}

bool markup( double signal1, double signal2 )
{
	if ( (signal2+0.01) / (signal1+0.01) > 2 )
		return true;
	else
		return false;
}
bool markdown( double signal1, double signal2 )
{
	if ( markup( signal2, signal1 ) )
		return true;
	else
		return false;
}

bool consregionmore( int rc1, int rc2 )
{
	if ( (rc2 + 0.1) / (rc1 + 0.1) > 2 )
		return true;
	else
		return false;
	
}

bool consregionless( int rc1, int rc2 )
{
	return consregionmore( rc2, rc1 );
}

int matrixtrans( int stage, int elem, vector<int > &matrix )
{
	int n = 0;
	if ((int)matrix.size() != stage )
	{
		cout<<"error in marixtrans: "<<matrix.size()<<", "<<stage<<endl;
		exit(1);
	}
	for ( size_t i = 0; i < matrix.size(); ++i )
	{
		if ( matrix[i] >= elem )
		{
			cout<<"error in marixtrans: "<<matrix[i]<<", "<<elem<<endl;
			exit(1);
		}
		n += matrix[i] * (int)pow((double)elem, (double)i );
	}
	return n;
}

vector<int > revtransmatrix( int stage, int elem, int digit )
{
	vector<int> matrix;
	for ( int i = 0; i < stage; ++i )
		matrix.push_back(0);
	int maxn = (int)pow((double)elem, (double)stage );
	if ( digit >= maxn )
	{
		cout<<"error in revtransmatrix "<<stage<<","<<elem<<","<<digit<<endl;
		exit(1);
	}
	for (int i = 0; i < stage; ++i )
	{
		matrix[i] = digit % elem;
		digit -= matrix[i];
		digit /= elem;
	}

	return matrix;
}

void getposfrom( map<string, vector< pair<int, int> > > &region, 
				map<string, map<pair<int, int>, string > > &region_map )
{
	region.clear();
	for ( map<string, map<pair<int, int>, string > >::iterator ite = region_map.begin(); ite != region_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int>, string >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			region[chr].push_back(seci->first);
		}
	}
}

void getposfrom( map<string, vector< pair<int, int> > > &region, 
				map<string, map<pair<int, int>, Region_id > > &region_map )
{
	region.clear();
	for ( map<string, map<pair<int, int>, Region_id > >::iterator ite = region_map.begin(); ite != region_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int>, Region_id >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			region[chr].push_back(seci->first);
		}
	}
}

void filter_region_remain( map<string, vector< pair<int, int> > > &filtered, 
				   map<string, vector< pair<int, int > > > &ori_region,
				   map<string, vector< pair<int, int > > > &fr )
{
	filtered.clear();
	for ( map<string, vector< pair<int, int > > >::iterator ite = ori_region.begin(); 
		ite != ori_region.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector< pair<int, int > >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			if ( fr.find(chr) == fr.end() )
				continue;
			else
			{
				bool overlap = false;
				for ( size_t i = 0; i < fr[chr].size(); ++i )
				{
					if ( fr[chr][i].second < seci->first )
						continue;
					if ( fr[chr][i].first > seci->second )
						break;

					overlap = true;
					break;
				}
				if ( overlap )
					filtered[chr].push_back(*seci);
			}
		}
	}
}

void filter_region_minus( map<string, vector< pair<int, int> > > &filtered, 
				   map<string, vector< pair<int, int > > > &ori_region,
				   map<string, vector< pair<int, int > > > &fr )
{
	filtered.clear();
	for ( map<string, vector< pair<int, int > > >::iterator ite = ori_region.begin(); 
		ite != ori_region.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector< pair<int, int > >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			if ( fr.find(chr) == fr.end() )
				filtered[chr].push_back( *seci );
			else
			{
				bool overlap = false;
				for ( size_t i = 0; i < fr[chr].size(); ++i )
				{
					if ( fr[chr][i].second < seci->first )
						continue;
					if ( fr[chr][i].first > seci->second )
						break;

					overlap = true;
					break;
				}
				if ( !overlap )
					filtered[chr].push_back(*seci);
			}
		}
	}
}

void filter_region_minus_remain( map<string, vector< pair<int, int> > > &overlapped, 
				   map<string, vector< pair<int, int> > > &nonoverlapped,
				   map<string, vector< pair<int, int > > > &ori_region,
				   map<string, vector< pair<int, int > > > &fr )
{
	overlapped.clear();
	nonoverlapped.clear();
	for ( map<string, vector< pair<int, int > > >::iterator ite = ori_region.begin(); 
		ite != ori_region.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector< pair<int, int > >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			if ( fr.find(chr) == fr.end() )
				nonoverlapped[chr].push_back( *seci );
			else
			{
				bool overlap = false;
				for ( size_t i = 0; i < fr[chr].size(); ++i )
				{
					if ( fr[chr][i].second < seci->first )
						continue;
					if ( fr[chr][i].first > seci->second )
						break;

					overlap = true;
					break;
				}
				if ( !overlap )
					nonoverlapped[chr].push_back(*seci);
				else
					overlapped[chr].push_back(*seci);
			}
		}
	}
}

// center in region 1000 flanking
void find_overlap_bycenterinregionplus(map<string, vector< pair<int, int> > > &overlapped, 
	map<string, vector< pair<int, int> > > &ori_region,
	map<string, vector< pair<int, int > > > &base_region,
	int exten )
{
	overlapped.clear();
	for ( map<string, vector< pair<int, int> > >::iterator ite = ori_region.begin(); 
		ite != ori_region.end(); ++ite )
	{
		string chr = ite->first;
		if ( base_region.find( chr ) ==  base_region.end() )
			continue;
		for ( vector< pair<int, int > >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			int center = (seci->second - seci->first) / 2 + seci->first;
			bool overlap = false;
			for ( size_t i = 0; i < base_region[chr].size(); ++i )
			{
				if ( base_region[chr][i].second < center - exten )
					continue;
				if ( base_region[chr][i].first > center + exten )
					break;

				overlap = true;
				break;
			}
			if ( overlap )
				overlapped[chr].push_back(*seci);
		}
	}
}

void region_merge_naive( vector< map<string, vector< pair<int, int> > > > &reg_ve,
	map<string, vector< pair<int, int> > > &merged )
{
	map<string, set<pair<int, int > > > chr_posset;
	for ( size_t i = 0; i < reg_ve.size(); ++i )
	{
		for ( map<string, vector<pair<int, int> > >::iterator ite = reg_ve[i].begin(); ite != reg_ve[i].end(); ++ite )
		{	
			string chr = ite->first;
			for ( size_t k = 0; k < ite->second.size(); ++k )
			{
				chr_posset[chr].insert(ite->second[k]);
			}
		}
	}
	
	for ( map<string, set<pair<int, int > > >::iterator ite = chr_posset.begin(); ite != chr_posset.end(); ++ite )
	{
		vector< pair<int, int > > m;
		m.push_back( *(ite->second.begin()) );
		for ( set<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			if ( si->first <= m.back().second )
			{
				if ( si->second > m.back().second )
					m.back().second = si->second;
			} else
			{
				m.push_back( *si );
			}
		}
		merged.insert( make_pair(ite->first, m ) );
	}
}

void stitch_region( map<string, vector<vector<pair<int, int> > > > &stitched,
				   map<string, vector<pair<int, int > > > &region,
				   int distance )
{
	stitched.clear();
	for ( map<string, vector<pair<int, int > > >::iterator ite = region.begin(); ite != region.end(); ++ite )
	{
		string chr = ite->first;
		size_t i = 0;
		vector<pair<int, int > > clu;
		clu.push_back(ite->second[i]);
		i++;
		for ( ; i < ite->second.size(); ++i )
		{
			if ( ite->second[i].first - ite->second[i-1].second < distance )
			{
				clu.push_back( ite->second[i]);
			} else
			{
				stitched[chr].push_back( clu );
				clu.clear();
				clu.push_back( ite->second[i] );
			}
		}
		stitched[chr].push_back( clu );
	}

}

void stitch_region( map<string, vector<Stitched_Region > > &stitched,
				   map<string, vector<pair<int, int > > > &region,
				   int distance )
{
	stitched.clear();
	for ( map<string, vector<pair<int, int > > >::iterator ite = region.begin(); ite != region.end(); ++ite )
	{
		string chr = ite->first;
		size_t i = 0;
		Stitched_Region clu;
		stitched[chr].push_back(clu);
		stitched[chr].back().chr = chr;
		stitched[chr].back().regions.push_back(ite->second[i]);
		i++;
		for ( ; i < ite->second.size(); ++i )
		{
			if ( ite->second[i].first - ite->second[i-1].second < distance )
			{
				stitched[chr].back().regions.push_back( ite->second[i]);
			} else
			{
				stitched[chr].push_back( clu );
				stitched[chr].back().chr = chr;
				stitched[chr].back().regions.push_back( ite->second[i] );
			}
		}
	}

}

void cluster_region_CrSample( map<string, vector< vector<pair<int, int> > > > &clustered,
							 vector< map<string, vector<pair<int, int> > > > &posVec,
							 int overlap_extension )
{
	map<string, set<pair<int, int > > > posmap;
	for ( size_t i = 0; i < posVec.size(); ++i )
	{
		for ( map<string, vector<pair<int, int> > >::iterator ite = posVec[i].begin(); ite != posVec[i].end(); ++ite )
		{
			string chr = ite->first;
			for ( vector<pair<int, int > >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
			{
				posmap[chr].insert( *seci );
			}
		}
	}

	clustered.clear();
	
	for ( map<string, set<pair<int, int > > >::iterator ite = posmap.begin(); ite != posmap.end(); ++ite )
	{
		string chr = ite->first;
		vector< pair<int, int > > miniv;
		set<pair<int, int > >::iterator seci = ite->second.begin();
		pair<int, int > prepos = *seci;
		miniv.push_back( prepos );
		++seci;
		for ( ; seci != ite->second.end(); ++seci )
		{
			pair<int, int> nextpos = *seci;
			if ( overlaptest( prepos, nextpos, overlap_extension ) )
			{
				miniv.push_back( nextpos );
				if ( nextpos.second > prepos.second )
					prepos = nextpos;
			} else
			{
				clustered[chr].push_back( miniv );
				miniv.clear();
				miniv.push_back( nextpos );
				prepos = nextpos;
			}
		}
		clustered[chr].push_back( miniv );
	}
}

string inttostr(int i )
{
	string Res;
	ostringstream convert;
	convert << i;
	Res = convert.str();
	return Res;
}

void outputtable( ofstream &outf, vector<vector<int> > &mat, vector<string > &rowname, vector<string> &colname )
{
	
	for ( size_t i = 0; i < colname.size(); ++i )
	{
		outf<<"\t"<<colname[i];
	}
	outf<<endl;
	if ( rowname.size() != mat.size() )
	{
		cout<<"error output table rowname size != mat size "<<rowname.size()<<","<<mat.size()<<endl; exit(1);
	}
	for ( size_t i = 0; i < rowname.size(); ++i )
	{
		outf<<rowname[i];
		if ( mat[i].size() != colname.size() )
		{
			cout<<"error output table mat[i] size != colname size "<<colname.size()<<", "<<mat[i].size()<<endl; exit(1);
		}
		for ( size_t j = 0; j < mat[i].size(); ++j) 
		{
			outf<<"\t"<<mat[i][j];
		}
		outf<<endl;
	}
}

void outputtable( ofstream &outf, vector<vector<double> > &mat, vector<string > &rowname, vector<string> &colname )
{
	
	for ( size_t i = 0; i < colname.size(); ++i )
	{
		outf<<"\t"<<colname[i];
	}
	outf<<endl;
	if ( rowname.size() != mat.size() )
	{
		cout<<"error output table rowname size != mat size "<<rowname.size()<<","<<mat.size()<<endl; exit(1);
	}
	for ( size_t i = 0; i < rowname.size(); ++i )
	{
		outf<<rowname[i];
		if ( mat[i].size() != colname.size() )
		{
			cout<<"error output table mat[i] size != colname size "<<colname.size()<<", "<<mat[i].size()<<endl; exit(1);
		}
		for ( size_t j = 0; j < mat[i].size(); ++j) 
		{
			outf<<"\t"<<setprecision(3)<<mat[i][j];
		}
		outf<<endl;
	}
}

vector<vector<double> > transposmat( vector<vector<double> > &mat )
{
	if ( mat.empty() )
	{
		cout<<"error mat empty transposmat"<<endl; exit(1);
	}
	size_t il = mat.size();
	size_t jl = mat[0].size(); 
	if ( jl == 0 )
	{
		cout<<"error mat[0] empty transposmat"<<endl; exit(1);
	}
	for ( size_t i = 0; i < il; ++i )
	{
		if ( mat[i].size() != jl )
		{	
			cout<<"error mat i "<<i<<" != mat[0] "<<jl<<endl; exit(1);
		}
		
	}
	vector<vector<double > > res;
	for ( size_t j = 0; j < jl; ++j )
	{
		vector< double > subv;
		for ( size_t i = 0; i < il; ++i )
		{
			subv.push_back( mat[i][j] );
		}
		res.push_back( subv );
	}
	return res;
}

void normalizemat( vector<vector<double > > &mat, vector<double > &thr )
{
	
	for ( size_t i = 0; i < mat.size(); ++i )
	{
		if ( mat[i].size() != thr.size() )
		{
			cout<<"error in normalizemat: mat[i]size != thr size "<<mat[i].size()<<","<<thr.size()<<endl; exit(1);
		}
		for ( size_t j = 0; j < mat[i].size(); ++j )
		{
			if ( thr[j] > 0 )
			{
				mat[i][j] /= thr[j];
				if ( mat[i][j] > 1 )
					mat[i][j] = 1;
			}
		}
	}
}

int findnearTSSinchr( set<int > &tssset, pair<int, int > region )
{
	int distance = -1;
	int selectTSS = -1;
	for ( set<int>::iterator ite = tssset.begin(); ite != tssset.end(); ++ite )
	{
		int td = min(abs(region.first - *ite), abs(region.second - *ite) );
		if ( distance == -1 )
		{
			selectTSS = *ite;
			distance = td;
		} else
		{
			if ( td < distance )
			{
				distance = td;
				selectTSS = *ite;
			}
		}
		if ( *ite < region.first )
			continue;
		if ( *ite > region.second )
			break;
		
	}
	return selectTSS;
}

int findnearTSSinchr( set<int > &tssset, pair<int, int > region, vector< int > &containedTSS )
{
	int distance = -1;
	int selectTSS = -1;
	for ( set<int>::iterator ite = tssset.begin(); ite != tssset.end(); ++ite )
	{
		int td = min(abs(region.first - *ite), abs(region.second - *ite) );
		if ( distance == -1 )
		{
			selectTSS = *ite;
			distance = td;
		} else
		{
			if ( td < distance )
			{
				distance = td;
				selectTSS = *ite;
			}
		}
		if ( *ite < region.first )
			continue;
		if ( *ite > region.second )
			break;
		
		containedTSS.push_back( *ite );
	}
	return selectTSS;
}

string onofftagmeaning( int tag, int stage )
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

string updowntagmeaning( int tag, int stage )
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
			m += "Up";
		} else if ( matrix1[i] == 2 )
		{
			m += "Down";
		}
		if ( i != matrix1.size()-1)
		{
			m += '_';
		}
	}
	return m;
}

vector<double > getqualter( vector< vector<vector<double > > > &mat_vec )
{
	multiset<double > d;
	for ( size_t i = 0; i < mat_vec.size(); ++i )
	{
		for ( size_t j = 0; j < mat_vec[i].size(); ++j )
		{
			for ( size_t k = 0; k < mat_vec[i][j].size(); ++k )
			{
				if ( mat_vec[i][j][k] > 0 )
					d.insert( mat_vec[i][j][k] );
			}
		}
	}
	int total = (int)d.size();
	int qu = (int)(0.1*total);
	int i = 0;
	vector<double> res;
	res.push_back(*d.begin());
	for ( multiset<double >::iterator ite = d.begin(); ite != d.end(); ++ite )
	{
		i++;
		if ( i == qu || i == 2*qu || i == 3*qu || i == 4*qu || i == 5*qu || i == 6*qu || i == 7*qu || i == 8*qu || i == 9*qu )
			res.push_back(*ite);
		
	}
	cout<<d.size()<<endl;
	res.push_back( *d.rbegin() );
	for ( size_t i = 0; i < res.size(); ++i )
		cout<<res[i]<<" ";
	cout<<endl;
	return res;
}

vector<vector<double > > mergecol(vector< vector<vector<double > > > &mat_vec )
{
	vector<vector<double > > res;
	// check row line
	if ( mat_vec.empty() )
	{
		cout<<"error mat_vec empty in mergecol()"<<endl; exit(1);
	}
	size_t rs = mat_vec[0].size();
	if ( mat_vec.size() == 1 )
	{
		res = mat_vec[0];
	} else
	{
		for ( size_t i = 1; i < mat_vec.size(); ++i )
		{
			if ( mat_vec[i].size() != rs )
			{
				cout<<"row n not consistent! "<<mat_vec[i].size()<<", "<<rs<<endl; exit(1);
			} 
		}
		for ( size_t i = 0; i < rs; ++i )
		{
			vector<double > subv;
			for ( size_t j = 0; j < mat_vec.size(); ++j )
			{
				subv.insert( subv.end(), mat_vec[j][i].begin(), mat_vec[j][i].end() );
			} 
			res.push_back( subv );
		}
	}
	return res;
}

int intersection_num( vector<set<Region_id > > &sets )
{
	if ( sets.empty() )
	{
		cout<<"error sets empty!"<<endl; exit(1);
	}
	set<Region_id > iniset = sets[0];
	for ( size_t i = 1; i < sets.size(); ++i )
	{
		set<Region_id > alset;
		for ( set< Region_id >::iterator ite = iniset.begin(); ite != iniset.end(); ++ite )
		{
			if (sets[i].find(*ite) != sets[i].end() )
				alset.insert(*ite);
		}
		iniset = alset;
	}
	return (int)iniset.size();
}

void venn_quintuple(ofstream &outf, vector<set<Region_id > > &sets, vector<string > &names )
{
	if ( (int)sets.size() != 5 )
	{
		cout<<"error sets size != 5 in venn_quituple "<<sets.size()<<endl; exit(1);
	}
	outf <<"venn.plot <- draw.quintuple.venn("<<endl;
	for ( size_t i = 0; i < sets.size(); ++i )
	{
		outf<<"area"<<(i+1)<<" = "<<sets[i].size()<<", ";
	}
	outf<<endl;
	for ( size_t i = 0; i < sets.size(); ++i )
	{
		for ( size_t j = i+1; j < sets.size(); ++j )
		{
			vector<set<Region_id > > asets;
			asets.push_back( sets[i] );
			asets.push_back( sets[j] );
			int inter = intersection_num( asets );
			outf<<"n"<<(i+1)<<(j+1)<<" = "<<inter<<", ";
		}
	}
	outf<<endl;
	for ( size_t i = 0; i < sets.size(); ++i )
	{
		for ( size_t j = i+1; j < sets.size(); ++j )
		{
			for ( size_t k = j +1; k < sets.size(); ++k )
			{
				vector<set<Region_id > > asets;
				asets.push_back( sets[i] );
				asets.push_back( sets[j] );
				asets.push_back( sets[k] );
				int inter = intersection_num( asets );
				outf<<"n"<<(i+1)<<(j+1)<<(k+1)<<" = "<<inter<<", ";
			}
		}
	}
	outf<<endl;
	for ( size_t i = 0; i < sets.size(); ++i )
	{
		for ( size_t j = i+1; j < sets.size(); ++j )
		{
			for ( size_t k = j +1; k < sets.size(); ++k )
			{
				for ( size_t l = k+1; l < sets.size(); ++l )
				{
					vector<set<Region_id > > asets;
					asets.push_back( sets[i] );
					asets.push_back( sets[j] );
					asets.push_back( sets[k] );
					asets.push_back( sets[l] );
					int inter = intersection_num( asets );
					outf<<"n"<<(i+1)<<(j+1)<<(k+1)<<(l+1)<<" = "<<inter<<", ";
				}
			}
		}
	}
	outf<<endl;
	int inter = intersection_num(sets);
	outf<<"n12345 = "<<inter<<","<<endl;
	outf<<"category = c(";
	for ( size_t i = 0; i < names.size(); ++i )
	{
		outf<<"\""<<names[i]<<"\"";
		if ( i != names.size()-1)
			outf<<",";
		else
			outf<<"), "<<endl;
	} 
	outf<<"fill=c(\"dodgerblue\",\"goldenrod1\",\"darkorange1\", \"seagreen3\", \"orchid3\"),"<<endl;
	outf<<"cat.col = c(\"dodgerblue\", \"goldenrod1\", \"darkorange1\", \"seagreen3\", \"orchid3\"),"<<endl;
	outf<<"cat.cex = 2, margin = 0.1, cex=c(1.5, 1.5, 1.5, 1.5, 1.5, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8, 1, 0.8,"<<endl;
	outf<<"1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 0.55, 1, 1, 1, 1, 1, 1.5),"<<endl;
	outf<<"ind=TRUE"<<endl;
	outf<<");"<<endl;
}

void venn_triple(ofstream &outf, vector<set<Region_id > > &sets, vector<string > &names, string outname )
{
	if ( (int)sets.size() != 3 )
	{
		cout<<"error sets size != 3 in venn_triple "<<sets.size()<<endl; exit(1);
	}
	outf <<"library(\"grid\")"<<endl;
	outf <<"library(\"VennDiagram\")"<<endl;
	outf <<"venn.plot <- draw.triple.venn("<<endl;
	for ( size_t i = 0; i < sets.size(); ++i )
	{
		outf<<"area"<<(i+1)<<" = "<<sets[i].size()<<", ";
	}
	outf<<endl;
	for ( size_t i = 0; i < sets.size(); ++i )
	{
		for ( size_t j = i+1; j < sets.size(); ++j )
		{
			vector<set<Region_id > > asets;
			asets.push_back( sets[i] );
			asets.push_back( sets[j] );
			int inter = intersection_num( asets );
			outf<<"n"<<(i+1)<<(j+1)<<" = "<<inter<<", ";
		}
	}
	outf<<endl;
	for ( size_t i = 0; i < sets.size(); ++i )
	{
		for ( size_t j = i+1; j < sets.size(); ++j )
		{
			for ( size_t k = j +1; k < sets.size(); ++k )
			{
				vector<set<Region_id > > asets;
				asets.push_back( sets[i] );
				asets.push_back( sets[j] );
				asets.push_back( sets[k] );
				int inter = intersection_num( asets );
				outf<<"n"<<(i+1)<<(j+1)<<(k+1)<<" = "<<inter<<", ";
			}
		}
	}
	outf<<endl;
	outf<<"category = c(";
	for ( size_t i = 0; i < names.size(); ++i )
	{
		outf<<"\""<<names[i]<<"\"";
		if ( i != names.size()-1)
			outf<<",";
		else
			outf<<"), "<<endl;
	} 
	outf<<"fill = c(\"blue\", \"red\", \"green\"),"<<endl;
	outf<<"lty=\"blank\", margin = 0.1, cex=3, cat.cex=3, cat.col=c(\"blue\", \"red\", \"green\")"<<endl;
	outf<<");"<<endl;
	outf<<"pdf(\""<<outname<<"\");"<<endl;
	outf<<"grid.draw(venn.plot);"<<endl;
	outf<<"dev.off()"<<endl;
}

void venn_threeset( vector<Region_id > &set1, vector<Region_id >& set2, 
	vector<Region_id > &set1_exc, vector<Region_id > &set2_exc, vector<Region_id >& overlap )
{
	set1_exc.clear();
	set2_exc.clear();
	overlap.clear();
	
	set<Region_id > set1s;
	set<Region_id > overs;
	for ( size_t i = 0; i < set1.size(); ++i )
	{
		set1s.insert( set1[i] );
	}
	for ( size_t i = 0; i < set2.size(); ++i )
	{
		if ( set1s.find(set2[i] ) != set1s.end() )
		{
			overs.insert( set2[i] );
			overlap.push_back( set2[i] );
		} else
		{
			set2_exc.push_back( set2[i] );
		}
		
	}
	for ( size_t i = 0; i < set1.size(); ++i )
	{
		if ( overs.find(set1[i]) == overs.end() )
		{
			set1_exc.push_back( set1[i] );
		}
	}
}


double match_rank( int rc, multimap<int, Region_id > &sig_r_map )
{
	int t = (int)sig_r_map.size();
	if ( sig_r_map.empty() )
	{
		cout<<"error in match_rank sig_r_map empty() "<<endl; exit(1);
	}
	int i = 0;
	for ( multimap<int, Region_id >::reverse_iterator ri = sig_r_map.rbegin(); ri != sig_r_map.rend(); ++ri )
	{
		i++;
		multimap<int, Region_id >::reverse_iterator nri = ri;
		++nri;
		
		if ( nri == sig_r_map.rend() )
		{
			return (double)i/t;
		}
		
		if ( ri->first <= rc && nri->first >= rc )
		{
			return (double)i/t;
		}
		
	}
	return 0;
}

double match_rank( double sig, multimap<double, Region_id > &sig_r_map )
{
	int t = (int)sig_r_map.size();
	if ( sig_r_map.empty() )
	{
		cout<<"error in match_rank sig_r_map empty() "<<endl; exit(1);
	}
	int i = 0;
	for ( multimap<double, Region_id >::reverse_iterator ri = sig_r_map.rbegin(); ri != sig_r_map.rend(); ++ri )
	{
		i++;
		multimap<double, Region_id >::reverse_iterator nri = ri;
		++nri;
		
		if ( nri == sig_r_map.rend() )
		{
			return (double)i/t;
		}
		
		if ( ri->first <= sig && nri->first >= sig )
		{
			return (double)i/t;
		}
		
	}
	return 0;
}

set<Region_id > vetoset(vector<Region_id > &ve )
{
	set<Region_id > rs;
	for ( size_t i = 0; i < ve.size(); ++i )
		rs.insert( ve[i]);
	return rs;
}

vector<pair<int, int > > mergeragions( set<pair<int, int > > &r )
{
	vector<pair<int, int > > mr;
	if ( r.empty() )
		return mr;
	
	mr.push_back( *r.begin() );
	for ( set<pair<int, int > >::iterator ite = r.begin(); ite != r.end(); ++ite )
	{
		if ( mr.back().second >= ite->first )
		{
			mr.back().second = max( ite->second, mr.back().second );
		} else
		{
			mr.push_back( *ite );
		}
	}
	return mr;
}

int shift_tagpos( int start, int end, char strand, int seg_len )
{
	if ( seg_len == 0 )
	{
		return start + (end-start)/2;
	}
	int p = 0;
	if ( strand == '+' )
	{
		return start + seg_len/2;
	} else if ( strand == '-')
		return end - 1 - seg_len/2;
	else
	{
		cout<<"error strand "<<strand<<endl; 
		exit(1);
		return 0;
	}
}

void assigntagpos( map<string, vector<int > > &tagpos, string chr, int start, int end, int win, int total, 
	vector<int > &counts, vector<double > &rpkm )
{
	vector<pair<int, int > > cut_reg_ve; 
	counts.clear();
	rpkm.clear();
	while ( start < end )
	{
		int tend = start + win - 1;
		if ( tend > end )
			tend = end;
		cut_reg_ve.push_back( make_pair(start, tend ) );
		counts.push_back( 0 );
		rpkm.push_back( 0 );
		start = tend + 1;
	}
	
	for ( vector<int >::iterator ite = tagpos[chr].begin(); ite != tagpos[chr].end(); ++ite )
	{
		for ( size_t i = 0; i < cut_reg_ve.size(); ++i )
		{
			if ( cut_reg_ve[i].second < *ite )
				continue;
			if ( cut_reg_ve[i].first > *ite )
				break;
			counts[i] += 1;
			break;
		}
	}
	
	for ( size_t i = 0; i < cut_reg_ve.size(); ++i )
	{
		int w = cut_reg_ve[i].second - cut_reg_ve[i].first + 1;
		double nc = ( counts[i] * 1.0 * 1000 * 1000000 ) / ( w * total );
		rpkm[i] = nc;
	}
}

void assigntagpos( int tagpos, vector<pair<int, int > > &cut_reg, vector<int > &counts )
{
	
	
	for ( size_t i = 0; i < cut_reg.size(); ++i )
	{
		if ( cut_reg[i].second < tagpos )
			continue;
		if ( cut_reg[i].first > tagpos )
			break;
		counts[i] += 1;
		break;
	}
	
	
	
}


int gettotaltagnumber( map<string, vector<int > > &tagpos )
{
	int n = 0;
	for ( map<string, vector<int > >::iterator ite = tagpos.begin(); ite != tagpos.end(); ++ite )
	{
		n += (int)ite->second.size();
	}
	return n;
}

void averagetablecolumn( vector< double > &ave_ve, vector<vector<double > > & rpkm_table )
{
	if ( rpkm_table.empty() )
	{
		cout<<"error in averagetablecolumn rpkm_table empty "<<endl; exit(1);
	} 
	size_t ncol = rpkm_table[0].size();
	for ( size_t i = 0; i < ncol; ++i )
	{
		ave_ve.push_back(0);
	}
	for ( size_t i = 0; i < rpkm_table.size(); ++i )
	{
		if ( rpkm_table[i].size() != ncol )
		{
			cout<<"error rpkm_table col i size "<<i<<" "<<rpkm_table[i].size()<<"!="<<ncol<<endl;
			exit(1);
		}
		for ( size_t j = 0; j < ncol; ++j )
		{
			ave_ve[j] += rpkm_table[i][j];
		}
	}
	int nrow = (int)rpkm_table.size();
	for ( size_t i = 0; i < ncol; ++i )
	{
		ave_ve[i] = ave_ve[i]*1.0 / nrow;
	}
}

void assigncountpool( map<string, vector<int > > &tagpos, Count_pool &cpool )
{
	int L = 100000;
	map<string, map<int, set<pair<int, int > > > > chr_index_region;
	for ( map<string, map<pair<int, int >, char > >::iterator ite = cpool.region_strand_map.begin();
		ite != cpool.region_strand_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, char >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int start = si->first.first;
			int end = si->first.second;
			int index1 = start / L;
			int index2 = end / L;
			chr_index_region[chr][index1].insert(  si->first  );
			if ( index2 > index1 )
			{
				for ( int i = index1+1; i <= index2; ++i )
					chr_index_region[chr][i].insert( si->first );
			}
		}
	} 

	int i = 0;
	
	for ( map<string, vector<int > >::iterator ite = tagpos.begin(); ite != tagpos.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			i++;
			if ( i % 100000 == 0 )
				cout<<i<<endl;
			
			int index = *si / L;
			if ( chr_index_region[chr].find(index) != chr_index_region[chr].end() )
			{
				for ( set<pair<int, int > >::iterator ci = chr_index_region[chr][index].begin();
					ci != chr_index_region[chr][index].end(); ++ci )
				{
					if ( ci->second < *si )
						continue;
					if ( ci->first > *si )
						break;
					size_t id = cpool.region_id_map[chr][*ci];
					assigntagpos( *si, cpool.cut_reg_ve[id], cpool.counts_table[id] );
				}
			}
				
		/*	for ( map<pair<int, int >, char >::iterator pi = cpool.region_strand_map[chr].begin(); 
				pi != cpool.region_strand_map[chr].end(); ++pi )
			{
				if ( pi->first.second < *si )
					continue;
				if ( pi->first.first > *si )
					break;
				size_t id = cpool.region_id_map[chr][pi->first];
				assigntagpos( *si, cpool.cut_reg_ve[id], cpool.counts_table[id]);
			}  */
		}
	}
}







