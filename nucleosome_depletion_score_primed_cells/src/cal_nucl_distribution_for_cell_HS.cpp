#include "Nucl_space.h"
#include "operation.h"


void readinHS( string infile, map<string, set<int> > &HS )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
		
	string line;
	while(!inf.eof())
	{
		getline(inf, line);
		if ( line.empty() )
			break;
		vector<string> ps = parse_string( line );
		string chr = ps[0];
		int pos = atoi(ps[1].c_str() );
		HS[chr].insert( pos );
	}
	inf.close();
}

void readinHS2( string infile, map<string, set<int> > &HS )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
		
	string line;
	getline(inf, line);
	string chr = line;
	while(!inf.eof())
	{
		getline(inf, line);
		if ( line.empty() )
			break;
		if ( line.substr(0,3) == "chr" )
			chr = line;
		else
		{
			int pos = atoi(line.c_str() );
			HS[chr].insert( pos );
		}
	}
	inf.close();
}

void mapping( map<string, set<int> > &HS, 
	map<string, map<string, set<int > > > &cell_nucl_pos,
	map< string, map<string, map<int, vector< int > > > > &cell_HS_relpos )
{
	int L = 100000;
	map<string, map<int, set<int > > > index_HS;
	for ( map<string, set<int> >::iterator ite = HS.begin(); ite != HS.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int index = *si / L;
			index_HS[chr][index].insert( *si );
		}
			
	}
	
	for ( map<string, map<string, set<int > > >::iterator ite = cell_nucl_pos.begin();
		ite != cell_nucl_pos.end(); ++ite )
	{
		string cell = ite->first;
		for ( map<string, set<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			string chr = si->first;
			for ( set<int >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int index = *ti / L;
				if ( index_HS[chr].find( index ) != index_HS[chr].end() )
				{
					for ( set<int >::iterator ci = index_HS[chr][index].begin();
						ci != index_HS[chr][index].end(); ++ci )
					{
						if ( *ci - 200 > *ti )
							break;
						if ( *ci -200 <= *ti && *ci+200 >= *ti )
						{
							int dis = (*ti - *ci);
							cell_HS_relpos[cell][chr][*ci].push_back( dis );
						}
					}
				}
			}
		}
	}
	
	
}

void distribute_for_cell( map< string, map<string, map<int, vector< int > > > > &cell_HS_relpos,
	map<string, vector<double > > &cell_distri,
	map<string, double > &cell_ratio,
	map<string, pair<int, int > > &cell_nr_tr )
{
	for ( map< string, map<string, map<int, vector< int > > > >::iterator ite = cell_HS_relpos.begin();
		ite != cell_HS_relpos.end(); ++ite )
	{
		string cell = ite->first;
		vector<int > relpos_ve;
		int val_c = 0;
		for ( map<string, map<int, vector< int > > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			string chr = si->first;
			
			for ( map<int, vector<int > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int HS = ti->first;
				if ( ti->second.size() > 4 )
					continue;
				val_c += 1;
				relpos_ve.insert( relpos_ve.end(), ti->second.begin(), ti->second.end() );
			}
			
		}
		if ( val_c < 100 )
			continue;
		
		int total = (int)relpos_ve.size();
		
		int nc = 0;
		int fc = 0;
		vector< int > c_ve;
		for ( int i = 0; i < 40; ++i )
			c_ve.push_back(0);
		int ref_total = 0;
		for ( size_t i = 0; i < relpos_ve.size(); ++i )
		{
			int index = (relpos_ve[i]+200) / 10;
			if ( index == 40 )
				continue;
			if ( index > 40 )
			{
				cout<<"error index > 40 "<<index<<endl;
				exit(1);
			}
			c_ve[index] += 1;
			ref_total += 1;
			if ( abs( relpos_ve[i] ) <= 100 )
				nc += 1;
			else
				fc += 1;
		}
		
		// smooth
		vector<double > sm_c_ve;
		for ( size_t i = 0; i < c_ve.size(); ++i )
		{
			if ( i == 0 )
			{
				double s = (c_ve[i] + c_ve[i+1]+c_ve[i+2]) / 3;
				sm_c_ve.push_back( s );
			} else if ( i == 1)
			{
				double s = (c_ve[i-1]+c_ve[i] + c_ve[i+1]+c_ve[i+2]) / 4;
				sm_c_ve.push_back( s );
			} else if ( i == c_ve.size()-1)
			{
				double s = (c_ve[i-2]+c_ve[i-1]+c_ve[i]) / 3;
				sm_c_ve.push_back( s );
			} else if ( i == c_ve.size()-2)
			{
				double s = (c_ve[i-2]+c_ve[i-1]+c_ve[i]+c_ve[i+1]) / 4;
				sm_c_ve.push_back( s );
			} else 
			{
				double s = (c_ve[i-2]+c_ve[i-1]+c_ve[i]+c_ve[i+1]+c_ve[i+2]) / 5;
				sm_c_ve.push_back( s );
			}
		}
		
		vector<double > val;
		for ( size_t i = 0; i < sm_c_ve.size(); ++i )
		{
			double v = sm_c_ve[i]*1.0 / ref_total;
			val.push_back( v );
		}
		cell_distri[cell] = val;
		
		double ratio = nc*1.0/(fc+nc);
		cell_ratio[cell] = ratio;
		cell_nr_tr[cell] = make_pair(nc, fc+nc );
	}
}

void distribute_for_HS( map< string, map<string, map<int, vector< int > > > > &cell_HS_relpos,
	map<string, map<int, vector<double > > > &HS_distri,
	map<string, map<int, double > > &HS_ratio )
{
	map<string, map<int, vector<int > > > HS_pos;
	for ( map< string, map<string, map<int, vector< int > > > >::iterator ite = cell_HS_relpos.begin();
		ite != cell_HS_relpos.end(); ++ite )
	{
		string cell = ite->first;
		
		for ( map<string, map<int, vector< int > > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			string chr = si->first;
			for ( map<int, vector<int > >::iterator ti = si->second.begin(); ti != si->second.end(); ++ti )
			{
				int HS = ti->first;
				HS_pos[chr][HS].insert( HS_pos[chr][HS].end(), ti->second.begin(), ti->second.end() );
			}
		}
	}
	
	for ( map<string, map<int, vector<int > > >::iterator ite = HS_pos.begin(); ite != HS_pos.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<int, vector<int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int HS = si->first;
			vector< int > c_ve;
			int nc = 0;
			int fc = 0;
			for ( int i = 0; i < 20; ++i )
				c_ve.push_back(0);
			int ref_total = 0;
			for ( size_t i = 0; i < si->second.size(); ++i )
			{
				int index = (si->second[i]+200) / 20;
				if ( index == 20 )
					continue;
				if ( index > 20 )
				{
					cout<<"error index > 20 "<<endl;
					exit(1);
				}
				c_ve[index] += 1;
				ref_total += 1;
				if ( abs( si->second[i] ) <= 80 )
					nc += 1;
				else
					fc += 1;
			}
			
			if ( ref_total < 20 )
				continue;
			
			vector<double > val;
			for ( size_t i = 0; i < c_ve.size(); ++i )
			{
				double v = c_ve[i]*1.0 / ref_total;
				val.push_back( v );
			}
			HS_distri[chr][HS] = val;
			double ratio = nc*1.0/(fc+nc);
			HS_ratio[chr][HS] = ratio;
		}
	}
	
}

void output_cell_distr( string prefix,
	map<string, vector<double > > &cell_distri,
	map<string, double > &cell_ratio,
	map<string, pair<int, int > > &cell_nr_tr )
{
	string outfile1 = prefix + ".cell_nuclpos_distri";
	string outfile2 = prefix + ".cell_near_far_ratio";
	ofstream outf1( outfile1.data() );
	ofstream outf2( outfile2.data() );
	for ( int i = 0; i < 60; ++i )
	{
		outf1<<"\tV"<<i+1;
	}
	outf1<<endl;
	for ( map<string, vector<double > >::iterator ite = cell_distri.begin();
		ite != cell_distri.end(); ++ite )
	{
		string cell = ite->first;
		outf1<<cell;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			outf1<<"\t"<<ite->second[i];
		}
		outf1<<endl;
	}
	outf1.close();
	
	for ( map<string, double >::iterator ite = cell_ratio.begin();
		ite != cell_ratio.end(); ++ite )
	{
		string cell = ite->first;
		outf2<<ite->first<<"\t"<<ite->second<<"\t"<<cell_nr_tr[cell].first<<"\t"<<cell_nr_tr[cell].second<<endl;
	}	
	outf2.close();
}

void output_HS_distr( string prefix,
	map<string, map<int, vector<double > > > &HS_distri,
	map<string, map<int, double > > &HS_ratio )
{
	string outfile1 = prefix + ".HS_nuclpos_distri";
	string outfile2 = prefix + ".HS_near_far_ratio";
	ofstream outf1( outfile1.data() );
	ofstream outf2( outfile2.data() );
	for ( int i = 0; i < 20; ++i )
	{
		outf1<<"\tV"<<i+1;
	}
	outf1<<endl;
	for ( map<string, map<int, vector<double > > > ::iterator ite = HS_distri.begin();
		ite != HS_distri.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<int, vector<double > >::iterator si = ite->second.begin(); 
			si != ite->second.end(); ++si )
		{
			outf1<<chr<<"_"<<si->first;
			for ( size_t i = 0; i < si->second.size(); ++i )
			{
				outf1<<"\t"<<si->second[i];
			}
			outf1<<endl;
		}
	}
	outf1.close();
	
	for ( map<string, map<int, double > >::iterator ite = HS_ratio.begin();
		ite != HS_ratio.end(); ++ite )
	{
		for ( map<int, double >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			outf2<<ite->first<<"_"<<si->first<<"\t"<<si->second<<endl;
		}
	}	
	outf2.close();
}


int main( int argc, char* argv[] )
{
	if ( argc == 1 )
	{
		cout<<"cal cell HS nucl pos distribution"<<endl;
		cout<<"Usage: prog bedfile(batch) HS_center/peak_center HS_type[0: twocol| 1: onecol] prefix"<<endl;
		exit(1);
		
	}
	
	string bedfile = argv[1];
	string HSfile = argv[2];
	int HStype=atoi(argv[3]);
	string prefix = argv[4];
	
	map<string, set<int> > HS;
	if ( HStype==0 )
		readinHS( HSfile, HS );
	else if ( HStype == 1)
		readinHS2( HSfile, HS );
	
	Cell_Nucl_Space_Bank bank;
	bank.readinnucl_bunch( bedfile );
	
	cout<<"mapping"<<endl;
	map< string, map<string, map<int, vector< int > > > > cell_HS_relpos;
	mapping( HS, bank.cell_nucl_pos, cell_HS_relpos );
	
	cout<<"distriute for cell"<<endl;
	map<string, vector<double > > cell_distri;
	map<string, double > cell_ratio;
	map<string, pair<int, int > > cell_nr_tr;
	distribute_for_cell( cell_HS_relpos, cell_distri, cell_ratio, cell_nr_tr );
	
	cout<<"output_cell_distr"<<endl;
	output_cell_distr( prefix, cell_distri, cell_ratio, cell_nr_tr );
	
/*	cout<<"distribute_for_HS"<<endl;
	map<string, map<int, vector<double > > > HS_distri;
	map<string, map<int, double > > HS_ratio;
	distribute_for_HS( cell_HS_relpos, HS_distri, HS_ratio );
	
	cout<<"output_HS_distr"<<endl;
	output_HS_distr( prefix, HS_distri, HS_ratio ); */
	
	return 0;
}


 




