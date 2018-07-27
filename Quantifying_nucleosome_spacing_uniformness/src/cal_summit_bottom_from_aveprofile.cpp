#include <iostream>
#include <string>
#include <map>
#include <list>
#include <vector>
#include <set>
#include <sstream>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <stdlib.h>
#include <stdio.h>

using namespace std;

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

void Gaussian_weight( double sigma, int ext, map<int, double > &w )
{
	for ( int j = (-1)*ext; j <= ext; ++j )
	{
		double e = exp( (-1)*(j*1.0/sigma)*(j*1.0/sigma)/2 );
		w.insert(make_pair(j, e) );
	}
}

double GW_score( map<int, double > &GW, int ext, int k, map<int, double > &score_map )
{
	double s = 0;
	for ( int j = (-1)*ext; j <= ext; ++j )
	{
		int i = j+k;
		double si = 0;
		if ( score_map.find(i) != score_map.end() )
		{
			si = score_map[i] * GW[j];
		}
		s += si;
	}	
	return s;
}

void GW_score( map<int, double > &GW, int ext, map<int, double > &score_map )
{
	map<int, double > gw_score_map;
	for ( map<int, double >::iterator ite = score_map.begin(); ite != score_map.end(); ++ite )
	{
		int k = ite->first;
		double gws = GW_score( GW, ext, k, score_map );
		gw_score_map.insert( make_pair(k, gws) );
	}
	score_map = gw_score_map;
}

void readin( string infile, map<int, double> &m )
{
	ifstream inf( infile.data() );
	
	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	
	cout<<"read in file "<<infile<<endl;
	string line;
	
	
	while(!inf.eof())
	{
		getline(inf, line );
		if ( line.empty() )
			break;
			
		vector<string > parsed_items = parse_string( line );
		int index = atoi(parsed_items[0].c_str());
		double sc = atof(parsed_items[1].c_str());
		m.insert( make_pair(index, sc ) );
		
	}
	inf.close();
}

void cal_summits_from_map( map<int, double> &m, map<int, double > &summits )
{
	map<int, double >::iterator ite = m.begin();
	int tsite = ite->first;
	int tstate = 1;
	double tsc = ite->second;
	++ite;
	for ( ; ite != m.end(); ++ite )
	{
		int site = ite->first;
		double sc = ite->second;
		int state = tstate;
			
		if ( sc > tsc )
		{
			state = 1;
		} else if ( sc < tsc )
		{
			state = -1;
		}
		
		if ( state == -1 && tstate == 1 )
		{
			summits.insert( make_pair( tsite, tsc ) );
		}
		
		tstate = state;
		tsite = site;
		tsc = sc;
	}
	
}

void cal_bottoms_from_map( map<int, double> &m, map<int, double > &bottoms )
{
	map<int, double >::iterator ite = m.begin();
	int tsite = ite->first;
	int tstate = 1;
	double tsc = ite->second;
	++ite;
	for ( ; ite != m.end(); ++ite )
	{
		int site = ite->first;
		double sc = ite->second;
		int state = tstate;
			
		if ( sc > tsc )
		{
			state = 1;
		} else if ( sc < tsc )
		{
			state = -1;
		}
		
		if ( state == 1 && tstate == -1 )
		{
			bottoms.insert( make_pair( tsite, tsc ) );
		}
		
		tstate = state;
		tsite = site;
		tsc = sc;
	}
	
}

void output( string outfile, map<int, double > &summits, map<int, double > &bottoms )
{
	ofstream outf( outfile.data() );
	outf<<"summits:";
	for ( map<int, double >::iterator ite = summits.begin(); ite != summits.end(); ++ite )
	{
		outf<<"\t"<<ite->first<<":"<<ite->second;
	}
	outf<<endl;
	outf<<"bottoms:";
	for ( map<int, double >::iterator ite = bottoms.begin(); ite != bottoms.end(); ++ite )
	{
		outf<<"\t"<<ite->first<<":"<<ite->second;
	}
	outf<<endl;
}

void outgwprofile( string outfile, map<int, double> &m )
{
	ofstream outf(outfile.data() );
	for ( map<int, double >::iterator ite = m.begin(); ite != m.end(); ++ite )
	{
		outf<<ite->first<<"\t"<<ite->second<<endl;
	}
	outf.close();
}

void exit_with_help()
{
	cerr <<"Get summits bottom from a profile vector"<<endl;
	cerr <<"Usage:	prog [OPTION1] [VALUE1] [[OPTION2] [VALUE2] ...]" <<endl;
	cerr <<"Options:" <<endl;
	
	cerr <<"-f		input profile file " <<endl;
	cerr <<"-w		[bool] Gaussian weight smooth (Optional default:0 NO)"<<endl;
	cerr <<"-s		sigma (for GW smooth)"<<endl;
	cerr <<"-e		extent (for GW smooth)" <<endl;
	cerr <<"-o		outfile (out summit bottom file)"<<endl;
	cerr <<"-g		output gw profile [for GW smooth only]"<<endl;
	exit(1);
}

int main( int argc, char* argv[] )
{
	
	
	string infile = "";
	int gw = 0;
	int sigma = 20;
	int extent = 73;
	string outfile = "";
	string gwoutfile = "";
	
	if ( argc == 1 )
		exit_with_help();
	
	for(int i=1; i<argc; i++)
	{
		if(argv[i][0] != '-')
			exit_with_help();

		if(argv[i][2] != '\0')
			exit_with_help();
		int option = argv[i][1];

		i++;

		switch(option)
		{
		
			
		case 's':
			sigma = atoi(argv[i]);
			break;
		case 'e':
			extent = atoi(argv[i]);
			break;
		case 'w':
			gw = atoi(argv[i]);
			break;
		case 'f':
			infile = argv[i];
			break;
		case 'o':
			outfile = argv[i];
			break;		
		case 'g':
			gwoutfile = argv[i];
			break;
		default:
			exit_with_help();
		}
	}
	
	map<int, double> m;
	readin( infile, m );
	
	if ( gw == 1)
	{
		map<int, double > w;
		Gaussian_weight( sigma, extent, w );
		GW_score( w, extent, m );
		outgwprofile( gwoutfile, m );
	}
	
	map<int, double > summits;
	cal_summits_from_map( m, summits );
	
	map<int, double > bottoms;
	cal_bottoms_from_map( m, bottoms );
	
	output( outfile, summits, bottoms );
	
	return 1;
	
}
