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
#include <algorithm>
#include <stdlib.h>
#include <stdio.h>
#include <ctime>
#include <cstdlib>

using namespace std;

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

void readinSB( string infile, map<int, double > &sm, map<int, double > &bt )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline(inf, line);
	vector<string > ps = parse_string( line );
	for ( size_t i = 1; i < ps.size();  ++i )
	{
		vector<string > ps2 = parse_string( ps[i], ':');
		int pos = atoi(ps2[0].c_str() );
		double v = atof(ps2[1].c_str() );
		sm[pos] = v;
	}
	getline(inf, line);
	ps = parse_string( line );
	for ( size_t i = 1; i < ps.size();  ++i )
	{
		vector<string > ps2 = parse_string( ps[i], ':');
		int pos = atoi(ps2[0].c_str() );
		double v = atof(ps2[1].c_str() );
		bt[pos] = v;
	}
	inf.close();
}

double find_best_match( int pos, map<int, double > &m, int type )
{
	map<int, double >::iterator ite = m.begin();
	int left_pos = ite->first;
	double left_v = ite->second;
	if ( left_pos > pos )
	{
		cout<<"error match "<<pos<<endl;
		exit(1);
	}
	while ( left_pos <= pos )
	{
		left_v = ite->second;
		++ite;
		if ( ite == m.end() )
			break;
		left_pos = ite->first;
	}
	if ( ite == m.end() )
	{
		cout<<"error match 2 "<<pos<<endl;
		exit(1);
	}
	double right_v = ite->second;
	
	if ( type == 1 )
	{
		return max(left_v, right_v );
	} else 
	{
		return min(left_v, right_v );
	}
}

void get_rh( map<int, double > &gw_sm, map<int, double > &gw_bt,
	map<int, double > &sm, map<int, double > &bt,
	string outfile )
{
	ofstream outf( outfile.data() );
	for ( map<int, double >::iterator ite = gw_sm.begin(); ite != gw_sm.end(); ++ite )
	{
		int smpos = ite->first;
		
		map<int, double >::iterator bi = gw_bt.begin();
		int left_bt = bi->first;
	//	cout<<left_bt<<endl;
		if ( left_bt >= smpos )
			continue;
		int f_left_bt = left_bt;
		while ( left_bt < smpos )
		{
			f_left_bt = left_bt;
			++bi;
			if ( bi == gw_bt.end() )
				break;
			left_bt = bi->first;
		}
		
		if ( bi == gw_bt.end() )
			continue;
		int f_right_bt = bi->first;
	//	cout<<smpos<<" "<<f_left_bt<<" "<<f_right_bt<<endl;
		double left_bt_v = find_best_match(f_left_bt, bt, 0);
		double right_bt_v = find_best_match(f_right_bt, bt, 0);
		double top_v = find_best_match(smpos, sm, 1 );
		double rh = top_v - (left_bt_v+right_bt_v) / 2;
		
		outf<<smpos<<":"<<top_v<<"\t"<<f_left_bt<<":"<<left_bt_v<<"\t"<<f_right_bt<<":"<<right_bt_v<<"\t"<<rh<<endl;
	}
	outf.close();
}

int main(int argc, char* argv[])
{
	cout<<"calculate relative height"<<endl;
	if ( argc == 1 )
	{
		cout<<"Usage: prog sm_gw_SBfile sm_SBfile outfile"<<endl;
		exit(1);
	}
	
	string smgwSBfile = argv[1];
	string smSBfile = argv[2];
	string outfile = argv[3];
	
	map<int, double > gw_sm;
	map<int, double > gw_bt;
	readinSB( smgwSBfile, gw_sm, gw_bt );
	
	map<int, double > sm;
	map<int, double > bt;
	readinSB( smSBfile, sm, bt );
	
	get_rh( gw_sm, gw_bt, sm, bt, outfile );
	
	return 0;
	
}




