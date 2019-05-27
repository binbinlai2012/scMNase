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

void rm_RD_beds( string infile, string outfile )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	
	ofstream outf( outfile.data() );
	map<string, set<pair< int, int > > > reads;
	int L = 100000;
	int i = 0;
	int rdc = 0;
	while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
		i += 1;
		if ( i % 1000000 == 0)
			cout<<i<<endl;
		vector<string> ps = parse_string( line );
		string chr = ps[0];
		int start = atoi(ps[1].c_str() );
		int end = atoi(ps[2].c_str() );
		
		if ( reads[chr].find(make_pair(start, end) ) == reads[chr].end() )
		{
			reads[chr].insert( make_pair(start, end ) );
			outf<<line<<endl;
		} else
		{
			rdc += 1;
		}
		
	}
	
	cout<<i<<" "<<rdc<<endl;
	inf.close();
	outf.close();
	
}

int main( int argc, char* argv[] )
{
	if ( argc == 1 )
	{
		cout<<"filter redundant bed"<<endl;
		cout<<"Usage: prog inbed output"<<endl;
		exit(1);
	}
	string infile = argv[1];
	string outfile = argv[2];
	rm_RD_beds( infile, outfile );
	
	return 0;
}







