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

map<int, double > cal_den( vector<int > &ve, int outl )
{
	map< int, int > d_c;
	d_c[outl] = 0;
	for ( size_t i = 0; i < ve.size(); ++i )
	{
		int d = ve[i];
		if ( d < outl )
		{
			if ( d_c.find(d) == d_c.end() )
			{
				d_c[d] = 1; 
			} else
			{
				d_c[d] += 1;
			}
		} else
		{
			d_c[outl] += 1;
		}
	}
	
	int n = (int)ve.size();
	map< int, double > d_d;
	for ( map<int, int >::iterator ite = d_c.begin(); ite != d_c.end(); ++ite )
	{
		int d = ite->first;
		double den = ite->second*1.0/n;
		d_d.insert( make_pair(d, den) );
	}
	
	return d_d;
	
}

void readinve( string infile, vector<int > &ve )
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
		vector<string > ps = parse_string(line);
		int b = atof( ps[0].c_str() );
	//	cout<<line<<" "<<b<<endl;
	//	exit(1);
		ve.push_back( b );
		
	}
	inf.close();
	
}




void output( string outfile, map<int, double > &d_den )
{
	ofstream outf( outfile.data() );
	for ( map<int, double >::iterator ite = d_den.begin(); ite != d_den.end(); ++ite )
	{
		outf<<ite->first<<"\t"<<ite->second<<endl;
	}
	outf.close();
	
		
}

int main( int argc, char* argv[] )
{
	if ( argc == 1 )
	{
		cout<<"Calculate density for integ value with outlier"<<endl;
		cout<<"Usage: prog invalue[1 column] outlier[int] outfile"<<endl;
		exit(1);
		
	}
	
	string infile = argv[1];
	int outl = atoi(argv[2]);
	string outfile = argv[3];
	
	vector<int > ve;
	readinve( infile, ve );
	
	map< int, double > den = cal_den( ve, outl );
	
	output( outfile, den );
	
	return 1;
	
}



