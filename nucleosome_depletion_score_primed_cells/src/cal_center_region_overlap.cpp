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

void getcenters( string infile, map<string, set<int > > &tfcenters )
{
	
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	map< string, set<int > > chr_pos_map;
	string line;
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
	//	if ( line[0] == '#' )
	//		continue;
		
		vector<string > parseditem = parse_string( line );
	/*	if ( parseditem.size() != 10 )
		{
			cout<<"error ucsc line: "<<line<<endl;
			exit(1);
		}
		if ( parseditem[0] == "REGION_ID")
			continue; */
		string chr = parseditem[0];
		
		int center = atoi(parseditem[1].c_str());
		tfcenters[chr].insert(center);
		
	}
	inf.close();
	
}

void getregion( string infile, map<string, vector< pair<int, int > > > &regions )
{
	
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	map<string, set< pair<int, int > > > chr_pos_map;
	string line;
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
	//	if ( line[0] == '#' )
	//		continue;
		
		vector<string > parseditem = parse_string( line );
	/*	if ( parseditem.size() != 10 )
		{
			cout<<"error ucsc line: "<<line<<endl;
			exit(1);
		}
		if ( parseditem[0] == "REGION_ID")
			continue; */
		string chr = parseditem[0];
		
		int start = atoi(parseditem[1].c_str() );
		int end = atoi(parseditem[2].c_str() );
		chr_pos_map[chr].insert(make_pair(start, end ) );
		
	}
	
	inf.close();
	
	regions.clear();
	for ( map<string, set< pair<int, int > > >::iterator ite = chr_pos_map.begin(); ite != chr_pos_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<pair<int, int > >::iterator subi = ite->second.begin(); subi != ite->second.end(); ++subi )
		{
			regions[chr].push_back( *subi );
		} 
	}
	
	
}

void cal_overlap( map<string, set<int > > &tfcenters, 
	map<string, vector< pair<int, int > > > &regions,
	string outovcfile,
	string outnovcfile,
	string outovrfile,
	string outnovrfile)
{
	map<string, set<int > > overlapped_centers;
	map<string, set<int > > nonoverlapped_centers;
	int ovc_num = 0;
	int novc_num = 0;
	for ( map<string, set<int > >::iterator ite = tfcenters.begin(); ite != tfcenters.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			bool overlap = false;
			for ( vector<pair<int, int > >::iterator pi = regions[chr].begin(); pi != regions[chr].end(); ++pi )
			{
				if ( pi->second < *si )
					continue;
				if ( pi->first > *si )
					break;
				overlap = true;
				break;
			}
			if ( overlap )
			{
				overlapped_centers[chr].insert( *si );
				ovc_num += 1;
			} else
			{
				nonoverlapped_centers[chr].insert( *si );
				novc_num += 1;
			}
			
		}
	}
	
	map<string, vector< pair<int, int > > > overlapped_regions;
	map<string, vector< pair<int, int > > > nonoverlapped_regions;
	int ovr_num = 0;
	int novr_num = 0;
	for ( map<string, vector< pair<int, int > > >::iterator ite = regions.begin(); ite != regions.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector< pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{	
			bool overlap = false;
			for ( set<int >::iterator pi = tfcenters[chr].begin(); pi != tfcenters[chr].end(); ++pi )
			{
				if ( *pi < si->first )
					continue;
				if ( *pi > si->second )
					break;
				overlap = true;
				break;
			}
			if ( overlap )
			{
				overlapped_regions[chr].push_back( *si );
				ovr_num += 1;
			} else
			{
				nonoverlapped_regions[chr].push_back( *si );
				novr_num += 1;
			}
		}
	}
	
	ofstream outovc( outovcfile.data() );
	for ( map<string, set<int > >::iterator ite = overlapped_centers.begin(); ite != overlapped_centers.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			outovc<<chr<<"\t"<<*si<<endl;
		}
	}
	outovc.close();
	ofstream outnovc( outnovcfile.data() );
	for ( map<string, set<int > >::iterator ite = nonoverlapped_centers.begin(); ite != nonoverlapped_centers.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			outnovc<<chr<<"\t"<<*si<<endl;
		}
	}
	outnovc.close();
	
	ofstream outovr( outovrfile.data() );
	for ( map<string, vector< pair<int, int > > >::iterator ite = overlapped_regions.begin(); ite != overlapped_regions.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector< pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{	
			outovr<<chr<<"\t"<<si->first<<"\t"<<si->second<<endl;
		}
	}
	outovr.close();
	ofstream outnovr( outnovrfile.data() );
	for ( map<string, vector< pair<int, int > > >::iterator ite = nonoverlapped_regions.begin(); ite != nonoverlapped_regions.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector< pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{	
			outnovr<<chr<<"\t"<<si->first<<"\t"<<si->second<<endl;
		}
	}
	outnovr.close();
	
	cout<<"total_centers\t"<<ovc_num+novc_num<<endl;
	cout<<"overlapped_centers\t"<<ovc_num<<endl;
	cout<<"nonoverlapped_centers\t"<<novc_num<<endl;
	cout<<"total_regions\t"<<ovr_num+novr_num<<endl;
	cout<<"overlapped_regions\t"<<ovr_num<<endl;
	cout<<"nonoverlapped_regions\t"<<novr_num<<endl;
	
}

void exit_with_help()
{
	cerr <<"calculate overlaps between centers and regions"<<endl;
	cerr <<"Usage:	prog [OPTION1] [VALUE1] [[OPTION2] [VALUE2] ...]" <<endl;
	cerr <<"Options:" <<endl;
	
	cerr <<"-C		Center file input" <<endl;
	cerr <<"-R		Region file input"<<endl;
	cerr <<"-1		output overlapped center"<<endl;
	cerr <<"-2		output nonoverlapped center"<<endl;
	cerr <<"-3		output overlapped regions" <<endl;
	cerr <<"-4		output nonoverlapped regions"<<endl;
	
	exit(1);
}

void exit_with_help(const char error[])
{
	cerr <<"Error:	" <<error <<endl;
	exit_with_help();

	exit(1);
}


int main( int argc, char* argv[] )
{
	string incenterfile;
	string inregionfile;
	string outovcfile = "overlapped_centers.txt";
	string outnovcfile = "nonoverlapped_centers.txt";
	string outovrfile = "overlapped_regions.txt";
	string outnovrfile = "nonoverlapped_regions.txt";
	
	if (argc == 1)
	{
		exit_with_help();
	}
	
	for(int i=1; i<argc; i++)
	{
		if(argv[i][0] != '-')
			exit_with_help("Options must start with \'-\'.");

		if(argv[i][2] != '\0')
			exit_with_help("The option should be exactly one letter.");
		int option = argv[i][1];

		i++;

		switch(option)
		{
		
			
		
		case 'C':
			incenterfile = argv[i];
			break;
		case 'R':
			inregionfile = argv[i];
			break;
		case '1':
			outovcfile = argv[i];
			break;
		case '2':
			outnovcfile = argv[i];
			break;
		case '3':
			outovrfile = argv[i];
			break;
		case '4':
			outnovrfile = argv[i];
			break;
		
		default:
			exit_with_help();
		}
	}
	
	if ( incenterfile.empty() )
	{
		exit_with_help("Error please assign centerfile!");
	}
	if ( inregionfile.empty() )
	{
		exit_with_help("Error please assign regionfile!");
	}
	
	map<string, set<int > > centers;
	getcenters( incenterfile, centers );
	
	map<string, vector< pair<int, int > > > regions;
	getregion( inregionfile, regions );
	
	cal_overlap( centers, regions, outovcfile, outnovcfile, outovrfile, outnovrfile);
	
	return 1;
	
}






