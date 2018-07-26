#include "operation.h"
#include "gene.h"

using namespace std;

void readinscnuclpos( string infile, map<string, set<int > > &nucl_pos, int l, int u )
{
	ifstream inf(infile.data() );
	
	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	
	
	while ( !inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string> ps = parse_string(line);
		string chr = ps[0];
		if ( chr == "chrM" )
			continue;
		int start = atoi(ps[1].c_str() );
		int end = atoi(ps[2].c_str() );
		int len = end - start + 1;
		if ( len >= l && len <= u )
		{
			int c = start + len / 2;
			if ( c < 0 )
			{
				cout<<chr<<" "<<start<<" "<<end<<endl; exit(1);
			}
			nucl_pos[chr].insert( c );
		}
	}
}

void readinscnuclpos( string infile, map<string, set<int > > &nucl_pos, int l, int u, string CHR )
{
	ifstream inf(infile.data() );
	
	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	
	
	while ( !inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string> ps = parse_string(line);
		string chr = ps[0];
		if ( chr == "chrM" )
			continue;
		if ( chr != CHR )
			continue;
		int start = atoi(ps[1].c_str() );
		int end = atoi(ps[2].c_str() );
		int len = end - start + 1;
		if ( len >= l && len <= u )
		{
			int c = start + len / 2;
			if ( c < 0 )
			{
				cout<<chr<<" "<<start<<" "<<end<<endl; exit(1);
			}
			nucl_pos[chr].insert( c );
		}
	}
}

void readinscnuclpos_exclude( string infile, map<string, set<int > > &nucl_pos, int l, int u, string excludechr )
{
	ifstream inf(infile.data() );
	
	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	
	vector<string > excr = parse_string( excludechr, ',' );
	set<string > excset;
	for ( size_t i = 0; i < excr.size(); ++i )
		excset.insert( excr[i]);
	
	while ( !inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string> ps = parse_string(line);
		string chr = ps[0];
		if ( chr == "chrM" )
			continue;
		if ( excset.find( chr ) != excset.end() )
			continue;
		int start = atoi(ps[1].c_str() );
		int end = atoi(ps[2].c_str() );
		int len = end - start + 1;
		if ( len >= l && len <= u )
		{
			int c = start + len / 2;
			if ( c < 0 )
			{
				cout<<chr<<" "<<start<<" "<<end<<endl; exit(1);
			}
			nucl_pos[chr].insert( c );
		}
	}
}

void readinregion( string infile, map<string, set<pair<int, int > > > &region )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	
	string line;
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string > parseditem = parse_string( line );
		string chr = parseditem[0];
		int start = atoi( parseditem[1].c_str() );
		int end = atoi( parseditem[2].c_str() );
		region[chr].insert( make_pair(start, end) );
	}
	inf.close();
}

void readingenelist( string infile, vector<string > &genenames )
{
	ifstream inf(infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	
	
	while(!inf.eof())
	{
		getline(inf,line);
		if ( line.empty() )
			break;
		
		vector<string> ps = parse_string( line );
		string gn = ps[0];
		genenames.push_back(gn);
		
		
	}
	inf.close();
}

void getlendis_singleanchorregion( map<string, set<int > > &nucl_pos, string prefix, map<string, set<pair< int, int > > > &region, int lenlower, int lenupper )
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
			for ( int index = index1; index <= index2; ++index)
				chr_index_site[chr][index].insert( *si );
			
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
		set<int > sites_in_region;
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
					sites_in_region.insert( n1 );
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
				if ( n2 - n1 <= lenlower )
					continue;
				if ( n2 - n1 > lenupper )
					break;
				outf<<n2-n1<<endl;
			} 
		}
		
	/*	for ( set<int>::iterator si = sites_in_region.begin();
			si != sites_in_region.end(); ++si )
		{
			num += 1;
			if ( num % 10000 == 0)
				cout<<num<<endl;
			if ( num >= 1000000 )
				break;
			int n1 = *si;
			set<int>::iterator nsi = si;
			++nsi;
			for ( ; nsi != sites_in_region.end(); ++nsi )
			{
				int n2 = *nsi;
				if ( n2 - n1 <= 30 )
					continue;
				if ( n2 - n1 > 2000 )
					break;
				outf<<n2-n1<<endl;
			}
		}  */
		
	}
	outf.close();
	
//	cout<<"nucl size " <<nucl_num<<endl;
	
	
}

void getlendis_bothanchorregion( map<string, set<int > > &nucl_pos, string prefix, map<string, set<pair< int, int > > > &region, int lenlower, int lenupper )
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
			for ( int index = index1; index <= index2; ++index)
				chr_index_site[chr][index].insert( *si );
			
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
		set<int > sites_in_region;
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
					sites_in_region.insert( n1 );
					break;
				}
			}
			if ( !fd )
				continue;
			
		/*	num += 1;
			if ( num % 10000 == 0)
				cout<<num<<endl;
			if ( num >= 1000000 )
				break;
			
			set<int>::iterator nsi = si;
			++nsi;
			for ( ; nsi != ite->second.end(); ++nsi )
			{
				int n2 = *nsi;
				if ( n2 - n1 <= lenlower )
					continue;
				if ( n2 - n1 > lenupper )
					break;
				outf<<n2-n1<<endl;
			} */
		}
		
		for ( set<int>::iterator si = sites_in_region.begin();
			si != sites_in_region.end(); ++si )
		{
			num += 1;
			if ( num % 10000 == 0)
				cout<<num<<endl;
			if ( num >= 1000000 )
				break;
			int n1 = *si;
			set<int>::iterator nsi = si;
			++nsi;
			for ( ; nsi != sites_in_region.end(); ++nsi )
			{
				int n2 = *nsi;
				if ( n2 - n1 <= lenlower )
					continue;
				if ( n2 - n1 > lenupper )
					break;
				outf<<n2-n1<<endl;
			}
		}  
		
	}
	outf.close();
	
//	cout<<"nucl size " <<nucl_num<<endl;
	
	
}

void getlendis_all( map<string, set<int > > &nucl_pos, string prefix, int lenlower, int lenupper )
{
	string outfile = prefix+".link_len.txt";
	ofstream outf( outfile.data() );
	int num = 0;
	int nucl_num = 0;
	for ( map<string, set<int > >::iterator ite = nucl_pos.begin(); ite != nucl_pos.end(); ++ite )
	{
		string chr = ite->first;
		nucl_num += ite->second.size();
		
		for ( set<int>::iterator si = ite->second.begin();
			si != ite->second.end(); ++si )
		{
			num += 1;
			if ( num % 10000 == 0)
				cout<<num<<endl;
			if ( num >= 1000000 )
				break;
			int n1 = *si;
			set<int>::iterator nsi = si;
			++nsi;
			for ( ; nsi != ite->second.end(); ++nsi )
			{
				int n2 = *nsi;
				if ( n2 - n1 <= lenlower )
					continue;
				if ( n2 - n1 > lenupper )
					break;
				outf<<n2-n1<<endl;
			}
		}  
		
	}
	outf.close();
}

void getlendis_gene( Genome &genome, vector<string> &genenames, map<string, set<int > > &nucl_pos, string prefix, int up, int down, int lenlower, int lenupper )
{
	int L = 100000;
	map<string, map<int, set<int> > > chr_index_sites;
	for ( map<string, set<int > >::iterator ite = nucl_pos.begin(); ite != nucl_pos.end(); ++ite )
	{
		string chr = ite->first;
		
		for ( set<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			int index = *si / L;
			chr_index_sites[chr][index].insert( *si );
		}
	}
	
	string outfile = prefix+".link_len.txt";
	ofstream outf( outfile.data() );
	
	for ( size_t i = 0; i < genenames.size(); ++i )
	{
	//	if ( i % 1000 == 0 )
	//		cout<<"gene"<<i<<endl;
		string gene = genenames[i];
		if ( genome.name_Transcript_map.find(gene) == genome.name_Transcript_map.end() )
		{
			cout<<"warning in cal_average_sc_on_tss_geneset(): cannot find gene "<<gene<<" in genome"<<endl;
			continue;
		}
		
		
		int tss = genome.name_Transcript_map[gene].getTSS();
		string chr = genome.name_Transcript_map[gene].chr;
		char strand = genome.name_Transcript_map[gene].strand;
		int index = tss / L;
		
		if ( strand == '+' )
		{
			int start = tss - up;
			int end = tss + down;
			for ( set<int>::iterator ci = chr_index_sites[chr][index].begin(); 
				ci != chr_index_sites[chr][index].end(); ++ci )
			{
				if ( *ci < start )
					continue;
				if ( *ci > end )
					break;
					
				int n1 = *ci;
				set<int>::iterator nci = ci;
				++nci;
				for ( ; nci != chr_index_sites[chr][index].end(); ++nci )
				{
					int n2 = *nci;
					if ( n2 - n1 <= lenlower )
						continue;
					if ( n2 - n1 > lenupper )
						break;
					outf<<n2-n1<<endl;
				}
			}
		} else if ( strand == '-' )
		{
			int start = tss - down; 
			int end = tss + up;
			for ( set<int>::iterator ci = chr_index_sites[chr][index].begin(); 
				ci != chr_index_sites[chr][index].end(); ++ci )
			{
				if ( *ci < start )
					continue;
				if ( *ci > end )
					break;
				
				int n1 = *ci;
				set<int>::iterator nci = ci;
				while( nci != chr_index_sites[chr][index].begin() )
				{
					--nci;
					int n2 = *nci;
					if ( n1 - n2 <= lenlower )
						continue;
					if ( n1 - n2 > lenupper )
						break;
					outf<<n1 - n2<<endl;
				}
			}
		}
		
	}
	
	
}

void exit_with_help()
{
	cerr <<"get nucleosome-nucleosome distance from sc nucl map"<<endl;
	cerr <<"Usage:	prog [OPTION1] [VALUE1] [[OPTION2] [VALUE2] ...]" <<endl;
	cerr <<"Options:" <<endl;
	
	cerr <<"-U		Ucsc annotation file input" <<endl;
	cerr <<"-f		input bed file [nucleosome bed]" <<endl;
	cerr <<"-R		input region file"<<endl;
	cerr <<"-t		type [1:around gene TSS; 2:within regions single nucl; 3:within regions both nucl; 4:all chromosome]"<<endl;
	cerr <<"-g		genelist (refseq ID, e.g. NM_***)"<<endl;
	cerr <<"-p		outprefix"<<endl;
	cerr <<"-u		upstream of TSS"<<endl;
	cerr <<"-d		downstream of TSS"<<endl;
	cerr <<"-L		lower bound of distance considered"<<endl;
	cerr <<"-E		upper bound of distance considered"<<endl;
	cerr <<"-i		min length of fragment for nucl"<<endl;
	cerr <<"-a		max length of fregment for nucl"<<endl;
	cerr <<"-C		Chromsome only [none]"<<endl;
	cerr <<"-c		Exclude chromsomes, seperate by comer [none]"<<endl;
	exit(1);
}

int main(int argc, char* argv[] )
{
	string ucscfile = "";
	string inbedfile = "";
	
	string prefix = "out";
	
	string genelistfile = "";
	string regionfile = "";
	int minl = 140;
	int maxl = 183;
	
	int upstream = 0;
	int downstream = 500;
	int lenlower = 30;
	int lenupper = 2000;
	int type = 1;
	
	string chr = "";
	string excludechr = "";

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
		
			
		case 'U':
			ucscfile = argv[i];
			break;
		case 'f':
			inbedfile = argv[i];
			break;
		case 'i':
			minl = atoi(argv[i]);
			break;
		case 'a':
			maxl = atoi(argv[i]);
			break;
		case 'p':
			prefix = argv[i];
			break;
		case 't':
			type = atoi(argv[i]);
			break;
		case 'g':
			genelistfile = argv[i];
			break;
		case 'R':
			regionfile = argv[i];
			break;
		case 'L':
			lenlower = atoi(argv[i]);
			break;
		case 'E':
			lenupper = atoi(argv[i]);
			break;
		case 'u':
			upstream = atoi(argv[i]);
			break;
		case 'd':
			downstream = atoi(argv[i]);
			break;
		case 'C':
			chr = argv[i];
			break;
		case 'c':
			excludechr = argv[i];
			break;
		default:
			exit_with_help();
		}
	}
	
	cout<<"read in nucl pos map"<<endl;
	
	map<string, set<int > > nucl_pos;
	if ( chr != "" )
		readinscnuclpos( inbedfile, nucl_pos, minl, maxl, chr );
	else if ( excludechr != "" )
		readinscnuclpos_exclude( inbedfile, nucl_pos, minl, maxl, excludechr );
	else 
		readinscnuclpos( inbedfile, nucl_pos, minl, maxl );
	
		
	if ( type == 1 )
	{
		Genome genome;
		genome.readtranscriptfromucsc( ucscfile );
		
		cout<<"read in gene list "<<endl;
		vector<string > Ts;
		readingenelist( genelistfile, Ts );
		
		cout<<"get len dis gene set"<<endl;
		getlendis_gene( genome, Ts, nucl_pos, prefix, upstream, downstream, lenlower, lenupper ); 
		return 0;
		
	} else if ( type == 2 || type == 3 )
	{
		cout<<"read in region"<<endl;
		map<string, set<pair<int, int > > > region;
		readinregion(  regionfile, region );
		
		if ( type == 2 )
		{
			cout<<"get len dis single anchor region"<<endl;
			getlendis_singleanchorregion( nucl_pos, prefix, region, lenlower, lenupper );
		} else if ( type == 3 )
		{
			cout<<"get len dis both anchor region"<<endl;
			getlendis_bothanchorregion( nucl_pos, prefix, region, lenlower, lenupper );
	
		}
		return 0;
		
	} else if ( type == 4 )
	{
		cout<<"get len dis all"<<endl;
		getlendis_all( nucl_pos, prefix, lenlower, lenupper );
	}
		
	return 0;	
}




