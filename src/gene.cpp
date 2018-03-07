#include "gene.h"


void Transcript::addTranscript(string inname, string inchr, char instrand, int instart, int inend, int incdsStart, int incdsEnd, int inexonCount, vector<int> &inexonStarts,
		vector<int> &inexonEnds, string ingeneName)
{
	start = instart;
	end = inend;
	chr = inchr;
	strand = instrand;
	name = inname;
	geneName = ingeneName;
	cdsStart = incdsStart;
	cdsEnd = incdsEnd;
	exonCount = inexonCount;
	exonStarts = inexonStarts;
	exonEnds = inexonEnds;
}

void Transcript::addTranscript(string inname, string inchr, char instrand, int instart, int inend, string ingeneName)
{
	start = instart;
	end = inend;
	chr = inchr;
	strand = instrand;
	name = inname;
	geneName = ingeneName;
}

int Transcript::getTSS()
{
	if ( strand == '+' )
		return start;
	else if ( strand == '-' )
		return end;
	else
	{
		cout<<"error getTSS unrecognized strand: "<<strand<<" "<<name<<endl;
		exit(1);
	}
}

int Transcript::getTES()
{
	if ( strand == '+' )
		return end;
	else if ( strand == '-' )
		return start;
	else
	{
		cout<<"error getTES unrecognized strand: "<<strand<<" "<<name<<endl;
		exit(1);
	}
}

pair<int, int > Transcript::getPromoter( int extension )
{
	int tss = getTSS();
	return make_pair( max(1, tss-extension), tss+extension-1 );
}

pair<int, int > Transcript::getPromoter( int up, int down )
{
	int tss = getTSS();
	if ( strand == '+' )
		return make_pair( max(1, tss-up), tss+down-1 );
	else if ( strand == '-' )
		return make_pair( max(1, tss-down), tss+up-1 );
	else
	{
		cout<<"error unknown strand "<<strand<<endl; 
		exit(1);
	}
}

pair<int, int > Transcript::getgenebody( int up, int down )
{
	int tss = getTSS();
	int tes = getTES();
	if ( strand == '+' )
		return make_pair( max(1, tss-up), tes+down-1 );
	else if ( strand == '-' )
		return make_pair( max(1, tes-down), tss+up-1 );
	else
	{
		cout<<"error unknown strand "<<strand<<endl; 
		exit(1);
	}
}

pair<int, int> Transcript::getPromoter()
{
	return getPromoter( 1000 );

}

void Gene::addTranscript(std::string name)
{
	TranscriptName.insert( name );
}

string Transcript::getgenetype()
{
	return name.substr(0,2);
}

map<string, map<pair<int, int>, string > > Genome::get_chr_promoter_TranscriptName_map( int up, int down )
{
	map<string, map<pair<int, int>, string > > pt;
/*	for ( map<string, map<pair<int, int >, string > >::iterator ite = chr_pos_TranscriptName_map.begin();
		ite != chr_pos_TranscriptName_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, string >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			pair<int, int > tpromoter = name_Transcript_map[seci->second].getPromoter();
			pt[chr][tpromoter] = seci->second;
		}
	} */
	for ( map<string, map<pair<int, int >, vector<string> > >::iterator ite = chr_pos_TranscriptNameS_map.begin();
		ite != chr_pos_TranscriptNameS_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, vector<string> >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{	
			string repren = seci->second.front();
			pair<int, int > tpromoter = name_Transcript_map[repren].getPromoter( up, down );
			pt[chr][tpromoter] = repren;
		}
	}
	return pt;
}

void Genome::getintergenicregion( )
{
	cout<<"pppp"<<endl;
//	intergenicregion.clear();
	cout<<"ddd"<<endl;
	map<string, vector< pair<int, int> > > coll;
	cout<<chr_pos_TranscriptNameS_map.size()<<endl;
	for ( map<string, map<pair<int, int >, vector<string> > >::iterator ite = chr_pos_TranscriptNameS_map.begin(); ite != chr_pos_TranscriptNameS_map.end(); ++ite )
	{
	//	cout<<ite->first<<endl;
		for ( map<pair<int, int >, vector<string> >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			coll[ite->first].push_back( si->first );
		}
	}
	cout<<"1"<<endl;
	vector< map<string, vector< pair<int, int> > > > reg_ve;
	reg_ve.push_back( coll );
	map<string, vector< pair<int, int> > > merged;
	cout<<"s"<<endl;
	region_merge_naive( reg_ve, merged );
	cout<<"2"<<endl;
	cout<<merged.size()<<endl;
	map<string, vector<pair<int, int > > > intergen;
	for ( map<string, vector< pair<int, int> > >::iterator ite = merged.begin(); ite != merged.end(); ++ite )
	{
	//	cout<<"m "<<ite->first<<"\t "<<ite->second.size()<<endl;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			if ( i == ite->second.size()-1 )
				break;
			int s = ite->second[i].second+1;
			int e = ite->second[i+1].first-1;
			intergen[ite->first].push_back( make_pair(s, e) );
			if ( i < 5 )
			{
			//	cout<<ite->first<<"\t"<<s<<"\t"<<e<<endl;
			}
		}
	}
	cout<<"dd"<<endl;
	intergenicregion = intergen;
	cout<<"kk"<<endl;
//	exit(1);
}

map<string, set<pair<int, int > > > Genome::get_chr_promoter( int up, int down )
{
	map<string, set<pair<int, int > > > pt;
	for ( map<string, map<pair<int, int >, vector<string> > >::iterator ite = chr_pos_TranscriptNameS_map.begin();
		ite != chr_pos_TranscriptNameS_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, vector<string> >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{	
			string repren = seci->second.front();
			pair<int, int > tpromoter = name_Transcript_map[repren].getPromoter( up, down );
			pt[chr].insert(tpromoter);
		}
	}
	return pt;
}

map<string, set<pair<int, int > > > Genome::get_genebody( int up, int down )
{
	map<string, set<pair<int, int > > > pt;
	for ( map<string, map<pair<int, int >, vector<string> > >::iterator ite = chr_pos_TranscriptNameS_map.begin();
		ite != chr_pos_TranscriptNameS_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, vector<string> >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{	
			string repren = seci->second.front();
			pair<int, int > tpgenebody = name_Transcript_map[repren].getgenebody( up, down );
			pt[chr].insert(tpgenebody);
		}
	}
	return pt;
}

map<string, set<int > > Genome::get_chr_TSS()
{
	map<string, set<int > > pt;
	for ( map<string, map<pair<int, int >, string > >::iterator ite = chr_pos_TranscriptName_map.begin();
		ite != chr_pos_TranscriptName_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( map<pair<int, int >, string >::iterator seci = ite->second.begin(); seci != ite->second.end(); ++seci )
		{
			int tss = name_Transcript_map[seci->second].getTSS();
			pt[chr].insert(tss);
		}
	}
	return pt;
}

vector<int > Genome::get_TSS_from_Gene( string gene)
{
	vector<int > res;
	if ( name_Gene_map.find(gene) == name_Gene_map.end() )
	{	
		cout<<"error name_Gene_map not find gene "<<gene<<endl; exit(1);
	}
	set<string > transname = name_Gene_map[gene].TranscriptName;
	for ( set<string >::iterator ite = name_Gene_map[gene].TranscriptName.begin(); ite != name_Gene_map[gene].TranscriptName.end(); ++ite )
	{
		if ( name_Transcript_map.find( *ite ) ==  name_Transcript_map.end() )
		{	
			cout<<"error name_Transcript_map not find trans "<<*ite<<endl; exit(1);
		}
		int tss = name_Transcript_map[*ite].getTSS();
		res.push_back(tss);
	}
	
	return res;
}

void Genome::transcripttogene()
{
	map< string, vector<string > > gene_transcript_map;
	for ( map<string, Transcript >::iterator ite = name_Transcript_map.begin(); ite != name_Transcript_map.end(); ++ite )
	{
		string genename = ite->second.geneName;
		gene_transcript_map[genename].push_back( ite->first );
	}
	for ( map< string, vector<string > >::iterator ite = gene_transcript_map.begin(); ite != gene_transcript_map.end(); ++ite )
	{
		string genename = ite->first;
		Gene gn;
		int start = -1;
		int end = -1;
		for ( size_t i = 0; i < ite->second.size(); ++i )
		{
			string trans = ite->second[i];
			gn.addTranscript( trans );
			gn.chr = name_Transcript_map[trans].chr;
			gn.strand = name_Transcript_map[trans].strand;
			gn.type = name_Transcript_map[trans].getgenetype();
			if ( start == -1 )
			{
				start = name_Transcript_map[trans].start;
				end = name_Transcript_map[trans].end;
			} else
			{
				if ( start > name_Transcript_map[trans].start )
					start = name_Transcript_map[trans].start;
				if ( end < name_Transcript_map[trans].end )
					end = name_Transcript_map[trans].end;
			}
		}
		gn.start = start;
		gn.end = end;
		name_Gene_map.insert( make_pair( genename, gn ) );
		chr_pos_GeneName_map[gn.chr][make_pair(start, end )] = genename;
	}
	cout<<"trans gene "<<name_Gene_map.begin()->first<<endl; 
}

void Genome::addgenomeseq(string infile)
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}

	string line;
	getline( inf, line );
	if (line[0] != '>' )
	{
		cout<<"error wrong format fasta"<<endl; exit(1);
	}
	string header = line.substr(1);
	string seq = "";
	while (!inf.eof() )
	{
		getline(inf, line );
		if ( line.empty() )
			break;
		while( line[0] != '>' )
		{
			seq += line;
			getline(inf, line );
			if ( inf.eof() )
				break;
			if ( line.empty() )
				break;
		}
		if ( !seq.empty() )
		{
			genomeseq.insert( make_pair( header, seq ) );
		}
		if ( line.empty() )
			break;
		if ( line[0] == '>' )
		{
			header = line.substr(1);
			seq = "";
		} 
	}
	inf.close();
}

string Genome::getsubseq(string chr, int pos, int len )
{
	if ( genomeseq.find( chr ) == genomeseq.end() )
	{
		cout<<"error not find chr in genome: "<<chr<<endl; exit(1);
		
	} 
	if ( pos+len > (int)genomeseq[chr].size() )
	{
		cout<<"error coord exceeds chrom "<<chr<<": "<<genomeseq[chr].size()<<" "<<pos<<" "<<len<<endl; exit(1);
		
	}
	return genomeseq[chr].substr(pos, len );
}

void Genome::readtranscriptfromucsc(string &infile )
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
	/*	if ( parseditem.size() != 11 )
		{
			cout<<"error ucsc line: "<<line<<endl;
			exit(1);
		} */
		vector<string> exonStarts = parse_string( parseditem[8], ',');
		vector<string> exonEnds = parse_string( parseditem[9], ',');
		vector<int > t_exonStarts;
		vector<int > t_exonEnds;
		if ( exonStarts.size() != exonEnds.size() )
		{
			cout<<"error exonStarts size != exonEnds size: "<<line<<endl;
			exit(1);
		}
		for ( size_t i = 0; i < exonStarts.size(); ++i )
		{
			t_exonStarts.push_back( atoi(exonStarts[i].c_str() ) );
			t_exonEnds.push_back( atoi(exonEnds[i].c_str() ) );
		}
	//	Transcript tob();
		Transcript tob( parseditem[0], parseditem[1], parseditem[2][0], atoi(parseditem[3].c_str()), atoi(parseditem[4].c_str()),
			atoi( parseditem[5].c_str()), atoi( parseditem[6].c_str()), atoi( parseditem[7].c_str() ), t_exonStarts, t_exonEnds, parseditem[10] );
		string chr = parseditem[1];
		string name = parseditem[0];
		int start = atoi( parseditem[3].c_str() );
		int end = atoi( parseditem[4].c_str() );
		name_Transcript_map.insert( make_pair(name, tob) );
		chr_pos_TranscriptName_map[chr][make_pair(start, end)] = name;
		chr_pos_TranscriptNameS_map[chr][make_pair(start, end)].push_back( name );
	}
	inf.close();
}


pair<string, int> getnearestgene( string chr, int start, int end, Genome &gn)
{
	string res = "";
	int dis = -1;
	if ( gn.chr_pos_TranscriptName_map.find(chr) == gn.chr_pos_TranscriptName_map.end() )
	{
		cout<<"error chr pos transcript not find "<<chr<<endl; exit(1);
	}
	for ( map<pair<int, int >, string >::iterator ite = gn.chr_pos_TranscriptName_map[chr].begin(); ite != gn.chr_pos_TranscriptName_map[chr].end(); ++ite )
	{
		string name = ite->second;
		int tss = gn.name_Transcript_map[name].getTSS();
		int d = min( abs(tss- start), abs(tss-end) );
		if ( dis == -1 )
		{
			dis = d;
			res = gn.name_Transcript_map[name].geneName;
		} else if ( d < dis )
		{
			dis = d;
			res = gn.name_Transcript_map[name].geneName;
		}
		
		if ( tss > end )
			break; 
	}
	return make_pair( res, dis );
}







