#include "operation.h"
#include "gene.h"
#include "Histonmark.h"
#include "enhancer.h"
#include "TF.h"
#include "databank.h"
#include "readfile.h"
#include "tfbinding_dynamic.h"
#include "main_tf_dynamics.h"
#include "pet.h"

using namespace std;

void gettsscpool( vector<string > &genelist, bool T, Genome &genome, Count_pool &cpool )
{
	set<string > geneset;
	for ( size_t i = 0; i < genelist.size(); ++i )
		geneset.insert( genelist[i] );
	if ( !T )
	{
		for ( map<string, Transcript >::iterator ite = genome.name_Transcript_map.begin(); 
			ite != genome.name_Transcript_map.end(); ++ite )
		{
			if ( geneset.find( ite->second.geneName ) != geneset.end() )
			{
				string chr = ite->second.chr;
				int tss = ite->second.getTSS();
				char strand = ite->second.strand;
				string sid = ite->first;
				if ( cpool.tss_id_map[chr].find( tss ) != cpool.tss_id_map[chr].end() )
					continue;
				cpool.id_ve.push_back( sid );
				vector<int > ini_ve;
				cpool.counts_table.push_back( ini_ve );
				vector<double > ini_rve;
				cpool.rpkm_table.push_back( ini_rve );
				cpool.tss_strand_map[chr][tss] = strand;
				cpool.tss_id_map[chr][tss] = cpool.id_ve.size()-1;
			}
		}
	} else
	{
		for ( vector<string >::iterator ite = genelist.begin(); ite != genelist.end(); ++ite )
		{
			int tss = genome.name_Transcript_map[*ite].getTSS();
			char strand = genome.name_Transcript_map[*ite].strand;
			string chr = genome.name_Transcript_map[*ite].chr;
			string sid = *ite;
			if ( cpool.tss_id_map[chr].find( tss ) != cpool.tss_id_map[chr].end() )
				continue;
			cpool.id_ve.push_back( sid );
			vector<int > ini_ve;
			cpool.counts_table.push_back( ini_ve );
			vector<double > ini_rve;
			cpool.rpkm_table.push_back( ini_rve );
			cpool.tss_strand_map[chr][tss] = strand;
			cpool.tss_id_map[chr][tss] = cpool.id_ve.size()-1;
		}
	}
}

void getTFCcpool( map<string, vector<int > > &tfcenter, Count_pool &cpool )
{
	for ( map<string, vector<int > >::iterator ite = tfcenter.begin(); ite != tfcenter.end(); ++ite )
	{
		for ( vector<int >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			string sid = ite->first+"_"+inttostr(*si);
			cpool.id_ve.push_back( sid );
			vector<int > ini_ve;
			cpool.counts_table.push_back( ini_ve );
			vector<double > ini_rve;
			cpool.rpkm_table.push_back( ini_rve );
			cpool.tss_strand_map[ite->first][*si] = '+';
			cpool.tss_id_map[ite->first][*si] = cpool.id_ve.size()-1;
		}
	}
}

void getTFcpool_direction( map<string, map<int, char > > &tfcenter_direction, Count_pool &cpool )
{
	for ( map<string, map<int, char > >::iterator ite = tfcenter_direction.begin(); ite != tfcenter_direction.end(); ++ite )
	{
		for ( map<int, char >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			string sid = ite->first+"_"+inttostr(si->first);
			cpool.id_ve.push_back( sid );
			vector<int > ini_ve;
			cpool.counts_table.push_back( ini_ve );
			vector<double > ini_rve;
			cpool.rpkm_table.push_back( ini_rve );
			cpool.tss_strand_map[ite->first][si->first] = si->second;
			cpool.tss_id_map[ite->first][si->first] = cpool.id_ve.size()-1;
		}
	}
}

void getbodycpool( vector<string > &genelist, bool T, Genome &genome, Count_pool &cpool )
{
	set<string > geneset;
	for ( size_t i = 0; i < genelist.size(); ++i )
		geneset.insert( genelist[i] );

	if ( !T )
	{	
		set<string > used;
		for ( map<string, Transcript >::iterator ite = genome.name_Transcript_map.begin(); 
			ite != genome.name_Transcript_map.end(); ++ite )
		{
			if ( geneset.find( ite->second.geneName ) != geneset.end() && used.find( ite->second.geneName ) == used.end() )
			{
				string chr = ite->second.chr;
				int start = ite->second.start;
				int end = ite->second.end;
				if ( end - start + 1 < 400 )
					continue;
				char strand = ite->second.strand;
				string sid = ite->first;
				if ( cpool.body_id_map[chr].find( make_pair(start, end) ) != cpool.body_id_map[chr].end() )
					continue;
				cpool.id_ve.push_back( sid );
				vector<int > ini_ve;
				cpool.counts_table.push_back( ini_ve );
				vector<double > ini_rve;
				cpool.rpkm_table.push_back( ini_rve );
				cpool.body_strand_map[chr][make_pair(start, end)] = strand;
				cpool.body_id_map[chr][make_pair(start, end)] = cpool.id_ve.size()-1;
				used.insert( ite->second.geneName );
			}
		}
	} else
	{
		for ( vector<string >::iterator ite = genelist.begin(); ite != genelist.end(); ++ite )
		{
			int start = genome.name_Transcript_map[*ite].start;
			int end = genome.name_Transcript_map[*ite].end;
			if ( start == end )
			{
				cout<<"error start == end "<<start<<" "<<end<<endl; exit(1);
			}
			if ( end - start + 1 < 400 )
				continue;
			char strand = genome.name_Transcript_map[*ite].strand;
			string chr = genome.name_Transcript_map[*ite].chr;
			string sid = *ite;
			if ( cpool.body_id_map[chr].find( make_pair(start, end) ) != cpool.body_id_map[chr].end() )
				continue;
			cpool.id_ve.push_back( sid );
			vector<int > ini_ve;
			cpool.counts_table.push_back( ini_ve );
			vector<double > ini_rve;
			cpool.rpkm_table.push_back( ini_rve );
			cpool.body_strand_map[chr][make_pair(start, end)] = strand;
			cpool.body_id_map[chr][make_pair(start, end)] = cpool.id_ve.size()-1;
		}
	}
}

void getRegioncpool( map<string, vector<pair<int, int > > > &regions, Count_pool &cpool )
{
	for ( map<string, vector<pair<int, int > > >::iterator ite = regions.begin(); ite != regions.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			string sid = ite->first+"_"+inttostr(si->first)+"_"+inttostr(si->second);
			int start = si->first;
			int end = si->second;
			cpool.id_ve.push_back(sid);
			vector<int > ini_ve;
			cpool.counts_table.push_back( ini_ve );
			vector<double > ini_rve;
			cpool.rpkm_table.push_back( ini_rve );
			cpool.body_strand_map[chr][make_pair(start, end)] = '+';
			cpool.body_id_map[chr][make_pair(start, end)] = cpool.id_ve.size()-1;
		}
	}
}

void assigncountpool_tss( Count_pool &cpool, map<string, vector<int > > &tagpos, 
	int ups, int downs, int win, int smooth_str, int total )
{
	cout<<"getprotect_tss"<<endl;
	cpool.gettssprotectmap();
	cout<<"getregionmap_tss"<<endl;
	cpool.getregionidmap_tss( ups, downs, win );
	
	cout<<"assigncountpool"<<endl;
	assigncountpool( tagpos, cpool );
	cout<<"norm rpkm"<<endl;
	cpool.normrpkm( total );
	cout<<"smooth rpkm"<<endl;
	cpool.smooth_ave_rpkm( smooth_str );
	cout<<"reverseminusstrand"<<endl;
	cpool.reverseminusstrand();
}

void assigncountpool_body( Count_pool &cpool, map<string, vector<int > > &tagpos, 
	int ups, int downs, int step, int smooth_str, int total )
{
	cout<<"getregionmap_body"<<endl;
	cpool.getregionidmap_body( ups, downs, step );
	cout<<"assigncountpool"<<endl;
	assigncountpool( tagpos, cpool );
	cout<<"norm rpkm"<<endl;
	cpool.normrpkm( total );
	cout<<"smooth rpkm"<<endl;
	cpool.smooth_ave_rpkm( smooth_str );
	cout<<"reverseminusstrand"<<endl;
	cpool.reverseminusstrand();
}

void assigncountpool_body_sp( Count_pool &cpool, map<string, vector<int > > &tagpos, 
	int ups, int downs, int win, int step, int smooth_str, int total )
{
	cout<<"getregionmap_body_sp"<<endl;
	cpool.getregionidmap_body_sp( ups, downs, win, step );
	cout<<"assigncountpool"<<endl;
	assigncountpool( tagpos, cpool );
	cout<<"norm rpkm"<<endl;
	cpool.normrpkm( total );
	cout<<"smooth rpkm"<<endl;
	cpool.smooth_ave_rpkm( smooth_str );
	cout<<"reverseminusstrand"<<endl;
	cpool.reverseminusstrand();
}

void assigncountpool_body_sp2( Count_pool &cpool, map<string, vector<int > > &tagpos, 
	int ups, int downs, int win, int step, int smooth_str, int total )
{
	cout<<"getregionmap_body_sp2"<<endl;
	cpool.getregionidmap_body_sp2( ups, downs, win, step );
	cout<<"assigncountpool"<<endl;
	assigncountpool( tagpos, cpool );
	cout<<"norm rpkm"<<endl;
	cpool.normrpkm( total );
	cout<<"smooth rpkm"<<endl;
	cpool.smooth_ave_rpkm( smooth_str );
	cout<<"reverseminusstrand"<<endl;
	cpool.reverseminusstrand();
}

void gettagcoverposfrompet( vector< PET > &pet_ve, map<string, vector<int > > &tagpos, int win )
{
	for ( size_t i = 0; i < pet_ve.size(); ++i )
	{
		if ( pet_ve[i].chr1 != pet_ve[i].chr2 )
			continue;
		string chr = pet_ve[i].chr1;
		int c1 = pet_ve[i].start1+(pet_ve[i].end1-pet_ve[i].start1)/2;
		int c2 = pet_ve[i].start2+(pet_ve[i].end2-pet_ve[i].start2)/2;
		int index1 = c1 / win;
		int index2 = c2 / win;
		if ( index2 > index1 )
		{
			for ( int index = index1; index <= index2; ++index )
			{
				int pos = index*win+1;
				tagpos[chr].push_back( pos );
			}
		}
	}
}

void outputcount( string prefix, vector< vector<double > > &rpkm_table, vector< string > &id_ve, vector< double > &ave_rpkm, vector< vector<int > > &counts_table, int type )
{
	if ( type == 0 || type == 1 )
	{
		string tablefile = prefix + ".counttable.txt";
		ofstream outf(tablefile.data() );
		for ( size_t i = 0; i < rpkm_table.size(); ++i )
		{
			outf<<id_ve[i];
			for ( size_t j = 0; j < rpkm_table[i].size(); ++j )
			{
				outf<<"\t"<<rpkm_table[i][j];
			}
			outf<<endl;
		} 
		outf.close();
	}
	
	if ( type == 0 || type == 2 )
	{
		string ctablefile = prefix + ".rawcounttable.txt";
		ofstream outf2(ctablefile.data() );
		for ( size_t i = 0; i < counts_table.size(); ++i )
		{
			outf2<<id_ve[i];
			for ( size_t j = 0; j < counts_table[i].size(); ++j )
			{
				outf2<<"\t"<<counts_table[i][j];
			}
			outf2<<endl;
		} 
		outf2.close();
	}
	
	if ( type == 0 || type == 3 )
	{
		string vefile = prefix + ".countaveprofile.txt";
		ofstream voutf( vefile.data() );
		for ( size_t i = 0; i < ave_rpkm.size(); ++i )
		{
			voutf<<(i+1)<<"\t"<<ave_rpkm[i]<<endl;
		}
		voutf.close();
	}
}

void exit_with_help()
{
	cerr <<"read count around tss and generate matrix"<<endl;
	cerr <<"Usage:	prog [OPTION1] [VALUE1] [[OPTION2] [VALUE2] ...]" <<endl;
	cerr <<"Options:" <<endl;
	
	cerr <<"-U		Ucsc file input" <<endl;
	cerr <<"-b		tag alignment bed file (first three column: chr start end) || bedpe file" <<endl;
	cerr <<"-u		Upstream (Optional default:2000)"<<endl;
	cerr <<"-d		Downstream (Optional default:2000)"<<endl;
	cerr <<"-g		genelist file"<<endl;
	cerr <<"-p		output prefix file input" <<endl;
	cerr <<"-f		fragment length (default:150)"<<endl;
	cerr <<"-w		window unit"<<endl;
	cerr <<"-T		1/0; Indicate whether gene list are transcripts or not. (default 0)"<<endl;
	cerr <<"-t		1 tss; 2 genebody; 3 TFcenter; 4 regionbody; 5 tss_genebody; 6 genebody_up_down; 7 regionbody_2; 8 TFcenter_direction; 9 TFcenter_&&_pet_coverage (default: 1)"<<endl;
	cerr <<"-s		window steps (for genebody only, default: 100)"<<endl;
	cerr <<"-m		smooth stretch steps (default: 0)"<<endl;
	cerr <<"-R		Region file with first three columns: chr start end"<<endl;
	cerr <<"-C		TFcenter file with first two columns: chr center"<<endl;
	cerr <<"-0		Output all the results: counttable; rawcounttable; countaveprofile"<<endl;
	cerr <<"-1		Output rpkmcounttable"<<endl;
	cerr <<"-2		Output rawcounttable"<<endl;
	cerr <<"-3		Output countaveprofile"<<endl;
	
	exit(1);
}

void exit_with_help(const char error[])
{
	cerr <<"Error:	" <<error <<endl;
	exit_with_help();

	exit(1);
}

int main(int argc, char* argv[])
{
	
//	string ucscfile = "/lincRNA/annotation/mm9-ensGene.txt";
//	string ucscfile = "/Other_analysis/annotation/backup/refFlat_mm9_EntrezID_filtered.ucsc";
	string ucscfile = "";
	string bedfile = "";
	int upstream = 2000;
	int downstream = 2000;
	string prefix = "";
	int seg_len = 150;
	int win = 20;
	int smooth_str = 0;
	string genefile = "";
	bool T = false;
	int type = 1;
	int step = 100;
	string TFCfile = "";
	string Regionfile = "";
	bool rawcout = false;
	bool rpkmout = false;
	bool avgcout = false;
	bool allout = false;
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
		
		case 'g':
			genefile = argv[i];
			break;
		
		case 'U':
			ucscfile = argv[i];
			break;
		case 'u':
			upstream = atoi(argv[i]);
			break;
		case 'd':
			downstream = atoi(argv[i]);
			break;
		case 'b':
			bedfile = argv[i];
			break;
		
		case 'p':
			prefix = argv[i];
			break;
		case 'f':
			seg_len = atoi(argv[i]);
			break;
		case 'w':
			win = atoi(argv[i]);
			break;		
		case 'm':
			smooth_str = atoi(argv[i]);
			break;
		case 'T':
			if ( atoi(argv[i]) == 1 )
				T = true;
			break;
		case 't':
			type = atoi(argv[i]);
			break;
		case 's':
			step = atoi(argv[i]);
			break;
		case 'C':
			TFCfile = argv[i];
			break;
		case 'R':
			Regionfile = argv[i];
			break;
		case '0':
			allout = true;
			break;
		case '1':
			rpkmout = true;
			break;
		case '2':
			rawcout = true;
			break;
		case '3':
			avgcout = true;
			break;
		default:
			exit_with_help();
		}
	}
/*	if ( genefile.empty() )
	{
		exit_with_help("Please assign gene file");
	} */
	if ( !rpkmout && !rawcout && !avgcout )
		allout = true;
	
	
	if ( type == 1 )
	{
		cout<<"read ucscfile"<<endl;
		Genome genome;
	//	readtranscriptfromEnsmbl( ucscfile, genome );
	//	readtranscriptfromucsc( ucscfile, genome );
		genome.readtranscriptfromucsc( ucscfile );
		
		vector< string > genelist;
		readgenelist( genefile, genelist );
	
		cout<<"gettsspool"<<endl;
		Count_pool cpool;
		gettsscpool( genelist, T, genome, cpool );
	
		cout<<"read tag"<<endl;
		map<string, vector<int > > tagposmap;
		readaligntagfile( bedfile, tagposmap, seg_len );
	
		int totaltag = gettotaltagnumber( tagposmap );
	
		cout<<"assign count"<<endl;
		assigncountpool_tss( cpool, tagposmap, upstream, downstream, win, smooth_str, totaltag );
	
		cout<<"cal avepro"<<endl;
		vector< double > ave_rpkm;
		averagetablecolumn( ave_rpkm, cpool.rpkm_table );
	
		cout<<"output"<<endl;
		if ( allout )
			outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 0 );
		else
		{
			if ( rpkmout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 1 );
			if ( rawcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 2 );
			if ( avgcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 3 );
			
		}
	} else if ( type == 2 )
	{
		
		Genome genome;
	//	readtranscriptfromEnsmbl( ucscfile, genome );
	//	readtranscriptfromucsc( ucscfile, genome );
		genome.readtranscriptfromucsc( ucscfile );
		
		vector< string > genelist;
		readgenelist( genefile, genelist );

		cout<<"gettsspool"<<endl;
		Count_pool cpool;
		getbodycpool( genelist, T, genome, cpool );
		
		cout<<"read tag"<<endl;
		map<string, vector<int > > tagposmap;
		readaligntagfile( bedfile, tagposmap, seg_len );
	
		int totaltag = gettotaltagnumber( tagposmap );
		
		cout<<"assign count"<<endl;
		assigncountpool_body( cpool, tagposmap, upstream, downstream, step, smooth_str, totaltag );
	
		cout<<"cal avepro"<<endl;
		vector< double > ave_rpkm;
		averagetablecolumn( ave_rpkm, cpool.rpkm_table );
	
		cout<<"output"<<endl;
		if ( allout )
			outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 0 );
		else
		{
			if ( rpkmout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 1 );
			if ( rawcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 2 );
			if ( avgcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 3 );
			
		}
	} else if ( type == 3 )
	{
		cout<<"TFCenter"<<endl;
		cout<<"Upstream "<<upstream<<endl;
		cout<<"Downstream "<<downstream<<endl;
		cout<<"win "<<win<<endl;
		map<string, vector<int > > tfcenters;
		readregioncenter( TFCfile, tfcenters );
		
		cout<<"getTFCpool"<<endl;
		Count_pool cpool;
		getTFCcpool( tfcenters, cpool );
		
		cout<<"read tag"<<endl;
		map<string, vector<int > > tagposmap;
		readaligntagfile( bedfile, tagposmap, seg_len );
	
		int totaltag = gettotaltagnumber( tagposmap );
	
		cout<<"assign count"<<endl;
		assigncountpool_tss( cpool, tagposmap, upstream, downstream, win, smooth_str, totaltag );
	
		cout<<"cal avepro"<<endl;
		vector< double > ave_rpkm;
		averagetablecolumn( ave_rpkm, cpool.rpkm_table );
	
		cout<<"output"<<endl;
		if ( allout )
			outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 0 );
		else
		{
			if ( rpkmout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 1 );
			if ( rawcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 2 );
			if ( avgcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 3 );
			
		}
		
	} else if ( type == 4 )
	{
		cout<<"Regionbody"<<endl;
		cout<<"Upstream "<<upstream<<endl;
		cout<<"Downstream "<<downstream<<endl;
		cout<<"step "<<step<<endl;
		cout<<"smooth_str " <<smooth_str<<endl;
		
		map<string, vector< pair<int, int > > > regions;
		readregionfile(Regionfile, regions );
		
		cout<<"getTFCpool"<<endl;
		Count_pool cpool;
		getRegioncpool( regions, cpool );
	
		cout<<"read tag"<<endl;
		map<string, vector<int > > tagposmap;
		readaligntagfile( bedfile, tagposmap, seg_len );
	
		int totaltag = gettotaltagnumber( tagposmap );
		
		cout<<"assign count"<<endl;
		assigncountpool_body( cpool, tagposmap, upstream, downstream, step, smooth_str, totaltag );
	
		cout<<"cal avepro"<<endl;
		vector< double > ave_rpkm;
		averagetablecolumn( ave_rpkm, cpool.rpkm_table );
			
		cout<<"output"<<endl;
		if ( allout )
			outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 0 );
		else
		{
			if ( rpkmout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 1 );
			if ( rawcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 2 );
			if ( avgcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 3 );
			
		}
		
	} else if ( type == 5 )
	{
		Genome genome;
	//	readtranscriptfromEnsmbl( ucscfile, genome );
	//	readtranscriptfromucsc( ucscfile, genome );
		genome.readtranscriptfromucsc( ucscfile );
	
		vector< string > genelist;
		readgenelist( genefile, genelist );

		cout<<"gettsspool"<<endl;
		Count_pool cpool;
		getbodycpool( genelist, T, genome, cpool );
		
		cout<<"read tag"<<endl;
		map<string, vector<int > > tagposmap;
		readaligntagfile( bedfile, tagposmap, seg_len );
	
		int totaltag = gettotaltagnumber( tagposmap );
		
		cout<<"assign count"<<endl;
		assigncountpool_body_sp( cpool, tagposmap, upstream, downstream, win, step, smooth_str, totaltag );
	
		cout<<"cal avepro"<<endl;
		vector< double > ave_rpkm;
		averagetablecolumn( ave_rpkm, cpool.rpkm_table );
	
		cout<<"output"<<endl;
		if ( allout )
			outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 0 );
		else
		{
			if ( rpkmout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 1 );
			if ( rawcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 2 );
			if ( avgcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 3 );
			
		}
		
	} else if ( type == 6 )
	{
		Genome genome;
	//	readtranscriptfromEnsmbl( ucscfile, genome );
	//	readtranscriptfromucsc( ucscfile, genome );
		genome.readtranscriptfromucsc( ucscfile );
	
		vector< string > genelist;
		readgenelist( genefile, genelist );

		cout<<"gettsspool"<<endl;
		Count_pool cpool;
		getbodycpool( genelist, T, genome, cpool );
		
		cout<<"read tag"<<endl;
		map<string, vector<int > > tagposmap;
		readaligntagfile( bedfile, tagposmap, seg_len );
	
		int totaltag = gettotaltagnumber( tagposmap );
		
		cout<<"assign count"<<endl;
		assigncountpool_body_sp2( cpool, tagposmap, upstream, downstream, win, step, smooth_str, totaltag );
	
		cout<<"cal avepro"<<endl;
		vector< double > ave_rpkm;
		averagetablecolumn( ave_rpkm, cpool.rpkm_table );
	
		cout<<"output"<<endl;
		if ( allout )
			outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 0 );
		else
		{
			if ( rpkmout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 1 );
			if ( rawcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 2 );
			if ( avgcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 3 );
			
		}
	} else if ( type == 7 )
	{
		cout<<"Regionbody"<<endl;
		cout<<"step "<<step<<endl;
		cout<<"Upstream "<<upstream<<endl;
		cout<<"Downstream "<<downstream<<endl;
		cout<<"win "<<win<<endl;
		cout<<"smooth_str " <<smooth_str<<endl;
		
		map<string, vector< pair<int, int > > > regions;
		readregionfile(Regionfile, regions );
		
		cout<<"getTFCpool"<<endl;
		Count_pool cpool;
		getRegioncpool( regions, cpool );
		
		cout<<"read tag"<<endl;
		map<string, vector<int > > tagposmap;
		readaligntagfile( bedfile, tagposmap, seg_len );
	
		int totaltag = gettotaltagnumber( tagposmap );
		
		cout<<"assign count"<<endl;
		assigncountpool_body_sp2( cpool, tagposmap, upstream, downstream, win, step, smooth_str, totaltag );
		
		cout<<"cal avepro"<<endl;
		vector< double > ave_rpkm;
		averagetablecolumn( ave_rpkm, cpool.rpkm_table );
	
		cout<<"output"<<endl;
		if ( allout )
			outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 0 );
		else
		{
			if ( rpkmout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 1 );
			if ( rawcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 2 );
			if ( avgcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 3 );
			
		}
		
	} else if ( type == 8 )
	{
		cout<<"TFCenter"<<endl;
		cout<<"Upstream "<<upstream<<endl;
		cout<<"Downstream "<<downstream<<endl;
		cout<<"win "<<win<<endl;
		map<string, map<int, char > > tfcenter_direction;
		readregioncenter_direction( TFCfile, tfcenter_direction );
		
		cout<<"getTFCpool"<<endl;
		Count_pool cpool;
		getTFcpool_direction( tfcenter_direction, cpool );
		
		cout<<"read tag"<<endl;
		map<string, vector<int > > tagposmap;
		readaligntagfile( bedfile, tagposmap, seg_len );
	
		int totaltag = gettotaltagnumber( tagposmap );
	
		cout<<"assign count"<<endl;
		assigncountpool_tss( cpool, tagposmap, upstream, downstream, win, smooth_str, totaltag );
	
		cout<<"cal avepro"<<endl;
		vector< double > ave_rpkm;
		averagetablecolumn( ave_rpkm, cpool.rpkm_table );
	
		cout<<"output"<<endl;
		if ( allout )
			outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 0 );
		else
		{
			if ( rpkmout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 1 );
			if ( rawcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 2 );
			if ( avgcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 3 );
			
		}
		
	} else if ( type==9)
	{
		cout<<"TFCenter"<<endl;
		cout<<"Upstream "<<upstream<<endl;
		cout<<"Downstream "<<downstream<<endl;
		cout<<"win "<<win<<endl;
		map<string, vector<int > > tfcenters;
		readregioncenter( TFCfile, tfcenters );
		
		cout<<"getTFCpool"<<endl;
		Count_pool cpool;
		getTFCcpool( tfcenters, cpool );
		
		cout<<"read pet "<<endl;
		PET_bank pet_bank;
		pet_bank.readinPET( bedfile );
		
		cout<<"assign pet cover tag pos"<<endl;
		map<string, vector<int > > tagposmap;
		gettagcoverposfrompet( pet_bank.pet_ve, tagposmap, win );
		
		int totaltag = (int)pet_bank.pet_ve.size();
		
		cout<<"assign count"<<endl;
		assigncountpool_tss( cpool, tagposmap, upstream, downstream, win, smooth_str, totaltag );
	
		cout<<"cal avepro"<<endl;
		vector< double > ave_rpkm;
		averagetablecolumn( ave_rpkm, cpool.rpkm_table );
	
		cout<<"output"<<endl;
		if ( allout )
			outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 0 );
		else
		{
			if ( rpkmout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 1 );
			if ( rawcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 2 );
			if ( avgcout )
				outputcount( prefix, cpool.rpkm_table, cpool.id_ve, ave_rpkm, cpool.counts_table, 3 );
			
		}
	}
	
	return 1;
	
}



