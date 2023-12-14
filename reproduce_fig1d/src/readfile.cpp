#include "readfile.h"

void readtranscriptfromucsc(string &infile, Bank &tbank )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}

	string line;
	getline( inf, line );
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
		char strand = parseditem[2][0];
		if ( strand != '+' && strand != '-' )
		{
			cout<<"error unrecognized strand "<<strand<<" "<<line<<endl; exit(1);
		}
		int start = atoi( parseditem[3].c_str() );
		int end = atoi( parseditem[4].c_str() );
		tbank.genome.name_Transcript_map.insert( make_pair(name, tob) );
		tbank.genome.chr_pos_TranscriptName_map[chr][make_pair(start, end)] = name;
		tbank.genome.chr_pos_TranscriptNameS_map[chr][make_pair(start, end)].push_back( name );
	}
}

void readtranscriptfromucsc(string &infile, Genome &genome )
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
		genome.name_Transcript_map.insert( make_pair(name, tob) );
		genome.chr_pos_TranscriptName_map[chr][make_pair(start, end)] = name;
		genome.chr_pos_TranscriptNameS_map[chr][make_pair(start, end)].push_back( name );
	}
}

void readtranscriptfromEnsmbl(string &infile, Genome &genome )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline(inf, line);
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		vector<string > parseditem = parse_string( line );
		if ( parseditem.size() < 16 )
		{
			cout<<"error ucsc line: "<<line<<endl;
			exit(1);
		}
		vector<string> exonStarts = parse_string( parseditem[9], ',');
		vector<string> exonEnds = parse_string( parseditem[10], ',');
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
		Transcript tob( parseditem[1], parseditem[2], parseditem[3][0], atoi(parseditem[4].c_str()), atoi(parseditem[5].c_str()),
			atoi( parseditem[6].c_str()), atoi( parseditem[7].c_str()), atoi( parseditem[8].c_str() ), t_exonStarts, t_exonEnds );
		string chr = parseditem[2];
		string name = parseditem[1];
		int start = atoi( parseditem[4].c_str() );
		int end = atoi( parseditem[5].c_str() );
		if ( tob.start == tob.end )
		{
			cout<<"error start == end "<<line<<endl; exit(1);
		}
		string Genename = parseditem[12];
		tob.geneName = Genename;
		genome.name_Transcript_map.insert( make_pair(name, tob) );
	//	genome.chr_pos_TranscriptName_map[chr][make_pair(start, end)] = name;
		genome.chr_pos_TranscriptNameS_map[chr][make_pair(start, end)].push_back( name );
		
		
	}
}

void readtranscriptfromucsc_Lite(string &infile, Bank &tbank )
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
		if ( parseditem.size() != 11 )
		{
			cout<<"error ucsc line: "<<line<<endl;
			exit(1);
		}
		
	//	Transcript tob;
		Transcript tob( parseditem[0], parseditem[1], parseditem[2][0], atoi(parseditem[3].c_str()), atoi(parseditem[4].c_str()), parseditem[10] );
		string chr = parseditem[1];
		if ( parseditem[2][0] != '+' && parseditem[2][0] != '-')
		{
			cout<<"error strand "<<line<<endl; exit(1);
		}
		string name = parseditem[0];
		int start = atoi( parseditem[2].c_str() );
		int end = atoi( parseditem[3].c_str() );
		tbank.genome.name_Transcript_map.insert( make_pair(name, tob) );
		tbank.genome.chr_pos_TranscriptName_map[chr][make_pair(start, end)] = name;
	}
}

void readhistonmark( string &infile, Histonmark_GW &mark_gw )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read "<<infile<<endl;
	string line;
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;

		vector<string > parseditem = parse_string( line );
		if ( parseditem.size() < 3 )
		{
			cout<<"error histonmark line: "<<line<<endl;
			exit(1);
		}

		string chr = parseditem[0];
		int start = atoi(parseditem[1].c_str() );
		int end = atoi(parseditem[2].c_str() );
		double signal = 0;
		double rcount = 0;
		if ( parseditem.size() == 4 )
		{
			signal = atof(parseditem[3].c_str() );
			rcount = signal;
		} else if ( parseditem.size() == 8 )
		{
			signal = atof(parseditem[6].c_str() );
			rcount = atoi(parseditem[3].c_str() )*(1-1.0/signal);
		}

		Histonmark mkb( start, end, chr );
		mkb.addsignal( signal );
		mkb.addrcount( rcount );
		mark_gw.histonmark_Vec.push_back(mkb);
		Region_id id = mark_gw.histonmark_Vec.size()-1;
		mark_gw.chr_pos_histonmark_map[chr][make_pair(start, end)] = id;
	}
	inf.close();
	//cout<<"done"<<endl;
}

void readbindingsignal(string &infile, TFBinding_GW &binding_gw )
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
		if ( parseditem.size() != 8 )
		{
			cout<<"error tfbinding line: "<<line<<endl;
			exit(1);
		}
		string chr = parseditem[0];
		int start = atoi(parseditem[1].c_str() );
		int end = atoi(parseditem[2].c_str() );
		
		double foldchg = atof(parseditem[6].c_str() );
		int count = (int)( atoi(parseditem[3].c_str() )*(1-1.0/foldchg) );
		foldchg = log(foldchg) / log(2.0);
		TFBinding tfb( start, end, chr );
		tfb.addfoldchange( foldchg );
		tfb.addrcount( count );
		binding_gw.binding_Vec.push_back(tfb);
		Region_id id = binding_gw.binding_Vec.size()-1;
		binding_gw.chr_pos_binding_map[chr][make_pair(start, end)] = id;
	}
}

void readbindingsignal_dbdec(string &infile, TFBinding_GW &binding_gw )
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
		if ( parseditem.size() != 13 )
		{
			cout<<"error ucsc line: "<<line<<endl;
			exit(1);
		}
		string chr = parseditem[0];
		int start = atoi(parseditem[1].c_str() );
		int end = atoi(parseditem[2].c_str() );
		
		double foldchg = atof(parseditem[7].c_str() );
		int count = (int)( atoi(parseditem[5].c_str() )*(1.0-foldchg) );
		foldchg = (-1)* ( log(foldchg) / log(2.0) );
		TFBinding tfb( start, end, chr );
		tfb.addfoldchange( foldchg );
		tfb.addrcount( count );
		binding_gw.binding_Vec.push_back(tfb);
		Region_id id = binding_gw.binding_Vec.size()-1;
		binding_gw.chr_pos_binding_map[chr][make_pair(start, end)] = id;
	}
	
	inf.close();
}

void readbindingsignal_rb(string &infile, TFBinding_GW &binding_gw )
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
		if ( parseditem.size() != 4 )
		{
			cout<<"error ucsc line: "<<line<<endl;
			exit(1);
		}
		string chr = parseditem[0];
		int start = atoi(parseditem[1].c_str() );
		int end = atoi(parseditem[2].c_str() );
		
		double foldchg = 0;
		int count = (int)( atoi(parseditem[3].c_str() ) );
	//	foldchg = (-1)* ( log(foldchg) / log(2.0) );
		TFBinding tfb( start, end, chr );
		tfb.addfoldchange( foldchg );
		tfb.addrcount( count );
		binding_gw.binding_Vec.push_back(tfb);
		Region_id id = binding_gw.binding_Vec.size()-1;
		binding_gw.chr_pos_binding_map[chr][make_pair(start, end)] = id;
	}
	
	inf.close();
}

void readbindingsiteGPS( string &infile, TFBindingGPS_GW &binding_gw )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}

	string line;
	getline( inf, line );
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		
		vector<string > parseditem = parse_string( line );
		string posi = parseditem[0];
		vector<string > spl = parse_string( posi, ':' );
		string chrp = spl[0];
		string chr = "chr"+chrp;
		int site =  atoi(spl[1].c_str());
		float qvalue = atof(parseditem[5].c_str());
		float pvalue = atof(parseditem[6].c_str());
		TFBindingGPS tfb(chr, site, qvalue, pvalue);
		binding_gw.binding_Vec.push_back(tfb);
		Region_id id = binding_gw.binding_Vec.size()-1;
		binding_gw.chr_pos_binding_map[chr][site] = id;
	}
	
	inf.close();
	
}

void readbindingsiteMACSsummit( string &infile, TFBindingMACS_GW &binding_gw )
{
	double qth = 1.0;
	ifstream inf( infile.data() );
	
	if ( !inf.good() )
	{
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	
	string line;
	while( !inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
			
		vector<string > parseditem = parse_string( line );
		string chr = parseditem[0];
		int summit = atoi(parseditem[1].c_str());
		string name = parseditem[3];
		float qvalue = atof(parseditem[4].c_str());
		
		if (qvalue < qth)
			continue;
		
		TFBindingMACS tfb(chr, summit, qvalue, name);
		binding_gw.binding_Vec.push_back(tfb);
		Region_id id = binding_gw.binding_Vec.size() - 1;
		binding_gw.chr_pos_binding_map[chr][summit] = id;
		binding_gw.peakname_binding_map[name] = id;
	}
	inf.close();
	
}

void readbindingsiteMACSpeak( string &infile, TFBindingMACS_GW &binding_gw )
{
	ifstream inf( infile.data() );
	
	if ( !inf.good() )
	{
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	
	string line;
	while ( !inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
			
		vector<string > parseditem = parse_string( line );
		
		int start = atoi(parseditem[1].c_str() );
		int end = atoi(parseditem[2].c_str() );
		string name = parseditem[3];
		if ( binding_gw.peakname_binding_map.find( name ) != binding_gw.peakname_binding_map.end() )
		{
			Region_id id = binding_gw.peakname_binding_map[name];
			binding_gw.binding_Vec[id].start = start;
			binding_gw.binding_Vec[id].end = end;
			binding_gw.binding_Vec[id].peak = true;
		}
	}
	inf.close();
}

void readmotifhitfromFIMO( string &infile, MotifHit_GW &motifhit_gw )
{
	ifstream inf( infile.data() );
	
	if ( !inf.good() )
	{
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	
	string line;
	getline( inf, line );
	while ( !inf.eof() )
	{
		getline(inf, line );
		if ( line.empty() )
			break;
		vector<string > parseditem = parse_string( line );
		
		string pos = parseditem[1];
		vector<string > spl = parse_string(pos, '_' );
		string chr = spl[0];
		
		int start = atoi(spl[1].c_str() );
		int end = atoi(spl[2].c_str() );
		int shifts = atoi( parseditem[2].c_str() );
		int shifte = atoi( parseditem[3].c_str() );
		int shiftm = (shifts+shifte) / 2;
		int mid = start + shiftm;
		char strand = parseditem[4][0];
		
		MotifHit moht(chr, mid, strand);
		motifhit_gw.motifhit_Vec.push_back( moht );
		Region_id id = motifhit_gw.motifhit_Vec.size()-1;
		motifhit_gw.chr_pos_motifhit_map[chr][mid] = id;
		
	
	}
	
	inf.close();
}

void readformatedfile_signal( string &infile, Bank &tbank )
{
	cout<<"read file "<<infile<<endl;
	vector<string > pit = parse_string( infile, '/');
	string sname = pit.back();
	vector<string> items = parse_string( sname, '.' );
	if ( items[0] != "tfsignal" || items.size() < 4 )
	{
		cout<<"error file name for tfsignal: "<<infile<<endl; exit(1);
	}
	string tfname = items[1];
	string cell = items[2];
/*	if ( cell != "WAT" && cell != "BAT" )
	{
		cout<<"error file name for tfsignal: cell error: "<<cell<<endl; exit(1);
	}*/
	string time_s = items[3];
	int time = atoi( time_s.substr(1).c_str() );
/*	if ( time != -3 && time != 0 && time != 2 && time != 7 )
	{
		cout<<"error file name for tfsignal: time error: "<<time<<endl; exit(1);
	}*/
	if ( tbank.cell_time_sample_map.find(cell) == tbank.cell_time_sample_map.end() )
	{
		Sample_Encode sample_encode;
		sample_encode.day = time;
		sample_encode.type = cell;
		tbank.encode_Vec.push_back( sample_encode );
		Sample_id id = tbank.encode_Vec.size()-1;
		tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
	} else
	{
		if ( tbank.cell_time_sample_map[cell].find( time ) == tbank.cell_time_sample_map[cell].end() )
		{
			Sample_Encode sample_encode;
			sample_encode.day = time;
			sample_encode.type = cell;
			tbank.encode_Vec.push_back( sample_encode );
			Sample_id id = tbank.encode_Vec.size()-1;
			tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
		}
	}
	Sample_id id = tbank.cell_time_sample_map[cell][time];
	TFBinding_GW tfb_gw;
	tbank.encode_Vec[id].name_TFBinding_GW_map.insert( make_pair(tfname, tfb_gw ) );

	if ( (tfname.find("UTX") != tfname.npos || tfname.find("MLL4") != tfname.npos || tfname.find("Brd4") != tfname.npos ) && cell == "BAT" )
		readbindingsignal_dbdec( infile, tbank.encode_Vec[id].name_TFBinding_GW_map[tfname]);
	else if ( tfname.find("FAIRE") != tfname.npos )
		readbindingsignal_rb( infile, tbank.encode_Vec[id].name_TFBinding_GW_map[tfname]);
	else
		readbindingsignal( infile, tbank.encode_Vec[id].name_TFBinding_GW_map[tfname]);
/*	if ( tfname != "MLL4" && tfname != "UTX" )
	{
		readbindingsignal( infile, tbank.encode_Vec[id].name_TFBinding_GW_map[tfname]);
	} else
	{
		
		readbindingsignal_dbdec( infile, tbank.encode_Vec[id].name_TFBinding_GW_map[tfname]);
	}*/
}

void readformatedfile_gps( string &infile, Bank &tbank )
{
	cout<<"read file "<<infile<<endl;
	vector<string > pit = parse_string( infile, '/');
	string sname = pit.back();
	vector<string> items = parse_string( sname, '.' );
	if ( items[0] != "gpstf" || items.size() < 4 )
	{
		cout<<"error file name for gpstf: "<<infile<<endl; exit(1);
	}
	string tfname = items[1];
	string cell = items[2];

	string time_s = items[3];
	int time = atoi( time_s.substr(1).c_str() );

	if ( tbank.cell_time_sample_map.find(cell) == tbank.cell_time_sample_map.end() )
	{
		Sample_Encode sample_encode;
		sample_encode.day = time;
		sample_encode.type = cell;
		tbank.encode_Vec.push_back( sample_encode );
		Sample_id id = tbank.encode_Vec.size()-1;
		tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
	} else
	{
		if ( tbank.cell_time_sample_map[cell].find( time ) == tbank.cell_time_sample_map[cell].end() )
		{
			Sample_Encode sample_encode;
			sample_encode.day = time;
			sample_encode.type = cell;
			tbank.encode_Vec.push_back( sample_encode );
			Sample_id id = tbank.encode_Vec.size()-1;
			tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
		}
	}
	Sample_id id = tbank.cell_time_sample_map[cell][time];
	TFBindingGPS_GW tfb_gw;
	tbank.encode_Vec[id].name_TFBindingGPS_GW_map.insert( make_pair(tfname, tfb_gw ) );

	readbindingsiteGPS( infile, tbank.encode_Vec[id].name_TFBindingGPS_GW_map[tfname] );
}

void readformatedfile_macs( string &infiles, Bank &tbank )
{
	vector<string > files = parse_string(infiles, ',');
	string infile = files[0];
	cout<<"read file "<<infile<<endl;
	vector<string > pit = parse_string( infile, '/');
	string sname = pit.back();
	vector<string> items = parse_string( sname, '.' );
	if ( items[0] != "macstf" || items.size() < 4 )
	{
		cout<<"error file name for macstf: "<<infile<<endl; exit(1);
	}
	string tfname = items[1];
	string cell = items[2];

	string time_s = items[3];
	int time = atoi( time_s.substr(1).c_str() );

	if ( tbank.cell_time_sample_map.find(cell) == tbank.cell_time_sample_map.end() )
	{
		Sample_Encode sample_encode;
		sample_encode.day = time;
		sample_encode.type = cell;
		tbank.encode_Vec.push_back( sample_encode );
		Sample_id id = tbank.encode_Vec.size()-1;
		tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
	} else
	{
		if ( tbank.cell_time_sample_map[cell].find( time ) == tbank.cell_time_sample_map[cell].end() )
		{
			Sample_Encode sample_encode;
			sample_encode.day = time;
			sample_encode.type = cell;
			tbank.encode_Vec.push_back( sample_encode );
			Sample_id id = tbank.encode_Vec.size()-1;
			tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
		}
	}
	Sample_id id = tbank.cell_time_sample_map[cell][time];
	TFBindingMACS_GW tfb_gw;
	tbank.encode_Vec[id].name_TFBindingMACS_GW_map.insert( make_pair(tfname, tfb_gw ) );

	readbindingsiteMACSsummit( infile, tbank.encode_Vec[id].name_TFBindingMACS_GW_map[tfname] );
	
	if ( files.size() == 2 )
	{	
		string file2 = files[1];
		readbindingsiteMACSpeak( infile, tbank.encode_Vec[id].name_TFBindingMACS_GW_map[tfname] );
	}
}

void readformatedfile_mark( string &infile, Bank &tbank )
{

	cout<<"read file "<<infile<<endl;
	vector<string > pit = parse_string( infile, '/');
	string sname = pit.back();
	vector<string> items = parse_string( sname, '.' );
	if ( items[0] != "hismark" || items.size() < 4 )
	{
		cout<<"error file name for hismark: "<<infile<<endl; exit(1);
	}
	string hmname = items[1];
	string cell = items[2];
/*	if ( cell != "WAT" && cell != "BAT"  )
	{
		cout<<"error file name for tfsignal: cell error: "<<cell<<endl; exit(1);
	}*/
	string time_s = items[3];
	int time = atoi( time_s.substr(1).c_str() );
/*	if ( time != -3 && time != 0 && time != 2 && time != 7 )
	{
		cout<<"error file name for tfsignal: time error: "<<time<<endl; exit(1);
	} */
	if ( tbank.cell_time_sample_map.find(cell) == tbank.cell_time_sample_map.end() )
	{
		Sample_Encode sample_encode;
		sample_encode.day = time;
		sample_encode.type = cell;
		tbank.encode_Vec.push_back( sample_encode );
		Sample_id id = tbank.encode_Vec.size()-1;
		tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
	} else
	{
		if ( tbank.cell_time_sample_map[cell].find( time ) == tbank.cell_time_sample_map[cell].end() )
		{
			Sample_Encode sample_encode;
			sample_encode.day = time;
			sample_encode.type = cell;
			tbank.encode_Vec.push_back( sample_encode );
			Sample_id id = tbank.encode_Vec.size()-1;
			tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
		}
	}
	Sample_id id = tbank.cell_time_sample_map[cell][time];
	Histonmark_GW hm_gw;
	tbank.encode_Vec[id].name_Histonmark_GW_map.insert( make_pair(hmname, hm_gw ) );

	readhistonmark( infile, tbank.encode_Vec[id].name_Histonmark_GW_map[hmname] );
	
}

void readtranscriptexp( string &infile, TranscriptExpression_GW &exp_gw )
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
		if ( line[0] == '#' )
			continue;

		vector<string > parseditem = parse_string( line );
		if ( parseditem.size() != 3 )
		{
			cout<<"error geneexpression line: "<<line<<endl;
			exit(1);
		}
		string name = parseditem[0];
		double rpkm = atof(parseditem[2].c_str());
		exp_gw.name_exp_map.insert( make_pair(name, rpkm ) );
	}
	inf.close();
}

void readformatedfile_gene( string &infile, Bank &tbank )
{
	cout<<"read file "<<infile<<endl;
	vector<string > pit = parse_string( infile, '/');
	string sname = pit.back();
	vector<string> items = parse_string( sname, '.' );
	if ( items[0] != "gene" || items.size() < 3 )
	{
		cout<<"error file name for geneexp: "<<infile<<endl; exit(1);
	}
	
	string cell = items[1];
/*	if ( cell != "WAT" && cell != "BAT"  )
	{
		cout<<"error file name for tfsignal: cell error: "<<cell<<endl; exit(1);
	}*/
	string time_s = items[2];
	int time = atoi( time_s.substr(1).c_str() );
/*	if ( time != -3 && time != 0 && time != 2 && time != 7 )
	{
		cout<<"error file name for tfsignal: time error: "<<time<<endl; exit(1);
	} */
	if ( tbank.cell_time_sample_map.find(cell) == tbank.cell_time_sample_map.end() )
	{
		Sample_Encode sample_encode;
		sample_encode.day = time;
		sample_encode.type = cell;
		tbank.encode_Vec.push_back( sample_encode );
		Sample_id id = tbank.encode_Vec.size()-1;
		tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
	} else
	{
		if ( tbank.cell_time_sample_map[cell].find( time ) == tbank.cell_time_sample_map[cell].end() )
		{
			Sample_Encode sample_encode;
			sample_encode.day = time;
			sample_encode.type = cell;
			tbank.encode_Vec.push_back( sample_encode );
			Sample_id id = tbank.encode_Vec.size()-1;
			tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
		}
	}
	Sample_id id = tbank.cell_time_sample_map[cell][time];
	readtranscriptexp( infile, tbank.encode_Vec[id].transcriptexpression );
}

void readdatatobank(Bank &bank, string &factorfile, string &hisfile, string &ucscfile, string &libfile )
{
	vector< string > facfiles;
	vector< string > hisfiles;
	if ( !factorfile.empty() )
	{
		ifstream inf( factorfile.data() );

	
		if( !inf.good() )
		{
			std::cout<<"file "<<factorfile<<" not found!"<<std::endl;
			exit(1);
		}

		string line;
	
		while (!inf.eof() )
	
		{
			getline(inf, line );
			if ( line.empty() )
				break;
			facfiles.push_back( line );
		}
		inf.close();
		
	}
	
	if ( !hisfile.empty() )
	{
		ifstream inf2( hisfile.data() );
		if( !inf2.good() )
		{
			std::cout<<"file "<<hisfile<<" not found!"<<std::endl;
			exit(1);
		}
		string line;
		while(!inf2.eof() )
	
		{
			getline(inf2, line );
			if ( line.empty() )
				break;
			hisfiles.push_back( line );
		}
		inf2.close();
	}
	
	for ( size_t i = 0; i < facfiles.size(); ++i )
	{
		readformatedfile_signal( facfiles[i], bank );
	}
	for ( size_t i = 0; i < hisfiles.size(); ++i )
	{
		readformatedfile_mark( hisfiles[i], bank );
	}
	if ( !libfile.empty() )
		readlibsizefile( libfile, bank );
		
	if ( !ucscfile.empty() )
	{
		readtranscriptfromucsc( ucscfile, bank );
		bank.genome.transcripttogene();
	}
		
	// normalize rank
	cout<<"normaize rank"<<endl;
	for ( size_t i = 0; i < bank.encode_Vec.size(); ++i )
	{
		for( map<string, TFBinding_GW >::iterator ite = bank.encode_Vec[i].name_TFBinding_GW_map.begin(); 
			ite != bank.encode_Vec[i].name_TFBinding_GW_map.end(); ++ite )
		{
			ite->second.ranknormalize_rc();
		}
		for ( map<string, Histonmark_GW >::iterator ite = bank.encode_Vec[i].name_Histonmark_GW_map.begin(); 
			ite != bank.encode_Vec[i].name_Histonmark_GW_map.end(); ++ite )
		{
			ite->second.ranknormalize_signal();
		}
	}
	
}

void readexpressionfile( string &infile, Bank &tbank )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}

	map< string, map< int, size_t > > sample_index_map;
	sample_index_map["BAT"][-3] = 5;
	sample_index_map["BAT"][0] = 6;
	sample_index_map["BAT"][2] = 7;
	sample_index_map["BAT"][7] = 8;
	sample_index_map["WAT"][-3] = 9;
	sample_index_map["WAT"][0] = 10;
	sample_index_map["WAT"][2] = 11;
	sample_index_map["WAT"][7] = 12;
	
	string line;
	getline(inf, line);
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;

		vector<string > parseditem = parse_string( line );
		if ( parseditem.size() != 13 )
		{
			cout<<"error expression line: "<<line<<endl;
			exit(1);
		}
		string genename = parseditem[0];
		for ( map< string, map<int, Sample_id > >::iterator ite = tbank.cell_time_sample_map.begin(); ite != tbank.cell_time_sample_map.end(); ++ite )
		{
			string cell = ite->first;
			for ( map<int, Sample_id >::iterator subi = ite->second.begin(); subi != ite->second.end(); ++subi )
			{
				int day = subi->first;
				Sample_id sample = subi->second;
				size_t index = sample_index_map[cell][day];
				tbank.encode_Vec[sample].geneexpression.name_exp_map.insert( make_pair( genename, atof(parseditem[index].c_str()) ) );
			}
		}
	}
	inf.close();
}

void readexpressionfile2( string &infile, Bank &tbank )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}

	map<size_t, pair<string, int > > index_sample_map;
	index_sample_map[5] = make_pair("BAT", -3);
	index_sample_map[6] = make_pair("BAT", 0);
	index_sample_map[7] = make_pair("BAT", 2);
	index_sample_map[8] = make_pair("BAT", 7);
	index_sample_map[9] = make_pair("WAT", -3);
	index_sample_map[10] = make_pair("WAT", 0);
	index_sample_map[11] = make_pair("WAT", 2);
	index_sample_map[12] = make_pair("WAT", 7);
	
	string line;
	getline(inf, line);
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;

		vector<string > parseditem = parse_string( line );
		if ( parseditem.size() != 13 )
		{
			cout<<"error expression line: "<<line<<endl;
			exit(1);
		}
		string genename = parseditem[0];
		for ( size_t i = 5; i < 13; ++i )
		{
			pair<string, int > sp = index_sample_map[i];
			string cell = sp.first;
			int time = sp.second;
			if ( tbank.cell_time_sample_map.find(cell) == tbank.cell_time_sample_map.end() )
			{
				Sample_Encode sample_encode;
				sample_encode.day = time;
				sample_encode.type = cell;
				tbank.encode_Vec.push_back( sample_encode );
				Sample_id id = tbank.encode_Vec.size()-1;
				tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
			} else
			{
				if ( tbank.cell_time_sample_map[cell].find( time ) == tbank.cell_time_sample_map[cell].end() )
				{
					Sample_Encode sample_encode;
					sample_encode.day = time;
					sample_encode.type = cell;
					tbank.encode_Vec.push_back( sample_encode );
					Sample_id id = tbank.encode_Vec.size()-1;
					tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
				}
			}
			
			Sample_id sample = tbank.cell_time_sample_map[cell][time];
			tbank.encode_Vec[sample].geneexpression.name_exp_map.insert( make_pair( genename, atof(parseditem[i].c_str()) ) );
				
		}
	}
	inf.close();
}

void readexpressionfile3( string &infile, Genome &genome )
{
	ifstream inf( infile.data() );
	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	
	string line;
	getline( inf, line );
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
			
		vector<string > parseditem = parse_string( line );
		
		string geneName = parseditem[0];
		string chr = parseditem[2];
		string ucscname = parseditem[3];
		vector<string > name_ve = parse_string( ucscname, ',' );
		vector<string > tss_ve = parse_string( parseditem[4], ',');
		map<Sample_id, double > sp_ex_map;
		for ( size_t i = 0; i < 8; ++i )
		{
			sp_ex_map.insert(make_pair( i, atof(parseditem[i+5].c_str())));
		}
		for ( size_t i = 0; i < name_ve.size(); ++i )
		{
			string gn = name_ve[i];
		//	int tss = atoi(tss_ve[i].c_str());
			if ( genome.name_Transcript_map.find( gn ) != genome.name_Transcript_map.end() )
			{
				genome.name_Transcript_map[gn].expression_map = sp_ex_map;
			} else
			{
				cout<<"warning "<<gn<<" not found"<<endl;
			}
		}
	}
	inf.close();
}

void readexpressionfile4( string &infile, Genome &genome )
{
	ifstream inf( infile.data() );
	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	
	string line;
	getline( inf, line );
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
			
		vector<string > parseditem = parse_string( line );
		
		string geneName = parseditem[0];
		
		map<Sample_id, double > sp_ex_map;
		for ( size_t i = 0; i < 8; ++i )
		{
			sp_ex_map.insert(make_pair( i, atof(parseditem[i+1].c_str())));
		}
		
		if ( genome.name_Transcript_map.find( geneName ) != genome.name_Transcript_map.end() )
		{
			genome.name_Transcript_map[geneName].expression_map = sp_ex_map;
		} else
		{
			cout<<"warning "<<geneName<<" not found"<<endl;
		}
		
	}
	inf.close();
}

void readexpressionfile_gene( string &infile, Genome &genome )
{
	ifstream inf( infile.data() );
	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	cout<<"read file "<<infile<<endl;
	
	string line;
	getline( inf, line );
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
			
		vector<string > parseditem = parse_string( line );
		
		string geneName = parseditem[0];
		
		map<Sample_id, double > sp_ex_map;
		for ( size_t i = 0; i < 8; ++i )
		{
			sp_ex_map.insert(make_pair( i, atof(parseditem[i+1].c_str())));
		}
		
		if ( genome.name_Gene_map.find( geneName ) != genome.name_Gene_map.end() )
		{
			genome.name_Gene_map[geneName].expression_map = sp_ex_map;
		} else
		{
			cout<<"warning "<<geneName<<" not found"<<endl; exit(1);
		}
		
	}
	inf.close();
}

void readregionfile( string &infile, map<string, vector< pair<int, int > > > &regions )
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

void readregioncenter( string &infile, map<string, vector<int > > &regioncenter )
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
		int c = atoi(parseditem[1].c_str());
		chr_pos_map[chr].insert(c);
	}
	
	inf.close();
	
	regioncenter.clear();
	for ( map< string, set<int > >::iterator ite = chr_pos_map.begin(); ite != chr_pos_map.end(); ++ite )
	{
		string chr = ite->first;
		for ( set<int>::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			regioncenter[chr].push_back(*si);
		}
	}
}

void readregioncenter_direction( string &infile, map<string, map<int, char > > &regioncenter_direction )
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
		int c = atoi(parseditem[1].c_str());
		char d = parseditem[2][0];
		regioncenter_direction[chr].insert( make_pair( c, d ));
	}
	
	inf.close();
	
	
}


void readlibsizefile( string &infile, Bank &tbank )
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
		vector<string> items = parse_string( parseditem[0], '.' );
		string tfname = items[0];
		string cell = items[1];
		int day = atoi( items[2].substr(1).c_str() );
		int libsize = atoi(parseditem[1].c_str() );
		bool fd = false;
		if ( tbank.cell_time_sample_map.find( cell ) != tbank.cell_time_sample_map.end() )
		{
			if ( tbank.cell_time_sample_map[cell].find( day ) != tbank.cell_time_sample_map[cell].end() )
			{
				Sample_id sample = tbank.cell_time_sample_map[cell][day];
				if ( tbank.encode_Vec[sample].name_TFBinding_GW_map.find( tfname ) != tbank.encode_Vec[sample].name_TFBinding_GW_map.end() )
				{
					tbank.encode_Vec[sample].name_TFBinding_GW_map[tfname].libsize = libsize;
					fd = true;
					tbank.encode_Vec[sample].name_TFBinding_GW_map[tfname].normalize_rc();
				}
				if ( tbank.encode_Vec[sample].name_Histonmark_GW_map.find( tfname ) != tbank.encode_Vec[sample].name_Histonmark_GW_map.end() )
				{
					tbank.encode_Vec[sample].name_Histonmark_GW_map[tfname].libsize = libsize;
					fd = true;
					tbank.encode_Vec[sample].name_Histonmark_GW_map[tfname].normalize_rc();
				}
			}
		}
		if ( !fd )
		{
			cout<<"warning lib "<<cell<<" "<<day<<" "<<tfname<<" not find databank"<<endl;
		}
	}
}

void readPIQcall( string &infile, PIQcall_GW &piqcall_gw )
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	
	string line;
	getline(inf, line );
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;

		vector<string > parseditem = parse_string( line, ',' );
		if ( (int)parseditem.size() != 7 )
		{
			cout<<"error PIQcall line "<<line<<endl; exit(1);
		}
		int id = atoi( parseditem[0].substr(1, parseditem[0].size()-2 ).c_str() );
		string chr = parseditem[1].substr(1, parseditem[1].size()-2 );
		int coord = atoi( parseditem[2].c_str());
		double purity = atof( parseditem[6].c_str() );
	//	cout<<id<<" "<<chr<<" "<<coord<<" "<<purity<<endl; exit(1);
		PIQcall piq(id, chr, coord, purity );
		piqcall_gw.callve.push_back( piq );
		piqcall_gw.chr_coord_region_map[chr][coord] = piqcall_gw.callve.size()-1;
		
	}
	
}


void readlincRNAs(string &infile, Bank &tbank)
{
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}
	
	map<size_t, pair<string, int > > index_sample_map;
	index_sample_map[5] = make_pair("BAT", -3);
	index_sample_map[6] = make_pair("BAT", 0);
	index_sample_map[7] = make_pair("BAT", 2);
	index_sample_map[8] = make_pair("BAT", 7);
	index_sample_map[9] = make_pair("WAT", -3);
	index_sample_map[10] = make_pair("WAT", 0);
	index_sample_map[11] = make_pair("WAT", 2);
	index_sample_map[12] = make_pair("WAT", 7);
	
	string line;
	getline( inf, line);
	while (!inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;

		vector<string > parseditem = parse_string( line );
		int id = atoi(parseditem[0].substr(14).c_str());
		string chr = parseditem[1];
		int start = atoi(parseditem[2].c_str());
		int end = atoi(parseditem[3].c_str());
		int win = atoi(parseditem[4].c_str());
		
		LincRNAs lob(id, chr, start, end, win);
		tbank.genome.id_LincRNAs_map.insert(make_pair(id, lob) );
		tbank.genome.chr_pos_LincRNAsId_map[chr][make_pair(start, end)] = id;
		
		for ( size_t i = 5; i < 13; ++i )
		{
			pair<string, int > sp = index_sample_map[i];
			string cell = sp.first;
			int time = sp.second;
			if ( tbank.cell_time_sample_map.find(cell) == tbank.cell_time_sample_map.end() )
			{
				Sample_Encode sample_encode;
				sample_encode.day = time;
				sample_encode.type = cell;
				tbank.encode_Vec.push_back( sample_encode );
				Sample_id id = tbank.encode_Vec.size()-1;
				tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
			} else
			{
				if ( tbank.cell_time_sample_map[cell].find( time ) == tbank.cell_time_sample_map[cell].end() )
				{
					Sample_Encode sample_encode;
					sample_encode.day = time;
					sample_encode.type = cell;
					tbank.encode_Vec.push_back( sample_encode );
					Sample_id id = tbank.encode_Vec.size()-1;
					tbank.cell_time_sample_map[cell].insert(make_pair(time, id ) );
				}
			}
			
			Sample_id sample = tbank.cell_time_sample_map[cell][time];
			tbank.encode_Vec[sample].lincrnaexpression.id_exp_map.insert( make_pair( id, atof(parseditem[i].c_str()) ) );
				
		}
		
	}
	
	inf.close();
}

void readenhancerpeak( string &npfile, string &bpfile, Enhancer_GW &enh_gw )
{
	ifstream npf( npfile.data() );
	
	if ( !npf.good() )
	{
		cout<<"file "<<npfile<<" not found!"<<endl; exit(1);
	}
	cout<<"read file "<<npfile<<endl;
	
	string line;
	while ( !npf.eof() )
	{
		getline( npf, line );
		if ( line.empty() )
			break;
			
		vector<string > parseditem = parse_string( line );
		string chr = parseditem[0];
		int start = atoi(parseditem[1].c_str());
		int end = atoi(parseditem[2].c_str());
		double inten = atof(parseditem.back().c_str() ) * 1000 / (end-start);
		Narrowpeak np(chr, start, end, inten );
		vector<string> sms = parse_string( parseditem[4], ',' );
		for ( size_t i = 0; i < sms.size(); ++i )
		{
			np.summits.insert( atoi(sms[i].c_str()) );
		}
		enh_gw.narrowpeak_Vec.push_back( np );
		enh_gw.chr_pos_narrowpeak_map[chr][make_pair(start, end)] = enh_gw.narrowpeak_Vec.size()-1;
		
	}
	npf.close();
	
	ifstream bpf( bpfile.data() );
	
	if ( !bpf.good() )
	{
		cout<<"file "<<bpfile<<" not found!"<<endl; exit(1);
	}
	
	cout<<"read file "<<bpfile<<endl;
	
	while( !bpf.eof() )
	{
		getline( bpf, line );
		if ( line.empty() )
			break;
			
		vector<string > parseditem = parse_string( line );
		string chr = parseditem[1];
		int start = atoi(parseditem[2].c_str() );
		int end = atoi(parseditem[3].c_str() );
		Broadpeak bp( chr, start, end );
		vector<string > sms = parse_string( parseditem[4], ',' );
		for ( size_t i = 0; i < sms.size(); ++i )
		{
			bp.summits.insert( atoi(sms[i].c_str()) );
		}
		enh_gw.broadpeak_Vec.push_back( bp );
		enh_gw.chr_pos_broadpeak_map[chr][make_pair(start, end )] = enh_gw.broadpeak_Vec.size()-1;
		
	}
	bpf.close();
	
}

void readenhancerpeak( string infile, Peak_Bank &p_bank )
{
	ifstream inf( infile.data() );
	if ( !inf.good() )
	{
		cout<<"file "<<infile<<" not found!"<<endl; exit(1);
	}
	
	string line;
	while ( !inf.eof() )
	{
		getline(inf, line );
		if ( line.empty() )
			break;
		
		vector<string > parseditem = parse_string( line );
		string type = parseditem[0];
		string cell = parseditem[1];
		int time = atoi( parseditem[2].c_str() );
		string npfile = parseditem[3];
		string bpfile = parseditem[4];
		Enhancer_GW enh_gw;
		
		readenhancerpeak( npfile, bpfile, enh_gw );
		if ( type == "k4me1")
		{
			p_bank.k4me1_GW_Vec.push_back( enh_gw );
			p_bank.cell_time_me1_idmap[cell][time] = p_bank.k4me1_GW_Vec.size()-1;
		} else if ( type == "k27ac")
		{
			p_bank.k27ac_GW_Vec.push_back( enh_gw );
			p_bank.cell_time_ac_idmap[cell][time] = p_bank.k27ac_GW_Vec.size()-1;
		} else
		{
			cout<<"error type "<<type<<endl; exit(1);
		}
	}
	inf.close();
}

void readmergedenhancerpeak( string infile, Peak_Bank &p_bank )
{
	p_bank.M_narrowpeak_Vec.clear();
	ifstream inf( infile.data() );
	if ( !inf.good() )
	{
		cout<<"file "<<infile<<" not found!"<<endl; exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline( inf, line );
	getline( inf, line );
//	int k = 0;
	while( !inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
	//	k+=1;
	//	cout<<k<<endl;
	//	cout<<line<<endl;
		vector<string > parseditem = parse_string( line );
		string chr = parseditem[0];
		int start = atoi(parseditem[1].c_str());
		int end = atoi( parseditem[2].c_str() );
		if ( parseditem.size() != 27 )
		{
			cout<<"error line "<<line<<endl;
		}
		Narrowpeak np;
		np.chr = chr;
		np.start = start;
		np.end = end;
		for ( size_t i = 0; i < 8; ++i )
		{
			
			int k = 3+(int)i*3;
			if ( k+2 >= (int)parseditem.size() )
			{
				cout<<"error line "<<line<<", "<<k<<endl; exit(1);
			}
			if ( parseditem[k] != "N" )
			{
				vector<string > sms = parse_string( parseditem[k], ',' );
				
				for ( size_t s = 0; s < sms.size(); ++s )
				{
					Summit st;
					st.summit = atoi(sms[s].c_str() );
					np.merged_summits[i]["k4me1"].push_back( st );
				} 
			}
			if ( parseditem[k+1] != "N" )
			{
				vector<string > sms = parse_string( parseditem[k+1], ',' );
				for ( size_t s = 0; s < sms.size(); ++s )
				{
					Summit st;
					st.summit = atoi(sms[s].c_str() );
					np.merged_summits[i]["k27ac"].push_back( st );
				} 
			}
			np.type_map.insert( make_pair( i, atoi(parseditem[k+2].c_str() ) ) );
		} 
		
		p_bank.M_narrowpeak_Vec.push_back( np );
	//	cout<<p_bank.M_narrowpeak_Vec.size()<<endl;
		if ( p_bank.chr_pos_M_narrowpeak_map.empty() )
		{
			map< pair<int, int >, Region_id > tm;
			p_bank.chr_pos_M_narrowpeak_map.insert( make_pair(chr, tm ) );
		} 	
		p_bank.chr_pos_M_narrowpeak_map[chr].insert( make_pair( make_pair(start, end ), p_bank.M_narrowpeak_Vec.size()-1 ) );
		
		
	}
	inf.close();
//	exit(1);
}

void readmergedenhancercount( string ctfile, string ictfile, Peak_Bank &p_bank, Sample_id sid, string type )
{
	map<string, map<pair<int, int >, double > > ctmap;
	map<string, map<pair<int, int >, double > > ictmap;
	
	ifstream infct( ctfile.data() );
	if ( !infct.good() )
	{
		cout<<"file "<<ctfile<<" not found!"<<endl; exit(1);
	}
	cout<<"read file "<<ctfile<<endl;
	string line;
	while( !infct.eof() )
	{
		getline( infct, line );
		if ( line.empty() )
		{
			break;
		}
		
		vector<string > parseditem = parse_string( line );
		string chr = parseditem[0];
		int start = atoi( parseditem[1].c_str() );
		int end = atoi( parseditem[2].c_str() );
		double ct = atof( parseditem[4].c_str() );
		ctmap[chr][make_pair(start, end )] = ct;
	}
	infct.close();
	
	ifstream infict( ictfile.data() );
	if ( !infict.good() )
	{
		cout<<"file "<<ictfile<<" not found!"<<endl; exit(1);
	}
	cout<<"read file "<<ictfile<<endl;
	while( !infict.eof() )
	{
		getline( infict, line );
		if ( line.empty() )
		{
			break;
		}
		
		vector<string > parseditem = parse_string( line );
		string chr = parseditem[0];
		int start = atoi( parseditem[1].c_str() );
		int end = atoi( parseditem[2].c_str() );
		double ct = atof( parseditem[4].c_str() );
		ictmap[chr][make_pair(start, end )] = ct;
	}
	infict.close();
	
	map<string, map<pair<int, int >, double > > dctmap;
	for ( map<string, map<pair<int, int >, double > >::iterator ite = ctmap.begin(); ite != ctmap.end(); ++ite )
	{
		for ( map<pair<int, int >, double >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			
			if ( ictmap[ite->first].find( si->first ) == ictmap[ite->first].end() )
			{
				cout<<"error ragion match "<<si->first.first<<","<<si->first.second<<","<<ite->first<<endl; exit(1);
			}
			double d = 0;
			if ( si->second > ictmap[ite->first][si->first] )
				d = si->second - ictmap[ite->first][si->first];
			dctmap[ite->first][si->first] = d;
			
		}
	} 
//	cout<<dctmap.begin()->second.begin()->first.first<<" "<<dctmap.begin()->second.begin()->first.second<<endl;
	
	for ( map< string, map< pair<int, int >, Region_id > >::iterator ite = p_bank.chr_pos_M_narrowpeak_map.begin(); 
		ite != p_bank.chr_pos_M_narrowpeak_map.end(); ++ite )
	{
		for ( map< pair<int, int >, Region_id >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			if ( dctmap[ite->first].find( si->first ) == dctmap[ite->first].end() )
			{
				cout<<"error npregion match "<<si->first.first<<","<<si->first.second<<","<<ite->first<<endl; exit(1);
			}
			p_bank.M_narrowpeak_Vec[si->second].merged_intensity[sid][type] = dctmap[ite->first][si->first];
		}
	}
}

void readmergedenhancercount( string infile, Peak_Bank &p_bank )
{
	ifstream inf( infile.data() );
	if ( !inf.good() )
	{
		cout<<"file "<<infile<<" not found!"<<endl; exit(1);
	}
	
	int i = 0;
	string line;
	while ( !inf.eof() )
	{
		getline(inf, line );
		if ( line.empty() )
			break;
		
		int s = i / 2;
		int j = i % 2;
		string type = "";
		if ( j == 0 )
		{
			type = "k4me1";
		} else
			type = "k27ac";
		
		vector< string > parseditem = parse_string( line );
		string ctfile = parseditem[0];
		string ictfile = parseditem[1];
		Region_id rid = s;
		readmergedenhancercount( ctfile, ictfile, p_bank, rid, type );
		i++; 
			
		
	}
	inf.close();
}

void readk4me3file( string infile, Genome &genome, Sample_id sid )
{
	ifstream inf( infile.data() );
	if ( !inf.good() )
	{	
		cout<<"file "<<infile<<" not found "<<endl; exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline( inf, line );
	while( !inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
		{
			break;
		}
		
		vector< string > parseditem = parse_string( line );
		string name = parseditem[0];
		string chr = parseditem[1];
		int tss = atoi( parseditem[2].c_str() );
		double dc = atof( parseditem[7].c_str() );
		if ( genome.chr_pos_actpro_map[chr].find( tss ) != genome.chr_pos_actpro_map[chr].end() )
		{
			Region_id rid = genome.chr_pos_actpro_map[chr][tss];
			genome.actpro_Vec[rid].k4me3_map.insert( make_pair(sid, dc ) );
		} else
		{
			Activepromoter ap;
			ap.chr = chr;
			ap.TSS = tss;
			ap.name = name;
			ap.k4me3_map.insert( make_pair(sid, dc ) );
			genome.actpro_Vec.push_back( ap );
			genome.chr_pos_actpro_map[chr][tss] = genome.actpro_Vec.size()-1;
		}
		
	}
	inf.close();
}

void readk4me3files( string infile, Genome &genome )
{
	ifstream inf( infile.data() );
	if ( !inf.good() )
	{	
		cout<<"file "<<infile<<" not found "<<endl; exit(1);
	}
	cout<<"read file "<<infile<<endl; 
	string line;
	vector< string > files;
	while( !inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
		{
			break;
		}
		files.push_back( line ); 
		
	}
	inf.close();
	
	for ( size_t i = 0; i < files.size(); ++i )
	{
		Sample_id sid = i;
		if ( i >= 4 )
			sid = i + 1;
		readk4me3file( files[i], genome, sid );
	}
}

void readepufile( string infile, vector<Position > &epus )
{
	ifstream inf( infile.data() );
	if ( !inf.good() )
	{
		cout<<"file "<<infile<<" not found "<<endl; exit(1);
	}
	cout<<"read file "<<infile<<endl;
	string line;
	getline( inf, line );
	while( !inf.eof() )
	{
		getline( inf, line );
		if ( line.empty() )
			break;
		
		vector< string > parseditem = parse_string( line );
		string chr = parseditem[0];
		int start = atoi( parseditem[5].c_str() );
		int end = atoi( parseditem[6].c_str() );
		Position pos( start, end, chr );
		epus.push_back( pos );
	}
	inf.close();
		
}

void readaligntagfile( string infile, map<string, vector< pair<int, int > > > &regions, 
	map<string, map<pair<int, int >, int > > &region_tag_map, int seg_len )
{
	region_tag_map.clear();
	for ( map<string, vector< pair<int, int > > >::iterator ite = regions.begin(); ite != regions.end(); ++ite )
	{
		string chr = ite->first;
		for ( vector<pair<int, int > >::iterator si = ite->second.begin(); si != ite->second.end(); ++si )
		{
			region_tag_map[chr][*si] = 0;
		}
	}
	
	ifstream inf( infile.data() );
	if ( !inf.good() )
	{
		cout<<"file "<<infile<<" not found "<<endl; exit(1);
	}
	cout<<"read file "<<infile<<endl;
	
	while( !inf.eof() )
	{
		string line;
		getline( inf, line );
		if ( line.empty() )
			break;
		
		
		vector< string > parseditem = parse_string( line );
		if ( parseditem.size() != 6 )
		{
			cout<<"error line: "<<line<<endl; exit(1);
		}
		string chr = parseditem[0];
		int start = atoi( parseditem[1].c_str() );
		int end = atoi( parseditem[2].c_str() );
		char strand = parseditem[5][0];
		if ( strand != '+' && strand != '-' )
		{
			cout<<"error strand "<<strand<<endl;
			exit(1);
		}
		int shiftpos = shift_tagpos( start, end, strand, seg_len );
		
		for ( map<pair<int, int >, int >::iterator ite = region_tag_map[chr].begin(); ite != region_tag_map[chr].end(); ++ite )
		{
			if ( ite->first.second < shiftpos )
				continue;
			if ( ite->first.first > shiftpos )
				break;
			ite->second += 1;
			break;
		}
	}
	inf.close();
	
}

void readaligntagfile( string infile, map<string, vector<int > > &tagpos, int seg_len )
{
	ifstream inf( infile.data() );
	if ( !inf.good() )
	{
		cout<<"file "<<infile<<" not found "<<endl; exit(1);
	}
	cout<<"read file "<<infile<<endl;
	
	int i = 0;
	while( !inf.eof() )
	{
		string line;
		getline( inf, line );
		if ( line.empty() )
			break;
		
		i += 1;
		if ( i % 100000 == 0 )
		{
			cout<<i<<endl;
		}
		
		vector< string > parseditem = parse_string( line );
		if ( parseditem.size() != 6 && seg_len > 0 )
		{
			cout<<"error line: "<<line<<endl; 
			cout<<"By default, you need to use BED6 format with the 6th column indicating strand."<<endl;
			cout<<"Otherwise, set -f to 0 if you use BED3 format"<<endl;
			exit(1);
		}
		string chr = parseditem[0];
		int start = atoi( parseditem[1].c_str() );
		int end = atoi( parseditem[2].c_str() );
		char strand = parseditem[5][0];
		if ( strand != '+' && strand != '-' )
		{
			cout<<"error strand "<<strand<<endl;
			exit(1);
		}
		int shiftpos = shift_tagpos( start, end, strand, seg_len );
		tagpos[chr].push_back( shiftpos );
	}
	
	inf.close();
}

void readgenelist( string infile, vector<string > &genelist )
{
	ifstream inf( infile.data() );
	if ( !inf.good() )
	{
		cout<<"file "<<infile<<" not found "<<endl; exit(1);
	}
	cout<<"read file "<<infile<<endl;
	
	while( !inf.eof() )
	{
		string line;
		getline( inf, line );
		if ( line.empty() )
			break;
			
		if ( line[0] == '#' )
			continue;
		genelist.push_back(line);
		
	}
	inf.close();
}


































