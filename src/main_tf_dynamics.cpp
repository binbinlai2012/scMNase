#include "main_tf_dynamics.h"


void main_tf_binding_dynamics( string infile, string libfile )
{

	// read in file 
	ifstream inf( infile.data() );

	
	if( !inf.good() )
	{
		
		std::cout<<"file "<<infile<<" not found!"<<std::endl;

		
		exit(1);
	
		
	}

	
	vector< string > files;
	
	string line;
	
	while (!inf.eof() )
	
	{
		getline(inf, line );
		if ( line.empty() )
			break;
		files.push_back( line );
	}

	Bank tbank;
	for ( size_t i = 0; i < files.size(); ++i )
	{
		readformatedfile_signal( files[i], tbank );
	}
	
	readlibsizefile( libfile, tbank );
	
	string cell = "BAT";
	int day1 = 0;
	int day2 = 2;
	int day3 = 7;
	Sample_id id1 = tbank.cell_time_sample_map[cell][day1];
	Sample_id id2 = tbank.cell_time_sample_map[cell][day2];
	Sample_id id3 = tbank.cell_time_sample_map[cell][day3];
	vector< pair<Sample_id, Sample_id > > stages;
	stages.push_back( make_pair(id1, id2) );
	stages.push_back( make_pair(id2, id3) );

	// get tfname
	vector< string > tfnameVec;
	for ( map<string, TFBinding_GW >::iterator ite = tbank.encode_Vec[0].name_TFBinding_GW_map.begin();
		ite != tbank.encode_Vec[0].name_TFBinding_GW_map.end(); ++ite )
	{
		tfnameVec.push_back( ite->first );
	}

	for ( size_t i = 0; i < tfnameVec.size(); ++i )
	{
		
		string tfname = tfnameVec[i];
		cout<<"tfname "<<tfname<<endl;
		Clustered_TFBinding_GW cluster_gw;
		map< int, set< Region_id > > onofftag_cluster_map;
		map< int, set< Region_id > > updowntag_cluster_map;
		cout<<"get_tfb_cluster_gw"<<endl;
		get_tfb_cluster_gw( tbank.encode_Vec, tfname, 50, cluster_gw );
		cout<<"get_tfb_dynamics_gw"<<endl;
		get_tfb_dynamics_gw( cluster_gw, stages, onofftag_cluster_map, updowntag_cluster_map );
		cout<<"output"<<endl;
		string outfile = tfname + ".D0-D2-D7.dynamics.summary";
		ofstream outf(outfile.data() );
		for ( map< int, set< Region_id > >::iterator ite = onofftag_cluster_map.begin(); ite != onofftag_cluster_map.end(); ++ite )
		{
			int n = (int)ite->second.size();
			string mean = onofftagmeaning( ite->first, 2 );
			outf<<mean<<"\t"<<n<<endl;
		}
		outf<<"-------"<<endl;
		for ( map< int, set< Region_id > >::iterator ite = updowntag_cluster_map.begin(); ite != updowntag_cluster_map.end(); ++ite )
		{
			int n = (int)ite->second.size();
			string mean = updowntagmeaning( ite->first, 2 );
			outf<<mean<<"\t"<<n<<endl;
		}
		outf<<"-------"<<endl;
	//	outf.close();
		int addrcount1 = 0;
		int addrcount2 = 0;
		int addrcount3 = 0;
		double addfdch1 = 0;
		double addfdch2 = 0;
		double addfdch3 = 0;
		int totalregion1 = 0;
		int totalregion2 = 0;
		int totalregion3 = 0;		
		string outfile2 = tfname + ".D0-D2-D7.dynamics.each";
		ofstream outf2(outfile2.data() );
		map<string, map<pair<int, int>, Region_id > >::iterator ite = cluster_gw.chr_pos_cluster_map.begin();
		for ( ; ite != cluster_gw.chr_pos_cluster_map.end(); ++ite )
		{
			string chr = ite->first;
			for ( map<pair<int, int>, Region_id >::iterator subi = ite->second.begin(); subi != ite->second.end(); ++subi )
			{

				int start = subi->first.first;
				int end = subi->first.second;
				int rcount1 = cluster_gw.clusters[subi->second].getrcount( id1 );
				int rcount2 = cluster_gw.clusters[subi->second].getrcount( id2 );
				int rcount3 = cluster_gw.clusters[subi->second].getrcount( id3 );
				double fdch1 = cluster_gw.clusters[subi->second].getfoldchange( id1 );
				double fdch2 = cluster_gw.clusters[subi->second].getfoldchange( id2 );
				double fdch3 = cluster_gw.clusters[subi->second].getfoldchange( id3 );
				if ( rcount1 > 0 )
				{
					addrcount1 += rcount1;
					addfdch1 += fdch1;
					totalregion1 += 1;
				}
				if ( rcount2 > 0 )
				{
					addrcount2 += rcount2;
					addfdch2 += fdch2;
					totalregion2 += 1;
				}
				if ( rcount3 > 0 )
				{
					addrcount3 += rcount3;
					addfdch3 += fdch3;
					totalregion3 += 1;
				}
				string onoffmean = onofftagmeaning( cluster_gw.clusters[subi->second].onoff_tag, 2 );
				string updownmean = updowntagmeaning( cluster_gw.clusters[subi->second].updown_tag, 2 );
				outf2 <<chr<<"\t"<<start<<"\t"<<end<<"\t"<<onoffmean<<"\t"<<updownmean<<"\t"<<rcount1<<"\t"<<rcount2<<"\t"<<rcount3<<"\t"<<fdch1<<"\t"<<fdch2<<"\t"<<fdch3<<endl;
			}
		}
		outf2.close();
		
		int avercount1 = addrcount1 / totalregion1;
		int avercount2 = addrcount2 / totalregion2;
		int avercount3 = addrcount3 / totalregion3;
		double avefdch1 = addfdch1 / totalregion1;
		double avefdch2 = addfdch2 / totalregion2;
		double avefdch3 = addfdch3 / totalregion3;
		outf<<"D0:\t"<<"total_rcount:"<<addrcount1<<"\ttotal_foldchange:"<<addfdch1<<"\ttotal_regions:"<<totalregion1<<"\tavercount:"<<avercount1<<"\tavefoldchange:"<<avefdch1<<endl;
		outf<<"D2:\t"<<"total_rcount:"<<addrcount2<<"\ttotal_foldchange:"<<addfdch2<<"\ttotal_regions:"<<totalregion2<<"\tavercount:"<<avercount2<<"\tavefoldchange:"<<avefdch2<<endl;
		outf<<"D7:\t"<<"total_rcount:"<<addrcount3<<"\ttotal_foldchange:"<<addfdch3<<"\ttotal_regions:"<<totalregion3<<"\tavercount:"<<avercount3<<"\tavefoldchange:"<<avefdch3<<endl;
		outf.close();
	}

}

void main_hm_binding_dynamics( string infile )
{

	// read in file 
	ifstream inf( infile.data() );

	if( !inf.good() ){
		std::cout<<"file "<<infile<<" not found!"<<std::endl;
		exit(1);
	}

	vector< string > files;
	string line;
	while (!inf.eof() )
	{
		getline(inf, line );
		if ( line.empty() )
			break;
		files.push_back( line );
	}

	Bank tbank;
	for ( size_t i = 0; i < files.size(); ++i )
	{
		readformatedfile_mark( files[i], tbank );
	}
	string cell = "BAT";
	int day1 = 0;
	int day2 = 2;
	int day3 = 7;
	Sample_id id1 = tbank.cell_time_sample_map[cell][day1];
	Sample_id id2 = tbank.cell_time_sample_map[cell][day2];
	Sample_id id3 = tbank.cell_time_sample_map[cell][day3];
	vector< pair<Sample_id, Sample_id > > stages;
	stages.push_back( make_pair(id1, id2) );
	stages.push_back( make_pair(id2, id3) );

	// get tfname
	vector< string > hmnameVec;
	for ( map<string, Histonmark_GW >::iterator ite = tbank.encode_Vec[0].name_Histonmark_GW_map.begin();
		ite != tbank.encode_Vec[0].name_Histonmark_GW_map.end(); ++ite )
	{
		hmnameVec.push_back( ite->first );
	}

	for ( size_t i = 0; i < hmnameVec.size(); ++i )
	{
		
		string hmname = hmnameVec[i];
		cout<<"hmname "<<hmname<<endl;
		Clustered_HMark_GW cluster_gw;
		map< int, set< Region_id > > onofftag_cluster_map;
		map< int, set< Region_id > > widenarrowtag_cluster_map;
		map< int, set< Region_id > > updowntag_cluster_map;
		cout<<"get_hm_cluster_gw"<<endl;
		get_hm_cluster_gw( tbank.encode_Vec, hmname, cluster_gw );
		cout<<"get_hm_dynamics_gw"<<endl;
		get_hm_dynamics_gw( cluster_gw, stages, onofftag_cluster_map, widenarrowtag_cluster_map, updowntag_cluster_map );
		cout<<"output"<<endl;
		string outfile = hmname + ".D0-D2-D7.dynamics.summary";
		ofstream outf(outfile.data() );
		for ( map< int, set< Region_id > >::iterator ite = onofftag_cluster_map.begin(); ite != onofftag_cluster_map.end(); ++ite )
		{
			int n = (int)ite->second.size();
			string mean = hmonofftagmeaning( ite->first, 2 );
			outf<<mean<<"\t"<<n<<endl;
		}
		outf<<"------"<<endl;
		for ( map< int, set< Region_id > >::iterator ite = widenarrowtag_cluster_map.begin(); ite != widenarrowtag_cluster_map.end(); ++ite )
		{
			int n = (int)ite->second.size();
			string mean = hmwidenarrowtagmeaning( ite->first, 2 );
			outf<<mean<<"\t"<<n<<endl;
		}
		outf<<"-------"<<endl;
		for ( map< int, set< Region_id > >::iterator ite = updowntag_cluster_map.begin(); ite != updowntag_cluster_map.end(); ++ite )
		{
			int n = (int)ite->second.size();
			string mean = hmupdowntagmeaning( ite->first, 2 );
			outf<<mean<<"\t"<<n<<endl;
		}
		outf<<"-------"<<endl;
	//	outf.close();

		string outfile2 = hmname + ".D0-D2-D7.dynamics.each";
		ofstream outf2(outfile2.data() );
		int addlength1 = 0;
		int addlength2 = 0;
		int addlength3 = 0;
		double addsignal1 = 0;
		double addsignal2 = 0;
		double addsignal3 = 0;
		int totalregion1 = 0;
		int totalregion2 = 0;
		int totalregion3 = 0;
		map<string, map<pair<int, int>, Region_id > >::iterator ite = cluster_gw.chr_pos_cluster_map.begin();
		for ( ; ite != cluster_gw.chr_pos_cluster_map.end(); ++ite )
		{
			string chr = ite->first;
			for ( map<pair<int, int>, Region_id >::iterator subi = ite->second.begin(); subi != ite->second.end(); ++subi )
			{

				int start = subi->first.first;
				int end = subi->first.second;
				int length1 = cluster_gw.clusters[subi->second].getlength( id1 );
				int length2 = cluster_gw.clusters[subi->second].getlength( id2 );
				int length3 = cluster_gw.clusters[subi->second].getlength( id3 );
				double signal1 = cluster_gw.clusters[subi->second].getsignal( id1 );
				double signal2 = cluster_gw.clusters[subi->second].getsignal( id2 );
				double signal3 = cluster_gw.clusters[subi->second].getsignal( id3 );
				addlength1 += length1;
				addlength2 += length2;
				addlength3 += length3;
				addsignal1 += signal1;
				addsignal2 += signal2;
				addsignal3 += signal3;
				if ( length1 != 0 )
					totalregion1 += 1;
				if ( length2 != 0 )
					totalregion2 += 1;
				if ( length3 != 0 )
					totalregion3 += 1;
				string onoffmean = hmonofftagmeaning( cluster_gw.clusters[subi->second].onoff_tag, 2 );
				string widenarrowmean = hmwidenarrowtagmeaning( cluster_gw.clusters[subi->second].widenarrow_tag, 2 );
				outf2 <<chr<<"\t"<<start<<"\t"<<end<<"\t"<<onoffmean<<"\t"<<widenarrowmean<<"\t"<<length1<<"\t"<<length2<<"\t"<<length3<<endl;
			}
		}
		outf2.close();

		int avelength1 = addlength1 / totalregion1;
		int avelength2 = addlength2 / totalregion2;
		int avelength3 = addlength3 / totalregion3;
		double avesignal1 = addsignal1 / totalregion1;
		double avesignal2 = addsignal2 / totalregion2;
		double avesignal3 = addsignal3 / totalregion3;
		outf<<"D0:\t"<<"total_length:"<<addlength1<<"\ttotal_signal:"<<addsignal1<<"\ttotal_regions:"<<totalregion1<<"\tavelength:"<<avelength1<<"\tavesignal:"<<avesignal1<<endl;
		outf<<"D2:\t"<<"total_length:"<<addlength2<<"\ttotal_signal:"<<addsignal2<<"\ttotal_regions:"<<totalregion2<<"\tavelength:"<<avelength2<<"\tavesignal:"<<avesignal2<<endl;
		outf<<"D7:\t"<<"total_length:"<<addlength3<<"\ttotal_signal:"<<addsignal3<<"\ttotal_regions:"<<totalregion3<<"\tavelength:"<<avelength3<<"\tavesignal:"<<avesignal3<<endl;;
		outf.close();
	}
}
