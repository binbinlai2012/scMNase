#ifndef READFILE_H
#define READFILE_H

#include "operation.h"
#include "databank.h"
#include "motif.h"
#include <fstream>

using namespace std;

void readtranscriptfromucsc(string &infile, Bank &tbank );
void readtranscriptfromucsc(string &infile, Genome &genome );
void readtranscriptfromEnsmbl(string &infile, Genome &genome );

void readtranscriptfromucsc_Lite(string &infile, Bank &tbank );

void readhistonmark( string &infile, Histonmark_GW &mark_gw );

void readtranscriptexp( string &infile, TranscriptExpression_GW &exp_gw );

void readbindingsignal(string &infile, TFBinding_GW &binding_gw );

void readbindingsignal_dbdec(string &infile, TFBinding_GW &binding_gw );

void readbindingsignal_rb(string &infile, TFBinding_GW &binding_gw );

void readbindingsiteGPS( string &infile, TFBindingGPS_GW &binding_gw );

void readbindingsiteMACSsummit( string &infile, TFBindingMACS_GW &binding_gw );

void readbindingsiteMACSpeak( string &infile, TFBindingMACS_GW &binding_gw );

void readformatedfile_signal( string &infile, Bank &tbank );

void readformatedfile_gps( string &infile, Bank &tbank );

void readformatedfile_macs( string &infile, Bank &tbank );

void readformatedfile_mark( string &infile, Bank &tbank );

void readformatedfile_gene( string &infile, Bank &tbank );

void readdatatobank(Bank &bank, string &facfile, string &hisfile, string &ucscfile, string &libfile);

void readexpressionfile( string &infile, Bank &tbank );
void readexpressionfile2( string &infile, Bank &tbank );
void readexpressionfile3( string &infile, Genome &genome );
void readexpressionfile4( string &infile, Genome &genome );
void readexpressionfile_gene( string &infile, Genome &genome );

void readregionfile( string &infile, map<string, vector< pair<int, int > > > &regions );
void readregioncenter( string &infile, map<string, vector<int > > &regioncenter );
void readregioncenter_direction( string &infile, map<string, map<int, char > > &regioncenter_direction );

void readlibsizefile( string &infile, Bank &tbank );

void readPIQcall( string &infile, PIQcall_GW &piqcall_gw);

void readlincRNAs(string &infile, Bank &tbank);

void readmotifhitfromFIMO( string &infile, MotifHit_GW &motifhit_gw );

void readenhancerpeak( string &npfile, string &bpfile, Enhancer_GW &enh_gw );

void readenhancerpeak( string infile, Peak_Bank &p_bank );

void readmergedenhancerpeak( string infile, Peak_Bank &p_bank );

void readmergedenhancercount( string ctfile, string ictfile, Peak_Bank &p_bank, Sample_id sid, string type );

void readmergedenhancercount( string infile, Peak_Bank &p_bank );

void readk4me3file( string infile, Genome &genome, Sample_id sid );

void readk4me3files( string infile, Genome &genome );

void readepufile( string infile, vector<Position > &epus );

void readaligntagfile( string infile, map<string, vector< pair<int, int > > > &regions, 
	map<string, map<pair<int, int >, int > > &region_tag_map, int seg_len );
	
void readaligntagfile( string infile, map<string, vector<int > > &tagpos, int seg_len );

void readgenelist( string infile, vector<string > &genelist );

#endif

