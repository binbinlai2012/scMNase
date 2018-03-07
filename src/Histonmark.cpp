#include "Histonmark.h"

int Histonmark::start()
{
	return pos.start;
}

int Histonmark::end()
{
	return pos.end;
}

string Histonmark::chr()
{
	return pos.chr;
}

int Histonmark::win()
{
	return pos.end-pos.start+1;
}

void Histonmark::addsignal(double insignal)
{
	signal = insignal;
}

void Histonmark::addrcount(double inr)
{
	rcount = inr;
}

int Histonmark::getlength()
{
	return pos.end - pos.start;
}

void Histonmark_GW::normalize_rc( )
{
	if ( !(libsize > 0) )
	{
		cout<<"error in normaize_rc libsize !>0 "<<libsize<<endl; exit(1);
	}
	for ( size_t i = 0; i < histonmark_Vec.size(); ++i )
	{
		histonmark_Vec[i].rcount = ( histonmark_Vec[i].rcount / ( ( (double)libsize / 1000000 ) * ( (double)histonmark_Vec[i].win() / 1000 ) ) );
	}
}

void Histonmark_GW::ranknormalize_signal()
{
	multimap<double, Region_id > rcount_region_map;
	for ( size_t i = 0; i < histonmark_Vec.size(); ++i )
	{
		rcount_region_map.insert(make_pair(histonmark_Vec[i].rcount, i ) );
	}
	sig_region_map = rcount_region_map;
	
	int totalsize = (int)rcount_region_map.size();
	if (rcount_region_map.empty() )
	{
		cout<<"error rcount_region_map empty "<<endl; exit(1);
		return;
	}
	int i = 0;
	for ( multimap<double, Region_id >::reverse_iterator ri = rcount_region_map.rbegin(); ri != rcount_region_map.rend(); ++ri )
	{
		i += 1;
		double r = (double)i / totalsize;
		histonmark_Vec[ri->second].rank = r;
	}
	
}

