#include "TF.h"

int TFBinding::start()
{
	return pos.start;
}

int TFBinding::end()
{
	return pos.end;
}

string TFBinding::chr()
{
	return pos.chr;
}

int TFBinding::win()
{
	return pos.end-pos.start+1;
}

void TFBinding::addsignal(int insignal)
{
	signal = insignal;
}

void TFBinding::addrcount( int inrc )
{
	rcount = inrc;
}

void TFBinding::addfoldchange( double fdcg )
{
	foldchange = fdcg;
}

void TFBinding_GW::normalize_rc( )
{
	if ( !(libsize > 0) )
	{
		cout<<"error in normaize_rc libsize !>0 "<<libsize<<endl; exit(1);
	}
	for ( size_t i = 0; i < binding_Vec.size(); ++i )
	{
		binding_Vec[i].rcount = (int)( (double)binding_Vec[i].rcount / ( ( (double)libsize / 1000000 ) * ( (double)binding_Vec[i].win() / 1000 ) ) );
	}
}

void TFBinding_GW::ranknormalize_rc()
{
	
	for ( size_t i = 0; i < binding_Vec.size(); ++i )
	{
		rc_region_map.insert(make_pair( binding_Vec[i].rcount, i ) );
	}
	
	int totalsize = (int)rc_region_map.size();
	if (rc_region_map.empty() )
	{
		cout<<"error rc_region_map empty "<<endl; exit(1);
		return;
	}
	int i = 0;
	for ( multimap<int, Region_id >::reverse_iterator ri = rc_region_map.rbegin(); ri != rc_region_map.rend(); ++ri )
	{
		i += 1;
		double r = (double)i / totalsize;
		binding_Vec[ri->second].rank = r;
	}
	
}

