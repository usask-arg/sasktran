#include <skclimatology21.h>
//#include "skosirismoderadius.h"

/*-----------------------------------------------------------------------------
 *					skModeRadius_FileNameDecoder::DecodeMjdFromName		2011-7-28*/
/** osmoderadiusv1_????_??_??.nc
  * 0123456789012345678901234
  *           111111111122222
 **/
/*---------------------------------------------------------------------------*/

double skModeRadius_FileNameDecoder::DecodeMjdFromName( nxFileSpec& spec ) const
{
	nxString	name;
	int			year;
	int			month;
	int			day;
	nxTimeStamp	mjd;
	bool		ok;

	name = spec.Name();
	name.MakeLower();
	ok   =  ( name.Mid(0,15) == "osmoderadiusv1_");
	if (ok)
	{
		year  = atoi( name.Mid(15,4) );
		month = atoi( name.Mid(20,2) );
		day   = atoi( name.Mid(23,2) );
		mjd.SetToUTC( day, month, year, 0, 0, 0, 0.0);
		ok = (mjd.MJD() > 45000.0) && (mjd.MJD() < 60000);
	}
	if (!ok) mjd = 0.0;
	return mjd.MJD();
}



/*-----------------------------------------------------------------------------
 *					skModeRadius_FileNameDecoder::EncodeYearMonthSubdirectory		2011-10-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool skModeRadius_FileNameDecoder::EncodeYearSubdirectory( int year, nxString* subdirname)
{
	subdirname->sprintf("%04d/",(int)year);
	return true;
}

/*-----------------------------------------------------------------------------
 *					skModeRadius_FileNameDecoder::DirectoryScanString		2011-7-28*/
/** **/
/*---------------------------------------------------------------------------*/

const char*	skModeRadius_FileNameDecoder::DirectoryScanString() const
{
	return "osmoderadiusv1_????_??_??.nc";
}

/*-----------------------------------------------------------------------------
 *					skModeRadius_FileNameDecoder::EncodeFileNamefromMjd		2011-7-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool skModeRadius_FileNameDecoder::EncodeFileNamefromMjd( double usermjd, const char* basedirectory, nxString* fullname, bool IsAscendingNode ) const
{
	int				day,month,year, hour, mins, secs;
	double			ticks;
	nxTimeStamp		mjd(usermjd);
	

	mjd.GetUTC( &day, &month, &year, &hour, &mins, &secs, &ticks );
	if (IsAscendingNode)
		fullname->sprintf("%s/PM/%04d/osmoderadiusv1_%04d_%02d_%02d.nc",(const char*)basedirectory,(int)year, (int)year, (int)month, (int)day);
	else
		fullname->sprintf("%s/AM/%04d/osmoderadiusv1_%04d_%02d_%02d.nc",(const char*)basedirectory,(int)year, (int)year, (int)month, (int)day);

	return true;
}

