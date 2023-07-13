#include <skclimatology21.h>
//#include "skosirismoderadius.h"



/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator								2005-7-21*/
/** **/
/*---------------------------------------------------------------------------*/

skModeRadius_FileLocator::skModeRadius_FileLocator()
: m_config( "USask-ARG", "/Climatology/OSIRISModeRadius/Storage/",nxRegistryConfiguration::GLOBAL_INI, true)

{
	m_startmjd = 0.0;
	m_endmjd   = 0.0;
	m_filemjd.reserve(1000);						// The file caches 3 months of ECMWF which is normally less than 400 points
	m_filenamedecoder = new skModeRadius_FileNameDecoder;
}


/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator::~skModeRadius_FileLocator		2005-7-21*/
/** **/
/*---------------------------------------------------------------------------*/

skModeRadius_FileLocator::~skModeRadius_FileLocator()
{

	if (m_filenamedecoder != NULL) delete m_filenamedecoder;
}



/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator::CheckBaseDirectory		2005-7-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skModeRadius_FileLocator::CheckBaseDirectory( )
{
	bool	ok;

	ok = !m_basedirectory.IsEmpty();
	if (!ok)
	{
		ok = m_config.LocateDirectoryFromKey(m_filenamedecoder->RegistryKeyName(), &m_basedirectory, nxTRUE, nxTRUE, "Browse for the root directory of the OSIRIS Mode Radius data files. e.g. \\\\OSIRUS\\ecmwf");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator::YearMonthCombo		2005-7-21*/
/** Get the year month combo as an integer yyyymm. yyyy = 4 digit year, mm = 2 digit month**/
/*---------------------------------------------------------------------------*/

int skModeRadius_FileLocator::Year( double usermjd )
{
	nxTimeStamp		mjd(usermjd);
	int				year,month,day;

	mjd.GetDate( &day, &month, &year );
	return year;
}


/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator::NextYearMonth		2005-7-21*/
/** Get the next year month combo (ie directory name of ECMWF storage)**/
/*---------------------------------------------------------------------------*/

int skModeRadius_FileLocator::NextYear( int year )
{
	year++;
	return year;
}


/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator::PrevYearMonth		2005-7-21*/
/** Get the previous year month combo (ie directory name of ECMWF storage)**/
/*---------------------------------------------------------------------------*/

int skModeRadius_FileLocator::PrevYear( int year )
{
	year--;
	return year;
}

/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator::EraseDirectoryListings		2005-7-21*/
/** **/
/*---------------------------------------------------------------------------*/

void skModeRadius_FileLocator::EraseDirectoryListings()
{
	m_filemjd.clear();
	m_startmjd = 0.0;
	m_endmjd   = 0.0;
}


/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator::operator()		2005-7-21*/
/** An operator overload called by nxDirectoryScan to parse filenames	**/
/*---------------------------------------------------------------------------*/

void skModeRadius_FileLocator::operator() (const char* fullfilename, bool isAdirectory )
{
	nxFileSpec		spec(fullfilename);
	double			mjd;

	if (!isAdirectory)
	{
		nxFileSpec	spec(fullfilename);
		mjd = m_filenamedecoder->DecodeMjdFromName( spec);
		if (mjd > 0.0)
		{
			m_filemjd.push_back(mjd);
		}
	}
}


/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator::InsertDirectoryListing		2005-7-21*/
/** **/
/*---------------------------------------------------------------------------*/

void skModeRadius_FileLocator::InsertDirectoryListing( int year, bool IsAscendingNode )
{
	nxString	dirname;
	nxString	subdirname;
	nxString	dirscan;

	m_filenamedecoder->EncodeYearSubdirectory( year, &subdirname) ;

	if (IsAscendingNode)
		dirname.sprintf("%s/PM/%s",(const char *)m_basedirectory,(const char*)subdirname);
	else
		dirname.sprintf("%s/AM/%s",(const char *)m_basedirectory,(const char*)subdirname);

	dirscan =  m_filenamedecoder->DirectoryScanString();
	DirectoryScan( dirname, *this, nxFALSE, dirscan.DangerousTypecast() );
}

/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator::UpdateListOfFiles		2005-7-21*/
/** Loads in the list of files available for the given year **/
/*---------------------------------------------------------------------------*/

bool skModeRadius_FileLocator::UpdateListOfFiles( double usermjd, bool IsAscendingNode )
{
	int		yearmonth;
	int		prev;
	int		next;
	bool	ok;

	ok = (usermjd > m_startmjd) &&  (usermjd < m_endmjd );
	if (!ok)
	{
		EraseDirectoryListings();
		CheckBaseDirectory();
		yearmonth = Year( usermjd );
		prev      = PrevYear ( yearmonth );
		next      = NextYear ( yearmonth );
		InsertDirectoryListing( prev, IsAscendingNode );
		InsertDirectoryListing( yearmonth, IsAscendingNode );
		InsertDirectoryListing( next, IsAscendingNode );
		std::sort( m_filemjd.begin(), m_filemjd.end() );
		ok = (m_filemjd.size() > 1);
		if (ok)
		{
			m_startmjd = m_filemjd.front();
			m_endmjd   = m_filemjd.back();
		}
		if (!ok) EraseDirectoryListings();

	}
	return ok;
}
	

/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator::LocateBoundingFileMjd		2005-7-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skModeRadius_FileLocator::LocateBoundingFileMjd( double usermjd, double* beforefilemjd, double* afterfilemjd, bool IsAscendingNode )
{
	iterator	iter;
	iterator	finish;
	iterator	start;

	bool		ok;
	
	NXASSERT((m_filenamedecoder != NULL));
	ok = UpdateListOfFiles( usermjd, IsAscendingNode );
	if (ok)
	{
		finish = m_filemjd.end();
		start  = m_filemjd.begin();
		iter = std::upper_bound( start, finish, usermjd );
		ok = ( !(iter == start) && !(iter == finish) );
		if (ok)
		{
			*afterfilemjd  = *iter--;
			*beforefilemjd = *iter;
		}
	}
	if (!ok)
	{
		nxLog::Verbose(NXLOG_WARNING, "skModeRadius_FileLocator::LocateBoundingFileMjd, Osiris Mode Radius files do not straddle/bound the requeste mjd");
		*beforefilemjd = 0;
		*afterfilemjd  = 0;
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator::GetECMWFFileName		2005-7-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skModeRadius_FileLocator::GetModeRadiusFileName( double usermjd, nxString* fullname, bool IsAscendingNode)
{
	bool	ok;

	CheckBaseDirectory();

	ok =        m_filenamedecoder->EncodeFileNamefromMjd( usermjd, m_basedirectory, fullname, IsAscendingNode );
	ok = ok &&  nxDirectory::FileExists(*fullname);
	return ok;

}

/*-----------------------------------------------------------------------------
 *					skModeRadius_FileLocator::LocateBoundingECMWFFiles		2005-7-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skModeRadius_FileLocator::LocateBoundingModeRadiusFiles( double usermjd,  nxString* beforename, nxString* aftername, bool IsAscendingNode )
{
	bool	ok;
	double	startmjd;
	double	finishmjd;

	ok = LocateBoundingFileMjd( usermjd, &startmjd, &finishmjd, IsAscendingNode );
	if (ok)
	{
		ok =       GetModeRadiusFileName( startmjd, beforename, IsAscendingNode );
		ok = ok && GetModeRadiusFileName( finishmjd, aftername, IsAscendingNode );
	}
	if (!ok)
	{
		beforename->Empty(nxFALSE);
		aftername->Empty(nxFALSE);
	}
	return ok;
}

