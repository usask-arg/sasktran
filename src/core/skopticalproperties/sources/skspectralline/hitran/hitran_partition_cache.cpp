#include <skopticalproperties21.h>

/*---------------------------------------------------------------------------
 *       HitranPartitioTableCache::HitranPartitioTableCache       2019-11-08 */
/** **/
/*---------------------------------------------------------------------------*/

HitranPartitionTableCache::HitranPartitionTableCache( const skHitranPartitionTableEntry*	parent)
{
	m_parent = parent;
	m_Tstart = -9999.0;
	m_Tend   = -9999.0;
	m_deltaT = 0.015625;
	m_versionid = 0x12345678;
}


/*---------------------------------------------------------------------------
 *  HitranPartitionTableCache::LoadBaseDirectoryNameFromRegistry  2019-11-08 */
/** **/
/*---------------------------------------------------------------------------*/

bool HitranPartitionTableCache::LoadBaseDirectoryNameFromRegistry( nxString* filename )
{
	nxRegistryConfiguration		config( "USask-ARG", "skOpticalProperties/Hitran/",nxRegistryConfiguration::GLOBAL_INI, true);
	bool						ok;

	ok = config.LocateDirectoryFromKey( "SpectralLineCacheDir", filename, true, true, "Enter location of HITRAN spectral line cache directory");
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"HitranLineStructCache::LoadBaseDirectoryNameFromRegistry, error loading directory name from registry");
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *                  HitranLineStructCache::FindFile                   2019-11-07 */
/** **/
/*---------------------------------------------------------------------------*/

bool HitranPartitionTableCache::FindFile( std::string* filename, bool* exists)
{
	nxString	basedir;
	nxString	name;
	bool		ok;
	int			molnum;
	int			isotopeid;

	molnum    = (int)m_parent->MoleculeNumber();
	isotopeid = (int)m_parent->IsotopeNumber();

	ok = LoadBaseDirectoryNameFromRegistry( &basedir );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"HitranLineStructCache::FindFile, error fetching base directory for the spectral line cache for molecule id (%d)", (int)molnum);
		filename->clear();
		*exists = false;
	}
	else
	{
		basedir.EnsureLastCharIsDirectoryChar();
		if (!nxDirectory::FileExists(basedir) && *exists)
		{
			nxDirectory::CreateADirectory(basedir);
		}
		name.sprintf("%sqcache_mol%03d_iso%03d.bin",(const char*)basedir, (int)molnum, (int)isotopeid);
		name.MakeDirectorySeparatorsOSConsistent();
		*filename = (const char*)name;
		*exists   = nxDirectory::FileExists(name);
	
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *              HitranPartitionTableCache::LoadCache              2019-11-08 */
/** **/
/*---------------------------------------------------------------------------*/

bool HitranPartitionTableCache::LoadCache()
{
	bool			ok;
	bool			exists = false;
	std::string		filename;
	FILE*			f;
	int				numpoints;
	int				versionid;

	ok = FindFile(&filename, &exists);
	ok = ok && exists;
	if (ok)
	{
		f = fopen( filename.c_str(), "rb");
		ok = (f != nullptr);
		if (ok)
		{

			ok = ok && (fread( &versionid,   sizeof(int),    1, f) == 1);			// Read in the file version
			ok = ok && (fread( &m_Tstart,    sizeof(double), 1, f) == 1);			// Read in the start temperature
			ok = ok && (fread( &m_deltaT,    sizeof(double), 1, f) == 1);			// Read in the delta T resolution.
			ok = ok && (fread( &numpoints,   sizeof(int),    1, f) == 1);			// Read in the delta T resolution.

			ok = ok && (versionid == m_versionid);								// Check that the versionid is as expected.
			ok = ok && ( numpoints > 10) && (numpoints < 1000000);
			ok = ok && m_Q.SetSize( numpoints );
			ok = ok && (fread( m_Q.UnsafeArrayBasePtr(), sizeof(double),  numpoints, f) == numpoints);			// Read in the delta T resolution.
			m_Tend = m_Tstart + (numpoints-1)*m_deltaT;
			fclose(f);

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"HitranPartitionTableCache::LoadCache, Error loading partition table cache. You should delete file [%s] so it can be regenerated", (const char*)filename.c_str());
			}
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"HitranPartitionTableCache::LoadCache, Error loading partition table cache. Regenerating partition table using software");
			CreateTable();
		}
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *             HitranPartitionTableCache::CreateCache             2019-11-08 */
/** **/
/*---------------------------------------------------------------------------*/

bool HitranPartitionTableCache::CreateCache()
{
	bool			ok;
	bool			exists = true;
	std::string		filename;
	FILE*			f;
	int				numpoints;
	int				versionid = m_versionid;

	CreateTable();
	ok = FindFile(&filename, &exists);
	if (ok)
	{
		numpoints = (int)m_Q.size();
		f = fopen( filename.c_str(), "wb");
		ok = (f != nullptr);
		if (ok)
		{

			ok = ok && (fwrite( &versionid,   sizeof(int),    1, f) == 1);			// Read in the file version
			ok = ok && (fwrite( &m_Tstart,    sizeof(double), 1, f) == 1);			// Read in the start temperature
			ok = ok && (fwrite( &m_deltaT,    sizeof(double), 1, f) == 1);			// Read in the delta T resolution.
			ok = ok && (fwrite( &numpoints,    sizeof(int),    1, f) == 1);			// Read in the delta T resolution.
			ok = ok && (fwrite( m_Q.UnsafeArrayBasePtr(), sizeof(double),  numpoints, f) == numpoints);			// Read in the delta T resolution.
			fclose(f);
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"HitranPartitionTableCache::WriteCache, Error writing partition table cache to file [%s] ", (const char*)filename.c_str());
	}
	return ok;
}




/*---------------------------------------------------------------------------
 *             HitranPartitionTableCache::CreateTable             2019-11-08 */
/** **/
/*---------------------------------------------------------------------------*/

void HitranPartitionTableCache::CreateTable()
{
	size_t	n;
	double	T;

	m_Tstart = 70.0;
	m_Tend   = 450.0;
	m_deltaT = 0.015625;
	n = (size_t)((m_Tend - m_Tstart)/m_deltaT) + 1;
	m_Q.SetSize( n );
	for (size_t i = 0; i < n; i++)
	{
		T = m_Tstart + i*m_deltaT;
		m_Q.at(i) = m_parent->InternalPartition(T);
	}

}


/*---------------------------------------------------------------------------
 *          HitranPartitionTableCache::InterpolateTable           2019-11-08 */
/** **/
/*---------------------------------------------------------------------------*/
double HitranPartitionTableCache::InterpolateTable(double T)
{
	size_t idx;
	double T0;
	double T1;
	double Q[2];
	double q;

	idx = (size_t)((T - m_Tstart)/m_deltaT);
	T0  = m_Tstart + idx*m_deltaT;
	T1 = T0 + m_deltaT;
	Q[0] = m_Q[idx++];
	Q[1] = m_Q[idx];
	q =  nxLinearInterpolate::FromTwoPoints(T, T0, T1, Q);
	return q;
}


/*---------------------------------------------------------------------------
 *          HitranPartitionTableCache::CheckAndLoadCache          2019-11-08 */
/** **/
/*---------------------------------------------------------------------------*/

bool HitranPartitionTableCache::CheckAndLoadCache()
{
	bool	ok;

	ok = (m_Tstart > 0.0);
	if (!ok)
	{
		ok = LoadCache();
		if (!ok)
		{
			ok = CreateCache();
			ok = ok && LoadCache();
		}
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *          HitranPartitioTableCache::InternalPartition           2019-11-08 */
/** **/
/*---------------------------------------------------------------------------*/

double	HitranPartitionTableCache::InternalPartition( double T )
{
	double value;

	CheckAndLoadCache();
	value = ( ( T < m_Tend) && (T >= m_Tstart)) ? InterpolateTable(T): m_parent->InternalPartition(T);
	return value;
}

