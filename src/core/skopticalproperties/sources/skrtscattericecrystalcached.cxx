#include <skopticalproperties21.h>



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystalCached::ThreadData::ThreadData		2014-2-22*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_IceCrystalCached::ThreadData::ThreadData()
{
	m_isdirtycached		  = true;
	m_wavenumber          = -99999.0;
}

/*-----------------------------------------------------------------------------
 *			skOpticalProperties_IceCrystalCached::skOpticalProperties_IceCrystalCached		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_IceCrystalCached::skOpticalProperties_IceCrystalCached()
{
	nxString	dirname;
	m_cachebasedir = LoadDirectoryNameFromRegistry();
	m_cachebasedir.EnsureLastCharIsDirectoryChar();
	m_numangles = 541;
}


/*-----------------------------------------------------------------------------
 *				skOpticalProperties_IceCrystalCached::~skOpticalProperties_IceCrystalCached		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_IceCrystalCached::~skOpticalProperties_IceCrystalCached()
{
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystalCached::SetRefractiveIndex		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystalCached::Set_RefractiveIndex( skRTRefractiveIndex* ri )
{
	SetDirtyCached();
	return skOpticalProperties_IceCrystal::Set_RefractiveIndex(ri);
}

/*-----------------------------------------------------------------------------
 *			skOpticalProperties_IceCrystalCached::Set_ParticleDistribution		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystalCached::Set_ParticleDistribution( skRTParticleDist* newdistrib )
{
	SetDirtyCached();
	return skOpticalProperties_IceCrystal::Set_ParticleDistribution( newdistrib );
}

/*-----------------------------------------------------------------------------
 *			skOpticalProperties_IceCrystalCached::Set_ScatterAlgorithm			2009-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystalCached::Set_ScatterAlgorithm( sk_NonsphericalParticle* newscat )
{
	SetDirtyCached();
	return skOpticalProperties_IceCrystal::Set_ScatterAlgorithm( newscat );
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystalCached::SetNumScatteringAngles		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystalCached::SetNumScatteringAngles( size_t numangles )
{
	m_numangles = numangles;

	for (auto& entry: m_data )
	{
		entry.second.m_phasematrix.resize(numangles);
	}
	SetDirtyCached();
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystalCached::SetDirtyCached		2014-2-22*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_IceCrystalCached::SetDirtyCached()
{
	for (auto& entry : m_data)
	{
		entry.second.m_isdirtycached = true;
	}
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystalCached::CheckDirtyAndUpdate		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystalCached::CheckDirtyAndUpdate(double wavenumber, ThreadData* data) 
{
	bool	ok;
	ok = (!data->m_isdirtycached && (wavenumber == data->m_wavenumber));
	if (!ok)
	{
		data->m_wavenumber = wavenumber;
		ok = UpdateTables( data );
		data->m_isdirtycached = !ok;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *				skOpticalProperties_IceCrystalCached::UpdateTables		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystalCached::UpdateTables( ThreadData* data) 
{
	nxString	filename;
	bool		ok;
	static std::mutex		cachemutex;

	std::lock_guard<std::mutex>		lock( cachemutex);				// Lock the cache mutex so no other thread is reading/or writing at this time

	filename = FullCacheName( data);				
	ok       = ReadCacheFile( filename,  data);			
	if (!ok)
	{
		ok =       CreateTables  ( data);
		ok = ok && WriteCacheFile( filename, data );
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *				skOpticalProperties_IceCrystalCached::UpdateTables		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystalCached::CreateTables(ThreadData* data) 
{
	size_t	idx;
	double	dmu;
	double	mu;
	bool	ok = true;
	size_t	numangles;
	

	NXTRACE_ONCEONLY(firsttime,("**** CHECK *****, skOpticalProperties_IceCrystalCached::CreateTables, Do the scattering angles defined in skOpticalProperties_IceCrystalCached, get through to class sk_NonsphericalParticle\n"));
	NXTRACE_ONCEONLY(firsttimeb,("***** CHECK **** skOpticalProperties_IceCrystalCached::CreateTables, It is implicitly assumed that the tables are linearly spaced. This is not the default case in sk_TmatrixRandomWrapper\n"));
	numangles  = data->NumScatteringAngles();
	dmu        = 2.0/(numangles-1);
	ok = skOpticalProperties_IceCrystal::CalculateCrossSections( data->m_wavenumber, &data->m_absxs, &data->m_extxs, &data->m_scattxs);
	if (ok)
	{
		for (idx = 0; idx < numangles; idx++ )
		{
			mu = -1.0 + (idx*dmu);
			ok = ok && skOpticalProperties_IceCrystal::CalculatePhaseMatrix( data->m_wavenumber, mu, &data->m_phasematrix.at(idx) );
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystalCached::FullCacheName		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skOpticalProperties_IceCrystalCached::FullCacheName( ThreadData* data ) 
{
	nxString	dirname;
	nxString	basename;
	nxString	fullfilename;
	size_t		wavelenint;
	double		wavelennm;

	wavelennm = 1.0E7/data->m_wavenumber;
	wavelenint = (size_t)(wavelennm*10);

	basename.sprintf("icecrystal/%s/%s%s_%0Iu.dat", (const char*)(skOpticalProperties_IceCrystal::Get_RefractiveIndex()->ChemicalName()),
												    (const char*)(skOpticalProperties_IceCrystal::Get_Distribution()->CachingDescriptor()), 
												    (const char*) Get_DescriptiveParticleString(),
												    (size_t)wavelenint );
	
	fullfilename = m_cachebasedir + basename;

	return fullfilename;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystalCached::LoadParametersFromRegistrySetting		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skOpticalProperties_IceCrystalCached::LoadDirectoryNameFromRegistry( )
{
	nxRegistryConfiguration		config( "USask-ARG", "skOpticalProperties/NonSpherical_IceCrystals/",nxRegistryConfiguration::GLOBAL_INI, true);
	nxString					filename;
	bool						ok;

	ok = config.LocateDirectoryFromKey( "Aerosol_Cache_Directory", &filename, true, true, "Enter location for caching ice crystal parameters");
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_IceCrystalCached::LoadDirectoryNameFromRegistry, error loading aersol cache directory name from registry, using TEMP");
		filename = getenv("TEMP");
	}
	return filename;
}

/*-----------------------------------------------------------------------------
 *				skOpticalProperties_IceCrystalCached::ReadCacheFile		2009-6-4*/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystalCached::ReadCacheFile( const char* filename,  ThreadData* data )
{
	bool	ok;
	nxFile	f;
	nxDWORD	numscat = 0;
	double	wavenumber;
	double	kabs;
	double	kext;
	double	kscat;

	f.Open( filename, "rb");
	ok = f.IsOpen();
	if (ok)
	{
		ok =       (f.Read( &wavenumber,          sizeof(wavenumber),      1 ) == 1);
		ok =       (f.Read( &kabs,				  sizeof(kabs),              1 ) == 1);
		ok =       (f.Read( &kext,				  sizeof(kext),              1 ) == 1);
		ok =       (f.Read( &kscat,			      sizeof(kscat),             1 ) == 1);
		ok = ok && (f.Read( &numscat,             sizeof(numscat),           1 ) == 1);
		data->m_absxs   = kabs;
		data->m_extxs   = kext;
		data->m_scattxs = kscat;
		ok = ok && data->SetNumScatteringAngles( numscat );
		ok = ok && (f.Read( &data->m_phasematrix.front(),   sizeof(data->m_phasematrix.front() ),    numscat) == numscat);
		f.Close();
		data->m_isdirtycached = !ok;
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_IceCrystalCached::ReadCacheFile, Error reading ice crystal file cache <%s>. Thats a problem!",(const char*)filename);
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *				skOpticalProperties_IceCrystalCached::ReadCacheFile		2009-6-4*/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystalCached::WriteCacheFile( const char* filename, ThreadData* data ) const
{
	bool	ok;
	nxFile	f;
	nxDWORD	numscat;
	nxFileSpec	spec(filename);
	double	wavenumber;
	double	kabs;
	double	kext;
	double	kscat;

	nxDirectory::CreateADirectory( spec.FullDirSpec() );
	f.Open( filename, "wb");
	ok = f.IsOpen();
	if (ok)
	{
		numscat = (nxDWORD)data->NumScatteringAngles();
		wavenumber = data->m_wavenumber;
		kabs       = data->m_absxs;
		kext       = data->m_extxs;
		kscat      = data->m_scattxs;

		ok =       (f.Write( &wavenumber,          sizeof(wavenumber),        1 ) == 1);
		ok =       (f.Write( &kabs,				   sizeof(kabs),              1 ) == 1);
		ok =       (f.Write( &kext,				   sizeof(kext),              1 ) == 1);
		ok =       (f.Write( &kscat,			   sizeof(kscat),             1 ) == 1);
		ok = ok && (f.Write( &numscat,             sizeof(numscat),           1 ) == 1);
		ok = ok && (f.Write( &data->m_phasematrix.front(),        sizeof(data->m_phasematrix.front()),    numscat) == numscat);
		f.Close();
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_IceCrystalCached::WriteCacheFile, Error writing ice crystal cache to file <%s>", (const char*)filename );
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::LookupUpThreadData		 2014- 10- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystalCached::LookupUpThreadData(ThreadData** data )
{
	size_t					threadnum;
	iterator				iter;
	bool					ok;
	static std::mutex		lock;							// Add a lock to make sure map.insert is thread safe.

	threadnum = nxWorkerThreadManager::GetCurrentThreadIdCode();
	lock.lock();
	iter      = m_data.find(threadnum);
	ok = (iter != m_data.end() );
	if (!ok)
	{
		ThreadData									blank;
		std::pair<iterator, bool>					result;
		std::pair<size_t, ThreadData>				newentry(threadnum, blank);

		result = m_data.insert(  newentry);
		ok = result.second;
		iter = result.first;
		(*iter).second.m_isdirtycached = true;
		(*iter).second.m_phasematrix.resize(m_numangles);
	}
	lock.unlock();
	if (!ok)
	{
		*data = nullptr;
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_Convolved::LookupUpThreadData, error looking/creating thread data for thread %d", (int)threadnum);
	}
	else
	{
		*data = &(*iter).second;
	}
	return ok;

}
/*-----------------------------------------------------------------------------
 *				skOpticalProperties_IceCrystalCached::CalculateCrossSections		2009-6-4*/
/*---------------------------------------------------------------------------*/

bool  skOpticalProperties_IceCrystalCached::CalculateCrossSections( double wavenumber, 	double*	absxs, double* extxs,  double* scattxs )
{
	bool ok;
	ThreadData* data;

	ok       =       LookupUpThreadData(&data);
	ok       = ok && CheckDirtyAndUpdate( wavenumber, data);
	if (ok)
	{
		*absxs   = data->m_absxs;
		*extxs   = data->m_extxs;
		*scattxs = data->m_scattxs;
	}
	else
	{
		*absxs   = std::numeric_limits<double>::quiet_NaN();
		*extxs   = std::numeric_limits<double>::quiet_NaN();
		*scattxs = std::numeric_limits<double>::quiet_NaN();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *				skOpticalProperties_IceCrystalCached::CalculatePhaseMatrix		2009-6-4*/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystalCached::CalculatePhaseMatrix( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix )
{
	bool		ok;
	ThreadData* data;

	ok  =       LookupUpThreadData(&data);
	ok  = ok && CheckDirtyAndUpdate( wavenumber, data);
	ok  = ok && GetPhaseMatrix( data, cosscatterangle, phasematrix);
	return ok;
}

/*-----------------------------------------------------------------------------
 *				skOpticalProperties_IceCrystalCached::GetPhaseMatrix		2009-6-4*/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystalCached::GetPhaseMatrix( ThreadData* data, double cosscatterangle, skRTPhaseMatrix* phasematrix )
{
	bool	ok;
	size_t	idx;
	size_t  numscat = data->NumScatteringAngles();
	double	dmu = 2.0/(numscat-1);

	NXASSERT(( (cosscatterangle >= -1.0) && (cosscatterangle <= 1.0) ));
	ok = (cosscatterangle >= -1.0) && (cosscatterangle <= 1.0);
	if (!ok)
	{
		NXTRACE_ONCEONLY(firsttime,("**** 2009-5-4 ***** skOpticalProperties_IceCrystalCached::GetPhaseMatrix, cosscatterangle is outside range -1 to + 1\n"));
		if (cosscatterangle < -1.0) cosscatterangle = -1.0;
		if (cosscatterangle >  1.0) cosscatterangle =  1.0;
	}
	idx = (size_t)((cosscatterangle + 1.0)/dmu + 0.5);
	if (idx >= numscat) idx = numscat -1;
	*phasematrix = data->m_phasematrix.at(idx);
	return ok;
}

