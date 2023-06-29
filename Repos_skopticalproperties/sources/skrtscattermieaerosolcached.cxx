#include <skopticalproperties21.h>
#include <nxbase_threads.h>
#include <omp.h>


static nxString g_aerosolcachedirectory;
static std::mutex g_mutex;


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::skOpticalProperties_MieAerosolCached		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_MieAerosolCached::skOpticalProperties_MieAerosolCached()
{
	SetNumScatteringAngles(5400);
	m_nummoments = 64;
	if ( g_aerosolcachedirectory.IsEmpty() )								// Load the registry directory just once
	{																		// to avoid excessive I/O if we create lots of objects.
		std::lock_guard<std::mutex>	lock(g_mutex);
		g_aerosolcachedirectory = LoadDirectoryNameFromRegistry();
	}
	m_cachebasedir = g_aerosolcachedirectory;
	m_cachebasedir.EnsureLastCharIsDirectoryChar();

}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystalCached::ThreadData::ThreadData		2014-2-22*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_MieAerosolCached::ThreadData::ThreadData()
{
	m_isdirtycached		  = true;
	m_wavenumber          = -99999.0;
//	m_phasematrix.resize( 1080);
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::~skOpticalProperties_MieAerosolCached		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_MieAerosolCached::~skOpticalProperties_MieAerosolCached()
{
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::SetRefractiveIndex		2008-4-15*/
/** Sets the refractive index of the mie aerosol. It is assumed this is
 *	only set once during initialization. Frequently changing the refractive
 *	index object may incur significant overhead as it always flushes the internal
 *	caches even if the refractive obejct is the same
 *
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosolCached::Set_RefractiveIndex( skRTRefractiveIndex* ri )
{
	bool	ok;

	ok = ( Get_RefractiveIndex() == ri );
	if (!ok)
	{
		SetDirty();
		ok = skOpticalProperties_MieAerosol::Set_RefractiveIndex(ri);
	}
	return ok; 
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::Set_ParticleDistribution		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosolCached::Set_ParticleDistribution( skRTParticleDist* distribution )
{
	bool				ok;
	skRTParticleDist*	currentdist;

	currentdist = Get_Distribution();													// Get the current particle distribution
	ok =       (currentdist != NULL) && (distribution != NULL);							// If neither distribution is NULL
	ok = ok && currentdist->IsSameDistributionAs( distribution );						// Then see if they are identical
	if (!ok)																			// if they are then we dont need to do anything
	{																					// otherwise
		SetDirtyCached();																		// Then set the optical proeprties cache as dirty
		ok = skOpticalProperties_MieAerosol::Set_ParticleDistribution( distribution );	// And set the new particle distribution
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::SetNumScatteringAngles		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosolCached::SetNumScatteringAngles( size_t numangles )
{
	m_numangles = numangles;
	SetDirtyCached();
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::SetDirtyCached		2014-2-22*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_MieAerosolCached::SetDirtyCached()
{
	for (auto& entry : m_data)
	{
		entry.second.m_isdirtycached = true;
	}
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::CheckDirtyAndUpdate		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosolCached::CheckDirtyAndUpdate(double wavenumber, ThreadData* data)
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
 *					skOpticalProperties_MieAerosolCached::UpdateTables		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosolCached::UpdateTables(ThreadData* data)
{
	nxString	filename;
	bool		ok;
	static std::mutex		cachemutex;

	std::lock_guard<std::mutex>		lock( cachemutex);				// Lock the cache mutex so no other thread is reading/or writing at this time

	filename = FullCacheName( data );
	ok       = ReadCacheFile( filename, data );
	if (!ok)
	{
		#if defined (NXDEBUG)
			nxLog::Record(NXLOG_INFO, "skOpticalProperties_MieAerosolCached::UpdateTables, creating table for <%s>", (const char*) filename );
		#endif
			ok = CreateTables( data);
			ok = ok && WriteCacheFile( filename, data );
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::UpdateTables		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosolCached::CreateTables(ThreadData* data)
{
	size_t	idx, idy;
	double	dmu;
	double	mu;
	bool	ok = true;
	bool	ok1;
	bool	iszero;
	size_t	numangle = data->NumScatteringAngles();
	size_t  nummoment = data->NumLegendreMoments();
	int     opticalmaxmoment;

	dmu        = 2.0/(numangle-1);
	skOpticalProperties_MieAerosol::CalculateCrossSections(data->m_wavenumber, &data->m_absxs, &data->m_extxs, &data->m_scattxs);
	iszero = (data->m_scattxs == 0.0);

	for (idx = 0; idx < numangle; idx++ )
	{
		if (iszero)
		{
			data->m_phasematrix[idx].SetTo(0.0);
		}
		else
		{
			mu = -1.0 + (idx*dmu);
			ok1 = skOpticalProperties_MieAerosol::CalculatePhaseMatrix( data->m_wavenumber, mu, &data->m_phasematrix.at(idx) );
			ok = ok && ok1;
		}
	}
	for (idx = 0; idx < nummoment; idx++)
	{
		if (iszero)
		{
			data->m_legendrep11[idx] = 0.0;
		}
		else
		{
			ok = ok && skOpticalProperties_MieAerosol::LegendreCoefficientsP11(data->m_wavenumber, data->m_legendrep11.data(), nummoment, opticalmaxmoment);
			// Just ensure that values that aren't calculated are set to 0, even though this probably should never happen
			for (idy = opticalmaxmoment; idy < nummoment; idy++)
			{
				data->m_legendrep11[idy] = 0.0;
			}
		}
	}

	return ok;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::FullCacheName		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skOpticalProperties_MieAerosolCached::FullCacheName( ThreadData* data )
{
	nxString	basename;
	nxString	fullfilename;
	size_t		wavelenint;
	double		wavelennm;
	std::complex<double> r;


	wavelennm = 1.0E7/data->m_wavenumber;
	wavelenint = (size_t)(wavelennm*10);

	r = Get_RefractiveIndex()->RefractiveIndex(data->m_wavenumber);

	basename.sprintf("mie_aerosol/%s/%s_mie_%05u_%08u_%08u_v12.dat", (const char*)(Get_RefractiveIndex()->ChemicalName()),
													  (const char*)(skOpticalProperties_MieAerosol::Get_Distribution()->CachingDescriptor()),
													  (unsigned int)wavelenint,
													  (unsigned int)(r.real() * 1e7),
		                                              (unsigned int)(r.imag() * 1e7));
	fullfilename = m_cachebasedir + basename;
	fullfilename.MakeDirectorySeparatorsOSConsistent();
	return fullfilename;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::LoadParametersFromRegistrySetting		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skOpticalProperties_MieAerosolCached::LoadDirectoryNameFromRegistry( )
{
	nxRegistryConfiguration		config( "USask-ARG", "SkOpticalProperties/MieAerosol/",nxRegistryConfiguration::GLOBAL_INI, true);
	nxString					filename;
	bool						ok;

	ok = config.LocateDirectoryFromKey( "Aerosol_Cache_Directory", &filename, true, true, "Enter location for caching aerosol parameters");
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MieAerosolCached::LoadDirectoryNameFromRegistry, error loading aersol cache directory name from registry, using TEMP");
		filename = getenv("TEMP");
		throw("Corrupted registry");
	}
	return filename;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::ReadCacheFile		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosolCached::ReadCacheFile( const char* filename, ThreadData* data )
{
	bool	ok;
	nxFile	f;
	nxDWORD	numscat = 0;
	nxDWORD nummoment = 0;
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
		ok = ok && (f.Read(&nummoment, sizeof(nummoment), 1) == 1);
		ok = ok && data->SetNumLegendreMoments(nummoment);
		ok = ok && (f.Read(&data->m_legendrep11.front(), sizeof(data->m_legendrep11.front()), nummoment) == nummoment);
		ok = ok && (f.Read( &numscat,             sizeof(numscat),           1 ) == 1);
		data->m_absxs   = kabs;
		data->m_extxs   = kext;
		data->m_scattxs = kscat;
		ok = ok && data->SetNumScatteringAngles( numscat );
		ok = ok && (f.Read( &data->m_phasematrix.front(),  sizeof(data->m_phasematrix.front()),    numscat ) == numscat );
		f.Close();
		data->m_isdirtycached = !ok;
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MieAerosolCached::ReadCacheFile, Error reading aerosol file cache <%s>. Thats a problem!",(const char*)filename);
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::ReadCacheFile		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosolCached::WriteCacheFile( const char* filename, ThreadData* data )
{
	bool	ok;
	nxFile	f;
	nxDWORD	numscat, nummoment;
	double	kscat;
	double	kext;
	double	kabs;
	double	wave;
	nxFileSpec	spec(filename);

	nxDirectory::CreateADirectory( spec.FullDirSpec() );
	f.Open( filename, "wb");
	ok = f.IsOpen();
	if (ok)
	{
		numscat = (nxDWORD)data->NumScatteringAngles();
		nummoment = (nxDWORD)data->NumLegendreMoments();
		kabs    = data->m_absxs;
		kext    = data->m_extxs;
		kscat   = data->m_scattxs;
		wave    = data->m_wavenumber;
		ok =       (f.Write( &wave,				   sizeof(wave),              1 ) == 1);
		ok =       (f.Write( &kabs,				   sizeof(kabs),              1 ) == 1);
		ok =       (f.Write( &kext,				   sizeof(kext),              1 ) == 1);
		ok =       (f.Write( &kscat,			   sizeof(kscat),             1 ) == 1);
		ok = ok && (f.Write(&nummoment, sizeof(nummoment), 1) == 1);
		ok = ok && (f.Write(&data->m_legendrep11.front(), sizeof(data->m_legendrep11.front()), nummoment) == nummoment);
		ok = ok && (f.Write( &numscat,             sizeof(numscat),           1 ) == 1);
		ok = ok && (f.Write( &data->m_phasematrix.front(),        sizeof(data->m_phasematrix.front()),    numscat ) == numscat);
		f.Close();
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MieAerosolCached::WriteCacheFile, Error writing aerosol cache to file <%s>", (const char*)filename );
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::LookupUpThreadData		 2014- 10- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosolCached::LookupUpThreadData(ThreadData** data )
{
	size_t					threadnum;
	iterator				iter;
	bool					ok;
	static std::mutex		lock;							// Add a lock to make sure map.insert is thread safe.

	if (omp_in_parallel())
	{
		lock.lock();
	}

	threadnum = nxWorkerThreadManager::GetCurrentThreadIdCode();
	iter      = m_data.find(threadnum);
	ok = (iter != m_data.end() );
	if (!ok)
	{
		ThreadData						blank;
		std::pair<iterator, bool>		result;
		std::pair<size_t, ThreadData>	newentry(threadnum, blank);

		result = m_data.insert(  newentry);
		ok = result.second;
		iter = result.first;
		if (ok)
		{
			(*iter).second.m_phasematrix.resize(m_numangles);
			(*iter).second.m_legendrep11.resize(m_nummoments);
			(*iter).second.m_isdirtycached = true;

		}
	}
	if (omp_in_parallel())
	{
		lock.unlock();
	}
	if (!ok)
	{
		*data = nullptr;
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MieAerosolCached::LookupUpThreadData, error looking/creating thread data for OMP thread %d", (int)threadnum);
	}
	else
	{
		*data = &(*iter).second;
	}
	return ok;

}
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::CalculateCrossSections		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool  skOpticalProperties_MieAerosolCached::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs)
{
	bool				ok;
	ThreadData*			threaddata;

	ok =       LookupUpThreadData( &threaddata );
	ok = ok && CheckDirtyAndUpdate(wavenumber, threaddata );
	if (ok)
	{
		*absxs   = threaddata->m_absxs;
		*extxs   = threaddata->m_extxs;
		*scattxs = threaddata->m_scattxs;
	}
	else
	{
		*absxs   = std::numeric_limits<double>::quiet_NaN();
		*extxs   = std::numeric_limits<double>::quiet_NaN();
		*scattxs = std::numeric_limits<double>::quiet_NaN();
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MieAerosolCached::CalculateCrossSections, Error calculating cross-sections thats not good, setting crosssections to NaN");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::CalculatePhaseMatrix		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosolCached::CalculatePhaseMatrix( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix)
{
	bool		ok;
	ThreadData* data;

	ok =       LookupUpThreadData( &data );
	ok = ok && CheckDirtyAndUpdate(wavenumber, data );
	ok = ok && GetPhaseMatrix( data, cosscatterangle, phasematrix);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosolCached::GetPhaseMatrix		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosolCached::GetPhaseMatrix( ThreadData* data, double cosscatterangle, skRTPhaseMatrix* phasematrix )
{
	bool	ok;
	size_t	idx;
	size_t	numscat = data->NumScatteringAngles();
	double	dmu = 2.0/(numscat-1);

	NXASSERT(( (cosscatterangle >= -1.0) && (cosscatterangle <= 1.0) ));
	ok = (cosscatterangle >= -1.0) && (cosscatterangle <= 1.0);
	if (!ok)
	{
		NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** skOpticalProperties_MieAerosolCached::GetPhaseMatrix, cosscatterangle is outside range -1 to + 1\n"));
		if (cosscatterangle < -1.0) cosscatterangle = -1.0;
		if (cosscatterangle >  1.0) cosscatterangle =  1.0;
	}
	idx = (size_t)((cosscatterangle + 1.0)/dmu + 0.5);
	if (idx >= numscat) idx = numscat-1;
	*phasematrix = data->m_phasematrix.at(idx);
	return ok;
}

bool skOpticalProperties_MieAerosolCached::LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff) 
{
	bool		ok;
	ThreadData* data;

	ok = LookupUpThreadData(&data);
	ok = ok && CheckDirtyAndUpdate(wavenumber, data);

	opticalmaxcoeff = std::min((size_t)usermaxcoeff, data->NumLegendreMoments());

	for (int i = 0; i < opticalmaxcoeff; ++i)
	{
		coeff[i] = data->m_legendrep11[i];
	}

	return true;
}