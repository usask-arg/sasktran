#include "../sasktranv21_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_RayTracing_V21::SKTRAN_SpecsInternal_RayTracing_V21		2010-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsInternal_RayTracing_V21::SKTRAN_SpecsInternal_RayTracing_V21()
{
	m_raytracingshells.reset(new SKTRAN_GridDefRayTracingShells_V21);
	m_raytracingshells->SetStatic();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_RayTracing_V21::~SKTRAN_SpecsInternal_RayTracing_V21		2010-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsInternal_RayTracing_V21::~SKTRAN_SpecsInternal_RayTracing_V21()
{
//	if (m_raytracingshells  != NULL) m_raytracingshells->Release();
//	m_raytracingshells = NULL;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_RayTracing_V21::MaxShellsAlongRay		2010-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_SpecsInternal_RayTracing_V21::MaxShellsAlongRay() const
{
	size_t	  nshells;
	nshells = m_raytracingshells->NumShells();
	return  (nshells*2 + nshells/2 + 2);
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_RayTracing_V21::ConfigureRayTracingShellAlts		2008-1-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsInternal_RayTracing_V21::ConfigureRayTracingShellAlts( const double*  alts, size_t npts )
{
	bool					ok;

	ok = m_raytracingshells->ConfigureHeights( alts, npts );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsInternal_RayTracing_V21::ConfigureRayTracingShellAlts, Error creating internal ray tracing specs. Thats a problem");
	}
	return (ok);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_RayTracing_V21::SKTRAN_SpecsUser_RayTracing_V21		2010-5-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsUser_RayTracing_V21::SKTRAN_SpecsUser_RayTracing_V21()
{
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_RayTracing_V21::ConfigureRayTracingShellAlts		2008-1-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_RayTracing_V21::ConfigureRayTracingShellAlts( const double*  alts, size_t npts )
{
	size_t					i;
	bool					ok = true;
	double					lastalt = -6500000.0;

	m_raytracingshells.resize(npts);
	for (i = 0; i < npts; i++ )
	{
		ok = ok && (lastalt < alts[i]);
		NXASSERT((ok));
		lastalt = alts[i];
		m_raytracingshells[i] = alts[i];
	}
	if (ok)
	{
		ok  = (m_raytracingshells.back() > 990);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsInternal_RayTracing_V21::ConfigureRayTracingShellAlts, The ray tracing altitudes appear to be in kilometers, they should be expressed in meters, Please Check");
		}

	}
	return (ok);
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_RayTracing_V21::CreateInternalSpecs		2010-5-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_RayTracing_V21::CreateInternalSpecs( const SKTRAN_SpecsInternal_RayTracing_V21** userinternalraytracing) const
{
	SKTRAN_SpecsInternal_RayTracing_V21*	internalraytracing;
	bool									ok;

	internalraytracing = new SKTRAN_SpecsInternal_RayTracing_V21;
	ok = (internalraytracing != NULL );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_RayTracing_V21::CreateInternalSpecs, Error allocating memory");
	}
	else
	{
		internalraytracing->AddRef();
		ok = internalraytracing->ConfigureRayTracingShellAlts( &(m_raytracingshells[0]), m_raytracingshells.size() );
	}
	*userinternalraytracing = internalraytracing;
	return ok;
}
