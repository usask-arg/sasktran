#include "include/sktran_hr_internals.h"


bool SKTRAN_HR_RayTracingRegionManager::UpdateUndefinedParametersFromLinesOfSight	( const SKTRAN_LineOfSightArray_V21& linesofsight )
{
	bool ok = true;

	ok = ok && SKTRAN_RayTracingRegionManager::UpdateUndefinedParametersFromLinesOfSight( linesofsight );
	ok = ok && UpdateBoundingReferences( linesofsight );
	ok = ok && UpdateLOSScatteringAngles( linesofsight );

	return ok;
}

bool SKTRAN_HR_RayTracingRegionManager::UpdateBoundingReferences( const SKTRAN_LineOfSightArray_V21& linesofsight )
{
	bool ok = true;

	std::vector<nxVector> points;
	//double dotp;
	double maxdotp = 5;
	nxVector avgstart;
	nxVector avgend;
	size_t numpointsin = 0;
	size_t numpointsout = 0;

	const SKTRAN_LineOfSightEntry_V2* entry;
	
	points.resize( linesofsight.NumRays()*2 );
	for( size_t losidx = 0; losidx < linesofsight.NumRays(); losidx++ )
	{
		ok = ok && linesofsight.GetRay( losidx, &entry );
		// Sometimes looking perfectly at the horizon this can fail due to some rounding errors, but if it does the end points are just set to
		// invalid so ignore them
		GetRayEndpoints( entry->Observer(), entry->Look(), &points[2*losidx], &points[2*losidx+1] );
		if( points[2*losidx].IsValid() && !points[2*losidx].IsZero() )
		{
			avgstart += points[2*losidx];
			numpointsin++;
		}
		if( points[2*losidx+1].IsValid() && !points[2*losidx+1].IsZero() )
		{
			avgend += points[2*losidx+1];
			numpointsout++;
		}
	}

	/*
	for( size_t startidx = 0; startidx < points.size(); startidx++ )
	{
		for( size_t endidx = startidx+1; endidx < points.size(); endidx++ )
		{
			dotp = points[startidx].UnitVector().Dot( points[endidx].UnitVector() );
			if( dotp < maxdotp )
			{
				m_inreferencepoint = points[startidx];
				m_outreferencepoint = points[endidx];
				maxdotp = dotp;
			}
		}
	}
	*/
	m_inreferencepoint = avgstart * (1.0/double(numpointsin));
	m_outreferencepoint = avgend * (1.0/double(numpointsout));

	return ok;
}

bool SKTRAN_HR_RayTracingRegionManager::UpdateLOSScatteringAngles( const SKTRAN_LineOfSightArray_V21& linesofsight )
{
	bool ok = true;
	double ssa;
	m_losminssa = 180.0;
	m_losmaxssa = 0.0;

	const SKTRAN_LineOfSightEntry_V2* entry;

	nxVector sun;
	this->GetSun(&sun);

	for(int losidx = 0; losidx < linesofsight.NumRays(); losidx++)
	{
		ok = ok && linesofsight.GetRay( losidx, &entry );

		ssa = entry->Look().AngleTo(sun);
		if( ssa < m_losminssa )
		{
			m_losminssa = ssa;
		}
		if( ssa > m_losmaxssa )
		{
			m_losmaxssa = ssa;
		}
	}
	return ok;
}


bool SKTRAN_HR_RayTracingRegionManager::GetBoundingReferences( nxVector& in, nxVector& out ) const
{
	bool ok = true;

	in = m_inreferencepoint;
	out = m_outreferencepoint;

	return ok;
}

bool SKTRAN_HR_RayTracingRegionManager::GetBoundingLOSScatteringAngles( double& minssa, double& maxssa ) const
{
	minssa = m_losminssa;
	maxssa = m_losmaxssa;

	return true;
}

bool SKTRAN_HR_RayTracingRegionManager::GetRayEndpointsObserverOutside( const nxVector& observer, const nxVector& look, nxVector* startpt, nxVector* endpt)
{
	double		tpheight;
	nxVector	tplocation;
	nxVector	dummy;
	bool		ok;

	Geoid().FromTangentPointLocation( observer, look );													// Get the tangent point	
	tplocation = Geoid().Location();																	// Get the location
	tpheight = Geoid().Height();																		// Get its altitude
	if (tpheight < GroundAlt())																			// If the observer is looking at the ground
	{																									// Then get the two intercepts, one with upper shell one with ground
		ok =      Geoid().GetShellHeightLocation( MaxAlt(),       observer, look, startpt, &dummy, &tplocation, tpheight );	// start of ray is at upper shell
		ok = Geoid().GetShellHeightLocation( GroundAlt(),    observer, look, endpt,   &dummy, &tplocation, tpheight );
	}																									// and that is that
	else																								// the observer is looking above ground
	{																									// so get the entry and exit points.
		ok =  Geoid().GetShellHeightLocation( MaxAlt(), observer, look, startpt, endpt, &tplocation, tpheight );
	}
	return ok;
}