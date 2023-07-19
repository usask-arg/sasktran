/**
 * SASKTRAN TIR Ray Tracing Region Manager
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_RayTracingRegionManager::UpdateUndefinedParametersFromLinesOfSight
 * 2018-09-13
 */
bool SKTRAN_TIR_RayTracingRegionManager::UpdateUndefinedParametersFromLinesOfSight(
	const SKTRAN_LineOfSightArray_V21& linesofsight)
{
	bool ok = true;

	ok = ok && SKTRAN_RayTracingRegionManager::UpdateUndefinedParametersFromLinesOfSight(linesofsight);
	ok = ok && UpdateBoundingReferences(linesofsight);

	return ok;
}

/**
 * SKTRAN_TIR_RayTracingRegionManager::UpdateBoundingReferences
 * 2018-09-13
 */
bool SKTRAN_TIR_RayTracingRegionManager::UpdateBoundingReferences(
	const SKTRAN_LineOfSightArray_V21& linesofsight)
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

	points.resize(linesofsight.NumRays() * 2);
	for (size_t losidx = 0; losidx < linesofsight.NumRays(); losidx++)
	{
		ok = ok && linesofsight.GetRay(losidx, &entry);
		// Sometimes looking perfectly at the horizon this can fail due to some rounding errors, but if it does the end points are just set to
		// invalid so ignore them
		GetRayEndpoints(entry->Observer(), entry->Look(), &points[2 * losidx], &points[2 * losidx + 1]);
		if (points[2 * losidx].IsValid() && !points[2 * losidx].IsZero())
		{
			avgstart += points[2 * losidx];
			numpointsin++;
		}
		if (points[2 * losidx + 1].IsValid() && !points[2 * losidx + 1].IsZero())
		{
			avgend += points[2 * losidx + 1];
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
	m_inreferencepoint = avgstart * (1.0 / double(numpointsin));
	m_outreferencepoint = avgend * (1.0 / double(numpointsout));

	return ok;
}

/**
 * SKTRAN_TIR_RayTracingRegionManager::GetBoundingReferences
 * 2018-09-13
 */
bool SKTRAN_TIR_RayTracingRegionManager::GetBoundingReferences(
	nxVector& in,
	nxVector& out) const
{
	bool ok = true;

	in = m_inreferencepoint;
	out = m_outreferencepoint;

	return ok;
}

/**
 * SKTRAN_TIR_RayTracingRegionManager::GetRayEndpointsObserverOutside
 * 2018-09-13
 */
bool SKTRAN_TIR_RayTracingRegionManager::GetRayEndpointsObserverOutside(
	const nxVector& observer,
	const nxVector& look,
	nxVector* startpt,
	nxVector* endpt)
{
	double		tpheight;
	nxVector	tplocation;
	nxVector	dummy;
	bool		ok;

	Geoid().FromTangentPointLocation(observer, look);													// Get the tangent point	
	tplocation = Geoid().Location();																	// Get the location
	tpheight = Geoid().Height();																		// Get its altitude
	if (tpheight < GroundAlt())																			// If the observer is looking at the ground
	{																									// Then get the two intercepts, one with upper shell one with ground
		ok = Geoid().GetShellHeightLocation(MaxAlt(), observer, look, startpt,
											&dummy, &tplocation, tpheight);	// start of ray is at upper shell
		ok = Geoid().GetShellHeightLocation(GroundAlt(), observer, look, endpt,
											&dummy, &tplocation, tpheight);
	}																									// and that is that
	else																								// the observer is looking above ground
	{																									// so get the entry and exit points.
		ok = Geoid().GetShellHeightLocation(MaxAlt(), observer, look, startpt,
											endpt, &tplocation, tpheight);
	}
	return ok;
}

/**
 * SKTRAN_TIR_RayTracingRegionManager::SetEarthRadius
 * 2019-06-05
 */
bool SKTRAN_TIR_RayTracingRegionManager::SetEarthRadius(
	double earthradius_meters)
{
	bool ok;

	SKTRAN_RayTracingRegionManager::SetEarthRadius(earthradius_meters);

	ok = (earthradius_meters > 100000.0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_RayTracingRegionManager::SetEarthRadius, it looks like your Earth radius (%f) is not in meters, you should enter a value similar to 6371000.000 for Earth. Nothing has been changed.", (double)earthradius_meters);
	}
	else
	{
		m_fixedearthradius = earthradius_meters;
	}

	return ok;
}

/**
 * SKTRAN_TIR_RayTracingRegionManager::Clear
 * 2019-06-05
 */
void SKTRAN_TIR_RayTracingRegionManager::Clear()
{
	SKTRAN_RayTracingRegionManager::Clear();

	m_minalt = 0.0;
}

/**
 * SKTRAN_TIR_RayTracingRegionManager::SetLowerBoundAltitude
 * 2019-06-05
 */
bool SKTRAN_TIR_RayTracingRegionManager::SetLowerBoundAltitude(
	double lowerheight_meters)
{
	bool ok = SKTRAN_RayTracingRegionManager::SetLowerBoundAltitude(lowerheight_meters);
	m_minalt = lowerheight_meters;
	return ok;
}

/**
 * SKTRAN_TIR_RayTracingRegionManager::CheckParameters
 * 2019-06-05
 */
bool SKTRAN_TIR_RayTracingRegionManager::CheckParameters() const
{
	bool	ok = true;
	bool	ok1;
	bool	ok2;

	double mjd;
	nxVector sun;

	ok = ok && GetMJD(&mjd);
	ok = ok && GetSun(&sun);

	ok = (MaxAlt() > 1000) && (MaxAlt() > m_minalt) && (MaxAlt() > GroundAlt());
	if (!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_RayTracingRegionManager::CheckParameters, The user specified altitude ranges are not acceptable");
	ok2 = (mjd >= 1000); if (!ok2) nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_RayTracingRegionManager::CheckParameters, The user has not defined the mjd");
	ok1 = (sun.IsValid()); if (!ok1) nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_RayTracingRegionManager::CheckParameters, The sun has not been defined");
	ok = ok && ok1 && ok2;
	return ok;
}

/**
 * SKTRAN_TIR_RayTracingRegionManager::MakeCoordinateSystem
 * 2019-06-05
 */
bool SKTRAN_TIR_RayTracingRegionManager::MakeCoordinateSystem(
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>* usercoordinates,
	double groundaltitude,
	double toa_altitude,
	nxGeodetic::GEOID_MODEL geoidmodel,
	bool userdefinedgeoidmodel) const
{
	bool ok;
	double latitude;
	double longitude;
	double mjd;
	nxVector sun;
	std::unique_ptr<SKTRAN_CoordinateTransform_V2> coordinates(new SKTRAN_CoordinateTransform_V2);

	if (userdefinedgeoidmodel)
	{
		const double EQUATORIAL_RADIUS = 6378100.0;

		const_cast<nxGeodetic&>(coordinates->TrueGeoid()).SelectGeoid(geoidmodel);

		if (geoidmodel == nxGeodetic::GEOID_SPHERE)
		{
			const_cast<nxGeodetic&>(coordinates->TrueGeoid()).SetTrueSphere(EQUATORIAL_RADIUS);
		}
	}

	ok = IsProperlyDefined();
	ok = ok && CheckParameters();

	if (!ok)
	{
		coordinates.release();
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_RayTracingRegionManager::MakeCoordinateSystem, Cannot make coordinates as the ray tracing region is nor properly defined. Thats a problem");
	}
	else
	{
		ok = GetReferencePoint(&latitude, &longitude);
		ok = ok && GetSun(&sun);
		ok = ok && GetMJD(&mjd);
		ok = ok && coordinates->ConfigureCoordinates(latitude, longitude, mjd, sun);
		ok = ok && coordinates->SetAtmosphereAltitudeBounds(groundaltitude, toa_altitude);
		if (NXFINITE(m_fixedearthradius)) ok = ok && coordinates->ManuallySetOsculatingSphereRadius(m_fixedearthradius);
		coordinates->SetStatic();
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_LineOfSightArray_V21::MakeCoordinateSystem, Error making coordinate system");
		}
	}
	*usercoordinates = std::move(coordinates);
	return ok;
}
