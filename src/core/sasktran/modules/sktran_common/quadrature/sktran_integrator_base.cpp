#include "../sktran_common.h"
#include <cmath>


bool SKTRAN_OpticalPropertiesIntegrator_Base::RoundScatterPosition(const HELIODETIC_POINT & start, const HELIODETIC_UNITVECTOR & look, const double & scatterDistance, const SKTRAN_CoordinateTransform_V2* coords, HELIODETIC_VECTOR & scatterVector, HELIODETIC_POINT & scatterPosition) const
{
	bool ok = true;
	double roundedScatterDistance = 0.0;

	scatterVector = start.Vector() + HELIODETIC_VECTOR(look, scatterDistance);

	ok = ok && coords->HelioVectorToHelioPoint(scatterVector, &scatterPosition);

	if (scatterPosition.Altitude() > coords->TOAAltitude() - 1.0)  // dangerous scatter point (close to the top of atmosphere)
	{
		double toaradius = coords->AltitudeToRadius(coords->TOAAltitude());
		double targetradius = toaradius - 1.0;
		double t1 = -(start.Vector() & look); // distance to tangent point
		double d = t1 * t1 + targetradius * targetradius - start.Radius() * start.Radius(); // discriminant for distance to intersection with target radius

		if (d >= 0.0) // the given line of sight intersects with the target radius (1m below toa)
		{
			if ((scatterVector & look) > 0.0) // look vector is exiting the atmosphere at the dangerous scatter point; we want the larger root
			{
				roundedScatterDistance = t1 + sqrt(d);
			}
			else if ((scatterVector & look) < 0.0) // look vector is entering the atmosphere at the dangerous scatter point; we want the smaller root
			{
				roundedScatterDistance = t1 - sqrt(d);
			}
			else // look vector is tangent at the dangerous scatter point
				 // this shouldn't happen because d >= 0.0 requieres that the altitude be <= target radius,
				 // but scatterPosition.Altitude() > coords->TOAAltitude() - 1.0 requires that the altitude be > target radius
				 // if it does happen, the tangent point should be very near the target altitude, so just use it
			{
				roundedScatterDistance = t1;
			}
		}
		else // the given line of sight does not intersect with the target radius
		{
			// the only way this can happen is if the start point is outside, and the line of sight misses the atmosphere or just brushes it (only interacting with the top 1m)
			// this shouldn't happen unless the initial los skims the top 1m, but in this case the radiance is pretty much 0 anyways
			// or unless a previous scatter point is outside the target radius, but the previous use of this function should avoid this
			ok = false;
		}
		scatterVector = start.Vector() + HELIODETIC_VECTOR(look, roundedScatterDistance);

		ok = ok && coords->HelioVectorToHelioPoint(scatterVector, &scatterPosition);
		ok = ok && scatterPosition.Altitude() <= coords->TOAAltitude();
	}

	return ok;
}

SKTRAN_OpticalPropertiesIntegrator_Base::SKTRAN_OpticalPropertiesIntegrator_Base( )
{
	m_opticalprops=NULL;
}

SKTRAN_OpticalPropertiesIntegrator_Base::~SKTRAN_OpticalPropertiesIntegrator_Base( )
{
	ReleaseResources();
}

void SKTRAN_OpticalPropertiesIntegrator_Base::ReleaseResources( )
{
	if(NULL!=m_opticalprops) m_opticalprops->Release(); 
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::SetOpticalProps		2013-06-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_OpticalPropertiesIntegrator_Base::SetOpticalProps( const SKTRAN_TableOpticalProperties_Base* optprop )
{
	bool ok = true;

	ok = ok && nullptr!=optprop;

	if( ok )
	{
		optprop->AddRef();
		ReleaseResources();
		m_opticalprops = optprop;
	} else	{
		nxLog::Record( NXLOG_WARNING, "Error, optical properties table is NULL in SKTRAN_OpticalPropertiesIntegrator_Straight::SetOpticalProps" );
		ok = false;
	}


	return ok;
}


SKTRAN_SourceTermIntegrator_Base::SKTRAN_SourceTermIntegrator_Base( )
{
	m_opticalprops = nullptr;
	m_fidx = 0;
}

SKTRAN_SourceTermIntegrator_Base::~SKTRAN_SourceTermIntegrator_Base( )
{
	ReleaseResources( );
}

void SKTRAN_SourceTermIntegrator_Base::ReleaseResources( )
{
	if(nullptr!=m_opticalprops) m_opticalprops->Release(); m_opticalprops=nullptr; 
}


bool SKTRAN_SourceTermIntegrator_Base::SetOpticalProps( const SKTRAN_TableOpticalProperties_Base* optprop )
{
	bool ok = true;
	
	ok = ok && optprop!=nullptr;

	if( ok )
	{
		optprop->AddRef();
		ReleaseResources();
		m_opticalprops = optprop;
	} else{
		nxLog::Record( NXLOG_WARNING, "Error, optical properties table is nullptr in SKTRAN_SourceTermIntegrator_Base::SetOpticalProps" );
	}

	return ok;
}
