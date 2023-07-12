#include "include/sktran_hr_internals.h"

bool SKTRAN_HR_Perturbation_Absorption::ExtinctionPerturbation( const HELIODETIC_POINT& location,
														         bool& isperturbation,
																 double& value ) const
{
	double alt = location.Altitude();

	if( alt <= m_upperbound && alt >= m_lowerbound )
	{
		isperturbation = true;
		value = GetPertVal();
	}
	else
	{
		isperturbation = false;
		value = 0.0;
	}
	return true;
}


bool SKTRAN_HR_Perturbation_Absorption::Initialize( double upperheight,
													double lowerheight,
													double pertval )
{
	m_upperbound = upperheight;
	m_lowerbound = lowerheight;
	SetPertVal(pertval);
	return true;
}

HELIODETIC_POINT SKTRAN_HR_Perturbation_Absorption::PerturbationLocation( const SKTRAN_CoordinateTransform_V2& coords ) const
{
	return coords.ReferencePoint( (m_upperbound + m_lowerbound)/2);
}

std::unique_ptr< SKTRAN_GeometryObject > SKTRAN_HR_Perturbation_Absorption::BoundingGeometryObject( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, size_t idx ) const
{
	std::unique_ptr<SKTRAN_GeometryObject> ret;

	double groundradius = coords->AltitudeToRadius( 0.0 );
	if( idx == 0 )
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject> ( new SKTRAN_GeometryObject_Sphere( m_lowerbound + groundradius ) );
	}
	else if ( idx == 1 )
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject> ( new SKTRAN_GeometryObject_Sphere( (m_lowerbound + m_upperbound)/2 + groundradius ) );
	}
	else if ( idx == 2 )
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject> ( new SKTRAN_GeometryObject_Sphere( m_upperbound + groundradius ) );
	}
	else
	{
		ret = nullptr;
	}
	return ret;
}



bool SKTRAN_HR_Perturbation_Absorption_Linear::Initialize( double center, double distancetozerobelow, double distancetozeroabove, double pertval )
{
	SetPertVal(pertval);
	m_center = center;
	m_distancetozerobelow = distancetozerobelow;
	m_distancetozeroabove = distancetozeroabove;
	
	return true;
}

bool SKTRAN_HR_Perturbation_Absorption_Linear::ExtinctionPerturbation( const HELIODETIC_POINT& location,
																	   bool& isperturbation,
																	   double& value ) const
{
	bool ok = true;

	double alt = location.Altitude();
	double disttocenter = alt - m_center;

	double deltadistance = disttocenter > 0.0 ? m_distancetozeroabove : m_distancetozerobelow;


	if( abs(disttocenter) < deltadistance )
	{
		isperturbation = true;
		value = GetPertVal() * (1-(abs(disttocenter) / deltadistance));
	}
	else
	{
		isperturbation = false;
		value = 0.0;
	}
	
	return ok;
}

HELIODETIC_POINT SKTRAN_HR_Perturbation_Absorption_Linear::PerturbationLocation( const SKTRAN_CoordinateTransform_V2& coords ) const
{
	return coords.ReferencePoint( (m_center) );
}

std::unique_ptr< SKTRAN_GeometryObject > SKTRAN_HR_Perturbation_Absorption_Linear::BoundingGeometryObject( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, size_t idx ) const
{
	std::unique_ptr<SKTRAN_GeometryObject> ret;

	double groundradius = coords->AltitudeToRadius( 0.0 );
	if( idx == 0 )
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject> ( new SKTRAN_GeometryObject_Sphere( m_center - m_distancetozerobelow + groundradius ) );
	}
	else if ( idx == 1 )
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject> ( new SKTRAN_GeometryObject_Sphere( m_center + groundradius ) );
	}
	else if ( idx == 2)
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject> ( new SKTRAN_GeometryObject_Sphere( m_center + m_distancetozeroabove + groundradius ) );
	}
	else
	{
		ret = nullptr;
	}
	return ret;
}
bool SKTRAN_HR_Perturbation_Absorption_Box::Initialize( double center, double distancetozero, const nxVector& normal1, const nxVector& normal2, double pertval )
{
	SetPertVal(pertval);
	m_centerradius = center;
	m_radiusdistancetozero = distancetozero;

	m_normal1 = normal1;
	m_normal2 = normal2;

	m_anglebetweennormals = abs(nxmath::acosd(m_normal1 & m_normal2));
	
	return true;
}

HELIODETIC_POINT SKTRAN_HR_Perturbation_Absorption_Box::PerturbationLocation( const SKTRAN_CoordinateTransform_V2& coords ) const
{
	return coords.ReferencePoint( m_centerradius );
}

/*
bool SKTRAN_HR_Perturbation_Absorption_Box::ExtinctionPerturbation( const HELIODETIC_POINT& location,
																	   bool& isperturbation,
																	   double& value ) const
{
	bool ok = true;

	double alt = location.Altitude();
	double disttocenter = abs(alt - m_centerradius);

	if( disttocenter < m_radiusdistancetozero )
	{
		const HELIODETIC_UNITVECTOR& heliounit = location.Vector().UnitVector();
		nxVector locunit( heliounit.X(), heliounit.Y(), heliounit.Z() );

		double angleToNormal1 = abs(abs(nxmath::acosd( locunit & m_normal1 ))-90);
		double angleToNormal2 = abs(abs(nxmath::acosd( locunit & m_normal2 )-90));
		double angleBetweenNormals = abs(nxmath::acosd( m_normal1 & m_normal2 ));

		if( angleToNormal1 <= angleBetweenNormals && angleToNormal2 <= angleBetweenNormals )
		{
			isperturbation = true;
			value = GetPertVal() * (1-(disttocenter / m_radiusdistancetozero));
			// linear interpolation in angle from the middle
			value *= (angleBetweenNormals/2 - abs(angleToNormal1-angleBetweenNormals/2))/(angleBetweenNormals/2);
		}
		else
		{
			isperturbation = false;
			value = 0.0;
		}
	}
	else
	{
		isperturbation = false;
		value = 0.0;
	}
	
	return ok;
}
*/

std::unique_ptr< SKTRAN_GeometryObject > SKTRAN_HR_Perturbation_Absorption_Box::BoundingGeometryObject( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, size_t idx ) const
{
	std::unique_ptr<SKTRAN_GeometryObject> ret;

	double groundradius = coords->AltitudeToRadius( 0.0 );
	if( idx == 0 )
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject> ( new SKTRAN_GeometryObject_Sphere( m_centerradius - m_radiusdistancetozero + groundradius ) );
	}
	else if ( idx == 1 )
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject> ( new SKTRAN_GeometryObject_Sphere( m_centerradius + m_radiusdistancetozero + groundradius ) );
	}
	else if ( idx == 2 )
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject> ( new SKTRAN_GeometryObject_Plane( m_normal1 ) );
	}
	else if ( idx == 3 )
	{
		ret = std::unique_ptr<SKTRAN_GeometryObject> ( new SKTRAN_GeometryObject_Plane( m_normal2 ) );
	}
	else
	{
		ret = nullptr;
	}
	return ret;
}