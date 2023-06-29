#include "../sasktranv21_internals.h"
/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileGeometry_V21::SKTRAN_DiffuseProfileGeometry_V21		2007-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_DiffuseProfileGeometry_V21::SKTRAN_DiffuseProfileGeometry_V21()
{
	m_parent             = NULL;
//	m_diffusepoints      = NULL;
//	m_numpoints          = 0;
	m_profileidx         = (SKTRAN_GridIndex)(-1);

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileGeometry_V21::~SKTRAN_DiffuseProfileGeometry_V21		2007-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_DiffuseProfileGeometry_V21::~SKTRAN_DiffuseProfileGeometry_V21()
{
	ReleaseResources();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileGeometry_V21::ReleaseResources		2007-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_DiffuseProfileGeometry_V21::ReleaseResources()
{
	m_diffusepoints.clear();
	//if (m_diffusepoints    != NULL) delete [] m_diffusepoints;
	//	m_diffusepoints    = NULL;
	//m_numpoints        = 0;
	m_profileidx         = (SKTRAN_GridIndex)(-1);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileGeometry_V21::CreateOpticalProfile		2010-2-23*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffuseProfileGeometry_V21::CreateOpticalProfile( SKTRAN_DiffuseProfileOptical_V21** optprofile ) 
{
	SKTRAN_DiffuseProfileOptical_V21*	op;
	bool								ok;

	op = new SKTRAN_DiffuseProfileOptical_V21(this);
	ok = (op != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_DiffuseProfileGeometry_V21::CreateOpticalProfile, error allocating memory for profile. Thats a problem");
	}
	*optprofile = op;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileGeometry_V21::AllocatePoints		2007-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffuseProfileGeometry_V21::AllocatePoints(size_t numradii)
{
	bool	ok;
	bool   ok1;

	m_diffusepoints.clear();
	ok   = (numradii == 0);
	if (!ok)
	{
		m_diffusepoints.resize( numradii);
		ok = (m_diffusepoints.size() == numradii);
		if (ok)
		{
			for (size_t i = 0; i < m_diffusepoints.size(); i++)
			{
				ok1 = m_diffusepoints.at(i).CreateOpticalPoint();
				ok = ok && ok1;
			}
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_DiffuseProfileGeometry_V21::AllocatePoints, Error allocating memeory for %Iu diffuse points", (size_t)numradii);
			m_diffusepoints.clear();
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_DiffuseProfileGeometry_V21::ConfigureGeometry		2007-12-12*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_DiffuseProfileGeometry_V21::ConfigureGeometry_Stage1 (	const  HELIODETIC_POINT&								location,
																	SKTRAN_GridIndex										profileidx,
																	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	coords, 
																	std::shared_ptr<const SKTRAN_RayFactory_Base>			diffuserayfactory,
																	const  SKTRAN_SpecsInternal_Diffuse_V21*				diffusespecs,
																	const  SKTRANSO_TableDiffusePoints*					parent
																 )
{
	double											radius;
	HELIODETIC_POINT								point;
	bool											ok;
	SKTRAN_GridIndex								pointidx;
	const SKTRAN_GridDefDiffuseHeights_V21*			diffuseheights;

	ReleaseResources();																				// release the current resources
	m_parent = parent;																				// Set the parent of this object
	m_location   = location;																		// Configure the cosine of the solar zenith angle of this profile
	m_profileidx = profileidx;
	NXASSERT(( (m_location.CosSZA() >= -1.0) && (m_location.CosSZA() <= 1.0) ));					// Make sure we are not being silly

	diffuseheights  = diffusespecs->DiffuseHeights ( profileidx );									// Get the heights of the diffuse points
	ok              = AllocatePoints( diffuseheights->NumAltitudes() );								// Allocate the points we need for this diffuse profile
	if (ok)																							// If that worked then
	{																								// proceed to initialize
		for (pointidx=0; pointidx < NumPoints(); pointidx++ )										// Now configure all of the points in the diffuse profile
		{																							// for each point
			radius = coords->AltitudeToRadius( diffuseheights->At(pointidx ) );						// Get the height of the diffuse altitude
			point.Initialize( m_location.LocalZenith(), radius, coords.get() );							// And set this diffuse point up so its at its own location but at the height of the scattering point
			ok = ok && m_diffusepoints[pointidx].ConfigureGeometry_Stage1(  diffuserayfactory,
																			diffusespecs,
																			point,
																			profileidx,
																			pointidx );
		}																							// do all of the points
	}	
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "SKTRAN_DiffuseProfileGeometry_V21::ConfigureGeometry, Error configuring geometry");
		ReleaseResources();
	}
	return ok;
}




