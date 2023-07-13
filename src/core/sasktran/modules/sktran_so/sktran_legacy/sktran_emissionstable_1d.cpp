#include "../sasktranv21_internals.h"
#include "../sasktran_legacy21.h"
#include "sktran_legacy_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::SKTRANSO_InternalEmissionPropertiesTable_1D_Height		2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_InternalEmissionPropertiesTable_1D_Height::SKTRANSO_InternalEmissionPropertiesTable_1D_Height()
{
	m_isempty               = true;
	m_wavelen               = 0.0;
	m_mjd                   = 0.0;
	m_altitudegrid          = NULL;	
	m_heighttoindextable    = NULL;
	m_minheight             = 9999999.0;	
	m_numheighttoindex      = 0;
	m_heightindexresolution = 200.0;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::~SKTRANSO_InternalEmissionPropertiesTable_1D_Height		2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_InternalEmissionPropertiesTable_1D_Height::~SKTRANSO_InternalEmissionPropertiesTable_1D_Height()
{
	ReleaseResources();
	ReleaseObjects();
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ReleaseResources		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ReleaseObjects()
{
	if ( m_altitudegrid     != NULL) m_altitudegrid->Release();
	m_altitudegrid     = NULL;
	m_radiance.clear();
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ReleaseResources		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ReleaseResources()
{
	if (m_heighttoindextable != NULL) delete [] m_heighttoindextable;
	m_numheighttoindex   = 0;
	m_heighttoindextable = NULL;
	m_radiance.clear();
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::Allocate		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_InternalEmissionPropertiesTable_1D_Height::Allocate( size_t numshells )
{
	bool	ok;

	ReleaseResources();
	ok =      (numshells > 0);
	if (ok)
	{
		m_radiance.resize ( numshells );
		ok = ok && ConfigureAltitudeToIndexTable();
	}
	if (!ok) ReleaseResources();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ConfigureAltitudeToIndexTable		2008-3-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ConfigureAltitudeToIndexTable()
{
	double	maxalt;
	double	nexth;
	double	h;
	size_t	hidx;
	size_t	idx;
	bool	ok;

	NXASSERT(( m_heighttoindextable == NULL ));

	maxalt                = m_altitudegrid->back();
	m_minheight           = m_altitudegrid->front();
	m_numheighttoindex    = size_t( (maxalt-m_minheight)/m_heightindexresolution ) + 1;
	m_heighttoindextable  = new size_t[m_numheighttoindex];
	h     = m_minheight;
	idx   = 0;
	for (hidx = 1; hidx < m_altitudegrid->NumAltitudes(); hidx++ )
	{
		nexth = m_altitudegrid->At(hidx);
		while ( h <= nexth )
		{
			m_heighttoindextable[idx++]  = hidx;
			h                           += m_heightindexresolution;
		}
	}
	NXASSERT(( idx <= m_numheighttoindex ));
	ok = (idx <= m_numheighttoindex);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ConfigureAltitudeToIndexTable, There is a logic error filling out the AltitudeToIndex table. It will need to be fixed");
	}
	hidx = m_altitudegrid->NumAltitudes();
	for (; idx < m_numheighttoindex; idx++) m_heighttoindextable[idx]  = hidx;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::IndexOfPointBelowOrEqual		2008-3-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_InternalEmissionPropertiesTable_1D_Height::IndexOfPointBelowOrEqual( double h0, size_t* lowindex ) const
{
	size_t	indexabove;
	size_t	i;
	bool	ok;


	i  =  (size_t)( (h0-m_minheight)/m_heightindexresolution);
	ok = (i <= m_numheighttoindex);
	if (ok)
	{
		indexabove = m_heighttoindextable[i];								// Get the index of the shell above this ray
		if ( m_altitudegrid->At(indexabove) != h0) indexabove--;			// If we dont have an exact match then shell below is one less
		*lowindex  = indexabove;											// copy it over
		NXASSERT(( indexabove < m_altitudegrid->NumAltitudes() ));			// and make sure we are good
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::IndexOfPointEqualOrAbove		2008-3-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_InternalEmissionPropertiesTable_1D_Height::IndexOfPointEqualOrAbove( double h0, size_t* hihindex ) const
{
	size_t	indexabove;
	size_t	i;
	bool	ok;


	i =  (size_t)( (h0-m_minheight)/m_heightindexresolution);
	ok = (i < m_numheighttoindex);
	if (ok)
	{
		indexabove = m_heighttoindextable[i];
		ok         = (indexabove > 0);
		*hihindex  = indexabove;
		NXASSERT(( indexabove < m_altitudegrid->NumAltitudes() )); 
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::GetIsotropicRadiance		2008-1-31*/
/** Get the effective extinction per meter as a linear function of altitude between
 *	the two altitude radii, r0 and r1.  The effective extinctions has an adjustment to the
 *	extinction cross-section that accounts for any strong delta function forward scatter.
 *	Note that the code finds the
 *	the radii in the optical properties radial altitiude grid just above and below
 *	these two radii. These point mat not necessarily be contiguous in the optical properties
 *	radial altitude grid.
 *
 *	returns the value of sigma(r) such that simga(r) = sigmaK + r*sigmaF
**/
/*---------------------------------------------------------------------------*/

double SKTRANSO_InternalEmissionPropertiesTable_1D_Height::GetIsotropicRadianceInAtmosphere( const HELIODETIC_POINT& location )const
{
	size_t	lowindex;
	size_t	hihindex;
	double	h;
	double	radiance = 0.0;
	bool	ok;
	bool	ok1;
	bool	ok2;
	double	hd, hu;
	double	r[2];

	h   = location.Altitude();
	ok1 = m_altitudegrid->IndexOfPointBelowOrEqual( h, &lowindex );		// Find the grid point below this radius
	ok2 = m_altitudegrid->IndexOfPointEqualOrAbove( h, &hihindex );		// find the grid point above this radius
	ok  = ok1 && ok2;
	if (!ok)
	{
		radiance = 0.0;
	}
	else
	{
		if (lowindex == hihindex)
		{
			radiance = m_radiance[lowindex];			// Convert extinction per CM to extenction per meter.
		}
		else
		{
			hd       = m_altitudegrid->At(lowindex);
			hu       = m_altitudegrid->At(hihindex);
			r[0]     = m_radiance[lowindex];			// Convert extinction per CM to extenction per meter.
			r[1]     = m_radiance[hihindex];
			radiance = nxLinearInterpolate::FromTwoPoints(h, hd,hu, r );
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::GetIsotropicGroundRadiance		 2015- 3- 3*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRANSO_InternalEmissionPropertiesTable_1D_Height::GetIsotropicGroundRadiance( const HELIODETIC_POINT& location )const
{
	return m_groundradiance;
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::Configure		2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ConfigureGeometry( const SKTRAN_SpecsInternal_Base* aspecs )
{
	bool								ok;
	size_t								numshells;
	const SKTRAN_SpecsInternal_V21*		specs;

	m_wavelen = 0.0;
	specs     = dynamic_cast<const SKTRAN_SpecsInternal_V21*>(aspecs);	
	ReleaseObjects();
	SetCoords( specs->CoordinateSystemObject() );
	m_location = specs->CoordinateSystemPtr()->ReferencePoint(0.0);
	m_mjd      = specs->CoordinateSystemPtr()->ReferencePointMJD();
	m_altitudegrid  = specs->OpticalTableSpecs()->OpticalPropertiesGrid();
	ok =    (m_altitudegrid     != NULL);
	if (ok)
	{
		m_altitudegrid->AddRef();
		numshells  = m_altitudegrid->NumAltitudes();
		ok        =  Allocate( numshells );
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ConfigureOptical		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ConfigureOptical( double wavelen, const SKTRAN_CoordinateTransform_V2* coords, SKTRAN_AtmosphericEmission* emission )
{
	SKTRAN_GridIndex		shellidx;
	GEODETIC_INSTANT		point;
	bool					ok;
	bool					ok1;

	ok = (emission == nullptr);
	if (ok)
	{
		SetEmpty(true);
	}
	else
	{
		SetEmpty(false);
		NXASSERT(( wavelen > 0.0 ));

		point.latitude  = coords->ReferencePtLatitude();
		point.longitude = coords->ReferencePtLongitude();
		point.heightm   = 0.0;
		point.mjd       = coords->ReferencePointMJD();			// Reference point MJD is now a weighted mean, weighted towards target altitude, it makes a tiny difference due to floating roundoff

		m_wavelen = wavelen;
		m_mjd     = point.mjd;
		NXASSERT(( m_mjd > 10000.0 ));

		point.heightm = m_altitudegrid->At( 0 );
		ok =  emission->SetWavelength (m_wavelen);										// Set the wavelength in the species
		ok1 = emission->SetTimeAndLocation( point, true, false );						// Calculate the isotropic emission radiance from the ground itself (thermal emissions have a discontinuity at the ground for example)
		m_groundradiance = emission->IsotropicRadiance();								// And save the ground emission
		m_radiance.resize(m_altitudegrid->NumAltitudes());
		
		for (shellidx = 0; shellidx < m_altitudegrid->NumAltitudes(); shellidx++)		// NOw do the atmopsheric emission signals.
		{
			ok1 = true;
			point.heightm = m_altitudegrid->At( shellidx );
			ok1 = emission->SetTimeAndLocation( point, false, false );	
			m_radiance.at(shellidx)  = emission->IsotropicRadiance();
			ok = ok && ok1;
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ConfigureOptical, Error configuring the optical state");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::InterpolateTable		 2015- 3- 3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRANSO_InternalEmissionPropertiesTable_1D_Height::InterpolateTable( const HELIODETIC_POINT&		location,				// Location where interpolation required
                                                         const HELIODETIC_UNITVECTOR&	look,					// Look direction where interpolation required
												         SKTRANSO_JIndex*				vertexdescriptortable, 
														 size_t							maxpts, 
														 size_t*						npts,
														 double							weight ) const
{
	SKTRAN_GridIndex	altindex	[2];	// The indices of the radii in each diffuse profile that bound the altitude
	double				altweight   [2];	// The weights for each of the radii
	size_t				numvertex;									
	size_t				idx;
	bool				isground;

	NXASSERT(( weight == 1.0 ));
	NXASSERT(( vertexdescriptortable != NULL ));
	NXASSERT(( npts != NULL ));
//	NXASSERT(( location.Altitude() < 250000.0 ));
	NXASSERT((m_altitudegrid != nullptr));

	idx      = 0;
	isground = (!NXFINITE(look.X())) && (!NXFINITE(look.Y())) && (!NXFINITE(look.Z()));
	if (isground)
	{
		vertexdescriptortable[idx++].ConfigureEmissionTableIndex( 0, 0, 1.0, 1.0, true);
	}
	else
	{
		numvertex = m_altitudegrid->FindingBoundingIndices(  location.Altitude(), SKTRAN_GridDefBase_V2::OUTOFBOUND_ZERO, altindex, altweight, 2);
		for (size_t i = 0; i < numvertex; i++)
		{
			if ( altweight[i] != 0.0 ) vertexdescriptortable[idx++].ConfigureEmissionTableIndex( 0, altindex[i], 1.0, altweight[i], false);
		}
	}
	*npts  = idx;
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ConvertJIndexToRadiancePtr		 2015- 3- 3*/
/** **/
/*---------------------------------------------------------------------------*/

const SKTRAN_StokesScalar*	SKTRANSO_InternalEmissionPropertiesTable_1D_Height::ConvertJIndexToRadiancePtr( const SKTRANSO_JIndex*		entry, 
																								ENUM_SKTRAN_JSOURCE			jsource ) const
{
	size_t						index;
	const SKTRAN_StokesScalar*	JPtr = nullptr;

	if ( entry->IsEmissionGround() )
	{
		JPtr  = &m_groundradiance;
	}
	else
	{
		index  = entry->HeightIndex();
		JPtr   = &m_radiance.at(index);
	}
	return JPtr;
}

