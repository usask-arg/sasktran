#include "../sasktranv21_internals.h"
#include "../sasktran_legacy21.h"



/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21		2010-4-29*/
/** The properties for the optical grid using teh legacy settings
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height : public SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21
{
	private:
	SKTRAN_GridDefOpticalPropertiesRadii_V21*					m_opticalpropheights;

	private:
		bool													ConfigureOpticalPropertyShells	( const std::vector<double>& heights );
	
	public:
																SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height( const std::vector<double>& heights );
		virtual												   ~SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height();

	public:
		virtual bool											CreateEmptyOpticalPropertiesTable( SKTRAN_TableOpticalProperties_V21**	 opticalpropertiestable ) const; 
		virtual const SKTRAN_GridDefOpticalPropertiesRadii_V21*	OpticalPropertiesGrid			 ()  const  { return m_opticalpropheights;}
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height		2007-11-20*/
/** A class used to cache the optical properties of the atmosphere. This class
 *	implements optical properties which are only a function of altitude.
 *	NOte the classes SKTRAN_QuadratureScatteringMatrixxxxxx are also affected by
 *	by this.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_TableOpticalProperties_1D_Height : public SKTRAN_TableOpticalProperties_V21
{
	private:
		double											m_wavelen;
		double											m_mjd;
		size_t											m_numheighttoindex;			//!< Number of elements in the height to index table (typically 10,000)
		size_t*											m_heighttoindextable;		//!< A table [m_numheighttoindex] that quickly converts a radius to an index with specified resolution
		double											m_heightindexresolution;	//!< The resolution o fthe radius to index rsolution, typically 10 meters
		double											m_minheight;

		std::vector<double>								m_scatextinction;			//!< JTW: Array [numshells] of scattering extinction .
		std::vector<double>								m_extinction;				//!< Array [numshells] of extinction .
//		nx1dArray<double>								m_effectiveextinction;		//!< Array [numshells] of effectiveextinction .
		nx2dArray<SKTRAN_PhaseMatrixScalar>				m_singleScatt;				//!< Array (numshells, numscatterangles) Indexed by shell and angle.
		const SKTRAN_GridDefOpticalPropertiesRadii_V21*	m_altitudegrid;				//!< The grid that defines the spatial grid for specifying optical properties
		const SKTRAN_GridDefScatterAngle_V21*					m_scatteranglegrid;			//!< The grid that defines the angular grid for specifying scattered rays.
		std::shared_ptr< const SKTRAN_CoordinateTransform_V2>	m_coordinates;				//!< 
		const skBRDF*								m_albedo;


	private:
		void										ReleaseResources			 ();
		void										ReleaseObjects				 ();
		bool										Allocate					 ( size_t numcells, size_t numangles );
		bool										ConfigureAltitudeToIndexTable();
		bool										IndexOfPointBelowOrEqual	 ( double h0, size_t* lowindex ) const;
		bool										IndexOfPointEqualOrAbove	 ( double h0, size_t* hihindex ) const;
		size_t										NumShells					 () const { return (m_altitudegrid     == NULL)? 0 : m_altitudegrid->NumAltitudes();}
		size_t										NumScatterAngles			 () const { return (m_scatteranglegrid == NULL)? 0 : m_scatteranglegrid->NumAngles();}
		double										GetExtinctionPerCM			 ( const HELIODETIC_POINT& location, const std::vector<double>& extinction  ) const;


	public:
													SKTRAN_TableOpticalProperties_1D_Height			();
		virtual 								   ~SKTRAN_TableOpticalProperties_1D_Height			();
		virtual bool								ConfigureGeometry					( const SKTRAN_SpecsInternal_Base* specs       ) override;
		virtual bool								ConfigureOptical					( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate ) override;
		virtual double								TotalExtinctionPerCM				( const HELIODETIC_POINT& point     ) const override;				//, double* shell_sigma_, double* shell_r, bool getshellbelow_radius ) const;
		virtual double								ScatteringExtinctionPerCM			( const HELIODETIC_POINT& location  )const override;
		virtual bool								GetScatteringCoefficientCM2			( const HELIODETIC_POINT& point, double cosangle,  SKTRAN_PhaseMatrixScalar* phasematrix ) const override;
		virtual bool								GetScatteringCoefficientCM2ForJindex( const SKTRANSO_JIndex* boundingpoints, size_t numpoints,  double cosangle, SKTRAN_PhaseMatrixScalar* coefficient ) const override;
		virtual bool								GetBoundingSpatialPoints			( const HELIODETIC_POINT& location, SKTRANSO_JIndex* interppoints, size_t* numpoints ) const override;
		virtual bool								GetBoundingScatteringPoints			( double cosscatteringangle, SKTRANSO_JIndex* interppoints, size_t* numpoints ) const override;
		virtual bool								GetBRDF								( const HELIODETIC_POINT& point, double mu_in, double mu_out, double cosdphi, double* brdf ) const override;
		virtual bool								GetBRDFGeodetic						( const GEODETIC_INSTANT& point, double mu_in, double mu_out, double cosdphi, double* brdf ) const override;

		virtual bool								GetEffectiveExtinctionPerCMWithHeight1( const HELIODETIC_POINT& point, double h0, double r0, double h1, double r1, double* sigma0, double* sigma1 ) const override;

		virtual bool								IsOptionTrue						( SKTRAN_TableOpticalProperties_V21::OPTIONSENUM options) const override;

//		virtual double								TotalExtinctionPerCM					( const SKTRAN_RayStorage_Base* r, size_t index ) const override                { nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_1D_Height, V3 interface not implemented in legacy."); return -9999.0;}
//		virtual bool								GetAlbedo								( const SKTRAN_RayStorage_Base* r, size_t index, double* albedo ) const override{ nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_1D_Height, V3 interface not implemented in legacy."); return false;}
//		virtual bool								GetEffectiveExtinctionPerCMWithHeight1	( const SKTRAN_RayStorage_Base* r, size_t startPtIndex, double* sigma0, double* sigma1 ) const override{ nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_1D_Height, V3 interface not implemented in legacy."); return false;}
//		virtual bool								GetScatteringCoefficientCM2				( const SKTRAN_RayStorage_Base* r, size_t index, double cosangle,  SKTRAN_PhaseMatrixScalar* phasematrix ) const override{ nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_1D_Height, V3 interface not implemented in legacy."); return false;}
//		virtual double								ScatteringExtinctionPerCM				( const SKTRAN_RayStorage_Base* r, size_t index ) const  override               { nxLog::Record(NXLOG_ERROR, "SKTRAN_TableOpticalProperties_1D_Height, V3 interface not implemented in legacy."); return -9999.0;}
        
        /* V3-only functions */
        bool                                        CreateInterpolationForPoint           ( const HELIODETIC_POINT& point, SKTRAN_TableOpticalProperties_PointCache& interpolator ) const override { nxLog::Record(NXLOG_WARNING, "SKTRAN_TableOpticalProperties_1D_Height::CreateInterpolationForPoint, V3 interface not implemented in legacy."); return false;}
		virtual double								TotalExtinctionPerCM					( double wavelength, const HELIODETIC_POINT& point) const override { return std::numeric_limits<double>::infinity(); }
		virtual bool								GetBRDF									( double wavelength, const HELIODETIC_POINT& point, double mu_in, double mu_out, double cosdphi, double* brdf) const override { return false; }
		virtual bool								GetBRDFGeodetic							( double wavelength, const GEODETIC_INSTANT& point, double mu_in, double mu_out, double cosdphi, double* brdf) const override { return false; }
		virtual bool								GetEffectiveExtinctionPerCMWithHeight1	( double wavelength, const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1) const override { return false; }
		virtual bool								GetEffectiveExtinctionPerCMWithHeight1	( double wavelength, const SKTRAN_RayStorage_Base* r, size_t startPtIndex, double* sigma0, double* sigma1) const override { return false; }
		virtual bool								GetScatteringCoefficientCM2				( double wavelength, const HELIODETIC_POINT& point, double cosangle, SKTRAN_PhaseMatrixScalar* scatcoeff) const override { return false; }
		virtual double								ScatteringExtinctionPerCM				( double wavelength, const HELIODETIC_POINT& point) const override { return std::numeric_limits<double>::infinity(); }
		virtual bool								GetScatteringMatrixCM2					( double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix) const override { return false; }
		virtual bool								GetResultOfUnpolarizedScatterCM2		( double wavelength, const HELIODETIC_POINT& point, double cosAngle, SKTRAN_Stokes_NC& stokesvec) const override { return false; }
		virtual bool								CreateInterpolationForPoint				( double wavelength, const HELIODETIC_POINT& point, SKTRAN_TableOpticalProperties_PointCache& interpolator) const override { return false; }


};




/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::SKTRAN_TableOpticalProperties_1D_Height		2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableOpticalProperties_1D_Height::SKTRAN_TableOpticalProperties_1D_Height()
{
	NXTRACE_ONCEONLY(firsttimeb,("**** 2008-12-23 ***** SKTRAN_TableOpticalProperties_1D_Height::SKTRAN_TableOpticalProperties_1D_Height, Needs work to get optical properties on Shell boundaries, CZR code is on Cell boundaries\n"));
	NXTRACE_ONCEONLY(firsttimec,("**** 2010-12-15 ***** SKTRAN_TableOpticalProperties_1D_Height::SKTRAN_TableOpticalProperties_1D_Height, NEED TO CHECK OUT ScatExtinctionPerCM, does it need to be virtual\n"));
	m_altitudegrid     = NULL;	
	m_scatteranglegrid = NULL;
	m_coordinates      = NULL;
	m_albedo           = NULL;
	m_mjd              = 0.0;
	m_minheight        = 9999999.0;
	
	m_numheighttoindex = 0;
	m_heighttoindextable = NULL;
	m_heightindexresolution = 10.0;


	NXTRACE(( "SKTRAN_TableOpticalProperties_1D_Height::Constructor, still need to implement the Update function\n"));
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::~SKTRAN_TableOpticalProperties_1D_Height		2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableOpticalProperties_1D_Height::~SKTRAN_TableOpticalProperties_1D_Height()
{
	NXTRACE(("SKTRAN_TableOpticalProperties_1D_Height::Destructor invoked\n"));
	ReleaseResources();
	ReleaseObjects();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::ReleaseResources		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_TableOpticalProperties_1D_Height::ReleaseObjects()
{
	if ( m_altitudegrid     != NULL) m_altitudegrid->Release();
	if ( m_scatteranglegrid != NULL) m_scatteranglegrid->Release();
	if ( m_albedo           != NULL) m_albedo->Release();
	m_altitudegrid     = NULL;
	m_scatteranglegrid = NULL;
	m_albedo           = NULL;
	m_extinction.clear();
//	m_effectiveextinction.erase();
	m_scatextinction.clear();
	m_singleScatt.erase();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::ReleaseResources		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_TableOpticalProperties_1D_Height::ReleaseResources()
{
	if (m_heighttoindextable != NULL) delete [] m_heighttoindextable;
	m_numheighttoindex   = 0;
	m_heighttoindextable = NULL;

	m_extinction.clear();
//	m_effectiveextinction.erase();
	m_scatextinction.clear();
	m_singleScatt.erase();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::Allocate		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height::Allocate( size_t numshells, size_t numangles )
{
	bool	ok;

	ReleaseResources();
	ok =      (numshells > 0) && (numangles > 0);
	/*ok = ok &&*/ m_extinction.resize         ( numshells );
//	ok = ok && m_effectiveextinction.SetSize( numshells);
	/*ok = ok &&*/ m_scatextinction.resize     ( numshells );
	ok = ok && m_singleScatt.SetSize        ( numshells, numangles );
	ok = ok && ConfigureAltitudeToIndexTable();
	if (!ok) ReleaseResources();
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::ConfigureAltitudeToIndexTable		2008-3-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height::ConfigureAltitudeToIndexTable()
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
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableOpticalProperties_1D_Height::ConfigureAltitudeToIndexTable, There is a logic error filling out the AltitudeToIndex table. It will need to be fixed");
	}
	hidx = m_altitudegrid->NumAltitudes();
	for (; idx < m_numheighttoindex; idx++) m_heighttoindextable[idx]  = hidx;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::IndexOfPointBelowOrEqual		2008-3-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height::IndexOfPointBelowOrEqual( double h0, size_t* lowindex ) const
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
 *					SKTRAN_TableOpticalProperties_1D_Height::IndexOfPointEqualOrAbove		2008-3-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height::IndexOfPointEqualOrAbove( double h0, size_t* hihindex ) const
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
 *					SKTRAN_TableOpticalProperties_1D_Height::GetExtinctionPerMeterWithHeight		2008-1-31*/
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

bool SKTRAN_TableOpticalProperties_1D_Height::GetEffectiveExtinctionPerCMWithHeight1( const HELIODETIC_POINT& location, double h0, double r0, double h1, double r1, double* sigma0, double* sigma1 ) const
{
	size_t	lowindex = 0;
	size_t	hihindex = 0;
	double	ru;
//	double	rd;
//	double	dr;
//	double	sigmad;
//	double	sigmau;
	bool	ok;
	bool	ok1;
	bool	ok2;

	if (r1 < r0)															// MAke sure r0 and r1 are in ascending order
	{
		ru = r0;
		r0 = r1;
		r1 = ru;

		ru = h0;
		h0 = h1;
		h1 = ru;
	}

	ok1 = IndexOfPointBelowOrEqual( h0, &lowindex );		// Get the optical properties shell below (or equal) the requested lower radius
	ok2 = IndexOfPointEqualOrAbove( h1, &hihindex );		// Get the optocal properties shell above (or equal) the requested upper radius
	ok  = ok1 && ok2;
	if (!ok)
	{
		*sigma0 = 0.0;
		*sigma1 = 0.0;
	}
	else
	{
		if (lowindex == hihindex)
		{
			*sigma0 = m_extinction[lowindex];			// Convert extinction per CM to extenction per meter.
			*sigma1  = *sigma0;
		}
		else
		{
//			rd     = m_coordinates->AltitudeToRadius( m_altitudegrid->At(lowindex) );
//			ru     = m_coordinates->AltitudeToRadius( m_altitudegrid->At(hihindex) );
//			dr     = ru - rd;
//			NXASSERT(( dr > 0.0)); 
			*sigma0 = m_extinction[lowindex];			// Convert extinction per CM to extenction per meter.
			*sigma1 = m_extinction[hihindex];
//			*sigmak = (sigmad*ru - sigmau*rd)/dr;
//			*sigmaf = (sigmau-sigmad)/dr;
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::GetExtinctionPerMeterWithHeight		2008-1-31*/
/** Getthe extinction per meter at altitude radius, r. This code performs log
 *	interpolation (as of April 25, 2014)
**/
/*---------------------------------------------------------------------------*/

double SKTRAN_TableOpticalProperties_1D_Height::GetExtinctionPerCM( const HELIODETIC_POINT& location, const std::vector<double>&	extinction  ) const
{
	size_t	lowindex;
	size_t	hihindex;
	double	ru;
	double	rd;
	double	sigmad;
	double	sigmau;
	double	sigma;
	bool	ok;
	bool	ok1;
	bool	ok2;
	double	r;
	;

	r = location.Altitude();
	ok1 = m_altitudegrid->IndexOfPointBelowOrEqual( r, &lowindex );						// Find the grid point below this radius
	ok2 = m_altitudegrid->IndexOfPointEqualOrAbove( r, &hihindex );						// find the grid point above this radius
	ok  = ok1 && ok2;													// Check that everythimg is good
	if (ok)																// it is
	{																	// so
		sigmad = extinction[lowindex];									// get extinction per CM
		if ( lowindex != hihindex)										// if we have different radii
		{																// then
			rd     = m_altitudegrid->At(lowindex);						// get the lower radius
			ru     = m_altitudegrid->At(hihindex);						// get the higher radius
			sigmau = extinction[hihindex];
			if ((sigmau <= 0.0) ||( sigmad <= 0.0))						// If either of our extinctions is identically 0 or less
			{
				sigma  = (sigmad + (sigmau-sigmad)*(r-rd)/(ru-rd));			// Do the linear interpolation
			}
			else															// Otherwise do log interpolation
			{
				sigmad = log(sigmad);
				sigmau = log(sigmau);							// Get extinction per CM
				sigma  = exp(sigmad + (sigmau-sigmad)*(r-rd)/(ru-rd));			// Do the linear interpolation
			}
		}																// and that is that
		else															// otherwise
		{																// upper and lower are the same index
			sigma = sigmad;												// so simpoly set it to the lower value
		}																// and that is that
	}																	// but
	else																// if we have an out of bounds condition
	{																	// then
		if (!ok1)														// if we failed to find a point below or equal
		{																// then
			sigma = 200.0;												// below ground so set optical depth per meter to a very high value		
			NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** SKTRAN_TableOpticalProperties_1D_Height::GetExtinctionPerMeter, point below ground\n"));
		}
		else															// otherwise
		{																// the point is above this grid
			sigma = 0.0;												// and signma is zero.
		}
	}
	return sigma;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::GetExtinctionPerMeterWithHeight		2008-1-31*/
/** Getthe extinction per meter at altitude radius, r.
**/
/*---------------------------------------------------------------------------*/

double SKTRAN_TableOpticalProperties_1D_Height::TotalExtinctionPerCM( const HELIODETIC_POINT& location  )const
{
	return GetExtinctionPerCM( location, m_extinction);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::ScatteringExtinctionPerCM		2010-5-31*/
/** Get the scattering extinction per centimeter at altitude radius, r.
**/
/*---------------------------------------------------------------------------*/

double SKTRAN_TableOpticalProperties_1D_Height::ScatteringExtinctionPerCM( const HELIODETIC_POINT& location  )const
{
	return GetExtinctionPerCM( location, m_scatextinction);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::Configure		2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height::ConfigureGeometry( const SKTRAN_SpecsInternal_Base* aspecs )
{

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::GetBoundingPoints		2009-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height::GetBoundingSpatialPoints( const HELIODETIC_POINT& location, SKTRANSO_JIndex* interppoints, size_t* numpoints ) const
{
	SKTRAN_GridIndex	uppercell;					// The upper shell that bounds radius
	SKTRAN_GridIndex	lowercell;					// The lower shell that bounds radius
	double				lowerw;						// The linear interpolating weight for the lower radius
	double				upperw;						// The linear interpolating weight for the upper radius
	bool				ok;
	size_t				idx;
	double				height;

	height = location.Altitude();
	ok  =  (*numpoints >= 2);
	ok  =  ok  && m_altitudegrid->FindBoundingIndices( height,   SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &lowercell,    &lowerw,   &uppercell,    &upperw );
	if (ok)
	{
		idx = 0;
		if (lowerw != 0.0) interppoints[idx++].ConfigureScatterMatrixTableIndex( 0, lowercell, 1.0, lowerw, 0, 1.0 );
		if (upperw != 0.0) interppoints[idx++].ConfigureScatterMatrixTableIndex( 0, uppercell, 1.0, upperw, 0, 1.0 );
		*numpoints = idx;
	}
	else
	{
		*numpoints = 0;
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableOpticalProperties_1D_Height::GetBoundingGridPoints, Error retrieving bounding grid points for height = %f", (double)height);
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::GetBoundingGridPoints		2009-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height::GetBoundingScatteringPoints( double cosscatteringangle, SKTRANSO_JIndex* interppoints, size_t* numpoints ) const
{
	SKTRAN_GridIndex	lowerscatidx;					// The lower scattering angle
	SKTRAN_GridIndex	upperscatidx;					// The upper scattering angle
	double				lowscatw;						// The linear interpolating weight for the lower radius
	double				upperscatw;						// The linear interpolating weight for the upper radius
	bool				ok;
	size_t				idx;

	ok  =  (*numpoints >= 2);
	ok  = ok && m_scatteranglegrid->FindBoundingIndices( cosscatteringangle, SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &lowerscatidx, &lowscatw, &upperscatidx, &upperscatw );
	if (ok)
	{
		idx = 0;
		if ( lowscatw   != 0.0) interppoints[idx++].ConfigureScatterMatrixTableIndex( 0, lowerscatidx, 1.0, lowscatw,   0, 1.0 );
		if ( upperscatw != 0.0) interppoints[idx++].ConfigureScatterMatrixTableIndex( 0, upperscatidx, 1.0, upperscatw, 0, 1.0 );
		*numpoints = idx;
	}
	else
	{
		*numpoints = 0;
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableOpticalProperties_1D_Height::GetBoundingScatteringPoints, Error retrieving bounding grid points for cos(scattering angle) = %f", (double)cosscatteringangle);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::GetScatteringCoefficient		2008-1-10*/
/** Calculates the scattering coefficient for a given altitude (expressed as a
 *	radius) and scattering angle (expressed as cos scatter angle).  The code
 *	performs a linear interpolation of the scattering array both in altitude
 *	and in cosangle space.
 *
 *	The function will return the coefficient as 0 if the code fails to properly
 *	interpolate the table. This normally occurs if the tables do not properly
 *	span the range of values used by the radiative transfer code.  In particular,
 *	the scattering angle table should probably have bounding values defined for scattering 
 *	angles of 0 degrees (cosscat = 1) and 180 degrees (cosscat = -1).
 **/
/*---------------------------------------------------------------------------*/


bool SKTRAN_TableOpticalProperties_1D_Height::GetScatteringCoefficientCM2( const HELIODETIC_POINT& location, double cosAngle, SKTRAN_PhaseMatrixScalar* coefficient ) const
{
	bool				ok;
	SKTRAN_GridIndex	uppercell;					// The upper shell that bounds radius
	SKTRAN_GridIndex	lowercell;					// The lower shell that bounds radius
	SKTRAN_GridIndex	upperscatidx;				// The upper scattering angle that bounds cosAngle
	SKTRAN_GridIndex	lowerscatidx;				// The lower scattering angle that bounds cosAngle
	double				lowradiusscat;				// Scattering on lower shell
	double				hihradiusscat;				// Scattering on upper shell
	double				lowerw;						// The linear interpolating weight for the lower radius
	double				upperw;						// The linear interpolating weight for the upper radius
	double				lowscatw;					// The linear interpolating weight for the lower scattering angle
	double				upperscatw;						// The linear interpolating weight for the upper scattering angle
	double				height;

	height = location.Altitude();
	ok  =       m_altitudegrid->FindBoundingIndices    ( height,   SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &lowercell,    &lowerw,   &uppercell,    &upperw );
	ok  = ok && m_scatteranglegrid->FindBoundingIndices( cosAngle, SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &lowerscatidx, &lowscatw, &upperscatidx, &upperscatw );
	if (ok)
	{
		lowradiusscat = m_singleScatt.At( lowercell, lowerscatidx)*lowscatw + m_singleScatt.At( lowercell, upperscatidx)*upperscatw;
		hihradiusscat = m_singleScatt.At( uppercell, lowerscatidx)*lowscatw + m_singleScatt.At( uppercell, upperscatidx)*upperscatw;
		*coefficient  = SKTRAN_DBL_TO_STOKES_SCALAR( (lowradiusscat*lowerw + hihradiusscat*upperw) );	
	}
	else
	{
		NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** SKTRAN_TableOpticalProperties_1D_Height::GetScatteringCoefficient, Bounding indices failed\n"));
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableOpticalProperties_1D_Height::GetScatteringCoefficient, failed finding index for height (%10g) meters, cosangle (%g).  Make sure the optical properties grid spans above the ray tracing and diffuse grids", (double) height, (double)cosAngle );
		*coefficient = SKTRAN_DBL_TO_STOKES_SCALAR(0);
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::GetScatteringCoefficientCM2		2009-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height::GetScatteringCoefficientCM2ForJindex( const SKTRANSO_JIndex* boundingpoints, size_t numpoints, double cosangle, SKTRAN_PhaseMatrixScalar* coefficient ) const
{
	size_t				pidx;
	size_t				cellidx;
	double				cellw;
	double				c;
	double				lowscatw;					// The linear interpolating weight for the lower scattering angle
	double				upperscatw;						// The linear interpolating weight for the upper scattering angle
	SKTRAN_GridIndex	lowerscatidx;
	SKTRAN_GridIndex	upperscatidx;
	bool				ok;


	ok      = m_scatteranglegrid->FindBoundingIndices( cosangle, SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &lowerscatidx, &lowscatw, &upperscatidx, &upperscatw );
	cellidx = boundingpoints->HeightIndex();
	cellw   = boundingpoints->VertexWeight();
	c       = (m_singleScatt.At(cellidx,lowerscatidx)*lowscatw + m_singleScatt.At(cellidx,upperscatidx)*upperscatw)*cellw;
	if (numpoints > 1)
	{
		for (pidx = 1; pidx < numpoints; pidx++ )
		{
			cellidx = boundingpoints[pidx].HeightIndex();
			cellw   = boundingpoints[pidx].VertexWeight();
			c      += (m_singleScatt.At(cellidx,lowerscatidx)*lowscatw + m_singleScatt.At(cellidx,upperscatidx)*upperscatw)*cellw;
		}
	}
	*coefficient = SKTRAN_DBL_TO_STOKES_SCALAR(c);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::ConfigureOptical		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height::ConfigureOptical( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate )
{
	SKTRAN_GridIndex		shellidx;
	SKTRAN_GridIndex		scatidx;
	GEODETIC_INSTANT		point;
	bool					ok;
	bool					ok1;
	bool					ok2;
	double					scatter;
	double					cosangle;
	skBRDF*				albedo;

	NXASSERT(( wavelen > 0.0 ));

	ok = opticalstate.GetAlbedoObject( &albedo );
	if (ok)
	{
		if (albedo   != NULL) albedo->AddRef(); 
		if (m_albedo != NULL) m_albedo->Release();
		m_albedo        = albedo;
		
		m_wavelen       = wavelen;
		point = opticalstate.GetTimeAndLocation();
		m_mjd = point.mjd;
		NXASSERT(( m_mjd > 10000.0 ));

		ok =  opticalstate.SetWavelength (m_wavelen);													// Set the wavelength in the species
		for (shellidx = 0; shellidx < m_altitudegrid->NumAltitudes(); shellidx++)
		{
			ok1 = true;
			point.heightm = m_altitudegrid->At( shellidx );
			ok1 = opticalstate.SetTimeAndLocation( point, false );	
			if (ok1)
			{
				m_extinction.at(shellidx)          = opticalstate.ExtinctionPercm();					
				m_scatextinction.at(shellidx)	   = opticalstate.ScatteringPercm();
				ok1  = ok1 && (m_extinction.at(shellidx) >= 0.0);
				ok1  = ok1 && (m_scatextinction.at(shellidx) >= 0.0);
				NXASSERT((m_scatextinction.at(shellidx)>= 0.0));
				NXASSERT((m_extinction.at(shellidx)>= 0.0));
//				m_effectiveextinction.At(shellidx) = opticalstate.EffectiveExtinctionPercm ();

				for (scatidx = 0; scatidx < m_scatteranglegrid->NumAngles(); scatidx++)
				{
					cosangle   = m_scatteranglegrid->At( scatidx );
					ok2        = opticalstate.ScalarScatteringCoefficient( cosangle, &scatter );
					m_singleScatt.At( shellidx, scatidx ) =  (SKTRAN_PhaseMatrixScalar)scatter;
					ok2     = ok2 && (scatter >= 0.0);
					ok1     = ok1 && ok2;
					NXASSERT((m_singleScatt.At( shellidx, scatidx ) >= 0.0));
				}
			}
			ok = ok && ok1;
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableOpticalProperties_1D_Height::ConfigureOptical, Error configuring the optical state");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height::GetAlbedo		2008-3-20*/
/** Returns the albedo at the specified point**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height::GetBRDF( const HELIODETIC_POINT& userpoint, double mu_in, double mu_out, double cosdphi, double* brdf ) const
{
	GEODETIC_INSTANT		point;
	bool					ok;

	if (m_albedo == NULL)
	{
		*brdf = 0.0;
		ok      = true;
	}
	else
	{
		point = m_coordinates->PointToGeodetic( userpoint, m_mjd );
		ok = m_albedo->BRDF( m_wavelen, point, mu_in, mu_out, cosdphi, brdf );
	}
	return ok;
}

bool SKTRAN_TableOpticalProperties_1D_Height::GetBRDFGeodetic(const GEODETIC_INSTANT& userpoint, double mu_in, double mu_out, double cosdphi, double* brdf) const
{
	bool					ok;

	if (m_albedo == NULL)
	{
		*brdf = 0.0;
		ok = true;
	}
	else
	{
		ok = m_albedo->BRDF(m_wavelen, userpoint, mu_in, mu_out, cosdphi, brdf);
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height::IsOptionTrue		2010-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height::IsOptionTrue( SKTRAN_TableOpticalProperties_V21::OPTIONSENUM options) const
{
	return (options == OPTION_IS_HORIZONTALLY_UNIFORM);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height::SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height		2010-4-29*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height::SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height( const std::vector<double>& heights  )
{
	m_opticalpropheights = new SKTRAN_GridDefOpticalPropertiesRadii_V21;		// create the blank ray tracing  grid
	m_opticalpropheights->AddRef();
	ConfigureOpticalPropertyShells( heights  );

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height::~SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height		2010-4-29*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height::~SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height()
{
	if (m_opticalpropheights != NULL) m_opticalpropheights->Release();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridSpecificationsLegacy_V2::ConfigureOpticalPropertyShells		2008-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height::ConfigureOpticalPropertyShells( const std::vector<double>& heights  )
{
	size_t  i;
	bool	ok;
	size_t	npts;

	npts = heights.size();
	ok = (m_opticalpropheights->AllocateGridArray(npts));
	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_GridSpecificationsLegacy_V2::MakeOpticalRadiiGrid, Error allocating memory for optical height buffer");
	}
	else
	{
		for (i = 0; i < npts; i++ )
		{
			m_opticalpropheights->AtVar(i) = heights[i];
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height::CreateOpticalPropertiesTable		2010-5-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height::CreateEmptyOpticalPropertiesTable ( SKTRAN_TableOpticalProperties_V21** opticalpropertiestable ) const
{
	SKTRAN_TableOpticalProperties_1D_Height*	table;
	bool										ok;

	table = new SKTRAN_TableOpticalProperties_1D_Height;
	ok    = (table != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height::CreateOpticalPropertiesTable, Error createing otpical properties table");
	}
	else
	{
		table->AddRef();
	}
	*opticalpropertiestable = table;
	return ok;
}
	

/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height::ConfigureOpticalPropertyShells		2010-5-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height::ConfigureOpticalPropertyShells( const double* alts, size_t npts)
{
	m_heights.resize(npts);
	std::copy( alts, alts+npts, m_heights.begin() );
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height::CreateInternalSpecs		2010-5-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height::CreateInternalSpecs( const SKTRAN_SpecsInternal_OpticalPropertiesGrid_V21** userinternalspecs ) const
{
	SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height*	internalspecs = NULL;
	bool													ok;	

	ok = (m_heights.size() > 0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height::CreateInternalSpecs, cannot create internal specs as the optical height grid is empty");
	}
	else
	{
		internalspecs = new SKTRAN_SpecsInternal_OpticalPropertiesGrid_1D_Height(m_heights);
		ok            = (internalspecs != NULL);
		if (!ok)
		{
			nxLog::Record( NXLOG_WARNING,"SKTRAN_SpecsUser_OpticalPropertiesGrid_1D_Height::CreateInternalSpecs, Error allocating memory");
		}
		else
		{
			internalspecs->AddRef();
		}
	}
	*userinternalspecs = internalspecs;
	return ok;

}





	

