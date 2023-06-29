#include "../sktran_common.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_Base::SKTRAN_TableOpticalProperties_Base		2013-10-16*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableOpticalProperties_Base::SKTRAN_TableOpticalProperties_Base()
{

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_Base::~SKTRAN_TableOpticalProperties_Base		2013-10-16*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableOpticalProperties_Base::~SKTRAN_TableOpticalProperties_Base()
{
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_Base::SetCoords		2013-10-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_Base::SetCoords( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coordinates)
{
	m_coordinates = coordinates;
	return true;
}


void SKTRAN_TableOpticalProperties_Base::SetPolarizationProperties ( std::unique_ptr< SKTRAN_PolarizationProperties_Base >& polprops )
{
    m_scatprops = std::move(polprops);
}


/*-----------------------------------------------------------------------------
 *					GroundBRDFAngles		 2016- 12- 22*/
/** Static member function to calculate the angles required for BRDF
 *	calculations.
 *
 *	\param point
 *		The location of the ground proint in Heliodetic coordinates
 *
 *	\param incomingray
 *		The unit vector of the incoming ray. This should be the direction away from the ground point.
 *
 *	\param outgoingray
 *		The unit vector o fthe outgoing ray. This should be the direction away from the ground point.
 *
 *	\param mu_in
 *		Returns the cosine of the  zenith angle of the incoming  ray.
 *
 *	\param mu_out
 *		Returns the cosine of the zenith angle of the outbound/scattered ray
 *
 *	\param cosdphi
 *		Returns the cosine of the azimuthal angle between the incoming and outbound ray.

 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_Base::GroundBRDFAngles( const HELIODETIC_POINT& point, const HELIODETIC_UNITVECTOR& incomingray, const HELIODETIC_UNITVECTOR& outgoingray, double* mu_in, double* mu_out, double* cosdphi )
{
	bool ok = true;
	HELIODETIC_UNITVECTOR		 localdirs[3];								// Gets local direction unit vectors: north, west, up
	const HELIODETIC_UNITVECTOR& upH    = point.LocalZenith();
	const HELIODETIC_UNITVECTOR& indirH  = incomingray;
	const HELIODETIC_UNITVECTOR& outdirH = outgoingray;

	point.LocalUnitVectors( localdirs, N_ELEMENTS(localdirs) );

	nxVector					 up		(     upH.X(),     upH.Y(),     upH.Z() );
	nxVector					 indir	(  indirH.X(),  indirH.Y(),  indirH.Z() );
	nxVector					 outdir ( outdirH.X(), outdirH.Y(), outdirH.Z() );	
	nxVector					 north  ( localdirs[0].X(), localdirs[0].Y(), localdirs[0].Z() );	// This not true geographic north but is the north which is towards the sun direction 
	nxVector					 west   ( localdirs[1].X(), localdirs[1].Y(), localdirs[1].Z() );	// NOt true geographic west but 90 degrees from the sun (ie either dawn or dusk not sure which one)


	*mu_in =  (up & indir);
	*mu_out = (up & outdir);
	if (*mu_in < 0.0) { 
		ok = false; 
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableOpticalProperties_Base::GroundBRDFAngles, Zenith angle of incoming ray is greater than 90 degrees (cos(zenith angle)=%e), i.e. it is \"coming into the point\". It should always be away from the point.", (double) *mu_in);
	}
	if (*mu_out < 0.0) { 
		ok = false; 
	nxLog::Record(NXLOG_WARNING,"SKTRAN_TableOpticalProperties_Base::GroundBRDFAngles, Zenith angle of outgoing ray is greater than 90 degrees (cos(zenith angle)=%e), i.e. it is \"coming into the point\". It should always be away from the point.", (double) *mu_out);
	}

	double inx    = indir & north;
	double iny    = indir & west;
	double inazi  = atan2( iny,inx);
	double outx   = outdir & north;
	double outy   = outdir & west;
	double outazi = atan2( outy, outx);
	
	*cosdphi = cos(outazi-inazi);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_Base::Get_AlbedoForDeprecatedLegacyCode		 2016- 12- 23*/
/** A little hack that enables us to use legacy code that still uses albedo instead 
 *	of the the full-up BRDF
 **/
/*---------------------------------------------------------------------------*/


bool SKTRAN_TableOpticalProperties_Base::Get_AlbedoForDeprecatedLegacyCode( const HELIODETIC_POINT& point, double* albedo ) const
{
	bool ok;

	ok = GetBRDF( point, 0, 0, 0,  albedo );
	*albedo *= nxmath::Pi;
	return ok;
}




/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::SKTRAN_TableOpticalProperties_1D_Height_V3		2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableOpticalProperties_1D_Height_V3::SKTRAN_TableOpticalProperties_1D_Height_V3()
{
	NXTRACE_ONCEONLY(firsttimeb,("**** 2008-12-23 ***** SKTRAN_TableOpticalProperties_1D_Height_V3::SKTRAN_TableOpticalProperties_1D_Height_V3, Needs work to get optical properties on Shell boundaries, CZR code is on Cell boundaries\n"));
	NXTRACE_ONCEONLY(firsttimec,("**** 2010-12-15 ***** SKTRAN_TableOpticalProperties_1D_Height_V3::SKTRAN_TableOpticalProperties_1D_Height_V3, NEED TO CHECK OUT ScatExtinctionPerCM, does it need to be virtual\n"));
	m_altitudegrid     = NULL;	
	m_scatteranglegrid = NULL;
	m_albedo           = NULL;
//	m_scatextinction   = new nx1dArray<double>*;
	m_scatextinction   = new std::vector<double>*;	// Should add some checks in here, and this should move to Allocate
	(*m_scatextinction)= NULL;
//	m_extinction	   = new nx1dArray<double>*;
	m_extinction	   = new std::vector<double>*;
	(*m_extinction)	   = NULL;
	m_mjd              = 0.0;
	m_minheight        = 9999999.0;
	M_GROUNDTOLERANCE  = -0.01;
	
	m_numheighttoindex = 0;
	m_heighttoindextable = NULL;
	m_heightindexresolution = 10.0;


	NXTRACE(( "SKTRAN_TableOpticalProperties_1D_Height_V3::Constructor, still need to implement the Update function\n"));
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::~SKTRAN_TableOpticalProperties_1D_Height_V3		2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableOpticalProperties_1D_Height_V3::~SKTRAN_TableOpticalProperties_1D_Height_V3()
{
	NXTRACE(("SKTRAN_TableOpticalProperties_1D_Height_V3::Destructor invoked\n"));
	ReleaseResources();
	ReleaseObjects();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::ReleaseResources		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_TableOpticalProperties_1D_Height_V3::ReleaseObjects()
{
	if ( m_altitudegrid     != NULL) m_altitudegrid->Release();
	if ( m_scatteranglegrid != NULL) m_scatteranglegrid->Release();
	if ( m_albedo           != NULL) m_albedo->Release();
	m_altitudegrid     = NULL;
	m_scatteranglegrid = NULL;
	m_albedo           = NULL;

	if(NULL!=m_extinction)		{ delete (*m_extinction);		(*m_extinction)=NULL;}
	if(NULL!=m_scatextinction)	{ delete (*m_scatextinction);	(*m_scatextinction)=NULL;}
	//if(NULL!=m_singleScatt)		{ delete (*m_singleScatt);		(*m_singleScatt)=NULL;}
	//std::vector<double>().swap(m_singleScatt);

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::ReleaseResources		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_TableOpticalProperties_1D_Height_V3::ReleaseResources()
{
	if (m_heighttoindextable != NULL) delete [] m_heighttoindextable;
	m_numheighttoindex   = 0;
	m_heighttoindextable = NULL;

	if(NULL!=m_extinction)	  { delete (*m_extinction); (*m_extinction)=NULL; delete m_extinction; m_extinction=NULL;}
	if(NULL!=m_scatextinction){ delete (*m_scatextinction); (*m_scatextinction)=NULL; delete m_scatextinction; m_scatextinction=NULL;}
	//if(NULL!=m_singleScatt)	  { delete (*m_singleScatt); (*m_singleScatt)=NULL; delete m_singleScatt; m_singleScatt=NULL;}
	//std::vector<double>().swap(m_singleScatt);

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::Allocate		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height_V3::Allocate( size_t numshells, size_t numangles )
{
	bool	ok;

	ReleaseResources();
	ok =      (numshells > 0) && (numangles > 0);

	// Allocate extinction tables
	//m_extinction = new nx1dArray<double>*; 
	//ok = ok && m_extinction!=NULL;
	//(*m_extinction)  = new nx1dArray<double>;
	//ok = ok && (*m_extinction) != NULL;
	//ok = ok && (*m_extinction)->SetSize( numshells );
	m_extinction = new std::vector<double>*; 
	ok = ok && m_extinction!=NULL;
	if(ok) (*m_extinction)  = new std::vector<double>;
	ok = ok && (*m_extinction) != NULL;
	if(ok) (*m_extinction)->resize( numshells );
	
	// Allocate scatter coefficient tables
	//m_scatextinction = new nx1dArray<double>*; ok = ok && m_scatextinction!=NULL;
	//(*m_scatextinction) = new nx1dArray<double>;
	//ok = ok && (*m_scatextinction) != NULL;
	//ok = ok && (*m_scatextinction)->SetSize     ( numshells );
	m_scatextinction = new std::vector<double>*;
	ok = ok && m_scatextinction!=NULL;
	if(ok) (*m_scatextinction) = new std::vector<double>;
	ok = ok && (*m_scatextinction) != NULL;
	if(ok) (*m_scatextinction)->resize ( numshells );

	// Allocate scatter pdf array
	//m_singleScatt = new nx2dArray<SKTRAN_PhaseMatrixScalar>*; ok = ok && m_singleScatt!=NULL;
	//(*m_singleScatt) = new nx2dArray<SKTRAN_PhaseMatrixScalar>;
	//ok = ok && (*m_singleScatt) != NULL;
	//ok = ok && (*m_singleScatt)->SetSize        ( numshells, numangles );
	//m_singleScatt.resize( numshells*numangles );

    m_scatprops->Allocate( numshells*numangles );

	// Check if everything's been done okay
	ok = ok && ConfigureAltitudeToIndexTable();
	if (!ok) ReleaseResources();

	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::ConfigureAltitudeToIndexTable		2008-3-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height_V3::ConfigureAltitudeToIndexTable()
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
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableOpticalProperties_1D_Height_V3::ConfigureAltitudeToIndexTable, There is a logic error filling out the AltitudeToIndex table. It will need to be fixed");
	}
	hidx = m_altitudegrid->NumAltitudes();
	for (; idx < m_numheighttoindex; idx++) m_heighttoindextable[idx]  = hidx;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::IndexOfPointBelowOrEqual		2008-3-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height_V3::IndexOfPointBelowOrEqual( double h0, size_t* lowindex ) const
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
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::IndexOfPointEqualOrAbove		2008-3-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height_V3::IndexOfPointEqualOrAbove( double h0, size_t* hihindex ) const
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


bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetAlbedo ( const SKTRAN_RayStorage_Base* r, size_t index, double* albedo ) const
{
	bool ok = true;

	nxLog::Record( NXLOG_ERROR, "SKTRAN_TableOpticalProperties_1D_Height_V3::Get albedo, Need to let minimum container know about albedo object type.");
	ok = false;

	return ok;
}


bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetEffectiveExtinctionPerCMWithHeight1( const SKTRAN_RayStorage_Base* r, size_t startPtIndex, double* sigma0, double* sigma1 ) const
{
	//double sigmad;
	//double sigmau;
	//double ru;
	//double rd;
	//double dr;

	*sigma0 = TotalExtinctionPerCM( r, startPtIndex   );
	*sigma1 = TotalExtinctionPerCM( r, startPtIndex+1 );
	NXASSERT( *sigma0 >= 0.0 );
	NXASSERT( *sigma1 >= 0.0 );

	//rd = r->RadiusOfPoint( startPtIndex   );
	//ru = r->RadiusOfPoint( startPtIndex+1 );

	//dr = ru - rd;
	//if( fabs(dr) > 0.001 )
	//{
	//	*sigmak = (sigmad*ru - sigmau*rd)/dr;
	//	*sigmaf = (sigmau - sigmad)/dr;
	//}
	//else
	//{
	//	*sigmak = sigmad;
	//	*sigmaf = 0;
	//}
	return true;

	//return GetEffectiveExtinctionPerCMWithHeightShell( startpoint, startpoint.Altitude(), startpoint.Radius(), endpoint.Altitude(), endpoint.Radius(), sigmak, sigmaf );
}


bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetScatteringCoefficientCM2 ( const SKTRAN_RayStorage_Base* r, size_t index, double cosangle,  SKTRAN_PhaseMatrixScalar* phasematrix ) const
{
	bool ok = true;

	ok = false;

	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::GetEffectiveExtinctionPerCMWithHeight		2013-05-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetEffectiveExtinctionPerCMWithHeight1( const HELIODETIC_POINT& startpoint, const HELIODETIC_POINT& endpoint, double* sigma0, double* sigma1 ) const
{
	//double sigmad;
	//double sigmau;
	//double ru;
	//double rd;
	//double dr;

	*sigma0 = TotalExtinctionPerCM( startpoint );
	*sigma1 = TotalExtinctionPerCM( endpoint   );
	NXASSERT( *sigma0 >= 0.0 );
	NXASSERT( *sigma1 >= 0.0 );

	//rd = startpoint.Radius();
	//ru = endpoint.Radius();

	//dr = ru - rd;
	//if( fabs(dr) > 0.001 )
	//{
	//	*sigmak = (sigmad*ru - sigmau*rd)/dr;
	//	*sigmaf = (sigmau - sigmad)/dr;
	//}
	//else
	//{
	//	*sigmak = sigmad;
	//	*sigmaf = 0;
	//}
	return true;

	//return GetEffectiveExtinctionPerCMWithHeightShell( startpoint, startpoint.Altitude(), startpoint.Radius(), endpoint.Altitude(), endpoint.Radius(), sigmak, sigmaf );
}


bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetLinearExtinctionPerCMVector( const std::vector< HELIODETIC_POINT>& quadpoints, std::vector<double>& sigmak, std::vector<double>& sigmaf, size_t numcells ) const
{
	bool ok = true;
	double sigmad;
	double sigmau;
	double rd;
	double ru;
	double dr;

	ok = ok && sigmak.size() >= numcells;
	ok = ok && sigmaf.size() >= numcells;
	ok = ok && quadpoints.size() > numcells;
	if( ok )
	{
		rd = quadpoints[0].Radius();
		sigmad = TotalExtinctionPerCM( quadpoints[0] );
		for( size_t pointidx = 0; pointidx < numcells; pointidx++ )
		{
			sigmau = TotalExtinctionPerCM( quadpoints[pointidx+1] );
			ru = quadpoints[pointidx+1].Radius();
			dr = ru - rd;
			if( fabs(dr) > 0.001 )
			{
				sigmak[pointidx] = 100.0*(sigmad*ru - sigmau*rd)/dr;	/* Watch that the multiplication by 100.0 doesn't screw you up -- this is done differently in the non-vectored version of this function */
				sigmaf[pointidx] = 100.0*(sigmau - sigmad)/dr;
			}
			else
			{
				sigmak[pointidx] = 100.0*sigmad;
				sigmaf[pointidx] = 100.0*0;
			}
			rd		= ru;
			sigmad	= sigmau;
		}
	}
	else
	{
		nxLog::Record( NXLOG_WARNING, "SKTRAN_TableOpticalProperties_1D_Height_V3::GetLinearExtinctionPerCMVector vectors not sized correctly" );
	}
	return ok;
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::GetExtinctionPerMeterWithHeight		2008-1-31*/
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

bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetEffectiveExtinctionPerCMWithHeightShell1( const HELIODETIC_POINT& location, double h0, double r0, double h1, double r1, double* sigma0, double* sigma1 ) const
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

	//nx1dArray<double>* extinction;
	std::vector<double>* extinction;

	extinction = *m_extinction;

	if (r1 < r0)															// MAke sure r0 and r1 are in ascending order
	{
		ru = r0;
		r0 = r1;
		r1 = ru;

		ru = h0;
		h0 = h1;
		h1 = ru;
	}

	NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** SKTRAN_TableOpticalProperties_1D_Height_V3::GetExtinctionPerCMWithHeight, check out indexing on shell boundaries is correct\n"));
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
			*sigma0 = (*extinction)[lowindex];			// Convert extinction per CM to extenction per meter.
			*sigma1  = *sigma0;
		}
		else
		{
			//rd     = CoordinatesPtr()->AltitudeToRadius( m_altitudegrid->At(lowindex) );
			//ru     = CoordinatesPtr()->AltitudeToRadius( m_altitudegrid->At(hihindex) );
			//dr     = ru - rd;
			//NXASSERT(( dr > 0.0)); 
			*sigma0 = (*extinction)[lowindex];			// Convert extinction per CM to extenction per meter.
			*sigma1 = (*extinction)[hihindex];
			//*sigmak = (sigmad*ru - sigmau*rd)/dr;
			//*sigmaf = (sigmau-sigmad)/dr;

		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::GetExtinctionPerMeterWithHeight		2008-1-31*/
/** Get the extinction per meter at altitude radius, r.
**/
/*---------------------------------------------------------------------------*/
//double SKTRAN_TableOpticalProperties_1D_Height_V3::GetExtinctionPerCM( const HELIODETIC_POINT& location, const nx1dArray<double>& extinction  ) const
//{
//	static bool track1DHeightMessage = true;
//	if (track1DHeightMessage){
//		nxLog::Record(NXLOG_WARNING, "SKTRAN_TableOpticalProperties_1D_Height_V3::GetExtinctionPerCM handle nx1dArray switchover.");
//		track1DHeightMessage = false;
//	}
//	return GetExtinctionPerCM(location,((nx1dArray<double>) extinction).STLVector());
//	//return 0.0;
//}

double SKTRAN_TableOpticalProperties_1D_Height_V3::GetExtinctionPerCM( double r, const std::vector<double>&	extinction  ) const
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
	

	ok1 = m_altitudegrid->IndexOfPointBelowOrEqual( r, &lowindex );						// Find the grid point below this radius
	ok2 = m_altitudegrid->IndexOfPointEqualOrAbove( r, &hihindex );						// find the grid point above this radius
	ok  = ok1 && ok2;													// Check that everythimg is good
	if (ok)																// it is
	{																	// so
		sigmad = extinction[lowindex];								// get extinction per CM
		if ( lowindex != hihindex)										// if we have different radii
		{																// then
			sigmau = extinction[hihindex];							// Get extinction per CM
			rd     = m_altitudegrid->At(lowindex);						// get the lower radius
			ru     = m_altitudegrid->At(hihindex);						// get the higher radius
			sigma  = sigmad + (sigmau-sigmad)*(r-rd)/(ru-rd);			// Do the linear interpolation
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
			if( r > M_GROUNDTOLERANCE ){									// If the point is *just* below the ground shell,
				sigma = extinction[lowindex];								// forgive some rounding error,
				NXTRACE_ONCEONLY(firsttime,("**** 2013-08-13 ***** SKTRAN_TableOpticalProperties_1D_Height_V3::GetExtinctionPerMeter, point below ground but treated as ground\n"));
			} else{															// otherwise, 
				sigma = 200.0;												// below ground so set optical depth per meter to a very high value		
				NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** SKTRAN_TableOpticalProperties_1D_Height_V3::GetExtinctionPerMeter, point below ground\n"));
			}
		}
		else															// otherwise
		{																// the point is above this grid
			sigma = 0.0;												// and signma is zero.
		}
	}
	return sigma;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::GetExtinctionPerMeterWithHeight		2008-1-31*/
/** Getthe extinction per meter at altitude radius, r.
**/
/*---------------------------------------------------------------------------*/

double SKTRAN_TableOpticalProperties_1D_Height_V3::TotalExtinctionPerCM( const HELIODETIC_POINT& location  )const
{
	//nx1dArray<double>* extinction = *m_extinction;
	std::vector<double>* extinction = *m_extinction;
	return GetExtinctionPerCM( location.Altitude(), *extinction );
}


double SKTRAN_TableOpticalProperties_1D_Height_V3::TotalExtinctionPerCM ( const SKTRAN_RayStorage_Base* r, size_t index ) const
{
	std::vector<double>* extinction = *m_extinction;
	return GetExtinctionPerCM( r->AltitudeOfPoint(index), *extinction );
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::ScatteringExtinctionPerCM		2010-5-31*/
/** Get the scattering extinction per centimeter at altitude radius, r.
**/
/*---------------------------------------------------------------------------*/

double SKTRAN_TableOpticalProperties_1D_Height_V3::ScatteringExtinctionPerCM( const HELIODETIC_POINT& location  )const
{
	//nx1dArray<double>* scatextinction = *m_scatextinction;
	std::vector<double>* scatextinction = *m_scatextinction;
	return GetExtinctionPerCM( location.Altitude(), *scatextinction );
}


double SKTRAN_TableOpticalProperties_1D_Height_V3::ScatteringExtinctionPerCM ( const SKTRAN_RayStorage_Base* r, size_t index ) const
{
		std::vector<double>* scatextinction = *m_scatextinction;
	return GetExtinctionPerCM( r->AltitudeOfPoint(index), *scatextinction );
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::Configure		2007-11-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height_V3::ConfigureGeometry( const SKTRAN_GridDefScatterAngle_V21& scatteranglegrid, const SKTRAN_GridDefOpticalPropertiesRadii_V21& altitudegrid )
{
	bool								ok;
	size_t								numshells;
	size_t								numangles;

	ReleaseObjects();
	m_wavelen = 0.0;
	ok = true;

	m_altitudegrid		= &altitudegrid;
	m_scatteranglegrid	= &scatteranglegrid;

	ok =    (m_altitudegrid     != NULL)
		 && (m_scatteranglegrid != NULL);

	if (ok)
	{
		m_altitudegrid->AddRef();
		m_scatteranglegrid->AddRef();
		numshells = m_altitudegrid->NumAltitudes();
		numangles = m_scatteranglegrid->NumAngles();
		ok        = Allocate( numshells, numangles );

		#if defined(NXDEBUG)
			if (ok)
			{
				nxdebugFunction();
			}
		#endif
	}
	return ok;
}


bool SKTRAN_TableOpticalProperties_1D_Height_V3::ConfigureGeometry( const SKTRAN_Specifications_Base* specs )
{
	bool ok = true;

	size_t								numshells;
	size_t								numangles;

	ReleaseObjects();
	m_wavelen = 0.0;
	ok = true;

	//SKTRAN_SpecsInternal_V21* specinternal = reinterpret_cast<const SKTRAN_SpecsInternal_V21*>(specs);

	//m_altitudegrid		= &altitudegrid;
	//m_scatteranglegrid	= &scatteranglegrid;

	ok =    (m_altitudegrid     != NULL)
		 && (m_scatteranglegrid != NULL);

	if (ok)
	{
		m_altitudegrid->AddRef();
		m_scatteranglegrid->AddRef();
		numshells = m_altitudegrid->NumAltitudes();
		numangles = m_scatteranglegrid->NumAngles();
		ok        = Allocate( numshells, numangles );

		#if defined(NXDEBUG)
			if (ok)
			{
				nxdebugFunction();
			}
		#endif
	}
	return ok;

}


void SKTRAN_TableOpticalProperties_1D_Height_V3::nxdebugFunction()
{
	//(*m_extinction)->SetTo(99999.0);			// Array [numcells] of extinction .
	//(*m_singleScatt)->SetTo(99999.0);			// Array (numcells, numangles) Indexed by cell and angle.
	(*m_extinction)->assign((*m_extinction)->size(),99999.0);			// Array [numcells] of extinction .
	//(*m_singleScatt)->SetTo(99999.0);			// Array (numcells, numangles) Indexed by cell and angle.
	//m_singleScatt.assign( m_singleScatt.size(), 99999.0 );
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::GetBoundingPoints		2009-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetBoundingSpatialPoints( const HELIODETIC_POINT& location, SKTRANSO_JIndex* interppoints, size_t* numpoints ) const
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
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableOpticalProperties_1D_Height_V3::GetBoundingGridPoints, Error retrieving bounding grid points for height = %f", (double)height);
	}
	return ok;
}

*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::GetBoundingGridPoints		2009-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetBoundingScatteringPoints( double cosscatteringangle, SKTRANSO_JIndex* interppoints, size_t* numpoints ) const
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
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableOpticalProperties_1D_Height_V3::GetBoundingScatteringPoints, Error retrieving bounding grid points for cos(scattering angle) = %f", (double)cosscatteringangle);
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::GetScatteringCoefficient		2008-1-10*/
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


bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetScatteringCoefficientCM2( const HELIODETIC_POINT& location, double cosAngle, SKTRAN_PhaseMatrixScalar* coefficient ) const
{
    bool ok = true;

    SKTRAN_GridIndex gridindices[4];
    double gridweights[4];
    size_t numNonZero;
    ok = ok && GetUniquePointWeights( location, cosAngle, gridindices, gridweights, numNonZero );
    ok = ok && m_scatprops->GetPhaseFunctionCM2( gridindices, gridweights, numNonZero, *coefficient );
    if(!ok){
        NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** SKTRAN_TableOpticalProperties_1D_Height_V3::GetScatteringCoefficient, Bounding indices failed\n"));
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableOpticalProperties_1D_Height_V3::GetScatteringCoefficient, failed finding index for height (%10g) meters, cosangle (%g).  Make sure the optical properties grid spans above the ray tracing and diffuse grids", (double) location.Altitude(), (double)cosAngle );
		*coefficient = SKTRAN_DBL_TO_STOKES_SCALAR(0);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::GetScatteringCoefficientCM2		2009-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetScatteringCoefficientCM2( const SKTRANSO_JIndex* boundingpoints, size_t numpoints, double cosangle, SKTRAN_PhaseMatrixScalar* coefficient ) const
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

	nx2dArray<double>*	singleScatt = m_singleScatt[0];

	ok      = m_scatteranglegrid->FindBoundingIndices( cosangle, SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &lowerscatidx, &lowscatw, &upperscatidx, &upperscatw );
	cellidx = boundingpoints->HeightIndex();
	cellw   = boundingpoints->VertexWeight();
	c       = (singleScatt->At(cellidx,lowerscatidx)*lowscatw + singleScatt->At(cellidx,upperscatidx)*upperscatw)*cellw;
	if (numpoints > 1)
	{
		for (pidx = 1; pidx < numpoints; pidx++ )
		{
			cellidx = boundingpoints[pidx].HeightIndex();
			cellw   = boundingpoints[pidx].VertexWeight();
			c      += (singleScatt->At(cellidx,lowerscatidx)*lowscatw + singleScatt->At(cellidx,upperscatidx)*upperscatw)*cellw;
		}
	}
	*coefficient = SKTRAN_DBL_TO_STOKES_SCALAR(c);
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::ConfigureOptical		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height_V3::ConfigureOptical( double wavelen, SKTRAN_AtmosphericOpticalState_V21& opticalstate )
{
	SKTRAN_GridIndex		shellidx;
	SKTRAN_GridIndex		scatidx;
	GEODETIC_INSTANT		point;
	bool					ok;
	bool					ok1;
	bool					ok2 = true;
//	double					scatter;
//	double					cosangle;
	skBRDF*					albedo;

	NXASSERT(( wavelen > 0.0 ));

	ok = opticalstate.GetAlbedoObject( &albedo );
	if (ok)
	{
		albedo->AddRef();
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

				(*m_extinction)->at(shellidx)          = opticalstate.ExtinctionPercm();	
				(*m_scatextinction)->at(shellidx)	   = opticalstate.ScatteringPercm();					
// Stuff below should be removed because this check is performed when the table are filled. It shouldn't be on the accessor to check the validity of table contents
//				ok1  = ok1 && (m_extinction.at(shellidx) >= 0.0);
//				ok1  = ok1 && (m_scatextinction.at(shellidx) >= 0.0);
//				NXASSERT((m_scatextinction.at(shellidx)>= 0.0));
//				NXASSERT((m_extinction.at(shellidx)>= 0.0));

				for (scatidx = 0; scatidx < m_scatteranglegrid->NumAngles(); scatidx++)
				{
					double cosangle;
					skRTPhaseMatrix pmatrix;

					cosangle	= m_scatteranglegrid->At( scatidx );
					ok          = opticalstate.VectorPhaseMatrix( cosangle, &pmatrix  );
                    m_scatprops->StorePolarizationPropsCM2( TableSubToInd(shellidx, scatidx), pmatrix, opticalstate );
					NXASSERT((m_scatprops->PhaseMatrixAccess(TableSubToInd(shellidx, scatidx)) >= 0.0));
				}
			}
			ok = ok && ok1;
		}
	}


	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableOpticalProperties_1D_Height_V3::ConfigureOptical, Error configuring the optical state");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableOpticalProperties_1D_Height_V3::GetAlbedo		2008-3-20*/
/** Returns the albedo at the specified point**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetBRDF( const HELIODETIC_POINT& userpoint, double mu_in, double mu_out, double cosdphi, double* brdf ) const

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
		point = CoordinatesPtr()->PointToGeodetic( userpoint, m_mjd );
		ok = m_albedo->BRDF( m_wavelen, point, mu_in, mu_out, cosdphi, brdf);
	}
	return ok;
}

/*-----------------------------------------------------------------------------
*					SKTRAN_TableOpticalProperties_1D_Height_V3::GetBRDFGeodetic		2017-05-01*/
/** Returns the albedo at the specified point**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetBRDFGeodetic(const GEODETIC_INSTANT& userpoint, double mu_in, double mu_out, double cosdphi, double* brdf) const

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

bool SKTRAN_TableOpticalProperties_1D_Height_V3::IsOptionTrue( SKTRAN_TableOpticalProperties_Base::OPTIONSENUM options) const
{
	return (options == OPTION_IS_HORIZONTALLY_UNIFORM);
}


bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetUniquePointWeights( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_GridIndex gridindices[4], double gridweights[4], size_t& numNonZero   ) const 
{
	bool ok = true;
	
	SKTRAN_GridIndex loAltIndex, hiAltIndex, loAngIndex, hiAngIndex;
	double loAltWeight, hiAltWeight, loAngWeight, hiAngWeight;
	const double minval = 0.0;

	ok = ok && m_altitudegrid->FindBoundingIndices     ( point.Altitude(), SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &loAltIndex, &loAltWeight, &hiAltIndex, &hiAltWeight );
	ok = ok && m_scatteranglegrid->FindBoundingIndices ( cosAngle,         SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &loAngIndex, &loAngWeight, &hiAngIndex, &hiAngWeight );
	
	numNonZero = 0;
	if( minval < loAltWeight*loAngWeight ){ gridweights[numNonZero]=loAltWeight*loAngWeight; gridindices[numNonZero]=TableSubToInd( loAltIndex, loAngIndex ); ++numNonZero;}
	if( minval < loAltWeight*hiAngWeight ){ gridweights[numNonZero]=loAltWeight*hiAngWeight; gridindices[numNonZero]=TableSubToInd( loAltIndex, hiAngIndex ); ++numNonZero;}
	if( minval < hiAltWeight*loAngWeight ){ gridweights[numNonZero]=hiAltWeight*loAngWeight; gridindices[numNonZero]=TableSubToInd( hiAltIndex, loAngIndex ); ++numNonZero;}
	if( minval < hiAltWeight*hiAngWeight ){ gridweights[numNonZero]=hiAltWeight*hiAngWeight; gridindices[numNonZero]=TableSubToInd( hiAltIndex, hiAngIndex ); ++numNonZero;}
    
	return ok;
}


SKTRAN_GridIndex SKTRAN_TableOpticalProperties_1D_Height_V3::TableSubToInd ( SKTRAN_GridIndex altidx, SKTRAN_GridIndex angidx ) const
{
	return altidx*m_scatteranglegrid->NumAngles() + angidx;
}


bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetScatteringMatrixCM2 ( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_ScatMat_MIMSNC&  pmatrix   ) const 
{
	bool ok = true;

	SKTRAN_GridIndex gridindices[4];
	double           gridweights[4];
	size_t           numNonZero;

	ok = ok && GetUniquePointWeights ( point, cosAngle, gridindices, gridweights, numNonZero );
	ok = ok && m_scatprops->GetScatteringMatrixCM2 ( gridindices, gridweights, numNonZero, pmatrix);

	return ok;
}


bool SKTRAN_TableOpticalProperties_1D_Height_V3::GetResultOfUnpolarizedScatterCM2 ( const HELIODETIC_POINT& point, double cosAngle, SKTRAN_Stokes_NC& stokesvec ) const
{
	bool ok = true;

	SKTRAN_GridIndex gridindices[4];
	double           gridweights[4];
	size_t           numNonZero;

	ok = ok && GetUniquePointWeights( point, cosAngle, gridindices, gridweights, numNonZero );
	ok && m_scatprops->GetResultOfUnpolarizedScatterCM2( gridindices, gridweights, numNonZero, stokesvec);

	return ok;
}


bool SKTRAN_TableOpticalProperties_1D_Height_V3::CreateInterpolationForPoint( const HELIODETIC_POINT& point, SKTRAN_TableOpticalProperties_PointCache& interpolator ) const
{
    bool ok = true;

    SKTRAN_GridIndex loAltIndex, hiAltIndex, loAngIndex, hiAngIndex;
    double loAltWeight, hiAltWeight, loAngWeight, hiAngWeight;

    ok = ok && m_altitudegrid->FindBoundingIndices     ( point.Altitude(), SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &loAltIndex, &loAltWeight, &hiAltIndex, &hiAltWeight );
    
    SKTRAN_GridIndex gridindices[4];
    double gridweights[4];
    size_t numNonZero;

    const size_t numAngles = m_scatteranglegrid->NumAngles();
    SKTRAN_ScatMat_MIMSNC pmatrix;
    std::vector<SKTRAN_ScatMat_MIMSNC> pmatrices;
    pmatrices.reserve( numAngles );
    const double minval = 0.0;

    for (size_t aidx = 0; aidx < numAngles; ++aidx) {
        numNonZero = 0;
        m_scatteranglegrid->FindBoundingIndices ( m_scatteranglegrid->At(aidx), SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &loAngIndex, &loAngWeight, &hiAngIndex, &hiAngWeight );
        if( minval < loAltWeight*loAngWeight ){ gridweights[numNonZero]=loAltWeight*loAngWeight; gridindices[numNonZero]=TableSubToInd( loAltIndex, loAngIndex ); ++numNonZero;}
        if( minval < loAltWeight*hiAngWeight ){ gridweights[numNonZero]=loAltWeight*hiAngWeight; gridindices[numNonZero]=TableSubToInd( loAltIndex, hiAngIndex ); ++numNonZero;}
        if( minval < hiAltWeight*loAngWeight ){ gridweights[numNonZero]=hiAltWeight*loAngWeight; gridindices[numNonZero]=TableSubToInd( hiAltIndex, loAngIndex ); ++numNonZero;}
        if( minval < hiAltWeight*hiAngWeight ){ gridweights[numNonZero]=hiAltWeight*hiAngWeight; gridindices[numNonZero]=TableSubToInd( hiAltIndex, hiAngIndex ); ++numNonZero;}

        ok = ok && m_scatprops->GetScatteringMatrixCM2( gridindices, gridweights, numNonZero, pmatrix);
        pmatrices.push_back( pmatrix*100.0 );
    }
    
    double scatteringExtinction = ScatteringExtinctionPerCM ( point );
    double totalExtinction      = TotalExtinctionPerCM      ( point );
    interpolator.InjectTable( pmatrices, m_scatteranglegrid, scatteringExtinction/totalExtinction, point );

    return ok;
}


