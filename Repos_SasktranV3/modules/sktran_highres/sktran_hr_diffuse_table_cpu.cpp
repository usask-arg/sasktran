#include "include/sktran_hr_internals.h"
#include <numeric>
#include <functional>
#include <iostream>


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::SKTRAN_HR_Diffuse_Table_CPU		2013-06-28*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Diffuse_Table_CPU::SKTRAN_HR_Diffuse_Table_CPU()
{
	m_optintegrator = nullptr;
	m_srcintegrator = nullptr;
	m_rayfactory     = nullptr;
	m_sources       = nullptr;
	m_opticaltable  = nullptr;
	m_radStorage    = nullptr;
	m_Avals         = nullptr;
	m_diffuselocationsphere = nullptr;
	m_numprofileinterp = 2;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::~SKTRAN_HR_Diffuse_Table_CPU		2013-06-28*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Diffuse_Table_CPU::~SKTRAN_HR_Diffuse_Table_CPU()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::ConfigureStorage		 2016- 12- 5*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::ConfigureStorage( std::unique_ptr< RadStore_Base >& radStorage, std::unique_ptr< Avals_Base >&   avals )
{
	m_radStorage = std::move( radStorage );
    m_Avals      = std::move( avals );
    m_Avals->SetOpticalTable( m_opticaltable );
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::PreSetup		 2016- 12- 5*/
/** Setup the scheme that interpolates in height by first dividing by extinction.
 *	This is much more accurate than linear interpolation  when the radiance and
 *	extinction scale exponentially.
 *	The albedo is used for ground points and we are quite happy to use the deprecated
 *	legacy code as it takes the albedo from the BRDF for straight down and straight up rays.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::PreSetup()
{
	for( size_t idx = 0; idx < m_kscatatpoints.size(); idx++ )
	{
        if(!m_diffusepoints[idx].IsGroundPoint())																	// IF this is not a ground point
		{																											// then 
            m_kscatatpoints[idx] = m_opticaltable->ScatteringExtinctionPerCM( m_diffusepoints[idx].Location() );	// Setup the interpolation table using the extinction
        }																											// and that is that
		else																										// otherwise this is a ground
		{																											// so
            double albedo;																							// setup the interpolation
            m_opticaltable->Get_AlbedoForDeprecatedLegacyCode( m_diffusepoints[idx].Location(), &albedo );			// using the albedo, the deprecated code is good fopr this
            m_kscatatpoints[idx] = albedo;																			
        }
	}
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::Initialize		 2016- 12- 5*/
/**  This is called prior to any real calculation, any sort of table setup/initialization that needs to be done should
 *   be done here
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::Initialize( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>&        coords,
											  const SKTRAN_OpticalPropertiesIntegrator_Base& optintegrator,
											  const SKTRAN_SourceTermIntegrator_Base&        srcintegrator,
											  std::shared_ptr<const SKTRAN_RayFactory_Base>   rayfactory,
											  const std::vector<SKTRAN_Source_Term*>&        sources,
											  const SKTRAN_TableOpticalProperties_Base&      opticaltable )
{
	bool ok = true;

	m_coords =        coords;
	m_optintegrator = &optintegrator;
	m_optintegrator-> AddRef();
	m_srcintegrator = &srcintegrator;
	m_srcintegrator-> AddRef();
	m_rayfactory    = rayfactory;
	m_sources =       &sources;
	m_opticaltable =  &opticaltable;
	m_opticaltable->  AddRef();

	if(nullptr!=m_radStorage){
		ok = ok && m_radStorage->Initialize( coords, optintegrator, srcintegrator, rayfactory, sources, opticaltable );
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Diffuse_Table_CPU::Initialize, radStorage should exist for initialization.");
	}

	if(nullptr!=m_Avals){
		m_Avals->SetOpticalTable( m_opticaltable );
	} else{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Diffuse_Table_CPU::Initialize, Avals should exist for initialization.");
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::ReleaseResources		2013-06-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::ReleaseResources()
{
	m_diffusepoints.clear();

	if( nullptr != m_optintegrator )  m_optintegrator-> Release();
	if( nullptr != m_srcintegrator )  m_srcintegrator-> Release();
	if( nullptr != m_opticaltable  )  m_opticaltable->  Release();

	m_optintegrator = nullptr;
	m_srcintegrator = nullptr;
	m_sources       = nullptr;
	m_opticaltable  = nullptr;
	m_radStorage    = nullptr;
	m_Avals         = nullptr;

	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::DumpOutGoingRadiances		2013-06-28*/
/** Dumps the current order outgoing radiances **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_HR_Diffuse_Table_CPU::DumpOutGoingRadiances( size_t order, double wlen )
{
	return m_radStorage->DumpOutgoingRadiances( this, order, wlen );
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::DumpIncomingRadiances		2013-06-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::DumpIncomingRadiances( size_t order, double wlen )
{
	return m_radStorage->DumpIncomingRadiances ( this, order, wlen );
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::DistanceWeight		2013-07-11*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_HR_Diffuse_Table_CPU::DistanceWeight( double distance ) const
{
	if( distance > 1E-5 )
	{
		return 1.0/( distance );
	}
	else
	{
		return 1E10;
	}
}

SKTRAN_Stokes_NC SKTRAN_HR_Diffuse_Table_CPU::IntegrateVectorIncoming(const SKTRAN_SourceTermQueryObject_Base & qobj,
																	  std::function<SKTRAN_ScatMat_MIMSNC(
																			  const nxVector &,
																			  const nxVector &)> func) const {
	std::array<size_t, 8> diffindex;
	std::array<SKTRAN_HR_WEIGHT_TYPE, 8> diffweight;
	size_t numpoint;

	const HELIODETIC_POINT& pt = qobj.GetPoint();
	const HELIODETIC_UNITVECTOR& heliolook = qobj.GetLookAway();



	numpoint = diffindex.max_size();
	ChooseDiffusePoints(pt, &diffindex[0], &diffweight[0], &numpoint);

	SKTRAN_Stokes_NC result;
	result.SetTo(0.0);
    std::unique_ptr< EtaCalculator_Base > etaCalculator;
    m_Avals->CreateEtaCalculator( etaCalculator );

	for( size_t pointidx = 0; pointidx < numpoint; pointidx++ )
	{
        // Away from observer
        nxVector look(heliolook.X(), heliolook.Y(), heliolook.Z());

		SKTRAN_Stokes_NC pointresult;
		pointresult.SetTo(0.0);
		const SKTRAN_HR_Diffuse_Point& diffusepoint = m_diffusepoints[diffindex[pointidx]];


        etaCalculator->SetPoint( diffusepoint );

        HELIODETIC_UNITVECTOR look_towards;

        RotateRayToDiffuse(qobj.GetPoint(), qobj.GetLookAway(), look);

        look_towards.SetCoords(-look.X(), -look.Y(), -look.Z());
        etaCalculator->SetOutgoingDirection(look_towards);

        for( size_t incomingidx = 0; incomingidx < diffusepoint.NumIncomingRays(); incomingidx++ )
		{
			SKTRAN_Stokes_NC incomingradiance = m_radStorage->IncomingVectorRadiance(diffusepoint.IncomingRadianceIdx(incomingidx));
			double incomingintegrationweight = diffusepoint.InCubatureWeight(incomingidx);
			//const HELIODETIC_UNITVECTOR& incominghelio = diffusepoint.IncomingRayGlobalCoords(incomingidx);
            //nxVector incomingunitray(incominghelio.X(), incominghelio.Y(), incominghelio.Z());

            nxVector incomingunitray = diffusepoint.IncomingUnitRayLocalCoords(incomingidx);
			// TODO: Check rotation of incoming unit ray
			incomingradiance *= incomingintegrationweight;
            etaCalculator->CalculateEtas( diffusepoint, incomingidx ); // Do geometry for this scatter
            etaCalculator->RefToScatt( incomingradiance );
			func(look, incomingunitray).LApplyTo(&incomingradiance);
            etaCalculator->ScattToRef( incomingradiance );
            pointresult += incomingradiance;
		}
		result += pointresult * diffweight[pointidx];
	}

	return result;
}

double SKTRAN_HR_Diffuse_Table_CPU::IntegrateScalarIncoming( const SKTRAN_SourceTermQueryObject_Base& qobj, std::function<double(const nxVector&, const nxVector&)> func ) const
{
	std::array<size_t, 8> diffindex;
	std::array<SKTRAN_HR_WEIGHT_TYPE, 8> diffweight;
	size_t numpoint;
	
	const HELIODETIC_POINT& pt = qobj.GetPoint();
	const HELIODETIC_UNITVECTOR& heliolook = qobj.GetLookAway();

	// Away from observer
	nxVector look(heliolook.X(), heliolook.Y(), heliolook.Z());

	numpoint = diffindex.max_size();
	ChooseDiffusePoints(pt, &diffindex[0], &diffweight[0], &numpoint);

	double result = 0;
	for( size_t pointidx = 0; pointidx < numpoint; pointidx++ )
	{
		double pointresult = 0;
		const SKTRAN_HR_Diffuse_Point& diffusepoint = m_diffusepoints[diffindex[pointidx]];

		for( size_t incomingidx = 0; incomingidx < diffusepoint.NumIncomingRays(); incomingidx++ )
		{
			double incomingradiance = m_radStorage->IncomingScalarRadiance(diffusepoint.IncomingRadianceIdx(incomingidx));
			double incomingintegrationweight = diffusepoint.InCubatureWeight(incomingidx);
			const HELIODETIC_UNITVECTOR& incominghelio = diffusepoint.IncomingRayGlobalCoords(incomingidx);
			nxVector incomingunitray(incominghelio.X(), incominghelio.Y(), incominghelio.Z());
			// TODO: Check rotation of incoming unit ray
			pointresult += incomingradiance * incomingintegrationweight * func(look, incomingunitray);
		}
		result += pointresult * diffweight[pointidx];
	}

	return result;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::AllocateDiffusePoints		2013-06-28*/
/** Allocates the required memory for diffuse points **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::AllocateDiffusePoints( size_t numdiffuse, size_t numground )
{
	bool ok = true;

	m_diffusepoints.resize( numdiffuse + numground );
	m_groundstartidx = numdiffuse;

	return ok;
}


bool PairLessThan( const std::pair<double,size_t>& l, const std::pair<double,size_t>& r )
{
	return l.first < r.first;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::AltWeightsForProfile		2013-07-11*/
/** For a given profile calculates which altitudes contribute to the source 
 *  and their weights
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::AltWeightsForProfile( double alt, size_t profileidx, SKTRAN_HR_WEIGHT_TYPE* altweights, size_t* altindex, size_t* numindex ) const
{
	bool								ok = true;
	std::vector<double>::const_iterator	low,up;
	double								dh;

	ok = (*numindex >= 2);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_HR_Diffuse_Table_CPU::AltWeightsForProfile, You must provide a buffer of at least 2 elements for diffuse height interpolation");
		*numindex = 0;
	}
	else
	{
		up = std::upper_bound( m_diffusealts[profileidx].begin(), m_diffusealts[profileidx].end(), alt );
		if( up != m_diffusealts[profileidx].end() )							// If we found our point below the upper most altitude
		{																	// then 
			NXASSERT( (alt <= *up ));
			if( up == m_diffusealts[profileidx].begin() )					// if the upper bound is before the first altitude
			{																// then the 
				*numindex = 1;												// set the number of points to 1
				altweights[0] = 1;											// set the weight to 1
				altindex[0] = 0;
			}																// otherwise
			else															// then lower bound
			{																// is the 
				low = up - 1;												// previous point
				*numindex = 2;
				dh    = std::max( *up - *low, 1.0E-08);
				NXASSERT( (dh > 0));
				altweights[0] = SKTRAN_HR_DBL_TO_WEIGHT((alt - *low) /dh);
				altweights[1] = SKTRAN_HR_DBL_TO_WEIGHT((*up - alt) / dh);
				altindex[0]   = up - m_diffusealts[profileidx].begin();
				altindex[1]   = low - m_diffusealts[profileidx].begin();
			}
		} 
		else if (m_diffusealts[profileidx].size() > 0)
		{
			altweights[0] = 1.0;
			altweights[1] = 0.0;
			altindex  [0] = m_diffusealts[profileidx].size()-1;
			altindex  [1] = 0;
			*numindex     = 1;
		} 
		else
		{
			altweights[0] = 0;
			altweights[1] = 0;
			altindex  [0] = 0;
			altindex  [1] = 0;
			*numindex = 0;
			//up = m_diffusealts[profileidx].end() - 1;
		};
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table::ChooseDiffusePoints		2013-07-11*/
/** Chooses which diffuse points contribute to the diffuse source at a specific point.
 *	The method finds the "N" closest diffuse profiles by calculating the distance
 *	of their unit vector from the point's unit vector.  For each diffuse profile
 *	calculate the set of weights and diffuse points indices that contribute.
 *
 *	\param pt
 *		The location of the point that will be interpolated
 *
 *	\param diffuseindex
 *		Pointer to an array that will store the indices of the diffuse points used in the interpolation.
 *		The user will guarantee the size of the array is at least the value indicated by maxnumpoints at the start of the code .
 *
 *	\param diffuseweights
 *		Pointer to an array of weights used for the interpolation. The user will guarantee this is at least maxnumpoints in size
 *		The user will guarantee the size of the array is at least the value indicated by maxnumpoints at the start of the code .
 *
 *	\param maxnumpoints
 *		On input this indicates the maximum size of the diffuseindex and diffuseweight arrays. On out it indicates the actual size of
 *		the arrays
 *	
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::ChooseDiffusePoints( const HELIODETIC_POINT& pt, size_t* diffuseindex, SKTRAN_HR_WEIGHT_TYPE* diffuseweights, size_t* maxnumpoints ) const
{
	bool									ok = true;
	bool									ok1;
	size_t									numprofile;
	double									alt;
	std::vector<double>::const_iterator		low,up;
	std::vector<std::pair<double,size_t> >	distance;
	SKTRAN_HR_WEIGHT_TYPE					altweights[4];
	size_t									altindex[4];
	SKTRAN_HR_WEIGHT_TYPE					profileweights[4];				// assume we dont have more than 4
	size_t									profileindex[4];
	SKTRAN_HR_WEIGHT_TYPE					totalprofileweight;
	double									x,y,z;
	size_t									numalt;
	HELIODETIC_UNITVECTOR					ptu;
	HELIODETIC_UNITVECTOR					profu;
	size_t									numpoints;

	numpoints = 0;
	if( m_diffuselocationsphere )
	{
		const HELIODETIC_UNITVECTOR unith = pt.Vector().UnitVector();
		const nxVector unit (unith.X(), unith.Y(), unith.Z());
		double weights[3];
		ok = ok && m_diffuselocationsphere->Triangulate( unit, profileindex, weights, 3 );
		numprofile = m_diffuselocationsphere->MaxNumInterpIndex();
		profileweights[0] = (SKTRAN_HR_WEIGHT_TYPE)weights[0];
		profileweights[1] = (SKTRAN_HR_WEIGHT_TYPE)weights[1];
		profileweights[2] = (SKTRAN_HR_WEIGHT_TYPE)weights[2];

		alt = pt.Altitude();
	}
	else
	{
		distance.resize( NumDiffuseProfiles() );											// Allocate space to find out which diffuse profile is closest to the point
		
		alt = pt.Altitude();
		ptu = pt.UnitVector();																// Get the unit vector toward this point
		for( size_t idx = 0; idx < NumDiffuseProfiles(); idx++ )															
		{
			profu = m_diffusepoints[m_profilestartidx[idx]].Location().UnitVector();		// Get the unit vector of each diffuse profile
			x = ptu.X() - profu.X();														// get the distance from our point
			y = ptu.Y() - profu.Y();														// To this diffuse profile
			z = ptu.Z() - profu.Z();
			distance[idx].first  = sqrt( x*x + y*y + z*z );									// and save the result
			distance[idx].second = idx;														// and the index of the profile
		}

		std::sort( distance.begin(), distance.end(), PairLessThan );						// Sort the distances
		numprofile = std::min( m_numprofileinterp, distance.size() );							// Restrict the number of profiles used to interpolate  to each requested interpolation number or the actial number of profiles

		totalprofileweight = 0;																		// reset  the total weight summation
		for( size_t idx = 0; idx < numprofile; idx++ )												// for each profile used for interpolation 
		{																							// get the
			profileindex[idx]   = distance[idx].second;												// index of the diffuse profile
			profileweights[idx] = SKTRAN_HR_DBL_TO_WEIGHT(DistanceWeight( distance[idx].first ));	// Get the non-normalized weight to be applied
			totalprofileweight += profileweights[idx];												// and sum the weights
		}
		
		NXASSERT(totalprofileweight != 0.0);

		for( size_t idx = 0; idx < numprofile; idx++ )											// renormalize the weights
		{																			
			profileweights[idx] /= totalprofileweight;
		}
	}

	for( size_t profileidx = 0; profileidx < numprofile; profileidx++ )							// for each diffuse profile that contributes to the interpolation
	{																							// get
		numalt = N_ELEMENTS(altweights);
		ok1 = ok && AltWeightsForProfile( alt, profileidx, altweights, altindex, &numalt );		// The indices and weights of points that contribute to the interpolation
		if (ok1)
		{
			ok1 = numpoints < *maxnumpoints;
			for( size_t altidx = 0; (altidx < numalt) && ok1; altidx++ )
			{
				diffuseindex[numpoints]   = m_profilestartidx[profileindex[profileidx]] + altindex[altidx];
				diffuseweights[numpoints] = profileweights[profileidx] * altweights[altidx];
				++numpoints;
				ok1 = numpoints < *maxnumpoints; // TODO if we need 6 points this will trigger an error
			}
			if (!ok1)
			{
				nxLog::Record(NXLOG_WARNING,"SKTRAN_HR_Diffuse_Table_CPU::ChooseDiffusePoints, Error generatinging interpolation tables, insuffucient storage passed in by user. Thats not good");
			}
		}
		ok = ok && ok1;
	}
	*maxnumpoints = numpoints;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table::ChooseGroundPoints		2013-06-27*/
/** For a given position, calculates the number of ground diffuse points
 *  that contribute to the ground radiance signal and their weights
 * Chooses which diffuse ground points contribute to the source at a ground point
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::ChooseGroundPoints( const HELIODETIC_POINT& pt, size_t* diffuseindex, SKTRAN_HR_WEIGHT_TYPE* diffuseweights, size_t& numpoints ) const
{
	bool ok = true;
	std::vector<std::pair<double,size_t> >	distance;
	size_t									profileindex[4];
	SKTRAN_HR_WEIGHT_TYPE					profileweights[4];
	double									x,y,z;
	SKTRAN_HR_WEIGHT_TYPE					totalprofileweight;
	numpoints = m_numprofileinterp;

	distance.resize( m_profilestartidx.size()-1 );
	for( size_t idx = 0; idx < m_profilestartidx.size()-1; idx++ )
	{
		x = pt.Vector().UnitVector().X() - m_diffusepoints[m_profilestartidx[idx]].Location().Vector().UnitVector().X();
		y = pt.Vector().UnitVector().Y() - m_diffusepoints[m_profilestartidx[idx]].Location().Vector().UnitVector().Y();
		z = pt.Vector().UnitVector().Z() - m_diffusepoints[m_profilestartidx[idx]].Location().Vector().UnitVector().Z();
		distance[idx].first = sqrt( x*x + y*y + z*z );
		distance[idx].second = idx;
	}
	std::sort( distance.begin(), distance.end(), PairLessThan );
	numpoints = std::min(numpoints, distance.size());
	totalprofileweight = 0;
	for( size_t idx = 0; idx < numpoints; idx++ )
	{
		profileindex[idx]   = distance[idx].second;
		profileweights[idx] = SKTRAN_HR_DBL_TO_WEIGHT(DistanceWeight( distance[idx].first ));
		totalprofileweight += profileweights[idx];
	}
	for( size_t idx = 0; idx < numpoints; idx++ )
	{
		profileweights[idx] /= totalprofileweight;
	}

	for( size_t idx = 0; idx < numpoints; idx++ )
	{
		diffuseweights[idx] = profileweights[idx];
		diffuseindex[idx] = m_groundstartidx + profileindex[idx];
	}

	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::CalcFirstOrderIncomingPoint		2013-07-11*/
/** Computes the first order incoming radiance for a specific point, a thread specific ray may be
 *  specified
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::CalcFirstOrderIncomingPoint( size_t pointidx, SKTRAN_RayOptical_Base* ray )
{
	bool ok = true;
	size_t					numinrays;
	const SKTRAN_HR_Diffuse_Point& diffusepoint = DiffusePointAt( pointidx );
	numinrays= diffusepoint.NumIncomingRays();

	// calculate and store the radiance from each incoming ray
	for( size_t rayidx = 0; rayidx < numinrays; rayidx++ )
	{
		ok = ok && CalcFirstOrderIncomingRay(pointidx, rayidx, ray);
	}

	if( !ok )
	{
		nxLog::Record( NXLOG_WARNING, "Problem in SKTRAN_HR_Diffuse_Table_CPU::CalcFirstOrderIncomingPoint" );
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::CalcFirstOrderIncomingPoint		2020-06-24*/
/** Computes the first order incoming radiance for a specific point and ray,
 *  put over rays so we can parallelize over rays instead of points.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::CalcFirstOrderIncomingRay( size_t pointidx, size_t rayidx, SKTRAN_RayOptical_Base* ray )
{
	bool ok = true;
	HELIODETIC_VECTOR		observer;
	HELIODETIC_UNITVECTOR	look;
	const SKTRAN_HR_Diffuse_Point& diffusepoint = DiffusePointAt( pointidx );

	observer = diffusepoint.Location().Vector();

	look = diffusepoint.IncomingRayGlobalCoords( rayidx );
	ok = ok && ray->MoveObserver( observer, look);
	ok = ok && ray->TraceRay_NewMethod();
	ok = ok && m_optintegrator->CalculateRayScalarTransmission_withMinContainer( ray, nullptr, false, false );
	m_radStorage->StoreFirstOrderIncoming( m_srcintegrator, m_sources, ray, diffusepoint.IncomingRadianceIdx(rayidx) );
	ok = ok && CreateDiffuseIndexesForRay( ray, pointidx, rayidx );

	if( !ok )
	{
		nxLog::Record( NXLOG_WARNING, "Problem in SKTRAN_HR_Diffuse_Table_CPU::CalcFirstOrderIncomingRay" );
	}

	return ok;
}




/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::CalcScatteringMatrixPoint		2013-07-11*/
/** Calculates the scattering matrix for a given point
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::CalcScatteringMatrixPoint	( size_t pointidx )
{
	bool ok = true;

	ok = ok && m_Avals->CalcScatteringMatrixPoint ( m_diffusepoints[pointidx] );

	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::ScatterPoint		2013-07-11*/
/** Scatters the current incoming signal to get the next order outgoing signal
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::ScatterPoint( size_t pointidx )
{
	bool ok = true;
	//const size_t							numops		= 0;
	const SKTRAN_HR_Diffuse_Point&	point	= DiffusePointAt( pointidx );
	
	ok = ok && m_radStorage->ScatterPoint( point, *m_Avals );

	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::ComputeNextOrderPoint		2013-07-11*/
/** Computes the next order incoming signal **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::ComputeNextOrderPoint( size_t pointidx )
{
	bool ok = true;
	const SKTRAN_HR_Diffuse_Point&	point		= DiffusePointAt( pointidx );

	ok = ok && m_radStorage->ComputeNextOrderPoint( point );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::ComputeMultipliersAndAdjustScatterArray		2013-06-25*/
/** Truitts photon conservation code modified to work inside the sparse matrix 
 *  format used here
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::ComputeMultipliersAndAdjustScatterArray( size_t pointidx )
{
	return m_Avals->ComputeMultipliersAndAdjustScatterArray( m_diffusepoints[pointidx] );
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::DiffuseSource		2013-06-27*/
/** Returns the source term from the interpolation of the diffuses points to 
 *	the location and in the opposite direction to the look vector stored in qobj. The answer is 
 *	returned as either a scalar or polarized radiance.
 *
 *  Takes a quadrature point and a look direction and calculates the diffuse 
 *  source(radiance/m) in the reverse look direction. The code provides a much 
 *	more accurate interpolation of source terms by dividing out the
 *	scattering cross-section before performing the linear interpolation. This works very well
 *	for exponential height variation and is acceptable for low order polynomials, see equation below for demo
 *	theory.
 *
 *	 \f{eqnarray*}{
 *      y &=& I_{0}e^{-kh}\\
 *		z &=& \frac{y}{e^{-kh}} = I_{0}\\
 *		z_0 &=& \frac{y_0}{e^{-kh_0}} \approx I_{0}\\
 *		z_1 &=& \frac{y_1}{e^{-kh_1}} \approx I_{0}\\	
 *		z &\approx& w_0z_0+w_1z_1\\
 *		y &\approx& e^{-kh}z 			
 * }
 *
 *	This technique only works properly if the underlying optical table has cross-sections
 *	stored at substantially higher spatial resolution than the source function tables but it
 *	degrades gracefully when that condition  is not met and works  much better when it is.
 **/
/*---------------------------------------------------------------------------*/
template< typename Radtype >
bool SKTRAN_HR_Diffuse_Table_CPU::DiffuseSource_impl( const SKTRAN_SourceTermQueryObject_Base& qobj, Radtype& source ) const
{
    bool								ok;
    size_t								diffuseindex     [6];
    SKTRAN_HR_WEIGHT_TYPE			    diffuseweight    [6];
    size_t								numpoints = N_ELEMENTS(diffuseindex);
	double								kscaatpoint;
	double								kscat;
    Radtype								source_temp;
    nxVector							inlook;
    //	HELIODETIC_BASIS					basis(pt,look);

    ok = ChooseDiffusePoints(qobj.GetPoint(), diffuseindex, diffuseweight, &numpoints );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_HR_Diffuse_Table_CPU::DiffuseSource_impl, Error calling ChooseDiffusePoints, This is a serious problem. Dont trust answers.");
	}
	else
	{
		kscat = m_opticaltable->ScatteringExtinctionPerCM( qobj.GetPoint() );
		skRTStokesVector::SetToZero(source);

		if ( kscat > 0)
		{
			for( size_t pointidx = 0; pointidx < numpoints; pointidx++ )
			{
				//printf("%zu: %e, %e, %e\n", diffuseindex[pointidx], m_diffusepoints[diffuseindex[pointidx]].Location().Vector().X(), m_diffusepoints[diffuseindex[pointidx]].Location().Vector().Y(), m_diffusepoints[diffuseindex[pointidx]].Location().Vector().Z() );
				ok = ok && RotateRayToDiffuse( qobj.GetPoint(), qobj.GetLookAway(), inlook ); // Should have the same effect as querPtLocalCoords( qobj.GetLook(), inlook );
				ok = ok && m_radStorage->DiffuseSource( qobj, m_diffusepoints[diffuseindex[pointidx]], inlook.UnitVector(), source_temp );
				kscaatpoint            = m_kscatatpoints[diffuseindex[pointidx]];							// Get the extinction at the interpolation ponts
				if (kscaatpoint > 0)																		// IF the extinction is above zero
				{																				// then 
					source_temp *= diffuseweight[pointidx];										// we can do the interpolation
					source_temp *= kscat /kscaatpoint;											// Perform the linear interpolation using the ratio of the two cross
					source += source_temp;
				}
				else																			// Otherwise extinction is zero at leat one of the interpolation points
				{																				// so
					 skRTStokesVector::SetToZero(source);										// Set the source to zero
					 break;																		// and break out of the loop
				}
			}
		}
	}
    return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::GroundSource_impl		 2016- 11- 28*/
/** Returns the source term from interpolation of the ground points at the location 
 *	and in the  opposite direction to  the look vector stored in qobj. The answer is 
 *	returned as either a scalar or polarized radiance.
 **/
/*---------------------------------------------------------------------------*/

template< typename Radtype >
bool SKTRAN_HR_Diffuse_Table_CPU::GroundSource_impl( const SKTRAN_SourceTermQueryObject_Base& qobj, Radtype& source ) const
{
    bool ok = true;
    size_t								diffuseindex     [6];
    SKTRAN_HR_WEIGHT_TYPE			    diffuseweight    [6];

    Radtype source_temp;
    size_t numpoints;
    nxVector inlook;
    //	HELIODETIC_BASIS					basis(pt,look);

    ok = ok && ChooseGroundPoints(qobj.GetPoint(), diffuseindex, diffuseweight, numpoints );
    skRTStokesVector::SetToZero(source);
    double coszen = 1.0/nxmath::Pi;



    for( size_t pointidx = 0; pointidx < numpoints; pointidx++ )
    {
		// radstore_scalar just gets the look from the qobj, but for consistency set both to be the same
		inlook = -1 * nxVector(qobj.GetLookAway().X(), qobj.GetLookAway().Y(), qobj.GetLookAway().Z());

        ok = ok && m_radStorage->GroundSource( qobj, m_diffusepoints[diffuseindex[pointidx]], inlook.UnitVector(), source_temp );
        source_temp *= diffuseweight[pointidx];
		source += source_temp;
    }

    return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::DiffuseSource		 2016- 12- 5*/
/** Calculates the diffuse source ( radiance/m ) at a specified (atmospheric) location
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::DiffuseSource( const SKTRAN_SourceTermQueryObject_Base& qNoTransform,  double& source ) const
{
    return DiffuseSource_impl( qNoTransform,  source );
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::DiffuseSource		 2016- 12- 5*/
/** Calculates the diffuse source ( radiance/m ) at a specified (atmospheric) location
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::DiffuseSource( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, SKTRAN_Stokes_NC& source ) const
{
    return DiffuseSource_impl( qNoTransform,  source );
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::GroundSource		 2016- 12- 5*/
/** Calculates the diffuse source ( radiance/m ) at a specified (ground) location 
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::GroundSource( const SKTRAN_SourceTermQueryObject_Base& qNoTransform,  double& source ) const
{
    return GroundSource_impl( qNoTransform,  source );
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::GroundSource		 2016- 12- 5*/
/** Calculates the diffuse source ( radiance/m ) at a specified (ground) location
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::GroundSource( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, SKTRAN_Stokes_NC& source ) const
{
    return GroundSource_impl( qNoTransform,  source );
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::RotateRayToDiffuse		2013-07-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::RotateRayToDiffuse( const HELIODETIC_POINT& quadlocation, const HELIODETIC_UNITVECTOR& look, nxVector& rotatedray ) const
{
	bool ok = true;

	HELIODETIC_UNITVECTOR						localunitvectors[3];
	HELIODETIC_UNITVECTOR						locallook;

	ok = ok && quadlocation.LocalUnitVectors( localunitvectors, N_ELEMENTS(localunitvectors) );
	locallook = quadlocation.TransformToLocalZenithCoords( look, localunitvectors );
	rotatedray.SetCoords( -1.0*locallook.X(), -1.0*locallook.Y(), -1.0*locallook.Z() );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::ConfigureIndexes		2013-06-21*/
/** Called anytime after the diffuse points have been created, here we create
 *  the necessary internal indexes for the diffuse points, and allocate the
 *  sparse scattering matrix
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::ConfigureIndexes()
{
	bool ok = true;
	size_t incomingcounter = 0;
	size_t outgoingcounter = 0;
	size_t scattvalcounter = 0;
//	size_t Browidx         = 0;
	size_t numpoints = NumDiffusePoints();
	size_t numground = NumGroundPoints();

	for( size_t pointidx = 0; pointidx < numpoints; pointidx++ )
	{
        m_diffusepoints[pointidx].SetPointIndices( pointidx, incomingcounter, outgoingcounter, scattvalcounter, m_diffusepoints[pointidx].NumIncomingRays() );
		//m_diffusepoints[pointidx].m_isgroundpoint = false;
		incomingcounter += m_diffusepoints[pointidx].NumIncomingRays();
		outgoingcounter += m_diffusepoints[pointidx].NumOutGoingRays();
		scattvalcounter += m_diffusepoints[pointidx].NumIncomingRays() * m_diffusepoints[pointidx].NumOutGoingRays();
	}

	for( size_t groundidx = 0; groundidx < numground; groundidx++ )
	{
        m_diffusepoints[m_groundstartidx + groundidx].SetPointIndices( (m_groundstartidx + groundidx), incomingcounter, outgoingcounter, scattvalcounter, m_diffusepoints[m_groundstartidx + groundidx].NumIncomingRays() );
		//m_diffusepoints[m_groundstartidx + groundidx].m_isgroundpoint = true;
		incomingcounter += m_diffusepoints[m_groundstartidx + groundidx].NumIncomingRays();
		outgoingcounter += m_diffusepoints[m_groundstartidx + groundidx].NumIncomingRays();		// ground point no saves the downward flux of each incoming ray
		scattvalcounter += m_diffusepoints[m_groundstartidx + groundidx].NumIncomingRays();
	}

	ok = ok && m_Avals     ->AllocateStorage( incomingcounter, outgoingcounter, scattvalcounter, numpoints+numground );
	ok = ok && m_radStorage->AllocateStorage( incomingcounter, outgoingcounter );

	m_kscatatpoints.resize( NumDiffusePoints() + NumGroundPoints() );
	
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::SetProfileStartIdx		2013-07-11*/
/**  Internally sets where profiles begin/end in terms of indexes.  The last
 *   element of profilestartidx should be the number of non-ground diffuse points
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::SetProfileStartIdx( const std::vector<size_t>& profilestartidx )
{
	bool ok = true;
	size_t		numheights;
	size_t		count = 0;

	m_diffusealts.resize( profilestartidx.size() - 1 );
	for( size_t profileidx = 0; profileidx < profilestartidx.size()-1; profileidx++ )
	{
		numheights = profilestartidx[profileidx+1] - profilestartidx[profileidx];
		m_diffusealts[profileidx].resize( numheights );
		for( size_t altidx = 0; altidx < numheights; altidx++ )
		{
			m_diffusealts[profileidx][altidx] = m_diffusepoints[count].Location().Altitude();
			++count;
		}
	}

	m_groundstartidx = profilestartidx[profilestartidx.size()-1];
	m_profilestartidx = profilestartidx;

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::MakeGroundSourceDiffuseIndexesForRay		 2016- 12- 23*/
/** Adds the ground source contributions to this ray. The ground sources from
 *	the last iteration are internally stored as discrete downward fluxes (I.cos(theta))
 *	To calculate the groun source reflected/scattered radiance in any given direction
 *	we mnust integrate the downward fluxes over the hemispherical 2.pi and apply
 *	the solid angle weight of each downward flux and the BRDF weight of each downward
 *	flux as it scatters into this ray.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::MakeGroundSourceDiffuseIndexesForRay( const SKTRAN_RayOptical_Base*	ray, 
																		size_t							pointidx,
																	    size_t							rayidx,
																		SKTRAN_HR_Diffuse_Index_Array*	indexarray_lo)
																		
{
	bool							ok = true;
	double							brdf;
	HELIODETIC_POINT				groundpoint;
	double							factor;
	double							mu_in;
	double							mu_out;
	double							cosdphi;
	double							transmissionfactor;
	size_t							diffuseindex[6];
	SKTRAN_HR_WEIGHT_TYPE			diffuseweight[6];
	size_t							numpoints;
	size_t							profidx = 0;
	HELIODETIC_UNITVECTOR			localcoords[3];
	HELIODETIC_UNITVECTOR			outgoingray_helio;
	HELIODETIC_UNITVECTOR			outgoingray_local;
	nxVector						outgoinghoriz;
	bool							ok1;
	SKTRAN_HR_Diffuse_Index			tempindex;


	transmissionfactor = exp(-1.0*ray->TotalOpticalDepth());
	ok = ok && ray->Storage()->LocationOfPoint( ray->Storage()->NumCells(), &groundpoint);
	ok = ok && ChooseGroundPoints( groundpoint, diffuseindex, diffuseweight, numpoints );

	const HELIODETIC_UNITVECTOR& localzenith_helio = groundpoint.LocalZenith();		// Get the local zenith at the point where the ray hits the ground
	outgoingray_helio       = ray->LookVector();									// Get the look direction of the ray (this is into the ground and we need out of the ground)
	outgoingray_helio.Negate();														// so negate the direction so it is away from the ground point
	mu_out = (localzenith_helio & outgoingray_helio);								// Get the cosine of zenith angle of outbound ray, take negative so the ray direction is out of the ground point, not into it.
	NXASSERT( (mu_out >= 0.0) );													// and make sure we have it right

	groundpoint.LocalUnitVectors(localcoords,3);														// Get North, West and Up at the current groundpoint in heliodetic corodinates
	outgoingray_local =  groundpoint.TransformToLocalZenithCoords( outgoingray_helio, localcoords );	// Transform the outgoing ray to the local coordinate system (in the local system, up is always (0,0,1)
	outgoinghoriz.SetCoords( outgoingray_local.X(), outgoingray_local.Y(), 0.0);						// Get just the horizontal component as we will need this later for the azimuth determination
	outgoinghoriz.UnitVector();
	
	if( SKTRAN_HR_DECOUPLE_DIFFUSE_PROFILES )
	{
		profidx = (std::upper_bound( std::begin( m_profilestartidx ), std::end( m_profilestartidx ), pointidx ) - std::begin( m_profilestartidx ) ) -1;
		for( size_t idx = 0; idx < numpoints; idx++ )
		{
			diffuseindex[idx] = profidx + m_groundstartidx;
		}
	}

		//transmissionfactor *= albedo/nxmath::Pi;

	for( size_t interppointidx = 0; interppointidx < numpoints; interppointidx++ )
	{
		SKTRAN_HR_Diffuse_Point& point = m_diffusepoints[diffuseindex[interppointidx]];
		const HELIODETIC_POINT&  location = point.Location();
		
		const GEODETIC_INSTANT geo_location = m_coords->PointToGeodetic(location, m_coords->ReferencePointMJD());

		for ( size_t inidx = 0;  inidx < point.NumGroundDownwardFlux(); inidx++)							// For each incoming downward flux
		{																									// Loop over incoming directions
			nxVector	incoming (point.IncomingUnitRayLocalCoords( inidx ));								// Get the full incoming ray unit vector (away from ground point) in local coords	
			nxVector	incominghoriz( incoming.X(), incoming.Y(), 0);										// Get the horizontal component of the look vector 
			incominghoriz = incominghoriz.UnitVector();														// and get the horizontal component as a unit vector

			mu_in    = incoming.Z();																		// cosine of incoming ray zenith angle
			NXASSERT(( mu_in >= 0.0));																		// make sure it is acceptable
			cosdphi  = outgoinghoriz & incominghoriz;														// get thecosine of the azimuth angle between rays
			ok1      = m_opticaltable->GetBRDFGeodetic( geo_location, mu_in, mu_out, cosdphi, &brdf );		// Do the BRDF calculation
			factor   = (ok1) ? transmissionfactor*brdf : 0;													// How much source function does this downward flux add to the ray we are considering 
			tempindex.index  = SKTRAN_HR_SIZET_TO_INDEX( point.GroundDownwardFluxIdx(inidx) );				// Get the index
			tempindex.weight = SKTRAN_HR_DBL_TO_WEIGHT ( factor * diffuseweight[interppointidx] );			// and interpolating weight of thsi point
			NXASSERT(( NXFINITE(tempindex.weight)));
			if( abs(tempindex.weight) > 1E-8 )																// and if the value is significant
			{																								// then add
				indexarray_lo->diffindex.push_back( tempindex );											// it to our list.
			}
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::CreateDiffuseIndexesForRay		2013-06-27*/
/** Takes an already traced ray, and indexes indicating which diffuse point it belongs
 *  to and which incoming ray it is, and fills the diffuse points internal index table
 *  The index table contains information about which diffuse points and which rays
 *  will contribute to the next order incoming signal
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::CreateDiffuseIndexesForRay( const SKTRAN_RayOptical_Base* ray, size_t pointidx, size_t rayidx )
{
	bool ok = true;

	size_t								diffuseindex[6];
	SKTRAN_HR_WEIGHT_TYPE				diffuseweight[6];
	HELIODETIC_POINT					centerpt;
	size_t								numpoints = N_ELEMENTS(diffuseindex);
	double				transmissionfactor;
	double				k;
	double				kscat;
	double				ds;
	size_t				unitsphereindex[3];
	double				unitsphereweight[3];
	double				opticaldepth;
	nxVector			inray;
	nxVector			temp;
	nxVector			rotationunit;
	size_t				numquad = ray->GetNumQuadraturePoints();
	SKTRAN_HR_Diffuse_Index tempindex;

	SKTRAN_HR_Diffuse_Index_Array indexarray_lo; 
	SKTRAN_HR_Diffuse_Index_Array indexarray_hi;
	indexarray_lo.diffindex.reserve( numquad * 6 );
	indexarray_hi.diffindex.reserve( numquad * 6 );

	size_t profidx = 0;

	if( SKTRAN_HR_DECOUPLE_DIFFUSE_PROFILES )
	{
		if( m_groundstartidx <= pointidx )
			profidx = m_profilestartidx[pointidx - m_groundstartidx];
		else
			profidx = *(std::upper_bound( std::begin( m_profilestartidx ), std::end( m_profilestartidx ), pointidx )-1);
	}


	// add in ground contributions
	if( ray->Storage()->GroundIsHit() )
	{
		ok = MakeGroundSourceDiffuseIndexesForRay( ray, pointidx, rayidx, &indexarray_lo);
	}

	for( size_t quadidx = 0; quadidx < ray->GetNumCells(); quadidx++ )
	{
		ds           = ray->Storage()->CellLength( quadidx );
		ok           = ok && ray->Storage()->CellMidPoint( quadidx, &centerpt );
		opticaldepth = ray->OpticalDepthArray().at(quadidx);
		numpoints    = N_ELEMENTS(diffuseindex);
		ok           = ok && ChooseDiffusePoints( centerpt, diffuseindex, diffuseweight, &numpoints );
		k            = m_opticaltable->TotalExtinctionPerCM( centerpt ) * 100.0;
		kscat        = m_opticaltable->ScatteringExtinctionPerCM( centerpt );
		transmissionfactor = (k > 0) ? exp(-1.0*opticaldepth) * (1-exp(-1.0*k*ds))/k : 0.0;

		// need to rotate the diffuse point onto the quadrature location, or equivalently rotate
		// the incoming look vector to the outgoing sphere
		RotateRayToDiffuse(centerpt, ray->LookVector(), inray);

		if( SKTRAN_HR_DECOUPLE_DIFFUSE_PROFILES )
		{
			for( size_t idx = 0; idx < numpoints; idx++ )
			{
				size_t oldprofidx = *(std::upper_bound( std::begin( m_profilestartidx ), std::end( m_profilestartidx), diffuseindex[idx] )-1);
				diffuseindex[idx] =(size_t) ( (int)diffuseindex[idx]  + ((int)profidx - (int)oldprofidx));
			}
		}
		
		for( size_t interppointidx = 0; interppointidx < numpoints; interppointidx++ )
		{
			SKTRAN_HR_Diffuse_Index_Array* indexarray;
			if( (numpoints/2) <= interppointidx ){
				indexarray = &indexarray_hi;
			} else{
				indexarray = &indexarray_lo;
			}

			ok = ok && m_diffusepoints[diffuseindex[interppointidx]].TriangulateOnOutgoing( inray, unitsphereindex, unitsphereweight , 3);

			double kscaatpoint =  m_kscatatpoints[diffuseindex[interppointidx]];
			double factor      = (kscaatpoint > 0) ?  transmissionfactor * kscat/kscaatpoint : 0;
						
			SKTRAN_HR_WEIGHT_TYPE weightfactor =  diffuseweight[interppointidx] * SKTRAN_HR_DBL_TO_WEIGHT(factor);
			for( size_t sphereidx = 0; sphereidx < 3; sphereidx++ )
			{
				tempindex.index = SKTRAN_HR_SIZET_TO_INDEX( m_diffusepoints[diffuseindex[interppointidx]].OutgoingRadianceIdx(unitsphereindex[sphereidx]) );
				tempindex.weight = weightfactor * SKTRAN_HR_DBL_TO_WEIGHT(unitsphereweight[sphereidx]);
				NXASSERT(( NXFINITE(tempindex.weight)));
				if( abs(tempindex.weight) > 1E-8 )
				{
					// don't add weights that are essentially 0
					indexarray->diffindex.push_back( tempindex );
				}
			}
		}
	}


    //std::sort(indexarray_lo.diffindex.begin(), indexarray_lo.diffindex.end(), []( const SKTRAN_HR_Diffuse_Index& d1, const SKTRAN_HR_Diffuse_Index& d2){ return d1.index < d2.index;} );
    //std::sort(indexarray_hi.diffindex.begin(), indexarray_hi.diffindex.end(), []( const SKTRAN_HR_Diffuse_Index& d1, const SKTRAN_HR_Diffuse_Index& d2){ return d1.index < d2.index;} );
    std::sort(indexarray_lo.diffindex.begin(), indexarray_lo.diffindex.end() );
    std::sort(indexarray_hi.diffindex.begin(), indexarray_hi.diffindex.end() );
	
	SKTRAN_HR_Diffuse_Index_Array& indexarray_internal = m_diffusepoints[pointidx].DiffuseIndicesVar()[rayidx];
	indexarray_internal.diffindex.resize( indexarray_lo.diffindex.size() + indexarray_hi.diffindex.size() );
	std::merge(indexarray_lo.diffindex.begin(), indexarray_lo.diffindex.end(), indexarray_hi.diffindex.begin(), indexarray_hi.diffindex.end(), indexarray_internal.diffindex.begin(), []( const SKTRAN_HR_Diffuse_Index& d1, const SKTRAN_HR_Diffuse_Index& d2){ return d1.index < d2.index;} );

	ok = ok && m_diffusepoints[pointidx].CorrectForIndexDuplicates_Sorted( rayidx );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::CleanDiffuseIndexes		2013-06-27*/
/** Called after all rays have been added to the next order matrix, cleans 
 *  up the excess memory used in the diffuse indexes.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::CleanDiffuseIndexes()
{
	bool ok = true;

	for( size_t pointidx = 0; pointidx <  m_diffusepoints.size() ; pointidx++ )
	{
		m_diffusepoints[pointidx].CleanDiffuseIndex();
	}

	m_radStorage->CleanDiffuseIndexes( );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::DeclareFirstOrderInitialized		 2016- 12- 5*/
/** Table needs to be told when all point indices have been processed.**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::DeclareFirstOrderInitialized ()
{
	return m_radStorage->DeclareFirstOrderInitialized();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::DeclareAllScattered		 2016- 12- 5*/
/** Table needs to be told when all point indices have been processed.**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::DeclareAllScattered ()
{
	return m_radStorage->DeclareAllScattered();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::DeclareAllIntegrated		 2016- 12- 5*/
/** Table needs to be told when all point indices have been processed.**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_CPU::DeclareAllIntegrated ()
{
	return m_radStorage->DeclareAllIntegrated();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::GetDiffuseSourceOrder		 2016- 12- 5*/
/** How many scatters of the radiance initialized in #CalcFirstOrderIncomingPoint 
 *	are accounted for in #DiffuseSource? 
 **/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_HR_Diffuse_Table_CPU::GetDiffuseSourceOrder () const
{
	return m_radStorage->GetDiffuseSourceOrder();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU::PrintMemReport		 2016- 12- 5*/
/** **/
/*---------------------------------------------------------------------------*/

void SKTRAN_HR_Diffuse_Table_CPU::PrintMemReport() const
{

	size_t numindices = 0;
	for(int didx=0; didx < (int)m_diffusepoints.size(); ++didx ){
		for(int ridx=0; ridx < (int)(m_diffusepoints[didx].DiffuseIndices().size()); ++ridx){
			numindices += m_diffusepoints[didx].DiffuseIndices()[ridx].diffindex.size();
		}
	}

	printf("num indices and weights: %u\n", (unsigned int)numindices);
	printf("Press any key to continue...\n"); cin.get();
}

//void SKTRAN_HR_Diffuse_Table_CPU::SetOrder( size_t order ) 
//{
//	m_radStorage->SetOrder( order ); 
//	return;
//}


