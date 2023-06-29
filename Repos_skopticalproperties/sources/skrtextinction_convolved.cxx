#include <skopticalproperties21.h>

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_FixedFWHM::skOpticalProperties_Convolved_FixedFWHM		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_Convolved_FixedFWHM::skOpticalProperties_Convolved_FixedFWHM	( )
{
	m_fwhm_nm = 0.0;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_FixedFWHM::skOpticalProperties_Convolved_FixedFWHM		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_Convolved_FixedFWHM::skOpticalProperties_Convolved_FixedFWHM	( skOpticalProperties* highresextinction, double fwhm_nm, double quad_resol_nm )
{
	m_fwhm_nm = fwhm_nm;
	SetHighResolutionCrossSectionObject( highresextinction );
	SetQuadratureResolution( quad_resol_nm );							// Set the resolution for the integral
	UseGaussianQuadrature(false);										// By default use trapezoidal quadrature
}


/*-----------------------------------------------------------------------------
 *					~skOpticalProperties_Convolved_FixedFWHM		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_Convolved_FixedFWHM::~skOpticalProperties_Convolved_FixedFWHM	( )
{
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_FixedFWHM::Get_FWHM		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved_FixedFWHM::Get_FWHM( double /*lambda_nm*/, double* fwhm_nm ) const
{
	NXASSERT((m_fwhm_nm > 0.0));
	*fwhm_nm = m_fwhm_nm;
	return (m_fwhm_nm > 0.0);
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_FixedFWHM::DeepCopy		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved_FixedFWHM::DeepCopy( const skOpticalProperties_Convolved_FixedFWHM& other )
{
	bool	ok;

	m_fwhm_nm = other.m_fwhm_nm;
	ok = skOpticalProperties_Convolved::DeepCopy( other );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_FixedFWHM::CreateClone		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_Convolved_FixedFWHM::CreateClone( skOpticalProperties** clone) const
{
	skOpticalProperties_Convolved_FixedFWHM*	other;
	bool						ok;

	other = new skOpticalProperties_Convolved_FixedFWHM;
	ok =  (other != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_Convolved::CreateClone, Error allocating memory for clone object");
	}
	else
	{
		other->AddRef();
		ok = other->DeepCopy(*this);
		if (!ok)
		{
			other->Release();
			other = NULL;
		}
	}
	*clone = other;
	return ok;
}
*/


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_GaussQuadrature::skOpticalProperties_Convolved_GaussQuadrature		2009-11-24*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_Convolved_GaussQuadrature::skOpticalProperties_Convolved_GaussQuadrature()
{
	m_lambda0       = 0;
	m_sigma         = 0;
	m_normalizer    = 0;
	m_numsteps      = 0;
	m_absorptioncm2 = 0;
	m_extinctioncm2 = 0;
	m_scattercm2    = 0;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_GaussQuadrature::GetIntegrationRange		2009-11-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved_GaussQuadrature::GetIntegrationRange( double* minwavelen_nm, double* maxwavelen_nm )
{
	*minwavelen_nm  = m_lambda0 - 4.0*m_sigma;					// Get the lower bound for integration in nm
	*maxwavelen_nm  = m_lambda0 + 4.0*m_sigma;					// get the upper bound for integration in nm
	return true;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::Initialize		2009-11-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved_GaussQuadrature::Initialize( double lambda0, double instpsf_fwhm, double quadrature_resolution_nm )
{
	double	w1;
	double	w2;


	m_lambda0    = lambda0;
	m_sigma      = 0.42466090014400952136075141705144 * instpsf_fwhm;			// 0.4246609001 = 1.0/(sqrt(8*log(2))  Convert fwhm to sigma in exp-( 0.5*(x/sigma)^2 )
	m_normalizer = 1.0/(m_sigma*sqrt(nxmath::TWOPI));							// Get the normalization constant
	GetIntegrationRange( &w1, &w2);
	m_numsteps   = (size_t)( (w2-w1)/quadrature_resolution_nm + 1);
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_GaussQuadrature::Distribution		2009-11-24*/
/** **/
/*---------------------------------------------------------------------------*/

double skOpticalProperties_Convolved_GaussQuadrature::Distribution( double lambda )
{
	double	f;
	double	distribution;

	f                = (lambda - m_lambda0)/m_sigma;									// Get the normalized distance of this point from the central wavelength
	distribution     = exp(-0.5*f*f);											// Get the point spread function gaussian distribution at this wavelength (un-normalized)
	return distribution;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_GaussQuadrature::UpdateCrossSectionSummation		2009-11-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved_GaussQuadrature::UpdateCrossSectionSummation( double lambda, double weight, skOpticalProperties* highresextinction, double* sumcheck)
{
	double		factor;
	double		wavenum;
	bool		ok1;
	double		absxs,extxs,scattxs;

	factor	         = weight*Distribution(lambda);								// Get the factor
	wavenum          = 1.0E7/lambda;												// Get the wavenumber of this quadrature point
	ok1              = highresextinction->CalculateCrossSections( wavenum, &absxs, &extxs, &scattxs );		// Get the cross-sections at this wavenumber
	*sumcheck       += factor; 
	m_absorptioncm2 += absxs*factor;
	m_extinctioncm2 += extxs*factor;
	m_scattercm2    += scattxs*factor;
	return ok1;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_GaussQuadrature::ConvolveTrapezoid		2009-11-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved_GaussQuadrature::ConvolveTrapezoid( skOpticalProperties* highresextinction, double* sumcheck)
{
	size_t		N1;
	size_t		i;
	double		h;
	double		x1;
	double		xi;
	double		x2;
	bool		ok;
	
	*sumcheck = 0.0;
	GetIntegrationRange( &x1, &x2 );
	ok    = (m_numsteps > 1);
	if (ok)
	{
		N1         = m_numsteps - 1;
		h		   = (x2-x1)/N1;
		ok = UpdateCrossSectionSummation( x1, h*0.5, highresextinction,  sumcheck);
		for (i = 1; i < N1; i++)
		{
			xi = x1 + i*h;
			ok = ok && UpdateCrossSectionSummation( xi, h, highresextinction, sumcheck);
		}
		ok = ok && UpdateCrossSectionSummation( x2, h*0.5, highresextinction, sumcheck);
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_GaussQuadrature::ConvolveGaussQuad		2009-11-24*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved_GaussQuadrature::ConvolveGaussQuad( skOpticalProperties* highresextinction, double* sumcheck)
{
	bool						ok;
	const nx1dArray<double>*	lambda;
	const nx1dArray<double>*	weights;
	size_t						i;
	double						minw0;
	double						maxw0;
	
	*sumcheck = 0.0;
	GetIntegrationRange( &minw0, &maxw0 );
	m_gauss.SetOrder   ( (int)m_numsteps );								// Set the number of points fo rthe gaussian integration
	m_gauss.SetRange   ( minw0, maxw0 );					// Set the range for the integration 
	ok         = (m_gauss.GetOrder() == (int)m_numsteps);
	lambda	   = &m_gauss.X();
	weights	   = &m_gauss.W();
	for (i = 0; i < m_numsteps; i++)
	{
		ok = ok && UpdateCrossSectionSummation( lambda->At(i), weights->At(i), highresextinction, sumcheck);	
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *'					nxGaussQuadrature::Integrate		2003-9-26
 *	This integrates over the specified range using a function object
 *	that has a function operator () of double ()( double x)
 *	The function object evaluates the mathematical function at the value
 *	X and returns the value Y.
 *-------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved_GaussQuadrature::Convolve(double lambda0, double instpsf_fwhm, double quadratureresolution_nm, bool usegaussianquadrature, skOpticalProperties* highresextinction)
{
	bool		ok;
	double		sumcheck;

	m_absorptioncm2 = 0;
	m_scattercm2    = 0;
	m_extinctioncm2 = 0;

	ok = Initialize( lambda0, instpsf_fwhm, quadratureresolution_nm );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_Convolved_GaussQuadrature::Convolve, Error initializing the Gaussian quadrature. Thats a problem");
	}
	else
	{
		if (usegaussianquadrature) ok = ConvolveGaussQuad( highresextinction, &sumcheck );
		else                       ok = ConvolveTrapezoid( highresextinction, &sumcheck );
		
		sumcheck        *= m_normalizer;
		m_absorptioncm2 *= m_normalizer;
		m_extinctioncm2 *= m_normalizer;
		m_scattercm2    *= m_normalizer;
		if  (( sumcheck < 0.97) || (sumcheck > 1.03))
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_Convolved_GaussQuadrature::Convolve, Low precision quadrature detected (%e, 1.000 is the perfect value). You should definitely use a higher quadrature resolution (or fix the bug if its really bad)", (double)sumcheck );
			ok = false;
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_Convolved_GaussQuadrature::Convolve, There was an error performing the convolution, I am setting cross-sections to 0.0");
			m_absorptioncm2 = 0;
			m_scattercm2    = 0;
			m_extinctioncm2 = 0;
		}
	}
	return ok;
}




/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved_GaussQuadrature::DeepCopy		2009-11-24*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_Convolved_GaussQuadrature::DeepCopy( const skOpticalProperties_Convolved_GaussQuadrature& other )
{
	bool	ok = true;

//	ok              =       m_gauss.DeepCopy      ( other.m_gauss );				// We dont need to copy
//	ok              = ok && m_trapezoidal.DeepCopy( other.m_trapezoidal );			// as every call to Convolve initializes these objects as needed
	m_numsteps      = other.m_numsteps;
	m_lambda0       = other.m_lambda0;
	m_sigma			= other.m_sigma;
	m_normalizer	= other.m_normalizer;
	m_absorptioncm2	= other.m_absorptioncm2;
	m_extinctioncm2	= other.m_extinctioncm2;
	m_scattercm2	= other.m_scattercm2;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_Convolved_GaussQuadrature::DeepCopy, There was an error making a copy. Thats a problem,");
	}
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::skOpticalProperties_Convolved		2009-11-9*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_Convolved::skOpticalProperties_Convolved()
{
	CommonInit();
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::~skOpticalProperties_Convolved		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_Convolved::~skOpticalProperties_Convolved()
{
	if (m_highresextinction    != NULL) m_highresextinction->Release();
	if (m_backgroundatmosphere != NULL) m_backgroundatmosphere->Release();
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::CommonInit		2009-11-9*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_Convolved::CommonInit()
{
	m_quadratureresolution_nm = -99999.0;;								// By default use a value that creates problems. User must set this value.
	m_usegaussianquad		  = false;									// By default use trapezoidal quadrature
	m_highresextinction  = NULL;
	m_backgroundatmosphere = NULL;

}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::SetHighResolutionCrossSectionObject		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved::SetHighResolutionCrossSectionObject(  skOpticalProperties* extinctionobject )
{
	
	if (extinctionobject     != NULL) extinctionobject->AddRef();
	if ( m_highresextinction != NULL) m_highresextinction->Release();
	m_highresextinction = extinctionobject;
	if (( m_backgroundatmosphere !=NULL) && (m_highresextinction != NULL))
	{
		m_highresextinction->SetAtmosphericState( m_backgroundatmosphere );
	}
	return true;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::SetAtmosphericState		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved::SetAtmosphericState ( skClimatology* neutralatmosphere )
{
	if ( neutralatmosphere != nullptr) neutralatmosphere->AddRef();
	if ( m_backgroundatmosphere != nullptr) m_backgroundatmosphere->Release();
	m_backgroundatmosphere = neutralatmosphere;
	if (m_highresextinction != NULL)
	{
		m_highresextinction->SetAtmosphericState( neutralatmosphere);		// Set the atmospherc state for the cross-sections
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::SetAtmosphericState		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved::SetLocation ( const GEODETIC_INSTANT& pt,  bool* crosssections_changed )
{
	bool	haschanged = false;
	bool	ok;

	ok = (m_highresextinction != NULL) && (m_backgroundatmosphere != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_Convolved::SetAtmosphericLocation, You must a background atmosphere with a call to SetAtmosphericState and define a high resolution cross-section object with a call to SetHighResolutionCrossSectionObject");
		haschanged = false;
	}
	else
	{
		ok = m_highresextinction->SetLocation( pt, &haschanged );		// Set the atmospherc state for the cross-sections
	}
	if (crosssections_changed != NULL) *crosssections_changed = haschanged;	
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::InternalClimatology_UpdateCache		2011-8-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved::InternalClimatology_UpdateCache( const GEODETIC_INSTANT& pt)
{
	bool ok;

	ok = (m_highresextinction != NULL);
	ok = ok && m_highresextinction->InternalClimatology_UpdateCache( pt );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_Convolved::InternalClimatology_UpdateCache, error updating the climatology for the internal high resolution object");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::CalculateCrossSections		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs,  skOpticalProperties_Convolved_GaussQuadrature* quadrature ) const
{
	double	lambda;
	double	fwhm;
	bool	ok;

	ok = (m_quadratureresolution_nm	 > 0.0 );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_Convolved::CalculateCrossSections, The quadrature resolution has not been set. You must call SetQuadratureResolution before calling CalculateCrossSections");
	}
	else
	{
		lambda  = 1.0E7/wavenumber;
		ok      =       Get_FWHM( lambda, &fwhm );
		ok      = ok && quadrature->Convolve( lambda, fwhm, m_quadratureresolution_nm, m_usegaussianquad, m_highresextinction);
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_Convolved::CalculateCrossSections, Error performing the convolution, thats a problem");
		*extxs   = 0.0;
		*absxs   = 0.0;
		*scattxs = 0.0;
	}
	else
	{
		*extxs   = quadrature->Extinctioncm2();
		*absxs   = quadrature->Absorptioncm2(); 
		*scattxs = quadrature->Scattercm2();
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::LookupUpThreadData		 2014- 10- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved::LookupUpThreadData(skOpticalProperties_Convolved_GaussQuadrature** data )
{
	size_t					threadnum;
	iterator				iter;
	bool					ok;
	static std::mutex		lock;							// Add a lock to make sure map.insert is thread safe.

	threadnum = nxWorkerThreadManager::GetCurrentThreadIdCode();
	lock.lock();
	iter      = m_quadrature.find(threadnum);
	ok = (iter != m_quadrature.end() );
	if (!ok)
	{
		skOpticalProperties_Convolved_GaussQuadrature						blank;
		std::pair<iterator, bool>											result;
		std::pair<size_t, skOpticalProperties_Convolved_GaussQuadrature>	newentry(threadnum, blank);

		result = m_quadrature.insert(  newentry);
		ok = result.second;
		iter = result.first;

	}
	lock.unlock();
	if (!ok)
	{
		*data = nullptr;
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_Convolved::LookupUpThreadData, error looking/creating thread data for thread %d", (int)threadnum);
	}
	else
	{
		*data = &(*iter).second;
	}
	return ok;

}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::CalculateCrossSections		2014-2-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs )
{
	bool											ok;
	skOpticalProperties_Convolved_GaussQuadrature*	quadrature;

	ok = LookupUpThreadData( &quadrature );
	ok = ok && CalculateCrossSections( wavenumber, absxs, extxs, scattxs, quadrature );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_Convolved::CalculateCrossSections, Error calculating convolved cross-sections returning NaNs");
		*absxs   = std::numeric_limits<double>::quiet_NaN();
		*extxs   = std::numeric_limits<double>::quiet_NaN();
		*scattxs = std::numeric_limits<double>::quiet_NaN();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::IsScatterer		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved::IsScatterer() const
{
	bool	isscatter = false;

	if (m_highresextinction != NULL) isscatter = m_highresextinction->IsScatterer();
	return isscatter;

}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::IsAbsorber		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved::IsAbsorber() const
{
	bool	isabsorber = false;

	if (m_highresextinction != NULL) isabsorber = m_highresextinction->IsAbsorber();
	return isabsorber;

}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::CalculatePhaseMatrix		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_Convolved::CalculatePhaseMatrix( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix)
{
	bool	ok;

	ok = (m_highresextinction != NULL);
	if (!ok)
	{
		NXTRACE_ONCEONLY(firsttime, ("***** skOpticalProperties_Convolved::CalculatePhaseMatrix, only calculatea phase matrix at one wavelength. There is no convolution calculation\n"));
		ok = m_highresextinction->CalculatePhaseMatrix( wavenumber, cosscatterangle, phasematrix);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_Convolved::CalculatePhaseMatrix, error calculating the phase matrix form the high resolution spectrum");
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_Convolved::DeepCopy		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_Convolved::DeepCopy( const skOpticalProperties_Convolved& other )
{
	bool	ok;

	if (m_highresextinction != NULL)
	{
		m_highresextinction->Release();
		m_highresextinction = NULL;
	}

	ok        =        skOpticalProperties::DeepCopy( other );			// Copy over the base class
	ok        = ok &&  m_quadrature.DeepCopy(other.m_quadrature );				// Copy over the quadrature object
	if (other.m_highresextinction != NULL)										// Copy over the high resolution spectral object
	{
		ok = ok && other.m_highresextinction->CreateClone( &m_highresextinction );	
	}
	m_isdirty = other.m_isdirty;												// Copy ove rthe dirty flag
	m_quadratureresolution_nm = other.m_quadratureresolution_nm;				// Copy over the integration resolution
	m_usegaussianquad		  = other.m_usegaussianquad;						// Copy over the gaussian quadrature flag

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_Convolved::DeepCopy, Error copying the object");
	}
	return ok;
}
*/

