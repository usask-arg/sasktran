#include <skopticalproperties21.h>
#include <boost/math/special_functions/gamma.hpp>



/*-----------------------------------------------------------------------------
 *					skRTParticleDist::ASA_To_N		 2015- 9- 25*/
/** Convert an aerosol surface area parameter in um2/cm3 to number of 
 *	particles per cm3. The conversion is distribution specific. 
 **/
/*---------------------------------------------------------------------------*/

double skRTParticleDist::ASA_To_N(double /*asa_um2percm3*/) const
{
	nxLog::Record(NXLOG_WARNING,"skRTParticleDist::ASA_To_N, This particle distribution does not yet support conversion between aerosol surface area and number density");
	return std::numeric_limits<double>::quiet_NaN();
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist::N_To_ASA		 2015- 9- 25*/
/** Convert an number density (particles/cm3) to aerosol surface area
 *	um2percm3. The conversion is distribution specific. 
 **/
/*---------------------------------------------------------------------------*/

double skRTParticleDist::N_To_ASA(double /*N_percm3*/) const
{
	nxLog::Record(NXLOG_WARNING,"skRTParticleDist::N_To_ASA, This particle distribution does not yet support conversion between number density and aerosol surface area");
	return std::numeric_limits<double>::quiet_NaN();
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist::IsSameDistributionAs		2013-6-26*/
/** used to compare if two particle distributions are identical. This is used
 *	to see if two distributions are identical so we can avoid cache hits.
**/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist::IsSameDistributionAs( const skRTParticleDist* other ) const
{
	bool		ok;
	double		buffer1[20];
	double		buffer2[20];
	size_t		numparam1 = 0;
	size_t		numparam2 = 0;

//	ok = (other != NULL) && typeid(*this) == typeid(*other);
	ok = (other != NULL) && ( Type() == other->Type());
	if (ok)
	{
		ok =		(      NumDistributionParameters() <= N_ELEMENTS(buffer1));
		ok = ok &&  (other->NumDistributionParameters() <= N_ELEMENTS(buffer2));
		ok = ok &&        GetDistributionParameterArray( buffer1, N_ELEMENTS(buffer1), &numparam1 );
		ok = ok &&  other->GetDistributionParameterArray( buffer2, N_ELEMENTS(buffer2), &numparam2 );
		ok = ok && (numparam1 == numparam2);
		if (ok)
		{
			for (size_t idx=0; ok && (idx < numparam1); idx++ )
			{
				ok = ok && ( abs(buffer1[idx] - buffer2[idx]) < 1e-10 );
			}
		}
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					skRTParticleDist_2Gamma::skRTParticleDist_2Gamma		2003-10-8
 *-------------------------------------------------------------------------*/


skRTParticleDist_2Gamma::skRTParticleDist_2Gamma()
{
	m_A    = 0;
	m_B    = 0;
	m_AB   = 0;
	m_RBM3 = 0;
	m_NORM = 0.0;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_2Gamma::SetDistributionParameters		2013-1-28*/
/** This class implements the 2 parameter gamma distribution. See Eqn 2.56 in
 *	Hansen And Jensen. Space Sci. Rev., 16. 527-610, 1974.
 *
 *	In this formulation A is the effective radius of the gamma distribution and
 *	B is effective variance.
 **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_2Gamma::SetDistributionParameters( double A, double B, double )
{
	bool	ok;
	double	RB;
	double  RBM2;
	double	DLAB;
	double  DLGM;

	m_A = A;
	m_B = B;

	ok = (m_B < 0.5);
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "skRTParticleDist_2Gamma::Distribution, The B (%g)for the Gamma distribution must be less than 0.5", (double) m_B);
		m_AB   = 0.0;
		m_RBM3 = 0.0;
		m_NORM = 0.0;
	}
	else
	{
		RB     = 1.0 / m_B;
		RBM2   = RB - 2.0;
		m_RBM3 = RB - 3.0;
		m_AB   = RB  / m_A;
		DLAB   = RBM2 * log(m_AB);
        DLGM   = lgamma(RBM2);
		m_NORM = DLAB-DLGM;							// Get the normalization constant
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					skRTParticleDist_2Gamma::GetParameters		2009-6-24
 *-------------------------------------------------------------------------*/

void skRTParticleDist_2Gamma::GetDistributionParameters( double* A, double* B, double* C  )
{
	*A = m_A;
	*B = m_B;
	*C = 0;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_PowerLaw::GetDistributionParameterArray		2013-6-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_2Gamma::GetDistributionParameterArray( double* parameters, size_t maxparams, size_t* numparams ) const
{
	bool	ok;
	ok = (maxparams >= 2);
	if (ok)
	{
		parameters[0] = m_A;
		parameters[1] = m_B;
		*numparams    = 2;
	}
	else
	{
		*numparams = 0;
		nxLog::Record(NXLOG_WARNING,"skRTParticleDist_2Gamma::GetDistributionParameterArray, Error returning parameters as user buffer is smaller than 2 elements. Thats a problem");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::GetDistributionParameterSpeciesID		2013-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_2Gamma::GetDistributionParameterSpeciesID ( CLIMATOLOGY_HANDLE* paramids, size_t maxparams, size_t* numparams) const
{
	bool	ok;
	ok = (maxparams >= 2);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skRTParticleDist_2Gamma::GetDistributionParameterSpeciesID, you need to pass in at least a buffer of 2 points");
		*numparams = 0;
	}
	else
	{
		paramids[0] = SKCLIMATOLOGY_GAMMA_EFFECTIVERADIUS_MICRONS;				// Parameter A is Gamma Distribution Effective Radius in microns
		paramids[1] = SKCLIMATOLOGY_GAMMA_EFFECTIVEVARIANCE_PERMICRON;			// Parameter B is Gama distribution Effective variance per micron
		*numparams  = 2;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_2Gamma::CachingDescriptor		2008-4-16*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skRTParticleDist_2Gamma::CachingDescriptor() const
{
	nxString	descriptor;
	size_t		RE,GAM;
	RE	= (size_t)(10000.0*m_A + 0.5);		// Get mode radius in 0.1 nanometers
	GAM	= (size_t)(1000.0*m_B + 0.5);		// Get 1000 times gamma
	descriptor.sprintf("gamma2/g2_%05u_%05u",(unsigned int)RE, (unsigned int)GAM);
	return descriptor;
}

/*---------------------------------------------------------------------------
 *'					skRTParticleDist_2Gamma::Density		2003-10-8
 *	Returns the particle distribution density at radius R.
 *-------------------------------------------------------------------------*/

double skRTParticleDist_2Gamma::Distribution( double R )
{
	NXASSERT( (m_B > 0.0) && (m_B < 0.5) );
	return exp( m_NORM +  m_RBM3*log(R) - m_AB*R);
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_2Gamma::Copy		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_2Gamma::Copy( const skRTParticleDist_2Gamma& other)
{

	m_A     = other.m_A;
	m_B     = other.m_B;
	m_AB    = other.m_AB;
	m_NORM  = other.m_NORM;
	m_RBM3  = other.m_RBM3;
	return true;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_2Gamma::GetQuadratureRadii		2009-7-9*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_2Gamma::GetQuadratureRadii(double* minradius, double* maxradius) const
{
	bool	ok;

	ok = ( (m_B > 0.0) && (m_B < 0.5) );
	if (ok)
	{
		*minradius = 0;
		*maxradius = 5.0*m_A;
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"skRTParticleDist_2Gamma::GetQuadratureRadii, cannot do quadrature stuff until coefficients are set");
		*minradius = 0.0;
		*maxradius = 0.0;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_2Gamma::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/
bool skRTParticleDist_2Gamma::CreateClone(skRTParticleDist** userclone) const
{
	skRTParticleDist_2Gamma*	clone;
	bool						ok;

	clone = new skRTParticleDist_2Gamma;
	ok    = (clone != NULL);
	if (ok)
	{
		clone->AddRef();
		ok = clone->Copy(*this);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skRTParticleDist_2Gamma::CreateClone, Error copying to the cloned object");
		}
	}
	*userclone = clone;
	return ok;
}

/*---------------------------------------------------------------------------
 *'					skRTParticleDist_3Gamma::skRTParticleDist_2Gamma		2003-10-8
 *-------------------------------------------------------------------------*/

skRTParticleDist_3Gamma::skRTParticleDist_3Gamma()
{
	m_A    = 0;
	m_B    = 0;
	m_C    = 0;
	m_DLB  = 0;
}

/*---------------------------------------------------------------------------
 *'					skRTParticleDist_3Gamma::SetParameters		2003-10-8
 *-------------------------------------------------------------------------*/

bool skRTParticleDist_3Gamma::SetDistributionParameters( double A, double B, double C )
{
	bool	ok;
	double	RAC;
	double	RRAC;
	double	DLGN;

	m_A = A;
	m_B = B;
	m_C = C;
	ok = (m_C != 0.0);
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "skRTParticleDist_3Gamma::Distribution, The C (%g) parameter must not be zero", (double)C);
		m_DLB = 0.0;
	}
	else
	{
		RAC   = (m_A + 1.0) / m_C;
		RRAC  = RAC;
        DLGN  = lgamma(RRAC);
		m_DLB = RAC * log(m_B) - DLGN;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_3Gamma::GetDistributionParameterArray		2013-6-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_3Gamma::GetDistributionParameterArray( double* parameters, size_t maxparams, size_t* numparams ) const
{
	bool	ok;
	ok = (maxparams >= 3);
	if (ok)
	{
		parameters[0] = m_A;
		parameters[1] = m_B;
		parameters[2] = m_C;
		*numparams    = 2;
	}
	else
	{
		*numparams = 0;
		nxLog::Record(NXLOG_WARNING,"skRTParticleDist_3Gamma::GetDistributionParameterArray, Error returning parameters as user buffer is smaller than 2 elements. Thats a problem");
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					skRTParticleDist_3Gamma::GetParameters		2009-6-24
 *-------------------------------------------------------------------------*/

void skRTParticleDist_3Gamma::GetDistributionParameters( double* A, double* B, double* C )
{
	*A = m_A;
	*B = m_B;
	*C = m_C;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_3Gamma::GetDistributionParameterSpeciesID		2013-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_3Gamma::GetDistributionParameterSpeciesID ( CLIMATOLOGY_HANDLE* /*paramids*/, size_t /*maxparams*/, size_t* numparams) const
{
	nxLog::Record(NXLOG_ERROR,"skRTParticleDist_3Gamma::GetDistributionParameterSpeciesID, Needs to be implemented");
	*numparams = 0;
	return false;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_3Gamma::CachingDescriptor		2008-4-16*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skRTParticleDist_3Gamma::CachingDescriptor() const
{
	nxString descriptor;

	descriptor.sprintf("gamma3_%0#12g_%0#12g_%0#12g",(double)m_A, (double)m_B, (double)m_C);
	return descriptor;
}
/*---------------------------------------------------------------------------
 *'					skRTParticleDist_3Gamma::Density		2003-10-8
 *	Returns the particle distribution density at radius R.
 *-------------------------------------------------------------------------*/

double skRTParticleDist_3Gamma::Distribution( double R )
{
	return exp( m_DLB +  m_A *log(R) - pow( m_B*R, m_C) )*m_C;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_3Gamma::GetQuadratureRadii		2008-5-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_3Gamma::GetQuadratureRadii	(double* /*minradius*/, double* /*maxradius*/) const
{
	nxLog::Record( NXLOG_ERROR, "skRTParticleDist_3Gamma::GetQuadratureRadii, Not yet implemented");
	NXASSERT(( false ));
	return false;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_3Gamma::Copy		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_3Gamma::Copy( const skRTParticleDist_3Gamma& other)
{
	m_A    = other.m_A;
	m_B    = other.m_B;
	m_C    = other.m_C;
	m_DLB  = other.m_DLB;
	return true;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_3Gamma::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_3Gamma::CreateClone(skRTParticleDist** userclone) const
{
	skRTParticleDist_3Gamma*	clone;
	bool						ok;

	clone = new skRTParticleDist_3Gamma;
	ok = (clone != NULL);
	if (ok)
	{
		clone->AddRef();
		ok = clone->Copy(*this);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skRTParticleDist_3Gamma::CreateClone,Error copying object to clone");
		}
	}
	*userclone = clone;
	return ok;

}


/*---------------------------------------------------------------------------
 *'					skRTParticleDist_BimodalGamma::SetParameters		2003-10-9
 *-------------------------------------------------------------------------*/

bool skRTParticleDist_BimodalGamma::SetDistributionParameters( double A, double B, double C )
{
	double	RB;
	double	RBM2;
	double	DLGM;
	bool	ok;

	RB = 1.0/B;
	ok = (B >= 0.5);
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "skRTParticleDist_BimodalGamma::Distribution, The B (%g) parameter must not be zero", (double)B);
		m_A1B   = 0.0;
		m_A2B   = 0.0;
		m_DLA1B = 0.0;
		m_DLA2B = 0.0;
		m_RBM3  = 0.0;
	}
	RBM2    = RB - 2.0;
    DLGM    = lgamma(RBM2);
	m_A1B   = RB / A;
	m_A2B   = RB / C;
	m_DLA1B = RBM2 * log(m_A1B) - DLGM;
	m_DLA2B = RBM2 * log(m_A2B) - DLGM;
	m_RBM3 = RB - 3.0;
	return ok;
}


/*---------------------------------------------------------------------------
 *'				skRTParticleDist_BimodalGamma::GetParameters		2009-6-24
 *-------------------------------------------------------------------------*/

void skRTParticleDist_BimodalGamma::GetDistributionParameters( double* /*A*/, double* /*B*/, double* /*C*/ )
{
	nxLog::Record(NXLOG_WARNING,"skRTParticleDist_BimodalGamma::GetDistributionParameters, not supported . Use GetDistributionParameterArray");
	throw("skRTParticleDist_BimodalGamma::GetDistributionParameters Unsupported");
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_BimodalGamma::GetDistributionParameterArray		2013-6-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_BimodalGamma::GetDistributionParameterArray( double* parameters, size_t maxparams, size_t* numparams ) const
{
	bool	ok;
	ok = (maxparams >= 5);
	if (ok)
	{
		parameters[0] = m_A1B;
		parameters[1] = m_A2B;
		parameters[2] = m_DLA1B;
		parameters[3] = m_DLA2B;
		parameters[4] = m_RBM3;
		*numparams    = 5;
	}
	else
	{
		*numparams = 0;
		nxLog::Record(NXLOG_WARNING,"skRTParticleDist_BimodalGamma::GetDistributionParameterArray, Error returning parameters as user buffer is smaller than 2 elements. Thats a problem");
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skRTParticleDist_BimodalGamma::CachingDescriptor		2013-6-27*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skRTParticleDist_BimodalGamma::CachingDescriptor	() const
{
	double		A;
	double		B;
	double		C;
	size_t		Ai;
	size_t		Bi;
	size_t		Ci;
	nxString	descriptor;

	B = 1.0/( m_RBM3 + 3.0 );
	A = 1.0/( m_A1B*B );
	C = 1.0/( m_A2B*B );

	Ai = (size_t)(1000.0*A + 0.5);
	Bi = (size_t)(1000.0*B + 0.5);
	Ci = (size_t)(1000.0*C + 0.5);

	descriptor.sprintf("bimodal_%05u_%05u_%05u",(unsigned int)Ai, (unsigned int)Bi, (unsigned int)Ci );
	return descriptor;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_BimodalGamma::GetDistributionParameterSpeciesID		2013-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_BimodalGamma::GetDistributionParameterSpeciesID	( CLIMATOLOGY_HANDLE* /*paramids*/, size_t /*maxparams*/, size_t* numparams) const
{
	nxLog::Record(NXLOG_WARNING,"skRTParticleDist_BimodalGamma::GetDistributionParameterSpeciesID, not yet implemeted");
	*numparams = 0;
	return false;
}

/*---------------------------------------------------------------------------
 *'					skRTParticleDist_BimodalGamma::Distribution		2003-10-9
 *-------------------------------------------------------------------------*/

double skRTParticleDist_BimodalGamma::Distribution( double R )
{
	double	D1;
	double	D2;
	double	logR = log(R);
	D1 = exp( m_DLA1B + m_RBM3*logR - m_A1B*R);
	D2 = exp( m_DLA2B + m_RBM3*logR - m_A2B*R);
	return 0.5*(D1 + D2);
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_BimodalGamma::GetQuadratureRadii		2008-5-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_BimodalGamma::GetQuadratureRadii	(double* /*minradius*/, double* /*maxradius*/) const
{
	nxLog::Record( NXLOG_ERROR, "skRTParticleDist_BimodalGamma::GetQuadratureRadii, Not yet implemented");
	NXASSERT(( false ));
	return false;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_BimodalGamma::Copy		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_BimodalGamma::Copy( const skRTParticleDist_BimodalGamma& other)
{

	m_A1B    = other.m_A1B;
	m_A2B    = other.m_A2B;
	m_DLA1B  = other.m_DLA1B;
	m_DLA2B  = other.m_DLA2B;
	m_RBM3   = other.m_RBM3;

	return true;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_BimodalGamma::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_BimodalGamma::CreateClone(skRTParticleDist** userclone) const
{
	skRTParticleDist_BimodalGamma*	clone;
	bool							ok;

	clone = new skRTParticleDist_BimodalGamma;
	ok = (clone != NULL);
	if (ok)
	{
		clone->AddRef();
		clone->Copy(*this);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skRTParticleDist_BimodalGamma::CreateClone, error copying this object to clone");
		}
	}
	*userclone = clone;
	return ok;

}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::CachingDescriptor		2008-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skRTParticleDist_LogNormal::CachingDescriptor	() const
{
	size_t	RG;
	size_t	SIGMA;
	nxString descriptor;
	
	RG    = (size_t)(10000.0*exp(m_A2) + 0.5);					// Get mode radius in 0.1 nanometers
	SIGMA = (size_t)(1000.0* exp(sqrt(0.5/m_A3 )) + 0.5);		// Get 1000 times mode width 
	
	descriptor.sprintf("lognormal/ln_%05u_%05u",(unsigned int)RG, (unsigned int)SIGMA );
	return descriptor;
}

/*---------------------------------------------------------------------------
 *'					skRTParticleDist_LogNormal::SetParameters		2003-10-9
 *	Note RG comes in as mode radius
 *	and MODEWIDTH comes in as is modewidth (>= 1.0);
 *	A little clumsy but provides backward compatibility, so be careful if you
 *	decide to change it.
 *-------------------------------------------------------------------------*/

bool skRTParticleDist_LogNormal::SetDistributionParameters( double RG, double MODEWIDTH, double )
{
	double	SIGMA;

	if ((RG == 0) || (MODEWIDTH <= 1.0))
	{
		m_QTEMP = 0.0;
		m_A2    = 0;
		m_A3    = 0.0;
		if ((MODEWIDTH <= 1.0) && (RG != 0.0 ))
		{
			nxLog::Record( NXLOG_WARNING, "skRTParticleDist_LogNormal::SetDistributionParameters, The mode width must be 1.0 or greater, make sure you are not passing in the natural log of mode width (we dont want that)");
		}
	}
	else
	{
		SIGMA = log(MODEWIDTH);
		NXASSERT((SIGMA > 0.0));
		m_QTEMP =sqrt( nxmath::TWOPI )*SIGMA;
		m_A2    = log(RG);
		m_A3    = 1.0/(2.0*SIGMA*SIGMA);
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::ModeRadiusMicrons		 2015- 9- 25*/
/** Get the mode radius in microns**/
/*---------------------------------------------------------------------------*/

double	skRTParticleDist_LogNormal::ModeRadiusMicrons() const
{
	double RG;
	RG  = (m_A3 > 0.0) ? exp(m_A2) : std::numeric_limits<double>::quiet_NaN();
	return RG;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::ModeWidth		 2015- 9- 25*/
/** Get the mode width**/
/*---------------------------------------------------------------------------*/

double	skRTParticleDist_LogNormal::ModeWidth() const
{
	double MODEWIDTH;

	MODEWIDTH  = (m_A3 > 0.0) ? exp(sqrt(0.5/m_A3 )) : std::numeric_limits<double>::quiet_NaN();
	return MODEWIDTH;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::Reff		 2015- 9- 25*/
/** Return the effective radius o fthe lognormal distribution
 *	See equation 2.53 of
 *  @ARTICLE{1974SSRv...16..527H,
 *  author = {{Hansen}, J.~E. and {Travis}, L.~D.},
 *  title = "{Light scattering in planetary atmospheres}",
 *	journal = {\ssr},
 *	keywords = {Atmospheric Optics, Atmospheric Scattering, Light Scattering, Optical Reflection, Planetary Atmospheres, Albedo, Angular Distribution, Integral Equations, Invariant Imbeddings, Mie Scattering, Monte Carlo Method, Particle Size Distribution, Polarization Characteristics, Rayleigh Scattering, Refractivity},
 *  year = 1974,
 *  month = oct,
 *  volume = 16,
 *  pages = {527-610},
 *  doi = {10.1007/BF00168069},
 *  adsurl = {http://adsabs.harvard.edu/abs/1974SSRv...16..527H},
 *	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
 *	}
 *
 *  and
 *
 *	See equation 7a of
 *  C. A. McLinden, J. C. McConnell, C. T. McElroy, and E. Griffioen, 1999:
 *	Observations of Stratospheric Aerosol Using CPFM Polarized Limb Radiances. J. Atmos. Sci., 56, 233�240.
 *	doi: http://dx.doi.org/10.1175/1520-0469(1999)056<0233:OOSAUC>2.0.CO;2
 *
**/
/*---------------------------------------------------------------------------*/

double skRTParticleDist_LogNormal::Reff() const
{
	double rg = ModeRadiusMicrons();
	double sg = log(ModeWidth());
	double reff;

	reff = rg*exp( 5.0*sg*sg/2.0);
	return reff;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::Veff		 2015- 9- 25*/
/**	Return the effective variance of the distribution
 *	See equation 2.54 of
 *  @ARTICLE{1974SSRv...16..527H,
 *  author = {{Hansen}, J.~E. and {Travis}, L.~D.},
 *  title = "{Light scattering in planetary atmospheres}",
 *	journal = {\ssr},
 *	keywords = {Atmospheric Optics, Atmospheric Scattering, Light Scattering, Optical Reflection, Planetary Atmospheres, Albedo, Angular Distribution, Integral Equations, Invariant Imbeddings, Mie Scattering, Monte Carlo Method, Particle Size Distribution, Polarization Characteristics, Rayleigh Scattering, Refractivity},
 *  year = 1974,
 *  month = oct,
 *  volume = 16,
 *  pages = {527-610},
 *  doi = {10.1007/BF00168069},
 *  adsurl = {http://adsabs.harvard.edu/abs/1974SSRv...16..527H},
 *	adsnote = {Provided by the SAO/NASA Astrophysics Data System}
 *	}
 *
 *  and
 *
 *	See equation 7b of
 *  C. A. McLinden, J. C. McConnell, C. T. McElroy, and E. Griffioen, 1999:
 *	Observations of Stratospheric Aerosol Using CPFM Polarized Limb Radiances. J. Atmos. Sci., 56, 233�240.
 *	doi: http://dx.doi.org/10.1175/1520-0469(1999)056<0233:OOSAUC>2.0.CO;2
**/
/*---------------------------------------------------------------------------*/

double skRTParticleDist_LogNormal::Veff() const
{
	double sg = log(ModeWidth());
	double veff;

	veff = exp(sg*sg) - 1;
	return veff;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::ASA_To_N		 2015- 9- 25*/
/** Calculate the number of aerosol particles per cm3 given the aerosol
 *	surface area in um2 per cm3. 
 *
 *	See equation 11a of
 *  C. A. McLinden, J. C. McConnell, C. T. McElroy, and E. Griffioen, 1999:
 *	Observations of Stratospheric Aerosol Using CPFM Polarized Limb Radiances. J. Atmos. Sci., 56, 233�240.
 *	doi: http://dx.doi.org/10.1175/1520-0469(1999)056<0233:OOSAUC>2.0.CO;2
 **/
/*---------------------------------------------------------------------------*/

double skRTParticleDist_LogNormal::ASA_To_N(double asa_um2percm3) const
{
	double N_percm3;
	double	reff = Reff();
	double  veff = Veff();

	N_percm3 = asa_um2percm3*pow( (veff+1), 3)/( 4.0*nxmath::Pi*reff*reff);
	return N_percm3;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::N_To_ASA		 2015- 9- 25*/
/** Calculate the aerosol surface area in um2 per cm3 given the number of
 *	aerosol particles per cm3.
 *
 *	See equation 11a of
 *  C. A. McLinden, J. C. McConnell, C. T. McElroy, and E. Griffioen, 1999:
 *	Observations of Stratospheric Aerosol Using CPFM Polarized Limb Radiances. J. Atmos. Sci., 56, 233�240.
 *	doi: http://dx.doi.org/10.1175/1520-0469(1999)056<0233:OOSAUC>2.0.CO;2
**/
/*---------------------------------------------------------------------------*/

double skRTParticleDist_LogNormal::N_To_ASA(double N_percm3) const
{
	double asa_um2percm3;
	double	reff = Reff();
	double  veff = Veff();

	asa_um2percm3 = 4.0*nxmath::Pi*reff*reff*N_percm3/pow( (veff+1), 3);
	return asa_um2percm3;
}


/*---------------------------------------------------------------------------
 *'					skRTParticleDist_LogNormal::GetParameters		2009-6-24
 *	Return the mode radius in microns and and teh mode width
 *
 *	\param A
 *	returns the mode radius in microns
 *
 *	\param B
 *	Returns the mode width
 *-------------------------------------------------------------------------*/

void skRTParticleDist_LogNormal::GetDistributionParameters( double* RG, double* MODEWIDTH, double* C )
{
	*RG        = ModeRadiusMicrons();
	*MODEWIDTH = ModeWidth();
	*C         = 0;
}



/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::GetDistributionParameterArray		2013-6-27*/
/** Returns the A2 and A3 coefficients which are derived from mode radius
 *	and mode width. This function is often used to see if two particle distributions
 *	are identical. Therefore it makes sense not to convert back to mode radius and mode width
 **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_LogNormal::GetDistributionParameterArray( double* parameters, size_t maxparams, size_t* numparams ) const
{
	bool	ok;

	ok = (maxparams >= 2);
	if ( ok)
	{
		parameters[0] = m_A2;
		parameters[1] = m_A3;
		*numparams     = 2;
	}
	else
	{
		numparams = 0;
		nxLog::Record(NXLOG_WARNING,"skRTParticleDist_LogNormal::GetDistributionParameterArray, User buffer is too small to hold the 2 parameters. Thats a problem");
	}
	return ok;
}
	

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::GetDistributionParameterSpeciesID		2013-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_LogNormal::GetDistributionParameterSpeciesID ( CLIMATOLOGY_HANDLE* paramids, size_t maxparams, size_t* numparams) const
{
	bool	ok;
	ok = (maxparams >= 2);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skRTParticleDist_LogNormal::GetDistributionParameterSpeciesID, you need to pass in at least a buffer of 2 points");
		*numparams = 0;
	}
	else
	{
		paramids[0] = SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS;				// Parameter A is Mode Radius
		paramids[1] = SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH;						// Parameter B is mode width
		*numparams  = 2;
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					skRTParticleDist_LogNormal::Distribution		2003-10-9
 *-------------------------------------------------------------------------*/

double skRTParticleDist_LogNormal::Distribution( double R )
{
	double	A1;
	double d = 0.0;

	if ( m_QTEMP > 0)
	{
		A1 = log(R);
		d  = exp( - nxmath::sqr(A1 - m_A2)*m_A3) / (m_QTEMP*R);
	}
	//else
	//{
	//	d = 0.0;
	//}
	return d;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::GetQuadratureRadii		2008-5-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_LogNormal::GetQuadratureRadii	(double* minradius, double* maxradius) const
{
	bool	ok;
	double	y0 = m_A2;
	double  dy;
	double	sigma;

	ok = (m_A3 != 0.0);
	if (ok)
	{
		sigma = sqrt( 0.5/m_A3);
		dy = 8.0*sigma;
		*minradius = exp( y0 - dy);
		*maxradius = exp( y0 + dy);
	}
	else
	{
//		nxLog::Record(NXLOG_WARNING,"skRTParticleDist_LogNormal::GetQuadratureRadii, cannot do quadrature stuff until coefficients are set");
		*minradius = 0.0;
		*maxradius = 1.0;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::DeepCopy		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_LogNormal::Copy( const skRTParticleDist_LogNormal& other)
{
	m_QTEMP = other.m_QTEMP;
	m_A2    = other.m_A2;
	m_A3    = other.m_A3;
	return true;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_LogNormal::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_LogNormal::CreateClone(skRTParticleDist** userclone) const
{
	skRTParticleDist_LogNormal*	clone;
	bool						ok;

	clone = new skRTParticleDist_LogNormal;
	ok = (clone != NULL );
	if (ok)
	{
		clone->AddRef();
		ok = clone->Copy(*this);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING," skRTParticleDist_LogNormal::CreateClone, Error copying this object to the clone");
		}
	}
	*userclone = clone;
	return ok;

}


/*---------------------------------------------------------------------------
 *'					skRTParticleDist_PowerLaw::SetParameters		2003-10-9
!c-----------------------------------------------------------------------
!c THE FOLLOWING (NSD=5) IS THE POWER LAW DISTRIBUTION (2.61,
!c N(R)=CONST*R**(-A) FROM RMIN=B TO RMAX=C AND ZERO OTHERWISE
!c-----------------------------------------------------------------------
 *-------------------------------------------------------------------------*/

bool skRTParticleDist_PowerLaw::SetDistributionParameters( double A, double B, double C )
{
	bool	ok;
	double	A1;

	ok = (B > 0);
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "nxRTParticleDist_Powerlaw::SetParameters, The B (%g) parameter must be greater than zero ", (double)B);
		m_A     = 0.0;
		m_B     = 0.0;
		m_C     = 0.0;
		m_CONST = 0.0;
	}
	else
	{
		m_A = A;
		m_B = B;
		m_C = C;

		A1 = 1.0 - m_A;
		if ( fabs(A1) > 1.0E-15) m_CONST = A1 / ( pow(m_C,A1) - pow(m_B,A1));
		else                     m_CONST = 1.0 /log(m_C/m_B);
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					skRTParticleDist_PowerLaw::GetParameters		2009-6-24
 *-------------------------------------------------------------------------*/

void skRTParticleDist_PowerLaw::GetDistributionParameters( double* A, double* B, double* C )
{
	*A = m_A;
	*B = m_B;
	*C = m_C;
}


/*-----------------------------------------------------------------------------
 *					skRTParticleDist_PowerLaw::GetDistributionParameterArray		2013-6-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_PowerLaw::GetDistributionParameterArray( double* parameters, size_t maxparams, size_t* numparams ) const
{
	bool	ok;
	ok = (maxparams >= 3);
	if (ok)
	{
		parameters[0] = m_A;
		parameters[1] = m_B;
		parameters[2] = m_C;
		*numparams    = 3;
	}
	else
	{
		*numparams = 0;
		nxLog::Record(NXLOG_WARNING,"skRTParticleDist_PowerLaw::GetDistributionParameterArray, Error returning parameters as user buffer is smaller than 3 elements. Thats a problem");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_PowerLaw::GetDistributionParameterSpeciesID		2013-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_PowerLaw::GetDistributionParameterSpeciesID	( CLIMATOLOGY_HANDLE* /*paramids*/, size_t /*maxparams*/, size_t* numparams) const
{
	nxLog::Record(NXLOG_WARNING,"skRTParticleDist_PowerLaw::GetDistributionParameterSpeciesID, not yet implemeted");
	*numparams = 0;
	return false;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_PowerLaw::CachingDescriptor		2008-4-16*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skRTParticleDist_PowerLaw::CachingDescriptor() const
{
	nxString descriptor;

	descriptor.sprintf("powerlaw_%12g_%12g_%12g",(double)m_A, (double)m_B, (double)m_C);
	return descriptor;
}
/*---------------------------------------------------------------------------
 *'					skRTParticleDist_PowerLaw::Distribution		2003-10-9
 *-------------------------------------------------------------------------*/

double skRTParticleDist_PowerLaw::Distribution( double R )
{
	double	sizedis;

	if ((R < m_B) || (R > m_C))
	{
		sizedis = 0.0;
	}
	else
	{
		sizedis = m_CONST/( pow(R,m_A) );
		if (sizedis < 1.E-30) sizedis = 0.0;
	}
	return sizedis;
}
/*-----------------------------------------------------------------------------
 *					skRTParticleDist_PowerLaw::GetQuadratureRadii		2008-5-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_PowerLaw::GetQuadratureRadii	(double* /*minradius*/, double* /*maxradius*/) const
{
	nxLog::Record( NXLOG_ERROR, "skRTParticleDist_PowerLaw::GetQuadratureRadii, Not yet implemented");
	NXASSERT(( false ));
	return false;
}

/*-----------------------------------------------------------------------------
 *					skRTParticleDist_PowerLaw::DeepCopy		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_PowerLaw::Copy( const skRTParticleDist_PowerLaw& other)
{
	m_CONST = other.m_CONST;
	m_A     = other.m_A;
	m_B     = other.m_B;
	m_C     = other.m_C;

	return true;
}
/*-----------------------------------------------------------------------------
 *					skRTParticleDist_PowerLaw::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTParticleDist_PowerLaw::CreateClone(skRTParticleDist** userclone) const
{
	skRTParticleDist_PowerLaw*	clone;
	bool						ok;

	clone = new skRTParticleDist_PowerLaw;
	ok = (clone != NULL );
	if (ok)
	{
		clone->AddRef();
		ok = clone->Copy(*this);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skRTParticleDist_PowerLaw::CreateClone,Error copying this object to the clone");
		}
	}
	*userclone = clone;
	return ok;

}

