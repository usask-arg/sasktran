#include <skopticalproperties21.h>
#include <omp.h>

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_RayleighDryAir::skOpticalProperties_RayleighDryAir		2003-10-23 */
/**	Construct the Rayleigh, Dry Air scattering object.   
 *-------------------------------------------------------------------------*/

skOpticalProperties_RayleighDryAir::skOpticalProperties_RayleighDryAir()
{
	m_O2mix   = 0.209476;
	m_N2mix   = 0.78084;
	m_CO2mix  = 0.000314;
	m_Armix   = 0.00934;
	m_XXmix   = 1.0 -( m_O2mix + m_N2mix + m_CO2mix + m_Armix);
}

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_RayleighDryAir::~skOpticalProperties_RayleighDryAir		2003-11-3
 *-------------------------------------------------------------------------*/

skOpticalProperties_RayleighDryAir::~skOpticalProperties_RayleighDryAir()
{
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_RayleighDryAir::Depol		2014-2-20*/
/** **/
/*---------------------------------------------------------------------------*/

//double skOpticalProperties_RayleighDryAir::Depol(size_t threadindex) const	
//{
//	double delta = m_threadstate[threadindex].m_delta;
//
//	return (1-delta)/(1+0.5*delta);
//}
//

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_RayleighDryAir::SetAtmosphericState		2009-11-10*/
/** The Rayleigh cross-section in this implementation has no dependency upon
 *	atmospheric state, in the sense that we assume the atmosphere has the same
 *	proportions of O2, N2, CO2 and Ar.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_RayleighDryAir::SetAtmosphericState( skClimatology* /*neutralatmosphere*/ )
{
	return true;
}


bool skOpticalProperties_RayleighDryAir::SetLocation( const GEODETIC_INSTANT& /*pt*/, bool* crosssectionschanged  )
{
	if (crosssectionschanged != NULL) *crosssectionschanged = false;
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_RayleighDryAir::InternalClimatology_UpdateCache		2014-2-20*/
/** The Rayleigh cross-section has no dependency on climatology parameters
 *	such as pressure and temperature. Thus this method does nothing.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_RayleighDryAir::InternalClimatology_UpdateCache( const GEODETIC_INSTANT& /*pt*/)
{
	return true;
}


bool skOpticalProperties_RayleighDryAir::CalculateCrossSectionsArray(const double *wavenumber, int numwavenumber,
                                                                     double *absxs, double *extxs, double *scattxs) {
    static const double PiCubed   = 31.006276680299820175476315067101;		// This is PI**3
    static const double lochsmidt = 1.0/2.686763E19;						// Particles per cm3 in Ideal gas at 0C, 1 atmosphere.

    #pragma omp parallel for schedule(guided, 1)
    for(int idx = 0; idx < numwavenumber; idx++) {
        double  lamda   = 1.0E4/wavenumber[idx];		// get the wavelength in microns
        double  sigma   = wavenumber[idx]*1.0E-4;		// Get the wavenumber in micron-1
        double  sigma2  = sigma*sigma;			// Get the terms for the
        double  sigma4  = sigma2*sigma2;		// polynomial expansions
        double  sigma6  = sigma2*sigma4;
        double  sigma8  = sigma4*sigma4;
        double	nr;
        double  O2xsect;
        double  N2xsect;
        double  Arxsect;
        double  CO2xsect;
        double  XXxsect;
        double	O2Fk;				// King Correction factor for O2
        double	N2Fk;				// King Correction factor for N2
        double	CO2Fk;				// King Correction factor for CO2
        double	ArFk;				// King Correction factor for Argon.
        double  XXFk;
        double	nO2;				// Refractivity of O2 at STP
        double	nN2;				// Refractivity of N2 at STP
        double	nAr;				// Refractivity of Argon at STP
        double	nCO2;				// Refractivity of CO2 at STP
        double	nXX;				// Refractivity of miscellaneous components at STP.

        // ---- Calculate Refractivities of individual components at STP
        // ---- formulae for refractivities are taken from Ref 1) Bates 1984
        // ---- and Ref 2)and Ref 3)
        //
        // ---- Note that more upto date values for the refractivity of air are
        // ---- available but not for the refractivity of the individual components.
        // ---- The biggest error in all of this is the lack of account for moist air.

        if      ( lamda  < 0.254 ) nN2 = 1.0E-8*( 6998.749 + 3233582.0/(144-sigma2));
        else if ( lamda  < 0.468 ) nN2 = 1.0E-8*( 5989.242 + 3363266.3/(144-sigma2));
        else                       nN2 = 1.0E-8*( 6855.200 + 3243157.0/(144-sigma2));

        if      ( lamda < 0.221 )  nO2 = 1.0E-8*( 23796.7  + 168988.4/(40.9-sigma2));
        else if ( lamda < 0.288 )  nO2 = 1.0E-8*( 22120.4  + 203187.6/(40.9-sigma2));
        else if ( lamda < 0.546 )  nO2 = 1.0E-8*( 20564.8  + 248089.9/(40.9-sigma2));
        else                       nO2 = 1.0E-8*( 21351.1  + 218567.0/(40.9-sigma2));

        nAr  = (5.547E-4*0.5)*(  1.0
                                 + 5.15E-3*sigma2
                                 + 4.19E-5*sigma4
                                 + 4.09E-7*sigma6
                                 + 4.32E-9*sigma8);

        nCO2 = 1.0E-5*( 1205.5*(5.79925/(166.175-sigma2))		// Refractivity of CO2 at STP
                        + 0.12005/(79.609-sigma2)
                        + 0.0053334/(56.3064-sigma2)
                        + 0.0043244/(46.0196-sigma2)
                        + 0.0001218145/(0.0584738 - sigma2)
        );
        nXX  = nAr;								// Refractivity of Miscellaneous constituents at STP (close)

        // ---- Calculate King Correction factors for different dry air constituents
        // ---- formulae for King Correction factors are taken from Ref 1) Bates 1984.

        O2Fk  = 1.096 + 1.385E-3*sigma2 + 1.448E-4*sigma4;
        N2Fk  = 1.034 + 3.17E-4*sigma2;
        CO2Fk = 1.15;
        ArFk  = 1.0;
        XXFk  = 1.0;

        O2xsect  = m_O2mix *nO2 *nO2 *O2Fk;			// (Unnormalized) fraction of cross-section from O2
        N2xsect  = m_N2mix *nN2 *nN2 *N2Fk;			// (Unnormalized) fraction of cross-section from N2
        Arxsect  = m_Armix *nAr *nAr *ArFk;			// (Unnormalized) fraction of cross-section from Ar
        CO2xsect = m_CO2mix*nCO2*nCO2*CO2Fk;			// (Unnormalized) fraction of cross-section from CO2
        XXxsect  = m_XXmix *nXX *nXX *XXFk;			// (Unnormalized) fraction of cross-section from Miscellaneous species
        nr       = O2xsect + N2xsect + Arxsect + CO2xsect + XXxsect;	// Total Unnormalized cross-section from everything

        extxs[idx] = 32.0*PiCubed/3.0*nr*nxmath::sqr(wavenumber[idx]*wavenumber[idx]*lochsmidt);	// Get the Rayleigh cross-section in cm2.
        scattxs[idx] = extxs[idx];
        absxs[idx] = 0.0;
    }
    return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_RayleighDryAir::CalculateCrossSections		2014-2-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_RayleighDryAir::CalculateCrossSections( double wavenum, double* absxs, double* extxs, double* scattxs, skOpticalProperties_RayleighDryAir::RayleighWavelength_TLS* state) const 
{
	static const double PiCubed   = 31.006276680299820175476315067101;		// This is PI**3
	static const double lochsmidt = 1.0/2.686763E19;						// Particles per cm3 in Ideal gas at 0C, 1 atmosphere.
	bool				ok;

	ok = (state->m_wavenumber == wavenum);
	if (!ok)
	{

		double  lamda   = 1.0E4/wavenum;		// get the wavelength in microns
		double  sigma   = wavenum*1.0E-4;		// Get the wavenumber in micron-1
		double  sigma2  = sigma*sigma;			// Get the terms for the 
		double  sigma4  = sigma2*sigma2;		// polynomial expansions
		double  sigma6  = sigma2*sigma4;
		double  sigma8  = sigma4*sigma4;
		double	nr;
		double  xsection;
		double  O2xsect;
		double  N2xsect;
		double  Arxsect;
		double  CO2xsect;
		double  XXxsect;
		double	O2Fk;				// King Correction factor for O2
		double	N2Fk;				// King Correction factor for N2
		double	CO2Fk;				// King Correction factor for CO2
		double	ArFk;				// King Correction factor for Argon.
		double  XXFk;
		double	nO2;				// Refractivity of O2 at STP
		double	nN2;				// Refractivity of N2 at STP
		double	nAr;				// Refractivity of Argon at STP
		double	nCO2;				// Refractivity of CO2 at STP
		double	nXX;				// Refractivity of miscellaneous components at STP.

		// ---- Calculate Refractivities of individual components at STP
		// ---- formulae for refractivities are taken from Ref 1) Bates 1984
		// ---- and Ref 2)and Ref 3)
		//
		// ---- Note that more upto date values for the refractivity of air are
		// ---- available but not for the refractivity of the individual components.
		// ---- The biggest error in all of this is the lack of account for moist air.

		if      ( lamda  < 0.254 ) nN2 = 1.0E-8*( 6998.749 + 3233582.0/(144-sigma2));
		else if ( lamda  < 0.468 ) nN2 = 1.0E-8*( 5989.242 + 3363266.3/(144-sigma2));
		else                       nN2 = 1.0E-8*( 6855.200 + 3243157.0/(144-sigma2));

		if      ( lamda < 0.221 )  nO2 = 1.0E-8*( 23796.7  + 168988.4/(40.9-sigma2));
		else if ( lamda < 0.288 )  nO2 = 1.0E-8*( 22120.4  + 203187.6/(40.9-sigma2));
		else if ( lamda < 0.546 )  nO2 = 1.0E-8*( 20564.8  + 248089.9/(40.9-sigma2));
		else                       nO2 = 1.0E-8*( 21351.1  + 218567.0/(40.9-sigma2));

		nAr  = (5.547E-4*0.5)*(  1.0  
								+ 5.15E-3*sigma2
								+ 4.19E-5*sigma4 
								+ 4.09E-7*sigma6 
								+ 4.32E-9*sigma8);

		nCO2 = 1.0E-5*( 1205.5*(5.79925/(166.175-sigma2))		// Refractivity of CO2 at STP
						+ 0.12005/(79.609-sigma2)
						+ 0.0053334/(56.3064-sigma2)
						+ 0.0043244/(46.0196-sigma2)
						+ 0.0001218145/(0.0584738 - sigma2)
						);                          
		nXX  = nAr;								// Refractivity of Miscellaneous constituents at STP (close)
	    
		// ---- Calculate King Correction factors for different dry air constituents
		// ---- formulae for King Correction factors are taken from Ref 1) Bates 1984.

		O2Fk  = 1.096 + 1.385E-3*sigma2 + 1.448E-4*sigma4;
		N2Fk  = 1.034 + 3.17E-4*sigma2;
		CO2Fk = 1.15;
		ArFk  = 1.0;
		XXFk  = 1.0;

		O2xsect  = m_O2mix *nO2 *nO2 *O2Fk;			// (Unnormalized) fraction of cross-section from O2
		N2xsect  = m_N2mix *nN2 *nN2 *N2Fk;			// (Unnormalized) fraction of cross-section from N2
		Arxsect  = m_Armix *nAr *nAr *ArFk;			// (Unnormalized) fraction of cross-section from Ar
		CO2xsect = m_CO2mix*nCO2*nCO2*CO2Fk;			// (Unnormalized) fraction of cross-section from CO2
		XXxsect  = m_XXmix *nXX *nXX *XXFk;			// (Unnormalized) fraction of cross-section from Miscellaneous species
		nr       = O2xsect + N2xsect + Arxsect + CO2xsect + XXxsect;	// Total Unnormalized cross-section from everything

		xsection = 32.0*PiCubed/3.0*nr*nxmath::sqr(wavenum*wavenum*lochsmidt);	// Get the Rayleigh cross-section in cm2.

		state->m_delta	    = 0;
		state->m_deltaprime = 0;
		AddDepolarization( O2Fk,  O2xsect,  state  );
		AddDepolarization( N2Fk,  N2xsect,  state  );
		AddDepolarization( ArFk,  Arxsect,  state  );
		AddDepolarization( CO2Fk, CO2xsect, state );
		AddDepolarization( XXFk,  XXxsect,  state  );
		state->m_delta      /= nr;
		state->m_deltaprime /= nr;
		state->m_xsection    = xsection;
		state->m_wavenumber  = wavenum;
	}
	*absxs = 0.0;
	*extxs   = state->m_xsection;
	*scattxs = state->m_xsection;
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_RayleighDryAir::LookupUpThreadData		 2014- 10- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_RayleighDryAir::LookupUpThreadData(RayleighWavelength_TLS** data )
{
	size_t					threadnum;
	iterator				iter;
	bool					ok;
	static std::mutex		lock;							// Add a lock to make sure map.insert is thread safe.

	if (omp_in_parallel())
	{
		lock.lock();
	}

	threadnum = nxWorkerThreadManager::GetCurrentThreadIdCode();
	iter      = m_threadstate.find(threadnum);
	ok = (iter != m_threadstate.end() );
	if (!ok)
	{
		std::pair<iterator, bool>							result;
		std::pair<size_t, RayleighWavelength_TLS >	newentry(threadnum, RayleighWavelength_TLS() );

		result = m_threadstate.insert(  newentry);
		ok = result.second;
		iter = result.first;
		(*iter).second.m_wavenumber = 0.0;
		(*iter).second.m_xsection   = 0.0;
		(*iter).second.m_delta      = 0.0;
		(*iter).second.m_deltaprime = 0.0;
	}
	if (omp_in_parallel())
	{
		lock.unlock();
	}
	if (!ok)
	{
		*data = nullptr;
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_RayleighDryAir::LookupUpThreadData, error looking/creating thread data for thread %d", (int)threadnum);
	}
	else
	{
		*data = &(*iter).second;
	}
	return ok;

}

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_RayleighDryAir::CalculateCrossSections		2003-10-30 */
/**
 *	Calculates the Rayleigh cross-sections of Dry Air using the formulae
 *	given in Bates 1984. Note that I have compared this calculation to Bate's Rayleigh
 *	cross section calculations in Table 1 of his 1984 paper and found exact agreement (to
 *  the 4 significant digits he published).
 *
 *	The only difference between my calculation and Bates is that I make an account
 *  for the tiny fraction of gas that is not N2, O2, Ar or CO2.  I assume the small
 * 	residual gas is similar in properties to Argon.  Bates does not describe what
 *	he does with the small residual (it is implied that it is ignored in his paper).
 *
 *	This class does not pre-calculate any phase matrix tables as its efficient to
 *	generate them on the fly.
 *
 */
/*-------------------------------------------------------------------------*/

bool skOpticalProperties_RayleighDryAir::CalculateCrossSections( double wavenum, double* absxs, double* extxs, double* scattxs )
{
	bool	ok;
	RayleighWavelength_TLS* threaddata;

	ok =       LookupUpThreadData( &threaddata );
	ok = ok && CalculateCrossSections( wavenum, absxs, extxs, scattxs, threaddata);
	return ok;
}


/*---------------------------------------------------------------------------
 *'					skOpticalProperties_RayleighDryAir::AddDepolarization		2003-11-18 */
/**	Adds the depolarization term from this fraction of species in the
 *	atmosphere to the overall depolarization of "air".
 */
/*--------------------------------------------------------------------------*/

void skOpticalProperties_RayleighDryAir::AddDepolarization( double Fk, double fraction, skOpticalProperties_RayleighDryAir::RayleighWavelength_TLS* state) const
{
	double depol;

	depol                = 6.0*(Fk-1.0)/(3.0+ 7*Fk);					// Get the depolarization term from the King correction factor
	state->m_delta      += fraction*(1 - depol  )/(1+depol*0.5);		// Get the DELTA term eqn 2.16 Hansen and Travis
	state->m_deltaprime += fraction*(1 - 2*depol)/(1-depol);			// Get the DELTA' term eqn 2.16 Hansen and Travis 
}

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_RayleighDryAir::Interpolate_PhaseMatrixTables		2003-10-30 */
/**	Get the Phase matrix for the scattering angle = cos-1(mu). The phase
 *	matrix when integrated over 4 pi solid angle normalizes to 4Pi. This calculation
 *	is performed on the fly as it is relatively quick.
 *
 *  Note: equations for P33 and P44 were changed to match equation 2.15 in Hansen
 *  and Travis, Light Scattering in Planetary Atmosphere's, 1975.  I believe there
 *  were some typos - Tony Bathgate 2009/05/26
 */
/*-------------------------------------------------------------------------*/

bool skOpticalProperties_RayleighDryAir::Interpolate_PhaseMatrixTables( double mu, skRTPhaseMatrix* P, skOpticalProperties_RayleighDryAir::RayleighWavelength_TLS* state) const
{
	double cosalpha    = mu;
	double cosalpha2   = mu*mu;
	double sinalpha2   = 1 - cosalpha2;
	double cosalpha2p1 = 1 + cosalpha2;

	P->At(1,1) = (SKRTFLOAT)(  state->m_delta*0.75*cosalpha2p1);
	P->At(1,2) = (SKRTFLOAT)( -state->m_delta*0.75*sinalpha2 );
	P->At(1,3) = (SKRTFLOAT)( 0 );
	P->At(1,4) = (SKRTFLOAT)( 0 );
	
	P->At(2,1) = (SKRTFLOAT)( P->At(1,2));
	P->At(2,2) = (SKRTFLOAT)( P->At(1,1));
	P->At(2,3) = (SKRTFLOAT)( 0);
	P->At(2,4) = (SKRTFLOAT)( 0);

	P->At(3,1) = (SKRTFLOAT)( 0);
	P->At(3,2) = (SKRTFLOAT)( 0);
	P->At(3,3) = (SKRTFLOAT)( state->m_delta*1.5*cosalpha);
	P->At(3,4) = (SKRTFLOAT)( 0);

	P->At(4,1) = (SKRTFLOAT)( 0);
	P->At(4,2) = (SKRTFLOAT)( 0);
	P->At(4,3) = (SKRTFLOAT)( 0);
	P->At(4,4) = (SKRTFLOAT)( state->m_deltaprime*state->m_delta*1.5*cosalpha);

	P->At(1,1) += (SKRTFLOAT)( (1-state->m_delta) );
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_RayleighDryAir::CalculatePhaseMatrix		2008-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_RayleighDryAir::CalculatePhaseMatrix( double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix)
{
	bool	ok;
	double  absxs, extxs, scattxs ;
	RayleighWavelength_TLS* threaddata;

	ok =       LookupUpThreadData( &threaddata );
	ok = ok && CalculateCrossSections(wavenumber, &absxs, &extxs, &scattxs, threaddata);
	ok = ok && Interpolate_PhaseMatrixTables( cosscatterangle, phasematrix, threaddata);
	return ok;
}

bool skOpticalProperties_RayleighDryAir::CalculateP11(double wavenumber, std::pair<double, size_t> cosscatterandindex, double& p11)
{
	bool	ok;
	double  absxs, extxs, scattxs;
	RayleighWavelength_TLS* threaddata;

	ok = LookupUpThreadData(&threaddata);
	ok = ok && CalculateCrossSections(wavenumber, &absxs, &extxs, &scattxs, threaddata);

	double cosalpha2 = cosscatterandindex.first * cosscatterandindex.first;
	double cosalpha2p1 = 1 + cosalpha2;

	p11 = (SKRTFLOAT)(threaddata->m_delta*0.75*cosalpha2p1);

	p11 += (SKRTFLOAT)( (1-threaddata->m_delta) );

	return ok;
}


bool skOpticalProperties_RayleighDryAir::LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff) {
    double  absxs, extxs, scattxs;
	RayleighWavelength_TLS* threaddata;
	bool ok = LookupUpThreadData(&threaddata);

    ok = ok && CalculateCrossSections(wavenumber, &absxs, &extxs, &scattxs, threaddata);

    // Rayleigh scattering has 3 legendre moments, (1, 0, (1-delta)/(2+delta))
	opticalmaxcoeff = 3;

	coeff[0] = 1.0;
	coeff[1] = 0.0;
	coeff[2] = threaddata->m_delta / 2;

	return ok;
}

bool skOpticalProperties_RayleighDryAir::LegendreCoefficientsPolarized(double wavenumber, double *a1, double *a2, double *a3,
                                                                       double *a4, double *b1, double *b2, int usermaxcoeff,
                                                                       int &opticalmaxcoeff) {
    double  absxs, extxs, scattxs;
    RayleighWavelength_TLS* threaddata;
    bool ok = LookupUpThreadData(&threaddata);

    ok = ok && CalculateCrossSections(wavenumber, &absxs, &extxs, &scattxs, threaddata);

    opticalmaxcoeff = 3;
    // (1-delta) / (2+delta) = x / 2
    // 1 - delta = (2 + delta) * x / 2
    // 2 - 2*delta = 2x + delta x
    // 2 - 2 x = delta x + 2*delta
    // 2(1-x) = delta(2+x)
    // delta = 2(1-x) / (2+x)

    double depol = 2 * (1 - threaddata->m_delta) / (2 + threaddata->m_delta);

    a1[0] = 1.0;
    a1[1] = 0.0;
    a1[2] = (1 - depol) / (2 + depol);

    a2[0] = 0.0;
    a2[1] = 0.0;
    a2[2] = 6 * (1 - depol) / (2 + depol);

    a3[0] = 0.0;
    a3[1] = 0.0;
    a3[2] = 0.0;

    a4[0] = 0.0;
    a4[1] = 3 * (1 - 2*depol) / (2 + depol);
    a4[2] = 0.0;

    b1[0] = 0.0;
    b1[1] = 0.0;
    b1[2] = sqrt(6.0) * (1- depol) / (2+ depol);

    b2[0] = 0.0;
    b2[1] = 0.0;
    b2[2] = 0.0;

    return ok;
}
