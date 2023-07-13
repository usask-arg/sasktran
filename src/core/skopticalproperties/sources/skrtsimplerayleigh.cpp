#include <skopticalproperties21.h>


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_SimpleRayleigh::skOpticalProperties_SimpleRayleigh		2008-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_SimpleRayleigh::skOpticalProperties_SimpleRayleigh()
{
	static bool firsttime = true;
	if (firsttime)
	{
		nxLog::Record(NXLOG_INFO,"skOpticalProperties_SimpleRayleigh::skOpticalProperties_SimpleRayleigh, Be warned that you are using the deprecated Rayleigh scattering algorithms.");
		firsttime = false;
	}
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_SimpleRayleigh::CalculateCrossSections		2008-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_SimpleRayleigh::CalculateCrossSections(double wavenumber, double* absxs, double* extxs, double* scattxs)
{
	double Wvlen;

	Wvlen = 1.0E7/wavenumber/1e3;
	*scattxs  = (4.0e-28/(exp((3.916+0.074*Wvlen+0.05/Wvlen)*log(Wvlen))));
	*extxs    = *scattxs;
	*absxs    = 0.0;
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_SimpleRayleigh::CalculatePhaseMatrix		2008-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_SimpleRayleigh::CalculatePhaseMatrix( double /*wavenumber*/ , double cosscatterangle, skRTPhaseMatrix* P)
{
	double mu = cosscatterangle;
	double cosalpha = mu;
	double cosalpha2 = mu * mu;
	double sinalpha2 = 1 - cosalpha2;
	double cosalpha2p1 = 1 + cosalpha2;


	P->At(1, 1) = (SKRTFLOAT)(0.75*cosalpha2p1);
	P->At(1, 2) = (SKRTFLOAT)(-0.75*sinalpha2);
	P->At(1, 3) = (SKRTFLOAT)(0);
	P->At(1, 4) = (SKRTFLOAT)(0);

	P->At(2, 1) = (SKRTFLOAT)(P->At(1, 2));
	P->At(2, 2) = (SKRTFLOAT)(P->At(1, 1));
	P->At(2, 3) = (SKRTFLOAT)(0);
	P->At(2, 4) = (SKRTFLOAT)(0);

	P->At(3, 1) = (SKRTFLOAT)(0);
	P->At(3, 2) = (SKRTFLOAT)(0);
	P->At(3, 3) = (SKRTFLOAT)(1.5*cosalpha);
	P->At(3, 4) = (SKRTFLOAT)(0);

	P->At(4, 1) = (SKRTFLOAT)(0);
	P->At(4, 2) = (SKRTFLOAT)(0);
	P->At(4, 3) = (SKRTFLOAT)(0);
	P->At(4, 4) = (SKRTFLOAT)(1.5*cosalpha);

	return true;

}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_SimpleRayleigh::CreateClone		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_SimpleRayleigh::CreateClone( skOpticalProperties** clone) const
{
	*clone = new skOpticalProperties_SimpleRayleigh();
	(*clone)->AddRef();
	return true;
}
*/

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_SimpleRayleigh::SetAtmosphericState		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_SimpleRayleigh::SetAtmosphericState( skClimatology* /*neutralatmosphere*/ )
{
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_SimpleRayleigh::SetLocation		 2015- 11- 17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_SimpleRayleigh::SetLocation( const GEODETIC_INSTANT& /*pt*/, bool* crosssectionschanged  )
{
	if (crosssectionschanged != NULL) *crosssectionschanged = false;
	return true;
}


bool skOpticalProperties_SimpleRayleigh::LegendreCoefficientsPolarized(double wavenumber, double *a1, double *a2, double *a3,
                                                                       double *a4, double *b1, double *b2, int usermaxcoeff,
                                                                       int &opticalmaxcoeff) {

    opticalmaxcoeff = 3;


    double depol = 0;

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

    return true;
}

bool skOpticalProperties_SimpleRayleigh::LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff) {
    double  absxs, extxs, scattxs;

    // Rayleigh scattering has 3 legendre moments, (1, 0, (1-delta)/(2+delta))
    opticalmaxcoeff = 3;

    coeff[0] = 1.0;
    coeff[1] = 0.0;
    coeff[2] = 0.5;

    return true;
}