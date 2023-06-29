#include <skopticalproperties21.h>
#include <limits.h>
#include <float.h>
#include <boost/math/quadrature/tanh_sinh.hpp>
#include <boost/math/special_functions/legendre.hpp>
/*---------------------------------------------------------------------------
 *		skOpticalProperties::skOpticalProperties		2003-11-25 */
/**	Initialize the base class. Initialize all numeric elements to zero and
  *	flag the class as dirty.
*/
/*-------------------------------------------------------------------------*/

skOpticalProperties::skOpticalProperties()
{
//	m_wavenumber      = 0.0;
//	m_kextperparticle = 0.0;
//	m_kabsperparticle = 0.0;
//	m_kscaperparticle = 0.0;
//	m_deltafunctionforwardscatterfraction = 0.0;		// By default there is no delta function with a forward scatter
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties::DeepCopy		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties::DeepCopy( const skOpticalProperties& /*other*/ )
{
//	m_wavenumber      = other.m_wavenumber;
//	m_kextperparticle = other.m_kextperparticle;
//	m_kabsperparticle = other.m_kabsperparticle;
//	m_kscaperparticle = other.m_kscaperparticle;
//	m_deltafunctionforwardscatterfraction = other.m_deltafunctionforwardscatterfraction;	// By default there is no delta function with a forward scatter
	return true;
}


/*---------------------------------------------------------------------------
 *					sk_ExtinctionPhaseMatrix::SetCrossSections		2003-10-22 */
/**	Set the extinction, ansorption and scattering cross-sections.
*/
/*-------------------------------------------------------------------------*/

/*
void skOpticalProperties::SetCrossSections( double wavenum, double kext, double kabs, double ksca)
{
	m_wavenumber      = wavenum;
	m_kextperparticle = kext;
	m_kabsperparticle = kabs;
	m_kscaperparticle = ksca;
}
*/


/*-----------------------------------------------------------------------------
 *					skOpticalProperties::CalculatePhaseMatrix		2008-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties::CalculatePhaseMatrix( double /*wavenumber*/ , double /*cosscatterangle*/, skRTPhaseMatrix* phasematrix)
{
	phasematrix->SetTo(0.0);
	phasematrix->At(1,1) = (SKRTFLOAT)1.0;
	return true;
}

bool skOpticalProperties::CalculateP11(double wavenumber, std::pair<double, size_t> cosscatterandindex, double& p11)
{
	skRTPhaseMatrix phase;
	bool ok = CalculatePhaseMatrix(wavenumber, cosscatterandindex.first, &phase);

	p11 = phase.At(1, 1);

	return ok;
}


/*----------------------------------------------------------------------------
 *'					sk_ExtinctionPhaseMatrix::CheckCosineRange		2003-12-11 */
/**	Quick check to get rid of nasty round off issues in cos(zenith angle).	*/
/*--------------------------------------------------------------------------*/

void skOpticalProperties::CheckCosineRange( double *mu)
{
	if (*mu < -1.0) *mu = -1.0;
	if (*mu >  1.0) *mu = 1.0;
}

/*---------------------------------------------------------------------------
 *					sk_ExtinctionPhaseMatrix::Get_RotatedPhaseMatrix		2003-10-22*/
/**	Get the phase matrix when the light is scattered from incoming zenith angle 
 *	#muprime to outgoing zenith angle #mu and the change in azimuth
 *	is #dphi radians (outgoing azimuth - incoming azimuth).
 *	This applies all of the rotation matrices to ensure the Stokes vectors refer to the plane
 *	formed by the ray and the vertical. Dont ask about beams exactly in the vertical direction.
 *	I think the code fails.
 *  A description of what this code does can be found in Appendix A of 
 *  reference 1. 
 *
 *	@param mu
 *		cosine of the outgoing zenith angle
 *
 *	@param muprime
 *		cosine of the incoming zenith angle
 *
 *	@param dphi
 *		The change in azimuth in radians: outgoing azimuth minus incoming azimuth in radians
 *
 *	@param rotatedmatrix
 *		Returns the phase matrix rproperly rotated so the reference polarization is
 *		in the zenith = 0 direction.t
 *
 *	\par Reference
 *	-# \b McLinden \b et \b al. A vector radiative-transfer model for the Odin/Osiris project.
 *	\e Can. \e J. \e Phys., \b 90: 375-393, 2002. \c doi:10.1139/P01-156
 *	.
**/
/*-------------------------------------------------------------------------*/

bool skOpticalProperties::GetRotatedPhaseMatrix( double wavenum, double mu, double muprime, double dphi, skRTPhaseMatrix* rotatedmatrix)
{
	double				costheta;
	skRTPhaseMatrix		phasematrix;
	bool				ok;

	CheckCosineRange(&mu);
	CheckCosineRange(&muprime);

	costheta =  phasematrix.GetScatteringAngle( mu, muprime, dphi );						// Get the scattering angle
	ok       = CalculatePhaseMatrix( wavenum, costheta, &phasematrix);									// Get the phase matrix for this scattering angle
	if (!ok)																				// IF it didnt work
	{																								// Then we have a problem
		nxLog::Record(NXLOG_WARNING,"sk_ExtinctionPhaseMatrix::Get_RotatedPhaseMatrix, Error retrieving phasematrix for cos(scattering angle) = %g, Returning empty phase matrix ", (double)costheta);
		rotatedmatrix->SetTo(0.0);																	// So tell everyone
	}																								// and that is that
	else																							// otherwise
	{																								// the other nasty
		if ( ((1.0-fabs(mu)) < 10.0*DBL_EPSILON) || (( 1.0 - fabs(muprime)) < 10.0*DBL_EPSILON))	// is to make sure nothing is parallel to vertical
		{																							// as the Stokes vector reference is undefined
			nxLog::Record(NXLOG_WARNING, "sk_ExtinctionPhaseMatrix::Get_RotatedPhaseMatrix, The input (%g) or output zenith angle (%g) is in the vertical direction, this is not allowed. Returning unrotated phase matrix", (double)mu, (double)muprime);
			ok              = false;														// So flag the problem
			*rotatedmatrix  = phasematrix;													// and copy the matrix over in unrotated form
		}																					// and that is that
		else																				// everything is checked out
		{																					// so
			if (ok) phasematrix.ApplyStokesRotation( mu, muprime, dphi, rotatedmatrix );	// Do the stokes rotation
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties::CalculateCrossSectionsArray		2014-2-28*/
/** Calculate the cross-sections for an array of wavenumbers. The defualt
 *	implementation is a multi threaded omp loop to calculate the cross-section
 *	for each wavenumber separately.
 *	Optimized versions, for Hitran and Voigt lineshapes for example, override this
 *	member and provide optimized, multi-threaded versions.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties::CalculateCrossSectionsArray( const double* wavenumber, int numwave, double* absxs, double* extxs, double* scattxs)
{
	bool	ok = true;

	
	#pragma omp parallel for schedule(dynamic)
	for (int iw = 0; iw < numwave; iw++)
	{
		bool ok1;
		ok1 = CalculateCrossSections( wavenumber[iw], &absxs[iw], &extxs[iw], &scattxs[iw]);
		#pragma omp critical
		{
			ok = ok && ok1;
		}
	}
	return ok;
}

bool skOpticalProperties::LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff) {
	boost::math::quadrature::tanh_sinh<double> integrator;

	for (int l = 0; l < usermaxcoeff; ++l) {
		auto f = [this, wavenumber, l](double x) {
			skRTPhaseMatrix phase;
			bool ok = this->CalculatePhaseMatrix(wavenumber, x, &phase);

			return phase.At(1, 1) * boost::math::legendre_p(l, x);
		};

		coeff[l] = integrator.integrate(f, 1e-5);
		coeff[l] *= double(2 * l + 1) / 2.0;
	}
	// Normalize the result
	double norm = coeff[0];
	for (int l = 0; l < usermaxcoeff; ++l) {
		coeff[l] /= norm;
	}

	opticalmaxcoeff = usermaxcoeff;

	return true;
}

bool skOpticalProperties::LegendreCoefficientsPolarized(double wavenumber, double *a1, double *a2, double *a3,
                                                        double *a4, double *b1, double *b2, int usermaxcoeff,
                                                        int &opticalmaxcoeff) {
    nxLog::Record(NXLOG_WARNING, "Requesting LegendreCoefficientsPolarized on a species that has not implemented them");
    return false;
}