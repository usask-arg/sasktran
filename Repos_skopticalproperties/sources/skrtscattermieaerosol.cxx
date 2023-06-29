#include <skopticalproperties21.h>
#include <nxbase_threads.h>

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_MieAerosol::skOpticalProperties_MieAerosol		2003-10-9*/
/**	Default constructor initializes the Mie Aerosol code so it uses
 *	100 fourier harmonics.
**/
/*-------------------------------------------------------------------------*/

skOpticalProperties_MieAerosol::skOpticalProperties_MieAerosol()
{
	init();
}

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_MieAerosol::~nxRTMieAerosolScatter		2003-10-9 */
/**	default destructor														*/
/*--------------------------------------------------------------------------*/

skOpticalProperties_MieAerosol::~skOpticalProperties_MieAerosol()
{
	ReleaseDistribution();
	ReleaseRI();
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosol::operator=		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/


/*
bool skOpticalProperties_MieAerosol::DeepCopy( const skOpticalProperties_MieAerosol& other )
{
	bool	ok1;
	bool	ok2;
	bool	ok3 = true;
	bool	ok4 = true;
	bool	ok5;
	bool	ok6;

	ReleaseDistribution();
	ReleaseRI();

	ok1 = skOpticalProperties::DeepCopy( other );
	ok2 = m_mie.DeepCopy( other.m_mie );
	if (other.m_distribution    != NULL)  ok3 = other.m_distribution->CreateClone( &m_distribution );
	if (other.m_refractiveindex != NULL)  ok4 = other.m_refractiveindex->CreateClone( &m_refractiveindex);
	ok5 = m_quadrature.DeepCopy( other.m_quadrature );
	ok6 = m_Pij.DeepCopy( other.m_Pij );
	m_isdirty = true;
	return (ok1 && ok2 && ok3 && ok4 && ok5 && ok6);
}
*/

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosol::CreateClone		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_MieAerosol::CreateClone(skOpticalProperties** userclone) const
{
	skOpticalProperties_MieAerosol*	clone;
	bool						ok;

	clone = new skOpticalProperties_MieAerosol;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_MieAerosol::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		ok = clone->DeepCopy( *this );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_MieAerosol::CreateClone, Error creating copying this object to clone object, this might cause significant issues");
		}
	}
	*userclone = clone;
	return ok;
}
*/

/*----------------------------------------------------------------------------
 *'					skOpticalProperties_MieAerosol::ReleaseRI		2003-10-21		*/
/**	Release the current refractive index object								*/
/*--------------------------------------------------------------------------*/

void skOpticalProperties_MieAerosol::ReleaseRI()
{
	if (m_refractiveindex != NULL)
	{
		m_refractiveindex->Release();
		m_refractiveindex = NULL;
	}
}
/*---------------------------------------------------------------------------
 *'					skOpticalProperties_MieAerosol::Set_RefractiveIndex	2003-10-15 */
/**	Set the refractive index of the Mie Aerosol spheres to the new value.
 *
 *	@param newri
 *	Pointer to the new refractive index object.  The class will call AddRef()
 *	on the new refractive index to keep the object alive and it will call
 *	Release() when it is finished with the object.
 */
/*--------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosol::Set_RefractiveIndex( skRTRefractiveIndex* newri)
{
	bool	ok;

	ok = (newri == m_refractiveindex);
	if (!ok)
	{
		if (newri != NULL) newri->AddRef();
		ReleaseRI();
		m_refractiveindex  = newri;
		ok = (m_refractiveindex != NULL);
	}
	SetDirty();
	return ok;
}

/*---------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosol::ReleaseDistribution		2003-10-9 */
/**	Release the curent paticle distribution object.							*/
/*--------------------------------------------------------------------------*/

void skOpticalProperties_MieAerosol::ReleaseDistribution()
{
	if (m_distribution != NULL)
	{
		m_distribution->Release();
		m_distribution = NULL;
	}
}

/*---------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosol::Set_ParticleDistribution		2003-10-15 */
/**	Set #newdistribution as the new particle distribution for the Mie
 *	scatterers.  The #newdistribution is cloned into memory. The code checks to see
 *	if the current and new distributions are identical
 *
 *	@param	newdistribution
 *		The new particle distribution.
 */
/*--------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosol::Set_ParticleDistribution( skRTParticleDist* newdistribution )
{
	bool	ok;
	bool	issame;

	issame =	(newdistribution != NULL)
		     && (m_distribution  != NULL)
			 && (m_distribution->IsSameDistributionAs(newdistribution));
	ok = issame;
	if (!ok)
	{
		ReleaseDistribution();
		ok = (newdistribution == NULL);
		if (!ok) ok  = newdistribution->CreateClone(&m_distribution);
		SetDirty();
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *'					skOpticalProperties_MieAerosol::LoadDefaultStratosphere		2003-10-9 */
/**	Load up a default stratospheric aerosol. Normally used for testing purposes.
 *	NOT RECOMMENDED FOR GENERAL USE.
 */
/*-------------------------------------------------------------------------*/

void skOpticalProperties_MieAerosol::LoadDefaultStratosphere()
{
	skRTParticleDist_LogNormal* dist = new skRTParticleDist_LogNormal;
	skRTRefractiveIndex_H2SO4*  ri   = new skRTRefractiveIndex_H2SO4;

	dist->SetDistributionParameters( 0.17, 1.552707, 0.0);		// Stratosphere default Log Normal parameters. Mode radius, Mode width
//	dist->SetRadiusRanges	 ( 0.001, 3.0 );					// effective (zero) to infinity ranges for integration
	Set_ParticleDistribution (dist);
	Set_RefractiveIndex      (ri);
}

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_MieAerosol::LoadDefaultTroposphere		2003-10-9 */
/**	Load up a default stratospheric aerosol. Normally used for testing purposes.
 *	NOT RECOMMENDED FOR GENERAL USE.
 */
/*--------------------------------------------------------------------------*/

void skOpticalProperties_MieAerosol::LoadDefaultTroposphere()
{
	skRTParticleDist_LogNormal* dist = new skRTParticleDist_LogNormal;
	skRTRefractiveIndex_Water*  ri   = new skRTRefractiveIndex_Water;

	dist->SetDistributionParameters	( 50.0, 150.0, 0.0);			// Stratosphere default Log Normal parameters.
//	dist->SetRadiusRanges		( 0.01, 50.0 );						// effective (zero) to infinity ranges for integration
	Set_ParticleDistribution	(dist);
	Set_RefractiveIndex			(ri);
}

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_MieAerosol::init		2003-10-9 */
/**	Inititalize the Mie aerosol class. By default we use 100 points in our
 *	radius distribution quadrature and we use 0.1 degree resolution for
 *	0 to 5 degrees and 175-180 scattering angles (forward scatter) and 1
 *	degree resolution for all others.
 **/
/*--------------------------------------------------------------------------*/

void skOpticalProperties_MieAerosol::init()
{
	double	startangle[2] = {0.0,  5.0};			// do 0 t0 5  in steps of 0.1 degrees
	double  endangle[2]   = {5.0, 90.0};			// do 5 to 90 in steps of 1   degrees
	double  resolution[2] = {0.1,  1.0};			// And the angles are symmetric about 90 degrees
	m_numlegendre = 64;

	m_quadrature.SetOrder(1000);						// Default is to use 1000 points when doing integration over particle size.
	m_distribution    = NULL;
	m_refractiveindex = NULL;
	m_mie.Set_SymmetricScatterAngles( startangle, endangle, resolution, 2 );

	m_mie.Set_MaxLegendreMoment(m_numlegendre);
	SetDirty();
	m_internalwavenumber = -1;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosol::CalculateCrossSections		 2014- 10- 24*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosol::CalculateCrossSections( double wavenum, double* absxs, double* extxs, double* scattxs)
{
	bool			ok;

	m_mutex_xsectionslock.lock();
	ok = CalculateCrossSectionsInternal( wavenum, absxs, extxs, scattxs);
	m_mutex_xsectionslock.unlock();
	return ok;
}

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_MieAerosol::CalculateCrossSections		2003-10-15 */
/**	This code is a protected virtual interface call that sk_ExtinctionPhaseMatrix
 *	will call to generate new extinction, absorption and scattering cross sections. 
 *	This call must calculate the new cross-sections in cm**2 and call SetCrossSections
 *	to pass the information back to the base class. In addition the code should use this opportunity to update any internal tables
 *	that enable rapid evaluation of the phase matrix for arbitrary scattering angle.
 *
 *	This code actually generates a lookup table of the phase matrix for a range of angles.
 *	The phase matrix has been integrated over the particle distribution.  The user may
 *	derive the phase matrix at any angle by linear interpolation. I have implicitly
 *	assumed that the internal scattering angle tables are monotonically
 *	increasing and the user has enough common sense to place a sufficient number of points
 *	in the angular distribution. The integration formulae are taken from eqns 2.47 to
 *	2.49 of Hansen and Travis.
 **/
/*-------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosol::CalculateCrossSectionsInternal( double wavenum, double* absxs, double* extxs, double* scattxs)
{
	bool					ok;
	double					minradius;
	double					maxradius;
	size_t					numangles;
	double					nr;
	double					nsum;
	nx2dArray<double>		fij, bl;
	double					c;
	double					factor;
	size_t					i;
	double					ksca = 0.0;
	double					kext = 0.0;
	double					kabs = 0.0;
	double					lamdamicron = 1.0E4/wavenum;		// Get the wavelength in microns
	std::complex<double>	ri;

	if (!m_isdirty && (abs(m_internalwavenumber - wavenum) / wavenum) < 1e-12)
	{
		// Do not need to repeat the calculation
		return true;
	}

	ok = ( m_distribution != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MieAerosol, Cannot calculate the Mie cross sections as no particle distribution is defined (PLEASE DEFINE ONE with Set_ParticleDistribution etc.)");
	}
	else
	{
		ok = (m_refractiveindex != NULL);										// Make sure the refractive index object is not NULL
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MieAerosol, Cannot calculate the Mie cross sections as no refractive index is defined (PLEASE DEFINE ONE with Set_RefractiveIndex)");
		}
		else
		{
			ri  = m_refractiveindex->RefractiveIndex( wavenum );
			m_mie.Set_RefractiveIndex( ri.real(), ri.imag() );
			if (ri.imag() == 0.0) m_mie.Set_MimCut( ri.real() + 1.0);
			else                  m_mie.Set_MimCut( 0.0 );
			m_mie.Set_Wavelength( lamdamicron  );						// Set the MIe scattering wavelnth to microns (Compatible with particle distribution radius )

			ok =       m_distribution->GetQuadratureRadii( &minradius, &maxradius );
			if (!ok)
			{
				ok = (minradius == 0.0)  && (maxradius == 1.0);				// This condition occurs when moddewidth or radius == 0 and does not require a message
				if (!ok)													// but if it is not this condition
				{															// then we have something goofy going on.
					nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MieAerosol::CalculateCrossSections, There was an error selecting the range of radii for the quadrature");
				}
				ksca  = 0.0;												// Convert from square microns to square centimeters
				kext  = 0.0;												// Convert from square microns to square centimeters
				kabs  = 0.0;												// Convert from square microns to square centimeters
			}
			else
			{
				m_quadrature.SetRange( minradius, maxradius );				// Set the quadrature range of points
					
				const nx1dArray<double>&	R = m_quadrature.X();			// R are the quadrature points in the radius distribution
				const nx1dArray<double>&	W = m_quadrature.W();			// W are the desired weights   in the radius distribution quadrature

				numangles = m_mie.Get_NumAngles();							// Get the number of cosine angles in the Mie scattering calculation
				m_Pij.SetSize(4,numangles);									// The Fij scattering array calculated by m_mie.( 4, num_Cosangles, )
				m_Pij.SetTo(0);

				m_legendre.SetSize(4, m_numlegendre);						// The legendre moments calculated by m_mie
				m_legendre.SetTo(0);

				nsum    = 0.0;
				for (i = 0; i < R.size(); i++)								// For each radius in the radius quadrature
				{															// then
					m_mie.Set_Radius( R.At(i) );							// configure the Mie scattering radius
					m_mie.Get_ScatteringMatrixCoeffs( &fij );				// Get the Fij calculated by the mie code.
					m_mie.Get_LegendreCoefficients( &bl );
					nr   = W.At(i)*m_distribution->Distribution( R.At(i) );	// Get n(r).w for evaluating the quadrature
					nsum   += nr;
					fij    *= nr;											// Get Fij.n(r).w
					m_Pij  += fij;											// Fij = Sum( Fij.n(r).w )
					bl *= nr;
					m_legendre += bl;
					ksca += (nr*m_mie.Csca());								// ksca = Sum
					kext += (nr*m_mie.Cext());								// kext = Sum
				}
				ok = ok && ((nsum > 0.99) && (nsum < 1.01));
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"skOpticalProperties_MieAerosol::CalculateCrossSections, The quadrature of the particle distribution was not normalized to unity with good accuracy. we got (%e) instead of 1.0", (double)nsum);
					NXASSERT(( ok ));
				}
				c = 4.0*nxmath::Pi/nxmath::sqr(m_mie.Get_k());				// Get the normalization constant from the Mueller matrix
				for (i=0; i < numangles; i++)								// to the phase matrix.
				{															// and for each angle
					factor        = c/ksca;								// apply the normalization number
					m_Pij.At(0,i) *= factor;								// to each of the 4 scattering elements
					m_Pij.At(1,i) *= factor;	
					m_Pij.At(2,i) *= factor;
					m_Pij.At(3,i) *= factor;
				}
				// Apply the normalization to the legendre moments
				std::vector<double> norm{ m_legendre.At(0, 0), m_legendre.At(1, 0), m_legendre.At(2, 0), m_legendre.At(3, 0) };
				for (i = 0; i < m_numlegendre; i++)
				{
					double factor = 2.0*i + 1;

					m_legendre.At(0, i) *= factor / norm[0];
					m_legendre.At(1, i) *= factor / norm[1];
					m_legendre.At(2, i) *= factor / norm[2];
					m_legendre.At(3, i) *= factor / norm[3];
				}

				ksca *= 1.0E-8;												// Convert from square microns to square centimeters
				kext *= 1.0E-8;												// Convert from square microns to square centimeters
				kabs  = kext - ksca;										// Convert from square microns to square centimeters
			}
			*absxs   = kabs; 
			*extxs   = kext;
			*scattxs = ksca;
			m_isdirty = !ok;
			m_internalwavenumber = wavenum;
		}
	}
	if (!ok)
	{
		*absxs   = 0.0;
		*extxs   = 0.0;
		*scattxs = 0.0;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosol::SetAtmosphericState		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosol::SetAtmosphericState( skClimatology* /*neutralatmosphere*/)
{
	return true;
}
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosol::SetAtmosphericState		2009-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosol::SetLocation( const GEODETIC_INSTANT& /*pt*/, bool* crosssectionschanged  )
{
	if (crosssectionschanged != NULL) *crosssectionschanged = false;
	return true;
}

/*---------------------------------------------------------------------------
 *'					ReverseCompareLessThan		2003-12-5
 *-------------------------------------------------------------------------*/

static bool ReverseCompareLessThan( double mu1, double mu2 )
{
	return mu1 > mu2;			// Reverse the meaning of less than as the array is in reverse order
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosol::InterpolatePhaseMatrixTables		2008-2-27*/
/**	Calculate the phase matrix at arbitrary scattering angle  cos-1(mu).
 *	You may assume that the internal phase matrix tables have been pre-calculated
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosol::InterpolatePhaseMatrixTables( double mu, skRTPhaseMatrix* P)
{
	bool						ok;
	nxArrayIter<double>			begin;
	nxArrayIter<double>			end;
	nxArrayIter<double>			iter;
	const nx1dArray<double>*	mu_table;
	intptr_t					idx;
	intptr_t					idx1;
	double						x1;
	double						x2;
	double						y1;
	double						y2;
	double						f;
	double						pcoeffs[4];

	mu_table = &m_mie.Get_CosAngles();
	ok = (P != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_MieAerosol, Error sizing LapackMatrix to 4x4, pointer is NULL or thread safe memory allocation error");
	}
	else
	{
		P->SetTo(0.0);
		ok = (mu_table->size() > 1);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_MieAerosol, Insufficient scattering angles are defined for interpolation, setting phase matrix to zero");
		}
		else
		{
			NXTRACE_ONCEONLY(firsttime,("**** 2008-12-23 ***** skOpticalProperties_MieAerosol::InterpolatePhaseMatrixTables, We need to walk through this code\n"));
			begin = mu_table->begin();												// Get the beginning iterator
			end   = mu_table->end();												// Get the end iterator
			iter  = std::lower_bound( begin, end, mu, ReverseCompareLessThan );		// search for the bounding angle
			idx   = (iter - begin);													// get the distance of entry point from beginning
			idx1  = idx-1;															// first point is point before it
			if (idx1 < 0)															// Make sure we are in bounds
			{
				idx  = 1;
				idx1 = 0;
			}
			if ( idx >= (intptr_t)mu_table->size() )					// at the front an at the end
			{
				idx  = (intptr_t)mu_table->size()-1;
				idx1 = idx-1;
			}
			x1 = mu_table->At(idx1);						// now do a linear interpolation of the data
			x2 = mu_table->At(idx);
			f  = (mu - x1)/(x2-x1);
			for (int i=0; i < 4; i++)
			{
				y1          = m_Pij.At(i, idx1);
				y2          = m_Pij.At(i, idx );
				pcoeffs[i]  = y1 + f*(y2-y1);
			}
            P->At( 1,1 ) =  (SKRTFLOAT)(pcoeffs[0]);						// Get P11
            P->At( 2,2 ) =  (SKRTFLOAT)(pcoeffs[0]);						// Get P22
			P->At( 1,2 ) =  (SKRTFLOAT)(pcoeffs[1]);						// Get P12
			P->At( 2,1 ) =  (SKRTFLOAT)(pcoeffs[1]);						// Get P21
			P->At( 3,3 ) =  (SKRTFLOAT)(pcoeffs[2]);						// Get P33
			P->At( 4,4 ) =  (SKRTFLOAT)(pcoeffs[2]);						// Get P44
			P->At( 3,4 ) =  (SKRTFLOAT)(pcoeffs[3]);						// get P34
			P->At( 4,3 ) =  (SKRTFLOAT)(-pcoeffs[3]);						// Get P43
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_MieAerosol::CalculatePhaseMatrix		2008-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_MieAerosol::CalculatePhaseMatrix( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix)
{
	bool	ok;
	double	absxs, extxs, scattxs;

	m_mutex_xsectionslock.lock();
	ok =        CalculateCrossSectionsInternal( wavenumber, &absxs, &extxs, &scattxs);
	ok = ok && 	InterpolatePhaseMatrixTables( cosscatterangle,phasematrix);
	m_mutex_xsectionslock.unlock();
	return ok;
}

bool skOpticalProperties_MieAerosol::LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff) {
	if (usermaxcoeff > m_numlegendre)
	{
		nxLog::Record(NXLOG_ERROR, "skOpticalProperties_MieAerosol::LegendreCoefficientsP11 Requesting more legendre moments than were calculated by MIEV0");
		return false;
	}

	double temp1, temp2, temp3;
	m_mutex_xsectionslock.lock();
	CalculateCrossSectionsInternal(wavenumber, &temp1, &temp2, &temp3);
	m_mutex_xsectionslock.unlock();

	opticalmaxcoeff = std::min((size_t)usermaxcoeff, m_numlegendre);

	for (int i = 0; i < opticalmaxcoeff; ++i) {
		coeff[i] = m_legendre.At(0, i);
	}

	return true;
}