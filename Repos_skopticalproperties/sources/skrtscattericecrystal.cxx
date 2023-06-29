#include <skopticalproperties21.h>
#include <nxbase_threads.h>

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_IceCrystal::skOpticalProperties_IceCrystal		2009-6-4	*/
/*-------------------------------------------------------------------------*/

skOpticalProperties_IceCrystal::skOpticalProperties_IceCrystal()
{
	m_nonsphere       = NULL;
	m_distribution    = NULL;
	m_refractiveindex = NULL;
	init();
}

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_IceCrystal::~skOpticalProperties_IceCrystal		2009-6-4	 */
/*--------------------------------------------------------------------------*/

skOpticalProperties_IceCrystal::~skOpticalProperties_IceCrystal()
{
	ReleaseScatterAlgorithm();
	ReleaseDistribution();
	ReleaseRI();
}


/*---------------------------------------------------------------------------
 *						skOpticalProperties_IceCrystal::init				2009-6-4*/
/**	Inititalize the base scatterer class. By default we use 100 points in our
 *	radius distribution quadrature and we use 0.1 degree resolution for
 *	0 to 5 degrees and 175-180 scattering angles (forward scatter) and 1
 *	degree resolution for all others.
 */
/*--------------------------------------------------------------------------*/

void skOpticalProperties_IceCrystal::init()
{
	skRTParticleDist_LogNormal*	lognormal = new skRTParticleDist_LogNormal;

	m_quadrature.SetOrder(50);											// Default is to use 50 points when doing integration over particle size. (This is ignored by TMatrixRandom it is hard coded to 50)
	m_distribution		= NULL;
	m_refractiveindex	= NULL;
	lognormal->SetDistributionParameters( 0.1, 1.6, 0.0);
	Set_ScatterAlgorithm    ( new sk_TMatrixRandomWrapper );				// instantiate temporarily to setup extinction object
	Set_RefractiveIndex     ( new skRTRefractiveIndex_ICE );
	Set_ParticleDistribution( lognormal );
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystal::DeepCopy			2009-6-4		*/
/*---------------------------------------------------------------------------*/

/* 
bool skOpticalProperties_IceCrystal::DeepCopy( const skOpticalProperties_IceCrystal& other )
{
	bool	ok1;
	bool	ok2 = true;
	bool	ok3 = true;
	bool	ok4 = true;
	bool	ok5;

	ReleaseScatterAlgorithm();
	ReleaseDistribution();
	ReleaseRI();

	ok1 = skOpticalProperties::DeepCopy( other );
	if (other.m_nonsphere		!= NULL)  ok2 = other.m_nonsphere->CreateClone( &m_nonsphere );
	if (other.m_distribution    != NULL)  ok3 = other.m_distribution->CreateClone( &m_distribution );
	if (other.m_refractiveindex != NULL)  ok4 = other.m_refractiveindex->CreateClone( &m_refractiveindex);
	ok5   = m_quadrature.DeepCopy( other.m_quadrature );
	m_Pij = other.m_Pij;
	m_isdirtyproperties = true;
	return (ok1 && ok2 && ok3 && ok4 && ok5 );
}
*/

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystal::CreateClone				2009-6-4	*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skOpticalProperties_IceCrystal::CreateClone(skOpticalProperties** userclone) const
{
	skOpticalProperties_IceCrystal*	clone;
	bool						ok;

	clone = new skOpticalProperties_IceCrystal;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_IceCrystal::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
		ok = clone->DeepCopy( *this );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_IceCrystal::CreateClone, Error creating copying this object to clone object, this might cause significant issues");
		}
	}
	*userclone = clone;
	return ok;
}
*/

/*----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystal::Set_AspectRatio			2009-7-14 */
/**	Set the aspect ratio of the base nonspherical particle */
/*--------------------------------------------------------------------------*/
bool skOpticalProperties_IceCrystal::Set_AspectRatio( double aspectRatio )
{
	bool	ok;
	ok = m_nonsphere->Set_AspectRatio( aspectRatio );
	return ok;
}
/*----------------------------------------------------------------------------
 *'					skOpticalProperties_IceCrystal::ReleaseRI			2009-6-4		*/
/**	Release the current refractive index object								*/
/*--------------------------------------------------------------------------*/

void skOpticalProperties_IceCrystal::ReleaseRI()
{
	if (m_refractiveindex != NULL)
	{
		m_refractiveindex->Release();
		m_refractiveindex = NULL;
	}
}

/*---------------------------------------------------------------------------
 *'					skOpticalProperties_IceCrystal::Set_RefractiveIndex			2009-6-4 */
/**	Set the refractive index of the ice crystal to the new value.
 *
 *	@param newri
 *	Pointer to the new refractive index object.  The class will call AddRef()
 *	on the new refractive index to keep the object alive and it will call
 *	Release() when it is finished with the object.
 */
/*--------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystal::Set_RefractiveIndex( skRTRefractiveIndex* newri)
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
	return ok;
}


/*---------------------------------------------------------------------------
 *				skOpticalProperties_IceCrystal::ReleaseScatterAlgorithm		2009-6-4 */
/**	Release the curent scattering algorithm object.							*/
/*--------------------------------------------------------------------------*/

void skOpticalProperties_IceCrystal::ReleaseScatterAlgorithm()
{
	if (m_nonsphere != NULL)
	{
		m_nonsphere->Release();
		m_nonsphere = NULL;
	}
}

/*---------------------------------------------------------------------------
 *'				skOpticalProperties_IceCrystal::Set_ScatterAlgorithm		2009-6-4 */
/**	Set the scattering algorithm (T-matrix oriented, random, DDA) of the 
 *	ice crystal to the new value.
 *
 *	@param newscat
 *	Pointer to the new scattering object.  The class will call AddRef()
 *	on the new refractive index to keep the object alive and it will call
 *	Release() when it is finished with the object.
 */
/*--------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystal::Set_ScatterAlgorithm( sk_NonsphericalParticle* newscat )
{
	bool	ok;

	ok = (newscat == m_nonsphere);
	if (!ok)
	{
		if (newscat != NULL) newscat->AddRef();
		ReleaseScatterAlgorithm();
		m_nonsphere  = newscat;
		ok = (m_nonsphere != NULL);
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystal::ReleaseDistribution		2009-6-4 */
/**	Release the curent paticle distribution object.							*/
/*--------------------------------------------------------------------------*/

void skOpticalProperties_IceCrystal::ReleaseDistribution()
{
	if (m_distribution != NULL)
	{
		m_distribution->Release();
		m_distribution = NULL;
	}
}


/*---------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystal::Set_ParticleDistribution		2009-6-4 
 *	Set #newdistribution as the new particle distribution for the nonspherical
 *	scatterers.  The #newdistribution is held in memory by placing an AddRef()
 *	on the object and is freed when the class calls Release()
 *	If the object m_nonsphere is of type 'TMatrixRandom', then the PSD information
 *	is passed along to that class' members for calling the DLL appropriately.
 *
 *--------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystal::Set_ParticleDistribution( skRTParticleDist* newdistribution )
{
	bool	ok;
	
	ok	= (newdistribution == m_distribution);
	if (!ok)
	{
		if ( newdistribution != NULL ) newdistribution->AddRef();
		ReleaseDistribution();
		m_distribution  = newdistribution;
		ok = (m_distribution != NULL);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystal::IntegrateOverRandomTmatrixDistribution		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystal::IntegrateOverRandomTmatrixDistribution(double wavenum, double* absxs, double* extxs, double* scattxs, std::vector<skRTPhaseMatrix>* Pij  )
{
	static std::mutex			g_mutex_xsectionslock;
	sk_TMatrixRandomWrapper*	tmatrixrandom;
	std::complex<double>		ri;
	bool						ok;
	double						lamdamicron = 1.0E4/wavenum;		// Get the wavelength in microns

	ok =		( m_nonsphere       != NULL)
			&&  ( m_distribution    != NULL)
			&&  ( m_refractiveindex != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_IceCrystal::IntegrateOverRandomTmatrixDistribution, cannot call Tmatrix code as this object is not properly initialized");
	}
	else
	{
		std::unique_lock<std::mutex>	lock(g_mutex_xsectionslock);				// The code is mutexed to make it thread safe. This is "okay" as we only call this code when we need to make mie scattering tables
																						// So we dont have to worry about parallelism and speed execution for the moment.	

		tmatrixrandom = dynamic_cast<sk_TMatrixRandomWrapper*>(m_nonsphere);				// Cast the non spherical object back into a sk_TMatrixRandomWrapper object
		ri  =       m_refractiveindex->RefractiveIndex        ( wavenum );					// Get the refractive index
		ok  =       tmatrixrandom->Set_RefractiveIndex        ( ri.real(), ri.imag() );		// and set it
		ok  = ok &&	tmatrixrandom->Set_Wavelength             ( lamdamicron  );				// Set the scattering wavelnth to microns (Compatible with particle distribution radius )
		ok  = ok && tmatrixrandom->Set_ParticleDistribution	  ( m_distribution );			// Set the particle distribution
		ok  = ok && tmatrixrandom->CalculateScattering        ();
		ok  = ok &&	tmatrixrandom->Get_ScatteringMatrixCoeffs ( Pij);				// Update Scattering coeffs and get scattering matrix coeffs calculated by the T-matrix code.
		if (ok)
		{
			NXTRACE_ONCEONLY(firsttime,("***** CHECK Phase matrix is  phase and not Mueller Matrix in skOpticalProperties_IceCrystal::IntegrateOverRandomTmatrixDistribution\n"));
			*absxs   = m_nonsphere->Cabs();
			*extxs   = m_nonsphere->Cext();
			*scattxs = m_nonsphere->Csca();
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_IceCrystal::IntegrateOverRandomTmatrixDistribution, There was an error calculating cross-sections and phase matrix. This is a problem");
		*absxs =  0.0;
		*extxs =  0.0;
		*scattxs =  0.0;
		 Pij->clear();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystal::IntegrateOverSingleRadiusTmatrixCode		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystal::IntegrateOverSingleRadiusCode(double wavenum, double* absxs, double* extxs, double* scattxs, std::vector<skRTPhaseMatrix>* Pij   )
{
	static std::mutex			g_mutex_xsectionslock;
	std::complex<double>			ri;
	std::vector<skRTPhaseMatrix>	fij;
	bool							ok;
	double					lamdamicron = 1.0E4/wavenum;		// Get the wavelength in microns
	size_t					numangles;
	double					nr;
	double					nsum;
	double					c;
	double					minradius = 0.0;
	double					maxradius = 0.0;
	size_t					i;
	double					ksca;
	double					kabs;
	double					kext;
	size_t					iangle;

	std::unique_lock<std::mutex>	lock(g_mutex_xsectionslock);				// The code is mutexed to make it thread safe. This is "okay" as we only call this code when we need to make mie scattering tables
																						// So we dont have to worry about parallelism and speed execution for the moment.	


	ri =       m_refractiveindex->RefractiveIndex	( wavenum );					// Get the refractive index at this wavenumber
	ok =       m_nonsphere->Set_RefractiveIndex		( ri.real(), ri.imag() );		// Set the refractive index in the Tmatrix code
	ok = ok && m_nonsphere->Set_Wavelength			( lamdamicron  );				// Set the wavelength in the TMatrix code
	ok = ok && m_distribution->GetQuadratureRadii	( &minradius, &maxradius );		// Getthe range of radii from the distribution
	m_quadrature.SetRange							( minradius, maxradius );		// Set the quadrature range of points
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_IceCrystal::CalculateCrossSections, There was an error selecting the range of radii for the quadrature");
		*absxs   = 0.0; 
		*extxs   = 0.0;
		*scattxs = 0.0;
	}
	else
	{
		const nx1dArray<double>&	R = m_quadrature.X();			// R are the quadrature points in the radius distribution
		const nx1dArray<double>&	W = m_quadrature.W();			// W are the desired weights   in the radius distribution quadrature

		numangles = m_nonsphere->Get_NumAngles();
		Pij->resize(numangles );				// The Fij scattering array calculated by m_nonsphere.( 16, num_Cosangles, )
		nsum    = 0.0;
		ksca    = 0.0;
		kext    = 0.0;
		for (i = 0; i < R.size(); i++)								// For each radius in the radius quadrature
		{															// then
			m_nonsphere->Set_Radius( R.At(i) );						// configure the scattering radius
			m_nonsphere->Get_ScatteringMatrixCoeffs( &fij );		// Get the Fij calculated by the code.
			nr      = W.At(i)*m_distribution->Distribution( R.At(i) );	// Get n(r).w for evaluating the quadrature
			nsum   += nr;
			for (iangle = 0; iangle < numangles; iangle++)
			{
				fij[iangle]      *= nr;											// Get Fij.n(r).w
				Pij->at(iangle)  += fij[iangle];											// Fij = Sum( Fij.n(r).w )
			}
			ksca   += (nr*m_nonsphere->Csca());						// ksca = Sum
			kext   += (nr*m_nonsphere->Cext());						// kext = Sum
		}
		ok = ok && ((nsum > 0.99) && (nsum < 1.01));
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_IceCrystal::CalculateCrossSections, The quadrature of the particle distribution was not normalized to unity with good accuracy. we got (%e) instead of 1.0", (double)nsum);
			NXASSERT(( ok ));
		}
		c   = 4.0*nxmath::Pi/nxmath::sqr(m_nonsphere->Get_k());		// Get the normalization constant from the Mueller matrix
		for (iangle = 0; iangle < numangles; iangle++)
		{
			Pij->at(iangle) *= (c/ksca);									// and apply it to get the phase matrix
		}
		kabs  = kext - ksca;										// Convert from square microns to square centimeters
		*absxs   = kabs; 
		*extxs   = kext;
		*scattxs = ksca;
	}
	return ok;
}
/*---------------------------------------------------------------------------
 *'					skOpticalProperties_IceCrystal::CalculateCrossSections		2009-6-4 */
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
*/
/*-------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystal::CalculateCrossSectionsInternal( double wavenum, double* absxs, double* extxs, double* scattxs, std::vector<skRTPhaseMatrix>* Pij )
{
	bool					ok;

	ok = (m_refractiveindex != NULL);										// Make sure the refractive index object is not NULL
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_IceCrystal, Cannot calculate the cross sections as no refractive index is defined (PLEASE DEFINE ONE with Set_RefractiveIndex)");
	}
	else
	{
		switch ( m_nonsphere->Algorithm())
		{
			case sk_NonsphericalParticle::ALG_TMATRAND :	ok = IntegrateOverRandomTmatrixDistribution( wavenum, absxs, extxs, scattxs, Pij  );
															break;
			case  sk_NonsphericalParticle::ALG_TMATORIENT : ok = IntegrateOverSingleRadiusCode( wavenum, absxs, extxs, scattxs, Pij  );
															break;
			case  sk_NonsphericalParticle::ALG_DDA        :	ok = IntegrateOverSingleRadiusCode( wavenum, absxs, extxs, scattxs, Pij  );
															break;
			default:
															nxLog::Record(NXLOG_WARNING, "skOpticalProperties_IceCrystal::CalculateCrossSections, Unrecognized non spherical particle algorithm");
															ok = false;
															break;
		};
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystal::LookupUpThreadData		 2014- 10- 22*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystal::LookupUpThreadData(std::vector<skRTPhaseMatrix>** data )
{
	size_t					threadnum;
	iterator				iter;
	bool					ok;
	static std::mutex		lock;							// Add a lock to make sure map.insert is thread safe.

	threadnum = nxWorkerThreadManager::GetCurrentThreadIdCode();
	lock.lock();
	iter      = m_PijArray.find(threadnum);
	ok = (iter != m_PijArray.end() );
	if (!ok)
	{
		std::pair<iterator, bool>							result;
		std::pair<size_t, std::vector<skRTPhaseMatrix> >	newentry(threadnum, std::vector<skRTPhaseMatrix>() );
		result = m_PijArray.insert(  newentry);
		ok = result.second;
		iter = result.first;

	}
	lock.unlock();
	if (!ok)
	{
		*data = nullptr;
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_IceCrystal::LookupUpThreadData, error looking/creating thread data for thread %d", (int)threadnum);
	}
	else
	{
		*data = &(*iter).second;
	}
	return ok;

}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystal::CalculateCrossSections		2014-2-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystal::CalculateCrossSections( double wavenum, double* absxs, double* extxs, double* scattxs )
{
	bool ok;
	std::vector<skRTPhaseMatrix>*	threadpij;

	ok = LookupUpThreadData(&threadpij);
	ok = ok && CalculateCrossSectionsInternal( wavenum, absxs, extxs, scattxs, threadpij);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ReverseCompareLessThan		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

static bool ReverseCompareLessThan( double mu1, double mu2 )
{
	return mu1 > mu2;			// Reverse the meaning of less than as the array is in reverse order
}

/*-----------------------------------------------------------------------------
 *				skOpticalProperties_IceCrystal::InterpolatePhaseMatrixTables		2009-6-4*/
 /**	Calculate the phase matrix at arbitrary scattering angle  cos-1(mu).
 *	You may assume that the internal phase matrix tables have been pre-calculated
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystal::InterpolatePhaseMatrixTables( double mu, skRTPhaseMatrix* P, std::vector<skRTPhaseMatrix>* Pij ) const
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
	double						f1;
	double						f2;

	mu_table = &m_nonsphere->Get_CosAngles();
	ok = (P != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_IceCrystal, Error sizing LapackMatrix to 4x4, pointer is NULL or memory allocation error");
	}
	else
	{
		P->SetTo(0.0);
		ok = (mu_table->size() > 1);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_IceCrystal, Insufficient scattering angles are defined for interpolation, setting phase matrix to zero");
		}
		else
		{
			NXTRACE_ONCEONLY(firsttime,("**** 2009-6-4 ***** skOpticalProperties_IceCrystal::InterpolatePhaseMatrixTables, We need to walk through this code\n"));
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
			x1 = mu_table->At(idx1);									// now do a linear interpolation of the data
			x2 = mu_table->At(idx);
			f1  = (x2 - mu)/(x2-x1);
			f2  = (mu - x1)/(x2-x1);
			*P  = Pij->at(idx1)*f1 + Pij->at(idx)*f2;
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystal::CalculatePhaseMatrix		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystal::CalculatePhaseMatrix( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix)
{
	bool	ok;
	double 	absxs, extxs, scattxs;
	std::vector<skRTPhaseMatrix>*	threadpij;

	ok = LookupUpThreadData(&threadpij);
	ok = ok &&  CalculateCrossSectionsInternal  ( wavenumber, &absxs, &extxs, &scattxs, threadpij);
	ok = ok && 	InterpolatePhaseMatrixTables( cosscatterangle,phasematrix, threadpij);
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystal::SetAtmosphericState		2009-11-10*/
/** The ice crystal size distributions and associated cross-sections have no 
	*temperature, pressure or altitude dependency.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystal::SetAtmosphericState( skClimatology* /*neutralatmosphere*/)
{
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_IceCrystal::SetLocation		 2015- 11- 17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_IceCrystal::SetLocation( const GEODETIC_INSTANT& /*pt*/, bool* crosssectionschanged  )
{
	if (crosssectionschanged != NULL) *crosssectionschanged = false;
	return true;
}

