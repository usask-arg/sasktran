#include <skopticalproperties21.h>
#include <mutex>

//using namespace std;

static std::mutex	g_mutexlock;

#if !defined(NX_WINDOWS)
 	#define MIEV0 miev0_
#endif

extern "C" void  MIEV0  (	double*					XX,				// REAL
							complex<double>*		CREFIN,			// COMPLEX
							char*					PERFCT,			// LOGICAL
							double*					MIMCUT,			// REAL
							char*					ANYANG,			// LOGICAL
							int*					NUMANG,			// INTEGER
							double*					XMU,			// REAL(*)
							int*					NMOM,			// INTEGER
							int*					IPOLZN,			// INTEGER
							int*					MOMDIM,			// INTEGER
							char*					PRNT,			// LOGICAL(*)
							double*					QEXT,			// REAL
							double*					QSCA,			// REAL
							double*					GQSC,			// REAL
							double*					PMOM,			// REAL(0:MOMDIM,*)
							complex<double>*		SFORW,			// COMPLEX
							complex<double>*		SBACK,			// COMPLEX
							complex<double>*		S1,				// COMPLEX (NUMANG))
							complex<double>*		S2,				// COMPLEX (NUMANG)
							complex<double>*		TFORW,			// COMPLEX (2)
							complex<double>*		TBACK,			// COMPLEX (2)
							double*					SPIKE );		// REAL




/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::sk_MieSphericalWiscombeWrapper		2003-10-2
 *-------------------------------------------------------------------------*/

sk_MieSphericalWiscombeWrapper::sk_MieSphericalWiscombeWrapper()
{
	m_radius  = 1.0E-04;						// Radius      = 1 micron by default = 10-4 cms
	m_lambda   = 557.0E-07;						// wavelength  = 557 nm
	UpdateXX();									// Update the m_xx field
	m_refr    = complex<double>(1.33,0.000);	// The real part of the refractive index (real,imaginary)
	m_perfect = false;;						// TRUE if the refractive index is infinite Wiscombe's PERFCT
	m_mimcut  = 0.0;							// Value below which imaginary refractive index is 0.
	m_anyang  = false;						// If TRUE then any angle can be placed in m_xmu. If FALSE then angles are mirror symmetric about 90
	m_nmom    = -1;								// The highest Legendre moment to calculate
	m_ipolzn  = +4321;							// Wicombe IPOLZN variable
	m_qext    = 0.0;							// Output QEXT from Wiscombes code
	m_qsca    = 0.0;							// Output QSCA from Wiscombes code
	m_gqsc    = 0.0;							// Output GSQC from Wiscombes code
	m_spike   = 1.0;
	m_Cext    = 0.0;							// Extinction cross-section (in units of radius**2)
	m_Cabs    = 0.0;							// Absorption cross-section (in units of radius**2)
	m_Csca    = 0.0;							// Scattering cross-section (in units of radius**2)
	Set_MaxLegendreMoment(0);
	m_isdirty = true;
}



/*-----------------------------------------------------------------------------
 *					sk_MieSphericalWiscombeWrapper::operator=		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::DeepCopy( const sk_MieSphericalWiscombeWrapper& other )
{
	bool	ok;
	bool	ok1;
	bool	ok2;
	bool	ok3;
	bool	ok4;

	m_radius   = other.m_radius;
	m_lambda   = other.m_lambda;
	m_k	       = other.m_k;
	m_xx	   = other.m_xx;
	m_refr	   = other.m_refr;
	m_perfect  = other.m_perfect;
	m_mimcut   = other.m_mimcut;
	m_anyang   = other.m_anyang;
	ok1        = m_xmu.DeepCopy( other.m_xmu );
	m_nmom	   = other.m_nmom;
	m_ipolzn   = other.m_ipolzn;
	ok2        = m_pmom.DeepCopy( other.m_pmom );
	m_qext	   = other.m_qext;
	m_qsca	   = other.m_qsca;
	m_qabs	   = other.m_qabs;
	m_Cext	   = other.m_Cext;
	m_Cabs	   = other.m_Cabs;
	m_Csca	   = other.m_Csca;
	m_gqsc	   = other.m_gqsc;
	ok3        = m_S1.DeepCopy( other.m_S1 );
	ok4        = m_S2.DeepCopy( other.m_S2 );
	m_Sforw	   = other.m_Sforw;
	m_Sback	   = other.m_Sback;
	m_Tforw[0] = other.m_Tforw[0];
	m_Tforw[1] = other.m_Tforw[1];
	m_Tback[0] = other.m_Tback[0];
	m_Tback[1] = other.m_Tback[1];
	m_spike	   = other.m_spike;
	m_isdirty  = other.m_isdirty;
	ok = (ok1 && ok2 && ok3 && ok4);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_MieSphericalWiscombeWrapper::DeepCopy, Error copying contents, this object is in an ill-defined state and may produce strange results");
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::SetDefaultAngles		2003-10-2
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::SetDefaultAngles()
{
	double	startangle[2] = {0.0,  10.0};
	double  endangle[2]   = {10.0, 90.0};
	double  resolution[2] = {0.1,   0.25};
    bool ok;

	ok =  Set_SymmetricScatterAngles( startangle, endangle, resolution, 2 );
	return ok;
}



/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::CalculateScattering		2003-10-2
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::CalculateScattering()
{
	bool	ok;

	ok = (!m_isdirty);
	if (!ok)
	{
		ok = (m_xmu.size() > 0);
		if (!ok)
		{
			nxLog::Verbose(NXLOG_WARNING, "sk_MieSphericalWiscombeWrapper::CalculateScattering, No scatterings defined so creating defaults");
			ok = SetDefaultAngles();
		}
		if (ok)
		{
			Wiscombe_Miev0();
		}
		m_isdirty = !ok;
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Wiscombe_Miev0		2003-10-1
 *-------------------------------------------------------------------------*/

void sk_MieSphericalWiscombeWrapper::Wiscombe_Miev0()
{
	char	PRNT[8]    = {0,0,0,0,0,0,0,0};
	char	PERFECT[4];
	char	ANYANG[4];
	int		NUMANG     = (int)m_xmu.size();
//	int		IPOLZN     = +1234;
	int		MOMDIM	   = m_nmom + 1;
	char	boolval;


	boolval = m_perfect ? -1 : 0;
	memset( PERFECT, boolval, sizeof(PERFECT));
	boolval = m_anyang ? -1 : 0;
	memset( ANYANG, boolval, sizeof(ANYANG));

	std::lock_guard<std::mutex>	lock(g_mutexlock);

	MIEV0(	&m_xx,						// double*	XX,				// REAL
			&m_refr,					// double*	CREFIN,			// COMPLEX
			&PERFECT[0],				// char*	PERFCT,			// LOGICAL
			&m_mimcut,					// double*	MIMCUT,			// REAL
			&ANYANG[0],					// char*	ANYANG,			// LOGICAL
			&NUMANG,					// INTEGER
			m_xmu.UnsafeArrayBasePtr(),	// double*	XMU,			// REAL(*)
			&m_nmom,					// int*		NMOM,			// INTEGER
			&m_ipolzn,					// int*		IPOLZN,			// INTEGER
			&MOMDIM,					// int*		MOMDIM,			// INTEGER
			&PRNT[0],					// LOGICAL(*)
			&m_qext,					// double*	QEXT,			// REAL
			&m_qsca,					// double*	QSCA,			// REAL
			&m_gqsc,					// double*	GQSC,			// REAL
			m_pmom.UnsafeArrayBasePtr(),// double*	PMOM,			// REAL(0:MOMDIM,*)
			&m_Sforw,					// double*	SFORW,			// COMPLEX
			&m_Sback,					// double*	SBACK,			// COMPLEX
			m_S1.UnsafeArrayBasePtr(),		// double*	S1,				// COMPLEX (NUMANG))
			m_S2.UnsafeArrayBasePtr(),		// double*	S2,				// COMPLEX (NUMANG)
			&m_Tforw[0],				// double*	TFORW,			// COMPLEX(*)
			&m_Tback[0],				// double*	TBACK,			// COMPLEX(*)
			&m_spike					// double* SPIKE );			// REAL
			);

	double area = nxmath::Pi*m_radius*m_radius;

	m_qabs = m_qext - m_qsca;				// Absorption efficiency
	m_Cext = m_qext*area;					// Extinction cross-section (in units of radius**2)
	m_Cabs = m_qabs*area;					// Absorption cross-section (in units of radius**2)
	m_Csca = m_qsca*area;					// Scattering cross-section (in units of radius**2)

}

/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::AllocateArrays		2003-10-1
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::AllocateArrays( size_t numangles )
{
	bool	ok;

	ok =        m_xmu.SetSize( numangles    ) && (numangles < 1001);
	ok = ok &&   m_S1.SetSize( numangles );
	ok = ok &&   m_S2.SetSize( numangles );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "sk_MieSphericalWiscombeWrapper::AllocateArrays, Error allocating arrays for %d angles (maximum number of angles is 1001)", (int)numangles);
		m_xmu.erase();
		m_S1.erase();
		m_S2.erase();
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::ConfigureSymmetricScatterAngles		2003-10-1
 *	Define the scattering angles which must be calculated. The user can provide
 *	a list of angular sections with different resoutions in each section.
 *	This allows the forward scatter peak near 0 degrees to be accurately
 *	modeled at high resolution while the other angles can use a much coarser
 *	resolution.
 *
 *	PARAMETERS
 *	startangle dblarr[numsteps]
 *		Specifies the start angle of eachregion in degrees. Must be between
 *		0 and less than 90.  I implicitly assume that the startangle of one
 *		section is greater than the endangle of the last section
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::Set_SymmetricScatterAngles( double *startangle,
														 double *endangle,
													     double *resolution,
													     int	 numrange )
{
	double		lastend = 0.0;
	int			numangles;
	int			nsteps;
	int			totalangles;
	bool		ok = true;
	int			i;
	int			j;
	double		theta;
	int			fidx;
	int			bidx;

	m_anyang = false;
	numangles = 0;
	for (i=0; i < numrange && ok; i++)
	{
		ok = (startangle[i] >= lastend) && (endangle[i] <= 90.0);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "sk_MieSphericalWiscombeWrapper::ConfigureScatterAngles, The startangle (%g) must be less than lastangle (%g) and endangle (%g) must be <= 90 degrees", (double)startangle[i], (double)lastend, (double)endangle[i]);
		}
		else
		{
			nsteps     = (int)((endangle[i] - startangle[i])/resolution[i] + 0.5);
			numangles += nsteps;
			lastend    = startangle[i] + (nsteps-1)*resolution[i];;
		}
	}
	if (ok)
	{
		totalangles = numangles *2;								// Account for 90 to 180 degree sector
		totalangles++;											// last angle is not exactly 90 (this is always true)
		ok = AllocateArrays(totalangles);						// Allocate the space for the cosine of scattering angles

		numangles = 0;
		fidx      = 0;
		bidx      = totalangles-1;
		for (i=0; i < numrange && ok; i++)
		{
			nsteps   = (int)((endangle[i] - startangle[i])/resolution[i] + 0.5);
			for (j = 0; j < nsteps; j++)
			{
				theta = nxmath::cosd(startangle[i] + j*resolution[i]);
				m_xmu.At(fidx++) =  theta;
				m_xmu.At(bidx--) = -theta;
			}
		}
		m_xmu.At(fidx) = 0;
	}
	if (!ok) m_xmu.erase();
	m_isdirty = true;
	return ok;
}

/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::ConfigureAnyScatterAngles		2003-10-1
 *	Define the scattering angles which must be calculated. The user just specifies
 *	the scattering angles.  They are not necessarily symmetric
 *
 *	PARAMETERS
 *	startangle dblarr[numsteps]
 *		Specifies the start angle of eachregion in degrees. Must be between
 *		0 and less than 90.  I implicitly assume that the startangle of one
 *		section is greater than the endangle of the last section
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::Set_AnyScatterAngles( nx1dArray<double>& degrees)
{
	size_t			numangles;
	bool		ok;

	m_anyang  = true;
	numangles = degrees.size();
	ok        = AllocateArrays(numangles);						// Allocate the space for the cosine of scattering angles
	if (ok)
	{
		for (size_t i=0; i < numangles; i++)
		{
			m_xmu.At(i) = nxmath::cosd( degrees[i] );
		}
	}
	else m_xmu.erase();
	m_isdirty = true;
	return ok;
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Set_Wavelength		2003-10-7
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::Set_Wavelength( double lamda )
{
	SetDirty (m_lambda != lamda );
	m_lambda = lamda;
	UpdateXX();
	return true;
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Set_Radius		2003-10-7
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::Set_Radius	( double radius )
{
	SetDirty (m_radius != radius );
	m_radius = radius;
	UpdateXX();
	return true;
}

/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Set_IsInfiniteRefractiveIndex		2003-10-7
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::Set_IsInfiniteRefractiveIndex	( bool val )
{
	SetDirty (m_perfect != val );
	m_perfect = val;
	return true;
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Set_RefractiveIndex		2003-10-7
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::Set_RefractiveIndex( double ri_real, double ri_imag )
{
	complex<double>	val( ri_real, ri_imag);

	SetDirty( m_refr != val );
	m_refr = val;
	return true;
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::SetMaxLegendreMoment		2003-10-7
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::Set_MaxLegendreMoment(int nmom )
{
	bool	ok;

	ok = (m_nmom == nmom);
	if (!ok)
	{
		m_isdirty = true;
		m_nmom = nxmax(0,nmom);
		ok = m_pmom.SetSize( m_nmom+2, 4 );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "sk_MieSphericalWiscombeWrapper::SetMaxLegendreMoment, Error allocating memory for Legendre moments");
			m_nmom = -1;
		}
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Set_IPolzn		2003-10-7
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::Set_IPolzn( int ipolzn )
{
	SetDirty( ipolzn != m_ipolzn );
	m_ipolzn = ipolzn;
	return true;
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Set_MimCut		2003-10-7
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::Set_MimCut( double mimcut )
{
	SetDirty( mimcut != m_mimcut );
	m_mimcut = mimcut;
	return true;
}

/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Qext		2003-10-7
 *-------------------------------------------------------------------------*/

double sk_MieSphericalWiscombeWrapper::Qext()
{
	CalculateScattering();
	return m_qext;
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Qsca		2003-10-7
 *-------------------------------------------------------------------------*/

double sk_MieSphericalWiscombeWrapper::Qsca()
{
	CalculateScattering();
	return m_qsca;
}

double sk_MieSphericalWiscombeWrapper::Qabs()
{
	CalculateScattering();
	return m_qabs;
}

/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Csca		2003-10-7
 *-------------------------------------------------------------------------*/

double sk_MieSphericalWiscombeWrapper::Csca()
{
	CalculateScattering();
	return m_Csca;
}

/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Cabs		2003-10-7
 *-------------------------------------------------------------------------*/

double sk_MieSphericalWiscombeWrapper::Cabs()
{
	CalculateScattering();
	return m_Cabs;
}

/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Cext		2003-10-7
 *-------------------------------------------------------------------------*/

double sk_MieSphericalWiscombeWrapper::Cext()
{
	CalculateScattering();
	return m_Cext;
}
/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Gqsc		2003-10-7
 *-------------------------------------------------------------------------*/

double sk_MieSphericalWiscombeWrapper::Gqsc()
{
	CalculateScattering();
	return m_gqsc;
}



/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::PMom		2003-10-7
 *-------------------------------------------------------------------------*/

nx2dArray<double>* sk_MieSphericalWiscombeWrapper::PMom()
{
	CalculateScattering();
	return &m_pmom;
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::S1		2003-10-7
 *-------------------------------------------------------------------------*/

nx1dArray< std::complex<double> >*	sk_MieSphericalWiscombeWrapper::S1()
{
	CalculateScattering();
	return &m_S1;
}

/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::S2		2003-10-7
 *-------------------------------------------------------------------------*/

nx1dArray< std::complex<double> >*	sk_MieSphericalWiscombeWrapper::S2()
{
	CalculateScattering();
	return &m_S2;

}

/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::SForward		2003-10-7
 *-------------------------------------------------------------------------*/

std::complex<double> sk_MieSphericalWiscombeWrapper::SForward	()
{
	CalculateScattering();
	return m_Sforw;
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::SBackward		2003-10-7
 *-------------------------------------------------------------------------*/

std::complex<double> sk_MieSphericalWiscombeWrapper::SBackward()
{
	CalculateScattering();
	return m_Sback;
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::TForward		2003-10-7
 *-------------------------------------------------------------------------*/

std::complex<double> sk_MieSphericalWiscombeWrapper::TForward	(int i)
{
	if ((i < 0) || (i > 1))
	{
		nxLog::Record( NXLOG_WARNING, "sk_MieSphericalWiscombeWrapper::TForward, index I is out of range (0..1) using value 0");
		i = 0;
	}
	CalculateScattering();
	return m_Tforw[i];
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::TBackward		2003-10-7
 *-------------------------------------------------------------------------*/

std::complex<double> sk_MieSphericalWiscombeWrapper::TBackward(int i)
{
	if ((i < 0) || (i > 1))
	{
		nxLog::Record( NXLOG_WARNING, "sk_MieSphericalWiscombeWrapper::TBackward, index I is out of range (0..1) using value 0");
		i = 0;
	}
	CalculateScattering();
	return m_Tback[i];
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Spike		2003-10-7
 *-------------------------------------------------------------------------*/

double sk_MieSphericalWiscombeWrapper::Spike()
{
	CalculateScattering();
	return m_spike;
}


/*---------------------------------------------------------------------------
 *'					sk_MieSphericalWiscombeWrapper::Get_ScatteringMatrixCeoffs		2003-10-9
 *	Get the array of scattering matrices for this mie scattering particle.
 *	Given that our Mie scattering code works out the scattering amplitude it is a
 *	simple matter to work out the Mueller scattering matrix.
 *
 *	This code follows eqn 4.77 of Bohren and Huffman: Absorption and
 *	Scattering of light by small particles, which is the same as eqn 2.36
 *	of Hansen and Travis (1974)
 *
 *	Note to convert the scattering matrix to the phase matrix you need
 *	a normalization constant. The constant is given by
 *	eqns 2.38 to 2.40 in Hansen and Travis 1974.
 *
 *	However if you are working with particle size distributions then the
 *	normalization constant is done after integrating over the scattering matrices
 *	for different particle sizes (see  eqn 2.49 in Hansen and Travis 1974).
 *
 *	matrixcoeff(0,*) = S11		*= all angles in m_xmu (see Get_CosAngles())
 *	matrixcoeff(1,*) = S12		*= all angles in m_xmu
 *	matrixcoeff(2,*) = S33		*= all angles in m_xmu
 *  matrixcoeff(3,*) = S34		*= all angles in m_xmu
 *
 *	Phase matrix can be reconstructed using
 *
 *		( S11   S12     0    0 )
 *		( S12   S11     0    0 )
 *		(   0     0   S33  S34 )
 *		(   0     0  -S34  S33 )
 *
 *
 *-------------------------------------------------------------------------*/

bool sk_MieSphericalWiscombeWrapper::Get_ScatteringMatrixCoeffs( nx2dArray<double>* coeffmatrix )
{
	bool								ok;
	complex<double>						S1;
	complex<double>						S2;
	double								s1sqr;
	double								s2sqr;
	double								s1x;
	double								s1y;
	double								s2x;
	double								s2y;
	size_t									i;
//	double								c;

	CalculateScattering();								// Update Cache so output parameters reflect latest input parameters
	ok = (m_xmu.size() > 0);
//	c   = (4.0*nxmath::Pi)/(m_k*m_k*m_Csca);			// Get the normalization constant reciprocal of equation 2.40 in Hansen And Travis
	coeffmatrix->SetSize(4, m_xmu.size());
	for (i=0; i < m_xmu.size(); i++)
	{
		S1    = m_S1[i];
		S2    = m_S2[i];
		s1sqr = norm(S1);
		s2sqr = norm(S2);
		s1x   = S1.real();
		s1y   = S1.imag();
		s2x   = S2.real();
		s2y   = S2.imag();

		coeffmatrix->At(0,i) = 0.5*( s1sqr   + s2sqr   );		// S11
		coeffmatrix->At(1,i) = 0.5*( s2sqr   - s1sqr   );		// S12
		coeffmatrix->At(2,i) =     ( s2x*s1x + s2y*s1y );		// S33
		coeffmatrix->At(3,i) =     ( s1x*s2y - s1y*s2x );		// S34
	}
	return true;
}

bool sk_MieSphericalWiscombeWrapper::Get_LegendreCoefficients(nx2dArray<double>* legendrecoeff)
{
	CalculateScattering();								// Update Cache so output parameters reflect latest input parameters
	legendrecoeff->SetSize(4, Get_MaxLegendreMoment());
	for (int i = 0; i < Get_MaxLegendreMoment(); i++)
	{
		legendrecoeff->At(0, i) = 0.5 * (m_pmom.At(i, 0) + m_pmom.At(i, 1));

		// TODO: Implement polarized legendre moments
		legendrecoeff->At(1, i) = 0.0;
		legendrecoeff->At(2, i) = 0.0;
		legendrecoeff->At(3, i) = 0.0;
	}
	return true;
}