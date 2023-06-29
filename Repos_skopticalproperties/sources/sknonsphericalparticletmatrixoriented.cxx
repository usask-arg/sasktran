#include <skopticalproperties21.h>
#include <nxbase_threads.h>

using namespace std;

static std::mutex	g_mutexlock;
//static nxMutex 		g_mutexlock;				// A mutex to avoid calling fortran code from multiple threads


   extern "C" void  TMATRIXORIENTED(	double*				AXI,		// REAL*8
										double*				RAT,		// REAL*8
										double*				LAM,		// REAL*8
										double*				MRR,		// REAL*8
										double*				MRI,		// REAL*8
										double*				EPS,		// REAL*8
										int*				NP,			// INTEGER                                             
										double*				DDELT,		// REAL*8
										int*				NDGS,		// INTEGER
										double*				ALPHA,		// REAL*8
										double*				BETA,		// REAL*8
										double*				THETA0,		// REAL*8
										double*				THET,		// REAL*8
										double*				PHI0,		// REAL*8
										double*				PHI,		// REAL*8
										int*				NUMSCAT,	// INTEGER
										double*				XMU,		// REAL*8
										complex<double>*	S11,		// COMPLEX*16
										complex<double>*	S12,		// COMPLEX*16
										complex<double>*	S21,		// COMPLEX*16
										complex<double>*	S22,		// COMPLEX*16
										double*				QSCA,
										double*				QEXT,
										double*				WALB
									);		



/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::sk_TMatrixOrientedWrapper		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

sk_TMatrixOrientedWrapper::sk_TMatrixOrientedWrapper()
{
	m_rat			= 1.0;							// equal-volume sphere by default
	m_np			= -2;							// cylinder
	m_ndgs			= 2;
	m_alpha			= 90;
	m_beta			= 0;
	m_thetaInc		= 56;
	m_thetaSca		= 65;
	m_phiInc		= 0;
	m_phiSca		= 128;
	m_ssalb         = -9999;
	SetDefaultAngles();
}


/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::~sk_TMatrixOrientedWrapper		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

sk_TMatrixOrientedWrapper::~sk_TMatrixOrientedWrapper()
{
	m_S11.erase();			
	m_S12.erase();			
	m_S21.erase();			
	m_S22.erase();			
}

/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::SetDefaultAngles		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixOrientedWrapper::SetDefaultAngles()
{
	double	startangle[2] = {0.0,  10.0};
	double  endangle[2]   = {10.0, 90.0};
	double  resolution[2] = {0.1,  0.25};
    bool ok;

	ok =  Set_SymmetricScatterAngles( startangle, endangle, resolution, 2 );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::DeepCopy		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixOrientedWrapper::DeepCopy( const sk_TMatrixOrientedWrapper & other )
{
	bool	ok;

	ok = sk_NonsphericalParticle::DeepCopy( other );
	m_np			= other.m_np;			
	m_ndgs			= other.m_ndgs;			
	m_rat			= other.m_rat;			
	m_alpha			= other.m_alpha;
	m_beta			= other.m_beta;
	m_thetaInc		= other.m_thetaInc;
	m_thetaSca		= other.m_thetaSca;
	m_phiInc		= other.m_phiInc;
	m_phiSca		= other.m_phiSca;
	m_ssalb         = other.m_ssalb;
	ok				= ok && m_S11.DeepCopy( other.m_S11 );
	ok				= ok && m_S12.DeepCopy( other.m_S12 );
	ok				= ok && m_S21.DeepCopy( other.m_S21 );
	ok				= ok && m_S22.DeepCopy( other.m_S22 );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_TMatrixOrientedWrapper::DeepCopy, Error copying contents, this object is in an ill-defined state and may produce strange results");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::CreateClone		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixOrientedWrapper::CreateClone(sk_NonsphericalParticle** userclone) const
{
	sk_TMatrixOrientedWrapper*	clone;
	bool						ok;

	clone = new sk_TMatrixOrientedWrapper;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_TMatrixOrientedWrapper::CreateClone, Error allocating memeory for newly cloned object");
	}
	else
	{
		clone->AddRef();
		ok = clone->DeepCopy(*this);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"sk_TMatrixOrientedWrapper::CreateClone, Error deep copying the contents of this to the newly created object");
			clone->Release();
			clone = NULL;
		}
	}
	*userclone = clone;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::AllocateScatteringAngleMatrixElements		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixOrientedWrapper::AllocateScatteringAngleMatrixElements( size_t numangles )
{
	bool	ok;

	ok =    m_S11.SetSize( numangles )
		 && m_S12.SetSize( numangles )		
		 && m_S21.SetSize( numangles )		
		 && m_S22.SetSize( numangles );		
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "sk_TMatrixOrientedWrapper::AllocateScatteringAngleMatrixElements, could not allocate scattering matrix elements");
		m_S11.erase();
		m_S12.erase();
		m_S21.erase();
		m_S22.erase();
	}
	return ok;
}




/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::Set_ShapeCylinder		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixOrientedWrapper::Set_ShapeCylinder()
{
    m_np				= -2;
	return true;
}

/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::Set_ShapeSpheroid		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixOrientedWrapper::Set_ShapeSpheroid()
{
    m_np		   = -1;
	return true;
}


/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::PartTypeString		2009-10-16*/
/** **/
/*---------------------------------------------------------------------------*/

const char* sk_TMatrixOrientedWrapper::PartTypeString() const
{
	const char* str;

	switch (m_np)
	{
	case -1 : str = "spher"; break;
	case -2 : str = "cyl";  break;
	default : str = "Unknown"; break;
	}
	return str;
}


/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::Set_ParticleOrientation		2009-10-15*/
/**Set the orientation of the particle in the lab frame, according to the Euler 
 *	angles 
 *		alpha: rotation of xy-plane about lab z axis
 *		beta: rotation z to z' about line of nodes
 *	Note that this code is implemented for rotationally-symmetric particles,
 *	so no third angle is necessary.
 **/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixOrientedWrapper::Set_ParticleOrientation( double alpha,double beta )
{
	SetDirty( m_alpha != alpha );
	SetDirty( m_beta != beta );
	m_alpha		= alpha;
	m_beta		= beta;
	return true;
}

/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::Set_SunDirection		2009-10-15*/
/** Set the direction of the sun: input solar zenith angle.**/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixOrientedWrapper::Set_SunDirection( double theta0 )
{
	SetDirty( m_thetaInc != theta0 );
	SetDirty( m_phiInc != 0. );
	m_thetaInc	= theta0;
	m_phiInc	= 0.;
	return true;
}


/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::Get_DescriptiveParticleString		2009-10-15*/
/** Return a string that defines a particle according to its parameters, which
 * is used for caching the computed particle's scattering properties
 **/
/*---------------------------------------------------------------------------*/

nxString sk_TMatrixOrientedWrapper::Get_DescriptiveParticleString() const
{
	nxString	basename;

	basename.sprintf("%s_%3.2f_orient_%.3f_%.3f_%.3f",	(const char*)PartTypeString(),  (double)Get_AspectRatio(),  (double)m_alpha, (double)m_beta, (double)m_thetaInc );
	return ( basename ); 
}


/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::CalculateScattering		2009-10-15*/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixOrientedWrapper::CalculateScattering()
{
	bool	ok;

	ok = !IsDirty();
	if (!ok)
	{
		ok = (Get_NumAngles() > 0);
		ok = ok && Mishchenko_TMatrix();
		if (ok) ClearDirty();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::Mishchenko_TMatrix		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixOrientedWrapper::Mishchenko_TMatrix()
{
	bool		ok;
	double		ri_real;
	double		ri_imag;
	double		radius;
	double		aspectratio;
	int			numscatang;
	double		delta;
	double		lambda;
	double		qext;
	double		qsca;
	double		csca;
	double		cext;
	double		area;


	ri_real	    = Get_RefractiveIndex().real();
	ri_imag	    = Get_RefractiveIndex().imag();
	radius      = Get_Radius();
	aspectratio = Get_AspectRatio();
	lambda      = Get_Lambda();
	delta       = Get_Delta();
	numscatang  = (int)Get_NumAngles();

	m_ssalb		= -9999.0; 
    cext		= 0.0;
	csca		= 0.0;
	qext		= 0.0;
	qsca	    = 0.0;


	{
	std::unique_lock<std::mutex>	lock(g_mutexlock);
//	g_mutexlock.Lock();																// grab the mutex, Not sure if Wiscombe is thread safe

	TMATRIXORIENTED(	&radius,								// double*			AXI,		// REAL*8
						&m_rat,									// double*			RAT,		// REAL*8
						&lambda,								// double*			LAM,		// REAL*8
						&ri_real,								// double*			MRR,		// REAL*8
						&ri_imag,								// double*			MRI,		// REAL*8
						&aspectratio,							// double*			EPS,		// REAL*8
						&m_np,									// int*				NP,			// INTEGER
						&delta,									// double*			DDELT,		// REAL*8
						&m_ndgs,								// int*				NDGS,		// INTEGER	
						&m_alpha,								// double*			ALPHA,		// REAL*8
						&m_beta,								// double*			BETA,		// REAL*8
						&m_thetaInc,							// double*			THETA0,		// REAL*8
						&m_thetaSca,							// double*			THET,		// REAL*8
						&m_phiInc,								// double*			PHI0,		// REAL*8
						&m_phiSca,								// double*			PHI,		// REAL*8
						&numscatang,
						 Get_CosAngles().UnsafeArrayBasePtr(),
						 m_S11.UnsafeArrayBasePtr(),
						 m_S12.UnsafeArrayBasePtr(),
						 m_S21.UnsafeArrayBasePtr(),
						 m_S22.UnsafeArrayBasePtr(),
						&qsca,
						&qext,
						&m_ssalb
						);
     }
//	g_mutexlock.Unlock();														// Release the mutex.
	area = nxmath::Pi*radius*radius;
	cext = qext*area;
	csca = qsca*area;
	ok = SetCrossSectionsFromFortran( cext, csca );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::Get_ScatteringMatrixCoeffs		2009-10-15*/
/**	Get the scattering matrix elements (as a function of scattering angle) for 
 *	this nonspherical scattering particle.
 *	Scattering matrix elements (Bohren and Huffman, 1983, eqn 3.16; a.k.a. 
 *	transformation matrix elements: Hansen and Travis, 1974, eqn 2.36; Liou, 
 *	2002, 5.2.105) are set in this routine. These are the following:
 *
 *	 F = 	( F11	F12		F13		F14	)
 *			( F21	F22		F23		F24	)
 *			( F31	F32		F33		F34	)
 *			( F41	F42		F43		F44	)
 *	 (at each value of cos(\Theta) = m_xmu)
 *	Note that the normalization
 *		F(\Theta)*(4*pi/k^2/Csca) = P(\Theta)
 *	to obtain the phase matrix is done within the calling routine of the 
 *	extinction class, skOpticalProperties_Ice::CalculateCrossSections.
 *	The arrays are set up as follows:
 *	coeffmatrix(0,0,*) = F11		*= all angles in m_xmu (see Get_CosAngles())
 *	coeffmatrix(1,0,*) = F21		*= all angles in m_xmu
 *	coeffmatrix(0,1,*) = F12		*= all angles in m_xmu
 *  ...
 *  coeffmatrix(3,2,*) = F43		*= all angles in m_xmu
 *  coeffmatrix(3,3,*) = F44		*= all angles in m_xmu,
 **/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixOrientedWrapper::Get_ScatteringMatrixCoeffs( std::vector<skRTPhaseMatrix>*	coeffmatrix )
{
	bool								ok;
	complex<double>						S11;
	complex<double>						S22;
	complex<double>						S12;
	complex<double>						S21;
	size_t								i;
	size_t								numangles;

	CalculateScattering();								// Update Cache so output parameters reflect latest input parameters
	numangles = Get_NumAngles();
	ok = (numangles > 0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_TMatrixOrientedWrapper::Get_ScatteringMatrixCoeffs, Cannot initialize matrix until scattering angles are defined");
	}
	else
	{

		try
		{

			coeffmatrix->resize(numangles);
			for (i=0; i < numangles; i++)
			{
				S11    = m_S11[i];
				S12    = m_S12[i];
				S21    = m_S21[i];
				S22    = m_S22[i];
				(*coeffmatrix)[i].At(1,1) =  0.5*( norm(S11) + norm(S12) + norm(S21) + norm(S22) );			// Z11
				(*coeffmatrix)[i].At(1,2) =  0.5*( norm(S11) - norm(S12) + norm(S21) - norm(S22) );			// Z12
				(*coeffmatrix)[i].At(1,3) = -real(         S11*conj(S12) + S22*conj(S21)			);		// Z13
				(*coeffmatrix)[i].At(1,4) = -imag(         S11*conj(S12) - S22*conj(S21)			);		// Z13

				(*coeffmatrix)[i].At(2,1) =  0.5*( norm(S11) + norm(S12) - norm(S21) - norm(S22) );			// Z21
				(*coeffmatrix)[i].At(2,2) =  0.5*( norm(S11) - norm(S12) - norm(S21) + norm(S22) );			// Z22
				(*coeffmatrix)[i].At(2,3) = -real(         S11*conj(S12) - S22*conj(S21)			);		// Z23
				(*coeffmatrix)[i].At(2,4) = -imag(         S11*conj(S12) + S22*conj(S21)			);		// Z24

				(*coeffmatrix)[i].At(3,1) = -real(         S11*conj(S21) + S22*conj(S12)			);		// Z31
				(*coeffmatrix)[i].At(3,2) = -real(         S11*conj(S21) - S22*conj(S12)			);		// Z32
				(*coeffmatrix)[i].At(3,3) =  real(         S11*conj(S22) + S12*conj(S21)			);		// Z33
				(*coeffmatrix)[i].At(3,4) =  imag(         S11*conj(S22) + S21*conj(S12)			);		// Z34

				(*coeffmatrix)[i].At(4,1) = -imag(         S21*conj(S11) + S22*conj(S12)			);		// Z41
				(*coeffmatrix)[i].At(4,2) = -imag(         S21*conj(S11) - S22*conj(S12)			);		// Z42
				(*coeffmatrix)[i].At(4,3) =  imag(         S22*conj(S11) - S12*conj(S21)			);		// Z43
				(*coeffmatrix)[i].At(4,4) =  real(         S22*conj(S11) - S12*conj(S21)			);		// Z44
			}
		}
		catch (...)
		{
			nxLog::Record(NXLOG_WARNING,"sk_TMatrixOrientedWrapper::Get_ScatteringMatrixCoeffs, error allocating coefficient matrix ");
			ok = false;
		}

	}
	if (!ok) coeffmatrix->clear();
	return ok;
}
