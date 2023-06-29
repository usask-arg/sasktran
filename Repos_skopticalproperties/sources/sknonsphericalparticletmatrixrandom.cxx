#include <skopticalproperties21.h>
#include <boost/thread.hpp>
//using namespace std;


static std::mutex	g_mutexlock;

#if !defined(NX_WINDOWS)
#define TMATRIXRANDOMFUNC tmatrixrandom_
#else
#define TMATRIXRANDOMFUNC TMATRIXRANDOMEP
#endif

   extern "C" void  TMATRIXRANDOMFUNC(	double*			REFF,		// REAL*8
									double*			VEFF,		// REAL*8
									double*			RAT,		// REAL*8
									int*			NDISTR,		// INTEGER
									double*			AXMAX,		// REAL*8
									double*			R1,			// REAL*8
									double*			R2,			// REAL*8
									int*			NPNAX,		// INTEGER
									double*			B,			// REAL*8
									double*			GAM,		// REAL*8
									int*			NKMAX,		// INTEGER
									double*			EPS,		// REAL*8
									int*			NP,			// INTEGER                                             
									double*			LAM,		// REAL*8
									double*			MRR,		// REAL*8
									double*			MRI,		// REAL*8
									double*			DDELT,		// REAL*8
									int*			NPNA,		// INTEGER
									int*			NDGS,		// INTEGER
									double*			CEXT,		// REAL*8
									double*			CSCA,		// REAL*8
									double*			W,			// REAL*8
									double*			ASYMM,		// REAL*8
									double*			XMU,		// REAL*8(0:NPNA-1)
									double*			F11,		// REAL*8(0:NPNA-1)
									double*			F22,		// REAL*8(0:NPNA-1)
									double*			F33,		// REAL*8(0:NPNA-1)
									double*			F44,		// REAL*8(0:NPNA-1)
									double*			F12,		// REAL*8(0:NPNA-1)
									double*			F34         // REAL*8(0:NPNA-1)
									);		


/*---------------------------------------------------------------------------
 *'			sk_TMatrixRandomWrapper::sk_TMatrixRandomWrapper		2009-6-3
 *-------------------------------------------------------------------------*/

sk_TMatrixRandomWrapper::sk_TMatrixRandomWrapper()
{
	m_rat			= 1.0;							// equal-volume sphere by default
	m_np			= -2;							// 
	m_ndgs			= 4;
	Set_AspectRatio ( 1.0 );
	Set_ShapeCylinder( );
	SetDefaultAngles();
	sk_NonsphericalParticle::SetDirty( true );
}


/*---------------------------------------------------------------------------
 *'		sk_TMatrixRandomWrapper::~sk_TMatrixRandomWrapper		2009-6-3
 *-------------------------------------------------------------------------*/

sk_TMatrixRandomWrapper::~sk_TMatrixRandomWrapper()
{
	m_f11.erase();			
	m_f22.erase();			
	m_f33.erase();			
	m_f44.erase();			
	m_f12.erase();			
	m_f34.erase();			
}

/*-----------------------------------------------------------------------------
 *					sk_NonsphericalParticle::SetDefaultAngles		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixRandomWrapper::SetDefaultAngles()
{
	double	startangle[2] = {0.0,  10.0};
	double  endangle[2]   = {10.0, 90.0};
	double  resolution[2] = {0.1,  0.25};
    bool ok;

	ok =  Set_SymmetricScatterAngles( startangle, endangle, resolution, 2 );
	return ok;
}

/*-----------------------------------------------------------------------------
 *				sk_TMatrixRandomWrapper::operator=			2009-6-3
 *---------------------------------------------------------------------------*/

bool sk_TMatrixRandomWrapper::DeepCopy( const sk_TMatrixRandomWrapper& other )
{
	bool	ok;

	ok              = sk_NonsphericalParticleDistributedSizes::DeepCopy( other );
	m_np			= other.m_np;			
	m_ndgs			= other.m_ndgs;			
	m_rat			= other.m_rat;			
	ok				= ok && m_f11.DeepCopy( other.m_f11 );
	ok				= ok && m_f22.DeepCopy( other.m_f22 );
	ok				= ok && m_f33.DeepCopy( other.m_f33 );
	ok				= ok && m_f44.DeepCopy( other.m_f44 );
	ok				= ok && m_f12.DeepCopy( other.m_f12 );
	ok				= ok && m_f34.DeepCopy( other.m_f34 );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_TMatrixRandomWrapper::DeepCopy, Error copying contents, this object is in an ill-defined state and may produce strange results");
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'		sk_TMatrixRandomWrapper::CreateClone				2009-6-3
 * 
 *-------------------------------------------------------------------------*/
bool sk_TMatrixRandomWrapper::CreateClone(sk_NonsphericalParticle** userclone) const
{
	sk_TMatrixRandomWrapper*	clone;
	bool						ok;

	clone = new sk_TMatrixRandomWrapper;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sk_TMatrixRandomWrapper::CreateClone, Error creating clone object");
	}
	else
	{
		clone->AddRef();
	}
	*userclone = clone;
	return ok;
}


/*---------------------------------------------------------------------------
 *'		sk_NonsphericalParticleDistributedSizes::sk_NonsphericalParticleDistributedSizes		2009-6-8
 *-------------------------------------------------------------------------*/
sk_NonsphericalParticleDistributedSizes::sk_NonsphericalParticleDistributedSizes( )
{
	m_b			= -9999;						// effective width
	m_npnax		= 1;							// one size distribution
	m_ndistr	= 2;							// the lognormal
	m_gam		= 0.5e0;						// IGNORED
	m_nkmax		= 48;							// NKMAX+2: num points over size distribution
	m_r1		= -9998;
	m_r2		= -9999;
	m_ssalb     = -9999;
	m_asym      = -9999;
}

/*---------------------------------------------------------------------------
 *'				sk_NonsphericalParticleDistributedSizes::DeepCopy				2009-6-8
 *-------------------------------------------------------------------------*/
bool sk_NonsphericalParticleDistributedSizes::DeepCopy( const sk_NonsphericalParticleDistributedSizes& other )
{
	bool	ok;

	ok = sk_NonsphericalParticle::DeepCopy( other );
	m_r1			= other.m_r1;
	m_r2			= other.m_r2;
	m_b				= other.m_b;			
	m_gam			= other.m_gam;			
	m_npnax			= other.m_npnax;		
	m_ndistr		= other.m_ndistr;		
	m_nkmax			= other.m_nkmax;	
	m_ssalb		    = other.m_ssalb;
	m_asym          = other.m_asym;	
	return ok;
}


/*---------------------------------------------------------------------------
 *'		sk_NonsphericalParticleDistributedSizes::Set_SizeDistMonodisperse		2009-6-8
 * 
 *-------------------------------------------------------------------------*/
void sk_NonsphericalParticleDistributedSizes::Set_SizeDistMonodisperse( double radius )
{
	Set_SizeDistGamma(  radius,1e-1 );
	m_npnax		= 1;
	Set_Radius(radius);
	m_nkmax		= -1;
	m_r1 = 0.9999999*radius;
	m_r2 = 1.0000001*radius;
}

/*---------------------------------------------------------------------------
 *'		sk_NonsphericalParticleDistributedSizes::Set_NumPartSizePoints		2009-6-8
 * 
 *-------------------------------------------------------------------------*/
void sk_NonsphericalParticleDistributedSizes::Set_NumPartSizePoints( int num_pts )
{
	m_nkmax	= num_pts;
}

/*---------------------------------------------------------------------------
 *'		sk_NonsphericalParticleDistributedSizes::Set_SizeDistGamma			2009-6-8
 * 
 *-------------------------------------------------------------------------*/
//void sk_NonsphericalParticleDistributedSizes::Set_SizeDistGamma( double radius,double width )
void sk_NonsphericalParticleDistributedSizes::Set_SizeDistGamma( double lambda,double mu )
{
	double		scale,shape;
	double		a,b;

	shape	= mu + 1;
	scale	= 1/lambda;
	a		= (mu + 3)/lambda;
	b		= 1/(mu+3);
	Set_Radius(a);
	m_b				= b;

	m_ndistr		= 4;		
	//Set_Radius(radius);
	//m_b				= width;
	m_gam			= 0;
	Set_SizeDistIntegrationLimits();
}

/*---------------------------------------------------------------------------
 *'		sk_NonsphericalParticleDistributedSizes::Set_SizeDistModGamma		2009-6-8
 * 
 *-------------------------------------------------------------------------*/
void sk_NonsphericalParticleDistributedSizes::Set_SizeDistModGamma( double /*rc*/, double /*alpha*/, double /*gamma*/ )
{
//	Set_SizeDistIntegrationLimits( rc,alpha,gamma );
}

/*---------------------------------------------------------------------------
 *'		sk_NonsphericalParticleDistributedSizes::Set_SizeDistLogNormal			2009-6-8
 * 
 *	\param radius
 *	The mode radius in microns
 *
 *	\param sigma
 *	The natural log of the mode width
 *-------------------------------------------------------------------------*/
void sk_NonsphericalParticleDistributedSizes::Set_SizeDistLogNormal( double radius,double sigma )
{
	bool	ischanged;
	double	newb;

	newb      = sigma*sigma;
	ischanged = ( m_b != newb);
	m_ndistr  = 2;				
	m_b		  = newb;								// Mischenko requires square of log(mode width)
	Set_Radius(radius);								// Mischenko requires mode radius in microns					
	SetDirty(ischanged);
	Set_SizeDistIntegrationLimits( );
}

/*---------------------------------------------------------------------------
 *'		sk_NonsphericalParticleDistributedSizes::Set_SizeDistPowerLaw			2009-6-8
 * 
 *-------------------------------------------------------------------------*/
void sk_NonsphericalParticleDistributedSizes::Set_SizeDistPowerLaw( double radius, double width )
{
	m_ndistr	= 3;							// the power law dist.
	m_b			= width;
	m_gam		= 0;
	Set_Radius(radius);
	Set_SizeDistIntegrationLimits( );
}

/*---------------------------------------------------------------------------
 *'		sk_NonsphericalParticleDistributedSizes::Set_SizeDistIntegrationLimits		2009-7-17
 * 
 *-------------------------------------------------------------------------*/
bool sk_NonsphericalParticleDistributedSizes::Set_SizeDistIntegrationLimits()
{
	double		widthNum;
	double		w;
	double		radius = Get_Radius();
	double		factor;

	switch(m_ndistr)
	{
		case 1:			// mod gamma
			nxLog::Record( NXLOG_ERROR, "sk_TMatrixRandomWrapper::Set_SizeDistIntegrationLimits, Not yet implemented for modified gamma distribution");
			break;
		case 2:			// lognormal
			widthNum	= 4.0;
			w           = exp( sqrt(m_b)  );	// Get the mode width from the square of log(mode radius)
			factor      = pow( w, widthNum );
			m_r1 = radius/factor;	
			m_r2 = radius*factor;
			break;
		case 3:			// power-law
			nxLog::Record( NXLOG_ERROR, "sk_TMatrixRandomWrapper::Set_SizeDistIntegrationLimits, Not yet implemented for power law distribution");
			break;
		case 4:			// gamma
			m_r1 = radius/20;
			m_r2 = radius*5.0;
			break;
	}
	SetDirty(true); 
	return true;
}

/*---------------------------------------------------------------------------
 *'		sk_TMatrixRandomWrapper::Set_ParticleDistribution		2009-6-24
 *	Make newdistribution->Type() trigger an update in the necessary parameters
 *	to make the call to DLL TMATRIXRANDOM.dll.
 *-------------------------------------------------------------------------*/

bool sk_TMatrixRandomWrapper::Set_ParticleDistribution( skRTParticleDist* newdistribution )
{
	bool	ok;
	double	A,B,C;
	double	moderadius;
	double	modewidth;
	ok = true;

	switch( newdistribution->Type() )
	{
		case skRTParticleDist::PDIST_2GAMMA:
			newdistribution->GetDistributionParameters( &A,&B,&C );
			Set_SizeDistGamma( A,B );
			break;

		case skRTParticleDist::PDIST_3GAMMA:	//not properly implemented
			newdistribution->GetDistributionParameters( &A,&B,&C );
			Set_SizeDistModGamma( A,B,C );
			break;

		case skRTParticleDist::PDIST_LOGNORMAL:
			newdistribution->GetDistributionParameters( &moderadius, &modewidth, &C );
			Set_SizeDistLogNormal( moderadius, log(modewidth) );
			break;

		case skRTParticleDist::PDIST_POWERLAW:	// not properly implemented
			newdistribution->GetDistributionParameters( &A,&B,&C );
			Set_SizeDistPowerLaw( A,B );
			break;

		case skRTParticleDist::PDIST_MONODISPERSE:	
			newdistribution->GetDistributionParameters( &A,&B,&C );
			Set_SizeDistMonodisperse( A );
			break;

		default :
			nxLog::Record(NXLOG_WARNING, "sk_TMatrixRandomWrapper::Set_ParticleDistribution, We do not currently support the requested particle distrubtion");
			ok = false;
			break;
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *'		sk_TMatrixRandomWrapper::Set_ShapeCylinder			2009-7-9
 * 
 *-------------------------------------------------------------------------*/
bool sk_TMatrixRandomWrapper::Set_ShapeCylinder()
{
    m_np		  = -2;
	return true;
}

/*---------------------------------------------------------------------------
 *'		sk_TMatrixRandomWrapper::Set_ShapeSpheroid			2009-7-9
 * 
 *-------------------------------------------------------------------------*/
bool sk_TMatrixRandomWrapper::Set_ShapeSpheroid()
{
    m_np				= -1;
	return true;
}

/*-----------------------------------------------------------------------------
 *					sk_TMatrixOrientedWrapper::PartTypeString		2009-10-16*/
/** **/
/*---------------------------------------------------------------------------*/

const char* sk_TMatrixRandomWrapper::PartTypeString() const
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

/*---------------------------------------------------------------------------
 *'		sk_TMatrixRandomWrapper::Get_DescriptiveParticleString		2009-6-3
 * Return a string that defines a particle according to its parameters, which
 * is used for caching the computed particle's scattering properties
 *-------------------------------------------------------------------------*/
nxString sk_TMatrixRandomWrapper::Get_DescriptiveParticleString() const
{
	nxString	basename;

	basename.sprintf("%s_%3.2f_rand",	(const char*)PartTypeString(),  (double)Get_AspectRatio() );
	return ( basename ); 
}


/*-----------------------------------------------------------------------------
 *			sk_TMatrixRandomWrapper::AllocateScatteringAngleMatrixElements		2009-10-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool sk_TMatrixRandomWrapper::AllocateScatteringAngleMatrixElements(size_t numangles)
{
	bool	ok;

	ok =      m_f11.SetSize( numangles );		
	ok = ok & m_f22.SetSize( numangles );		
	ok = ok & m_f33.SetSize( numangles );		
	ok = ok & m_f44.SetSize( numangles );		
	ok = ok & m_f12.SetSize( numangles );	
	ok = ok & m_f34.SetSize( numangles );		

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "sk_TMatrixRandomWrapper::AllocateScatteringMatrixElements, could not allocate scattering matrix elements");
		m_f11.erase();
		m_f22.erase();
		m_f33.erase();
		m_f44.erase();
		m_f12.erase();
		m_f34.erase();
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *'				sk_TMatrixRandomWrapper::CalculateScattering		2009-6-3
 *-------------------------------------------------------------------------*/

bool sk_TMatrixRandomWrapper::CalculateScattering()
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

/*---------------------------------------------------------------------------
 *'				sk_TMatrixRandomWrapper::Mishchenko_TMatrix			2009-6-3
 *-------------------------------------------------------------------------*/

bool sk_TMatrixRandomWrapper::Mishchenko_TMatrix()
{
	bool		ok = true;
	double		reff;
	double		veff;
	double		ri_real;
	double		ri_imag;
	double		radius;
	double		aspectratio;
	int			numscatang;
	double		delta;
	double		lambda;
	double		cext;
	double		csca;
	double		qext;
	double		qsca;
	double		ssalb;
	double		asym;


	ri_real	    = Get_RefractiveIndex().real();
	ri_imag	    = Get_RefractiveIndex().imag();
	radius      = Get_Radius();
	aspectratio = Get_AspectRatio();
	lambda      = Get_Lambda();
	delta       = Get_Delta();
	numscatang  = (int)Get_NumAngles();

	ssalb		= 0.0; //Get_SSALB();
	asym        = 0.0; //Get_ASYMP();
    cext		= 0.0;
	csca		= 0.0;
	qext		= 0.0;
	qsca	    = 0.0;


	{
	std::unique_lock<std::mutex>	lock(g_mutexlock);
	TMATRIXRANDOMFUNC(	&reff,									// double*			REFF		// REAL*8	
					&veff,									// double*			VEFF		// REAL*8	
 					&m_rat,									// double*			RAT,		// INPUT REAL*8
					&m_ndistr,								// int*				NDISTR,		// INPUT INTEGER	
					&radius,								// double*			AXMAX,		// INPUT REAL*8   particle radius
					&m_r1,									// double*			R1,			// INPUT REAL*8   min radius: size dist int.
					&m_r2,									// double*			R2,			// INPUT REAL*8   max radius: size dist int.
					&m_npnax,								// int*				NPNAX,		// INPUT INTEGER	number of size distributions of spec. type
					&m_b,									// double*			B,			// INPUT REAL*8
					&m_gam,									// double*			GAM,		// INPUT REAL*8
					&m_nkmax,								// int*				NKMAX,		// INPUT INTEGER
					&aspectratio,							// double*			EPS,		// REAL*8
					&m_np,									// int*				NP,			// INTEGER
					&lambda,								// double*			LAM,		// REAL*8
					&ri_real,								// double*			MRR,		// REAL*8
					&ri_imag,								// double*			MRI,		// REAL*8
					&delta,									// double*			DDELT,		// REAL*8
					&numscatang,							// int*				NPNA,		// INTEGER
					&m_ndgs,								// int*				NDGS,		// INTEGER	
					&cext,									// double*			CSCA,		// REAL*8
					&csca,									// double*			CEXT		// REAL*8
					&ssalb,									// double*			W,			// REAL*8
					&asym,									// double*			ASYMM,		// REAL*8
					Get_CosAngles().UnsafeArrayBasePtr(),	// double*			XMU,		// REAL*8(1..NPNA)
					m_f11.UnsafeArrayBasePtr(),				// double*			F11,		// REAL*8(1..NPNA)
					m_f22.UnsafeArrayBasePtr(),				// double*			F22,		// REAL*8(1..NPNA)
					m_f33.UnsafeArrayBasePtr(),				// double*			F33,		// REAL*8(1..NPNA)
					m_f44.UnsafeArrayBasePtr(),				// double*			F44,		// REAL*8(1..NPNA)
					m_f12.UnsafeArrayBasePtr(),				// double*			F12,		// REAL*8(1..NPNA)
					m_f34.UnsafeArrayBasePtr()				// double*			F34,		// REAL*8(1..NPNA)
					);

    }
	if (abs(csca+1e6)<1e-3)												// Code returns -1.0E6 if there is a problem
	{
		nxLog::Verbose(NXLOG_WARNING, "sk_TMatrixRandomWrapper::Mishchenko_TMatrix, computations for reff=%.2d did not converge.",(double)reff );
		ok	= false;
	}
	else
	{
		Set_SSALB( ssalb);
		Set_ASYM(  asym );
		SetCrossSectionsFromFortran	( cext, csca );
	}

	return ok;
}





/*---------------------------------------------------------------------------
 *'				sk_TMatrixRandomWrapper::Get_ScatteringMatrixCeoffs		2009-6-3
 *	Get the scattering matrix elements (as a function of scattering angle) for 
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
 *	to obtain the phase matrix is already done within the DLL code, so there is 
 *	no need to perform this within the extinction class calling routine, 
 *	skOpticalProperties_Ice::CalculateCrossSections.
 *
 *	The arrays are set up as follows:
 *		coeffmatrix(0,*) = a1		* = all angles in m_xmu (see Get_CosAngles())
 *		coeffmatrix(1,*) = a2		
 *		coeffmatrix(2,*) = a3		
 *		coeffmatrix(3,*) = a4		
 *		coeffmatrix(4,*) = b1		
 *		coeffmatrix(5,*) = b2		
 *	where the elements of the normalized scattering matrix (Mishchenko, eqn 4.51; 
 *	a.k.a. phase matrix) are set in this routine.  These are the following:
 *
 *		( a1(mu)	b1(mu)  0		0		) 
 *		( b1(mu)	a2(mu)  0		0		)
 *		(  0		0		a3(mu)	b2(mu)	)
 *		(  0		0		-b2(mu)	a4(mu)	)
 *	which is block-diagonal due to the random orientation of the axissymmetric 
 *	scattering particles.
 *	The normalization of the phase matrix elements is performed such that
 *				\int_{0}^{4\pi}d\Omega[a1(\theta)] = 4\pi
 *-------------------------------------------------------------------------*/

bool sk_TMatrixRandomWrapper::Get_ScatteringMatrixCoeffs( std::vector<skRTPhaseMatrix>*	coeffmatrix )
{
	bool			ok;
	size_t			i;
	size_t			numangles;
	double			a1,a2,a3,a4;
	double			b1,b2;


	numangles = Get_NumAngles();
	ok = CalculateScattering();					// Update Cache so output parameters reflect latest input parameters
	ok = ok && (numangles > 0);
	if (ok)
	{
		try
		{
			coeffmatrix->resize(numangles);
			for (i=0; i < numangles; i++)
			{
				a1 = m_f11.At(i);			// a1
				a2 = m_f22.At(i);			// a2
				a3 = m_f33.At(i);			// a3
				a4 = m_f44.At(i);			// a4
				b1 = m_f12.At(i);			// b1
				b2 = m_f34.At(i);			// b2

				(*coeffmatrix)[i].At(1,1) = a1;
				(*coeffmatrix)[i].At(1,2) = b1;
				(*coeffmatrix)[i].At(1,3) = 0;
				(*coeffmatrix)[i].At(1,4) = 0;

				(*coeffmatrix)[i].At(2,1) = b1;
				(*coeffmatrix)[i].At(2,2) = a2;
				(*coeffmatrix)[i].At(2,3) = 0;
				(*coeffmatrix)[i].At(2,4) = 0;

				(*coeffmatrix)[i].At(3,1) = 0;
				(*coeffmatrix)[i].At(3,2) = 0;
				(*coeffmatrix)[i].At(3,3) = a3;
				(*coeffmatrix)[i].At(3,4) = b2;

				(*coeffmatrix)[i].At(4,1) =  0;
				(*coeffmatrix)[i].At(4,2) =  0;
				(*coeffmatrix)[i].At(4,3) = -b2;
				(*coeffmatrix)[i].At(4,4) =  a4;
			}
		}
		catch (...)
		{
			nxLog::Record(NXLOG_WARNING,"sk_TMatrixRandomWrapper::Get_ScatteringMatrixCoeffs, Error allocating phase matrix array. Thats a problem");
			ok = false;
		}
	}

	return ok;
}

///*---------------------------------------------------------------------------
// *'				sk_TMatrixRandomWrapper::Set_UncommonParameters		2009-10-9
// *	Function to allow user to directly set parameters directly (testing purposes
// *	only).  These members are typically set by standard interface functions.
// *-------------------------------------------------------------------------*/
//
//bool sk_TMatrixRandomWrapper::Set_UncommonParameters( double RAT,double R1,double R2,int NP,int NDISTR,int NPNAX, int NKMAX, int NDGS )
//{
//	sk_NonsphericalParticleDistributedSizes::m_r1		= R1;
//	sk_NonsphericalParticleDistributedSizes::m_r2		= R2;
//	sk_NonsphericalParticleDistributedSizes::m_npnax	= NPNAX;
//	sk_NonsphericalParticleDistributedSizes::m_ndistr	= NDISTR;
//	sk_NonsphericalParticleDistributedSizes::m_nkmax	= NKMAX;
//	m_rat				= RAT;	
//	m_np				= NP;		
//	m_ndgs				= NDGS;		
//	SetDirty(  sk_NonsphericalParticleDistributedSizes::m_r1	!= R1
//           && sk_NonsphericalParticleDistributedSizes::m_r2		!= R2
//           && sk_NonsphericalParticleDistributedSizes::m_npnax	!= NPNAX
//           && sk_NonsphericalParticleDistributedSizes::m_ndistr	!= NDISTR
//           && sk_NonsphericalParticleDistributedSizes::m_nkmax	!= NKMAX
//           && m_rat != RAT && m_np != NP && m_ndgs != NDGS ); 
//	return true;
//}

