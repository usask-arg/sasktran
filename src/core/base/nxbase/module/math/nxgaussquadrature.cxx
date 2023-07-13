#include "nxbase_math.h"


/*---------------------------------------------------------------------------
 *'					GAUSST											2003-9-25
 *	THIS SUBROUTINE CALCULATES THE GAUSSIAN QUADRATURE POINTS AND WEIGHTS
 *	FOR AN INTEGRAL FROM X1 TO X2 USING NG POINTS.  This code was ported from
 *	Chris McLindens Fortran Code.
 *
 *	However it does seem to work. I have checked it against various tables
 *	and it does appear to give the same result.
 *
 * PARAMETERS:
 *	NG
 *		TOTAL NUMBER OF GAUSS POINTS INCLUDING BOTH SIDES
 *
 *	X1
 *		START
 *
 *	X2
 *		END
 *
 *	XP	dblarr[NG]
 *		POINTS
 *
 *	WT  dblarr[NG]
 *		WEIGHTS
 *-------------------------------------------------------------------------*/

static void GAUSST( int NG, double X1, double X2, double* XP, double* WT)
{

	const double PS  = 1.013211836423378E-01;
	const double DXL = 1.E-16 ;
	double	XMID = (X2+X1) / 2.0;
	double	XHAF = (X2-X1) / 2.0;
	double	DNG = NG;
	int		NN = NG / 2;
	int		N2 = NN * 2;
	double	PN;
	int		N;
	int		I;
	double	C;
	double	DI;
	double	Z;
	double	ZZ;
	double	DN;
	double	DM;
	double	X;
	double	PNI;
	double	PNJ;
	double	PNK;
	double	DX;
	int		J;

	--XP;								// Place XP and WT on a 1-based
	--WT;								// indexing scheme
	if ( N2 != NG)						// If the number of points is od
	{
		XP[NN+1] = XMID;				// 
		WT[NN+1] = 1.0;
		if (NG < 2)
		{
			WT[NN+1] = 2.0;
			return;
		}
		WT[NN+1] = 1.0;
		PN = 1.0;
		N = 0;

		 do
		{
			N  = N + 2;
			DN = N;
			DM = DN - 1.0;
			PN = PN * (DM / DN);
		}
		while (N < N2);
		WT[NN+1] = 2.0 * XHAF / nxmath::sqr(DNG * PN);
	}

    I  = 0;
	C  = nxmath::Pi/sqrt(DNG * (DNG + 1.0) + 0.50 - PS)/105.0;
	do
	{
	    I  = I + 1;
		DI = I;
		Z  = PS / nxmath::sqr(4.0 * DI - 1.0);
		ZZ = (105.0 + Z * (210.0 - Z * (2170.0 - Z * (105812.0 - 12554474.0 * Z))));
		X = cos(ZZ * C * (DI - 0.25));
	    do
		{
			N = 1;
			DM = 1.0;
			PNI = 1.0;
			PNJ = X;
			do
			{
				N = N + 1;
				DN = N;
				PNK = ((DM + DN) * X * PNJ - DM * PNI) / DN;
				PNI = PNJ;
				PNJ = PNK;
				DM = DN;
			} while (N < NG);

			DX = PNJ * (1.0 - X * X) / DNG / (PNI - X * PNJ);
			X  = X - DX;
		} while ( fabs(DX) > DXL);

		J     = NG + 1 - I;
		XP[I] = XMID - XHAF * X;
		XP[J] = XMID + XHAF * X;
		WT[I] = 2.0 * XHAF * (1.0 - X * X) / nxmath::sqr(DNG * PNI);
		WT[J] = WT[I];
	} while ( I < NN );
}

/*-----------------------------------------------------------------------------
 *					nxGaussQuadratureBase::nxGaussQuadratureBase	2003-9-26*/
/** Default constructor**/
/*---------------------------------------------------------------------------*/

nxGaussQuadratureBase::nxGaussQuadratureBase()
{
	m_N  = -1;				// The "order" of the quadrature.
	m_x1 = -1;				// The value of X at the start of the interval
	m_x2 =  1;				// The value of X at the end of the interval
	m_isdirty = nxFALSE;
}



/*-----------------------------------------------------------------------------
 *					nxGaussQuadratureBase::operator=		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxGaussQuadratureBase::DeepCopy( const nxGaussQuadratureBase& other )
{
	bool		ok;
	ok =              m_x.DeepCopy      (other.m_x);
	ok        = ok && m_weights.DeepCopy(other.m_weights);
	m_N       = other.m_N;
	m_x1      = other.m_x1;				// The value of X at the start of the interval
	m_x2      = other.m_x2;				// The value of X at the end of the interval
	m_isdirty = other.m_isdirty;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxGaussQuadratureBase::ConfigureMultiRange		2004-3-1 */
/** This allows the gauss quadrature to be broken into \b N individual intervals.
 *	Each interval is actually treated as its own seperate gaussian quadrature
 *	interval. This was initially developed to break the zenith quadrature integration in
 *	plane parallel radiative transfer into a downeard and upward hemisphere as
 *	there is a discontinuity at the 90 degree zenith angle between the upward
 *	and downward fluxes.  A single range of integration can be set by calling
 *	#SetRange
 *
 *	\param startx
 *		 The starting x coordinate of each interval.\c double[numranges]
 *
 *	\param 	endx  
 *			The ending x coordinate of each interval, \c double[numranges]
 *
 *	\param 	ptsperrange
 *			The number of quadrature points in each interval. \c int[numranges]
 *
 *	\param 	numranges
 *		The number of intervals.
 *
 *	\return
 *	Returns nxTRUE if success.
**/
/*---------------------------------------------------------------------------*/

nxBOOL nxGaussQuadratureBase::ConfigureMultiRange( double* startx, double* endx, int* ptsperrange, int numranges)
{
	int		i;
	int		N;
	int		offset;
	double*	xptr;
	double*	wptr;
	nxBOOL	ok;

	N = 0;
	for (i=0; i < numranges; i++)
	{
		NXASSERT(ptsperrange[i] > 0);
		N += ptsperrange[i];
	}
	m_N = N;
	ok  = m_x.SetSize(m_N);							// The location of X coordinates
	ok  = ok && m_weights.SetSize(m_N);				// The weights for each coordinate
	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR, "nxGaussQuadratureBase::ConfigureMultiRange, memory allocation error");
		m_x.erase();
		m_weights.erase();
		m_isdirty = nxFALSE;
		m_N       = 0;

	}
	else
	{
		m_x1 = startx[0];
		m_x2 = endx[numranges-1];
		xptr = m_x.UnsafeArrayBasePtr();
		wptr =  m_weights.UnsafeArrayBasePtr();
		offset = 0;
		for (i=0; i < numranges; i++)
		{
			GAUSST(ptsperrange[i], startx[i], endx[i], xptr+offset, wptr + offset);
			offset += ptsperrange[i];
		}
		m_isdirty = nxFALSE;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxGaussQuadratureBase::CheckDirty				2003-9-26*/
/** Check to see if the configuration is dirty. If it is then calculate a new
 *	set of quadrature points.
 **/
/*---------------------------------------------------------------------------*/

void nxGaussQuadratureBase::CheckDirty()
{
	NXASSERT(m_N > 0);
	if (m_isdirty)
	{
		m_isdirty = !ConfigureMultiRange( &m_x1, &m_x2, &m_N, 1);
	}
}


/*---------------------------------------------------------------------------
 *					nxGaussQuadratureBase::SetRange				2003-9-26 */
/** Sets the range of integration from X1 to X2.  
**/
/*-------------------------------------------------------------------------*/

void nxGaussQuadratureBase::SetRange( double x1, double x2 )
{
	NXASSERT( x1 < x2);
	m_isdirty = m_isdirty || (x1 != m_x1) || (x2 != m_x2);
	m_x1 = x1;
	m_x2 = x2;
}

/*-----------------------------------------------------------------------------
 *					nxGaussQuadratureBase::SetOrder					2003-9-26*/
/** Sets the order of integration to the specified value **/
/*---------------------------------------------------------------------------*/

void nxGaussQuadratureBase::SetOrder( int N )
{
	NXASSERT( N > 0 );
	m_isdirty = m_isdirty || (N != m_N);
	m_N = N;
}

/*-----------------------------------------------------------------------------
 *					nxGaussQuadratureBase::Integrate		2003-9-26		 */
/** Integrates tabulated versions of Y.  The values of Y must correspond to
 *	the function eveluated at the quadrature points specified in #X.
 **/
/*---------------------------------------------------------------------------*/

double nxGaussQuadratureBase::Integrate( nx1dArray<double>& Y )
{
	nxArrayIter<double>	yi = Y.begin();
	nxArrayIter<double>	ye = Y.end();
	nxArrayIter<double> wi;
	double				sum;
	nxBOOL				ok;
	
	CheckDirty();
	ok = ((int)Y.size() == m_N);
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "nxGaussQuadratureBase::Integrate, The size of Y (%d) does not match the quadrature order (%d), the integral is returned as 1.0E300", (int)Y.size(), (int)m_N);
		sum = 1.0E300;
	}
	else
	{
		sum = 0.0;
		wi  = m_weights.begin();
		while( yi != ye)
		{
			sum += (*wi)*(*yi);
			++wi;
			++yi;
		}
	}
	return sum;
}


/*-----------------------------------------------------------------------------
 *					nxTrapezoidalQuadratureBase::nxTrapezoidalQuadratureBase		2004-8-8*/
/** **/
/*---------------------------------------------------------------------------*/

nxTrapezoidalQuadratureBase::nxTrapezoidalQuadratureBase()
{
	m_N = 0;									// The "order" of the quadrature.
	m_x1 = 0.0;									// The value of X at the start of the interval
	m_x2 = 0.0;									// The value of X at the end of the interval
}


/*-----------------------------------------------------------------------------
 *					nxTrapezoidalQuadratureBase::SetRange		2004-8-8*/
/** Set the range of integration from x1 to x2
**/
/*---------------------------------------------------------------------------*/

void nxTrapezoidalQuadratureBase::SetRange	( double x1, double x2 )
{
	m_x1 = x1;
	m_x2 = x2;
}


/*-----------------------------------------------------------------------------
 *					nxTrapezoidalQuadratureBase::SetOrder		2004-8-8*/
/** Set the number of quadrature points in the integration.  The number of
 *	must be at least 2. It is implicitly assumed that the points are evenly
 *	spaced between the start and end of integration. 
 **/
/*---------------------------------------------------------------------------*/

void  nxTrapezoidalQuadratureBase::SetOrder( int N )
{
	if (N < 2) N = 2;
	m_N = N;
}


/*-----------------------------------------------------------------------------
 *					nxTrapezoidalQuadratureBase::Integrate		2004-8-8*/
/** Integrate the array Y. It is implicitly assumed that Y is evaluated at the
 *	evenly spaced quadrature points.  There must eb at least 2 points in the
 *	quadrature (evaluated at the start and end of integration)
 **/
/*---------------------------------------------------------------------------*/

double nxTrapezoidalQuadratureBase::Integrate( nx1dArray<double>& Y )
{
	double				h;
	nxArrayIter<double>	yi = Y.begin();
	nxArrayIter<double>	ye = Y.end();
	double				sum;
	nxBOOL				ok;
	int					N1;
	
	ok = ((int)Y.size() == m_N) && ( m_N >= 2);
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "nxTrapezoidalQuadratureBase::Integrate, The size of Y (%d) does not match the quadrature order (%d), the integral is returned as 1.0E300", (int)Y.size(), (int)m_N);
		sum = 0;
	}
	else
	{
		N1 = m_N-1;
		h   = (m_x2-m_x1)/N1;
		sum = 0.5*(*yi);
		++yi;
		for (int i = 1; i < N1; i++)
		{
			sum += *yi;
			++yi;
		}
		sum += 0.5*(*yi);
		sum *= h;
	}
	return sum;
}
