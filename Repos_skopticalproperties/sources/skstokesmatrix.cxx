#include <skopticalproperties21.h>

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::skRTStokesVector		2003-10-22
 *-------------------------------------------------------------------------*/

skRTStokesVector::skRTStokesVector()
{
	NXASSERT( sizeof(skRTStokesVector) == 4*sizeof(m_phasematrixstorage[0]) );			// **** Some code relies upon fact its exactly 4 doubles 
	SetTo(0.0);
}


/*---------------------------------------------------------------------------
 *'					:skRTStokesVector::skRTStokesVector		2003-12-18
 *-------------------------------------------------------------------------*/

skRTStokesVector::skRTStokesVector( const skRTStokesVector& other )
{
	*this = other;
}		

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::operator =		2003-12-18
 *-------------------------------------------------------------------------*/

skRTStokesVector&	skRTStokesVector::operator = ( const skRTStokesVector& other )
{
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage);i++) m_phasematrixstorage[i] = other.m_phasematrixstorage[i];
	return *this;
}

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::operator*=		2003-12-18
 *-------------------------------------------------------------------------*/

skRTStokesVector&	skRTStokesVector::operator*= ( float val )
{
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage);i++) m_phasematrixstorage[i] *= (SKRTFLOAT)val;
	return *this;
}

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::operator*=		2003-12-18
 *-------------------------------------------------------------------------*/

skRTStokesVector&	skRTStokesVector::operator*= ( double val )
{
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage);i++) m_phasematrixstorage[i] *= (SKRTFLOAT)val;
	return *this;
}

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::operator +=		2004-2-5
 *-------------------------------------------------------------------------*/

skRTStokesVector&	skRTStokesVector::operator += (const skRTStokesVector& other)
{
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage);i++) m_phasematrixstorage[i] += other.m_phasematrixstorage[i];
	return *this;
}

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::operator +		2004-2-5
 *-------------------------------------------------------------------------*/

skRTStokesVector	skRTStokesVector::operator + (const skRTStokesVector& other) const
{
	skRTStokesVector	answer;
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage);i++)
	{
		answer.m_phasematrixstorage[i] = m_phasematrixstorage[i] + other.m_phasematrixstorage[i];
	}
	return answer;
}
/*---------------------------------------------------------------------------
 *'					skRTStokesVector::operator -		2004-2-5
 *-------------------------------------------------------------------------*/

skRTStokesVector	skRTStokesVector::operator - (const skRTStokesVector& other) const
{
	skRTStokesVector	answer;
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage);i++)
	{
		answer.m_phasematrixstorage[i] = m_phasematrixstorage[i] - other.m_phasematrixstorage[i];
	}
	return answer;
}
/*---------------------------------------------------------------------------
 *'					skRTStokesVector::operator -		2004-2-5
 *-------------------------------------------------------------------------*/

skRTStokesVector	skRTStokesVector::operator * (double value) const
{
	skRTStokesVector	answer;
	SKRTFLOAT			val = (SKRTFLOAT)value;

	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage);i++)
	{
		answer.m_phasematrixstorage[i] = m_phasematrixstorage[i]*val;
	}
	return answer;
}

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::operator -		2004-2-5
 *-------------------------------------------------------------------------*/

skRTStokesVector	skRTStokesVector::operator * (float value) const
{
	skRTStokesVector	answer;
	SKRTFLOAT			val = (SKRTFLOAT)value;

	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage);i++)
	{
		answer.m_phasematrixstorage[i] = m_phasematrixstorage[i]*val;
	}
	return answer;
}

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::SetTo		2004-1-30
 *-------------------------------------------------------------------------*/

void skRTStokesVector::SetTo ( float val )
{
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage);i++) m_phasematrixstorage[i] = (SKRTFLOAT)val;
}

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::SetTo		2004-1-30
 *-------------------------------------------------------------------------*/

void skRTStokesVector::SetTo ( double val )
{
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage);i++) m_phasematrixstorage[i] = (SKRTFLOAT)val;
}

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::SetTo		2004-1-30
 *-------------------------------------------------------------------------*/

void skRTStokesVector::SetTo (  double I, double Q, double U, double V )
{
	At(1) = (SKRTFLOAT)I;
	At(2) = (SKRTFLOAT)Q;
	At(3) = (SKRTFLOAT)U;
	At(4) = (SKRTFLOAT)V;
}

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::SetTo		2004-1-30
 *-------------------------------------------------------------------------*/

void skRTStokesVector::SetTo (  double* IQUV)
{
	At(1) = (SKRTFLOAT)IQUV[0];
	At(2) = (SKRTFLOAT)IQUV[1];
	At(3) = (SKRTFLOAT)IQUV[2];
	At(4) = (SKRTFLOAT)IQUV[3];
}

/* Left multiply 
 * [ 1           0            0  0
 *   0  cos(2*eta)  -sin(2*eta)  0
 *   0  sin(2*eta)   cos(2*eta)  0
 *   0           0            0  1 ]
 * (Mishchenko 2002)
 * onto this vector.
 */
void skRTStokesVector::RotatePolarPlaneThru( double cosEta, double sinEta )
{
	double cte, ste, a, b;
	
	cte = cosEta*cosEta - sinEta*sinEta;
	ste = 2.0*cosEta*sinEta;

	a = cte*m_phasematrixstorage[1] - ste*m_phasematrixstorage[2];
	b = ste*m_phasematrixstorage[1] + cte*m_phasematrixstorage[2];

	m_phasematrixstorage[1] = a;
	m_phasematrixstorage[2] = b;
}

#if defined(NXDEBUG)

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::At		2004-1-30
 *-------------------------------------------------------------------------*/

SKRTFLOAT& skRTStokesVector::At(int idx)
{
	bool	ok;

	ok =    (idx > 0) && (idx < 5);
	if (!ok)
	{
		NXTRACE(("skRTStokesVector::At, index (%d) is out of range (1-4)\n", (int)idx));
		NXASSERT( false );
		idx = 1;
	}
	return m_phasematrixstorage[ --idx ];
}

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::At		2004-1-30
 *-------------------------------------------------------------------------*/

SKRTFLOAT skRTStokesVector::At(int idx) const
{
	bool	ok;

	ok =    (idx > 0) && (idx < 5);
	if (!ok)
	{
		NXTRACE(("skRTStokesVector::At, index (%d) is out of range (1-4)\n", (int)idx));
		NXASSERT( false );
		idx = 1;
	}
	return m_phasematrixstorage[ --idx ];
}

#endif

/*---------------------------------------------------------------------------
 *'					skRTStokesVector::dgemv		2004-2-6
 *	Wrapper for the Intel Math Kernel Library dgemv (BLAS 2) which computes
 *	the product of a matrix and a vector:  y = alpha*a*x + beta*y
 *
 *	A "slow" version of the code is provided in the comments for those
 *	who want to 
 *-------------------------------------------------------------------------*/

//skRTStokesVector& skRTStokesVector::dgemv( double alpha, skRTPhaseMatrix& a, skRTStokesVector& x, double beta)
//{
//	cblas_dgemv( CblasColMajor, CblasNoTrans, 4, 4, alpha, a.ArrayBasePtr(), 4, x.ArrayBasePtr(), 1, beta, ArrayBasePtr(), 1);
/*	
	for (int k = 1; k <= 4; k++)			// Do the (4x4)#(4x1) matrix multiplication explicitly.
	{
		At(k)  = (   a.At(k,1)*x.At(1)
				   + a.At(k,2)*x.At(2)
		           + a.At(k,3)*x.At(3)
				   + a.At(k,4)*x.At(4))*alpha + beta*At(k);
	}
*/
//	return *this;
//}

template<>
void skRTStokesVector::SetToZero(skRTStokesVector& a){
	a.SetTo(0.0);
}

template<>
void skRTStokesVector::SetToZero(double& a){
	a = 0.0;
}

template<>
bool skRTStokesVector::IsNegative(const skRTStokesVector& a)
{
	return (a.At(1)<0.0);
}

template<>
bool skRTStokesVector::IsNegative(const double& a)
{
	return (a<0.0);
}

template<>
void skRTStokesVector::SetNansToZero(skRTStokesVector& a)
{
	if( a.At(1) != a.At(1) )
		SetToZero(a);
}

template<>
void skRTStokesVector::SetNansToZero(double& a)
{
	// if compiler is set to ffast math this can be optimized out
	// and cause problems
	if( a != a )
		SetToZero(a);
}

template<>
void skRTStokesVector::SetNegativesToZero(skRTStokesVector& a)
{
	if(IsNegative(a)) SetToZero(a);
}


template<>
void skRTStokesVector::SetNegativesToZero(double& a)
{
	if(IsNegative(a)) SetToZero(a);
}
