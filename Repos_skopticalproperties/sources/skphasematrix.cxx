#include <skopticalproperties21.h>
#include <limits.h>
#include <float.h>
/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::skRTPhaseMatrix		2003-10-22 */
/**	Construct the phase matrix and initialize all of the elements to zero.
 */
/*-------------------------------------------------------------------------*/

skRTPhaseMatrix::skRTPhaseMatrix()
{
	NXASSERT( sizeof(skRTPhaseMatrix) == 16*sizeof(m_phasematrixstorage[0]) );			// **** Some code relies upon fact its exactly 4 doubles 
	SetTo(0.0);
}


/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::skRTPhaseMatrix		2003-10-22 */
/** Copy the phase matrix.  Simply copy the entire contents of the other
  *	matrix
 */
/*-------------------------------------------------------------------------*/

skRTPhaseMatrix::skRTPhaseMatrix( const skRTPhaseMatrix& other )
{
	*this = other;
}
/*---------------------------------------------------------------------------
 *					skRTPhaseMatrix::At					2003-10-22      */
/** Return a reference to the element at the specified location. The code uses a 1 based
 *	indexing.  Since it returns a reference it can be modified by the user.
 *	This code has been specifically designed so the indexes passed
 *	in look like the indices used in the matrix formulation.  The DEBUG version
 *	of the code will check the index bounds but the RELEASE version will not.
 *
 *	@param row
 *		The vertical (1st) index of the 4x4 matrix
 *
 *	@param col
 *		The horizontal (2nd) index of the 4x4 matrix.
 *
 */
 /*-------------------------------------------------------------------------*/

SKRTFLOAT& skRTPhaseMatrix::At ( int row, int col)
{
#if defined(NXDEBUG)
	bool	ok;

	ok =    (col > 0) && (col < 5)
		 && (row > 0) && (row < 5);
	if (!ok)
	{
		NXTRACE(("skRTPhaseMatrix::At, indces are out of range, (%d, %d)\n", (int)row, (int)col));
		NXASSERT( false );
		col = 1;
		row = 1;
	}
#endif

	return m_phasematrixstorage[ ((--col)<< 2) | (--row) ];
}
/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::At					2003-10-22
 *-------------------------------------------------------------------------*/

SKRTFLOAT skRTPhaseMatrix::At ( int row, int col) const
{

#if defined(NXDEBUG)
	bool	ok;

	ok =    (col > 0) && (col < 5)
		 && (row > 0) && (row < 5);
	if (!ok)
	{
		NXTRACE(("skRTPhaseMatrix::At, indces are out of range, (%d, %d)\n", (int)row, (int)col));
		NXASSERT( false );
		col = 1;
		row = 1;
	}
#endif

	return m_phasematrixstorage[ ((--col)<< 2) | (--row) ];
}

/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::operator = 					2003-10-22
 *-------------------------------------------------------------------------*/

skRTPhaseMatrix& skRTPhaseMatrix::operator = ( const skRTPhaseMatrix& other )
{
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage); i++) m_phasematrixstorage[i] = other.m_phasematrixstorage[i];
	return *this;
}


/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::SetTo		2004-1-30
 *-------------------------------------------------------------------------*/

void skRTPhaseMatrix::SetTo( double val )
{
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage); i++) m_phasematrixstorage[i] = (SKRTFLOAT)val;
}

/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::SetTo		2004-1-30
 *-------------------------------------------------------------------------*/

void skRTPhaseMatrix::SetTo( float val )
{
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage); i++) m_phasematrixstorage[i] = (SKRTFLOAT)val;
}

/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::operator +=		2004-1-30
 *-------------------------------------------------------------------------*/

skRTPhaseMatrix& skRTPhaseMatrix::operator +=( const skRTPhaseMatrix& other)
{
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage); i++) m_phasematrixstorage[i] += other.m_phasematrixstorage[i];
	return *this;
}

/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::operator +		2004-1-30
 *-------------------------------------------------------------------------*/

skRTPhaseMatrix skRTPhaseMatrix::operator + ( const skRTPhaseMatrix& other) const
{
	skRTPhaseMatrix	answer;
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage); i++)
	{
		answer.m_phasematrixstorage[i] = m_phasematrixstorage[i] + other.m_phasematrixstorage[i];
	}
	return answer;
}
/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::operator -		2004-1-30
 *-------------------------------------------------------------------------*/

skRTPhaseMatrix skRTPhaseMatrix::operator - ( const skRTPhaseMatrix& other) const
{
	skRTPhaseMatrix	answer;
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage); i++)
	{
		answer.m_phasematrixstorage[i] = m_phasematrixstorage[i] - other.m_phasematrixstorage[i];
	}
	return answer;
}
/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::operator *		2004-1-30
 *-------------------------------------------------------------------------*/

skRTPhaseMatrix skRTPhaseMatrix::operator * (  double value) const
{
	skRTPhaseMatrix	answer;
	SKRTFLOAT			val=(SKRTFLOAT)value;

	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage); i++)
	{
		answer.m_phasematrixstorage[i] = m_phasematrixstorage[i]*val;
	}
	return answer;
}

/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::operator *		2004-1-30
 *-------------------------------------------------------------------------*/

skRTPhaseMatrix skRTPhaseMatrix::operator * (  float value) const
{
	skRTPhaseMatrix	answer;
	SKRTFLOAT			val=(SKRTFLOAT)value;

	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage); i++)
	{
		answer.m_phasematrixstorage[i] = m_phasematrixstorage[i]*val;
	}
	return answer;
}

/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::operator *=		2004-1-30
 *-------------------------------------------------------------------------*/

skRTPhaseMatrix& skRTPhaseMatrix::operator *= ( double value)
{
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage); i++) m_phasematrixstorage[i] *= (SKRTFLOAT)value;
	return *this;
}
/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::operator *=		2004-1-30
 *-------------------------------------------------------------------------*/

skRTPhaseMatrix& skRTPhaseMatrix::operator *= ( float value)
{
	for (size_t i=0; i < N_ELEMENTS(m_phasematrixstorage); i++) m_phasematrixstorage[i] *= (SKRTFLOAT)value;
	return *this;
}

/*---------------------------------------------------------------------------
 *'					 skRTPhaseMatrix::operator *		2004-4-27
 *-------------------------------------------------------------------------*/

skRTStokesVector skRTPhaseMatrix::operator *( const skRTStokesVector& x ) const
{
	skRTStokesVector	answer;
	for (int k = 1; k <= 4; k++)			// Do the (4x4)#(4x1) matrix multiplication explicitly.
	{
		answer.At(k)  = (     At(k,1)*x.At(1)
							+ At(k,2)*x.At(2)
							+ At(k,3)*x.At(3)
							+ At(k,4)*x.At(4));
	}
	return answer;
}

/*---------------------------------------------------------------------------
 *'					skRTPhaseMatrix::ApplyStokesRotation		2003-10-22 */
/**	Apply a Stokes rotation to the current phase matrix.  This code implements
 *	equation A.5 of McLinden et al. Can J. Phys. 80: 375-393, 2002.  Namely, given
 *	the input is a phase matrix representing scattering within the scattering plane then 
 *	rotate the phase matrix so it has the reference polarization
 *	in the vertical direction.  This requires rotating from the vertical to the incoming
 *	direction, applying the original scattering and then rotating from the scattering plane's outgoing
 *	direction back to the vertical.
 *
 *  \par Assumptions
 *	Note this code implicitly assumes that the phase matrix is only
 *	dependendent upon 4 paramaters. This is the case for the Rayleigh and
 *	spherical MIE. More advanced spherical scattering code will need to do the
 *	Rotation as the application/multiplication of two rotation matrices
 *
 *	
 *	@param	cosmu 
 *		    The cosine of outgoing zenith angle.  This should never be exactly -1 or +1
 *			as the rotation is undefined for these angles.
 *
 *	@param cosmuprime
 *			The cosine of the incoming zenith angle.  This should never be exactly -1 or +1
 *			as the rotation is undefined for these angles.
 *
 *	@param dphi
 *			The angle in radians of (outgoing azimuth - incoming azimuth ) 
 *
 *	@param  rotatedmatrix 
 *			returns the rotated phase matrix. Must not be NULL.
*/
/*-------------------------------------------------------------------------*/

bool skRTPhaseMatrix::ApplyStokesRotation( double cosmu, double cosmuprime, double dphi, skRTPhaseMatrix* rotatedmatrix )
{
	double sinmu       = sqrt( 1 - nxmath::sqr(cosmu));
	double sinmuprime  = sqrt( 1 - nxmath::sqr(cosmuprime));
	double cosdphi     = cos(dphi);
	double sindphi     = sin(dphi);
	double p11         = At(1,1);
	double p12         = At(1,2);
	double p33         = At(3,3);
	double p34         = At(3,4);
	double y1;
	double z1;
	double y2;
	double z2;
	double i1;
	double i2;
	double c1;
	double c2;
	double s1;
	double s2;


	y1 = sinmu*sindphi;
	z1 = cosmu*sinmuprime - sinmu*cosmuprime*cosdphi;
	y2 = sinmuprime*sindphi;
	z2 = sinmu*cosmuprime - sinmuprime*cosmu*cosdphi;
	i1 = -2*atan2( y1,z1);
	i2 = -2*atan2( y2,z2);

	c1 = cos(i1);
	c2 = cos( nxmath::TWOPI + i2 );
	s1 = sin( i1 );
	s2 = sin( nxmath::TWOPI + i2 );

	rotatedmatrix->At(1,1) = (SKRTFLOAT)(p11);
	rotatedmatrix->At(1,2) = (SKRTFLOAT)(p12*c1);
	rotatedmatrix->At(1,3) = (SKRTFLOAT)(p12*s1);
	rotatedmatrix->At(1,4) = (SKRTFLOAT)0;

	rotatedmatrix->At(2,1) = (SKRTFLOAT)(p12*c2);
	rotatedmatrix->At(2,2) = (SKRTFLOAT)(c2*p11*c1 - s2*p33*s1);
	rotatedmatrix->At(2,3) = (SKRTFLOAT)(c2*p11*s1 + s2*p33*c1);
	rotatedmatrix->At(2,4) = (SKRTFLOAT)(-p34*s2);

	rotatedmatrix->At(3,1) = (SKRTFLOAT)(-p12*s2);
	rotatedmatrix->At(3,2) = (SKRTFLOAT)(-s2*p11*c1 - c2*p33*s1);
	rotatedmatrix->At(3,3) = (SKRTFLOAT)(-s2*p11*s1 + c2*p33*c1);
	rotatedmatrix->At(3,4) = (SKRTFLOAT)(-p34*c2);

	rotatedmatrix->At(4,1) = (SKRTFLOAT)(0);
	rotatedmatrix->At(4,2) = (SKRTFLOAT)(-p34*s1);
	rotatedmatrix->At(4,3) = (SKRTFLOAT)( p34*c1);
	rotatedmatrix->At(4,4) = (SKRTFLOAT)( p33);

	return true;
}


skRTPhaseMatrix& skRTPhaseMatrix::RMultBy ( const skRTPhaseMatrix& other )
{
	SKRTFLOAT      r1, r2, r3, r4;
	SKRTFLOAT      l1,l2,l3,l4;
	const_iterator mr  = other.begin()-1;
	const_iterator ml  = begin();
	iterator       res = begin();
	
	if(&other==this){
		skRTPhaseMatrix copy = *this;
		RMultBy(copy);
	} else{
		for(int n=4; n>0; --n)
		{
			l1 = *(ml+0);
			l2 = *(ml+4);
			l3 = *(ml+8);
			l4 = *(ml+12);
			r1  = l1* *(++mr);
			r1 += l2* *(++mr);
			r1 += l3* *(++mr);
			r1 += l4* *(++mr);
			r2  = l1* *(++mr);
			r2 += l2* *(++mr);
			r2 += l3* *(++mr);
			r2 += l4* *(++mr);
			r3  = l1* *(++mr);
			r3 += l2* *(++mr);
			r3 += l3* *(++mr);
			r3 += l4* *(++mr);
			r4  = l1* *(++mr);
			r4 += l2* *(++mr);
			r4 += l3* *(++mr);
			r4 += l4* *(++mr);
			*(res+0)  = r1;
			*(res+4)  = r2;
			*(res+8)  = r3;
			*(res+12) = r4;
			++res;
			++ml;
			mr-=16;
		}
	}

	return *this;
}

skRTPhaseMatrix& skRTPhaseMatrix::LMultBy ( const skRTPhaseMatrix& other )
{
	skRTPhaseMatrix temp;
	temp = other;
	temp.RMultBy(*this);
	*(this) = temp;

	return *this;
}