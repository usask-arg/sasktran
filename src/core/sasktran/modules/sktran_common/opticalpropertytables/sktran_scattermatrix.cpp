#include "../sktran_common.h"


void SKTRAN_Stokes_NC::LMultBy(const SKTRAN_ScatMat_MIMSNC& s) 
{
	s.LApplyTo(this); 
}

void SKTRAN_Stokes_NC::Normalize( )
{
    double normfac( m_I>0.0?(1.0/m_I):0.0 );
    m_I *= normfac;
    m_Q *= normfac;
    m_U *= normfac;
}


//void SKTRAN_Stokes_RONC::RotatePolarPlaneThru( double cosEta, double sinEta )
//{
//	SKRTFLOAT cte, ste, a, b;
//
//	cte = cosEta*cosEta - sinEta*sinEta;
//	ste = 2.0*cosEta*sinEta;
//
//	a = cte*m_Q - ste*m_U;
//	b = ste*m_Q + cte*m_U;
//
//	m_Q = a;
//	m_U = b;
//}

SKTRAN_Stokes_NC operator*(const skRTPhaseMatrix&p, const SKTRAN_Stokes_NC& v){
	SKRTFLOAT I, Q, U, V;
	
	I = p.At(1,1)*v.I() + p.At(1,2)*v.Q() + p.At(1,3)*v.U() + p.At(1,4)*v.V();
	Q = p.At(2,1)*v.I() + p.At(2,2)*v.Q() + p.At(2,3)*v.U() + p.At(2,4)*v.V();
	U = p.At(3,1)*v.I() + p.At(3,2)*v.Q() + p.At(3,3)*v.U() + p.At(3,4)*v.V();
	V = p.At(4,1)*v.I() + p.At(4,2)*v.Q() + p.At(4,3)*v.U() + p.At(4,4)*v.V();

    SKTRAN_Stokes_NC result;
    result.SetTo(I, Q, U, V);

	return result;
}


template<>
void skRTStokesVector::SetToZero(SKTRAN_Stokes_NC& a){
	a.SetTo(0.0);
}

template<>
bool skRTStokesVector::IsNegative(const SKTRAN_Stokes_NC& a)
{
	return (a.I()<0.0);
}

template<>
void skRTStokesVector::SetNansToZero(SKTRAN_Stokes_NC& a)
{
	if( a.I() != a.I() )
		skRTStokesVector::SetToZero(a);
}

template<>
void skRTStokesVector::SetNegativesToZero(SKTRAN_Stokes_NC& a)
{
	if(IsNegative(a)) skRTStokesVector::SetToZero(a);
}


/*---------------------------------------------------------------------------
 *'      SKTRAN_ScatMat_RONC::SKTRAN_ScatMat_RONC		2015-06-04 */
/**	Construct the phase matrix and initialize all of the elements to zero.
 */
/*-------------------------------------------------------------------------*/

SKTRAN_ScatMat_MIMSNC::SKTRAN_ScatMat_MIMSNC()
{
	SetTo(0.0);
}


/*---------------------------------------------------------------------------
 *'     SKTRAN_ScatMat_RONC::SKTRAN_ScatMat_RONC		2015-06-04 */
/** Copy the phase matrix.  Simply copy the entire contents of the other
  *	matrix
 */
/*-------------------------------------------------------------------------*/

SKTRAN_ScatMat_MIMSNC::SKTRAN_ScatMat_MIMSNC( const SKTRAN_ScatMat_MIMSNC& other )
{
	*this = other;
}

SKTRAN_ScatMat_MIMSNC::SKTRAN_ScatMat_MIMSNC ( const skRTPhaseMatrix& p )
{
	m_p11 = p.At(1,1);
	m_p21 = p.At(2,1);
	m_p22 = p.At(2,2);
	m_p33 = p.At(3,3);
}

/*---------------------------------------------------------------------------
 *         SKTRAN_ScatMat_RONC::At					2015-06-04      */
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

void SKTRAN_ScatMat_MIMSNC::AssignAt ( int row, int col, SKRTFLOAT val )
{
#if defined(NXDEBUG)
	bool	ok;

	ok =    (col > 0) && (col < 5)
		 && (row > 0) && (row < 5);
	if (!ok)
	{
		NXTRACE(("SKTRAN_ScatMat_RONC::At, indces are out of range, (%d, %d)\n", (int)row, (int)col));
		NXASSERT( false );
		col = 1;
		row = 1;
	}
#endif

	int rxc = row^col;
	if( 0 == rxc ){
		// Reference is on the diagonal
		if( 1 == row ){
			m_p11 = val;
		} else if(2 == row){
			m_p22 = val;
		} else if(3 == row){
			m_p33 = val;
		} else{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_ScatMat_RONC::AssignAt, Assignment (%1.16e) ignored, element (%u,%u) is assumed to be zero.", val, row, col );
		}
	} else if ( 3 == rxc ){
		// Reference is (1,2) or (2,1)
		m_p21 = val;
	} else{
		// Reference is off-diagonal and not (1,2) or (2,1) -- warning, these are assumed to be zero 
		nxLog::Record(NXLOG_WARNING, "SKTRAN_ScatMat_RONC::AssignAt, Assignment (%1.16e) ignored, element (%u,%u) is assumed to be zero.", val, row, col );
	}
}

/*---------------------------------------------------------------------------
 *'					SKTRAN_ScatMat_RONC::At					2015-06-04
 *-------------------------------------------------------------------------*/

SKRTFLOAT SKTRAN_ScatMat_MIMSNC::At ( int row, int col) const
{

#if defined(NXDEBUG)
	bool	ok;

	ok =    (col > 0) && (col < 5)
		 && (row > 0) && (row < 5);
	if (!ok)
	{
		NXTRACE(("SKTRAN_ScatMat_RONC::At, indces are out of range, (%d, %d)\n", (int)row, (int)col));
		NXASSERT( false );
		col = 1;
		row = 1;
	}
#endif
	
	int rxc = row^col;
	if( 0 == rxc ){
		// Query is on the diagonal
		if( 1 == row ){
			return m_p11; // (1,1)
		} else if( 2 == row ){
			return m_p22; // (2,2)
		} else if( 3 == row ){
			return m_p33; // (3,3)
		} else{
			return 0.0; // (4,4)
		}
	} else if ( 3 == rxc ){
		// Querry is (1,2) or (2,1)
		return m_p21;
	} else{
		// Query is off-diagonal and not (1,2) or (2,1)
		return 0.0;
	}
}

/*---------------------------------------------------------------------------
 *'					SKTRAN_ScatMat_RONC::operator = 					2015-06-04
 *-------------------------------------------------------------------------*/

SKTRAN_ScatMat_MIMSNC& SKTRAN_ScatMat_MIMSNC::operator = ( const SKTRAN_ScatMat_MIMSNC& other )
{
	m_p11 = other.m_p11;
	m_p21 = other.m_p21;
	m_p22 = other.m_p22;
	m_p33 = other.m_p33;
	return *this;
}


/*---------------------------------------------------------------------------
 *'					SKTRAN_ScatMat_RONC::SetTo		2015-06-04
 *-------------------------------------------------------------------------*/

void SKTRAN_ScatMat_MIMSNC::SetTo( double value )
{
	SKRTFLOAT val = (SKRTFLOAT)value;
	m_p11 = val;
	m_p21 = val;
	m_p22 = val;
	m_p33 = val;
}

/*---------------------------------------------------------------------------
 *'					SKTRAN_ScatMat_RONC::SetTo		2015-06-04
 *-------------------------------------------------------------------------*/

void SKTRAN_ScatMat_MIMSNC::SetTo( float value )
{
	SKRTFLOAT val = (SKRTFLOAT)value;
	m_p11 = val;
	m_p21 = val;
	m_p22 = val;
	m_p33 = val;
}

void SKTRAN_ScatMat_MIMSNC::SetToIdentity ( )
{
	m_p11 = 1.0;
	m_p21 = 0.0;
	m_p22 = 1.0;
	m_p33 = 1.0;
}

void SKTRAN_ScatMat_MIMSNC::AddToThis ( const SKTRAN_ScatMat_MIMSNC& other, double weight )
{
	m_p11 += weight*other.m_p11;
	m_p21 += weight*other.m_p21;
	m_p22 += weight*other.m_p22;
	m_p33 += weight*other.m_p33;
}

/*---------------------------------------------------------------------------
 *'					SKTRAN_ScatMat_RONC::operator +=		2015-06-04
 *-------------------------------------------------------------------------*/

SKTRAN_ScatMat_MIMSNC& SKTRAN_ScatMat_MIMSNC::operator +=( const SKTRAN_ScatMat_MIMSNC& other)
{
	m_p11 += other.m_p11;
	m_p21 += other.m_p21;
	m_p22 += other.m_p22;
	m_p33 += other.m_p33;
	return *this;
}

/*---------------------------------------------------------------------------
 *'					SKTRAN_ScatMat_RONC::operator +		2015-06-04
 *-------------------------------------------------------------------------*/

SKTRAN_ScatMat_MIMSNC SKTRAN_ScatMat_MIMSNC::operator + ( const SKTRAN_ScatMat_MIMSNC& other) const
{
	SKTRAN_ScatMat_MIMSNC	answer;
	answer.m_p11 = m_p11 + other.m_p11;
	answer.m_p21 = m_p21 + other.m_p21;
	answer.m_p22 = m_p22 + other.m_p22;
	answer.m_p33 = m_p33 + other.m_p33;
	return answer;
}
/*---------------------------------------------------------------------------
 *'					SKTRAN_ScatMat_RONC::operator -		2015-06-04
 *-------------------------------------------------------------------------*/

SKTRAN_ScatMat_MIMSNC SKTRAN_ScatMat_MIMSNC::operator - ( const SKTRAN_ScatMat_MIMSNC& other) const
{
	SKTRAN_ScatMat_MIMSNC	answer;
	answer.m_p11 = m_p11 - other.m_p11;
	answer.m_p21 = m_p21 - other.m_p21;
	answer.m_p22 = m_p22 - other.m_p22;
	answer.m_p33 = m_p33 - other.m_p33;
	return answer;
}
/*---------------------------------------------------------------------------
 *'					SKTRAN_ScatMat_RONC::operator *		2015-06-04
 *-------------------------------------------------------------------------*/

SKTRAN_ScatMat_MIMSNC SKTRAN_ScatMat_MIMSNC::operator * (  double value) const
{
    SKTRAN_ScatMat_MIMSNC	answer;
	SKRTFLOAT val=(SKRTFLOAT)value;
	answer.m_p11 = m_p11 * val;
	answer.m_p21 = m_p21 * val;
	answer.m_p22 = m_p22 * val;
	answer.m_p33 = m_p33 * val;
	return answer;
}

/*---------------------------------------------------------------------------
 *'					SKTRAN_ScatMat_RONC::operator *		2015-06-04
 *-------------------------------------------------------------------------*/

SKTRAN_ScatMat_MIMSNC SKTRAN_ScatMat_MIMSNC::operator * (  float value) const
{
    SKTRAN_ScatMat_MIMSNC	answer;
	SKRTFLOAT val=(SKRTFLOAT)value;
	answer.m_p11 = m_p11 * val;
	answer.m_p21 = m_p21 * val;
	answer.m_p22 = m_p22 * val;
	answer.m_p33 = m_p33 * val;
	return answer;
}

/*---------------------------------------------------------------------------
 *'					SKTRAN_ScatMat_RONC::operator *=		2015-06-04
 *-------------------------------------------------------------------------*/

SKTRAN_ScatMat_MIMSNC& SKTRAN_ScatMat_MIMSNC::operator *= ( double value)
{
	SKRTFLOAT val=(SKRTFLOAT)value;
	m_p11 *= val;
	m_p21 *= val;
	m_p22 *= val;
	m_p33 *= val;
	return *this;
}

/*---------------------------------------------------------------------------
 *'					SKTRAN_ScatMat_RONC::operator *=		2015-06-04
 *-------------------------------------------------------------------------*/

SKTRAN_ScatMat_MIMSNC& SKTRAN_ScatMat_MIMSNC::operator *= ( float value)
{
	SKRTFLOAT val=(SKRTFLOAT)value;
	m_p11 *= val;
	m_p21 *= val;
	m_p22 *= val;
	m_p33 *= val;
	return *this;
}

/*---------------------------------------------------------------------------
 *'					 SKTRAN_ScatMat_RONC::operator *		2015-06-04
 *-------------------------------------------------------------------------*/

SKTRAN_Stokes_NC SKTRAN_ScatMat_MIMSNC::operator *( const SKTRAN_Stokes_NC& x ) const
{
	SKTRAN_Stokes_NC	result;
	result.Assign_I ( m_p11*x.I() + m_p21*x.Q() );
	result.Assign_Q ( m_p21*x.I() + m_p22*x.Q() );
	result.Assign_U ( m_p33*x.U() );
	result.Assign_V ( 0.0 );

	return result;
}

//
//SKTRAN_ScatMat_RONC& SKTRAN_ScatMat_RONC::RMultBy ( const SKTRAN_ScatMat_RONC& other )
//{
//	nxLog::Record(NXLOG_ERROR, "SKTRAN_ScatMat_RONC::RMultBy, Would need five elements to store product! Should remove this function.");
//	return *this;
//}
//
//SKTRAN_ScatMat_RONC& SKTRAN_ScatMat_RONC::LMultBy ( const SKTRAN_ScatMat_RONC& other )
//{
//	nxLog::Record(NXLOG_ERROR, "SKTRAN_ScatMat_RONC::LMultBy, Would need five elements to store product! Should remove this function.");
//	return *this;
//}

void SKTRAN_ScatMat_MIMSNC::LApplyTo( SKTRAN_Stokes_NC* v ) const
{
	//double r_I = m_p11*v.I() + m_p21*v.Q();
	//double r_Q = m_p21*v.I() + m_p22*v.Q();
	//double r_U = m_p33*v.U();
	v->SetTo( m_p11*v->I() + m_p21*v->Q()
		   ,  m_p21*v->I() + m_p22*v->Q()
		   ,  m_p33*v->U()
		   ,  0.0);
}


void SKTRAN_ScatMat_MIMSNC::Normalize ( )
{
    double normfac( m_p11>0.0 ? (1.0/m_p11) : 0.0 );
    m_p11 *= normfac;
    m_p21 *= normfac;
    m_p22 *= normfac;
    m_p33 *= normfac;

}


SKTRAN_ScatMat_Rot::SKTRAN_ScatMat_Rot ( ){
	m_cos2eta = 1.0;
	m_sin2eta = 0.0;
}

SKTRAN_ScatMat_Rot::SKTRAN_ScatMat_Rot ( double cosEta, double sinEta ){
	SetAngle( cosEta, sinEta );
}

void SKTRAN_ScatMat_Rot::SetAngle ( double cosEta, double sinEta ){
	m_cos2eta = cosEta*cosEta - sinEta*sinEta;
	m_sin2eta = 2.0*cosEta*sinEta;
}

SKTRAN_Stokes_NC SKTRAN_ScatMat_Rot::operator* ( const SKTRAN_Stokes_NC& v ) const{
	
	SKTRAN_Stokes_NC result(v);
	
	result.Assign_Q( m_cos2eta*v.Q() - m_sin2eta*v.U() );
	result.Assign_U( m_sin2eta*v.Q() + m_cos2eta*v.U() );

	return result;
}

double SKTRAN_ScatMat_Rot::At( int row, int col ) const{

	double result;

	if( row == col ){
		if( 2==row || 3==row ){
			result = m_cos2eta;
		} else{
			result = 1.0;
		}
	} else if( 1 == (row^col) ){
		if( 2 == row ){
			result =  m_sin2eta;
		} else{
			result = -m_sin2eta;
		}
	} else{
		result = 0.0;
	}

	return result;
}



#define SKTRAN_PhaseMat_MIMSNC_sub2ind(row, col) (3*(row-1)+col-1)

SKTRAN_PhaseMat_MIMSNC::SKTRAN_PhaseMat_MIMSNC()
{
	SetTo(0.0);
}

SKTRAN_PhaseMat_MIMSNC::SKTRAN_PhaseMat_MIMSNC( const SKTRAN_ScatMat_MIMSNC& other )
{
	*this = other;
}

SKTRAN_PhaseMat_MIMSNC::SKTRAN_PhaseMat_MIMSNC( const SKTRAN_ScatMat_Rot& other )
{
	*this = other;
}


SKTRAN_PhaseMat_MIMSNC::SKTRAN_PhaseMat_MIMSNC( const SKTRAN_PhaseMat_MIMSNC& other )
{
	*this = other;
}

SKTRAN_PhaseMat_MIMSNC::SKTRAN_PhaseMat_MIMSNC ( const skRTPhaseMatrix& p )
{
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] = p.At(1,1); // Remember SKTRAN_PhaseMat_MIMSNC is row major whereas skRTPhaseMatrix is column-major
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] = p.At(1,2);
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)] = p.At(1,3);
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] = p.At(2,1);
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] = p.At(2,2);
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] = p.At(2,3);
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)] = p.At(3,1);
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] = p.At(3,2);
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)] = p.At(3,3);
}

void SKTRAN_PhaseMat_MIMSNC::AssignAt ( int row, int col, SKRTFLOAT val )
{
#if defined(NXDEBUG)
	bool	ok;

	ok =    (col > 0) && (col < 5)
		 && (row > 0) && (row < 5);
	if (!ok)
	{
		NXTRACE(("SKTRAN_ScatMat_RONC::At, indces are out of range, (%d, %d)\n", (int)row, (int)col));
		NXASSERT( false );
		col = 1;
		row = 1;
	}
#endif

	 if( 4!=row && 4!=col) m_p[ SKTRAN_PhaseMat_MIMSNC_sub2ind(row,col) ] = val;
}


SKRTFLOAT SKTRAN_PhaseMat_MIMSNC::At ( int row, int col) const
{

#if defined(NXDEBUG)
	bool	ok;

	ok =    (col > 0) && (col < 5)
		 && (row > 0) && (row < 5);
	if (!ok)
	{
		NXTRACE(("SKTRAN_ScatMat_RONC::At, indces are out of range, (%d, %d)\n", (int)row, (int)col));
		NXASSERT( false );
		col = 1;
		row = 1;
	}
#endif
	
	
	 if( 4!=row && 4!=col){
		 return m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(row,col) ];
	 } else{
		 return 0.0;
	 }
}

void SKTRAN_PhaseMat_MIMSNC::SetTo( double value )
{
	m_p.fill( (SKRTFLOAT)value);
}

void SKTRAN_PhaseMat_MIMSNC::SetTo( float value )
{
	m_p.fill( (SKRTFLOAT)value);
}

void SKTRAN_PhaseMat_MIMSNC::SetToIdentity ( )
{
	SetTo( 0.0 );
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] = 1.0;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] = 1.0;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)] = 1.0;
}

SKTRAN_PhaseMat_MIMSNC&  SKTRAN_PhaseMat_MIMSNC::RMultBy ( const SKTRAN_ScatMat_MIMSNC& other ){
	SKRTFLOAT temp1, temp2;

	temp1 = other.m_p11*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] + other.m_p21*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)];
	temp2 = other.m_p21*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] + other.m_p22*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] = temp2;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)] *= other.m_p33;

	temp1 = other.m_p11*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] + other.m_p21*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)];
	temp2 = other.m_p21*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] + other.m_p22*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] = temp2;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] *= other.m_p33;

	temp1 = other.m_p11*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)] + other.m_p21*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)];
	temp2 = other.m_p21*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)] + other.m_p22*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] = temp2;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)] *= other.m_p33;

	return *this;
}



SKTRAN_PhaseMat_MIMSNC&  SKTRAN_PhaseMat_MIMSNC::LMultBy ( const SKTRAN_ScatMat_MIMSNC& other ){
	
	SKRTFLOAT p11, p12, p13, p21, p22, p23;
	p11 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)];
	p12 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)];
	p13 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)];
	p21 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)];
	p22 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)];
	p23 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)];
	
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] = other.m_p11*p11 + other.m_p21*p21;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] = other.m_p11*p12 + other.m_p21*p22;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)] = other.m_p11*p13 + other.m_p21*p23;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] = other.m_p21*p11 + other.m_p22*p21;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] = other.m_p21*p12 + other.m_p22*p22;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] = other.m_p21*p13 + other.m_p22*p23;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)] *= other.m_p33;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] *= other.m_p33;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)] *= other.m_p33;

	return *this;
}


SKTRAN_PhaseMat_MIMSNC&  SKTRAN_PhaseMat_MIMSNC::RMultBy ( const SKTRAN_ScatMat_Rot& other ){
	SKRTFLOAT temp1, temp2;
	
	temp1 =  other.m_cos2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] + other.m_sin2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)];
	temp2 = -other.m_sin2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] + other.m_cos2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)] = temp2;
	
	temp1 =  other.m_cos2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] + other.m_sin2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)];
	temp2 = -other.m_sin2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] + other.m_cos2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] = temp2;
	
	temp1 =  other.m_cos2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] + other.m_sin2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)];
	temp2 = -other.m_sin2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] + other.m_cos2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)] = temp2;
	
	return *this;
}



SKTRAN_PhaseMat_MIMSNC&  SKTRAN_PhaseMat_MIMSNC::LMultBy ( const SKTRAN_ScatMat_Rot& other ){
	SKRTFLOAT temp1, temp2;
	
	temp1 = other.m_cos2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] - other.m_sin2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)];
	temp2 = other.m_sin2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] + other.m_cos2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)] = temp2;
	
	temp1 = other.m_cos2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] - other.m_sin2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)];
	temp2 = other.m_sin2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] + other.m_cos2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] = temp2;
	
	temp1 = other.m_cos2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] - other.m_sin2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)];
	temp2 = other.m_sin2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] + other.m_cos2eta*m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)] = temp2;
	
	return *this;
}




SKTRAN_PhaseMat_MIMSNC&  SKTRAN_PhaseMat_MIMSNC::RMultBy ( const SKTRAN_PhaseMat_MIMSNC& other ){
	SKRTFLOAT temp1, temp2, temp3;
	
	temp1 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)];
	temp2 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)];
	temp3 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] = temp2;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)] = temp3;
	
	temp1 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)];
	temp2 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)];
	temp3 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] = temp2;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] = temp3;

	temp1 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)];
	temp2 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)];
	temp3 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] = temp2;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)] = temp3;

	return *this;
}

SKTRAN_PhaseMat_MIMSNC&  SKTRAN_PhaseMat_MIMSNC::LMultBy ( const SKTRAN_PhaseMat_MIMSNC& other ){
	SKRTFLOAT temp1, temp2, temp3;
	
	temp1 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)];
	temp2 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)];
	temp3 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] = temp2;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)] = temp3;
	
	temp1 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)];
	temp2 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)];
	temp3 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] = temp2;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] = temp3;
	
	temp1 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)];
	temp2 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)];
	temp3 = m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)]*other.m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)];
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)] = temp1;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] = temp2;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)] = temp3;

	return *this;
}





SKTRAN_PhaseMat_MIMSNC& SKTRAN_PhaseMat_MIMSNC::operator +=( const SKTRAN_PhaseMat_MIMSNC& other)
{
	auto oit = other.m_p.cbegin();
	std::for_each( m_p.begin(), m_p.end(), [&oit]( SKRTFLOAT& pelem ) { pelem += *oit; ++oit;} );
	return *this;
}


SKTRAN_PhaseMat_MIMSNC& SKTRAN_PhaseMat_MIMSNC::operator = ( const SKTRAN_PhaseMat_MIMSNC& other )
{
	m_p = other.m_p;
	return *this;
}

SKTRAN_PhaseMat_MIMSNC& SKTRAN_PhaseMat_MIMSNC::operator = ( const SKTRAN_ScatMat_Rot& other )
{
	SetTo( 0.0 );
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] =  1.0;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] =  other.m_cos2eta;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)] = -other.m_sin2eta;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)] =  other.m_sin2eta;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)] =  other.m_cos2eta;
	
	return *this;
}

SKTRAN_PhaseMat_MIMSNC& SKTRAN_PhaseMat_MIMSNC::operator = ( const SKTRAN_ScatMat_MIMSNC& other )
{
	SetTo( 0.0 );
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)] = other.m_p11;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)] = other.m_p21;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)] = other.m_p21;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)] = other.m_p22;
	m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)] = other.m_p33;
	return *this;
}

SKTRAN_PhaseMat_MIMSNC SKTRAN_PhaseMat_MIMSNC::operator + ( const SKTRAN_PhaseMat_MIMSNC& other) const
{
	SKTRAN_PhaseMat_MIMSNC	answer;
	auto pit =       m_p.cbegin();
	auto oit = other.m_p.cbegin();
	std::for_each( answer.m_p.begin(), answer.m_p.end(), [&pit, &oit]( SKRTFLOAT& pelem ) { pelem = *pit + *oit; ++pit, ++oit;} );
	return answer;
}

SKTRAN_PhaseMat_MIMSNC SKTRAN_PhaseMat_MIMSNC::operator - ( const SKTRAN_PhaseMat_MIMSNC& other) const
{
	SKTRAN_PhaseMat_MIMSNC	answer;
	auto pit =       m_p.cbegin();
	auto oit = other.m_p.cbegin();
	std::for_each( answer.m_p.begin(), answer.m_p.end(), [&pit, &oit]( SKRTFLOAT& pelem ) { pelem = *pit - *oit; ++pit, ++oit;} );
	return answer;
}


SKTRAN_PhaseMat_MIMSNC SKTRAN_PhaseMat_MIMSNC::operator * ( double value ) const
{
	SKTRAN_PhaseMat_MIMSNC	answer;
	SKRTFLOAT val = (SKRTFLOAT) value;
	auto pit =       m_p.cbegin();
	std::for_each( answer.m_p.begin(), answer.m_p.end(), [&pit, val]( SKRTFLOAT& pelem ) { pelem = *pit * val; ++pit;} );
	return answer;
}


SKTRAN_PhaseMat_MIMSNC SKTRAN_PhaseMat_MIMSNC::operator * ( float value ) const
{
	SKTRAN_PhaseMat_MIMSNC	answer;
	SKRTFLOAT val = (SKRTFLOAT) value;
	auto pit =       m_p.cbegin();
	std::for_each( answer.m_p.begin(), answer.m_p.end(), [&pit, val]( SKRTFLOAT& pelem ) { pelem = *pit * val;} );
	return answer;
}


SKTRAN_PhaseMat_MIMSNC& SKTRAN_PhaseMat_MIMSNC::operator *= ( double value)
{
	SKRTFLOAT val=(SKRTFLOAT)value;
	std::for_each( m_p.begin(), m_p.end(), [val]( SKRTFLOAT& pelem ) { pelem *= val;} );
	return *this;
}

SKTRAN_PhaseMat_MIMSNC& SKTRAN_PhaseMat_MIMSNC::operator *= ( float value)
{
	SKRTFLOAT val=(SKRTFLOAT)value;
	std::for_each( m_p.begin(), m_p.end(), [val]( SKRTFLOAT& pelem ) { pelem *= val;} );
	return *this;
}

/*---------------------------------------------------------------------------
 *'					 SKTRAN_ScatMat_RONC::operator *		2015-06-04
 *-------------------------------------------------------------------------*/

SKTRAN_Stokes_NC SKTRAN_PhaseMat_MIMSNC::operator *( const SKTRAN_Stokes_NC& x ) const
{
	SKTRAN_Stokes_NC	answer;
	answer.Assign_I ( m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,1)]*x.I() + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,2)]*x.Q() + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1,3)]*x.U() );
	answer.Assign_Q ( m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,1)]*x.I() + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,2)]*x.Q() + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2,3)]*x.U() );
	answer.Assign_U ( m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,1)]*x.I() + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,2)]*x.Q() + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3,3)]*x.U() );
	answer.Assign_V ( 0.0 );

	return answer;
}


void SKTRAN_PhaseMat_MIMSNC::LApplyTo(SKTRAN_Stokes_NC* v) const
{
	double I(v->I()), Q(v->Q()), U(v->U());
	double Iprime( m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1, 1)]*I + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1, 2)]*Q + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(1, 3)]*U );
	double Qprime( m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2, 1)]*I + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2, 2)]*Q + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(2, 3)]*U );
	double Uprime( m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3, 1)]*I + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3, 2)]*Q + m_p[SKTRAN_PhaseMat_MIMSNC_sub2ind(3, 3)]*U );
	
	v->SetTo( Iprime, Qprime, Uprime, 0.0 );
}
