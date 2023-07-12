#include <sasktran.h>
#include "./include/testscattermatrices.h"


bool TestScatterMatrices::Compare ( double d1, double d2 ) const{
	return abs(d2-d1) < m_tolerance;
}

bool TestScatterMatrices::Compare ( const SKTRAN_Stokes_NC& v1, const SKTRAN_Stokes_NC& v2 ) const{
	bool ok = true;
	ok = ok && Compare( v1.I(), v2.I() );
	ok = ok && Compare( v1.Q(), v2.Q() );
	ok = ok && Compare( v1.U(), v2.U() );
	ok = ok && Compare( v1.V(), v2.V() );
	return ok;
}

template<typename SCATMAT1, typename SCATMAT2>
bool TestScatterMatrices::Compare_Impl(const SCATMAT1& p, const SCATMAT2& s) const{
	bool ok = true;
	for( int xidx=1; xidx<5; ++xidx )
		for( int yidx=1; yidx<5; ++yidx )
			ok = ok && Compare( (double) p.At(xidx,yidx), (double) s.At(xidx,yidx) );
	return ok;
}


bool TestScatterMatrices::TestRONC ( ) const
{
	bool ok = true;

	SKTRAN_Stokes_NC v;
	v.SetTo( 0.79, 0.53, 0.43, 0.0 );

	skRTPhaseMatrix p_zero;
	p_zero.SetTo(0.0);

	skRTPhaseMatrix p_eye;
	p_eye.SetTo(0.0);
	p_eye.At(1,1) = 1.0;
	p_eye.At(2,2) = 1.0;
	p_eye.At(3,3) = 1.0;

	skRTPhaseMatrix p_ronc;
	p_ronc.SetTo(0.0);
	p_ronc.At(1,1) = 1.0;
	p_ronc.At(1,2) = 0.5;
	p_ronc.At(2,1) = 0.5;
	p_ronc.At(2,2) = 0.3;
	p_ronc.At(3,3) = 0.1;


	// Test SKTRAN_ScatMat_RONC
    SKTRAN_ScatMat_MIMSNC s_ronc_zero;
	ok = ok && Compare( p_zero, s_ronc_zero );
    SKTRAN_ScatMat_MIMSNC s_ronc(p_ronc);
	ok = ok && Compare( p_ronc, s_ronc );
	ok = ok && Compare( SKTRAN_ScatMat_MIMSNC(s_ronc), s_ronc );

    SKTRAN_ScatMat_MIMSNC s_ronc_mod;
	skRTPhaseMatrix     p_ronc_mod;
	s_ronc_mod = s_ronc;
	ok = ok && Compare( s_ronc_mod, s_ronc );
	p_ronc_mod = p_ronc;

	s_ronc_mod.AssignAt(1,1, 0.19);
	p_ronc_mod.At(1,1) = 0.19;
	ok = ok && Compare(p_ronc_mod, s_ronc_mod);
	s_ronc_mod.AssignAt(2,2, 0.19);
	p_ronc_mod.At(2,2) = 0.19;
	ok = ok && Compare(p_ronc_mod, s_ronc_mod);
	s_ronc_mod.AssignAt(3,3, 0.19);
	p_ronc_mod.At(3,3) = 0.19;
	ok = ok && Compare(p_ronc_mod, s_ronc_mod);
	s_ronc_mod.AssignAt(1,2, 0.19);
	p_ronc_mod.At(1,2) = 0.19;
	p_ronc_mod.At(2,1) = 0.19;
	ok = ok && Compare(p_ronc_mod, s_ronc_mod);
	s_ronc_mod.AssignAt(2,1, 0.79);
	p_ronc_mod.At(1,2) = 0.79;
	p_ronc_mod.At(2,1) = 0.79;
	ok = ok && Compare(p_ronc_mod, s_ronc_mod);

    SKTRAN_ScatMat_MIMSNC s_ronc_eye;
	s_ronc_eye.SetToIdentity();
	ok = ok && Compare( p_eye, s_ronc_eye );
	s_ronc_mod = s_ronc;
	s_ronc_mod.AddToThis( s_ronc_eye, 0.019 );
	ok = ok && Compare( p_ronc+(p_eye*0.019), s_ronc_mod );

	s_ronc_mod = s_ronc;
	p_ronc_mod = p_ronc;
	s_ronc_mod += s_ronc_eye;
	p_ronc_mod += p_eye;
	ok = ok && Compare( p_ronc_mod, s_ronc_mod );
	s_ronc_mod *= 0.19;
	p_ronc_mod *= 0.19;
	ok = ok && Compare( p_ronc_mod, s_ronc_mod );

	ok = ok && Compare( p_ronc_mod*v, s_ronc_mod*v );
	ok = ok && Compare( p_ronc_mod*0.19, s_ronc_mod*0.19 );
	ok = ok && Compare( p_ronc_mod-p_eye, s_ronc_mod-s_ronc_eye );
	ok = ok && Compare( p_ronc_mod+p_eye, s_ronc_mod+s_ronc_eye );

	return ok;
}

bool TestScatterMatrices::TestRot ( ) const{

	bool ok = true;

	SKTRAN_Stokes_NC v;
	v.SetTo( 0.79, 0.53, 0.43, 0.0 );
	
	const double theta = 2.3;
	const double cosTheta = std::cos( theta );
	const double sinTheta = std::sin( theta );

	skRTPhaseMatrix p;
	p.SetTo(0.0);
	p.At(1,1) =  1.0;
	p.At(2,2) =  1.0;
	p.At(3,3) =  1.0;
	SKTRAN_ScatMat_Rot s;
	ok = ok && Compare( p*v, s*v );

	p.At(2,2) =  cosTheta;
	p.At(3,3) =  cosTheta;
	p.At(2,3) = -sinTheta;
	p.At(3,2) =  sinTheta;
	s.SetAngle(cosTheta, sinTheta);
	ok = ok && Compare( p*v, s*v );

	return ok;
}


bool TestScatterMatrices::TestRONC_Composite ( ) const
{
	bool ok = true;

	skRTPhaseMatrix p_zero;
	p_zero.SetTo( 0.0 );

	SKTRAN_Stokes_NC v;
	v.SetTo( 0.79, 0.53, 0.43, 0.0 );

	skRTPhaseMatrix p;
	p.SetTo(0.0);
	p.At(1,1) = 1.00;
	p.At(1,2) = 0.19;
	p.At(1,3) = 0.17;
	p.At(2,1) = 0.13;
	p.At(2,2) = 0.11;
	p.At(2,3) = 0.07;
	p.At(3,1) = 0.05;
	p.At(3,2) = 0.03;
	p.At(3,3) = 0.02;

	skRTPhaseMatrix p_ronc;
	p_ronc.SetTo( 0.0 );
	p_ronc.At(1,1) = 0.37;
	p_ronc.At(1,2) = 0.31;
	p_ronc.At(2,1) = 0.31;
	p_ronc.At(2,2) = 0.29;
	p_ronc.At(3,3) = 0.23;
    SKTRAN_ScatMat_MIMSNC s_ronc(p_ronc);
	
	const double theta = 2.1;
	SKTRAN_ScatMat_Rot s_rot;
	s_rot.SetAngle( std::cos(theta), std::sin(theta) );
	skRTPhaseMatrix p_rot;
	p_rot.SetTo(0.0);
	p_rot.At(1,1) =  1.0;
	p_rot.At(2,2) =  std::cos(2.0*theta);
	p_rot.At(3,3) =  std::cos(2.0*theta);
	p_rot.At(2,3) = -std::sin(2.0*theta);
	p_rot.At(3,2) =  std::sin(2.0*theta);

	SKTRAN_PhaseMat_MIMSNC s(p);
	SKTRAN_PhaseMat_MIMSNC s_zero;
	ok = ok && Compare( p_zero, s_zero );
	ok = ok && Compare( p, s );
	ok = ok && Compare( p_rot, SKTRAN_PhaseMat_MIMSNC(s_rot) );
	ok = ok && Compare( p_ronc, SKTRAN_PhaseMat_MIMSNC(p_ronc) );
	ok = ok && Compare( s, SKTRAN_PhaseMat_MIMSNC(s) );

	for( int xidx=1; xidx<4; ++xidx ){
		for( int yidx=1; yidx<4; ++yidx ){
			skRTPhaseMatrix p_mod( p_zero );
			SKTRAN_PhaseMat_MIMSNC s_mod( s_zero );
			p_mod.At( xidx, yidx ) = 0.19;
			s_mod.AssignAt( xidx, yidx, 0.19 );
			ok = ok && Compare( p_mod, s_mod );
		}
	}
	
	skRTPhaseMatrix p_oneval; 
	p_oneval.SetTo( 0.19 );
	for( int xidx = 1; xidx<5; ++xidx ){
		p_oneval.At(xidx, 4) = 0.0;
		p_oneval.At(4, xidx) = 0.0;
	}
	SKTRAN_PhaseMat_MIMSNC s_oneval;
	s_oneval.SetTo( 0.19 );
	ok = ok && Compare( p_oneval, s_oneval );

	skRTPhaseMatrix p_eye;
	p_eye.SetTo( 0.0 );
	p_eye.At(1,1) = 1.0;
	p_eye.At(2,2) = 1.0;
	p_eye.At(3,3) = 1.0;
	SKTRAN_PhaseMat_MIMSNC s_eye;
	s_eye.SetToIdentity();
	ok = ok && Compare( p_eye, s_eye );

	SKTRAN_PhaseMat_MIMSNC s_mod;
	skRTPhaseMatrix p_mod;
	s_mod = s; p_mod = p;
	s_mod.RMultBy( SKTRAN_ScatMat_Rot( s_rot ) );
	p_mod.RMultBy( p_rot );
	ok = ok && Compare( p_mod, s_mod );
	s_mod = s; p_mod = p;
	s_mod.LMultBy( SKTRAN_ScatMat_Rot( s_rot ) );
	p_mod.LMultBy( p_rot );
	ok = ok && Compare( p_mod, s_mod );

	s_mod = s; p_mod = p;
	s_mod.RMultBy( SKTRAN_ScatMat_MIMSNC( s_ronc ) );
	p_mod.RMultBy( p_ronc );
	ok = ok && Compare( p_mod, s_mod );
	s_mod = s; p_mod = p;
	s_mod.LMultBy( SKTRAN_ScatMat_MIMSNC( s_ronc ) );
	p_mod.LMultBy( p_ronc );
	ok = ok && Compare( p_mod, s_mod );

	s_mod = s; p_mod = p;
	s_mod.RMultBy( SKTRAN_PhaseMat_MIMSNC( s ) );
	p_mod.RMultBy( p );
	ok = ok && Compare( p_mod, s_mod );
	s_mod = s; p_mod = p;
	s_mod.LMultBy( SKTRAN_PhaseMat_MIMSNC( s ) );
	p_mod.LMultBy( p );
	ok = ok && Compare( p_mod, s_mod );

	s_mod  = s;
	s_mod += s;
	p_mod  = p;
	p_mod += p;
	ok = ok && Compare( p_mod, s_mod );
	s_mod *= 0.19;
	p_mod *= 0.19;
	ok = ok && Compare( p_mod, s_mod );

	s_mod = s_rot;
	ok = ok && Compare( p_rot, s_mod );
	s_mod = s_ronc;
	ok = ok && Compare( p_ronc, s_mod );

	ok = ok && Compare( p*v, s*v );
	ok = ok && Compare( p * 0.19, s * 0.19 );

	ok = ok && Compare( p - p_eye, s - s_eye );
	ok = ok && Compare( p + p_eye, s + s_eye );

	return ok;
}

bool TestScatterMatrices::RunTest ( ) const
{	
	bool ok = true;

	ok = ok && TestRONC           ( );
	ok = ok && TestRot            ( );
	ok = ok && TestRONC_Composite ( );
	
    return ok;
}

