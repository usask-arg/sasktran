#include "include/sktran_hr_internals.h"


void EtaCalculator_DoRotation::CalculateEtas ( const SKTRAN_HR_Diffuse_Point& point, size_t inidx ) 
{
	size_t scatidx = sub2ind( inidx, m_outidx, point.NumIncomingRays() );

	// Compute rotation matrices from spherical geometry
	const HELIODETIC_UNITVECTOR inray(point.IncomingRayGlobalCoords(inidx));
	double nol1 = -point.IncomingRayGlobalCoords( inidx ).Z(); // "(ref dir [n]= +sun dir), dot prop_1=-ell_1"
	double l1l2 = -(inray & m_outray ); // prop_1=-ell_1, dot prop_2=ell_2 
	double invSinScattAngle ( (1.0-l1l2*l1l2)>m_denomtol ? 1.0 / sqrt(1.0-l1l2*l1l2) : 1.0 );
	double sinDeltaLon = fabs(m_incoming_inPlane[inidx].x*m_outgoing_inPlane.y - m_incoming_inPlane[inidx].y*m_outgoing_inPlane.x);
	double sinRuleRatio = sinDeltaLon*invSinScattAngle;

	// Use spherical law of cosines to get "rotation from reference plane into scattering plane" angle
	if(  1e-10 < m_sinInToSun[inidx] ){
		m_crts = (m_nol2 - nol1*l1l2) * invSinScattAngle * m_invSinInToSun[inidx];
		m_srts = m_sinOutToSun*sinRuleRatio;
	} else{
		m_crts = 1.0;
		m_srts = 0.0;
	}

	// and "rotation from scattering plane back to reference plane" angle
	if( 1e-10 < m_sinOutToSun ){
		m_cstr = (nol1 - m_nol2*l1l2) * invSinScattAngle * m_invSinOutToSun;
		m_sstr = -m_sinInToSun[inidx]*sinRuleRatio; // Picks up negative sign, Mishchenko p. 88
	} else{
		m_cstr = 1.0;
		m_sstr = 0.0;
	}

	// One rotation will be clockwise, the other CCW, unless in and out directions are at same "local latitude", in which case eta=0 anyway
	double discr = -inray.X()*m_outray.Y() + inray.Y()*m_outray.X();
	if( 0.0 < discr ){
		m_crts = -m_crts;
	} else{
		m_cstr = -m_cstr;
	}

}


void EtaCalculator_DoRotation::UpdateOutgoingIndex( const SKTRAN_HR_Diffuse_Point& point, size_t outidx ) 
{
	const double dtol = 1.0; //-1e-10;
	const double denomtol = 1e-10;

	m_outidx = outidx;

	nxVector outray_unrotated;
	point.OutgoingRayLocalCoords( outidx, outray_unrotated );
	m_outray = (HELIODETIC_VECTOR(m_qPtVecs[0],outray_unrotated.X()) + HELIODETIC_VECTOR(m_qPtVecs[1],outray_unrotated.Y()) + HELIODETIC_VECTOR(m_qPtVecs[2],outray_unrotated.Z())).UnitVector(); // Get local incoming LOOK direction

	m_nol2 = m_outray.Z(); // "prop_2=ell_2, dot (ref dir [n]= +sun dir)"
	m_sinOutToSun = fabs(m_nol2)<dtol ? sqrt( 1.0 - m_nol2*m_nol2 ) : 0.0;
	m_invSinOutToSun = m_sinOutToSun>denomtol ? 1.0/m_sinOutToSun : 0.0;

	double norm = fabs(m_sinOutToSun)>1e-20 ? 1.0/m_sinOutToSun : 0.0;
	m_outgoing_inPlane = inPlaneVec{ m_outray.X()*norm, m_outray.Y()*norm };
}

void EtaCalculator_DoRotation::SetOutgoingDirection(const HELIODETIC_UNITVECTOR &out) {
    const double dtol = 1.0; //-1e-10;
    const double denomtol = 1e-10;

    m_outray = out;

    m_nol2 = m_outray.Z(); // "prop_2=ell_2, dot (ref dir [n]= +sun dir)"
    m_sinOutToSun = fabs(m_nol2)<dtol ? sqrt( 1.0 - m_nol2*m_nol2 ) : 0.0;
    m_invSinOutToSun = m_sinOutToSun>denomtol ? 1.0/m_sinOutToSun : 0.0;

    double norm = fabs(m_sinOutToSun)>1e-20 ? 1.0/m_sinOutToSun : 0.0;
    m_outgoing_inPlane = inPlaneVec{ m_outray.X()*norm, m_outray.Y()*norm };
}


void EtaCalculator_DoRotation::SetPoint( const SKTRAN_HR_Diffuse_Point& point ) 
{
	const size_t          numincoming = point.NumIncomingRays();

	const double dtol = 1.0; //-1e-10;
	const double denomtol = 1e-10;

	m_incoming_inPlane .resize ( numincoming ); 
	m_sinInToSun       .resize( numincoming );
	m_invSinInToSun    .resize( numincoming );

	// Unit vectors for the "local" reference frame
	point.Location().LocalUnitVectors( m_qPtVecs, 3 ); 

	// Pre-calculate (some) eta terms for stokes vector rotation
	for( int inidx=0; inidx<numincoming; ++inidx ){
		const HELIODETIC_UNITVECTOR& in = point.IncomingRayGlobalCoords( inidx );
		double nol1 = -in.Z(); // "(ref dir [n]= +sun dir), dot prop_1=-ell_1"
		m_sinInToSun[inidx] = fabs(nol1)<dtol ? sqrt( 1.0 - nol1*nol1 ) : 0.0;
		m_invSinInToSun[inidx] = m_sinInToSun[inidx]>denomtol ? 1.0/m_sinInToSun[inidx] : 0.0;
		double norm = fabs(m_sinInToSun[inidx])>1e-20 ? 1.0/m_sinInToSun[inidx] : 0.0;
		m_incoming_inPlane[inidx] = inPlaneVec{ -in.X()*norm, -in.Y()*norm };
	}
}


void EtaCalculator_DoRotation::RefToScatt( SKTRAN_Stokes_NC& stokes ) const 
{
	//stokes.RotatePolarPlaneThru( m_currentEtas.cosEta_refToScatt, m_currentEtas.sinEta_refToScatt ); // Rotate first order incoming into scatter plane
	stokes.RotatePolarPlaneThru( m_crts, m_srts ); // Rotate first order incoming into scatter plane
};


void EtaCalculator_DoRotation::ScattToRef( SKTRAN_Stokes_NC& stokes ) const 
{
	//stokes.RotatePolarPlaneThru( m_currentEtas.cosEta_scattToRef, m_currentEtas.sinEta_scattToRef ); // Rotate into reference plane
	stokes.RotatePolarPlaneThru( m_cstr, m_sstr ); // Rotate into reference plane
};


void EtaCalculator_DoRotation::ScattMatToPhaseMat  ( const SKTRAN_ScatMat_MIMSNC& s, SKTRAN_PhaseMat_MIMSNC& p ) const
{
	SKTRAN_ScatMat_Rot r;

	// Load scatter matrix
	p = s;

	// Rotate to and from scattering frame 
	r.SetAngle( m_crts, m_srts );
	p.RMultBy( r );
	r.SetAngle( m_cstr, m_sstr );
	p.LMultBy( r );
};


