#include "include/sktran_hr_internals.h"
#include <numeric>
#include <functional>


size_t RadStore_PV1::GetDiffuseSourceOrder ( ) const
{
	return m_order_incoming;
}


bool RadStore_PV1::DiffuseSource ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& lookvec, double& source ) const  
{
	bool ok = true;
	
	SKTRAN_Stokes_NC vec;
	ok = ok && DiffuseSource( qNoTransform, point, lookvec, vec );
	source = vec.I();

	return ok;
}


bool RadStore_PV1::GetPolarizedComponent_basisImplicit ( size_t radIndex, const nxVector& look_transformed, SKTRAN_Stokes_NC& incomingvec, const HELIODETIC_POINT& queryPoint, const HELIODETIC_UNITVECTOR& inrayQCoords ) const
{
    incomingvec.SetTo( 0.0 );
    return true;
}


bool RadStore_PV1::DiffuseSource ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& negLookUnit_local, SKTRAN_Stokes_NC& source ) const  
{
    bool ok = true;

    const double           dtol(1e-10);  // Numerical error tolerance (double)
    nxVector               inray_local;
    SKTRAN_HR_WEIGHT_TYPE  incomingscalar;
    SKTRAN_Stokes_NC       incomingvec;
    SKTRAN_Stokes_NC       pointvec;
    SKTRAN_ScatMat_MIMSNC  scatMat;
    double                 cosEta_refToScatt, sinEta_refToScatt;
    double                 cosEta_scattToRef, sinEta_scattToRef;
    size_t                 radIndex;

    const size_t					numincoming = point.NumIncomingRays();

    // Fetch query point basis
    const HELIODETIC_BASIS& qPt_obsBasis_helio = qNoTransform.GetBasis();

    // Be able to get scattering frame at querry point
    const size_t numQPtVecs = 3;
    HELIODETIC_UNITVECTOR qPtVecs[numQPtVecs];
    qNoTransform.GetPoint().LocalUnitVectors(qPtVecs, numQPtVecs); 

    double nol2 = -qNoTransform.GetLookAway().Z(); // "prop_2=-ell_2, dot (ref dir [n] = +sun dir)"

    SKTRAN_TableOpticalProperties_PointCache scatMatInterpolator;
    ok = ok && m_opticaltable->CreateInterpolationForPoint( qNoTransform.GetPoint(), scatMatInterpolator );

    double p11running(0.0); // Also integrate the phase function so we can do photon conservation 

    pointvec.SetTo(0.0);
    for( size_t inidx = 0; inidx < numincoming; inidx++ )
    {
        radIndex = point.IncomingRadianceIdx( inidx );
        inray_local = point.IncomingUnitRayLocalCoords( inidx );

        // Rotate incomingray to qPt
        HELIODETIC_VECTOR qPtInrayHelio_notUnit = HELIODETIC_VECTOR(qPtVecs[0],inray_local.X()) + HELIODETIC_VECTOR(qPtVecs[1],inray_local.Y()) + HELIODETIC_VECTOR(qPtVecs[2],inray_local.Z()); // Get local incoming LOOK direction
        HELIODETIC_UNITVECTOR qPtInrayHelio;
        qPtInrayHelio.SetCoords( qPtInrayHelio_notUnit.X(), qPtInrayHelio_notUnit.Y(), qPtInrayHelio_notUnit.Z() );

        // Get incoming polarized radiance, defined for a basis with y in the incomingRay-observerRay (scattering) plane
        ok = ok && GetPolarizedComponent_basisImplicit( radIndex, negLookUnit_local, incomingvec, qNoTransform.GetPoint(), qPtInrayHelio );
        
        // Get scalar incoming from this direction and add as unpolarized 
        incomingscalar = m_totalincomingscalar[radIndex]; 
        incomingvec.Add_I( incomingscalar ); 
        
        // Rotate from solar (implicit) frame into scattering frame
        // Use spherical law of cosines to calculate rotation angles
        double nol1 = -qPtInrayHelio.Z(); // "prop1=-ell1, dot (ref dir [n] = +sun dir"
        double l1l2 = qPtInrayHelio & qNoTransform.GetLookAway(); // "prop1=-ell1, dot prop2=-ell2"
        double sinInToSun    = sqrt( 1.0 - nol1*nol1 );
        double sinOutToSun   = sqrt( 1.0 - nol2*nol2 );
        double sinScattAngle = sqrt( 1.0 - l1l2*l1l2 );
        if( (dtol < sinInToSun)  &&  (dtol < sinOutToSun)  &&  (dtol < sinScattAngle) ){
            cosEta_refToScatt = (nol2 - nol1*l1l2) / (sinInToSun*sinScattAngle);
            double sintemp = 1.0 - cosEta_refToScatt*cosEta_refToScatt;
            sinEta_refToScatt = 0.0 < sintemp ? sqrt( sintemp ) : 0.0;
            cosEta_scattToRef = (nol1 - nol2*l1l2) / (sinOutToSun*sinScattAngle);
            sintemp = 1.0 - cosEta_scattToRef*cosEta_scattToRef;
            sinEta_scattToRef = 0.0 < sintemp ? sqrt( sintemp ) : 0.0;
            double discr = qPtInrayHelio.X()*qNoTransform.GetLookAway().Y() - qPtInrayHelio.Y()*qNoTransform.GetLookAway().X(); // cross(prop1=-ell1, prop2=-ell2), dot (ref dir [n] = +sin dir)

            if( 0.0 < discr ){
                cosEta_refToScatt = -cosEta_refToScatt;
                sinEta_scattToRef = -sinEta_scattToRef;
            } else{
                cosEta_scattToRef = -cosEta_scattToRef;
                sinEta_scattToRef = -sinEta_scattToRef;
            }
        } else{
            cosEta_refToScatt = 1.0;
            sinEta_refToScatt = 0.0;
            cosEta_scattToRef = 1.0;
            sinEta_scattToRef = 0.0;
        }

        // Rotate the basis for the electric field components used in this vector through angle eta (Tony's thesis page 44), rotates polarization plane CW about look dir
        incomingvec.RotatePolarPlaneThru( cosEta_refToScatt, sinEta_refToScatt ); // Result of pmatrix acting on source vector is in scatter frame -- rotate polarization into my frame

        // Do the scatter (multiply the Mueller matrix)
        scatMatInterpolator.Interpolate_kscattPerM( l1l2, scatMat );
        incomingvec.LMultBy( scatMat );

        // Rotate the basis for the electric field components used in this vector through angle eta (Tony's thesis page 44), rotates polarization plane CW about look dir
        incomingvec.RotatePolarPlaneThru( cosEta_scattToRef, sinEta_scattToRef ); // Result of pmatrix acting on source vector is in scatter frame -- rotate polarization into my frame

        // Add contribution from this inoming direction to the total contribution at this point
        incomingvec *= point.InCubatureWeight( inidx );
        p11running += scatMat.p11() * point.InCubatureWeight( inidx ); // Add to the numerical integral of the phase matrix 
        pointvec += incomingvec;
    }
    
    // Photon conservation. \int_{4\pi} d\Omega p_{11}(\Omega\cdot\textrm{look}) = 1
    // Normalize by the numerically calculated integral of the phase function to make this true
    p11running = 0<p11running ? p11running : 1.0; 
    pointvec *= 100.0 * m_opticaltable->ScatteringExtinctionPerCM( qNoTransform.GetPoint() ) / p11running;

    // Need to do this stupid step because SKTRAN_HR_Diffuse_Table_CPU::DiffuseSource_impl assumes the output has been multiplied by the scattering
    // coefficient at the *diffuse point* (this behaviour is grandfathered in from the scalar code, which interpolates outgoing radiance on an 
    // outgoing sphere, whereas here we scatter directly from the incoming sphere to the line of sight and look up the phase matrix at the query 
    // point directly from the optical properties table. 
    double extinctionRatio = m_opticaltable->ScatteringExtinctionPerCM( point.Location() ) / m_opticaltable->ScatteringExtinctionPerCM( qNoTransform.GetPoint() );
    pointvec *= isfinite(extinctionRatio) ? extinctionRatio : 0.0;

    source = pointvec;

    return ok;
}

