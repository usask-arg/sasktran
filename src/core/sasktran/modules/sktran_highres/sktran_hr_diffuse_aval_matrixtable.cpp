#include "include/sktran_hr_internals.h"


bool Avals_MatrixTable::AllocateStorage( size_t incomingcounter, size_t outgoingcounter, size_t scattvalcounter, size_t pointscounter )
{
    m_interpFunctions.resize( pointscounter );
    m_incomingWeightMultipliers .resize( incomingcounter     );

    return true;
}


bool Avals_MatrixTable::ComputeMultipliersAndAdjustScatterArray ( const SKTRAN_HR_Diffuse_Point& point )
{
    bool						ok = true;
    double						w;
    size_t						inidx;
    size_t						outidx;
    SKTRAN_StokesScalar			ratio;
    const SKTRAN_UnitSphere_V2* incomingsphere;
    size_t						numinscats;
    size_t						numoutscats;
    size_t                      inrayidx;
    size_t                      outrayidx;
    double						scatextinction;
    SKTRAN_ScatMat_MIMSNC         p;

    incomingsphere = point.IncomingUnitSphere();
    numinscats     = point.NumUniqueScatterIncoming();
    numoutscats    = point.NumUniqueScatterOutgoing();

    scatextinction = m_opticaltable->ScatteringExtinctionPerCM( point.Location() ) * 100;

    nxVector outray;
    for( inidx = 0; inidx < numinscats; inidx++ )
    {
        ratio = 0.0;
        inrayidx = point.UniqueScatterIncoming( inidx );
        nxVector inray(point.IncomingUnitRayLocalCoords(inrayidx));
        for( outidx = 0; outidx < numoutscats; outidx++ )
        {
            outrayidx = point.UniqueScatterOutgoing( outidx );
            w		= point.OutgoingCubatureWeight( outrayidx );
            point.OutgoingRayLocalCoords(outrayidx, outray);
            m_interpFunctions[point.UniquePointIdentifier()].Interpolate_kscattPerM( -(inray&outray), p);
            ratio += p.At(1,1)*w;
        }
        //ratio /= incomingWeight; // Didn't multiply this onto Avals, so don't need to remove it
        ratio *= scatextinction>0 ? 1/scatextinction : 0.0;

        m_incomingWeightMultipliers[point.IncomingRadianceIdx(inidx)] = (SKTRAN_HR_DBL_TO_WEIGHT(ratio)>1e-20)?SKTRAN_HR_DBL_TO_WEIGHT(ratio)*point.InCubatureWeight(inidx):0.0;
        if ( (ratio <0.177827941003892) || (ratio > 5.62341325190349) || (ratio!=ratio) )  // if ( (ratio < pow(10,-0.75)) || (ratio > pow(10,0.75)) )
        {
			static int counter = 0;
			if (counter < 25)
			{
				nxLog::Record(NXLOG_WARNING, "   Large Phase Function Adjustment: %7.4e  %7.4e  %7.4e", ratio, incomingsphere->CubatureWeightAt( inrayidx ), scatextinction);
				counter++;
			}
        }
    }

    return ok;
}


bool Avals_MatrixTable::CalcScatteringMatrixPoint( const SKTRAN_HR_Diffuse_Point& point )
{
    bool ok = true;

    if( ! point.IsGroundPoint() ){
        ok = ok && CalcScatteringMatrixPoint_Participating ( point );
    } else{
        ok = ok && CalcScatteringMatrixPoint_Boundary( point );
    }

    return ok;
}


bool Avals_MatrixTable::CalcScatteringMatrixPoint_Participating ( const SKTRAN_HR_Diffuse_Point& point )
{
    bool ok = true;
    
    ok = ok && m_opticaltable->CreateInterpolationForPoint( point.Location(), m_interpFunctions[point.UniquePointIdentifier()] );
    ok = ok && ComputeMultipliersAndAdjustScatterArray( point );

    return ok;
}


bool Avals_MatrixTable::CalcScatteringMatrixPoint_Boundary ( const SKTRAN_HR_Diffuse_Point& point )
{
    bool ok = true;
    size_t numinscats		= point.NumUniqueScatterIncoming();
    HELIODETIC_UNITVECTOR	ingrounddir;
    SKTRAN_ScatMat_MIMSNC     p;

    p.SetTo(0.0);

    std::vector< SKTRAN_ScatMat_MIMSNC > pmatrices;
    SKTRAN_GridDefScatterAngle_V21* angles = new SKTRAN_GridDefScatterAngle_V21;
    const size_t numAngles = 201;
    angles->Configure( 180.0 / (numAngles-1) );
    angles->SetGridSearchMode( SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM );
    pmatrices.reserve( angles->NumAngles() );
    p.SetTo( 0.0 );
    for (int aidx = 0; aidx < angles->NumAngles(); ++aidx) {
        p.AssignAt(1,1,angles->At(aidx));
        pmatrices.push_back( p );
    }

    double albedo;
    m_opticaltable->Get_AlbedoForDeprecatedLegacyCode( point.Location(), &albedo );
    m_interpFunctions[point.UniquePointIdentifier()].InjectTable( pmatrices, angles, albedo, point.Location() );

    // Need to multiply by point.InCubatureWeight ( point.UniqueScatterIncoming( inidx );)
    for(int inidx=0; inidx<point.NumIncomingRays(); ++inidx)
        m_incomingWeightMultipliers[point.IncomingRadianceIdx(inidx)] = point.InCubatureWeight(inidx);

    return ok;
}


void Avals_MatrixTable::GetValue ( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_ScatMat_MIMSNC*       Aval ) const 
{
    m_interpFunctions[point.UniquePointIdentifier()].Interpolate_kscattPerM( point.CosScatteringAngle(outidx, inidx), *Aval );
    *Aval *= m_incomingWeightMultipliers[point.IncomingRadianceIdx(inidx)];
}


void Avals_MatrixTable::ApplyValue ( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_Stokes_NC* radiance ) const 
{
	SKTRAN_ScatMat_MIMSNC s; 
	GetValue( point, inidx, outidx, &s );
	s.LApplyTo( radiance );
}


void Avals_MatrixTable::ApplyValue ( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_HR_WEIGHT_TYPE* radiance ) const 
{
	SKTRAN_ScatMat_MIMSNC s; 
	GetValue( point, inidx, outidx, &s );
	*radiance *= (SKTRAN_HR_WEIGHT_TYPE)(s.p11());
}


std::unique_ptr< Aval_ScalarIteratorManager_Base > Avals_MatrixTable::ScalarAvalIteratorManager ( const SKTRAN_HR_Diffuse_Point& point ) const
{
  
    SKTRAN_ScatMat_MIMSNC aval;

    const size_t numincoming ( point.NumIncomingRays( ) );
    const size_t numoutgoing ( point.NumOutGoingRays( ) );

    std::unique_ptr< Aval_CScalarIteratorManager > manager( new Aval_CScalarIteratorManager(numincoming,numoutgoing) );
	std::vector< SKTRAN_HR_WEIGHT_TYPE >::iterator iterator ( manager->m_vals.begin() );

    // By convention, scalar Avals are accessed in "incomingRayIndex-major" order. 
    // If this changes, it needs to be changed in all the iterator manager clases 
    // and in the scalar radstore classes
    for( size_t inidx=0; inidx<numincoming; ++inidx ){
        for( size_t outidx=0; outidx<numoutgoing; ++outidx ){
            GetValue( point, inidx, outidx, &aval );
            *iterator = (SKTRAN_HR_WEIGHT_TYPE) aval.p11();
            ++iterator;
        }
    }
    
    return std::unique_ptr< Aval_ScalarIteratorManager_Base > ( std::move(manager) );
}


