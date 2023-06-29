#include "include/sktran_hr_internals.h"


template< class MatrixType > 
bool Avals_MatrixStore< MatrixType >::AllocateStorage( size_t incomingcounter, size_t outgoingcounter, size_t scattvalcounter, size_t pointscounter )
{
    m_Avals   .resize( scattvalcounter     );
 //   m_Acolind .resize( scattvalcounter     );
 //   m_Arowptr .resize( outgoingcounter + 1 );

 //   m_Arowptr[m_Arowptr.size()-1] = static_cast<int>(m_Avals.size());

    return true;
}


template< class MatrixType > 
bool Avals_MatrixStore< MatrixType >::ComputeMultipliersAndAdjustScatterArray ( const SKTRAN_HR_Diffuse_Point& point )
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

    incomingsphere = point.IncomingUnitSphere();
    numinscats     = point.NumUniqueScatterIncoming();
    numoutscats    = point.NumUniqueScatterOutgoing();

    scatextinction = m_opticaltable->ScatteringExtinctionPerCM( point.Location() ) * 100;

    for( inidx = 0; inidx < numinscats; inidx++ )
    {
        ratio = 0.0;
        inrayidx = point.UniqueScatterIncoming( inidx );
        for( outidx = 0; outidx < numoutscats; outidx++ )
        {
            outrayidx = point.UniqueScatterOutgoing( outidx );
            w		= point.OutgoingCubatureWeight( outrayidx );
            ratio  += m_Avals[ point.ScatterPropertyIdx(outrayidx, inrayidx) ].At(1,1)*w;
        }
        ratio /= incomingsphere->CubatureWeightAt( inrayidx );
        ratio /= scatextinction;
        for( outidx = 0; outidx < numoutscats; outidx++ )
        {
            outrayidx = point.UniqueScatterOutgoing( outidx );
            m_Avals[ point.ScatterPropertyIdx(outrayidx, inrayidx) ] *= (SKTRAN_HR_DBL_TO_WEIGHT(ratio)>1e-20)?1.0/SKTRAN_HR_DBL_TO_WEIGHT(ratio):0.0;
            //if(outidx==3 && point.UniquePointIdentifier()==160){
            //    printf("\n0 %i %20.17e  %20.17e  %20.17e  %20.17e  0.0", (int)inidx, m_Avals[ point.ScatterPropertyIdx(outrayidx, inrayidx) ].At(1,1), m_Avals[ point.ScatterPropertyIdx(outrayidx, inrayidx) ].At(2,1), m_Avals[ point.ScatterPropertyIdx(outrayidx, inrayidx) ].At(2,2), m_Avals[ point.ScatterPropertyIdx(outrayidx, inrayidx) ].At(3,3) );
            //}
        }
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


template< class MatrixType > 
bool Avals_MatrixStore< MatrixType >::CalcScatteringMatrixPoint( const SKTRAN_HR_Diffuse_Point& point )
{
    bool ok = true;

    if( ! point.IsGroundPoint() ){
        ok = ok && CalcScatteringMatrixPoint_Participating ( point );
    } else{
        ok = ok && CalcScatteringMatrixPoint_Boundary( point );
    }

    return ok;
}


template< class MatrixType > 
bool Avals_MatrixStore< MatrixType >::CalcScatteringMatrixPoint_Participating ( const SKTRAN_HR_Diffuse_Point& point )
{
    bool ok = true;

    size_t     numinscats  ( point.NumUniqueScatterIncoming() );
    size_t     numoutscats ( point.NumUniqueScatterOutgoing() );
    size_t     outrayidx;
    size_t     inrayidx;

    nxVector   outdir;
    nxVector   indir;
	MatrixType s;

	std::unique_ptr< EtaCalculator_Base > etaCalculator;

	// Create object to calculate rotation angles
	CreateEtaCalculatorForStoragePhase( etaCalculator ); 
	etaCalculator->SetPoint( point );

    for( size_t outidx = 0; outidx < numoutscats; outidx++ )
    {
        outrayidx = point.UniqueScatterOutgoing( outidx );
        point.OutgoingRayLocalCoords( outrayidx, outdir );
		etaCalculator->UpdateOutgoingIndex( point, outidx );
        for( size_t inidx = 0; inidx < numinscats; inidx++ )
        {
            inrayidx = point.UniqueScatterIncoming( inidx );
            indir = point.IncomingUnitRayLocalCoords( inrayidx );
			etaCalculator->CalculateEtas( point, inidx );
            ok = ok && CalcScatteringMatrixPoint_Participating_Impl( point, indir, outdir, point.InCubatureWeight(inidx), s, etaCalculator );
            m_Avals  [ point.ScatterPropertyIdx(outrayidx, inrayidx) ] = s;
 //           m_Acolind[ point.ScatterPropertyIdx(outrayidx, inrayidx) ] = static_cast<int>(point.IncomingRadianceIdx(inrayidx));

        }
//        m_Arowptr[ point.OutgoingRadianceIdx(outrayidx) ] = static_cast<int>( point.ScatterPropertyIdx(outrayidx, 0) );
    }
    ok = ok && ComputeMultipliersAndAdjustScatterArray( point );

    return ok;
}


template< > 
bool Avals_MatrixStore< SKTRAN_ScatMat_MIMSNC >::CalcScatteringMatrixPoint_Participating_Impl( const SKTRAN_HR_Diffuse_Point& point, const nxVector& indir, const nxVector& outdir, double cubatureWeight, SKTRAN_ScatMat_MIMSNC& s, std::unique_ptr< EtaCalculator_Base >& )
{
	bool ok = true;

	// Don't need etaCalculator, could change implementation so it doesn't have to get passed in, but the way it is now makes multithreading easy
	double cosangle = -1 * (indir.X()*outdir.X() +
		indir.Y()*outdir.Y() +
		indir.Z()*outdir.Z());
	ok = ok && m_opticaltable->GetScatteringMatrixCM2(point.Location(), cosangle, s); 
	s *= (100.0*cubatureWeight);

	return ok;
}


template< > 
bool Avals_MatrixStore< SKTRAN_PhaseMat_MIMSNC >::CalcScatteringMatrixPoint_Participating_Impl( const SKTRAN_HR_Diffuse_Point& point, const nxVector& indir, const nxVector& outdir, double cubatureWeight, SKTRAN_PhaseMat_MIMSNC& p, std::unique_ptr< EtaCalculator_Base >& etaCalculator )
{
	bool ok = true;
	SKTRAN_ScatMat_MIMSNC s;

	double cosangle = -1 * (indir.X()*outdir.X() +
		indir.Y()*outdir.Y() +
		indir.Z()*outdir.Z());
	ok = ok && m_opticaltable->GetScatteringMatrixCM2(point.Location(), cosangle, s); 
	s *= (100.0*cubatureWeight);

	// Now take Eta from eta calculator and apply it to s to make p
	etaCalculator->ScattMatToPhaseMat( s, p );
		
	return ok;
}


template< class MatrixType > 
bool Avals_MatrixStore< MatrixType >::CalcScatteringMatrixPoint_Boundary ( const SKTRAN_HR_Diffuse_Point& point )
{
    bool ok = true;

    size_t numinscats		= point.NumUniqueScatterIncoming();

    HELIODETIC_UNITVECTOR	ingrounddir;
    double					cosangle;
    double					scattcoeff;
    SKTRAN_ScatMat_MIMSNC   p;
    size_t                  inrayidx;

    p.SetTo(0.0);
    // ground point
    // note that we dont take into account albedo in the scattering matrix
    // this is done so that we can use an albedo that varies based on position
    for( size_t inidx = 0; inidx < numinscats; inidx++ )
    {
        inrayidx = point.UniqueScatterIncoming( inidx );
        ingrounddir = point.IncomingRayGlobalCoords( inrayidx );
        cosangle = fabs(ingrounddir & point.Location().LocalZenith() );
        //printf("\ns: %e  %e", cosangle, ingrounddir & point.Location().LocalZenith());

        scattcoeff = cosangle;
        scattcoeff *= point.InCubatureWeight( inrayidx );
        p.AssignAt(1,1,scattcoeff); // Here is assumption that boundary is lambertian! 
        m_Avals  [ point.GroundDownwardFluxFactorsIdx( inrayidx) ] = p;
//        m_Acolind[ point.ScatterPropertyIdx( 0, inrayidx ) ] = static_cast<int>( point.IncomingRadianceIdx(inrayidx) );
    }
 //   m_Arowptr[ point.OutgoingRadianceIdx(0) ] = static_cast<int>( point.ScatterPropertyIdx(0,0) );

    return ok;
}

template< class MatrixType >
const MatrixType&	Avals_MatrixStore< MatrixType >::GetValueRef( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx ) const 
{
	return point.IsGroundPoint() ? m_Avals[ point.GroundDownwardFluxFactorsIdx(inidx) ] :  m_Avals[ point.ScatterPropertyIdx(outidx, inidx) ];
}


template< class MatrixType > 
void Avals_MatrixStore< MatrixType >::ApplyValue ( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_Stokes_NC* radiance ) const 
{
	GetValueRef( point, inidx, outidx ).LApplyTo( radiance );
}


template< class MatrixType > 
void Avals_MatrixStore< MatrixType >::ApplyValue ( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_HR_WEIGHT_TYPE* radiance ) const 
{
	*radiance *= (SKTRAN_HR_WEIGHT_TYPE)(GetValueRef( point, inidx, outidx ).p11());
}


template< class MatrixType >
std::unique_ptr< Aval_ScalarIteratorManager_Base > Avals_MatrixStore< MatrixType >::ScalarAvalIteratorManager ( const SKTRAN_HR_Diffuse_Point& point ) const
{
    std::unique_ptr< Aval_BScalarIteratorManager > manager( new Aval_BScalarIteratorManager );
  
    SKTRAN_ScatMat_MIMSNC aval;

    const size_t numincoming ( point.NumIncomingRays( ) );
    const size_t numoutgoing ( (point.IsGroundPoint() ? 1 : point.NumOutGoingRays( ))  );

    manager->m_vals.resize( numincoming*numoutgoing );
    std::vector< SKTRAN_HR_WEIGHT_TYPE >::iterator iterator ( manager->m_vals.begin() );
    typename std::vector< MatrixType >::const_iterator matrixIterator ( m_Avals.begin() + point.ScatterPropertyIdx( 0, 0 ) );

    // By convention, scalar Avals are accessed in "incomingRayIndex-major" order. 
    // If this changes, it needs to be changed in all the iterator manager clases 
    // and in the scalar radstore classes
    for( size_t inidx=0; inidx<numincoming; ++inidx ){
        for( size_t outidx=0; outidx<numoutgoing; ++outidx ){
            *iterator = (SKTRAN_HR_WEIGHT_TYPE) matrixIterator->p11();

            ++iterator;
            ++matrixIterator;
        }
    }

    manager->m_valIterator = manager->m_vals.begin();
    return std::unique_ptr< Aval_ScalarIteratorManager_Base > ( std::move(manager) );
}


template<>
void Avals_MatrixStore< SKTRAN_ScatMat_MIMSNC >::CreateEtaCalculator( std::unique_ptr< EtaCalculator_Base >& etaCalculator ) const
{
	etaCalculator= std::unique_ptr< EtaCalculator_Base >( new EtaCalculator_DoRotation ); 
}


template<>
void Avals_MatrixStore< SKTRAN_PhaseMat_MIMSNC >::CreateEtaCalculator( std::unique_ptr< EtaCalculator_Base >& etaCalculator ) const
{
	etaCalculator= std::unique_ptr< EtaCalculator_Base >( new EtaCalculator_NoRotation ); 
}


template<>
void Avals_MatrixStore< SKTRAN_ScatMat_MIMSNC >::CreateEtaCalculatorForStoragePhase( std::unique_ptr< EtaCalculator_Base >& etaCalculator ) const
{
	etaCalculator= std::unique_ptr< EtaCalculator_Base >( new EtaCalculator_NoRotation ); 
}


template<>
void Avals_MatrixStore< SKTRAN_PhaseMat_MIMSNC >::CreateEtaCalculatorForStoragePhase( std::unique_ptr< EtaCalculator_Base >& etaCalculator ) const
{
	etaCalculator= std::unique_ptr< EtaCalculator_Base >( new EtaCalculator_DoRotation ); 
}


// Explicitly instantiate classes so we can keep implementation in the .cpp 
template class Avals_MatrixStore< SKTRAN_ScatMat_MIMSNC >;
template class Avals_MatrixStore< SKTRAN_PhaseMat_MIMSNC >;


