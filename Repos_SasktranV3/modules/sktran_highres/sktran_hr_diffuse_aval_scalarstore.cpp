#include "include/sktran_hr_internals.h"


bool Avals_ScalarStore::AllocateStorage( size_t incomingcounter, size_t outgoingcounter, size_t scattvalcounter, size_t pointscounter )
{
	m_Avals   .resize( scattvalcounter     );
//	m_Acolind .resize( scattvalcounter     );
//	m_Arowptr .resize( outgoingcounter + 1 );

//	m_Arowptr[m_Arowptr.size()-1] = static_cast<int>(m_Avals.size());

	return true;
}


bool Avals_ScalarStore::ComputeMultipliersAndAdjustScatterArray( const SKTRAN_HR_Diffuse_Point& point )
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
			ratio  += m_Avals[ point.ScatterPropertyIdx(outrayidx, inrayidx) ]*w;
		}
		if ( scatextinction <= 0.0)					// It is occassionally possible to a get a point where there is no scatter for example
		{											// right at the top of the atmosphere
			if ( ratio != 0.0 )
			{
				nxLog::Record(NXLOG_WARNING, "Avals_ScalarStore::ComputeMultipliersAndAdjustScatterArray, Oh oh. we have a non-zero integral (Sum (A*w) = %e) even though the scatextinction value is zero or less", (double)ratio, (double)scatextinction);
				throw("Avals_ScalarStore::ComputeMultipliersAndAdjustScatterArray, inconsistent parameters");
			}
		}
		else
		{
			NXASSERT( (scatextinction > 0.0) );
			ratio /= incomingsphere->CubatureWeightAt( inrayidx );
			ratio /= scatextinction;
			NXASSERT( (NXFINITE(ratio)) );

			for( outidx = 0; outidx < numoutscats; outidx++ )
			{
				outrayidx = point.UniqueScatterOutgoing( outidx );
				m_Avals[ point.ScatterPropertyIdx(outrayidx, inrayidx) ] /= SKTRAN_HR_DBL_TO_WEIGHT(ratio);
			}
			if ( (ratio <0.177827941003892) || (ratio > 5.62341325190349) || ( !NXFINITE(ratio)) )  // if ( (ratio < pow(10,-0.75)) || (ratio > pow(10,0.75)) )
			{
				static int counter = 0;
				if (counter < 25)
				{
					nxLog::Record(NXLOG_WARNING, "   Large Phase Function Adjustment: %7.4e  %7.4e  %7.4e", ratio, incomingsphere->CubatureWeightAt( inrayidx ), scatextinction);
					counter++;
				}
			}
		}
	}
	
	return ok;
}


bool Avals_ScalarStore::CalcScatteringMatrixPoint( const SKTRAN_HR_Diffuse_Point& point )
{
	bool ok = true;

	if( ! point.IsGroundPoint() ){
		ok = ok && CalcScatteringMatrixPoint_Participating ( point );
	} else{
		ok = ok && CalcScatteringMatrixPoint_Boundary( point );
	}

	return ok;
}


bool Avals_ScalarStore::CalcScatteringMatrixPoint_Participating ( const SKTRAN_HR_Diffuse_Point& point )
{
	bool ok = true;

	size_t numinscats		= point.NumUniqueScatterIncoming();
	size_t numoutscats		= point.NumUniqueScatterOutgoing();
	size_t outrayidx;
	size_t inrayidx;

	nxVector				outdir;
	nxVector				indir;
	double					cosangle;
	double					scattcoeff;

	SKTRAN_TableOpticalProperties_PointCache cache;
	SKTRAN_ScatMat_MIMSNC temp;

	NXASSERT( (!point.IsGroundPoint()));

	m_opticaltable->CreateInterpolationForPoint(point.Location(), cache);

	for( size_t outidx = 0; outidx < numoutscats; outidx++ )
	{
		outrayidx = point.UniqueScatterOutgoing( outidx );
		point.OutgoingRayLocalCoords( outrayidx, outdir );
		for( size_t inidx = 0; inidx < numinscats; inidx++ )
		{
			inrayidx = point.UniqueScatterIncoming( inidx );
			indir = point.IncomingUnitRayLocalCoords( inrayidx );
			cosangle =	-1*(indir.X()*outdir.X() +
						indir.Y()*outdir.Y() +
						indir.Z()*outdir.Z());
			cache.Interpolate_kscattPerM(cosangle, temp);
			scattcoeff = (double) temp.At(1, 1);
			scattcoeff *= point.InCubatureWeight( inrayidx );
			m_Avals  [ point.ScatterPropertyIdx(outrayidx, inrayidx) ] = SKTRAN_HR_DBL_TO_WEIGHT(scattcoeff);
			
		}
	}
	ok = ok && ComputeMultipliersAndAdjustScatterArray( point );

	return ok;
}



/*-----------------------------------------------------------------------------
 *					Avals_ScalarStore::CalcScatteringMatrixPoint_Boundary		 2016- 11- 30*/
/** Sets up the A-matrix so it can calculate the downward flux of each incoming 
 *	radiance of a ground point. This stores the cos(Theta).DOmega of each incoming 
 *	radiance. THe function RadStore_Scalar::ScatterPoint_Scalar will use the
 *	A matrix values to generate the N incoming downward fluxes into the "outgoing" storage area
 *	of RadStore_Scalar.
 *
 *	Note: that we dont take into account albedo in the scattering matrix. This is done so 
 *	that we can use an albedo that varies based on position
 *
 **/
/*---------------------------------------------------------------------------*/

bool Avals_ScalarStore::CalcScatteringMatrixPoint_Boundary ( const SKTRAN_HR_Diffuse_Point& point )
{
	bool ok = true;

	size_t numinscats		= point.NumUniqueScatterIncoming();

	HELIODETIC_UNITVECTOR	ingrounddir;
	double					cosangle;
	double					scattcoeff;
	size_t                  inrayidx;

	NXASSERT( (point.IsGroundPoint()));													// make sure this is a ground point 
	ok =  point.IsGroundPoint();
	if (ok)
	{
		for( size_t inidx = 0; inidx < numinscats; inidx++ )								// For each incoming ray
		{																					// of the grounbd point
			inrayidx = point.UniqueScatterIncoming( inidx );								// Get the index of the ray from our table	
			ingrounddir = point.IncomingRayGlobalCoords( inrayidx );									// Get the incoming HELIODETIC_UNITVECTOR
			cosangle = fabs(ingrounddir & point.Location().LocalZenith() );					// and get the cosine of the zenith angle
			scattcoeff = cosangle;															// save the cos(theta) term
			scattcoeff *= point.InCubatureWeight( inrayidx );												// get the delta solid angle  (DOmega) of this incoming ray
			m_Avals  [ point.GroundDownwardFluxFactorsIdx( inrayidx ) ] = SKTRAN_HR_DBL_TO_WEIGHT(scattcoeff);	// and store the product in the A-Matrix 
		}
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"Avals_ScalarStore::CalcScatteringMatrixPoint_Boundary, the incoming point is not a ground point. Thats not good");
	}
	return ok;
}


void Avals_ScalarStore::GetValue ( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_ScatMat_MIMSNC* aval_matrix ) const
{


	if( ! point.IsGroundPoint() ){
		nxVector outdir;
		nxVector indir;
		double cosangle;

		point.OutgoingRayLocalCoords( outidx, outdir );
		indir = point.IncomingUnitRayLocalCoords( inidx );
		cosangle =	-1*(indir.X()*outdir.X() +
						indir.Y()*outdir.Y() +
						indir.Z()*outdir.Z());
		m_opticaltable->GetScatteringMatrixCM2( point.Location(), cosangle, *aval_matrix ); // Don't have to check ok because we know this point has been queried before
        double ratio = ( m_Avals[ point.ScatterPropertyIdx(outidx, inidx) ] / aval_matrix->p11());
        *aval_matrix *= ratio;
	} else{
        SKRTFLOAT scalarValue( m_Avals[ point.GroundDownwardFluxFactorsIdx( inidx ) ] ); // Get post-photon-conservation value
		aval_matrix->SetTo( 0.0 );
		aval_matrix->AssignAt( 1,1, scalarValue );
	}
    
}


void Avals_ScalarStore::ApplyValue ( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_Stokes_NC* radiance ) const 
{
	SKTRAN_ScatMat_MIMSNC s;
	GetValue( point, inidx, outidx, &s );
	s.LApplyTo( radiance );
}


void Avals_ScalarStore::ApplyValue ( const SKTRAN_HR_Diffuse_Point& point, size_t inidx, size_t outidx, SKTRAN_HR_WEIGHT_TYPE* radiance ) const 
{
	SKTRAN_ScatMat_MIMSNC s;
	GetValue( point, inidx, outidx, &s );
	*radiance *= (SKTRAN_HR_WEIGHT_TYPE)(s.p11());
}


std::unique_ptr< Aval_ScalarIteratorManager_Base >  Avals_ScalarStore::ScalarAvalIteratorManager ( const SKTRAN_HR_Diffuse_Point& point ) const
{
	const SKTRAN_HR_WEIGHT_TYPE*		parent;
    
	parent =  (point.IsGroundPoint() ) ? &m_Avals[point.GroundDownwardFluxFactorsIdx(0)]  : &m_Avals[point.ScatterPropertyIdx(0, 0)];
    return std::unique_ptr< Aval_ScalarIteratorManager_Base > ( new Aval_AScalarIteratorManager(parent) );
}


void Avals_ScalarStore::Aval_AScalarIteratorManager::Advance ( size_t numToAdvance )
{
    m_currentAval += numToAdvance;
}


SKTRAN_HR_WEIGHT_TYPE Avals_ScalarStore::Aval_AScalarIteratorManager::ApplyAndAdvance ( SKTRAN_HR_WEIGHT_TYPE radiance )
{
    return radiance * *(m_currentAval++);
}

