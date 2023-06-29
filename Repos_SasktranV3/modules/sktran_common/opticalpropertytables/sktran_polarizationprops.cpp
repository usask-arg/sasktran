#include "../sktran_common.h"

/**************
 * NoPolarized class
 *************/

SKTRAN_PolarizationProperties_NoPolarization::SKTRAN_PolarizationProperties_NoPolarization ()
{
}

SKTRAN_PolarizationProperties_NoPolarization::~SKTRAN_PolarizationProperties_NoPolarization ()
{
}


bool SKTRAN_PolarizationProperties_NoPolarization::Allocate	( size_t numcachepoints )
{
	bool	ok = true;
    m_phaseFunction.resize( numcachepoints );
	return ok;
}

bool SKTRAN_PolarizationProperties_NoPolarization::StorePolarizationPropsCM2( size_t index, const skRTPhaseMatrix& pmatrix, SKTRAN_AtmosphericOpticalState_V21& state )
{
	bool ok = true;

    NXASSERT( index < m_phaseFunction.size() );
    m_phaseFunction[index] = pmatrix.At(1,1);
    
	return ok;
}


bool SKTRAN_PolarizationProperties_NoPolarization::GetPhaseFunctionCM2 ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, double& phaseFunction ) const
{
    bool ok = true;

    phaseFunction = 0;

    for(size_t midx=0; midx<numels; ++midx){
        phaseFunction += m_phaseFunction[gridindex[midx]]*gridweight[midx];
    }

    return ok;
}

bool SKTRAN_PolarizationProperties_NoPolarization::GetScatteringMatrixCM2 ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_ScatMat_MIMSNC& pmatrix ) const
{
    bool ok = true;
    double p11; 

    pmatrix.SetTo( 0.0 );
    ok = ok && GetPhaseFunctionCM2( gridindex, gridweight, numels, p11 );
    pmatrix.AssignAt( 1,1, p11 );

    return ok;
}


bool SKTRAN_PolarizationProperties_NoPolarization::GetResultOfUnpolarizedScatterCM2( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_Stokes_NC& stokesvec ) const
{
	bool ok = true;
    double p11;

    ok = ok && GetPhaseFunctionCM2 ( gridindex, gridweight, numels, p11 );
	stokesvec.SetTo(0.0);
	stokesvec.Assign_I( p11 );

	return ok;
}


/**************
 * Polarized class
 *************/

SKTRAN_PolarizationProperties_Polarized::SKTRAN_PolarizationProperties_Polarized			()
{
}

SKTRAN_PolarizationProperties_Polarized::~SKTRAN_PolarizationProperties_Polarized			()
{
}

bool SKTRAN_PolarizationProperties_Polarized::Allocate ( size_t numcachepoints )
{
	bool	ok = true;

	m_muellerMatrices.resize( numcachepoints );
	
	return ok;
}


bool SKTRAN_PolarizationProperties_Polarized::StorePolarizationPropsCM2( size_t index, const skRTPhaseMatrix& pmatrix, SKTRAN_AtmosphericOpticalState_V21& state )
{
	bool ok = true;

	NXASSERT( index < m_muellerMatrices.size() );
	m_muellerMatrices[index] = ScatMatType(pmatrix);

	return ok;
}

bool SKTRAN_PolarizationProperties_Polarized::GetPhaseFunctionCM2( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, double& phaseFunction ) const
{
    bool ok = true;

    phaseFunction = 0.0;
    for(size_t midx=0; midx<numels; ++midx){
        phaseFunction += m_muellerMatrices[gridindex[midx]].p11() * gridweight[midx];
    }

    return ok;
}

bool SKTRAN_PolarizationProperties_Polarized::GetScatteringMatrixCM2( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_ScatMat_MIMSNC& matrix ) const
{
	bool ok = true;

	matrix.SetTo(0.0);

	for(size_t midx=0; midx<numels; ++midx){
		matrix.AddToThis( m_muellerMatrices[gridindex[midx]], gridweight[midx] );
	}
	
	return ok;
}

bool SKTRAN_PolarizationProperties_Polarized::GetResultOfUnpolarizedScatterCM2 ( const SKTRAN_GridIndex gridindex[], const double gridweight[], size_t numels, SKTRAN_Stokes_NC& stokesvec ) const
{
	bool ok = true;

	stokesvec.SetTo(0.0);
	for(size_t midx=0; midx<numels; ++midx){
        stokesvec.Add_I( m_muellerMatrices[gridindex[midx]].At(1,1) * gridweight[midx] );
        stokesvec.Add_Q( m_muellerMatrices[gridindex[midx]].At(2,1) * gridweight[midx] );
		stokesvec.Add_U( m_muellerMatrices[gridindex[midx]].At(3,1) * gridweight[midx] );
		stokesvec.Add_V( m_muellerMatrices[gridindex[midx]].At(4,1) * gridweight[midx] );
	}
	
	return ok;
}




/**************
 * Polarized with Eddington class
 *************/

SKTRAN_PolarizationProperties_Polarized_Eddington::SKTRAN_PolarizationProperties_Polarized_Eddington ()
{
}

SKTRAN_PolarizationProperties_Polarized_Eddington::~SKTRAN_PolarizationProperties_Polarized_Eddington ()
{
}

bool SKTRAN_PolarizationProperties_Polarized_Eddington::Allocate ( size_t numcachepoints )
{
	bool	ok = true;

	m_eddingtonExtinction.resize( numcachepoints );
	ok = ok && SKTRAN_PolarizationProperties_Polarized::Allocate ( numcachepoints );
	
	return ok;
}


bool SKTRAN_PolarizationProperties_Polarized_Eddington::StorePolarizationPropsCM2( size_t index, const skRTPhaseMatrix& pmatrix, SKTRAN_AtmosphericOpticalState_V21& state )
{
	bool ok = true;

	m_eddingtonExtinction[index] = state.DeltaForwardPercm();
	ok = ok && SKTRAN_PolarizationProperties_Polarized::StorePolarizationPropsCM2( index, pmatrix, state );

	return ok;
}
