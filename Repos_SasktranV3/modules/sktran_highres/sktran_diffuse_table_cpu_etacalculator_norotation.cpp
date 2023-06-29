#include "include/sktran_hr_internals.h"


void EtaCalculator_NoRotation::UpdateOutgoingIndex( const SKTRAN_HR_Diffuse_Point& point, size_t outidx ) 
{
	// Do nothing, nothing has been stored 
}


void EtaCalculator_NoRotation::RefToScatt( SKTRAN_Stokes_NC& stokes ) const 
{
	return; // eta terms are built into the phase matrix
};


void EtaCalculator_NoRotation::ScattToRef( SKTRAN_Stokes_NC& stokes ) const 
{
	return; // eta terms are built into the phase matrix
};


void EtaCalculator_NoRotation::CalculateEtas( const SKTRAN_HR_Diffuse_Point& point, size_t inidx ) 
{
	return; // this class does *not* calculate etas
};


void EtaCalculator_NoRotation::ScattMatToPhaseMat  ( const SKTRAN_ScatMat_MIMSNC& s, SKTRAN_PhaseMat_MIMSNC& p ) const 
{
	p = s; // No rotation is performed
};


void EtaCalculator_NoRotation::SetPoint( const SKTRAN_HR_Diffuse_Point& point ) 
{
	// Do nothing, we'll just look up the phase matrices soon
}


