#include "../sktran_common.h"


SKTRAN_Sun_Base::SKTRAN_Sun_Base( )
{

}

SKTRAN_Sun_Base::~SKTRAN_Sun_Base( )
{

}

/***********************
 * SKTRAN_Sun_Point
 ***********************/

SKTRAN_Sun_Point::SKTRAN_Sun_Point( )
{
	m_sunUnit.SetCoords(0.0, 0.0, 1.0);
}

SKTRAN_Sun_Point::~SKTRAN_Sun_Point( )
{

}

//bool SKTRAN_Sun_Point::MakeThreadsafeFor( size_t )
//{
//	return true;
//}

void SKTRAN_Sun_Point::UpdateSun( ) const
{
	
}

void SKTRAN_Sun_Point::SunUnitVector( HELIODETIC_UNITVECTOR* h ) const
{
	h->SetCoords(0.0, 0.0, 1.0);
}

//void SKTRAN_Sun_Point::SunUnitVector( HELIODETIC_VECTOR& h ) const
//{
//	h.SetCoords(0.0, 0.0, 1.0);
//}

const HELIODETIC_UNITVECTOR& SKTRAN_Sun_Point::GetSunUnit( ) const
{
	return m_sunUnit;
}

double SKTRAN_Sun_Point::CosAngleToSun( const HELIODETIC_UNITVECTOR&  h) const
{
	return h.Z();
}

//double SKTRAN_Sun_Point::CosAngleToSun( const HELIODETIC_VECTOR& h) const
//{
//	return h.UnitVector().Z();
//}

double SKTRAN_Sun_Point::ComponentToSun( const HELIODETIC_VECTOR& h) const
{
	return h.Z();
}
