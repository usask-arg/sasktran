#include "include/sktran_hr_internals.h"

SKTRAN_HR_OutgoingSphereObject_Base::SKTRAN_HR_OutgoingSphereObject_Base()
{
	m_outgoingsphere = nullptr;
}
SKTRAN_HR_OutgoingSphereObject_Base::~SKTRAN_HR_OutgoingSphereObject_Base()
{
	ReleaseResources ( );
}

void SKTRAN_HR_OutgoingSphereObject_Base::ReleaseResources ( )
{
	if(nullptr!=m_outgoingsphere) m_outgoingsphere->Release();
	m_outgoingsphere = nullptr;
}

void SKTRAN_HR_OutgoingSphereObject_Base::SetOutgoingSphere ( const SKTRAN_UnitSphere_V2* sphere )
{
	if( nullptr!=sphere ) sphere->AddRef();
	if(nullptr!=m_outgoingsphere) m_outgoingsphere->Release();
	m_outgoingsphere = sphere;
}

bool SKTRAN_HR_OutgoingSphereObject_Base::TriangulateOnOutgoing( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices ) const
{
	return m_outgoingsphere->Triangulate( unit, unit_indexptr, unit_weightptr, maxvertices );
}

size_t SKTRAN_HR_OutgoingSphereObject_Base::NumOutgoingRays() const{
	return m_outgoingsphere->NumUnitVectors();
}

void SKTRAN_HR_OutgoingSphereObject_Base::OutgoingRayLocalCoords( size_t idx, nxVector& outray ) const
{
	outray = m_outgoingsphere->UnitVectorAt( idx );
}

double SKTRAN_HR_OutgoingSphereObject_Base::OutgoingCubatureWeight( size_t outidx ) const
{
	return m_outgoingsphere->CubatureWeightAt( outidx );
}
