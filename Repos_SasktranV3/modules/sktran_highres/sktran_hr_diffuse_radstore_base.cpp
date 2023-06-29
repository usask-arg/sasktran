#include "include/sktran_hr_internals.h"
#include <numeric>
#include <functional>


RadStore_Base::RadStore_Base ( )
{
	m_opticaltable = nullptr;

	m_order_incoming = 0;
	m_order_outgoing = 0;

	m_numincomingrays = 0;
	m_numoutgoingrays = 0;
}


RadStore_Base::~RadStore_Base ( )
{
	ReleaseResources( );
}


bool RadStore_Base::AllocateStorage ( size_t incomingcounter, size_t outgoingcounter )
{
	m_numincomingrays = incomingcounter;
	m_numoutgoingrays = outgoingcounter;
	return true;
}


bool RadStore_Base::CleanDiffuseIndexes ( )
{
	m_order_incoming = 0;
	m_order_outgoing = 0;
	return true;
}


void RadStore_Base::ReleaseResources ( )
{
	if(nullptr!=m_opticaltable) m_opticaltable->Release( );
	m_opticaltable = nullptr;
}


bool RadStore_Base::Initialize( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>&        coords,
																const SKTRAN_OpticalPropertiesIntegrator_Base& optintegrator,
																const SKTRAN_SourceTermIntegrator_Base&        srcintegrator,
																std::shared_ptr<const SKTRAN_RayFactory_Base>   rayfactory,
																const std::vector<SKTRAN_Source_Term*>&        sources,
																const SKTRAN_TableOpticalProperties_Base&      opticaltable )
{
	bool ok = true;
	
	opticaltable.AddRef();
	if(nullptr!=m_opticaltable) m_opticaltable->Release( );
	m_opticaltable = &opticaltable;

	return ok;
}


