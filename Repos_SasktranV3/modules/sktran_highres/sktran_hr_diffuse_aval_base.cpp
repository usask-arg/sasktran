#include "include/sktran_hr_internals.h"


Avals_Base::Avals_Base ( )
{
	m_opticaltable = nullptr;
}


Avals_Base::~Avals_Base ( )
{
	ReleaseResources ( );
}


void Avals_Base::ReleaseResources ( )
{
	if(nullptr!=m_opticaltable) m_opticaltable->Release ( );
	m_opticaltable = nullptr;
}


void Avals_Base::SetOpticalTable ( const SKTRAN_TableOpticalProperties_Base* table )
{
	if(nullptr!=table) table->AddRef( );
	if(nullptr!=m_opticaltable) m_opticaltable->Release( );
	m_opticaltable = table;
}


