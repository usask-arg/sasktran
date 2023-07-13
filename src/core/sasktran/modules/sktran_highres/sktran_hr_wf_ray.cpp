#include "include/sktran_hr_internals.h"

bool SKTRAN_HR_WF_Ray::Allocate( size_t numperturb )
{
	m_weights.clear();
	m_addedopticaldepth.clear();

	m_weights.resize( numperturb );
	m_addedopticaldepth.resize( numperturb );

	return true;
}

bool SKTRAN_HR_WF_Ray_Polarized::Allocate( size_t numperturb )
{
    m_weights.clear();
    m_addedopticaldepth.clear();

    m_weights.resize( numperturb );
    m_addedopticaldepth.resize( numperturb );

    return true;
}