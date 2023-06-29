#include "include/sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Source::SourceTermAtPoint		2013-06-27*/
/** Calculates the diffuse source term(radiance/m)  at a given point
 *  in the atmosphere and for a given look direction
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Source::SourceTermAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, double* source ) const
{
	bool ok = true;

	ok = ok && m_diffusetable->DiffuseSource( qobj, *source );

	return ok;
}

bool SKTRAN_HR_Diffuse_Source::SourceTermAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source ) const
{
	bool ok = true;

	ok = ok && m_diffusetable->DiffuseSource( qobj, *source );

	return ok;
}


bool SKTRAN_HR_Diffuse_Source::GroundSourceAtPoint(const SKTRAN_SourceTermQueryObject_Base& qobj, double*             source) const
{
    bool ok = true;

    ok = ok && m_diffusetable->GroundSource( qobj, *source );

    return ok;
}

bool SKTRAN_HR_Diffuse_Source::GroundSourceAtPoint(const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source) const
{
    bool ok = true;

    ok = ok && m_diffusetable->GroundSource( qobj, *source );

    return ok;
}
