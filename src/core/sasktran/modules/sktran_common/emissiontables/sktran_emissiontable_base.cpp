#include "../sktran_common.h"


SKTRAN_EmissionTable_Base::~SKTRAN_EmissionTable_Base( )
{

}

bool SKTRAN_EmissionTable_Base::SourceTermAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source   ) const
{
    bool ok = true;
    double scalar;
    ok = ok && SourceTermAtPoint( qobj, &scalar );
    source->SetTo( 0.0 );
    source->Assign_I( scalar );
    return ok;
}


bool SKTRAN_EmissionTable_Base::GroundSourceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC* source   ) const
{
    bool ok = true;
    double scalar;
    ok = ok && GroundSourceAtPoint( qobj, &scalar );
    source->SetTo( 0.0 );
    source->Assign_I( scalar );
    return ok;
}


bool SKTRAN_EmissionTable_Base::MonteCarlo_SingleScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance ) const
{
    bool ok = true;
    double scalar;
    ok = ok && MonteCarlo_SingleScatteredRadianceAtPoint( qobj, scalar );
    radiance.SetTo( 0.0 );
    radiance.Assign_I( scalar );
    return ok;
}


bool SKTRAN_EmissionTable_Base::MonteCarlo_GroundScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC& radiance ) const
{
    bool ok = true;
    double scalar;
    ok = ok && MonteCarlo_GroundScatteredRadianceAtPoint( qobj, scalar );
    radiance.SetTo( 0.0 );
    radiance.Assign_I( scalar );
    return ok;
}


