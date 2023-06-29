#include "../sasktranv21_internals.h"
#include "../sasktran_legacy21.h"
#include "sktran_legacy_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableRayLOS_Legacy		 2014- 12- 9*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableRayLOS_Legacy::SKTRAN_TableRayLOS_Legacy   ( std::weak_ptr< const SKTRAN_RayFactory_Base> rayfactory)
	                     : SKTRANSO_TableRayLOS( rayfactory)
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableDiffusePoints_2D_Height_SZA::ConfigureGeometry		2010-5-28*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableRayLOS_Legacy::ConfigureGeometry( SKTRANSO_Quadrature_TLS_V21* quadrature,	 SKTRANSO_RayLOSGeometry_V21* losray )
{
	bool									ok = true;
	std::vector< SKTRAN_Distance>	locations;
	SKTRAN_Distance				startintercept;
	SKTRAN_Distance				exitintercept;
	SKTRAN_Distance				centerintercept;
	size_t									idx;
	size_t									idx1;
	size_t									numcells;

	numcells = losray->Storage()->NumCells();														
	locations.resize( numcells + 2 );

	startintercept =  losray->Storage()->DistanceOfPointFromOrigin( 0);
	locations[0] = startintercept;
	idx1 = 1;
	for (idx = 0; idx < numcells; idx++)
	{
		startintercept =  losray->Storage()->DistanceOfPointFromOrigin( idx);
		exitintercept  =  losray->Storage()->DistanceOfPointFromOrigin( idx+1);
		centerintercept  =  SKTRAN_DBL_TO_DISTANCE( 0.5*(startintercept + exitintercept) );					// get the center intercepts
		locations[idx1]  = centerintercept;
		++idx1;
	}
	locations[idx1] = exitintercept;
	ok = ok && ConfigureRayToSunLocations( locations, quadrature, losray );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TableRayLOS_Legacy::ConfigureGeometry, Error configure the line of sight ray-to-sun table. Thats a problem");
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_TableRayLOSFactory_Legacy::SKTRAN_TableRayLOSFactory_Legacy		2010-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableRayLOSFactory_Legacy::SKTRAN_TableRayLOSFactory_Legacy()
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableRayLOSFactory_Legacy::~SKTRAN_TableRayLOSFactory_Legacy		2010-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_TableRayLOSFactory_Legacy::~SKTRAN_TableRayLOSFactory_Legacy()
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_TableRayLOSFactory_Legacy::CreateInternalSingleScatterTable		2010-6-22*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_TableRayLOSFactory_Legacy::CreateInternalSingleScatterTable( SKTRANSO_TableRayLOS** internalsinglescattertable, std::weak_ptr< const SKTRAN_RayFactory_Base> rayfactory ) const
{
	bool	ok;

	*internalsinglescattertable = new SKTRAN_TableRayLOS_Legacy( rayfactory);

	ok = (*internalsinglescattertable != NULL);
	if (ok)
	{
		(*internalsinglescattertable)->AddRef();
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_TableRayLOSFactory_Legacy::CreateInternalSingleScatterTable, Error allocating memory for internal single scatter table object");
	}
	return ok;
}

