#include "../sasktranv21_internals.h"

/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayStorageBaseGeometry::SKTRANSO_RayStorageBaseGeometry		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_RayStorageBaseGeometry::SKTRANSO_RayStorageBaseGeometry(  std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords )
	                            :SKTRAN_RayStorage_Straight( coords)
{
}

/*-----------------------------------------------------------------------------
 *					SKTRANSO_RayStorageBaseGeometry::~SKTRANSO_RayStorageBaseGeometry		2007-11-16*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRANSO_RayStorageBaseGeometry::~SKTRANSO_RayStorageBaseGeometry()
{
}

