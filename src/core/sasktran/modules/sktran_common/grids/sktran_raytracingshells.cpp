#include "../sktran_common.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefRayTracingShells_V21::SKTRAN_GridDefRayTracingShells_V21		2007-11-9*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_GridDefRayTracingShells_V21::SKTRAN_GridDefRayTracingShells_V21()
{
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefRayTracingShells_V21::~SKTRAN_GridDefRayTracingShells_V21		2007-11-9*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_GridDefRayTracingShells_V21::~SKTRAN_GridDefRayTracingShells_V21()
{
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefRayTracingShells_V21::NumCells						2007-11-10*/
/** Return the number of cells
**/
/*---------------------------------------------------------------------------*/

size_t SKTRAN_GridDefRayTracingShells_V21::NumCells	() const
{
	size_t	n;

	n = NumShells();
	if (n > 0) --n;
	return n;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefRayTracingShells_V21::LowestShell		2007-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_GridDefRayTracingShells_V21::LowestShell() const
{
	NXASSERT(( (NumShells() > 0) ));
	return front();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefRayTracingShells_V21::HighestShell		2007-11-10*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_GridDefRayTracingShells_V21::HighestShell() const
{
	NXASSERT(( (NumShells() > 0) ));
	return back();
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefRayTracingShells_V21::ConfigureGeometry		2007-11-9*/
/** Configure the spatial, ray-tracing grid. The latitude and longitude
 *	are the coordinates **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GridDefRayTracingShells_V21::ConfigureHeights( const double* shellAlts, size_t numshells )
{
	bool				ok;

	ok  = CopyGridArray( shellAlts, numshells );				// Copy the shell altitudes over to the grid
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefRayTracingShells_V21::ConfigureHeights		2014-1-31*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GridDefRayTracingShells_V21::ConfigureHeights( const std::vector<double>& shellAlts )
{
	bool				ok;

	ok  = CopyGridArray( shellAlts );				// Copy the shell altitudes over to the grid
	return ok;
}
