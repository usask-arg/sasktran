#include "../sktran_common.h"

SKTRAN_UnitSphere_Dummy::SKTRAN_UnitSphere_Dummy( const nxVector& loc )
{
	AllocateVertices( 1 );
	UnitVectorAtVar( 0 ) = loc.UnitVector();
}

bool SKTRAN_UnitSphere_Dummy::Triangulate( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices, size_t& speedHelper) const
{
	return Triangulate( unit, unit_indexptr, unit_weightptr, maxvertices );
}

bool SKTRAN_UnitSphere_Dummy::Triangulate( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices ) const
{
	unit_indexptr[0] = 0;
	unit_weightptr[0] = 1;
	return true;
}