#include "../sktran_common.h"


/*-----------------------------------------------------------------------------
 *					100 point Minimum Potential enrgy point distribution	2008-2-10*/
/** Rob Womersley, http://web.maths.unsw.edu.au/~rsw/Sphere/
 *	These are the vertices of 100 points on the unit sphere using a minimum energy configuration
 *	For each point set the text file has four items per row: the xj, yj, and zj cartesian coordinates in [-1, 1],
 *	and the cubature weight wj for that point. The number of rows is equal to the number of points.
 *	All points are on the unit sphere so should have xj2 + yj2 + zj2 = 1.
 *	The file names have three components: point set (me), degree(09), number of points(100).
 **/
/*---------------------------------------------------------------------------*/

#include "me06_0049.hpp"
#include "me07_0064.hpp"
#include "me08_0081.hpp"
#include "me09_0100.hpp"
#include "me12_0169.hpp"
#include "me13_0196.hpp"
#include "me14_0225.hpp"
#include "me15_0256.hpp"
#include "me17_0324.hpp"
#include "me19_0400.hpp"
#include "me22_0529.hpp"
#include "me24_0625.hpp"
#include "me26_0729.hpp"
#include "me28_0841.hpp"
#include "me29_0900.hpp"
#include "me31_1024.hpp"
#include "me35_1296.hpp"
#include "me39_1600.hpp"
#include <climits>

/*-----------------------------------------------------------------------------
 *					ME_TABLE_ENTRY		2008-8-12*/
/** a quick structure for storing the Minimum energy sphere tables**/
/*---------------------------------------------------------------------------*/

struct ME_TABLE_ENTRY
{
	double*		array;
	size_t		numelements;
};


/*-----------------------------------------------------------------------------
 *					g_tablentries		2008-8-12*/
/** The array of minimum energy spheres supported by this class**/
/*---------------------------------------------------------------------------*/

static ME_TABLE_ENTRY g_tablentries[] = 
{
	{g_me06_0049,  N_ELEMENTS(g_me06_0049)/4},
	{g_me07_0064,  N_ELEMENTS(g_me07_0064)/4},
	{g_me08_0081,  N_ELEMENTS(g_me08_0081)/4},
	{g_me09_0100,  N_ELEMENTS(g_me09_0100)/4},
	{g_me12_0169,  N_ELEMENTS(g_me12_0169)/4},
	{g_me13_0196,  N_ELEMENTS(g_me13_0196)/4},
	{g_me14_0225,  N_ELEMENTS(g_me14_0225)/4},
	{g_me15_0256,  N_ELEMENTS(g_me15_0256)/4},
	{g_me17_0324,  N_ELEMENTS(g_me17_0324)/4},
	{g_me19_0400,  N_ELEMENTS(g_me19_0400)/4},
	{g_me22_0529,  N_ELEMENTS(g_me22_0529)/4},
	{g_me24_0625,  N_ELEMENTS(g_me24_0625)/4},
	{g_me26_0729,  N_ELEMENTS(g_me26_0729)/4},
	{g_me28_0841,  N_ELEMENTS(g_me28_0841)/4},
	{g_me29_0900,  N_ELEMENTS(g_me29_0900)/4},
	{g_me31_1024,  N_ELEMENTS(g_me31_1024)/4},
	{g_me35_1296,  N_ELEMENTS(g_me35_1296)/4},
	{g_me39_1600,  N_ELEMENTS(g_me39_1600)/4}
};


/*-----------------------------------------------------------------------------
 *					EntryClosestToRequestedNumber		2008-8-12*/
/** Find the minimum energy sphere closest to the requested number of vertices
**/
/*---------------------------------------------------------------------------*/

static const ME_TABLE_ENTRY* EntryClosestToRequestedNumber( size_t numvertexrequested )
{
	ME_TABLE_ENTRY*	bestentry    = NULL;
	int			    bestdistance = INT_MAX;
	int				distance;
	size_t			idx;

	for (idx = 0; idx < N_ELEMENTS(g_tablentries); idx++)
	{
		distance = abs( (int)numvertexrequested - (int)(g_tablentries[idx].numelements) );
		if (distance < bestdistance )
		{
			bestentry = &g_tablentries[idx];
			bestdistance = distance;
		}
	}
	NXASSERT(( bestentry != NULL ));
	return  bestentry;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefOutboundUnitSphereME100::SKTRAN_GridDefOutboundUnitSphereME100		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_UnitSphereME::SKTRAN_UnitSphereME( size_t numvertex )
{
	bool					ok;
	size_t					inidx;
	size_t					idx;
	const ME_TABLE_ENTRY*	entry;
//	size_t					numvertex = 100;

	entry = EntryClosestToRequestedNumber(numvertex);
	NXASSERT(( entry != NULL ));
	if (entry->numelements != numvertex)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_UnitSphereME: Configure, Using %u outgoing vertices as closest match to requested %u vertices", (size_t)entry->numelements, (size_t)numvertex );
	}

	ok     = AllocateVertices( entry->numelements );
	if (ok)
	{

		for (idx = 0; idx < entry->numelements; idx++ )
		{
			inidx = 4*idx;
			UnitVectorAtVar(idx).SetCoords( entry->array[inidx+0], entry->array[inidx+1], entry->array[inidx+2] );
		    CubatureWeightAtVar(idx) = entry->array[inidx+3];			
		}
		ok = ok && InitializeLookupTable();
	}
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereME::~SKTRAN_UnitSphereME		2008-2-10*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_UnitSphereME::~SKTRAN_UnitSphereME()
{
}

