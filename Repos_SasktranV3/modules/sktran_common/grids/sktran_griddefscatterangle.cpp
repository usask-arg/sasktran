#include "../sktran_common.h"

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefScatterAngle_V21::Configure		2009-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GridDefScatterAngle_V21::Configure( double degreeresolution, double minssa, double maxssa)
{
	SKTRAN_GridIndex	numangles;
	bool				ok;
	SKTRAN_GridIndex	idx;
	m_maxcosssa = nxmath::cosd(minssa);
	m_mincosssa = nxmath::cosd(maxssa);

	numangles   = (SKTRAN_GridIndex)((maxssa - minssa)/degreeresolution + 0.5);
	double deltastep = (m_maxcosssa - m_mincosssa)/numangles;										// Cosine varies from -1.0 to +1.0 over range
	ok = AllocateGridArray( numangles+1);								// allocate array, add on the last element for the final +1
	if (ok)	
	{
		for (idx = 0; idx < numangles; idx ++)
		{
			AtVar(idx) = m_mincosssa + idx*deltastep;
		}
		AtVar(numangles) = m_maxcosssa;
        m_invdeltastep = 1.0 / deltastep;
		ok = ok && SetGridSearchMode( SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM );
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_GridDefScatterAngle_V21::Configure, error configuring the scattering angle table. Thats not good");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefScatterAngle_V21::CopyGridArray		2009-2-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GridDefScatterAngle_V21::CopyGridArray(  const double* source, size_t numpoints )
{
	nxLog::Record(NXLOG_WARNING, "SKTRAN_GridDefScatterAngle_V21::CopyGridArray, do not call CopyGridArray for the scattering angles, call Configure");
	NXASSERT((false));
	return false;
}




/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefScatterAngle_V21::FindBoundingIndices		2009-2-4*/
/** An overloaded version that makes use of the fact that the scatttering angle
 *	table is uniformly spaced.**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_GridDefScatterAngle_V21::FindBoundingIndices( double x, ENUM_INTERPOLATIONMODE outrange, SKTRAN_GridIndex* lowercell, double* lowerweight, SKTRAN_GridIndex* uppercell, double* upperweight) const
{
	SKTRAN_GridIndex	idx1;
	SKTRAN_GridIndex	idx2;

	NXASSERT(( (x >= m_mincosssa - 0.0001) && ( x < m_maxcosssa + 0.0001) ));
	idx1 = (SKTRAN_GridIndex)((x - m_mincosssa)*m_invdeltastep);
	NXASSERT((idx1 >= 0 ));
	idx2 = idx1 + 1;
	if (idx2 >= NumAngles())
	{
		--idx2;
		--idx1;
	}
	*lowercell	= idx1;
	*uppercell  = idx2;
	*upperweight   = (x - At(idx1))*m_invdeltastep;		// Get the linear interpolation weight to be applied to the upper index
	*lowerweight   = 1.0 - *upperweight;						// Get the linear interpolation weight to be applied to the lower index
	return true;													// return ok.
}

