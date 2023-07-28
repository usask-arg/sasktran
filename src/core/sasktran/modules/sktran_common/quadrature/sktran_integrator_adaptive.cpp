#include "../sktran_common.h"
#include <cmath>


SKTRAN_OpticalPropertiesIntegrator_Adaptive::SKTRAN_OpticalPropertiesIntegrator_Adaptive()
{
	m_maxextinctiongradient = 0.95;
	m_maxopticaldepth = 100000;
	m_maxrayopticaldepthtosplit = 1000000;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Adaptive::CalculateRayTransmission_withMinContainer		 2014- 12- 2*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_OpticalPropertiesIntegrator_Adaptive::CalculateRayScalarTransmission_withMinContainer( SKTRAN_RayOptical_Base* baseray, double* transmission, bool totaltransmissiononly, bool usecachedtransmission  ) const
{
	SKTRAN_RayOptical_StraightQuadrature_Base*			straightray = dynamic_cast<SKTRAN_RayOptical_StraightQuadrature_Base*>(baseray);

	bool								ok = true;
	size_t								cellindex;
	double								opticaldepthcell;				// Optical Depth of the current cell size
	double								opticaldepthbuffer = 0;
	std::vector<double>*				odstorage  = straightray->OpticalDepthArrayVar();
	SKTRAN_RayStorage_Base*				raystorage = straightray->StorageVar();

	if( totaltransmissiononly )
	{
		odstorage->resize(1);
	}
	else
	{
		// start with the geometry-defined shells
		odstorage->reserve( raystorage->NumQuadraturePoints() * 3 );
		odstorage->resize(0);
		odstorage->push_back(0.0);
		raystorage->Reserve( odstorage->capacity() );
	}

    // Note: Cannot vectorize by simply making loop upper bound loop-invariant. NumQuadraturePoints may change between iterations. 
	for( size_t endquadpt = 1; endquadpt < raystorage->NumQuadraturePoints(); endquadpt++ )
	{
		// calculate the optical depth, if its too big we need to add a cell boundary
		cellindex        = endquadpt-1;
		opticaldepthcell = OpticalDepthOfCell_withMinCache( straightray, cellindex);

		// this will have put the relevant extinctions into the ray minimum container
		double kstart = raystorage->ExtinctionAtCellStart( cellindex );
		double kend   = raystorage->ExtinctionAtCellStart( endquadpt );

		while( opticaldepthcell > m_maxopticaldepth && std::min(kstart, kend) / std::max(kstart, kend) < m_maxextinctiongradient && opticaldepthbuffer < m_maxrayopticaldepthtosplit)
		{
			raystorage->SplitCell(cellindex);

			opticaldepthcell = OpticalDepthOfCell_withMinCache( straightray, cellindex );
			kend = raystorage->ExtinctionAtCellStart( cellindex+1 );
		}
		opticaldepthbuffer += opticaldepthcell;
		
		if(!totaltransmissiononly) odstorage->push_back( opticaldepthbuffer );
	}
	if( totaltransmissiononly )
	{
		odstorage->at(0) = opticaldepthbuffer;
	}

	return ok;
}


bool SKTRAN_OpticalPropertiesIntegrator_Adaptive::CalculateRayScalarTransmission( SKTRAN_RayOptical_Base* baseray, double* transmission, bool totaltransmissiononly, bool usecachedtransmission  ) const
{
	SKTRAN_RayOptical_StraightQuadrature_Base*			ray = dynamic_cast<SKTRAN_RayOptical_StraightQuadrature_Base*>(baseray);
	bool								ok = true;
	double								opticaldepthcell;				// Optical Depth of the current cell size
	double								kstart;							// extinction at the start of the cell
	double								kend;							// extinction at the end of the cell
	HELIODETIC_POINT					cache;
	HELIODETIC_POINT					cacheSP;
	HELIODETIC_POINT					cacheEP;
	double								opticaldepthbuffer = 0;
	std::vector<double>*				odstorage  = ray->OpticalDepthArrayVar();
	SKTRAN_RayStorage_Base*				raystorage = ray->StorageVar();
//	SKTRAN_RayStorage_Base*				raystorage = ray->StorageVar();//GRAY()->StraightStorageVar();

	if (raystorage->NumCells() > 0)						// initialize buffer points
	{
		raystorage->LocationOfPoint(0,&cacheSP);
		cacheEP = cacheSP;
	}

	if( totaltransmissiononly )
	{
		odstorage->resize(1);
	}
	else
	{
		odstorage->reserve( raystorage->NumQuadraturePoints() * 3 );
		odstorage->resize(0);							// start with the geometry-defined shells
		odstorage->push_back(0.0);
	}
	raystorage->Reserve( raystorage->NumQuadraturePoints() * 3 );

    // Note: Cannot vectorize by simply making loop upper bound loop-invariant. NumQuadraturePoints may change between iterations. 
	for( size_t endquadpt = 1; endquadpt < raystorage->NumQuadraturePoints(); endquadpt++ )
	{
		cache = cacheEP;
		opticaldepthcell = OpticalDepthOfCell_advanceCachePoints( ray, endquadpt-1, cacheSP, cacheEP, &kstart, &kend );

		size_t cellindex = endquadpt - 1;

		
		// we want to add a new point if the gradient of the extinction is large, and the optical depth
		// has an appreciable effect on radiance
		while( std::min(kstart, kend) / std::max(kstart, kend) < m_maxextinctiongradient && opticaldepthcell > m_maxopticaldepth && opticaldepthbuffer < m_maxrayopticaldepthtosplit )
		{
			raystorage->SplitCell(cellindex);
			cacheEP = cache;	// Need to start at previous segment's start point
			opticaldepthcell = OpticalDepthOfCell_advanceCachePoints( ray, endquadpt-1, cacheSP, cacheEP, &kstart, &kend );
		}
		opticaldepthbuffer += opticaldepthcell;
		if( !totaltransmissiononly ) odstorage->push_back( opticaldepthbuffer );
	}
	if( totaltransmissiononly )
	{
		odstorage->at(0) = opticaldepthbuffer;
	}
	return ok;
}

