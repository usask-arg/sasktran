#include "../sktran_common.h"

#include <vector>
#include <algorithm>


/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfCell		2013-06-11*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfCell( const SKTRAN_RayOptical_Base* baseray, size_t cellidx ) const
{
	const SKTRAN_RayOptical_StraightQuadrature_Base*	ray			= reinterpret_cast< const SKTRAN_RayOptical_StraightQuadrature_Base*>(baseray);
	double								opticaldepth;
	SKTRAN_Distance						startintercept;
	SKTRAN_Distance						endintercept;
	bool								ok = true;

	NXASSERT( cellidx < (ray->GetNumQuadraturePoints() - 1) );
	
	startintercept = ray->Storage()->DistanceOfPointFromOrigin( cellidx   );
	endintercept   = ray->Storage()->DistanceOfPointFromOrigin( cellidx+1 );
	
	opticaldepth = OpticalDepthOfSegment( cellidx, startintercept, endintercept, ray );
	return opticaldepth;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfCell_advanceCachePoints		2014-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfCell_advanceCachePoints( const SKTRAN_RayOptical_Base* baseray, size_t cellidx, HELIODETIC_POINT& cacheSP, HELIODETIC_POINT& cacheEP, double* kstart, double* kend ) const
{
	const SKTRAN_RayOptical_StraightQuadrature_Base*	ray			= reinterpret_cast< const SKTRAN_RayOptical_StraightQuadrature_Base*>(baseray);
	double			opticaldepth;
	SKTRAN_Distance startintercept;
	SKTRAN_Distance endintercept;
	bool			ok				= true;

	NXASSERT( cellidx < (ray->GetNumQuadraturePoints() - 1) );
	
	startintercept = ray->Storage()->DistanceOfPointFromOrigin( cellidx  );
	endintercept   = ray->Storage()->DistanceOfPointFromOrigin( cellidx+1);
	opticaldepth = OpticalDepthOfSegment_advanceCachePoints( cellidx, startintercept, endintercept, ray, cacheSP, cacheEP, kstart, kend );
	return opticaldepth;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfCell_withMinCache		2014-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfCell_withMinCache( const SKTRAN_RayOptical_Base* ray, size_t cellidx ) const
{
	double			opticaldepth;
	bool			ok				= true;

	NXASSERT( cellidx < (ray->GetNumQuadraturePoints() - 1) );
	
	opticaldepth = OpticalDepthOfSegment_withMinCache( cellidx, ray );
	NXASSERT(( (opticaldepth >= 0.0) && ( opticaldepth < 1.0E6) ));
	return opticaldepth;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfSegment		2013-06-12*/
/** Assume the input cacheEP is the new start point, input cacheSP is the previous (out of date) start point **/
/*---------------------------------------------------------------------------*/

double SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfSegment_withMinCache( size_t cellidx, const SKTRAN_RayOptical_Base* ray ) const
{
	double opticaldepth = -9999.0;

//	const SKTRAN_RayOptical_Straight*	ray			= reinterpret_cast< const SKTRAN_RayOptical_Straight*>(baseray);

	bool						ok = true;
	double						r0;
	double						r1;
	double						rt;
	double						t0;
	double						t1;
	double						sigma0 = -9999.0;
	double						sigma1 = -9999.0;
	//SKTRAN_OpticalDepthCalculator_LinearWithHeight	odcalculator;
	
//	r0 = cacheSP.Radius();

//	t0 = ray->Storage()->DistanceOfPointFromCellTangentPoint( cellidx, cellidx);
//	ok = ok && ray->GetQuadratureInterpParams_startPoint( cellidx+1, r1, t1, rt, cacheEP);
//	ok = ok && m_opticalprops->GetEffectiveExtinctionPerCMWithHeight( cacheSP, cacheEP, &sigmak, &sigmaf );

	// ok =       ray->Storage()->LocationOfPoint(  cellidx,   &location0);
	// ok = ok && ray->Storage()->LocationOfPoint(  cellidx+1, &location1);
	//ok = ok && m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1( ray->Storage(), cellidx, &sigma0, &sigma1 );

	//ok = ok && m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1(ray->GetWavelength(), ray->Storage(), cellidx, &sigma0, &sigma1);
	ok = ok && GetEffectiveExtinctionPerCMWithHeight1(ray, cellidx, &sigma0, &sigma1);
	
	NXASSERT( ( (sigma0 >= 0) && (sigma0 < 1.0E4)));
	NXASSERT( ( (sigma1 >= 0) && (sigma1 < 1.0E4)));

	//r0 = ray->GeometryRay()->Coordinates()->AltitudeToRadius( ray->StorageAccess().AltitudeOfQuadraturePoint( cellidx   ) );
	//r1 = ray->GeometryRay()->Coordinates()->AltitudeToRadius( ray->StorageAccess().AltitudeOfQuadraturePoint( cellidx+1 ) );
	r0 = ray->Storage()->RadiusOfPoint( cellidx   );
	r1 = ray->Storage()->RadiusOfPoint( cellidx+1 );
	t0 = ray->Storage()->DistanceOfPointFromCellTangentPoint( cellidx,   cellidx );
	t1 = ray->Storage()->DistanceOfPointFromCellTangentPoint( cellidx+1, cellidx );
	rt = ray->Storage()->RadiusOfCellTangentPoint(cellidx);

	//ok = ok && odcalculator.ConfigureQuadratureCoefficients( r0, r1, t0, t1, rt );
	if (ok)
	{
		//opticaldepth = odcalculator.OpticalDepthFromStartToEnd( sigma0, sigma1);
		ok = ok && GetOpticalDepthFromParams( r0, r1, t0, t1, rt, sigma0, sigma1, opticaldepth );
	}
	ok = ok && (opticaldepth >= 0.0);
	if (!ok)
	{
		if(1e-7 < fabs(opticaldepth) || fabs(r0-r1)>1.0)
		{
			// Only print warning if the error is really of a magnitude that wew care about
			nxLog::Record(NXLOG_WARNING,"SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfSegment_withMinCache, Error looking up optical depth of a segment, cellidx =%d, opticaldepth =%18.8e, sigma0= %18.8e, sigma1=%18.8e,  r0=%18.8e, r1=%18.8e, t0=%18.8e, t1=%18.8e", (int)cellidx, (double)opticaldepth, (double)sigma0, (double)sigma1, (double)r0, (double)r1, (double)t0, (double)t1);
			//printf("%1.20e\n%1.20e\n%1.20e\n%1.20e\n%1.20e\n%1.20e\n%1.20e\n%u\n\n", r0, r1, t0, t1, rt, sigma0, sigma1, cellidx);
			//printf("%1.20e\n%1.20e\n%1.20e\n%1.20e\n%1.20e\n%1.20e\n\n",ray->GRAY()->GetObserver().X(),ray->GRAY()->GetObserver().Y(),ray->GRAY()->GetObserver().Z(),ray->GRAY()->LookVector().X(),ray->GRAY()->LookVector().Y(),ray->GRAY()->LookVector().Z());
		}
		opticaldepth = 0;
	}
	return opticaldepth * ray->Storage()->CellCurvature(cellidx);

}

double SKTRAN_OpticalPropertiesIntegrator_Straight::TotalExtinctionPerCM(const SKTRAN_RayOptical_Base * ray, const HELIODETIC_POINT & point) const
{
	return m_opticalprops->TotalExtinctionPerCM(point);
}

bool SKTRAN_OpticalPropertiesIntegrator_Straight::GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base * ray, const HELIODETIC_POINT & startpoint, HELIODETIC_POINT & endpoint, double * sigma0, double * sigma1) const
{
	return m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1(startpoint, endpoint, sigma0, sigma1);
}

bool SKTRAN_OpticalPropertiesIntegrator_Straight::GetEffectiveExtinctionPerCMWithHeight1(const SKTRAN_RayOptical_Base * ray, size_t startPtIndex, double * sigma0, double * sigma1) const
{
	return m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1(ray->Storage(), startPtIndex, sigma0, sigma1);
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfSegment_advanceCachePoints		2014-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfSegment_advanceCachePoints( size_t cellidx, const SKTRAN_Distance& startintercept, const SKTRAN_Distance& endintercept, const SKTRAN_RayOptical_Base* baseray, HELIODETIC_POINT& cacheSP, HELIODETIC_POINT& cacheEP, double* kstart, double* kend ) const
{
	double opticaldepth;
	const SKTRAN_RayOptical_StraightQuadrature_Base*	ray			= reinterpret_cast< const SKTRAN_RayOptical_StraightQuadrature_Base*>(baseray);

	bool						ok = true;
	double						r0;
	double						r1;
	double						rt;
	double						t0;
	double						t1;
	double						sigma0;
	double						sigma1;
	//SKTRAN_OpticalDepthCalculator_LinearWithHeight	odcalculator;
	
	r0 = ray->Storage()->RadiusOfPoint(cellidx);
	t0 = ray->Storage()->DistanceOfPointFromCellTangentPoint( cellidx, cellidx);
	ok = ok && ray->GetQuadratureInterpParams_startPoint(   cellidx+1, cellidx,  &r1, &t1, &rt, &cacheEP);
	//ok = ok && m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1( cacheSP, cacheEP, &sigma0, &sigma1 );
	//ok = ok && GetEffectiveExtinctionPerCMWithHeight1(baseray, cacheSP, cacheEP, &sigma0, &sigma1);
	ok = ok && GetEffectiveExtinctionPerCMWithHeight1(baseray, cellidx, &sigma0, &sigma1);
	//ok = ok && odcalculator.ConfigureQuadratureCoefficients( r0,r1,t0,t1,rt);
	if (ok)
	{
		//opticaldepth = odcalculator.OpticalDepthFromStartToEnd( sigma0, sigma1);
		ok = ok && GetOpticalDepthFromParams(r0, r1, t0, t1, rt, sigma0, sigma1, opticaldepth);
	}
	ok = ok && (opticaldepth >= 0.0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfSegment_advanceCachePoints, Error looking up optical depth of a segment, cellidx = %d, opticaldepth =%18.8e, sigma0= %18.8e, sigma1=%18.8e,  r0=%18.8e, r1=%18.8e, t0=%18.8e, t1=%18.8e", (int)cellidx, (double)opticaldepth, (double)sigma0, (double)sigma1, (double)r0, (double)r1, (double)t0, (double)t1);
		opticaldepth = 0;
	}
	if( nullptr != kstart )
	{
		*kstart = sigma0;
	}
	if( nullptr != kend )
	{
		*kend = sigma1;
	}
	return opticaldepth * ray->Storage()->CellCurvature(cellidx);

}


/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfSegment		2014-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

double SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfSegment( size_t cellidx, const SKTRAN_Distance& startintercept, const SKTRAN_Distance& endintercept, const SKTRAN_RayOptical_Base* baseray ) const
{
	double opticaldepth = -99999.0;
	const SKTRAN_RayOptical_StraightQuadrature_Base*	ray			= reinterpret_cast< const SKTRAN_RayOptical_StraightQuadrature_Base*>(baseray);

	bool						ok = true;
	double						r0;
	double						r1;
	double						rt;
	double						t0;
	double						t1;
	double						sigma0;
	double						sigma1;
	HELIODETIC_POINT			startpoint;
	HELIODETIC_POINT			endpoint;
	//SKTRAN_OpticalDepthCalculator_LinearWithHeight	odcalculator;

	ok = ok && ray->GetQuadratureInterpParams( cellidx, &r0, &r1, &t0, &t1, &rt, &startpoint, &endpoint );
	//ok = ok && m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1( ray->GetWavelength(), startpoint, endpoint, &sigma0, &sigma1 );
	//ok = ok && GetEffectiveExtinctionPerCMWithHeight1(ray, startpoint, endpoint, &sigma0, &sigma1);
	ok = ok && GetEffectiveExtinctionPerCMWithHeight1(ray, cellidx, &sigma0, &sigma1);	
	//ok = ok && odcalculator.ConfigureQuadratureCoefficients( r0, r1, t0, t1, rt);
	if (ok)
	{
		//opticaldepth = odcalculator.OpticalDepthFromStartToEnd( sigma0, sigma1) ;
		ok = ok && GetOpticalDepthFromParams(r0, r1, t0, t1, rt, sigma0, sigma1, opticaldepth);
	}
	ok = ok && (opticaldepth >= 0.0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_OpticalPropertiesIntegrator_Straight::OpticalDepthOfSegment, Error looking up optical depth of a segment");
		opticaldepth = 0;
	}
	return opticaldepth * ray->Storage()->CellCurvature(cellidx);

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::CalculateRayTransmission		2013-06-11*/
/** **/
/*---------------------------------------------------------------------------*/

 bool SKTRAN_OpticalPropertiesIntegrator_Straight::CalculateRayScalarTransmission_withMinContainer   ( SKTRAN_RayOptical_Base* baseray, double* transmission, bool totaltransmissiononly, bool usecachedtransmission  ) const
 {
	bool					ok					= true;
	std::vector<double>*	odstorage			= baseray->OpticalDepthArrayVar();
	double					opticaldepthbuffer	= 0;
	size_t					npts;

	if( totaltransmissiononly )
	{
		odstorage->resize(1);
        const size_t numCells(baseray->GetNumCells()); // Loop invariant for MS auto-vectorizer
		for( size_t cellidx = 0; cellidx < numCells; cellidx++ )
		{
			opticaldepthbuffer += OpticalDepthOfCell_withMinCache( baseray, cellidx );
		}
		odstorage->at(0) = opticaldepthbuffer;
	}
	else
	{
		npts = baseray->GetNumQuadraturePoints();
		odstorage->resize(npts );
		if (npts > 0)
		{
			odstorage->at(0) = 0.0;
            const size_t numQuadPoints(baseray->GetNumQuadraturePoints()); // Loop invariant for MS auto-vectorizer
			for( size_t cellidx = 1; cellidx < numQuadPoints; cellidx++ )
			{
				opticaldepthbuffer += OpticalDepthOfCell_withMinCache( baseray, cellidx-1 );
				odstorage->at(cellidx) = opticaldepthbuffer;
			}
		}
	}

	return ok;
 }


/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::CalculateRayTransmission		2014-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_OpticalPropertiesIntegrator_Straight::CalculateRayScalarTransmission( SKTRAN_RayOptical_Base* baseray, double* transmission, bool totaltransmissiononly, bool usecachedtransmission  ) const
{
	bool					ok					= true;
	std::vector<double>*	odstorage			= baseray->OpticalDepthArrayVar();
	double					opticaldepthbuffer	= 0;

	if( totaltransmissiononly )
	{
		//straightray->m_quaddistances.resize(1);
		odstorage->resize(1);
        const size_t numCells(baseray->GetNumCells()); // Loop invariant for MS auto-vectorizer
		for( size_t cellidx = 0; cellidx < numCells; cellidx++ )
		{
			opticaldepthbuffer += OpticalDepthOfCell( baseray, cellidx );
		}
		odstorage->at(0) = opticaldepthbuffer;
	}
	else
	{
		odstorage->resize( baseray->GetNumQuadraturePoints() );
		odstorage->at(0) = 0.0;
        const size_t numQuadPoints(baseray->GetNumQuadraturePoints()); // Loop invariant for MS auto-vectorizer
		for( size_t cellidx = 1; cellidx < numQuadPoints; cellidx++ )
		{
			opticaldepthbuffer += OpticalDepthOfCell( baseray, cellidx-1 );
			odstorage->at(cellidx) = opticaldepthbuffer;
		}
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::CalculateRayTransmissionVector		2014-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_OpticalPropertiesIntegrator_Straight::CalculateRayScalarTransmissionVector( SKTRAN_RayOptical_Base* baseray, double* transmission, bool totaltransmissiononly, bool usecachedtransmission, std::vector<double>* sigmak, std::vector<double>* sigmaf  ) const
{
	nxLog::Record(NXLOG_ERROR, "SKTRAN_OpticalPropertiesIntegrator_Straight::CalculateRayScalarTransmissionVector is not implemented");
	return false;
}

bool SKTRAN_OpticalPropertiesIntegrator_Straight::GetOpticalDepthFromParams(double r0, double r1, double t0, double t1, double rt, double sigma0, double sigma1, double& opticaldepth) const
{
	SKTRAN_OpticalDepthCalculator_LinearWithHeight	odcalculator;

	odcalculator.ConfigureQuadratureCoefficients( r0,r1,t0,t1,rt);
	opticaldepth = odcalculator.OpticalDepthFromStartToEnd( sigma0, sigma1);

	return opticaldepth >= 0.0;

}
	




/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::SingleScatterRadiance		2014-1-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_OpticalPropertiesIntegrator_Straight::SingleScatterRadiance( SKTRAN_RayOptical_Base* baseray, double& radiance, const SKTRAN_Source_Term& source ) const
{
	bool						ok = true;
	HELIODETIC_POINT			centerpoint;
	double						sourcevalue;
	HELIODETIC_UNITVECTOR		look;
	double						ds;
	double						k;
	double						opticaldepth;
	std::vector<double>*		odstorage			= baseray->OpticalDepthArrayVar();

	radiance = 0;
	
	HELIODETIC_POINT obspt;
	baseray->Coordinates()->HelioVectorToHelioPoint( baseray->GetObserver(), &obspt);
	SKTRAN_SourceTermQueryObject_StraightPolarized qobj( obspt, baseray->LookVector() );
    const size_t numCells(baseray->GetNumCells()); // Loop invariant for MS auto-vectorizer 
	for( size_t cellidx = 0; cellidx < numCells; cellidx++ )
	{
		look = baseray->Storage()->AverageLookVectorAwayFromObserver( cellidx);
		ds   = baseray->Storage()->CellLength( cellidx);
		// check for very small path lengths
		if( ds < 1E-6 )
		{
			continue;
		}
		ok = ok && baseray->Storage()->CellMidPoint( cellidx, &centerpoint );
		qobj.UpdateQuery( centerpoint, look );
		ok = ok && source.SourceTermAtPoint( qobj, &sourcevalue );
		
		// this method gives identical to all decimal places to sasktranv21
		//k = m_opticalprops->TotalExtinctionPerCM( centerpoint ) * 100.0;
		k = TotalExtinctionPerCM(baseray, centerpoint) * 100.0;
		opticaldepth = odstorage->at( cellidx);
		radiance += (sourcevalue *(1 - exp(-1.0*k*ds)) / k) * exp(-1.0*opticaldepth);
	}

	// check for ground contributions
	if( baseray->Storage()->GroundIsHit() )
	{
		HELIODETIC_POINT groundpoint;
		
		baseray->Storage()->LocationOfPoint( baseray->Storage()->NumCells(), &groundpoint );

        qobj.UpdateQuery( groundpoint, baseray->Storage()->AverageLookVectorAwayFromObserver(numCells - 1));
		
		//coszen = -(groundpoint.LocalZenith() & straightray->GeometryRay()->LookVector());
        double groundsourcevalue;
        source.GroundSourceAtPoint( qobj, &groundsourcevalue ); 
        radiance += groundsourcevalue * exp(-1.0*baseray->TotalOpticalDepth());
	}
	
	return ok;
}




/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::FindNextScatterPosition		2013-06-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_OpticalPropertiesIntegrator_Straight::FindNextScatterPosition(	double									rand,
																			double									resolution,
																			SKTRAN_TableOpticalProperties_Base*		cachedOpticalProperties, 
																			SKTRAN_RayOptical_Base const*			originalRay_base,
																			HELIODETIC_VECTOR&						scatterVector,
																			double&									scatterProb, 
																			double&									targetOpticalDepth_toUser, 
																			double userImposedDistance ) const
{
	
	const SKTRAN_RayStorage_Base*			storage		= originalRay_base->Storage();
	const std::vector<double>&				odarray     = originalRay_base->OpticalDepthArray();
	bool									ok = true;
	double									scatterDistance;
	double									targetTau, tau;
	double									totalod;
	double									totalabsorption;
	std::vector<double>::const_iterator		exitod_iter;
	HELIODETIC_UNITVECTOR					stepLook;

	NXASSERT(( storage->NumQuadraturePoints() == odarray.size() ));

	if( storage->NumQuadraturePoints() == 0)
	{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_OpticalPropertiesIntegrator_Straight_MC::FindNextScatterDistance, transmission array has size zero. ");
		return false;
	}

	
	totalod          = odarray.back();
	totalabsorption  = 1.0 - exp(-totalod);
	targetTau        = -log(1 - rand*(totalabsorption));
	targetOpticalDepth_toUser = targetTau;

	//dbtarget = targetTau;
	
	//scatterProb = 1e-16 < originalRay->TotalOpticalDepth() ? (1.0-exp(-targetTau)) / (1.0 - exp(-originalRay->TotalOpticalDepth())) : 0.0; // Normalized pdf
	scatterProb = exp(-targetTau) / (totalabsorption); // Still needs to be multiplied by extinction at scatter point

	// perform binary search for index of shell previous to scattering point
	if(0.0==targetTau)
	{
		exitod_iter = odarray.begin()+1;
	}
	else if (targetTau >= totalod)
	{
		// A weird special case that can happen when the ray optical depth is close to machine precision and rounding errors occur when transforming
		// the probability space
		targetTau = totalod;
		exitod_iter = odarray.end() - 1;
		targetTau -= *(exitod_iter - 1);	// Calculate the optical depth increase seen inside this shell only
	}
	else
	{
		exitod_iter= std::lower_bound( odarray.begin(), odarray.end(), targetTau);
		targetTau -= *(exitod_iter-1);	// Calculate the optical depth increase seen inside this shell only
	}

	// Perform "binary search" to find location of scatter point
	double dt;
	double sigma0, sigma1;
	double t0, t1, r0, r1, rt;
	//SKTRAN_OpticalDepthCalculator_LinearWithHeight	odcalculator;

	//scatShell = 10*(exitInterceptOpticalDepth - originalRay->m_opticaldepth.begin());
	size_t endinterceptindex   = (exitod_iter - odarray.begin());
	size_t tangentIndex        = exitod_iter - odarray.begin()-1;

	if( endinterceptindex == 0 )
	{
		tangentIndex = 0;
		endinterceptindex = 1;
	}

	size_t startinterceptindex = endinterceptindex - 1;
	
	HELIODETIC_POINT			startpoint;
	HELIODETIC_POINT			endpoint;
	//originalRay->GeometryRay()->LocationAlongRayAsPoint( *startintercept, &startpoint );
	//originalRay->GeometryRay()->LocationAlongRayAsPoint( *endintercept, &endpoint );
	
	//originalRay->GetCellQuadraturePoint(exitInterceptOpticalDepth-originalRay->m_opticaldepth.cbegin()-1,startpoint);
	//originalRay->GetCellQuadraturePoint(exitInterceptOpticalDepth-originalRay->m_opticaldepth.cbegin(),endpoint);
	ok = ok && storage->LocationOfPoint( startinterceptindex, &startpoint);
	ok = ok && storage->LocationOfPoint( endinterceptindex ,  &endpoint);
	//ok = ok && m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1(originalRay_base->GetWavelength(), startpoint, endpoint, &sigma0, &sigma1 );
	ok = ok && GetEffectiveExtinctionPerCMWithHeight1(originalRay_base, startpoint, endpoint, &sigma0, &sigma1);


	rt = storage->RadiusOfCellTangentPoint(0);
	r0 = startpoint.Radius();
	r1 = endpoint.Radius();
	t0 = storage->DistanceOfPointFromCellTangentPoint( startinterceptindex, startinterceptindex);
	t1 = storage->DistanceOfPointFromCellTangentPoint( endinterceptindex,   startinterceptindex);
	//dt = t1 - t0;
	
	//odcalculator.ConfigureQuadratureCoefficients( r0,r1,t0,t1,rt);
	//while(fabs(dt)>resolution)
	//{
	//	dt /= 2.0;					// decrease step size
	//	t1 -= dt;					// step forwards along look direction
	//	r1 = sqrt(rt*rt + t1*t1);	// calculate radius at this point
	//	sigma1 = odcalculator.SigmaAtRadius(r1, sigma0, sigma1);				// Get the cross-section at this new point
	//	odcalculator.ConfigureQuadratureCoefficients( r0,r1,t0,t1,rt);			// reconfigure the optical depth calculation
	//	tau = odcalculator.OpticalDepthFromStartToEnd( sigma0, sigma1);			// And do the calculation

	//	//if( dt > 0.0 ){		// case 1, upward rays, or very short cell segments at the tangent point
	//	//	tau = sigmak*(t1-t0) + sigmaf*0.5*( (r1*t1-r0*t0) + rt*rt*log( (r1+t1)/(r0+t0) ) );
	//	//} else{				// case 2. downward ray, passes from upper shell to lower shell
	//	//	tau = sigmak*(t0-t1) + sigmaf*0.5*( (r0*t0-r1*t1) + rt*rt*log( (r0+t0)/(r1+t1) ) );
	//	//}
	//	if(targetTau > tau) t1 += dt;
	//}
	//// scatterProb *= sigmak + sqrt(rt*rt+t1*t1)*sigmaf; // Slope; old way of doing this

	//scatterDistance = fabs(t1-t0);

	ok = ok && CalculateDistanceToTargetOpticalDepth(r0, r1, t0, t1, rt, sigma0, sigma1, resolution, targetTau, scatterDistance);

	if( 0.0 <= userImposedDistance ){
		scatterDistance = userImposedDistance - storage->DistanceOfPointFromOrigin( startinterceptindex );
	}

	//scatShell += (size_t)(floor(10*fabs(t1-t0) / (*endintercept-*startintercept)));
	//dbscatdistance = scatterDistance;
	//scatterVector = originalRay->GeometryRay()->GetObserver() + HELIODETIC_VECTOR(originalRay->GeometryRay()->LookVector(), scatterDistance);
	stepLook      = storage->AverageLookVectorAwayFromObserver( tangentIndex );

	ok = ok && RoundScatterPosition(startpoint, stepLook, scatterDistance, storage->GetCoordsPtr(), scatterVector, endpoint); // move the scatter point back slightly if on the edge of the atmosphere		

	//scatterProb *= m_opticalprops->TotalExtinctionPerCM(originalRay_base->GetWavelength(), endpoint ) * 100.0;
	scatterProb *= TotalExtinctionPerCM(originalRay_base, endpoint) * 100.0;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_OpticalPropertiesIntegrator_Straight::FindNextScatterPosition, Error locating next scatter position. That might cause problems");
	}
	return ok;
}


bool SKTRAN_OpticalPropertiesIntegrator_Straight::CalculatePartialOpticalDepth(SKTRAN_RayOptical_Base const * ray, size_t cellindex, double fraction, double & opticaldepth) const
{
	bool ok = true;
	double sigma0, sigma1;
	const SKTRAN_RayStorage_Base* storage = ray->Storage();

	
	HELIODETIC_POINT pt0, pt1;
	ok = ok && storage->LocationOfPoint(cellindex, &pt0);
	ok = ok && storage->LocationOfPoint(cellindex + 1, &pt1);
	pt1.FromVector(HELIODETIC_VECTOR(pt0.UnitVector(), pt0.Radius() * (1.0 - fraction)) + HELIODETIC_VECTOR(pt1.UnitVector(), pt1.Radius() * fraction), m_opticalprops->CoordinatesPtr());

	double rt = storage->RadiusOfCellTangentPoint(cellindex);
	double r0 = pt0.Radius();
	double r1 = pt1.Radius();
	double t0 = storage->DistanceOfPointFromCellTangentPoint(cellindex, cellindex);
	double t1 = storage->DistanceOfPointFromCellTangentPoint(cellindex + 1, cellindex + 1);
	t1 = (1.0 - fraction) * t0 + fraction * t1;

	ok = ok && GetEffectiveExtinctionPerCMWithHeight1(ray, pt0, pt1, &sigma0, &sigma1);

	ok = ok && GetOpticalDepthFromParams(r0, r1, t0, t1, rt, sigma0, sigma1, opticaldepth);

	return ok;	
}


bool SKTRAN_OpticalPropertiesIntegrator_Straight::CalculateDistanceToTargetOpticalDepth(double r0, double r1, double t0, double t1, double rt, double sigma0, double sigma1, double resolution, double targetTau, double & scatterDistance) const
{
	bool ok = true;

	double tau;
	double 	dt = t1 - t0;
	SKTRAN_OpticalDepthCalculator_LinearWithHeight	odcalculator;
	ok = ok && odcalculator.ConfigureQuadratureCoefficients(r0, r1, t0, t1, rt);
	while (fabs(dt) > resolution)
	{
		dt /= 2.0;					// decrease step size
		t1 -= dt;					// step forwards along look direction
		r1 = sqrt(rt*rt + t1 * t1);	// calculate radius at this point
		sigma1 = odcalculator.SigmaAtRadius(r1, sigma0, sigma1);				// Get the cross-section at this new point
		ok = ok && odcalculator.ConfigureQuadratureCoefficients(r0, r1, t0, t1, rt);			// reconfigure the optical depth calculation
		tau = odcalculator.OpticalDepthFromStartToEnd(sigma0, sigma1);			// And do the calculation

		//if( dt > 0.0 ){		// case 1, upward rays, or very short cell segments at the tangent point
		//	tau = sigmak*(t1-t0) + sigmaf*0.5*( (r1*t1-r0*t0) + rt*rt*log( (r1+t1)/(r0+t0) ) );
		//} else{				// case 2. downward ray, passes from upper shell to lower shell
		//	tau = sigmak*(t0-t1) + sigmaf*0.5*( (r0*t0-r1*t1) + rt*rt*log( (r0+t0)/(r1+t1) ) );
		//}
		if (targetTau > tau) t1 += dt;
	}
	// scatterProb *= sigmak + sqrt(rt*rt+t1*t1)*sigmaf; // Slope; old way of doing this

	scatterDistance = fabs(t1 - t0);
	return ok;
}

bool SKTRAN_OpticalPropertiesIntegrator_ConstantLayers::GetOpticalDepthFromParams(double r0, double r1, double t0, double t1, double rt, double sigma0, double sigma1, double & opticaldepth) const
{
	double sigma = r0 > r1 ? sigma1 : sigma0; // the lower point should represent the state of the layer unless the radii are the same (which shouldn't happen since rays are split at tangent points)
	opticaldepth = 1e2 * fabs(t1 - t0) * sigma;
	return true;
}

bool SKTRAN_OpticalPropertiesIntegrator_ConstantLayers::CalculateDistanceToTargetOpticalDepth(double r0, double r1, double t0, double t1, double rt, double sigma0, double sigma1, double resolution, double targetTau, double & scatterDistance) const
{
	bool ok = true;

	double sigma = r0 > r1 ? sigma1 : sigma0;
	ok = ok && sigma > 0.0;
	scatterDistance = targetTau / (1e2 * sigma);
	ok = ok && scatterDistance < fabs(t1 - t0);
	return ok;
}

