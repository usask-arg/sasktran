#include "../sktran_common.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_SourceTermIntegrator_Order0::IntegrateSourceTerm		 2014- 12- 3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SourceTermIntegrator_Order0::IntegrateSourceTerm ( const SKTRAN_RayOptical_Base* ray, SKTRAN_Stokes_NC& radiance, const std::vector<SKTRAN_Source_Term*>& sources ) const 
{
	Integration_Impl<SKTRAN_Stokes_NC> temp;
    return temp.IntegrateSourceTerm(ray, radiance, sources, m_opticalprops); 
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SourceTermIntegrator_Order0::IntegrateSourceTerm		 2014- 12- 3*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_SourceTermIntegrator_Order0::IntegrateSourceTerm ( const SKTRAN_RayOptical_Base* ray, double& radiance, const std::vector<SKTRAN_Source_Term*>& sources ) const 
{
	Integration_Impl<double> temp; 
    return temp.IntegrateSourceTerm(ray, radiance, sources, m_opticalprops); 
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_SourceTermIntegrator_Order0::Impl<radtype>::IntegrateSourceTerm		 2014- 12- 3*/
/** **/
/*---------------------------------------------------------------------------*/

template<typename radtype>
bool  SKTRAN_SourceTermIntegrator_Order0::Integration_Impl<radtype>::IntegrateSourceTerm( const SKTRAN_RayOptical_Base* baseray, radtype& radiance, const std::vector<SKTRAN_Source_Term*>& sources, const SKTRAN_TableOpticalProperties_Base* optprops ) const
{
	bool ok = true;
	
	HELIODETIC_POINT			centerpoint;
	radtype						source;
	radtype						tempsource;
	HELIODETIC_UNITVECTOR		look;

	double						ds;

	double						k;
	double						opticaldepth;
	const SKTRAN_RayStorage_Base*			storage		= baseray->Storage();
	const std::vector<double>&				odarray     = baseray->OpticalDepthArray();
    const size_t numSources(sources.size()); // Loop invariant for MS auto-vectoriser 
	
	skRTStokesVector::SetToZero(radiance);
	skRTStokesVector::SetToZero(source);
	skRTStokesVector::SetToZero(tempsource);
	HELIODETIC_POINT obspt;
	baseray->Coordinates()->HelioVectorToHelioPoint( baseray->GetObserver(), &obspt );
	SKTRAN_SourceTermQueryObject_StraightPolarized qobj(obspt, baseray->LookVector() );
	for( size_t cellidx = 0; cellidx < baseray->GetNumCells(); cellidx++ )
	{
		look = storage->AverageLookVectorAwayFromObserver( cellidx );
		ds   = storage->CellLength( cellidx );
		//// check for very small path lengths
		//if( ds < 1E-6 )
		//{
		//	continue;
		//}
		storage->CellMidPoint( cellidx, &centerpoint );
		skRTStokesVector::SetToZero(source);
		qobj.UpdateQuery( centerpoint, look );
		for( size_t sourceidx = 0; sourceidx < numSources; sourceidx++ )
		{

			ok = ok && sources[sourceidx]->SourceTermAtPoint( qobj,  &tempsource );
			source += tempsource;
		}
		
		// this method gives identical to all decimal places to sasktranv21
		k = optprops->TotalExtinctionPerCM( centerpoint ) * 100;
		opticaldepth = odarray.at(cellidx);
		radiance += (source  * ((1 - exp(-1.0*k*ds)) / k)) * exp(-1.0*opticaldepth);
	}

	// add ground contributions
	if( baseray->Storage()->GroundIsHit() )
	{
		radtype groundsource;
		skRTStokesVector::SetToZero(groundsource);
		HELIODETIC_POINT groundpoint;

		baseray->Storage()->LocationOfPoint( baseray->Storage()->NumCells(), &groundpoint);
		qobj.UpdateQuery( groundpoint, look );
		for( size_t sourceidx = 0; sourceidx < numSources; sourceidx++ )
		{
			ok = ok && sources[sourceidx]->GroundSourceAtPoint( qobj, &tempsource );
			groundsource += tempsource;
		}

		radiance += groundsource *  exp(-1.0*baseray->TotalOpticalDepth());
	}
	return ok;
}





/* Order 2 code */

/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Straight::IntegrateSourceTerm		2014-1-25*/
/** Integrate source term assuming the source at the center of a segment is well-represented
    by a quadratic fit to the source intensity at the end- and mid-points
**/
/*---------------------------------------------------------------------------*/

template<typename radtype>
bool SKTRAN_SourceTermIntegrator_Order2::IntegrateSourceTerm_Impl( const SKTRAN_RayOptical_Base* baseray, radtype& radiance, const std::vector<SKTRAN_Source_Term*>& sources ) const
{
	bool								ok			= true;
	const SKTRAN_RayOptical_Straight*	ray			= reinterpret_cast< const SKTRAN_RayOptical_Straight*>(baseray);
	HELIODETIC_POINT					point;
	HELIODETIC_POINT					centerpoint;
	//HELIODETIC_POINT					startpoint;
	HELIODETIC_UNITVECTOR				look;
	double								ds;
	double								opticaldepth_start;
	double								opticaldepth_end;
	radtype								tempsource;
	radtype								sourcestart;
	radtype								sourcemid;
	radtype								sourceend;
	radtype                             radiance_cell;

	skRTStokesVector::SetToZero( radiance    );
	skRTStokesVector::SetToZero( tempsource  );
	skRTStokesVector::SetToZero( sourcestart );
	skRTStokesVector::SetToZero( sourcemid   );
	skRTStokesVector::SetToZero( sourceend   );

	HELIODETIC_POINT obspt;
	ray->Coordinates()->HelioVectorToHelioPoint(ray->GetObserver(), &obspt);
	SKTRAN_SourceTermQueryObject_StraightPolarized qPoint ( obspt, ray->LookVector() );
	SKTRAN_SourceTermQueryObject_StraightPolarized qCenter( obspt, ray->LookVector() );

	if( 0 < ray->GetNumQuadraturePoints() ){
		// need to initialize the start source
		ok   = ok && ray->Storage()->LocationOfPoint( 0, &point );
		look = ray->Storage()->AverageLookVectorAwayFromObserver( 0 );

		qPoint.UpdateQuery( point, look );

		for( size_t sourceidx = 0; sourceidx < sources.size(); sourceidx++ )
		{
			ok = ok && sources[sourceidx]->SourceTermAtPoint( qPoint, &tempsource );
			sourcestart += tempsource;
		}

		for( size_t cellidx = 0; cellidx < ray->GetNumCells(); cellidx++ )
		{
			look =  ray->Storage()->AverageLookVectorAwayFromObserver( cellidx );
			ds   =  ray->Storage()->CellLength( cellidx );
			// check for very small path lengths
			if( ds < 1E-6 )
			{
				continue;
			}
			opticaldepth_start = baseray->OpticalDepthArray().at( cellidx );	 //GetOpticalDepthAtCell
			opticaldepth_end   = baseray->OpticalDepthArray().at( cellidx+1 );
			//ok = ok && baseray->Storage()->LocationOfPoint( cellidx,   &startpoint );
			ok = ok && baseray->Storage()->LocationOfPoint( cellidx+1, &point );
			ok = ok && baseray->Storage()->CellMidPoint   ( cellidx,   &centerpoint );
			qPoint.UpdateQuery( point, look);
			qCenter.UpdateQuery( centerpoint, look );
			skRTStokesVector::SetToZero( sourceend );
			skRTStokesVector::SetToZero( sourcemid );
			for( size_t sourceidx = 0; sourceidx < sources.size(); sourceidx++ )
			{
				ok = ok && sources[sourceidx]->SourceTermAtPoint( qPoint, &tempsource );

				sourceend += tempsource;
				ok = ok && sources[sourceidx]->SourceTermAtPoint( qCenter,  &tempsource );

				sourcemid += tempsource;
			
			}
			// attenuate back to the observer
			radiance_cell = GetRadianceContrib( sourcestart, sourcemid, sourceend, ds, opticaldepth_start, opticaldepth_end ) * exp(-1.0*opticaldepth_start);
			radiance += radiance_cell;

			ok = ok && ( !skRTStokesVector::IsNegative(radiance) );
			if( !ok )
			{
				nxLog::Record(NXLOG_WARNING, "SKTRAN_OpticalPropertiesIntegrator_Adaptive::IntegrateSourceTerm Radiance Less than 0" );
				break;
			}
			sourcestart = sourceend;
		}
	}

	// check for ground contributions
	if( baseray->Storage()->GroundIsHit() )
	{
		radtype groundsource;
		skRTStokesVector::SetToZero( groundsource );

		HELIODETIC_POINT groundpoint;
		ray->Storage()->LocationOfPoint( ray->Storage()->NumCells(), &groundpoint);
		qPoint.UpdateQuery( groundpoint, look );
		for( size_t sourceidx = 0; sourceidx < sources.size(); sourceidx++ )
		{
			ok = ok && sources[sourceidx]->GroundSourceAtPoint( qPoint, &tempsource );
			groundsource += tempsource;
		}

		radiance += groundsource *  exp(-1.0*ray->TotalOpticalDepth());	
	}
	
	return ok;
}



bool SKTRAN_SourceTermIntegrator_Order2::IntegrateSourceTerm( const SKTRAN_RayOptical_Base* baseray, SKTRAN_Stokes_NC& radiance, const std::vector<SKTRAN_Source_Term*>& sources ) const
{
	return IntegrateSourceTerm_Impl( baseray, radiance, sources );
}

bool SKTRAN_SourceTermIntegrator_Order2::IntegrateSourceTerm( const SKTRAN_RayOptical_Base* baseray, double& radiance, const std::vector<SKTRAN_Source_Term*>& sources ) const
{
	return IntegrateSourceTerm_Impl( baseray, radiance, sources );
}


template<typename radtype>
bool SKTRAN_SourceTermIntegrator_Order2::GetQuadraticCoeff( const radtype& sourcestart, const radtype& sourcemid, const radtype& sourceend, double ds, radtype& b, radtype& c ) const
{
	bool ok = true;

	b = ( sourcestart*3.0 + sourcemid*(-4.0) + sourceend ) * ( -1.0 / ds );
	c = ( sourcestart + sourcemid*(-2.0) + sourceend ) * ( 2.0 / (ds*ds) );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_OpticalPropertiesIntegrator_Adaptive::GetRadianceContrib		2013-08-22*/
/** Calculates the (unattenuated) radiance contribution from a cell.  To get the actual 
 *  radiance this value needs to be attenuated by the transmission at the START of
 *  the cell.
 *
 *  ------------- OLD SASKTRANV21 METHOD ---------------
 *  Sample the source term at the mid point of the cell, J
 *  Sample the extinction at the mid point of the cell,  k
 *  Perform the integral over the cell
 *     int( exp(-k*s) * J ds )
 *
 *  -------------- NEW METHOD            ----------------
 *  Sample the source term at the start, mid, and end points of the cell
 *  Use these values to find J = a + b*s + c*s^2 , i.e. J as a function of
 *  distance along the cell. (s=0 at start of cell)
 *  Calculate k as the average extinction over the cell, i.e.
 *     k = (opticaldepthofcell)/ds
 *  Perform the integral over the cell
 *     int( exp(-k*s) * (a + b*s + c*s^2) )
 *
 *  Notes on implentation:
 *  In general the integral for s^n can be found analytically, however there
 *  are precision problems that can be run into.  Usually the exact form for
 *  these integrals involves terms similar to (1 - exp(-somethingcloseto0)) * somethignbig
 *  which causes precision errors on 64bit.  The s^0 and s^1 terms can be
 *  evaluated fine on 64bit precision for all reasonable values of k, however
 *  the s^2 term can be off by a very large margin (getting answers like
 *  10^-5 when they should be 10^-13)  To fix this, if the value of k*ds
 *  is too small (the region where the errors occur), we use 7 terms in
 *  the taylor expansion of the answer rather than the exact formula. 
 *  This fixes the precision errors.
 *
 *  Validity of Approximation:
 *  The approximation was tested inside Mathematica.  The exact formula was compiled
 *  as C code within mathematica, which exhibited the same numerical instability
 *  described above.  The exact answer to arbritrary precision and the taylor
 *  result were also both computed inside mathematica.
 **/
/*---------------------------------------------------------------------------*/

template<typename radtype>
radtype SKTRAN_SourceTermIntegrator_Order2::GetRadianceContrib( const radtype& sourcestart, const radtype& sourcemid, const radtype& sourceend, double ds, double sigmastart, double sigmaend ) const
{
	radtype a,b,c;
	radtype sourceconst;
	radtype sourcelinear;
	radtype sourcequad;
	radtype result;

	GetQuadraticCoeff( sourcestart, sourcemid, sourceend, ds, b, c );
	a = sourcestart;

	double k = (sigmaend - sigmastart)/ ds;
	sourceconst  = a * ( (1 - exp(-1.0*k*ds)) / k);

	// the exact formula causes precision errors, so use an approximation
	// if it is valid
	if( ds * k < 0.01 ){
        double ds1k1 ( ds*k        );
        double ds2k0 ( ds*ds       );
        double ds3k1 ( ds2k0*ds1k1 );
        double ds4k2 ( ds3k1*ds1k1 );
        double ds5k3 ( ds4k2*ds1k1 );
        double ds6k4 ( ds5k3*ds1k1 );
        double ds7k5 ( ds6k4*ds1k1 );
        sourcelinear = b * (    ds2k0/2.0  -     ds3k1/3.0  +     ds4k2/ 8.0  -    ds5k3/30.0  +    ds6k4/144.0  -    ds7k5/840.0 ); // +  (ds7k5*ds1k1)/5760.0 );
		sourcequad   = c * ( ds*ds2k0/3.0  -  ds*ds3k1/4.0  +  ds*ds4k2/10.0  - ds*ds5k3/36.0  + ds*ds6k4/180.0  - ds*ds7k5/720.0 );
    } else{
        sourcelinear = b * ( (1.0 - exp(-1.0*ds*k)*(ds*k+1)) / (k*k) );
		sourcequad   = c * ( (2.0 + exp(-1.0*ds*k)*( -2.0 - ds*k*(2.0 + ds*k)))/ (k*k*k) );
    }
    
    result = sourceconst + sourcelinear + sourcequad;
	skRTStokesVector::SetNegativesToZero(result);
	skRTStokesVector::SetNansToZero(result);

	return result;
}

