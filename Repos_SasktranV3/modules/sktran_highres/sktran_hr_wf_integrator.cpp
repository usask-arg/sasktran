#include "include/sktran_hr_internals.h"

//bool SKTRAN_HR_WF_Integrator::CalculateWeightingFunctions_DirectEffectsOnly( const SKTRAN_RayOptical_Base& raybase,
//														   const SKTRAN_HR_WF_Store& perturbations,
//														   const std::vector< SKTRAN_Source_Term* >& sources,
//																 SKTRAN_HR_WF_Ray& wf ) const
//{
//	bool ok = true;
//	double kperturb;
//	const SKTRAN_RayOptical_Straight&		ray			= dynamic_cast< const SKTRAN_RayOptical_Straight&>(raybase);
//
//	double ds;
//	HELIODETIC_POINT midpoint;
//	HELIODETIC_UNITVECTOR look;
//	double source, solarsource, diffusesource;
//	double singlescatt;
//	double tempsource;
//	double k;
//	double opticaldepth;
//	const SKTRAN_HR_Diffuse_Source* difftable = nullptr;
//	if (sources.size() > 1)
//	{
//		difftable = dynamic_cast<const SKTRAN_HR_Diffuse_Source*> (sources[1]);
//	}
//	const SKTRAN_SolarTransmission_Base* solartable = dynamic_cast<const SKTRAN_SolarTransmission_Base*> (sources[0]);
//
//    HELIODETIC_POINT obspt;
//    ray.Coordinates()->HelioVectorToHelioPoint( ray.GetObserver(), &obspt );
//    SKTRAN_SourceTermQueryObject_Simple qMidpoint ( obspt, ray.LookVector() );
//
//	std::vector<double> pertvals;
//	std::vector<double> pertalts;
//	std::vector<double> pertwidth;
//	std::vector<HELIODETIC_POINT> pertlocations;
//	pertalts.resize(perturbations.StoreSize());
//	perturbations.PerturbationLocation(*ray.Coordinates(), pertlocations);
//	perturbations.PerturbationAltitudeWidth(pertwidth);
//	for( size_t idx = 0; idx < pertlocations.size(); idx++ )
//	{
//		pertalts[idx] = pertlocations[idx].Altitude();
//	}
//	wf.Allocate( perturbations.StoreSize() );
//
//	for( size_t cellidx = 0; cellidx < ray.GetNumCells(); cellidx++ )
//	{
//		look = ray.Storage()->AverageLookVectorAwayFromObserver( cellidx);
//		ds   = ray.Storage()->CellLength( cellidx);
//		if( ds < 1E-6 )
//			continue;
//		ok = ok && ray.Storage()->CellMidPoint( cellidx, &midpoint );
//
//		source = 0;
//		solarsource = 0;
//		diffusesource = 0;
//        qMidpoint.UpdateQuery( midpoint, look );
//		for( size_t sourceidx = 0; sourceidx < sources.size(); sourceidx++ )
//		{
//			ok = ok && sources[sourceidx]->SourceTermAtPoint( qMidpoint, &tempsource );
//			source += tempsource;
//			if( sourceidx == 0)
//				solarsource = tempsource;
//			if( sourceidx == 1)
//				diffusesource = tempsource;
//		}
//		
//		opticaldepth = ray.OpticalDepthArray().at(cellidx);
//		k = GetOpticalProps()->TotalExtinctionPerCM( midpoint ) * 100;
//		GetOpticalProps()->GetScatteringCoefficientCM2( midpoint, look.Z(), &singlescatt );
//		singlescatt = singlescatt * 100;
//		double cellalt = midpoint.Altitude();
//		perturbations.ExtinctionPerturbation(midpoint, pertvals);
//
//		// Calculate effictive cell length ( if opticall thin, this is ds, if optically thick this decreases)
//		double solartransmission;
//		double effectivecelllength;
//		double effectivecelllength_secondorder;
//		if (k > 0)
//		{
//			effectivecelllength = (1 - exp(-1.0*k*ds)) / k;
//			effectivecelllength_secondorder = (1 - exp(-ds*k) * (1 + ds*k)) / (k*k);
//			solartransmission = solarsource / singlescatt;
//		}
//		else
//		{
//			// This can happen if there is nothing in the cell
//			effectivecelllength = ds;
//			effectivecelllength_secondorder = ds*ds;
//			solartable->TransmissionAtPoint(midpoint, solartransmission);
//		}
//
//		for( size_t pertidx = 0; pertidx < perturbations.StoreSize(); pertidx++ )
//		{
//			double pertalt = pertalts[pertidx];		
//			kperturb = pertvals[pertidx];
//			wf.WeightAt( pertidx ) += -1.0*(wf.AddedOpticalDepth( pertidx )) * (source * effectivecelllength) * exp(-1.0*opticaldepth);
//			if( kperturb > 0 )
//			{
//				if( m_isscatterer )
//				{
//					double phasefunction = PhaseFunction( pertidx, look.Z() );
//					double integrationresult = 0;
//					if( difftable )
//					{
//						std::function<double(const nxVector&, const nxVector&)> phasefp = [&](const nxVector& a, const nxVector& b) {
//							double cosangle = a.Dot(b);
//							return PhaseFunction(pertidx, cosangle);
//						};
//						integrationresult = difftable->IntegrateScalarIncoming(qMidpoint, phasefp);
//						if (integrationresult != integrationresult)
//						{
//							// This can happen due to some divisions when k=0
//							integrationresult = 0.0;
//						}
//					}
//					//wf.WeightAt( pertidx ) += kperturb * ( solarsource * phasefunction + (1 - pertalt/100000) * diffusesource /(4.0*nxmath::Pi)) / singlescatt * (1-exp(-1.0*k*ds))/k * exp(-1.0*opticaldepth);
//					wf.WeightAt( pertidx ) += kperturb * (( solartransmission * phasefunction) + integrationresult) * effectivecelllength * exp(-1.0*opticaldepth);
//				
//				}
//				wf.WeightAt( pertidx ) += source * (-1.0* kperturb) * effectivecelllength_secondorder * exp(-1.0*opticaldepth);
//				wf.AddedOpticalDepth( pertidx ) += kperturb * ds;
//			}
//		}
//		
//	}
//
//	if (ray.Storage()->GroundIsHit())
//	{
//		double totalgroundsource = 0.0;
//		double totalopticaldepth = ray.TotalOpticalDepth();
//
//		HELIODETIC_POINT groundpt;
//		ray.Storage()->LocationOfPoint(ray.Storage()->NumQuadraturePoints() - 1, &groundpt);
//
//		SKTRAN_SourceTermQueryObject_Simple qground(groundpt, ray.LookVector());
//		for (size_t sourceidx = 0; sourceidx < sources.size(); sourceidx++)
//		{
//			ok = ok && sources[sourceidx]->GroundSourceAtPoint(qground, &tempsource);
//			totalgroundsource += tempsource;
//		}
//
//		for (size_t pertidx = 0; pertidx < perturbations.StoreSize(); pertidx++)
//		{
//			double addedoptical = wf.AddedOpticalDepth(pertidx);
//
//			wf.WeightAt(pertidx) += -1 * addedoptical * totalgroundsource * exp(-1.0*totalopticaldepth);
//		}
//	}
//
//	return ok;
//}


bool SKTRAN_HR_WF_Integrator::CachePerturbationInformation(const SKTRAN_CoordinateTransform_V2&								coords,
														   const SKTRAN_HR_WF_Store&	perturbations)
{
	std::vector<HELIODETIC_POINT> pertlocations;
	m_pertalts.resize(perturbations.StoreSize());
	perturbations.PerturbationLocation(coords, pertlocations);
	perturbations.PerturbationAltitudeWidth(m_pertwidths);
	for (size_t idx = 0; idx < pertlocations.size(); idx++)
	{
		m_pertalts[idx] = pertlocations[idx].Altitude();
	}
	return true;
}



bool SKTRAN_HR_WF_Integrator::CalculateWeightingFunctions( const SKTRAN_RayOptical_Base& raybase,
														   const SKTRAN_HR_WF_Store& perturbations,
														   const std::vector< SKTRAN_Source_Term* >& sources,
	                                                       const std::vector<std::unique_ptr<SKTRAN_HR_WF_SpeciesInformationBase>>& wfinfo,
	                                                       std::vector<SKTRAN_HR_WF_Ray>& wfs,
													       bool directeffectsonly ) const
{
	bool ok = true;
	double kperturb;
	const SKTRAN_RayOptical_StraightQuadrature_Base&		ray			= dynamic_cast< const SKTRAN_RayOptical_StraightQuadrature_Base&>(raybase);
	std::unique_ptr<SKTRAN_RayOptical_Base> solarray;
	std::unique_ptr<SKTRAN_RayOptical_Base> groundray;

	HELIODETIC_UNITVECTOR sunlook, groundlook;
	sunlook.SetCoords(0,0,1);

	std::vector<double> addedsolaropticaldepth;
	if (directeffectsonly == false)
	{
		addedsolaropticaldepth.resize(perturbations.StoreSize());
	}

	m_rayfactory->CreateRayObject( &solarray );
	m_rayfactory->CreateRayObject( &groundray );

	std::vector<double> pertvals;
	double ds;
	HELIODETIC_POINT midpoint;
	HELIODETIC_UNITVECTOR look;
	double source, solarsource, diffusesource;
	double singlescatt;
	double tempsource;
	double k;
	double opticaldepth;

	// TODO: go through the sources and try to find these rather than assuming the order
	const SKTRAN_HR_Diffuse_Source* difftable = nullptr;
	if (sources.size() > 1)
	{
		difftable = dynamic_cast<const SKTRAN_HR_Diffuse_Source*> (sources[1]);
	}
	const SKTRAN_SolarTransmission_Base* solartable = dynamic_cast<const SKTRAN_SolarTransmission_Base*> (sources[0]);

    HELIODETIC_POINT obspt;
    ray.Coordinates()->HelioVectorToHelioPoint( ray.GetObserver(), &obspt );

	SKTRAN_SourceTermQueryObject_StraightPolarized qMidpoint(obspt, ray.LookVector());

	for (int idspecies = 0; idspecies < wfs.size(); idspecies++)
	{
		wfs[idspecies].Allocate(perturbations.StoreSize());
	}

	for( size_t cellidx = 0; cellidx < ray.GetNumCells(); cellidx++ )
	{
		look = ray.Storage()->AverageLookVectorAwayFromObserver( cellidx);
		ds   = ray.Storage()->CellLength( cellidx);
		if( ds < 1E-6 )
			continue;
		ok = ok && ray.Storage()->CellMidPoint(cellidx, &midpoint);

		if (directeffectsonly == false)
		{
			ok = ok && solarray->MoveObserver(midpoint.Vector(), sunlook);
			ok = ok && solarray->TraceRay_NewMethod();

			AddedOpticalDepth(perturbations, addedsolaropticaldepth, *solarray);
		}
		source = 0;
		solarsource = 0;
		diffusesource = 0;
        qMidpoint.UpdateQuery( midpoint, look );
		for( size_t sourceidx = 0; sourceidx < sources.size(); sourceidx++ )
		{
			ok = ok && sources[sourceidx]->SourceTermAtPoint( qMidpoint, &tempsource );
			source += tempsource;
			if( sourceidx == 0)
				solarsource = tempsource;
			if( sourceidx == 1)
				diffusesource = tempsource;
		}
		
		opticaldepth = ray.OpticalDepthArray().at(cellidx);
		k = GetOpticalProps()->TotalExtinctionPerCM( midpoint ) * 100;
		double kscat = GetOpticalProps()->ScatteringExtinctionPerCM( midpoint ) * 100;
		GetOpticalProps()->GetScatteringCoefficientCM2( midpoint, look.Z(), &singlescatt );
		singlescatt = singlescatt * 100;
		double cellalt = midpoint.Altitude();
		perturbations.ExtinctionPerturbation(midpoint, pertvals);

		// Calculate effictive cell length ( if opticall thin, this is ds, if optically thick this decreases)
		// Also get the solar transmission
		double solartransmission;
		double effectivecelllength;
		double effectivecelllength_secondorder;
		if (k > 0)
		{
			effectivecelllength = (1 - exp(-1.0*k*ds)) / k;
			effectivecelllength_secondorder = (1 - exp(-ds*k) * (1 + ds*k)) / (k*k);
			solartransmission = solarsource / singlescatt;
		}
		else
		{
			// This can happen if there is nothing in the cell
			effectivecelllength = ds;
			effectivecelllength_secondorder = ds*ds;
			solartable->TransmissionAtPoint(midpoint, solartransmission);
		}


		for (size_t pertidx = 0; pertidx < perturbations.StoreSize(); pertidx++)
		{
			double pertalt = m_pertalts[pertidx];
			kperturb = pertvals[pertidx];
			for (int idspecies = 0; idspecies < wfs.size(); idspecies++)
			{
				if( addedsolaropticaldepth.size() > 0 )
				{
					wfs[idspecies].WeightAt( pertidx ) += -1.0*(addedsolaropticaldepth[pertidx] ) * (solarsource * effectivecelllength) * exp(-1.0*opticaldepth);
					if( pertalt > cellalt )
						wfs[idspecies].WeightAt( pertidx ) += -1.0*(addedsolaropticaldepth[pertidx] ) * (diffusesource * effectivecelllength) * exp(-1.0*opticaldepth);

				}
				wfs[idspecies].WeightAt( pertidx ) += -1.0*(wfs[idspecies].AddedOpticalDepth( pertidx )) * (source * effectivecelllength) * exp(-1.0*opticaldepth);
				if (kperturb > 0)
				{
					if (wfinfo[idspecies]->IsScatterer())
					{
						double phasefunction = PhaseFunction(*wfinfo[idspecies], pertidx, look.Z());
						double integrationresult = 0;
						if (difftable)
						{
							std::function<double(const nxVector&, const nxVector&)> phasefp = [&](const nxVector& a, const nxVector& b) {
								double cosangle = a.Dot(b);
								return PhaseFunction(*wfinfo[idspecies], pertidx, cosangle);
							};
							integrationresult = difftable->IntegrateScalarIncoming(qMidpoint, phasefp);
							if (integrationresult != integrationresult)
							{
								// This can happen due to some divisions when k=0
								integrationresult = 0.0;
							}
						}
						wfs[idspecies].WeightAt(pertidx) += kperturb * ((solartransmission * phasefunction) + integrationresult)* effectivecelllength * exp(-1.0*opticaldepth);

					}
					wfs[idspecies].WeightAt(pertidx) += source * (-1.0* kperturb) * effectivecelllength_secondorder * exp(-1.0*opticaldepth);
					// TODO: do real optical depth calculation from k on cell bounds
					wfs[idspecies].AddedOpticalDepth(pertidx) += kperturb * ds;
				}
			}
		}
		
	}

	if (ray.Storage()->GroundIsHit())
	{
		double totalgroundsource = 0.0;
		double totalopticaldepth = ray.TotalOpticalDepth();

		HELIODETIC_POINT groundpt;
		ray.Storage()->LocationOfPoint(ray.Storage()->NumQuadraturePoints() - 1, &groundpt);

		SKTRAN_SourceTermQueryObject_Simple qground(groundpt, ray.LookVector());
		double solartransmission = 0.0;
		for (size_t sourceidx = 0; sourceidx < sources.size(); sourceidx++)
		{
			ok = ok && sources[sourceidx]->GroundSourceAtPoint(qground, &tempsource);
			totalgroundsource += tempsource;
			if (sourceidx == 0)
				solarsource = tempsource;
			if (sourceidx == 1)
				diffusesource = tempsource;

			auto solarsourceterm = dynamic_cast<SKTRAN_SolarTransmission_Base*>(sources[sourceidx]);
			if(solarsourceterm != nullptr)
			{
				solarsourceterm->TransmissionAtPoint(groundpt, solartransmission);
			}
		}

		// radiance is solartrans * mu_in * raytrans, but scale this by the ratio of total ground source to direct solar source to put in some 
		// multiple scattering effects
		double mu_in = groundpt.CosSZA();
		wfs[0].brdfwf() += solartransmission * mu_in * exp(-1.0*totalopticaldepth) * (solarsource + diffusesource) / solarsource;

		for (size_t pertidx = 0; pertidx < perturbations.StoreSize(); pertidx++)
		{
			for (int idspecies = 0; idspecies < wfs.size(); idspecies++)
			{
				double addedoptical = wfs[idspecies].AddedOpticalDepth(pertidx);

				wfs[idspecies].WeightAt(pertidx) += -1.0 * addedoptical * totalgroundsource * exp(-1.0*totalopticaldepth);
				if (addedsolaropticaldepth.size() > 0) {
					wfs[idspecies].WeightAt(pertidx) += -1.0 * addedsolaropticaldepth[pertidx] * solarsource * exp(-1.0*totalopticaldepth);
				}
			}

		}
	}

	return ok;
}



bool SKTRAN_HR_WF_Integrator::CalculateWeightingFunctionsPolarized(const SKTRAN_RayOptical_Base &raybase,
                                                                   const SKTRAN_HR_WF_Store &perturbations,
                                                                   const std::vector<SKTRAN_Source_Term *> &sources,
                                                                   const std::vector<std::unique_ptr<SKTRAN_HR_WF_SpeciesInformationBase>> &wfinfo,
                                                                   std::vector<SKTRAN_HR_WF_Ray_Polarized> &wfs,
                                                                   bool directeffectsonly) const {
    bool ok = true;
    double kperturb;
    const SKTRAN_RayOptical_StraightQuadrature_Base&		ray			= dynamic_cast< const SKTRAN_RayOptical_StraightQuadrature_Base&>(raybase);
    std::unique_ptr<SKTRAN_RayOptical_Base> solarray;
    std::unique_ptr<SKTRAN_RayOptical_Base> groundray;

    HELIODETIC_UNITVECTOR sunlook, groundlook;
    sunlook.SetCoords(0,0,1);

    std::vector<double> addedsolaropticaldepth;
    if (directeffectsonly == false)
    {
        addedsolaropticaldepth.resize(perturbations.StoreSize());
    }

    m_rayfactory->CreateRayObject( &solarray );
    m_rayfactory->CreateRayObject( &groundray );

    std::vector<double> pertvals;
    double ds;
    HELIODETIC_POINT midpoint;
    HELIODETIC_UNITVECTOR look;
    SKTRAN_Stokes_NC source, solarsource, diffusesource;
    double singlescatt;
    SKTRAN_Stokes_NC tempsource;
    double k;
    double opticaldepth;

    // TODO: go through the sources and try to find these rather than assuming the order
    const SKTRAN_HR_Diffuse_Source* difftable = nullptr;
    if (sources.size() > 1)
    {
        difftable = dynamic_cast<const SKTRAN_HR_Diffuse_Source*> (sources[1]);
    }
    const SKTRAN_SolarTransmission_Base* solartable = dynamic_cast<const SKTRAN_SolarTransmission_Base*> (sources[0]);

    HELIODETIC_POINT obspt;
    ray.Coordinates()->HelioVectorToHelioPoint( ray.GetObserver(), &obspt );

    SKTRAN_SourceTermQueryObject_StraightPolarized qMidpoint(obspt, ray.LookVector());

    for (int idspecies = 0; idspecies < wfs.size(); idspecies++)
    {
        wfs[idspecies].Allocate(perturbations.StoreSize());
    }

    for( size_t cellidx = 0; cellidx < ray.GetNumCells(); cellidx++ )
    {
        look = ray.Storage()->AverageLookVectorAwayFromObserver( cellidx);
        ds   = ray.Storage()->CellLength( cellidx);
        if( ds < 1E-6 )
            continue;
        ok = ok && ray.Storage()->CellMidPoint(cellidx, &midpoint);

        if (directeffectsonly == false)
        {
            ok = ok && solarray->MoveObserver(midpoint.Vector(), sunlook);
            ok = ok && solarray->TraceRay_NewMethod();

            AddedOpticalDepth(perturbations, addedsolaropticaldepth, *solarray);
        }
        source.SetTo(0);
        solarsource.SetTo(0);
        diffusesource.SetTo(0);
        qMidpoint.UpdateQuery( midpoint, look );
        for( size_t sourceidx = 0; sourceidx < sources.size(); sourceidx++ )
        {
            ok = ok && sources[sourceidx]->SourceTermAtPoint( qMidpoint, &tempsource );
            source += tempsource;
            if( sourceidx == 0)
                solarsource = tempsource;
            if( sourceidx == 1)
                diffusesource = tempsource;
        }

        opticaldepth = ray.OpticalDepthArray().at(cellidx);
        k = GetOpticalProps()->TotalExtinctionPerCM( midpoint ) * 100;
        double kscat = GetOpticalProps()->ScatteringExtinctionPerCM( midpoint ) * 100;
        GetOpticalProps()->GetScatteringCoefficientCM2( midpoint, look.Z(), &singlescatt );
        singlescatt = singlescatt * 100;
        double cellalt = midpoint.Altitude();
        perturbations.ExtinctionPerturbation(midpoint, pertvals);

        // Calculate effictive cell length ( if opticall thin, this is ds, if optically thick this decreases)
        // Also get the solar transmission
        SKTRAN_Stokes_NC solartransmission;
        double effectivecelllength;
        double effectivecelllength_secondorder;
        if (k > 0)
        {
            effectivecelllength = (1 - exp(-1.0*k*ds)) / k;
            effectivecelllength_secondorder = (1 - exp(-ds*k) * (1 + ds*k)) / (k*k);
            solartransmission.SetTo(solarsource.I() * (1/singlescatt), 0, 0, 0);
        }
        else
        {
            // This can happen if there is nothing in the cell
            double temptrans;
            effectivecelllength = ds;
            effectivecelllength_secondorder = ds*ds;
            solartable->TransmissionAtPoint(midpoint, temptrans);
            solartransmission.SetTo(temptrans, 0, 0, 0);
        }


        for (size_t pertidx = 0; pertidx < perturbations.StoreSize(); pertidx++)
        {
            double pertalt = m_pertalts[pertidx];
            kperturb = pertvals[pertidx];
            for (int idspecies = 0; idspecies < wfs.size(); idspecies++)
            {
                if( addedsolaropticaldepth.size() > 0 )
                {
                    wfs[idspecies].WeightAt( pertidx ) += solarsource * (-1.0*(addedsolaropticaldepth[pertidx] ) * (effectivecelllength) * exp(-1.0*opticaldepth));
                    if( pertalt > cellalt )
                        wfs[idspecies].WeightAt( pertidx ) += diffusesource * (-1.0*(addedsolaropticaldepth[pertidx] ) * (effectivecelllength) * exp(-1.0*opticaldepth));

                }
                wfs[idspecies].WeightAt( pertidx ) += source*(-1.0*(wfs[idspecies].AddedOpticalDepth( pertidx )) * (effectivecelllength) * exp(-1.0*opticaldepth));
                if (kperturb > 0)
                {
                    if (wfinfo[idspecies]->IsScatterer())
                    {
                        SKTRAN_ScatMat_MIMSNC phasefunction = PhaseMatrix(*wfinfo[idspecies], pertidx, look.Z());
                        SKTRAN_Stokes_NC integrationresult;
                        integrationresult.SetTo(0.0);
                        if (difftable)
                        {

                            std::function<SKTRAN_ScatMat_MIMSNC(const nxVector&, const nxVector&)> phasefp = [&](const nxVector& a, const nxVector& b) {
                                double cosangle = a.Dot(b);
                                return PhaseMatrix(*wfinfo[idspecies], pertidx, cosangle);
                            };
                            integrationresult = difftable->IntegrateVectorIncoming(qMidpoint, phasefp);
                            if (integrationresult.I() != integrationresult.I())
                            {
                                // This can happen due to some divisions when k=0
                                integrationresult.SetTo(0.0);
                            }

                        }
                        SKTRAN_Stokes_NC temp_stokes = solartransmission;
                        phasefunction.LApplyTo(&temp_stokes);
                        wfs[idspecies].WeightAt(pertidx) += (temp_stokes + integrationresult) * (kperturb * effectivecelllength * exp(-1.0*opticaldepth));

                    }
                    wfs[idspecies].WeightAt(pertidx) += source * (-1.0* kperturb) * effectivecelllength_secondorder * exp(-1.0*opticaldepth);
                    // TODO: do real optical depth calculation from k on cell bounds
                    wfs[idspecies].AddedOpticalDepth(pertidx) += kperturb * ds;
                }
            }
        }

    }
    return ok;
}

double SKTRAN_HR_WF_Integrator::PhaseFunction(const SKTRAN_HR_WF_SpeciesInformationBase& info, size_t pertindex, double cosangle ) const
{
	double result = 0.0;

	std::array<double,2> weight;
	std::array<size_t,2> index;
	info.PhaseAngleGrid().FindBoundingIndices( cosangle, SKTRAN_GridDefBase_V2::OUTOFBOUND_TRUNCATE, &index[0], &weight[0], &index[1], &weight[1] );

	for( size_t idx = 0; idx < 2; idx++ )
	{
		result += info.PhaseFunction().At( index[idx], pertindex ).At(1, 1) * weight[idx];
	}

	return result;
}

SKTRAN_ScatMat_MIMSNC SKTRAN_HR_WF_Integrator::PhaseMatrix(const SKTRAN_HR_WF_SpeciesInformationBase& info, size_t pertindex, double cosangle ) const
{
    SKTRAN_ScatMat_MIMSNC result;
    result.SetTo(double(0));

    std::array<double,2> weight;
    std::array<size_t,2> index;
    info.PhaseAngleGrid().FindBoundingIndices( cosangle, SKTRAN_GridDefBase_V2::OUTOFBOUND_TRUNCATE, &index[0], &weight[0], &index[1], &weight[1] );

    for( size_t idx = 0; idx < 2; idx++ )
    {
        result += info.PhaseFunction().At( index[idx], pertindex ) * weight[idx];
    }

    return result;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_WF_Integrator::AddedOpticalDepth		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_WF_Integrator::AddedOpticalDepth( const SKTRAN_HR_WF_Store& perturbations,
													   std::vector<double>& addedopticaldepth,
													   const SKTRAN_RayOptical_Base& ray
													   ) const
{
	//return true;
	bool ok = true;
	//double kperturb;
	HELIODETIC_POINT midpoint;
	double ds;
	if( addedopticaldepth.size() != perturbations.StoreSize() )
		addedopticaldepth.resize( perturbations.StoreSize() );

	std::vector<double> pertvals;
	std::fill( std::begin(addedopticaldepth), std::end(addedopticaldepth), 0.0 );
	for( size_t cellidx = 0; cellidx < ray.GetNumCells(); cellidx++ )
	{
		ok = ok && ray.Storage()->CellMidPoint( cellidx, &midpoint );
		ds = ray.Storage()->CellLength( cellidx);
		if( ds < 1E-6 )
		{
			continue;
		}
		perturbations.AddExtinctionPerturbation(midpoint, addedopticaldepth, ds);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_WF_Integrator::CreateSolarRayFactory		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_WF_Integrator::CreateSolarRayFactory( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
													 const SKTRAN_HR_WF_Store& perturbation, 
	                                                 double toaHeight )
{
	bool ok = true;

	SKTRAN_RayFactory<  SKTRAN_RayOptical_Straight,
						SKTRAN_RayTracer_Straight_Generic,
						SKTRAN_RayStorage_Straight_HR >*			factory = new SKTRAN_RayFactory<  SKTRAN_RayOptical_Straight,
																								  SKTRAN_RayTracer_Straight_Generic,
																								  SKTRAN_RayStorage_Straight_HR> (coords);

	
	factory->RayTracer()->SetEarthRadius( coords->AltitudeToRadius( 0.0 ) );
	factory->RayTracer()->SetUpperAtmoRadius( coords->AltitudeToRadius( 0.0 ) + toaHeight );

	perturbation.AddGeometryToRayTracer(coords, *factory->RayTracer());

	m_rayfactory.reset(factory);
//	factory->SetThisFactory( m_rayfactory );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_WF_Integrator::SetSolarRayFactory		 2014- 12- 4*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_WF_Integrator::SetSolarRayFactory( SKTRAN_RayFactory_Base*	rayfactory )
{
	m_rayfactory.reset(rayfactory);
	//rayfactory->SetThisFactory(m_rayfactory);
	return true;
}

bool SKTRAN_HR_WF_Integrator::CreatePhaseFunctionTableParticleSize(skClimatology& neutral, SKTRAN_AtmosphericOpticalState_V21& opticalstate, skOpticalProperties& optprop, double wavel, SKTRAN_HR_WF_Store& perturbations, const SKTRAN_CoordinateTransform_V2& coords,
	double moderadiuschange, double modewidthchange)
{
	bool ok = true;

	auto aero_optprop = dynamic_cast<skOpticalProperties_AerosolProfile*>(&optprop);

	if (!aero_optprop)
	{
		nxLog::Record(NXLOG_WARNING, "Error calculating particle size wf");
	}

	skClimatology* aero_clim;

	opticalstate.GetSpeciesClimatology(SKCLIMATOLOGY_AEROSOL_CM3, &aero_clim);

	skOpticalProperties_AerosolProfileH2SO4 mie_aerosol;

	std::vector<double> particlesize_heights;
	particlesize_heights.push_back(0);
	particlesize_heights.push_back(100000);

	std::vector<double> moderadius;
	std::vector<double> modewidth;

	moderadius.push_back(0.08);
	moderadius.push_back(0.08);

	modewidth.push_back(1.6);
	modewidth.push_back(1.6);

	double pert_moderadius;
	double pert_modewidth;
	
	m_phaseanglegrid.Configure(0.5);
	m_phasefunction.SetSize(m_phaseanglegrid.NumAngles(), perturbations.StoreSize());

	for (size_t pertidx = 0; pertidx < perturbations.StoreSize(); pertidx++)
	{
		HELIODETIC_POINT loc = perturbations.RawAccess(pertidx)->PerturbationLocation(coords);
		GEODETIC_INSTANT geoloc = coords.PointToGeodetic(loc, coords.ReferencePointMJD());
		bool haschanged = true;
		optprop.SetAtmosphericState(&neutral);
		optprop.SetLocation(geoloc, &haschanged);

		double aero_n;
		aero_clim->GetParameter(SKCLIMATOLOGY_AEROSOL_CM3, geoloc, &aero_n, true);

		haschanged = true;
		mie_aerosol.SetAtmosphericState(&neutral);
		mie_aerosol.SetLocation(geoloc, &haschanged);

		aero_optprop->GetDistributionParameter(SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS, geoloc, &pert_moderadius);
		aero_optprop->GetDistributionParameter(SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH, geoloc, &pert_modewidth);

		skRTPhaseMatrix pmatrix;
		skRTPhaseMatrix pmatrix_pert;
		double absxs, extxs, scattxs;
		optprop.CalculateCrossSections(1E7 / wavel, &absxs, &extxs, &scattxs);

		double dmr = pert_moderadius * moderadiuschange;
		double dmw = pert_modewidth * modewidthchange;

		moderadius[0] = pert_moderadius + dmr;
		moderadius[1] = pert_moderadius + dmr;

		modewidth[0] = pert_modewidth + dmw;
		modewidth[1] = pert_modewidth + dmw;

		mie_aerosol.SetLogNormalProfileClimatology(&particlesize_heights[0] , &moderadius[0], &modewidth[0], 2);
		double absxs_pert, extxs_pert, scattxs_pert;
		mie_aerosol.CalculateCrossSections(1E7 / wavel, &absxs_pert, &extxs_pert, &scattxs_pert);
		
		// Derivative of extinction wrt to particle size change
		// Will be multiplied by extxs*100 LATER
		double ext_deriv = (extxs_pert - extxs) / (dmr + dmw) / extxs * aero_n;
		perturbations.RawAccess(pertidx)->SetPertVal(ext_deriv);

		for (size_t angleidx = 0; angleidx < m_phaseanglegrid.NumAngles(); angleidx++)
		{
			optprop.CalculatePhaseMatrix(1E7 / wavel, m_phaseanglegrid.At(angleidx), &pmatrix);
			mie_aerosol.CalculatePhaseMatrix(1E7 / wavel, m_phaseanglegrid.At(angleidx), &pmatrix_pert);

			// Derivative of phase function wrt to particle size change
			double phase_deriv = (pmatrix_pert.At(1, 1) - pmatrix.At(1, 1)) / (dmr + dmw);

			phase_deriv *= aero_n;
			phase_deriv /= ext_deriv;

			m_phasefunction.At(angleidx, pertidx) = (pmatrix.At(1, 1) + phase_deriv) * (scattxs / extxs) / (4.0*nxmath::Pi); // scale by single scattering albedo
		}
	}
	m_isscatterer = true;
	return ok;
}