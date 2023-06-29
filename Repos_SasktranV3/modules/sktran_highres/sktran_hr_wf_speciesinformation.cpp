#include "include/sktran_hr_internals.h"

SKTRAN_HR_WF_SpeciesInformationBase::SKTRAN_HR_WF_SpeciesInformationBase()
{
	m_isscatterer = false;
	m_phaseanglegrid = nullptr;
	m_phasefunction = nullptr;
	m_xs = nullptr;
}

void SKTRAN_HR_WF_SpeciesInformationBase::Initialize(std::unique_ptr<SKTRAN_GridDefScatterAngle_V21>&& phaseanglegrid,
	std::unique_ptr<nx2dArray<SKTRAN_ScatMat_MIMSNC>>&& phasefunction,
	bool isscatterer,
	std::unique_ptr<std::vector<double>>&& xs)
{
	m_phaseanglegrid = std::move(phaseanglegrid);
	m_phasefunction = std::move(phasefunction);
	m_isscatterer = isscatterer;
	m_xs = std::move(xs);
}


SKTRAN_HR_WF_SpeciesInformationStandard::SKTRAN_HR_WF_SpeciesInformationStandard(skClimatology& backgroundstate, skOpticalProperties& opticalproperty, const SKTRAN_HR_WF_Store& wfstore, const SKTRAN_CoordinateTransform_V2& coordinates,
	double wavelength, const SKTRAN_HR_Specs_Internal_OpticalPropertiesTable& opticalspecs)
{
	bool isscatterer = opticalproperty.IsScatterer();
	double f;

    skRTPhaseMatrix temp;

	auto xs = std::unique_ptr<std::vector<double>>(new std::vector<double>);
	std::unique_ptr<SKTRAN_GridDefScatterAngle_V21> phaseanglegrid = nullptr;
	std::unique_ptr<nx2dArray<SKTRAN_ScatMat_MIMSNC>> phasefunction = nullptr;

	xs->resize(wfstore.StoreSize());

	for (int idx = 0; idx < wfstore.StoreSize(); idx++)
	{
		HELIODETIC_POINT pertpoint = wfstore.RawAccess(idx)->PerturbationLocation(coordinates);
		GEODETIC_INSTANT pertgeo = coordinates.PointToGeodetic(pertpoint, coordinates.ReferencePointMJD());
		bool crosssectionchanged = true;
		opticalproperty.SetLocation(pertgeo, &crosssectionchanged);

		double absxs, scattxs, extxs;
		opticalproperty.CalculateCrossSections(1E7 / wavelength, &absxs, &extxs, &scattxs);

		if (opticalproperty.IsDeltaFunctionForwardScatter())
		{
			f = opticalproperty.DeltaFunctionForwardScatterFraction();
		}
		else
		{
			f = 0.0;
		}

		xs->at(idx) = (1 - (scattxs / extxs) * f) * extxs;
	}

	if (isscatterer)
	{
		phaseanglegrid = std::unique_ptr<SKTRAN_GridDefScatterAngle_V21>(new SKTRAN_GridDefScatterAngle_V21);
		phasefunction = std::unique_ptr<nx2dArray<SKTRAN_ScatMat_MIMSNC>>(new nx2dArray<SKTRAN_ScatMat_MIMSNC>);

		opticalspecs.MakeScatterAngleGrid(*phaseanglegrid);
		phasefunction->SetSize(phaseanglegrid->NumAngles(), wfstore.StoreSize());

		for (size_t pertidx = 0; pertidx < wfstore.StoreSize(); pertidx++)
		{
			HELIODETIC_POINT loc = wfstore.RawAccess(pertidx)->PerturbationLocation(coordinates);
			GEODETIC_INSTANT geoloc = coordinates.PointToGeodetic(loc, coordinates.ReferencePointMJD());
			bool haschanged = true;
			opticalproperty.SetAtmosphericState(&backgroundstate);
			opticalproperty.SetLocation(geoloc, &haschanged);
			double p11, norm;

			double absxs, extxs, scattxs;
			opticalproperty.CalculateCrossSections(1E7 / wavelength, &absxs, &extxs, &scattxs);
			if (opticalproperty.IsDeltaFunctionForwardScatter())
			{
				f = opticalproperty.DeltaFunctionForwardScatterFraction();
				norm = 4.0*nxmath::Pi * (1 - f);
			}
			else
			{
				f = 0.0;
				norm = 4.0*nxmath::Pi;
			}

			double adjustedkscat = (1 - f) * scattxs;
			double adjustedkext = (1 - scattxs / extxs * f) * extxs;
			double ssa = adjustedkscat / adjustedkext;

			for (size_t angleidx = 0; angleidx < phaseanglegrid->NumAngles(); angleidx++)
			{
				// opticalproperty.CalculateP11(1E7 / wavelength, std::make_pair(phaseanglegrid->At(angleidx), angleidx), temp);
                opticalproperty.CalculatePhaseMatrix(1E7 / wavelength, phaseanglegrid->At(angleidx), &temp);
                temp *= (ssa) / (norm); // scale by single scattering albedo
				phasefunction->At(angleidx, pertidx) = SKTRAN_ScatMat_MIMSNC(temp);
			}
		}
	}

	this->Initialize(std::move(phaseanglegrid), std::move(phasefunction), isscatterer, std::move(xs));
}


SKTRAN_HR_WF_SpeciesInformationAerosolLogNormal::SKTRAN_HR_WF_SpeciesInformationAerosolLogNormal(skClimatology& backgroundstate, skOpticalProperties& optprop, const SKTRAN_HR_WF_Store& perturbations, const SKTRAN_CoordinateTransform_V2& coordinates,
	double wavel, SKTRAN_AtmosphericOpticalState_V21& opticalstate, skClimatology* aero_clim, const CLIMATOLOGY_HANDLE& aerohandle, double moderadiuschange, double modewidthchange, const SKTRAN_HR_Specs_Internal_OpticalPropertiesTable& opticalspecs)
{
	bool isscatterer = true;

	auto xs = std::unique_ptr<std::vector<double>>(new std::vector<double>);
	std::unique_ptr<SKTRAN_GridDefScatterAngle_V21> phaseanglegrid = nullptr;
	std::unique_ptr<nx2dArray<SKTRAN_ScatMat_MIMSNC>> phasefunction = nullptr;

	phaseanglegrid = std::unique_ptr<SKTRAN_GridDefScatterAngle_V21>(new SKTRAN_GridDefScatterAngle_V21);
	phasefunction = std::unique_ptr<nx2dArray<SKTRAN_ScatMat_MIMSNC>>(new nx2dArray<SKTRAN_ScatMat_MIMSNC>);

	xs->resize(perturbations.StoreSize());

	bool ok = true;

	auto aero_optprop = dynamic_cast<skOpticalProperties_AerosolProfile*>(&optprop);

	if (!aero_optprop)
	{
		nxLog::Record(NXLOG_WARNING, "Error calculating particle size wf");
	}

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

	opticalspecs.MakeScatterAngleGrid(*phaseanglegrid);
	phasefunction->SetSize(phaseanglegrid->NumAngles(), perturbations.StoreSize());

	for (size_t pertidx = 0; pertidx < perturbations.StoreSize(); pertidx++)
	{
		HELIODETIC_POINT loc = perturbations.RawAccess(pertidx)->PerturbationLocation(coordinates);
		GEODETIC_INSTANT geoloc = coordinates.PointToGeodetic(loc, coordinates.ReferencePointMJD());
		bool haschanged = true;
		optprop.SetAtmosphericState(&backgroundstate);
		optprop.SetLocation(geoloc, &haschanged);

		double aero_n;
		aero_clim->GetParameter(aerohandle, geoloc, &aero_n, true);

		haschanged = true;
		mie_aerosol.SetAtmosphericState(&backgroundstate);
		mie_aerosol.SetLocation(geoloc, &haschanged);

		aero_optprop->GetDistributionParameter(SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS, geoloc, &pert_moderadius);
		aero_optprop->GetDistributionParameter(SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH, geoloc, &pert_modewidth);

		skRTPhaseMatrix pmatrix;
		skRTPhaseMatrix pmatrix_pert;
        skRTPhaseMatrix pmatrix_temp;
		double absxs, extxs, scattxs;
		optprop.CalculateCrossSections(1E7 / wavel, &absxs, &extxs, &scattxs);

		double dmr = pert_moderadius * moderadiuschange;
		double dmw = pert_modewidth * modewidthchange;

		moderadius[0] = pert_moderadius + dmr;
		moderadius[1] = pert_moderadius + dmr;

		modewidth[0] = pert_modewidth + dmw;
		modewidth[1] = pert_modewidth + dmw;

		mie_aerosol.SetLogNormalProfileClimatology(&particlesize_heights[0], &moderadius[0], &modewidth[0], 2);
		double absxs_pert, extxs_pert, scattxs_pert;
		haschanged = true;
		mie_aerosol.SetLocation(geoloc, &haschanged);
		mie_aerosol.CalculateCrossSections(1E7 / wavel, &absxs_pert, &extxs_pert, &scattxs_pert);

		// Derivative of extinction wrt to particle size change
		double ext_deriv = (extxs_pert - extxs) / (dmr + dmw) * aero_n;
		xs->at(pertidx) = ext_deriv;

		for (size_t angleidx = 0; angleidx < phaseanglegrid->NumAngles(); angleidx++)
		{
			optprop.CalculatePhaseMatrix(1E7 / wavel, phaseanglegrid->At(angleidx), &pmatrix);
			mie_aerosol.CalculatePhaseMatrix(1E7 / wavel, phaseanglegrid->At(angleidx), &pmatrix_pert);

			// Derivative of phase function wrt to particle size change

            pmatrix_temp = (pmatrix_pert - pmatrix) * (aero_n * extxs/(dmr + dmw)/ext_deriv);

            phasefunction->At(angleidx, pertidx) = SKTRAN_ScatMat_MIMSNC((pmatrix + pmatrix_temp) * (scattxs / extxs) * (1/(4.0*nxmath::Pi)));
		}
	}
	this->Initialize(std::move(phaseanglegrid), std::move(phasefunction), isscatterer, std::move(xs));
}