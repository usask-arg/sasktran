/**
 * SASKTRAN TIR Integrator
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_Integrator::SKTRAN_TIR_Integrator
 * 2018-09-13
 */
SKTRAN_TIR_Integrator::SKTRAN_TIR_Integrator()
{
	m_opticalprops = nullptr;
	m_maxopticaldepth = 100000.0;
	m_minextinctionratio = 0.95;
	m_optinttype = OpticalPropertiesIntegratorTypeTIR::adaptive;
	m_srcorder = SourceTermOrderTIR::order0;
	m_extinctiontype = LayerExtinctionTypeTIR::linearwithheight;
	m_wfunit = SpeciesWFUnitTIR::numberdensity;
}

/**
 * SKTRAN_TIR_Integrator::~SKTRAN_TIR_Integrator
 * 2018-09-13
 */
SKTRAN_TIR_Integrator::~SKTRAN_TIR_Integrator()
{
	ReleaseResources();
}

/**
 * SKTRAN_TIR_Integrator::ReleaseResources
 * 2018-09-13
 */
void SKTRAN_TIR_Integrator::ReleaseResources() {
	if (nullptr != m_opticalprops) m_opticalprops->Release(); m_opticalprops = nullptr;
}

/**
 * SKTRAN_TIR_Integrator::SetOpticalProps
 * 2018-09-13
 *
 * Sets the optical properties table for this integrator to use in calculations. Saves a reference to the provided table,
 * so any changes made to the optical properties table after calling this function will impact the integration calculation.
 *
 * @param[in] optprop Pointer to the optical properties table to use
 *
 * @post If the optical properties were previously set, that table is released before setting the new one.
 */
bool SKTRAN_TIR_Integrator::SetOpticalProps(
	const SKTRAN_TIR_TableOpticalProperties* optprop)
{
	bool ok = true;

	ok = ok && optprop != nullptr;

	if (ok)
	{
		optprop->AddRef();
		ReleaseResources();
		m_opticalprops = optprop;
	}
	else
	{
		nxLog::Record(NXLOG_WARNING, "Error, optical properties table is nullptr in SKTRAN_TIR_Integrator::SetOpticalProps");
		ok = false;
	}

	return ok;
}

/**
 * SKTRAN_TIR_Integrator::OpticalDepthOfCell
 * 2018-09-13
 *
 * Calculates the optical depth across a single cell (between two quadrature points). Calls the appropriate method
 * based on configuration settings.
 *
 * @param[in,out] ray Optical ray which has been traced, after calling this function the ray's storage will be
 *                    modified to contain the extinction values at the start and end points of the cell. Note that
 *                    despite being const the ray's storage can be modified because the members of the storage class
 *                    are mutable.
 * @param[in] cellidx Index of this cell in the ray's storage
 * @param[in] wavelidx Index of this wavelength in the optical properties table
 */
double SKTRAN_TIR_Integrator::OpticalDepthOfCell(
	const SKTRAN_RayOptical_Base* ray,
	size_t cellidx,
	size_t wavelidx) const
{
	double			opticaldepth;
	bool			ok = true;

	NXASSERT(cellidx < (ray->GetNumQuadraturePoints() - 1));

	if (m_extinctiontype == LayerExtinctionTypeTIR::linearwithheight)
	{
		opticaldepth = OpticalDepthOfSegment_LinearWithHeight(cellidx, wavelidx, ray);
	}
	else  // m_extinctiontype == LayerExtinctionTypeTIR::constant
	{
		opticaldepth = OpticalDepthOfSegment_Constant(cellidx, wavelidx, ray);
	}

	NXASSERT(((opticaldepth >= 0.0) && (opticaldepth < 1.0E6)));
	return opticaldepth;
}

/**
 * SKTRAN_TIR_Integrator::OpticalDepthOfSegment_Constant
 * 2018-09-13
 *
 * Calculates the optical depth across a single cell, assuming a constant extinction within the cell. The cell extinction
 * is computed as the average of the extinction values at the start and end points of the cell.
 *
 * @param[in,out] ray Optical ray which has been traced, after calling this function the ray's storage will be
 *                    modified to contain the extinction values at the start and end points of the cell
 * @param[in] cellidx Index of this cell in the ray's storage
 * @param[in] wavelidx Index of this wavelength in the optical properties table
 */
double SKTRAN_TIR_Integrator::OpticalDepthOfSegment_Constant(
	size_t cellidx,
	size_t wavelidx,
	const SKTRAN_RayOptical_Base* ray) const
{
	double opticaldepth = -9999.0;

	bool						ok = true;
	double						sigma0 = -9999.0;
	double						sigma1 = -9999.0;
	double						ds;

	ok = ok && m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1(ray->Storage(), cellidx, &sigma0, &sigma1, wavelidx);

	NXASSERT(((sigma0 >= 0) && (sigma0 < 1.0E4)));
	NXASSERT(((sigma1 >= 0) && (sigma1 < 1.0E4)));

	if (ok)
	{
		// calculate optical depth along the cell pathlength
		ds = ray->Storage()->CellLength(cellidx);
		opticaldepth = 100.0 * ds * (sigma0 + sigma1) / 2.0;
	}
	ok = ok && (opticaldepth >= 0.0);
	if (!ok)
	{
		if (1e-7 < fabs(opticaldepth))
		{
			// Only print warning if the error is really of a magnitude that we care about
			nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_Integrator::OpticalDepthOfSegment_Constant, Error looking up optical depth of a segment, cellidx =%d, opticaldepth =%18.8e, sigma0= %18.8e, sigma1=%18.8e", (int)cellidx, (double)opticaldepth, (double)sigma0, (double)sigma1);
		}
		opticaldepth = 0;
	}
	// Unlike for OpticalDepthOfSegment_LinearWithHeight, we do not need to multiply by the cell curvature because the
	// path length (ds) obtained from the ray storage will be the length of the curved path if refraction is considered
	return opticaldepth;
}

/**
 * SKTRAN_TIR_Integrator::OpticalDepthOfSegment_LinearWithHeight
 * 2018-09-13
 *
 * Calculates the optical depth across a single cell, where the extinction is taken to vary across the cell,
 * linearly with height. The SKTRAN_OpticalDepthCalculator_LinearWithHeight class from the SASKTRAN common
 * code is utilized for this calculation.
 *
 * @param[in,out] ray Optical ray which has been traced, after calling this function the ray's storage will be
 *                    modified to contain the extinction values at the start and end points of the cell
 * @param[in] cellidx Index of this cell in the ray's storage
 * @param[in] wavelidx Index of this wavelength in the optical properties table
 */
double SKTRAN_TIR_Integrator::OpticalDepthOfSegment_LinearWithHeight(
	size_t cellidx,
	size_t wavelidx,
	const SKTRAN_RayOptical_Base* ray) const
{
	double opticaldepth = -9999.0;

	bool						ok = true;
	double						r0;
	double						r1;
	double						rt;
	double						t0;
	double						t1;
	double						sigma0 = -9999.0;
	double						sigma1 = -9999.0;
	SKTRAN_OpticalDepthCalculator_LinearWithHeight	odcalculator;

	ok = ok && m_opticalprops->GetEffectiveExtinctionPerCMWithHeight1(ray->Storage(), cellidx, &sigma0, &sigma1, wavelidx);

	NXASSERT(((sigma0 >= 0) && (sigma0 < 1.0E4)));
	NXASSERT(((sigma1 >= 0) && (sigma1 < 1.0E4)));

	r0 = ray->Storage()->RadiusOfPoint(cellidx);
	r1 = ray->Storage()->RadiusOfPoint(cellidx + 1);
	t0 = ray->Storage()->DistanceOfPointFromCellTangentPoint(cellidx, cellidx);
	t1 = ray->Storage()->DistanceOfPointFromCellTangentPoint(cellidx + 1, cellidx);
	rt = ray->Storage()->RadiusOfCellTangentPoint(cellidx);

	ok = ok && odcalculator.ConfigureQuadratureCoefficients(r0, r1, t0, t1, rt);
	if (ok)
	{
		// calculate optical depth for a straight ray
		opticaldepth = odcalculator.OpticalDepthFromStartToEnd(sigma0, sigma1);
	}
	ok = ok && (opticaldepth >= 0.0);
	if (!ok)
	{
		if (1e-7 < fabs(opticaldepth) || fabs(r0 - r1)>1.0)
		{
			// Only print warning if the error is really of a magnitude that wew care about
			nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_Integrator::OpticalDepthOfSegment_LinearWithHeight, Error looking up optical depth of a segment, cellidx =%d, opticaldepth =%18.8e, sigma0= %18.8e, sigma1=%18.8e,  r0=%18.8e, r1=%18.8e, t0=%18.8e, t1=%18.8e", (int)cellidx, (double)opticaldepth, (double)sigma0, (double)sigma1, (double)r0, (double)r1, (double)t0, (double)t1);
		}
		opticaldepth = 0;
	}
	// correct optical depth for the curved ray path length
	return opticaldepth * ray->Storage()->CellCurvature(cellidx);
}

/**
 * SKTRAN_TIR_Integrator::GetSecondOrderTerms
 * 2019-02-04
 *
 * Computes terms required if the source term in a cell is fitted to a second
 * order function.
 * 
 * @param[in] sourcestart The source term evaluated at the beginning of the
 *            cell, which is the point furthest from the observer.
 * @param[in] sourcemid The source term evaluated at the midpoint of the cell.
 * @param[in] sourceend the source term evaluated at the end of the cell, which
 *            is the point nearest the observer.
 * @param[in] ds The path length of the current cell.
 * @param[in] opticaldepthcell The optical depth across the current cell.
 * @param[in] transmissioncell The transmission across the current cell.
 * @param[in] rad The radiance at the start of the cell.
 * @param[out] sourcecell The order 2 source term value computed for this cell.
 * @param[out] derirad A term in the radiance calculation which accounts for the
 *             difference between the incoming radiance and the source term in
 *             the cell.
 * @param[out] drad_dopt The derivative of radiance with respect to optical
 *             depth, used to calculate analytic weighting functions.
 */
bool SKTRAN_TIR_Integrator::GetSecondOrderTerms(
	const double& sourcestart,
	const double& sourcemid,
	const double& sourceend,
	double ds,
	double opticaldepthcell,
	double transmissioncell,
	double rad,
	double& sourcecell,
	double& derirad,
	double& drad_dopt)
{
	bool ok = true;
	double a, b, c;

	ok = ok && GetQuadraticCoeff(sourcestart, sourcemid, sourceend, ds, a, b, c);

	double k = opticaldepthcell / ds;

	// the exact formula causes precision errors, so use an approximation
	// if it is valid
	if (ds * k < 0.01)
	{
		double ds1k1(ds * k);
		double ds2k1(ds * ds1k1);
		double ds3k2(ds2k1 * ds1k1);
		double ds4k3(ds3k2 * ds1k1);
		double ds5k4(ds4k3 * ds1k1);
		double ds6k5(ds5k4 * ds1k1);
		double ds7k6(ds6k5 * ds1k1);
		derirad = rad * transmissioncell + a * (-transmissioncell) +
			b * (-ds + ds2k1 / 2.0 - ds3k2 / 6.0 + ds4k3 / 24.0 - ds5k4 / 120.0 + ds6k5 / 720.0 - ds7k6 / 5040.0) +
			c * (-(ds * ds) + ds * ds2k1 / 3.0 - ds * ds3k2 / 12.0 + ds * ds4k3 / 60.0 - ds * ds5k4 / 360.0 + ds * ds6k5 / 2520.0 - ds * ds7k6 / 20160.0);
		sourcecell = a + b * ds + c * ds * ds;
		// The below approximation still results in precision errors, so use order 0 form for now
		//drad_dopt = -rad * transmissioncell + a * transmissioncell +
		//	b * ((2.0 + ds * k) / (2.0 * k) - ds2k1 / 3.0 + ds3k2 / 8.0 - ds4k3 / 30.0 + ds5k4 / 144.0 - ds6k5 / 840.0 + ds7k6 / 5760.0) +
		//	2 * c * ((-6.0 + 6.0 * ds * k + ds * ds * k * k) / (6.0 * k * k) - ds * ds2k1 / 12.0 + ds * ds3k2 / 40.0 - ds * ds4k3 / 180.0 + ds * ds5k4 / 1008.0 - ds * ds6k5 / 6720.0 + ds * ds7k6 / 51840.0);
		drad_dopt = transmissioncell * (sourcemid - rad);
	}
	else
	{
		derirad = rad * transmissioncell + a * (-transmissioncell) + b * ((transmissioncell - 1) / k) + c * (2 / (k * k) * (1 - transmissioncell) - 2 * ds / k);
		sourcecell = a + b * ds + c * ds * ds;
		// Results in precision errors, so use order 0 term for now
		//drad_dopt = -rad * transmissioncell +
		//	a * transmissioncell +
		//	b * ((1 - transmissioncell) * (1 / k + 1 / (k * k * ds))) +
		//	2 * c * (ds / k + 1 / (k * k) - (1 - transmissioncell) * (2 / (k * k * k * ds) + 1 / (k * k)));
		drad_dopt = transmissioncell * (sourcemid - rad);	// using order 0 weighting function because order 2 blows up
	}

	if (sourcecell < 0) sourcecell = 0.0;
	if (std::isnan(sourcecell)) sourcecell = 0.0;

	return ok;
}

/**
 * SKTRAN_TIR_Integrator::GetQuadraticCoeff
 * 2018-09-24
 *
 * Computes the quadratic coefficients, b and c, for a 2nd order representation of the source term within a cell,
 * where the source term, J(s), is a quadratic function of distance given by
 *   J(s) = a + b * s + c * s^2
 * where s is the distance measured from the start point of the cell.
 *
 * @param[in] sourcestart Source term at start point of cell
 * @param[in] sourcemid Source term at midpoint of cell
 * @param[in] sourceend Source term at endpoint of cell
 * @param[in] ds radiative transfer path length through the cell
 * @param[out] a constant term in J(s)
 * @param[out] b coefficient of the linear part of J(s)
 * @param[out] c coefficient of the 2nd order part of J(s)
 */
bool SKTRAN_TIR_Integrator::GetQuadraticCoeff(
	const double& sourcestart,
	const double& sourcemid,
	const double& sourceend,
	double ds,
	double& a,
	double& b,
	double& c)
{
	bool ok = true;

	a = sourcestart;
	b = (sourcestart * 3.0 + sourcemid * (-4.0) + sourceend) * (-1.0 / ds);
	c = (sourcestart + sourcemid * (-2.0) + sourceend) * (2.0 / (ds * ds));

	return ok;
}

/**
 * SKTRAN_TIR_Integrator::IntegrateRay
 * 2018-09-24
 *
 * Performs the radiative transfer integration along a line of sight for a given wavelength.
 *
 * @param[in,out] ray The optical ray to perform the integration along. The ray must be traced before calling this
 *                    this function. This function will modify the ray's internal storage. In particular, the extinction
 *                    and weighting function storage will be altered.
 * @param[out] radiance The calculated radiance which would be observed from the start point of this ray, looking along
 *                      the line of sight to the end point of this ray
 * @param[in] wf_handles The climatology handles of molecular species whose weighting functions are being computed.
 *                       The cross sections of these species must be included in the optical properties table which was
 *                       assigned to this integrator use SetOpticalProps
 * @param[in] wavelidx The wavelength index
 */
bool SKTRAN_TIR_Integrator::IntegrateRay(
	SKTRAN_RayOptical_Base* ray,
	double& radiance,
	const std::vector<CLIMATOLOGY_HANDLE>& wf_handles,
	size_t wavelidx) const
{
	bool ok = true;
	bool calcwf;
	SKTRAN_RayStorage_Base* raystorage = dynamic_cast<SKTRAN_RayStorage_Base*>(ray->StorageVar());
	std::vector<double>* odstorage = ray->OpticalDepthArrayVar();
	std::map<CLIMATOLOGY_HANDLE, std::vector<double>> wfcellstorage;  // storage for WF contributions from each cell before they are interpolated to quadrature points
	double sourcestart, sourcemid, sourceend, sourcecell, derirad, rad, opticaldepthcell, transmissioncell, ds, drad_dopt, opticaldepthtoobserver;

	HELIODETIC_POINT startpoint, midpoint, endpoint;
	HELIODETIC_UNITVECTOR look;

	HELIODETIC_POINT obspt;
	ray->Coordinates()->HelioVectorToHelioPoint(ray->GetObserver(), &obspt);

	// are weighting functions being calculated?
	calcwf = wf_handles.size() > 0;

	for (auto const& species : wf_handles)
	{
		raystorage->AddWFSpecies(species);
		wfcellstorage.emplace(species, std::vector<double>());
		wfcellstorage[species].reserve(raystorage->NumQuadraturePoints() * 3);
		wfcellstorage[species].resize(0);
	}

	// Check end LOS
	rad = 0.0;
	if (raystorage->GroundIsHit())
	{
		// LOS ends at ground
		skRTStokesVector::SetToZero(rad);
		HELIODETIC_POINT groundpoint;
		raystorage->LocationOfPoint(raystorage->NumQuadraturePoints() - 1, &groundpoint);
		look = raystorage->AverageLookVectorAwayFromObserver(raystorage->NumCells() - 1);
		rad = m_opticalprops->GroundSourceAtPointAndWavel(groundpoint, wavelidx);
	}

	// allocate optical depth storage
	if (calcwf)
	{
		odstorage->reserve(raystorage->NumQuadraturePoints() * 3);
		odstorage->resize(0);
	}

	ok = ok && raystorage->LocationOfPoint(raystorage->NumQuadraturePoints() - 1, &startpoint);
	look = raystorage->AverageLookVectorAwayFromObserver(raystorage->NumCells() - 1);
	sourcestart = m_opticalprops->SourceTermAtPointAndWavel(startpoint, wavelidx);
	// Begin with layer furthest from observer
	for (size_t startquadpt = raystorage->NumQuadraturePoints() - 1; startquadpt > 0; startquadpt--)
	{
		size_t cellidx = startquadpt - 1;

		opticaldepthcell = OpticalDepthOfCell(ray, cellidx, wavelidx);

		double kstart = raystorage->ExtinctionAtCellStart(startquadpt);
		double kend = raystorage->ExtinctionAtCellStart(cellidx);

		ds = raystorage->CellLength(cellidx);

		// Adaptive integration
		if (m_optinttype == OpticalPropertiesIntegratorTypeTIR::adaptive)
		{
			while (opticaldepthcell > m_maxopticaldepth && std::min(kstart, kend) / std::max(kstart, kend) < m_minextinctionratio)
			{
				raystorage->SplitCell(cellidx);

				cellidx++;
				startquadpt++;

				opticaldepthcell = OpticalDepthOfCell(ray, cellidx, wavelidx);

				kend = raystorage->ExtinctionAtCellStart(cellidx);

				ds = raystorage->CellLength(cellidx);
			}
		}

		if (calcwf)
		{
			odstorage->insert(odstorage->begin(), opticaldepthcell);
		}
		transmissioncell = exp(-1.0 * opticaldepthcell);

		look = raystorage->AverageLookVectorAwayFromObserver(cellidx);

		// take average source term
		ok = ok && raystorage->LocationOfPoint(startquadpt, &startpoint);
		ok = ok && raystorage->LocationOfPoint(cellidx, &endpoint);
		sourceend = m_opticalprops->SourceTermAtPointAndWavel(endpoint, wavelidx);

		// Source term in middle of cell
		ok = ok && raystorage->CellMidPoint(cellidx, &midpoint);
		sourcemid = m_opticalprops->SourceTermAtPointAndWavel(midpoint, wavelidx);

		// Is the source term 2nd or 0th order
		if (m_srcorder == SourceTermOrderTIR::order2)
		{
			// Source term at end of cell
			ok = ok && raystorage->LocationOfPoint(cellidx, &endpoint);
			ok = ok && GetSecondOrderTerms(sourcestart, sourcemid, sourceend, ds, opticaldepthcell, transmissioncell, rad, sourcecell, derirad, drad_dopt);
			sourcestart = sourceend;
		}
		else // m_srcorder == SourceTermOrderTIR::order0
		{
			sourcecell = (sourcestart + sourceend) / 2.0;
			derirad = transmissioncell * (rad - sourcecell);
			drad_dopt = -derirad;
			sourcestart = sourceend;
		}

		// computes dI/dn and/or dI/dT
		for (auto const& species : wf_handles)
		{
			if (species == SKCLIMATOLOGY_TEMPERATURE_K)
			{
				double dkdTcell = (m_opticalprops->AbsorptionTemperatureDerivativeAtPointAndWavel(startpoint, wavelidx)
								   + m_opticalprops->AbsorptionTemperatureDerivativeAtPointAndWavel(endpoint, wavelidx)) / 2.0;
				double dBdTcell = (m_opticalprops->PlanckFunctionTemperatureDerivativeAtPointAndWavel(startpoint, wavelidx)
								   + m_opticalprops->PlanckFunctionTemperatureDerivativeAtPointAndWavel(endpoint, wavelidx)) / 2.0;
				wfcellstorage[species].insert(wfcellstorage[species].begin(), (sourcecell - rad) * transmissioncell * ds * 100.0 * dkdTcell + (1 - transmissioncell) * dBdTcell);
			}
			else if (m_wfunit == SpeciesWFUnitTIR::numberdensity)
			{
				double cellxs = (m_opticalprops->SpeciesCrossSectionCM2AtWavel(startpoint, species, wavelidx)
								 + m_opticalprops->SpeciesCrossSectionCM2AtWavel(endpoint, species, wavelidx)) / 2.0;
				wfcellstorage[species].insert(wfcellstorage[species].begin(), cellxs * ds * 100.0 * (sourcecell - rad) * transmissioncell);
			}
			else  // m_wfunit == SpeciesWFUnitTIR::vmr
			{
				double cellxs = (m_opticalprops->SpeciesCrossSectionCM2AtWavel(startpoint, species, wavelidx)
								 + m_opticalprops->SpeciesCrossSectionCM2AtWavel(endpoint, species, wavelidx)) / 2.0;
				double cellairdensity = (m_opticalprops->AirNumberDensityAtPoint(startpoint)
										 + m_opticalprops->AirNumberDensityAtPoint(endpoint)) / 2.0;
				wfcellstorage[species].insert(wfcellstorage[species].begin(), cellxs * ds * 100.0 * (sourcecell - rad) * cellairdensity * transmissioncell);
			}
		}

		rad = derirad + sourcecell;

	}

	radiance = rad;

	// if weighting functions need to be computed
	if (calcwf)
	{
		opticaldepthtoobserver = 0.0;

		// loop backwards from observer (skip first point, observer location, because the optical depth from that point is 0)
		for (size_t quadpt = 1; quadpt < raystorage->NumQuadraturePoints() - 1; quadpt++)
		{
			// multiply weighting function by transmission to the observer
			for (auto const& species : wf_handles)
			{
				raystorage->SetWFAtPoint(species, quadpt,
					(wfcellstorage[species][quadpt] * exp(-odstorage->at(quadpt - 1)) + wfcellstorage[species][quadpt - 1]) / 2.0 * exp(-opticaldepthtoobserver));
			}
			opticaldepthtoobserver += odstorage->at(quadpt - 1);
		}

		// set weighting function at end points to be equal to the value in the middle of the cell on either end
		for (auto const& species : wf_handles)
		{
			raystorage->SetWFAtPoint(species, 0, wfcellstorage[species][0]);
			raystorage->SetWFAtPoint(species, raystorage->NumQuadraturePoints() - 1, wfcellstorage[species][raystorage->NumQuadraturePoints() - 2] * exp(-opticaldepthtoobserver));
		}
	}

	return ok;
}
