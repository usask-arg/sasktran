#include "include/sktran_montecarlo_internals.h"

SKTRAN_GridDefAirMassFactorShells::SKTRAN_GridDefAirMassFactorShells()
{
	m_extendedToGround = false;
	m_extendedToTOA = false;
}

bool SKTRAN_GridDefAirMassFactorShells::ConfigureHeights(const double* shellAlts, size_t numshells)
{
	return SKTRAN_GridDefRayTracingShells_V21::ConfigureHeights(shellAlts, numshells);
}

bool SKTRAN_GridDefAirMassFactorShells::ConfigureHeights(const std::vector<double>& shellAlts)
{
	return SKTRAN_GridDefRayTracingShells_V21::ConfigureHeights(shellAlts);
}

bool SKTRAN_GridDefAirMassFactorShells::ConfigureHeights(const double* shellAlts, size_t numshells, double surfaceHeight, double toaHeight)
{
	std::vector<double> newShellAlts(shellAlts, shellAlts + numshells);
	if (newShellAlts.front() > surfaceHeight) m_extendedToGround = true; newShellAlts.insert(newShellAlts.begin(), surfaceHeight);
	if (newShellAlts.back() < toaHeight) m_extendedToTOA = true; newShellAlts.push_back(toaHeight);
	return SKTRAN_GridDefRayTracingShells_V21::ConfigureHeights(newShellAlts);
}

bool SKTRAN_GridDefAirMassFactorShells::ConfigureHeights(const std::vector<double>& shellAlts, double surfaceHeight, double toaHeight)
{
	std::vector<double> newShellAlts(shellAlts.begin(), shellAlts.end());
	if (newShellAlts.front() > surfaceHeight) { m_extendedToGround = true; newShellAlts.insert(newShellAlts.begin(), surfaceHeight); }
	if (newShellAlts.back() < toaHeight) { m_extendedToTOA = true; newShellAlts.push_back(toaHeight); }
	return SKTRAN_GridDefRayTracingShells_V21::ConfigureHeights(newShellAlts);
}

SKTRAN_MCAirMassFactorCalculator_Base::SKTRAN_MCAirMassFactorCalculator_Base()
{
	m_coords = nullptr;
	m_amfopticalpropertiestable = nullptr; 
	m_amfopticalpropsintegrator = nullptr; 
	m_amfRayFactory_los = nullptr;
	m_amfRayFactory_solar = nullptr;
	m_amfRayFactory_secondary = nullptr;

	m_loscurved = false;
	m_solarcurved = false;
	m_secondarycurved = false;
}

SKTRAN_MCAirMassFactorCalculator_Base::~SKTRAN_MCAirMassFactorCalculator_Base()
{
	ReleaseResources();
}

bool SKTRAN_MCAirMassFactorCalculator_Base::ReleaseResources()
{
	if (NULL != m_amfopticalpropertiestable) m_amfopticalpropertiestable->Release(); m_amfopticalpropertiestable = NULL;
	if (NULL != m_amfopticalpropsintegrator) m_amfopticalpropsintegrator->Release(); m_amfopticalpropsintegrator = NULL;
	return true;
}

bool SKTRAN_MCAirMassFactorCalculator_Base::TraceRay(const HELIODETIC_VECTOR& observer, const HELIODETIC_UNITVECTOR& look, bool curved, bool optical, SKTRAN_RayOptical_Base* ray) const
{
	bool ok = true;
	ok = ok && ray->MoveObserver(observer, look);
	ok = ok && ray->TraceRay_NewMethod();

	if (optical)
	{
		if (curved)
		{
			std::vector<double> sigmak, sigmaf;
			sigmak.resize(ray->GetNumQuadraturePoints());
			sigmaf.resize(ray->GetNumQuadraturePoints());
			ok = ok && m_amfopticalpropsintegrator->CalculateRayScalarTransmissionVector(ray, NULL, false, true, &sigmak, &sigmaf);
		}
		else
		{
			ok = ok && m_amfopticalpropsintegrator->CalculateRayScalarTransmission_withMinContainer(ray, NULL, false, true);
		}
	}
	return ok;
}


bool SKTRAN_MCAirMassFactorCalculator_Base::FindScatterPoint(const SKTRAN_RayOptical_Base* ray, const HELIODETIC_POINT& scatterPoint, size_t& numCells, double& finalCellWeight) const
{
	bool ok = true;

	const SKTRAN_RayStorage_Base* storage = ray->Storage();

	HELIODETIC_VECTOR scatterVector = scatterPoint.Vector();

	size_t					idx;			// index of the current ray-tracing cell
	HELIODETIC_POINT		point;			// start point of the current ray-tracing cell
	HELIODETIC_VECTOR		vector;			// corresponding vector
	HELIODETIC_UNITVECTOR	look;			// look vector of the current ray-tracing cell
	double					cellLength;		// length of the current ray-tracing cell
	HELIODETIC_VECTOR		diffVector;		// vector from the start of the current ray-tracing cell to the scatter point;
	double					diffProjected;	// diffVector projected in the look direction

	// first check if the scatter point is at the end (should happen quite regularly, whenever it reflects off the ground)
	ok = ok && storage->LocationOfPoint(storage->NumCells(), &point);
	vector = point.Vector();
	if ((scatterVector - vector).Magnitude() < 1.0)
	{
		numCells = storage->NumCells();
		finalCellWeight = 1.0;
		return ok;
	}

	// binary search
	size_t minIdx = 0;
	size_t maxIdx = storage->NumCells() - 1;
	idx = 0;
	while (minIdx <= maxIdx)
	{
		idx = (minIdx + maxIdx) / 2;
		ok = ok && storage->LocationOfPoint(idx, &point);
		look = storage->AverageLookVectorAwayFromObserver(idx);
		cellLength = storage->CellLength(idx);

		vector = point.Vector();
		diffVector = scatterVector - vector;
		diffProjected = diffVector.X() * look.X() + diffVector.Y() * look.Y() + diffVector.Z() * look.Z();

		if (diffProjected < -1)
		{
			maxIdx = idx - 1;
		}
		else if (diffProjected < cellLength)
		{
			break;
		}
		else
		{
			minIdx = idx + 1;
		}
	}

	// make sure the scatter point actually lies on the ray-tracing cell (within a meter)
	ok = ok && diffProjected > -1.0 && diffProjected < cellLength + 1;
	ok = ok && (scatterVector - vector - HELIODETIC_VECTOR(look, diffProjected)).Magnitude() < 1.0;
	ok = ok && idx < storage->NumCells();

	if (ok) numCells = idx + 1;
	if (ok) finalCellWeight = std::max(0.0, std::min(1.0, diffProjected / cellLength));
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_MCAirMassFactorCalculator_Base::ScatterPointIndex, Error finding the scatter point on the ray.");
	return ok;
}

bool SKTRAN_MCAirMassFactorCalculator_Base::IndexRayCells(const SKTRAN_RayOptical_Base* ray, size_t numRayCells, std::vector<size_t>& amfCellIndices) const
{
	bool ok = true;

	const SKTRAN_RayStorage_Base* storage = ray->Storage();

	size_t numAmfCells = m_shellGrid_amf->NumCells();
	ok = ok && numRayCells <= storage->NumCells();

	HELIODETIC_POINT startPoint, endPoint;				// start and end of the current ray-tracing cell (the following parameters correspond to these points)
	size_t startAmfCellIdx, endAmfCellIdx;				// index of the amf shell below (or at) each point (equivalently the index of the amf cell that contains the point)
	double startAltitude, endAltitude;					// altitude of each point
	double startAltAbove;								// altitude of the shell above the start point
	double startAltBelowOrEqual;						// altitude of the shell below or equal the start point
	size_t startRayQuadIdx, endRayQuadIdx;				// ray-tracing quadrature index of each point

	size_t amfCellIdx;	// index of amf cell that contains the current ray-tracing cell

	// find the start of the ray (observer) in the AMF grid
	startRayQuadIdx = 0;
	ok = ok && storage->LocationOfPoint(startRayQuadIdx, &startPoint);
	startAltitude = startPoint.Altitude();
	m_shellGrid_amf->IndexOfPointBelowOrEqual(startAltitude, &startAmfCellIdx);

	// exception: if a point is directly on the highest shell, assign it to the cell below
	// it shouldn't break the logic below since if it points down (which it should) it will land in case 3
	if (startAmfCellIdx == numAmfCells) startAmfCellIdx--;

	// loop through every ray-tracing cell and determine which amf cell it belongs to
	amfCellIndices.resize(numRayCells);
	for (size_t rayCellIdx = 0; ok && rayCellIdx < numRayCells; rayCellIdx++)
	{
		// assume that ray-tracing cells do not cross amf shells (since these are the shells that were used to trace the rays)
		// do not assume that each ray-tracing cell spans an entire amf cell

		// find the shells surrounding the start point
		startAltBelowOrEqual = m_shellGrid_amf->At(startAmfCellIdx);
		startAltAbove = m_shellGrid_amf->At(startAmfCellIdx + 1);

		// find the altitude of the end point
		endRayQuadIdx = rayCellIdx + 1;
		ok = ok && storage->LocationOfPoint(endRayQuadIdx, &endPoint);
		endAltitude = endPoint.Altitude();

		// find which amf cell contains the ray-tracing cell
		if (endAltitude == startAltAbove)				// case 1: ray-tracing cell ends on the upper shell
		{
			endAmfCellIdx = startAmfCellIdx + 1;
			amfCellIdx = startAmfCellIdx; 
		}
		else if (endAltitude < startAltBelowOrEqual)	// case 2: ray-tracing cell starts on the lower shell and then points down 
		{
			endAmfCellIdx = startAmfCellIdx - 1;
			amfCellIdx = endAmfCellIdx; 
		}
		else											// case 3: ray-tracing cell starts and ends within the same amf cell
		{
			amfCellIdx = endAmfCellIdx = startAmfCellIdx;
		}
		ok = ok && amfCellIdx < numAmfCells;

		amfCellIndices[rayCellIdx] = amfCellIdx;

		startAltitude = endAltitude;
		startAmfCellIdx = endAmfCellIdx;
	}

	
	return ok;
}

size_t SKTRAN_MCAirMassFactorCalculator_Base::NumAMFCells() const
{
	size_t numcells = m_shellGrid_amf->NumCells();
	if (m_shellGrid_amf->ExtendedToGround()) numcells -= 1;
	if (m_shellGrid_amf->ExtendedToTOA()) numcells -= 1;
	return numcells;
}

std::vector<double> SKTRAN_MCAirMassFactorCalculator_Base::AMFShellHeights() const
{
	std::vector<double> heights = m_shellGrid_amf->ShellHeight();
	if (m_shellGrid_amf->ExtendedToGround()) heights.erase(heights.begin());
	if (m_shellGrid_amf->ExtendedToTOA()) heights.erase(--heights.end());
	return heights;
}

SKTRAN_MCAirMassFactorCalculator_Length::SKTRAN_MCAirMassFactorCalculator_Length()
{

}

bool SKTRAN_MCAirMassFactorCalculator_Length::TraceMotherRay(SKTRAN_RayOptical_Base const* motherRay)
{
	bool ok = true;
	std::unique_ptr<SKTRAN_RayOptical_Base> temp;
	ok = ok && m_amfRayFactory_los->CreateRayObject(&temp);
	ok = ok && TraceRay(motherRay->GetObserver(), motherRay->LookVector(), m_loscurved, false, temp.get());
	m_losray = std::move(temp);
	return ok;
}

bool SKTRAN_MCAirMassFactorCalculator_Length::AllocateRayOptical(size_t numthreads)
{
	m_solarray.resize(numthreads);
	m_secondaryray.resize(numthreads);
	for (size_t idx = 0; idx < numthreads; idx++)
	{
		m_amfRayFactory_solar->CreateRayObject(&m_solarray[idx]);
		m_amfRayFactory_secondary->CreateRayObject(&m_secondaryray[idx]);
	}
	return true;
}

bool SKTRAN_MCAirMassFactorCalculator_Length::CalculateSlantContribution(size_t order, SKTRAN_RayOptical_Base const* scatterRay, const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* mcphoton, size_t threadid)
{
	bool ok = true;

	HELIODETIC_UNITVECTOR sun;
	const SKTRAN_RayOptical_Base* ray;
	const SKTRAN_RayStorage_Base* storage;

	size_t numCellsUpToScatter;		// number of ray-tracing cells up to and including the cell containing the scatter point
	double fraction;				// the fraction of the length of the final cell that occurs before the scatter point

	// if the user-specified amf shells do not reach the ground and/or toa, the amf shells are extended for ray-tracing purposes
	// therefore there are two sets of indices: 
	//	extended (where 0 refers to the cell added between the lowest amf cell and the ground) 
	//  unextended (where 0 refers to the lowest amf cell)
	size_t rayCellIdx;						// index for ray-tracing cells
	size_t amfCellIdx;						// index of amf cells
	std::vector<size_t> amfCellIndices;		// vector of amf cell indices corresponding to each ray-tracing cell (extended)
	size_t numAmfCells = NumAMFCells();		// number of amf cells (unextended)
	size_t ground = m_shellGrid_amf->ExtendedToGround() ? 1 : 0; // used to decrement extended amf cell indices to make unextended amf cell indices

	// find slant path contributions from the ray between the previous scatter point (observer) and the current scatter point
	if (order == 1)
	{
		ray = m_losray.get();
	}
	else
	{
		ok = ok && TraceRay(scatterRay->GetObserver(), scatterRay->LookVector(), m_secondarycurved, false, m_secondaryray[threadid].get());
		ray = m_secondaryray[threadid].get();
	}

	ok = ok && FindScatterPoint(ray, scatterPoint, numCellsUpToScatter, fraction);
	ok = ok && IndexRayCells(ray, numCellsUpToScatter, amfCellIndices);

	if (ok && numCellsUpToScatter > 0)
	{
		storage = ray->Storage();
		for (rayCellIdx = 0; rayCellIdx < numCellsUpToScatter - 1; rayCellIdx++)
		{
			amfCellIdx = amfCellIndices[rayCellIdx] - ground;
			if (amfCellIdx < numAmfCells) mcphoton->m_scatterSlantColumns[amfCellIdx] += storage->CellLength(rayCellIdx);
		}
		amfCellIdx = amfCellIndices[numCellsUpToScatter - 1] - ground;
		if (amfCellIdx < numAmfCells) mcphoton->m_scatterSlantColumns[amfCellIdx] += storage->CellLength(numCellsUpToScatter - 1) * fraction;
	}

	// find slant path contributions from the ray from the sun to the current scatter point
	sun.SetCoords(0.0, 0.0, 1.0);
	ok = ok && TraceRay(scatterPoint.Vector(), sun, m_solarcurved, false, m_solarray[threadid].get());
	std::fill(mcphoton->m_solarSlantColumns.begin(), mcphoton->m_solarSlantColumns.end(), 0.0);

	ray = m_solarray[threadid].get();
	storage = ray->Storage();

	if (ok && storage->NumCells() > 0)
	{
		ok = ok && IndexRayCells(ray, storage->NumCells(), amfCellIndices);
		if (ok)
		{			
			for (rayCellIdx = 0; rayCellIdx < storage->NumCells(); rayCellIdx++)
			{
				amfCellIdx = amfCellIndices[rayCellIdx] - ground;
				if (amfCellIdx < numAmfCells) mcphoton->m_solarSlantColumns[amfCellIdx] += storage->CellLength(rayCellIdx);
			}
		}
	}

	return ok;
}

bool SKTRAN_MCAirMassFactorCalculator_Length::AllocatePhotons(std::vector<std::unique_ptr<SKTRAN_MCPhoton_Base>>& mcphotons) const
{
	bool ok = true;
	size_t n = NumAMFCells();
	for (auto it = mcphotons.begin(); it != mcphotons.end(); it++)
	{
		(*it)->m_solarSlantColumns.resize(n);
		(*it)->m_scatterSlantColumns.resize(n);
		ok = ok && ClearPhoton((*it).get());
	}
	return ok;
}

bool SKTRAN_MCAirMassFactorCalculator_Length::ClearPhoton(SKTRAN_MCPhoton_Base* mcphoton) const
{
	std::fill(mcphoton->m_solarSlantColumns.begin(), mcphoton->m_solarSlantColumns.end(), 0.0);
	std::fill(mcphoton->m_scatterSlantColumns.begin(), mcphoton->m_scatterSlantColumns.end(), 0.0);
	return true;
}


bool SKTRAN_MCAirMassFactorCalculator_Length::InitializeLogger(SKTRAN_MCAirMassFactorLogger* logger) const
{
	bool ok = true;
	size_t numcells = NumAMFCells();
	std::vector<double> heights = AMFShellHeights();
	ok = ok && numcells > 0 && numcells == heights.size() - 1;
	if (ok)
	{
		std::vector<double> layerThickness(heights.size() - 1);
		for (size_t idx = 0; idx < numcells; idx++)
			layerThickness[idx] = heights[idx + 1] - heights[idx];
		logger->Initialize(numcells, layerThickness);
	}
	return ok;
}

SKTRAN_MCAirMassFactorCalculator_OpticalDepth::SKTRAN_MCAirMassFactorCalculator_OpticalDepth()
{

}

bool SKTRAN_MCAirMassFactorCalculator_OpticalDepth::CalculateOpticalPropertiesTable(double wavelen, SKTRAN_AtmosphericOpticalState_V21* opticalstate, bool userupdateclimatology)
{
	bool ok = true;
	SKTRAN_AtmosphericOpticalState_V21 amfopticalstate;
	skClimatology* amfclimatology;
	skOpticalProperties* amfopticalproperties;
	GEODETIC_INSTANT point;

	// construct the AMF optical state (containing only the AMF species)
	opticalstate->GetSpeciesClimatology(m_amfspecieshandle, &amfclimatology);
	opticalstate->GetSpeciesOpticalProperties(m_amfspecieshandle, &amfopticalproperties);
	amfopticalstate.AddSpecies(m_amfspecieshandle, amfclimatology, amfopticalproperties);

	// Get reference point for filling optical property tables
	const SKTRAN_CoordinateTransform_V2* coords = Coordinates();
	point.latitude = coords->ReferencePtLatitude();
	point.longitude = coords->ReferencePtLongitude();
	point.heightm = 0.0;
	point.mjd = coords->ReferencePointMJD();
	NXASSERT((point.mjd > 10000.0));
	if (point.mjd < 10000.0)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_MCAirMassFactorCalculator_OpticalDepth::CalculateOpticalPropertiesTable the mjd being used for the climatologies is probably out of range. Its value is %e", (double)point.mjd);
	}

	// Fill optical property tables
	ok = ok && amfopticalstate.SetTimeAndLocation(point, userupdateclimatology);
	ok = ok && m_amfopticalpropertiestable->ConfigureOptical(wavelen, amfopticalstate);

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_MCAirMassFactorCalculator_OpticalDepth::CalculateRadiance, Error calculating the AMF optical properties table. Thats not good");
	}

	return ok;
}

bool SKTRAN_MCAirMassFactorCalculator_OpticalDepth::TraceMotherRay(SKTRAN_RayOptical_Base const* motherRay)
{
	bool ok = true;
	std::unique_ptr<SKTRAN_RayOptical_Base> temp;
	ok = ok && m_amfRayFactory_los->CreateRayObject(&temp);
	ok = ok && TraceRay(motherRay->GetObserver(), motherRay->LookVector(), m_loscurved, true, temp.get());
	m_losray = std::move(temp);
	return ok;
}

bool SKTRAN_MCAirMassFactorCalculator_OpticalDepth::AllocateRayOptical(size_t numthreads)
{
	m_solarray.resize(numthreads);
	m_secondaryray.resize(numthreads);
	for (size_t idx = 0; idx < numthreads; idx++)
	{
		m_amfRayFactory_solar->CreateRayObject(&m_solarray[idx]);
		m_amfRayFactory_secondary->CreateRayObject(&m_secondaryray[idx]);
	}
	return true;
}

bool SKTRAN_MCAirMassFactorCalculator_OpticalDepth::PartialOpticalDepth(SKTRAN_RayOptical_Base const* ray, size_t cellindex, double fraction, double& opticaldepth) const
{
	bool ok = true;

	if (fraction <= 0)
	{
		opticaldepth = 0;
		return ok;
	}
	else if (fraction >= 1)
	{
		std::vector<double> odarray = ray->OpticalDepthArray();
		ok = ok && cellindex < odarray.size() - 1;
		if (ok) opticaldepth = odarray[cellindex + 1] - odarray[cellindex];
		return ok;
	}
	else
	{
		ok = ok && m_amfopticalpropsintegrator->CalculatePartialOpticalDepth(ray, cellindex, fraction, opticaldepth);
		return ok;
	}
}

bool SKTRAN_MCAirMassFactorCalculator_OpticalDepth::CalculateSlantContribution(size_t order, SKTRAN_RayOptical_Base const* scatterRay, const HELIODETIC_POINT& scatterPoint, SKTRAN_MCPhoton_Base* mcphoton, size_t threadid)
{
	bool ok = true;

	HELIODETIC_UNITVECTOR sun;	
	const SKTRAN_RayOptical_Base* ray;
	
	size_t numCellsUpToScatter;		// number of ray-tracing cells up to and including the cell containing the scatter point
	double fraction;				// the fraction of the length of the final cell that occurs before the scatter point
	double partialod;				// the optical depth of the final cell up to the scatter point

	// if the user-specified amf shells do not reach the ground and/or toa, the amf shells are extended for ray-tracing purposes
	// therefore there are two sets of indices: 
	//	extended (where 0 refers to the cell added between the lowest amf cell and the ground) 
	//  unextended (where 0 refers to the lowest amf cell)
	size_t rayCellIdx;						// index for ray-tracing cells
	size_t amfCellIdx;						// index of amf cells
	std::vector<size_t> amfCellIndices;		// vector of amf cell indices corresponding to each ray-tracing cell (extended)
	size_t numAmfCells = NumAMFCells();		// number of amf cells (unextended)
	size_t ground = m_shellGrid_amf->ExtendedToGround() ? 1 : 0; // used to decrement extended amf cell indices to make unextended amf cell indices

	// find slant path contributions from the ray between the previous scatter point (observer) and the current scatter point
	if (order == 1)
	{
		ray = m_losray.get();
	}
	else
	{
		ok = ok && TraceRay(scatterRay->GetObserver(), scatterRay->LookVector(), m_secondarycurved, true, m_secondaryray[threadid].get());
		ray = m_secondaryray[threadid].get();
	}

	ok = ok && FindScatterPoint(ray, scatterPoint, numCellsUpToScatter, fraction);
	ok = ok && IndexRayCells(ray, numCellsUpToScatter, amfCellIndices);

	if (ok && numCellsUpToScatter > 0) 
	{
		const std::vector<double>& odarray = ray->OpticalDepthArray(); // cumulative optical depth
		// scatterSlantColumns contains previous scatter segments - do not reset here

		for (rayCellIdx = 0; rayCellIdx < numCellsUpToScatter - 1; rayCellIdx++)
		{
			amfCellIdx = amfCellIndices[rayCellIdx] - ground;
			if (amfCellIdx < numAmfCells) mcphoton->m_scatterSlantColumns[amfCellIdx] += odarray[rayCellIdx + 1] - odarray[rayCellIdx];
		}

		rayCellIdx = numCellsUpToScatter - 1;
		amfCellIdx = amfCellIndices[rayCellIdx] - ground;
		ok = ok && PartialOpticalDepth(ray, rayCellIdx, fraction, partialod);
		partialod = fraction * (odarray[rayCellIdx + 1] - odarray[rayCellIdx]);
		if (amfCellIdx < numAmfCells) mcphoton->m_scatterSlantColumns[amfCellIdx] += partialod;
	}

	// find slant path contributions from the ray from the sun to the current scatter point
	sun.SetCoords(0.0, 0.0, 1.0);
	ok = ok && TraceRay(scatterPoint.Vector(), sun, m_solarcurved, true, m_solarray[threadid].get());

	ray = m_solarray[threadid].get();
	const std::vector<double>& odarray = ray->OpticalDepthArray();

	ok = ok && IndexRayCells(ray, ray->Storage()->NumCells(), amfCellIndices);

	if (ok && odarray.size() > 0)
	{
		const std::vector<double>& odarray = ray->OpticalDepthArray(); // cumulative optical depth
		// solarSlantColumns only holds the current segment, so we need to reset it here
		std::fill(mcphoton->m_solarSlantColumns.begin(), mcphoton->m_solarSlantColumns.end(), 0.0);

		// optical depth must be calculated by subtracting consecutive cumulative optical depths
		for (rayCellIdx = 0; rayCellIdx < ray->Storage()->NumCells(); rayCellIdx++)
		{
			amfCellIdx = amfCellIndices[rayCellIdx] - ground;
			if (amfCellIdx < numAmfCells) mcphoton->m_solarSlantColumns[amfCellIdx] += odarray[rayCellIdx + 1] - odarray[rayCellIdx];
		}
	}

	return ok;
}

bool SKTRAN_MCAirMassFactorCalculator_OpticalDepth::AllocatePhotons(std::vector<std::unique_ptr<SKTRAN_MCPhoton_Base>>& mcphotons) const
{
	bool ok = true;
	size_t n = NumAMFCells();
	for (auto it = mcphotons.begin(); it != mcphotons.end(); it++)
	{
		(*it)->m_solarSlantColumns.resize(n);
		(*it)->m_scatterSlantColumns.resize(n);
		ok = ok && ClearPhoton((*it).get());
	}
	return ok;
}

bool SKTRAN_MCAirMassFactorCalculator_OpticalDepth::ClearPhoton(SKTRAN_MCPhoton_Base* mcphoton) const
{
	std::fill(mcphoton->m_solarSlantColumns.begin(), mcphoton->m_solarSlantColumns.end(), 0.0);
	std::fill(mcphoton->m_scatterSlantColumns.begin(), mcphoton->m_scatterSlantColumns.end(), 0.0);
	return true;
}

bool SKTRAN_MCAirMassFactorCalculator_OpticalDepth::InitializeLogger(SKTRAN_MCAirMassFactorLogger* logger) const
{
	bool ok = true;
	std::vector<double> vod;
	size_t numcells = NumAMFCells();

	ok = ok && numcells > 0;

	ok = ok && VerticalOpticalDepth(vod);
	if (ok) logger->Initialize(numcells, vod);
	return ok;
}

bool SKTRAN_MCAirMassFactorCalculator_OpticalDepth::VerticalOpticalDepth(std::vector<double>& vod) const
{
	// TODO there is some overlap of this function with CalculateSlantContribution which should be condensed
	bool ok = true;
	std::vector<double> heights;
	std::vector<double> opticaldepths;
	HELIODETIC_POINT pt;
	double ext0, ext1;

	size_t ground = m_shellGrid_amf->ExtendedToGround() ? 1 : 0; // used to decrement extended amf cell indices to make unextended amf cell indices
	size_t numAmfCells = NumAMFCells();		// number of amf cells (unextended)

	HELIODETIC_POINT ref = m_amfopticalpropertiestable->CoordinatesPtr()->ReferencePoint(0.0);

	std::unique_ptr<SKTRAN_RayOptical_Base> ray;
	m_amfRayFactory_secondary->CreateRayObject(&ray);
	ok = ok && TraceRay(ref.Vector(), ref.LocalZenith(), m_secondarycurved, true, ray.get());
	const SKTRAN_RayStorage_Base* storage = ray->Storage();

	std::vector<size_t> amfCellIndices;
	ok = ok && IndexRayCells(ray.get(), storage->NumCells(), amfCellIndices);

	size_t amfCellIdx;
	vod.resize(numAmfCells, 0.0);
	const std::vector<double>& odarray = ray->OpticalDepthArray();
	for (size_t rayCellIdx = 0; rayCellIdx < storage->NumCells(); rayCellIdx++)
	{
		amfCellIdx = amfCellIndices[rayCellIdx] - ground;
		if (amfCellIdx < numAmfCells) vod[amfCellIdx] += odarray[rayCellIdx + 1] - odarray[rayCellIdx];
	}

	return ok;

	//heights = AMFShellHeights();
	//opticaldepths.resize(NumAMFCells(), 0.0);

	//pt = Coordinates()->ReferencePoint(heights[0]);
	//ext0 = m_amfopticalpropertiestable->TotalExtinctionPerCM(pt);
	//for (size_t i = 0; i < NumAMFCells(); i++) {
	//	pt = Coordinates()->ReferencePoint(heights[i+1]);
	//	ext1 = m_amfopticalpropertiestable->TotalExtinctionPerCM(pt);
	//	opticaldepths[i] = 50.0 * (ext0 + ext1) * (heights[i + 1] - heights[i]);
	//	ext0 = ext1;
	//}
	//return opticaldepths;
}

bool SKTRAN_MCAirMassFactorCalculator_DoNothing::InitializeLogger(SKTRAN_MCAirMassFactorLogger* logger) const
{
	bool ok = true;
	size_t numcells = 0;
	std::vector<double> vod;
	if (ok) logger->Initialize(numcells, vod);
	return ok;
}