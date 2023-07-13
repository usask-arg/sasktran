#include "include/sktran_montecarlo_internals.h"

SKTRAN_MCPhoton_RadInfo::SKTRAN_MCPhoton_RadInfo()
{
	m_vec.SetTo(0.0);
}


SKTRAN_MCPhoton_Base::SKTRAN_MCPhoton_Base ( )
{
	Initialize();
}

SKTRAN_MCPhoton_Base::SKTRAN_MCPhoton_Base( const SKTRAN_MCPhoton_Base& other )
{
	//if ( other.m_photonOptical.get() != nullptr)
	//{
	//	nxLog::Record(NXLOG_WARNING,"SKTRAN_MCPhoton::SKTRAN_MCPhoton, Copy constructor cannot move non-null other.m_photonOptical. Code needs to be re-worked");
	//	// m_photonOptical = std::move(other.m_photonOptical);
	//}
	Configure(other);
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_MCPhoton::operator=		 2014- 11- 15*/
/** The assignmenmt operator is really treated as a move operator
**/
/*---------------------------------------------------------------------------*/

SKTRAN_MCPhoton_Base& SKTRAN_MCPhoton_Base::operator= ( SKTRAN_MCPhoton_Base& other )
{
	m_photonOptical = std::move(other.m_photonOptical);
	
	Configure(other);

	return *this;
}

SKTRAN_MCPhoton_Base::SKTRAN_MCPhoton_Base( std::unique_ptr<SKTRAN_RayOptical_Base> r ) 
{
	Initialize( );
	SetOpticalRay( std::move(r));
}

bool SKTRAN_MCPhoton_Base::FindScatterPointCellIndex( const SKTRAN_RayStorage_Base* storage, const std::vector<double>& odarray, const double & targetTau, size_t& scatterCellIndex, HELIODETIC_POINT& scatterCellStartPoint ) const
{
	bool ok = true;

	std::vector<double>::const_iterator		exitod_iter;

	double totalod = odarray.back();
	double target = targetTau;

	if (0.0 == target)
	{
		exitod_iter = odarray.begin() + 1;
	}
	else if (target >= totalod)
	{
		// A weird special case that can happen when the ray optical depth is close to machine precision and rounding errors occur when transforming
		// the probability space
		target = totalod;
		exitod_iter = odarray.end() - 1;
	}
	else
	{
		exitod_iter = std::lower_bound(odarray.begin(), odarray.end(), target);
	}

	size_t endinterceptindex = (exitod_iter - odarray.begin());
	size_t tangentIndex = exitod_iter - odarray.begin() - 1;

	if (endinterceptindex == 0)
	{
		tangentIndex = 0;
		endinterceptindex = 1;
	}

	//size_t startinterceptindex = endinterceptindex - 1;
	scatterCellIndex = endinterceptindex - 1;
	ok = ok && storage->LocationOfPoint(scatterCellIndex, &scatterCellStartPoint);

	return ok;
}

bool SKTRAN_MCPhoton_Base::ConfigureQuadratureCoefficients(const SKTRAN_RayStorage_Base* storage, const HELIODETIC_POINT& scatterPoint, const size_t& scatterCellIndex, const HELIODETIC_POINT& scatterCellStartPoint, SKTRAN_OpticalDepthCalculator_LinearWithHeight& odcalculator) const
{
	bool ok = true;

	size_t startinterceptindex = scatterCellIndex;
	size_t endinterceptindex = scatterCellIndex + 1;

	// distance from startintercept to scatterPoint
	double finalcelldistance = (scatterPoint.Vector() - scatterCellStartPoint.Vector()).Magnitude();

	double rt = storage->RadiusOfCellTangentPoint(startinterceptindex);
	double r0 = scatterCellStartPoint.Radius();
	double r1 = scatterPoint.Radius();
	double t0 = storage->DistanceOfPointFromCellTangentPoint(startinterceptindex, startinterceptindex);
	double t2 = storage->DistanceOfPointFromCellTangentPoint(endinterceptindex, startinterceptindex);
	double t1 = t0 > t2 ? t0 - finalcelldistance : t0 + finalcelldistance;
	odcalculator.ConfigureQuadratureCoefficients(r0, r1, t0, t1, rt);

	return ok;
}


SKTRAN_MCPhoton_Base::~SKTRAN_MCPhoton_Base ( )
{ 
//	ReleaseResources( );
}

void SKTRAN_MCPhoton_Base::Initialize( )
{
	m_photonOptical = nullptr;
	//m_photonRadiance.resize(1);
	//m_photonSource.resize(1);
	//m_scatterWeight.resize(1);
	//m_weightFactor.resize(1);
	//m_albedo.resize(1);
	//m_transmission.resize(1);
	//m_opticaldepth.resize(1);
	//m_scatterFactor.resize(1);
	//m_finalWavelength.resize(1);
	//m_currentWavelength.resize(1);

	//m_eScatterWeight.resize(1);
	//m_eWeightFactor.resize(1);
	//m_eAlbedo.resize(1);
	//m_eTransmission.resize(1);
	//m_eOpticalDepth.resize(1);

	ResetRadiance(); 
	m_basis.x.Clear();
	m_basis.y.Clear();
	m_basis.z.Clear();
	m_distanceProb  = 1.0;
	m_scatterVector.SetCoords(-9999.9, -9999.9, -9999.9);
	m_targetTau     = 0.0;
}

//void SKTRAN_MCPhoton::ReleaseResources( )
//{
//	if(m_manageOpticalLifetime)	m_photonOptical.release();
//}

SKTRAN_RayOptical_Base* SKTRAN_MCPhoton_Base::photonOptical ( )
{
	return m_photonOptical.get(); 
}

const SKTRAN_RayOptical_Base* SKTRAN_MCPhoton_Base::photonOptical ( ) const
{
	return m_photonOptical.get(); 
}

SKTRAN_MCPhoton_RadInfo& SKTRAN_MCPhoton_Base::photonRadiance  ( )      
{
	return *(m_photonRadiance.begin() + m_primaryWavelengthIndex);
}

const SKTRAN_MCPhoton_RadInfo& SKTRAN_MCPhoton_Base::photonRadiance ( ) const
{
	return *(m_photonRadiance.begin() + m_primaryWavelengthIndex);
}

std::vector<SKTRAN_MCPhoton_RadInfo>& SKTRAN_MCPhoton_Base::photonRadiances( )
{
	return m_photonRadiance;
}

const std::vector<SKTRAN_MCPhoton_RadInfo>& SKTRAN_MCPhoton_Base::photonRadiances( ) const
{
	return m_photonRadiance;
}

SKTRAN_MCPhoton_RadInfo& SKTRAN_MCPhoton_Base::photonSource()
{
	return *(m_photonSource.begin() + m_primaryWavelengthIndex);
}

const SKTRAN_MCPhoton_RadInfo& SKTRAN_MCPhoton_Base::photonSource() const
{
	return *(m_photonSource.begin() + m_primaryWavelengthIndex);
}

std::vector<SKTRAN_MCPhoton_RadInfo>& SKTRAN_MCPhoton_Base::photonSources()
{
	return m_photonSource;
}

const std::vector<SKTRAN_MCPhoton_RadInfo>& SKTRAN_MCPhoton_Base::photonSources() const
{
	return m_photonSource;
}

void SKTRAN_MCPhoton_Base::ResetRadiance ( )

{
	// Start with no radiance contribution
	//m_photonRadiance.SetVector(0.0);

	ResetFactors();

	for (auto&& it : m_photonRadiance) it.SetVector(0.0);
	std::fill(m_scatterWeight.begin(), m_scatterWeight.end(), 1.0);	
	
	auto finalwl = m_finalWavelength.begin();
	for (auto&& currentwl : m_currentWavelength) currentwl = *(finalwl++);
	
	// start witht the ray set to the primary wavelength
	if (m_photonOptical != nullptr) m_photonOptical->SetWavelength(m_primaryWavelength);

	// Start with an non-polarizing path in the original basis
	m_toObserverOp.SetToIdentity();

	// Start with the manual scatter selection factor set to 1.0
	m_optFactor = 1.0;
}

void SKTRAN_MCPhoton_Base::ResetFactors()
{
	for (auto&& it : m_photonSource) it.SetVector(0.0);
	std::fill(m_weightFactor.begin(), m_weightFactor.end(), 1.0);
	std::fill(m_albedo.begin(), m_albedo.end(), 1.0);
	std::fill(m_transmission.begin(), m_transmission.end(), 1.0);
	std::fill(m_scatterFactor.begin(), m_scatterFactor.end(), 1.0);
	for (auto&& it : m_opticaldepth) std::fill(it.begin(), it.end(), 0.0);
}

bool SKTRAN_MCPhoton_Base::CalculateAlbedo(const SKTRAN_TableOpticalProperties_Base * opticalpropertiestable, const SKTRAN_TableOpticalProperties_MCBase * mcOptTable, const HELIODETIC_POINT & scatterPoint)
{
	bool ok = true;
	double wl = m_currentWavelength[m_primaryWavelengthIndex];
	double albedo = opticalpropertiestable->ScatteringExtinctionPerCM(wl, scatterPoint);
	albedo += mcOptTable->InelasticProperties()->InelasticExtinctionPerCM(wl, scatterPoint);
	albedo /= opticalpropertiestable->TotalExtinctionPerCM(wl, scatterPoint);		// Can't use the simplification we used in the single LOS code because we need #albedo w/o scattCoeff later
	std::fill(m_albedo.begin(), m_albedo.end(), albedo);
	m_isGroundScatter = false;
	return ok;
}

bool SKTRAN_MCPhoton_Base::SetOpticalRay ( std::unique_ptr<SKTRAN_RayOptical_Base> r ) 
{
//		ReleaseResources(); 
    	m_photonOptical = std::move(r);
		if (m_photonOptical) m_photonOptical->SetWavelength(m_primaryWavelength);
		return true;
}

bool SKTRAN_MCPhoton_Base::DefineRayBasis( )
{
    bool ok = true;
	HELIODETIC_POINT pt;

	ok = ok && nullptr!=m_photonOptical;
	NXASSERT((m_photonOptical->Coordinates() != nullptr));
	ok = ok && m_photonOptical->Coordinates()->HelioVectorToHelioPoint(m_photonOptical->GetObserver(), &pt );

	if(ok){
		m_basis.ProduceBasis( pt, m_photonOptical->LookVector() );
	}

    return ok;
}


const SKTRAN_MCBasis& SKTRAN_MCPhoton_Base::GetBasis ( ) const
{
    return m_basis;
}

SKTRAN_MCBasis& SKTRAN_MCPhoton_Base::GetBasisVar ( )
{
    return m_basis;
}

void SKTRAN_MCPhoton_Base::AddPhaseOpInPath ( const SKTRAN_ScatMat_MIMSNC& op )
{
	m_toObserverOp.RMultBy(op);
}

void SKTRAN_MCPhoton_Base::AddRotateOpInPath ( const SKTRAN_ScatMat_Rot& op )
{
	m_toObserverOp.RMultBy(op);
}

void SKTRAN_MCPhoton_Base::LeftApplyPathTo  ( SKTRAN_Stokes_NC& source ) const
{
	source = m_toObserverOp * source;
}


const SKTRAN_PhaseMat_MIMSNC&  SKTRAN_MCPhoton_Base::GetPathOp ( ) const
{
	return m_toObserverOp;
}

SKTRAN_PhaseMat_MIMSNC&        SKTRAN_MCPhoton_Base::GetPathOp ( )
{
	return m_toObserverOp;
}

//bool SKTRAN_MCPhoton_Base::SetCurrentWavelength(double wavelength)
//{
//	bool ok = true;
//
//	if (m_numWavelengths > 1)
//	{
//		auto it = std::find(m_finalWavelength.begin(), m_finalWavelength.end(), wavelength);
//		ok = ok && it != m_finalWavelength.end();
//		if (ok)
//		{
//			m_primaryWavelengthIndex = it - m_finalWavelength.begin();
//			m_primaryWavelength = wavelength;
//			if (m_photonOptical) m_photonOptical->SetWavelength(wavelength);
//		}
//	}
//	else
//	{
//		m_primaryWavelengthIndex = 0;
//		m_primaryWavelength = wavelength;
//		if (m_photonOptical) m_photonOptical->SetWavelength(wavelength);
//		m_currentWavelength = { wavelength };
//		m_finalWavelength = { wavelength };
//	}
//	return ok;
//}

bool SKTRAN_MCPhoton_Base::Configure(const SKTRAN_MCPhoton_Base & other)
{
	bool ok = true;

	m_numWavelengths = other.m_numWavelengths;
	m_primaryWavelengthIndex = other.m_primaryWavelengthIndex;
	m_primaryWavelength = other.m_primaryWavelength;
	m_manualScatter = other.m_manualScatter;

	m_photonRadiance = other.m_photonRadiance;
	m_photonSource = other.m_photonSource;
	m_currentWavelength = other.m_currentWavelength;
	m_finalWavelength = other.m_finalWavelength;
	m_scatterWeight = other.m_scatterWeight;
	m_weightFactor = other.m_weightFactor;
	m_albedo = other.m_albedo;
	m_transmission = other.m_transmission;
	m_scatterFactor = other.m_scatterFactor;
	m_opticaldepth = other.m_opticaldepth;
	m_optFactor = other.m_optFactor;

	m_basis.x = other.m_basis.x;
	m_basis.y = other.m_basis.y;
	m_basis.z = other.m_basis.z;

	m_solarSlantColumns = other.m_solarSlantColumns;
	m_scatterSlantColumns = other.m_scatterSlantColumns;

	m_distanceProb = other.m_distanceProb;
	m_scatterVector.SetCoords(other.m_scatterVector.X(), other.m_scatterVector.Y(), other.m_scatterVector.Z());
	m_toObserverOp = other.m_toObserverOp;
	
	return ok;
}

bool SKTRAN_MCPhoton::TraceRays(const SKTRAN_OpticalPropertiesIntegrator_Base* integrator, bool curved)
{
	bool ok = true;
	
	ok = ok && m_photonOptical->TraceRay_NewMethod();

	if (curved)
	{
		std::vector<double> sigmak(m_photonOptical->GetNumQuadraturePoints()); // should allocate persistent vectors if this is ever used
		std::vector<double> sigmaf(m_photonOptical->GetNumQuadraturePoints());

		for (size_t wlidx = 0; wlidx < m_numWavelengths; wlidx++)
		{
			if (wlidx != m_primaryWavelengthIndex)
			{
				m_photonOptical->SetWavelength(m_currentWavelength[wlidx]);
				ok = ok && integrator->CalculateRayScalarTransmissionVector(m_photonOptical.get(), NULL, false, true, &sigmak, &sigmaf);
				m_opticaldepth[wlidx] = m_photonOptical->OpticalDepthArray();
			}
		}
		m_photonOptical->SetWavelength(m_primaryWavelength);
		ok = ok && integrator->CalculateRayScalarTransmissionVector(m_photonOptical.get(), NULL, false, true, &sigmak, &sigmaf); // integrate the primary wavelength last so that photonOptical remembers the optical depth for each cell
		m_opticaldepth[m_primaryWavelengthIndex] = m_photonOptical->OpticalDepthArray();
	}
	else
	{
		for (size_t wlidx = 0; wlidx < m_numWavelengths; wlidx++)
		{
			if (wlidx != m_primaryWavelengthIndex)
			{
				m_photonOptical->SetWavelength(m_currentWavelength[wlidx]);
				ok = ok && integrator->CalculateRayScalarTransmission_withMinContainer(m_photonOptical.get(), NULL, false, true);
				m_opticaldepth[wlidx] = m_photonOptical->OpticalDepthArray();
			}
		}
		m_photonOptical->SetWavelength(m_currentWavelength[m_primaryWavelengthIndex]);
		ok = ok && integrator->CalculateRayScalarTransmission_withMinContainer(m_photonOptical.get(), NULL, false, true); // integrate the primary wavelength last so that photonOptical remembers the optical depth for each cell
		m_opticaldepth[m_primaryWavelengthIndex] = m_photonOptical->OpticalDepthArray();
	}

	return ok;
}

bool SKTRAN_MCPhoton::UpdateScatterWeight(double forcedScatterCorrFactor)
{
	bool ok = true;

	double weightFactor;
	for (size_t idx = 0; idx < m_numWavelengths; idx++)
	{
		weightFactor = m_albedo[idx];													// single scatter albedo for atmo scatter, surface albedo for ground scatter
		weightFactor *= forcedScatterCorrFactor;										// forced scatter correction factor
		//weightFactor *= m_transmission[idx] / m_transmission[m_primaryWavelengthIndex]; // simultaneous-wavelength factor: from scatter direction and/or wavelength-shift distributions
		//weightFactor *= m_scatterFactor[idx];											// simultaneous-wavelength factor: from picking a scatter point from the transmission function
		weightFactor *= m_optFactor;													// optimization factor: from forcing an elastic or inelastic scatter
		m_weightFactor[idx] = weightFactor;
		m_scatterWeight[idx] *= weightFactor;
	}

	return ok;
}

bool SKTRAN_MCPhoton::SetWavelengths(const std::vector<double>& wavelengths)
{
	bool ok = true;

	m_numWavelengths = 1;

	m_currentWavelength.resize(m_numWavelengths);
	m_finalWavelength.resize(m_numWavelengths);
	m_photonRadiance.resize(m_numWavelengths);
	m_photonSource.resize(m_numWavelengths);
	m_currentWavelength.resize(m_numWavelengths);
	m_scatterWeight.resize(m_numWavelengths);
	m_weightFactor.resize(m_numWavelengths);
	m_albedo.resize(m_numWavelengths);
	m_transmission.resize(m_numWavelengths);
	m_scatterFactor.resize(m_numWavelengths);
	m_opticaldepth.resize(m_numWavelengths);

	ResetRadiance();

	return ok;
}

bool SKTRAN_MCPhoton::SetCurrentWavelength(double wavelength)
{
	bool ok = true;
	m_primaryWavelengthIndex = 0;
	m_primaryWavelength = wavelength;
	if (m_photonOptical) m_photonOptical->SetWavelength(wavelength);
	m_currentWavelength = { wavelength };
	m_finalWavelength = { wavelength };
	return ok;
}


SKTRAN_MCPhoton_Ring::SKTRAN_MCPhoton_Ring(const SKTRAN_MCPhoton_Ring & other) : SKTRAN_MCPhoton_Base(other)
{
	m_ePhotonRadiance = other.m_ePhotonRadiance;
	m_ePhotonSource = other.m_ePhotonSource;
	m_eScatterWeight = other.m_eScatterWeight;
	m_eWeightFactor = other.m_eWeightFactor;
	m_eAlbedo = other.m_eAlbedo;
	m_eTransmission = other.m_eTransmission;
	m_eOpticalDepth = other.m_eOpticalDepth;
}

SKTRAN_MCPhoton_RadInfo& SKTRAN_MCPhoton_Ring::photonRadiance(bool elasticRaman)
{
	return elasticRaman ? *(m_ePhotonRadiance.begin() + m_primaryWavelengthIndex) : SKTRAN_MCPhoton_Base::photonRadiance();
}

const SKTRAN_MCPhoton_RadInfo& SKTRAN_MCPhoton_Ring::photonRadiance(bool elasticRaman) const
{
	return elasticRaman ? *(m_ePhotonRadiance.begin() + m_primaryWavelengthIndex) : SKTRAN_MCPhoton_Base::photonRadiance();
}

std::vector<SKTRAN_MCPhoton_RadInfo>& SKTRAN_MCPhoton_Ring::photonRadiances(bool elasticRaman)
{
	return elasticRaman ? m_ePhotonRadiance : SKTRAN_MCPhoton_Base::photonRadiances();
}

const std::vector<SKTRAN_MCPhoton_RadInfo>& SKTRAN_MCPhoton_Ring::photonRadiances(bool elasticRaman) const
{
	return elasticRaman ? m_ePhotonRadiance : SKTRAN_MCPhoton_Base::photonRadiances();
}

SKTRAN_MCPhoton_RadInfo& SKTRAN_MCPhoton_Ring::photonSource(bool elasticRaman)
{
	return elasticRaman ? *(m_ePhotonSource.begin() + m_primaryWavelengthIndex) : SKTRAN_MCPhoton_Base::photonSource();
}

const SKTRAN_MCPhoton_RadInfo& SKTRAN_MCPhoton_Ring::photonSource(bool elasticRaman) const
{
	return elasticRaman ? *(m_ePhotonSource.begin() + m_primaryWavelengthIndex) : SKTRAN_MCPhoton_Base::photonSource();
}

std::vector<SKTRAN_MCPhoton_RadInfo>& SKTRAN_MCPhoton_Ring::photonSources(bool elasticRaman)
{
	return elasticRaman ? m_ePhotonSource : SKTRAN_MCPhoton_Base::photonSources();
}

const std::vector<SKTRAN_MCPhoton_RadInfo>& SKTRAN_MCPhoton_Ring::photonSources(bool elasticRaman) const
{
	return elasticRaman ? m_ePhotonSource : SKTRAN_MCPhoton_Base::photonSources();
}

void SKTRAN_MCPhoton_Ring::ResetRadiance()
{
	SKTRAN_MCPhoton_Base::ResetRadiance();
	for (auto&& it : m_ePhotonRadiance) it.SetVector(0.0);
	std::fill(m_eScatterWeight.begin(), m_eScatterWeight.end(), 1.0);
}

void SKTRAN_MCPhoton_Ring::ResetFactors()
{
	SKTRAN_MCPhoton_Base::ResetFactors();
	for (auto&& it : m_ePhotonSource) it.SetVector(0.0);
	std::fill(m_eWeightFactor.begin(), m_eWeightFactor.end(), 1.0);
	std::fill(m_eAlbedo.begin(), m_eAlbedo.end(), 1.0);
	std::fill(m_eTransmission.begin(), m_eTransmission.end(), 1.0);
	for (auto&& it : m_eOpticalDepth) std::fill(it.begin(), it.end(), 0.0);
}

bool SKTRAN_MCPhoton_Ring::SetWavelengths(const std::vector<double>& wavelengths)
{
	bool ok = true;

	m_finalWavelength = wavelengths; 
	m_currentWavelength = wavelengths;
	m_numWavelengths = 1;

	m_currentWavelength.resize(m_numWavelengths);
	m_finalWavelength.resize(m_numWavelengths);
	m_photonRadiance.resize(m_numWavelengths);
	m_photonSource.resize(m_numWavelengths);
	m_scatterWeight.resize(m_numWavelengths);
	m_weightFactor.resize(m_numWavelengths);
	m_albedo.resize(m_numWavelengths);
	m_transmission.resize(m_numWavelengths);
	m_scatterFactor.resize(m_numWavelengths);
	m_opticaldepth.resize(m_numWavelengths);

	m_ePhotonRadiance.resize(m_numWavelengths);
	m_ePhotonSource.resize(m_numWavelengths);
	m_eScatterWeight.resize(m_numWavelengths);
	m_eWeightFactor.resize(m_numWavelengths);
	m_eAlbedo.resize(m_numWavelengths);
	m_eTransmission.resize(m_numWavelengths);
	m_eOpticalDepth.resize(m_numWavelengths);

	ResetRadiance();
	
	return ok;
}

bool SKTRAN_MCPhoton_Ring::SetCurrentWavelength(double wavelength)
{
	bool ok = true;

	m_primaryWavelengthIndex = 0;
	m_primaryWavelength = wavelength;
	if (m_photonOptical) m_photonOptical->SetWavelength(wavelength);
	m_currentWavelength = { wavelength };
	m_finalWavelength = { wavelength };
	
	return ok;
}

bool SKTRAN_MCPhoton_Ring::TraceRays(const SKTRAN_OpticalPropertiesIntegrator_Base * integrator, bool curved)
{
	bool ok = true;

	ok = ok && m_photonOptical->TraceRay_NewMethod();

	if (curved)
	{
		std::vector<double> sigmak(m_photonOptical->GetNumQuadraturePoints()); // should allocate persistent vectors if this is ever used
		std::vector<double> sigmaf(m_photonOptical->GetNumQuadraturePoints());

		if (m_currentWavelength[0] != m_finalWavelength[0])
		{
			m_photonOptical->SetWavelength(m_finalWavelength[0]);
			ok = ok && integrator->CalculateRayScalarTransmissionVector(m_photonOptical.get(), NULL, false, true, &sigmak, &sigmaf);
			m_eOpticalDepth[0] = m_photonOptical->OpticalDepthArray();
		}
		m_photonOptical->SetWavelength(m_currentWavelength[0]);
		ok = ok && integrator->CalculateRayScalarTransmissionVector(m_photonOptical.get(), NULL, false, true, &sigmak, &sigmaf); // integrate the primary wavelength last so that photonOptical remembers the optical depth for each cell
		m_opticaldepth[0] = m_photonOptical->OpticalDepthArray();
	}
	else
	{
		if (m_currentWavelength[0] != m_finalWavelength[0])
		{
			m_photonOptical->SetWavelength(m_finalWavelength[0]);
			ok = ok && integrator->CalculateRayScalarTransmission_withMinContainer(m_photonOptical.get(), NULL, false, true);
			m_eOpticalDepth[0] = m_photonOptical->OpticalDepthArray();
		}
		m_photonOptical->SetWavelength(m_currentWavelength[0]);
		ok = ok && integrator->CalculateRayScalarTransmission_withMinContainer(m_photonOptical.get(), NULL, false, true); // integrate the primary wavelength last so that photonOptical remembers the optical depth for each cell
		m_opticaldepth[0] = m_photonOptical->OpticalDepthArray();
	}

	return ok;
}

bool SKTRAN_MCPhoton_Ring::CalculateAlbedo(const SKTRAN_TableOpticalProperties_Base * opticalpropertiestable, const SKTRAN_TableOpticalProperties_MCBase* mcOptTable, const HELIODETIC_POINT & scatterPoint)
{
	bool ok = SKTRAN_MCPhoton_Base::CalculateAlbedo(opticalpropertiestable, mcOptTable, scatterPoint);

	m_eAlbedo[0] = m_albedo[0];
	if (m_currentWavelength[0] != m_finalWavelength[0])
	{
		double kcurrent = mcOptTable->InelasticProperties()->InelasticExtinctionPerCM(m_currentWavelength[0], scatterPoint);
		double kfinal = mcOptTable->InelasticProperties()->InelasticExtinctionPerCM(m_finalWavelength[0], scatterPoint);
		m_eAlbedo[0] *= kfinal / kcurrent;
	}

	return ok;
}

bool SKTRAN_MCPhoton_Ring::UpdateScatterWeight(double forcedScatterCorrFactor)
{
	bool ok = true;

	double weightFactor;

	weightFactor = m_albedo[0];					// single scatter albedo for atmo scatter, surface albedo for ground scatter
	weightFactor *= forcedScatterCorrFactor;	// forced scatter correction factor
	weightFactor *= m_optFactor;				// optimization factor: from forcing an elastic or inelastic scatter - should be 1.0 when not in optimized mode
	m_weightFactor[0] = weightFactor;
	m_scatterWeight[0] *= weightFactor;

	weightFactor = m_eAlbedo[0];
	weightFactor *= forcedScatterCorrFactor;
	weightFactor *= m_optFactor;
	weightFactor *= m_eTransmission[0] / m_transmission[0];
	m_eWeightFactor[0] = weightFactor;
	m_eScatterWeight[0] *= weightFactor;

	return ok;
}

bool SKTRAN_MCPhoton_Ring::CalculateTransmissionsAtmoScatter(SKTRAN_TableOpticalProperties_Base const * opticalpropertiestable, SKTRAN_MCPhoton_Base const * incomingPhoton, const HELIODETIC_POINT & scatterPoint)
{
	bool ok = true;
	size_t scatterCellIndex;
	HELIODETIC_POINT scatterCellStartPoint;
	SKTRAN_OpticalDepthCalculator_LinearWithHeight odcalculator;
	const std::vector<double>& odarray = incomingPhoton->photonOptical()->OpticalDepthArray();
	const SKTRAN_RayStorage_Base* storage = incomingPhoton->photonOptical()->Storage();

	// start by finding (again) the cell containing the scatter point (should be passed from integrator instead at some point)
	ok = ok && FindScatterPointCellIndex(storage, odarray, m_targetTau, scatterCellIndex, scatterCellStartPoint);
	// now configure the integrator to integrate from the start of this cell to the scatterPoint
	ok = ok && ConfigureQuadratureCoefficients(storage, scatterPoint, scatterCellIndex, scatterCellStartPoint, odcalculator);

	// calculate and store the required transmissions
	double sigma0, sigma1;
	ok = ok && opticalpropertiestable->GetEffectiveExtinctionPerCMWithHeight1(incomingPhoton->CurrentWavelength(), scatterCellStartPoint, scatterPoint, &sigma0, &sigma1);
	m_transmission[0] = exp(-incomingPhoton->OpticalDepth(false)[scatterCellIndex] - odcalculator.OpticalDepthFromStartToEnd(sigma0, sigma1));
	if (m_currentWavelength[0] != m_finalWavelength[0])
	{
		ok = ok && opticalpropertiestable->GetEffectiveExtinctionPerCMWithHeight1(incomingPhoton->FinalWavelength(), scatterCellStartPoint, scatterPoint, &sigma0, &sigma1);
		m_eTransmission[0] = exp(-incomingPhoton->OpticalDepth(true)[scatterCellIndex] - odcalculator.OpticalDepthFromStartToEnd(sigma0, sigma1));
	}
	else
	{
		m_eTransmission[0] = m_transmission[0];
	}

	return ok;
}

bool SKTRAN_MCPhoton_Ring::CalculateTransmissionsGroundScatter(SKTRAN_MCPhoton_Base const * incomingPhoton, const HELIODETIC_POINT & scatterPoint)
{
	bool ok = true;

	m_transmission[0] = exp(-incomingPhoton->OpticalDepth(false).back());
	if (m_currentWavelength[0] != m_finalWavelength[0]) m_eTransmission[0] = exp(-incomingPhoton->OpticalDepth(true).back());
	else m_eTransmission[0] = m_transmission[0];

	return ok;
}

bool SKTRAN_MCPhoton_Ring::Configure(const SKTRAN_MCPhoton_Ring& other)
{
	bool ok = SKTRAN_MCPhoton_Base::Configure(other);

	m_eScatterWeight = other.m_eScatterWeight;
	m_eWeightFactor = other.m_eWeightFactor;
	m_eAlbedo = other.m_eAlbedo;
	m_eTransmission = other.m_eTransmission;
	m_eOpticalDepth = other.m_eOpticalDepth;

	return ok;
}

bool SKTRAN_MCPhoton_Simultaneous::TraceRays(const SKTRAN_OpticalPropertiesIntegrator_Base * integrator, bool curved)
{
	bool ok = true;

	ok = ok && m_photonOptical->TraceRay_NewMethod();

	if (curved)
	{
		std::vector<double> sigmak(m_photonOptical->GetNumQuadraturePoints()); // should allocate persistent vectors if this is ever used
		std::vector<double> sigmaf(m_photonOptical->GetNumQuadraturePoints());

		for (size_t wlidx = 0; wlidx < m_numWavelengths; wlidx++)
		{
			if (wlidx != m_primaryWavelengthIndex)
			{
				m_photonOptical->SetWavelength(m_currentWavelength[wlidx]);
				ok = ok && integrator->CalculateRayScalarTransmissionVector(m_photonOptical.get(), NULL, false, true, &sigmak, &sigmaf);
				m_opticaldepth[wlidx] = m_photonOptical->OpticalDepthArray();
			}
		}
		m_photonOptical->SetWavelength(m_primaryWavelength);
		ok = ok && integrator->CalculateRayScalarTransmissionVector(m_photonOptical.get(), NULL, false, true, &sigmak, &sigmaf); // integrate the primary wavelength last so that photonOptical remembers the optical depth for each cell
		m_opticaldepth[m_primaryWavelengthIndex] = m_photonOptical->OpticalDepthArray();
	}
	else
	{
		for (size_t wlidx = 0; wlidx < m_numWavelengths; wlidx++)
		{
			if (wlidx != m_primaryWavelengthIndex)
			{
				m_photonOptical->SetWavelength(m_currentWavelength[wlidx]);
				ok = ok && integrator->CalculateRayScalarTransmission_withMinContainer(m_photonOptical.get(), NULL, false, true);
				m_opticaldepth[wlidx] = m_photonOptical->OpticalDepthArray();
			}
		}
		m_photonOptical->SetWavelength(m_currentWavelength[m_primaryWavelengthIndex]);
		ok = ok && integrator->CalculateRayScalarTransmission_withMinContainer(m_photonOptical.get(), NULL, false, true); // integrate the primary wavelength last so that photonOptical remembers the optical depth for each cell
		m_opticaldepth[m_primaryWavelengthIndex] = m_photonOptical->OpticalDepthArray();
	}

	return ok;
}

bool SKTRAN_MCPhoton_Simultaneous::UpdateScatterWeight(double forcedScatterCorrFactor)
{
	bool ok = true;

	double weightFactor;
	for (size_t idx = 0; idx < m_numWavelengths; idx++)
	{
		weightFactor = m_albedo[idx];													// single scatter albedo for atmo scatter, surface albedo for ground scatter
		weightFactor *= forcedScatterCorrFactor;										// forced scatter correction factor
		weightFactor *= m_transmission[idx] / m_transmission[m_primaryWavelengthIndex]; // simultaneous-wavelength factor: from picking a scatter point from the transmission function
		weightFactor *= m_scatterFactor[idx];											// simultaneous-wavelength factor: from scatter direction and/or wavelength-shift distributions
		weightFactor *= m_optFactor;													// optimization factor: from forcing an elastic or inelastic scatter
		m_weightFactor[idx] = weightFactor;
		m_scatterWeight[idx] *= weightFactor;
	}

	return ok;
}

bool SKTRAN_MCPhoton_Simultaneous::CalculateTransmissionsAtmoScatter(const SKTRAN_TableOpticalProperties_Base* opticalpropertiestable, SKTRAN_MCPhoton_Base const * incomingPhoton, const HELIODETIC_POINT & scatterPoint)
{
	bool ok = true;
	size_t scatterCellIndex;
	HELIODETIC_POINT scatterCellStartPoint;
	SKTRAN_OpticalDepthCalculator_LinearWithHeight odcalculator;
	const std::vector<double>& odarray = incomingPhoton->photonOptical()->OpticalDepthArray();
	const SKTRAN_RayStorage_Base* storage = incomingPhoton->photonOptical()->Storage();

	// start by finding (again) the cell containing the scatter point (should be passed from integrator instead at some point)
	ok = ok && FindScatterPointCellIndex(storage, odarray, m_targetTau, scatterCellIndex, scatterCellStartPoint);
	// now configure the integrator to integrate from the start of this cell to the scatterPoint
	ok = ok && ConfigureQuadratureCoefficients(storage, scatterPoint, scatterCellIndex, scatterCellStartPoint, odcalculator);

	// calculate and store the required transmissions
	double sigma0, sigma1;
	auto wl = incomingPhoton->CurrentWavelengths().cbegin();
	auto od = incomingPhoton->OpticalDepths().cbegin();
	for (auto tr = m_transmission.begin(); tr != m_transmission.end(); tr++, wl++, od++)
	{
		ok = ok && opticalpropertiestable->GetEffectiveExtinctionPerCMWithHeight1(*wl, scatterCellStartPoint, scatterPoint, &sigma0, &sigma1);
		*tr = exp(-(*od)[scatterCellIndex] - odcalculator.OpticalDepthFromStartToEnd(sigma0, sigma1));
	}

	return ok;
}

bool SKTRAN_MCPhoton_Simultaneous::CalculateTransmissionsGroundScatter(SKTRAN_MCPhoton_Base const * incomingPhoton, const HELIODETIC_POINT & scatterPoint)
{
	bool ok = true;

	auto od = incomingPhoton->OpticalDepths().cbegin();
	for (auto tr = m_transmission.begin(); tr != m_transmission.end(); tr++, od++)
	{
		*tr = exp(-(*od).back());
	}
	return ok;
}

bool SKTRAN_MCPhoton_Simultaneous::SetWavelengths(const std::vector<double>& wavelengths)
{
	bool ok = true;

	m_numWavelengths = wavelengths.size();
	ok = ok && m_numWavelengths > 1;

	if (ok)
	{
		m_finalWavelength = wavelengths;
		m_currentWavelength = wavelengths;
		m_photonRadiance.resize(m_numWavelengths);
		m_photonSource.resize(m_numWavelengths);
		m_scatterWeight.resize(m_numWavelengths);
		m_weightFactor.resize(m_numWavelengths);
		m_albedo.resize(m_numWavelengths);
		m_transmission.resize(m_numWavelengths);
		m_scatterFactor.resize(m_numWavelengths);
		m_opticaldepth.resize(m_numWavelengths);

		ResetRadiance();
	}

	return ok;
}

bool SKTRAN_MCPhoton_Simultaneous::SetCurrentWavelength(double wavelength)
{
	bool ok = true;
	
	auto it = std::find(m_finalWavelength.begin(), m_finalWavelength.end(), wavelength);
	ok = ok && it != m_finalWavelength.end();
	if (ok)
	{
		m_primaryWavelengthIndex = it - m_finalWavelength.begin();
		m_primaryWavelength = wavelength;
		if (m_photonOptical) m_photonOptical->SetWavelength(wavelength);
	}
	
	return ok;
}


bool SKTRAN_MCPhoton_SimultaneousRing::CalculateTransmissionsAtmoScatter(SKTRAN_TableOpticalProperties_Base const * opticalpropertiestable, SKTRAN_MCPhoton_Base const * incomingPhoton, const HELIODETIC_POINT & scatterPoint)
{
	bool ok = true;
	size_t scatterCellIndex;
	HELIODETIC_POINT scatterCellStartPoint;
	SKTRAN_OpticalDepthCalculator_LinearWithHeight odcalculator;
	const std::vector<double>& odarray = incomingPhoton->photonOptical()->OpticalDepthArray();
	const SKTRAN_RayStorage_Base* storage = incomingPhoton->photonOptical()->Storage();

	// start by finding (again) the cell containing the scatter point (should be passed from integrator instead at some point)
	ok = ok && FindScatterPointCellIndex(storage, odarray, m_targetTau, scatterCellIndex, scatterCellStartPoint);
	// now configure the integrator to integrate from the start of this cell to the scatterPoint
	ok = ok && ConfigureQuadratureCoefficients(storage, scatterPoint, scatterCellIndex, scatterCellStartPoint, odcalculator);

	// calculate and store the required transmissions
	double sigma0, sigma1;
	for (size_t wlidx = 0; wlidx < m_numWavelengths; wlidx++)
	{
		ok = ok && opticalpropertiestable->GetEffectiveExtinctionPerCMWithHeight1(incomingPhoton->CurrentWavelengths()[wlidx], scatterCellStartPoint, scatterPoint, &sigma0, &sigma1);
		m_transmission[wlidx] = exp(-incomingPhoton->OpticalDepths(false)[wlidx][scatterCellIndex] - odcalculator.OpticalDepthFromStartToEnd(sigma0, sigma1));
		if (m_currentWavelength[wlidx] != m_finalWavelength[wlidx])
		{
			ok = ok && opticalpropertiestable->GetEffectiveExtinctionPerCMWithHeight1(incomingPhoton->FinalWavelengths()[wlidx], scatterCellStartPoint, scatterPoint, &sigma0, &sigma1);
			m_eTransmission[wlidx] = exp(-incomingPhoton->OpticalDepths(true)[wlidx][scatterCellIndex] - odcalculator.OpticalDepthFromStartToEnd(sigma0, sigma1));
		}
		else
		{
			m_eTransmission[wlidx] = m_transmission[wlidx];
		}
	}

	return ok;
}

bool SKTRAN_MCPhoton_SimultaneousRing::CalculateTransmissionsGroundScatter(SKTRAN_MCPhoton_Base const * incomingPhoton, const HELIODETIC_POINT & scatterPoint)
{
	bool ok = true;

	for (size_t wlidx = 0; wlidx < m_numWavelengths; wlidx++)
	{
		m_transmission[wlidx] = exp(-incomingPhoton->OpticalDepths(false)[wlidx].back());
		if (m_currentWavelength[wlidx] != m_finalWavelength[wlidx]) m_eTransmission[wlidx] = exp(-incomingPhoton->OpticalDepths(true)[wlidx].back());
		else m_eTransmission[wlidx] = m_transmission[wlidx];
	}

	return ok;
}

bool SKTRAN_MCPhoton_SimultaneousRing::TraceRays(const SKTRAN_OpticalPropertiesIntegrator_Base * integrator, bool curved)
{
	bool ok = true;

	ok = ok && m_photonOptical->TraceRay_NewMethod();

	if (curved)
	{
		std::vector<double> sigmak(m_photonOptical->GetNumQuadraturePoints()); // should allocate persistent vectors if this is ever used
		std::vector<double> sigmaf(m_photonOptical->GetNumQuadraturePoints());

		for (size_t wlidx = 0; wlidx < m_numWavelengths; wlidx++)
		{
			if (wlidx != m_primaryWavelengthIndex)
			{
				m_photonOptical->SetWavelength(m_currentWavelength[wlidx]);
				ok = ok && integrator->CalculateRayScalarTransmissionVector(m_photonOptical.get(), NULL, false, true, &sigmak, &sigmaf);
				m_opticaldepth[wlidx] = m_photonOptical->OpticalDepthArray();
			}
		}
		m_photonOptical->SetWavelength(m_primaryWavelength);
		ok = ok && integrator->CalculateRayScalarTransmissionVector(m_photonOptical.get(), NULL, false, true, &sigmak, &sigmaf); // integrate the primary wavelength last so that photonOptical remembers the optical depth for each cell
		m_opticaldepth[m_primaryWavelengthIndex] = m_photonOptical->OpticalDepthArray();
	}
	else
	{
		for (size_t wlidx = 0; wlidx < m_numWavelengths; wlidx++)
		{
			if (wlidx != m_primaryWavelengthIndex)
			{
				m_photonOptical->SetWavelength(m_currentWavelength[wlidx]);
				ok = ok && integrator->CalculateRayScalarTransmission_withMinContainer(m_photonOptical.get(), NULL, false, true);
				m_opticaldepth[wlidx] = m_photonOptical->OpticalDepthArray();
			}
		}
		m_photonOptical->SetWavelength(m_currentWavelength[m_primaryWavelengthIndex]);
		ok = ok && integrator->CalculateRayScalarTransmission_withMinContainer(m_photonOptical.get(), NULL, false, true); // integrate the primary wavelength last so that photonOptical remembers the optical depth for each cell
		m_opticaldepth[m_primaryWavelengthIndex] = m_photonOptical->OpticalDepthArray();
	}

	return ok;
}

bool SKTRAN_MCPhoton_SimultaneousRing::CalculateAlbedo(const SKTRAN_TableOpticalProperties_Base * opticalpropertiestable, const SKTRAN_TableOpticalProperties_MCBase * mcOptTable, const HELIODETIC_POINT & scatterPoint)
{


	bool ok = SKTRAN_MCPhoton_Base::CalculateAlbedo(opticalpropertiestable, mcOptTable, scatterPoint);

	for (size_t wlidx = 0; wlidx < m_numWavelengths; wlidx++)
	{
		m_eAlbedo[wlidx] = m_albedo[wlidx];
		if (m_currentWavelength[wlidx] != m_finalWavelength[wlidx])
		{
			double kcurrent = mcOptTable->InelasticProperties()->InelasticExtinctionPerCM(m_currentWavelength[wlidx], scatterPoint);
			double kfinal = mcOptTable->InelasticProperties()->InelasticExtinctionPerCM(m_finalWavelength[wlidx], scatterPoint);
			m_eAlbedo[wlidx] *= kfinal / kcurrent;
		}
	}
	return ok;
}

bool SKTRAN_MCPhoton_SimultaneousRing::UpdateScatterWeight(double forcedScatterCorrFactor)
{
	bool ok = true;

	double weightFactor;
	for (size_t idx = 0; idx < m_numWavelengths; idx++)
	{
		weightFactor = m_albedo[idx];													// single scatter albedo for atmo scatter, surface albedo for ground scatter
		weightFactor *= forcedScatterCorrFactor;										// forced scatter correction factor
		weightFactor *= m_transmission[idx] / m_transmission[m_primaryWavelengthIndex]; // simultaneous-wavelength factor: from scatter direction and/or wavelength-shift distributions
		weightFactor *= m_scatterFactor[idx];											// simultaneous-wavelength factor: from picking a scatter point from the transmission function
		weightFactor *= m_optFactor;													// optimization factor: from forcing an elastic or inelastic scatter - should be 1.0 when not in optimized mode
		m_weightFactor[idx] = weightFactor;
		m_scatterWeight[idx] *= weightFactor;

		weightFactor = m_eAlbedo[idx];
		weightFactor *= forcedScatterCorrFactor;
		weightFactor *= m_optFactor;
		weightFactor *= m_eTransmission[idx] / m_transmission[m_primaryWavelengthIndex];
		weightFactor *= m_scatterFactor[idx];
		m_eWeightFactor[idx] = weightFactor;
		m_eScatterWeight[idx] *= weightFactor;
	}
	return ok;
}

bool SKTRAN_MCPhoton_SimultaneousRing::SetWavelengths(const std::vector<double>& wavelengths)
{
	bool ok = true;

	m_numWavelengths = wavelengths.size();
	ok = ok && m_numWavelengths > 1;

	if (ok)
	{
		m_finalWavelength = wavelengths;
		m_currentWavelength = wavelengths;
		m_photonRadiance.resize(m_numWavelengths);
		m_photonSource.resize(m_numWavelengths);
		m_scatterWeight.resize(m_numWavelengths);
		m_weightFactor.resize(m_numWavelengths);
		m_albedo.resize(m_numWavelengths);
		m_transmission.resize(m_numWavelengths);
		m_scatterFactor.resize(m_numWavelengths);
		m_opticaldepth.resize(m_numWavelengths);

		m_ePhotonRadiance.resize(m_numWavelengths);
		m_ePhotonSource.resize(m_numWavelengths);
		m_eScatterWeight.resize(m_numWavelengths);
		m_eWeightFactor.resize(m_numWavelengths);
		m_eAlbedo.resize(m_numWavelengths);
		m_eTransmission.resize(m_numWavelengths);
		m_eOpticalDepth.resize(m_numWavelengths);

		ResetRadiance();
	}

	return ok;
}

bool SKTRAN_MCPhoton_SimultaneousRing::SetCurrentWavelength(double wavelength)
{
	bool ok = true;

	auto it = std::find(m_finalWavelength.begin(), m_finalWavelength.end(), wavelength);
	ok = ok && it != m_finalWavelength.end();
	if (ok)
	{
		m_primaryWavelengthIndex = it - m_finalWavelength.begin();
		m_primaryWavelength = wavelength;
		if (m_photonOptical) m_photonOptical->SetWavelength(wavelength);
	}
	
	return ok;
}
