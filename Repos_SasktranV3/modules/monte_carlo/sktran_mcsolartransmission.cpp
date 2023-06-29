#include "include/sktran_montecarlo_internals.h"



bool SKTRAN_SolarTransmission_NoTable_reuseRays_MC::MonteCarlo_GroundScatteredRadianceAtPoint(const SKTRAN_SourceTermQueryObject_Base & qobj, SKTRAN_MCPhoton_Base* photon) const
{
	bool ok = true;
	double scattCoeff_solar;

	// use scatter factors as temporary vectors (it is not needed at this point in the calculation)
	std::vector<double>& scattCoeffs = photon->ScatterFactors(); // this vector will store the solar scatter coefficient (BRDF evaluated at solar angles)

	Sun()->UpdateSun();
	scattCoeff_solar = Sun()->CosAngleToSun(qobj.GetPoint().UnitVector()) / nxmath::Pi;
	std::fill(scattCoeffs.begin(), scattCoeffs.end(), scattCoeff_solar);

	ok = ok && MonteCarlo_SourceAtPoint(qobj, photon, photon->CurrentWavelengths(), scattCoeffs, photon->photonSources());

	return ok;
}

bool SKTRAN_SolarTransmission_NoTable_reuseRays_MC::MonteCarlo_SingleScatteredRadianceAtPoint(const SKTRAN_SourceTermQueryObject_Base & qobj, SKTRAN_MCPhoton_Base* photon) const
{
	bool ok = true;

	// use transmissions and scatter factors as temporary vectors (they are not needed at this point in the calculation)
	std::vector<double>& wavelengths = photon->Transmissions();
	std::vector<double>& scattCoeffs = photon->ScatterFactors(); // this vector will store the solar scatter coefficient (phase function at ssa) times the scatter xs ratio

	ok = ok && ScatterCoefficient(qobj, photon, photon->CurrentWavelengths(), wavelengths, scattCoeffs);
	ok = ok && MonteCarlo_SourceAtPoint(qobj, photon, wavelengths, scattCoeffs, photon->photonSources());

	return ok;
}

bool SKTRAN_SolarTransmission_NoTable_reuseRays_MC::ScatterCoefficient(const SKTRAN_SourceTermQueryObject_Base & qobj, const SKTRAN_MCPhoton_Base* photon, const std::vector<double>& outgoingwavelengths, std::vector<double>& incomingwavelengths, std::vector<double>& scattCoeffs) const
{
	bool ok = true;

	double scattCoeff_solar, kscatt;
	kscatt = Integrator()->GetOpticalProps()->ScatteringExtinctionPerCM(outgoingwavelengths[photon->m_primaryWavelengthIndex], qobj.GetPoint()); // k_s* = scattering extinction at primary wavelength

	auto cwl = outgoingwavelengths.cbegin();
	auto nwl = incomingwavelengths.begin();
	for (auto sc = scattCoeffs.begin(); sc != scattCoeffs.end(); sc++, nwl++, cwl++)
	{
		ok = ok && Integrator()->GetOpticalProps()->GetScatteringCoefficientCM2(*cwl, qobj.GetPoint(), Sun()->CosAngleToSun(qobj.GetLookAway()), &scattCoeff_solar); // scattCoeff_solar now contains k_s * p, what we want is p * (xs ratio) = p * (k_s / k_s*), therefore divide by primary xs
		*sc = 0.0 == kscatt ? 0.0 : scattCoeff_solar / kscatt;
		*nwl = *cwl; // no wavelength change: assign current wavelengths to "new" wavelengths
	}

	return ok;
}

bool SKTRAN_SolarTransmission_NoTable_reuseRays_MC::MonteCarlo_SourceAtPoint(const SKTRAN_SourceTermQueryObject_Base & qobj, const SKTRAN_MCPhoton_Base* photon, const std::vector<double> wavelengths, const std::vector<double>& scattCoeffs, std::vector<SKTRAN_MCPhoton_RadInfo>& sources) const
{
	bool ok = true;

	double transmission;

	Sun()->UpdateSun();

	auto wl = wavelengths.cbegin();
	auto sc = scattCoeffs.cbegin();
	for (auto source = sources.begin(); source != sources.end(); source++, wl++, sc++)
	{
		ok = ok && TransmissionAtPoint(*wl, qobj.GetPoint(), transmission);
		source->AddToVector(transmission * (*sc));
	}

	return ok;
}

bool SKTRAN_SolarTransmission_NoTable_reuseRays_SolarSpectrum_MC::SetSolarSpectrum(const std::vector<double>& solarSpectrum)
{
	bool ok = true;
	m_solarSpectrum = solarSpectrum;
	return ok;
}

bool SKTRAN_SolarTransmission_NoTable_reuseRays_SolarSpectrum_MC::MonteCarlo_SourceAtPoint(const SKTRAN_SourceTermQueryObject_Base & qobj, const SKTRAN_MCPhoton_Base* photon, const std::vector<double> wavelengths, const std::vector<double>& scattCoeffs, std::vector<SKTRAN_MCPhoton_RadInfo>& sources) const
{
	bool ok = true;

	double irradiance, transmission;
	SKTRAN_GridIndex lowercell, uppercell;
	double lowerweight, upperweight;

	Sun()->UpdateSun();

	auto wl = wavelengths.cbegin();
	auto sc = scattCoeffs.cbegin();
	for (auto source = sources.begin(); source != sources.end(); source++, wl++, sc++)
	{
		ok = ok && TransmissionAtPoint(*wl, qobj.GetPoint(), transmission);
		ok = ok && m_wavelengthgrid->FindBoundingIndices(*wl, SKTRAN_GridDefBase_V2::ENUM_INTERPOLATIONMODE::OUTOFBOUND_ERROR, &lowercell, &lowerweight, &uppercell, &upperweight);
		irradiance = m_solarSpectrum.size() > 0 ? m_solarSpectrum[lowercell] * lowerweight + m_solarSpectrum[uppercell] * upperweight : 1.0;
		source->AddToVector(irradiance * transmission * (*sc) * photon->ManualScatterFactor());
	}

	return ok;
}

SKTRAN_SolarTransmission_Inelastic_MC::SKTRAN_SolarTransmission_Inelastic_MC()
{
	m_mcOptProps = nullptr;
}

SKTRAN_SolarTransmission_Inelastic_MC::~SKTRAN_SolarTransmission_Inelastic_MC()
{
}

//bool SKTRAN_SolarTransmission_Inelastic_MC::MonteCarlo_SingleScatteredRadianceAtPoint(const double& wavelength, const SKTRAN_SourceTermQueryObject_Base & qobj, double & radiance) const
//{
//	bool ok = true;
//
//	size_t threadid = omp_get_thread_num();
//	double incomingwavelength;
//	double r = (*m_rng)[threadid]();
//
//	double transmission;
//	double scattCoeff_solar;
//	double kscatt;
//
//	double kelastic = Integrator()->GetOpticalProps()->ScatteringExtinctionPerCM(wavelength, qobj.GetPoint());
//	double kinelastic = m_mcOptProps->InelasticProperties()->InelasticExtinctionPerCM(wavelength, qobj.GetPoint());
//	double k = r * (kelastic + kinelastic);
//	if (k < kinelastic)
//	{
//		r = k / kinelastic;
//		m_mcOptProps->InelasticProperties()->GetIncomingWavelength(wavelength, qobj.GetPoint(), r, incomingwavelength);
//		ok = ok && m_mcOptProps->InelasticProperties()->GetScatteringCoefficientCM2(incomingwavelength, qobj.GetPoint(), Sun()->CosAngleToSun(qobj.GetLookAway()), &scattCoeff_solar);	// New look is heliodetic --> cosScatAngle = dot(-zhat,newLook) = -newLook.z
//		kscatt = m_mcOptProps->InelasticProperties()->InelasticExtinctionPerCM(incomingwavelength, qobj.GetPoint());
//	}
//	else
//	{
//		incomingwavelength = wavelength;
//		ok = ok && Integrator()->GetOpticalProps()->GetScatteringCoefficientCM2(incomingwavelength, qobj.GetPoint(), Sun()->CosAngleToSun(qobj.GetLookAway()), &scattCoeff_solar);
//		kscatt = Integrator()->GetOpticalProps()->ScatteringExtinctionPerCM(incomingwavelength, qobj.GetPoint());
//	}
//
//	Sun()->UpdateSun();
//	ok = ok && TransmissionAtPoint(incomingwavelength, qobj.GetPoint(), transmission);
//	scattCoeff_solar = 0.0 == kscatt ? 0.0 : scattCoeff_solar / kscatt;			// GetScatteringCoefficientCM2 returns P_normalized*kscatt
//	radiance = transmission * scattCoeff_solar;
//
//	SKTRAN_GridIndex lowercell, uppercell;
//	double lowerweight, upperweight;
//	ok = ok && m_wavelengthgrid->FindBoundingIndices(wavelength, SKTRAN_GridDefBase_V2::ENUM_INTERPOLATIONMODE::OUTOFBOUND_ERROR, &lowercell, &lowerweight, &uppercell, &upperweight);
//	radiance *= m_solarSpectrum[lowercell] * lowerweight + m_solarSpectrum[uppercell] * upperweight;
//
//	return ok;
//}
//
//bool SKTRAN_SolarTransmission_Inelastic_MC::MonteCarlo_GroundScatteredRadianceAtPoint(const double& wavelength, const SKTRAN_SourceTermQueryObject_Base & qobj, double & radiance) const
//{
//	bool ok = true;
//
//	ok = ok && SKTRAN_SolarTransmission_NoTable::MonteCarlo_GroundScatteredRadianceAtPoint(wavelength, qobj, radiance);
//
//	SKTRAN_GridIndex lowercell, uppercell;
//	double lowerweight, upperweight;
//	ok = ok && m_wavelengthgrid->FindBoundingIndices(wavelength, SKTRAN_GridDefBase_V2::ENUM_INTERPOLATIONMODE::OUTOFBOUND_ERROR, &lowercell, &lowerweight, &uppercell, &upperweight);
//	radiance *= m_solarSpectrum[lowercell] * lowerweight + m_solarSpectrum[uppercell] * upperweight;
//
//	return ok;
//}

bool SKTRAN_SolarTransmission_Inelastic_MC::MonteCarlo_SingleScatteredRadianceAtPoint(const SKTRAN_SourceTermQueryObject_Base & qobj, SKTRAN_MCPhoton_Base* photon) const
{
	bool ok = true;
	
	size_t threadid = omp_get_thread_num();
	double r = 0.0;
	bool elastic = true;

	const double wavelength = photon->CurrentWavelength();
	const double kelastic = Integrator()->GetOpticalProps()->ScatteringExtinctionPerCM(wavelength, qobj.GetPoint());
	const double kinelastic = m_mcOptProps->InelasticProperties()->InelasticExtinctionPerCM(wavelength, qobj.GetPoint());

	if (photon->m_manualScatter)
	{
		elastic = photon->m_elasticScatter;
		r = photon->m_randNum;
		photon->ManualScatterFactor() = elastic ? kelastic / (kelastic + kinelastic) : kinelastic / (kelastic + kinelastic);
	}
	else
	{
		double k = photon->m_randNum * (kelastic + kinelastic);
		elastic = k > kinelastic;
		if (!elastic) r = k / kinelastic;
	}
	
	// use transmissions and scatter factors as temporary vectors (they are not needed at this point in the calculation)
	std::vector<double>& wavelengths = photon->Transmissions();
	std::vector<double>& scattCoeffs = photon->ScatterFactors(); // this vector will store the solar scatter coefficient (phase function at ssa) times the scatter xs ratio

	if (elastic)
	{		
		ok = ok && ScatterCoefficient(qobj, photon, photon->CurrentWavelengths(), wavelengths, scattCoeffs);
	}
	else
	{
		ok = ok && InelasticScatterCoefficient(qobj, photon, r, photon->CurrentWavelengths(), wavelengths, scattCoeffs);
	}
	
	ok = ok && MonteCarlo_SourceAtPoint(qobj, photon, wavelengths, scattCoeffs, photon->photonSources());

	return ok;
}

bool SKTRAN_SolarTransmission_Inelastic_MC::FillTable_ClassSpecific()
{
	bool ok = true;

	// no table to fill but this is a good time to make sure we have a mc optical property object
	m_mcOptProps = dynamic_cast<const SKTRAN_TableOpticalProperties_MCBase*>(Integrator()->GetOpticalProps());
	ok = ok && m_mcOptProps != nullptr;

	// and that the inelastic optical properties are defined
	ok = ok && m_mcOptProps->InelasticProperties() != nullptr;

	// and that the solar spectrum has been set
	ok = ok && m_wavelengthgrid != nullptr;
	ok = ok && m_wavelengthgrid->NumWavelengths() == m_solarSpectrum.size();

	return ok;
}

bool SKTRAN_SolarTransmission_Inelastic_MC::InelasticScatterCoefficient(const SKTRAN_SourceTermQueryObject_Base & qobj, const SKTRAN_MCPhoton_Base* photon, const double& randNum, const std::vector<double>& outgoingwavelengths, std::vector<double>& incomingwavelengths, std::vector<double>& scattCoeffs) const
{
	bool ok = true;

	double scattCoeff_solar, kscatt;

 	m_mcOptProps->InelasticProperties()->GetIncomingWavelength(outgoingwavelengths, photon->m_primaryWavelengthIndex, qobj.GetPoint(), randNum, incomingwavelengths, scattCoeffs); // choose new wavelengths, calculate xs ratios

	// Raman phase function is always wavelength-independent
	ok = ok && m_mcOptProps->InelasticProperties()->GetScatteringCoefficientCM2(photon->CurrentWavelength(), qobj.GetPoint(), Sun()->CosAngleToSun(qobj.GetLookAway()), &scattCoeff_solar);	// New look is heliodetic --> cosScatAngle = dot(-zhat,newLook) = -newLook.z
	kscatt = m_mcOptProps->InelasticProperties()->InelasticExtinctionPerCM(photon->CurrentWavelength(), qobj.GetPoint());
	scattCoeff_solar = 0.0 == kscatt ? 0.0 : scattCoeff_solar / kscatt;	// GetScatteringCoefficientCM2 returns P_normalized*kscatt | scattCoeff_solar should be wavelength-independent
	 
	for (auto&& sc : scattCoeffs) sc *= scattCoeff_solar;

	return ok;
}

bool SKTRAN_SolarTransmission_Ring_MC::MonteCarlo_SingleScatteredRadianceAtPoint(const SKTRAN_SourceTermQueryObject_Base & qobj, SKTRAN_MCPhoton_Base * photon) const
{
	bool ok = true;

	size_t threadid = omp_get_thread_num();
	double r = 0.0;
	bool elastic = true;

	const double wavelength = photon->CurrentWavelength();
	const double kelastic = Integrator()->GetOpticalProps()->ScatteringExtinctionPerCM(wavelength, qobj.GetPoint());
	const double kinelastic = m_mcOptProps->InelasticProperties()->InelasticExtinctionPerCM(wavelength, qobj.GetPoint());

	if (photon->m_manualScatter)
	{
		elastic = photon->m_elasticScatter;
		r = photon->m_randNum;
		photon->ManualScatterFactor() = elastic ? kelastic / (kelastic + kinelastic) : kinelastic / (kelastic + kinelastic);
	}
	else
	{
		double k = photon->m_randNum * (kelastic + kinelastic);
		elastic = k > kinelastic;
		if (!elastic) r = k / kinelastic;
	}

	// use transmissions and scatter factors as temporary vectors (they are not needed at this point in the calculation)
	std::vector<double>& wavelengths = photon->Transmissions();
	std::vector<double>& scattCoeffs = photon->ScatterFactors(); // this vector will store the solar scatter coefficient (phase function at ssa) times the scatter xs ratio

	if (elastic)
	{
		ok = ok && ScatterCoefficient(qobj, photon, photon->CurrentWavelengths() , wavelengths, scattCoeffs);
		ok = ok && MonteCarlo_SourceAtPoint(qobj, photon, wavelengths, scattCoeffs, photon->photonSources(false));
		if (photon->CurrentWavelength() != photon->FinalWavelength()) // calculate equivalent source term at final wavelength, if different from current wavelength
		{
			ok = ok && ScatterCoefficient(qobj, photon, photon->FinalWavelengths(), wavelengths, scattCoeffs);
			ok = ok && MonteCarlo_SourceAtPoint(qobj, photon, wavelengths, scattCoeffs, photon->photonSources(true));
		}
		else // if final and current wavelengths are the same, store the same source term
		{
			auto ps = photon->photonSources(false).cbegin();
			for (auto eps = photon->photonSources(true).begin(); eps != photon->photonSources(true).end(); ps++, eps++) eps->AddToVector(ps->GetRecentContribVec());
		}

	}
	else
	{
		ok = ok && InelasticScatterCoefficient(qobj, photon, r, photon->CurrentWavelengths(), wavelengths, scattCoeffs);
		ok = ok && MonteCarlo_SourceAtPoint(qobj, photon, wavelengths, scattCoeffs, photon->photonSources(false));
		ok = ok && InelasticScatterCoefficient(qobj, photon, r, photon->FinalWavelengths(), wavelengths, scattCoeffs);
		ok = ok && MonteCarlo_SourceAtPoint(qobj, photon, photon->FinalWavelengths(), scattCoeffs, photon->photonSources(true));
	}

	return ok;
}

bool SKTRAN_SolarTransmission_Ring_MC::MonteCarlo_GroundScatteredRadianceAtPoint(const SKTRAN_SourceTermQueryObject_Base & qobj, SKTRAN_MCPhoton_Base* photon) const
{
	bool ok = true;

	double scattCoeff_solar;

	// use scatter factors as temporary vectors (it is not needed at this point in the calculation)
	std::vector<double>& scattCoeffs = photon->ScatterFactors(); // this vector will store the solar scatter coefficient (BRDF evaluated at solar angles)

	Sun()->UpdateSun();
	scattCoeff_solar = Sun()->CosAngleToSun(qobj.GetPoint().UnitVector()) / nxmath::Pi;
	std::fill(scattCoeffs.begin(), scattCoeffs.end(), scattCoeff_solar);

	ok = ok && MonteCarlo_SourceAtPoint(qobj, photon, photon->CurrentWavelengths(), scattCoeffs, photon->photonSources(false));

	if (photon->CurrentWavelength() != photon->FinalWavelength())
	{
		ok = ok && MonteCarlo_SourceAtPoint(qobj, photon, photon->FinalWavelengths(), scattCoeffs, photon->photonSources(true));
	}
	else
	{
		auto source = photon->photonSources(false).cbegin();
		for (auto&& esource : photon->photonSources(true))
		{
			esource.AddToVector((*source++).GetRecentContribVec());
		}
	}

	return ok;
}