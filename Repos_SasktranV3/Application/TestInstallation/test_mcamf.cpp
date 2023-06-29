#include <sasktran.h>
#include "include/test_mcamf.h"
#include <float.h>

SKTRAN_MCAMF_Test::SKTRAN_MCAMF_Test(bool verbose, double tolerance)
{
	m_hardcodefolder = "../../../../Application/TestInstallation/hardcodetestvalues/";
	m_verbose = verbose;
	m_tol = tolerance;
}

void SKTRAN_MCAMF_Test::MakeAtmosphericState()
{
	skOpticalProperties_RayleighDryAir*		rayleigh;					// Optical properties of one air molecule
	skClimatology_MSIS90*					msis90;
	skClimatology_LabowOzoneVMR*			o3numberdensity;
	skOpticalProperties_O3_OSIRISRes*		o3_opticalprops;			// optical properties of one O3 molecule
	bool									ok = true;

	m_opticalstate = std::unique_ptr<SKTRAN_AtmosphericOpticalState_V21>(new SKTRAN_AtmosphericOpticalState_V21);

	rayleigh = new skOpticalProperties_RayleighDryAir;
	msis90 = new skClimatology_MSIS90;
	o3numberdensity = new skClimatology_LabowOzoneVMR;
	o3_opticalprops = new skOpticalProperties_O3_OSIRISRes;

	ok = (rayleigh != NULL) && (msis90 != NULL) && (o3_opticalprops != NULL) && (o3numberdensity != NULL);

	ok = ok && m_opticalstate->AddSpecies(SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, msis90, rayleigh);
	ok = ok && m_opticalstate->AddSpecies(SKCLIMATOLOGY_O3_CM3, o3numberdensity, o3_opticalprops);
	ok = ok && m_opticalstate->SetAlbedo(1.0);

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "MakeSpeciesList, there was an error making the species list");
	}
}

void SKTRAN_MCAMF_Test::MakeLinesOfSight()
{
	double	mjd = 58330.75;

	nxGeodetic observer = nxGeodetic();
	nxGeodetic ground = nxGeodetic();
	nxVector look;

	observer.FromGeodetic(0., 260., 35.786e6); // geostationary over 100W - approximate TEMPO location
	ground.FromGeodetic(50., 260., 0.); // approximate center of TEMPO's Canadian field of regard
	look = (ground.Location() - observer.Location()).UnitVector();

	m_linesofsight = std::unique_ptr<SKTRAN_LineOfSightArray_V21>(new SKTRAN_LineOfSightArray_V21);
	m_linesofsight->AddLineOfSight(observer.Location(), look, mjd);
}

bool SKTRAN_MCAMF_Test::RunCurrentTest()
{
	nx1dArray<double> radiance;
	nx1dArray<double> variance;
	nx2dArray<double> airmassfactor;
	nx2dArray<double> airmassfactorvariance;

	std::vector<SKTRAN_StokesScalar>	radiance_temp;
	std::vector<double>					amf_temp;
	std::vector<double>					amf_var_temp;
	
	radiance.SetSize(m_linesofsight->NumRays());
	variance.SetSize(m_linesofsight->NumRays());
	airmassfactor.SetSize(m_linesofsight->NumRays(), m_amfheights.size() > 0 ? m_amfheights.size() - 1 : 0);
	airmassfactorvariance.SetSize(m_linesofsight->NumRays(), m_amfheights.size() > 0 ? m_amfheights.size() - 1 : 0);

	//m_specs->SetSunGeographicPosition(m_sun);
	m_engine->ConfigureModel(*m_specs, *m_linesofsight, 1);
	std::cout << "Running test " << m_testname << ":";

	m_engine->CalculateRadiance(&radiance_temp, m_wavel[0], m_scatterorder, m_opticalstate.get());
	for (size_t j = 0; j < radiance_temp.size(); j++)
	{
		radiance.At(j) = radiance_temp[j];
		m_engine->GetMeasurementVariance(j, variance.At(j));
		m_engine->GetAirMassFactors(j, amf_temp);
		m_engine->GetAirMassFactorVariance(j, amf_var_temp);
		for (size_t k = 0; k < amf_temp.size(); k++)
		{
			airmassfactor.At(j, k)         = amf_temp[k];
			airmassfactorvariance.At(j, k) = amf_var_temp[k];
		}
	}
	
	bool writeVals = false;
	if (writeVals) {
		radiance.WriteToTextFile((m_hardcodefolder + m_testname).c_str(), true, "%1.16e");
		variance.WriteToTextFile((m_hardcodefolder + m_testname + "_var").c_str(), true, "%1.16e");
		airmassfactor.WriteToTextFile((m_hardcodefolder + m_testname + "_amf").c_str(), true, "%1.16e");
		airmassfactorvariance.WriteToTextFile((m_hardcodefolder + m_testname + "_amf_var").c_str(), true, "%1.16e");
	}

	if (CheckSuccess(radiance, variance, airmassfactor, airmassfactorvariance))
	{
		std::cout << " Passed" << std::endl;
		if (m_verbose)
			DisplayFullTestInfo(radiance, variance, airmassfactor, airmassfactorvariance);
		return true;
	}
	else
	{
		std::cout << " Failed" << std::endl;
		DisplayFullTestInfo(radiance, variance, airmassfactor, airmassfactorvariance);
		return false;
	}
}



void SKTRAN_MCAMF_Test::SetUpDefault()
{
	m_engine = std::unique_ptr<SKTRAN_Engine_MC_V21>(new SKTRAN_Engine_MC_V21);
	m_specs = std::unique_ptr<SKTRAN_Specifications_MC>(new SKTRAN_Specifications_MC);

	MakeAtmosphericState();
	MakeLinesOfSight();

	m_specs->ConfigureDefaults();
	m_specs->SetSunGeographicPosition(m_sun);
	m_specs->SetNumPhotonsPerLOS(500);
	m_specs->SetPrecisionMC(0.0);
	m_specs->SetMinimumRelPathWeight(0.0);
	m_specs->SetMinFractionHigherOrder(1.0); // optimization turned off for AMF calculation for now 
	m_specs->SetNumRayTracingAlts(100 + 1);
	m_specs->SetLOSRayTracerType(SKTRAN_Specifications_MC::RayTracerType::shell);
	m_specs->SetMSRayTracerType(SKTRAN_Specifications_MC::RayTracerType::shell);
	m_specs->SetSolarRayTracerType(SKTRAN_Specifications_MC::RayTracerType::shell);
	m_specs->SetSunType(SKTRAN_Specifications_MC::SunType::point);
	m_specs->SetSolarTableType(SKTRAN_Specifications_MC::SolarTableType::noTable);
	m_specs->SetOptPropIntType(SKTRAN_Specifications_MC::OptPropIntType::straight);
	m_specs->SetSunType(SKTRAN_Specifications_MC::SunType::point);
	m_specs->SetOptTableType(SKTRAN_Specifications_MC::OptTableType::dim1);
	m_specs->SetPolType(SKTRAN_Specifications_MC::PolType::none);

	m_specs->SetAllowDynamicThreads(false);
	//m_specs->SetChunkSize             ( 0 );
	m_specs->SetRngSeed(1234);
	m_specs->SetScatterPositionRes(50.0);

	m_specs->SetSecondaryOutput(SKTRAN_Specifications_MC::SecondaryOutput::lengthAMF);

	m_amfheights.resize(11);
	for (size_t i = 0; i < 11; i++)
		m_amfheights[i] = 5e3 * i;
	m_specs->SetManualAMFHeights(m_amfheights);

	m_scatterorder = 50;

	m_wavel.clear();
	m_wavel.push_back(440.);
}

void SKTRAN_MCAMF_Test::SetUpShortTest()
{
	m_testname = "MCAMF_ShortTest";

	SetUpDefault();
	m_specs->FinalizeSpecs();

	LoadHardcodedValues(1, 10);
}

void SKTRAN_MCAMF_Test::SetUpLengthTest()
{
	m_testname = "MCAMF_LengthTest";

	SetUpDefault();
	m_specs->SetNumPhotonsPerLOS(5000);
	m_amfheights.resize(51);
	for (size_t i = 0; i < 51; i++)
		m_amfheights[i] = 1e3 * i;
	m_specs->SetManualAMFHeights(m_amfheights);
	m_specs->SetSecondaryOutput(SKTRAN_Specifications_MC::SecondaryOutput::lengthAMF);
	m_specs->FinalizeSpecs();

	LoadHardcodedValues(1, 50);
}

void SKTRAN_MCAMF_Test::SetUpOpticalDepthTest()
{
	m_testname = "MCAMF_OpticalDepthTest";

	SetUpDefault();
	m_specs->SetNumPhotonsPerLOS(5000);
	m_amfheights.resize(51);
	for (size_t i = 0; i < 51; i++)
		m_amfheights[i] = 1e3 * i;
	m_specs->SetManualAMFHeights(m_amfheights);
	m_specs->SetSecondaryOutput(SKTRAN_Specifications_MC::SecondaryOutput::opticalDepthAMF);
	m_specs->SetAMFSpeciesHandle(SKCLIMATOLOGY_O3_CM3);
	m_specs->FinalizeSpecs();

	LoadHardcodedValues(1, 50);
}


void SKTRAN_MCAMF_Test::LoadHardcodedValues(int numLOS, int numAMF)
{
	bool ok;

	m_hardcoderadiance.SetSize(numLOS);
	ok = m_hardcoderadiance.InputColumnMajorText((m_hardcodefolder + m_testname).c_str());
	if (!ok) { m_hardcoderadiance.SetSize(numLOS); m_hardcoderadiance.SetTo(std::numeric_limits<double>::infinity()); }
	m_hardcodevariance.SetSize(numLOS);
	ok = m_hardcodevariance.InputColumnMajorText((m_hardcodefolder + m_testname + "_var").c_str());
	if (!ok) { m_hardcodevariance.SetSize(numLOS); m_hardcodevariance.SetTo(std::numeric_limits<double>::infinity()); }
	ok = m_hardcodeamf.InputColumnMajorText((m_hardcodefolder + m_testname + "_amf").c_str());
	if (!ok) { m_hardcodeamf.SetSize(numLOS, numAMF); m_hardcodeamf.SetTo(std::numeric_limits<double>::infinity()); }
	ok = m_hardcodeamfvariance.InputColumnMajorText((m_hardcodefolder + m_testname + "_amf_var").c_str());
	if (!ok) { m_hardcodeamfvariance.SetSize(numLOS, numAMF); m_hardcodeamfvariance.SetTo(std::numeric_limits<double>::infinity()); }

}


void SKTRAN_MCAMF_Test::DisplayFullTestInfo(const nx1dArray<double>& radiance, const nx1dArray<double>& variance, const nx2dArray<double>& airmassfactor, const nx2dArray<double>& airmassfactorvariance)
{
	double raddiff, vardiff, amfdiff, amfvardiff;

	if (m_verbose)
		printf("------------------------------- Values ---------------------------------\n");
	else
		printf("--------------------------- Failed Values ------------------------------\n");
	
	for (size_t losidx = 0; losidx < m_linesofsight->NumRays(); losidx++)
	{
		raddiff = (radiance.At(losidx) - m_hardcoderadiance.At(losidx)) / m_hardcoderadiance.At(losidx) * 100;
		vardiff = (variance.At(losidx) - m_hardcodevariance.At(losidx)) / m_hardcodevariance.At(losidx) * 100;
		if (!(raddiff == raddiff && vardiff == vardiff) || abs(raddiff) > m_tol || abs(vardiff) > m_tol || m_verbose)
		{
			printf("Wavelength %f nm, losidx %i:\n    r:%22.15e v:%22.15e (%%)\n", m_wavel[0], (unsigned int)losidx, raddiff, vardiff);
		}
		for (size_t amfidx = 0; amfidx < m_amfheights.size() - 1; amfidx++)
		{
			amfdiff = (airmassfactor.At(losidx, amfidx) / m_hardcodeamf.At(losidx, amfidx) - 1) * 100;
			amfvardiff = (airmassfactorvariance.At(losidx, amfidx) / m_hardcodeamfvariance.At(losidx, amfidx) - 1) * 100;
			if (!(amfdiff == amfdiff && amfvardiff == amfvardiff) || abs(amfdiff) > m_tol || abs(amfvardiff) > m_tol || m_verbose)
			{
				printf("  amfidx %i:\n    a:%22.15e v:%22.15e (%%)\n", (unsigned int)amfidx, amfdiff, amfvardiff);
			}
		}
	}
	
	printf("------------------------------------------------------------------------\n");
}

bool SKTRAN_MCAMF_Test::CheckSuccess(const nx1dArray<double>& radiance, const nx1dArray<double>& variance, const nx2dArray<double>& airmassfactor, const nx2dArray<double>& airmassfactorvariance) 
{
	double raddiff, vardiff, amfdiff, amfvardiff;
	for (size_t losidx = 0; losidx < m_linesofsight->NumRays(); losidx++) {
		raddiff = (radiance.At(losidx) - m_hardcoderadiance.At(losidx)) / m_hardcoderadiance.At(losidx) * 100;
		vardiff = (variance.At(losidx) - m_hardcodevariance.At(losidx)) / m_hardcodevariance.At(losidx) * 100;
		if (abs(raddiff) > m_tol || abs(vardiff) > m_tol)
			return false;
		for (size_t amfidx = 0; amfidx < m_amfheights.size() - 1; amfidx++)
		{
			amfdiff = (airmassfactor.At(losidx, amfidx) / m_hardcodeamf.At(losidx, amfidx) - 1) * 100;
			amfvardiff = (airmassfactorvariance.At(losidx, amfidx) / m_hardcodeamfvariance.At(losidx, amfidx) - 1) * 100;
			if (abs(amfdiff) > m_tol || abs(amfvardiff) > m_tol)
				return false;
		}
	}
	return true;

}

bool SKTRAN_MCAMF_Test::RunShortTests()
{
	bool ok1;

	SetUpShortTest();
	ok1 = RunCurrentTest();

	return ok1;
}

bool SKTRAN_MCAMF_Test::RunAllTests()
{
	bool ok1 = true, ok2 = true, ok3 = true;

	SetUpShortTest();
	ok1 = RunCurrentTest();

	SetUpLengthTest();
	ok2 = RunCurrentTest();

	SetUpOpticalDepthTest();
	ok3 = RunCurrentTest();

	return ok1 && ok2 && ok3;
}