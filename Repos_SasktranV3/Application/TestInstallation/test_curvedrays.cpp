#include <sasktran.h>
#include "include/test_curvedrays.h"

bool TestCurvedRays()
{
	SKTRAN_CurvedRayTests tester(true, 1e-6);

	tester.RunAllTests();

	return true;
}

SKTRAN_CurvedRayTests::SKTRAN_CurvedRayTests(bool verbose, double tolerance)
{
	m_hardcodefolder = "../../../../Application/TestInstallation/hardcodetestvalues/";
	m_verbose = verbose;
	m_tol = tolerance;
}

void SKTRAN_CurvedRayTests::MakeAtmosphericState()
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

void SKTRAN_CurvedRayTests::MakeLinesOfSight()
{
	double	mjd = 54832.0;
	double  lat = -57.5;
	double	lng = 70;
	double 	sza = 60;
	double  saa = 157.5;
	double  rayazi = 0;
	double tanheights_meters[40];
	const int numHeights = 40;

	double deltaHeight = 100.0;
	double lowestHeight = 500.0;

	m_linesofsight = std::unique_ptr<SKTRAN_LineOfSightArray_V21>(new SKTRAN_LineOfSightArray_V21);

	// make the lines of sight
	for (int i = 0; i < numHeights; i++)
	{
		tanheights_meters[i] = i * deltaHeight + lowestHeight;
	}
	m_linesofsight->SetRaysFromTangentHeightArray(mjd, lat, lng, sza, saa, rayazi, tanheights_meters, numHeights, 600000.0, &m_sun);
}

void SKTRAN_CurvedRayTests::MakeObserverInsideLinesOfSight()
{
	m_linesofsight = std::unique_ptr<SKTRAN_LineOfSightArray_V21>(new SKTRAN_LineOfSightArray_V21);

	// Observer inside, looking down at ground ray
	nxVector observer(6398140.0, 0.0, -801444.94508357);
	nxVector look(-0.02235022, 0.0, 0.9997502);
	double mjd = 54372.0;

	m_linesofsight->AddLineOfSight(observer, look, mjd);

	// Observer inside, looking through limb ray
	look = nxVector(-0.01417615, 0.0, 0.99989951);
	m_linesofsight->AddLineOfSight(observer, look, mjd);

	// Observer inside, looking up ray
	look = nxVector(0.06830655, 0.0, 0.99766438);
	m_linesofsight->AddLineOfSight(observer, look, mjd);


	m_sun = nxVector(6.12323400e-17, 8.66025404e-01, 5.00000000e-01);

}

bool SKTRAN_CurvedRayTests::RunCurrentTest()
{
	nx2dArray<double> radiance;
	std::vector<SKTRAN_StokesScalar>	radiance_temp;
	radiance.SetSize(m_linesofsight->NumRays(), m_wavel.size());

	m_engine->SetSun(m_sun);
	m_engine->ConfigureModel(*m_specs, *m_linesofsight, 0);
	std::cout << "Running test " << m_testname << ":";
	for (size_t i = 0; i < m_wavel.size(); i++)
	{
		m_engine->CalculateRadiance(&radiance_temp, m_wavel[i], m_scatterorder, m_opticalstate.get());
		for (size_t j = 0; j < radiance_temp.size(); j++)
		{
			radiance.At(j, i) = radiance_temp[j];
		}
	}
	// uncomment to remake the hardcoded values
	// radiance.WriteToTextFile( (m_hardcodefolder + m_testname).c_str(), true, "%.15e");
	for (size_t i = 0; i < m_wavel.size(); i++)
	{
		for (size_t j = 0; j < radiance_temp.size(); j++)
		{
			double percdiff;
			percdiff = (m_hardcodevalues.At(j, i) - radiance.At(j, i)) / m_hardcodevalues.At(j, i) * 100;
			if (abs(percdiff) > m_tol)
			{
				std::cout << " Failed" << std::endl;
				DisplayFullTestInfo(radiance);
				return false;
			}
		}
	}
	std::cout << " Passed" << std::endl;
	if (m_verbose)
		DisplayFullTestInfo(radiance);
	return true;
}

void SKTRAN_CurvedRayTests::SetupObserverOutsideTest()
{
	m_engine = std::unique_ptr<SKTRAN_HR_Engine>(new SKTRAN_HR_Engine);
	m_specs = std::unique_ptr<SKTRAN_HR_Specs_User>(new SKTRAN_HR_Specs_User);

	m_specs->RayTracingSpecs().SetLinesOfSightType(SKTRAN_HR_RayTracer_Type::SKTRAN_HR_RayTracer_Curved);
	//m_specs->RayTracingSpecs().SetLinesOfSightType(SKTRAN_HR_RayTracer_Type::SKTRAN_HR_RayTracer_Shells);
	m_specs->IntegratorSpecs().SetMaxOpticalDepth(999999);

	MakeAtmosphericState();
	MakeLinesOfSight();
	m_scatterorder = 1;
	m_testname = "hr_curved_observeroutside";

	m_hardcodevalues.SetSize(4, 2);
	m_hardcodevalues.InputColumnMajorText((m_hardcodefolder + m_testname).c_str());
	m_wavel.clear();
	m_wavel.push_back(340);
	m_wavel.push_back(600);

}

void SKTRAN_CurvedRayTests::SetupObserverInsideTest()
{
	m_engine = std::unique_ptr<SKTRAN_HR_Engine>(new SKTRAN_HR_Engine);
	m_specs = std::unique_ptr<SKTRAN_HR_Specs_User>(new SKTRAN_HR_Specs_User);

	m_specs->RayTracingSpecs().SetLinesOfSightType(SKTRAN_HR_RayTracer_Type::SKTRAN_HR_RayTracer_Curved);
	//m_specs->RayTracingSpecs().SetLinesOfSightType(SKTRAN_HR_RayTracer_Type::SKTRAN_HR_RayTracer_Shells);
	m_specs->IntegratorSpecs().SetMaxOpticalDepth(999999);

	MakeAtmosphericState();
	MakeObserverInsideLinesOfSight();
	m_scatterorder = 1;
	m_testname = "hr_curved_observerinside";

	m_hardcodevalues.SetSize(4, 2);
	m_hardcodevalues.InputColumnMajorText((m_hardcodefolder + m_testname).c_str());
	m_wavel.clear();
	m_wavel.push_back(340);
	m_wavel.push_back(600);

}

void SKTRAN_CurvedRayTests::SetupSingleScatterTest()
{
	m_engine = std::unique_ptr<SKTRAN_HR_Engine>(new SKTRAN_HR_Engine);
	m_specs = std::unique_ptr<SKTRAN_HR_Specs_User>(new SKTRAN_HR_Specs_User);

	m_specs->RayTracingSpecs().SetLinesOfSightType(SKTRAN_HR_RayTracer_Type::SKTRAN_HR_RayTracer_Curved);


	MakeAtmosphericState();
	MakeLinesOfSight();
	m_scatterorder = 1;
	m_testname = "hr_singlescatter";

	m_wavel.clear();
	for (double wav = 340; wav < 800; wav += 5)
		m_wavel.push_back(wav);
	m_hardcodevalues.SetSize(4, m_wavel.size());
	m_hardcodevalues.InputColumnMajorText((m_hardcodefolder + m_testname).c_str());
}

void SKTRAN_CurvedRayTests::DisplayFullTestInfo(const nx2dArray<double>& radiance)
{
	double percdiff;
	if (m_verbose)
		printf("----------------------- Values -------------------------\n");
	else
		printf("------------------- Failed Values ----------------------\n");
	for (size_t waveidx = 0; waveidx < m_wavel.size(); waveidx++)
	{
		for (size_t altidx = 0; altidx < m_linesofsight->NumRays(); altidx++)
		{
			percdiff = (radiance.At(altidx, waveidx) - m_hardcodevalues.At(altidx, waveidx)) / m_hardcodevalues.At(altidx, waveidx) * 100;
			if (abs(percdiff) > m_tol || m_verbose)
				printf("Wavelength %f nm, losidx %i: %19.15f %%\n", m_wavel[waveidx], (unsigned int)altidx, percdiff);
		}
	}
	printf("--------------------------------------------------------\n");
}

void SKTRAN_CurvedRayTests::RunShortTests()
{
	SetupObserverOutsideTest();
	RunCurrentTest();
}

void SKTRAN_CurvedRayTests::RunAllTests()
{
	SetupObserverOutsideTest();
	RunCurrentTest();

	SetupObserverInsideTest();
	RunCurrentTest();
}