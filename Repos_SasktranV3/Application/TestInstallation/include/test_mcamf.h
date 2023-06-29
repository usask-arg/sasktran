#pragma once

#include <sasktran.h>
#include <memory>

class SKTRAN_MCAMF_Test
{
private:
	std::unique_ptr<SKTRAN_Engine_MC_V21>					m_engine;
	std::unique_ptr<SKTRAN_AtmosphericOpticalState_V21>		m_opticalstate;
	std::unique_ptr<SKTRAN_LineOfSightArray_V21>			m_linesofsight;
	std::unique_ptr<SKTRAN_Specifications_MC>				m_specs;
	std::vector<double>										m_wavel;
	std::vector<double>										m_amfheights;
	nx1dArray<double>										m_hardcoderadiance;
	nx1dArray<double>										m_hardcodevariance;
	nx2dArray<double>										m_hardcodeamf;
	nx2dArray<double>										m_hardcodeamfvariance;
	std::string												m_testname;
	nxVector												m_sun;
	int														m_scatterorder;
	std::string												m_hardcodefolder;
	bool													m_verbose;
	double													m_tol;

private:
	bool													RunCurrentTest();
	void													MakeAtmosphericState();
	void													MakeLinesOfSight();

	void                                                    LoadHardcodedValues(int numLOS, int numAMF);

	void                                                    SetUpDefault();
	void													SetUpShortTest();
	void													SetUpLengthTest();
	void													SetUpOpticalDepthTest();

	bool													CheckSuccess(const nx1dArray<double>& radiance, const nx1dArray<double>& variance, const nx2dArray<double>& airmassfactor, const nx2dArray<double>& airmassfactorvariance);
	void													DisplayFullTestInfo(const nx1dArray<double>& radiance, const nx1dArray<double>& variance, const nx2dArray<double>& airmassfactor, const nx2dArray<double>& airmassfactorvariance);

public:
	SKTRAN_MCAMF_Test(bool verbose, double tolerance);
	bool													RunShortTests();
	bool													RunAllTests();
};