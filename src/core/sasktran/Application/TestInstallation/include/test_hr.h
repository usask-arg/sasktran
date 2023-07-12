#pragma once

#include <sasktran.h>
#include <memory>

class SKTRAN_HR_Test
{
	private:
		std::unique_ptr<SKTRAN_HR_Engine>						m_engine;
		std::unique_ptr<SKTRAN_AtmosphericOpticalState_V21>		m_opticalstate;
		std::unique_ptr<SKTRAN_LineOfSightArray_V21>			m_linesofsight;
		std::unique_ptr<SKTRAN_HR_Specs_User>					m_specs;
		std::vector<double>										m_wavel;
		nx2dArray<double>										m_hardcodevalues;
		std::string												m_testname;
		nxVector												m_sun;
		int														m_scatterorder;
		std::string												m_hardcodefolder;
		bool													m_verbose;
		double													m_tol;

		bool													RunCurrentTest();
		void													MakeAtmosphericState();
		void													MakeLinesOfSight();

		void													SetupSingleScatterTest();
		void													SetupShortTest();
		void													SetupReplicateSOTest();
		void													SetupMultipleDiffuseProfileTest();
		void													SetupTwoDimTest();

		void													DisplayFullTestInfo( const nx2dArray<double>& radiance );
		
	public:
		SKTRAN_HR_Test( bool verbose, double tolerance );
		void													RunShortTests();
		void													RunAllTests();
};