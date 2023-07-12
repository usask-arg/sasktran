#pragma once

#include <sasktran.h>
#include <memory>

class SKTRAN_MC_Test
{
	private:
		std::unique_ptr<SKTRAN_Engine_MC_V21>					m_engine;
		std::unique_ptr<SKTRAN_AtmosphericOpticalState_V21>		m_opticalstate;
		std::unique_ptr<SKTRAN_LineOfSightArray_V21>			m_linesofsight;
		std::unique_ptr<SKTRAN_Specifications_MC>				m_specs;
		std::vector<double>										m_wavel;
		nx2dArray<double>										m_hardcoderadiance;
		nx2dArray<double>										m_hardcodevariance;
		nx2dArray<double>                                       m_hardcodeangles;
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
		
		void                                                    LoadHardcodedValues ( int numLOS, int numWavs );

		void                                                    SetUpDefault ( );
		void													SetUpShortTest();
		void 													SetUpSingleScattTest ();
		void 													SetUpSecondOrderTest ();
		void 													SetUpPseudoPolTest   ();
		void 													SetUpPolarizationTest();
		void 													SetUpHighPrecisionTest();

		void													DisplayFullTestInfo( const nx2dArray<double>& radiance, const nx2dArray<double>& variance, const nx2dArray<double>& angles );
		
	public:
																SKTRAN_MC_Test( bool verbose, double tolerance );
		void													RunShortTests();
		void													RunAllTests();
};