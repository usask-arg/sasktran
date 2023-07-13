#pragma once

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sasktran.h>


class SKTRAN_Polarization_Test
{
	public:

		double										m_tolerance;
		std::string									m_hardcodefolder;
		bool										m_writeToFile;
		bool										m_writeErrorsToScreen;
		std::vector<double>							m_wavs;
		std::vector<double>							m_alts;
		std::vector<double>							m_szas;
		std::vector<double>							m_saas;
		size_t										m_geometryIndex;
		int											m_seed;
		double										m_precisionMC;
		SKTRAN_Specifications_MC::SunType			m_sunType;
		SKTRAN_Specifications_MC::SolarTableType	m_sttype;
		size_t										m_numThreadsMC;
		size_t										m_numThreadsHR;
		size_t										m_maxOrderScatter;

	public:
													SKTRAN_Polarization_Test				(double errortolerance );
										 		   ~SKTRAN_Polarization_Test				( );  
		void										SetupForTest							( );
		bool										LoadReferenceValues						( std::string& filename, nx2dArray<skRTStokesVector>& referenceValues ) const; 
		bool										WriteReferenceValuesToFile				( std::string& filename, const nx2dArray<skRTStokesVector>& newReferenceValues ) const; 
		bool										CompareValues							( const nx2dArray<skRTStokesVector>& newValues, const nx2dArray<skRTStokesVector>& referenceValues, bool writeErrorsToScreen ) const; 
		bool										MakeSingleScatterComparisonLinesOfSight	( SKTRAN_LineOfSightArray_V21& linesofsight, nxVector& sunvec );
		bool										CreateOpticalState						( SKTRAN_AtmosphericOpticalState_V21* opticalstate );
		bool										runHR									(  int pvorder, SKTRAN_HR_PolHOType pvhotype );
		bool										runMC									( );
		int											main									( );
    
};


