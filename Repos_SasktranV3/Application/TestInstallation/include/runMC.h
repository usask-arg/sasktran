#pragma once

//#include <sasktran.h>

//#include <modules/monte_carlo/include/sktran_montecarlo.h>
//#include <boost/timer/timer.hpp>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <sasktran.h>

class ProfilerMC
{
public:
	std::string m_dirstr;

	int m_runctr;
	std::vector<double> m_wavs;
	std::vector<double> m_alts;
	std::vector<double> m_lats;
	std::vector<double> m_szas;
	std::vector<double> m_saas;
	std::vector<double> m_albs;
	std::vector<int>    m_userayleigh;
	std::vector<int>	m_useozone;
	std::vector<int>    m_useaero;
    std::vector<int>    m_useemission;
	std::vector<int>    m_useno2;
    bool                m_useo2aband;
	size_t m_wavidx;
	size_t m_altidx;
	size_t m_latidx;
	size_t m_szaidx;
	size_t m_saaidx;
	size_t m_albidx;
	size_t m_uraidx;
	size_t m_uo3idx;
	size_t m_uaeidx;
    size_t m_uemidx;

	int    m_seed;
	double m_precisionMC;
	std::string m_inputsdir;
	std::string m_outfilename;
    std::string m_aeroProfileFile;
    std::string m_emissionTableFile;
	SKTRAN_Specifications_MC::SunType m_sunType;
	SKTRAN_Specifications_MC::SolarTableType m_sttype;
	SKTRAN_Specifications_MC::PolType m_polType;
	SKTRAN_Specifications_MC::SymmetryType m_smType;

	size_t m_numThreads;
	size_t m_maxOrderScatter;

    std::vector<std::string> namestr;
    bool m_writeToFile;
    bool m_writeToScreen;
    bool m_ignoreExistingFile;

public:
	      ProfilerMC();
	     ~ProfilerMC();
    void  LoadParamsFromFile( nx1dArray<double>& wavsToDo, std::string& outdir, std::string& emissiontablefilepath, size_t& maxOrderOfScatter, double& albedo );
	void  SetDefaults();
	void  LoadFromInputFile( const char* filename );
	void  PrintOptions();
	void  SetSunType   (SKTRAN_Specifications_MC::SunType t) {m_sunType = t;}
	void  SetStType    (SKTRAN_Specifications_MC::SolarTableType t) {m_sttype = t;}
	void  SetPolType    (SKTRAN_Specifications_MC::PolType   t) {m_polType = t;}
	SKTRAN_Specifications_MC::PolType GetPolType() const {return m_polType; }
	void  SetSmType    (SKTRAN_Specifications_MC::SymmetryType   t) {m_smType = t;}
	void  SetSeed      (int s) {m_seed = s;}
	void  SetNumThreads (size_t s) {m_numThreads = s;}
	void  SetAeroProfileFile(std::string s) {m_aeroProfileFile = s;}
	void  SetPrecisionMC (double d ) { m_precisionMC = d;}
	void  SetMaxOrderScatter (size_t s) {m_maxOrderScatter=s;}
	bool  ConfigureStratAerosolOpticalProps( skOpticalProperties_AerosolProfileH2SO4* optaerosol, const char* climatologyDatabaseDir  );
	bool  MakeSingleScatterComparisonLinesOfSight( SKTRAN_LineOfSightArray_V21& linesofsight, nxVector& sunvec );
	bool  Make2dImageObservers                   ( SKTRAN_LineOfSightArray_V21& linesofsight, nxVector& sunvec );
	bool  CreateOpticalState( SKTRAN_AtmosphericOpticalState_V21* opticalstate );
	//bool  runHR( const SKTRAN_LineOfSightArray_V21& linesofsight, SKTRAN_AtmosphericOpticalState_V21& opticalstate, const nxVector& sun, nx2dArray<double>& radiance, bool usecurve, std::string filename, nx2dArray<double>& retrad, nx2dArray<skRTStokesVector>* retsvec );
	//int   runsimHR( nx2dArray<double>& retrad, nx2dArray<skRTStokesVector>* retsvec );
	bool runHR( nx2dArray<double>& retrad, nx2dArray<skRTStokesVector>* retsvec, int pvorder, SKTRAN_HR_PolHOType pvhotype );
	//bool  runMC( const SKTRAN_LineOfSightArray_V21& linesofsight, SKTRAN_AtmosphericOpticalState_V21& opticalstate, const nxVector& sun, nx2dArray<double>& radiance, bool usecurve, std::string filename, nx2dArray<double>& retrad, nx2dArray<skRTStokesVector>* retsvec );
	bool runMC( nx2dArray<double>& retrad, nx2dArray<skRTStokesVector>* retsvec, nx2dArray<double>& timing );
	//int   runsim( nx2dArray<double>& retrad, nx2dArray<skRTStokesVector>* retsvec );
	int   main( std::vector<nx2dArray<double>>& retrad, std::vector<nx2dArray<skRTStokesVector>>* retsvec, std::vector<std::string>& retstr );

	int SuggestNumDiffuseProfiles() const;
	void SetDirStr(std::string s) {m_dirstr = s;}

};


//
//int main(int argc, char* argv[])
//{	
//
//	{
//	std::vector<nx2dArray<double>> retrad;
//	std::vector<nx2dArray<skRTStokesVector>> retsvec;
//	std::vector<std::string> retstr;
//	stringstream sstr;
//
//	ProfilerMC m;
//	m.SetSunType(SKTRAN_Specifications_MC::sunType_point);
//	m.SetStType(SKTRAN_Specifications_MC::stType_2d);
//	m.SetStType(SKTRAN_Specifications_MC::stType_noTable);
//	m.SetNumThreads( 0 );
//	nx2dArray<double> toprint;
//
//	// Change to vectorized code
//	m.SetPrecisionMC(0.01);
//
//	m.SetSeed( 0 );
//	m.SetDirStr("C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/2014/trash/mc_");
//
//
//	// NO AEROSOL
//	m.useaero.at(0) = 0;
//
//	int* lkptr = new int;
//	delete lkptr;
//	//m.SetOtType(SKTRAN_Specifications_MC::otType_1D);
//	//m.SetMaxOrderScatter(50);
//	//m.SetDirStr("C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/2014/trash/mc_");
//	//m.main( retrad, NULL, retstr );
//	//retrad.clear(); retsvec.clear(); retstr.clear();
//	//std::printf("\n\n");
//
//	//m.SetOtType(SKTRAN_Specifications_MC::otType_1dPol);
//	//m.SetMaxOrderScatter(50);
//	//m.SetDirStr("C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/2014/poltesting/vector50_noAero/mc_");
//	//m.main( retrad, &retsvec, retstr );
//	//retrad.clear(); retsvec.clear(); retstr.clear();
//	//std::printf("\n\n");
//
//	//m.SetOtType(SKTRAN_Specifications_MC::otType_1dPol_pseudo);
//	//m.SetMaxOrderScatter(2);
//	//m.SetDirStr("C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/2014/poltesting/vector2_pseudovecs_noAero/mc_");
//	//m.main( retrad, &retsvec, retstr );
//	//retrad.clear(); retsvec.clear(); retstr.clear();
//	//std::printf("\n\n");
//
//	//m.SetOtType(SKTRAN_Specifications_MC::otType_1dPol);
//	//m.SetMaxOrderScatter(50);
//	////m.SetDirStr("C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/2014/poltesting/vector50_noAero/mc_");
//	//m.main( retrad, &retsvec, retstr );
//	//retrad.clear(); retsvec.clear(); retstr.clear();
//	//std::printf("\n\n");
//
//	m.SetPolType(SKTRAN_Specifications_MC::SKTRAN_Specifications_MC_PolType::polType_none);
//	
//	m.SetMaxOrderScatter(50);
//	//m.SetDirStr("C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/2014/convtesting/vardistr/"); // Include trailing /
//	m.SetDirStr(""); 
//	m.main( retrad, m.GetPolType()==SKTRAN_Specifications_MC::SKTRAN_Specifications_MC_PolType::polType_none? nullptr : &retsvec , retstr );
//	retrad.clear(); retsvec.clear(); retstr.clear();
//	std::printf("\n\n");
//
//	//m.SetOtType(SKTRAN_Specifications_MC::otType_1dPol);
//	//m.SetMaxOrderScatter(50);
//	//m.SetDirStr("C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/2014/poltesting/vector50_vecs/mc_");
//	//m.main( retrad, &retsvec, retstr );
//	//retrad.clear(); retsvec.clear(); retstr.clear();
//	//std::printf("\n\n");
//
//
//	}
//	
//	printf("done\n");
//	_CrtDumpMemoryLeaks();
//	system("pause");
//
//	return 0;
//
//}
//
//
//
