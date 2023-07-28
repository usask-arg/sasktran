#include<sasktran.h>
#include "include/runMC.h"
#include<algorithm>
#include <boost/timer/timer.hpp>

#include <stdlib.h>

 //   {
	//nxLogConsole	thelog;
	//std::timer::auto_cpu_timer t;
	////themain(argc, argv);
	//std::vector<nx2dArray<double>> retrad;
	//std::vector<nx2dArray<skRTStokesVector>> retsvec;
	//std::vector<std::string> retstr;
	//ProfilerMC pmc;
	//pmc.main( retrad, &retsvec, retstr );
 //   }

ProfilerMC::ProfilerMC()
{
	m_wavidx = 0;
	m_altidx = 0;
	m_latidx = 0;
	m_szaidx = 0;
	m_saaidx = 0;
	m_albidx = 0;
	m_uraidx = 0;
	m_uo3idx = 0;
    m_uaeidx = 0;
    m_uemidx = 0;
    m_precisionMC = 0.01;
    m_outfilename = "";
	m_inputsdir   = "C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/inputs";
    //m_dirstr    = "C:/ARGsoftware/pv-paper/hrmccomp/data/full_default/";
    //m_dirstr    = "C:/ARGsoftware/pv-paper/hrmccomp/data/trash/";
    m_dirstr    = "C:/ARGsoftware/a-band-spectra/data/data_single_saa30_emissionOnly/";
	//m_aeroProfileFile = "/Aerosol_6432012.txt";
	m_aeroProfileFile = "/AerosolProfileTropics.txt";
	m_sunType   = SKTRAN_Specifications_MC::SunType::point;
	m_polType = SKTRAN_Specifications_MC::PolType::pol;
    m_sttype = SKTRAN_Specifications_MC::SolarTableType::noTable;
	m_smType = SKTRAN_Specifications_MC::SymmetryType::none;

	m_maxOrderScatter = 10;
	m_seed = 0;

    m_useo2aband  = true;

    SetDefaults();

}

ProfilerMC::~ProfilerMC()
{

}

void ProfilerMC::LoadParamsFromFile( nx1dArray<double>& wavsToDo, std::string& outdir, std::string& emissiontablefilepath, size_t& order, double& albedo )
{
    std::ifstream istream;
    istream.open( "C:/ARGsoftware/a-band-spectra/sktranconfig.txt", std::ifstream::in );

    std::string wavsToDoPath;

    istream >> wavsToDoPath;
    istream >> outdir;
    istream >> emissiontablefilepath;
    istream >> order; 
    istream >> albedo;
    istream.close();
    
    //emissiontablefilepath = "C:/ARGsoftware/a-band-spectra/tabulatedEmission_nxFormat.txt";
	//order = 10;
	
    size_t numWavsToDo;
    istream.open( wavsToDoPath.c_str() );
    istream >> numWavsToDo;
    wavsToDo.SetSize( numWavsToDo );
    istream >> wavsToDo;
    istream.close();

    return;
}


void ProfilerMC::SetDefaults()
{

    nx1dArray<double> wavsToDo;
    double albedo;
    LoadParamsFromFile( wavsToDo, m_dirstr, m_emissionTableFile, m_maxOrderScatter, albedo );

    m_wavs.resize(0);
    //m_wavs.push_back(350.0);
    //m_wavs.push_back( 767.875 );
    //m_wavs.push_back( 290.0 );
    //m_wavs.push_back(600.0);
	//m_wavs.push_back(380.0);
	//m_wavs.push_back(400.0);
	//m_wavs.push_back(500.0);
	//m_wavs.push_back(550.0);
	//m_wavs.push_back(600.0);
	//m_wavs.push_back(650.0);
	//m_wavs.push_back(700.0);
	//m_wavs.push_back(750.0);
	//m_wavs.push_back(800.0);
	//m_wavs.push_back(850.0);
	//m_wavs.push_back(900.0);
	//m_wavs.push_back(950.0);

    //m_wavs.push_back(749.94);
    //m_wavs.push_back(602.39);
    //m_wavs.push_back(322.50);
    //m_wavs.push_back(350.31);
    
    //for (double wav = 763.0261; wav < 763.037; wav += 0.0002) {
    //for ( double wav=lowav; wav<hiwav; wav+=0.0002 ){
    //    m_wavs.push_back( wav );
    //}
    //for (int widx = 0; widx < wavsToDo.size(); ++widx) {
    //    m_wavs.push_back( wavsToDo.At(widx) );
    //}
    //m_wavs.push_back( 7.77986105168162e+02 );
    m_wavs.push_back(7.63313620759077e+02);
    m_wavs.push_back(7.63314203407206e+02);
    m_wavs.push_back(7.63315077381066e+02);
    

	m_lats.resize(0);
	m_lats.push_back(   75.0 );
    
	m_alts.resize(0);
    //for(int repidx=10000; repidx<120000; repidx+=1000 ){
    //    m_alts.push_back( repidx );
    //}   
    //m_alts.push_back( 115000.0 );
    m_alts.push_back(100000.0);
    //m_alts.push_back(100000.0);
    //m_alts.push_back( 10000.0 );
    ////m_alts.push_back( 20500.0 );
    ////m_alts.push_back( 30000.0 );
    //m_alts.push_back( 40000.0 );
    ////m_alts.push_back( 50000.0 );
    ////m_alts.push_back( 60000.0 );
    //m_alts.push_back( 70000.0 );
    ////m_alts.push_back( 80000.0 );
    ////m_alts.push_back( 90000.0 );
    
	m_szas.resize(0);
    //m_szas.push_back(88.0);
    //m_szas.push_back(88.0);
    //m_szas.push_back(30.0);
    //m_szas.push_back(20.0);
    m_szas.push_back(60.0);
    //m_szas.push_back(80.0);
    //m_szas.push_back(89.0);

	
    m_saas.resize(0);
    //m_saas.push_back(   0.0 );
    //m_saas.push_back(  45.0 );
    //m_saas.push_back(  90.0 );
    //m_saas.push_back( 135.0 );
    //m_saas.push_back( 90.0 );
    //m_saas.push_back( 150.0 );
    m_saas.push_back( 30.0 );

    m_albs.resize(0);
    m_albs.push_back(albedo);
    //m_albs.push_back(0.3);
    //m_albs.push_back(0.0);

    m_userayleigh.resize(0);
    //m_userayleigh.push_back(1);
	m_userayleigh.push_back(1);

    m_useozone.resize(0);
    m_useozone.push_back(1);

    m_useaero.resize(0);
    //m_useaero.push_back(0);
    m_useaero.push_back(1);

    m_useemission.resize(1);
    m_useemission[0] = 1;

	//m_useno2.resize(0);
	//m_useno2.push_back(1);

}

void ProfilerMC::LoadFromInputFile( const char* filename )
{
/*
    std::ifstream infile;
    infile.open(filename);
    std::string temp;
    double tempd;
    std::stringstream stream;

    infile >> temp; //
    infile >> sza;
    infile >> temp;
    infile >> saa;
    infile >> temp;
    getline(infile, temp);
    stream.str(temp);
    while( stream.good() )
    {
        stream >> tempd;
        alt.push_back( tempd );
    }
    stream.clear();
    infile >> temp;
    getline(infile, temp);
    stream.str(temp);
    while( stream.good() )
    {
        stream >> tempd;
        wavelen.push_back( tempd );
    }
    infile >> temp;
    infile >> numprofile;
    infile >> temp;
    infile >> albedo;
    infile >> temp;
    infile >> userayleigh;
    infile >> temp;
    infile >> useozone;
    infile >> temp;
    infile >> useaero;
    infile >> temp;
    infile >> outfilename;
*/

}

void ProfilerMC::PrintOptions()
{
    std::cout << "sza = "<< m_szas[m_szaidx] << std::endl;
    std::cout << "saa = " << m_saas[m_saaidx] << std::endl;
    std::cout << "albedo = " << m_albs[m_albidx] << std::endl;
    std::cout << "rayleigh = "<< m_userayleigh[m_uraidx] <<  std::endl;
    std::cout << "aero = " << m_useaero[m_uaeidx] << std::endl;
    std::cout << "ozone = " << m_useozone[m_uo3idx] << std::endl;
    std::cout << "emission = " << m_useemission[m_uemidx] << std::endl;

    for(size_t i = 0; i < m_alts.size(); i++ )
        std::cout << m_alts[i] << " ";
    std::cout << std::endl;
    for(size_t i = 0; i < m_wavs.size(); i++ )
        std::cout << m_wavs[i] << " ";
    std::cout << std::endl;
    std::cout << m_outfilename << std::endl;
}

bool ProfilerMC::ConfigureStratAerosolOpticalProps( skOpticalProperties_AerosolProfileH2SO4* optaerosol, const char* climatologyDatabaseDir  )
{
    nxString	        aersizeparamfilename;           // file where the aerosol parameters are stroed
    nx2dArray<double>       aerosolparams;                  // the actual array of aerosol parameters
    nx1dArray<double>       buffer;                         // temporary buffer
    nx1dArray<double>       ah;                             // altitudes
    nx1dArray<double>       moderad;                        // mode radius
    nx1dArray<double>       modewidth;                      // mode width
    bool                    ok;

    // get the actual file name + path of the aerosol parameters
    aersizeparamfilename.sprintf( "%s/StratAerosolSizeParameters.txt", climatologyDatabaseDir );

    // load in the data
    ok = aerosolparams.InputColumnMajorText( aersizeparamfilename, 3, 0 );
    if (ok)
    {
        // read in the data to the specified variables
        ah.SetSize       ( aerosolparams.XSize( ) );
        moderad.SetSize	 ( aerosolparams.XSize( ) );
        modewidth.SetSize( aerosolparams.XSize( ) );
        for (size_t altctr=0; altctr<aerosolparams.XSize(); altctr++ )
        {
            ah.At       (altctr) = 1000.0*aerosolparams.At(altctr,0);
            moderad.At  (altctr) = aerosolparams.At(altctr,1);
            modewidth.At(altctr) = aerosolparams.At(altctr,2);
        }
        // set the log normal climatology
        ok = optaerosol->SetLogNormalProfileClimatology( ah.UnsafeArrayBasePtr(), moderad.UnsafeArrayBasePtr(), modewidth.UnsafeArrayBasePtr(), ah.size() );
    }

    if (!ok)
    {
        nxLog::Record(NXLOG_WARNING,"ConfigureAerosolOpticalProps, Error configuring aerosol mie scattering height profile from file <%s>", (const char*)aersizeparamfilename );
    }
    return ok;
}


bool ProfilerMC::MakeSingleScatterComparisonLinesOfSight( SKTRAN_LineOfSightArray_V21& linesofsight, nxVector& sunvec )
{
	bool 	ok = true;
	double	mjd 	= 54832.0;
	 //std::vector<double> tanheights_meters;
        /*
        for( size_t idx = 0; idx < 60; idx++ )
        {
            tanheights_meters.push_back( 10000 + 1000*idx );
        }
        */

	
	//tanheights_meters[0] = 10000;
	//for (int i=0; i < 36; i++) tanheights_meters[i] = i*2000.0 + 10000;	
	//ok = linesofsight.SetRaysFromTangentHeightArray( mjd, lat, lng, sza, saa, rayazi, tanheights_meters, 36, 600000.0, NULL );
	//for( altidx = 0; altidx < alts.size(); altidx++ )
	
		//ok = linesofsight.AddEquatorialLineOfSight( szas[szaidx], saas[saaidx], alts[altidx], 600000.0, mjd );

	double* alist = new double[m_alts.size()];
	for(int aidx = 0; aidx< (int)m_alts.size(); aidx++){
		alist[aidx] = m_alts[aidx];
	}
	ok = ok && linesofsight.SetRaysFromTangentHeightArray(54000, m_lats[m_latidx], 0.0, m_szas[m_szaidx], m_saas[m_saaidx], 0.0, alist, (int) m_alts.size(), 600000.0, &sunvec);
	
    //nxVector obs(4306163.7358194590, 0.00000000000000000, 5474201.7747975970 );
    //nxVector look(-0.5810, 0.3420, -0.7386 );
    //sunvec.SetCoords(-0.017337588530253717, -0.43301270189221924, 0.90122106501343813);

    //ok = ok && linesofsight.AddLineOfSight( obs, look, 54000.0 );
	
    delete[] alist;

	return ok;
}


bool ProfilerMC::Make2dImageObservers( SKTRAN_LineOfSightArray_V21& linesofsight, nxVector& sunvec )
{
	bool ok = true;
	double mjd = 54832.0;

	double minUpAngle = 0.0;
	double maxUpAngle = 0.5; // Tangent point is at 100km for 1.954deg
	double minLRAngle = -5.0; // 500km at tangent point corresponds to 10deg
	double maxLRAngle =  5.0;
	int numUDangles = 20;
	int numLRAngles = 10;
	nxVector obs(6382140.0000000000+30000.0, 0.00000000000000000, -2813183.7949759066);
	nxVector look;
	const double rotBehindSun = -0.0* nxmath::Pi / 180.0;
	double R[9] = {cos(rotBehindSun), 0.0, -sin(rotBehindSun), 0.0, 1.0, 0.0, sin(rotBehindSun), 0.0, cos(rotBehindSun)};
	obs.SetCoords( R[0]*obs.X()+R[1]*obs.Y()+R[2]*obs.Z(), R[3]*obs.X()+R[4]*obs.Y()+R[5]*obs.Z(), R[6]*obs.X()+R[7]*obs.Y()+R[8]*obs.Z() );

	for(int udidx=0; udidx<numUDangles; ++udidx){
		for(int lridx = 0; lridx<numLRAngles; ++lridx){
			double udang = ( ((maxUpAngle-minUpAngle)/(numUDangles-1)) * udidx + minUpAngle ) * ( nxmath::Pi / 180.0 );
			double lrang = ( ((maxLRAngle-minLRAngle)/(numLRAngles-1)) * lridx + minLRAngle ) * ( nxmath::Pi / 180.0 );
			look.SetCoords( sin(udang)*cos(lrang), sin(lrang), cos(udang)*cos(lrang) );
			look.SetCoords( R[0]*look.X()+R[1]*look.Y()+R[2]*look.Z(), R[3]*look.X()+R[4]*look.Y()+R[5]*look.Z(), R[6]*look.X()+R[7]*look.Y()+R[8]*look.Z() );
			ok = ok && linesofsight.AddLineOfSight( obs, look.UnitVector(), mjd );
		}
	}

	sunvec.SetCoords(0.0, 0.0, 1.0);

	return ok;
}

/*-----------------------------------------------------------------------------
 *					MakeSpeciesList		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool ProfilerMC::CreateOpticalState( SKTRAN_AtmosphericOpticalState_V21* opticalstate )
{
	bool                                       ok = true;
	bool                                       ok_nonempty = false;

	skOpticalProperties_RayleighDryAir*        rayleigh;					// Optical properties of one air molecule
	skClimatology_MSIS90*                      msis90;	
	skClimatology_LabowOzoneVMR*               o3numberdensity;
	skOpticalProperties_O3_OSIRISRes*          o3_opticalprops;			// optical properties of one O3 molecule
    skClimatology_UserDefinedTable*            aerosoldensity;
    skOpticalProperties_AerosolProfileH2SO4*   aero_opticalprops;

	std::string indir;

	ok = ok && NULL!=opticalstate;
	if(ok) opticalstate->erase();
	
	msis90           = new skClimatology_MSIS90;
	rayleigh         = new skOpticalProperties_RayleighDryAir;
	ok = ok && (NULL!=msis90) && (NULL!=rayleigh);
	
	if( 0.1 < m_userayleigh[m_uraidx] ){
		ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, msis90, rayleigh );
		ok_nonempty = true;
	}
	if( 0.1 < m_useozone[m_uo3idx] ){
		o3numberdensity  = new skClimatology_LabowOzoneVMR;
		o3_opticalprops  = new skOpticalProperties_O3_OSIRISRes;
		ok = ok && (NULL!=o3numberdensity) && (NULL!=o3_opticalprops);
		ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_O3_CM3, o3numberdensity, o3_opticalprops);
		ok_nonempty = true;
	}
    if( 0.1 < m_useaero[m_uaeidx] ){
        indir = m_inputsdir + m_aeroProfileFile;
        aerosoldensity	        	= new skClimatology_UserDefinedTable(SKCLIMATOLOGY_AEROSOLH2SO4_CM3, indir.c_str());
        //skClimatology_Constant* aerosoldensity = new skClimatology_Constant( 10.0 ); 
        aero_opticalprops           = new skOpticalProperties_AerosolProfileH2SO4;
        ok = ok && (NULL!=aerosoldensity) && (NULL!=aero_opticalprops);
        ok = ok && ConfigureStratAerosolOpticalProps( aero_opticalprops, m_inputsdir.c_str()  );
        ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_AEROSOLH2SO4_CM3, aerosoldensity,  aero_opticalprops );
        ok_nonempty = true;
    }
    if( 0.1 < m_useemission[m_uemidx] ){    
        skEmission_Tabulated_HeightWavelength* emission = new skEmission_Tabulated_HeightWavelength; 
        ok = ok && emission->LoadHeightWavelengthProfileFromFile( m_emissionTableFile.c_str() );
        ok = ok && opticalstate->AddEmission( SKEMISSION_PHOTOCHEMICAL_0, emission );
        ok_nonempty = true;
    }

    if( m_useo2aband && (749.06<m_wavs[m_wavidx]) && (778.82>m_wavs[m_wavidx])  ){
        skClimatology_MSIS90*  msis90 = new skClimatology_MSIS90;	
        skOpticalProperties_HitranChemical* o2 = new skOpticalProperties_HitranChemical("O2");
		o2->SetWavenumberRange(12839.9373411058, 13350.0654153205 );
        ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_O2_CM3, msis90, o2 );
    }

	ok = ok && opticalstate->SetAlbedo( m_albs[m_albidx] );

	ok = ok && ok_nonempty;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"MakeSpeciesList, there was an error making the species list.");
		if(NULL!=opticalstate){
			opticalstate->erase();
		} else{
			delete rayleigh;
			delete msis90;	
			delete o3numberdensity;
			delete o3_opticalprops;	
			delete aerosoldensity;
			delete aero_opticalprops;
		}
	}

	return ok;
}

//mc_radiance, false, outfilename, retrad, retsvec 
bool ProfilerMC::runMC( nx2dArray<double>& retrad, nx2dArray<skRTStokesVector>* retsvec, nx2dArray<double>& timing )
{
	bool ok = true;
    nx2dArray<double> radiance;

    {
	SKTRAN_Engine_MC_V21				engine_mc;
	std::vector<SKTRAN_StokesScalar>    radiance_mc;
	std::vector<skRTStokesVector>     radiance_vec;
	SKTRAN_Specifications_MC			specs_mc;
	nx1dArray<skRTStokesVector>         tempsvec;

	SKTRAN_LineOfSightArray_V21			linesofsight;
	nx2dArray<double>					mc_radiance;
	SKTRAN_AtmosphericOpticalState_V21	opticalstate;
	nxVector							sun ( 1.0, 0, 0 );
	ok = ok && MakeSingleScatterComparisonLinesOfSight( linesofsight, sun );
    //printf("\n\n mcsun: %17.14e %17.14e %17.14e\n\n", sun.X(), sun.Y(), sun.Z());
	ok = ok && CreateOpticalState( &opticalstate );


    ok = ok && specs_mc.ConfigureDefaults        ( );
    ok = ok && specs_mc.SetSunGeographicPosition ( sun );
    ok = ok && specs_mc.SetSineSunApexAngle      ( nxmath::sind(0.268) ); // Angle from sun center to sun edge at Earth
    ok = ok && specs_mc.SetScatterPositionRes    ( 10.0 );
    ok = ok && specs_mc.SetAllowDynamicThreads   ( false );

    ok = ok && specs_mc.SetMinFractionHigherOrder( 0.1 );
    ok = ok && specs_mc.SetMinimumRelPathWeight  ( 0.0*m_precisionMC*1e-3 );
    //ok = ok && specs_mc.SetMinFractionHigherOrder( 0.1 );
    //ok = ok && specs_mc.SetMinimumRelPathWeight  ( m_precisionMC*1e-3 );
	ok = ok && specs_mc.SetPrecisionMC           ( m_precisionMC );
    ok = ok && specs_mc.SetNumPhotonsPerLOS      ( true ? (size_t)(3.0/(m_precisionMC*m_precisionMC)) : 1 );
    //ok = ok && specs_mc.SetChunkSize             ( 1 );
    ok = ok && specs_mc.SetRngSeed               ( m_seed );
    ok = ok && specs_mc.SetSunType               ( m_sunType );
    ok = ok && specs_mc.SetSolarTableType        ( m_sttype );
	//ok = ok && specs_mc.SetThermalEmissionType   ( SKTRAN_Specifications_MC::ThermalType::none );
	ok = ok && specs_mc.SetSolarRayTracerType    ( SKTRAN_Specifications_MC::RayTracerType::shell );
    ok = ok && specs_mc.SetMSRayTracerType       ( SKTRAN_Specifications_MC::RayTracerType::shell );
	ok = ok && specs_mc.SetLOSRayTracerType      ( SKTRAN_Specifications_MC::RayTracerType::shell );
	ok = ok && specs_mc.SetOptPropIntType        ( SKTRAN_Specifications_MC::OptPropIntType::straight );
	ok = ok && specs_mc.SetKernelType            ( SKTRAN_Specifications_MC::LogType::none );
    ok = ok && specs_mc.SetPolType				 ( SKTRAN_Specifications_MC::PolType::pol /*m_polType */);
    //ok = ok && specs_mc.SetPolType				 ( SKTRAN_Specifications_MC::PolType::none );
	ok = ok && specs_mc.SetSymmetryType          ( m_smType );
	ok = ok && specs_mc.SetAdaptOptDepthMax      ( 0.05 );
    //ok = ok && specs_mc.SetSolarTableType        ( SKTRAN_Specifications_MC::SolarTableType::doNothing );
    ok = ok && specs_mc.SetEmissionTableType     ( SKTRAN_Specifications_MC::EmissionTableType::doNothing );
    
    ok = ok && specs_mc.FinalizeSpecs            ( );


	ok = ok && engine_mc.ConfigureModel( specs_mc, linesofsight, m_numThreads);
	radiance_mc.resize( linesofsight.NumRays() );
	radiance_vec.resize( linesofsight.NumRays() );

	radiance.SetSize( linesofsight.NumRays(), 1/*wavs.size()*/ );
	radiance.SetTo(-1.0);

//	for( size_t wavidx = 0; wavidx < wavs.size(); wavidx++ )
	{
                //std::cout << wavelen[i] << std::endl;
		std::stringstream cloudstream;
		//cloudstream << "C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/2013/randomDisc/AerosolProfileTropics/output/scattcloud_point/saa" << saas[saaidx];
		//engine_mc.SetSaveDir(cloudstream.str());
		ok = ok && engine_mc.CalculateRadiance(&radiance_mc, m_wavs[m_wavidx], m_maxOrderScatter, &opticalstate, &radiance_vec );
        ok = ok && engine_mc.GetTimingData( timing );
		
		//hr_engine.CalculateSecondOrderRadiance(&hr_radiance_temp, wavelen[i],2,&opticalstate );
		for(size_t losidx = 0; losidx < std::min(radiance_mc.size(),radiance.XSize()); losidx++ )
		{
			radiance.At(losidx,0/*wavidx*/) = radiance_mc[losidx];
			//printf("%u:\t%1.10e\n", m_runctr++, radiance_mc[losidx]);
			//printf("%1.16e\n", radiance_mc[losidx]); m_runctr++;
		}
        
		if(NULL!=retsvec){
			//ok = ok && engine_mc.GetStokesVectors(radiance_vec);
			retsvec->SetSize(radiance_vec.size(),4);
			for(size_t losidx = 0; losidx < std::min(radiance_mc.size(),radiance.XSize()); ++losidx){
				retsvec->At(losidx,0) = radiance_vec[losidx];
			}
		}

	}
	//printf(";\n");
    }
        
    if(ok) retrad = radiance;
    
    return ok;
}


int ProfilerMC::SuggestNumDiffuseProfiles() const
{
	const int numdpszas = 12;
	const int numdpsaas = 11;
	const int numdps_hardmax = 15;

	int dpszas[numdpszas] = {/*15*/ 0, 45, 60, 70, 80, 85, 87, 88, 89, 90 ,91, 92};
	int dpsaas[numdpsaas] = {0, 20, 40, 60, 80, 90, 100, 120, 140, 160, 180};
	int dptable[numdpszas*numdpsaas] = {
	2,	2,	2,	2,	2,	2,	1,	2,	2,	2,
	2,	2,	2,	2,	2,	1,	2,	2,	2,	2,
	3,	3,	3,	2,	2,	1,	2,	3,	4,	4,
	4,	4,	3,	3,	2,	1,	2,	4,	4,	4,
	4,	4,	4,	3,	2,	1,	2,	4,	4,	4,
	6,	7,	6,	2,	2,	1,	3,	4,	4,	4,
	23,	20,	15,	9,	4,	1,	3,	4,	4,	4,
	46,	38,	26,	13,	4,	1,	3,	4,	4,	4,
	81,	71,	46,	19,	6,	1,	3,	4,	4,	5,
	91,	86,	71,	32,	8,	2,	4,	6,	6,	6,
	96,	91,	81,	46,	11,	3,	6,	13,	10,	9,
	96,	96,	86,	61,	14,	3,	9,	18,	18,	17
	};
	
	int dpszahiidx;
	for(dpszahiidx=0; dpszahiidx<numdpszas && dpszas[dpszahiidx]<=m_szas[m_szaidx]; ++dpszahiidx);
	int dpszaloidx = std::max(dpszahiidx-1, 0);
	int dpsaahiidx;
	for(dpsaahiidx=0; dpsaahiidx<numdpsaas && dpsaas[dpsaahiidx]<=m_saas[m_saaidx]; ++dpsaahiidx);
	int dpsaaloidx = std::max(dpsaahiidx-1, 0);

	int numdps = 1;
	numdps = std::max( numdps, dptable[ dpszaloidx*numdpsaas + dpsaaloidx ] );
	numdps = std::max( numdps, dptable[ dpszaloidx*numdpsaas + dpsaahiidx ] );
	numdps = std::max( numdps, dptable[ dpszahiidx*numdpsaas + dpsaaloidx ] );
	numdps = std::max( numdps, dptable[ dpszahiidx*numdpsaas + dpsaahiidx ] );
	if( (numdps%2) == 0) numdps += 1;

	numdps = std::min(numdps, numdps_hardmax);

	int retval = 3+0*numdps;
	//std::printf(" %d", retval );
	return retval;
}

bool ProfilerMC::runHR( nx2dArray<double>& retrad, nx2dArray<skRTStokesVector>* retsvec, int pvorder, SKTRAN_HR_PolHOType pvhotype )
{
	bool ok = true;
	SKTRAN_HR_Engine					engine_hr;
	std::vector<SKTRAN_StokesScalar>    radiance_hr;
	SKTRAN_HR_Specs_User                specs_hr;

	SKTRAN_LineOfSightArray_V21			linesofsight;
	nx2dArray<double>					mc_radiance;
	SKTRAN_AtmosphericOpticalState_V21	opticalstate;
	nxVector							sun ( 1.0, 0, 0 );
	ok = ok && MakeSingleScatterComparisonLinesOfSight( linesofsight, sun );
    //sun.SetCoords(-8.71557427476582e-02, 0.00000000000000e+00, 9.96194698091746e-01);

	ok = ok && CreateOpticalState( &opticalstate );


	specs_hr.OpticalPropertiesSpecs().SetMaxPolarizationOrder ( pvorder );
	specs_hr.OpticalPropertiesSpecs().SetPolarizationHigherOrderBehaviour( pvhotype );
    //specs_hr.OpticalPropertiesSpecs().SetPolarizationOrder( 0 );
    //specs_hr.OpticalPropertiesSpecs().SetPolarizationHigherOrderBehaviour( SKTRAN_HR_PolHOType::unpolarized );

	specs_hr.OpticalPropertiesSpecs().SetAtmosphereHasDelta( SKTRAN_HR_AtmosphereHasDelta::no );

	specs_hr.IntegratorSpecs().SetMaxOpticalDepth( 0.1 );
	specs_hr.IntegratorSpecs().UseLegacySasktran21Technique( false );
	//specs_hr.RayTracingSpecs().SetShellSpacing( 20000.0 );
	//specs_hr.RayTracingSpecs().SetSolarShellSpacing( 20000.0 );
	specs_hr.RayTracingSpecs().SetSolarShellSpacing( 500.0 );
	//specs_hr.IntegratorSpecs().SetMaxOpticalDepth( 1 );
	//specs_hr.IntegratorSpecs().SetMaxExtinctionGradient( 0.9 );
	//specs_hr.RayTracingSpecs().SetLinesOfSightType( SKTRAN_HR_RayTracer_Curved );
	specs_hr.DiffuseSpecs().SetNumProfiles ( 1 + 0*3 + 0*SuggestNumDiffuseProfiles() );
    
    specs_hr.OpticalPropertiesSpecs().SetPrecacheWavel( m_wavs );

    //specs_hr.DiffuseSpecs().SetNumOutgoing( 324 );

	// hrmccomp default: 
	//     m_numbeforehoriz = 6;
	//     m_numhoriz = 8;
	//     m_numafterhoriz = 10;
	//     m_numazi = 12;
    // hrmccomp hires 
        //specs_hr.DiffuseSpecs().SetNumAzi( 27 ); // This is the only thing that's different... probably leftover from testing
        //specs_hr.DiffuseSpecs().SetNumAfterHoriz( 21 );
        //specs_hr.DiffuseSpecs().SetNumBeforeHoriz( 21 );
        //specs_hr.DiffuseSpecs().SetNumHoriz( 16 );
    // Diffuse sphere plots
	    int mfac = 1;
	    specs_hr.DiffuseSpecs().SetNumBeforeHoriz( mfac* 6 );
	    specs_hr.DiffuseSpecs().SetNumHoriz      ( mfac* 8 );
	    specs_hr.DiffuseSpecs().SetNumAfterHoriz ( mfac*10 );
	    specs_hr.DiffuseSpecs().SetNumAzi        ( mfac*12 );
    //// debugging
        //specs_hr.DiffuseSpecs().SetNumBeforeHoriz( 3 );
        //specs_hr.DiffuseSpecs().SetNumHoriz      ( 3 );
        //specs_hr.DiffuseSpecs().SetNumAfterHoriz ( 3 );
        //specs_hr.DiffuseSpecs().SetNumAzi        ( 3 );
        //specs_hr.DiffuseSpecs().SetHeightRes     ( 60000.0 );

    specs_hr.IntegratorSpecs().SetUseEmissions( true );
    specs_hr.IntegratorSpecs().SetUseSolarTransmission( true );
    specs_hr.RayTracingSpecs().SetTOAHeight( 120000.0 );
    specs_hr.DiffuseSpecs()   .SetMaxDiffuseHeight( 120000.0 );
    
	ok = ok && engine_hr.SetSun( sun );
	ok = ok && engine_hr.ConfigureModel( specs_hr, linesofsight, m_numThreads);
	radiance_hr.resize( linesofsight.NumRays() );

	nx2dArray<double> radiance;
	radiance.SetSize( linesofsight.NumRays(), 1/*wavs.size()*/ );
	radiance.SetTo(-1.0);

	std::vector<skRTStokesVector>* vecptr = nullptr;
	if(nullptr!=retsvec){
		vecptr = new std::vector<skRTStokesVector>;
	}


	for( size_t wavidx = 0; wavidx < m_wavs.size(); wavidx++ )
	{
        std::stringstream tempstr;
        tempstr.str("");
        auto precisiondefault = tempstr.precision();
        tempstr.precision(8);
        tempstr << namestr.back() << "_wav" << m_wavs[wavidx];
        tempstr.precision(precisiondefault);
        namestr.push_back(tempstr.str());
        tempstr.str("");
        tempstr << m_dirstr << "hr_" << namestr.back() << "_vec5.txt";
        FILE* f = fopen(tempstr.str().c_str(),"r");
        if( m_ignoreExistingFile || nullptr==f ){

        //std::cout << wavelen[i] << std::endl;
		std::stringstream cloudstream;
		//cloudstream << "C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/2013/randomDisc/AerosolProfileTropics/output/scattcloud_point/saa" << saas[saaidx];
		//engine_mc.SetSaveDir(cloudstream.str());
        
        ok = ok && engine_hr.CalculateRadiance(&radiance_hr, m_wavs[wavidx], m_maxOrderScatter, &opticalstate, vecptr );
		
		//hr_engine.CalculateSecondOrderRadiance(&hr_radiance_temp, wavelen[i],2,&opticalstate );
		for(size_t losidx = 0; losidx < std::min(radiance_hr.size(),radiance.XSize()); losidx++ )
		{
			radiance.At(losidx,0/*wavidx*/) = radiance_hr[losidx];
			//printf("%u:\t%1.10e\n", m_runctr++, radiance_mc[losidx]);
			//printf("%1.5e ", radiance_hr[losidx]); m_runctr++;
		}
		if( ok && nullptr!=vecptr && nullptr!=retsvec )
        {
			retsvec->SetSize(vecptr->size(),4);
			for(size_t losidx = 0; losidx < std::min(radiance_hr.size(),radiance.XSize()); ++losidx){
				retsvec->At(losidx,0) = vecptr->at(losidx);
			}
		}


            printf("%s", tempstr.str().c_str() );
        //retVal = runHR( rettemp, &svectemp, 5, SKTRAN_HR_PolHOType::unpolarized );
        if( !m_dirstr.empty() ){
            //retrad.push_back( rettemp );
            //retstr.push_back(namestr.back());
            if( ok ){
                //retsvec->push_back( svectemp );
                nx2dArray<double> toprint;
                toprint.SetSize(4, vecptr->size() );
                for( int xidx=0; xidx<vecptr->size(); xidx++){
                    toprint.At(0,xidx) = vecptr->at(xidx).At(1);
                    toprint.At(1,xidx) = vecptr->at(xidx).At(2);
                    toprint.At(2,xidx) = vecptr->at(xidx).At(3);
                    toprint.At(3,xidx) = vecptr->at(xidx).At(4);
                }
                if( m_writeToFile   ) toprint.WriteToTextFile(tempstr.str().c_str(),false,"%1.16e");
                if( m_writeToScreen ){
                    printf("\n");
                    for( int xidx=0; xidx< (int)vecptr->size(); ++xidx){
                        printf("%1.16e\t %1.16e\t %1.16e\n", toprint.At(0,xidx), toprint.At(1,xidx), toprint.At(2,xidx) );
                    }
                }
                printf(" ]/\n");
            } else{
                printf(" :(\n");
            }
            }
        }
        if(nullptr!=f) fclose(f);
        namestr.pop_back();
	}

	if(ok) retrad = radiance;

	delete vecptr; vecptr = nullptr;

	return ok;
}




int ProfilerMC::main( std::vector<nx2dArray<double>>& retrad,  std::vector<nx2dArray<skRTStokesVector>>* retsvec, std::vector<std::string>& retstr )
{
    int retVal = 0;
	m_runctr = 0;

	bool runMCflag   = true;
	bool runAnyHR    = true;
	bool runPV0flag  = runAnyHR && false;
    bool runPV1flag  = runAnyHR && false;
    bool runPV2flag  = runAnyHR && false;
    bool runPV5flag  = runAnyHR && true;
	bool runPVfullflag = runAnyHR && false;

	m_ignoreExistingFile = true;
	bool deleteFileMode = false; // Goes through and deletes files without doing any simulation

	m_writeToFile   = false;
	m_writeToScreen = true;

    std::stringstream sstream;
    std::stringstream dirstrm;
    std::string       lockfile;
	
	//std::timer::auto_cpu_timer t;
	
	nx2dArray<double> rettemp;
	nx2dArray<skRTStokesVector> svectemp;
	nx2dArray<double> toprint;

	std::stringstream tempstr;

	SetNumThreads( 7 );

	if(deleteFileMode){
		printf("Confirm, twice, that you want to delete data.\n");
        system("pause");
        system("pause");
	}

    for( m_albidx = 0; m_albidx < m_albs.size(); m_albidx++ )
    {
		tempstr.str("");
		tempstr << "a" << m_albs[m_albidx];
		namestr.push_back(tempstr.str());
        for( m_uo3idx = 0; m_uo3idx < m_useozone.size(); m_uo3idx++ )
        {
			tempstr.str("");
			tempstr << namestr.back() << "_o" << m_useozone[m_uo3idx];
			namestr.push_back(tempstr.str());
        for( m_uaeidx = 0; m_uaeidx < m_useaero.size(); m_uaeidx++ )
        {
			tempstr.str("");
			tempstr << namestr.back() << "_ae" << m_useaero[m_uaeidx];
			namestr.push_back(tempstr.str());
		for( m_uraidx = 0; m_uraidx < m_userayleigh.size(); m_uraidx++)
		{
			tempstr.str("");
			tempstr << namestr.back() << "_r" << m_userayleigh[m_uraidx];
			namestr.push_back(tempstr.str());
		for( m_latidx = 0; m_latidx < m_lats.size(); m_latidx++ )
		{
			tempstr.str("");
			tempstr << namestr.back() << "_lat" << m_lats[m_latidx]; 
			namestr.push_back(tempstr.str());
		for( m_saaidx = 0; m_saaidx < m_saas.size(); m_saaidx++ )
		{
			tempstr.str("");
			tempstr << namestr.back() << "_saa" << m_saas[m_saaidx];
			namestr.push_back(tempstr.str());
        //for( m_szaidx = m_saaidx; m_szaidx == m_saaidx; m_szaidx++ )
        for( m_szaidx=0; m_szaidx<m_szas.size(); ++m_szaidx )
		{
			tempstr.str("");
			tempstr << namestr.back() << "_sza" << m_szas[m_szaidx];
			namestr.push_back(tempstr.str());
		//for( m_wavidx = 0; m_wavidx < m_wavs.size(); m_wavidx++)
		//{
			//tempstr.str("");
   //         auto precisiondefault = tempstr.precision();
   //         tempstr.precision(8);
			//tempstr << namestr.back() << "_wav" << m_wavs[m_wavidx];
   //         tempstr.precision(precisiondefault);
			//namestr.push_back(tempstr.str());
		for( m_altidx = 0; m_altidx < 1 /*alts.size()*/; m_altidx++)
		{
			tempstr.str("");
			//tempstr << namestr.back() << "_alt" << alts[altidx];
            //tempstr << m_dirstr << namestr.back() << ".txt";
			//namestr.push_back(tempstr.str());
			//printf("%s\n",namestr.back().c_str());


            // Finish off filename string
            //tempstr << m_dirstr << namestr.back() << ".txt";



			if( !deleteFileMode ){
			// RUN MC SIMULATION
			if( runMCflag ){// || !(m_wavs[m_wavidx]<341.0 && 80.5<m_szas[m_szaidx]) ){
				tempstr.str("");
				tempstr << m_dirstr << "mc_" << namestr.back() << "_vec.txt";
				FILE* fmc = fopen(tempstr.str().c_str(),"r");
				if( m_ignoreExistingFile || nullptr==fmc ){
					//printf("%s", tempstr.str().c_str() );
                    nx2dArray<double> timing;
					retVal = runMC( rettemp, &svectemp, timing );
					if( !m_dirstr.empty() ){
						//printf("[_)");
						//rettemp.WriteToTextFile(tempstr.str().c_str(),false,"%1.16e");
						retrad.push_back( rettemp );
						retstr.push_back(namestr.back());
						if(NULL!=retsvec){
							retsvec->push_back( svectemp );
							toprint.SetSize(4,retsvec->back().XSize());
							for( int xidx=0; xidx< (int)retrad.back().XSize(); xidx++){
								toprint.At(0,xidx) = retsvec->back().At(xidx,0).At(1);
								toprint.At(1,xidx) = retsvec->back().At(xidx,0).At(2);
								toprint.At(2,xidx) = retsvec->back().At(xidx,0).At(3);
								toprint.At(3,xidx) = retsvec->back().At(xidx,0).At(4);
								//for( int yidx=0; yidx<retrad[rridx].YSize(); yidx++){
								//	std::printf("%1.16e\n",retrad[rridx].At(xidx,yidx));
								//}
							}

							if( m_writeToFile   ){
                                toprint.WriteToTextFile(tempstr.str().c_str(),false,"%1.16e");
                                timing.WriteToTextFile((tempstr.str()+"_timing").c_str(), false, "%1.16e");
                            }
							if( m_writeToScreen ){
								printf("\n");
								for( int xidx=0; xidx< (int)retrad.back().XSize(); ++xidx){
									printf("%1.16e\t %1.16e\t %1.16e\n", toprint.At(0,xidx), toprint.At(1,xidx), toprint.At(2,xidx) );
								}
							}
							//printf(" ]/\n");
						} else{
							printf(" :(\n");
						}
					}
				} else{
					fclose(fmc);
				}

			}

			// RUN HR SIMULATION -- pv0
			if( runPV0flag ){
				//std::timer::auto_cpu_timer t;
				tempstr.str("");
				tempstr << m_dirstr << "hr_" << namestr.back() << "_vec0.txt";
				FILE* fpv0 = fopen(tempstr.str().c_str(),"r");
				if( m_ignoreExistingFile || nullptr==fpv0 ){
					printf("%s", tempstr.str().c_str() );
					retVal = runHR( rettemp, &svectemp, 0, SKTRAN_HR_PolHOType::unpolarized );
					//return 0;
					if( !m_dirstr.empty() ){
						retrad.push_back( rettemp );
						retstr.push_back(namestr.back());
						if(NULL!=retsvec){
							retsvec->push_back( svectemp );
							toprint.SetSize(4,retsvec->back().XSize());
							for( int xidx=0; xidx< (int)retrad.back().XSize(); xidx++){
								toprint.At(0,xidx) = retsvec->back().At(xidx,0).At(1);
								toprint.At(1,xidx) = retsvec->back().At(xidx,0).At(2);
								toprint.At(2,xidx) = retsvec->back().At(xidx,0).At(3);
								toprint.At(3,xidx) = retsvec->back().At(xidx,0).At(4);
							}
							if( m_writeToFile   ) toprint.WriteToTextFile(tempstr.str().c_str(),false,"%1.16e");
							if( m_writeToScreen ){
								printf("\n");
								for( int xidx=0; xidx< (int)retrad.back().XSize(); ++xidx){
									printf("%1.16e\t %1.16e\t %1.16e\n", toprint.At(0,xidx), toprint.At(1,xidx), toprint.At(2,xidx) );
								}
							}
							printf(" ]/\n");
						} else{
							printf(" :(\n");
						}
					}
					if(nullptr!=fpv0){
						fclose( fpv0 );
					}
				} else{
					fclose(fpv0);
				}
			}


            // RUN HR SIMULATION -- pv1
            if( runPV1flag ){
                //std::timer::auto_cpu_timer t;
                tempstr.str("");
                tempstr << m_dirstr << "hr_" << namestr.back() << "_vec1.txt";
                FILE* fpv1 = fopen(tempstr.str().c_str(),"r");
                if( m_ignoreExistingFile || nullptr==fpv1 ){
                    printf("%s", tempstr.str().c_str() );
                    retVal = runHR( rettemp, &svectemp, 1, SKTRAN_HR_PolHOType::unpolarized );
                    //return 0;
                    if( !m_dirstr.empty() ){
                        retrad.push_back( rettemp );
                        retstr.push_back(namestr.back());
                        if(NULL!=retsvec){
                            retsvec->push_back( svectemp );
                            toprint.SetSize(4,retsvec->back().XSize());
                            for( int xidx=0; xidx< (int)retrad.back().XSize(); xidx++){
                                toprint.At(0,xidx) = retsvec->back().At(xidx,0).At(1);
                                toprint.At(1,xidx) = retsvec->back().At(xidx,0).At(2);
                                toprint.At(2,xidx) = retsvec->back().At(xidx,0).At(3);
                                toprint.At(3,xidx) = retsvec->back().At(xidx,0).At(4);
                            }
                            if( m_writeToFile   ) toprint.WriteToTextFile(tempstr.str().c_str(),false,"%1.16e");
                            if( m_writeToScreen ){
                                printf("\n");
                                for( int xidx=0; xidx< (int)retrad.back().XSize(); ++xidx){
                                    printf("%1.16e\t %1.16e\t %1.16e\n", toprint.At(0,xidx), toprint.At(1,xidx), toprint.At(2,xidx) );
                                }
                            }
                            printf(" ]/\n");
                        } else{
                            printf(" :(\n");
                        }
                    }
                    if(nullptr!=fpv1){
                        fclose( fpv1 );
                    }
                } else{
                    fclose(fpv1);
                }
            }

            // RUN HR SIMULATION -- pv2
            if( runPV2flag ){
                //std::timer::auto_cpu_timer t;
                tempstr.str("");
                tempstr << m_dirstr << "hr_" << namestr.back() << "_vec2.txt";
                FILE* fpv2 = fopen(tempstr.str().c_str(),"r");
                if( m_ignoreExistingFile || nullptr==fpv2 ){
                    printf("%s", tempstr.str().c_str() );
                    retVal = runHR( rettemp, &svectemp, 2, SKTRAN_HR_PolHOType::unpolarized );
                    //return 0;
                    if( !m_dirstr.empty() ){
                        retrad.push_back( rettemp );
                        retstr.push_back(namestr.back());
                        if(NULL!=retsvec){
                            retsvec->push_back( svectemp );
                            toprint.SetSize(4,retsvec->back().XSize());
                            for( int xidx=0; xidx< (int)retrad.back().XSize(); xidx++){
                                toprint.At(0,xidx) = retsvec->back().At(xidx,0).At(1);
                                toprint.At(1,xidx) = retsvec->back().At(xidx,0).At(2);
                                toprint.At(2,xidx) = retsvec->back().At(xidx,0).At(3);
                                toprint.At(3,xidx) = retsvec->back().At(xidx,0).At(4);
                            }
                            if( m_writeToFile   ) toprint.WriteToTextFile(tempstr.str().c_str(),false,"%1.16e");
                            if( m_writeToScreen ){
                                printf("\n");
                                for( int xidx=0; xidx< (int)retrad.back().XSize(); ++xidx){
                                    printf("%1.16e\t %1.16e\t %1.16e\n", toprint.At(0,xidx), toprint.At(1,xidx), toprint.At(2,xidx) );
                                }
                            }
                            printf(" ]/\n");
                        } else{
                            printf(" :(\n");
                        }
                    }
                    if(nullptr!=fpv2){
                        fclose( fpv2 );
                    }
                } else{
                    fclose(fpv2);
                }
            }



            // RUN HR SIMULATION -- pv5
            if( runPV5flag ){
                //std::timer::auto_cpu_timer t;
                //if( m_ignoreExistingFile || nullptr==fpv5 ){
                    retVal = runHR( rettemp, &svectemp, 5, SKTRAN_HR_PolHOType::unpolarized );
            //        
            //        if(nullptr!=fpv5){
            //            fclose( fpv5 );
            //        }
            //    } else{
            //        fclose(fpv5);
            //    }
            }



			if( runPVfullflag ){
				//std::timer::auto_cpu_timer t;
				// RUN HR SIMULATION -- pvfull
				tempstr.str("");
				tempstr << m_dirstr << "hr_" << namestr.back() << "_vecfull.txt";
				FILE* fpvfull = fopen(tempstr.str().c_str(),"r");
				if( m_ignoreExistingFile || nullptr==fpvfull ){
					printf("%s", tempstr.str().c_str() );
					retVal = runHR( rettemp, &svectemp, (int)m_maxOrderScatter, SKTRAN_HR_PolHOType::unpolarized );
					//return 0;
					if( !m_dirstr.empty() ){
						retrad.push_back( rettemp );
						retstr.push_back(namestr.back());
						if(NULL!=retsvec){
							retsvec->push_back( svectemp );
							toprint.SetSize(4,retsvec->back().XSize());
							for( int xidx=0; xidx< (int)retrad.back().XSize(); xidx++){
								toprint.At(0,xidx) = retsvec->back().At(xidx,0).At(1);
								toprint.At(1,xidx) = retsvec->back().At(xidx,0).At(2);
								toprint.At(2,xidx) = retsvec->back().At(xidx,0).At(3);
								toprint.At(3,xidx) = retsvec->back().At(xidx,0).At(4);
							}
							if( m_writeToFile   ) toprint.WriteToTextFile(tempstr.str().c_str(),false,"%1.16e");
							if( m_writeToScreen ){
								printf("\n");
								for( int xidx=0; xidx< (int)retrad.back().XSize(); ++xidx){
									printf("%1.16e\t %1.16e\t %1.16e\n", toprint.At(0,xidx), toprint.At(1,xidx), toprint.At(2,xidx) );
								}
							}
							printf(" ]/\n");
						} else{
							printf(" :(\n");
						}
					}
					if(nullptr!=fpvfull){
						fclose( fpvfull );
					}
				} else{
					fclose(fpvfull);
				}
			}
			} else{
				
				if( runMCflag ){
					tempstr.str("");
					tempstr << m_dirstr << "mc_" << namestr.back() << "_vec.txt";
					printf("%s deleted\n", tempstr.str().c_str());
					remove( tempstr.str().c_str() );
				}
				if( runPV0flag ){
					tempstr.str("");
					tempstr << m_dirstr << "hr_" << namestr.back() << "_vec0.txt";
					printf("%s deleted\n", tempstr.str().c_str());
					remove( tempstr.str().c_str() );
				}
				if( runPV1flag ){
					tempstr.str("");
					tempstr << m_dirstr << "hr_" << namestr.back() << "_vec1.txt";
					printf("%s deleted\n", tempstr.str().c_str());
					remove( tempstr.str().c_str() );
				}
				if( runPVfullflag ){
					tempstr.str("");
					tempstr << m_dirstr << "hr_" << namestr.back() << "_vecfull.txt";
					printf("%s deleted\n", tempstr.str().c_str());
					remove( tempstr.str().c_str() );
				}
			}
		} //namestr.pop_back();
		//}
        namestr.pop_back();
		} namestr.pop_back();
		} namestr.pop_back();
		} namestr.pop_back();
        } namestr.pop_back();
        } namestr.pop_back();
        } namestr.pop_back();
    }
    //PrintOptions();
    
	namestr.clear();

    printf("\n\n");
	//system("pause");
    return 0;
}
