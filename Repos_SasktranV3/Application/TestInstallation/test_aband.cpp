#include <sasktran.h>
#if defined (NX_WINDOWS)
#include "include/test_aband.h"
#include <algorithm>
#include <boost/timer/timer.hpp>


SKTRAN_ABand_Test::SKTRAN_ABand_Test()
{
    SetupForTest();
}


SKTRAN_ABand_Test::~SKTRAN_ABand_Test()
{

}

void SKTRAN_ABand_Test::SetupForTest()
{
   	m_hardcodefolder   = "../../../../Application/TestInstallation/hardcodetestvalues/";
    m_emissionFileName = "abandTestVER";

    m_tolerance = 1e-10;
    m_toaHeight = 119000.0;

    m_writeToFile = false;
    m_writeErrorsToScreen = true;

    m_geometryIndex = 0;
    m_precisionMC = 0.01;

    // Open-MP as implemented in Visual Studio Version 14.0.25123.00 Update 2, does not 
    // support runtime-generated chunk sizes. These are needed to make sure multithreaded 
    // MC runs a known number of ray histories on each thread. The multithreaded-ness of 
    // MC can't be tested quickly until this is fixed. 
    m_numThreadsMC = 1;
    m_numThreadsHR = 7; 

    m_maxOrderScatter = 10;
    m_seed = 1234;

    m_wavs.resize(0);
    m_wavs.push_back(7.63313620759077e+02);
    m_wavs.push_back(7.63314203407206e+02);
    m_wavs.push_back(7.63315077381066e+02);

    m_alts.resize(0);
    //for(int repidx=11000; repidx<120000; repidx+=2000 ){
    //    m_alts.push_back( repidx );
    //}   
    m_alts.push_back(  20000.0 );
    m_alts.push_back(  50000.0 );
    m_alts.push_back(  80000.0 );
    m_alts.push_back( 110000.0 );
    
    m_szas.resize(0);
    m_szas.push_back(60.0);
    m_szas.push_back(89.0);

    m_saas.resize(0);
    m_saas.push_back( 30.0 );
    m_saas.push_back( 90.0 );
}


bool SKTRAN_ABand_Test::MakeSingleScatterComparisonLinesOfSight( SKTRAN_LineOfSightArray_V21& linesofsight, nxVector& sunvec )
{
    bool     ok = true;
    double    mjd     = 54832.0;

    double* alist = new double[m_alts.size()];
    for(int aidx = 0; aidx< (int)m_alts.size(); aidx++){
        alist[aidx] = m_alts[aidx];
    }
    ok = ok && linesofsight.SetRaysFromTangentHeightArray(54000, 75.0, 0.0, m_szas[m_geometryIndex], m_saas[m_geometryIndex], 0.0, alist, (int) m_alts.size(), 600000.0, &sunvec);
    
    delete[] alist;

    return ok;
}


bool SKTRAN_ABand_Test::CreateOpticalState( SKTRAN_AtmosphericOpticalState_V21* opticalstate )
{
    bool                                       ok = true;

    skOpticalProperties_RayleighDryAir*        rayleigh;                    // Optical properties of one air molecule
    skClimatology_MSIS90*                      msis90;    
    skClimatology_LabowOzoneVMR*               o3numberdensity;
    skOpticalProperties_O3_OSIRISRes*          o3_opticalprops;            // optical properties of one O3 molecule

    ok = ok && NULL!=opticalstate;
    if(ok) opticalstate->erase();
    
    msis90           = new skClimatology_MSIS90;
    rayleigh         = new skOpticalProperties_RayleighDryAir;
    ok = ok && (NULL!=msis90) && (NULL!=rayleigh);
    ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, msis90, rayleigh );
    
    o3numberdensity  = new skClimatology_LabowOzoneVMR;
    o3_opticalprops  = new skOpticalProperties_O3_OSIRISRes;
    ok = ok && (NULL!=o3numberdensity) && (NULL!=o3_opticalprops);
    ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_O3_CM3, o3numberdensity, o3_opticalprops);
    
    skOpticalProperties_HitranChemical* o2 = new skOpticalProperties_HitranChemical("O2");
	o2->SetWavenumberRange(12839.9373411058, 13350.0654153205 );
    ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_O2_CM3, msis90, o2 );
    
    skEmission_Tabulated_HeightWavelength* emission = new skEmission_Tabulated_HeightWavelength; 
    ok = ok && emission->LoadHeightWavelengthProfileFromFile( (m_hardcodefolder + m_emissionFileName).c_str() );
    ok = ok && opticalstate->AddEmission( SKEMISSION_PHOTOCHEMICAL_0, emission );
    
    ok = ok && opticalstate->SetAlbedo( 0.3 );

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
        }
    }

    return ok;
}


bool SKTRAN_ABand_Test::runMC( )
{
    bool ok = true;

    SKTRAN_Engine_MC_V21               engine_mc;
    std::vector<SKTRAN_StokesScalar>   radiance_mc;
	std::vector<skRTStokesVector>      vradiance_mc;
	nx2dArray<skRTStokesVector>        referenceValues;
    SKTRAN_Specifications_MC           specs_mc;
    nx1dArray<skRTStokesVector>        tempsvec;

    SKTRAN_LineOfSightArray_V21        linesofsight;
    nx2dArray<double>                  mc_radiance;
    SKTRAN_AtmosphericOpticalState_V21 opticalstate;
    nxVector                           sun ( 1.0, 0, 0 );

	bool errorIsFound = false;


    ok = ok && MakeSingleScatterComparisonLinesOfSight ( linesofsight, sun );
    ok = ok && CreateOpticalState                      ( &opticalstate );

    ok = ok && specs_mc.ConfigureDefaults        ( );
    ok = ok && specs_mc.SetSunGeographicPosition ( sun );
    ok = ok && specs_mc.SetSineSunApexAngle      ( nxmath::sind(0.268) ); // Angle from sun center to sun edge at Earth
    ok = ok && specs_mc.SetScatterPositionRes    ( 10.0 );
    ok = ok && specs_mc.SetAllowDynamicThreads   ( false );

    ok = ok && specs_mc.SetMinFractionHigherOrder( 0.1 );
    ok = ok && specs_mc.SetMinimumRelPathWeight  ( 0.0*m_precisionMC*1e-3 );
    
    ok = ok && specs_mc.SetPrecisionMC           ( m_precisionMC );
    ok = ok && specs_mc.SetNumPhotonsPerLOS      ( true ? (size_t)(3.0/(m_precisionMC*m_precisionMC)) : 1 );
    //ok = ok && specs_mc.SetChunkSize             ( 1 );
    ok = ok && specs_mc.SetRngSeed               ( m_seed );
    ok = ok && specs_mc.SetSunType               ( SKTRAN_Specifications_MC::SunType::point );
    ok = ok && specs_mc.SetSolarRayTracerType    ( SKTRAN_Specifications_MC::RayTracerType::shell );
    ok = ok && specs_mc.SetMSRayTracerType       ( SKTRAN_Specifications_MC::RayTracerType::shell );
    ok = ok && specs_mc.SetLOSRayTracerType      ( SKTRAN_Specifications_MC::RayTracerType::shell );
    ok = ok && specs_mc.SetOptPropIntType        ( SKTRAN_Specifications_MC::OptPropIntType::straight );
    ok = ok && specs_mc.SetKernelType            ( SKTRAN_Specifications_MC::LogType::none );
    ok = ok && specs_mc.SetPolType                 ( SKTRAN_Specifications_MC::PolType::pol /*m_polType */);
    
    ok = ok && specs_mc.SetSolarTableType        ( SKTRAN_Specifications_MC::SolarTableType::doNothing );
    ok = ok && specs_mc.SetEmissionTableType     ( SKTRAN_Specifications_MC::EmissionTableType::dim1 );
    ok = ok && specs_mc.SetTOAHeight( m_toaHeight );

    ok = ok && specs_mc.SetAdaptOptDepthMax      ( 0.05 );
    ok = ok && specs_mc.FinalizeSpecs            ( );

    ok = ok && engine_mc.ConfigureModel( specs_mc, linesofsight, m_numThreadsMC);
    radiance_mc.resize( linesofsight.NumRays() );
	vradiance_mc.resize( linesofsight.NumRays() );

	nx2dArray<double> radiance;
	radiance.SetSize( linesofsight.NumRays(), m_wavs.size() );
	radiance.SetTo(-1.0);
	nx2dArray<skRTStokesVector> vradiance;
	vradiance.SetSize( linesofsight.NumRays(), m_wavs.size() );

    for( size_t wavidx = 0; wavidx < m_wavs.size(); wavidx++ )
    {
		ok = ok && engine_mc.CalculateRadiance(&radiance_mc, m_wavs[wavidx], m_maxOrderScatter, &opticalstate, &vradiance_mc );
		for(size_t losidx = 0; losidx < min(radiance_mc.size(),radiance.XSize()); losidx++ )
		{
			vradiance.At( losidx, wavidx ) = vradiance_mc[losidx];
		}
    }

    for(size_t losidx = 0; losidx < min(radiance_mc.size(),radiance.XSize()); losidx++ )
    {
        radiance.At(losidx,0/*wavidx*/) = radiance_mc[losidx];
    }


	std::stringstream hardcodefilename;
	hardcodefilename << m_hardcodefolder << "abandTest_MC_geo" << m_geometryIndex;

	referenceValues.SetSize( linesofsight.NumRays(), m_wavs.size() );
	if( m_writeToFile ) ok = ok && WriteReferenceValuesToFile( hardcodefilename.str(), vradiance );
	ok = ok && LoadReferenceValues( hardcodefilename.str(), referenceValues ); 
	errorIsFound = errorIsFound | CompareValues( vradiance, referenceValues, m_writeErrorsToScreen );

    if(ok){
	    if(errorIsFound) 
		    printf("\nAt least one MC-calculated stokes vector had a relative difference of at least %e compared to the reference value.\n\n\n", m_tolerance );
    } else{
        printf("\nAn error occurred while running SKTRAN_ABand_Test::runMC. This was caught by a bad 'ok' flag.");
    }
   
    return errorIsFound;
}


bool SKTRAN_ABand_Test::runHR( int pvorder, SKTRAN_HR_PolHOType pvhotype )
{
	bool ok = true;

	SKTRAN_HR_Engine                   engine_hr;
	std::vector<SKTRAN_StokesScalar>   radiance_hr;
	std::vector<skRTStokesVector>      vradiance_hr;
	nx2dArray<skRTStokesVector>        referenceValues;
	SKTRAN_HR_Specs_User               specs_hr;

	SKTRAN_LineOfSightArray_V21        linesofsight;
	SKTRAN_AtmosphericOpticalState_V21 opticalstate;
	nxVector                           sun ( 1.0, 0, 0 );

	bool errorIsFound = false;

	ok = ok && MakeSingleScatterComparisonLinesOfSight( linesofsight, sun );
	ok = ok && CreateOpticalState( &opticalstate );

	specs_hr.DiffuseSpecs().SetNumProfiles    ( 3 );
	specs_hr.DiffuseSpecs().SetNumBeforeHoriz ( 6 );
	specs_hr.DiffuseSpecs().SetNumHoriz       ( 8 );
	specs_hr.DiffuseSpecs().SetNumAfterHoriz  ( 10 );
	specs_hr.DiffuseSpecs().SetNumAzi         ( 12 );
    
	specs_hr.OpticalPropertiesSpecs().SetAtmosphereHasDelta( SKTRAN_HR_AtmosphereHasDelta::no );
	specs_hr.IntegratorSpecs().UseLegacySasktran21Technique( false );
	specs_hr.IntegratorSpecs().SetUseEmissions( true );
	specs_hr.IntegratorSpecs().SetUseSolarTransmission( false );
    specs_hr.IntegratorSpecs().SetMaxOpticalDepth( 0.02 );
    specs_hr.RayTracingSpecs().SetTOAHeight( m_toaHeight );
    specs_hr.DiffuseSpecs()   .SetMaxDiffuseHeight( m_toaHeight );

	specs_hr.OpticalPropertiesSpecs().SetMaxPolarizationOrder ( pvorder );
	specs_hr.OpticalPropertiesSpecs().SetPolarizationHigherOrderBehaviour( pvhotype );

	ok = ok && engine_hr.SetSun( sun );
	ok = ok && engine_hr.ConfigureModel( specs_hr, linesofsight, m_numThreadsHR);
	radiance_hr.resize( linesofsight.NumRays() );

	nx2dArray<double> radiance;
	radiance.SetSize( linesofsight.NumRays(), m_wavs.size() );
	radiance.SetTo(-1.0);
	nx2dArray<skRTStokesVector> vradiance;
	vradiance.SetSize( linesofsight.NumRays(), m_wavs.size() );

	for( size_t wavidx = 0; wavidx < m_wavs.size(); wavidx++ )
	{
		ok = ok && engine_hr.CalculateRadiance(&radiance_hr, m_wavs[wavidx], m_maxOrderScatter, &opticalstate, &vradiance_hr );
		for(size_t losidx = 0; losidx < min(radiance_hr.size(),radiance.XSize()); losidx++ )
		{
			vradiance.At( losidx, wavidx ) = vradiance_hr[losidx];
		}
	}


	std::stringstream hardcodefilename;
	hardcodefilename << m_hardcodefolder << "abandTest_HR_pv" << pvorder << "_geo" << m_geometryIndex;

	referenceValues.SetSize( linesofsight.NumRays(), m_wavs.size() );
	if( m_writeToFile ) ok = ok && WriteReferenceValuesToFile( hardcodefilename.str(), vradiance );
	ok = ok && LoadReferenceValues( hardcodefilename.str(), referenceValues ); 
	errorIsFound = errorIsFound | CompareValues( vradiance, referenceValues, m_writeErrorsToScreen );

    if(ok){
	    if(errorIsFound) 
		    printf("\nAt least one HR-calculated stokes vector had a relative difference of at least %e compared to the reference value for polarization order %d.\n\n\n", m_tolerance, pvorder);
    } else{
        printf("\nAn error occurred while running SKTRAN_ABand_Test::runHR. This was caught by a bad 'ok' flag.");
    }
	
	return errorIsFound;
}


bool SKTRAN_ABand_Test::LoadReferenceValues( std::string& filename, nx2dArray<skRTStokesVector>& referenceValues ) const 
{
    bool ok = true;

    ifstream fstream;
    nx2dArray<double> tempI, tempQ, tempU, tempV;
    tempI.SetSize( referenceValues.XSize(), referenceValues.YSize() );
    tempQ.SetSize( referenceValues.XSize(), referenceValues.YSize() );
    tempU.SetSize( referenceValues.XSize(), referenceValues.YSize() );
    tempV.SetSize( referenceValues.XSize(), referenceValues.YSize() );
    
    fstream.open((filename+"I").c_str(), ios_base::in);
    ok = ok && fstream.is_open();
    if(ok) fstream >> tempI;
    fstream.close();
    
    fstream.open((filename+"Q").c_str(), ios_base::in);
    ok = ok && fstream.is_open();
    if(ok) fstream >> tempQ;
    fstream.close();
    
    fstream.open((filename+"U").c_str(), ios_base::in);
    ok = ok && fstream.is_open();
    if(ok) fstream >> tempU;
    fstream.close();
    
    fstream.open((filename+"V").c_str(), ios_base::in);
    ok = ok && fstream.is_open();
    if(ok) fstream >> tempV;
    fstream.close();

    if(ok){
        for( int yidx=0; yidx<referenceValues.YSize(); ++yidx )
            for( int xidx=0; xidx<referenceValues.XSize(); ++xidx )
                referenceValues.At(xidx,yidx).SetTo( tempI.At(xidx,yidx), tempQ.At(xidx,yidx), tempU.At(xidx,yidx), tempV.At(xidx,yidx) );
    } else{
        referenceValues.SetSize(0,0);
        nxLog::Record(NXLOG_WARNING, "SKTRAN_ABand_Test::LoadReferenceValues, Could not read reference values from file.");
    }

    return ok;
}


bool SKTRAN_ABand_Test::WriteReferenceValuesToFile( std::string& filename, const nx2dArray<skRTStokesVector>& newReferenceValues ) const
{
    bool ok = true;

    nx2dArray<double> tempI, tempQ, tempU, tempV;
    tempI.SetSize( newReferenceValues.XSize(), newReferenceValues.YSize() );
    tempQ.SetSize( newReferenceValues.XSize(), newReferenceValues.YSize() );
    tempU.SetSize( newReferenceValues.XSize(), newReferenceValues.YSize() );
    tempV.SetSize( newReferenceValues.XSize(), newReferenceValues.YSize() );
    
    for( int yidx=0; yidx<newReferenceValues.YSize(); ++yidx ){
        for( int xidx=0; xidx<newReferenceValues.XSize(); ++xidx ){
			tempI.At(xidx,yidx) = newReferenceValues.At(xidx,yidx).At(1);
			tempQ.At(xidx,yidx) = newReferenceValues.At(xidx,yidx).At(2);
			tempU.At(xidx,yidx) = newReferenceValues.At(xidx,yidx).At(3);
			tempV.At(xidx,yidx) = newReferenceValues.At(xidx,yidx).At(4);
        }
    }
    
    ok = ok && tempI.WriteToTextFile( (filename+"I").c_str(), false, "%1.17e" );
    ok = ok && tempQ.WriteToTextFile( (filename+"Q").c_str(), false, "%1.17e" );
    ok = ok && tempU.WriteToTextFile( (filename+"U").c_str(), false, "%1.17e" );
    ok = ok && tempV.WriteToTextFile( (filename+"V").c_str(), false, "%1.17e" );

    if(!ok) nxLog::Record(NXLOG_WARNING, "SKTRAN_ABand_Test::WriteReferenceValues, Could not write new reference values to file.");

    return ok;
}


bool SKTRAN_ABand_Test::CompareValues( const nx2dArray<skRTStokesVector>& newValues, const nx2dArray<skRTStokesVector>& referenceValues, bool writeErrorsToScreen ) const
{
    bool errorIsFound = false;

    errorIsFound = errorIsFound || (newValues.XSize() != referenceValues.XSize());
    errorIsFound = errorIsFound || (newValues.YSize() != referenceValues.YSize());

	if( !errorIsFound ){
        for( int yidx=0; yidx<(int)newValues.YSize(); ++yidx){
            for( int xidx=0; xidx<(int)newValues.XSize(); ++xidx){
                skRTStokesVector difference = referenceValues.At(xidx,yidx) - newValues.At(xidx,yidx);
                skRTStokesVector mean = (referenceValues.At(xidx,yidx) + newValues.At(xidx,yidx)) * 0.5;
                double relErrorI = abs(mean.At(1))>0.0 ? abs(difference.At(1)/mean.At(1)) : 0.0;
                double relErrorQ = abs(mean.At(2))>0.0 ? abs(difference.At(2)/mean.At(2)) : 0.0;
                double relErrorU = abs(mean.At(3))>0.0 ? abs(difference.At(3)/mean.At(3)) : 0.0;
                double relErrorV = abs(mean.At(4))>0.0 ? abs(difference.At(4)/mean.At(4)) : 0.0;
                double sumErrors = relErrorI + relErrorQ + relErrorU + relErrorV;
                bool iserror = sumErrors > m_tolerance;
	            if(iserror){
					printf("Stokes vec rel error: %1.16e %1.16e %1.16e %1.16e\n", difference.At(1), difference.At(2), difference.At(3), difference.At(4) );
					errorIsFound = true;
                }
            }
        }
    }

    return errorIsFound;
}


int SKTRAN_ABand_Test::main(  )
{
	bool errorDetected = false;
    
    printf("\n------------------------------------*------------------------------------");
    printf("\nThis program will test the emission components of the MC and HR engines.");
    printf("\nIt will take about a minute.");
    printf("\n-------------------------------------------------------------------------");
    printf("\n");


    for( m_geometryIndex=0; m_geometryIndex<m_saas.size(); ++m_geometryIndex )
    {
        errorDetected = errorDetected | runMC( );
		errorDetected = errorDetected | runHR( 5, SKTRAN_HR_PolHOType::unpolarized );
    }

    
    printf("\n\n");
    printf("&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
    printf("$                        $\n");
	if( errorDetected ){
    printf("&  Emission test FAILED  &\n");
	}else {
    printf("&  Emission test PASSED  &\n");
    }
    printf("$                        $\n");
    printf("&&&&&&&&&&&&&&&&&&&&&&&&&&\n");
    printf("\n\n");


    return 0;
}
#endif
