#include <sasktran.h>

static bool MakeLinesOfSightFromArray( SKTRAN_LineOfSightArray_V21& linesofsight )
{
	bool 	ok;
	double	mjd 	= 54832.0;
	double  lat 	= 0.0;
	double	lng 	= 0.0;
	double 	sza 	= 60.0;
	double  saa 	= 0.0;
	double  rayazi 	= 270.0;
	double	tanheights_meters[36];
	
	for (int i=0; i < 36; i++) tanheights_meters[i] = i*2000.0 + 10000.0;	
	ok = linesofsight.SetRaysFromTangentHeightArray( mjd, lat, lng, sza, saa, rayazi, tanheights_meters, 36, 600000.0, NULL );

	return ok;
}



/*-----------------------------------------------------------------------------
 *					MakeSingleScatterLinesOfSight		2010-7-28*/
/** I generated this test data to diagnose/optimize the single scatter calculations.
 *	The OSIRIS processing uses single scatter during the Flat-field absolute
 *	calibration. I have taken 5 lines of sight from Odin-OSiris scan 34746056.
 *	The 5 lines of sight correspond to tangent altitudes 45, 50, 55, 60, and 65 km.
 *	The sun is at a sza angle around 97.3 degrees. 
 **/
/*---------------------------------------------------------------------------*/
static bool MakeSingleScatterLinesOfSight(	SKTRAN_LineOfSightArray_V21*	linesofsight )
{

	// Odin-osiris scannumber = 34746056, 5 lines of sight for 45,50,55,60,65  km

	double observerarray [15] = 
	{
	     -2520207.2578727026,      -6301346.9729181193,      -1549481.7901466270,
		 -2505836.1230591196,      -6294401.9933427088,      -1600406.6631476590,
	     -2491181.3525196696,      -6286866.0002852920,      -1652039.6527833764,
	     -2476510.4623296079,      -6279155.8510839343,      -1702856.0304066506,
	     -2462084.9387370660,      -6271351.8349050386,      -1752170.4555878886
	};

	double  lookarray [15]= 
	{
	   -0.054815521705618,        0.207906915766870,        0.976611474925647,
	   -0.058830891772998,        0.199408914700985,        0.978148767269564,
	   -0.062918988864250,        0.190767320555488,        0.979616777239130,
	   -0.066960426444567,        0.182221675478060,        0.980974802059834,
	   -0.070704457668028,        0.173768449263020,        0.982245084338219
	};
	
	double  mjdarray [5]=  { 54290.0990619432, 54290.0989812359, 54290.0988992517, 54290.0988184148, 54290.0987398143 };

	nxVector	observer;
	nxVector	look;
	double		mjd;
	size_t		i;
	size_t		ix;
	size_t		iy;
	size_t		iz;
	bool		ok = true;

	for (i = 0; i < 5; i++)
	{
		ix = 3*i;
		iy = ix + 1;
		iz = ix + 2;
		observer.SetCoords(  observerarray[ix], observerarray[iy], observerarray[iz]  );
		look.SetCoords    (  lookarray[ix],     lookarray[iy],     lookarray[iz]      );
		mjd = mjdarray[i];
		ok = ok && linesofsight->AddLineOfSight( observer, look, mjd );						// For the test we add just one line of sight.
	}
	return ok;
}





/*-----------------------------------------------------------------------------
 *					MakeOffOrbitLinesOfSight		2010-7-28*/
/** **/
/*---------------------------------------------------------------------------*/

static void MakeOffOrbitLinesOfSight( nxVector* satpos, nxVector* lookv, nxVector* sun, double tangenth)
{
	nxGeodetic	geoid;

	nxVector	west;
	nxVector	south;
	nxVector	up;
	nxVector	horiz;
	double		offorbitangle= 17; 

	sun->SetCoords					( 0.0,  -nxmath::cosd(15.0), nxmath::cosd(15.0) );
	geoid.FromGeodetic				( 70.0, 0.0, 600000.0 );				// Get the satellite at 82 degrees north
	*satpos = geoid.Location();												// and at 600 km altitude
	geoid.GetGeodeticWestSouthUp	( &west, &south, &up);					// Get the topocentric unit vectors at the current location

	horiz = (-nxmath::sind(offorbitangle-7.0))*west  + (-nxmath::cosd(offorbitangle-7.0))*south;	// get an angle east of due north corresponding to the off "orbit" angle (orbit inclination = 97.0 )

	geoid.FromTangentAltitude			( tangenth, *satpos, horiz, lookv);							// Get the look vector at off orbit angle bearing that looks at 23 km altitude.
	printf("%f  %f %f\n", (double)geoid.GeodeticLatitude(), geoid.GeodeticLongitude(), geoid.Height() );
}
	


/*-----------------------------------------------------------------------------
 *					MakeSpeciesList		2008-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

static bool MakeOpticalState( SKTRAN_AtmosphericOpticalState_V21* opticalstate )
{
	skOpticalProperties_RayleighDryAir*			rayleigh;					// Optical properties of one air molecule
	skClimatology_MSIS90*					msis90;	
	skClimatology_LabowOzoneVMR*			o3numberdensity;
	skOpticalProperties_O3_OSIRISRes*			o3_opticalprops;			// optical properties of one O3 molecule
//	skClimatology_Ecmwf*					ecmwf;	
//	skRTExtinction_MieAerosol*				mieAerosol;					// optical properties of one mie aerosol molecule
//	skRTExtinction_AerosolProfile*			aerosolProfile;	
	bool									ok;

	rayleigh         = new skOpticalProperties_RayleighDryAir;
	msis90           = new skClimatology_MSIS90;
	o3numberdensity  = new skClimatology_LabowOzoneVMR;
	o3_opticalprops  = new skOpticalProperties_O3_OSIRISRes;
//	mieAerosol		 = new skRTExtinction_MieAerosol;
//	aerosolProfile   = new skRTExtinction_AerosolProfile;
//	ecmwf			 = new skClimatology_Ecmwf;

	opticalstate->erase();
	ok =       (rayleigh != NULL) && (msis90 != NULL) && (o3_opticalprops != NULL) && (o3numberdensity != NULL);

	ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3,	msis90,				rayleigh );
	ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_PRESSURE_PA,				msis90,				rayleigh );
	ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_TEMPERATURE_K,			msis90,				rayleigh );
	ok = ok && opticalstate->AddSpecies( SKCLIMATOLOGY_O3_CM3,					o3numberdensity,	o3_opticalprops);
	ok = ok && opticalstate->SetAlbedo( 0.2 );		// set the albedo to 0.2

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"MakeSpeciesList, there was an error making the species list");
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					ConfigureSpecifications		2011-07-12*/
/** 
*	This function sets up all the geometry necessary for use in the
*	Monte Carlo ray-tracing class, MC_RayGeometry_Array.  The sun and
*	mjd match up with those used in MakeLinesOfSightFromArray
**/
/*---------------------------------------------------------------------------*/
static bool ConfigureSpecifications( SKTRANSO_SpecificationsUser_Legacy& specs, const SKTRAN_LineOfSightArray_V21& linesofsight )
{
	SKTRANSO_Engine		engine;
	bool					ok;
	double					mjd = 54832.0;
	nxVector				sun ( 4.999999998815835e-001, 0.000000000000000e+000, 8.660254038528065e-001 );
	
	ok = 	   engine.ConfigureModel( specs, linesofsight, linesofsight.NumRays() );
	ok = ok && specs.UpdateUndefinedParametersFromLinesOfSight( linesofsight );
	ok = ok && specs.RayTracingRegionManagerVar()->SetSun( sun );	// explicitly set the sun unit vector
	
	return ok;
}



/*-----------------------------------------------------------------------------
 *					TestSingleScatter		2010-7-28*/
/** **/
/*---------------------------------------------------------------------------*/

void TestSingleScatter()
{

	SKTRANSO_Engine						sasktranengine;
	SKTRAN_AtmosphericOpticalState_V21		opticalstate;
	SKTRANSO_SpecificationsUser_Legacy		specs;
	SKTRAN_LineOfSightArray_V21				linesofsight;


	nx2dArray<double>						radiance;
	nxTimeStamp								tstart;
	nxTimeStamp								tfinish;
	std::vector<SKTRAN_StokesScalar>		losradiance;
	size_t									numwavelen = 1353;
	bool									ok;
	double									wavelen;


	tstart.FromSystem();
	ok =        MakeOpticalState( &opticalstate );
//	ok  = ok && MakeSingleScatterLinesOfSight( &linesofsight );
	ok  = ok && MakeLinesOfSightFromArray( linesofsight );
	ok  = ok && ConfigureSpecifications( specs, linesofsight );
	ok  = ok && sasktranengine.ConfigureModel( specs, linesofsight, linesofsight.NumRays() );
	ok  = ok && radiance.SetSize( linesofsight.NumRays(), numwavelen);
	for (size_t iw = 0; iw < numwavelen ; iw++)
	{
		wavelen = 274 + 0.39*iw;
		ok = ok &&  sasktranengine.CalculateRadiance ( &losradiance, wavelen,  1, &opticalstate );
		for ( size_t ip = 0; ip < linesofsight.NumRays(); ip++)
		{
			radiance.At(ip,iw) = losradiance[ip];
		}
	}
	if (!ok)																	// IF we failed to do the calculation properly
	{																			// then log a messaage
		nxLog::Record(NXLOG_WARNING, "TestSingleScatter, There was an error calculating radiances in the sasktran engine");
	}		// otherwise

	const char* filename = "C:/ARGsoftware/Repos_SasktranV21_MonteCarlo/output/sasktran_singlescatter.txt";
	radiance.WriteToTextFile( filename, true, "%e" );
}



/*-----------------------------------------------------------------------------
 *					TestLatLonUnitSphere		2011-2-3*/
/** **/
/*---------------------------------------------------------------------------*/

void TestLatLonUnitSphere()
{
	SKTRAN_UnitSphereLatLonGrid*	latlongrid;
	bool							ok;
	size_t							i;
	double							sum;
	double							w;

	double zenith  [19] = {5, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120, 130, 140, 150, 160, 170, 180};
	double azimuth [12] = {0, 30, 60, 90, 120, 150, 180, 210, 240, 270, 300, 330};

	SKTRAN_GridDefDiffuseIncomingZenith_V21		zengrid;
	SKTRAN_GridDefDiffuseIncomingAzimuth_V21		azigrid;

	zengrid.SetStatic();
	azigrid.SetStatic();

	zengrid.CopyGridArray( zenith,  N_ELEMENTS(zenith)) ;
	azigrid.CopyGridArray( azimuth, N_ELEMENTS(azimuth)) ;

	latlongrid = new SKTRAN_UnitSphereLatLonGrid;
	latlongrid->AddRef();

	ok = latlongrid->DefineGrid( &zengrid, &azigrid );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"TestLatLonUnitSphere, There was an error defining the Lat/Lon grid"); 
	}
	else
	{
		sum = 0.0;
		for (i = 0; i < latlongrid->NumUnitVectors(); i++)
		{
			w = latlongrid->CubatureWeightAt( i );
			if (w <= 0.0)
			{
				nxLog::Record(NXLOG_WARNING, "TestLatLonUnitSphere(), Vertex %u has negative or zero cubature weight (%e). Thats not very good", (unsigned int)i, (double)w);
			}
			sum += w;;
		}
		ok 	= (sum > 12.56637) && (sum < 12.56638);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"TestLatLonUnitSphere(), The cubature weights on the lat/lon grid do not add up to 4 pi. Thats not good");
		}
	}

	size_t	izen;
	size_t	iazi;
	double	zen;
	double	azi;
	nxVector	unit;
	size_t		vertexindex[3];
	double		weights[3];
	double		newweights[3];
	double		maxdiff;
	double		thediff;
	double		diff;

	maxdiff = 0.0;
	for (izen = 0; izen <= 180*16; izen++)
	{
		zen = izen/16.0;
		for (iazi = 0; iazi <= 360*16; iazi++)
		{
			azi = iazi/16.0;
			unit.FromLatLong( 90.0-zen, azi );
			ok = latlongrid->Triangulate( unit, vertexindex, weights, 3);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "TestLatLonUnitSphere, Error looking up interpolation for ($f,%f)", (double)zen, (double)azi );
			}
			else
			{
	//			if ((zen > 6) && (zen < 168))			// Dont test the fans around the poles
				{
					ok = latlongrid->InterpolateTriangle( unit, vertexindex, newweights);
					if (ok)
					{
						thediff = 0.0;
						for (i = 0; i < 3; i++)
						{
							ok = ok && (newweights[i] < 1.2) && (newweights[i] > -0.2);
							diff = fabs( newweights[i] - weights[i]);
							if (diff > thediff) thediff = diff;
						}
						ok = (thediff < 0.2);
						if (thediff > maxdiff) maxdiff = thediff;
					}
					if (!ok)
					{
						nxLog::Record(NXLOG_WARNING, "TestLatLonUnitSphere(), The vertices returned by Triangulate did not work with INterpolateTriangle");
					}
				}
			}
		}
	}
	printf("The maximum difference between interpolation schemes is %f\n", (double)maxdiff);



	latlongrid->Release();
}
