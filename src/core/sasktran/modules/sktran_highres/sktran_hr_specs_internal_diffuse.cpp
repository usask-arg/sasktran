#include "include/sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Diffuse::SKTRAN_HR_Specs_Internal_Diffuse		2013-06-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Specs_Internal_Diffuse::SKTRAN_HR_Specs_Internal_Diffuse()
{
	m_outgoingsphereobj  = nullptr;
	m_numabove           = 0;
	m_numhoriz           = 0;
	m_numbelow           = 0;
	m_numazi             = 0;
	m_numoutgoing        = 0;
	m_uselegacyv21sphere = false;
	m_horizonsize        = 0.0;
	m_incomingspheretype = SKTRAN_HR_DiffuseIncomingSphereType::SKTRAN_HR_DiffuseIncomingType_Default;
    m_diffuseMatrixStorageMethod = SKTRAN_HR_DiffuseMatrixStorageMethod::scalar;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Diffuse::~SKTRAN_HR_Specs_Internal_Diffuse		2013-06-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Specs_Internal_Diffuse::~SKTRAN_HR_Specs_Internal_Diffuse()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Diffuse::Initialize		2013-06-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_Diffuse::Initialize( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords )
{
	m_coords = coords;
	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Diffuse::ConfigureIncomingUnitSphere		2013-06-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_Diffuse::ConfigureIncomingUnitSphere( const SKTRAN_UnitSphere_V2** unitsphere, const double& altitude, bool isground )
{
	bool ok = true;

	SKTRAN_UnitSphereLatLonGrid* unitsphere_temp = new SKTRAN_UnitSphereLatLonGrid;

	SKTRAN_GridDefDiffuseIncomingAzimuth_V21 azigrid;
	SKTRAN_GridDefDiffuseIncomingZenith_V21  zengrid;

	if( m_incomingspheretype == SKTRAN_HR_DiffuseIncomingType_Default )
	{
		ok = ok && MakeDefaultIncomingZenGrid( zengrid, altitude, isground );
	}
	else if ( m_incomingspheretype == SKTRAN_HR_DiffuseIncomingType_HardCode )
	{
		ok = ok && MakeSasktran21IncomingZenGrid( zengrid, altitude, isground );
	}
	else if ( m_incomingspheretype == SKTRAN_HR_DiffuseIncomingType_v21Shifted )
	{
		ok = ok && MakeSasktran21HorizonShiftZenGrid( zengrid, altitude, isground );
	}
	else
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_Internal_Diffuse::ConfigureIncomingUnitSphere Shouldn't be here" );
	}
	ok = ok && MakeDefaultIncomingAziGrid( azigrid, altitude );

	unitsphere_temp->DefineGrid( &zengrid, &azigrid );
	unitsphere_temp->AddRef();

	*unitsphere = unitsphere_temp;

	return ok;
}

bool SKTRAN_HR_Specs_Internal_Diffuse::ConfigureOutgoingUnitSphere()
{
	bool ok = true;

	SKTRAN_UnitSphereME* sphere = new SKTRAN_UnitSphereME( m_numoutgoing );
	if(nullptr!=m_outgoingsphereobj) m_outgoingsphereobj->Release();
	m_outgoingsphereobj = new SKTRAN_HR_OutgoingSphereObject_Base;
	if( nullptr!=m_outgoingsphereobj ){
		m_outgoingsphereobj->AddRef();
		m_outgoingsphereobj->SetOutgoingSphere( sphere );
	} else{
		delete sphere;
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Diffuse::ReleaseResources		2013-06-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_Diffuse::ReleaseResources()
{
	if( nullptr!=m_outgoingsphereobj ) m_outgoingsphereobj->Release();
	m_outgoingsphereobj = nullptr;
	return true;
}


bool SKTRAN_HR_Specs_Internal_Diffuse::MakeSasktran21IncomingZenGrid( SKTRAN_GridDefDiffuseIncomingZenith_V21& zengrid, const double& altitude, const bool& isgroundpoint )
{
	bool ok = true;
	static double g_defaultincomingzenith[] = { 15, 30, 40, 50, 60, 70, 75, 80, 85, 87, 89, 90, 91, 93, 95, 100, 105, 110, 120, 130, 140, 150, 165, 180};
	double lastzen;
	std::vector<double> definedzenith;
	std::vector<double>::const_iterator iter;
	definedzenith.resize(24);
	size_t numpoints;


	lastzen = 0.0;
	for( size_t idx = 0; idx < 24; idx++ )
	{
		definedzenith[idx] = (g_defaultincomingzenith[idx] + lastzen) * 0.5;
		lastzen = g_defaultincomingzenith[idx];
	}

	if( isgroundpoint )
	{
		iter = std::lower_bound( definedzenith.begin(), definedzenith.end(), 90.0+0.00001 );
		numpoints = iter - definedzenith.begin();
	}
	else
	{
		numpoints = 24;
	}

	zengrid.AllocateGridArray( numpoints );

	for( size_t idx = 0; idx < numpoints; idx++ )
	{
		zengrid.AtVar(idx) = definedzenith[idx];
	}
	zengrid.SetIsGroundPoint( isgroundpoint );


	return ok;
}

bool SKTRAN_HR_Specs_Internal_Diffuse::MakeSasktran21HorizonShiftZenGrid( SKTRAN_GridDefDiffuseIncomingZenith_V21& zengrid, const double& altitude, const bool& isgroundpoint )
{
	double		zen_eighty;
	double		atmosDelta;
	double		horizonDelta;
	double		groundDelta;
	double		horizon;
	size_t		idx;
	size_t		numpoints;
	size_t		zenithCtr;
	size_t		horizonCtr;
	size_t		groundCtr;
	bool		ok;
	size_t		numatmos;
	size_t		numhorizon;
	size_t		numground;
//	bool		isgroundpoint;
	std::vector<double>	zenithbuffer;
	double				earthradius = 6360000.0;

	
	numatmos         = m_numabove;				// Number of zenith points in upward regions above 80 degrees zenith angle
	numhorizon       = m_numhoriz;			// Number of zenith points between 80 degrees zenith and the horizon
	numground        = m_numbelow;				// Number of zenith points below the horizon.

//		isgroundpoint    = IsGroundPoint( altitude );							// See is this is a ground point
	earthradius      =  m_coords->AltitudeToRadius( 0.0 );
	zengrid.SetIsGroundPoint(isgroundpoint);							//
	zen_eighty       = 80.0;												// begin horizon zeniths at 80 degrees 

	if (isgroundpoint )														// If this is for a point that hits the ground
	{																		// then only get half the number of points
		numhorizon   = numhorizon/2;										// For some reason halves the number of points for ground points
		numground    = 0;													// There are no points below the ground
		horizon      = 90.0;
		atmosDelta   = zen_eighty/numatmos;									// Atmosphere region between zenith angle 0 and 80 degrees.
		horizonDelta = ( horizon  - zen_eighty )/numhorizon;				// Horizon is SZA between 80 and the horizon.
		groundDelta  = 0.0;
	}
	else
	{
		horizon      = 90.0 + 180.0*acos( (earthradius)/(altitude + earthradius) )/nxmath::Pi;	// Find the horizon.
		atmosDelta   = zen_eighty/numatmos;									// Atmosphere region between zenith angle 0 and 80 degrees.
		horizonDelta = ( horizon	- zen_eighty )/numhorizon;				// Horizon is SZA between 80 and the horizon.
		groundDelta  = ( 180.0		- horizon    )/numground;				// Ground is SZA between horizon and 180 degrees.
	}

	numpoints   = numatmos + numhorizon + numground + 1;
	zenithbuffer.assign(numpoints, 0.0);
	ok          = zengrid.AllocateGridArray( numpoints-1 );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_GridSpecificationsLegacy_V2::ConfigureIncomingZenithGrid, Error allocating memory for %Iu elements", (size_t)numpoints );
	}
	else
	{
		idx = 0;

		for( zenithCtr = 0; zenithCtr < numatmos; ++zenithCtr )
		{
			zenithbuffer[idx++] = zenithCtr * atmosDelta;
		}

		for( horizonCtr = 0; horizonCtr < numhorizon+1; ++horizonCtr )
		{
			zenithbuffer[idx++]  = zen_eighty  + horizonCtr * horizonDelta;
		}

		if (numground > 0)
		{
			for( groundCtr = numground; groundCtr > 0; --groundCtr )
			{
				zenithbuffer[idx++] = 180 - (groundCtr-1)* groundDelta;
			}
		}
		NXASSERT(( idx == numpoints ));
		for (idx=0; idx < numpoints-1; ++idx )
		{
			zengrid.AtVar(idx) = 0.5*(zenithbuffer[idx] + zenithbuffer[idx+1]);
		}
	}
	
	return ok;
}
/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Diffuse::MakeDefaultIncomingZenGrid		2013-06-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_Diffuse::MakeDefaultIncomingZenGrid( SKTRAN_GridDefDiffuseIncomingZenith_V21& zengrid, const double& altitude, const bool& isgroundpoint )
{
	double		atmosDelta;
	double		horizonDelta;
	double		groundDelta;
	double		horizon;
	double		horizonend;
	double		horizonstart;
	double		horizonsize;
	size_t		idx;
	size_t		numpoints;
	size_t		zenithCtr;
	size_t		horizonCtr;
	size_t		groundCtr;
	bool		ok;
	size_t		numatmos;
	size_t		numhorizon;
	size_t		numground;
	std::vector<double>	zenithbuffer;
	double				earthradius;
	{

		numatmos         = m_numabove;				// Number of zenith points in upward regions above 80 degrees zenith angle
		numhorizon       = m_numhoriz;				// Number of zenith points between 80 degrees zenith and the horizon
		numground        = m_numbelow;				// Number of zenith points below the horizon.

		earthradius      = m_coords->AltitudeToRadius( 0.0 );
		zengrid.SetIsGroundPoint(isgroundpoint);							//
		horizonsize		 = m_horizonsize;

		if (isgroundpoint )														// If this is for a point that hits the ground
		{																		// then only get half the number of points
			horizonstart = 90 - horizonsize/2;
			numhorizon   = numhorizon/2;										// For some reason halves the number of points for ground points
			numground    = 0;													// There are no points below the ground
			horizon      = 90.0;
			atmosDelta   = horizonstart/numatmos;									// Atmosphere region between zenith angle 0 and 80 degrees.
			horizonDelta = ( horizon  - horizonstart )/numhorizon;				// Horizon is SZA between 80 and the horizon.
			groundDelta  = 0.0;
		}
		else
		{
			
			horizon      = 90.0 + 180.0*acos( (earthradius)/(altitude + earthradius) )/nxmath::Pi;	// Find the horizon.
			horizonend   = horizon + horizonsize/2;
			horizonstart = horizon - horizonsize/2;
			atmosDelta   = horizonstart/numatmos;									// Atmosphere region between zenith angle 0 and 80 degrees.
			horizonDelta = ( horizonend	- horizonstart )/numhorizon;				// Horizon is SZA between 80 and the horizon.
			groundDelta  = ( 180.0		- horizonend    )/numground;				// Ground is SZA between horizon and 180 degrees.
		}

		numpoints   = numatmos + numhorizon + numground + 1;
		zenithbuffer.assign(numpoints, 0.0);
		ok          = zengrid.AllocateGridArray( numpoints-1 );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_Internal_Diffuse::MakeDefaultZenithGrid, Error allocating memory for %Iu elements", (size_t)numpoints );
		}
		else
		{
			idx = 0;

			for( zenithCtr = 0; zenithCtr < numatmos; ++zenithCtr )
			{
				zenithbuffer[idx++] = zenithCtr * atmosDelta;
			}

			for( horizonCtr = 0; horizonCtr < numhorizon+1; ++horizonCtr )
			{
				zenithbuffer[idx++]  = horizonstart  + horizonCtr * horizonDelta;
			}

			if (numground > 0)
			{
				for( groundCtr = numground; groundCtr > 0; --groundCtr )
				{
					zenithbuffer[idx++] = 180 - (groundCtr-1)* groundDelta;
				}
			}
			NXASSERT(( idx == numpoints ));
			for (idx=0; idx < numpoints-1; ++idx )
			{
				zengrid.AtVar(idx) = 0.5*(zenithbuffer[idx] + zenithbuffer[idx+1]);
			}
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Diffuse::MakeDefaultIncomingAziGrid		2013-06-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_Diffuse::MakeDefaultIncomingAziGrid( SKTRAN_GridDefDiffuseIncomingAzimuth_V21& azigrid, const double& altitude )
{
	bool										ok = true;
	size_t										aziidx;
	azigrid.AllocateGridArray( m_numazi );

	for (aziidx = 0; aziidx < m_numazi; aziidx++)
	{
		azigrid.AtVar(aziidx) = double(aziidx)/double(m_numazi) * 360.0;
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Diffuse::MakeDiffusePoint		2013-06-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_Diffuse::MakeDiffusePoint( SKTRAN_HR_Diffuse_Point& diffusepoint, const HELIODETIC_POINT& location, bool isground )
{
	const SKTRAN_UnitSphere_V2* incomingsphere = nullptr;
	bool						ok = true;
	
	ok = ok && ConfigureIncomingUnitSphere( &incomingsphere, location.Altitude(), isground );
	ok = ok && diffusepoint.ConfigureSpheres( incomingsphere, m_outgoingsphereobj, location, isground);
	if (incomingsphere != nullptr)incomingsphere->Release();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Diffuse::Configure		2013-06-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_Diffuse::Configure( const SKTRAN_HR_Specs_User& specs )
{
	bool ok = true;
	m_numoutgoing = specs.DiffuseSpecsConst().GetNumOutgoing();
	ok = ok && ConfigureDefaults();

	m_numabove = specs.DiffuseSpecsConst().GetNumBeforeHoriz();
	m_numhoriz = specs.DiffuseSpecsConst().GetNumHoriz();
	m_numbelow = specs.DiffuseSpecsConst().GetNumAfterHoriz();
	m_numazi   = specs.DiffuseSpecsConst().GetNumAzi();
	m_horizonsize = specs.DiffuseSpecsConst().m_horizonsize;
	m_uselegacyv21sphere = specs.DiffuseSpecsConst().m_forcev21incomingsphere;
	m_incomingspheretype = specs.DiffuseSpecsConst().m_incomingspheretype;
    m_diffuseMatrixStorageMethod = specs.DiffuseSpecsConst().GetScatteringMatrixStorageMethod();

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Diffuse::ConfigureDefaults		2013-06-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_Diffuse::ConfigureDefaults()
{
	bool ok = true;

	//m_numabove = 6;
	//m_numhoriz = 8;
	//m_numbelow = 10;
	//m_numazi = 12;

	ok = ok && ConfigureOutgoingUnitSphere();

	return ok;
}

bool SKTRAN_HR_Specs_Internal_Diffuse::MakeSecondOrderSource( SKTRAN_HR_Diffuse_Second_Order_Source** source )
{
	bool ok = true;

	const SKTRAN_UnitSphere_V2*			unitsphere;
	SKTRAN_HR_Diffuse_Second_Order_Source*		source_local = new SKTRAN_HR_Diffuse_Second_Order_Source;
	ok = ok && ConfigureIncomingUnitSphere( &unitsphere, 500, false );
	
	source_local->SetUnitSphere( *unitsphere );

	*source = source_local;
	source_local->AddRef();

	return ok;
}
