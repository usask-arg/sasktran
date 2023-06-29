#include "include/sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_RayTracer::SKTRAN_HR_Specs_Internal_RayTracer		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Specs_Internal_RayTracer::SKTRAN_HR_Specs_Internal_RayTracer()
{
	m_linesofsighttype		= SKTRAN_HR_RayTracer_Type::SKTRAN_HR_RayTracer_Straight_Generic;
	m_diffusetype			= SKTRAN_HR_RayTracer_Type::SKTRAN_HR_RayTracer_Straight_Generic;
	m_solartype				= SKTRAN_HR_RayTracer_Type::SKTRAN_HR_RayTracer_Straight_Generic;
	m_shellspacing			= 0.0;
	m_manualshells			= { 0.0 };
	m_setmanualshells		= false;
	m_solarshellspacing		= 0.0;
	m_manualsolarshells		= { 0.0 };
	m_setmanualsolarshells	= false;
	m_curvedseparation		= 0.0;
	m_usecurve				= false;
	m_parentspecs			= nullptr;
	m_groundshiftalt		= 0.0;

}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_RayTracer::~SKTRAN_HR_Specs_Internal_RayTracer		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Specs_Internal_RayTracer::~SKTRAN_HR_Specs_Internal_RayTracer()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_RayTracer::CreateRayTracer		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_RayTracer::CreateRayFactory( std::shared_ptr<SKTRAN_RayFactory_Base>&			userrayfactory,
														   std::shared_ptr< const SKTRAN_CoordinateTransform_V2>&	coords,
														   SKTRAN_HR_RayTracer_Type type, bool islos )
{
	bool ok = true;
	std::shared_ptr<SKTRAN_RayFactory_Base>		rayfactory;

	if( SKTRAN_HR_RayTracer_Shells == type )
	{
		ok = ok && CreateShellRayFactory( rayfactory, coords );
	}
	else if( SKTRAN_HR_RayTracer_Curved == type )
	{
		m_usecurve = true;
		ok = ok && CreateCurvedRayFactory( rayfactory, coords );
	}
	else if(SKTRAN_HR_RayTracer_Curved_NoCurve == type )
	{
		m_usecurve = false;
		ok = ok && CreateCurvedRayFactory( rayfactory, coords );
	}
	else if( SKTRAN_HR_RayTracer_Straight_Generic == type )
	{
		ok = ok && CreateGenericShellRayFactory( rayfactory, coords, islos );
	}
	else
	{
		ok = false;
		nxLog::Record(NXLOG_WARNING, "Error Shouldnt be here in SKTRAN_HR_Specs_Internal_RayTracer::CreateRayTracer");
	}
//	rayfactory->SetThisFactory( rayfactory);
	userrayfactory = rayfactory;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_RayTracer::CreateCurvedRayTracer		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_RayTracer::CreateCurvedRayFactory( std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords )
{
	bool ok = true;

	std::unique_ptr< SKTRAN_RayFactory< SKTRAN_RayOptical_Curved,
										SKTRAN_RayTracer_Curved_Shells,
										SKTRAN_RayStorage_CurvedPiecewise_HR
									  > >										raytracer_local ( new SKTRAN_RayFactory< SKTRAN_RayOptical_Curved,
																														 SKTRAN_RayTracer_Curved_Shells,
																														 SKTRAN_RayStorage_CurvedPiecewise_HR>( coords)
																								);

	std::unique_ptr< skRTRefractiveIndex_Profile>								n_profile( new skRTRefractiveIndex_Profile );

	shellsptr linesofsightshells;
	if (UseManualShells())
	{
		ok = ok && CreateManualShells(linesofsightshells, m_manualshells);
	}
	else
	{
		ok = ok && CreateEvenlySpacedShells(linesofsightshells, m_shellspacing);
	}
	

	raytracer_local->RayTracer()->Initialize( linesofsightshells, std::move(n_profile));

	raytracer = std::move( raytracer_local );
//	raytracer->SetThisFactory( raytracer );
//	raytracer->AddRef();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_RayTracer::CreateGenericShellRayTracer		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_RayTracer::CreateGenericShellRayFactory( std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, bool islos )
{
	bool ok = true;

	std::unique_ptr< SKTRAN_RayFactory< SKTRAN_RayOptical_Straight,
										SKTRAN_RayTracer_Straight_Generic,
										SKTRAN_RayStorage_Straight_HR
									  > >										raytracer_local ( new SKTRAN_RayFactory< SKTRAN_RayOptical_Straight,
																														 SKTRAN_RayTracer_Straight_Generic,
																														 SKTRAN_RayStorage_Straight_HR> (coords)
																								);




	shellsptr raytracingshells;
	if (UseManualShells())
	{
		ok = ok && CreateManualShells(raytracingshells, m_manualshells);
	}
	else
	{
		ok = ok && CreateEvenlySpacedShells(raytracingshells, m_shellspacing);
	}	

	double groundradius = coords->AltitudeToRadius(0.0);
	for( int i = 1; i < (int)raytracingshells->NumGridPoints(); i++ )
	{
		std::unique_ptr<SKTRAN_GeometryObject_Sphere> tobj ( new SKTRAN_GeometryObject_Sphere( raytracingshells->At( i ) + groundradius ) );
		raytracer_local->RayTracer()->AddGeometryObject( std::move( tobj ) );
	}
	raytracer_local->RayTracer()->SetEarthRadius( groundradius     + m_groundshiftalt);
	raytracer_local->RayTracer()->SetUpperAtmoRadius( groundradius + TOAHeight() );

	if( m_parentspecs->CalcWf() )
	{
		m_parentspecs->OpticalPropertiesSpecs().AddOpticalInformationToRayTracer( *(raytracer_local->RayTracer()), coords.get() );
		if( islos )
		{
			m_parentspecs->WeightingFunctionSpecs().AddWfGeometryToRayTracer( *(raytracer_local->RayTracer()), coords );
		}
	}
	else
	{
		m_parentspecs->OpticalPropertiesSpecs().AddOpticalInformationToRayTracer( *(raytracer_local->RayTracer()), coords.get() );
	}

	/* add in a cylinders for the terminator */
	double cylinderspacing = 100;
	int    numcylinders = 10;
	for( int i = 0; i < numcylinders; i++ )
	{
		std::unique_ptr<SKTRAN_GeometryObject_Cylinder> tobj ( new SKTRAN_GeometryObject_Cylinder( nxVector(0,0,1), m_groundshiftalt + groundradius + (i-(numcylinders-1)/2)*cylinderspacing ) );
		raytracer_local->RayTracer()->AddGeometryObject( std::move( tobj ) );
	}

	raytracer = std::move( raytracer_local );
//	raytracer->SetThisFactory( raytracer );
	return true;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_RayTracer::CreateLineOfSightRayFactory		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_RayTracer::CreateLineOfSightRayFactory( std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords )
{
	bool ok = true;

	ok = ok && CreateRayFactory( raytracer, coords, m_linesofsighttype, true );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_RayTracer::CreateDiffuseRayFactory		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_RayTracer::CreateDiffuseRayFactory(  std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords)
{
	bool ok = true;

	ok = ok && CreateRayFactory( raytracer, coords, m_diffusetype, false );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_RayTracer::CreateSolarRayTracer		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_RayTracer::CreateSolarRayFactory(  std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords )
{
	bool ok = true;
	
	double oldspacing = m_shellspacing;
	bool oldsetspacing = m_setshellspacing;
	std::vector<double> oldcustomshells = m_manualshells;
	bool oldsetcustomshells = m_setmanualshells;
	double oldgroundshiftalt = m_groundshiftalt;

	m_shellspacing = m_solarshellspacing;
	m_setshellspacing = m_setsolarshellspacing;
	m_manualshells = m_manualsolarshells;
	m_setmanualshells = m_setmanualsolarshells;
	m_groundshiftalt = m_groundshiftalt - 0.1;
	ok = ok && CreateRayFactory( raytracer, coords, m_solartype, false );

	m_shellspacing = oldspacing;
	m_setshellspacing = oldsetspacing;
	m_manualshells = oldcustomshells;
	m_setmanualshells = oldsetcustomshells;
	m_groundshiftalt = oldgroundshiftalt;

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_RayTracer::CreateShellRayTracer		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_RayTracer::CreateShellRayFactory(  std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords  )
{
	bool ok = true;

	std::unique_ptr< SKTRAN_RayFactory< SKTRAN_RayOptical_Straight,
										SKTRAN_RayTracer_Shells,
										SKTRAN_RayStorage_Straight_HR
									  > >								   raytracer_local ( new SKTRAN_RayFactory< SKTRAN_RayOptical_Straight,
																														SKTRAN_RayTracer_Shells,
																														SKTRAN_RayStorage_Straight_HR> (coords)
																							);


	shellsptr raytracingshells;
	if (UseManualShells())
	{
		ok = ok && CreateManualShells(raytracingshells, m_manualshells);
	}
	else
	{
		ok = ok && CreateEvenlySpacedShells(raytracingshells, m_shellspacing);
	}

	raytracer_local->RayTracer()->Initialize( raytracingshells );
	raytracer = std::move( raytracer_local );
//	raytracer->SetThisFactory(raytracer);
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_RayTracer::ReleaseResources		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_RayTracer::ReleaseResources()
{
	bool ok = true;
	
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_RayTracer::ConfigureDefaults		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_RayTracer::ConfigureDefaults()
{
	bool ok = true;


	return ok;
}

bool SKTRAN_HR_Specs_Internal_RayTracer::CreateEvenlySpacedShells( shellsptr& rayshells, double shellspacing )
{
	bool ok = true;
	size_t numshells;

	rayshells = shellsptr (new SKTRAN_GridDefRayTracingShells_V21);
	numshells = static_cast<size_t>(ceil(TOAHeight()/shellspacing)) + 1;

	std::vector<double> shellheights;
	shellheights.resize( numshells );

	for( size_t idx = 0; idx < numshells; idx++ )
	{
		shellheights[idx] = idx*shellspacing;
	}

	ok = ok && rayshells->ConfigureHeights( &shellheights[0], numshells );
	rayshells->SetStatic();
	rayshells->AddRef();
	return ok;
}
/*-----------------------------------------------------------------------------
*					SKTRAN_HR_Specs_Internal_RayTracer::CreateCustomShells		 2018-05-18 */
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_HR_Specs_Internal_RayTracer::CreateManualShells(shellsptr& rayshells, std::vector<double> customshells)
{
	bool ok = true;
	size_t numshells;

	rayshells = shellsptr(new SKTRAN_GridDefRayTracingShells_V21);
	numshells = static_cast<size_t>(customshells.size());

	ok = ok && rayshells->ConfigureHeights(&customshells[0], numshells);
	rayshells->SetStatic();
	rayshells->AddRef();
	return ok;
}

bool SKTRAN_HR_Specs_Internal_RayTracer::Configure(const SKTRAN_HR_Specs_User& specs, const SKTRAN_HR_Specs_Internal_Core* parentspecs)
{
	bool ok = true;

	const SKTRAN_HR_Specs_User_RayTracer& rayspecs = specs.RayTracingSpecsConst();
	rayspecs.CheckShellParameters();

	m_linesofsighttype		= rayspecs.m_linesofsighttype;
	m_diffusetype			= rayspecs.m_diffusetype;
	m_solartype				= rayspecs.m_solartype;

	m_shellspacing			= rayspecs.m_shellspacing;
	m_setshellspacing		= rayspecs.m_setshellspacing;
	m_manualshells			= rayspecs.m_manualshells;
	m_setmanualshells		= rayspecs.m_setmanualshells;
	m_solarshellspacing		= rayspecs.m_solarshellspacing;
	m_setsolarshellspacing	= rayspecs.m_setsolarshellspacing;
	m_manualsolarshells		= rayspecs.m_manualsolarshells;
	m_setmanualsolarshells	= rayspecs.m_setmanualsolarshells;
	m_curvedseparation		= rayspecs.m_curvedseparation;
	m_usecurve				= rayspecs.m_usecurve;
	m_groundshiftalt		= rayspecs.m_groundshiftalt;
	m_parentspecs			= parentspecs;

    m_toaHeight				= specs.RayTracingSpecsConst().m_toaHeight;
        
    return ok;
}

std::vector<double> SKTRAN_HR_Specs_Internal_RayTracer::SolarShellHeights() const
{
	if (UseManualSolarShells())
	{
		return m_manualsolarshells;
	}
	else
	{
		bool ok = true;
		size_t numshells;

		double surface = SurfaceHeight();
		numshells = static_cast<size_t>(ceil((TOAHeight() -surface) / m_solarshellspacing)) + 1;

		std::vector<double> shellheights;
		shellheights.resize(numshells);

		for (size_t idx = 0; idx < numshells; idx++)
		{
			shellheights[idx] = idx*m_solarshellspacing + surface;
		}
		return shellheights;
	}
}