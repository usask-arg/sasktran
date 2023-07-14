/**
 * SASKTRAN TIR Internal Ray Tracer Specs
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_Specs_Internal_RayTracer::SKTRAN_TIR_Specs_Internal_RayTracer
 * 2019-04-22
 */
SKTRAN_TIR_Specs_Internal_RayTracer::SKTRAN_TIR_Specs_Internal_RayTracer()
{
	m_linesofsighttype = RayTracerTypeTIR::straight;
	m_shellspacing = 0.0;
	m_manualshells = { 0.0 };
	m_setmanualshells = false;
	m_usecurve = false;
	m_parentspecs = nullptr;
	m_groundshiftalt = 0.0;
	m_geoidmodel = nxGeodetic::WGS84;
	m_setmanualgeoidmodel = false;
}

/**
 * SKTRAN_TIR_Specs_Internal_RayTracer::~SKTRAN_TIR_Specs_Internal_RayTracer
 * 2019-04-22
 */
SKTRAN_TIR_Specs_Internal_RayTracer::~SKTRAN_TIR_Specs_Internal_RayTracer()
{
	ReleaseResources();
}


/**
 * SKTRAN_TIR_Specs_Internal_RayTracer::CreateRayFactory
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_RayTracer::CreateRayFactory(
	std::shared_ptr<SKTRAN_RayFactory_Base>& userrayfactory,
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
	RayTracerTypeTIR type,
	bool islos)
{
	bool ok = true;
	std::shared_ptr<SKTRAN_RayFactory_Base>		rayfactory;

	if (RayTracerTypeTIR::curved == type)
	{
		m_usecurve = true;
		ok = ok && CreateCurvedRayFactory(rayfactory, coords);
	}
	else if (RayTracerTypeTIR::straight == type)
	{
		m_usecurve = false;
		ok = ok && CreateGenericShellRayFactory(rayfactory, coords, islos);
	}
	else
	{
		ok = false;
		nxLog::Record(NXLOG_WARNING, "Error Shouldnt be here in SKTRAN_TIR_Specs_Internal_RayTracer::CreateRayTracer");
	}
	userrayfactory = rayfactory;
	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_RayTracer::CreateCurvedRayFactory
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_RayTracer::CreateCurvedRayFactory(
	std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer,
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords)
{
	bool ok = true;

	std::unique_ptr<SKTRAN_RayFactory<SKTRAN_RayOptical_Curved,
		                              SKTRAN_RayTracer_Curved_Shells,
		                              SKTRAN_RayStorage_CurvedPiecewise_TIR>>
		raytracer_local(new SKTRAN_RayFactory<SKTRAN_RayOptical_Curved,
						                      SKTRAN_RayTracer_Curved_Shells,
						                      SKTRAN_RayStorage_CurvedPiecewise_TIR>(coords));

	std::unique_ptr< skRTRefractiveIndex_Profile_TIR> n_profile(new skRTRefractiveIndex_Profile_TIR);

	shellsptr linesofsightshells;
	if (UseManualShells())
	{
		ok = ok && CreateManualShells(linesofsightshells, m_manualshells);
	}
	else
	{
		ok = ok && CreateEvenlySpacedShells(linesofsightshells, m_shellspacing);
	}

	raytracer_local->RayTracer()->Initialize(linesofsightshells, std::move(n_profile));

	raytracer = std::move(raytracer_local);
	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_RayTracer::CreateGenericShellRayFactory
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_RayTracer::CreateGenericShellRayFactory(
	std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer,
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
	bool islos)
{
	bool ok = true;

	std::unique_ptr<SKTRAN_RayFactory<SKTRAN_RayOptical_Straight,
		                              SKTRAN_RayTracer_Straight_Generic,
		                              SKTRAN_RayStorage_Straight_TIR>>
		raytracer_local(new SKTRAN_RayFactory<SKTRAN_RayOptical_Straight,
						                      SKTRAN_RayTracer_Straight_Generic,
						                      SKTRAN_RayStorage_Straight_TIR>(coords));

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
	for (int i = 1; i < (int)raytracingshells->NumGridPoints(); i++)
	{
		std::unique_ptr<SKTRAN_GeometryObject_Sphere> tobj(new SKTRAN_GeometryObject_Sphere(raytracingshells->At(i) + groundradius));
		raytracer_local->RayTracer()->AddGeometryObject(std::move(tobj));
	}
	raytracer_local->RayTracer()->SetEarthRadius(groundradius + m_groundshiftalt);
	raytracer_local->RayTracer()->SetUpperAtmoRadius(groundradius + TOAHeight());

	if (m_parentspecs->CalcWf())
	{
		m_parentspecs->OpticalPropertiesSpecs().AddOpticalInformationToRayTracer(*(raytracer_local->RayTracer()), coords.get());
		if (islos)
		{
			m_parentspecs->WeightingFunctionSpecs().AddWfGeometryToRayTracer(*(raytracer_local->RayTracer()), coords);
		}
	}
	else
	{
		m_parentspecs->OpticalPropertiesSpecs().AddOpticalInformationToRayTracer(*(raytracer_local->RayTracer()), coords.get());
	}

	/* add in a cylinders for the terminator */
	double cylinderspacing = 100;
	int    numcylinders = 10;
	for (int i = 0; i < numcylinders; i++)
	{
		std::unique_ptr<SKTRAN_GeometryObject_Cylinder> tobj(new SKTRAN_GeometryObject_Cylinder(nxVector(0, 0, 1), m_groundshiftalt + groundradius + (i - (numcylinders - 1) / 2)*cylinderspacing));
		raytracer_local->RayTracer()->AddGeometryObject(std::move(tobj));
	}

	raytracer = std::move(raytracer_local);
	return true;
}



/**
 * SKTRAN_TIR_Specs_Internal_RayTracer::CreateLineOfSightRayFactory
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_RayTracer::CreateLineOfSightRayFactory(
	std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer,
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords)
{
	bool ok = true;

	ok = ok && CreateRayFactory(raytracer, coords, m_linesofsighttype, true);

	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_RayTracer::CreateShellRayFactory
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_RayTracer::CreateShellRayFactory(
	std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer,
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords)
{
	bool ok = true;

	std::unique_ptr<SKTRAN_RayFactory<SKTRAN_RayOptical_Straight,
		                              SKTRAN_RayTracer_Shells,
		                              SKTRAN_RayStorage_Straight_TIR>>
		raytracer_local(new SKTRAN_RayFactory<SKTRAN_RayOptical_Straight,
						                      SKTRAN_RayTracer_Shells,
						                      SKTRAN_RayStorage_Straight_TIR>(coords));

	shellsptr raytracingshells;
	if (UseManualShells())
	{
		ok = ok && CreateManualShells(raytracingshells, m_manualshells);
	}
	else
	{
		ok = ok && CreateEvenlySpacedShells(raytracingshells, m_shellspacing);
	}

	raytracer_local->RayTracer()->Initialize(raytracingshells);
	raytracer = std::move(raytracer_local);
	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_RayTracer::ReleaseResources
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_RayTracer::ReleaseResources()
{
	bool ok = true;

	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_RayTracer::ConfigureDefaults
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_RayTracer::ConfigureDefaults()
{
	bool ok = true;


	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_RayTracer::CreateEvenlySpacedShells
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_RayTracer::CreateEvenlySpacedShells(
	shellsptr& rayshells,
	double shellspacing)
{
	bool ok = true;
	size_t numshells;

	rayshells = shellsptr(new SKTRAN_GridDefRayTracingShells_V21);
	numshells = static_cast<size_t>(ceil(TOAHeight() / shellspacing)) + 1;

	std::vector<double> shellheights;
	shellheights.resize(numshells);

	for (size_t idx = 0; idx < numshells; idx++)
	{
		shellheights[idx] = idx * shellspacing;
	}

	ok = ok && rayshells->ConfigureHeights(&shellheights[0], numshells);
	rayshells->SetStatic();
	rayshells->AddRef();
	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_RayTracer::CreateManualShells
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_RayTracer::CreateManualShells(
	shellsptr& rayshells,
	std::vector<double> customshells)
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

/**
 * SKTRAN_TIR_Specs_Internal_RayTracer::Configure
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_RayTracer::Configure(
	const SKTRAN_TIR_Specs_User& specs,
	const SKTRAN_TIR_Specs_Internal_Core* parentspecs)
{
	bool ok = true;

	const SKTRAN_TIR_Specs_User_RayTracer& rayspecs = specs.RayTracingSpecsConst();
	rayspecs.CheckShellParameters();

	m_linesofsighttype = rayspecs.m_linesofsighttype;

	m_shellspacing = rayspecs.m_shellspacing;
	m_setshellspacing = rayspecs.m_setshellspacing;
	m_manualshells = rayspecs.m_manualshells;
	m_setmanualshells = rayspecs.m_setmanualshells;
	m_usecurve = rayspecs.m_usecurve;
	m_groundshiftalt = rayspecs.m_groundshiftalt;
	m_parentspecs = parentspecs;

	m_toaHeight = specs.RayTracingSpecsConst().m_toaHeight;

	m_geoidmodel = rayspecs.m_geoidmodel;
	m_setmanualgeoidmodel = rayspecs.m_setmanualgeoidmodel;

	return ok;
}
