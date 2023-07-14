/**
 * SASKTRAN TIR Internal Optical Properties Specifications
 */

#include "include/sktran_tir_internals.h"

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable
 * 2019-04-22
 */
SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable()
{
	m_heightspacing = 0.0;
	m_numprofiles = 0;
	m_forcecacheupdates = false;
	m_raymanager = nullptr;
}

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::~SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable
 * 2019-04-22
 */
SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::~SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable()
{

}

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::Configure
 * 2019-04-22
 *
 * Copies user settings to this internal specifications object.
 *
 * @param[in] specs
 * @param[in] raymanager
 */
bool SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::Configure(
	const SKTRAN_TIR_Specs_User& specs,
	const SKTRAN_TIR_RayTracingRegionManager* raymanager)
{
	bool ok = true;

	const SKTRAN_TIR_Specs_User_OpticalPropertiesTable& optspecs = specs.OpticalPropertiesSpecsConst();

	m_numprofiles = optspecs.m_numprofiles;
	m_tabledim = optspecs.m_tabledim;

	m_heightspacing = optspecs.m_heightres;
	m_raymanager = raymanager;
	m_forcecacheupdates = optspecs.m_forcecacheupdates;


	if (optspecs.m_normal.IsValid() && optspecs.m_reference.IsValid())
	{
		// user specified normal and reference
		m_normal = optspecs.m_normal;
		m_reference = optspecs.m_reference;
	}
	else
	{
		// Calculate normal and reference from the LOS
		MakeNormalAndReferenceFromLOS();
	}

	if (optspecs.m_anglegrid.size() != 0)
	{
		// user specified angle grid
		m_anglegrid = optspecs.m_anglegrid;
	}
	else
	{
		// create default angle grid
		MakeDefaultAngleGrid();
	}
	m_heightgrid = optspecs.m_heightgrid;
	//ok = ok && ConfigureDefaults();

	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::MakeNormalAndReferenceFromLOS
 * 2019-04-22
 */
void SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::MakeNormalAndReferenceFromLOS()
{
	double lat, lon;
	m_raymanager->GetReferencePoint(&lat, &lon);

	nxGeodetic geo = nxGeodetic(lat, lon);
	m_reference = geo.Location();
	m_reference = m_reference.UnitVector();

	nxVector in, out;
	m_raymanager->GetBoundingReferences(in, out);

	m_normal = in.Cross(m_reference);
	m_normal = -1.0*m_normal.UnitVector();

	// if 'in' is aligned with 'm_reference' the cross product will be [0, 0, 0]
	// we do not want this, so use 'out' instead
	if (m_normal.IsZero())
	{
		m_normal = out.Cross(m_reference);
		m_normal = m_normal.UnitVector();
	}
}

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::AddOpticalInformationToRayTracer
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::AddOpticalInformationToRayTracer(
	SKTRAN_RayTracer_Straight_Generic& raytracer,
	const SKTRAN_CoordinateTransform_V2* coords) const
{
	if (m_tabledim == AtmosphereDimensionTIR::dim1)
	{
		return true;
	}
	else if (m_tabledim == AtmosphereDimensionTIR::dim2)
	{
		return AddPlaneInformationToRayTracer(raytracer, *coords);
	}
	return false;
}

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::AddPlaneInformationToRayTracer
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::AddPlaneInformationToRayTracer(
	SKTRAN_RayTracer_Straight_Generic& raytracer,
	const SKTRAN_CoordinateTransform_V2& coords) const
{
	nxVector x, y;
	HELIODETIC_VECTOR xh, yh;
	xh = coords.GeographicToHelio(m_reference);
	yh = coords.GeographicToHelio(m_reference.Cross(m_normal));
	x.SetCoords(xh.X(), xh.Y(), xh.Z());
	y.SetCoords(yh.X(), yh.Y(), yh.Z());
	for (size_t idx = 0; idx < m_anglegrid.size(); idx++)
	{
		nxVector normal;
		normal = nxmath::cosd(m_anglegrid[idx] + 90.0) * x + nxmath::sind(m_anglegrid[idx] + 90.0) * y;
		raytracer.AddGeometryObject(std::unique_ptr<SKTRAN_GeometryObject_Plane>(new SKTRAN_GeometryObject_Plane(normal)));
	}
	return true;
}

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::MakeDefaultAngleGrid
 * 2019-04-22
 */
void SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::MakeDefaultAngleGrid()
{
	// default is 15 degrees on either side of the tangent point
	double spacing = 30.0 / double(m_numprofiles - 1);
	m_anglegrid.resize(m_numprofiles);
	for (size_t idx = 0; idx < m_numprofiles; idx++)
	{
		m_anglegrid[idx] = -15.0 + spacing * idx;
	}
}

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::ConfigureDefaults
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::ConfigureDefaults()
{
	bool ok = true;

	m_tabledim = AtmosphereDimensionTIR::dim1;
	m_heightspacing = 500;

	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::Create1dTable
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::Create1dTable(
	OpticalTablePtrTIR& table,
	const SKTRAN_CoordinateTransform_V2& coords,
	double toaHeight)
{
	bool ok = true;

	std::unique_ptr<SKTRAN_TIR_TableOpticalProperties> opttable(new SKTRAN_TIR_TableOpticalProperties);
	SKTRAN_GridDefOpticalPropertiesRadii_V21* heightgrid = new SKTRAN_GridDefOpticalPropertiesRadii_V21;

	nxVector referenceHELIO;
	referenceHELIO.SetCoords(coords.ReferencePoint(0.0).UnitVector().X(), coords.ReferencePoint(0.0).UnitVector().Y(), coords.ReferencePoint(0.0).UnitVector().Z());

	SKTRAN_UnitSphere_V2* unitsphere = new SKTRAN_UnitSphere_Dummy(referenceHELIO);

	ok = ok && MakeHeightGrid(*heightgrid, toaHeight);

	opttable->SetAltitudes(*heightgrid);
	opttable->SetUnitSphere(*unitsphere);
	table = std::move(opttable);
	table->AddRef();

	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::MakeHeightGrid
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::MakeHeightGrid(
	SKTRAN_GridDefOpticalPropertiesRadii_V21& heightgrid,
	double toaHeight)
{
	bool ok = true;

	if (m_heightgrid.size() != 0)
	{
		return heightgrid.ConfigureAltitudes(&m_heightgrid[0], m_heightgrid.size());
	}

	size_t numalts = static_cast<size_t>(ceil((toaHeight) / m_heightspacing)) + 1;

	std::vector<double> altgrid;
	altgrid.resize(numalts);
	for (size_t i = 0; i < numalts; i++)
	{
		altgrid[i] = i * m_heightspacing;
	}
	ok = ok && heightgrid.ConfigureAltitudes(&altgrid[0], numalts);
	heightgrid.SetGridSearchMode(SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM);
	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::MakeLOSPlaneSphere
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::MakeLOSPlaneSphere(
	SKTRAN_UnitSphere_V2** sphere,
	const SKTRAN_CoordinateTransform_V2& coords)
{
	bool ok = true;
	nxVector look;
	std::vector<nxVector> locations;

	SKTRAN_UnitSphere_Plane* spherelocal = new SKTRAN_UnitSphere_Plane;

	nxVector x, y;
	HELIODETIC_VECTOR xh, yh;
	xh = coords.GeographicToHelio(m_reference);
	yh = coords.GeographicToHelio(m_reference.Cross(m_normal));
	x.SetCoords(xh.X(), xh.Y(), xh.Z());
	y.SetCoords(yh.X(), yh.Y(), yh.Z());

	for (const auto& th : m_anglegrid)
	{
		nxVector loc = nxmath::cosd(th) * x + nxmath::sind(th) * y;
		locations.emplace_back(loc);
	}

	spherelocal->ConstructPlane(locations);
	*sphere = spherelocal;

	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::CreateOpticalTable
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::CreateOpticalTable(
	OpticalTablePtrTIR& table,
	const SKTRAN_CoordinateTransform_V2& coords,
	double toaHeight)
{
	bool ok = true;

	if (m_tabledim == AtmosphereDimensionTIR::dim1)
	{
		ok = ok && Create1dTable(table, coords, toaHeight);
	}
	else if (m_tabledim == AtmosphereDimensionTIR::dim2)
	{
		ok = ok && Create3DUnitSphereTable(table, coords, toaHeight);
	}
	else
	{
		ok = false;
		nxLog::Record(NXLOG_WARNING, "Error Should not be here, SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::CreateOpticalTable");
	}
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::CreateOpticalTable, Could not create optical properties table.");

	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::Create3DUnitSphereTable
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::Create3DUnitSphereTable(
	OpticalTablePtrTIR& table,
	const SKTRAN_CoordinateTransform_V2& coords,
	double toaHeight)
{
	bool ok = true;

	std::unique_ptr<SKTRAN_TIR_TableOpticalProperties> opttable(new SKTRAN_TIR_TableOpticalProperties);
	SKTRAN_GridDefOpticalPropertiesRadii_V21* heightgrid = new SKTRAN_GridDefOpticalPropertiesRadii_V21;
	SKTRAN_UnitSphere_V2* unitsphere;

	ok = ok && MakeHeightGrid(*heightgrid, toaHeight);
	if (AtmosphereDimensionTIR::dim2 == m_tabledim)
	{
		ok = ok && MakeLOSPlaneSphere(&unitsphere, coords);
	}
	else
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable::Create3DUnitSphereTable, Invalid atmosphere dimensionality. This function should only be called for 2D atmospheres, 3D atmospheres are not supported by the TIR engine");
		ok = false;
	}

	opttable->SetAltitudes(*heightgrid);
	opttable->SetUnitSphere(*unitsphere);
	opttable->SetForceCacheUpdates(m_forcecacheupdates);

	table = std::move(opttable);
	table->AddRef();

	return ok;
}
