/**
 * SASKTRAN TIR Internal Core Specs
 */

#include "include/sktran_tir_internals.h"


/**
 * SKTRAN_TIR_Specs_Internal_Core::SKTRAN_TIR_Specs_Internal_Core
 * 2019-04-22
 */
SKTRAN_TIR_Specs_Internal_Core::SKTRAN_TIR_Specs_Internal_Core()
{
	m_linesofsight = nullptr;
	m_in3dmode = false;
	m_calcwf = false;
}

/**
 * SKTRAN_TIR_Specs_Internal_Core::~SKTRAN_TIR_Specs_Internal_Core
 * 2019-04-22
 */
SKTRAN_TIR_Specs_Internal_Core::~SKTRAN_TIR_Specs_Internal_Core()
{
	ReleaseResources();
}

/**
 * SKTRAN_TIR_Specs_Internal_Core::Configure
 * 2019-04-22
 *
 * Configures the TIR internal specifications from the provided user specifications and lines of sight. Must be called
 * before using the internal specifications object to create optical properties tables or perturbation lists.
 *
 * @param[in] specs The SKTRAN_TIR_Specs_User object containing user-defined settings for the calculation.
 * @param[in] linesofsight Lines of sight to perform the radiative transfer calculation along.
 */
bool SKTRAN_TIR_Specs_Internal_Core::Configure(
	const SKTRAN_Specifications_Base& specs,
	const SKTRAN_LineOfSightArray_V21& linesofsight)
{
	bool ok = true;
	const SKTRAN_TIR_Specs_User userspecs = dynamic_cast<const SKTRAN_TIR_Specs_User&>(specs);
	m_linesofsight = &linesofsight;

	ok = ok && m_raytracerspecs.Configure(userspecs, (this));
	ok = ok && m_integratorspecs.Configure(userspecs);
	ok = ok && m_opttablespecs.Configure(userspecs, &m_raymanager);
	ok = ok && m_wfspecs.Configure(userspecs);

	m_in3dmode = m_opttablespecs.In3dMode();
	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_Core::CreateCoordinates
 * 2019-04-22
 *
 * Creates the coordinate system for use in radiative transfer calculations.
 *
 * @param[out] coords
 * @param[in] linesofsight
 * @param[in] toaHeight
 * @param[in] geoidmodel
 */
bool SKTRAN_TIR_Specs_Internal_Core::CreateCoordinates(
	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>* coords,
	const SKTRAN_LineOfSightArray_V21& linesofsight,
	double toaHeight,
	nxGeodetic::GEOID_MODEL geoidmodel,
	bool userdefinedgeoidmodel)
{
	bool ok = true;

	ok = ok && m_raymanager.SetUpperBoundAltitude(toaHeight);
	ok = ok && m_raymanager.UpdateUndefinedParametersFromLinesOfSight(linesofsight);
	ok = ok && m_raymanager.MakeCoordinateSystem(coords, 0.0, toaHeight, geoidmodel, userdefinedgeoidmodel); // m_rayTracerSpecs is not yet configured 

	m_coords = *coords;

	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_TIR_Specs_Internal_Core::CreateCoordinates, Error creating coordinates.");
	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_Core::ReleaseResources
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_Core::ReleaseResources()
{
	return true;
}

/**
 * SKTRAN_TIR_Specs_Internal_Core::RotateStartAndEnd
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_Core::RotateStartAndEnd(
	const HELIODETIC_UNITVECTOR& look,
	const HELIODETIC_UNITVECTOR& in,
	HELIODETIC_UNITVECTOR& out,
	double theta)
{
	bool ok = true;

	double costheta = nxmath::cosd(theta);
	double sintheta = nxmath::sind(theta);

	double ux = look.X();
	double uy = look.Y();
	double uz = look.Z();

	out.SetCoords((costheta + ux * ux*(1 - costheta)) * in.X() + (ux*uy*(1 - costheta) - uz * sintheta)*in.Y() + (ux*uz*(1 - costheta) + uy * sintheta)*in.Z(),
		(uy*ux*(1 - costheta) + uz * sintheta) * in.X() + (costheta + uy * uy*(1 - costheta))*in.Y() + (uy*uz*(1 - costheta) - ux * sintheta)*in.Z(),
				  (uz*ux*(1 - costheta) - uy * sintheta) * in.X() + (uz*uy*(1 - costheta) + ux * sintheta)*in.Y() + (costheta + uz * uz*(1 - costheta))*in.Z());

	return ok;
}

/**
 * SKTRAN_TIR_Specs_Internal_Core::CalcReferencePoints
 * 2019-04-22
 */
bool SKTRAN_TIR_Specs_Internal_Core::CalcReferencePoints(
	nxVector& in,
	nxVector& out)
{
	bool ok = true;

	ok = ok && m_raymanager.GetBoundingReferences(in, out);

	return ok;
}
