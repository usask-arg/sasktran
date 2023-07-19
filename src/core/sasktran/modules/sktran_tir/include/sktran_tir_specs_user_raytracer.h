/**
 * SASKTRAN TIR User Ray Tracer Specifications
 */

class SKTRAN_TIR_Specs_Internal_Core;

/**
 * SKTRAN_TIR_Specs_User_RayTracer
 *
 * Allows the user to configure properties of the TIR ray tracer.
 * Available settings are:
 *
 *   LinesOfSightType - Sets whether lines of sight are straight or curved.
 *   ShellSpacing - Sets the vertical spacing between the spherical shells
 *     which are used to segment the rays.
 *   ManualShells - A vector of heights at which to create spherical shells.
 *     Typically left blank, but if specified allows the ray to be segmented
 *     at custom heights.
 *   CurvedSeparation
 *   GroundShiftAlt - Allows the ground to be set at a height other than 0,
 *     specified in meters.
 *   TOAHeight - Sets the height of the top of atmosphere in meters.
 *   GeoidModel - Sets the reference ellipsoid used. Default is IAU1976
 */
class SKTRAN_TIR_Specs_User_RayTracer
{
private:
	RayTracerTypeTIR					m_linesofsighttype;
	double								m_shellspacing;
	bool								m_setshellspacing;
	std::vector<double>					m_manualshells;
	bool								m_setmanualshells;
	bool								m_usecurve;
	double								m_groundshiftalt;
	double								m_toaHeight;
	nxGeodetic::GEOID_MODEL				m_geoidmodel;
	bool								m_setmanualgeoidmodel;

private:
	bool ConfigureDefaults();
	bool CheckShellParameters() const;

public:
	SKTRAN_TIR_Specs_User_RayTracer();
	~SKTRAN_TIR_Specs_User_RayTracer();

	bool SetLinesOfSightType(RayTracerTypeTIR type) { m_linesofsighttype = type; return true; }
	bool SetShellSpacing(double spacing) { m_shellspacing = spacing; m_setshellspacing = true; return true; }
	bool SetManualShells(std::vector<double> customshells);
	bool SetGroundShiftAlt(double alt) { m_groundshiftalt = alt; return true; }
	bool SetTOAHeight(double height) { m_toaHeight = height; return true; }
	bool SetGeoidModel(nxGeodetic::GEOID_MODEL model) { m_setmanualgeoidmodel = true; m_geoidmodel = model; return true; }
	const bool UseManualGeoidModel()  const { return m_setmanualgeoidmodel; }
	double GetTOAHeight() const { return m_toaHeight; }
	double GetGroundShiftAlt() const { return m_groundshiftalt; }
	double GetShellSpacing() const { return m_shellspacing; }
	const bool GetUseCurve() const { return m_usecurve; }
	const nxGeodetic::GEOID_MODEL GetGeoidModel() const { return m_geoidmodel; }
	friend class SKTRAN_TIR_Specs_Internal_RayTracer;
	friend class SKTRAN_TIR_Specs_Internal_Core;
};
