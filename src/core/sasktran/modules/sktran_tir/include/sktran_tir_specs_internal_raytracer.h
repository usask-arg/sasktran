/**
 * SASKTRAN TIR Internal Ray Tracer Specifications
 * 2019-04-16
 */

/**
 * SKTRAN_TIR_Specs_Internal_RayTracer
 *
 * Factory class for creating ray tracers in the TIR model. See the TIR User
 * RayTracing Specs header file for a detailed description of the available
 * settings. Supports both straight and curved rays.
 */
class SKTRAN_TIR_Specs_Internal_RayTracer
{
private:
	typedef std::shared_ptr<SKTRAN_GridDefRayTracingShells_V21> shellsptr;
	
	RayTracerTypeTIR						m_linesofsighttype;
	double									m_shellspacing;
	bool									m_setshellspacing;
	std::vector<double>						m_manualshells;
	bool									m_setmanualshells;
	bool									m_usecurve;
	const SKTRAN_TIR_Specs_Internal_Core*	m_parentspecs;
	double									m_groundshiftalt;
	double									m_toaHeight;
	nxGeodetic::GEOID_MODEL					m_geoidmodel;
	bool									m_setmanualgeoidmodel;

private:
	virtual bool CreateShellRayFactory(std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer,
									   std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords);
	virtual bool CreateGenericShellRayFactory(std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer,
											  std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
											  bool islos);
	virtual bool CreateCurvedRayFactory(std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer,
										std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords);
	virtual bool CreateEvenlySpacedShells(shellsptr& rayshells,
										  double shellspacing);
	virtual bool CreateManualShells(shellsptr& rayshells,
									std::vector<double> shellboundaries);
	virtual bool ConfigureDefaults();
	virtual bool ReleaseResources();
	virtual bool CreateRayFactory(std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer,
								  std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
								  RayTracerTypeTIR type,
								  bool islos);
	bool UseManualShells() const { return (m_setmanualshells) && (!m_setshellspacing); };

public:
	SKTRAN_TIR_Specs_Internal_RayTracer();
	virtual ~SKTRAN_TIR_Specs_Internal_RayTracer();
	virtual bool Configure(const SKTRAN_TIR_Specs_User& specs,
						   const SKTRAN_TIR_Specs_Internal_Core* parentspecs);
	virtual bool CreateLineOfSightRayFactory(std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer,
											 std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords);

	double TOAHeight() const { return m_toaHeight; }
	double SurfaceHeight() const { return m_groundshiftalt; }

	bool UseCurvedRays() const { return m_usecurve; }

	const nxGeodetic::GEOID_MODEL GetGeoidModel() const { return m_geoidmodel; }
};
