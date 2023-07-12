//#include "sktran_hr_internals.h"

class SKTRAN_HR_Specs_Internal_Core;


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_RayTracer		2014-10-30*/
/** Factory class for creating ray tracers used in the HR model.
 *
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Specs_Internal_RayTracer
{
	private:
		typedef										std::shared_ptr<SKTRAN_GridDefRayTracingShells_V21> shellsptr;
		SKTRAN_HR_RayTracer_Type					m_linesofsighttype;
		SKTRAN_HR_RayTracer_Type					m_diffusetype;
		SKTRAN_HR_RayTracer_Type					m_solartype;
		double										m_shellspacing;
		bool										m_setshellspacing;
		std::vector<double>							m_manualshells;
		bool										m_setmanualshells;
		double										m_solarshellspacing;
		bool										m_setsolarshellspacing;
		std::vector<double>							m_manualsolarshells;
		bool										m_setmanualsolarshells;
		double										m_curvedseparation;
		bool										m_usecurve;
		const SKTRAN_HR_Specs_Internal_Core*		m_parentspecs;
		double										m_groundshiftalt;
        double                                      m_toaHeight;

	private:
		virtual bool				CreateShellRayFactory		( std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords );
		virtual bool				CreateGenericShellRayFactory( std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, bool islos );
		virtual bool				CreateCurvedRayFactory		( std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords );
		virtual bool				CreateEvenlySpacedShells	( shellsptr& rayshells, double shellspacing );
		virtual bool				CreateManualShells			( shellsptr& rayshells, std::vector<double> shellboundaries);
		virtual bool				ConfigureDefaults			();
		virtual bool				ReleaseResources			();
		virtual bool				CreateRayFactory				( std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_HR_RayTracer_Type type, bool islos );
		bool						UseManualShells				() const { return (m_setmanualshells) && (!m_setshellspacing); };
	public:
									SKTRAN_HR_Specs_Internal_RayTracer	();
		virtual					   ~SKTRAN_HR_Specs_Internal_RayTracer	();
		virtual bool				Configure							( const SKTRAN_HR_Specs_User& specs, const SKTRAN_HR_Specs_Internal_Core* parentspecs );
		virtual bool				CreateLineOfSightRayFactory			( std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords);
		virtual bool				CreateDiffuseRayFactory				( std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords);
		virtual bool				CreateSolarRayFactory				( std::shared_ptr<SKTRAN_RayFactory_Base>& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords);

		double						SolarShellSpacing() const { return m_solarshellspacing; }
		std::vector<double>			SolarShellHeights() const;
        double                      TOAHeight		 () const { return m_toaHeight;}
        double                      SurfaceHeight    () const { return m_groundshiftalt;}
		bool						UseManualSolarShells() const { return (m_setmanualsolarshells) && (!m_setsolarshellspacing); };
		bool 						LineOfSightRefractionEnabled() const {return (m_linesofsighttype == SKTRAN_HR_RayTracer_Type::SKTRAN_HR_RayTracer_Curved);}
		bool 						SolarRefractionEnabled() const {return (m_solartype == SKTRAN_HR_RayTracer_Type::SKTRAN_HR_RayTracer_Curved);}
		bool 						DiffuseRefractionEnabled() const {return (m_diffusetype == SKTRAN_HR_RayTracer_Type::SKTRAN_HR_RayTracer_Curved);}


};
