//#include "sktran_hr_internals.h"

class SKTRAN_HR_Specs_User_RayTracer
{
	private:
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
		double										m_groundshiftalt;
        double                                      m_toaHeight;
		bool										m_nadir_referencepoint_onground;

	private:
		bool										ConfigureDefaults();
		bool										CheckShellParameters() const;
	public:
		SKTRAN_HR_Specs_User_RayTracer();
		~SKTRAN_HR_Specs_User_RayTracer();

		bool										SetLinesOfSightType( SKTRAN_HR_RayTracer_Type type ) { m_linesofsighttype = type; return true;}
		bool										SetDiffuseType( SKTRAN_HR_RayTracer_Type type ) { m_diffusetype = type; return true;}
		bool										SetSolarType( SKTRAN_HR_RayTracer_Type type ) { m_solartype = type; return true;}
		bool										SetShellSpacing( double spacing ) { m_shellspacing = spacing; m_setshellspacing = true; return true; }
		bool										SetManualShells( std::vector<double> customshells );
		bool										SetCurvedSeparation( double separation ) { m_curvedseparation = separation; return true;}
		bool										SetUseCurve	( bool usecurve )	{ m_usecurve = usecurve; return true; }
		bool										SetSolarShellSpacing( double spacing ) { m_solarshellspacing = spacing; m_setsolarshellspacing = true; return true; }
		bool										SetManualSolarShells( std::vector<double> customshells );
		bool										SetGroundShiftAlt(double alt) { m_groundshiftalt = alt; return true; }
		double										GetGroundShiftAlt ( ) const { return m_groundshiftalt; }
        bool                                        SetTOAHeight (double height) { m_toaHeight=height; return true; }
        double                                      GetTOAHeight ( ) const { return m_toaHeight; }
		bool										SetNadirReferencePointOnGround ( bool onground ) { m_nadir_referencepoint_onground = onground; return true; }
		bool										GetNadirReferencePointOnGround ( ) const { return m_nadir_referencepoint_onground; }
		friend class SKTRAN_HR_Specs_Internal_RayTracer;
		friend class SKTRAN_HR_Specs_Internal_Core;
		
};
