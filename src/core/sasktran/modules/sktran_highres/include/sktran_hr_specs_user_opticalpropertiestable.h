//#include "sktran_hr_internals.h"

class SKTRAN_HR_Specs_User_OpticalPropertiesTable
{
	private:
		SKTRAN_HR_OpticalPropertiesTableType			m_optproptype;
		size_t											m_numprofiles;
		double											m_heightres;
		double											m_scatres;
		
		size_t                                          m_maxPolarizationOrder;
		double 											m_polarizationHigherOrderFraction;
		SKTRAN_HR_PolHOType                             m_polHigherOrderBehaviour;
		SKTRAN_HR_AtmosphereHasDelta                    m_atmosphereHasDelta;

		double											m_coneanglesep;
		size_t											m_numcones;
		size_t											m_profilepercone;
		
		nxVector										m_normal;			
		nxVector										m_reference;
		std::vector<double>								m_anglegrid;
		std::vector<double>								m_heightgrid;
		std::vector<double>								m_precachewavel;

		bool											m_forcecacheupdates;
	private:
		bool											ConfigureDefaults();
	public:
		SKTRAN_HR_Specs_User_OpticalPropertiesTable() { ConfigureDefaults(); }
		~SKTRAN_HR_Specs_User_OpticalPropertiesTable() {};
		void											SetOpticalPropertiesType( SKTRAN_HR_OpticalPropertiesTableType type ) { m_optproptype = type; }
		void                                            SetMaxPolarizationOrder     ( size_t order ) { m_maxPolarizationOrder = order; }
		void                                            SetPolarizationHigherOrderBehaviour ( SKTRAN_HR_PolHOType type ){ m_polHigherOrderBehaviour = type;}
		void                                            SetAtmosphereHasDelta   ( SKTRAN_HR_AtmosphereHasDelta         val  ) { m_atmosphereHasDelta = val; }
		void											SetNumProfiles( size_t numprofiles ) { m_numprofiles = numprofiles; }
		void											SetHeightRes( double h ) { m_heightres = h; }
		void											SetNormalAndReference( const nxVector& normal, const nxVector& reference ) { m_normal = normal; m_reference = reference; }
		void											SetAngleGrid( const std::vector<double>& anglegrid ) { m_anglegrid = anglegrid; }
		void											SetHeightGrid( const std::vector<double>& heightgrid ) { m_heightgrid = heightgrid; }
		void											SetDelaunayParam( double coneanglesep, size_t numcones, size_t profilepercone ) { m_coneanglesep = coneanglesep; m_numcones = numcones; m_profilepercone = profilepercone; }
		size_t                                          GetMaxPolarizationOrder ( ) const { return m_maxPolarizationOrder;}
		SKTRAN_HR_PolHOType                             GetPolarizationHigherOrderBehaviour ( ) const { return m_polHigherOrderBehaviour;}
		double 											GetPolarizationHigherOrderFraction () const { return m_polarizationHigherOrderFraction;}
		void											SetPrecacheWavel( const std::vector<double>& wavel ) { m_precachewavel = wavel; }
		void											SetScatterResolution(double res) { m_scatres = res; }

		void											SetForceCacheUpdates( bool update ) { m_forcecacheupdates = update; }
		void											SetPolarizationHigherOrderFraction( double fraction ) { m_polarizationHigherOrderFraction = fraction; }

		friend class SKTRAN_HR_Specs_Internal_OpticalPropertiesTable;
};
