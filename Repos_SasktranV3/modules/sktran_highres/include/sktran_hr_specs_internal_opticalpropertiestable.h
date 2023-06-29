//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_OpticalPropertiesTable		2014-10-30*/
/** Factory class for creating optical property tables within the HR engine.
 *
 *  The class is also responsible for adding additional information to the
 *  raytracer which may be dependant on how the optical property table is
 *  set up
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Specs_Internal_OpticalPropertiesTable
{
	private:
		SKTRAN_HR_OpticalPropertiesTableType    m_optproptype;
		size_t                                  m_maxPolarizationOrder;
		SKTRAN_HR_PolHOType                     m_polHigherOrderBehaviour;
		double 									m_polarizationHigherOrderFraction;
		double									m_scatgridresolution;
		double									m_heightspacing;
		size_t									m_numprofiles;
		std::vector<double>						m_anglegrid;
		std::vector<double>						m_heightgrid;
		nxVector								m_normal;
		nxVector								m_reference;
		double									m_coneanglesep;
		size_t									m_numcones;
		size_t									m_profilepercone;

		std::vector<double>						m_precachewavel;

		bool									m_forcecacheupdates;

		const SKTRAN_HR_RayTracingRegionManager*   m_raymanager;

		SKTRAN_HR_AtmosphereHasDelta            m_atmosphereHasDelta;
		bool 									m_norefractionsinglescatter;

	private:
		virtual bool                            CreatePolarizationObject( std::unique_ptr< SKTRAN_PolarizationProperties_Base >& polobject ) const;
		virtual bool							Create1dTable			( OpticalTablePtr& table, const SKTRAN_CoordinateTransform_V2& coords, double toaHeight );
		virtual bool							Created3DUnitSphereTable( OpticalTablePtr&, const SKTRAN_CoordinateTransform_V2& coords, double toaHeight );
		void									MakeDefaultAngleGrid    ();
		void									MakeNormalAndReferenceFromLOS();
		virtual bool							MakeHeightGrid			( SKTRAN_GridDefOpticalPropertiesRadii_V21& heightgrid, double toaHeight );
		virtual bool							MakeUnitSphere			( SKTRAN_UnitSphere_V2** unitsphere );
		virtual bool							MakeDelaunaySphere		( SKTRAN_UnitSphere_V2** unitsphere, const SKTRAN_CoordinateTransform_V2& coords );
		virtual bool							MakeLOSPlaneSphere      ( SKTRAN_UnitSphere_V2** unitsphere, const SKTRAN_CoordinateTransform_V2& coords );
		virtual nxVector						CalcRotatedVector		( const nxVector& tangent, double zenith, double azimuth ) const;

		bool									AddPlaneInformationToRayTracer( SKTRAN_RayTracer_Straight_Generic& raytracer, const SKTRAN_CoordinateTransform_V2& coords ) const;
		bool									AddDelaunayInformationToRayTracer( SKTRAN_RayTracer_Straight_Generic& raytracer, const SKTRAN_CoordinateTransform_V2& coords ) const;

	public:
												SKTRAN_HR_Specs_Internal_OpticalPropertiesTable();
		virtual								   ~SKTRAN_HR_Specs_Internal_OpticalPropertiesTable();

		virtual bool							Configure				( const SKTRAN_HR_Specs_User& specs, const SKTRAN_HR_RayTracingRegionManager* raymanager );
		virtual bool							CreateOpticalTable		( OpticalTablePtr& table, const SKTRAN_CoordinateTransform_V2& coords, double toaHeight );
		virtual bool							ConfigureDefaults		();
		virtual bool							In3dMode				() const { return m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_3D_UnitSphere ||
																					m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_LOSPlane ||
																					m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_SZA; }
		size_t                                  GetMaxPolarizationOrder    ( ) const { return m_maxPolarizationOrder;}
        double                                  GetHeightSpacing        ( ) const { return m_heightspacing; }

		SKTRAN_GridDefSLON_V21					MakeSLONGrid			( const SKTRAN_CoordinateTransform_V2& coords ) const;
		virtual bool							MakeScatterAngleGrid(SKTRAN_GridDefScatterAngle_V21& scattergrid) const;

		const std::vector<double>&				PrecacheWavel			() const { return m_precachewavel; }

		bool									AddOpticalInformationToRayTracer( SKTRAN_RayTracer_Straight_Generic& raytracer, const SKTRAN_CoordinateTransform_V2* coords ) const;
		friend class SKTRAN_HR_Specs_Internal_wf;
		void									SetNoRefractionSingleScatter( bool norefrac ) { m_norefractionsinglescatter = norefrac;}
};
