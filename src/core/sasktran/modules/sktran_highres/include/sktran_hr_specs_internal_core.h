//#include "sktran_hr_internals.h"

class SKTRAN_HR_Diffuse_Table;
class SKTRAN_HR_Diffuse_Table_CPU;
class SKTRAN_HR_Diffuse_Table_Base;

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Core		2013-08-23*/
/** 
 *  Factory class is responsible for creating and configuring objects required
 *  by the engine.  a SKTRAN_HR_Specs_User object is passed in first to
 *  initialize any user settings, then base class pointers are passed in and
 *  set to created derived objects.
 *
 *  Class is also responsible for creating the diffuse table, as the diffuse
 *  table creation requires information from many sources it was placed here
 *  rather than the internal_specs_diffuse class
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Specs_Internal_Core
{
	private:
		size_t														m_numprofiles;
		size_t														m_numofflook;
		size_t														m_numinterp;
		double														m_angleofflook;
		double														m_diffuseheightres;				// Diffuse height resolution, value copied from  SKTRAN_HR_Specs_User_Diffuse
		double														m_diffusemaxheight;				// Diffuse maximum height,    This value is usually NAN . It is only provided so we can get agreement with old SASKTRAN code
		std::vector<double>											m_manualdiffuseheights;			// Users can specify there own diffuse heights if they wish. They should be in meters and should be ascending 
		std::vector<nxVector>										m_diffuseprofilelocations;      // User specified manual diffuse profile locations (by xyz)
		std::vector<double>											m_diffuseprofilelatlons;        // User specified manual diffuse profile locations (by xyz)
		std::vector<double>											m_diffuseprofileszas;           // User specified diffuse profile SZAs
		std::vector<double>											m_diffuseprofilelospositions;   // User specified diffuse profile locations: linear scale such that 0.0 puts it under the entry point, 1.0 puts it under/over the exit/ground point (extrapolation permissible)
		nxVector													m_diffuseplanereference;		// Reference vector for user specified diffuse profile locations in a plane
		nxVector													m_diffuseplanenormal;			// Normal vector for user specified diffuse profile locations in a plane
		std::vector<double>											m_diffuseplaneangles;			// Angles from reference vector in degress for user specified diffuse profile locations in a plane

		enum class InternalDiffusePlacementType						{ ManualLocation, ManualLatLon, ManualSZA, ManualLOS, ManualPlane, LinearLOS2D, LinearLOS1D, LinearLOS1DExit, LinearLOS1DEntry, LinearLOS1DEntryExit, LinearSZA, LinearSZAForceTP, SmartSZA, Undefined }; // changes made here should be updated in DiffuseLocationTypeChar
		InternalDiffusePlacementType								m_internaldiffuseplacementtype;

		std::vector<size_t>											m_trackedscatords;
		std::vector<size_t>											m_trackeddiffprofs;
		bool														m_track;
		SKTRAN_HR_RayTracingRegionManager							m_raymanager;
		const SKTRAN_LineOfSightArray_V21*							m_linesofsight;
		std::shared_ptr< const SKTRAN_CoordinateTransform_V2>		m_coords;
		SKTRAN_HR_DiffuseProfilePlacementType						m_userdiffuseplacementtype;
		bool														m_in3dmode;
		bool														m_calcwf;
		int 														m_scatterorder;
		size_t                                                      m_maxPolarizationOrder;
		SKTRAN_HR_PolHOType                                         m_polHigherOrderBehaviour; 
		double 														m_polarizationHigherOrderFraction;
		SKTRAN_HR_Specs_Internal_RayTracer							m_raytracerspecs;
		SKTRAN_HR_Specs_Internal_Integrator							m_integratorspecs;
		SKTRAN_HR_Specs_Internal_OpticalPropertiesTable				m_opttablespecs;
		SKTRAN_HR_Specs_Internal_Diffuse							m_diffusespecs;
		SKTRAN_HR_Specs_Internal_wf									m_wfspecs;

	private:												
		const char*													DiffuseLocationTypeChar			( InternalDiffusePlacementType type ) const;
		bool														DiffusePlacementWarning			( InternalDiffusePlacementType selected, InternalDiffusePlacementType overridden ) const;

	private:
		bool														ReleaseResources				();
		bool														RotateStartAndEnd				( const HELIODETIC_UNITVECTOR& look, const HELIODETIC_UNITVECTOR& in, HELIODETIC_UNITVECTOR& out, double theta );
		bool														CreateDiffusePoints				( SKTRAN_HR_Diffuse_Table_CPU& diffusetable );
		bool														CreateDiffuseLocations          ( std::vector<HELIODETIC_POINT>& locations );
		bool														CreateDiffuseHeights			( std::vector<double>* heights );
		bool														CreateLinearDiffuseLocations    ( std::vector<HELIODETIC_POINT>& locations );
		bool														CreateOffLOSDiffuseLocations	( std::vector<HELIODETIC_POINT>& locations );
		bool														CreateLinearSZADiffuseLocations	( std::vector<HELIODETIC_POINT>& locations );
		bool														CreateManualDiffuseLocations	( std::vector<HELIODETIC_POINT>& locations );
		bool														CreateLatLonDiffuseLocations	( std::vector<HELIODETIC_POINT>& locations );
		bool														CreateManualLOSDiffuseLocations	( std::vector<HELIODETIC_POINT>& locations );
		bool														CreatePlaneDiffuseLocations		( std::vector<HELIODETIC_POINT>& locations );
		bool														CreateDiffuseTableType			( std::unique_ptr<SKTRAN_HR_Diffuse_Table_CPU>& diffusetable, bool usecache );
        bool                                                        CreateDiffuseTableAValues       ( std::unique_ptr< SKTRAN_HR_Diffuse_Table_CPU >& diffusetable );
		bool														CreateSZAs						( std::vector<double>& szas ) const;
		bool														CalcReferencePoints				( nxVector& in, nxVector& out );
		const SKTRAN_CoordinateTransform_V2*						CoordinatePtr					() const { return m_coords.get();}

	public:
																	SKTRAN_HR_Specs_Internal_Core	();
		virtual													   ~SKTRAN_HR_Specs_Internal_Core	();
		SKTRAN_HR_Specs_Internal_RayTracer&							RayTracerSpecs					()			{ return m_raytracerspecs; }
		SKTRAN_HR_Specs_Internal_Integrator&						IntegratorSpecs					()			{ return m_integratorspecs; }
		SKTRAN_HR_Specs_Internal_OpticalPropertiesTable&			OpticalPropertiesSpecs			()			{ return m_opttablespecs; }
		const SKTRAN_HR_Specs_Internal_OpticalPropertiesTable&		OpticalPropertiesSpecs			() const	{ return m_opttablespecs; }
		SKTRAN_HR_Specs_Internal_Diffuse&							DiffuseSpecs					()			{ return m_diffusespecs; }
		SKTRAN_HR_Specs_Internal_wf&								WeightingFunctionSpecs			()			{ return m_wfspecs; }
		const SKTRAN_HR_Specs_Internal_wf&							WeightingFunctionSpecs			() const	{ return m_wfspecs; }
		bool														Is3dMode						()			{ return m_in3dmode; }
		bool														CalcWf							() const	{ return m_wfspecs.DoWfCalculation(); }
		SKTRAN_HR_RayTracingRegionManager&							RayManager						()			{ return m_raymanager; }
		int 														ScatterOrder					() const    { return m_scatterorder; }
        int                                                         MaxPolarizationOrder            () const    { return m_maxPolarizationOrder; }
//		double														GroundShift						() const    { return m_groundshiftalt; }

	public:
		/*virtual*/ bool											Configure						( const SKTRAN_Specifications_Base& specs, const SKTRAN_LineOfSightArray_V21& linesofsight );
		/*virtual*/ bool											CreateCoordinates				( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>* coords, nxVector& sun, const SKTRAN_LineOfSightArray_V21& linesofsight, double surfaceHeight, double toaHeight, bool nadir_referencepoint_onground );
		/*virtual*/ bool											CreateDiffuseTable				( std::unique_ptr<SKTRAN_HR_Diffuse_Table_CPU>& diffusetable );
};
