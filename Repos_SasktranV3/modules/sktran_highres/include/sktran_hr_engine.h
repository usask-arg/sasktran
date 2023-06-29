//#pragma once

//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Engine		2013-06-07*/
/** \par Pointer convention
 *  I have tried to follow a specific pointer convention for the engine internally
 *  which is the 'recommended' practice for c++11, of course, stack allocated
 *  objects should be preferred whenever possible
 *  - unique_ptr
 *      Implies sole ownership, should be the default when creating a new object
 *  - shared_ptr
 *      Shared ownership, not currently used anywhere.  Should only be used if
 *      a unique_ptr cannot be used.
 *  - raw pointer
 *      Treated as non-owning, if an object receives a raw pointer and holds on to it
 *      the pointer is assumed to stay valid for the lifetime of the object
 *  - weak_ptr
 *      Used when the method does not require the object to exist, but will
 *      change behavior depending on if it exists or not.  Not currently
 *      used anywhere
 *  In addition, checking for nullptr with smart pointers should be done
 *  as if( smartptr ) { }  rather than if ( nullptr != smartpr ) { },
 *  the conversion to bool is an overloaded operator in the smartptr classes.
 *  Null checking for raw pointers should be done if( nullptr != rawptr ) { },
 *  nullptr is a new literal in c++11 to replace NULL or 0 as checks for
 *  null pointers.
 *
 *  \par Passing by reference or passing a pointer
 *  Currently both methods are used without any apparent convention.  In the
 *  future I will try to move the codebase over to the following,
 *  - Pass by Reference
 *      Used when the method does not hold a copy of the object, and is thus
 *      only used locally in the method
 *  - Pass by pointer
 *      Used when the method holds on to a copy of the object, the ownership
 *      is determined from the pointer conventions outlined above.  Internally,
 *      unless explicitly stated we will assume that passed pointer are not nullptr,
 *      
 *  Any future additions will try to follow this convention
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Engine : public SKTRAN_Engine_Base
{
	private:
		OpticalTablePtr											m_opticalpropertiestable;
		SKTRAN_HR_Specs_Internal_Core							m_internalspecs;
		RayFactoryPtr											m_diffuserayfactory;
		RayFactoryPtr											m_linesofsightrayfactory;
		RayFactoryPtr											m_solarrayfactory;
		OptIntegratorPtr										m_optintegrator;
		SrcIntegratorPtr										m_srcintegrator;
		
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	m_coordinates;
		SKTRAN_RayTracingRegionManager							m_raymanager;
		SolarTablePtr											m_solartable;
        EmissionTablePtr                                        m_emissiontable; 
        SolarTablePtr											m_diffusesolartable;
		SKTRAN_HR_LinesOfSightTable								m_linesofsighttable;
		nxVector												m_sun;
        std::unique_ptr<SKTRAN_HR_Diffuse_Table_CPU>			m_diffusetable;
		SKTRAN_HR_Thread_Manager								m_threadmanager;
		std::vector<SKTRAN_Source_Term*>						m_diffusesources;
		bool													m_calcwf;
		std::vector<std::vector<SKTRAN_HR_WF_Ray>>				m_wf;
        std::vector<std::vector<SKTRAN_HR_WF_Ray_Polarized>>	m_wf_pol;

        std::unique_ptr<SKTRAN_HR_WF_Extinction_Table>			m_wf_extinction;

		size_t													m_numcossza;
		bool 													m_prefillsolartransmission;
		bool 													m_usesolartableforsinglescattering;

		bool													m_enablelosscattering;

	private:
		virtual void								ReleaseResources						();
		virtual bool								CreateSolarTable						();
        virtual bool                                CreateEmissionsTables                   ();
		virtual bool								CreateDiffuseSolarTable					();
		virtual bool								CreateDiffuseSources					();

	public:
													SKTRAN_HR_Engine();
		virtual									   ~SKTRAN_HR_Engine();

		virtual bool								ConfigureModel							(       SKTRAN_SpecsUser_Base&			modelspecifications,
																							  const SKTRAN_LineOfSightArray_V21&	linesofsight,
																							  size_t								numthreads );

		virtual bool								CalculateRadiance						( std::vector<SKTRAN_StokesScalar>*		losradiance,
																							  double								wavelen,
																							  size_t								numordersofscatter,
																							  SKTRAN_AtmosphericOpticalState_V21*	opticalstate,
																							  std::vector<skRTStokesVector>*        losvector=nullptr,
																							  bool									updateclimatology = false,
																							  SKTRAN_DiagnosticInterface*			diag = NULL) override;

		virtual bool								SetSun									( nxVector sun ) { m_sun = sun; return true; }
		nxVector									GetSun									() const { return m_sun; }
		const SKTRAN_HR_WF_Ray&						WFForRay								( size_t speciesidx, size_t rayidx ) const { return m_wf[rayidx][speciesidx]; }
        const SKTRAN_HR_WF_Ray_Polarized&			PolarizedWFForRay						( size_t speciesidx, size_t rayidx ) const { return m_wf_pol[rayidx][speciesidx]; }


    bool										CalculateWeightingFunctions				( double wavelen, const std::vector< SKTRAN_Source_Term* > sources, SKTRAN_AtmosphericOpticalState_V21& opticalstate );
        bool										CalculateWeightingFunctionsPolarized	( double wavelen, const std::vector< SKTRAN_Source_Term* > sources, SKTRAN_AtmosphericOpticalState_V21& opticalstate );


    virtual GEODETIC_INSTANT					ReferencePoint							( ) const;

		const SKTRAN_HR_Specs_Internal_Core&		InternalSpecs							() const { return m_internalspecs; };
		SKTRAN_HR_Specs_Internal_Core*				InternalSpecsVar						() { return &m_internalspecs; };
		const SKTRAN_CoordinateTransform_V2*		Coordinates								()		 { return m_coordinates.get();}
		SKTRAN_HR_LinesOfSightTable&				LinesOfSight							() { return m_linesofsighttable; }
		SKTRAN_RayTracingRegionManager*				RayTracingRegionManager					() { return &m_raymanager;}	

		void										SetNumCosSza							(size_t num) { m_numcossza = num; }
		void 										SetPrefillSolarTransmission 			(bool prefill) { m_prefillsolartransmission = prefill; }
		void 										SetUseSolarTableForSingleScattering     (bool use) { m_usesolartableforsinglescattering = use; }
		void 										SetNumSolarTableSZA						(int numsza) { m_numcossza = numsza; }

		void										SetUseLOSScattering						(bool use) { m_enablelosscattering = use; }

};
