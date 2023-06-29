#pragma once

/*-----------------------------------------------------------------------------
 *					class SKTRANSO_Quadrature_TLS_V2_Legacy		2010-3-8*/
/** A specific quadrature instance for use in the Sasktran radiative transfer model.
 *	This class provides single scatter using interpolation of the
 *	Solar transmission Table. Thus it is not efficient for single scatter calculations
 *	as it requires the entire solar transmission table to be initialized.
 *
 *	The code calculates optical depth along a ray using a linear interpolation of the
 *	extinction cross-section with altitude. However it integrates the source function
 *	terms across a ray segment using a fixed source function and a fixed extinction
 *	cross-section (we should look into improving this).
 *
 *		#- Atmospheric Single scatter from interpolation of Solar transmisison table
 *		#- Ground Single scatter from interpolation of the ????
 *		#- Atmospheric diffuse signal from interpolation of the Diffuse Points Table
 *		#- Ground diffuse scatter from interpolation of the Ground Diffuse Points Table
 *
 **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_Quadrature_TLS_V2_Legacy : public SKTRANSO_Quadrature_TLS_V21
{	
	private:
		SKTRANSO_JIndexArray								m_jindexworkspace;				//!< A "large" work space to pre-calculate each Jindex table for any ray before copying to specific instance.
		double*												m_cellfactors;					//!< Reserved for use by method FillCellFactorBuffer, A work buffer [maxshells] for storing cell factor terms.
		double*												m_shelltransmissionsbuffer;		//!< Reserved for method CalculateRayTransmission.    A work buffer [maxshells] for storing shell transmissions, /f$e^{-\tau_0}/f$, from edge of cell to observer
		double*												m_scratchbuffer;				//!< A temporary work buffer [maxquadpoints] for storing quadrature point weights along a ray.
		size_t												m_maxshellsalongray;			//!< Current size of the m_shelltransmissionsbuffer work buffer.
		size_t												m_maxquadpointsalongray;		//!< Current size of the m_scratchbuffer work buffers.
		const SKTRAN_AlbedoBRDF_V2*							m_albedobrdf;
		const SKTRAN_EngineDiffuseTables*					m_enginetables;
		const SKTRAN_TableOpticalProperties_V21*			m_opticalprops;					//!< The optical properties table;
		std::shared_ptr< const SKTRAN_CoordinateTransform_V2>	m_coordinates;


	private:
		size_t							MaxQuadraturePointsAlongRay							( const SKTRAN_SpecsInternal_RayTracing_V21* raytracingspecs);
		size_t							MaxQuadraturePointsPerCell							() const { return 1;}
		void							ReleaseTemporaryWorkAreas							();
		bool							FillCellFactorBuffer								( const SKTRANSO_RayInternalOptical* ray, size_t* numweights, bool usecachedvalue );
		double							CellTransmissionFactor								( const HELIODETIC_POINT& point, double celllengthmeters ) ;
		double							OpticalDepthOfSegment								( size_t cellidx, const SKTRAN_RayStorage_Base* storage) ;
		double							OpticalDepthOfCell									( const SKTRAN_RayStorage_Base* ray, size_t cellidx );
		bool							CreateJIndexTableForQuadraturePointsAlongRay		( const SKTRANSO_JindexTableBase* tablesource, const SKTRANSO_RayInternalGeometry* ray, SKTRANSO_JIndexArray* jindex);
		bool							CreateJIndexTableForGroundPoint						( const SKTRANSO_JindexTableBase* tablesource, const SKTRANSO_RayInternalGeometry* ray, SKTRANSO_JIndexArray* jindex, bool isSingleScatterTerm);
		const SKTRANSO_JindexTableBase*		RaySolarTransmissionTable							( const SKTRANSO_RayInternalGeometry* ray);

		bool							QuadratureWeightsForAtmosphericEmission				( const SKTRANSO_RayInternalOptical*		ray,
																							  double*									weights,
																							  size_t*									numweights);

		bool							QuadratureWeightsForMSGroundPoint					( const SKTRANSO_RayInternalOptical*		ray,
																							  double*									weights,
																							  size_t*									numweights
																							 );

		bool							QuadratureWeightsForSSGroundPoint					( const SKTRANSO_RayInternalOptical*			ray, 
																							  const SKTRAN_TableOpticalProperties_V21*	opticalprops,
																							  double*									weights,
																							  size_t*									numweights);

		bool							GetCellQuadraturePointAltitudeSzaAndLookTowardsObserver	(const SKTRAN_RayStorage_Base* ray, size_t cellidx, HELIODETIC_POINT* point, HELIODETIC_UNITVECTOR* look, double* celllengthmeters ) const;

		bool							QuadratureWeightsForAtmosphericSS					( const SKTRANSO_RayInternalOptical* ray, double* weights, size_t* numweights             ) ;

		bool							Initialize											( const SKTRAN_SpecsInternal_V21* modelspecifications );

	public:
										SKTRANSO_Quadrature_TLS_V2_Legacy						( const SKTRAN_SpecsInternal_V21*		modelspecifications,
																							  const SKTRAN_EngineDiffuseTables*		modeltables);

		virtual						   ~SKTRANSO_Quadrature_TLS_V2_Legacy					();
		
		virtual bool					SetOpticalProps										( const SKTRAN_TableOpticalProperties_V21*	optprop)override;


		virtual bool					CalculateRayTransmission							( SKTRANSO_RayInternalOptical*			ray,
																							  double*							transmisison,
																							  bool								totaltransmissiononly,
																							  bool								usecachedtransmission
																							)override;

		virtual bool					CreateJIndexTable_AtmosphericSingleScatter			( const SKTRANSO_RayInternalGeometry*	ray,
																									SKTRANSO_JIndexArray*		jindex
																							)override;

		virtual bool					CreateJValueTable_AtmosphericSingleScatter			(       SKTRAN_RayInternalDiffuseOptical_V2*			ray,
																							  const SKTRANSO_JIndexArray&					jindex,
																							        SKTRAN_JValueTable_V21*					jvalue,
																							  bool											usecachedtransmission,
																							  bool											usecachedcellfactors 
																							)override;

		virtual bool					CreateJIndexTable_AtmosphericDiffuseScatter			( const SKTRANSO_RayInternalGeometry*	ray,
																					        SKTRANSO_JIndexArray*		jindex
																							)override;

		virtual bool					CreateJValueTable_AtmosphericDiffuseScatter			(       SKTRANSO_RayInternalOptical*	ray,
																							  const SKTRANSO_JIndexArray&		jindex,
																								    SKTRAN_JValueTable_V21*		jvalue,
																							  bool								usecachedtransmission,
																							  bool                              usecachedcell_transmissions
																							)override;

		virtual bool					CreateJIndexTable_InterpolateGroundDiffuseScatter	( const SKTRANSO_RayInternalGeometry*			ray,
																							        SKTRANSO_JIndexArray*					jindex
																							)override;

		virtual bool					CreateJValueTable_InterpolateGroundDiffuseScatter	(       SKTRAN_RayInternalDiffuseOptical_V2*	ray,
																							  const SKTRANSO_JIndexArray&					jindex,
																							        SKTRAN_JValueTable_V21*					jvalue,
																							  bool											usecachedtransmission
																							)override;

		virtual bool					CreateJIndexTable_GroundSingleScatter				( const SKTRANSO_RayInternalGeometry*			ray,
																							        SKTRANSO_JIndexArray*					jindex
																							)override;

		virtual bool					CreateJValueTable_GroundSingleScatter				(       SKTRAN_RayInternalDiffuseOptical_V2*	ray,
																							  const SKTRANSO_JIndexArray&					jindex,
																							        SKTRAN_JValueTable_V21*					jvalue,
																							  bool											usecachedtransmission
																							)override;

		
		virtual bool					CreateJIndexTable_AtmosphericEmissions				( const SKTRANSO_RayInternalGeometry*			ray,
																									SKTRANSO_JIndexArray*					jindex
																							) override;

		virtual bool					CreateJValueTable_AtmosphericEmissions				(		SKTRAN_RayInternalDiffuseOptical_V2*	ray,
																							  const SKTRANSO_JIndexArray&					jindex,
																									SKTRAN_JValueTable_V21*					jvalue,
																									bool									usecachedtransmission,
																									bool									usecachedcellfactors
																							) override;

		virtual bool					CreateJIndexTable_GroundEmissions					( const SKTRANSO_RayInternalGeometry*			ray,
																									SKTRANSO_JIndexArray*					jindex
																							) override;

		virtual bool					CreateJValueTable_GroundEmissions					(		SKTRAN_RayInternalDiffuseOptical_V2*	ray,
																							  const SKTRANSO_JIndexArray&					jindex,
																									SKTRAN_JValueTable_V21*					jvalue,
																									bool									usecachedtransmission
																							) override;


		virtual bool					CreateJIndexTable_GroundPointIncomingRaysUpwardFlux	( SKTRAN_GroundPointDiffuseGeometry_V21* groundpoint )override;

		virtual bool					CreateJValueTable_GroundPointIncomingRaysUpwardFlux	( SKTRANSO_GroundPointDiffuseOptical*	groundpoint,
																							  const SKTRAN_TableOpticalProperties_V21* optprop
																							)override;

};



