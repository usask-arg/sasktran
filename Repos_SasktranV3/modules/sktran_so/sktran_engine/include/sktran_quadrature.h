
class SKTRAN_GroundPointDiffuseGeometry_V21;
class SKTRANSO_GroundPointDiffuseOptical;

/*-----------------------------------------------------------------------------
 *					SKTRANSO_Quadrature_TLS_V21		2010-3-8*/
/** A class that is the workhorse for all of the quadrature/integration calculations
 *	used in the Sasktran model. The class is designed to be replaced by user
 *	specific implementations although Sasktran provides two implementations.
 *
 *	The Sasktran model creates one instance of this class for each processing
 *	thread. It is therefore important that each quadrature implementation not
 *	use static or global variables of any kind (without mutex protection at least) as
 *	it is extremely likely that several instances will be running simultaneously
 *	inside any given model run.
 *
 *	This class is responsible for all quadrature based calculations within the model.
 *	This version of Sasktran considers the following quadrature terms,
 *
 *	#-	Optical Depth calculations along a ray. 
 *	#-	Multiple scatter from the atmosphere.
 *	#-  Multiple scatter from the ground.
 *	#-	Single scatter from the atmosphere.
 *	#-  Single scatter from the ground.
 *
 *	\section Ray Transmission
 *	\section Diffuse Atmospheric Scatter
 *	The quadrature object provides the following methods associated with the
 *	diffuse scatter off the atmosphere,
 *		- #CreateJIndexTable_AtmosphericDiffuseScatter
 *		- #CreateJValueTable_AtmosphericDiffuseScatter
 *
 *	\section Diffuse Ground Points
 *	The quadrature object provides the following methods associated with the
 *	diffuse ground points,
 *		- #CreateJIndexTable_InterpolateGroundDiffuseScatter
 *		- #CreateJValueTable_InterpolateGroundDiffuseScatter
 *		- #CreateJIndexTable_GroundPointDiffuseUpwardFlux
 *		- #CreateJValueTable_GroundPointDiffuseUpwardFlux
 *
 *	\section Single Scatter from the Ground
 *		- #CreateJIndexTable_GroundSingleScatter
 *		- #CreateJValueTable_GroundSingleScatter

 *	The quadrature object provides the following methods associated with the
 *	single scatter off the ground,
 *	\par Caching Strategies
 *	The class provides a few hooks that allow the callers some flexibility in calculation
 *	strategies. All caching strategies must be thread safe! The following strategies can
 *	be controlled by the caller,
 *
 *	#-	Only calculate the optical depth of a ray. This strategy is beneficial for the
 *		rays in the Solar Transmission Table that only require the total transmission of the ray.
 *		The strategy avoids the unnecessary conversion of the optical depth of a shell
 *		boundary to a transmission factor, avoids an exponential function call. 
 *
 *	#-	Cache the optical depth of the start of each segment (expressed as a transmission) from 
 *		within an internal buffer. This is important as all quadratures along the ray require
 *		the optical depth of the ray from the start of the ray. Thus it is important that the
 *		option is properly implemented by both the derived class and the caller
 *		when processing all of the Jindex ->JValue tables associated with a given ray. However,
 *		the user can extend the optimization to multiple rays at the same point and solar zenith
 *		angle in a homogeneous atmosphere
 *		
 *	#-	Cache the cell factors for each segment of the ray. The same factors can often be
 *		applied to single scatter or diffuse source functions and do not need recalculating
 *		when processing the Jindex ->JValue tables associated with a given ray. The user
 *		can extend this optimization to multiple rays at the same point and solar zenith
 *		angle in a homogeneous atmosphere
 **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_Quadrature_TLS_V21 : public SKTRAN_Quadrature_TLS_Base
{
	private:
		size_t					m_numpointsprocessed;					// Stores the number of points processed in the last thread action (used for tuning performance

	public:
								SKTRANSO_Quadrature_TLS_V21(){}
		virtual				   ~SKTRANSO_Quadrature_TLS_V21(){}
		void					SetNumPointsProcessed( size_t numpoints)  { m_numpointsprocessed = numpoints;}
		size_t					NumPointsProcessed() const { return m_numpointsprocessed;}


		/*-----------------------------------------------------------------------------
		 *						SetOpticalProps								2010-2-19*/
		/**
		 *
		**/
		/*-----------------------------------------------------------------------------*/
		virtual bool			SetOpticalProps										( const SKTRAN_TableOpticalProperties_V21*	optprop) = 0;

		/*-----------------------------------------------------------------------------
		 *						CalculateRayTransmission					2010-2-19*/
		/**
		 *
		**/
		/*-----------------------------------------------------------------------------*/
		virtual bool			CalculateRayTransmission							( SKTRANSO_RayInternalOptical*			ray,
																					  double*							totaltransmission,
																					  bool								totaltransmissiononly,
																					  bool								usecachedtransmission
																					) = 0;

		/*-----------------------------------------------------------------------------
		 *						CreateJIndexTable_AtmosphericSingleScatter	2010-2-19*/
		/**
		 *
		**/
		/*-----------------------------------------------------------------------------*/
		virtual bool			CreateJIndexTable_AtmosphericSingleScatter			( const SKTRANSO_RayInternalGeometry*			ray,
																					        SKTRANSO_JIndexArray*					jindex
																					) = 0;

		/*-----------------------------------------------------------------------------
		 *																	2010-2-19*/
		/**
		 *
		**/
		/*-----------------------------------------------------------------------------*/
		virtual bool			CreateJValueTable_AtmosphericSingleScatter			(       SKTRAN_RayInternalDiffuseOptical_V2*	ray,
																					  const SKTRANSO_JIndexArray&					jindex,
																					        SKTRAN_JValueTable_V21*					jvalue,
																					  bool											usecachedtransmission,
																					  bool											usecachedcellfactors
																					) = 0;


		/*-----------------------------------------------------------------------------
		 *																	2010-2-19*/
		/**
		 *
		**/
		/*-----------------------------------------------------------------------------*/

		virtual bool			CreateJIndexTable_AtmosphericDiffuseScatter			( const SKTRANSO_RayInternalGeometry*	ray,
																					        SKTRANSO_JIndexArray*		jindex
																					) = 0;

		/*-----------------------------------------------------------------------------
		 *																	2010-2-19*/
		/**
		 *
		**/
		/*-----------------------------------------------------------------------------*/
		virtual bool			CreateJValueTable_AtmosphericDiffuseScatter			(       SKTRANSO_RayInternalOptical*	ray,
																					  const SKTRANSO_JIndexArray&		jindex,
																					        SKTRAN_JValueTable_V21*		jvalue,
																					  bool								usecachedtransmission,
																					  bool                              usecachedcell_transmissions
																					) = 0;

		/*-----------------------------------------------------------------------------
		 *																	2010-2-19*/
		/**
		 *
		**/
		/*-----------------------------------------------------------------------------*/
/*
	virtual bool			CreateJIndexTable_LOSSingleScatter					( SKTRANSO_RayLOSGeometry_V21*	ray,
																					  SKTRANSO_JIndexArray*		jindex,
																					  SKTRANSO_JIndexArray*		groundjindex 
																					 ) = 0;

*/
		/*-----------------------------------------------------------------------------
		 *																	2010-2-19*/
		/**
		 *
		**/
		/*-----------------------------------------------------------------------------*/
/*
	virtual bool			CreateJValueTable_LOSSingleScatter					(       SKTRANSO_RayInternalOptical*	ray,
																					  const SKTRANSO_JIndexArray&		jindex,
																					        SKTRAN_JValueTable_V21*		jvalue
																					) = 0;
*/

		/*-----------------------------------------------------------------------------
		 *						CreateJIndexTable_InterpolateGroundDiffuseScatter	 */
		/**	Used by rays during the geometry section to interpolate the
		 *	diffuse upward flux from the diffuse ground points table. The interpolation
		 *	is returned as a JIndexTable.  The table is converted to a JValueTable for a
		 *	specific wavelength by the matching function #CreateJIndexTable_InterpolateGroundDiffuseScatter
		 *
		 *	\param ray
		 *	The ray that requires the upward flux calculation.
		 *
		 *	\param jindex
		 *	returns the JIndex table for the upward flux. The table may also include the
		 *	angular response of the albedo BRDF (as we only support a wavelength independent BRDF)
		**/
		/*-----------------------------------------------------------------------------*/

		virtual bool			CreateJIndexTable_InterpolateGroundDiffuseScatter	( const SKTRANSO_RayInternalGeometry*	ray,
																					        SKTRANSO_JIndexArray*					jindex
																					) = 0;

		/*-----------------------------------------------------------------------------
		 *						CreateJValueTable_GroundDiffuseScatter				 */
		/**	Converts the Jindex table created by #CreateJIndexTable_InterpolateGroundDiffuseScatter
		 *	into a JValue table for a specific wavelength. The quadrature object will
		 *	return the JValue table so provides the radiance due to diffuse ground
		 *	scatter at the observer's end of the ray. It will include the transmission
		 *	of the ray.
		**/
		/*-----------------------------------------------------------------------------*/

		virtual bool			CreateJValueTable_InterpolateGroundDiffuseScatter	(       SKTRAN_RayInternalDiffuseOptical_V2*	ray,
																					  const SKTRANSO_JIndexArray&					jindex,
																					        SKTRAN_JValueTable_V21*					jvalue,
																					  bool											usecachedtransmission
																					) = 0;

		/*-----------------------------------------------------------------------------
		 *						CreateJIndexTable_GroundPointDiffuseUpwardFlux	     */
		/**	Used by the diffuse ground points table during the geometry section of
		 *	code to creates the diffuse upward flux at each ground point. This will
		 *	normally evaluate all of the angular areas over the unit hemisphere of
		 *	the incoming rays as well as convert the incoming ray to a downward flux. This
		 *	method is matched with #CreateJValueTable_GroundPointDiffuseUpwardFlux which
		 *	converts the Jindex table to a wavelength specific JValue table.
		 *
		 *	\param groundpoint
		 *	The ground point that will be processed.
		 *
		 *	\param jindex
		 *	Returns the Jindex Table for the diffuse upward flux calculation. This code is
		 *	wavelength independent.
		**/
		/*-----------------------------------------------------------------------------*/

		virtual bool			CreateJIndexTable_GroundPointIncomingRaysUpwardFlux	( SKTRAN_GroundPointDiffuseGeometry_V21* groundpoint 
																					) = 0;

		/*-----------------------------------------------------------------------------
		 *																			 */
		/**
		 *
		**/
		/*-----------------------------------------------------------------------------*/
		virtual bool			CreateJValueTable_GroundPointIncomingRaysUpwardFlux	( SKTRANSO_GroundPointDiffuseOptical*	groundpoint,
																					  const SKTRAN_TableOpticalProperties_V21* optprop
																					) = 0;

		/*-----------------------------------------------------------------------------
		 *																	2010-2-19*/
		/**
		 *
		**/
		/*-----------------------------------------------------------------------------*/

		virtual bool			CreateJIndexTable_GroundSingleScatter				( const SKTRANSO_RayInternalGeometry*			ray,
																					        SKTRANSO_JIndexArray*					jindex
																					) = 0;

		/*-----------------------------------------------------------------------------
		 *																	2010-2-19*/
		/**
		 *
		**/
		/*-----------------------------------------------------------------------------*/
		virtual bool			CreateJValueTable_GroundSingleScatter				(       SKTRAN_RayInternalDiffuseOptical_V2*	ray,
																					  const SKTRANSO_JIndexArray&					jindex,
																					        SKTRAN_JValueTable_V21*					jvalue,
																					  bool											usecachedtransmission
																					) = 0;


		virtual bool			CreateJIndexTable_AtmosphericEmissions				( const SKTRANSO_RayInternalGeometry*			ray,
																					        SKTRANSO_JIndexArray*					jindex
																					) = 0;

		virtual bool			CreateJValueTable_AtmosphericEmissions				(     SKTRAN_RayInternalDiffuseOptical_V2*		ray,
																					  const SKTRANSO_JIndexArray&					jindex,
																					        SKTRAN_JValueTable_V21*					jvalue,
																					  bool											usecachedtransmission,
																					  bool											usecachedcellfactors
																					) = 0;

		virtual bool			CreateJIndexTable_GroundEmissions					( const SKTRANSO_RayInternalGeometry*			ray,
																					        SKTRANSO_JIndexArray*					jindex
																					) = 0;

		virtual bool			CreateJValueTable_GroundEmissions					(     SKTRAN_RayInternalDiffuseOptical_V2*		ray,
																					  const SKTRANSO_JIndexArray&					jindex,
																					        SKTRAN_JValueTable_V21*					jvalue,
																					  bool											usecachedtransmission
																					) = 0;


};


/*-----------------------------------------------------------------------------
 *					class SKTRAN_Quadrature_Factory_V21				2010-3-22*/
/** The class used to create the quadrature instances used by the
 *	multiple threads in Sasktran. This factory class will enter the 
 *	Sasktran model through the specifications. The class is small and simple. The
 *	derived factory class should be developed alongside the speciifc derived class instances
 *	of SKTRANSO_Quadrature_TLS_V21.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_Quadrature_Factory_V21 : public nxUnknown 
{

	public:
							SKTRAN_Quadrature_Factory_V21(){}
		virtual			   ~SKTRAN_Quadrature_Factory_V21(){}


		/** Create a new instance of a quadrature object which can be used by a
		 *	single thread in the Sasktran model.
		 *
		 *	\param modelspecifications
		 *	The specifications used for this model calculation. The quadrature object
		 *	can use this parameter to configure itself.
		 *
		 *	\param modeltables
		 *	A pointer to the internal tables used by the engine. The quadrature object
		 *	may use the tables to perform its calculations but cannot modify the tables.
		 *	For example, the quadrature object will normally choose to use the internal
		 *	diffuse points table to calculate the diffuse signal along any ray.
		 *
		 *	\param instance
		 *	Returns a pointer to the the new thread specific quadrature object. The
		 *	object, if not null, has a reference count already added. The user must
		 *	call Release on this object to free all of the internal resources
		 *	used by the object. The Sasktran code will do most of the subsequent initialization
		 *	by calling instance->ConfigureGeometry and instance->ConfigureOptical
		 **/

		virtual bool		CreateThreadInstance( const SKTRAN_SpecsInternal_V21*		modelspecifications,
												  const SKTRAN_EngineDiffuseTables*		modeltables,
												  SKTRANSO_Quadrature_TLS_V21**			instance ) const = 0; 
};
	
