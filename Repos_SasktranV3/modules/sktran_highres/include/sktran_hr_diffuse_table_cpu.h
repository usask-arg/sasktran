//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_CPU		2014-10-30*/
/** The primary implementation of the diffuse table for the HR engine.
 *
 *  Diffuse points are placed in profiles at arbritrary locations on the
 *  surface of the Earth, interpolation is then done between n (by default 3)
 *  nearest profiles.  
 * 
 *  Incoming/outgoing radiances are stored as single vectors for the entire
 *  diffuse field to force contiguous memory, and to allow the problem to
 *  be formulated as successive matrix multiplications.  
 *
 *  Implementation Notes:
 *    - The resulting diffuse field, J(x),  is proportional to Kscat(x),
 *      we compute J = A(x_p)*Kscat(x_p) where x_p is the location of a diffuse point,
 *      and then multiply by Kscat(x) / Kscat(x_p) when DiffuseSource is called. 
 *
 *   Base interface for the diffuse table in the HR engine.  The diffuse table
 *   has two primary tasks, the first is to provide diffuse source information
 *   to the engine, and the second is to compute the diffuse source information
 *   The only requirement is that diffuse points must be indexable internally by
 *   a single index.
 *   
 *   Set up of the diffuse table is done by the thread manager in a multithreaded
 *   environment, the order of function calls done is as follows, functions
 *   marked TS must be thread safe
 *
 *   AllocateDiffusePoints()
 *   ConfigureIndexes()
 *	 ForEach Point
 *		CalcFirstOrderIncomingPoint()	TS
 *		CalcDiffuseIndexesPoint()		TS
 *		CalcScatteringMatrixPoint()		TS
 *   End ForEach
 *   ForEach Point
 *		PostProcessPoint()				TS
 *	 End ForEach
 *	 Scatter()
 *   For ScatteringOrders
 *		ComputeNextOrder()
 *		Scatter()
 *	 End For
 *	 ReleaseResources()
 *		

 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Diffuse_Table_CPU  //: public SKTRAN_HR_Diffuse_Table_Base
{
    //friend class SKTRAN_HR_Specs_Internal_Core; // Specs need to create member objects for this class 

	private:	
		std::unique_ptr< RadStore_Base >                m_radStorage;
        std::unique_ptr< Avals_Base >                   m_Avals;
		std::vector<SKTRAN_HR_Diffuse_Point>			m_diffusepoints;
		size_t											m_groundstartidx;
		std::vector<size_t>								m_profilestartidx;						// Stores where diffuse profiles begin/end in terms of indexes.  The last element of profilestartidx should be the number of non-ground diffuse points
		std::vector< std::vector<double> >				m_diffusealts;
		size_t											m_numprofileinterp;					// The number of diffuse profiles to interpolate ( eg 1, 2 or 3)
		std::unique_ptr<const SKTRAN_UnitSphere_V2>		m_diffuselocationsphere;		// As an option diffuse profiles can be specified on a Unit sphere.
		std::vector<double>								m_kscatatpoints;				// An array of scattering cross-sections at each point. Used exclusively to weight linear interpolations between point. Note ground points are treated differently.
		std::vector<double>								m_szas;
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2> m_coords;
		const SKTRAN_OpticalPropertiesIntegrator_Base*	m_optintegrator;
		const SKTRAN_SourceTermIntegrator_Base*         m_srcintegrator;
		std::shared_ptr<const SKTRAN_RayFactory_Base>	m_rayfactory;
		const std::vector<SKTRAN_Source_Term*>*			m_sources;
		const SKTRAN_TableOpticalProperties_Base*		m_opticaltable;
		std::vector<size_t>								m_profindices;
		std::vector<size_t>								m_scatords;

	private:
		size_t											NumDiffuseProfiles						() const { return m_profilestartidx.size() -1;}
        template< typename Radtype > bool				DiffuseSource_impl						( const SKTRAN_SourceTermQueryObject_Base& qobj,  Radtype& source ) const;
        template< typename Radtype > bool				GroundSource_impl						( const SKTRAN_SourceTermQueryObject_Base& qobj,  Radtype& source ) const;
		bool											ReleaseResources						();
		bool											RotateRayToDiffuse						(	const HELIODETIC_POINT& quadlocation,
																									const HELIODETIC_UNITVECTOR& look,
																									nxVector& rotatedray ) const;
		double											DistanceWeight							( double distance ) const;
		bool											ComputeMultipliersAndAdjustScatterArray	( size_t pointdx );
		bool											CreateDiffuseIndexesForRay				(	const SKTRAN_RayOptical_Base* ray, size_t pointidx,  size_t rayidx );
		const SKTRAN_CoordinateTransform_V2*			CoordinatesPtr							() const { return m_coords.get(); }
		const SKTRAN_RayFactory_Base*					RayFactoryPtr							() const { return m_rayfactory.get();}
		bool											MakeGroundSourceDiffuseIndexesForRay	( const SKTRAN_RayOptical_Base*		ray,
																								  size_t							pointidx, 
																								  size_t							rayidx,
																								  SKTRAN_HR_Diffuse_Index_Array*	indexarray_lo); 
	protected:
		bool											AltWeightsForProfile		( double alt, size_t profileidx, SKTRAN_HR_WEIGHT_TYPE* altweights, size_t* altindex, size_t* numindex ) const;
		const std::vector<size_t>&						ProfileStartIndex			() const { return m_profilestartidx;}
		size_t											GroundStartIndex			() const { return m_groundstartidx;}
		const std::vector<SKTRAN_HR_Diffuse_Point>&		DiffusePoints				() const { return m_diffusepoints;}

	public:
														SKTRAN_HR_Diffuse_Table_CPU();
		virtual										   ~SKTRAN_HR_Diffuse_Table_CPU();

		bool											Initialize					( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>&     coords,
																					  const SKTRAN_OpticalPropertiesIntegrator_Base&			integrator,
																					  const SKTRAN_SourceTermIntegrator_Base&					srcintegrator,
																					  std::shared_ptr<const SKTRAN_RayFactory_Base>				rayfactory,
																					  const std::vector<SKTRAN_Source_Term*>&					sources,
																					  const SKTRAN_TableOpticalProperties_Base&					opticaltable );

		bool											ConfigureStorage			( std::unique_ptr< RadStore_Base >& radStorage, std::unique_ptr< Avals_Base >&   avals );
		bool											DiffuseSource				( const SKTRAN_SourceTermQueryObject_Base& qobj, double&            source ) const; 
		bool											DiffuseSource               ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC&  losvec ) const ;
		bool											GroundSource				( const SKTRAN_SourceTermQueryObject_Base& qobj, double&			source ) const ;
        bool											GroundSource                ( const SKTRAN_SourceTermQueryObject_Base& qobj, SKTRAN_Stokes_NC&	losvec ) const; 
		bool											DeclareFirstOrderInitialized ( ) ;
		bool											DeclareAllScattered		     ( ) ;
		bool											DeclareAllIntegrated         ( ) ; 
		size_t											GetDiffuseSourceOrder        ( ) const ;
		void											SetGroundStartIdx			( size_t idx )				{ m_groundstartidx = idx; }
		bool											SetProfileStartIdx			( const std::vector<size_t>& profilestartidx );
		bool											AllocateDiffusePoints		( size_t numdiffuse, size_t numground );
		size_t											NumDiffusePoints			() const					{ return m_groundstartidx; } //!<Returns the number of non-ground diffuse points
		size_t											NumGroundPoints				() const					{ return m_diffusepoints.size() - m_groundstartidx; } //!< The number of ground diffuse points only
		const SKTRAN_HR_Diffuse_Point&					DiffusePointAt				( size_t idx ) const		{ return m_diffusepoints[idx]; }
		const SKTRAN_HR_Diffuse_Point&					GroundPointAt               ( size_t idx ) const		{ return m_diffusepoints[idx+m_groundstartidx];}
		SKTRAN_HR_Diffuse_Point&						DiffusePointAtVar			( size_t idx )				{ return m_diffusepoints[idx]; }
		bool											ConfigureIndexes			();
		bool											CleanDiffuseIndexes			();
		double											IntegrateScalarIncoming		( const SKTRAN_SourceTermQueryObject_Base&, std::function<double(const nxVector&, const nxVector&)> func ) const;
        SKTRAN_Stokes_NC                                IntegrateVectorIncoming     ( const SKTRAN_SourceTermQueryObject_Base&, std::function<SKTRAN_ScatMat_MIMSNC(const nxVector&, const nxVector&)> func ) const;
		bool											PreSetup					();
		bool											CalcFirstOrderIncomingPoint	( size_t pointidx,  SKTRAN_RayOptical_Base* ray = NULL );
		bool 											CalcFirstOrderIncomingRay   ( size_t pointidx, size_t rayidx, SKTRAN_RayOptical_Base* ray = NULL);
		bool											CalcScatteringMatrixPoint	( size_t pointidx ) ;
		bool											PostProcessPoint			( size_t pointidx ) { return true; }
		bool											ScatterPoint		     	( size_t pointidx ); 
		bool											ComputeNextOrderPoint		( size_t pointidx ); 
		void											SetDiffuseLocationSphere    ( std::unique_ptr<const SKTRAN_UnitSphere_V2> sphere ) { m_diffuselocationsphere = std::move( sphere ); }
		void											PrintMemReport				() const;
		void											SetDiagnostic				( std::vector<size_t> idxs, std::vector<size_t> scatords ) { m_profindices = idxs; m_scatords = scatords; }
		bool											IsSetDiagnostic				( size_t order ) {return ( std::find(m_scatords.begin(), m_scatords.end(), order) != m_scatords.end() ); }
		const std::vector<size_t>&						DiagnosticProfIdxs			() const { return m_profindices; }
		size_t											DiagnosticProfAlts			() const { return m_diffusealts[0].size(); }

	public:
		virtual void									SetNumProfileInterp			( size_t numinterp )		{ m_numprofileinterp = numinterp; }

		virtual bool									DumpIncomingRadiances		( size_t order, double wlen );
		
		virtual bool									DumpOutGoingRadiances		( size_t order, double wlen );
		
		virtual bool									ChooseDiffusePoints			( const HELIODETIC_POINT&				pt, 
																					  size_t*								diffuseindex, 
																					  SKTRAN_HR_WEIGHT_TYPE*				diffuseweights, 
																					  size_t*								numpoints ) const;

		virtual bool									ChooseGroundPoints			( const HELIODETIC_POINT&				pt,
																					  size_t*								diffuseindex,
																					  SKTRAN_HR_WEIGHT_TYPE*				diffuseweights,
																					  size_t&								numpoints ) const;

};


