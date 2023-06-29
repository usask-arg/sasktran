#pragma once



/*-----------------------------------------------------------------------------
 *					RadStore_Base		 2016- 12- 4*/
/** Store and manipulate radiances in either the scalar, pseudo-vector (order 1) 
 *	or pseudo-vector (order > 1) paradigm. 
**/
/*---------------------------------------------------------------------------*/

class RadStore_Base
{
	private:
		size_t										m_numincomingrays;
		size_t										m_numoutgoingrays;

	protected:
		size_t										m_order_incoming;
		size_t										m_order_outgoing;
		const SKTRAN_TableOpticalProperties_Base*	m_opticaltable;

	private:
		void										ReleaseResources			( );

	public:
													RadStore_Base               ( );
		virtual									   ~RadStore_Base               ( );
		size_t										GetNumIncoming              ( ) const { return m_numincomingrays;}
		size_t										GetNumOutgoing				( ) const { return m_numoutgoingrays;}

		virtual bool								Initialize					(	std::shared_ptr<const SKTRAN_CoordinateTransform_V2>&       coords,
																					const SKTRAN_OpticalPropertiesIntegrator_Base&				optintegrator,
																					const SKTRAN_SourceTermIntegrator_Base&						srcintegrator,
																					std::shared_ptr<const SKTRAN_RayFactory_Base>				rayfactory,
																					const std::vector<SKTRAN_Source_Term*>&						sources,
																					const SKTRAN_TableOpticalProperties_Base&					opticaltable );


	public:
		virtual bool								AllocateStorage             ( size_t incomingcounter, size_t outgoingcounter );
		virtual bool								CleanDiffuseIndexes         ( );
		virtual bool								DumpIncomingRadiances       ( const SKTRAN_HR_Diffuse_Table_CPU* table, size_t order, double wlen ) = 0;
		virtual bool								DumpOutgoingRadiances       ( const SKTRAN_HR_Diffuse_Table_CPU* table, size_t order, double wlen ) = 0;
		virtual double								IncomingScalarRadiance		( size_t idx ) const = 0;
        virtual SKTRAN_Stokes_NC                    IncomingVectorRadiance      ( size_t idx ) const = 0;
		virtual bool                                DeclareFirstOrderInitialized( ) = 0;
		virtual bool                                DeclareAllScattered		    ( ) = 0;
		virtual bool                                DeclareAllIntegrated        ( ) = 0; 
		virtual size_t								GetDiffuseSourceOrder       ( ) const = 0;
		virtual bool								StoreFirstOrderIncoming     ( const SKTRAN_SourceTermIntegrator_Base* srcintegrator, const std::vector<SKTRAN_Source_Term*>* sources, const SKTRAN_RayOptical_Base* ray, size_t index ) = 0;
		virtual bool								ComputeNextOrderPoint       ( const SKTRAN_HR_Diffuse_Point& point ) = 0;
		virtual bool								ScatterPoint                ( const SKTRAN_HR_Diffuse_Point& point, const Avals_Base& Avals ) = 0;
		virtual bool								DiffuseSource               ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& look, double&           source ) const = 0;
		virtual bool								DiffuseSource               ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& look, SKTRAN_Stokes_NC& source ) const = 0;
		virtual bool								GroundSource                ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& look, double&           source ) const = 0;
		virtual bool								GroundSource                ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& look, SKTRAN_Stokes_NC& source ) const = 0;
};



/*-----------------------------------------------------------------------------
 *					RadStore_Scalar		 2016- 12- 4*/
/** **/
/*---------------------------------------------------------------------------*/

class RadStore_Scalar : public RadStore_Base
{
	protected:
		std::vector<SKTRAN_HR_WEIGHT_TYPE>      m_totalincomingscalar;
		std::vector<SKTRAN_HR_WEIGHT_TYPE>      m_incomingscalar;
		std::vector<SKTRAN_HR_WEIGHT_TYPE>      m_outgoingscalar;
		std::vector<double>						m_totaloutgoingscalar;

	protected:
		bool									ScatterPoint_Scalar							( const SKTRAN_HR_Diffuse_Point& point, const Avals_Base& Avals );
		virtual double							DumpIncomingRadiances_GetIncomingRadiance	( size_t incomingRadIdx ) const { return m_incomingscalar[incomingRadIdx];}

	public:
												RadStore_Scalar             ( );
		virtual								   ~RadStore_Scalar             ( ) {;}
		virtual bool							AllocateStorage             ( size_t incomingcounter, size_t outgoingcounter ) override;
		virtual bool							CleanDiffuseIndexes         ( ) override;
		virtual bool							DeclareFirstOrderInitialized( ) override;
		virtual bool							DeclareAllScattered			( ) override;
		virtual bool							DeclareAllIntegrated        ( ) override; 
		virtual size_t							GetDiffuseSourceOrder		( ) const override;
		virtual bool							DumpIncomingRadiances       ( const SKTRAN_HR_Diffuse_Table_CPU* table, size_t order, double wlen ) override;
		virtual bool							DumpOutgoingRadiances       ( const SKTRAN_HR_Diffuse_Table_CPU* table, size_t order, double wlen ) override;
		virtual double							IncomingScalarRadiance		( size_t idx ) const override { return m_totalincomingscalar[idx]; }
        virtual SKTRAN_Stokes_NC                IncomingVectorRadiance      ( size_t idx ) const override { SKTRAN_Stokes_NC ret; ret.SetTo(m_totalincomingscalar[idx], 0, 0, 0); return ret; }
		virtual bool							StoreFirstOrderIncoming     ( const SKTRAN_SourceTermIntegrator_Base* srcintegrator, const std::vector<SKTRAN_Source_Term*>* sources, const SKTRAN_RayOptical_Base* ray, size_t index ) override;
		virtual bool							ComputeNextOrderPoint       ( const SKTRAN_HR_Diffuse_Point& point ) override;
		virtual bool							ScatterPoint                ( const SKTRAN_HR_Diffuse_Point& point, const Avals_Base& Avals ) override;
		virtual bool							DiffuseSource               ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& look, double&           source ) const override;
		virtual bool							DiffuseSource               ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& look, SKTRAN_Stokes_NC& source ) const override;
		virtual bool							GroundSource                ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& look, double&           source ) const override;
		virtual bool							GroundSource                ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& look, SKTRAN_Stokes_NC& source ) const override;

};


/*-----------------------------------------------------------------------------
 *					RadStore_PV1		 2016- 12- 4*/
/** **/
/*---------------------------------------------------------------------------*/

class RadStore_PV1 : public RadStore_Scalar
{	
	protected:
		virtual bool							GetPolarizedComponent_basisImplicit ( size_t radIndex, const nxVector& look_transformed, SKTRAN_Stokes_NC& incomingvec, const HELIODETIC_POINT& queryPoint, const HELIODETIC_UNITVECTOR& inrayQCoords ) const;

	public:
		virtual								   ~RadStore_PV1				( ) {;}
		virtual size_t							GetDiffuseSourceOrder		( ) const override;
		virtual bool							DiffuseSource               ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& look, double&           source ) const override;
		virtual bool							DiffuseSource               ( const SKTRAN_SourceTermQueryObject_Base& qNoTransform, const SKTRAN_HR_Diffuse_Point& point, const nxVector& look, SKTRAN_Stokes_NC& source ) const override;

};

/*-----------------------------------------------------------------------------
 *					RadStore_Polarized		 2016- 12- 4*/
/** **/
/*---------------------------------------------------------------------------*/

class RadStore_Polarized : public RadStore_PV1
{
	protected:
		struct linPolComponents
		{ 
			double Q, U;
			linPolComponents ( ) {;}
			linPolComponents (double q, double u) : Q(q), U(u) {;}
		};

	protected:
		std::vector<SKTRAN_Stokes_NC>            m_totalincomingvector;
		std::vector<SKTRAN_Stokes_NC>            m_incomingvector;
		std::vector<linPolComponents>            m_outgoingLinPols; // This could be moved to a PV3 class, but the memory hit is small
		std::vector<linPolComponents>			 m_outgoingLinPols_lastOrder;

		size_t                                   m_maxPolarizedOrder_integration;
		size_t                                   m_maxPolarizedOrder_scattering;
		size_t 									 m_numDiffusePolarizedOrder;
		bool 									 m_modec;
		double 									 m_adaptivefraction;

	protected:		
		bool									ScatterPoint_Participating					( const SKTRAN_HR_Diffuse_Point& point, const Avals_Base& Avals );
		bool									ScatterPoint_Boundary						( const SKTRAN_HR_Diffuse_Point& point, const Avals_Base& Avals );
		//virtual void							ScatterPoint_PolarizedToScalar_StoreData    ( size_t outidx, const skRTStokesVector& value );
		void									ScatterPoint_StoreData						( size_t outidx, const SKTRAN_Stokes_NC& value );
		bool									ComputeNextOrderPoint_IntegrateScalars		( const SKTRAN_HR_Diffuse_Point& point );
		bool									ComputeNextOrderPoint_IntegrateVectors		( const SKTRAN_HR_Diffuse_Point& point );
		virtual bool							GetPolarizedComponent_basisImplicit			( size_t radIndex, const nxVector& look_transformed, SKTRAN_Stokes_NC& incomingvec, const HELIODETIC_POINT& queryPoint, const HELIODETIC_UNITVECTOR& inrayQCoords ) const override;
		virtual double							DumpIncomingRadiances_GetIncomingRadiance	( size_t incomingRadIdx ) const { return m_totalincomingvector[incomingRadIdx].I();}

	public:
												RadStore_Polarized							( );
		virtual								   ~RadStore_Polarized							( ) {;}
		void									ConfigurePolarizedScattering				( size_t numPolarizedDiffusions, bool mode_c, double adaptiveFraction=0.0);
		virtual	bool							AllocateStorage								( size_t incomingcounter, size_t outgoingcounter ) override;
		virtual bool							CleanDiffuseIndexes							( ) override;
		virtual bool							DumpIncomingRadiances						( const SKTRAN_HR_Diffuse_Table_CPU* table, size_t order, double wlen ) override;
		virtual bool							StoreFirstOrderIncoming						( const SKTRAN_SourceTermIntegrator_Base* srcintegrator, const std::vector<SKTRAN_Source_Term*>* sources, const SKTRAN_RayOptical_Base* ray, size_t index ) override;
		virtual bool							ScatterPoint								( const SKTRAN_HR_Diffuse_Point& point, const Avals_Base& Avals ) override;
		virtual bool							ComputeNextOrderPoint						( const SKTRAN_HR_Diffuse_Point& point ) override;

		// Necessary for the weighting function calculation
		virtual double							IncomingScalarRadiance(size_t idx) const override { return m_totalincomingvector[idx].I(); }
        virtual SKTRAN_Stokes_NC                IncomingVectorRadiance(size_t idx) const override { SKTRAN_Stokes_NC ret; ret.SetTo(m_totalincomingvector[idx].I(),0,0,0); return ret; }

		virtual bool							DeclareAllScattered			( ) override;

};

