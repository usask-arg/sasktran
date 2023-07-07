#include "sktran_montecarlo_internals.h" 

class SKTRAN_TableOpticalProperties_MCBase;

/*
 * Right now this is a class. It would be better if it were
 * a union so we wouldn't have to align the next member, but then
 * we'd have to make the stokes vector's constructor trivial, which 
 * it basically is anyway...
 */
class SKTRAN_MCPhoton_RadInfo
{
	SKTRAN_Stokes_NC  m_vec;
	SKTRAN_Stokes_NC  m_debugVec;
	SKTRAN_Stokes_NC  m_recentContrib;

	//SKTRAN_Stokes_NC  m_vec2;
	//SKTRAN_Stokes_NC  m_debugVec2;
	//SKTRAN_Stokes_NC  m_recentContrib2;
	
public:
	SKTRAN_MCPhoton_RadInfo ( );

	void SetVector          ( SKTRAN_StokesScalar          s ) { m_vec.SetTo(s,0.0,0.0,0.0); m_recentContrib.SetTo(s,0.0,0.0,0.0);}
	void AddToVector        ( SKTRAN_StokesScalar          s ) { m_vec.Add_I(s);             m_recentContrib.SetTo(s,0.0,0.0,0.0);}
	void SetVector          ( const SKTRAN_Stokes_NC&    v ) { m_vec = v;                  m_recentContrib = v;}
	void AddToVector        ( const SKTRAN_Stokes_NC&    v ) { m_vec+= v;                  m_recentContrib = v;}
	void SetDebugVector     ( const SKTRAN_Stokes_NC& v ) { m_debugVec = v;}
	void ScaleVector		( double s ) { m_vec.SetTo(m_vec.I() * s, m_vec.Q() * s, m_vec.U() * s, m_vec.V() * s); }

	      SKTRAN_StokesScalar  GetScalar           ( ) const {return m_vec.I();}
	const SKTRAN_Stokes_NC&  GetVector           ( ) const {return m_vec;}
	const SKTRAN_Stokes_NC&  GetDebugVector      ( ) const {return m_debugVec;}
	      SKTRAN_StokesScalar  GetRecentContribSca ( ) const { return m_recentContrib.I();}
	const SKTRAN_Stokes_NC&  GetRecentContribVec ( ) const { return m_recentContrib;}

	//void SetVector          ( SKTRAN_StokesScalar          s, bool secondary ) { if (secondary) { m_vec2.SetTo(s,0.0,0.0,0.0); m_recentContrib2.SetTo(s,0.0,0.0,0.0); } else SetVector(s); }
	//void AddToVector        ( SKTRAN_StokesScalar          s, bool secondary ) { if (secondary) { m_vec2.Add_I(s); m_recentContrib2.SetTo(s,0.0,0.0,0.0); } else AddToVector(s); }
	//void SetVector          ( const SKTRAN_Stokes_NC&    v, bool secondary) { if (secondary) { m_vec2 = v; m_recentContrib2 = v; } else SetVector(v); }
	//void AddToVector        ( const SKTRAN_Stokes_NC&    v, bool secondary) { if (secondary) { m_vec2+= v; m_recentContrib2 = v; } else AddToVector(v); }
	//void SetDebugVector     ( const SKTRAN_Stokes_NC& v, bool secondary) { if (secondary) { m_debugVec2 = v; } else SetDebugVector(v); }
	//void ScaleVector		( double s, bool secondary) { if (secondary) { m_vec2.SetTo(m_vec2.I() * s, m_vec2.Q() * s, m_vec2.U() * s, m_vec2.V() * s); } else ScaleVector(s); }

	//      SKTRAN_StokesScalar  GetScalar           ( bool secondary ) const { return secondary ? m_vec2.I() : m_vec.I(); }
	//const SKTRAN_Stokes_NC&  GetVector           ( bool secondary ) const { return secondary ? m_vec2 : m_vec; }
	//const SKTRAN_Stokes_NC&  GetDebugVector      ( bool secondary ) const { return secondary ? m_debugVec2 : m_debugVec; }
	//      SKTRAN_StokesScalar  GetRecentContribSca ( bool secondary ) const { return secondary ? m_recentContrib2.I() : m_recentContrib.I(); }
	//const SKTRAN_Stokes_NC&  GetRecentContribVec ( bool secondary ) const { return secondary ? m_recentContrib2 : m_recentContrib; }

}; 

//class SKTRAN_MCBasis{
//	public: 
//		HELIODETIC_UNITVECTOR x;
//		HELIODETIC_UNITVECTOR y;
//		HELIODETIC_UNITVECTOR z;
//
//		SKTRAN_MCBasis( ){ }
//		SKTRAN_MCBasis(const SKTRAN_MCBasis& other)
//		{
//			x = other.x;
//			y = other.y;
//			z = other.z;
//		}
//		~SKTRAN_MCBasis( ){ }
//};

typedef HELIODETIC_BASIS SKTRAN_MCBasis;

class SKTRAN_MCPhoton_Base
{
	protected:
		std::unique_ptr<SKTRAN_RayOptical_Base>		m_photonOptical;
		std::vector<SKTRAN_MCPhoton_RadInfo>		m_photonRadiance;
		std::vector<SKTRAN_MCPhoton_RadInfo>		m_photonSource;
		
		std::vector<double> m_finalWavelength;
		std::vector<double> m_currentWavelength;
		std::vector<double> m_scatterWeight;		// 
		std::vector<double> m_weightFactor;		//
		std::vector<double> m_albedo;
		std::vector<double> m_transmission;
		std::vector<double> m_scatterFactor;
		std::vector<std::vector<double>> m_opticaldepth;

		double m_optFactor;

	public:
        SKTRAN_MCBasis                  m_basis;           // {look, theta, phi} as in Tony's thesis
		bool                            m_isGroundScatter;
		//double                          m_scatterWeight;
		//double                          m_weightFactor;		// #m_scatterWeight = m_weightFactor*prevScatterWeight
		std::vector<double>				m_solarSlantColumns; // integrated number density OR length (within each atmospheric cell) of direct solar ray
		std::vector<double>				m_scatterSlantColumns; // integrated number density OR length (within each atmospheric cell) of the ray history so far, excluding direct solar rays
		//double                          m_albedo;
		double                          m_distanceProb;
		double                          m_targetTau;
		HELIODETIC_VECTOR               m_scatterVector;

		size_t							m_numWavelengths;
		size_t							m_primaryWavelengthIndex;
		double							m_primaryWavelength;

		double							m_sourcesRadiance;
		double							m_sourcesRadiance2;

		double							m_randNum;
		bool							m_manualScatter;
		bool							m_elasticScatter;

	private:
		SKTRAN_PhaseMat_MIMSNC                 m_toObserverOp;

	private:
		void                            Initialize          ( );
//										SKTRAN_MCPhoton		( const SKTRAN_MCPhoton& other );					// Not defined, not allowed
										SKTRAN_MCPhoton_Base		( std::unique_ptr<SKTRAN_RayOptical_Base> r );		// Defined but deprecated.

	protected: // helper functions for CalculateTransmissionsAtmoScatter and CalculateTransmissionsGroundScatter
		bool							FindScatterPointCellIndex( const SKTRAN_RayStorage_Base* storage, const std::vector<double>& odarray, const double & targetTau, size_t& scatterCellIndex, HELIODETIC_POINT& scatterCellStartPoint ) const;
		bool							ConfigureQuadratureCoefficients( const SKTRAN_RayStorage_Base* storage, const HELIODETIC_POINT& scatterPoint, const size_t& scatterCellIndex, const HELIODETIC_POINT& scatterCellStartPoint, SKTRAN_OpticalDepthCalculator_LinearWithHeight& odcalculator ) const;

	public:
		                                SKTRAN_MCPhoton_Base   ( ); 
		virtual                        ~SKTRAN_MCPhoton_Base   ( );
	//	                                SKTRAN_MCPhoton		( SKTRAN_MCPhoton&& other );	// The 
										SKTRAN_MCPhoton_Base		( const SKTRAN_MCPhoton_Base& other );		// The copy constructor is really a move operator
	//									SKTRAN_MCPhoton		( const SKTRAN_MCPhoton& other );		// Copy constructor used by std::vector::resize under g++
		SKTRAN_MCPhoton_Base&				operator =			( SKTRAN_MCPhoton_Base& other );		// The assignment operator is really a move operator, its reqiuired for vector storage
		virtual SKTRAN_MCPhoton_Base*	Clone() const = 0;

		virtual bool					SimultaneousRaman() const = 0;
									   
		virtual SKTRAN_RayOptical_Base*         photonOptical     ( );
		virtual const SKTRAN_RayOptical_Base*   photonOptical     ( ) const;
		virtual SKTRAN_MCPhoton_RadInfo&  photonRadiance    ( );
		virtual SKTRAN_MCPhoton_RadInfo&  photonRadiance    ( bool elasticRaman ) { return photonRadiance(); }
		const virtual SKTRAN_MCPhoton_RadInfo&  photonRadiance    ( ) const;
		const virtual SKTRAN_MCPhoton_RadInfo&  photonRadiance    ( bool elasticRaman ) const { return photonRadiance(); }
		virtual std::vector<SKTRAN_MCPhoton_RadInfo>&  photonRadiances   ( );
		virtual std::vector<SKTRAN_MCPhoton_RadInfo>&  photonRadiances   ( bool elasticRaman ) { return photonRadiances(); }
		virtual const std::vector<SKTRAN_MCPhoton_RadInfo>&  photonRadiances   ( ) const;
		virtual const std::vector<SKTRAN_MCPhoton_RadInfo>&  photonRadiances   ( bool elasticRaman ) const { return photonRadiances(); }
		virtual SKTRAN_MCPhoton_RadInfo&  photonSource( );
		virtual SKTRAN_MCPhoton_RadInfo&  photonSource( bool elasticRaman ) { return photonSource(); }
		virtual const SKTRAN_MCPhoton_RadInfo&  photonSource( ) const;
		virtual const SKTRAN_MCPhoton_RadInfo&  photonSource( bool elasticRaman ) const { return photonSource(); }
		virtual std::vector<SKTRAN_MCPhoton_RadInfo>&  photonSources( );
		virtual std::vector<SKTRAN_MCPhoton_RadInfo>&  photonSources( bool elasticRaman ) { return photonSources(); }
		virtual const std::vector<SKTRAN_MCPhoton_RadInfo>&  photonSources( ) const;
		virtual const std::vector<SKTRAN_MCPhoton_RadInfo>&  photonSources( bool elasticRaman ) const { return photonSources(); }


		virtual void					ResetRadiance     ( );
		virtual void					ResetFactors	  ( );
		virtual bool					CalculateAlbedo	  ( const SKTRAN_TableOpticalProperties_Base* opticalpropertiestable, const SKTRAN_TableOpticalProperties_MCBase* mcOptTable, const HELIODETIC_POINT & scatterPoint );
		virtual bool					UpdateScatterWeight( double forcedScatterCorrFactor ) = 0;
		virtual bool					PreGenerateRandNum ( SKTRAN_RNG& rng ) = 0;

		bool                            SetOpticalRay     ( std::unique_ptr<SKTRAN_RayOptical_Base>  r ); // #controlLifetimeManagement true if mcphoton has only reference to r
		bool                            DefineRayBasis    ( );
		const SKTRAN_MCBasis&           GetBasis          ( ) const;
		      SKTRAN_MCBasis&           GetBasisVar       ( );
	    const SKTRAN_PhaseMat_MIMSNC&          GetPathOp         ( ) const;
		      SKTRAN_PhaseMat_MIMSNC&          GetPathOp         ( );
			  
        void                            AddPhaseOpInPath  ( const SKTRAN_ScatMat_MIMSNC& op );
        void                            AddRotateOpInPath ( const SKTRAN_ScatMat_Rot&  op );
		void                            LeftApplyPathTo   ( SKTRAN_Stokes_NC& source ) const;

		virtual bool					SetWavelengths    ( const std::vector<double>& wavelengths ) = 0;
		virtual bool					SetCurrentWavelength( double wavelength ) = 0;
		virtual bool					Configure		  ( const SKTRAN_MCPhoton_Base& other );
		virtual bool					TraceRays		  ( const SKTRAN_OpticalPropertiesIntegrator_Base* integrator, bool curved ) = 0;
		//bool							ResetScatterWeights( );
		//bool							CalculateAlbedos		  ( const HELIODETIC_POINT& scatterPoint, bool isgroundscatter, const SKTRAN_TableOpticalProperties_Base* opttable, const SKTRAN_TableOpticalProperties_MCBase* mcopttable );
		//bool							SetZeroContribution( );

		virtual	bool					CalculateTransmissionsAtmoScatter( SKTRAN_TableOpticalProperties_Base const* opticalpropertiestable, SKTRAN_MCPhoton_Base const* incomingPhoton, const HELIODETIC_POINT& scatterPoint ) = 0;
		virtual bool					CalculateTransmissionsGroundScatter( SKTRAN_MCPhoton_Base const * incomingPhoton, const HELIODETIC_POINT & scatterPoint ) = 0;

		virtual double&							FinalWavelength( ) { return m_finalWavelength[m_primaryWavelengthIndex]; }
		virtual std::vector<double>&			FinalWavelengths( ) { return m_finalWavelength; }
		virtual double&							CurrentWavelength( ) { return m_currentWavelength[m_primaryWavelengthIndex]; }
		virtual std::vector<double>&			CurrentWavelengths( ) { return m_currentWavelength; }
		virtual double&							ScatterWeight() { return m_scatterWeight[m_primaryWavelengthIndex]; }
		virtual double&							ScatterWeight( bool elasticRaman ) { return m_scatterWeight[m_primaryWavelengthIndex]; }
		virtual std::vector<double>&			ScatterWeights() { return m_scatterWeight; }
		virtual std::vector<double>&			ScatterWeights( bool elasticRaman ) { return m_scatterWeight; }
		virtual double&							WeightFactor() { return m_weightFactor[m_primaryWavelengthIndex]; }
		virtual double&							WeightFactor( bool elasticRaman ) { return m_weightFactor[m_primaryWavelengthIndex]; }
		virtual std::vector<double>&			WeightFactors() { return m_weightFactor; }
		virtual std::vector<double>&			WeightFactors( bool elasticRaman ) { return m_weightFactor; }
		virtual double&							Albedo() { return m_albedo[m_primaryWavelengthIndex]; }
		virtual double&							Albedo( bool elasticRaman ) { return m_albedo[m_primaryWavelengthIndex]; }
		virtual std::vector<double>&			Albedos() { return m_albedo; }
		virtual std::vector<double>&			Albedos( bool elasticRaman ) { return m_albedo; }
		virtual double&							Transmission() { return m_transmission[m_primaryWavelengthIndex]; }
		virtual double&							Transmission( bool elasticRaman ) { return m_transmission[m_primaryWavelengthIndex]; }
		virtual std::vector<double>&			Transmissions() { return m_transmission; }
		virtual std::vector<double>&			Transmissions( bool elasticRaman ) { return m_transmission; }
		virtual double&							ScatterFactor() { return m_scatterFactor[m_primaryWavelengthIndex]; }
		virtual std::vector<double>&			ScatterFactors() { return m_scatterFactor; }
		virtual std::vector<double>&			OpticalDepth() { return m_opticaldepth[m_primaryWavelengthIndex]; }
		virtual std::vector<double>&			OpticalDepth( bool elasticRaman ) { return m_opticaldepth[m_primaryWavelengthIndex]; }
		virtual std::vector<std::vector<double>>& OpticalDepths() { return m_opticaldepth; }
		virtual std::vector<std::vector<double>>& OpticalDepths( bool elasticRaman ) { return m_opticaldepth; }
		virtual double&							ManualScatterFactor() { return m_optFactor; }

		virtual const double&					FinalWavelength() const { return m_finalWavelength[m_primaryWavelengthIndex]; }
		virtual const std::vector<double>&		FinalWavelengths() const { return m_finalWavelength; }
		virtual const double&					CurrentWavelength() const { return m_currentWavelength[m_primaryWavelengthIndex]; }
		virtual const std::vector<double>&		CurrentWavelengths() const { return m_currentWavelength; }
		virtual const double&					ScatterWeight() const { return m_scatterWeight[m_primaryWavelengthIndex]; }
		virtual const double&					ScatterWeight( bool elasticRaman ) const { return m_scatterWeight[m_primaryWavelengthIndex]; }
		virtual const std::vector<double>&		ScatterWeights() const { return m_scatterWeight; }
		virtual const std::vector<double>&		ScatterWeights( bool elasticRaman ) const { return m_scatterWeight; }
		virtual const double&					WeightFactor() const { return m_weightFactor[m_primaryWavelengthIndex]; }
		virtual const double&					WeightFactor( bool elasticRaman ) const { return m_weightFactor[m_primaryWavelengthIndex]; }
		virtual const std::vector<double>&		WeightFactors() const { return m_weightFactor; }
		virtual const std::vector<double>&		WeightFactors( bool elasticRaman ) const { return m_weightFactor; }
		virtual const double&					Albedo() const { return m_albedo[m_primaryWavelengthIndex]; }
		virtual const double&					Albedo( bool elasticRaman ) const { return m_albedo[m_primaryWavelengthIndex]; }
		virtual const std::vector<double>&		Albedos() const { return m_albedo; }
		virtual const std::vector<double>&		Albedos( bool elasticRaman ) const { return m_albedo; }
		virtual const double&					Transmission() const { return m_transmission[m_primaryWavelengthIndex]; }
		virtual const double&					Transmission( bool elasticRaman ) const { return m_transmission[m_primaryWavelengthIndex]; }
		virtual const std::vector<double>&		Transmissions() const { return m_transmission; }
		virtual const std::vector<double>&		Transmissions( bool elasticRaman ) const { return m_transmission; }
		virtual const double&					ScatterFactor() const { return m_scatterFactor[m_primaryWavelengthIndex]; }
		virtual const std::vector<double>&		ScatterFactors() const { return m_scatterFactor; }
		virtual const std::vector<double>&		OpticalDepth() const { return m_opticaldepth[m_primaryWavelengthIndex]; }
		virtual const std::vector<double>&		OpticalDepth( bool elasticRaman ) const { return m_opticaldepth[m_primaryWavelengthIndex]; }
		virtual const std::vector<std::vector<double>>&	OpticalDepths() const { return m_opticaldepth; }
		virtual const std::vector<std::vector<double>>&	OpticalDepths( bool elasticRaman ) const { return m_opticaldepth; }
		virtual const double&					ManualScatterFactor() const { return m_optFactor; }
};

class SKTRAN_MCPhoton : public SKTRAN_MCPhoton_Base
{
	public:
		SKTRAN_MCPhoton() {}
		~SKTRAN_MCPhoton() {}

		SKTRAN_MCPhoton* Clone() const { return new SKTRAN_MCPhoton(*this);  }

		virtual bool					SimultaneousRaman() const override { return false; }

		virtual	bool					CalculateTransmissionsAtmoScatter( SKTRAN_TableOpticalProperties_Base const* opticalpropertiestable, SKTRAN_MCPhoton_Base const* incomingPhoton, const HELIODETIC_POINT& scatterPoint ) override { return true; }
		virtual bool					CalculateTransmissionsGroundScatter( SKTRAN_MCPhoton_Base const * incomingPhoton, const HELIODETIC_POINT & scatterPoint ) override { return true; }
		virtual bool					TraceRays(const SKTRAN_OpticalPropertiesIntegrator_Base* integrator, bool curved) override;

		virtual bool					UpdateScatterWeight( double forcedScatterCorrFactor ) override;
		virtual bool					PreGenerateRandNum(SKTRAN_RNG& rng) override { return true; }

		virtual bool					SetWavelengths(const std::vector<double>& wavelengths) override;
		virtual bool					SetCurrentWavelength(double wavelength) override;
};

class SKTRAN_MCPhoton_Simultaneous : public SKTRAN_MCPhoton_Base
{
	public:
		SKTRAN_MCPhoton_Simultaneous() {}
		~SKTRAN_MCPhoton_Simultaneous() {}

		SKTRAN_MCPhoton_Simultaneous* Clone() const { return new SKTRAN_MCPhoton_Simultaneous(*this); }

		virtual bool					SimultaneousRaman() const override { return false; }

		virtual	bool					CalculateTransmissionsAtmoScatter( SKTRAN_TableOpticalProperties_Base const* opticalpropertiestable, SKTRAN_MCPhoton_Base const* incomingPhoton, const HELIODETIC_POINT& scatterPoint );
		virtual bool					CalculateTransmissionsGroundScatter( SKTRAN_MCPhoton_Base const * incomingPhoton, const HELIODETIC_POINT & scatterPoint );
		virtual bool					TraceRays(const SKTRAN_OpticalPropertiesIntegrator_Base* integrator, bool curved) override;

		virtual bool					UpdateScatterWeight(double forcedScatterCorrFactor) override;
		virtual bool					PreGenerateRandNum(SKTRAN_RNG& rng) override { return true; }

		virtual bool					SetWavelengths(const std::vector<double>& wavelengths) override;
		virtual bool					SetCurrentWavelength(double wavelength) override;

};

class SKTRAN_MCPhoton_Inelastic : public SKTRAN_MCPhoton
{
	public:
		SKTRAN_MCPhoton_Inelastic() {}
		~SKTRAN_MCPhoton_Inelastic() {}

		SKTRAN_MCPhoton* Clone() const { return new SKTRAN_MCPhoton_Inelastic(*this); }

		virtual bool					PreGenerateRandNum(SKTRAN_RNG& rng) override { m_randNum = rng(); return true; }

};

class SKTRAN_MCPhoton_SimultaneousInelastic : public SKTRAN_MCPhoton_Simultaneous
{
	public:
		SKTRAN_MCPhoton_SimultaneousInelastic() {}
		~SKTRAN_MCPhoton_SimultaneousInelastic() {}

		SKTRAN_MCPhoton_SimultaneousInelastic* Clone() const { return new SKTRAN_MCPhoton_SimultaneousInelastic(*this); }

		virtual bool					PreGenerateRandNum(SKTRAN_RNG& rng) override { m_randNum = rng(); return true; }
};

class SKTRAN_MCPhoton_Ring : public SKTRAN_MCPhoton_Base
{	
	protected: 
		std::vector<SKTRAN_MCPhoton_RadInfo> m_ePhotonRadiance;
		std::vector<SKTRAN_MCPhoton_RadInfo> m_ePhotonSource;
		std::vector<double> m_eScatterWeight;
		std::vector<double> m_eWeightFactor;
		std::vector<double> m_eAlbedo;
		std::vector<double> m_eTransmission;
		std::vector<std::vector<double>> m_eOpticalDepth;

	public:
		SKTRAN_MCPhoton_Ring() {}
		SKTRAN_MCPhoton_Ring(const SKTRAN_MCPhoton_Ring& other);
		~SKTRAN_MCPhoton_Ring() {}

		SKTRAN_MCPhoton_Ring* Clone() const { return new SKTRAN_MCPhoton_Ring(*this); }
		virtual bool						Configure(const SKTRAN_MCPhoton_Ring& other);

		virtual bool						SimultaneousRaman() const override { return true; }

		virtual void						ResetRadiance() override;
		virtual void						ResetFactors() override;

		virtual	bool						CalculateTransmissionsAtmoScatter( SKTRAN_TableOpticalProperties_Base const* opticalpropertiestable, SKTRAN_MCPhoton_Base const* incomingPhoton, const HELIODETIC_POINT& scatterPoint );
		virtual bool						CalculateTransmissionsGroundScatter( SKTRAN_MCPhoton_Base const * incomingPhoton, const HELIODETIC_POINT & scatterPoint );
		virtual bool						TraceRays( const SKTRAN_OpticalPropertiesIntegrator_Base* integrator, bool curved ) override;
		
		virtual bool						CalculateAlbedo(const SKTRAN_TableOpticalProperties_Base* opticalpropertiestable, const SKTRAN_TableOpticalProperties_MCBase* mcOptTable, const HELIODETIC_POINT & scatterPoint);
		virtual bool						UpdateScatterWeight(double forcedScatterCorrFactor) override;
		virtual bool						PreGenerateRandNum(SKTRAN_RNG& rng) override { m_randNum = rng(); return true; }

		virtual bool						SetWavelengths(const std::vector<double>& wavelengths);
		virtual bool						SetCurrentWavelength(double wavelength) override;

		virtual SKTRAN_MCPhoton_RadInfo&					photonRadiance(bool elasticRaman) override;
		virtual const SKTRAN_MCPhoton_RadInfo&				photonRadiance(bool elasticRaman) const override;
		virtual std::vector<SKTRAN_MCPhoton_RadInfo>&		photonRadiances(bool elasticRaman) override;
		virtual const std::vector<SKTRAN_MCPhoton_RadInfo>& photonRadiances(bool elasticRaman) const override;
		virtual SKTRAN_MCPhoton_RadInfo&					photonSource(bool elasticRaman) override;
		virtual const SKTRAN_MCPhoton_RadInfo&				photonSource(bool elasticRaman) const override;
		virtual std::vector<SKTRAN_MCPhoton_RadInfo>&		photonSources(bool elasticRaman) override;
		virtual const std::vector<SKTRAN_MCPhoton_RadInfo>& photonSources(bool elasticRaman) const override;

		virtual double&										ScatterWeight(bool elasticRaman) override { return elasticRaman ? m_eScatterWeight[m_primaryWavelengthIndex] : SKTRAN_MCPhoton_Base::ScatterWeight(); }
		virtual std::vector<double>&						ScatterWeights(bool elasticRaman) override { return elasticRaman ? m_eScatterWeight : SKTRAN_MCPhoton_Base::ScatterWeights(); }
		virtual double&										WeightFactor(bool elasticRaman) override { return elasticRaman ? m_weightFactor[m_primaryWavelengthIndex] : SKTRAN_MCPhoton_Base::WeightFactor(); }
		virtual std::vector<double>&						WeightFactors(bool elasticRaman) override { return elasticRaman ? m_weightFactor : SKTRAN_MCPhoton_Base::WeightFactors(); }
		virtual double&										Albedo(bool elasticRaman) override { return elasticRaman ? m_albedo[m_primaryWavelengthIndex] : SKTRAN_MCPhoton_Base::Albedo(); }
		virtual std::vector<double>&						Albedos(bool elasticRaman) override { return elasticRaman ? m_albedo : SKTRAN_MCPhoton_Base::Albedos(); }
		virtual double&										Transmission(bool elasticRaman) override { return elasticRaman ? m_transmission[m_primaryWavelengthIndex] : SKTRAN_MCPhoton_Base::Transmission(); }
		virtual std::vector<double>&						Transmissions(bool elasticRaman) override { return elasticRaman ? m_transmission : SKTRAN_MCPhoton_Base::Transmissions(); }
		virtual std::vector<double>&						OpticalDepth(bool elasticRaman) override { return elasticRaman ? m_opticaldepth[m_primaryWavelengthIndex] : SKTRAN_MCPhoton_Base::OpticalDepth(); }
		virtual std::vector<std::vector<double>>&			OpticalDepths(bool elasticRaman) override { return elasticRaman ? m_opticaldepth : SKTRAN_MCPhoton_Base::OpticalDepths(); }

		virtual const double&								ScatterWeight(bool elasticRaman) const override { return elasticRaman ? m_scatterWeight[m_primaryWavelengthIndex] : SKTRAN_MCPhoton_Base::ScatterWeight(); }
		virtual const std::vector<double>&					ScatterWeights(bool elasticRaman) const override { return elasticRaman ? m_scatterWeight : SKTRAN_MCPhoton_Base::ScatterWeights(); }
		virtual const double&								WeightFactor(bool elasticRaman) const override { return elasticRaman ? m_weightFactor[m_primaryWavelengthIndex] : SKTRAN_MCPhoton_Base::WeightFactor(); }
		virtual const std::vector<double>&					WeightFactors(bool elasticRaman) const override { return elasticRaman ? m_weightFactor : SKTRAN_MCPhoton_Base::WeightFactors(); }
		virtual const double&								Albedo(bool elasticRaman) const override { return elasticRaman ? m_albedo[m_primaryWavelengthIndex] : SKTRAN_MCPhoton_Base::Albedo(); }
		virtual const std::vector<double>&					Albedos(bool elasticRaman) const override { return elasticRaman ? m_albedo : SKTRAN_MCPhoton_Base::Albedos(); }
		virtual const double&								Transmission(bool elasticRaman) const override { return elasticRaman ? m_transmission[m_primaryWavelengthIndex] : SKTRAN_MCPhoton_Base::Transmission(); }
		virtual const std::vector<double>&					Transmissions(bool elasticRaman) const override { return elasticRaman ? m_transmission : SKTRAN_MCPhoton_Base::Transmissions(); }
		virtual const std::vector<double>&					OpticalDepth(bool elasticRaman) const override { return elasticRaman ? m_opticaldepth[m_primaryWavelengthIndex] : SKTRAN_MCPhoton_Base::OpticalDepth(); }
		virtual const std::vector<std::vector<double>>&		OpticalDepths(bool elasticRaman) const override { return elasticRaman ? m_opticaldepth : SKTRAN_MCPhoton_Base::OpticalDepths(); }
};

class SKTRAN_MCPhoton_SimultaneousRing : public SKTRAN_MCPhoton_Ring
{
	public:
		SKTRAN_MCPhoton_SimultaneousRing() {}
		~SKTRAN_MCPhoton_SimultaneousRing() {}
		SKTRAN_MCPhoton_SimultaneousRing* Clone() const { return new SKTRAN_MCPhoton_SimultaneousRing(*this); }

		virtual bool					SimultaneousRaman() const override { return true; }


		virtual	bool					CalculateTransmissionsAtmoScatter( SKTRAN_TableOpticalProperties_Base const* opticalpropertiestable, SKTRAN_MCPhoton_Base const* incomingPhoton, const HELIODETIC_POINT& scatterPoint );
		virtual bool					CalculateTransmissionsGroundScatter( SKTRAN_MCPhoton_Base const * incomingPhoton, const HELIODETIC_POINT & scatterPoint );
		virtual bool					TraceRays(const SKTRAN_OpticalPropertiesIntegrator_Base* integrator, bool curved) override;

		virtual bool					CalculateAlbedo(const SKTRAN_TableOpticalProperties_Base* opticalpropertiestable, const SKTRAN_TableOpticalProperties_MCBase* mcOptTable, const HELIODETIC_POINT & scatterPoint);
		virtual bool					UpdateScatterWeight(double forcedScatterCorrFactor) override;
		virtual bool					PreGenerateRandNum(SKTRAN_RNG& rng) override { m_randNum = rng(); return true; }

		virtual bool					SetWavelengths(const std::vector<double>& wavelengths);
		virtual bool					SetCurrentWavelength(double wavelength) override;
};

/* Container for a set of vertices that define a ray (photon path) in R^{3*N} space */
typedef std::vector< std::unique_ptr<SKTRAN_MCPhoton_Base> > SKTRAN_RaySpaceElement;
