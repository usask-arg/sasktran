#include <unordered_map>

/*-----------------------------------------------------------------------------
 *					skBaumIceCrystals_DataBase		2014-4-14*/
/** **/
/*---------------------------------------------------------------------------*/

class skBaumIceCrystals_DataBase
{
	public:
		enum ICE_CRYSTAL_SHAPE {		BAUM_AGGREGATE_SOLID_COLUMNS_SEVERELYROUGH,
										BAUM_GENERAL_HABIT_MIXTURE_SEVERELYROUGH,
										BAUM_SOLID_COLUMNS_SEVERELYROUGH,
										BAUM_UNDEFINED_ICECRYSTAL
							   };

	private:
		struct CurrentIndex												// A small helper class to help with the linear interpolation
		{																// of the wavelen, De and phase angle arrays.
			private:
				const nx1dArray<double>&	m_valuearray;
				size_t						m_index_0;
				size_t						m_index_1;
				double						m_value_0;
				double						m_value_1;
				double						m_currentvalue;
			public:
											CurrentIndex	( const nx1dArray<double>&	values);						//!< The wavelength in nanometers
				bool						UpdateIndices	( double v);
				size_t						Idx0			() const { return m_index_0;}
				size_t						Idx1			() const { return m_index_1;}
				double						V0				() const { return m_value_0;}
				double						V1				() const { return m_value_1;}
				double						V				() const { return m_currentvalue;}
		};

	private:
		double							m_thetaCutoff;					//!< The cutoff angle in degrees for determining the forward scatter contribution in p11
		bool							m_thetacutoff_auto;				//!< IF true then use an automatic algorithm to set the Theta Cutoff angle
		bool							m_legendrecached;				//!< If true then the file contains cached legendre moments
		CurrentIndex					m_Wavelenindex;
		CurrentIndex					m_Deindex;
		CurrentIndex					m_PhaseAngleindex;
		nxRegistryConfiguration			m_config;						//!< The object used to access the registry.
		nxString						m_basedir;						//!< String to hold the directory of the BAUM netcdf files
		bool							m_isdirty;						//!< True if the local variables are not loaded or out of sync with requested ice crystal shape.
		ICE_CRYSTAL_SHAPE				m_icecrystalshape;				//!< The requested ic crystal database. BAUM 2014 offers three ice-crystal shape databases			.
		nx1dArray<double>				m_wavelen;						//!< Array( wavelength ).	The wavelength in nanometers
		nx1dArray<double>				m_De;							//!< Array( De ).			The effective diameter in microns
		nx1dArray<double>				m_phaseangles;					//!< Array(angles).			The scattering angles in degrees (0 to 180.0) for the phase matrix.
		nx1dArray<double>				m_phasemoments;					//!< Array( moments ).      The phase moments (0, 1, 2, 3 ...) for the phase matrix
		nx2dArray<double>				m_ext_crosssection;				//!< Array(wavelength, De). The extinction cross-section, cm2
		nx2dArray<double>				m_abs_crosssection;				//!< Array(wavelength, De).	The absorption cross-section, cm2, indexed by wavelength and De
		nx2dArray<double>				m_scatt_crosssection;			//!< Array(wavelength, De).	The scattering cross-section, cm2, indexed 
		nx2dArray<double>				m_deltafraction;				//!< Array(wavelength, De).	The forward scattering component of the phase matrix
		nx3dArray<double>				m_p11;							//!< Array(wavelength, De, numagles). Phase matrix element p11 (adjusted for forward scatter component)
		nx3dArray<double>				m_p21;							//!< Array(wavelength, De, numagles). Phase matrix element p21 and p12.
		nx3dArray<double>				m_p22;							//!< Array(wavelength, De, numagles). Phase matrix element p22
		nx3dArray<double>				m_p33;							//!< Array(wavelength, De, numagles). Phase matrix element p33
		nx3dArray<double>				m_p43;							//!< Array(wavelength, De, numagles). Phase matrix element p43 and p34
		nx3dArray<double>				m_p44;							//!< Array(wavelength, De, numagles). Phase matrix element p44
		nx3dArray<double>				m_p11_legendre;					//!< Array(wavelength, De, nummoment).  Phase moments for P11
	
	private:
		bool							FetchFilename					( nxString* filename );
		bool							FetchBaseDirectory				( );
		bool							InterpolateCrossSectionArray    ( const nx2dArray<double>& xs, double* value );
		bool							InterpolatePhaseArrayAtIndex    ( const nx3dArray<double>& xs, size_t phaseidx, double* value );
		bool							InterpolateP11					( double w, double de, double angle, nx3dArray<double>& p, double* p11);
		bool							TruncateAndComputeDeltaFraction	( nx3dArray<double>& p11, nx2dArray<double>* deltafraction );
		bool							TruncateAndComputeDeltaFraction_NoTruncation	( nx3dArray<double>& p11, nx2dArray<double>* deltafraction );
		double							ForwardScatterCutoffAngle		( double de );

	public:
										skBaumIceCrystals_DataBase		();
									   ~skBaumIceCrystals_DataBase		();
		bool							SetForwardScatterCutoffAngle	( double cutoff_degrees );
		bool							LoadDatabase					( bool useDeltaEddingtonApproximation );
		bool							IsLoaded						() const {return m_wavelen.size() > 0;};
		bool							IsLegendreCached                () const { return m_legendrecached; };
		bool							InterpolateCrossSections		( double w, double de, double* absxs, double* extxs, double* scattxs);
		bool							InterpolatePhaseMatrix			( double w, double de, double angle, skRTPhaseMatrix* phasematrix);
		bool							InterpolatePhaseScalar          ( double w, double de, double angle, double& p11 );
		bool							InterpolateForwardScatter		( double w, double de, double* forwardscatter);
		bool							InterpolateLegendre             ( double w, double de, std::vector<double>& legendremoments );

};



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_BaumIceCrystals2014		2014-4-15*/
/** **/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_BaumIceCrystals2014 : public skOpticalProperties
{
	private:
		skClimatology*						m_backgroundatmosphere;
		skBaumIceCrystals_DataBase			m_icecrystalsdb;
		skClimatology*						m_effectivesize;
		bool                                m_useEddingtonApproximation;

		double								m_current_wavenumber;
		double								m_current_De;
		double								m_current_absxs;
		double								m_current_extxs;
		double								m_current_scattxs;
		double								m_current_forwardscatter;

		const std::vector<double>*			m_phasegridhint;
		std::vector<double>					m_p11cache;
		bool								m_p11cachevalid;

	private:
		void								ResetCurrentValues							( double de);

	public:
											skOpticalProperties_BaumIceCrystals2014		();
		virtual							   ~skOpticalProperties_BaumIceCrystals2014		() override;
		bool								SetEffectiveSizeClimatology					( skClimatology* effectivesize);
		bool								SetForwardScatterCutoffAngle				( double cutoff_degrees );
		void                                SetUseDeltaEddingtonApproximation           ( bool b ) { m_useEddingtonApproximation = b; }

	public:
		virtual bool						SetAtmosphericState					( skClimatology* neutralatmosphere) override;
		virtual bool						SetLocation							( const GEODETIC_INSTANT& pt, bool* crosssectionschanged ) override;
		virtual bool						InternalClimatology_UpdateCache			( const GEODETIC_INSTANT& pt) override;
		virtual bool						CalculateCrossSections				( double wavenumber, double* absxs, double* extxs, double* scattxs) override;										//!< Calculate cross-sections at the specified wave-number.
		virtual bool						CalculatePhaseMatrix				( double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix) override;	//!< Calculate the phase matrix at the specified wavenumber and scattering angle
		virtual bool						CalculateP11						( double wavenumber, std::pair<double, size_t> cosscatterandindex, double& p11 ) override;	//!< Calculate the phase matrix at the specified wavenumber and scattering angle
		virtual bool						IsScatterer							() const override;			//!< Returns true if this particles scatters radiation
		virtual bool						IsAbsorber							() const override;			//!< Returns true if this particles absorbs radiation radiation
		virtual double						DeltaFunctionForwardScatterFraction () const override;
		virtual bool						IsHeightDependent() const { return false; }
		virtual bool						PhaseGridHint(const std::vector<double>& cosscatterangles) { m_phasegridhint = &cosscatterangles; return true; }
		virtual bool						LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff) override;
};

