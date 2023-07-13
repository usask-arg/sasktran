/*---------------------------------------------------------------------------
 *					class skEmission_Tabulated_HeightWavelength						2003-11-28*/
/**	\ingroup skopticalpropmisc
 *	A class for tabulating isotropic emission as a function of altitude
 *	and wavelength. 
 **/
/*-------------------------------------------------------------------------*/

class skEmission_Tabulated_HeightWavelength: public skEmission
{
	private:
		nx1dArray<double>					m_heights;
		nx1dArray<double>					m_wavelennm;
		nx2dArray<double>					m_emission;
		bool								m_isground;
		size_t								m_idxh1;
		size_t								m_idxh2;
		double								m_h1;
		double								m_h2;

	private:
		void								ReleaseResources			();
		bool								LookupIndicesAndWeights		( const nx1dArray<double>& h, double value, double* w1, size_t* idx1, double* w2, size_t* idx2 ) const;
											skEmission_Tabulated_HeightWavelength( const skEmission_Tabulated_HeightWavelength& other );	// Dont allow copy constructor
		skEmission_Tabulated_HeightWavelength& operator =						( const skEmission_Tabulated_HeightWavelength& other );	// Dont allow assignment operator

	public:
											skEmission_Tabulated_HeightWavelength	();
		virtual							   ~skEmission_Tabulated_HeightWavelength	() override {};
		bool								SetEmissionTable					( const nx2dArray<double>& extinction, const nx1dArray<double>& wavelens, const nx1dArray<double>& heights_meters);
		bool								LoadHeightWavelengthProfileFromFile	( const char* filename );


	public:
		virtual bool						UpdateCache					( const GEODETIC_INSTANT& pt)  override { return true;}			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						UpdateLocation				( const GEODETIC_INSTANT& pt, bool isground) override;				//!< Allows the emission object to update any internal caches to this time and location
		virtual bool						IsotropicEmission			( double wavenumber, double* isotropicradiance) override;				//!< Calculate the isotropic emission as a radiance at the specified wave-number at the location specified in last call to SetAtmosphericState
};

