
/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulatedLUT		2014-2-12*/
/** Implements an internal Lookup table o calculate the Voigt function.
**/
/* --------------------------------------------------------------------------- */

class skSpectralLineShape_VoigtTabulatedLUT : public nxUnknown
{
	private:
		std::vector< std::vector<double> >		m_table;						// A patchy 2-D array. The first index indexes Voigt Y. The second index is for Voigt X
		double									m_xresolution;
		double									m_yresolution;
		size_t									m_numx;
		size_t									m_numy;

	public:
												skSpectralLineShape_VoigtTabulatedLUT();
		size_t									Yindex							( double Y ) const;
		size_t									Xindex							( double X ) const;
		size_t									NumY							( )			 const { return m_numy;}
		size_t									NumX							( )			 const { return m_numx;}
		bool									CalculateXEntries				( std::vector<double>& K, size_t yindex) const;
		nxString								LoadDirectoryNameFromRegistry	( );
		nxString								FullCacheName					( size_t yindex);
		bool									ReadCacheFile					( std::vector<double>&	K, const char* filename );
		bool									WriteCacheFile					( const std::vector<double>&	K, const char* filename );
		bool									CheckTableEntry					( size_t yindex);
		double									At								( size_t iy, size_t ix) const { return m_table.at(iy).at(ix);}
};


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtTabulated		2014-2-12*/
/** **/
/*---------------------------------------------------------------------------*/

class skSpectralLineShapeStorageBuffer_VoigtTabulated
{
	private:
		double												m_nu00;			//!< The wave number of the spectral line
		double												m_aD;			//!< The doppler halh width
		double												m_aL;			//!< The Lorentz half width.
		double												m_Y;			//!< The current value of Y set from the last call to SetLineParams
		static skSpectralLineShape_VoigtTabulatedLUT		m_lookuptable;	//!< The Voigt lookup table.

	public:
		double												Nu00	() const { return m_nu00;}
		double												Y		() const { return m_Y;}
		double												aD		() const { return m_aD;}
		double												aL		() const { return m_aL;}
		const skSpectralLineShape_VoigtTabulatedLUT&		LUT		() const { return   m_lookuptable;}
		skSpectralLineShape_VoigtTabulatedLUT*				LUTVar	()		 { return  &m_lookuptable;}


	public:
						skSpectralLineShapeStorageBuffer_VoigtTabulated();
		virtual		   ~skSpectralLineShapeStorageBuffer_VoigtTabulated();
		void			ResetLineParams		();															// Reset the line parameters to "invalid"
		bool			SetLineParams		( double nu00,												// Set line parameters
											  double pressure,
											  double partialpressure,
											  double temperature,
											  double tref,
											  double tempcoeff,
											  double mass,
											  double airhalfwidth,
											  double selfhalfwidth);
		bool			IsConfigured		( ) const					{ return (NXFINITE(m_Y));}	// See if the user has properly called SetLineParams

};

/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtTabulated		2013-3-14*/
/** \ingroup spectralline
 *	**/
/*---------------------------------------------------------------------------*/

class skSpectralLineShape_VoigtTabulated : public skSpectralLineShape
{
	private:
	
	private:
		bool									NormalizedShape						( double X, double Y, double *W, skSpectralLineShape_VoigtTabulatedLUT* lut) const;

	public:
												skSpectralLineShape_VoigtTabulated	();
		virtual bool							LineShapeFunction					( double nu, double* uservalue, const skSpectralLine* spectralline, skSpectralLineShapeStorageBuffer* storagebuffer ) const;
		virtual bool							ConfigureLineParameters				( const skSpectralLine* spectralline,  double	temperature, double	pressure, const GEODETIC_INSTANT& geopt, skClimatology* atmopshericstate, skSpectralLineShapeStorageBuffer*	storagebuffer );
		virtual bool							CreateStorageBuffer					( skSpectralLineShapeStorageBuffer** storagebuffer );
};
