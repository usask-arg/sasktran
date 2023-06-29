
/*-----------------------------------------------------------------------------
 *					class skEmission					 */
/**	\ingroup skemission
 *	\par Overview
 *	The skEmission class is a base class that provides an interface for calculating emissions within the atmosphere.
 *	The class is currently designed to simulate thermal and photo-chemical emissions. The skEmission interface allows the user to specify
 *	a location and then returns a radiance at the requested wavenumber.
 */
/*-------------------------------------------------------------------------*/

class skEmission : public nxUnknown
{
	private:
		skEmission&							operator =							( const skEmission& other );  // =delete; Visual Studio 2012 does not like this yet			// Dont allow assignment
											skEmission							( const skEmission& other );  // =delete;			// Dont allow copy constructor

	public:
											skEmission							(){};
		virtual							   ~skEmission							(){}

	public:
		virtual bool						UpdateCache							( const GEODETIC_INSTANT& pt) = 0;														//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						UpdateLocation						( const GEODETIC_INSTANT& pt, bool isground) = 0;										//!< Changes the location of the emission object.
		virtual bool						IsotropicEmission					( double wavenumber, double* isotropicradiance)      = 0;								//!< Calculate the isotropic emission as a radiance at the specified wave-number. Photons/cm2/sec/steradian/nm 
		virtual bool						IsotropicEmissionArray				( const std::vector<double>& wavenumber, std::vector<double>* isotropicradiance);		//!< Calculate the isotropic emissions for an array of wavelengths. Photons/cm2/sec/steradian/nm 
};

/*-----------------------------------------------------------------------------
 *					skEmission_Constant		 2015- 3- 2 */
/** **/
/*---------------------------------------------------------------------------*/

class skEmission_Constant : public skEmission
{
	private:
		double								m_radiance; 

	public:
											skEmission_Constant( double constantvalue = 0.0 )			{ m_radiance = constantvalue;}
		virtual							   ~skEmission_Constant() override		{ }
		void								SetValue							( double value)	{ m_radiance = value;}
		virtual bool						UpdateCache							( const GEODETIC_INSTANT& /*pt*/ ) override				{ return true;};	//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						UpdateLocation						( const GEODETIC_INSTANT& /*pt*/, bool /* isground*/) override  { return true;}														//!< Allows the emission object to update any internal caches to this time and location
		virtual bool						IsotropicEmission					( double /*wavenumber*/, double* isotropicradiance)  override { *isotropicradiance = m_radiance; return true;}					//!< Calculate the isotropic emission as a radiance at the specified wave-number at the location specified in last call to SetAtmosphericState

};
