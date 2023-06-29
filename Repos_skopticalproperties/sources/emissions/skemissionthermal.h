/*---------------------------------------------------------------------------
 *					class skEmission_Thermal					2017-08-23*/
/**	\ingroup skopticalpropmisc
 *	A class for adding thermal emissions with isotropic emission given by
 *	the Planck function. 
 **/
/*-------------------------------------------------------------------------*/

class skEmission_Thermal : public skEmission
{
	private:
		skClimatology*						m_atmospheric_state;
		double								m_groundemissivity;
		double								m_temperatureK;
		double								m_absorptionperm;
		bool								m_isground;
		bool								m_absorptionisset;

	private:
		double								PlanckBlackbody				( double wavelen_nm, double temperature_K );
											skEmission_Thermal			( const skEmission_Thermal& other );	// Dont allow copy constructor
		skEmission_Thermal&					operator =					( const skEmission_Thermal& other );	// Dont allow assignment operator

	public:
											skEmission_Thermal			();
											skEmission_Thermal			( double groundemissivity );
		virtual							   ~skEmission_Thermal			();
		void								SetGroundEmissivity			( double value ) { m_groundemissivity = value; }
		void								SetAbsorptionPerM			( double value );																			// Sets the absorption coefficient for the current point in the atmosphere; must be called after SetAtmosphericState and before IsotropicEmission
		bool								SetAtmosphericState			( skClimatology* neutralatmosphere ) ;	//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						UpdateCache					( const GEODETIC_INSTANT& pt) override;
		virtual bool						UpdateLocation				( const GEODETIC_INSTANT& pt, bool isground ) override;										//!< Allows the emission object to update any internal caches to this time and location
		virtual bool						IsotropicEmission			( double wavenumber, double* isotropicradiance ) override;									//!< Calculate the isotropic emission as a radiance at the specified wave-number at the location specified in last call to SetAtmosphericState
};

