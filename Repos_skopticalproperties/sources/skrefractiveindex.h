
#include <complex>


/*---------------------------------------------------------------------------
 *					class skRTRefractiveIndex						2003-10-10 */
/**	\ingroup skrefindex
 *	Base class for calculating refractive index or differnce media.
 *	This is typically used for calculating the refractive index of various
 *	particles for Mie scattering calculations.
 */
/*-------------------------------------------------------------------------*/

class skRTRefractiveIndex : public nxUnknown
{
	public:
		virtual						   ~skRTRefractiveIndex()	{};
		virtual std::complex<double>	RefractiveIndex         ( double wavenum )                  const = 0;
		virtual const char*				ChemicalName			()                                  const = 0;
};


/*---------------------------------------------------------------------------
 *					class skRTRefractiveIndex_MoistAir		2003-10-23	*/
 /** \ingroup skrefindex
  *	Calculates the refractive index of moist air.*/
 /*-------------------------------------------------------------------------*/

class skRTRefractiveIndex_MoistAir : public skRTRefractiveIndex
{
	private:
		double	m_T;					// Temperature in kelvins
		double	m_P;					// Total Pressure in pascals
		double	m_WaterVapPres;			// Water Vapour Partial Pressure in pascals
		double	m_ppmCO2;				// ppm of CO2 in the dry air mixture
		double	m_rho_axs;				// Density of dry air at 15C, 101325 Pa, (CO2 present but no water)
		double	m_rho_ws;				// density of pure water vapour at 20C, 1333 Pa.


	private:
		double							MoistAirDensity				( double xw, double* rhoair = NULL, double* rhowater = NULL) const;
		double							Compressibility				( double xw ) const;										// Compressibility of moist air
		double							MolarMassDryAir				() const;													// Depends only upon CO2 density in ppm
		double							MolarFractionWaterVapour	() const;
		double							EnhancementFactor			() const;
		double							SaturationVapourPressure	() const;
		void							UpdateStandardDensities		();
		bool							DeepCopy					( const skRTRefractiveIndex_MoistAir& other );

	public:
		std::complex<double>			Refractivity					( double wavenum ) const;
		bool							Set_CO2ppm						( double co2ppm  );
		bool							Set_Temperature					( double kelvins );
		bool							Set_WaterVapourPressure			( double wvp);
		bool							Set_TotalPressure				( double totalpressure);
		bool							Set_PTandRH						( double totalpressure, double kelvins, double relativehumidity);
		double							DryAirRefractivityAtSTP			( double vacuumwavenum ) const;
		double							WaterVapourRefractivityAtSTP	( double vacuumwavenum ) const;

	public:
										skRTRefractiveIndex_MoistAir	( );
		virtual						   ~skRTRefractiveIndex_MoistAir	( ){};
		virtual std::complex<double>	RefractiveIndex					( double wavenum )  const override { return Refractivity(wavenum) + std::complex<double>(1,0);} 
		virtual const char*				ChemicalName					( ) 				const override {return "moistair";}
};

/*---------------------------------------------------------------------------
 *					class skRTRefractiveIndex_Tabulated			2003-10-10*/
/**	\ingroup skrefindex
 *	Calculates the refractive index of a medium by interpolating tables
 */
 /*-------------------------------------------------------------------------*/

class skRTRefractiveIndex_Tabulated : public skRTRefractiveIndex
{
	private:
		nx1dArray<double>				m_wavelen;				// wavelength in nanometers, monotonically ascending
		nx1dArray<double>				m_real;					// Real part of Refractive Index
		nx1dArray<double>				m_imag;					// Imaginary part of Refractive index
		
	protected:
		void							operator=		( const skRTRefractiveIndex_Tabulated& other );


	public:
		virtual						   ~skRTRefractiveIndex_Tabulated(){};
		bool							AttachToStatic(double* data, size_t nelements);
		virtual std::complex<double>	RefractiveIndex( double wavenum ) const override;
		virtual const char*				ChemicalName() const override { return "usertable"; }
};


/*---------------------------------------------------------------------------
 *					class skRTRefractiveIndex_ICE				2003-10-10 */
/**	\ingroup skrefindex
 *	Tabulates the complex refractive index of ICE. Used in MIE scattering
 *	calculations.
 */
 /*-------------------------------------------------------------------------*/

class skRTRefractiveIndex_ICE : public skRTRefractiveIndex_Tabulated
{
	public:
										skRTRefractiveIndex_ICE();
		virtual						   ~skRTRefractiveIndex_ICE(){};
		virtual const char*				ChemicalName			() const {return "ice";};
};

/*---------------------------------------------------------------------------
 *					class skRTRefractiveIndex_H2SO4		2003-10-10			*/
/**	\ingroup skrefindex
 *	Tabulates the complex refractive index of sulphate particles. Used in
 *	MIE scattering calculations.
 */
 /*-------------------------------------------------------------------------*/

class skRTRefractiveIndex_H2SO4 : public skRTRefractiveIndex_Tabulated
{
	public:
										skRTRefractiveIndex_H2SO4();
		virtual						   ~skRTRefractiveIndex_H2SO4(){};
		virtual const char*				ChemicalName			() const {return "h2so4";}
};

/*---------------------------------------------------------------------------
 *					class skRTRefractiveIndex_Dust		2003-10-10			*/
/**	\ingroup skrefindex
 *	Tabulates the complex refractive index of dust particles. Used in
 *	MIE scattering calculations.
 */
/*-------------------------------------------------------------------------*/

class skRTRefractiveIndex_Dust : public skRTRefractiveIndex_Tabulated
{
	public:
										skRTRefractiveIndex_Dust();
		virtual						   ~skRTRefractiveIndex_Dust(){};
		virtual const char*				ChemicalName			() const {return "dust";};
};

/*---------------------------------------------------------------------------
 *					class skRTRefractiveIndex_Water		2003-10-10			*/
/**	\ingroup skrefindex
 *	Tabulates the complex refractive index of water particles. Used in
 *	MIE scattering calculations.
 */
/*-------------------------------------------------------------------------*/

class skRTRefractiveIndex_Water : public skRTRefractiveIndex_Tabulated
{
	public:
										skRTRefractiveIndex_Water();
		virtual						   ~skRTRefractiveIndex_Water(){};
		virtual const char*				ChemicalName			() const {return "water";}
};


