
/*-----------------------------------------------------------------------------
 *					class nxMeteorology		 2016- 5- 5*/
/** **/
/*---------------------------------------------------------------------------*/

class nxMeteorology
{
	private:
		double		m_latitude;
		double		m_phi;
		double		m_cos2phi;
		double		m_cosphi;
		double		m_sinphi;

	private:
		double			Latitude								() const { return m_latitude;}

	public:
		static double	GasConstantDryAir						();
		static double	GasConstantWaterVapour					();
						nxMeteorology							( double latitude);
		virtual 	   ~nxMeteorology							() {}
		virtual void	SetLatitude								( double latitude);						
		double			LocalGravityAtSeaLevel					( ) const;
		double			LocalEarthRadius						( ) const;
		double			BelowSurfaceLapseRate					() const { return 0.0065;}
		double			BelowSurfaceAlpha						(double *usertstar,double T0, double phisurf) const;
		double			ExtrapolateTemperatureBelowLowestLevel	( double T0, double H0, double h) const;
		double			ExtrapolatePressueBelowLowestLevel		( double P0, double T0, double H0, double h) const;
		double			GeometricHeightToGeopotential			( double z_in_meters ) const;
		double			GeopotentialToGeometricHeight			( double Phi ) const;
		double			GeopotentialToGeopotentialHeight		( double geopotential ) const;
		double			GeopotentialHeightToGeopotential		( double geopotentialheight ) const;
		double			GeopotentialHeightToGeometricHeight		( double geopotentialheight ) const;
		double			GeopotentialHeightFromGeometric			( double z_in_meters) const;
		double			GeometricFromGeopotentialHeight			( double geopotentialheight) const;
		double			PotentialTemperature					( double P, double T) const;
		double			VirtualTemperature						( double T, double q) const;
		double			PotentialTemperatureToTemperature		( double P, double Theta) const;
};


/*-----------------------------------------------------------------------------
 *				nxMeterology_ICAOSimple			 2016- 5- 24*/
/** Implements a simplified version of the International Civil Aviation Organization atmospheric
 *	model as described in the vertical interpolation section 1.3 of ECMWF IFS documentation,
 *	IFS Documentation - Cy40r1, PartI: Observations, November 2013. The simplified model has linearly
 *	decreasing temperature profile across the tropopause and a fixed temperature at all heights
 *	above the tropopause. The model is used to accurately interpolate the ECMWF pressure and geopotential data. 
 *
 * The temperature profile of the ICAO atmosphere below the tropopause is given by 
 * \f[
 * T_{ICAO}=T_{0}+\lambda\phi_{ICAO}
 * \f]
 * where \f$ \phi= \f$ geopotential above sea level, \f$ \lambda = -\frac{\Lambda}{g} \f$ , \f$ \Lambda=0.0065 Km^{-1} \f$ in the troposphere and 0 above the stratosphere, \f$g = \f$ gravity at sea level and \f$ T_0 = 288 \f$. The hydrostatic and ideal gas law equations yields
 * \f[
 * \frac{dp}{p}= -\frac{d\phi}{RT}
 * \f] 
 * substituting for T below the tropopause gives
 * \f[
 * \frac{dp}{p}= -\frac{d\phi}{R\left( T_{0} + \lambda\phi\right) }
 * \f]
 * Analytically integrating gives an expression allowing \f$p\f$ to be determined given \f$\phi\f$ 
 * \f[
 * \frac{p}{p_{0}}= \left(  1 + \frac{\lambda}{T_0}\phi\right) ^{-\frac{1}{\lambda R}}
 * \f] 
 * where \f$p_0 =\f$ 101325 Pa and is the pressure at sea level. The geopotential can be calculated given the pressure below the tropopause
 * \f[
 * \phi = \left[ \left( \frac{p}{p_0}\right) ^{-\lambda R} - 1\right] \frac{T_0}{\lambda}
 * \f]
 * The ICAO tropopause is given by \f$T=216.5K\f$. Thus,
 * \f[
 * \phi_{trop} = \left[ \frac{216.5 - T_0}{\lambda}\right]
 * \f]
 * and
 * \f[
 * \frac{p_{trop}}{p_{0}}= \left(  1 + \frac{\lambda}{T_0}\phi_{trop}\right) ^{-\frac{1}{\lambda R}}
 * \f]
 * above the tropopause, \f$p\f$ is given
 * \f[
 * \frac{p}{p_{trop}}= \exp -\left( \frac{\phi - \phi_{trop}}{RT_{trop}} \right)
 * \f]
 * and geopotential by
 * \f[
 * \phi = \phi_{trop} - \ln \left( \frac{p}{p_{trop}}\right)RT_{trop}
 * \f]
 **/
/*---------------------------------------------------------------------------*/

class nxMeterology_ICAOSimple: public nxMeteorology
{
	private:
		double			m_T0;
		double			m_P0;
		double			m_Ttrop;
		double			m_PHItrop;
		double			m_Ptrop;
		double			m_lambda;
		double			m_lambdaR;

	private:
		double			AnalyticallyIntegratePhiInTroposphere (double p) const;
		double			AnalyticallyIntegratePhiInStratosphere(double p) const;

	public:
						nxMeterology_ICAOSimple	( double latitude);
		double			T_ICAO					( double phi_icao ) const;
		double			Phi_ICAO				( double p) const;
};


/*-----------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_Heightentry		 2016- 5- 31*/
/**  Stores the pressure and temperature at a given altitude as interpolated
 *	from the ECMWF grib data by class nxMeteorology_EcwmfIFS_ObservationOperator
 **/
/*---------------------------------------------------------------------------*/

class nxMeteorology_EcwmfIFS_HeightEntry
{
	public:
		double	heightm;
		double	geopotential;
		double	pressure_pa;
		double	temperature_K;
		bool	isvalid;

	public:
		nxMeteorology_EcwmfIFS_HeightEntry( ) { isvalid = false; heightm = std::numeric_limits<double>::quiet_NaN(); }
};

/*-----------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator		 2016- 5- 30*/
/** A class that implements the ECMWF Observation Operator described in IFS DOcumentation Cy40r1
 *	Operational Implementation 22 November 2013, Part I: Observations. The purpose of the class
 *	is to provide vertical interpolation of the ecmwf model fields so we can
 *	extract pressure and temperature at any altitude within an skClimatology framework. This class
 *	will calculate pressure and temperature for any geometric altitude between sea-level and the
 *	upper altitude of the ECMWF. It does not provide extrapolation beyond the upper altitude.
 *
 *	\par Interpolation  between surface and upper altitude
 *	The ECMWF surface (approximately) follows the topographic surface of Earth and we use the
 *	algorithms described in the IFS documentation to derive pressure and temperature at any
 *	geometric altitude contained by the ECMWF data. Note that interpolation between levels includes
 *	a small correction to account for the changes in temperature between levels (by referencing the ICAO atmosphere).
 *	Apparently this is somewhat better than a simple linear interpolation in Log(p) space.

 * \par Below Surface Extrapolation
 *	The class also implements the atmospheric below surface extrapolation outlined in the IFS observation documentation 
 *	and described in the references:
 *	Briefly, we use a fixed lapse rate for all altitudes below the surface and perform adiabatic extrapolation
 * \f[ \alpha = 0.0065\frac{R_d}{g} \f]
 *	thus the temperature at a given geopotential below the surface is given by
 *	\f[ T = T_* + \frac{\alpha}{R_d}\left( \phi_s - \phi\right) \f]
 *  where \f$ T_*\f$ is the surface temperature
 *	The geopotential below the surface at a given pressure is given by,
 *	\f[ \phi = \phi_s - \frac{R_{d}T_*}{\alpha} \left[  \left( \frac{p}{p_s}\right)^\alpha - 1\right]   \f]
 *  and the pressure at a given geopotential below the surface is given by
 *	\f[ p = p_s \left[ 1 + \frac{\alpha\left(\phi_s-\phi\right) }{R_{d}T_*}\right]^{1/\alpha} \f]
 * which for small \f$ \alpha \f$ can be approximated as 
 *	\f[ p \approx p_{s}\exp \left[ \frac{\phi_s}{R_{d}T_*}\left(1 - \frac{1}{2}x + \frac{1}{3}x^2  \right) \right] \f]
 *	where,
 *	\f[ x = \frac{\alpha\phi_s}{R_{d}T_*} \f]
 *
 * \par  References
 *	- Vertical Interpolation and Truncation of Model-Coordinate Data, Kevin E. Trenberth, Jeffery C. Berry, Lawrence E. Buja, December 1993.
 *	- IFS Documentation-Cy38r1, June 2012, Part III: Dynamics and Numerical Procedures
 *	- Ryab El Khatib, July 2002, Full-Pos Scientists Guide in Arpege IFS cycle, CY24T1 and ALADIN cycle AL15 with ARPEGE/IFS cycle CY24T1.
 *	- An Energy and Angular Momentum Conserving Vertical Finite-Difference Scheme and Hybrid Vertical Coordinate, Simmons and Burridge, 1981
 *	- The Calculation of Geopotential and the Pressure Gradient in the ECMWF Atmospheric Model: Influence on the simulation of the
 *	  polar atmosphere and on temperature analyses: Simmons and Chen, 1991. (also in Technical Report No. 66 from ECMWF Reserach Dept).
 *
**/
/*---------------------------------------------------------------------------*/

class nxMeteorology_EcwmfIFS_ObservationOperator : public nxMeteorology
{
	private:
		double						m_longitude;
		nxMeterology_ICAOSimple		m_icao;							//!< The ICAO model, used for vertical interpolation
		double						m_psurf;						//!< Surface Pressure
		double						m_phisurf;						//!< Surface geopotential
		double						m_phisurf_icao;					//!< Surface geopotential of ICAO model
		double						m_Tstar;						//!< Surface temperature in K
		double						m_T0;							//!< Mean sea level temperature
		std::vector<double>			m_lnp;							//!< Log pressure on full model levels
		std::vector<double>			m_p;							//!< Pressure on full model levels
		std::vector<double>			m_T;							//!< Temperature on full model levels
		std::vector<double>			m_Tv;							//!< Virtual temperature on full model levels, used for vertical interpolation
		std::vector<double>			m_Ticao;						//!< ICAO model temperature of full model levels
		std::vector<double>			m_phi;							//!< Geopotential on full model levels.
		std::vector<double>			m_deltaphik;					//!< Delta Phi on full model levels. Difference between ECMWF and ICAO geopotentials from surface, used for vertical interpolation

	private:
		size_t						NumLevels										() const { return m_lnp.size();}
		double						Interpolate_Temperature							( double P0, double T0, double P1, double T1, double P) const;
		bool						GeopotentialToPressureAndTemperature			( double phi, double * p, double * T ) const;
		bool						GeopotentialToBelowGroundPressureAndTemperature	( double phi, double * p, double * T ) const;
		double						ModelLevelPressure								( double Pkminushalf, double Pkplushalf ) const;
		double						AboveModelPressure								(double geopotential) const;
		double						AboveModelTemperature							( double geopotential) const;


	// -- Methods to support below surface extrapolation
	private:
		void						SetSurfaceTemperature							( double Tnl, double Pnl );
		double						BelowSurfacePressure							( double geopotential ) const;
		double						BelowSurfaceGeopotential						( double pressure     ) const;
		double						BelowSurfaceTemperature							( double geopotential ) const;

	public:
									nxMeteorology_EcwmfIFS_ObservationOperator		();
		virtual					   ~nxMeteorology_EcwmfIFS_ObservationOperator		();

		bool						Configure										(	double						latitude, 
																						double					    longitude,
																						double						lnsurfacepressure, 
																						double						surfacegeopotential, 
																						const nx1dArray<double>&	T,
																						const nx1dArray<double>*	phalflevel = nullptr,	// You must specify pressure either on half levels
																						const nx1dArray<double>*	pfulllevel = nullptr,	// or on full levels but not both.
																						const nx1dArray<double>*	q          = nullptr);	// humidity is optional. It is 0.0 if not defined.

		bool						LogPressureToGeopotential						( double lnp,   double * phi ) const;
		bool						UpdateHeightEntry								( double heightm, nxMeteorology_EcwmfIFS_HeightEntry* entry ) const;
		double						UpperAltitude									() const { return GeopotentialToGeometricHeight(m_phi.front());}
		double						Longitude										() const { return m_longitude;}
		const std::vector<double>&	LnPArray										() const { return m_lnp;}
		const std::vector<double>&	PArray											() const { return m_p;}									
		const std::vector<double>&	PhiArray										() const { return m_phi;}									
		const std::vector<double>&	TArray											() const { return m_T;}									
};

