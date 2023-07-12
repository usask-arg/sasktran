#include "nxbase_geodesy.h"
#include "module/math/zbrent.h"
#include <algorithm>


/*-----------------------------------------------------------------------------
 *					nxMeteorology::GasConstantDryAir		 2016- 5- 4*/
/** **/
/*---------------------------------------------------------------------------*/

double  nxMeteorology::GasConstantDryAir()
{	
	return 287.058;		// (R/Mr * 1000), where R = Na.k, Mr = Molecular weight of air
}


/*-----------------------------------------------------------------------------
 *					nxMeteorology::GasConstantWaterVapour		 2016- 5- 4*/
/** **/
/*---------------------------------------------------------------------------*/

double  nxMeteorology::GasConstantWaterVapour()
{
	return 461.53;
}


/*-----------------------------------------------------------------------------
 *					nxMeteorology::nxMeteorology		 2016- 5- 5*/
/** **/
/*---------------------------------------------------------------------------*/

nxMeteorology::nxMeteorology( double latitude)
{
	SetLatitude(latitude);
}


/*-----------------------------------------------------------------------------
 *					nxMeteorology::SetLatitude		 2016- 5- 31*/
/** **/
/*---------------------------------------------------------------------------*/

void nxMeteorology::SetLatitude( double latitude)
{
	m_latitude = latitude;
	m_phi      = nxmath::DegreesToRadians(latitude);
	m_sinphi   = sin(m_phi);
	m_cosphi   = cos(m_phi);
	m_cos2phi  = cos(2*m_phi);
}


/*-----------------------------------------------------------------------------
 *					nxMeteorology::LocalGravityAtSeaLevel		 2016- 5- 4*/
/** **/
/*---------------------------------------------------------------------------*/

double nxMeteorology::LocalGravityAtSeaLevel( ) const
{
	double			G;

	G       = 9.80616*( 1 + ( -0.002637 + 0.0000059*m_cos2phi)*m_cos2phi);		// Gravity at ground as function of latitude
	return G;
}


/*-----------------------------------------------------------------------------
 *					nxMeteorology::LocalEarthRadius		 2016- 5- 4*/
/** **/
/*---------------------------------------------------------------------------*/

double nxMeteorology::LocalEarthRadius() const
{
	const double	Rmax = 6378137.0;			// Semi major axis in meters
	const double	Rmin = 6356752.0;			// Semi minor eaxis in meters
	double			Re;

	Re      = sqrt( nxmath::sqr(Rmax*m_cosphi) + nxmath::sqr(Rmin*m_sinphi));	// "Radius of Earth" as a function latitude
	return Re;
}

/*---------------------------------------------------------------------------
 *					GeometricHeightToGeopotentialHeight				*/
/**
 *	@param double z_in_meters
 *	The geometric height in meters above sea level
 *
 *	@param double latitude
 *	The latitude in degrees of the conversion
 *
 *	@par RETURN
 *	Returns the geopotential of the geometric height. ( As an FYI This is approximately the geometric height multiplied by 10).
 *
 *	@par HISTORY
 * 	2016-5-2. Copied from the skClimatology ECMWF code and placed into centralised location.
 **/
/*------------------------------------------------------------------------*/

double nxMeteorology::GeometricHeightToGeopotential( double z_in_meters ) const
{
	double			G;
	double			Re;
	double			H;

	Re = LocalEarthRadius();						// "Radius of Earth" as a function latitude
	G  = LocalGravityAtSeaLevel( );					// Gravity at ground as function of latitude
	H  = z_in_meters*G*Re/( Re + z_in_meters );		// Conversion from geometric to geopotential
	return H;
}


/*---------------------------------------------------------------------------
 *					nxMeteorology::GeopotentialToGeometricHeight				*/
/**
 *	@param double H
 *	Geopotential
 *	@param double latitude
 *	latitude
 *	@par RETURN
 *	height in meters
 *	@par HISTORY
 * 	2016-5-6
 **/
/*------------------------------------------------------------------------*/

double nxMeteorology::GeopotentialToGeometricHeight( double H ) const
{
	double			G;
	double			Re;
	double			z;

	Re = LocalEarthRadius();				// "Radius of Eart" as a function latitude
	G  = LocalGravityAtSeaLevel( );			// Gravity at ground as function of latitude
	z  = H*Re/( Re*G - H);					// Conversion from geometric to geopotential
	return z;
}

/*---------------------------------------------------------------------------
 *					GeopotentialToGeopotentialHeight				*/
/**
 *	@param double geopotential
 *		The incoming geopotential (typically ~1/10'th the geometric height)
 *
 *	@par RETURN
 *	the geopotential height in meters. This is not the same as the geometric height.
 *	@par HISTORY
 * 	2016-5-2
 **/
/*------------------------------------------------------------------------*/

double nxMeteorology::GeopotentialToGeopotentialHeight ( double geopotential ) const
{
	const double	g0   = 9.80665;				// Accel due to gravity used fro potential heighty calcs

	return geopotential/g0; 
}



/*---------------------------------------------------------------------------
 *                              nxMeteorology::GeopotentialHeightToGeopotential                             2019-10-15 */
/** **/
/*---------------------------------------------------------------------------*/

double nxMeteorology::GeopotentialHeightToGeopotential ( double geopotentialheight ) const
{
	const double	g0   = 9.80665;				// Accel due to gravity used fro potential heighty calcs

	return geopotentialheight*g0; 
}

/*---------------------------------------------------------------------------
 *       nxMeteorology::GeopotentialHeightToGeometricHeight       2019-12-13 */
/** **/
/*---------------------------------------------------------------------------*/

double nxMeteorology::GeopotentialHeightToGeometricHeight( double geopotentialheight ) const
{
	double phi = GeopotentialHeightToGeopotential (geopotentialheight);
	double h   = GeopotentialToGeometricHeight(phi);
	return h;
}


/*-----------------------------------------------------------------------------
 *					nxMeteorology::VirtualTemperature		 2016- 5- 25*/
/** Calculates the virtual temperature given the actual temperature T and
 *	specific humidity, q.
 *	\f[
 *	T_{v}=T\left[ 1+\left( \frac{R_{vap}}{R_{dry}} - 1\right)q\right]
 *   \f]
 *	@param double T
 *	The temperature of the air in Kelvins
 *
 *	@param double q
 *	The specific humidity
 *
 *	@par RETURN
 *	double, the virtual temperature in Kelvins
**/
/*---------------------------------------------------------------------------*/

double	nxMeteorology::VirtualTemperature( double T, double q) const
{
	double Tv;
	double Rd	= GasConstantDryAir();
	double Rv   = GasConstantWaterVapour();
	double f    = (Rv/Rd - 1);

	Tv = T*(1.0 + f*q);
	return Tv;
}

/*---------------------------------------------------------------------------
 *					nxMeteorology::PotentialTemperature				*/
/** Calculates the potential temperature assuming dry air of the given 
 *	pressure and temperature.
 *	\f[
 *	\theta=T\left( \frac{100000}{p}\right)^{{R}/{c_p}}
 *	\f]
 *
 *	@param double P
 *	The pressure of the air in pascals
 *
 *	@param double T
 *	The temperature of the air in Kelvins
 *
 *	@par RETURN
 *	double
 *	@par HISTORY
 * 	2016-5-5
 **/
/*------------------------------------------------------------------------*/

double nxMeteorology::PotentialTemperature( double P, double T) const
{
	double			theta;
	const double	kappa = 2.0/7.0;			// R/Cp  = 2/7 for dry air

	theta = T*pow(100000.0/P, kappa);
	return theta;
}



/*---------------------------------------------------------------------------
 *					nxMeteorology::PotentialTemperatureToTemperature				*/
/**	Calculates the temperature given the pressure and potential temperature.
 *	\f[
 *	T=\theta\left( \frac{100000}{p}\right)^{-{R}/{c_p}}
 *	\f]
 
 *	@param double P
 *	@param double Theta
 *	@par RETURN
 *	double
 *	@par HISTORY
 * 	2016-5-5
 **/
/*------------------------------------------------------------------------*/

double nxMeteorology::PotentialTemperatureToTemperature( double P, double Theta) const
{
	double			T;
	const double	kappa = 2.0/7.0;			// R/Cp  = 2/7 for dry air

	T = Theta*pow(100000.0/P, -kappa);
	return T;
}


/*-----------------------------------------------------------------------------
 *					nxMeteorology::GeopotentialHeightFromGeometric		2005-7-21*/
/** Get the geopotential height from the specified geometric height (in meters) at the
  *	the current latitude and longitude
 **/
/*---------------------------------------------------------------------------*/

double nxMeteorology::GeopotentialHeightFromGeometric( double z_in_meters ) const
{
	double			phi;
	double			H;
	phi = GeometricHeightToGeopotential(z_in_meters ) ;
	H   = GeopotentialToGeopotentialHeight( phi );
	return H;
}


/*---------------------------------------------------------------------------
 *       nxMeteorology::GeometricFromGeopotentialHeight			2019-10-15 */
/** **/
/*---------------------------------------------------------------------------*/
double nxMeteorology::GeometricFromGeopotentialHeight( double geopotentialheight) const
{
	double			phi;
	double			H;
	phi = GeopotentialHeightToGeopotential(geopotentialheight ) ;
	H   = GeopotentialToGeometricHeight( phi );
	return H;
}

/*-----------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::BelowSurfaceAlpha		 2016- 6- 8*/
/** Claculates the \f$ \alpha \f$ defined by Trenberth et al. 1993 for geopotential
 *	and pressure calculations below the surface. This implements equation (6) and follows
 *	the adjustments presented in equation (13) through (14.3).  The IFS Observations code
 *	seems to follow Trenberth even though the IF observations document is not quite as clear
 *	in this area.
 *
 *	Note the \f$\alpha \f$calculated here is not the same as the \f$ \alpha \f$ used for
 *	temperature calculations below the surface. 
 *
 **/
/*---------------------------------------------------------------------------*/

double nxMeteorology::BelowSurfaceAlpha(double *usertstar, double T0, double phisurf) const
{
	double	tstar;
	double  g         = LocalGravityAtSeaLevel();
	double	Rd        = GasConstantDryAir();
	double	lapserate = BelowSurfaceLapseRate();
	double	alpha;

	tstar    = *usertstar;											// Get a copy o fthe surface temperature
	alpha    = lapserate*Rd/g;									// and use a standard lapse rate (0.0065 K/m)
	if ( (tstar <= 290.5) && (T0 > 290.5) )
	{
		alpha = Rd*(290.5 - tstar)/phisurf;
	}
	if ( ( tstar > 290.5) && ( T0  > 290.5 ))
	{
		tstar = 0.5*( 290.5 + tstar);
	}
	else if ( tstar < 255.0)
	{
		tstar = 0.5*( 255.0 + tstar);
	}
	*usertstar = tstar;
	return alpha;
}


/*---------------------------------------------------------------------------
 *     nxMeteorology::ExtrapolateTemperatureBelowLowestLevel      2019-10-21 */
/** Calculates the temperature at geometric height "h" below geometric height 
 *	H0 which is at temperature T0. Assumes a fixed tropospheric lapse rate
 **/
/*---------------------------------------------------------------------------*/

double nxMeteorology::ExtrapolateTemperatureBelowLowestLevel( double T0, double H0, double h) const
{
	return T0 + (H0-h)*BelowSurfaceLapseRate();
}



/*---------------------------------------------------------------------------
 *                             const                              2019-10-21 */
/*** Calculates the pressure at geometric height "h" below geometric height 
 *	H0 which is at temperature T0 and pressure P0. Assumes hydrostatic equilibrium
	and a fixed tropospheric lapse rate
 **/
/*---------------------------------------------------------------------------*/

double nxMeteorology::ExtrapolatePressueBelowLowestLevel( double P0, double T0, double H0, double h) const
{
	double  g      = LocalGravityAtSeaLevel();
	double	Rd     = GasConstantDryAir();
	double	lambda = BelowSurfaceLapseRate();
	double  offset;

	offset = log(  1.0 + (lambda*(H0-h)/T0))/(Rd*lambda);
	if (h > H0) offset = - offset;
	return P0*exp(offset);
}



/*-----------------------------------------------------------------------------
 *					nxMeterology_ICAOSimple::nxMeterology_ICAOSimple		 2016- 5- 24*/
/** **/
/*---------------------------------------------------------------------------*/

nxMeterology_ICAOSimple::nxMeterology_ICAOSimple( double latitude) : nxMeteorology(latitude)
{
	double	g;

	g       = LocalGravityAtSeaLevel();
	m_lambda  = -0.0065/g;
	m_lambdaR = m_lambda*GasConstantDryAir();
	m_T0      = 288.0;
	m_P0      = 101325.0;
	m_Ttrop   = 216.5;
	m_PHItrop = (m_Ttrop-m_T0)/m_lambda;
	m_Ptrop   = m_P0*pow( 1 + m_PHItrop*m_lambda/m_T0, -1.0/(m_lambdaR) );
}

/*-----------------------------------------------------------------------------
 *					nxMeterology_ICAOSimple::T_ICAO		 2016- 5- 24*/
/** Calculates the temperature, \f$ T_{ICAO}\f$, given geopotential \f$\phi_{ICAO} \f$
 *	Uses the ICAO troposphere temperature profile but is constant above the troposphere.
**/
/*---------------------------------------------------------------------------*/

double nxMeterology_ICAOSimple::T_ICAO( double phi_icao ) const
{
	double T;
	T = (phi_icao <= m_PHItrop) ? m_T0 + m_lambda*phi_icao 
		                        : m_Ttrop;
	return T;
}


/*-----------------------------------------------------------------------------
 *					nxMeterology_ICAOSimple::IntegratePhiInTroposphere		 2016- 6- 7*/
/** **/
/*---------------------------------------------------------------------------*/

double nxMeterology_ICAOSimple::AnalyticallyIntegratePhiInTroposphere(double p) const
{
	double phi;

	phi = (pow( p/m_P0, -m_lambdaR) - 1.0 )*m_T0/m_lambda ;
	return phi;
}


/*-----------------------------------------------------------------------------
 *					nxMeterology_ICAOSimple::IntegratePhiInStratosphere		 2016- 6- 7*/
/** **/
/*---------------------------------------------------------------------------*/

double nxMeterology_ICAOSimple::AnalyticallyIntegratePhiInStratosphere(double p) const
{
	double phi;

	phi =  m_PHItrop - log( p/m_Ptrop)*GasConstantDryAir()*m_Ttrop;
	return phi;
}

/*-----------------------------------------------------------------------------
 *					nxMeterology_ICAOSimple::Phi_ICAO		 2016- 5- 24*/
/** Calculates the geopotential  \f$ \phi_{ICAO}\f$ given geopotential \f$T_{ICAO} \f$
 *	Uses the ICAO troposphere profile but is constant above the troposphere.
**/
/*---------------------------------------------------------------------------*/

double nxMeterology_ICAOSimple::Phi_ICAO( double p) const
{
	double phi;


	phi  = (p >= m_Ptrop)  ? AnalyticallyIntegratePhiInTroposphere (p)
		                   : AnalyticallyIntegratePhiInStratosphere(p);
	return phi;
}



/*---------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::Constructor				*/
/** *
 **/
/*------------------------------------------------------------------------*/

nxMeteorology_EcwmfIFS_ObservationOperator::nxMeteorology_EcwmfIFS_ObservationOperator()
										   : m_icao(0.0), nxMeteorology(0.0)
{
	double nan     = std::numeric_limits<double>::quiet_NaN();
	m_longitude    = nan;
	m_psurf        = nan;						// Surface Pressure
	m_phisurf      = nan;						// Surface geopotential
	m_phisurf_icao = nan;						// Surface geopotential of ICAO model
	m_Tstar        = nan;						// Surface temperature in K
	m_T0           = nan;						// Mean sea level temperature
}


/*---------------------------------------------------------------------------
 * nxMeteorology_EcwmfIFS_ObservationOperator::~nxMeteorology_EcwmfIFS_ObservationOperator2019-10-16 */
/** **/
/*---------------------------------------------------------------------------*/
nxMeteorology_EcwmfIFS_ObservationOperator::~nxMeteorology_EcwmfIFS_ObservationOperator()
{
}

/*---------------------------------------------------------------------------
 *     nxMeteorology_EcwmfIFS_ObservationOperator::Configure      2019-10-16 */
/** Construct the ecmwf observation operator for this vertical profile. It is normally
 *	initialized with array data directly extracted from the ECMWF GRIB files at the ECMWF 
 *	indexed locations.
 *
 *	@param double latitude
 *		The latitude of this vertical profile in degrees.
 *
 *	@param double longitude
 *		The longitude of this vertical profile in degrees
 *
 *	@param const  nx1dArray<double>& T
 *		An array temperature on (typically 60) full pressure levels. There is one less element than in the half level pressures array.
 *		The highest index is the lowest altitude and is just above the surface (typically a few meters).
 *
 *	@param const nx1dArray<double>* phalflevel
 *		An array of ECMWF pressure on (typically 61) half level pressures.
 *		The highest index is the highest pressure and corresponds to the surface pressure.
 *
 *	@param const nx1dArray<double>* fulllevel
 *		An array of ECMWF pressure on (typically 60) full level pressures.
 *		The highest index is the highest pressure and corresponds to the surface pressure.
 *
 *	@param const  nx1dArray<double>* q
 *		An array of specific humdify on (typically 60) full pressure levels. 
 *
 *	@par HISTORY
 * 	2016-5-30
**/
/*---------------------------------------------------------------------------*/

bool nxMeteorology_EcwmfIFS_ObservationOperator::Configure(	double						latitude, 
															double					    longitude,
															double						lnsurfacepressure, 
															double						surfacegeopotential, 
															const nx1dArray<double>&	T,
															const nx1dArray<double>*	phalflevel,
															const nx1dArray<double>*	pfulllevel, 
															const nx1dArray<double>*	q)
{

	bool	ok, ok1, ok2, ok3, ok4;
	size_t	L;
	size_t	j;
	size_t	k;
	double	Rd = GasConstantDryAir();
	double	pplus;
	double	pminus;
	double	lnpp;
	double	p;
	double	phiicao;
	double	phihalflevel;
	double  dt;
	double	dp1;
	double	dp2;
	double  alpha;
	double	tv;

	m_icao.SetLatitude(latitude);
	this->SetLatitude(latitude);
	m_longitude = longitude;

	L   = T.size();
	ok1 = (q          == nullptr || (q->size()          == L));
	ok2 = (phalflevel == nullptr || (phalflevel->size() == (L+1)));
	ok3 = (pfulllevel == nullptr || (pfulllevel->size() == L));
	ok4 = (phalflevel == nullptr) ^ (pfulllevel == nullptr);
	ok = ok1 && ok4; 
	if (!ok4)
	{
		nxLog::Record(NXLOG_WARNING,"nxMeteorology_EcwmfIFS_ObservationOperator::Configure, you must specify one and only one of half pressure levels or full pressure levels.");
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxMeteorology_EcwmfIFS_ObservationOperator::Configure, the array sizes are not as expected. Thats not good");
	}
	else
	{
		m_psurf        = exp(lnsurfacepressure);						// save the surface pressure (pressure on lowest model half level)
		m_phisurf      = surfacegeopotential,							// Save the geopotential on the surface
		m_phisurf_icao = m_icao.Phi_ICAO( m_psurf);						// Get the surface potential from the ICAO model using the surface pressure
																	
		m_Tv.resize   (L);												// Allocate the internal storage arrays.
		m_Ticao.resize(L);												// These arrays are all on the model full level
		m_phi.resize(L);
		m_deltaphik.resize(L);
		m_lnp.resize(L);
		m_p.resize(L);
		m_T.resize(L);

		phihalflevel      = m_phisurf;									// Reset summation counter of equation 2.21 from IFS part III Dynamics and Numerical Procedures, Chapter 2 basic equations
		k = L;
		for (j = 0; j < L; j++)											// Loop from the groun to the upper altitude
		{																// so start at opposite end of array
			k--;														// and step down.
			if (phalflevel != nullptr)
			{
				pplus         = phalflevel->At(k+1);					// Get the "plus"  model half pressure level (Pk+1/2)
				pminus        = phalflevel->At(k);						// Get the "minus" model half pressure level (Pk-1/2)
				lnpp          = log( pplus/pminus );					// get the intermediate log term
				p			  = ModelLevelPressure( pminus, pplus);		// Calculate the model full pressure value
			}
			else
			{
				p = pfulllevel->At(k);
			}
			m_p.at(k)     = p;
			m_lnp.at(k)   = log(p);										// and take the logs.
			phiicao       = m_icao.Phi_ICAO(p);							// Get the ICAO geopotential of the model full-level pressure 
			m_Ticao.at(k) = m_icao.T_ICAO(phiicao);							// Get the ICAO temperature  at the model full-level pressure
			tv            = VirtualTemperature( T.At(k), (q == nullptr) ? 0.0 : q->At(k) );		// Get the virtual temperature at this model full level
			m_Tv.at(k)    = tv;											// And assign it to this model level.
			m_T.at(k)     = T.At(k);
			dt            = tv - m_Ticao.at(k);

			if ( k == 0) alpha = log(2);								// evaluate the alpha term, from either equation 2.23 IFS part III or equation 1.3 IFS Part 1
			else         alpha = 1.0 - pminus/(pplus-pminus)*lnpp;

			m_phi.at(k)        = phihalflevel + Rd*tv*alpha;			// Get phi on full model levels by evaluating equation 2.22 of IFS part III.
			phihalflevel      += Rd*tv*lnpp;
			dp1                = (m_phi.at(k) - m_phisurf);
			dp2                = (phiicao - m_phisurf_icao);
			m_deltaphik.at(k)  =  dp1 - dp2;							// Get Delta phi on full model levels by evaluating "the spirit" of equation 1.3 of IFS part 1. (Their explicit formula is actually not appropriate for our software).
		}
		SetSurfaceTemperature( m_T.at(L-2), m_p.at(L-2) );		// Set surface temperature using the model level above the lowest model level
	}

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxMeteorology_EcwmfIFS_ObservationOperator::Constructor, Error during initialization");
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::ModelLevelPressure				*/
/** Calculate the pressure on the ECWMF model level ('k') given the pressure on the two half levels
 *	surrounding the model level. We use the formula given at the top of page 7 in the 
 *	IFS Cy38r1 documentation June 2012, Part III: Dynamics and Numerical Procedures. Note the equation
 *	is unlabelled and is in the text between equation 2.11 and 2.12.
 *
 *	Note that we use a simple aritmetic mean to calculate the full-level pressure. This is different to 
 *	the formula used in Trenberth et al 1993 (their table 1 is not a simple mean) and to the interpolation formula  
 *	given in equation 1 of Ryab El Khatib, July 2002, Full-Pos Scientists Guide in Arpege IFS cycle
 *	CY24T1 and ALADIN cycle AL15 with ARPEGE/IFS cycle CY24T1. 
 *
 *	@param double Pkminushalf
 *	Pressure on the half level k-1/2. This will be the smaller pressure at the higher altitude.
 *
 *	@param double Pkplushalf
 *	Pressure on the half level, k+1/2. This will be the higher pressure at the lower altitude
 *
 *	@par RETURN
 *	The pressure on the model level.
 *
 *	@par HISTORY
 * 	2016-5-5
 **/
/*------------------------------------------------------------------------*/

double nxMeteorology_EcwmfIFS_ObservationOperator::ModelLevelPressure( double Pkminushalf, double Pkplushalf ) const
{
	return 0.5*(Pkminushalf + Pkplushalf);
}

/*---------------------------------------------------------------------------
 *					nxMeteorology::ECMWF_SurfaceTemperature				*/
/** Calculate the ECMWF lowest 1/2 level surface temperature, \f$T^{*}\f$,  and mean sea level temperature
 *	using the formulae given in equations 1.6 to 1.9 of IFS Cy40R1, Part I Observatons. Also see equation 5 of
 *	Trenberth, Berry and Buja, Vertical Interpolation and Truncation of Model-Coordinate Data (Dec 1993).
 *
 *	@param double Tnl
 *		Temperature at the lowest (integer) model level 
 *	@param double Pnl
 *		Pressure at the lowest 
 *	@par RETURN
 *	double
 *	@par HISTORY
 * 	2016-5-5
 **/
/*------------------------------------------------------------------------*/

void nxMeteorology_EcwmfIFS_ObservationOperator::SetSurfaceTemperature( double Tnl, double Pnl )
{
	double	tcorr;
	double	Tx = 290.5;
	double	Ty = 255.0;
	double Rd = GasConstantDryAir();
	double g  = LocalGravityAtSeaLevel();
	double lapserate = BelowSurfaceLapseRate();


	m_Tstar = Tnl*( 1.0 + lapserate*Rd/g*log(m_psurf/Pnl));				// Get the surface temperature from lapse rate calculation
	tcorr   = std::max( Ty, std::min(Tx, m_Tstar));						// do a correction for excessively hot or cold surfaces
	m_Tstar = 0.5*( m_Tstar + tcorr );									// and average the surface temperature
	m_T0    = m_Tstar + lapserate*m_phisurf/LocalGravityAtSeaLevel();	// Get the mean sea level temperature
}


	
/*---------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::BelowSurfacePressure				*/
/** Calculates the pressure at a geopotential below the ECMWF surface. This uses the a rearrangement of
 *	equation 1.10 in IFS Cy40R1, Part I Observations
 *
 *	@param double geopotential
 *		the geopotential at which the pressure is desired. This geopotential must be less than the surfacepotential. 
 *
 *	@par RETURN
 *	Returns the pressure at the requested geopotential height 
 *
 *	@par HISTORY
 * 	2016-5-5
 **/
/*------------------------------------------------------------------------*/

double nxMeteorology_EcwmfIFS_ObservationOperator::BelowSurfacePressure( double geopotential ) const
{
	double	Rd        = GasConstantDryAir();
	double	P         = std::numeric_limits<double>::quiet_NaN();
	double	deltaphi;
	double	x;
	double	tstar;
	double  alpha;

	NXASSERT(( geopotential <= m_phi.back() ));
	if ( geopotential >= m_phisurf)									// If we are above the surface then we should be between the first model level and the surface
	{																// so given this is usually quite close (< 100 meters) to the surface
		double	lnp0 = log(m_psurf);								// just perform linear interpolation in log pressure and geoptential height
		double  dy = m_lnp.back() - lnp0;
		double	dx = m_phi.back() - m_phisurf;
		P   =  lnp0 + dy/dx*(geopotential - m_phisurf);
		P   = exp(P);
	}
	else															// if we are below the surface then
	{																// perform the adiabtaic extrapolation downwards from the surface 
		tstar    = m_Tstar;											// Get a copy of the surface temperature
		alpha    = BelowSurfaceAlpha(&tstar, m_T0, m_phisurf);						// and calculate the (possibly adjusted) surface temperature and lapse rate/alpha coefficient
		deltaphi = (m_phisurf - geopotential);
		x        = alpha*deltaphi/(Rd*tstar)  + 1;
		P        = m_psurf*pow( x, 1.0/alpha);
	}
	NXASSERT(( NXFINITE(P)  ));
	return P;
}


/*-----------------------------------------------------------------------------
 *					QuadraticExtrapolation		 2016- 6- 9*/
/** **/
/*---------------------------------------------------------------------------*/

static double QuadraticExtrapolation( double x, double X[3], double Y[3] )
{
	double	y;
	double	x0 = X[0];
	double	x1 = X[1];
	double	x2 = X[2];
	double	y0 = Y[0];
	double	y1 = Y[1];
	double	y2 = Y[2];

	double dy0 =  (y0-y1)/(x0-x1);
	double dy1 =  (y1-y2)/(x1-x2);
	double a, b, c;

	c = (dy0 - dy1)/(x0-x2);
	b = dy0 - c*(x0+x1);
	a = y0 - (b + c*x0)*x0;
	
	y =  a + (b + c*x)*x;
	return y;

}
/*-----------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::AboveModelGeopotential		 2016- 6- 9*/
/** **/
/*---------------------------------------------------------------------------*/

double nxMeteorology_EcwmfIFS_ObservationOperator::AboveModelPressure( double geopotential) const
{
	double	X[3] = { m_phi.at(0), m_phi.at(1), m_phi.at(2)};
	double	Y[3] = { m_lnp.at(0), m_lnp.at(1), m_lnp.at(2)};
	double	lnp;

	lnp = QuadraticExtrapolation( geopotential, X, Y);
	return exp(lnp);
}


/*-----------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::AboveModelTemperature		 2016- 6- 9*/
/** **/
/*---------------------------------------------------------------------------*/

double nxMeteorology_EcwmfIFS_ObservationOperator::AboveModelTemperature( double /*geopotential*/) const
{
	return m_T.front();
}
/*---------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::BelowSurfaceGeopotential				*/
/** Calculates the geopotential at a (higher) pressure  level below the ECMWF model surface.
 *	Uses the adiabatic method described  in equation 15 of Trenberth et al. (dec 1993), slightly
 *	modified to agree with equation 1.10 in IFS Cy40R1, Part I Observations
 *
 *	@param double pressure
 *		the pressure at which the geopotential is desired. This pressure must be greater than the surfacepressure.
 *	@par RETURN
 *	Returns the geopotential at the requested pressure. Returns NaN if the requested pressure is less than the surface pressure.
 *	double
 *	@par HISTORY
 * 	2016-5-5
 **/
/*------------------------------------------------------------------------*/

double nxMeteorology_EcwmfIFS_ObservationOperator::BelowSurfaceGeopotential( double pressure ) const
{
	double	Rd    = GasConstantDryAir();
	double	tstar;
	double  alpha;
	double	phi;

	NXASSERT(( pressure <= m_p.back() ));
	tstar = m_Tstar;
	alpha = BelowSurfaceAlpha(&tstar, m_T0, m_phisurf);						// Calculate the (possibly adjusted) surface temperature and lapse rate/alpha coefficient
	phi   = m_phisurf -  (Rd*tstar/alpha)*( pow( pressure/m_psurf, alpha) - 1);		// equation 10 from Trenberth et al 1993, equation 1.10 from IFS observations Cy40r1 Part I
	return phi;
}

/*---------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::BelowSurfaceTemperature				*/
/** Calculates the temperature at a (higher) pressure level below the ECMWF surface. Uses 
 *	equation 1.15 to 1.18 of IFS Cy40r1 Part I: Observations. The below surface  code is copied from
 *	the IFS fortran in ppt_old.f90. Its a little opaque to follow but it seems to agree with the IFS documentation.
 *	I also tried the code in ppt.F90 and it makes no difference.
 *
 *	@param double geopotential
 *		The geopotential at which the temperature is required.
 *	@par RETURN
 *	Returns the below surface temperature at the requested geopotential.
 *	@par HISTORY
 * 	2016-5-5
 **/
/*------------------------------------------------------------------------*/

double nxMeteorology_EcwmfIFS_ObservationOperator::BelowSurfaceTemperature( double geopotential) const
{
	double  g         = LocalGravityAtSeaLevel();
//	double	Rd        = GasConstantDryAir();
	double  lapserate = BelowSurfaceLapseRate();
	double	alpha     = std::numeric_limits<double>::quiet_NaN();
	double	T         = std::numeric_limits<double>::quiet_NaN();
//	double	hsurf;
//	double	T0prime;
	double	deltaphi;

	NXASSERT(( geopotential <= m_phi.back() ));
	if ( geopotential > m_phisurf)
	{
		NXASSERT( (geopotential <= m_phi.back()));
		T  = m_Tstar + ( geopotential - m_phisurf)*( m_T.back()-m_Tstar)/(m_phi.back() - m_phisurf);
	}
	else
	{

		double TXX    = 298.0;
		double phi2000 = 2000.0*g;
		double phi2500 = 2500.0*g;
		double USDFI  = 1.0/(phi2500 - phi2000);
		double alpha0  = lapserate/g;
		double Tpl     = std::min(m_T0,TXX);
		double ZCOEFPL = std::min(1.0, std::max(0.0,(m_phisurf-phi2000)*USDFI));
		alpha  = std::max( 0.0 , alpha0 + ZCOEFPL*(Tpl - m_T0)/std::max(m_phisurf, phi2000) );

		//hsurf = m_phisurf/g;
		//if (hsurf < 2000.0)
		//{
		//	alpha = lapserate/g;
		//}
		//else 
		//{
		//	if ( hsurf > 2500.0)
		//	{
		//		T0prime = std::min( m_T0, 298.0);
		//	}
		//	else
		//	{
		//		T0prime = 0.002*( (2500.0 - hsurf)*m_T0 + (hsurf-2000.0)*std::min(m_T0, 298.0) );
		//	}			
		//	alpha = (T0prime-m_Tstar)/m_phisurf;
		//	if (T0prime < m_Tstar) alpha = 0.0;
		//}
		deltaphi = (m_phisurf - geopotential);
		T = m_Tstar + alpha*deltaphi;
	}
	NXASSERT(( NXFINITE(T)  ));
	return T;
}

/*-----------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::PressureToGeoPotential		 2016- 5- 30*/
/** The operation that finds the geopotential of a given pressure. This follows
 *	section 1.31. Geopotential Height of ECMWF IFS Documentation Cy40r1, Part I Observations. 
 *	This involves linearly interpolating in \f$\log p\f$ the \f$\Delta\phi_{k}\f$ generated
 *	from the difference in the ICAO and ECWMF model. Note that we have not implemented the quadratic
 *	interpolation for altitudes the second model level.
**/
/*---------------------------------------------------------------------------*/

bool nxMeteorology_EcwmfIFS_ObservationOperator::LogPressureToGeopotential( double lnp, double * phi ) const
{
	double	deltaphi;
	bool	ok;

	NXASSERT( (lnp < 20.0));
	deltaphi = nxLinearInterpolate::EvaluateYatX(	lnp, m_lnp, m_deltaphik, nxLinearInterpolate::ENUM_MISSINGVALUE, std::numeric_limits<double>::quiet_NaN() );
	*phi     =  deltaphi + ( m_icao.Phi_ICAO(exp(lnp))- m_phisurf_icao) + m_phisurf;
	ok = NXFINITE( *phi );
	NXASSERT((ok));
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::Interpolate_Temperature		 2016- 5- 31*/
/** Interpolates temperature between pressure levels in the model. This follows the first
 *	paragraph of section 1.3.4 Temperature from ECMWF IFS part I: Observations Cy40r1. 
 *	"Temperature is interpolated linearly in pressure (PPINTP)... ". The extrapolation outside
 *	the bounds of the model levels are done elsewhere.
**/
/*---------------------------------------------------------------------------*/

double nxMeteorology_EcwmfIFS_ObservationOperator::Interpolate_Temperature( double P0, double T0, double P1, double T1, double P) const
{
	double	Tvals[2] = { T0, T1 };
	double	z0     = P0;
	double  z1	   = P1;
	double	z      = P;
	double	T;

	T = nxLinearInterpolate::FromTwoPoints( z, z0, z1, Tvals);
	return T;
}

/*-----------------------------------------------------------------------------
 *					GeopotentialHeightZeroFinder		 2016- 5- 31*/
/** A small class used to find the pressure of a given geopotential height. This
 *	is the inverse operation  of method nxMeteorology_EcwmfIFS_ObservationOperator::LogPressureToGeopotential.
 *	It is applied as an iterative technique that adjusts the pressure between two two upper limits to 
 *	until the calculated geopotential equals the required geopotential. This class is used as the functor for
 *	a zero finding algorithm .
**/
/*---------------------------------------------------------------------------*/

class GeopotentialHeightZeroFinder
{
	private:
		double												m_targetphi;
		const nxMeteorology_EcwmfIFS_ObservationOperator&	m_ecmwfifs;

	public:
					GeopotentialHeightZeroFinder(  double targetphi,  const nxMeteorology_EcwmfIFS_ObservationOperator&	ecmwfifs) : m_ecmwfifs( ecmwfifs), m_targetphi(targetphi) {}
		double		operator()(double p) { double phi; m_ecmwfifs.LogPressureToGeopotential(p, &phi); return (phi-m_targetphi);}
};

/*-----------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::GeopotentialToPressure		 2016- 5- 31*/
/** Retrieve the pressure an temperature at a given geopotential for all altitudes 
 *	fully surrounded by the model levels.
 **/
/*---------------------------------------------------------------------------*/

bool nxMeteorology_EcwmfIFS_ObservationOperator::GeopotentialToPressureAndTemperature( double phi, double * p, double *T ) const
{
	bool	ok;
	size_t	idx0 = 0;
	size_t	idx1 = 0;
	double	lnp0;
	double	lnp1;
	double	phi0;
	double	phi1;
	int		status;

	ok = (phi > m_phi.back()) && (phi <= m_phi.front());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxMeteorology_EcwmfIFS_ObservationOperator::GeopotentialToPressure, requested geopotential %10.4f is out of ECMWF standard model range %10.4f to %10.4f", (double)phi, (double)m_phi.back(), (double)m_phi.front());
		*p = std::numeric_limits<double>::quiet_NaN();
		*T = std::numeric_limits<double>::quiet_NaN();
	}
	else
	{
		ok = ok && nxLinearInterpolate::FindBoundingIndicesAscending<double, std::vector<double>::const_reverse_iterator>( m_phi.rbegin(), m_phi.rend(), phi, &idx0, &idx1, &phi0, &phi1);
		if (ok)
		{
			idx0 = m_phi.size() - idx0 -1;
			idx1 = m_phi.size() - idx1 -1;

			lnp0 = std::min( m_lnp.at(idx0)+0.05, m_lnp.back()  );			// Set nominal bounds for pressure search
			lnp1 = std::max( m_lnp.at(idx1)-0.05, m_lnp.front() );			// a little bit bigger to avoid any round off issues near the edges, but dont go out of bounds
			ok   = zbrent2<GeopotentialHeightZeroFinder> ( GeopotentialHeightZeroFinder( phi, *this), lnp1, lnp0, 0.000001, p, &status, 15);

		#if defined(NXDEBUG)
			double phitest; 
			LogPressureToGeopotential( *p, &phitest ); 
			NXASSERT( (fabs(phitest-phi) < 1.0) );
		#endif
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"nxMeteorology_EcwmfIFS_ObservationOperator::GeopotentialToPressureAndTemperature, failed to find zero match for geopotential height");
			}
			else
			{
				*p = exp(*p);
				*T   = Interpolate_Temperature( m_p.at(idx0), m_T.at(idx0), m_p.at(idx1), m_T.at(idx1), *p);
			}
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"nxMeteorology_EcwmfIFS_ObservationOperator::GeopotentialToPressure, Cannot find a pressure for geopotential %e", (double)phi);
			*p = std::numeric_limits<double>::quiet_NaN();
			*T = std::numeric_limits<double>::quiet_NaN();
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::GeopotentialToBelowGroundPressureAndTemperature		 2016- 6- 1*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxMeteorology_EcwmfIFS_ObservationOperator::GeopotentialToBelowGroundPressureAndTemperature( double phi, double * p, double *T ) const
{
	bool	ok;

	ok = (phi <= m_phi.back());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxMeteorology_EcwmfIFS_ObservationOperator::GeopotentialToBelowGroundPressureAndTemperature, requested geopotential %10.4f is not below the model lowest bound %10.4f", (double)phi, (double)m_phi.back() );
		*p = std::numeric_limits<double>::quiet_NaN();
		*T = std::numeric_limits<double>::quiet_NaN();
	}
	else
	{
		*p     = BelowSurfacePressure   ( phi);
		*T     = BelowSurfaceTemperature( phi);
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxMeteorology_EcwmfIFS_ObservationOperator::UpdateHeightEntry		 2016- 5- 31*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxMeteorology_EcwmfIFS_ObservationOperator::UpdateHeightEntry( double heightm, nxMeteorology_EcwmfIFS_HeightEntry* entry ) const
{
	bool ok;

	ok = (entry->heightm == heightm) && entry->isvalid;
	if (!ok)
	{
		entry->heightm      = heightm;
		entry->geopotential = GeometricHeightToGeopotential( heightm);
		if ( (entry->geopotential > m_phi.back() ) &&  (entry->geopotential <= m_phi.front() ))
		{
			entry->isvalid  = GeopotentialToPressureAndTemperature( entry->geopotential, &entry->pressure_pa, &entry->temperature_K);
		}
		else if (entry->geopotential <= m_phi.back() )
		{
			entry->isvalid  = GeopotentialToBelowGroundPressureAndTemperature( entry->geopotential, &entry->pressure_pa, &entry->temperature_K);
		}
		else if (entry->geopotential >= m_phi.front())
		{
			entry->pressure_pa    = AboveModelPressure  ( entry->geopotential );
			entry->temperature_K  = AboveModelTemperature(  entry->geopotential );
			entry->isvalid        = true;
		}
		ok = entry->isvalid;
	}
	return ok;
}
