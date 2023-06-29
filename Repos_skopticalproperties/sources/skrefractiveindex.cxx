#include <skopticalproperties21.h>


/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndex_MoistAir::operator=		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

void  skRTRefractiveIndex_Tabulated::operator=( const skRTRefractiveIndex_Tabulated& other )
{
	m_wavelen = other.m_wavelen;
	m_real    = other.m_real;
	m_imag    = other.m_imag;					// Imaginary part of Refractive index
}


/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_Tabulated::AttachToStatic		2003-10-10
 *	Attach to a linear array of refractive index data where tge data are org
 *-------------------------------------------------------------------------*/

bool skRTRefractiveIndex_Tabulated::AttachToStatic( double* data, size_t nelements)
{
	bool		ok;
	size_t		strides = 3*sizeof(*data);

	ok = (nelements >0 ) && (nelements%3 == 0);
	if (ok)
	{
		nelements = nelements/3;
		ok = ok && m_wavelen.nxArrayLinear<double>::Attach( 1, &nelements, data, NULL, &strides );
		ok = ok && m_real.nxArrayLinear<double>::Attach( 1, &nelements, data+1, NULL, &strides );
		ok = ok && m_imag.nxArrayLinear<double>::Attach( 1, &nelements, data+2, NULL, &strides );
	}
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING, "skRTRefractiveIndex_Tabulated::AttachToStatic, Error attaching arrays to specified array");
	}
	return ok;
}



/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::skRTRefractiveIndex_MoistAir 2003-10-23
 *	This code follows the upto date evaluation of the refractive index of air
 *	given by Ciddor.
 *
 *	Reference 1
 *	------------
 *	Refractive Index of Air: new equations for the visible and near infrared.
 *  Philip. E. Ciddor, Appl. Opt., 35, 9, 1566-1575, 1996.
 *-------------------------------------------------------------------------*/

skRTRefractiveIndex_MoistAir::skRTRefractiveIndex_MoistAir()
{
	m_T            = 293.15;			// default is 20 C
	m_P            = 101325.0;			// 1 atmosphere of pressure
	m_WaterVapPres = 0.0;				// Water Vapour Partial Pressure in pascals
	m_ppmCO2       = 450.0;				// Default value of CO2 density in dry air in parts per million

	UpdateStandardDensities();
}


/*-----------------------------------------------------------------------------
 *					skRTRefractiveIndex_MoistAir::operator=		2008-3-4*/
/** **/
/*---------------------------------------------------------------------------*/

bool skRTRefractiveIndex_MoistAir::DeepCopy( const skRTRefractiveIndex_MoistAir& other )
{
	m_T            = other.m_T;
	m_P            = other.m_P;
	m_WaterVapPres = other.m_WaterVapPres;
	m_ppmCO2       = other.m_ppmCO2;
	m_rho_axs      = other.m_rho_axs;
	m_rho_ws       = other.m_rho_ws;
	return true;
}


/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::UpdateStandardDensities		2003-10-24
 *	Updates the denisty of dry air at 15C 101325 Pa and density of pure
 *	water vapour at 15C 1333Pa using directions in paragraph following eqn 5 
  *	of P.E. Ciddor, Applied Optics, Vol. 35, No. 9., 1566-1573, 1996.
*-------------------------------------------------------------------------*/

void skRTRefractiveIndex_MoistAir::UpdateStandardDensities( )
{
	double T = m_T;							// Store the value of T and P
	double P = m_P;

	m_T = (273.15 + 20.0);					// Temperature for compressibility of pure waterin kelvins
	m_P = 1333.0;							// Total Pressure in pascals
	m_rho_ws = MoistAirDensity(1.0);		// density of pure water vapour at 20C, 1333 Pa.
 

	m_T = (273.15 + 15.0);					// Temperature for compressibility of dry air in Kelvins
	m_P = 101325.0;							// Total Pressure in pascals
	m_rho_axs = MoistAirDensity	( 0.0);		// Density of dry air at 15C, 101325 Pa, (CO2 present but no water)

	m_T = T;
	m_P = P;
}

/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::MolarMassDryAir		2003-10-24
 *	Calculates the Molar mass of dry air in kg/mol.  Given just below equation 4
 *	of P.E. Ciddor, Applied Optics, Vol. 35, No. 9., 1566-1573, 1996.
 *
 *-------------------------------------------------------------------------*/

double skRTRefractiveIndex_MoistAir::MolarMassDryAir() const
{
	double Ma;

	Ma  = (0.0289635 + 12.011E-9*(m_ppmCO2 - 400.0));		// Molar Mass of dry air with m_xc ppm of CO2 kg/mol
	return Ma;
}

/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::EnhancementFactor		2003-10-24
 *	Calculates the enhancement factor of water vapour in air using equation
 *	in lower middle of paragraph folloowing equation 4 
 *	of P.E. Ciddor, Applied Optics, Vol. 35, No. 9., 1566-1573, 1996..
 *-------------------------------------------------------------------------*/

double skRTRefractiveIndex_MoistAir::EnhancementFactor() const
{
	static const double ALPHA = 1.00062;
	static const double BETA  = 3.14E-8;
	static const double GAMMA = 5.6E-7;
	double  t                 = m_T - 273.15;

	return (ALPHA + BETA*m_P + GAMMA*t*t);				// enhancement factor of water vapour in air
}


/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::SaturationVapourPressure		2003-10-24
 *	Calculates the saturation vapour pressure of water vapour in air using equation
 *	in middle of paragraph following equation 4. Returns the answer in Pascals
 *	of P.E. Ciddor, Applied Optics, Vol. 35, No. 9., 1566-1573, 1996..
 *-------------------------------------------------------------------------*/

double skRTRefractiveIndex_MoistAir::SaturationVapourPressure() const
{
	static const double A	  =  1.2378847e-5;
	static const double B	  = -1.9121316e-2;
	static const double C	  = 33.93711047;
	static const double D	  = -6.3431645e3;
	
	return exp( (A*m_T +B)*m_T +C + D/m_T);			// Water saturation Vapour Pressure in Pascals (over liquid)
}

/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::MolarFractionWaterVapour		2003-10-24
 *	Calculates the Molar Fraction of water vapour using equation 4
 *	of P.E. Ciddor, Applied Optics, Vol. 35, No. 9., 1566-1573, 1996. MKS units.
 *-------------------------------------------------------------------------*/

double skRTRefractiveIndex_MoistAir::MolarFractionWaterVapour() const
{
	return EnhancementFactor()*m_WaterVapPres/m_P;					// Molar fraction of water = Enhancement*Partial Pressure Water/Total Pressure
}

/*---------------------------------------------------------------------------
 *'					MoistAirDensity		2003-10-23
 *	This subroutine calculates the density of moist air following equation 4
 *	of P.E. Ciddor, Applied Optics, Vol. 35, No. 9., 1566-1573, 1996. MKS units.
 *
 *	This follows the BIPM 1981/91 equation
 *
 *	Calculates density of dry air component in kg/m3
 *  and density of water vapor component
 *  m_T = temperature in Kelvin
 *  m_P = total pressure
 *  m_WaterVapourPressure = partial pressure of water vapor
 *  Fractional Humidity   = PW/SVP
 *-------------------------------------------------------------------------*/

double skRTRefractiveIndex_MoistAir::MoistAirDensity( double xw, double* rhoair, double* rhowater ) const
{
	static const double R     = 8.314510;					// Gas constant J.mol-1.K-1
	static const double	Mw    = 0.0180150;					// Molar mass of Water kg/mol
	double	density;
	double  zrt;
	double  rhoa;
	double  rhow;
	double  Ma;
	double  Z;

	Ma       = MolarMassDryAir();							// Molar Mass of dry air with m_xc ppm of CO2 kg/mol
	Z        = Compressibility(xw);							// Compressibility of dry air and pure water vapour.
	zrt      = m_P/( Z*R*m_T);
	rhoa     = zrt*Ma*(1-xw);
	rhow     = zrt*Mw*xw;
	density  = rhoa+rhow;
	if (rhoair   != NULL) *rhoair   = rhoa;
	if (rhowater != NULL) *rhowater = rhow;
	return density;
}

/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::Compressibility		2003-10-23
 *	This subroutine calculates the compressibility of moist air
 *	as described by P.E. Ciddor, Applied Optics, Vol. 35,
 *	No. 9., March 20, 1996. MKS units.
 *
 * T  = temperature in Kelvin
 * P  = total pressure (Pascals)
 * Xw = molar fraction of water vapor
 *-------------------------------------------------------------------------*/

double skRTRefractiveIndex_MoistAir::Compressibility( double Xw) const
{
	static const double a0 =  1.58123E-6;
	static const double a1 = -2.9331E-8;
	static const double a2 =  1.1043E-10;
	static const double b0  =  5.707E-6;
	static const double b1 = -2.051E-8;
	static const double c0  =  1.9898E-4;
	static const double c1 = -2.376E-6;
	static const double d  =  1.83E-11;
	static const double e  = -0.765E-8;

	double t = m_T - 273.15;
	double f = m_P/m_T;
	double Z;
	double Xw2 = Xw*Xw;

	Z = 1 - f*(a0 + a1*t + a2*t*t + (b0 + b1*t)*Xw + (c0 + c1*t)*Xw2) + f*f*(d + e*Xw2);
	return Z;
}


/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::DryAirRefractivityAtSTP		2003-10-23
 *	Get the refractivity of dry air at 15C, 101325 Pa.
 *	This includes contributions from CO2. See equation 1 
 *	of P.E. Ciddor, Applied Optics, Vol. 35, No. 9., 1566-1573, 1996. MKS units.
 *-------------------------------------------------------------------------*/

double skRTRefractiveIndex_MoistAir::DryAirRefractivityAtSTP( double vacuumwavenum ) const
{
	static const double k0 = 238.0185;			// micron m-2
	static const double k1 = 5792105.0;			// micron m-2
	static const double k2 = 57.362;			// micron m-2
	static const double k3 = 167917.0;			// micron m-2;
	double  sigma  = vacuumwavenum*1.0E-04;		// convert cm-1 to micron-1
	double	sigma2 = sigma*sigma;
	double	nasm1;
	double  naxs;


	nasm1 = 1.0E-8*( k1/( k0 - sigma2) + k3/(k2-sigma2) );			// Equation 1 of reference 1
	naxs  =  nasm1*( 1.0 + 0.534E-06*(m_ppmCO2 - 450.0) );			// Correction for CO2 density
	return naxs;
};

/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::WaterVapourRefractivityAtSTP		2003-10-24
 *	Get the refractivity of watervapour at 20C, 1333 Pa. See equation 3 
 *	of P.E. Ciddor, Applied Optics, Vol. 35, No. 9., 1566-1573, 1996. MKS units.
 *-------------------------------------------------------------------------*/

double skRTRefractiveIndex_MoistAir::WaterVapourRefractivityAtSTP( double vacuumwavenum ) const
{
	static const double w0 = 294.235;			// micron m2
	static const double w1 = 2.6422;			// micron m2
	static const double w2 = -0.032380;			// micron m2
	static const double w3 = 0.004028;			// micron m2;
	double  sigma  = vacuumwavenum*1.0E-04;		// convert cm-1 to micron-1
	double	sigma2 = sigma*sigma;
	double  nws;

	nws = 1.022E-8*( w0 + (w1 + (w2 + w3*sigma2)*sigma2)*sigma2 );
	return nws;
};


/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::Refractivity		2003-10-24
 *	Get the refractivity of moist air at curenttemperature and pressure.
 *	Uses equation 5 which is almost as accurate as equations 6,7 8 of
 *	of P.E. Ciddor, Applied Optics, Vol. 35, No. 9., 1566-1573, 1996. MKS units.
 *-------------------------------------------------------------------------*/

std::complex<double> skRTRefractiveIndex_MoistAir::Refractivity( double wavenum ) const
{
	double	Nws;		// Refractivity of water vapour at 20C 1333 Pa
	double	Naxs;		// Refractivity of dry air 	15C 101325 Pa.
	double	rhoa;		// density of air at current pressure and temperature
	double  rhow;		// density of water vapour at current pressure and temperature
	double  Na;			// Refractivity of dry air;
	double	Nw;			// Refractivity of water vapour
	double  xw;

	xw   = MolarFractionWaterVapour    ();
	Nws  = WaterVapourRefractivityAtSTP( wavenum );		// get refractivity of water vapour at standard temp and pressure
	Naxs = DryAirRefractivityAtSTP     ( wavenum );		// get refractivity of dry air at standard temp and pressure
	MoistAirDensity( xw, &rhoa, &rhow );				// Get the density of dry air and the density of water vapour
	Na   = rhoa/m_rho_axs*Naxs;							// refractivity of dry air at current temperature and pressure
	Nw   = rhow/m_rho_ws*Nws;							// refractivity of water vapour at current temperature and pressure
	return Na+Nw;
};


/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::Set_CO2ppm		2003-10-24
 *	set the CO2 content of the air in parts per million. The default value
 *	is 450.
 *-------------------------------------------------------------------------*/

bool skRTRefractiveIndex_MoistAir::Set_CO2ppm( double co2ppm)
{
	m_ppmCO2 = co2ppm;
	UpdateStandardDensities();
	return true;
}

/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::Set_TotalPressure		2003-10-24
 *	Set the total pressure (Pascals), default is 101325.0
 *-------------------------------------------------------------------------*/

bool skRTRefractiveIndex_MoistAir::Set_TotalPressure( double totalpressure)
{
	m_P = totalpressure;
	return true;
}
/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::Set_Temperature		2003-10-24
 *	Set the temperature (Kelvins) , default is 293.15
 *-------------------------------------------------------------------------*/

bool skRTRefractiveIndex_MoistAir::Set_Temperature( double kelvins)
{
	bool	ok;
	m_T = kelvins;
	ok = (m_T > 30);
	if (!ok) nxLog::Record( NXLOG_WARNING, "skRTRefractiveIndex_MoistAir::Set_PTW, The temperature entered (%g) looks like it is in Centigrade. You must use kelvins", (double)m_T);
	return ok;
}

/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::Set_WaterVapourPressure 2003-10-24
 *	Set the water vapour pressure (Pascals).  Default is 0.0
 *-------------------------------------------------------------------------*/

bool skRTRefractiveIndex_MoistAir::Set_WaterVapourPressure( double wvp)
{
	m_WaterVapPres = wvp;
	return true;
}


/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_MoistAir::Set_PTandRH		2003-10-24
 *-------------------------------------------------------------------------*/

bool skRTRefractiveIndex_MoistAir::Set_PTandRH( double totalpressure, double kelvins, double fractionalhumidity)
{
	Set_TotalPressure( totalpressure );
	Set_Temperature  ( kelvins );
	m_WaterVapPres = fractionalhumidity*SaturationVapourPressure();
	return true;
}


/*---------------------------------------------------------------------------
 *'					skRTRefractiveIndex_Tabulated::RefractiveIndex		2003-10-10
 *-------------------------------------------------------------------------*/

std::complex<double> skRTRefractiveIndex_Tabulated::RefractiveIndex( double wavenum ) const
{
	std::complex<double>		answer;
	double						real_ri;
	double						imag_ri;
	nxArrayIter<double>			iter;
	intptr_t					idx2;
	intptr_t					idx1;
	intptr_t					npts   = (intptr_t)m_wavelen.size();
	nxArrayIter<double>			wbegin(m_wavelen.begin());
	nxArrayIter<double>			wend (m_wavelen.end());
	double						nanometer;
	double						w1;
	double						g;
	double						gm1;

	NXASSERT(( wavenum != 0.0 ));
	nanometer = 1.0E7/wavenum;
	iter = std::lower_bound( wbegin, wend, nanometer );
	idx2  = (iter - wbegin);
	idx1  = idx2-1;
	if ( (idx2 <= 0) || ( idx2 >= npts))
	{
		if (idx2 >= npts) idx2 = npts-1;
		if (idx2  < 0)    idx2 = 0;
		real_ri = m_real[idx2];
		imag_ri = m_imag[idx2];
	}
	else
	{
		w1      = m_wavelen[idx1];
		g       = (nanometer - w1)/(m_wavelen[idx2] - w1);
		gm1     = 1 - g;
		real_ri = gm1*m_real[idx1] + g*m_real[idx2];
		imag_ri = gm1*m_imag[idx1] + g*m_imag[idx2];
	}
	return std::complex<double>(real_ri, imag_ri);
}
