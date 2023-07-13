#include <skopticalproperties21.h>
#include <limits>
/*-----------------------------------------------------------------------------
 *					struct skSpectralLineShapeStorageBuffer_VoigtHumlik		2013-3-7*/
/** \ingroup spectrallineinternals
 *	[DEPRECATED] This code has been replaced by the Kuntz+errata algorithm.
 *	The storage buffer for each spectral line used in the Voigt Humlik
 *	calculation. 
 *
 *	The class implements a quick and accurate Voigt profile based on
 *	the HUMLIK.FOR code derived from the paper: 
 *	"Rapid Approximation to the Voigt/Faddeeva function and its derivatives",
 *	R. J Wells, JQSRT 62, 1999, 29-48

 **/
/*---------------------------------------------------------------------------*/

class skSpectralLineShapeStorageBuffer_VoigtHumlik
{
	public:
		double	m_nu00;			// The wave number of the spectral line
		double	m_aD;			// The doppler halh width
		double	m_aL;			// The Lorenz half width.

		double	Y;				// Current value of Y
		double	YQ;				// Y*Y
		double	YRRTPI;			// y/SQRT(pi)

		double	XLIM0;			// |x| on region boundaries
		double	XLIM1;
		double	XLIM2;
		double	XLIM3;
		double	XLIM4;

		double	A0;				// Region 1 y dependents
		double	D0;
		double	D2;

		double	H0;				// Region 2 Y-dependents
		double	H2;
		double	H4;
		double	H6;
		double	E0;
		double	E2;
		double	E4;

		double	Z0;				// Region 3 Y-dependents
		double	Z2;
		double	Z4;
		double	Z6;
		double	Z8;
		double	P0;
		double	P2;
		double	P4;
		double	P6;
		double	P8;

	private:
		bool			SetY				( double newY);
		double			NormalizedVoigt		( double X   );

	public:
						skSpectralLineShapeStorageBuffer_VoigtHumlik();

		double			Voigt				( double nu  );												// Get the Voigt profile
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
		bool			IsConfigured		( ) const					{ return (NXFINITE(Y));}	// See if the user has properly called SetLineParams
		double			LorentzHWHH			( ) const					{ return m_aL;}					// Get the Lorentz HWHH in wavenumbers
		double			DopplerHWHH			( ) const					{ return m_aD;}					// Get the Doppler HWHH in wavenumbers
};



using namespace::nxcgs;			// Use the CGS constant
static	const double  R0         = 146.7;											// Region boundaries
static	const double  R1         = 14.67;											// for R=4
static	const double  RRTPI      = 0.56418958354775628694807945156077;				// 1/SQRT(pi) 
static 	const double  Y0         = 1.5;												// for CPF12 algorithm
static 	const double  Y0PY0      = 3.0;												// Y0+Y0;								
static 	const double  Y0Q        = 2.25;											// Y0*Y0;
static  const double  twoln2     = 1.3862943611198906188344642429164;				// 2*log(2);
static  const double  sqrtln2	 = 0.8325546111576977563531646448952;				// SQRT( LN(2) )
static  const double  sqrtln2opi = 0.46971863934982566688617016420509;				// SQRT((log(2)/pi)=0.46971863934983
static  const double  oneAMU     = 1.66053886E-24;									// One atomic mass unit in grams 

/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtHumlik::skSpectralLineShapeStorageBuffer_VoigtHumlik		2013-3-7*/
/** Initialize the Voigt Humlik buffer so the first valid call will reset all
 *	of the interbally cached variables.
 **/
/*---------------------------------------------------------------------------*/

skSpectralLineShapeStorageBuffer_VoigtHumlik::skSpectralLineShapeStorageBuffer_VoigtHumlik()
{
	m_nu00 = 10000;					// default central wavenumber to 1 micron
 	m_aL   = 0.05;					// The Lorentz half width
    m_aD   = 0.01;					// The Doppler half width
	ResetLineParams();
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_HumlicekWells::SetY								2002-9-17*/
/**	Set the Y parameter of the Voigt function.  This basically configures
 *	a whole bunch of parameters that can be pre-calculated.
 *	Y is normally SQRT(LN(2))*LorentzGamma/DopplerGamma
 *	and hence is normally set by pressure, temperature and mass of species.
 **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShapeStorageBuffer_VoigtHumlik::SetY( double newY )
{
	if ( newY != Y)											// make sure we do need to update cache
	{																	// we do so
		Y  = newY;
		YQ = Y*Y;
		Y       = Y;
		YQ      = YQ;										// y^2
		YRRTPI  = Y*RRTPI;									// y/SQRT(pi)

		XLIM0 = R0 - Y; 
		XLIM1 = R1 - Y;
		XLIM3 = 3.097*Y - 0.45;
		XLIM2 = 6.8 - Y;
		XLIM4 = 18.1*Y + 1.65;
		if ( Y <= 0.000001 )													// When y<10^-6
		{
			XLIM1 = XLIM0;														// avoid W4 algorithm
			XLIM2 = XLIM0;
		}

		// ---- Region 1  y-dependents

		A0 = YQ + 0.5;                                             // Region 1 y-dependents
		D0 = A0*A0;
		D2 = YQ + YQ - 1.0;

		// ---- Region 2 y-dependents

		H0 =  0.5625 + YQ*(4.5 + YQ*(10.5 + YQ*(6.0 + YQ)));				// Region 2 y-dependents
		H2 = -4.5    + YQ*(9.0 + YQ*( 6.0 + YQ* 4.0));
		H4 = 10.5    - YQ*(6.0 - YQ*  6.0);
		H6 = -6.0    + YQ* 4.0;
		E0 =  1.875  + YQ*(8.25 + YQ*(5.5 + YQ));
		E2 =  5.25   + YQ*(1.0  + YQ* 3.0);
		E4 =  0.75*H6;

		// ---- region 3  y-dependents

		Z0 = 272.1014     + Y*(1280.829 + Y*(2802.870 + Y*(3764.966	+ Y*(3447.629 + Y*(2256.981 + Y*(1074.409 + Y*(369.1989 + Y*(88.26741 + Y*(13.39880 + Y)))))))));
		Z2 = 211.678      + Y*(902.3066 + Y*(1758.336 + Y*(2037.310 + Y*(1549.675 + Y*(793.4273 + Y*(266.2987 + Y*(53.59518 + Y*5.0)))))));
		Z4 = 78.86585     + Y*(308.1852 + Y*(497.3014 + Y*(479.2576 + Y*(269.2916 + Y*(80.39278 + Y*10.0)))));
		Z6 = 22.03523     + Y*(55.02933 + Y*(92.75679 + Y*(53.59518 + Y*10.0)));
		Z8 = 1.496460     + Y*(13.39880 + Y*5.0);

		P0 = 153.5168     + Y*(549.3954 + Y*(919.4955 + Y*(946.8970 + Y*(662.8097 + Y*(328.2151 + Y*(115.3772 + Y*(27.93941 + Y*(4.264678 + Y*0.3183291))))))));
		P2 = -34.16955    + Y*(-1.322256+ Y*(124.5975 + Y*(189.7730 + Y*(139.4665 + Y*(56.81652 + Y*(12.79458 + Y*1.2733163))))));
		P4 = 2.584042     + Y*(10.46332 + Y*(24.01655 + Y*(29.81482 + Y*(12.79568 + Y*1.9099744))));
		P6 = -0.07272979  + Y*(0.9377051+ Y*(4.266322 + Y*1.273316));
		P8 = 0.0005480304 + Y*0.3183291;
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtHumlik::ResetLineParams		2013-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

void skSpectralLineShapeStorageBuffer_VoigtHumlik::ResetLineParams()
{
	Y = std::numeric_limits<double>::quiet_NaN();				// Current value of Y
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_HumlicekWells::NormalizedVoigt		2013-3-6*/
/**	Calculate the Voigt profile at the specified Voigt Normalized Coordinate
 *	X.
 *
 *	X is normally SQRT(LN(2))*( NU-NU00)/DopplerGamma
 */
/*---------------------------------------------------------------------------*/

double skSpectralLineShapeStorageBuffer_VoigtHumlik::NormalizedVoigt( double X )
{
	static  const double C[6]   = {1.0117281, -0.75197147,  0.012557727, 0.010022008,  -0.00024206814,   0.00000050084806 };
	static  const double S[6]   = {1.393237,   0.23115241, -0.15535147,  0.0062183662,  0.000091908299, -0.00000062752596 };
	static  const double T[6]   = {0.31424038, 0.94778839,  1.5976826,   2.2795071,     3.0206370,       3.8897249 };

// Local variables

	int    J;																	// Loop variables
	double K;
	double ABX, XQ;																// |x|, x^2
	double XP[6], XM[6], YP[6], YM[6];											// CPF12 temporary values
	double MQ[6], PQ[6], MF[6], PF[6];
	double D, YF, YPY0, YPY0Q;  
	bool	ok;

	ok = (NXFINITE(Y) );
	if (!ok)
	{
		K = 0;
	}
	else
	{
		ABX = (double)fabs( X );																// |x|
		XQ  = ABX*ABX;																	// x^2
		if ( ABX > XLIM0 )																// Region 0 algorithm
		{
			K = YRRTPI / (XQ + YQ);
		}
		else if ( ABX > XLIM1 )															// Humlicek W4 Region 1
		{
			D = RRTPI / (D0 + XQ*(D2 + XQ));
			K = D*Y   * (A0 + XQ);
		}
		else if ( ABX > XLIM2 )															// Humlicek W4 Region 2 
		{
			D = RRTPI / (H0 + XQ*(H2 + XQ*(H4 + XQ*(H6 + XQ))));
			K = D*Y    *(E0 + XQ*(E2 + XQ*(E4 + XQ)));
		}

		else if ( ABX < XLIM3 )															// Humlicek W4 Region 3
		{
			D = 1.7724538 / (Z0 + XQ*(Z2 + XQ*(Z4 + XQ*(Z6 + XQ*(Z8+XQ)))));
			K = D*(P0 + XQ*(P2 + XQ*(P4 + XQ*(P6 + XQ*P8))));
		}
		else
		{																				// Humlicek CPF12 algorithm
				                                                         
			YPY0  = Y + Y0;
			YPY0Q = YPY0*YPY0;
			K     = 0.0;
			for (J=0; J < 6; J++)
			{
				D     = X - T[J];
				MQ[J] = D*D;
				MF[J] = 1.0 / (MQ[J] + YPY0Q);
				XM[J] = MF[J]*D;
				YM[J] = MF[J]*YPY0;
				D     = X + T[J];
				PQ[J] = D*D;
				PF[J] = 1.0 / (PQ[J] + YPY0Q);
				XP[J] = PF[J]*D;
				YP[J] = PF[J]*YPY0;
			}

			if ( ABX <= XLIM4 )															// Humlicek CPF12 Region I
			{
				for (J=0; J < 6; J++)
				{
					K += C[J]*(YM[J]+YP[J]) - S[J]*(XM[J]-XP[J]);
				}

			}
			else
			{				                                                            // Humlicek CPF12 Region II

				YF   = Y + Y0PY0;
				for (J=0; J < 6; J++)
				{
					K    +=   (C[J]*(MQ[J]*MF[J]-Y0*YM[J]) + S[J]*YF*XM[J]) / (MQ[J]+Y0Q)
							+ (C[J]*(PQ[J]*PF[J]-Y0*YP[J]) - S[J]*YF*XP[J]) / (PQ[J]+Y0Q);
				}
				K = Y*K + exp( -XQ );
			}
		}
	}
	return K;
}


/*---------------------------------------------------------------------------
 *'					VoigtProfile::Voigt                              2002-9-17
 *-------------------------------------------------------------------------*/

double skSpectralLineShapeStorageBuffer_VoigtHumlik::Voigt( double nu )
{

	double x;
	double voigt;

	x     = (nu - m_nu00)*sqrtln2/m_aD;
	voigt = (sqrtln2opi/m_aD)*NormalizedVoigt( x );
	return voigt;
}


/*---------------------------------------------------------------------------
 *'					VoigtProfile::SetLineParams                      2002-9-17
 *
 *	nu00              = lines central wavelength in wavenumbers
 *	pressure          = air atmospheric pressure in dynes, cgs 
 *	partialpressure   = partial pressure of the species in dynes 
 *	temperature       = temperature in Kelvin
 *	mass              = molecular mass in grams
 *	LorentzWidthMult  = Lorentz Half Width at 1 atmosphere.
 *
 *	Calculates the Half width Half max of the Doppler (m_aD) and the 
 *	Half width half max of the Lorentz width (m_aL)
 *
 *	Doppler HWHM from Page 92, Exploration of the Solar System by Infra-red
 *	remote sensing.  Note their formula gives full width Half Max and we
 *`	want Half width half Max.
 *
 *	Lorentz pressure broadening from "The Hitran Molecular Spectroscopic Database
 *	and HAWKS (Hitran Atmospheric Workstation): 1996 Edition, Appendix equation A12
 *
 *-------------------------------------------------------------------------*/

bool skSpectralLineShapeStorageBuffer_VoigtHumlik::SetLineParams( double nu00,
																  double pressure,
																  double partialpressure,
																  double temperature,
																  double tref,
																  double tempcoeff,
																  double mass,
																  double airhalfwidth,
																  double selfhalfwidth)
{
	bool		  ok;							 
	double		  newy;
	const double  pref    = 1013250.0;									// Reference pressure in Dynes

	m_nu00 = nu00;
	m_aD   = nu00*sqrt(twoln2*KBOLTZMAN*temperature/(mass))/CLIGHT;
	m_aL   = pow(tref/temperature,tempcoeff)*( (airhalfwidth*(pressure-partialpressure)/pref)+(selfhalfwidth*partialpressure/pref) );
	newy   = sqrtln2*m_aL/m_aD;
	ok     = SetY(newy);
	if (!ok)
	{
		ResetLineParams();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ExtractStorageBuffer		2013-3-7*/
/** "Typecast" the baseclass storage buffer object so m_cache points to the
 *	structure we want. Again it looks big and ugly but its actually quite
 *	easy. Its way nicer to debug this than  using void* as type info is
 *	still preserved.
 **/
/*---------------------------------------------------------------------------*/

static skSpectralLineShapeStorageBuffer_VoigtHumlik*	ExtractStorageBuffer	( skSpectralLineShapeStorageBuffer* storagebuffer)
{
	return skSpectralLineShapeStorageBuffer_Type< skSpectralLineShapeStorageBuffer_VoigtHumlik >::Cache( storagebuffer);
}

/*---------------------------------------------------------------------------
 *'					skSpectralLineShape_HumlicekWells::skSpectralLineShape_HumlicekWells                        2002-9-17
 *  To calculate the Faddeeva function with relative error less than 10^(-R).
 *  R0=1.51*EXP(1.144*R) and R1=1.60*EXP(0.554*R) can be set by the the user
 *  subject to the constraints 14.88<R0<460.4 and 4.85<R1<25.5
 *-------------------------------------------------------------------------*/

skSpectralLineShape_HumlicekWells::skSpectralLineShape_HumlicekWells()
{
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_HumlicekWells::~skSpectralLineShape_HumlicekWells		2013-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineShape_HumlicekWells::~skSpectralLineShape_HumlicekWells()
{
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_HumlicekWells::CreateStorageBuffer		2013-3-7*/
/** Create an instance of the stroage buffer. It looks big and ugly
 *	but its actually quite simple
 **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_HumlicekWells::CreateStorageBuffer( skSpectralLineShapeStorageBuffer**	storagebuffer )
{
	skSpectralLineShapeStorageBuffer_Type< skSpectralLineShapeStorageBuffer_VoigtHumlik >*		buffer;
	bool																						ok;

	buffer = new skSpectralLineShapeStorageBuffer_Type< skSpectralLineShapeStorageBuffer_VoigtHumlik >;			// Allocate the storage object
	ok = (buffer != NULL);																						// 
	if ( ok)
	{
		buffer->AddRef();
	}
	*storagebuffer = dynamic_cast<skSpectralLineShapeStorageBuffer*>( buffer );
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_HumlicekWells::ConfigureLineParameters		2013-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_HumlicekWells::ConfigureLineParameters( const skSpectralLine*						spectralline,
															   double								temperature,
															   double								pressure,
															   const GEODETIC_INSTANT&				geopt,
															   skClimatology*						/*atmosphericstate*/,
															   skSpectralLineShapeStorageBuffer*	storagebuffer )
{
	skSpectralLineShapeStorageBuffer_VoigtHumlik*	cache;
	bool			ok;
	double			tref;
	double			nu00;
	double			partialpressure;
	double			tempcoeff;
	double			mass;
	double			airhalfwidth;
	double			selfhalfwidth;
	const double	pref    = 1013250.0;									// Reference pressure in Dynes
	double			patm;

	cache = ExtractStorageBuffer( storagebuffer );					// set m_cache to point to storage buffer
	ok =       (cache != NULL);
//	ok = ok && atmosphericstate->GetParameter( SKCLIMATOLOGY_TEMPERATURE_K, geopt, &temperature, false );
//	ok = ok && atmosphericstate->GetParameter( SKCLIMATOLOGY_PRESSURE_PA,   geopt, &pressure,    false );

	if (ok)
	{
		tref            = spectralline->Tref();				// reference temperature for temp coeff
		nu00            = spectralline->Nu();				// Get the wave number of the spectral line
		tempcoeff       = spectralline->Nair();				// Get the temperature coeffiecient
		partialpressure = spectralline->ParentMolecule()->PartialPressure( geopt, pressure, temperature );
		mass            = spectralline->ParentMolecule()->MassAMU();
		airhalfwidth    = spectralline->GammaAir();
		selfhalfwidth   = spectralline->GammaSelf();
		pressure        = pressure*10.0;			// Convert Pascals to Dynes
		partialpressure = partialpressure*10.0;		// Convert Pascals to Dynes.
		mass            = mass*oneAMU;				// Convert the AMU to grams
		patm            = pressure/pref;
		nu00            += spectralline->Deltaair()*patm;
		ok = cache->SetLineParams(  nu00, pressure, partialpressure, temperature, tref, tempcoeff, mass, airhalfwidth, selfhalfwidth );
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineShape_HumlicekWells::ConfigureLineParameters, There was an error configuring the line parameters for the Voigt Humlik code");
		cache->ResetLineParams();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_HumlicekWells::LineShapeFunction		2013-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_HumlicekWells::LineShapeFunction( double nu, double* uservalue, const skSpectralLine* /*spectralline*/, double	/*maxlinestrength*/, skSpectralLineShapeStorageBuffer* storagebuffer ) const
{
	skSpectralLineShapeStorageBuffer_VoigtHumlik*	cache;
	bool	ok;

	cache      =       ExtractStorageBuffer( storagebuffer );					// set m_cache to point to storage buffer
	ok         =       (cache != NULL);
	ok         = ok && cache->IsConfigured();
	*uservalue = ok ?  cache->Voigt(nu) : 0.0;
	return ok;
}
