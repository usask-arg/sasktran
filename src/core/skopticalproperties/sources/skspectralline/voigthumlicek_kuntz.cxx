#include <skopticalproperties21.h>
#include <limits>




using namespace::nxcgs;			// Use the CGS constant

static  const double  twoln2     = 1.3862943611198906188344642429164;				// 2*log(2);
static  const double  sqrtln2	 = 0.8325546111576977563531646448952;				// SQRT( LN(2) )
static  const double  sqrtln2opi = 0.46971863934982566688617016420509;				// SQRT((log(2)/pi)=0.46971863934983
//static  const double  oneAMU     = 1.66053886E-24;									// One atomic mass unit in grams 


/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtKuntz::skSpectralLineShape_VoigtKuntz		2014-2-18*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineShape_VoigtKuntz::skSpectralLineShape_VoigtKuntz		( )
{
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtKuntz::~skSpectralLineShape_VoigtKuntz		2014-2-18*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineShape_VoigtKuntz::~skSpectralLineShape_VoigtKuntz		( )
{
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtKuntz::LineShapeFunction		2014-2-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtKuntz::LineShapeFunction( double nu, double* uservalue, const skSpectralLine* spectralline, skSpectralLineShapeStorageBuffer* storagebuffer )
{
	skSpectralLineShapeStorageBuffer_VoigtKuntz*	cache;

	cache = dynamic_cast < skSpectralLineShapeStorageBuffer_VoigtKuntz*>( storagebuffer );					// set m_cache to point to storage buffer
	*uservalue = cache->Voigt( nu, spectralline->LineIntensity() );
	return true;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtKuntz::AddLineShapeFunctionArray		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtKuntz::AddLineShapeFunctionArray( const std::vector<double>& nu, std::vector<double>* uservalue, const skSpectralLine* spectralline, skSpectralLineShapeStorageBuffer* storagebuffer )
{
	skSpectralLineShapeStorageBuffer_VoigtKuntz*	cache;

	cache = dynamic_cast < skSpectralLineShapeStorageBuffer_VoigtKuntz*>( storagebuffer );					// set m_cache to point to storage buffer
	return cache->AddVoigt( nu, uservalue, spectralline->LineIntensity() );
}



/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtKuntz::SetParentMaxLineStrength		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtKuntz::SetParentMaxLineStrength( double parentmaxlinestrength, const skSpectralLine* spectralline, skSpectralLineShapeStorageBuffer* storagebuffer ) const
{
	skSpectralLineShapeStorageBuffer_VoigtKuntz*	cache;

	cache = dynamic_cast < skSpectralLineShapeStorageBuffer_VoigtKuntz*>( storagebuffer );
	cache->SetLimitsFromParentMaxLineStrength( parentmaxlinestrength, spectralline->LineIntensity() );
	return true;
}

bool skSpectralLineShape_VoigtKuntz::SetTolerance( double tolerance, const skSpectralLine* /*spectralline*/, skSpectralLineShapeStorageBuffer* storagebuffer ) 
{
	skSpectralLineShapeStorageBuffer_VoigtKuntz*	cache;

	cache = dynamic_cast < skSpectralLineShapeStorageBuffer_VoigtKuntz*>( storagebuffer );
	cache->SetTolerance( tolerance );
	return true;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtKuntz::ConfigureLineParameters		2014-2-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtKuntz::ConfigureLineParameters( const skSpectralLine*						spectralline,
															   double								temperature,
															   double								pressure,
															   const GEODETIC_INSTANT&				geopt,
															   skClimatology*						/*atmosphericstate*/,
															   skSpectralLineShapeStorageBuffer*	storagebuffer )
{
	skSpectralLineShapeStorageBuffer_VoigtKuntz*	cache;
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

	cache = dynamic_cast < skSpectralLineShapeStorageBuffer_VoigtKuntz*>( storagebuffer );					// set m_cache to point to storage buffer
	ok =       (cache != NULL);
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
		mass            = mass*AMU;					// Convert the AMU to grams
		patm            = pressure/pref;
		nu00            += spectralline->Deltaair()*patm;
		ok = cache->SetLineParams(  nu00, pressure, partialpressure, temperature, tref, tempcoeff, mass, airhalfwidth, selfhalfwidth );
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineShape_VoigtKuntz::ConfigureLineParameters, There was an error configuring the line parameters for the Voigt Kuntz code");
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtKuntz::CreateStorageBuffer		2014-2-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape_VoigtKuntz::CreateStorageBuffer( skSpectralLineShapeStorageBuffer** storagebuffer )
{
	skSpectralLineShapeStorageBuffer_VoigtKuntz*	ptr;

	ptr = new skSpectralLineShapeStorageBuffer_VoigtKuntz;
	*storagebuffer = ptr;
	ptr->AddRef();
	return true;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz		2014-2-18*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineShapeStorageBuffer_VoigtKuntz::skSpectralLineShapeStorageBuffer_VoigtKuntz()
{
	m_tolerance = 1.0E-09;			// Dont bother with Voigt profiles that contribute that contribute less than this amount of signal
	m_y         = -9.0E20;
	ResetCoefficients();
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::Newy1		2014-2-13*/
/** **/
/*---------------------------------------------------------------------------*/

void skSpectralLineShapeStorageBuffer_VoigtKuntz::ConfigureRegion1()
{
    m_yy = m_y*m_y;
    a1   = m_y*(0.2820948 + 0.5641896*m_yy);
    b1   = 0.5641896*m_y;
//	a2   =   0.25 + m_yy + m_yy*m_yy;			// Kuntz original paper
	a2   =   0.5 + m_yy + m_yy*m_yy;			// Proper form with errata correction from Ruyten
    b2   = - 1.0 + 2.0*m_yy ;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::Humlicek1		2014-2-13*/
/** **/
/*---------------------------------------------------------------------------*/

double skSpectralLineShapeStorageBuffer_VoigtKuntz::K1( double x ) 
{
	double xx;
	double	voigt;

	if ( a2 == 0.0 ) ConfigureRegion1();
	xx = x*x;
    voigt = (a1 + b1*xx)/(a2 + xx*(b2 + xx));
	return voigt;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::Newy2		2014-2-13*/
/** **/
/*---------------------------------------------------------------------------*/

void skSpectralLineShapeStorageBuffer_VoigtKuntz::ConfigureRegion2()
{
	m_yy = m_y*m_y;
	a3 =  m_y*(1.05786  +m_yy*(4.65456 + m_yy*(3.10304 + m_yy*0.56419)));
	b3 =  m_y*(2.9620   +m_yy*(0.56419 + m_yy*1.69257));
//	c3 =  m_y*(1.69257-m_yy*2.53885);	// Kuntz original paper
	c3 =  m_y*(-2.53885 +m_yy*1.69257);		// Proper form with errata correction from Ruyt	d3 =  m_y*0.5641896;
	d3 =  m_y*0.56419;
	
	a4 = 0.5625 +  m_yy*(4.5  + m_yy*(10.5 + m_yy*(6.0 + m_yy)));
	b4 = -4.5   +  m_yy*(9.0  + m_yy*(6.0  + m_yy*4.0)) ;
	c4 = 10.5   +  m_yy*(-6.0 + 6.0*m_yy );
	d4 = -6.0   +  m_yy*4.0 ;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz:Humlicek2		2014-2-13*/
/** **/
/*---------------------------------------------------------------------------*/

double skSpectralLineShapeStorageBuffer_VoigtKuntz::K2( double x) 
{
      double	xx;
	  double	voigt;

	  if ( a4 == 0.0 ) ConfigureRegion2();
      xx = x*x;
      voigt = ( a3 + xx*(b3 + xx*(c3 + xx*d3))) / ( a4 + xx*(b4 + xx*(c4 + xx*(d4+xx))));
      return voigt;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::Newy3		2014-2-13*/
/** **/
/*---------------------------------------------------------------------------*/

void skSpectralLineShapeStorageBuffer_VoigtKuntz::ConfigureRegion3( )
{

	a5 = 272.102     + m_y*(973.778  +m_y*(1629.76 +m_y*(1678.33 +m_y*(1174.8  +m_y*(581.746 +m_y*(204.501 +m_y*(49.5213 +m_y*(7.55895 +m_y*0.564224))))))));	// Corrected to follow errata
    b5 = -60.5644    + m_y*(-2.34403 +m_y*(220.843 +m_y*(336.364 +m_y*(247.198 +m_y*(100.705 +m_y*(22.6778 +m_y*2.25689))))));
    c5 = 4.58029     + m_y*(18.546   +m_y*(42.5683 +m_y*(52.8454 +m_y*(22.6798 +m_y*3.38534))));
    d5 = -0.128922   + m_y*(1.66203  +m_y*(7.56186 +m_y*2.25689)) ;
	e5 = 0.000971457 + 0.564224*m_y;

    a6  = 272.102 +m_y*(1280.83 +m_y*(2802.87 +m_y*(3764.97 +m_y*(3447.63 +m_y*(2256.98 +m_y*(1074.41 +m_y*(369.199 +m_y*(88.2674 +m_y*(13.3988 + m_y)))))))));
    b6  = 211.678 +m_y*(902.306 +m_y*(1758.34 +m_y*(2037.31 +m_y*(1549.68 +m_y*(793.427 +m_y*(266.299 +m_y*(53.5952 +5.0*m_y)))))));
    c6  = 78.866  +m_y*(308.186 +m_y*(497.302 +m_y*(479.258 +m_y*(269.292 +m_y*(80.3928 +m_y*10.0))))); 
    d6  = 22.0353 +m_y*(55.0293 +m_y*(92.7568 +m_y*(53.5952 +m_y*10.0)));		// Corrected to follow errata
    e6  = 1.49645 +m_y*(13.3988 +m_y*5.0);
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz:Humlicek3		2014-2-13*/
/** **/
/*---------------------------------------------------------------------------*/

double skSpectralLineShapeStorageBuffer_VoigtKuntz::K3( double x) 
{
	double	xx;
	double	voigt;

	if (a5 == 0.0) ConfigureRegion3();
	xx    = x*x;
	voigt = (a5 + xx*(b5+xx*(c5+xx*(d5 + xx*e5))))/(a6+xx*(b6+ xx*(c6+xx*(d6+xx*(e6+xx)))));
	return voigt;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::Newy4		2014-2-13*/
/** **/
/*---------------------------------------------------------------------------*/

void skSpectralLineShapeStorageBuffer_VoigtKuntz::ConfigureRegion4( )
{
	m_yy = m_y*m_y;
	m_y2 = 2.0*m_y;
	a7 = m_y*( 1.16028e9 +m_yy*(-9.86604e8 +m_yy*( 4.56662e8 +m_yy*(-1.53575e8 +m_yy*( 4.08168e7 +m_yy*(-9.69463e6 +m_yy*( 1.6841e6  +m_yy*(-320772.0 +m_yy*( 40649.2 +m_yy*(-5860.68  +m_yy*( 571.687 +m_yy*(-72.9359 +m_yy*(2.35944-m_yy*0.56419)))))))))))));
	b7 = m_y*(-5.60505e8 +m_yy*(-9.85386e8 +m_yy*( 8.06985e8 +m_yy*(-2.91876e8 +m_yy*( 8.64829e7 +m_yy*(-7.72359e6 +m_yy*( 3.59915e6 +m_yy*(-234417.0 +m_yy*( 45251.3 +m_yy*(-2269.19  +m_yy*(-234.143 +m_yy*( 23.0312 -m_yy*7.33447))))))))))));
	c7 = m_y*(-6.51523e8 +m_yy*( 2.47157e8 +m_yy*( 2.94262e8 +m_yy*(-2.04467e8 +m_yy*( 2.29303e7 +m_yy*(-2.38180e7 +m_yy*( 576054.0  +m_yy*(  98079.1 +m_yy*(-25338.3 +m_yy*( 1097.77  +m_yy*(97.6203  -m_yy*44.0068)))))))))));
	d7 = m_y*(-2.63894e8 +m_yy*( 2.70167e8 +m_yy*(-9.96224e7 +m_yy*(-4.15013e7 +m_yy*( 3.83112e7 +m_yy*( 2.24040e6 +m_yy*(-303569.0  +m_yy*( -66431.2 +m_yy*( 8381.97 +m_yy*(  228.563 -m_yy*161.358))))))))));
	e7 = m_y*(-6.31771e7 +m_yy*( 1.40677e8 +m_yy*( 5.56965e6 +m_yy*( 2.46201e7 +m_yy*( 468142.0  +m_yy*(-1.00300e6 +m_yy*(-66212.1   +m_yy*(  23507.6 +m_yy*(  296.38 -m_yy*403.396)))))))));
	f7 = m_y*(-1.69846e7 +m_yy*( 4.07382e6 +m_yy*(-3.32896e7 +m_yy*(-1.93114e6 +m_yy*(-934717.0  +m_yy*( 8820.94   +m_yy*( 37544.8   +m_yy*(  125.591 -m_yy*726.112))))))));
	g7 = m_y*(-1.23165e6 +m_yy*( 7.52883e6 +m_yy*(-900010.0  +m_yy*(-186682.0  +m_yy*( 79902.5   +m_yy*( 37371.9   +m_yy*(-260.198   -m_yy*968.15)))))));
	h7 = m_y*(-610622.   +m_yy*( 86407.6   +m_yy*( 153468.   +m_yy*( 72520.9   +m_yy*( 23137.1   +m_yy*(-571.645   -m_yy*968.15))))));
	o7 = m_y*(-23586.5   +m_yy*( 49883.8   +m_yy*( 26538.5   +m_yy*( 8073.15   +m_yy*(-575.164   -m_yy*726.112)))));
	p7 = m_y*(-8009.1    +m_yy*( 2198.86   +m_yy*( 953.655   +m_yy*(-352.467   -m_yy*403.396))));
	q7 = m_y*(-622.056   +m_yy*(-271.202   +m_yy*(-134.792   -m_yy*161.358)));
	r7 = m_y*(-77.0535   +m_yy*(-29.7896   -m_yy*44.0068));
	s7 = m_y*(-2.92264   -m_yy*7.33447);
	t7 = -0.56419*m_y;

	a8 = 1.02827e9 +m_yy*(-1.5599e9  +m_yy*( 1.17022e9 +m_yy*(-5.79099e8 +m_yy*( 2.11107e8 +m_yy*(-6.11148e7 +m_yy*( 1.44647e7 +m_yy*(-2.85721e6 +m_yy*( 483737.0 +m_yy*(-70946.1 +m_yy*( 9504.65 +m_yy*(-955.194 +m_yy*(126.532  +m_yy*(-3.68288+m_yy)))))))))))));
	b8 = 1.5599e9  +m_yy*(-2.28855e9 +m_yy*( 1.66421e9 +m_yy*(-7.53828e8 +m_yy*( 2.89676e8 +m_yy*(-7.01358e7 +m_yy*( 1.39465e7 +m_yy*(-2.84954e6 +m_yy*( 498334.0 +m_yy*(-55600.0 +m_yy*( 3058.26 +m_yy*( 533.254 +m_yy*(-40.5117 +m_yy*14.))))))))))));
	c8 = 1.17022e9 +m_yy*(-1.66421e9 +m_yy*( 1.06002e9 +m_yy*(-6.60078e8 +m_yy*( 6.33496e7 +m_yy*(-4.60396e7 +m_yy*( 1.4841e7  +m_yy*(-1.06352e6 +m_yy*(-217801.0 +m_yy*( 48153.3 +m_yy*(-1500.17 +m_yy*(-198.876 +m_yy*91.)))))))))));
	d8 = 5.79099e8 +m_yy*(-7.53828e8 +m_yy*( 6.60078e8 +m_yy*( 5.40367e7 +m_yy*( 1.99846e8 +m_yy*(-6.87656e6 +m_yy*(-6.89002e6 +m_yy*( 280428.0  +m_yy*( 161461.0 +m_yy*(-16493.7 +m_yy*(-567.163 +m_yy*364.0))))))))));
	e8 = 2.11107e8 +m_yy*(-2.89676e8 +m_yy*( 6.33496e7 +m_yy*(-1.99846e8 +m_yy*(-5.01017e7 +m_yy*(-5.25722e6 +m_yy*( 1.9547e6  +m_yy*( 240373.0  +m_yy*(-55582.0  +m_yy*(-1012.79 +m_yy*1001.0)))))))));
	f8 = 6.11148e7 +m_yy*(-7.01358e7 +m_yy*( 4.60396e7 +m_yy*(-6.87656e6 +m_yy*( 5.25722e6 +m_yy*( 3.04316e6 +m_yy*( 123052.0  +m_yy*(-106663.0  +m_yy*(-1093.82  +m_yy*2002.0))))))));
	g8 = 1.44647e7 +m_yy*(-1.39465e7 +m_yy*( 1.4841e7  +m_yy*( 6.89002e6 +m_yy*( 1.9547e6  +m_yy*(-123052.0  +m_yy*(-131337.0  +m_yy*(-486.14    +m_yy*3003.0)))))));
	h8 = 2.85721e6 +m_yy*(-2.84954e6 +m_yy*( 1.06352e6 +m_yy*( 280428.0  +m_yy*(-240373.0  +m_yy*(-106663.0  +m_yy*( 486.14    +m_yy*3432.0))))));
	o8 = 483737.0  +m_yy*(-498334.0  +m_yy*(-217801.0  +m_yy*(-161461.0  +m_yy*(-55582.0   +m_yy*(1093.82    +m_yy*3003.0)))));
	p8 = 70946.1   +m_yy*(-55600.0   +m_yy*(-48153.3   +m_yy*(-16493.7   +m_yy*(1012.79    +m_yy*2002.0))));
	q8 = 9504.65   +m_yy*(-3058.26   +m_yy*(-1500.17   +m_yy*( 567.163   +m_yy*1001.0)));
	r8 = 955.194   +m_yy*( 533.254   +m_yy*( 198.875   +m_yy*364.0));
	s8 = 126.532   +m_yy*( 40.5117   +m_yy*91.0);
	t8 = 3.68288   +m_yy*14.0;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz:Humlicek4		2014-2-13*/
/** **/
/*---------------------------------------------------------------------------*/

double skSpectralLineShapeStorageBuffer_VoigtKuntz::K4( double x) 
{
	double	xx;
	double	voigt;

	if (a8 == 0.0) ConfigureRegion4();

	xx = x*x;
	voigt = (a7 + xx*(b7+xx*(c7+xx*(d7+xx*(e7+xx*(f7+xx*(g7+xx*(h7+xx*(o7+xx*(p7+xx*(q7+xx*(r7+ xx*(s7+xx*t7)))))))))))))
			/
			(a8 + xx*(b8+xx*(c8+xx*(d8+xx*(e8+xx*(f8+xx*(g8+xx*(h8+xx*(o8+xx*(p8+xx*(q8+xx*(r8+xx*(s8+xx*(t8+xx))))))))))))));

	voigt = exp(m_yy-xx)*cos(m_y2*x) - voigt;
	return voigt;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::ResetCoefficients		2014-2-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShapeStorageBuffer_VoigtKuntz::ResetCoefficients()
{
	a2 = 0.0;			// Use the value of 0.0 to indicate that Y has changed and has not been refreshed
	a4 = 0.0;
	a5 = 0.0;
	a8 = 0.0;
	return true;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::ConfigureRegionBoundaries		2014-2-18*/
/** Find the x values for each of the 4 regions defined by Kuntz for a given
 *	value of Y **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShapeStorageBuffer_VoigtKuntz::ConfigureRegionBoundaries()
{
	m_maxnu = 1.0E030;
	m_minnu = 0.0;

	m_xlim1 = 15.0 - m_y;								// This is where region 1 finishes for abs(x). Its a line from (0,15), to (15,0) 
	m_xlim2 =  5.5 - m_y;								// This is where region 2 finishes
	if (m_xlim1 < 0.0) m_xlim1 = 0.0;
	if (m_xlim2 < 0.0) m_xlim2 = 0.0;
	m_xlim3 =  0.0;										// Region 3 always finishes at the Y axis
	if ( m_y < 0.75) m_xlim4 = (m_y + 0.176)/0.195;		// See if region 4 is valid
	else             m_xlim4 = 9.0E20;					// Set it to an unreachable value if it is invalid.
	return true;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::SetY		2014-2-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShapeStorageBuffer_VoigtKuntz::SetY( double userY)
{
	bool	ok;
	double	Y = fabs(userY);

	ok = (m_y == Y);
	if (!ok)
	{
		m_y = Y;
		ok  =		 ConfigureRegionBoundaries();
		ok = ok &&   ResetCoefficients();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::K		2014-2-18*/
/** Retyrn the value of K(x,y) using the technique outlined by Kuntz
 **/
/*---------------------------------------------------------------------------*/

double skSpectralLineShapeStorageBuffer_VoigtKuntz::K( double userx) 
{
	double	x = fabs(userx);
	double	v;

	if      ( x >= m_xlim1) v = K1(x);
	else if ( x >= m_xlim2) v = K2(x);
	else if ( x >= m_xlim4) v = K4(x);
	else                    v = K3(x);
	return v;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::Voigt		2014-2-18*/
/** **/
/*---------------------------------------------------------------------------*/

double skSpectralLineShapeStorageBuffer_VoigtKuntz::Voigt( double nu, double Snm ) 
{

	double x;
	double voigt;

	if ( (nu > m_minnu)  & ( nu < m_maxnu))
	{
		x     = fabs(nu - m_nu00)*m_xfactor;					//*sqrtln2/m_aD;
		voigt = Snm*K( x )*m_normfactor;
	}
	else
	{
		voigt = 0.0;
	}
	return voigt;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::AddVoigt		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShapeStorageBuffer_VoigtKuntz::AddVoigt( const std::vector<double>& nu, std::vector<double>* uservalue, double Snm ) 
{
	std::vector<double>::const_iterator	loweriter;
	std::vector<double>::const_iterator	upperiter;
	size_t								il,iu;
	double								x;
	double								f;

	NXASSERT( uservalue->size() == nu.size() );
	f = m_normfactor*Snm;
	loweriter = std::lower_bound( nu.begin(), nu.end(), m_minnu );
	upperiter = std::lower_bound( loweriter,  nu.end(), m_maxnu );
	il = loweriter - nu.begin();
	iu = upperiter - nu.begin();
	for ( size_t i  = il; i < iu;  i++)
	{
		x                 = fabs(nu.at(i) - m_nu00)*m_xfactor;
		uservalue->at(i) += K( x )*f;
	}
	return true;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::SetTolerance		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShapeStorageBuffer_VoigtKuntz::SetTolerance( double tolerance )
{
	m_tolerance = tolerance;
	return true;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz::SetLimitsFromParentMaxLineStrength		2014-2-27*/
/** Sets limits on the acceptable wavenumber range based upon the parent maximum line strength
 *	This function is normally called after SetLineParams has been called, usually in response
 *	to a "setAtmopshericState" somewehere else in the system. The idea is that the user in charge
 *	of  a micro-window finds the "brightest" line in that window and we only calculate signals to
 *	fraction of the brightest signal. For large X and Y in the Voigt function we can use the approximation,
 * \f[
 *	K(x,y)=frac{1}{\sqrt{\pi}}\frac{y}{x^2 +y^2} 
 *	\f]
 *	we then use the tolerance criteria that
 *	\f[
 *	S_{max}T = S_{line}K(x,y)
 *	where T is the overall tolerance of the micro-window.  This can be solved for x.
 *
 *	By default we will always go out to at least 15 ( used to be xlim2), which is approx 5 standard deviations
 *	from the line centre, even if there is no indication we need that level of accuracy.
 **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShapeStorageBuffer_VoigtKuntz::SetLimitsFromParentMaxLineStrength( double Smax, double Sline )
{
	double	x;
	double	f;
	const double rootinvpi = 0.56418958354775628694807945156077;		// 1/sqrt(pi)
	bool	ok;
	double	maxx;
	double	xhalf;

	NXASSERT(( Smax > 0) && (m_tolerance >= 0) && (m_xfactor > 0.0) );
	ok = (Smax > 0.0) ;
	if (ok && (m_tolerance > 0))
	{
		f = (Sline*rootinvpi/(Smax*m_tolerance)  - m_y)*m_y;
		x = (f > 0) ? sqrt(f) : 0.0;
	    xhalf = 5*0.5*( m_y + sqrt( m_y*m_y + 2.772588722239781));		// Get 5 times the empirical half width of the Voigt formula
		maxx = nxmax( xhalf, x);										// We want to make sure lines add to at least ~5 times the half width of the profile, so even weak lines add something if they are close.
//		NXTRACE(("%10.3f, %10.3f, %10.3f\n", (double)maxx, (double)x, (double)xhalf) );
	}
	else
	{
		maxx = 300000.0;												// This should ensure that every line goes through the system and makes a far wing contribution.
	}
	m_maxnu = m_nu00 + maxx/m_xfactor;
	m_minnu = m_nu00 - maxx/m_xfactor;
	return ok;
}

/*---------------------------------------------------------------------------
 *'					skSpectralLineShapeStorageBuffer_VoigtKuntz::SetLineParams                      2002-9-17
 *
 *	nu00              = lines central wavelength in wavenumbers
 *	pressure          = air atmospheric pressure in dynes, cgs 
 *	partialpressure   = partial pressure of the species in dynes 
 *	temperature       = temperature in Kelvin
 *	mass              = molecular mass in grams
 *	LorentzWidthMult  = Lorentz Half Width at 1 atmosphere.
 *
 *	Calculates the Half width Half max of the Doppler (m_aD) and the 
skSpectralLineShapeStorageBuffer_VoigtHumlik *	Half width half max of the Lorentz width (m_aL)
 *
 *	Doppler HWHM from Page 92, Exploration of the Solar System by Infra-red
 *	remote sensing.  Note their formula gives full width Half Max and we
 *`	want Half width half Max.
 *
 *	Lorentz pressure broadening from "The Hitran Molecular Spectroscopic Database
 *	and HAWKS (Hitran Atmospheric Workstation): 1996 Edition, Appendix equation A12
 *
 *-------------------------------------------------------------------------*/

bool skSpectralLineShapeStorageBuffer_VoigtKuntz::SetLineParams( double nu00,
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
	const double  pref    = 1013250.0;					// Reference pressure in Dynes

	m_nu00       = nu00;
	m_aD         = nu00*sqrt(twoln2*KBOLTZMAN*temperature/(mass))/CLIGHT;
	m_aL         = pow(tref/temperature,tempcoeff)*( (airhalfwidth*(pressure-partialpressure)/pref)+(selfhalfwidth*partialpressure/pref) );
	m_xfactor    = sqrtln2/m_aD;
	newy         = sqrtln2*m_aL/m_aD;
	m_normfactor = sqrtln2opi/m_aD;
	ok           = SetY(newy);
	return ok;
}

