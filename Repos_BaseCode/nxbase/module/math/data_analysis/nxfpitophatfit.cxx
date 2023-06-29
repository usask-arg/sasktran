/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_math.h"
using namespace::nxmath;
#include "lsqanalytic.h"

//---------------------------------------------------------------------------
//					nxFpiTopHatFit::FP
//	Calculates the standard Fabry-Perot formula given a phase and
//	reflective finesse.
//
//	BRACK = Phase of airy function in radians.
//	F     = Finesse
//
//---------------------------------------------------------------------------

double nxFpiTopHatFit::FP( double BRACK ) 
{ 
   return 1.0 / (1.0 + m_F * sqr(sin(BRACK))); 
} 

//---------------------------------------------------------------------------
//						nxFpiTopHatFit::TOPHAT
//	Calculates the convolution of a tophat spectral line
//	with the Airy function plus a constant background.
//	It also saves the phases and terms used to generate
//	the convolution. This is needed when the derivatives
//	are to be calculated in CHI_DERIVS.
//
//---------------------------------------------------------------------------

const nxdblarr& nxFpiTopHatFit::TOPHAT( const nxdblarr& parameter, nxdblarr* yfit, int lim1, int lim2 )
{ 

	if (parameter.N_Elements() < NUMBER_OF_GAUSS_FRINGE_PARAMETERS )
	{
		nxLog::Record( NXLOG_WARNING, "nxFpiTopHatFit::Tophat, # of parameters os too small");
		return m_yfit;
	}
	double		Y0    = parameter [(int)PEAK_HEIGHT];
	double		R0    = parameter [(int)PEAK_POSITION];
	double		DTOP  = parameter [(int)PEAK_WIDTH];
	double		BACK  = parameter [(int)DC_BACKGROUND];

	int         I;
	double      AATAN1;												// temporary buffer for m_zplus[I]
	double      AATAN2;												// temporary buffer for m_zneg [I]
	double      NORM;												// Fitted data normalising factor
	double      SQRT1F;												// SQRT (1 + F)
	double      PIDFSR;												// PI divided by FSR

	if (DTOP < 1.0E-05) DTOP = 1.0E-05;								// avoid divide by zeroes
	SQRT1F  = sqrt(1 + m_F);										// get Sqrt of 1+F
	PIDFSR  = Pi / m_FSR;											// get PI divided by FSR
	m_zero    = PIDFSR * DTOP;										// get Z(o) phase
	m_z0      = atan2(SQRT1F * sin(m_zero), cos(m_zero));			// get m_z0
	NORM    = Y0 / (2 * m_z0);										// get a normalising factor
	for (I = lim1; I <= lim2; I++)								// for the range of data
	{																// get the fitted data
		m_plus[I] = PIDFSR * (I - R0 + DTOP);						// get the Z(+) phase
		m_negs[I] = PIDFSR * (I - R0 - DTOP);						// get the Z(-) phase
		AATAN1 = atan2(SQRT1F * sin(m_plus[I]), cos(m_plus[I])); 
		AATAN2 = atan2(SQRT1F * sin(m_negs[I]), cos(m_negs[I])); 
		if (AATAN2 > AATAN1) AATAN2 = AATAN2 - TWOPI; 
		m_zplus[I]   = AATAN1;										// save the Z(+) function
		m_zneg[I]    = AATAN2;										// save the Z(-) function
		m_yfit[I]    = NORM * (AATAN1 - AATAN2) + BACK;				// save the Y fitted
		yfit->At(I)  = m_yfit[I];
   } 
	return m_yfit;
} 

void nxFpiTopHatFit::operator() ( const nxdblarr& olda, nxdblarr* yfit, int lim1, int lim2)
{
	TOPHAT( olda, yfit, lim1, lim2 );
}


//---------------------------------------------------------------------------
//						chi_derivs 
//   Gets the  derivatives of chi-squared with respect to
//   Y0, R0, and DTOP. It will use the values of m_plus,m_negs
//   m_zero,m_zplus,m_zneg,m_z0 and m_yfit as calculated in procedure TOPHAT.
//---------------------------------------------------------------------------

#if defined(__BORLANDC__)					/* is this a Borland C++ compiler */
#pragma argsused							// I dont use lim1 and lim2
#endif


void nxFpiTopHatFit::CHI_DERIVS( const nxdblarr& parameters, nx2dArray<double>* derivs, int, int )
{ 
	double      FPP;
	double      FPN;
	double      NORM;
	double		BACK;
	double		Y0;
	double		FP0;

	Y0   = parameters[(int)PEAK_HEIGHT];
	BACK = parameters[(int)DC_BACKGROUND];
	FP0  = FP(m_zero);

	NORM        = Y0 * Pi * sqrt(1 + m_F) / (2.0 * m_FSR * m_z0);			// get normalising number
	for (int I = m_lim1; I <= m_lim2; I++)
	{
		FPP							   = FP(m_plus[I]);									// get the FP transmission
		FPN							   = FP(m_negs[I]);									// for m_plus and m_negs phases
		derivs->At( PEAK_HEIGHT, I)   = (m_yfit[I] - BACK) / Y0;								// get Y0, R0
		derivs->At( PEAK_POSITION, I) = NORM * (FPN - FPP);									// and DTOP derivatives 
		derivs->At( PEAK_WIDTH, I )   = NORM * (FPP + FPN - FP0 * (m_zplus[I] - m_zneg[I]) / m_z0); 
		derivs->At( DC_BACKGROUND,I)  = 1.0; 
	}
}


//---------------------------------------------------------------------------
//
//						nxFpiTopHatFit::Topper
//						----------------------
//	Least square fits the convolution of a tophat with the Airy
//	function to radius squared data.
//
//	Y0  := Peak height.
//	R0  := Peak position.
//	DTOP:= 1/2 width of top hat in bins.
//	BACK:= background level.
//
//	From my simple investigation it seems better to start with a
//	relatively large value of DTOP, eg 3 bins rather than a small
//	value. At small values of DTOP the routine is very close to a
//	perfect FP profile. The routine finds it hard to get the true
//	fit values, ie it homes into a local trough in ch-squared
//	rather than the true minimum. I would recommend at least 10
//	iterative loops before deciding to quit analysis.
//
//---------------------------------------------------------------------------

nxBOOL nxFpiTopHatFit::Topper( const nxdblarr&  y,
							   nxdblarr*        parameters,
						       int                     lim1,
							   int                     lim2,
							   double                  f,
							   double                  fsr,
							   int                     maxloop,
							   nxdblarr*        fittedline)
{

	int		npts      = y.N_Elements();
	nxBOOL	ok = nxFALSE;

	m_yfit.SetSize(npts);										// Get the fitting function storage space
	m_yfit.SetTo(0.0);											// and clear it all to zero

	if (parameters->N_Elements() != NUMBER_OF_GAUSS_FRINGE_PARAMETERS )
	{
		nxLog::Record(NXLOG_WARNING, "nxFpiTopHatFit::topper, TOPPER variable 'nxdblarr parameters' is not the correct size");
		return ok;
	}
	m_plus.SetSize(npts); 										// Holds the plus phases
	m_negs.SetSize(npts); 										// Holds the negative phases
	m_zplus.SetSize(npts);										// The Z(+) function
	m_zneg.SetSize(npts); 										// The Z(-) function
	m_derivs.SetSize(NUMBER_OF_GAUSS_FRINGE_PARAMETERS, npts);						// Holds the derivatives.

	m_F    = f;													// copy over the finess
	m_FSR  = fsr;												// the free spectral range in R**2 bins
	m_lim1 = lim1;												// the lower limit for the fit
	m_lim2 = lim2;												// and the upper limit for the fit

	parameters->At(PEAK_POSITION) -= 0.5;						// and shift peaks from centre of pixel
	ok = LsqAnalyticFit ( *this, y, parameters, lim1, lim2, maxloop, fittedline );
	parameters->At(PEAK_POSITION) += 0.5;
	return ok;
}


const nxdblarr&	nxFpiTopHatFit::SyntheticProfile( const nxdblarr& param, int npts, double fsr, double f )
{
	nxdblarr parameters( param);
	m_yfit.SetSize(npts);							// Get the fitting function storage space
	m_yfit.SetTo(0.0);							// and clear it all to zero

	if (parameters.N_Elements() < NUMBER_OF_GAUSS_FRINGE_PARAMETERS )
	{
		nxLog::Record(NXLOG_WARNING, "nxFpiTopHatFit::SyntheticProfile, TOPPER variable 'nxdblarr parameters' is not the correct size");
		return m_yfit;
	}

	m_plus.SetSize(npts); 								// Holds the plus phases
	m_negs.SetSize(npts); 								// Holds the negative phases
	m_zplus.SetSize(npts);								// The Z(+) function
	m_zneg.SetSize(npts); 								// The Z(-) function
	m_derivs.SetSize(NUMBER_OF_GAUSS_FRINGE_PARAMETERS, npts);						// Holds the derivatives.

	m_F    = f;									// copy over the finess
	m_FSR  = fsr;								// the free spectral range in R**2 bins
	m_lim1 = 0;									// the lower limit for the fit
	m_lim2 = npts-1;							// and the upper limit for the fit

	parameters[(int)PEAK_POSITION] -= 0.5;
	return TOPHAT( parameters,&m_yfit, m_lim1, m_lim2 );
}

 
