
/*-----------------------------------------------------------------------------
 *					skSpectralLineShape_VoigtKuntz		2014-2-18*/
/** Implements the Voigt profile using the techniques outlined by Kuntz, 1997 and
*	corrected for the errata published by Ruyten(2004). I have intentionally kept
*	all of the coefficients identical to the published values even though more accurate
*	versions may be out there. The paper claims to have systematic fractional errors
*	of less than 10^-6 for the Voigt function over its entire domain .
*
*	\par Multi-Threading.
*	This line shape object supports multithreading simultaneous spectral lines where 
*	each spectral line uses a single trhead thread to calculate a voigt profile for
*	an array of wavelengths.
*
*	\par Line thresholding
*	The class and its associated buffer #skSpectralLineShapeStorageBuffer_VoigtKuntz
*	use a thresholding technique to reduce the number of wavenumbers considered for
*	any one spectral line. The basic premise is that we are tring to model the spectrum across
*	a micro-window at thousands of wavenumbers with hundreds or thousands of spectral lines. Higher level
*	code calculates the strongest line intensity from all of the spectral lines that contribute to the window
*	and then we only consider wavenumber regions that contribute more than fraction m_tolerance (see 
*	class skSpectralLineShapeStorageBuffer_VoigtKuntz) to the signal. This allows us to quickly limit the upper and lower
*	limits of X in the Voigt profile for a given Y. This technique often allows thousands of calculations to be eliminated
*	and provides massive speedups of code with no significant loss of precision.
*
*	\par References
*	-# M. Kuntz: A New Implementation of the Humlicek Algorithm For The Calculation Of The Voigt Profile Function, JQSRT, 57, 819, 1997.
*	-# W. Ruyten: Comment "A new implementation of the Humlicek...", JQSRT, 86, 231, 2003, doi:10.1016/j.jqsrt.2003.12.027
**/
/*---------------------------------------------------------------------------*/

class skSpectralLineShape_VoigtKuntz : public skSpectralLineShape
{
	public:
							skSpectralLineShape_VoigtKuntz		( );
		virtual			   ~skSpectralLineShape_VoigtKuntz		( ) override;
		virtual bool		LineShapeFunction					(                     double nu,            double*   uservalue, const skSpectralLine* spectralline, skSpectralLineShapeStorageBuffer* storagebuffer ) override;
		virtual bool		AddLineShapeFunctionArray			( const std::vector<double>& nu, std::vector<double>* uservalue, const skSpectralLine* spectralline, skSpectralLineShapeStorageBuffer* storagebuffer ) override;
		virtual bool		ConfigureLineParameters				( const skSpectralLine* spectralline,  double	temperature, double	pressure, const GEODETIC_INSTANT& geopt, skClimatology* atmopshericstate, skSpectralLineShapeStorageBuffer*	storagebuffer ) override;
		virtual bool		CreateStorageBuffer					( skSpectralLineShapeStorageBuffer** storagebuffer ) override;
		virtual bool		SetParentMaxLineStrength			( double parentmaxlinestrength, const skSpectralLine* spectralline, skSpectralLineShapeStorageBuffer* storagebuffer ) const override;
		virtual bool		SetTolerance						( double tolerance, const skSpectralLine* spectralline, skSpectralLineShapeStorageBuffer* storagebuffer ) override;
};


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_VoigtKuntz		2014-2-18*/
/** Implements the Voigt profile using the techniques outlined by Kuntz, 1997.
*	This is an exact implementation of the published Kuntz algorithm but has
*	been corrected for the errata published by Ruyten(2004). I have intentionally kept
*	all of the coefficients identical to the published values. The paper claims to
*	have systematic fractional errors of less than 10^-6 for the Voigt function
*	over its entire domain 
*
*	\par References
*	-# M. Kuntz: A New Implementation of the Humlicek Algorithm For The Calculation Of The Voigt Profile Function, JQSRT, 57, 819, 1997.
*	-# W. Ruyten: Comment "A new implementation of the Humlicek...", JQSRT, 86, 231, 2003, doi:10.1016/j.jqsrt.2003.12.027
*/
/*---------------------------------------------------------------------------*/

class skSpectralLineShapeStorageBuffer_VoigtKuntz : public skSpectralLineShapeStorageBuffer
{
	public:
		double	m_tolerance;		// The overall tolerance of the Voigt function, used by SetLimitsfromParentMaxLineStrength
		double	m_nu00;				// The wave number of the spectral line
		double	m_aD;				// The doppler halh width
		double	m_aL;				// The Lorenz half width.
		double	m_xfactor;
		double  m_normfactor;

		double	m_y;				// Current value of Y
		double  m_yy;				// Current value of y-squared (used by region 4)
		double	m_y2;				// Current value of 2.0*y     (used by region 4)

		double m_maxnu;				// Upper limit in wavenumber for this line and still exceed tolerance
		double m_minnu;				// Lower limit in wavenumber for this line and still exceed tolerance
		double	m_xlim1,m_xlim2,m_xlim3,m_xlim4;
		double	a1,b1;
		double	a2,b2;
		double	a3,b3,c3,d3;
		double	a4,b4,c4,d4;
		double	a5,b5,c5,d5,e5;
		double	a6,b6,c6,d6,e6;
		double  a7,b7,c7,d7,e7,f7,g7,h7,o7,p7,q7,r7,s7,t7;
		double  a8,b8,c8,d8,e8,f8,g8,h8,o8,p8,q8,r8,s8,t8;

	private:
		bool			ResetCoefficients			();
		bool			ConfigureRegionBoundaries	();
		double			K1							( double x ) ;
		double			K2							( double x ) ;
		double			K3							( double x ) ;
		double			K4							( double x ) ;
		void			ConfigureRegion1			( );
		void			ConfigureRegion2			( );
		void			ConfigureRegion3			( );
		void			ConfigureRegion4			( );
//		bool			SetY						( double Y);
//		double			K							( double userx);

	public:
						skSpectralLineShapeStorageBuffer_VoigtKuntz	();
		virtual		   ~skSpectralLineShapeStorageBuffer_VoigtKuntz	(){}
		bool			SetY										( double Y);
		double			K											( double userx) ;
		double			Voigt										( double nu, double Snm );
		bool			AddVoigt									( const std::vector<double>& nu, std::vector<double>* uservalue, double Snm );
		bool			SetLineParams								( double nu00, double pressure, double partialpressure, double temperature, double tref, double tempcoeff, double mass, double airhalfwidth, double selfhalfwidth);
		bool			SetLimitsFromParentMaxLineStrength			( double Smax, double Sline);
		bool			SetTolerance								( double tolerance);
};
