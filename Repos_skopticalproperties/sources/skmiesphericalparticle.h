#include <complex>


/*-----------------------------------------------------------------------------
 *					class sk_MieSphericalWiscombeWrapper		2003-10-1*/
/** \ingroup miescat
 *	This is MIE scattering for a single spherical particle. The code
 *	provides a wrapper for the Mie scattering code written by Warren
 *	Wiscombe.  Wiscombe's Fortran code has been compiled "as is" only adding
 *	one line to make the function MIEV0 exportable from a DLL.
 *
 *	The Wiscombe Fortran code is compiled under the Fortran_Wiscombemie directory
 *	of the nxlib subproject directory.
 *
 *	Reference:
 *	W.J. Wiscombe, Improved Mie scattering algorithms. Appl Opt., 19,
 *	1505-1509, 1980.
 *
 *	Also see Warren Wiscombes FTP site
 *  Dr. Warren J. Wiscombe (wiscombe@climate.gsfc.nasa.gov)
 *        NASA Goddard Space Flight Center
 *        Code 913
 *        Greenbelt, MD 20771
 *
 *  FTP availability of the MIEV0 code: 
 *  The entire package is available by anonymous ftp
 *  from Internet site climate.gsfc.nasa.gov in subdirectory
 *  pub/wiscombe.  (ftp to 'climate', login as 'anonymous', give
 *  your e-mail address as password, then 'cd' to pub/wiscombe.)
 *
 *
 *	Basic operation.
 *	----------------
 *	The class computes the Mie scattering for a range of scattering angles.
 *	Wiscombe's code accepts either an array of angles symmetric about 90
 *	or arbitrary angles.  I have provided two functions that accept 
 *	eitehr configuration.
 *
 *  Most applications should call:
 *
 *	Set_Wavelength							???
 *  Set_Radius								???
 *  Set_SymmetricScatterAngles				???
 *  Set_RefractiveIndex						Default (1.33,0)
 *
 *  and may choose to call
 *
 *	Set_IsInfiniteRefractiveIndex			Default "Not Infinite"
 *  Set_MaxLegendreMoment					Default 0
 *	Set_IPolzn								Default +4321
 *	Set_MimCut								Default 0.0
 *
 *
 *
 *	I have written wrapper funtions for all of the Fortran variables
 *	used in Wiscombe's code (apart from the ones that are not required)
 *
 *	Note the following:
 *	Extinction cross-section, Cext = Qext . pi.a^2		a = radius of particle
 *  Scattering cross-section, Csca = Qsca x pi.a^2
 *	Absorbtion cross-section = (Cext - Csca)
 *
**/
/*---------------------------------------------------------------------------*/

class sk_MieSphericalWiscombeWrapper
{
	private:
		double								m_radius;		// Radius of particle
		double								m_lambda;		// Wavelength in same units as m_radius
		double								m_k;			// 2Pi/lambda
		double								m_xx;			// The Mie Size parameter (2*i*radius/wavelength)
		std::complex<double>				m_refr;			// The complex refractive index (real,imaginary)
		bool								m_perfect;		// TRUE if the refractive index is infinite Wiscombe's PERFCT
		double								m_mimcut;		// Value below which imaginary refractive index is 0.
		bool								m_anyang;		// If TRUE then any angle can be placed in m_xmu. If FALSE then angles are mirror symmetric about 90
		nx1dArray< double>					m_xmu;			// The array of cos(scattering angles) (0-Pi radians, always, monotone odd and always includes pi/2)
		int									m_nmom;			// The highest Legendre moment to calculate
		int									m_ipolzn;		// Wiscombe IPOLZN variable
		nx2dArray<double>					m_pmom;			// Wiscombe PMOM array
		double								m_qext;			// Output QEXT from Wiscombes code
		double								m_qsca;			// Output QSCA from Wiscombes code
		double								m_qabs;			// Absorption efficiency
		double								m_Cext;			// Extinction cross-section (in units of radius**2)
		double								m_Cabs;			// Absorption cross-section (in units of radius**2)
		double								m_Csca;			// Scattering cross-section (in units of radius**2)
		double								m_gqsc;			// Output GSQC from Wiscombes code
		nx1dArray< std::complex<double> >	m_S1;			// The Complex array of S1 Mie Scattering amplitudes
		nx1dArray< std::complex<double> >	m_S2;			// The Complex array of S1 Mie Scattering amplitudes
		std::complex<double>				m_Sforw;
		std::complex<double>				m_Sback;
		std::complex<double>				m_Tforw[2];
		std::complex<double>				m_Tback[2];
		double								m_spike;
		bool								m_isdirty;


	private:
		sk_MieSphericalWiscombeWrapper&		operator=       ( const sk_MieSphericalWiscombeWrapper& other );

		bool					SetDefaultAngles				( );
		bool					AllocateArrays					( size_t numangles );
		void					Wiscombe_Miev0					();
		void					SetDirty						(bool val ){ m_isdirty = m_isdirty || val;}
		bool					CalculateScattering				();
		void					UpdateXX						() { m_k = nxmath::TWOPI/m_lambda; m_xx = m_k*m_radius;}

	public:
								sk_MieSphericalWiscombeWrapper			();

	public:	// ---- Frequently called input parameters

		bool					Set_Wavelength					( double lamda );
		bool					Set_Radius						( double radius );
		bool					Set_SymmetricScatterAngles		( double *startangle, double *endangle, double *resolution, int numsteps);
		bool					Set_RefractiveIndex				( double ri_real, double ri_imag );
		bool					Set_AnyScatterAngles			( nx1dArray<double>& degrees);

		// ---- Less frequently called input parameters

		bool					Set_IsInfiniteRefractiveIndex	( bool val );
		bool					Set_MaxLegendreMoment			( int nmom );
		bool					Set_IPolzn						( int ipolzn );
		bool					Set_MimCut						( double mimcut );

		// ---- Get the curent status of input parameters
	

		double					Get_Wavelength					( ) const { return m_lambda;}
		double					Get_Radius						( ) const { return m_radius;}
		double					Get_MieSizeParameter			( )	const { return m_xx;}						
		bool					Get_IsInfiniteRefractiveIndex	( )	const { return m_perfect;}
		std::complex<double>	Get_RefractiveIndex				( ) const { return m_refr;}
		int						Get_MaxLegendreMoment			( )	const { return m_nmom;}
		int						Get_IPolzn						( )	const { return m_ipolzn;}
		double					Get_k							( ) const { return m_k;}


	public: // ---- Get outputs of MIE scattering calculation
		bool								DeepCopy		( const sk_MieSphericalWiscombeWrapper& other );
		const nx1dArray< double>&			Get_CosAngles	() const { return m_xmu;}
		size_t								Get_NumAngles   () const { return m_xmu.size();}
		double								Qext			();			// Get the extinction efficiency
		double								Qsca		();			// Get the scattering efficiency
		double								Qabs		();			// Get the absoptionefficiency
		double								Cext		();			// Get the extinction cross-section (area in units of radius**2)
		double								Csca		();			// Get the scattering cross-section (area in units of radius**2)
		double								Cabs		();			// Get the absoption  cross-section (area in units of radius**2)
		double								Gqsc		();			// Get the asymmetry factor
		nx1dArray< std::complex<double> >*	S1			();			// The Complex array of S1 Mie Scattering amplitudes
		nx1dArray< std::complex<double> >*	S2			();			// The Complex array of S1 Mie Scattering amplitudes
		nx2dArray<double>*					PMom		();			// The array of Legendre Polynomials (0:nmom, 4 )
		std::complex<double>				SForward	();
		std::complex<double>				SBackward	();
		std::complex<double>				TForward	(int i);
		std::complex<double>				TBackward	(int i);
		double								Spike		();
		bool								Get_ScatteringMatrixCoeffs( nx2dArray<double>* coeffmatrix );
		bool								Get_LegendreCoefficients(nx2dArray<double>* legendrecoeff);
};
