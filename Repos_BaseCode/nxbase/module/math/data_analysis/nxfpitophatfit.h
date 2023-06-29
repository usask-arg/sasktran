
/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

/*---------------------------------------------------------------------------
 *						class nxdblarr										*/
/** \ingroup nxLinearArray_Internal
 *	\deprecated since 2001-01-01
 *	A class to represent double precision arrays. Used by the Fabry-Perot
 *	software and not much else.
**/
/*--------------------------------------------------------------------------*/

class  nxdblarr : public nx1dArray<double>
{
	private:
		nxBOOL			AdjustParabolaFittingWindow( int maxbin, double minvalue, int *lim1, int *lim2, int maxrad);
		nxBOOL 			QuadraticCentredOnPeak( int lim1, int lim2, double* peak, double* peakerr, double* intensity, double* A, double* B, double* C);

	public:
						nxdblarr(){}												// default constructor
						nxdblarr( int nx ):nx1dArray<double> ( (size_t)nx) {}				// constructor that makes an array
						nxdblarr( const nxdblarr& a) : nx1dArray<double> (a) {}
		nxdblarr&		operator= ( const nx1dArray<double>& a ) { nx1dArray<double>::operator=(a); return *this;}

		nxBOOL			LeastSquareQuadratic  ( double *A, double *B, double *C, double *peak, double *intensity, double *peakerr, int lim1, int lim2);
		nxBOOL			ParabolaPeakFitter    ( double *peak, double *peakerr, double *pkHeight, int* lim1, int* lim2, int maxrad );
		int				N_Elements				() const { return (int)size();}
		int				LimitBoundOfIndex		( int idx ) const { if (idx < 0) idx = 0;
																	if (idx >= N_Elements()) idx = N_Elements() - 1;
																	return idx;}
		double			MaxValueBelowThreshold(  double threshold, int lim1, int lim2, int* maxbin); 
		double			Min( int i1, int i2) const;
		double			Max( int i1, int i2) const;


};

/*-----------------------------------------------------------------------------
 *					nxFpiTopHatFit		2004-11-23*/
/** \ingroup Math_fitting
 *	Fits an approximate airy function to data. This was used for many years
 *	to analyse Fabry-Perot atmospheric data.  Its glorious days are past.  **/
/*---------------------------------------------------------------------------*/

class  nxFpiTopHatFit
{
	public:
	enum  { PEAK_HEIGHT	  = 0,
			PEAK_POSITION = 1,
			PEAK_WIDTH    = 2,
			DC_BACKGROUND = 3,
			NUMBER_OF_GAUSS_FRINGE_PARAMETERS = 4};

	private:
		int					m_lim1;
		int					m_lim2;
		double				m_F;
		double				m_FSR;
		nxdblarr			m_plus; 					// Holds the plus phases
		nxdblarr			m_negs; 					// Holds the negative phases
		nxdblarr			m_zplus;					// The Z(+) function
		nxdblarr			m_zneg; 					// The Z(-) function
		nxdblarr			m_yfit;						// Holds the fitted function (used to evalutae derivs)
		nx2dArray<double>	m_derivs;
		double				m_z0;
		double				m_zero; 					// Holds the zero phases

	private:
		double					FP              ( double BRACK );
		const nxdblarr&			TOPHAT          ( const nxdblarr& parameters, nxdblarr* yfit, int lim1, int lim2 );

	public:
		void					operator()      ( const nxdblarr& parameters, nxdblarr* yfit, int lim1, int lim2);								// get theoretical curve from
		void					CHI_DERIVS      ( const nxdblarr& parameters, nx2dArray<double>* derivs, int lim1, int lim2 );
		nxBOOL					Topper          ( const nxdblarr&  y, nxdblarr* parameters, int lim1, int lim2, double f, double fsr, int  maxloop = 10, nxdblarr* fittedline = NULL);
		const nxdblarr&			SyntheticProfile( const nxdblarr& parameters, int npts, double fsr, double f );


};
