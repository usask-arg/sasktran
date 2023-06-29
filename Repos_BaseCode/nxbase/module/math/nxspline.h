/*-----------------------------------------------------------------------------
 *					class nxLinearInterpolate	2005-7-28*/
/** \ingroup Math_fitting
 *	A class for providing linear interpolation. **/
/*---------------------------------------------------------------------------*/

class nxLinearInterpolate
{
	public:
		enum EnumOutOfBoundAction	{ ENUM_INTERPOLATE, ENUM_TRUNCATE, ENUM_MISSINGVALUE};
	public:
		static double FromTwoPoints						(		double					x,
																double					x0,
																double					x1,
																double*					v
														);

		static double FromSquare						(		double					x,
																double					y,
																double					x0,
																double					x1,
																double					y0,
																double					y1,
																double*					v
														);

		static double EvaluateYatX						(		double					x,
																const double*			xarray,
																const double*			yarray,
																size_t					numpoints,
																EnumOutOfBoundAction	action,
																double					missingval,
																double					maxgapinx = -1.0
														);

		static double EvaluateYatX						(		double					x,
														  const std::vector<double>&	xarray,
														  const std::vector<double>&	yarray,
																EnumOutOfBoundAction	action,
																double					missingval,
																double					maxgapinx = -1.0
														);

		static double LogInterpolateYatX				(		double					x,
																const double*			xarray,
																const double*			yarray,
																size_t					numpoints,
																EnumOutOfBoundAction	action,
																double					missingvalue,
																double					maxgapinx = -1.0
														);
		
		static double LogInterpolateYatX				(		double					x,
														  const std::vector<double>&	xarray,
														  const std::vector<double>&	yarray,
																EnumOutOfBoundAction	action,
																double					missingvalue,
																double					maxgapinx = -1.0
														);

		static double LogInterpolateYatX				(		double					x,
																nxArrayIter<double>&	xbegin, 
																nxArrayIter<double>&	ybegin,
																size_t					numpoints,
																EnumOutOfBoundAction	action,
																double					missingvalue,
																double maxgapinx
														);

		static double EvaluateYatX						(		double					x,
																nxArrayIter<double>&	xarray,
																nxArrayIter<double>&	yarray,
																size_t					numpoints,
																EnumOutOfBoundAction	action,
																double					missingvalue,
																double					maxgapinx= -1.0
														);

		static bool	  FindBoundingIndicesDescending		( const	std::vector<double>&	location, 
																double					x,
																size_t*					lowercell,
																size_t*					uppercell,
																double*					lowerx,
																double*					upperx
														);

		template <class T, class CONSTITERATOR>
		static bool	 FindBoundingIndicesAscending		( CONSTITERATOR					start,
														  CONSTITERATOR					finish,
														  T								x,
														  size_t*						lowercell,
														  size_t*						uppercell,
														  T*							lowerx,
														  T*							upperx  );

		static bool	  FindBoundingIndicesAscending		( const	std::vector<double>&	location, 
																double					x,
																size_t*					lowercell,
																size_t*					uppercell,
																double*					lowerx,
																double*					upperx
														);
		static bool	FindBoundingIndicesAscendingCyclic	( const	std::vector<double>&	longitude,
																double					x,
																double					cyclicspan, 
																size_t*					lowercell,
																size_t*					uppercell,
																double*					lowerx,
																double*					upperx
														);

};
#include "nxlinearinterpolate_templates.hpp"


/*-----------------------------------------------------------------------------
 *					nxSpline		2005-3-17*/
/** \ingroup Math_fitting
 *	A class for fitting splines to data and interpolating and integrating.
 **/
/*---------------------------------------------------------------------------*/

class nxSpline
{
	private:
		int					m_numcoefs;
		double*				m_x;						// The X values of the tabulated data points
		double*				m_y;						// The Y values of the tabulated data points
		double*				m_y2;						// The second derivative of the interpolating function at the tabulated points
		double*				m_u;

	private:
		nxBOOL				Allocate( int ncoefs );
		void				ReleaseResources();

	public:
							nxSpline();
						   ~nxSpline();
		void				Clear		() { ReleaseResources();}
		bool				IsDefined	()  const { return m_numcoefs > 0;}
		bool				Configure	( const double* x,  const double* y, size_t n, double yp0, double ypn);
		bool				Configure	( const nx1dArray<double>&   x,  const nx1dArray<double>&    y, double yp0 = 1.1E30, double ypn = 1.1E30);
		bool				Configure	( const std::vector<double>& x,  const std::vector<double>&  y, double yp0 = 1.1E30, double ypn = 1.1E30);
		double				Interpolate	( double x, double badvalue ) const;
		double				Integrate	( nx1dArray<double>* integral ) const;
		bool				DeepCopy	( const nxSpline& other );



};


/*-----------------------------------------------------------------------------
 *					nxSpline		2005-3-17*/
/** \ingroup Math_fitting
 *	A class for fitting splines to data and interpolating and integrating.
 **/
/*---------------------------------------------------------------------------*/

class nxSpline2
{
	private:
		double				m_padvalue;
		size_t				m_startsplineindex;			// remove a sequence of constant values from the front of the spline array. Start spline at  this index. Typically this will return "begin() of the array if there is no leading padding
		size_t				m_endsplineindex;			// remove a sequence of constant values from the end   of the spline array. Finish spline at this index. Typically this will return "end()" of the array if there is no following padding
		nx1dArray<double>	m_x;						// The X values of the tabulated data points
		nx1dArray<double>	m_y;						// The Y values of the tabulated data points
		nx1dArray<double>	m_y2;						// The second derivative of the interpolating function at the tabulated points
		nx1dArray<double>	m_u;

	private:
		nxBOOL				Allocate				( size_t ncoefs );
		void				ReleaseResources		();
		bool				CheckBounds				(const double* x,  const double* y, size_t n);
		bool				CheckIsAscending		(const double* x,  size_t n );
		bool				CheckNotBadValues		(const double* y,  size_t n );
		bool				GetStartSplineIndex		(const double* y,  size_t n);
		bool				DeepCopy				( const nxSpline2& other );
		bool				GetEndSplineIndex		(const double* y,  size_t n);
		bool				ConfigureSplineSegment	(double* x,  double* y, double* y2, double* u, size_t n, double yp0, double ypn);
		bool				CopyXY					( const double*x, const double* y,  size_t n );

	public:
							nxSpline2			();
							nxSpline2			( const nxSpline2& other);
						   ~nxSpline2			();
		void				Clear				() { ReleaseResources();}
		bool				IsDefined			()  const { return m_x.size()> 0;}
		bool				SetPadValue			( double badvalue) { m_padvalue = badvalue; return true;}
		bool				Configure			( const double* x,  const double* y, size_t n,                  double yp0, double ypn);
		bool				Configure			( const nx1dArray<double>&   x,  const nx1dArray<double>&    y, double yp0 = 1.1E30, double ypn = 1.1E30);
		bool				Configure			( const std::vector<double>& x,  const std::vector<double>&  y, double yp0 = 1.1E30, double ypn = 1.1E30);
		double				Interpolate			( double x ) const;




};


/*-----------------------------------------------------------------------------
 *					nxPiecewiseLinear		2005-3-17*/
/** \ingroup Math_fitting
 *	A class for piecewise linear curves.
 **/
/*---------------------------------------------------------------------------*/

class nxPiecewiseLinear
{
	private:

		nxLinearInterpolate::EnumOutOfBoundAction	m_outofboundsaction;
		double										m_missingval;
		size_t										m_numcoefs;
		double*										m_x;						// The X values of the tabulated data points
		double*										m_y;						// The Y values of the tabulated data points

	private:
		nxBOOL				Allocate				( size_t ncoefs );
		void				ReleaseResources		();
		bool				DeepCopy				( const nxPiecewiseLinear& other );

	public:
							nxPiecewiseLinear		();
							nxPiecewiseLinear		(const nxPiecewiseLinear& other);
						   ~nxPiecewiseLinear		();
		bool				Configure				( const double* x,  const double* y, size_t n);
		bool				Configure				( const nx1dArray<double>&   x,  const nx1dArray<double>&    y);
		bool				Configure				( const std::vector<double>& x,  const std::vector<double>&  y);
		double				Interpolate				( double x, double badvalue ) const;
		void				SetExtrapolationMethod	( nxLinearInterpolate::EnumOutOfBoundAction action )		{ m_outofboundsaction = action; }
		void				SetMissingValue			( double value )											{ m_missingval        = value; }
		void				Clear					()															{ ReleaseResources();}
		bool				IsDefined				()  const													{ return m_numcoefs > 0;}
};


/*-----------------------------------------------------------------------------
 *					nxPolynomial		2007-7-23*/
/** \ingroup Math_fitting
 *	A class for handling polynomials functions.
 **/
/*---------------------------------------------------------------------------*/

class nxPolynomial
{
	private:
		nx1dArray<double>			m_coeffs;				// !< The array of coefficients, highest power of x first

	public:
									nxPolynomial			();
								   ~nxPolynomial			();
	   bool							SetCoeffs				( nxArrayLinear<double>& coeffs);
	   bool							SetLinear				( double m, double c );
	   bool							SetQuadratic			( double a, double b, double c );
	   bool							SetCoeffsFromPolyFit	( const nxArrayLinear<double>& user_x, const nxArrayLinear<double>& user_y, int norder );
	   const nx1dArray<double>&		Coeffs					() { return m_coeffs;}
	   bool							Evaluate				( nxArrayLinear<double>& x, nxArrayLinear<double>*  y );
	   bool							Evaluate				( double x, double* y );
	   void							Erase					() { m_coeffs.erase();} 
	   size_t						NumCoeffs				() { return m_coeffs.size();}

};
