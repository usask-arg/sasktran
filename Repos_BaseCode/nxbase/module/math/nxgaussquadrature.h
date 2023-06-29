
typedef double (*nxQuadrature_funcptr)(double);

/*---------------------------------------------------------------------------
 *					class nxTrapezoidalQuadratureBase				2003-9-26*/
 /** \ingroup Math_quadrature_internals
 *	This is the base class for the trapezoidal quadrature.  It has all of the
 *	functionality apart from the templated Function Object.
 **/
 /*-------------------------------------------------------------------------*/

class nxTrapezoidalQuadratureBase
{
	protected:
		int					m_N;				// The "order" of the quadrature.
		double				m_x1;				// The value of X at the start of the interval
		double				m_x2;				// The value of X at the end of the interval

	public:
									nxTrapezoidalQuadratureBase	();						//!< The default constructor
		virtual					   ~nxTrapezoidalQuadratureBase	(){};					//!< The default destructor
		void						SetRange	( double x1, double x2 );				//!< Set the integration range
		void						SetOrder	( int N );								//!< Set the number of points in the range
		int							GetOrder	() const	{ return m_N;}				//!< The "order" of the quadrature.
		double						StartPoint	() const	{ return m_x1;}				//!< The value of X at the start of the interval
		double						EndPoint	() const	{ return m_x2;}				//!< The value of X at the end of the interval
		double						Integrate	( nx1dArray<double>& Y );				//!< Integrate table implicitly evaluated at X() points
//		double						Integrate	( double (* Yfunc)(double x)  );		// Integrate when a function can evaluate Y at location X
};


/*---------------------------------------------------------------------------
 *					class nxTrapezoidalQuadrature					2003-10-8*/
/**	\ingroup Math_quadrature
 *	A class for trapezoidal quadrature. Most of the configuration methods are in the
 *	base class #nxTrapezoidalQuadratureBase. This templated class allows the quadrature to be used with
 *	a function object or regular function.  The function object is a class that has an operator
 *	of the form \c double \c operator \c ()(double x). The operator evaluates the mathematical
 *	function at the point x.
 *	
 *	\par Simple Usage
 *	Users who wish to use a regular C function like cos(x) or sqrt(x) will
 *	create an instance of nxTrapezoidalQuadrature using the default template
 *  \c nxTrapezoidalQuadrature<>
 *
 *	\par Example
 *
 *	\code
 *	nxTrapezoidalQuadrature<>			quadrature;
 *	double								answer;
 *
 *	quadrature.SetRange	( 0, nxmath::TWOPI double x2 );
 *	quadrature.SetOrder	( 360 );
 *	answer = quadrature.Integrate( cos );
 *
 *	\endcode
**/
/*-------------------------------------------------------------------------*/

template <class FuncObj = nxQuadrature_funcptr>
class nxTrapezoidalQuadrature : public nxTrapezoidalQuadratureBase
{
	public:
		double						Integrate	( FuncObj function  );					//!< Integrate when a function can evaluate Y at location X
};

/*---------------------------------------------------------------------------
 *					class nxGaussQuadratureBase					  2003-9-26 */
/**	\ingroup Math_quadrature_internals
 *	This is the base class for the Guassian Quadrature.  It has all of the
 *	functionality apart from the templated Function Object.
 *
**/
/*-------------------------------------------------------------------------*/

class nxGaussQuadratureBase
{
	protected:
		nx1dArray<double>	m_x;				// The location of X coordinates
		nx1dArray<double>	m_weights;			// The weights for each coordinate
		int					m_N;				// The "order" of the quadrature.
		double				m_x1;				// The value of X at the start of the interval
		double				m_x2;				// The value of X at the end of the interval
		nxBOOL				m_isdirty;

	protected:
		void						CheckDirty();

	public:
									nxGaussQuadratureBase	();
		nxBOOL						ConfigureMultiRange		( double* startx, double* endx, int* ptsperrange, int numranges);

		void						SetRange	( double x1, double x2 );						//!< Set the range of integration
		void						SetOrder	( int N );										//!< Set the order
		const nx1dArray<double>&	X			()			{ CheckDirty(); return m_x;}		//!< Get the array of X quadrature coordinates, used to evaluate the function
		const nx1dArray<double>&	W			()			{ CheckDirty(); return m_weights;}	//!< Get the array of quadrature weights used in the integration
		int							GetOrder	() const	{ return m_N;}						//!< The "order" of the quadrature.
		double						StartPoint	() const	{ return m_x1;}						//!< The value of X at the start of the interval
		double						EndPoint	() const	{ return m_x2;}						//!< The value of X at the end of the interval
		double						Integrate	( nx1dArray<double>& Y );						//!< Integrate table implicitly evaluated at X() points
		bool						DeepCopy	( const nxGaussQuadratureBase& other );

};

/*---------------------------------------------------------------------------
 *					class nxGaussQuadrature							2003-10-8*/
/**	\ingroup Math_quadrature
 *	A class for Guassian quadrature. Gaussian quadrature implicitly allows the algorithm
 *	to select the points at which the function is evaluated (it uses the zeroes of the Legendre polynomial).
 *	This class allows the integration range to be broken into multiple sub-units, see #ConfigureMultiRange. This feature was
 *	introduced so an integral in radiative transfer over zenith angle could be broken into a upward and
 *	downward components to ovaoid the discontinuity at horizontal directions.
 *	The configuration methods are in the base class #nxGaussQuadratureBase. This templated class
 *	allows the quadrature to be used with a function object or regular function.  The function object is a class that has an operator
 *	of the form \c double \c operator \c ()(double x). The operator evaluates the mathematical
 *	function at the point x.
 *	
 *	\par Simple Usage
 *	Users who wish to use a regular C function like cos(x) or sqrt(x) will
 *	create an instance of nxGaussQuadrature using the default template
 *  \c nxGaussQuadrature<>
 *
 *	\par Example
 *
 *	\code
 *	nxGaussQuadrature<>			quadrature;
 *	double						answer;
 *
 *	quadrature.SetRange	( 0, nxmath::TWOPI double x2 );
 *	quadrature.SetOrder	( 10 );
 *	answer = quadrature.Integrate( cos );
 *
 *	\endcode
**/
/*-------------------------------------------------------------------------*/

template <class FuncObj = nxQuadrature_funcptr>
class nxGaussQuadrature : public nxGaussQuadratureBase
{
	public:
		double						Integrate	( FuncObj function  );					// Integrate when a function can evaluate Y at location X
		bool						DeepCopy   ( const nxGaussQuadrature& other )		{ return nxGaussQuadratureBase::DeepCopy(other);}
};

#include "nxgaussquadrature.hpp"


