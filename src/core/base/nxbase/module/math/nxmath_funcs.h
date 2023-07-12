/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

#define NX_NXMATH_H			1

/** \ingroup MATH_functions General Functions */
/** @{*/

/** The nxmath namespace */
#include <math.h>
namespace nxmath
{
	extern const double TWOPI;										//!< \f$2\pi\f$ as a double precision
	extern const double TwoPi;										//!< \f$\pi\f$ as a double precision
	extern const double Pi;											//!< \f$\pi\f$ as a double precision
	extern const double PiOver2;									//!< \f$\pi/2\f$ as a double precision
	extern const double ONE_DEGREE;									//!< One degree in radians.
	extern const double ONE_RADIAN;									//!< One radian in degrees
	double  TRUNC ( double x);										//!< returns complete integer less than X
	double  sqr   ( double x);										//!< return \f$x^{2}\f$;
	double  acosd ( double x);										//!< return \f$\cos^{-1}(x)\f$ in degrees
	double  cosd  ( double x);										//!< return \f$\cos(x)\f$ x in degrees
	double  asind ( double x);										//!< return \f$\sin^{-1}(x)\f$ in degree
	double  sind  ( double x);										//!< return \f$\sin(x)\f$  x in degrees
	double  atan2d( double y, double x);							//!< return \f$\tan^{-1}(\frac{y}{x})\f$ in degrees in range -180 to +180
	double  DegreesToRadians( double x);							//!< Convert degrees to radians
	double  RadiansToDegrees( double x);							//!< Convert radians to degrees
	double  inrange (double x, double r);							//!< return "x mod r"
	short   sign    (double x);										//!< -1 if x is negative, +1 if x is positive or zero
	double	round( double x );										//!< returns floor(x+0.5)
	double	PowerOf10Magnitude( double x, int *firstdigit );		//!< returns the power of 10 magnitude of X
};
/** @}*/


