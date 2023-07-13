#include "nxbase_math.h"
/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

namespace nxmath
{

const double TWOPI      	   = 6.283185307179586476925287;      	//@define Value of 2.PI to 24 decimal places
const double TwoPi      	   = 6.283185307179586476925287;      	//@define Value of 2.PI to 24 decimal places
const double Pi         	   = 3.141592653589793238462643;      	//@define Value of PI to 24 decimal places
const double PiOver2           = 1.5707963267948966192313215;      	//@define Value of PI to 24 decimal places
const double ONE_DEGREE 	   = 0.017453292519943295769237;      	//@define One degree expressed as radians to 24 decomal places
const double ONE_RADIAN        = 57.295779513082320876798155;     	//@define One radian expressed in degrees to 24 decimal places
//const double ARC_SECONDS_RAD   = 206264.8062						//@define Converts radians to arc-seconds to 9 decimal places

double TRUNC( double x)
{
	return floor(x);
}												//!< returns complete integer less than X

double sqr( double x)
{
   return x*x;
}

double acosd ( double x)
{
	if (x >  1.0) x = 1.0;
	if (x < -1.0) x = -1.0;
	return acos(x)*ONE_RADIAN;
}

double cosd( double x)
{
   return  cos(x*ONE_DEGREE);
}

double asind ( double x)
{
	if (x >  1.0) x =  1.0;
	if (x < -1.0) x = -1.0;
	return asin(x)*ONE_RADIAN;
}

double sind  ( double x)
{
   return  sin(x*ONE_DEGREE);
}

double atan2d( double y, double x)
{
	if ((y ==0) && ( x == 0)) return 0;
	return atan2(y,x)*ONE_RADIAN;
}

 double DegreesToRadians( double x)
{
   return x*ONE_DEGREE;
}

 double RadiansToDegrees( double x)
{
   return x*ONE_RADIAN;
}

double inrange (double x, double r)
{
	if (r == 0) x = 0;
	else x = fmod(x,r);
	if (x < 0.0) x+=r;
	return x;
}

short  sign(double x)
{
   if (x < 0.0) return -1; else return 1;
}

double round( double x )
{
	return floor( x + 0.5 );
}

/*<--------------------------------------------------------------------------
 *'					PowerOf10Magnitude
 *	Return the power of 10 that makes (absolute value of) the first digit of X fall in
 *	the range 1.0-9.999999..
 *
 *''PARAMETERS:
 *`		double x
 *			The value for which the power of 10 is required
 *
 *`		int *firstdigit
 *			returns the first power of 10 digit it will be in the range 1-9			
 *
 *''RETURNS:
 *	double, the power of 10, ie 10**N as opposed to N
 *
 *'	HISTORY:
 * 2002-10-9
 *>------------------------------------------------------------------------*/

double PowerOf10Magnitude( double x, int *firstdigit )
{
	double p10;
	int    digit;

	if (x == 0)
	{
		p10 = 1.0E-100;
		if (firstdigit != NULL) *firstdigit = 1;
	}
	else
	{
		double ax = fabs(x);
		p10 = pow(10.0, floor(log10(ax)));		// Get the power of 10 of the data range
		digit = (int)(ax/p10);					// get the first digit.
		if (digit == 10)						// Numerical round off might make 9.99999999999 equal to 10
		{										// if it does
			p10++;								// then increase the power of 10
			digit = 1;							// and set the digit equal to 0
		}
		if (digit == 0)							// Numerical round off may truncate 0.99999999999 to zero
		{										// if it does 
			digit = 1;							// then just set it to 1.
		}
		if (firstdigit != NULL)	*firstdigit = digit;
	}
	return p10;
}

};
