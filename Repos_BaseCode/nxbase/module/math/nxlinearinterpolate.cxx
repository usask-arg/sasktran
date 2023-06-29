#include "nxbase_math.h"



/*-----------------------------------------------------------------------------
 *				nxLinearInterpolate::fromTwoPoints							2005-7-28*/
/** Linearly interpolates within a square.
 *
 *	\param x
 *	The x coordinate where a value is required. Best results occur if x is between x0 and x1.
 *
 *	\param x0
 *		The lower bounding value of x. Typically the x-coordinate of the lower left corner of the square
 *
 *	\param x1
 *		The upper bounding value of x. Typically the x-coordinate of the upper right corner of the square
 *
 *	\param v
 *		An array of 2 values, One value for each end of the linear interpolation where:
 *		-# v[0] = v(x0)
 *		-# v[1] = v(x1)
 *
 *	\return
 *		The interpolated value.
 **/
/*---------------------------------------------------------------------------*/

double nxLinearInterpolate::FromTwoPoints( double x, double x0, double x1, double* v )
{
	double	f;
	double	dx;
	double	y;

	dx = x1-x0;

	if (dx == 0.0)
	{
		y = *v;
	}
	else
	{
		f = (x1-x)/dx;
		y = f*v[0] + (1-f)*v[1];
	}
	return y;
}

/*-----------------------------------------------------------------------------
 *				nxLinearInterpolate::EvaluateYatX							2005-7-28*/
/** Linearly interpolates within a square.
 *
 *	\param x
 *	The x coordinate where a value is required. Best results occur if x is between x0 and x1.
 *
 *	\param x0
 *		The lower bounding value of x. Typically the x-coordinate of the lower left corner of the square
 *
 *	\param x1
 *		The upper bounding value of x. Typically the x-coordinate of the upper right corner of the square
 *
 *	\param v
 *		An array of 2 values, One value for each end of the linear interpolation where:
 *		-# v[0] = v(x0)
 *		-# v[1] = v(x1)
 *
 *	\return
 *		The interpolated value.
 **/
/*---------------------------------------------------------------------------*/

double nxLinearInterpolate::EvaluateYatX( double x, nxArrayIter<double>& xarray, nxArrayIter<double>& yarray, size_t numpoints, EnumOutOfBoundAction action, double missingvalue, double maxgapinx  )
{
	double					y=0.0;
	bool					ok;
	nxArrayIter<double>		hiter;
	nxArrayIter<double>		start;
	nxArrayIter<double>		finish;
	double			yval[2];
	size_t			idx;
	double			gap;
	double			x0;
	double			x1;

	ok     = (numpoints > 1);
	if (!ok)
	{
		y = missingvalue;
	}
	else
	{
		finish   = xarray + numpoints;
		start    = xarray;
		hiter    = std::upper_bound( start, finish, x );
		ok       = (hiter > start) && (hiter < finish);
		if (!ok)
		{
			switch ( action)
			{
			case ENUM_TRUNCATE		:	if (hiter == start) y = *yarray;
										else                y = *(yarray + (numpoints-1) );
										break;

			case ENUM_MISSINGVALUE  :	y = missingvalue;
										break;

			case ENUM_INTERPOLATE	:	if (hiter == finish) hiter--;
										if (hiter == start) hiter++;
										ok = true;
										break;

			default					:	y = missingvalue;
										break;
			};
		}
		if (ok)
		{
			idx = (hiter - start) - 1;
			x0      = *(xarray + idx);
			x1      = *(xarray + idx + 1);
			gap     = x1-x0;
			if ( (maxgapinx < 0) || (gap < maxgapinx) )
			{
				yval[0] = *( yarray + idx  );
				yval[1] = *( yarray + idx + 1);
				y       = FromTwoPoints( x, *(xarray + idx), *(xarray + idx + 1), &yval[0] );
			}
			else
			{
				y = missingvalue;
			}
		}
	}
	return y;
}

/*-----------------------------------------------------------------------------
 *				nxLinearInterpolate::EvaluateYatX							2005-7-28*/
/** Linearly interpolates within a square.
 *
 *	\param x
 *	The x coordinate where a value is required. Best results occur if x is between x0 and x1.
 *
 *	\param xarray
 *		The array[numpoints] of x values. The array must be montonically increasing
 *
 *	\param yarray
 *		The array[numpoints] of y values that will be interpolated
 *
 *	\param numpoints
 *		The number of points in the xarray and yarray:
 *
 *	\param action
 *		The action to be taken when extrapolating outside the x array bounds.
 *		-# nxLinearInterpolate::ENUM_INTERPOLATE extrapolate from last two (nearest) end points
 *		-# nxLinearInterpolate::ENUM_TRUNCATE truncate to nearest end point
 *		-# nxLinearInterpolate::ENUM_MISSINGVALUE insert missing value
 *
 *	\param missingvalue
 *		The missing value to be inserted for out of range interpolation or maximum gap points this is often the quiet NaN
 *		value from std::numeric_limits<double>::quiet_NaN().
 *
 *	\param maxgapinx
 *		if the distance between 2 interpolating X values is greater than this value then dont interpolate but replace with missingvalue. If this value is 0 or less then
 *		ignore the gap test.
 *
 *	\

 *	\return
 *		The interpolated value.
 **/
/*---------------------------------------------------------------------------*/

double nxLinearInterpolate::EvaluateYatX( double x, const double* xarray, const double* yarray, size_t numpoints, EnumOutOfBoundAction action, double missingvalue, double maxgapinx  )
{
	double			y= 0.0;
	bool			ok;
	const double*	hiter;
	const double*	start;
	const double*	finish;
	double			yval[2];
	size_t			idx;
	double			gap;

	ok     = (numpoints > 1);
	if (!ok)
	{
		if  (     (numpoints == 1)
			  && ( (action == ENUM_TRUNCATE) || (action == ENUM_INTERPOLATE) )
			 ) y = yarray[0];
		else y = missingvalue;
	}
	else
	{
		finish   = xarray + numpoints;
		start    = xarray;
		hiter    = std::upper_bound( start, finish, x );
		ok       = (hiter > start) && (hiter < finish);
		if (!ok)
		{
			if ( action == ENUM_MISSINGVALUE)
			{
				if (hiter == finish)
				{
					double	xfrac;
					xfrac = ( x != 0.0) ? fabs( (x - xarray[numpoints-1])/x) : fabs(xarray[numpoints-1]);
					if (xfrac < 1.0E-05) action = ENUM_INTERPOLATE;
				}
				else if (hiter == start)
				{
					double	xfrac;
					xfrac = ( x != 0.0) ? fabs( (x - xarray[0])/x) : fabs(xarray[0]);
					if (xfrac < 1.0E-05) action = ENUM_INTERPOLATE;
				}
			}

			switch ( action)
			{
			case ENUM_TRUNCATE		:	if (hiter == start) y = yarray[0];
										else                y = yarray[numpoints-1];
										break;

			case ENUM_MISSINGVALUE  :	y = missingvalue;
										break;

			case ENUM_INTERPOLATE	:	if (hiter == finish) hiter--;
										if (hiter == start) hiter++;
										ok = true;
										break;

			default					:	y = missingvalue;
										break;
			};
		}
		if (ok)
		{
			idx = (hiter - start) - 1;
			gap = (xarray[idx+1] - xarray[idx]);
			ok = (maxgapinx < 0) || (gap <= maxgapinx);					
			if (ok)
			{
				yval[0] = yarray[idx];
				yval[1] = yarray[idx+1];
				y       = FromTwoPoints( x, xarray[idx], xarray[idx+1], &yval[0] );
			}
			else
			{
				y = missingvalue;
			}
		}
	}
	return y;
}


/*-----------------------------------------------------------------------------
 *					nxLinearInterpolate::EvaluateYatX		2012-11-6*/
/** **/
/*---------------------------------------------------------------------------*/

double nxLinearInterpolate::EvaluateYatX( double x, const std::vector<double>& xarray, const std::vector<double>& yarray, EnumOutOfBoundAction action, double missingvalue, double maxgapinx  )
{
	bool	ok;
	double	value;

	ok = ( xarray.size() == yarray.size());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "nxLinearInterpolate::EvaluateYatX, The two arrays are not the same size that is not good. Returning missing value");
		value = missingvalue;
	}
	else
	{
		if (xarray.size() > 0)
		{
			value = nxLinearInterpolate::EvaluateYatX( x, &xarray.front(), &yarray.front(), xarray.size(), action, missingvalue, maxgapinx  );
		}
		else
		{
			value = missingvalue;
		}
	}
	return value;

}



/*-----------------------------------------------------------------------------
 *				nxLinearInterpolate::EvaluateYatX							2005-7-28*/
/** Linearly interpolates within a square.
 *
 *	\param x
 *	The x coordinate where a value is required. Best results occur if x is between x0 and x1.
 *
 *	\param x0
 *		The lower bounding value of x. Typically the x-coordinate of the lower left corner of the square
 *
 *	\param x1
 *		The upper bounding value of x. Typically the x-coordinate of the upper right corner of the square
 *
 *	\param v
 *		An array of 2 values, One value for each end of the linear interpolation where:
 *		-# v[0] = v(x0)
 *		-# v[1] = v(x1)
 *
 *	\return
 *		The interpolated value.
 **/
/*---------------------------------------------------------------------------*/

double nxLinearInterpolate::LogInterpolateYatX( double x, const double* xarray, const double* yarray, size_t numpoints, EnumOutOfBoundAction action, double missingvalue, double maxgapinx  )
{
	double			y= 0.0;
	bool			ok;
	const double*	hiter;
	const double*	start;
	const double*	finish;
	double			yval[2];
	size_t			idx;
	double			gap;

	ok     = (numpoints > 1);
	if (!ok)
	{
		if  (     (numpoints == 1)
			  && ( (action == ENUM_TRUNCATE) || (action == ENUM_INTERPOLATE) )
			 ) y = yarray[0];
		else y = missingvalue;
	}
	else
	{
		finish   = xarray + numpoints;
		start    = xarray;
		hiter    = std::upper_bound( start, finish, x );
		ok       = (hiter > start) && (hiter < finish);
		if (!ok)
		{
			switch ( action)
			{
			case ENUM_TRUNCATE		:	if (hiter == start) y = yarray[0];
										else                y = yarray[numpoints-1];
										break;

			case ENUM_MISSINGVALUE  :	y = missingvalue;
										break;

			case ENUM_INTERPOLATE	:	if (hiter == finish) hiter--;
										if (hiter == start) hiter++;
										ok = true;
										break;

			default					:	y = missingvalue;
										break;
			};
		}
		if (ok)
		{
			idx = (hiter - start) - 1;
			gap = xarray[idx+1] - xarray[idx];
			if ( (maxgapinx < 0) || (gap <= maxgapinx))
			{
				yval[0] = yarray[idx];
				yval[1] = yarray[idx+1];
				if ((yval[0] <= 0.0) || (yval[1] <= 0.0))								// Dont do log interpolation on points less than or equal to 0.0
				{
					y  = FromTwoPoints( x, xarray[idx], xarray[idx+1], &yval[0] );
				}
				else
				{
					yval[0] = log( yval[0] );
					yval[1] = log( yval[1] );
					y       = FromTwoPoints( x, xarray[idx], xarray[idx+1], &yval[0] );
					y       = exp(y);
				}
			}
			else
			{
				y = missingvalue;
			}
		}		
	}
	return y;
}


/*-----------------------------------------------------------------------------
 *				nxLinearInterpolate::EvaluateYatX							2005-7-28*/
/** Linearly interpolates within a square.
 *
 *	\param x
 *	The x coordinate where a value is required. Best results occur if x is between x0 and x1.
 *
 *	\param x0
 *		The lower bounding value of x. Typically the x-coordinate of the lower left corner of the square
 *
 *	\param x1
 *		The upper bounding value of x. Typically the x-coordinate of the upper right corner of the square
 *
 *	\param v
 *		An array of 2 values, One value for each end of the linear interpolation where:
 *		-# v[0] = v(x0)
 *		-# v[1] = v(x1)
 *
 *	\return
 *		The interpolated value.
 **/
/*---------------------------------------------------------------------------*/

double nxLinearInterpolate::LogInterpolateYatX( double x, nxArrayIter<double>& xbegin, nxArrayIter<double>& ybegin, size_t numpoints, EnumOutOfBoundAction action, double missingvalue, double maxgapinx  )
{
	double				y= 0.0;
	bool				ok;
	nxArrayIter<double>	hiter;
	nxArrayIter<double>	start;
	nxArrayIter<double>	finish;
	double				yval[2];
	size_t				idx;
	double				gap;
	double				x0,x1;

	ok     = (numpoints > 1);
	if (!ok)
	{
		if  (     (numpoints == 1)
			  && ( (action == ENUM_TRUNCATE) || (action == ENUM_INTERPOLATE) )
			 ) y = *ybegin;
		else y = missingvalue;
	}
	else
	{
		finish   = xbegin + numpoints;
		start    = xbegin;
		hiter    = std::upper_bound( start, finish, x );
		ok       = (hiter > start) && (hiter < finish);
		if (!ok)
		{
			switch ( action)
			{
			case ENUM_TRUNCATE		:	if (hiter == start) y = *ybegin;
										else                y = *(ybegin +(numpoints - 1));
										break;

			case ENUM_MISSINGVALUE  :	y = missingvalue;
										break;

			case ENUM_INTERPOLATE	:	if (hiter == finish) hiter--;
										if (hiter == start) hiter++;
										ok = true;
										break;

			default					:	y = missingvalue;
										break;
			};
		}
		if (ok)
		{
			idx = (hiter - start) - 1;
			x1 =  *(xbegin + (idx+1));
			x0 = *(xbegin + idx);
			gap = x1 - x0;
			if ( (maxgapinx < 0) || (gap <= maxgapinx))
			{
				yval[0] = log( *(ybegin + idx)   );
				yval[1] = log( *(ybegin + idx+1) );
				y       = FromTwoPoints( x, x0, x1, &yval[0] );
				y       = exp(y);
			}
			else
			{
				y = missingvalue;
			}
		}		
	}
	return y;
}

/*-----------------------------------------------------------------------------
 *					nxLinearInterpolate::LogInterpolateYatX		2012-5-5*/
/** **/
/*---------------------------------------------------------------------------*/

double nxLinearInterpolate::LogInterpolateYatX( double x, const std::vector<double>& xarray, const std::vector<double>& yarray, EnumOutOfBoundAction action, double missingvalue, double maxgapinx  )
{
	bool	ok;
	double	value;

	ok = ( xarray.size() == yarray.size());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "nxLinearInterpolate::LogInterpolateYatX, The two arrays are not the same size that is not good. Returning missing value");
		value = missingvalue;
	}
	else
	{
		if (xarray.size() > 0)
		{
			value = nxLinearInterpolate::LogInterpolateYatX( x, &xarray.front(), &yarray.front(), xarray.size(), action, missingvalue , maxgapinx );
		}
		else
		{
			value = missingvalue;
		}
	}
	return value;

}

/*-----------------------------------------------------------------------------
 *				nxLinearInterpolate::FromSquare							2005-7-28*/
/** Linearly interpolates within a square.
 *
 *	\param x
 *	The x coordinate where a value is required. Best results occur if x is between x0 and x1.
 *
 *	\param y
 *	The y coordinate where a value is required. Best results occur if y is between y0 and y1.
 *
 *	\param x0
 *		The lower bounding value of x. Typically the x-coordinate of the lower left corner of the square
 *
 *	\param y0
 *		The lower bounding value of y. Typically the y-coordinate of the lower left corner of the square
 *
 *	\param x1
 *		The upper bounding value of x. Typically the x-coordinate of the upper right corner of the square
 *
 *	\param y1
 *		The upper bounding value of y. Typically the y-coordinate of the upper right corner of the square
 *	
 *	\param v
 *		An array of 4 corner values, One value for each of the four corners of the square where
 *		-# v[0] = v(x0, y0)
 *		-# v[1] = v(x0, y1)
 *		-# v[2] = v(x1, y1)
 *		-# v[3] = v(x1, y0)
 *	Typically the 4 corners are specified in a clockwise fashion starting from the bottom left corner.
 *
 *	\return
 *		The interpolated value.
 **/
/*---------------------------------------------------------------------------*/

double nxLinearInterpolate::FromSquare( double x, double y, double x0, double x1, double y0, double y1, double* v)
{
	double				dy;
	double				dx;
	double				dxdy;
	double				d00,d10,d01,d11;
	double				dy1,dx1,dy0,dx0;
	double				answer;

	dx  = (x1-x0);
	dy  = (y1-y0);
	dxdy = dx*dy;
	if (dxdy == 0.0)														// IF dx or dy is zero
	{																		// then its a special case
		if (dx == 0.0)														// if dx is zero
		{																	// then
			answer = nxLinearInterpolate::FromTwoPoints( y, y0,y1, v );		// linearly interpolate in Y using v[0] and v[1]
		}																	// otherwise
		else																// dy must be zero
		{																	// so
			answer = nxLinearInterpolate::FromTwoPoints( x, x0,x1, &v[1] );	// linearly interpolate in X using v[1] and v[2]
		}
	}
	else																	// otherise our corner points are sensibly seperated 
	{																		// so
		dy1    = (y1-y);
		dy0    = (y-y0);
		dx1    = (x1-x);
		dx0    = (x-x0);

		d00    = dy1*dx1;							// get all of the linear combinations 
		d10    = dy1*dx0;
		d01    = dy0*dx1;
		d11    = dy0*dx0;
													// and evaluate the answer.
		answer = (d00*v[0] + d10*v[3] + d01*v[1] + d11*v[2])/dxdy;
	}
	return answer;
}
