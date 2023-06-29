#include "nxbase_math.h"

/*-----------------------------------------------------------------------------
 *					FindBoundingIndicesDescending		2011-7-25*/
/**  This is a simple interpolating helper function program that is useful for linear
 *	interpolation of data set in descending order. This code does not return the actual
 *	interpolation but instead returns the two bounding points for the interpolation. It is assumed that
 *	the user has another array of the same dimension as myarray that they will interpolate.
 *
 *	\param location
 *		The array of values that specify the location of the interpolation points placed in descending
 *		order.
 *
 *	\param x
 *		The value of the location where we need the interpolation. 
 *
 *	\param lowercell
 *	Returns the index in the location array of the lower bound point for point "x"
 *
 *	\param uppercell
 *	Returns the index in the location array of the upper bound point for point "x"
 *
 *	\param lowerx
 *	Returns the location of the lower bound point. This is the same as location[*lowercell]
 *
 *	\param upperx
 *	Returns the location of the upper bound point. This is the same as location[*uppercell]
 *
 *	\returns
 *	True if success.
**/
/*---------------------------------------------------------------------------*/

bool nxLinearInterpolate::FindBoundingIndicesDescending( const	std::vector<double>&	location,
																double					x,
																size_t*					lowercell,
																size_t*					uppercell,
																double*					lowerx,
																double*					upperx  )
{
	std::vector<double>::const_reverse_iterator		x1;							// The value just above our value of x
	std::vector<double>::const_reverse_iterator		x0;							// The value just below our value of x
	std::vector<double>::const_reverse_iterator		start;
	std::vector<double>::const_reverse_iterator		finish;
	bool											ok;
	bool											outofbounds;
	size_t											offset;

	NXASSERT( location.front() > location.back() );							// make sure the array is descending
	ok = (location.size() > 1 ) && (location.front() > location.back());		// Make sure we have a meaningful array
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"FindBoundingIndicesDescending, the array passed in is only 1 or zero elements. Thats too small to interpolate!");
	}
	else
	{
		start  = location.rbegin();										// get the start of the altitude grid
		finish = location.rend();											// Get the end of the altitude grid
		x1     = std::upper_bound( start, finish, x );				// Find the pointer to the value greater than x

		outofbounds = (x1 == start) || (x1 == finish);			// this is OK as long as it is not out of bounds
		if (outofbounds)										// if it is out of bounds
		{														// then
			if (x1 == start  ) x1 = start  + 1;					// adjust values before the start of the table
			if (x1 == finish ) x1 = finish - 1;					// adjust values beyond the end of the table
		}														// we now havethe upper value as sensible;
		x0           = x1 - 1;									// get the lower value
		offset       = location.size() - 1;
		*uppercell   = offset - (x1 - start);					// Get the index of the upper cell
		*lowercell   = offset - (x0 - start);					// Get the index of the lower cell
		*lowerx      = location[*lowercell];
		*upperx      = location[*uppercell];
	}
	if (!ok)
	{
		*uppercell = 1999999999;		// the first point is the last point in our cyclic series
		*lowercell = 1999999999;		// the second point is the first point in our array
		*lowerx    = 0;					// Get the value of the last point and adjust for cyclic behaviour.
		*upperx    = 0;
	}

	return ok;													// return ok.
}

/*-----------------------------------------------------------------------------
 *					FindBoundingIndicesDescending		2011-7-25*/
/**  This is a simple interpolating helper function program that is useful for linear
 *	interpolation of data set in descending order. This code does not return the actual
 *	interpolation but instead returns the two bounding points for the interpolation. It is assumed that
 *	the user has another array of the same dimension as myarray that they will interpolate.
 *
 *	\param location
 *		The array of values that specify the location of the interpolation points placed in descending
 *		order.
 *
 *	\param x
 *		The value of the location where we need the interpolation. 
 *
 *	\param lowercell
 *	Returns the index in the location array of the lower bound point for point "x"
 *
 *	\param uppercell
 *	Returns the index in the location array of the upper bound point for point "x"
 *
 *	\param lowerx
 *	Returns the location of the lower bound point. This is the same as location[*lowercell]
 *
 *	\param upperx
 *	Returns the location of the upper bound point. This is the same as location[*uppercell]
 *
 *	\returns
 *	True if success.
**/
/*---------------------------------------------------------------------------*/

bool nxLinearInterpolate::FindBoundingIndicesAscending( const	std::vector<double>&	location,
																double					x,
																size_t*					lowercell,
																size_t*					uppercell,
																double*					lowerx,
																double*					upperx  )
{
	bool											ok;

	ok = FindBoundingIndicesAscending<double, std::vector<double>::const_iterator>( location.begin(), location.end(), x, lowercell, uppercell, lowerx, upperx);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"FindBoundingIndicesDescending, the array passed in is only 1 or zero elements. Thats too small to interpolate!");
		*uppercell = 1999999999;		// the first point is the last point in our cyclic series
		*lowercell = 1999999999;		// the second point is the first point in our array
		*lowerx    = 0;					// Get the value of the last point and adjust for cyclic behaviour.
		*upperx    = 0;
	}
	return ok;													// return ok.
}


/*-----------------------------------------------------------------------------
 *					FindBoundingIndicesAscendingCyclic				2011-7-25*/
/** This is a simple interpolating helper function program that is useful for linear
 *	interpolation in a cyclic parameter such as longitude. The cyclic aspect
 *	is important when dealing with linear interpolation of the end points. Eg
 *	interpolation before the first point is really interpolation of the first and 
 *	last point wher the last points "angular" value has been reduced by the cyclic
 *	value.
 *
 *	This code does not return the actual interpolation but instead returns the two
 *	bounding points for the cyclic interpolation, which for longitude may be less than 0 or more than 360.
 *	In addition it returns the indices of the cells in longitude which should be used. It is assumed that
 *	the user has another array of the same dimension as longitude that they will interpolate.
 *
 *	\param longitude
 *		The array of longitudes (or cyclic values) of the interpolation points placed in ascending
 *		(longitude) order. The array is cyclic in that a value at longitude of 0 is the same as a
 *		value at longitude of 360.I have only tested code that starts at 0 and finishes at 360. The
 *		code will not work for non-zero starting ranges such as -180 to +180.
 *
 *	\param x
 *		The longitude (or cyclic variable) where we need the interpolation. will be enterpolated to
 *
 *	\param cyclicspan
 *	The lebngth of the cycle. Typically 360.0 for longitude
 *
 *	\param lowercell
 *	Returns the index in the longitude array of the lower bound point for point "x"
 *
 *	\param uppercell
 *	Returns the index in the longitude array of the upper bound point for point "x"
 *
 *	\param lowerx
 *	Returns the "longitude" of the lower bound point. This is usually the same as longitude[*lowercell]
 *	but can different by "cyclicspan" on occassion.
 *
 *	\param upperx
 *	Returns the "longitude" of the upper bound point. This is usually the same as longitude[*uppercell]
 *	but can different by "cyclicspan" on occassion.
 *
 *	\returns
 *	True if success.
 **/
/*---------------------------------------------------------------------------*/

bool nxLinearInterpolate::FindBoundingIndicesAscendingCyclic( const	std::vector<double>&		longitude,
																	double						x,
																	double						cyclicspan, 
																	size_t*						lowercell,
																	size_t*						uppercell,
																	double*						lowerx,
																	double*						upperx )
{
	std::vector<double>::const_iterator		x1;							// The value just above our value of x
	std::vector<double>::const_iterator		x0;							// The value just below our value of x
	std::vector<double>::const_iterator		start;
	std::vector<double>::const_iterator		finish;
	size_t									lastidx;
	bool									ok = true;
	bool									outofbounds;

	NXASSERT( longitude.front() < longitude.back() );							// make sure the array is Ascending
	ok = (longitude.size() > 1 ) && (longitude.front() < longitude.back());		// Make sure we have a meaningful array
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxLinearInterpolate::FindBoundingIndicesAscendingCyclic, the array passed in is only 1 or zero elements. Thats too small to interpolate!");
	}
	else
	{
		start  = longitude.begin();									// get the start of the longitude grid
		finish = longitude.end();									// Get the end of the longitude grid
		x1     = std::upper_bound( start, finish, x );				// Find the pointer to the value greater than x
		outofbounds = (x1 == start) || (x1 == finish);				// this is OK as long as it is not out of bounds
		if (outofbounds)											// if it is out of bounds
		{															// then
			lastidx = longitude.size() - 1;							// Get the index of the last array
			if (x1 == finish )										// if our value is greater than the last value (eg > 359 but less than 360)
			{														// then do the cyclic wrap
				*uppercell = 0;										// the 2nd point is the first index in the array
				*lowercell = lastidx;								// the 
				*lowerx    = longitude[lastidx];
				*upperx    = longitude[0] + cyclicspan;
			}
			else if (x1 == start)									// if the point is less than the first value (eg > 0 but less than 1.0
			{														// then
				*uppercell = lastidx;								// the first point is the last point in our cyclic series
				*lowercell = 0;										// the second point is the first point in our array
				*lowerx    = longitude[lastidx] - cyclicspan;		// Get the value of the last point and adjust for cyclic behaviour.
				*upperx    = longitude[0];
			}
			else
			{
				nxLog::Record(NXLOG_WARNING,"nxLinearInterpolate::FindBoundingIndicesAscendingCyclic, bad behaviour. This code loop should not happen but its amazing how often the impossible is not.");
				ok = false;
			}
		}
		else
		{	x0           = x1 - 1;									// we are not on the edges so ;
			*uppercell   = (x1 - start);							// Get the index of the upper cell
			*lowercell   = (x0 - start);							// Get the index of the lower cell
			*lowerx      = *x0;
			*upperx      = *x1;
		}
	}
	if (!ok)
	{
		*uppercell = 1999999999;		// the first point is the last point in our cyclic series
		*lowercell = 1999999999;		// the second point is the first point in our array
		*lowerx    = 0;					// Get the value of the last point and adjust for cyclic behaviour.
		*upperx    = 0;
	}
	return ok;													// return ok.
}


