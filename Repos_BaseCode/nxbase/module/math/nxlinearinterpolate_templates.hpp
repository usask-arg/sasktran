/*-----------------------------------------------------------------------------
 *					FindBoundingIndicesDescending		2011-7-25*/
/**  This is a simple interpolating helper function program that is useful for linear
 *	interpolation of data set in descending order. This code does not return the actual
 *	interpolation but instead returns the two bounding points for the interpolation. It is assumed that
 *	the user has another array of the same dimension as myarray that they will interpolate.
 *
 *	\param start
 *		The starting iterator. Usually the beginning of the array
 *		order.
 *
 *	\param finsih
 *		The end iterator. Usually the end of the array
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

template <class T, class CONSTITERATOR>
bool nxLinearInterpolate::FindBoundingIndicesAscending( CONSTITERATOR					start,
													    CONSTITERATOR					finish,
																T						x,
																size_t*					lowercell,
																size_t*					uppercell,
																T*						lowerx,
																T*						upperx  )
{
	CONSTITERATOR									x1;							// The value just above our value of x
	CONSTITERATOR									x0;							// The value just below our value of x
	bool											ok;
	bool											outofbounds;

	ok = ((finish-start) > 1) && ( (*start) <= (*(finish-1)) );	// make sure we have a meaningful array where the first point is less than the last point
	if (ok)
	{
		x1     = std::upper_bound( start, finish, x );			// Find the pointer to the value greater than x
		outofbounds = (x1 == start) || (x1 == finish);			// this is OK as long as it is not out of bounds
		if (outofbounds)										// if it is out of bounds
		{														// then
			if (x1 == start  ) x1 = start  + 1;					// adjust values before the start of the table
			if (x1 == finish ) x1 = finish - 1;					// adjust values beyond the end of the table
		}														// we now havethe upper value as sensible;
		x0           = x1 - 1;									// get the lower value
		*uppercell   = (x1 - start);							// Get the index of the upper cell
		*lowercell   = (x0 - start);							// Get the index of the lower cell
		*lowerx      = *x0;
		*upperx      = *x1;
	}
	return ok;													// return ok.
}

