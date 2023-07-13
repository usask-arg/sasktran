



/*-----------------------------------------------------------------------------
 *					nxarray::SwapElements		2009-7-13*/
/** \ingroup nxLinearArray_Algorithms
 *	Swaps two elements. This is really a worker function for
 *	function Median
 **/
/*---------------------------------------------------------------------------*/

template< class T >
void SwapElements( T& a, T&b )
{
	T	t(a);
	a = b;
	b = t;
}

/*-----------------------------------------------------------------------------
 *					Median		2009-7-13*/
/** \ingroup nxLinearArray_Algorithms
 *	Finds the median of an array using a quick select algorithm. The algorithm
 *	makes a copy of the array as it modifie sthe arry it searches. The allocation
 *	is probably the slow part. The rest is nice and fast.
 *  This Quickselect routine is based on the algorithm described in
 *  "Numerical recipes in C", Second Edition,
 *  Cambridge University Press, 1992, Section 8.5, ISBN 0-521-43108-5
 *  This code by Nicolas Devillard - 1998. Public domain
**/
/*---------------------------------------------------------------------------*/

template < class T >
T Median( const nxArrayLinear<T>& userarray ) 
{
	T				medianvalue;
    size_t			low;
	size_t			high;
    size_t			median;
    size_t			middle;
	size_t			ll;
	size_t			hh;
	nx1dArray<T>	arr;
	bool			ok;
	size_t			n;

	ok     = arr.DeepCopy( userarray );					// Make a copy of the users array
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "nxarray::Median, Error allocating a copy of the array. We cannot find the median, returning 0.");
		medianvalue = 0.0;
	}
	else
	{
		n      = arr.size();								// Get the size of the array
		low    = 0;											// reset the indices
		high   = n-1;										// for the start and end points
		median = (low + high) / 2;							// and the mid-point

		for (;;)
		{
			if (high <= low) break;							// Exit if we have repartitioned to one element only
			if (high == low + 1)							// Exit if we have repartitioned to two elements only
			{
				if (arr[low] > arr[high]) SwapElements(arr[low], arr[high]) ;
				break;
			}
			/* Find median of low, middle and high items; swap into position low */
	    
			middle = (low + high) / 2;
			if ( arr[middle] > arr[high]) SwapElements(arr[middle], arr[high]) ;
			if ( arr[low]    > arr[high]) SwapElements(arr[low],    arr[high]) ;
			if ( arr[middle] > arr[low] ) SwapElements(arr[middle], arr[low] ) ;

			/* Swap low item (now in position middle) into position (low+1) */
			SwapElements(arr[middle], arr[low+1]) ;

			/* Nibble from each end towards middle, swapping items when stuck */
			ll = low + 1;
			hh = high;
			for (;;)
			{
				do ll++; while (arr[low] > arr[ll]) ;
				do hh--; while (arr[hh]  > arr[low]) ;
				if (hh < ll) break;
				SwapElements(arr[ll], arr[hh]) ;
			}

			/* Swap middle item (in position low) back into correct position */
			SwapElements(arr[low], arr[hh]) ;

			/* Re-set active partition */
			if (hh <= median) low  = ll;
			if (hh >= median) high = hh - 1;
		}
		medianvalue = arr[median];
	}
	return medianvalue;
}

