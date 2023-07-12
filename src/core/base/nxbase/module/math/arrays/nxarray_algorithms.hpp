namespace nxarray
{

/*-----------------------------------------------------------------------------
 *					Transpose		2009-7-31*/
/** \ingroup nxLinearArray_Algorithms
 *	returns the transpose of a in array t.
 **/
/*---------------------------------------------------------------------------*/

template <class T>
bool Transpose( const nx2dArray<T>& a, nx2dArray<T>* t  )
{
	bool			ok;
	size_t			numrows = a.YSize();
	size_t			numcols = a.XSize();

	ok = t->SetSize( numrows, numcols );
	if (ok)
	{
		for ( size_t i=0; i < numrows; i++ )
		{
			for ( size_t j=0; j < numcols; j++ )
			{
				t->At(i,j) = a.At(j,i);
			}
		}
	}
	else
	{
		t->erase();
		nxLog::Record( NXLOG_WARNING, "nxarray::Transpose, Error allocating memory for temporary array." );
	}
	return( ok );
}

/*--------------------------------------------------------------------------
 *					nxArray::Min                                2002-10-30*/
/** \ingroup nxLinearArray_Algorithms
 *	Returns the minimum value of an array.
 **/
/*------------------------------------------------------------------------*/

template <class T>
T  Min( const nxArrayLinear<T>& a )
{
	T   minvalue = 0;
	T	v;
	nxArrayIter<T>	aiter      = a.begin();
	nxArrayIter<T>	aend       = a.end();

	if (aiter < aend) minvalue = *aiter;
	while (aiter < aend)
	{
		v	     = *aiter;
		if ( v < minvalue) minvalue = v;
		++aiter;
	}
	return minvalue;
}

/*--------------------------------------------------------------------------
 *					nxArray::LesserOf                            2002-10-30*/
/** \ingroup nxLinearArray_Algorithms
 *	nxarray::LesserOf is in the nxarray namespace. Returns a copy of the array a
 *	where each element is less than or equal to #maxv. It is analogous to
 *	the IDL "<" operator
 **/
/*------------------------------------------------------------------------*/

template <class T>
nxArrayLinear<T>  LesserOf( const nxArrayLinear<T>& a, const T& maxv )
{
	nxArrayTemporary<T>		temp;
	nxBOOL					ok;

	temp.DeepCopy(a, nxFALSE);
	ok = (temp.size() == a.size());
	if (!ok)  nxLog::Record(NXLOG_WARNING, "nxarray::LesserOf, temp array is not correct size, operation not performed");
	else
	{
		nxArrayIter<T>	titer	   = temp.begin();
		nxArrayIter<T>	aiter      = a.begin();
		nxArrayIter<T>	aend       = a.end();
		while ( aiter != aend ) { *titer = (*aiter < maxv) ? *aiter : maxv; ++titer; ++aiter;}
	}
	return temp;
}

/*--------------------------------------------------------------------------
 *					nxArray::GreaterOf                           2002-10-30*/
/** \ingroup nxLinearArray_Algorithms
 *	Returns a copy of the array \e a  where each element is greater than or equal to
 *	#minv. It is analagous to the IDL ">" operator.
**/
/*------------------------------------------------------------------------*/

template <class T>
 nxArrayLinear<T>  GreaterOf( const nxArrayLinear<T>& a, const T& minv )
{
	nxArrayTemporary<T>		temp;
	nxBOOL					ok;

	temp.DeepCopy(a, nxFALSE);
	ok = (temp.size() == a.size());
	if (!ok)  nxLog::Record(NXLOG_WARNING, "nxarray::GreaterOf, temp array is not correct size, operation not performed");\
	else
	{
		nxArrayIter<T>	titer	   = temp.begin();
		nxArrayIter<T>	aiter      = a.begin();
		nxArrayIter<T>	aend       = a.end();
		while ( aiter != aend )
		{
			*titer = (*aiter > minv) ? *aiter : minv;
			++titer;
			++aiter;
		}
	}
	return temp;
}

/*--------------------------------------------------------------------------
 *					nxArray::Max                                   2002-10-30*/
/** \ingroup nxLinearArray_Algorithms
 *	Returns the maximum value of an array
 **/
/*------------------------------------------------------------------------*/
template <class T>
T  Max( const nxArrayLinear<T>& a )
{
	T   maxvalue = 0;
	T	v;
	nxArrayIter<T>	aiter      = a.begin();
	nxArrayIter<T>	aend       = a.end();

	if (aiter < aend) maxvalue = *aiter;
	while (aiter < aend)
	{
		v	     = *aiter;
		if ( v > maxvalue) maxvalue = v;
		++aiter;
	}
	return maxvalue;
}


/*-----------------------------------------------------------------------------
 *					nxarray::Total		2009-7-9*/
/** \ingroup nxLinearArray_Algorithms
 *	Returns the total of an array
 **/
/*---------------------------------------------------------------------------*/

template <class T>
T  Total( const nxArrayLinear<T>& a )
{
	T totalsig        = 0;

	nxArrayIter<T>	aiter      = a.begin();
	nxArrayIter<T>	aend       = a.end();

	while (aiter < aend)
	{
		totalsig += (*aiter++);
	}
	return totalsig;
}

/*-----------------------------------------------------------------------------
 *					nxarray::Total		2009-7-9*/
/** \ingroup nxLinearArray_Algorithms
 *	Returns the average of an array
 **/
/*---------------------------------------------------------------------------*/

template <class T>
T  Average( const nxArrayLinear<T>& a )
{
	T		avg;
	double npts;

	npts = (double)a.size();
	avg = (npts > 0) ? avg = (T)(Total( a )/npts) : 0;
	return avg;
}

/*-----------------------------------------------------------------------------
 *					Interpolate		2006-3-15*/
/** \ingroup nxLinearArray_Algorithms
 *	nxarray::Interpolate() is in the nxarray namespace. It linearly interpolates
 *	an array #y, specified at locations #xn to location #x and returns the value.
 *
 *	\param y
 *		The array of y ordinate values that will be interpolated.
 *
 *	\param xn
 *		The x-abscissa values corresponding to the array of y values. The user must ensure that
 *		xn and y have the same number of elements.
 *
 *	\param x
 *		The scalar abscissa value at which attay #y will be interpolated.
 **/
/*---------------------------------------------------------------------------*/

template <class T1, class T2>
T1 Interpolate( const nxArrayLinear<T1>& y, const nxArrayLinear<T2>& xn, T2 x)
{
	nxArrayIter<T2>	iter;
	nxArrayIter<T2>	first = xn.begin();
	nxArrayIter<T2>	last  = xn.end();
	double				f;
	intptr_t			idx0;
	intptr_t			idx1;
	double				x0;
	double				dx;
	double				x1;
	T1					v;

	iter = std::upper_bound( xn.begin(), xn.end(), x );
	if (iter == first)
	{
		f    = 1.0;
		idx0 = 0;
		idx1 = 0;
	}
	else if (iter == last)
	{
		f = 1.0;
		idx0 = (intptr_t)xn.size() - 1;
		idx1 = idx0;
	}
	else
	{
		idx1 = (iter - first);
		idx0 = idx1 - 1;
		x0   = xn.At(idx0 );
		x1   = xn.At(idx1 );
		dx   = (x1-x0);
		if (::fabs(dx) != 0.0) f = (x1-x)/dx;
		else                   f = 1.0;
	}
	v = y.At(idx0)*f + y.At(idx1)*(1.0 - f);
	return v;
}
/*-----------------------------------------------------------------------------
 *					nxarray::Interpolate		2006-3-15*/
/** \ingroup nxLinearArray_Algorithms
 *	This function is in the nxarray namespace.
 *	Linearly interpolates array #y, specified at locations #xn to locations #x
 *  and returns the results in #values.
 **/
/*---------------------------------------------------------------------------*/

template <class T1, class T2>
bool Interpolate( const nxArrayLinear<T1>& y, const nxArrayLinear<T2>& xn,  const nxArrayLinear<T2>& x, nxArrayLinear<T2> values)
{
	nxArrayIter<T2>	last  = x.end();
	nxArrayIter<T2>	viter;
	nxArrayIter<T2> xiter;
	double				f;
	intptr_t			idx0;
	intptr_t			idx1;
	T2					x0;
	T2					dx;
	T2					x1;
	T1					v;
	bool				ok;

	ok = values->SetSize( x.size() );
	if (ok)
	{
		viter = values->begin();
		for( xiter = x.begin(); !(xiter == last); ++xiter)
		{
			*viter = Interpolate( y, xn, *xiter );
			viter++;
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					WHERE		2004-12-13*/
/** \ingroup nxLinearArray_Algorithms
 *	This is a function that emulates the WHERE function IDL for two  It creates an
 *	array that indicates where the two arrays compare TRUE. The comparison
 *	is done with the templated NXCOMPARE_OBJECT object function. This function or
 *	object must overload the function operator such that it has
 *	bool comp( const T& a, const T& b)
**/
/*---------------------------------------------------------------------------*/
template <class T, class NXCOMPARE_OBJECT>
size_t where( nxArrayLinear<T>& array1, NXCOMPARE_OBJECT& comp, nxArrayLinear<T>& array2,  nxArrayLinear<size_t>* index )
{
	nxBOOL					ok;
	size_t					n = array1.size();
	nxArrayIter<T>			i1;
	nxArrayIter<T>			i2;
	size_t					logicalindex;
	size_t					ngood;

	logicalindex = 0;
	ngood        = 0;
	ok = (n == array2.size());
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "nxArrayLinear<T>::where, Cannot do comparison as the two arrays are not the same size");
	}
	if (!ok || n < 1)
	{
		index->erase();
	}
	else
	{
		index->SetSize( 1, &n);
		i1 = array1.begin();
		i2 = array2.begin();
		while (logicalindex < n)
		{
			if (comp( *i1, *i2))
			{
				index->At(&ngood) = logicalindex;
				ngood++;
			}
			++i1;
			++i2;
			++logicalindex;
		}
		index->TrimSize(ngood);
	}
	return ngood;
}

/*-----------------------------------------------------------------------------
 *					nxarray::where		2004-12-13*/
/** \ingroup nxLinearArray_Algorithms
 *	This is a function that emulates the WHERE function IDL for two  It creates an
 *	array that indicates where the two arrays compare TRUE. The comparison
 *	is done with the templated NXCOMPARE_OBJECT object function. This function or
 *	object must overfload the function operator such that it has
 *	bool comp( const T& a, const T& b)
**/
/*---------------------------------------------------------------------------*/
template <class T, class NXCOMPARE_OBJECT>
size_t where( nxArrayLinear<T>& array1, const NXCOMPARE_OBJECT& comp, T value,  nxArrayLinear<size_t>* index )
{
	size_t					    n = array1.size();
	nxArrayIter<T>			i1;
	size_t					logicalindex;
	size_t						ngood;

	logicalindex = 0;
	ngood        = 0;
	if (n < 1)
	{
		index->erase();
	}
	else
	{
		index->SetSize( 1, &n);
		i1 = array1.begin();
		while (logicalindex < (size_t)n)
		{
			if (comp(*i1, value))
			{
				index->At(&ngood) = logicalindex;
				ngood++;
			}
			++i1;
			++logicalindex;
		}
		index->TrimSize(ngood);
	}
	return ngood;
}

#include "median_quickselect.hpp"

// template<class A>  nxArrayLinear<A>	conjg ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( ::conjg ) }

/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	cos   ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( ::cos   ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	cosd  ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( nxmath::cosd  ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	acos  ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( ::acos  ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	acosd ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( nxmath::acosd ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	cosh  ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( ::cosh  ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	exp   ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( ::exp   ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	log   ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( ::log   ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	log10 ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( ::log10 ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	sin   ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( ::sin   ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	sind  ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( nxmath::sind  ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	asin  ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( ::asin  ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	asind ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( nxmath::asind ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	sinh  ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( ::sinh  ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	sqrt  ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( ::sqrt  ) }
/** \ingroup nxLinearArray_Algorithms */ template<class A>  nxArrayLinear<A>	fabs  ( const nxArrayLinear<A>& a ) { nxArrayLinearUnaryFunction( ::fabs  ) }

};
