
/*---------------------------------------------------------------------------
 *'					nxArrayTemporary::Copy Constructor				2002-9-9
 *	fast copy constructor for the nxArraTemporary.  Since the parent
 *	is guaranteed to be temporary then this will also be temporary
 *-------------------------------------------------------------------------*/

template <class T>
nxArrayTemporary<T>::nxArrayTemporary( const nxArrayTemporary<T>& other)
                 :nxArrayLinear<T>( other )
{
	NXTRACE(("nxArrayTemporary<T>:: COPY CONSTRUCTOR  invoked\n"));
	this->SetIsTemporary();
}

/*---------------------------------------------------------------------------
 *'					nxArrayTemporary<T>::PrepareBinaryOperator       2002-9-8
 *-------------------------------------------------------------------------*/

template <class T>
bool nxArrayTemporary<T>::PrepareBinaryOperator( const nxArrayLinear<T>& a, const nxArrayLinear<T>& b)
{
	RankSpecification			localrank;
	const RankSpecification*	rank;
	const RankSpecification*	arank = a.ArrayRankSpecs();
	const RankSpecification*	brank = b.ArrayRankSpecs();
//	InxMemoryManager<T>*		manager = NULL;
	nxBOOL						ok;

	ok = (a.size() == b.size()) ;										// operator only works if arrays have same number of elements
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "nxArrayTemporary<T>::PrepareBinaryOperator, the two arrays are not the same size (%d and %d)", (int)a.size(), (int)b.size());
	}
	else																	// Arrays are the same size
	{																		// so check for zero sized arrays
		ok = (a.size() > 0);												// Zero-sized arrays need no preparartion as there is nothing to do
		if (ok)																// its a finite size
		{																	// so

			// ---- Get a pointer to the highest dimension rank.
			// ---- Make sure the stride is contiguous for temporaries.

			rank = (arank->Rank() >= brank->Rank()) ? arank : brank;		// Get the highest dimension rank;
			if (!rank->IsContiguous())													// If the stride of this RankSpecification is not 1 (its not from a nxArrayTemporary
			{																			// so then make a copy of the rank specification
				localrank.Configure( rank->Rank(), rank->Dims(), sizeof(T), NULL );		// Make the copy
				rank = &localrank;														// and point to the copy instead.
			}																			// now rank points to a stride 1 configuration

			// ---- look for a temporary array, use it if available

			if ( a.IsTemporary() )											// If a is a temporary object
			{																// then
				ok = this->SetTemporaryArrayParams( a.ArrayManager(), rank, a.UnsafeArrayBasePtr() );
			}
			else if ( b.IsTemporary() )										// otherwise if b is temporary
			{																// then
				ok = this->SetTemporaryArrayParams( b.ArrayManager(), rank, b.UnsafeArrayBasePtr() );
			}
			else															// otherwise no option but to
			{																// allocate memory
				ok = this->SetSize( rank->Rank(), rank->Dims() );
			}

			// ---- Now configure this temporary memory allocation

			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "nxArrayTemporary<T>::PrepareBinaryOperator, Error allocating/locating a proper memory manager");
				this->Detach();
			}
		}
	}
	return ok;
}


