


/*---------------------------------------------------------------------------
 *'					nxArrayIterVariableStride::nxArrayIterVariableStride    2003-6-6
 *-------------------------------------------------------------------------*/

template <class T>
nxArrayIterVariableStride<T>::nxArrayIterVariableStride( const nxArrayIterVariableStride<T>& other )
                             
{
	this->m_pointer		  = other.m_pointer;
	m_basepointer     = other.m_basepointer;
	m_rankspec		  = other.m_rankspec;
	m_logicalposition = other.m_logicalposition;
}

/*---------------------------------------------------------------------------
 *'					nxArrayIterVariableStride::operator=             2003-6-6
 *-------------------------------------------------------------------------*/

template <class T>
nxArrayIterVariableStride<T>& nxArrayIterVariableStride<T>::operator=( const nxArrayIterVariableStride<T>& other )
{
	this->m_pointer		  = other.m_pointer;
	m_basepointer     = other.m_basepointer;
	m_rankspec		  = other.m_rankspec;
	m_logicalposition = other.m_logicalposition;
	return *this;
}

/*---------------------------------------------------------------------------
 *'					nxArrayIterVariableStride::Configure             2003-6-6
 *-------------------------------------------------------------------------*/

template <class T>
void nxArrayIterVariableStride<T>::Configure(nxBYTE* basepointer, const RankSpecification* rankspec)
{
	m_basepointer     = basepointer;
	this->m_pointer         = (T*)basepointer;
	m_rankspec        = rankspec;
	m_logicalposition = 0;
}


/*---------------------------------------------------------------------------
 *'					nxArrayIter<T>::operator=                        2003-6-6
 *-------------------------------------------------------------------------*/

template <class T> nxArrayIter< T>& nxArrayIter<T>::operator=( const nxArrayIter<T>& other)
{
	if  (other.m_iterimp == (nxArrayIterCore<T>*)&other.m_contiguous )
	{
		m_contiguous = (other.m_contiguous);
		m_iterimp = &m_contiguous;
	}
	else if (other.m_iterimp == (nxArrayIterCore<T>*)&other.m_fixedstride)
	{
		m_fixedstride = other.m_fixedstride;
		m_iterimp     = &m_fixedstride;
	}
	else
	{
		m_variablestride = other.m_variablestride;
		m_iterimp        = &m_variablestride;
	}
	return *this;
}


/*---------------------------------------------------------------------------
 *'					nxArrayIter<T>::nxArrayIter                        2003-6-6
 *-------------------------------------------------------------------------*/

template <class T>
nxArrayIter<T>::nxArrayIter( nxBYTE* basepointer, const RankSpecification* rank)
{
	if      (rank->IsContiguous  ()) m_iterimp = &m_contiguous;
	else if (rank->IsFixedStride ()) m_iterimp = &m_fixedstride;
	else                             m_iterimp = &m_variablestride;
	m_iterimp->Configure( basepointer, rank );
}
