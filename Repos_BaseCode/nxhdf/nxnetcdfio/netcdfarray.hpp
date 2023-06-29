template <class nxNetcdfEntityType>
const nxNetcdfEntityType*	nxNetcdfEntityArray<nxNetcdfEntityType>::At( const char* name ) const
{
	std::string					dummy(name);
//	nxNetcdfEntityType			dummy;
	const_iterator				iter;
	const nxNetcdfEntityType*	ptr;

//	dummy.SetNCID(0);
//	dummy.SetName(name);

	iter = m_array.find( dummy );
	ptr  = (iter != m_array.end()) ? &(iter->second) : NULL;
	return ptr;
}



/*-----------------------------------------------------------------------------
 *					_findncid		2011-7-19*/
/** **/
/*---------------------------------------------------------------------------*/

template <class nxNetcdfEntityType>
class _findncid
{
	private:
		int		m_ncid;

	public:
					_findncid( int ncid )					{ m_ncid = ncid;}
		bool		operator () ( const typename nxNetcdfEntityArray<nxNetcdfEntityType>::value_type& a )	{ return a.second.NCID() == m_ncid;} 
};


/*-----------------------------------------------------------------------------
 *					nxNetcdfEntityArray<nxNetcdfEntityType>::Find		2011-7-19*/
/** **/
/*---------------------------------------------------------------------------*/

template <class nxNetcdfEntityType>
const nxNetcdfEntityType* nxNetcdfEntityArray<nxNetcdfEntityType>::Find( int ncid ) const
{
	const_iterator					iter;
	const nxNetcdfEntityType*	ptr;


	iter = std::find_if( m_array.begin(), m_array.end(), _findncid<nxNetcdfEntityType>(ncid ));
	ptr  = (iter != m_array.end()) ? &(*iter).second : NULL;
	return ptr;
}

/*-----------------------------------------------------------------------------
 *					nxNetcdfEntityArray<nxNetcdfEntityType>::AddBlankEntry		2011-7-15*/
/** Adds a blank entry to our array of Netcdfentitries. It returnms a pointer
 *	to the new entity if it worked, otherwise it returns NULL**/
/*---------------------------------------------------------------------------*/

/*template <class nxNetcdfEntityType>
const nxNetcdfEntityType* nxNetcdfEntityArray<nxNetcdfEntityType>::AddBlankEntry( int ncid, const char* name) const
{
	nxNetcdfEntityType				dummy;
	std::pair <iterator, bool>		p;
	const nxNetcdfEntityType*		ptr;


	dummy.SetNCID(ncid);
	dummy.SetName(name);

	p = m_array.insert( dummy );
	ptr = (p.second) ? &(*p.first) : NULL;
	return ptr;
}
*/
