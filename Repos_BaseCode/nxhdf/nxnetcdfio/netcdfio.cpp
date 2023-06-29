//#include "makeecmwf.h"
#include "nxnetcdfio.h"

/*-----------------------------------------------------------------------------
 *					nxNetcdfEntity::DeepCopy		2011-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfEntity::DeepCopy		( const nxNetcdfEntity& other)
{
	m_ncid = other.m_ncid;
	m_name = other.m_name;
	m_parent = other.m_parent;
	return true;
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfEntity::ParentGroup		2011-7-19*/
/** **/
/*---------------------------------------------------------------------------*/

const nxNetcdfGroup* nxNetcdfEntity::ParentGroup() const
{
	return dynamic_cast<const nxNetcdfGroup*>(m_parent);
}

/*-----------------------------------------------------------------------------
 *					nxNetcdfEntity::RootGroup		2011-7-19*/
/** **/
/*---------------------------------------------------------------------------*/

const nxNetcdfGroup* nxNetcdfEntity::RootGroup() const
{
	const nxNetcdfEntity*	parent;


	parent = this;
	while ( parent->m_parent != NULL)
	{
		parent = parent->m_parent;
	}
	return dynamic_cast<const nxNetcdfGroup*>(parent);
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfVar::LoadAttributes		2011-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfVar::LoadAttributes( )
{
	bool	ok;

	ok = m_atts.LoadAttributes( this );
	return ok;
}



/*-----------------------------------------------------------------------------
 *					nxNetcdfVar::AdjustForScaleAndOffset		2011-7-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfVar::AdjustForScaleAndOffset( nxArrayLinear<double>* vararray, double usermissingvalue ) const
{
	double					scale = 1.0;
	double					offset = 0.0;
	double					missingvalue = HUGE_VAL;
	double					v;
	bool					gotfill, gotscale, gotoffset, ok;
	nxArrayIter<double>		iter;
	nxArrayIter<double>		enditer;

	gotfill   = AttributeDouble( "_FillValue",   &missingvalue);
	gotscale  = AttributeDouble( "scale_factor", &scale);
	gotoffset = AttributeDouble( "add_offset",   &offset);
	ok        = gotfill || gotscale || gotoffset;

	if (ok)
	{
		enditer = vararray->end();
		iter    = vararray->begin();
		while (!(iter == enditer))
		{
			v = *iter;
			if (v == missingvalue)
			{
				v = usermissingvalue;
			}
			else
			{
				v = v*scale + offset;
			}
			*iter = v;
			++iter;
		}
	}
	return true;
}

/*-----------------------------------------------------------------------------
 *					nxNetcdfVar::LoadData		2011-7-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfVar::LoadData( nxArrayLinear<double>* vararray, double usermissingvalue ) const
{
	bool					ok;
	nx1dArray<size_t>		rankarray;

	ok =       LoadRankSpecs( &rankarray );
	ok = ok && vararray->SetSize( (int)rankarray.size(), rankarray.UnsafeArrayBasePtr() );
	vararray->SetTo(-999999.0);
	ok = ok && ( nc_get_var_double(ParentNCID(), NCID(), vararray->UnsafeArrayBasePtr() ) == NC_NOERR);
	ok = ok && AdjustForScaleAndOffset( vararray, usermissingvalue );

	if (!ok )
	{
		vararray->erase();
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *                   nxNetcdfVar::LoadDataSlice                   2019-10-03 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfVar::LoadDataSlice( nxArrayLinear<double>*	 vararray, const nxNetcdfHyperSlabDefn& hyperslab, double usermissingvalue) const
{
	bool	ok;
	int		err;

	ok = hyperslab.AllocateArrayToHoldHyperSlab( vararray );
	if (ok)
	{
//		vararray->SetTo(std::numeric_limits<double>::quiet_NaN() );
		const size_t* start=hyperslab.Start();
		const size_t* count= hyperslab.Count();
		err = nc_get_vara_double	( ParentNCID(), NCID(),   start, count , vararray->UnsafeArrayBasePtr());
		ok = ok && (err == NC_NOERR);
		ok = ok && AdjustForScaleAndOffset( vararray, usermissingvalue );
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxNetcdfVar::LoadVariable		2011-7-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfVar::LoadVariable( nxArrayLinear<double>* vararray, std::vector<nxNetcdfDimStruct>* dimsinfo, double usermissingvalue  )const
{
	bool				ok;

	ok =        LoadData    ( vararray, usermissingvalue );
	ok = ok &&  LoadDimInfo ( dimsinfo );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfVar::AttributeString		2011-7-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfVar::AttributeString( const char* attributename, std::string* value ) const
{
	bool	    ok;
	size_t		len;
	char*		strbuffer;

	ok =       ( nc_inq_attlen(ParentNCID(), NCID(), attributename, &len ) == NC_NOERR);
	if (ok)
	{
		strbuffer = new char[len+1];
		strbuffer[len] = '\0';
		ok = (nc_get_att_text(ParentNCID(), NCID(), attributename, strbuffer ) == NC_NOERR);
		if (ok)
		{
			*value = strbuffer;
		}
		delete [] strbuffer;
	}
	if (!ok) *value = "";
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfVar::AttributeDouble		2011-7-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfVar::AttributeDouble			( const char* attributename, double* value ) const
{
	bool	ok;

	ok = (nc_get_att_double(ParentNCID(), NCID(), attributename, value) == NC_NOERR);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxNetcdfVar::LoadRankSpecs		2011-7-19*/
/** Loads the rank specifications of this variable**/
/*---------------------------------------------------------------------------*/

bool nxNetcdfVar::LoadRankSpecs( nx1dArray<size_t>* rankarray ) const
{
	bool				ok;
	int					ndims;
	std::vector<int>	dims;
	int					i;
	int					idx;
	size_t				len;
	const nxNetcdfGroup*		group;
	const nxNetcdfDim*			dim;		


	ok =  (nc_inq_varndims (ParentNCID(), NCID(), &ndims) == NC_NOERR);								// get the number of dimensions
	if (ok )																						// If we ARE GOOD
	{																								// then see if we have a real array or a scalar
		if ( ndims > 0)																				// if we have an array
		{																							// then
			dims.resize(ndims);																		// allocate temporary (int) space to hold the size of each dimension
			ok    =  rankarray->SetSize( ndims );													// allocate final (size_t) space to hold the size of each dimension 
			ok    =  ok && (nc_inq_vardimid ( ParentNCID(), NCID(), &dims.front() ) == NC_NOERR);	// get the dimension netcdf id codes of netcdf dimension entities
			group = ParentGroup();																	// get the parent group which we shall use to look up the dimension entities 
			ok    = ok && (group != NULL);															// make sure the parent group is sensible
			if (ok)																					// if it is
			{																						// then
				idx = ndims;
				for (i = 0; i < ndims; i++)															// for each dimension entity id
				{																					// we are going to find the netcdf entity 
					idx--;																			// but we shall go in reverse to match how data is loaded in.
					dim              = group->DimAt(dims[idx]);										// the netcdf dimension entity
					len              = (dim != NULL) ? dim->Length() : 0;							// and get the Length of the dimension, If we have no dimension set it to 0.
					rankarray->At(i) = len;															// save the length of this dimension
					ok = ok && (len > 0);															// make sure we have nothing that is zero or less
				}																					// and do all of the dimensions
			}																						// and that is that for an real array
		}																							// if we have a scalar
		else																						// then
		{																							// emulate
			rankarray->SetSize( 1 );																// with an array of 1 dimension
			rankarray->At( 0 ) = 1;																	// of length 1
		}
	}
	if (!ok )
	{
		rankarray->erase();
		nxLog::Record(NXLOG_WARNING,"nxNetcdfVar::LoadRankSpecs, There was an error loading teh rank specifications for variable %s", (const char*)Name().c_str() );
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfVar::getVarAttributes		2011-7-20*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfVar::GetVarAttributes( nxNetcdfDimStruct* info ) const
{
	AttributeString("longname", &info->name );
	AttributeString("units",    &info->units );
	info->dimensionpoints.clear();
	return true;
}

/*-----------------------------------------------------------------------------
 *					nxNetcdfVar::LoadRankSpecs		2011-7-19*/
/** Loads the rank specifications of this variable**/
/*---------------------------------------------------------------------------*/

bool nxNetcdfVar::LoadDimInfo(  std::vector<nxNetcdfDimStruct>* dimsinfo ) const
{
	bool				ok;
	bool				ok1;
	int					ndims;
	std::vector<int>	dims;
	int					i;
	int					idx;
	const nxNetcdfGroup*		group;
	const nxNetcdfDim*			dim;		
	nxNetcdfDimStruct		dummy;

	ok =  (nc_inq_varndims (ParentNCID(), NCID(), &ndims) == NC_NOERR);								// get the number of dimensions
	ok =  ok && ( ndims > 0);																		// we shoudl have at least one dimension
	if (ok )																						// If we are good
	{																								// then see if we have a real array or a scalar
		dims.resize      (ndims+1);																	// allocate temporary (int) space to hold the size of each dimension
		dimsinfo->resize (ndims+1);																	// allocate space to hold the returned info
		ok    =  (nc_inq_vardimid ( ParentNCID(), NCID(), &dims.front() ) == NC_NOERR);				// get the dimension netcdf id codes of netcdf dimension entities
		group = ParentGroup();																		// get the parent group which we shall use to look up the dimension entities 
		ok    = ok && (group != NULL);																// make sure the parent group is sensible
		if (ok)																						// if it is
		{																							// then
			idx = ndims;
			for (i = 0; i < ndims; i++)																// for each dimension entity id
			{																						// we are going to find the dimension entity
				idx--;																				// but go backwards as LoadData returns last dim as the fastest/first dim
				dim  = group->DimAt( dims[idx] );															// the dimension entity
				ok1  = dim->GetDimensionInfo( &dimsinfo->at(i) );									// and get the dimension information into our structure
				ok   = ok && ok1;																	// check it is onk
			}																						// and do all of the dimensions
		}
		ok = ok && GetVarAttributes(  &dimsinfo->at(ndims) );
	}
	if (!ok )
	{
		dimsinfo->clear();
		nxLog::Record(NXLOG_WARNING,"nxNetcdfVar::LoadDimInfo, There was an error loading the dimension info for variable %s", (const char*)Name().c_str() );
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxNetcdfDim::GetDimensionInfo		2011-7-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfDim::GetDimensionInfo( nxNetcdfDimStruct* info ) const
{
	bool					ok;
	size_t					len;
	const nxNetcdfVar*		var;
	nx1dArray<double>		axisdata;
	std::string				unitstr;


	len   = Length();														// Get the length of this variable
	var   = RootGroup()->CoordinateVariable( Name().c_str(), false );		// Can we find a coordinate variable with the same name as this dimesnion
	ok    = (var != NULL);
	if (ok)
	{
		var->AttributeString("units", &unitstr);
		ok    = var->LoadData( &axisdata , HUGE_VAL );
		if (!ok) axisdata = axisdata.Indgen(len);
		ok    = ok && (axisdata.size() == len);
		if (ok)
		{
			info->name = Name();
			info->units = unitstr;
			info->dimensionpoints.assign( axisdata.begin(), axisdata.end() );
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING," nxNetcdfVar::GetDimensionInfo, error retrieving dimension info.");
		info->name = Name();
		info->units = "";
		info->dimensionpoints.clear();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					nxNetcdfGroup::nxNetcdfGroup		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

nxNetcdfGroup::nxNetcdfGroup()
{
	m_groups = new nxNetcdfGroups;

}


/*-----------------------------------------------------------------------------
 *					nxNetcdfGroup::nxNetcdfGroup		2011-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

nxNetcdfGroup::nxNetcdfGroup( const nxNetcdfGroup& other )
{
	nxNetcdfEntity::DeepCopy( other);
	m_groups  = new nxNetcdfGroups;
	*m_groups = *other.m_groups;
	m_vars    = other.m_vars;
	m_dims    = other.m_dims;
}



/*-----------------------------------------------------------------------------
 *					nxNetcdfGroup::~nxNetcdfGroup		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

nxNetcdfGroup::~nxNetcdfGroup()
{
	ReleaseResources();
	delete m_groups;
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfGroup:ReleaseResources		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

void nxNetcdfGroup::ReleaseResources()
{
	m_vars.Clear();
	m_dims.Clear();
	if (m_groups != NULL) m_groups->ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfGroup::At		2011-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

const nxNetcdfGroup* nxNetcdfGroup::At( const char* groupname) const
{
	return m_groups->At(groupname);
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfGroup::VarAt		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

const nxNetcdfVar* nxNetcdfGroup::VarAt( const char* varname, bool recursegroups) const
{
	const nxNetcdfVar*	varptr;

	varptr = m_vars.At( varname );
	if ((varptr == NULL) && (recursegroups ))
	{
		varptr = m_groups->VarAt( varname, recursegroups );
	}
	return varptr;
}

/*-----------------------------------------------------------------------------
 *					nxNetcdfGroup::Load		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfGroup::Load( )
{
	bool	ok;

	ok =       m_groups->Load       ( this );
	ok = ok && m_dims.Load          ( this );
	ok = ok && m_vars.Load          ( this );
	ok = ok && LoadVariableAttributes();
	ok = ok && LoadSubGroups();
	return ok;
}

	

/*-----------------------------------------------------------------------------
 *					nxNetcdfGroup::LoadVariableAttributes		2011-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfGroup::LoadVariableAttributes()
{
	nxNetcdfEntityArray<nxNetcdfVar>::iterator	variter;
	bool									ok = true;
	bool									ok1;

	for (variter = m_vars.Array().begin();!(variter == m_vars.Array().end()); ++variter)
	{
		ok1 = variter->second.LoadAttributes();
		ok = ok && ok1;
	}
	return ok;
}
/*-----------------------------------------------------------------------------
 *					nxNetcdfGroup::LoadVariableAttributes		2011-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfGroup::LoadSubGroups()
{
	nxNetcdfEntityArray<nxNetcdfGroup>::iterator	variter;
	bool									ok = true;
	bool									ok1;

	for (variter = m_groups->Groups().Array().begin();!(variter == m_groups->Groups().Array().end()); ++variter)
	{
		ok1 = variter->second.Load();
		ok = ok && ok1;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfGroup::CoordinateVariable		2011-7-19*/
/** **/
/*---------------------------------------------------------------------------*/

const nxNetcdfVar* nxNetcdfGroup::CoordinateVariable	( const char* varname, bool searchfromroot ) const
{
	const nxNetcdfVar*									var;
	nxNetcdfEntityArray<nxNetcdfGroup>::iterator	variter;

	if (searchfromroot)
	{
		var = RootGroup()->CoordinateVariable( varname, false );
	}
	else
	{
		var = VarAt(varname, false);
		if (var == NULL)
		{
			variter = m_groups->Groups().Array().begin();
			while ((var == NULL) && !(variter == m_groups->Groups().Array().end()) )
			{
				var = variter->second.CoordinateVariable( varname, false );
				++variter;
			}
		}
	}
	return var;
}




/*-----------------------------------------------------------------------------
 *					nxNetcdfGroups::Load		2011-7-15*/
/** Load all of the groups into an array of groups. This function is used for 
 *	Variables and dimensions
 **/
/*---------------------------------------------------------------------------*/

template <class nxNetcdfEntityType>
bool nxNetcdfEntityArray<nxNetcdfEntityType>::Load(  nxNetcdfGroup* parent )
{
	int												numvals;
	int												numvals2;
	std::vector<int>								id;
	char											name[NC_MAX_NAME+2];
	bool											ok;
	bool											ok1;
	size_t											i;
	iterator										iter;
	nxNetcdfEntityType								dummy;
	int												ncid;
	std::string										keyname;

	m_array.clear();
	dummy.SetParent( parent );
	ncid = dummy.ParentNCID();

	ok = (dummy.Handler_Inq_NumIDs( ncid, &numvals, NULL ) == NC_NOERR);			// Find out how many groups are children of this level
	if (ok && (numvals > 0))														// If we are good and we have entries
	{																				// then
		id.resize(numvals);															// resize the number of elements in the temporary array
		ok = (dummy.Handler_Inq_NumIDs( ncid, &numvals2, &id.front()) == NC_NOERR);	// And get the ncid of each child group
		ok = ok && (numvals2 == numvals);											// make sure we are still seeing the same number of groups
		if (ok)																		// if that worked ok.
		{																			// then we can create the new objects
			for (i = 0; i < (size_t) numvals; i++)									// so for each retrieved group
			{																		// get the group details
				ok1 = (dummy.Handler_Inq_Names(  ncid, id[i], name ) == NC_NOERR);	// from the netcdf file
				if (ok1)															// if that worked ok
				{																	// then
					dummy.SetNCID( id[i] );
					dummy.SetName( name );
					keyname = name;
					m_array.insert( value_type(keyname, dummy) );
				}
				ok = ok && ok1;
			}
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxNetcdfEntityArray<nxNetcdfEntityType>::Load, Error loading in entity. This is a problem as it probably indicates some sort of file corruption");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *							nxNetcdfEntityArray::LoadAttributes		2011-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

template <class nxNetcdfEntityType>
bool nxNetcdfEntityArray<nxNetcdfEntityType>::LoadAttributes(  nxNetcdfVar* varparent )
{
	int												numatts;
	char											name[NC_MAX_NAME+2];
	bool											ok;
	bool											ok1;
	int												i;
	nxNetcdfAtt										dummy;						// This code is opnly designed for instances of nxNetcdfAtt
	int												ncid;
	int												varid;
	std::string										keyname;

	m_array.clear();

	ncid  = varparent->ParentNCID();
	varid = varparent->NCID();
	ok = (nc_inq_varnatts( ncid, varid, &numatts ) == NC_NOERR);				// Find out how many groups are children of this level
	if (ok && (numatts > 0))														// If we are good and we have entries
	{																				// then
		for (i = 0; i <  numatts; i++)												// so for each retrieved group
		{																			// get the group details
			ok1 = (nc_inq_attname( ncid, varid, i , name) == NC_NOERR);
			if (ok1)																// if that worked ok
			{																		// then
				keyname = name;														// Copy the name over to a std::string
				dummy.SetParent(varparent);
				dummy.SetNCID( i );
				dummy.SetName( name );
				m_array.insert( value_type( keyname, dummy)  );
			}
			ok = ok && ok1;
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"nxNetcdfEntityArray<nxNetcdfEntityType>::LoadAttributes, Error loading in entity. This is a problem as it probably indicates some sort of file corruption");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfGroups::ReleaseResoucres		2011-7-18*/
/** **/
/*---------------------------------------------------------------------------*/

void nxNetcdfGroups::ReleaseResources( )
{
	m_groups.Clear();
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfGroups::VarAt		2011-7-25*/
/** recurse over all of the groups looking for the desired variable. **/
/*---------------------------------------------------------------------------*/

const nxNetcdfVar* nxNetcdfGroups::VarAt( const char* varname, bool recursegroups )
{
	nxNetcdfEntityArray<nxNetcdfGroup>::iterator	iter;
	const nxNetcdfVar* 									ptr = NULL;

	iter = m_groups.Array().begin();
	while ( (ptr == NULL) && !(iter ==  m_groups.Array().end()))
	{
		ptr = (*iter).second.VarAt(varname,recursegroups);
		++iter;
	}
	return ptr;
}

/*-----------------------------------------------------------------------------
 *					nxNetcdfDim::LoadData		2011-7-19*/
/** **/
/*---------------------------------------------------------------------------*/

size_t nxNetcdfDim::Length() const
{
	size_t	length;
	bool	ok;

	ok = (nc_inq_dimlen  ( ParentNCID(), NCID(), &length) == NC_NOERR);
	if (!ok) length = 0;
	return length;
}

/*-----------------------------------------------------------------------------
 *					nxNetcdfFile::nxNetcdfFile		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

nxNetcdfFile::nxNetcdfFile()
{
	m_isopen = false;
	m_topgroup.SetName("Root");
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfFile::~nxNetcdfFile		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

nxNetcdfFile::~nxNetcdfFile()
{
	Close();
	m_isopen = false;
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfFile::Close		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

void  nxNetcdfFile::Close()
{
	if (m_isopen)
	{
		nc_close(NCID());
		m_isopen = false;
		SetNCID(0);
		SetName("");
		m_topgroup.ReleaseResources();
	}
}


/*-----------------------------------------------------------------------------
 *					nxNetcdfFile::OpenRead		2011-7-15*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfFile::OpenRead( const char* filename )
{
	int		ncid;
	bool	ok;

	Close();						
	ok =  (nc_open( filename, NC_NOWRITE, &ncid) == NC_NOERR);                // open existing netCDF dataset
	if (ok)
	{
		m_isopen = true;
		SetNCID(ncid);
		SetName(filename);
		m_topgroup.SetNCID( ncid );
		ok = m_topgroup.Load( );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "nxNetcdfFile::OpenRead, Error reading metadata information from an existing NetCDF file. Thats not good");
		}
	}
	return ok;
}



/*---------------------------------------------------------------------------
 *          nxNetcdfHyperSlabDefn::nxNetcdfHyperSlabDefn          2019-10-03 */
/** **/
/*---------------------------------------------------------------------------*/
nxNetcdfHyperSlabDefn::nxNetcdfHyperSlabDefn()
{
	m_variable = nullptr;
}



/*---------------------------------------------------------------------------
 *            nxNetcdfHyperSlabDefn::AttachToVariable             2019-10-03 */
/** Attach this hyperslab object to the given variable and define a hyper slab
 *	that loads in the entire array
 **/
/*---------------------------------------------------------------------------*/
bool nxNetcdfHyperSlabDefn::AttachToVariable( const nxNetcdfVar* variable )
{
	bool ok;

	ok = (variable != nullptr);
	m_variable = variable;
	if (ok)
	{
		ok = m_variable->LoadRankSpecs( &m_rankspecs );
		if (ok)
		{
			size_t npts = m_rankspecs.size();
			m_start.resize (npts,0);
			m_count.resize (npts);
			m_stride.resize(npts,1);
			for (size_t i = 0; i < npts; i++)
			{
				m_count[i] = m_rankspecs[i];
			}
		}
		else
		{
			nxLog::Record(NXLOG_WARNING,"nxNetcdfHyperSlabDefn::AttachToVariable, Error loading the rank specifications from the variable.");
		}

	}
	if (!ok)
	{
		m_rankspecs.erase();
		m_start.resize(0);
		m_count.resize(0);
		m_stride.resize(0);
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *             nxNetcdfHyperSlabDefn::ClearSlabRange              2019-10-03 */
/** **/
/*---------------------------------------------------------------------------*/

void nxNetcdfHyperSlabDefn::ClearSlabRange()
{
	size_t	npts = m_rankspecs.size();

	m_start.assign(npts, 0);
	m_count.assign(npts, 0);
	m_stride.assign(npts,1);
}

/*---------------------------------------------------------------------------
 *              nxNetcdfHyperSlabDefn::DefineRanges               2019-10-03 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfHyperSlabDefn::DefineRanges( size_t npoints, size_t* start, size_t* count, size_t* stride)
{
	bool	ok;
	bool	ok1;
	size_t	n;
	size_t	n1;
	size_t	s;
	size_t	d;
	size_t	t;

	ok = m_rankspecs.size() == npoints;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "nxNetcdfHyperSlabDefn::DefineRanges, dimension size mismatch. Expected %d dimensions to be given but got %d.", (int) m_rankspecs.size(), (int)npoints);
	}
	else
	{
		n1 = npoints - 1;
		for (size_t i=0; i < npoints; i++)
		{
			n = m_rankspecs[i];
			s = start[i];
			d = count[i];
			if ( d == NXARRAY_STARSELECT)
			{
				d = m_rankspecs[i] - s;
			}
			t = (stride != nullptr)? stride[i] : 1;
			ok1 = ( s < n ) && ((s + d*t) <= n);
			ok = ok && ok1;
			if (ok1)
			{
				m_start[n1-i]  = s;
				m_count[n1-i]  = d;
				m_stride[n1-i] = t;
			}
			else
			{
				nxLog::Record(NXLOG_WARNING, "nxNetcdfHyperSlabDefn::DefineRanges, dimension (%d), requested start, count and stride (%d, %d, %d) do not match dimension size (%d).", (int)i, (int)s, (int)d, (int)t, (int)n);
			}
		}
	}
	if (!ok)
	{
		ClearSlabRange();
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *      nxNetcdfHyperSlabDefn::AllocateArrayToHoldHyperSlab       2019-10-03 */
/** **/
/*---------------------------------------------------------------------------*/

bool nxNetcdfHyperSlabDefn::AllocateArrayToHoldHyperSlab( nxArrayLinear<double>* slabarray) const
{
	bool ok;
	size_t ndims = m_rankspecs.size();
	size_t n1    = ndims - 1;
	std::vector<size_t> count;

	count.resize(ndims);
	for (size_t i = 0; i < ndims; i++)
	{
		count[i] = m_count[n1-i];
	}

	ok =       (ndims > 0);
	ok = ok && slabarray->SetSize( int(ndims), &count.front());
	if (!ok)
	{
		slabarray->erase();
	}
	return ok;
}
