#include <skclimatology21.h>
#include "sknetcdf4.h"

/*-----------------------------------------------------------------------------
 *					sknetcdfVariable::sknetcdfVariable		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

sknetcdfVariable::sknetcdfVariable()
{
	m_data = NULL;
}

/*-----------------------------------------------------------------------------
 *					sknetcdfVariable::sknetcdfVariable		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

sknetcdfVariable::sknetcdfVariable( const char* name)
{
	m_data = NULL;
	m_name = name;
}

/*-----------------------------------------------------------------------------
 *					sknetcdfVariable::sknetcdfVariable		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

sknetcdfVariable::sknetcdfVariable( const sknetcdfVariable& other )
{
	m_name     = other.m_name;
	m_diminfo  = other.m_diminfo;
	if (other.m_data == NULL)
	{
		m_data = NULL;
	}
	else
	{
		m_data = new nxArrayLinear<double>;
		m_data->DeepCopy( *other.m_data);
	}
}


/*-----------------------------------------------------------------------------
 *					sknetcdfVariable::~sknetcdfVariable		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

sknetcdfVariable::~sknetcdfVariable()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					sknetcdfVariable::ReleaseResources		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

void sknetcdfVariable::ReleaseResources()
{
	m_diminfo.clear();
	if (m_data != NULL)
	{
		delete m_data;
		m_data = NULL;
	}
}
		


/*-----------------------------------------------------------------------------
 *					sknetcdfVariable::Dimension		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

const nxNetcdfDimStruct* sknetcdfVariable::Dimension( const char* name) const
{
	size_t							i;
	const nxNetcdfDimStruct*		ptr = NULL;
	
	for ( i = 0; (( ptr == NULL) && (i < m_diminfo.size())); i++)
	{
		if (m_diminfo.at(i).name == name)
		{
			ptr  = &(m_diminfo.at(i));
		}
	}
	return ptr;
}


/*-----------------------------------------------------------------------------
 *					sknetcdfVariable::Dimension		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

const nxNetcdfDimStruct* sknetcdfVariable::Dimension( size_t idx ) const
{
	bool	ok;

	ok = (idx < m_diminfo.size() );
	return ok ? &(m_diminfo.at(idx)) : NULL;
}


/*-----------------------------------------------------------------------------
 *					sknetcdfVariable::LoadFromNetCDF		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool sknetcdfVariable::LoadFromNetCDF( nxNetcdfFile& ncfile, double missingvalue )
{
	const nxNetcdfVar*	var;
	bool				ok;

	ReleaseResources();
	ok = (m_name.size() > 0);
	if (ok)
	{
		var = ncfile.VarAt( m_name.c_str() );
		ok = (var != NULL);
		ok = ok && var->LoadDimInfo( &m_diminfo );
		if (ok)
		{
			switch (m_diminfo.size())
			{
				case 2 :	m_data = new nx1dArray<double>; break;
				case 3 :	m_data = new nx2dArray<double>; break;
				case 4 :	m_data = new nx3dArray<double>; break;
				default:	m_data = new nxArrayLinear<double>; break;
			};
			ok = ok && (m_data != NULL);
			ok = ok && var->LoadData( m_data, missingvalue );
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sknetcdfVariable::LoadFromNetCDF, error loading variable <%s> from netcdf file. Thats probably not good", (const char*)m_name.c_str());
		ReleaseResources();
	}
	return ok;
}







/*-----------------------------------------------------------------------------
 *					sknetcdf_IONetCDF::sknetcdf_IONetCDF		2011-7-21*/
/** **/
/*---------------------------------------------------------------------------*/

sknetcdf_IONetCDF::sknetcdf_IONetCDF()
{
	m_missingvalue = -99999;
}



/*-----------------------------------------------------------------------------
 *					sknetcdf_IONetCDF::~sknetcdf_IONetCDF		2011-7-21*/
/** **/
/*---------------------------------------------------------------------------*/

sknetcdf_IONetCDF::~sknetcdf_IONetCDF()
{
}


/*-----------------------------------------------------------------------------
 *					sknetcdf_IONetCDF::AddEmptyField		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool sknetcdf_IONetCDF::AddEmptyField( const char* name)
{
	sknetcdfVariable	empty(name);

	m_fields.insert(empty);
	return true;
}


/*-----------------------------------------------------------------------------
 *					sknetcdf_IONetCDF::Field		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

const sknetcdfVariable*	sknetcdf_IONetCDF::Field(const char* name )
{
	sknetcdfVariable			dummy(name);
	const sknetcdfVariable*		ptr;
	fielditerator			iter;

	iter = m_fields.find( dummy );
	ptr  = (iter != m_fields.end()) ? &(*iter) : NULL;
	return ptr;
}


/*-----------------------------------------------------------------------------
 *					sknetcdf_IONetCDF::LoadFields		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

bool sknetcdf_IONetCDF::LoadFields( const char* filename)
{
	fielditerator			iter;
	const sknetcdfVariable*	constptr;
	sknetcdfVariable*		ptr;
	bool					ok = true;
	bool					ok1;
	nxNetcdfFile			ncfile;

	ok = ncfile.OpenRead( filename );
	if (ok)
	{
		for (iter = m_fields.begin(); !( iter == m_fields.end()); ++iter)
		{
			constptr = &(*iter);									// ndl303, 2012-10-04. We are defetaing typechecking in the compiler as the code was originally written
			ptr      = (sknetcdfVariable*)((intptr_t)constptr);		// to work when <set> members were not forced to be constant. gcc in 2012 is demanding 
			ok1 = ptr->LoadFromNetCDF( ncfile, m_missingvalue);		// they be constant. The only other option other than this hack is a lot of reworking
			ok = ok && ok1;											// Knock yourself out if you object to the hack.
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					sknetcdf_IONetCDF::LoadFile		2011-7-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool sknetcdf_IONetCDF::LoadFile( const char* filename, const char* fieldlist, double missingvalue)
{
	bool				ok = true;
	nxStringArray		tokens;
	
	m_missingvalue = missingvalue;
	nxStrtok( fieldlist, &tokens, ";,:" );

	m_fields.clear();
	
	for (int i = 0; i <tokens.GetSize(); i++)
	{
		ok = ok && AddEmptyField( (const char*)(tokens.GetAt(i)) );
	}
	ok = ok && LoadFields  ( filename );
	return ok;
}




/*-----------------------------------------------------------------------------
 *					FindBoundingIndices		2011-7-25*/
/** **/
/*---------------------------------------------------------------------------*/

/*
static bool FindBoundingIndicesAscending( const std::vector<double>& myarray, double x, size_t* lowercell, size_t* uppercell )
{
	std::vector<double>::const_iterator		x1;							// The value just above our value of x
	std::vector<double>::const_iterator		x0;							// The value just below our value of x
	std::vector<double>::const_iterator		start;
	std::vector<double>::const_iterator		finish;
	bool									ok;
	bool									outofbounds;

	start  = myarray.begin();										// get the start of the altitude grid
	finish = myarray.end();											// Get the end of the altitude grid
	x1     = std::upper_bound( start, finish, x );				// Find the pointer to the value greater than x

	outofbounds = (x1 == start) || (x1 == finish);			// this is OK as long as it is not out of bounds
	if (outofbounds)										// if it is out of bounds
	{														// then
		if (x1 == start  ) x1 = start  + 1;					// adjust values before the start of the table
		if (x1 == finish ) x1 = finish - 1;					// adjust values beyond the end of the table
	}														// we now havethe upper value as sensible;
	x0           = x1 - 1;									// get the lower value
	*uppercell   = (x1 - start);							// Get the index of the upper cell
	*lowercell   = (x0 - start);							// Get the index of the lower cell

	return ok;													// return ok.
}
*/

/*-----------------------------------------------------------------------------
 *					skEcmwf_IO::InterpolateLatLongAllPressures		2005-6-23*/
/** Interpolates the primary value (Temperature or geopotential height) to the
 *	specified latitiude and longitude usning linear interpolation of the 4 lat/long
 *	points enclosing the required place. The interpolation is made at all
 *	of the pressure levels.
 **/
/*---------------------------------------------------------------------------*/

/*

bool sknetcdf_IONetCDF::VariableProfile_InterpolateLatLong	( const char* varname, double latitude, double longitude, std::vector<double>* profile)
{
	const sknetcdfVariable*		var;
	const nxNetcdfDimStruct*	latdim;
	const nxNetcdfDimStruct*	lngdim;
	const nxNetcdfDimStruct*	levdim;

	size_t						iy0, iy1; 
	size_t						ix0, ix1;
	double						y0,y1,y;
	double						x0,x1,x;
	double						vs[4];
	bool						ok;
	bool						ok1;
	size_t						pidx;
	size_t						nlevels;
	const nx3dArray<double>*	data;


	var = Field(varname );
	ok  = (var != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"sknetcdf_IONetCDF::VariableProfile_InterpolateLatLong, cannot find variable <%s>. Thats not good", (const char*)varname );
	}
	else
	{	
		lngdim = var->Dimension( (size_t) 0);										// get the longitude dimension. It should be the second
		latdim = var->Dimension(1);										// get the latitude dimension. It should be the third.
		levdim = var->Dimension(2);										// Get level dimension. It should be the first

		ok =    (latdim != NULL)										// Make sure we got all of the dimension info 
			 && (lngdim != NULL)
			 && (levdim != NULL)
			 && (levdim->name == "level")								// and make sure the dimension info makes sense.
			 && (lngdim->name == "lon")									// and is what we expect.
			 && (latdim->name == "lat");

		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"sknetcdf_IONetCDF::VariableProfile_InterpolateLatLong, Error finidng the lat, lon or levels dimension for variable <%s>, This is not going to work.", (const char*)varname );
		}
		else
		{
			if (longitude < 0.0) longitude += 360.0;							// get Longitude in range 0-360.0;
			NXASSERT(( (longitude >= 0.0)   && (longitude <= 360.0) ));			// Check input variables are in range
			NXASSERT(( (latitude  >= -90.0) && (latitude  <= 90.0)  ));

			y   = latitude;
			x   = nxmath::inrange(longitude, 360.0);

			nlevels = levdim->dimensionpoints.size();
			profile->resize( nlevels );
			data    = var->Data3D();
			ok      =       FindBoundingIndicesDescending     ( latdim->dimensionpoints, y,          &iy0, &iy1, &y0, &y1 );		// Interpolate latitude
			ok      = ok && FindBoundingIndicesAscendingCyclic( lngdim->dimensionpoints, x,  360.0,  &ix0, &ix1, &x0, &x1 );		// Interpolate longitude
			ok      = ok && (profile->size() == nlevels );
			ok      = ok && (data != NULL);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"VariableProfile_InterpolateLatLong, There was an error interpolating in latitude and longitude. Thats rather odd");
			}
			else
			{
				for (pidx = 0; pidx < nlevels; pidx++)
				{
					vs[0] = data->At( ix0, iy0, pidx );
					vs[1] = data->At( ix0, iy1, pidx );
					vs[2] = data->At( ix1, iy1, pidx );
					vs[3] = data->At( ix1, iy0, pidx );

					ok1 =	   (vs[0] != m_missingvalue )
							&& (vs[1] != m_missingvalue )
							&& (vs[2] != m_missingvalue )
							&& (vs[3] != m_missingvalue );
					if (!ok1)
					{
						profile->at(pidx) = m_missingvalue;
					}
					else
					{
						profile->at(pidx) = nxLinearInterpolate::FromSquare( x,y,x0,x1,y0,y1, vs);
					}
				}
			}
		}
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"VariableProfile_InterpolateLatLong, There was an error generating the requested profile. It is set to empty");
		profile->clear();
	}
	return ok;
}

*/
