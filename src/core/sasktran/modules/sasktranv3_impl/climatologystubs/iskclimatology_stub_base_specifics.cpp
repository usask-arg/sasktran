
#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_Constant::ISKClimatology_Stub_Constant		 2015- 6- 30*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_Constant::ISKClimatology_Stub_Constant( skClimatology_Constant* ptclimate)
	                         :ISKClimatology_Stub_Base(ptclimate)
{
	m_ptclimate = ptclimate;
	MakeSetFunctions();
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_Constant::~ISKClimatology_Stub_Constant		 2015- 6- 30*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_Constant::~ISKClimatology_Stub_Constant()
{
}
/*---------------------------------------------------------------------------
 *        ISKClimatology_Stub_Constant::MakeSetFunctions         2019-10-02 */
/** **/
/*---------------------------------------------------------------------------*/

void ISKClimatology_Stub_Constant::MakeSetFunctions()
{
	AddSetScalarFunction( 	"SetConstantValue", 
		[&, this](double d)
		{
			m_ptclimate->SetConstantValue(d);
			return true;
		}
	);
}

/*-----------------------------------------------------------------------------
 *					:ISKClimatology_Stub_Base		 2015- 1- 9*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_OnePressureTemp::ISKClimatology_Stub_OnePressureTemp( skClimatology_OneTemperatureAndPressure* ptclimate)
	                                :ISKClimatology_Stub_Base(ptclimate)
{
	m_ptclimate = ptclimate;
	MakeSetFunctions			();
}


/*---------------------------------------------------------------------------
 *     ISKClimatology_Stub_OnePressureTemp::MakeSetFunctions     2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

void ISKClimatology_Stub_OnePressureTemp::MakeSetFunctions()
{
	AddSetScalarFunction("SetTemperature",
		[&,this](double value)
		{
			m_ptclimate->SetTemperature(value);
			return true;
		}
	);
	AddSetScalarFunction( "SetPressure",
		[&,this](double value)
		{
			m_ptclimate->SetPressure( value);
			return true;
		}
	);

}


/*---------------------------------------------------------------------------
 *              ~ISKClimatology_Stub_OnePressureTemp              2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_OnePressureTemp:: ~ISKClimatology_Stub_OnePressureTemp()
{
}

/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefined		2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_UserDefined::ISKClimatology_Stub_UserDefined( skClimatology_UserTableSpline* climate)
	                            :ISKClimatology_Stub_Base(climate)
{
	m_userdefinedclimatology = climate;				// Note that the reference count is in the base class
	m_dologinterpolation     = false;
	m_dopiecewiselinear      = false;
	m_currentbadvalue        = 0.0;
	MakeSetFunctions();

}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefined::~ISKClimatology_Stub_UserDefined		2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_UserDefined::~ISKClimatology_Stub_UserDefined	()
{
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefined::Set		2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology_Stub_UserDefined::SetPropertyUserDefined( const CLIMATOLOGY_HANDLE& species, double* profile, int numpoints )
{
	std::vector<double>	p;
	bool				ok;

	ok = (m_currentheightarray.size() == (size_t)numpoints);
	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING,"ISKClimatology_Stub_UserDefined::SetPropertyUserDefined, the number of point in the height array (%z) set in a previous call to SetPropertyArray(\"Heights\", ...) does not equal the number of points in the profile (%d)", (size_t)m_currentheightarray.size(), (int)numpoints);
	}
	else
	{
		p.assign( profile,  profile + numpoints);
		ok = m_userdefinedclimatology->AddProfile( species, m_currentheightarray, p, m_dologinterpolation, m_dopiecewiselinear, m_currentbadvalue);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefined::SetPropertyScalar		2014-4-16*/
/** **/
/*---------------------------------------------------------------------------*/

void ISKClimatology_Stub_UserDefined::MakeSetFunctions()
{
	AddSetScalarFunction("dologinterpolation",
		[&,this](double value)
		{
			m_dologinterpolation = (value != 0.0);
			return true;
		}
	);

	AddSetScalarFunction( "dopiecewiselinear",
		[&, this] (double value)
		{
			m_dopiecewiselinear = (value != 0.0);
			return true;
		}
	);

	AddSetScalarFunction( "badvalue",
		[&, this] (double value)
		{
			m_currentbadvalue = value;
			return true;
		}
	);

	AddSetVectorFunction( "heights", 
		[&, this] ( const double* value, int numpoints)
		{
			m_currentheightarray.assign( value, value+numpoints);
			return true;
		}
	);

}



/*-----------------------------------------------------------------------------
 *					:ISKClimatology_Stub_UserDefinedTable		 2016- 9- 26*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_UserDefinedTable::ISKClimatology_Stub_UserDefinedTable( skClimatology_UserDefinedTable* climate)
	                                 :ISKClimatology_Stub_Base(climate)
{
	m_userdefinedclim = climate;
	m_currenthandle = SKCLIMATOLOGY_UNDEFINED;
	MakeSetFunctions();
}

/*---------------------------------------------------------------------------
 *    ISKClimatology_Stub_UserDefinedTable::MakeSetFunctions     2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

void ISKClimatology_Stub_UserDefinedTable::MakeSetFunctions()
{
	AddSetVectorFunction( "heights",
		[&, this](const double* value, int numpoints)
		{
			m_heights.assign( value, value+numpoints);
			return true;
		}
	);
	
	AddSetStringFunction( "fromtextfile",
		[&, this](const char *value )
		{
			bool ok;
			ok = !(m_currenthandle == SKCLIMATOLOGY_UNDEFINED);		// Note the != operator is not supported on Linux systems
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"ISKClimatology_Stub_UserDefined::SetPropertyString, cannot lad profile frm file as the current climatology handle is Undefined. Set Property \'climatology_handle\' first");
			}
			else
			{
				const char* filename = value;
				ok = m_userdefinedclim->LoadProfileFromTextFile( &m_currenthandle, 1, filename);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"USERDEFINED_TABLE, There were errors loading a climatology profile from file <%s>", (const char*)filename);
				}
			}
			return ok;
		}
	);

	AddSetStringFunction( "climatology_handle",
		[&, this](const char *value )
		{
			bool				ok;
			CLIMATOLOGY_HANDLE*	handleptr;
				
			handleptr = FindGlobalClimatologyHandle ( value );					// MAke sure this handle already exists. If it does not then ** IN THE FUTURE *** we should perhaps create the entry
			ok  = (handleptr != nullptr);
			m_currenthandle = ok ? *handleptr : SKCLIMATOLOGY_UNDEFINED;
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"USERDEFINED_TABLE, SetProperty(CLIMATOLOGY_HANDLE) failed as <%s> was not a recognised existing handle.", (const char*)value);
			}
			return ok;
		}
	);
}



/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefinedTable::SetPropertyUserDefined		 2016- 9- 26*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology_Stub_UserDefinedTable::SetPropertyUserDefined( const CLIMATOLOGY_HANDLE& species, double* profile, int numpoints )
{
	bool ok = true;

	nx2dArray<double> prof;

	prof.SetSize( numpoints, 2 );
	for( size_t idx = 0; (int) idx < numpoints; idx++ )
	{
		prof.At( idx, 0 ) = m_heights[idx];
		prof.At( idx, 1 ) = profile[idx];
	}
	m_userdefinedclim->LoadProfileFrom2dArray( &species, 1, prof );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefined		2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_UserDefined3D::ISKClimatology_Stub_UserDefined3D( skClimatology_UserDefined3D_LatLonHeight* climate)
	                              :ISKClimatology_Stub_Base(climate)
{
	m_userdefinedclimatology = climate;				// Note that the reference count is in the base class
	m_currentbadvalue        = std::numeric_limits<double>::quiet_NaN();
	MakeSetFunctions();
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefined::~ISKClimatology_Stub_UserDefined		2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_UserDefined3D::~ISKClimatology_Stub_UserDefined3D	()
{
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefined::Set		2014-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKClimatology_Stub_UserDefined3D::SetPropertyUserDefined( const CLIMATOLOGY_HANDLE& species, double* profile, int numpoints )
{
	nx3dArray<double>	p;
	bool				ok;
	size_t				np;

	np = m_currentheightarray.size() * m_currentlatarray.size() * m_currentlonarray.size();
	ok =    (np == (size_t)numpoints);

	if (!ok)
	{
		nxLog::Record( NXLOG_WARNING,"ISKClimatology_Stub_UserDefined3D::SetPropertyUserDefined, the total number of points in the height, latitude and longitude arrays (%z,%z %z = %z) set in a previous call to SetPropertyArraydoes not equal the number of points in this profile (%d)",
			                          (size_t)m_currentheightarray.size(), (size_t) m_currentlatarray.size(), (size_t)m_currentlonarray.size(), (size_t)np, (int)numpoints);
	}
	else
	{
		ok =       p.Attach( m_currentheightarray.size(),  m_currentlonarray.size(), m_currentlatarray.size(), profile);
		ok = ok && m_userdefinedclimatology->LoadProfile( species, m_currentheightarray, m_currentlonarray, m_currentlatarray, p, m_currentbadvalue);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_UserDefined::SetPropertyArray		2014-4-16*/
/** **/
/*---------------------------------------------------------------------------*/

void ISKClimatology_Stub_UserDefined3D::MakeSetFunctions( )
{

	AddSetVectorFunction( "heights", 
		[&, this](const double* value, int numpoints)
		{
			m_currentheightarray.assign( value, value+numpoints);
			return true;
		}
	);

	AddSetVectorFunction( "latitudes",
		[&, this](const double* value, int numpoints)
		{
			m_currentlatarray.assign( value, value+numpoints);
			return true;
		}
	);

	AddSetVectorFunction( "longitudes",
		[&, this](const double* value, int numpoints)
		{
			m_currentlonarray.assign( value, value+numpoints);
			return true;
		}
	);
	
	AddSetScalarFunction( "badvalue",
		[&,this] (double value)
		{
			m_currentbadvalue = value;
			return true;
		}
	);
}










/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600::Constructor		2014-3-31*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600::ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600( skClimatology_OsirisAerosolModeRadiusV600* climate)
													:ISKClimatology_Stub_Base(climate)
{
	m_climatologyv600 = climate;
	MakeSetFunctions();
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600:::ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600		2014-3-31*/
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600::~ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600()
{
}


/*-----------------------------------------------------------------------------
 *					ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600::MakeSetFunctions		2014-3-31*/
/** **/
/*---------------------------------------------------------------------------*/

void ISKClimatology_Stub_OSIRISL2_AEROSOLMODERADIUS_V600::MakeSetFunctions( )
{

	AddSetScalarFunction("setisascendingnode",
		[&,this](double value)
		{
			bool setascending;
			
			setascending = (value != 0.0);
			m_climatologyv600->SetIsAscendingNode(setascending);
			return true;
		}
	);
}

/*---------------------------------------------------------------------------
 * ISKClimatology_Stub_UserDefinedPlane::ISKClimatology_Stub_UserDefinedPlane2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

ISKClimatology_Stub_UserDefinedPlane::ISKClimatology_Stub_UserDefinedPlane( skClimatology_UserDefinedPlane* userclim )
	: ISKClimatology_Stub_Base( userclim )
{
	m_userdefinedclim = userclim;
	m_dolog = false;
	MakeSetFunctions();
}

/*---------------------------------------------------------------------------
 *    ISKClimatology_Stub_UserDefinedPlane::MakeSetFunctions     2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

void ISKClimatology_Stub_UserDefinedPlane::MakeSetFunctions()
{

	AddSetScalarFunction( "dologinterpolation",
		[&,this](double value )
		{
			m_dolog = value > 0.5;
			return true;
		}
	);

	AddSetVectorFunction( "heights",
		[&,this] (const double* value, int numpoints)
		{
			bool ok;

			ok = m_userdefinedclim->SetHeightGrid( std::vector<double>( value, value+numpoints) );
			m_numheights = numpoints;
			return ok;
		}
	);
	
	AddSetVectorFunction( "angles",
		[&,this] (const double* value, int numpoints)
		{
			bool ok;
			ok = m_userdefinedclim->SetAngleGrid( std::vector<double>( value, value+numpoints) );
			m_numangles = numpoints;
			return ok;
		}
	);

	AddSetVectorFunction( "normalandreference",
		[&,this] (const double* value, int numpoints)
		{
			bool ok = (numpoints == 6);
			if (ok)
			{
				nxVector normal   (value[0], value[1], value[2]);
				nxVector reference(value[3], value[4], value[5]);
				ok = m_userdefinedclim->SetPlane( normal, reference );
			}
			if (!ok)
			{
				nxLog::Record( NXLOG_WARNING,"ISKClimatology_Stub_UserDefinedPlane::there were errors executing property  normalandreference. Did you pass in exactly 6 points (%d) ?", (int)numpoints);
			}
			return ok;
		}
	);

}


/*---------------------------------------------------------------------------
 *  ISKClimatology_Stub_UserDefinedPlane::SetPropertyUserDefined  2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/
bool ISKClimatology_Stub_UserDefinedPlane::SetPropertyUserDefined( const CLIMATOLOGY_HANDLE& species,  double* profile, int numheights )
{
	bool ok;

	ok = (m_numheights*m_numangles == numheights );
	if( ok )
	{
		nx2dArray<double> prof(m_numheights, m_numangles, profile);
		m_userdefinedclim->AddSpecies( prof, species, m_dolog );
	}
	else
	{
		nxLog::Record(NXLOG_WARNING, "Error setting user defined property in ISKClimatology_Stub_UserDefinedPlane");
	}

	return ok;
}

