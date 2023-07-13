#include "../sktran_common.h"
#include <boost/algorithm/string.hpp>

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::ISKEngine_Stub_HR		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

ISKModuleBase_Stub::ISKModuleBase_Stub()
{
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::~ISKEngine_Stub_HR		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

ISKModuleBase_Stub::~ISKModuleBase_Stub()
{
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::ParseCommandAndIndex		 2015- 10- 29*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase_Stub::ParseCommandAndIndex( const char* inputchars,  std::string& cmdstr, int& index )
{
	nxStringArray	tokens;
	nxString		input( inputchars);
	nxString		cmd;
	int				numtoken;
	bool			ok;

	numtoken = nxStrtok( input, &tokens, "([]) ,:;");
	ok = (numtoken == 2);
	if (ok)
	{
		cmd = tokens.GetAt(0);
		index = atoi( tokens.GetAt(1) );
		cmd.MakeLower();
	}



	//// first check if we have an input with an index
	//std::regex e("^([^\\[\\(]*)[\\[\\(](\\d*)[\\]\\)]");			// i feel bad
	//std::smatch m;
	//std::string str = (const char*)input;

	//bool ok = true;
	//ok = ok && std::regex_search( str, m, e );
	//if( ok )
	//{
	//	// we have an input with an index
	//	ok = (m.size() == 3);
	//	if ( ok )
	//	{
	//		cmd = m[1].str().c_str();
	//		index = atoi(m[2].str().c_str());
	//	}
	//}
	if( !ok )
	{
		cmd = input;
		cmd.MakeLower();
		index = -1;
	}
	cmdstr = (const char*)cmd;
	return true;
}


/*---------------------------------------------------------------------------
 *            ISKModuleBase_Stub::AddSetScalarFunction            2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase_Stub::AddSetScalarFunction		( const char* aname, std::function<bool(double)> func)
{
	std::string	name(aname);
	boost::algorithm::to_lower(name);
	m_scalarsetfunctions[name] = func;
	return true;
}


/*---------------------------------------------------------------------------
 *            ISKModuleBase_Stub::AddSetStringFunction            2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase_Stub::AddSetStringFunction( const char * aname, std::function<bool(const char*)> func)
{
	std::string	name(aname);
	boost::algorithm::to_lower(name);
	m_stringsetfunctions[name] = func;
	return true;
}


/*---------------------------------------------------------------------------
 *            ISKModuleBase_Stub::AddSetVectorFunction            2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase_Stub::AddSetVectorFunction( const char* aname, std::function<bool(const double*, int)> func)
{
	std::string	name(aname);
	boost::algorithm::to_lower(name);
	m_vectorsetfunctions[name] = func;
	return true;
}


/*---------------------------------------------------------------------------
 *            ISKModuleBase_Stub::AddSetObjectFunction            2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/
bool ISKModuleBase_Stub::AddSetObjectFunction		( const char* aname, std::function<bool(nxUnknown*)> func)
{
	std::string	name(aname);
	boost::algorithm::to_lower(name);
	m_objectsetfunctions[name] = func;
	return true;
}


/*---------------------------------------------------------------------------
 *            ISKModuleBase_Stub::AddGetScalarFunction            2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/
bool ISKModuleBase_Stub::AddGetScalarFunction		( const char* aname, std::function<bool(double*)> func)
{
	std::string	name(aname);
	boost::algorithm::to_lower(name);
	m_scalargetfunctions[ name] = func;
	return true;
}


/*---------------------------------------------------------------------------
 *            ISKModuleBase_Stub::AddGetVectorFunction            2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase_Stub::AddGetVectorFunction		( const char* aname, std::function<bool(int)> func)
{
	std::string	name(aname);
	boost::algorithm::to_lower(name);
	m_vectorgetfunctions[ name ] = func;
	return true;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::Set		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase_Stub::SetPropertyScalar( const char* propertyname, double value )
{
	std::string		str(propertyname);
	boost::algorithm::to_lower(str);

	auto funciterator = m_scalarsetfunctions.find( str );
	if( funciterator == std::end(m_scalarsetfunctions) )
	{
		// 2015-11-20 ndl303, We dont want Set Scalar reporting an error as matlab checks it when passed 1 element arrays. We need it to fail quietly.
		//nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_HR::Set, this object does not support any scalar properties including [%s]\n", (const char*)propertyname); 
		return false;
	}
	else
	{
		return funciterator->second(value);
	}
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::Set		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase_Stub::SetPropertyString( const char* propertyname, const char* value )
{
	std::string		str(propertyname);
	boost::algorithm::to_lower(str);

	auto funciterator = m_stringsetfunctions.find( str );
	if( funciterator == std::end(m_stringsetfunctions) )
	{
		return false;
	}
	else
	{
		return funciterator->second(value);
	}
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::Set		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase_Stub::SetPropertyArray( const char* propertyname, const double* value, int numpoints )
{
	bool		ok;
	std::string str(propertyname);
	boost::algorithm::to_lower(str);

	auto funciterator = m_vectorsetfunctions.find( str );
	if( funciterator == std::end( m_vectorsetfunctions ) )
	{
		nxLog::Record(NXLOG_WARNING,"ISKModuleBase_Stub::SetPropertyArray, this object does not support array property [%s]\n", (const char*)str.c_str());
		return false;
	}
	else
	{
		return funciterator->second(value, numpoints);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::SetPropertyObject		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase_Stub::SetPropertyObject( const char* propertyname, nxUnknown* object )
{
	std::string str(propertyname);
	boost::algorithm::to_lower(str);

	auto funciterator = m_objectsetfunctions.find( str );
	if( funciterator == std::end(m_objectsetfunctions) )
	{
		nxLog::Record(NXLOG_WARNING,"ISKModuleBase_Stub::SetPropertyObject, this object does not support any object properties including [%s]\n", (const char*)str.c_str());
		return false;
	}
	else
	{
		return funciterator->second(object);
	};
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::GetPropertyScalar		 2014- 5- 16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase_Stub::GetPropertyScalar( const char* propertyname, double* value )
{
	std::string str(propertyname);
	boost::algorithm::to_lower(str);

	auto funciterator = m_scalargetfunctions.find( str );
	if( funciterator == std::end(m_scalargetfunctions) )
	{
		return false;
	}
	else
	{
		return funciterator->second(value);
	}
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::GetPropertyArray		 2014- 5- 16*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKModuleBase_Stub::GetProperty( const char* propertyname, const double** value, int* numpoints )
{
	std::string cmd;
	int index;
	double	scalarvalue;
	bool ok;

	ok = GetPropertyScalar( propertyname, &scalarvalue);
	if (ok)
	{
		m_getpropertybuffer.resize(1, scalarvalue);
		*numpoints = 0;
		*value = &m_getpropertybuffer[0];
	}
	else
	{

		ok = ParseCommandAndIndex( propertyname, cmd, index );
		auto funcinterator = m_vectorgetfunctions.find( cmd );
		if( funcinterator == std::end(m_vectorgetfunctions ) )
		{
			nxLog::Record(NXLOG_WARNING,"Sasktran Component, The HR engine not support property <%s>", (const char*) cmd.c_str() );
			*numpoints = 0;
			*value     = NULL;
			ok = false;
		}
		else
		{
			ok = funcinterator->second(index);
			*numpoints = (int)m_getpropertybuffer.size();
			*value = &m_getpropertybuffer[0];
		}
	}
	return ok;
}

