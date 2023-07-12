#include <skopticalproperties21.h>
#include "skbrdf_stubs.h"


/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Base::ISKBrdf_Stub_Base		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Base::ISKBrdf_Stub_Base	( skBRDF* optprop)
{
	m_brdf = optprop;
	if (m_brdf != NULL) m_brdf->AddRef();
}


/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Base::~ISKBrdf_Stub_Base		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

ISKBrdf_Stub_Base::~ISKBrdf_Stub_Base	()
{
	if (m_brdf != NULL) m_brdf->Release();
}


/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Base::SkOpticalPropertyPointer		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Base::CalculateCrossSectionsArray		2014-3-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_Base::BRDF( double wavelennm, const GEODETIC_INSTANT& pt, double MU_in, double MU_out, double COSDPHI, double* brdf)
{
	bool	ok = (m_brdf != nullptr);
	ok = ok && m_brdf->BRDF( wavelennm, pt, MU_in, MU_out, COSDPHI, brdf);
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Base::SetPropertyScalar		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_Base::SetPropertyScalar( const char* propertyname, double value )
{
	nxLog::Record(NXLOG_WARNING,"ISKBrdf_Stub_Base::SetPropertyScalar, this optical property does not support scalar property [%s]\n", (const char*)propertyname);
	return false;
}

/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Base::SetPropertyArray		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_Base::SetPropertyArray( const char* propertyname, const double* value, int numpoints )
{
	nxLog::Record(NXLOG_WARNING,"ISKBrdf_Stub_Base::SetPropertyArray, this optical property does not support array property [%s]\n", (const char*)propertyname);
	return false;
}


/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Base::SetPropertyObject		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_Base::SetPropertyObject( const char* propertyname, nxUnknown* object )
{
	nxLog::Record(NXLOG_WARNING,"ISKBrdf_Stub_Base::SetPropertyObject, this optical property does not support object property [%s]\n", (const char*)propertyname);
	return false;
}



/*-----------------------------------------------------------------------------
 *					ISKBrdf_Stub_Base::SetPropertyString		 2016- 9- 26*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKBrdf_Stub_Base::SetPropertyString( const char* propertyname, const char* str)

{
	nxLog::Record(NXLOG_WARNING,"ISKBrdf_Stub_Base::SetPropertyString, this optical property does not support string property [%s]\n", (const char*)propertyname);
	return false;
}

