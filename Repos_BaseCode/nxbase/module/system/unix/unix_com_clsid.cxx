/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

#include "nxbase_core.h"
//#if defined(NX_UNIX_VER)

#include "unix_com_emmulate.h"

//---------------------------------------------------------------------------
//						nxUnixCLSID::constructor
//	Open up the CLSID initialization file (and read it into memory).
//---------------------------------------------------------------------------

nxUnixCLSID::nxUnixCLSID()
{
}

nxUnixCLSID::~nxUnixCLSID()
{
}

//---------------------------------------------------------------------------
//						nxUnixCLSID::KeyAsString
//	Look up the CLSID key for the given ID
//---------------------------------------------------------------------------

bool nxUnixCLSID::KeyAsString( REFCLSID cls, nxString* str )
{
	str->sprintf( "%.8X-%.4X-%.4X-%.2X%.2x-%.2X%.2X%.2X%.2X%.2X%.2X",
				 (unsigned int)cls.Data1,
			     (unsigned int)cls.Data2,
				 (unsigned int)cls.Data3,
				 (unsigned int)cls.Data4[0],
				 (unsigned int)cls.Data4[1],
				 (unsigned int)cls.Data4[2],
				 (unsigned int)cls.Data4[3],
				 (unsigned int)cls.Data4[4],
				 (unsigned int)cls.Data4[5],
				 (unsigned int)cls.Data4[6],
				 (unsigned int)cls.Data4[7]);
	return true;
}

static const GUID CLSID_OnyxOsirisname   = { 0x80707760, 0x89D6, 0x11D0, {0xB7,0x5B,0x00,0x00,0xC0,0x54,0x85,0x54}};
static const GUID CLSID_OnyxCoreDatabase = { 0xdc8e98b0, 0x49c2, 0x11d2, {0xb8,0x76,0x00,0x00,0xc0,0x54,0x85,0x54}};

STDAPI OnyxCoreIsLoaded();
STDAPI OnyxOsirisIsLoaded();

//---------------------------------------------------------------------------
//						nxUnixCLSID::DllName
//	Look up the Dll (Shareable Image) name for the given ID
//
// November 20, 2017
// -----------------
// The ONYX code is in its final days and we now dont use a global registry
// and this is creating compatibility problems. ONYX is only used for onyxcore
// and osirisonyx so I have hacked solutions for these two objects. OnyxCore now exposes
//	function OnyxCoreIsLoaded which is called, i.e. linked into this file. This saves us
// having too look for the DLL object in the registry
//---------------------------------------------------------------------------

bool nxUnixCLSID::GetDllName( REFCLSID cls, nxString* dllname )
{
	bool						ok = true;

	if (cls == CLSID_OnyxOsirisname   )
	{
	//	ok = (OnyxOsirisIsLoaded() == 0);	// A little UNIX Hack in the final days of Onyx code . Make sure the linker loas in this object and then we dont have to find it on the system path"
		*dllname="onyxosiris.so";
	}
	else if (cls == CLSID_OnyxCoreDatabase )
	{
	//	ok = (OnyxCoreIsLoaded() == 0);
		*dllname="onyx.so";				// A little UNIX Hack in the final days of Onyx code . Make sure the linker loas in this object and then we dont have to find it on the system path"
	}
	else
	{

		nxString					str;
		nxRegistryConfiguration		key("USask-ARG/CLSID","",nxRegistryConfiguration::GLOBAL_INI, false);

		ok =       KeyAsString( cls, &str );
		ok = ok && key.SetFileKeyName( str );
		ok = ok && key.GetString( "InprocServer32", dllname );
	}
	return ok;
}

//---------------------------------------------------------------------------
//						nxUnixCLSID::CreateCLSIDEntry
//	Write the COM CLSID entry into the nxregistry system
//---------------------------------------------------------------------------

bool nxUnixCLSID::CreateCLSIDEntry( REFCLSID cls, nxString dllname )
{
	nxString					str;
	bool						ok;
    nxRegistryConfiguration		key("USask-ARG/CLSID","",nxRegistryConfiguration::GLOBAL_INI, false);


	ok =       KeyAsString( cls, &str );
	ok = ok && key.SetFileKeyName( str );
	ok = ok && key.SetString( "InprocServer32", dllname );
	return ok;
}




