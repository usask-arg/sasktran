#include "sasktranif_internals.h"

#if defined(NX_WINDOWS)
	#include <boost/uuid/uuid.hpp>
	#include <boost/uuid/uuid_generators.hpp>
	#include <boost/uuid/uuid_io.hpp>
#else
	#undef byte
	#include <boost/uuid/uuid.hpp>
	#include <boost/uuid/uuid_generators.hpp>
	#include <boost/uuid/uuid_io.hpp>
#endif




/*-----------------------------------------------------------------------------
 *					SasktranIF_Global_Climatology_Handles		 2015- 11- 18*/
/** **/
/*---------------------------------------------------------------------------*/

class SasktranIF_Global_Climatology_Handles
{
	private:
		std::map< nxString, CLIMATOLOGY_HANDLE>								m_globalhandletable;
		typedef std::map< nxString, CLIMATOLOGY_HANDLE>::value_type			value_type;

		std::map< nxString, CLIMATOLOGY_HANDLE>*                            m_parenthandletable;   // If not null then the handle table of the parent DLL

		bool											InitializeStandardHandles();
	public:
		std::map< nxString, CLIMATOLOGY_HANDLE>*		HandleTable();
		const std::map< nxString, CLIMATOLOGY_HANDLE>*	HandleTable() const;

	public:
														SasktranIF_Global_Climatology_Handles();
		bool											Add_Global_Climatology_Handle ( const char* name, const CLIMATOLOGY_HANDLE&  handle);
		CLIMATOLOGY_HANDLE&								Handle( const char* name, bool printerror = true);
		const char*										NameOfHandle( CLIMATOLOGY_HANDLE& handle);
		bool											HasKey( const char* name) const;
		void											SetParentHandleTable(std::map< nxString, CLIMATOLOGY_HANDLE>* table) { m_parenthandletable = table; }
};

SasktranIF_Global_Climatology_Handles	g_handles;

/*-----------------------------------------------------------------------------
 *					AddGlobalClimatologyHandle		 2015- 11- 19*/
/** **/
/*---------------------------------------------------------------------------*/

bool AddGlobalClimatologyHandle ( const char* name, const CLIMATOLOGY_HANDLE&  handle)
{
	return g_handles.Add_Global_Climatology_Handle( name, handle);
}

std::map<nxString, CLIMATOLOGY_HANDLE>* InternalGlobalClimatologyHandleTable()
{
	return g_handles.HandleTable();
}

std::map< nxString, CLIMATOLOGY_HANDLE>* SasktranIF_Global_Climatology_Handles::HandleTable()
{
	if (m_parenthandletable != nullptr)
	{
		return m_parenthandletable;
	}
	else 
	{
		return &m_globalhandletable;
	}
}

const std::map< nxString, CLIMATOLOGY_HANDLE>* SasktranIF_Global_Climatology_Handles::HandleTable() const
{
	if (m_parenthandletable != nullptr)
	{
		return m_parenthandletable;
	}
	else
	{
		return &m_globalhandletable;
	}
}


/*-----------------------------------------------------------------------------
 *					FindGlobalClimatologyHandle		 2015- 11- 19*/
/** **/
/*---------------------------------------------------------------------------*/

CLIMATOLOGY_HANDLE*  FindGlobalClimatologyHandle ( const char* name, bool printerror)
{
	return &g_handles.Handle( name, printerror);
}


/*---------------------------------------------------------------------------
 *               FindGlobalClimatologyNameOfHandle                2019-09-18 */
/** **/
/*---------------------------------------------------------------------------*/
const char* FindGlobalClimatologyNameOfHandle ( CLIMATOLOGY_HANDLE& handle)
{
	return g_handles.NameOfHandle( handle);
}

bool  AddGeneratedGlobalClimatologyHandleIfNotExists(const char* name)
{
	bool exists = HasKey_InGlobalClimatologyHandle(name);

	if (!exists)
	{
		boost::uuids::uuid uuid = boost::uuids::random_generator()();

		const CLIMATOLOGY_HANDLE* guid = reinterpret_cast<const CLIMATOLOGY_HANDLE*>(&uuid.data[0]);

		AddGlobalClimatologyHandle(name, *guid);
	}

	return true;
}

bool SetParentHandleTable(std::map<nxString, CLIMATOLOGY_HANDLE>* parenttable)
{
	g_handles.SetParentHandleTable(parenttable);
	
	return true;
}


/*-----------------------------------------------------------------------------
 *					HasKey_InGlobalClimatologyHandle		 2016- 7- 8*/
/** **/
/*---------------------------------------------------------------------------*/

bool HasKey_InGlobalClimatologyHandle ( const char* name)
{
	return g_handles.HasKey( name);
}

/*-----------------------------------------------------------------------------
 *					SasktranIF_Global_Climatology_Handles::SasktranIF_Global_Climatology_Handles		 2015- 11- 19*/
/** **/
/*---------------------------------------------------------------------------*/

SasktranIF_Global_Climatology_Handles::SasktranIF_Global_Climatology_Handles()
{
	m_parenthandletable = nullptr;
	InitializeStandardHandles();
}

/*-----------------------------------------------------------------------------
 *					SasktranIF_Global_Climatology_Handles::Add_Global_Climatology_Handle		 2015- 11- 19*/
/** **/
/*---------------------------------------------------------------------------*/

bool SasktranIF_Global_Climatology_Handles::Add_Global_Climatology_Handle ( const char* name, const CLIMATOLOGY_HANDLE&  handle)
{
	bool	ok = true;
	nxString	keyname(name);

	keyname.MakeUpper();

	auto iter = HandleTable()->find( keyname );
	if (iter == HandleTable()->end())
	{
		auto status = HandleTable()->insert( value_type( keyname, handle) );
		ok = status.second;
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SasktranIF CLIMATOLOGY_HANDLE, Error inserting handle for <%s>.", (const char*) name);
		}
	}
	else
	{
		ok = ((*iter).second == handle);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "SasktranIF CLIMATOLOGY_HANDLE, Cannot insert handle for <%s> as this entry already exists with a different value", (const char*) name);
		}
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					SasktranIF_Global_Climatology_Handles::Handle		 2015- 11- 19*/
/** **/
/*---------------------------------------------------------------------------*/

CLIMATOLOGY_HANDLE&	SasktranIF_Global_Climatology_Handles::Handle( const char* name, bool printerror)
{
	bool		ok = true;
	nxString	keyname(name);

	keyname.MakeUpper();
	auto iter = HandleTable()->find( keyname );
	ok = !(iter == HandleTable()->end());
	if (!ok && printerror)
	{
		nxLog::Record(NXLOG_WARNING,"SasktranIF CLIMATOLOGY_HANDLE, cannot find an entry in the global table for <%s>. Returning SKCLIMATOLOGY_UNDEFINED", (const char*)name );
	}
	return ok ? (*iter).second : SKCLIMATOLOGY_UNDEFINED;
}

/*---------------------------------------------------------------------------
 *      SasktranIF_Global_Climatology_Handles::NameOfHandle       2019-09-18 */
/** **/
/*---------------------------------------------------------------------------*/

const char*	SasktranIF_Global_Climatology_Handles::NameOfHandle( CLIMATOLOGY_HANDLE& handle)
{
	bool	ok = true;
	bool	found = false;
	const	char*	name = "UNKNOWN_CLIMATOLOGY";

	auto iter = HandleTable()->begin();
	while (!found && !(iter == HandleTable()->end() ))
	{
		found = (iter->second == handle);
		if (found)
		{
			name = (const char*)(iter->first);
			break;
		}
		++iter;
	}
	if (!found)
	{
		nxLog::Record(NXLOG_WARNING,"SasktranIF_Global_Climatology_Handles::NameOfHandle, could not find the requested climatology handle in the internal global table. Thats not good and may cause knock on issues");
	}
	return name;	
}

/*-----------------------------------------------------------------------------
 *					SasktranIF_Global_Climatology_Handles::Handle		 2015- 11- 19*/
/** **/
/*---------------------------------------------------------------------------*/

bool SasktranIF_Global_Climatology_Handles::HasKey( const char* name) const
{
	bool		ok = true;
	nxString	keyname(name);

	keyname.MakeUpper();
	auto iter = HandleTable()->find( keyname );
	ok = !(iter == HandleTable()->end());
	return ok;
}



/*-----------------------------------------------------------------------------
 *					MAKE_CLIMATOLOGY_HANDLE		 2015- 11- 18*/
/** **/
/*---------------------------------------------------------------------------*/

static CLIMATOLOGY_HANDLE MAKE_CLIMATOLOGY_HANDLE( nxDWORD  DATA1, unsigned short DATA2, unsigned short DATA3,  unsigned char  B0, unsigned char  B1, unsigned char  B2, unsigned char  B3,unsigned char  B4, unsigned char  B5, unsigned char  B6, unsigned char  B7)
{
	CLIMATOLOGY_HANDLE x = {  DATA1, DATA2, DATA3, { B0, B1, B2, B3, B4, B5, B6, B7 }  };
	return x;

}

#define ADD_STANDARD_ENTRY(NAME, DATA1, DATA2, DATA3, B0, B1, B2, B3, B4, B5, B6, B7) NAME=MAKE_CLIMATOLOGY_HANDLE( DATA1, DATA2, DATA3, B0, B1, B2, B3, B4, B5, B6, B7); Add_Global_Climatology_Handle( #NAME, NAME);	


/*-----------------------------------------------------------------------------
 *					SasktranIF_Global_Climatology_Handles::InitializeStandardHandles		 2015- 11- 18*/
/** **/
/*---------------------------------------------------------------------------*/

bool SasktranIF_Global_Climatology_Handles::InitializeStandardHandles()
{
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_AOA_DAYS,								0x9655be40, 0x69d7, 0x4d49,  0x9c, 0x11, 0xd8, 0x2e, 0x67, 0x62, 0xc9, 0xe9 );	// Age of Air in days
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_AEROSOLDUST_CM3,						0xd6b6cc7f, 0x5254, 0x4c00,  0x95, 0x76, 0x8c, 0xfb, 0x00, 0x81, 0xf6, 0x8e );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_AEROSOLH2SO4_CM3,						0x83c71d08, 0x52e6, 0x479b,  0xb4, 0x3e, 0x6a, 0x53, 0xca, 0x2c, 0x67, 0x91 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_AEROSOLICE_CM3,						0xd77f970e, 0xee72, 0x4915,  0x80, 0xe0, 0x8c, 0xd0, 0xb2, 0x63, 0x32, 0x67 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_AEROSOLSURFACEAREA_UM2PerCM3,			0x59ab91f1, 0x120e, 0x4bb7,  0x9c, 0x78, 0xb2, 0x9c, 0x94, 0x9e, 0x2c, 0xef );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_AEROSOLWATER_CM3,						0xd4142f12, 0xd975, 0x4f9f,  0xa4, 0x30, 0xd5, 0xfb, 0xb7, 0xf4, 0xd0, 0xd7 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_AEROSOL_CM3,							0x853c5f2b, 0x45a9, 0x42e7,  0x94, 0x9e, 0x17, 0x97, 0xff, 0x84, 0x33, 0xa1 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_AEROSOL_EXTINCTIONPERKM,				0x7cbec40a, 0x533e, 0x4cf3,  0x92, 0x11, 0x82, 0xad, 0x13, 0xba, 0x38, 0x0e );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3,					0xe7b90bbd, 0x1212, 0x4820,  0x81, 0x62, 0x69, 0x1a, 0x6a, 0x26, 0x90, 0xa4 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_ALBEDO,								0xb4498463, 0xe2e0, 0x4b2d,  0x96, 0x9f, 0x62, 0x49, 0x63, 0xdc, 0xc8, 0xef );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_Ar_CM3,								0xa5d739a1, 0x2985, 0x4645,  0x90, 0x40, 0xc0, 0xdc, 0x1d, 0xea, 0x49, 0x87 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_BCPI,									0x2580e2ba, 0x6fbe, 0x4768,  0x88, 0xe5, 0x71, 0x43, 0xae, 0xf1, 0x12, 0x4f );	// Black Carbon 
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_BCPO,									0xa9432fb9, 0x5fe3, 0x4d5a,  0x95, 0x6e, 0xef, 0x0d, 0x32, 0xcc, 0x81, 0x64 );	// Black Carbon 
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_BRCL_CM3,								0x0ac086dc, 0x56d2, 0x41b5,  0xb7, 0xb4, 0xb0, 0x4f, 0xe0, 0x3f, 0x57, 0x24 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_BRNO3_CM3,								0x002931f2, 0x35c7, 0x4a9f,  0x9a, 0x5f, 0x93, 0x9c, 0x73, 0x0e, 0x75, 0x90 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_BRO_CM3,								0xbd6480ce, 0xc70d, 0x489c,  0xab, 0x4c, 0x88, 0xbd, 0x61, 0xa3, 0x2f, 0x3a );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_BRO_VMR,								0xe773778a, 0xf999, 0x4ac6,  0x95, 0x70, 0x53, 0x24, 0xce, 0x27, 0x4f, 0xfb );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_BRX_CM3,								0xbdb5ed8e, 0x0214, 0x4a26,  0x95, 0xf4, 0xd0, 0x2e, 0x67, 0xf3, 0x1c, 0xdd );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_BRX_VMR,								0x68049a8c, 0x8d88, 0x438e,  0xba, 0xda, 0x0f, 0xf9, 0xc5, 0x87, 0x3d, 0xd9 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_BRY_CM3,								0xa9c7e1dc, 0x20f8, 0x4cdc,  0x87, 0x9c, 0xa2, 0x05, 0x9f, 0x73, 0xc6, 0xe2 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_BRY_VMR,								0x1ce710f4, 0xa5d7, 0x4d87,  0x97, 0xca, 0x06, 0x13, 0xac, 0x7d, 0x1c, 0x15 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_BR_CM3,								0x0df46c61, 0x44f8, 0x4d19,  0x9c, 0x5b, 0x75, 0xda, 0xdd, 0x5c, 0xdf, 0xb2 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_C2H2_CM3,								0x37bfe48a, 0xfe5f, 0x4d9b,  0x8a, 0x4a, 0xb2, 0xbb, 0xc5, 0xaa, 0xf2, 0x40 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_C2H4_CM3,								0x89b882a7, 0x494b, 0x47fe,  0xb2, 0x74, 0x6e, 0x11, 0x72, 0x17, 0x03, 0x61 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_C2H6_CM3,								0xd7360803, 0x0d69, 0x4a71,  0xa9, 0x73, 0x60, 0x49, 0x62, 0xae, 0xa5, 0xc6 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_C3H6O_CM3,								0xc402ca2b, 0x0fb1, 0x4277,  0xb4, 0x37, 0xf7, 0x74, 0x5b, 0xbb, 0x1b, 0xd4 );	// Acetone
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_C3H6O_VMR,								0xfb8d75b6, 0x1e09, 0x41ce,  0xa5, 0x10, 0x07, 0x98, 0x5d, 0x2d, 0xda, 0xb6 );	// Acetone
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_C5H8_CM3,								0xde12a207, 0xe9f3, 0x443a,  0xbc, 0x04, 0x02, 0x70, 0xe4, 0x48, 0x94, 0x5f );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_C5H8_VMR,								0x2f223c94, 0x74bd, 0x4c71,  0x93, 0xc5, 0x86, 0xd3, 0xdd, 0xdf, 0x32, 0x06 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CCL4_CM3,								0x46e70f90, 0x90dc, 0x487f,  0x8a, 0xb2, 0x21, 0x4c, 0x0e, 0x7f, 0x9a, 0xd4 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CCL4_VMR,								0x7aaee908, 0x5cd6, 0x4d49,  0x80, 0x3f, 0xaa, 0x90, 0x0d, 0x8f, 0x3c, 0x1c );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CF2CL2_CM3,							0x7ad2a63e, 0x1f39, 0x45fb,  0x8f, 0x07, 0x51, 0xf6, 0xab, 0xf4, 0x81, 0x56 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CF2CL2_VMR,							0x26f8bc13, 0xe36f, 0x48d5,  0xb0, 0x92, 0xf1, 0x80, 0xa5, 0x41, 0xe6, 0xb0 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CF4_CM3,								0xfe8fdfd7, 0x442e, 0x4680,  0xa9, 0xc6, 0x47, 0x91, 0xa2, 0x1d, 0x08, 0x7f );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CFCL3_CM3,								0x685c28e4, 0xe2d1, 0x4246,  0xba, 0x1f, 0xbb, 0x4a, 0x6f, 0xbb, 0x82, 0x03 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CFCL3_VMR,								0xc1012ef1, 0x5180, 0x4dcc,  0x86, 0x76, 0xd1, 0x6a, 0x00, 0xd2, 0x3d, 0xa5 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH2O_CM3,								0xff01a5f7, 0x21c8, 0x4a1b,  0xaa, 0x33, 0x0c, 0x14, 0x2a, 0xd9, 0x17, 0x74 );	// Formaldehyde
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH2O_VMR,								0x67f082ae, 0x18da, 0x46e4,  0x90, 0x90, 0x2d, 0x16, 0x84, 0x3a, 0x26, 0x0f );	// Formaldehyde
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH3BR_CM3,								0x7efb84e5, 0x7ddd, 0x4702,  0xb2, 0x76, 0x5e, 0x09, 0xe5, 0x50, 0xba, 0xad );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH3BR_VMR,								0x625f90cd, 0x64c2, 0x424d,  0xb3, 0x1e, 0xc5, 0xb8, 0x57, 0x08, 0x96, 0x12 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH3CL_CM3,								0xee929203, 0x38e6, 0x4cdd,  0x99, 0x56, 0x0d, 0x86, 0x33, 0x47, 0xb4, 0xf7 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH3CL_VMR,								0x4432ce62, 0x122d, 0x46dd,  0xb0, 0xf6, 0x7e, 0xe6, 0x0c, 0xe3, 0xd2, 0x81 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH3CN_CM3,								0xb53983c5, 0xfd51, 0x4eb6,  0x98, 0x09, 0x5b, 0xab, 0x3f, 0x66, 0x62, 0x2d );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH3I_CM3,								0x04338ce3, 0x6cf6, 0x4c0f,  0x93, 0x7b, 0xa4, 0x93, 0x71, 0x4b, 0x38, 0xbd );	// Methyl Iodide
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH3I_VMR,								0xff146136, 0xc1d5, 0x40d5,  0xb7, 0x58, 0x7b, 0x05, 0x33, 0xa8, 0x26, 0x03 );	// Methyl Iodide
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH3OH_CM3,								0x64999dc3, 0x3ea6, 0x4fd0,  0xa4, 0x00, 0x92, 0x5f, 0xdb, 0xd7, 0xb9, 0x49 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH4_CM3,								0x5e971864, 0xc45c, 0x4228,  0xad, 0x51, 0x0f, 0xac, 0xbc, 0xbf, 0x6b, 0xe3 );	// Methane
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH4_VMR,								0x2f9719b3, 0xc2c3, 0x4493,  0xbf, 0x70, 0x9c, 0x1e, 0x18, 0x55, 0xc6, 0x06 );	// Methane
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CL2O2_CM3,								0x30078ee6, 0xf1d2, 0x45b5,  0xaf, 0x83, 0x1f, 0xa3, 0x97, 0x4b, 0x4d, 0x2c );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CL2_CM3,								0x245c3e5c, 0x4c50, 0x4ff2,  0x82, 0xf1, 0xd2, 0x55, 0xdb, 0x86, 0x94, 0xab );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CLNO3_CM3,								0x9d688eb3, 0x52b5, 0x4eb0,  0x83, 0x9f, 0x23, 0xd2, 0xba, 0xba, 0x3b, 0x6e );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CLONO2_CM3,							0x1f0121b6, 0x1e2b, 0x435b,  0x92, 0x37, 0xb6, 0xd1, 0x54, 0xae, 0x2c, 0xa9 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CLOUD_FRACTION,						0xb241a613, 0xf103, 0x4f6c,  0x92, 0xa3, 0x5b, 0x0a, 0xfb, 0x0f, 0xbe, 0xe7 );	// Inorganic Nitrates
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CLO_CM3,								0x2d182b56, 0x0e44, 0x4c11,  0xb8, 0xc1, 0x26, 0x44, 0x54, 0x79, 0xe7, 0xbc );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CLY_CM3,								0x109a535f, 0x0c1d, 0x4b72,  0xa9, 0xc2, 0xb2, 0x51, 0xa9, 0x33, 0x9a, 0x9e );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CLY_VMR,								0xcad6adcf, 0x59ba, 0x4450,  0x97, 0xda, 0x1d, 0x9a, 0xd3, 0x53, 0x30, 0xaf );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CL_CM3,								0x9d98fc4f, 0x8453, 0x4853,  0x98, 0xf6, 0x32, 0x69, 0x5f, 0x16, 0x07, 0xbd );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CO2_CM3,								0xd72e084e, 0xcf7a, 0x4003,  0x9a, 0xd4, 0x38, 0x05, 0x3b, 0xe5, 0x41, 0x24 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CO2_VMR,								0x77226852, 0x7a55, 0x47d5,  0x86, 0xb8, 0xee, 0x34, 0xb3, 0xdb, 0x33, 0xae );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_COF2_CM3,								0x2819d7f9, 0x5a43, 0x4285,  0xb6, 0x92, 0xde, 0xd0, 0x44, 0x51, 0x61, 0xbe );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CO_CM3,								0xedd67f72, 0x4dac, 0x4ab2,  0xad, 0x78, 0x30, 0x42, 0x5f, 0x1b, 0xa8, 0xe1 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CO_VMR,								0xced33bca, 0x751a, 0x4883,  0x9d, 0x81, 0x85, 0x55, 0x15, 0xab, 0xf8, 0xa1 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_DUST_0p7mu,							0x36aba684, 0xc223, 0x4365,  0xbd, 0x6c, 0x96, 0xb9, 0x85, 0x4d, 0x6a, 0xb4 );	// Dust Aerosol, Reff 0.7 micron
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_DUST_1p4mu,							0x1501a55e, 0x02f4, 0x444f,  0xba, 0xcc, 0x60, 0xb9, 0xc6, 0x6a, 0x2e, 0x31 );	// Dust Aerosol, Reff 1.4 micron
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_DUST_2p4mu,							0x560560e6, 0xe20d, 0x4f28,  0x82, 0x17, 0xd4, 0x38, 0x29, 0xb8, 0x53, 0x7d );	// Dust Aerosol, Reff 1.4 micron
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_DUST_4p5mu,							0x07b3d22e, 0x7a28, 0x4633,  0xa4, 0xb6, 0x9a, 0x53, 0xf0, 0x84, 0x69, 0xe9 );	// Dust Aerosol, Reff 4.5 micron
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_EFFECTIVESIZE_MICRONS,					0x71062c28, 0x6e13, 0x4296,  0xa5, 0x26, 0x27, 0xfa, 0xc9, 0x6a, 0xe9, 0xcd );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_EPV,									0xd31ffad4, 0x427c, 0x4008,  0x9c, 0x03, 0x36, 0x01, 0x23, 0xb9, 0x06, 0x17 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_GEOMETRIC_HEIGHT,						0x14556024, 0x8d1b, 0x4f69,  0x9b, 0x8f, 0xfb, 0x0e, 0x19, 0x41, 0x6b, 0x9f );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_GEOPOTENTIAL_HEIGHT,					0xed808183, 0x20fd, 0x4214,  0x88, 0x49, 0xd7, 0xfc, 0x38, 0x6c, 0x3e, 0x72 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_H2CO_CM3,								0xd8aa034a, 0xaa47, 0x4771,  0x97, 0xcc, 0x5f, 0x59, 0xb8, 0x27, 0xd0, 0x67 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_H2O2_CM3,								0x5488ab42, 0x77c2, 0x4ad0,  0x88, 0x76, 0x3b, 0x37, 0x7f, 0x0d, 0x72, 0x17 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_H2O_CM3,								0xc91894c1, 0x0328, 0x47f3,  0x9a, 0x38, 0xa3, 0x77, 0x7c, 0x33, 0x6e, 0x02 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_H2O_VMR,								0x20b8db2e, 0x11f0, 0x436d,  0x8c, 0x94, 0x20, 0x79, 0xd3, 0xec, 0x71, 0xd7 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_H2S_CM3,								0x82a06da3, 0x4554, 0x4b23,  0xa6, 0x73, 0x83, 0x2a, 0x70, 0xbe, 0x54, 0xd5 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_H2_CM3,								0xdcfa09da, 0x131b, 0x49be,  0x89, 0x67, 0x81, 0xc5, 0xd9, 0x4c, 0xe1, 0xdb );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_H2_VMR,								0xd9266719, 0x40b1, 0x4e58,  0x8a, 0xb5, 0x52, 0xef, 0xb5, 0x8c, 0x2d, 0xaf );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HBR_CM3,								0xdf43ec49, 0x30db, 0x46e5,  0x8e, 0xdd, 0x2c, 0xa8, 0x07, 0x91, 0xca, 0x88 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HCL_CM3,								0x171bb589, 0x6314, 0x460f,  0xad, 0xad, 0x1a, 0x4a, 0xcd, 0x54, 0x61, 0x94 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HCN_CM3,								0x389e83ff, 0xde73, 0x41fe,  0xa2, 0x31, 0xa0, 0x13, 0x69, 0x73, 0xdb, 0xd0 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HCOOH_CM3,								0x58e4dfa4, 0x7731, 0x4b75,  0x91, 0x89, 0x3c, 0xad, 0xdf, 0x9d, 0x8b, 0xe7 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HF_CM3,								0xc78a42cc, 0x7aa2, 0x442b,  0xaf, 0xd8, 0x8a, 0xa9, 0xf6, 0xeb, 0x2b, 0x16 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HI_CM3,								0xb8891417, 0xdb8b, 0x460a,  0xae, 0x27, 0x90, 0x98, 0xb1, 0xb9, 0x29, 0xef );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HNO2_CM3,								0xac45ebea, 0x5ff3, 0x46e8,  0xa2, 0x09, 0x3b, 0x09, 0xc6, 0xc5, 0x08, 0x2b );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HNO2_VMR,								0xf79e74b6, 0xe4c1, 0x4269,  0x87, 0x45, 0x3c, 0x17, 0x9c, 0xb0, 0x3f, 0x33 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HNO3_CM3,								0x5341ba6d, 0xf6f7, 0x479d,  0x8c, 0x22, 0x3a, 0x32, 0xd8, 0x69, 0x3a, 0xb2 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HNO3_VMR,								0xb5cd58fb, 0x0bad, 0x4f57,  0xaa, 0x45, 0x70, 0xb0, 0x0b, 0x9f, 0x32, 0x58 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HNO4_CM3,								0x3a827321, 0xf89e, 0x47a8,  0x8c, 0x86, 0xa5, 0x66, 0x72, 0x49, 0x85, 0x3c );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HO2_CM3,								0x48fd60ea, 0x5eff, 0x44d9,  0x8e, 0xe3, 0x84, 0x77, 0xd3, 0x2c, 0xc9, 0x22 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HOBR_CM3,								0x430ee9ef, 0x45db, 0x4b57,  0x98, 0x60, 0x98, 0x46, 0x08, 0x24, 0x0e, 0x31 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_HOCL_CM3,								0x48271b4b, 0x9337, 0x454e,  0x9c, 0x61, 0x06, 0x09, 0x3f, 0x7d, 0x5b, 0x80 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_H_CM3,									0x31251479, 0x3d77, 0x49d9,  0x8d, 0x05, 0x58, 0xc2, 0x99, 0x64, 0x95, 0x98 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_He_CM3,								0x8e8f6d2f, 0x1158, 0x4bf0,  0x89, 0xb7, 0x26, 0x07, 0x2c, 0xfa, 0x90, 0xa7 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_ICE_CM3,								0x9526c08d, 0x5c4f, 0x4f49,  0xb4, 0x33, 0x89, 0x70, 0x50, 0x99, 0x34, 0x6d );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_JH2O,									0xd1bfb436, 0xc2c5, 0x460d,  0xb5, 0xf3, 0xed, 0x9d, 0x28, 0x19, 0xde, 0x05 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS,			0xdfa33ed1, 0x3e0b, 0x4b84,  0xb3, 0x88, 0xaa, 0xc3, 0xf9, 0x60, 0xb8, 0x70 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH,					0xe59f7ef0, 0xacdd, 0x473a,  0xaf, 0x55, 0x01, 0x76, 0x5b, 0x16, 0x5f, 0x2a );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_MECL_CM3,								0xa206cd04, 0x8f8e, 0x429f,  0x94, 0x88, 0x8d, 0x0d, 0x33, 0xdf, 0xb2, 0xb3 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_MECL_VMR,								0x6e9c6ebf, 0x400c, 0x47c3,  0x90, 0xbd, 0x88, 0xf4, 0x65, 0x09, 0xd1, 0x76 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_N2O5_CM3,								0xbcc2cd2d, 0xe357, 0x4e35,  0x84, 0x5d, 0x7a, 0x1e, 0x4c, 0x6b, 0x53, 0xa6 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_N2O_CM3,								0x80eff6a0, 0xaaac, 0x4141,  0x94, 0xf6, 0x2d, 0xfe, 0xd4, 0xbe, 0x7a, 0x62 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_N2O_VMR,								0x967938a2, 0x4dbe, 0x4fea,  0xb4, 0x10, 0x23, 0xb5, 0xe8, 0xbf, 0x1f, 0xd3 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_N2_CM3,								0x971af1d3, 0xce72, 0x4b0d,  0x8d, 0x41, 0xe6, 0x3b, 0x72, 0x75, 0x11, 0xaf );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_N2_VMR,								0x4d18ad9f, 0xf679, 0x47e5,  0x92, 0x80, 0x42, 0xbf, 0xa8, 0x96, 0x27, 0x61 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_NH3_CM3,								0xfe4bb454, 0x4236, 0x434e,  0x8c, 0x34, 0x3e, 0x63, 0x5c, 0x21, 0xdf, 0x8f );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_NH3_VMR,								0xd8bd3a31, 0x8b1c, 0x49f5,  0xb7, 0x12, 0x81, 0x7d, 0x84, 0x95, 0x06, 0x0d );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_NITS,									0x4805a300, 0xa650, 0x4655,  0x99, 0x3c, 0xd0, 0x45, 0x79, 0x05, 0x67, 0xe8 );	// Inorganic Nitrates
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_NO2_CM3,								0xdc60d0a2, 0x8a95, 0x4803,  0x99, 0xb4, 0x4b, 0x6d, 0xd3, 0xa4, 0xaf, 0xda );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_NO2_VMR,								0xb7dc3cbd, 0xd16e, 0x4b0f,  0xac, 0x4b, 0x85, 0x1c, 0xe8, 0xc8, 0x7a, 0x33 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_NO3_CM3,								0x804751f5, 0x508d, 0x4a95,  0xa3, 0x1b, 0x94, 0xe3, 0xbd, 0xb5, 0x6c, 0x41 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_NOPLUS_CM3,							0x2beaadc0, 0x44ed, 0x4a39,  0x88, 0x66, 0xa2, 0x00, 0x6c, 0x72, 0xc8, 0xd0 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_NOY_CM3,								0x55cffc62, 0x12d7, 0x48e5,  0x9c, 0xec, 0x60, 0xe7, 0xbc, 0xe3, 0x9a, 0x5d );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_NOY_VMR,								0xb31c77a7, 0x7b02, 0x464b,  0x88, 0x8f, 0xf2, 0x0e, 0x37, 0x20, 0xea, 0xb7 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_NO_CM3,								0xa6d247ef, 0xdd46, 0x47b6,  0x84, 0xe9, 0x0d, 0x26, 0x4b, 0x8a, 0x0d, 0x75 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_NO_VMR,								0xa81a1db0, 0x0c27, 0x4151,  0x98, 0x35, 0x9d, 0xaa, 0x64, 0xe0, 0x94, 0xca );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_N_CM3,                                 0x363f6e4f, 0xc08e, 0x480e,  0x96, 0x9f, 0xa6, 0xe1, 0x10, 0x1f, 0xa6, 0x93 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_O2_CM3,								0x05e3fb16, 0xe2e6, 0x4feb,  0xa7, 0x30, 0xa6, 0x74, 0xaf, 0x32, 0xad, 0x34 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_O2_O2_CM6,                             0xa7e3d865, 0x4ea8, 0x4ef0,  0x8f, 0xa9, 0x51, 0xfb, 0xf2, 0x54, 0x08, 0x7c );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_O2_VMR,								0x373f9b2c, 0xae05, 0x4fe1,  0x8e, 0x26, 0xe6, 0x71, 0x90, 0xf8, 0x56, 0x09 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_O3_CM3,								0xdfe865bc, 0xe8ff, 0x4603,  0xa5, 0xba, 0xf3, 0x9e, 0xef, 0x1f, 0x08, 0x2f );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_O3_VMR,								0x820b6ee7, 0x0f8b, 0x4e6f,  0x85, 0x29, 0x27, 0x18, 0xfd, 0x4e, 0xf2, 0x5d );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_OCLO_CM3,								0x9e8fd34c, 0x3c14, 0x42ad,  0x9c, 0x77, 0xd3, 0x27, 0x19, 0xec, 0x88, 0x28 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_OCPI,									0xb1163887, 0x6815, 0x4435,  0x8e, 0x14, 0x94, 0x90, 0x57, 0xcd, 0x2d, 0x38 );	// Organic Carbon Aerosol I
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_OCPO,									0xfb357aa7, 0x6662, 0x45e3,  0xb9, 0x23, 0xe4, 0x31, 0x59, 0xbc, 0xa7, 0xbf );	// Organic Carbon Aerosol O
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_OCS_CM3,								0xcc443e47, 0x7a87, 0x4a20,  0x9f, 0xd6, 0x40, 0xe1, 0x56, 0xb9, 0x3c, 0x01 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_OH_CM3,								0x1f3d66ee, 0x330f, 0x4d8c,  0x8d, 0x56, 0xbb, 0x55, 0x11, 0x66, 0xcf, 0x45 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_O_CM3,									0x2f907c2d, 0x72e1, 0x4e9e,  0x8e, 0xad, 0x02, 0x16, 0x1a, 0x48, 0x33, 0x26 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_PAN_CM3,								0x724b435e, 0x0bc7, 0x43e9,  0x9d, 0x12, 0xc8, 0xaa, 0x98, 0x41, 0x9a, 0x99 );	// Peroxy acetyl nitrate
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_PAN_VMR,								0x82aac576, 0xfd63, 0x4416,  0xa6, 0xfb, 0x98, 0xcd, 0xa2, 0x4a, 0x81, 0x3a );	// Peroxy acetyl nitrate
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_PH3_CM3,								0xe035b73b, 0xc4d2, 0x4ce8,  0x8c, 0xcc, 0xc8, 0x37, 0x80, 0x29, 0x19, 0xa4 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_POTENTIAL_TEMPERATURE_K,				0x2e9919d6, 0x38f4, 0x4700,  0x9e, 0x42, 0xdb, 0x54, 0x23, 0x39, 0xbf, 0xc8 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_PRESSURE_PA,							0x9d908f8e, 0xdb2d, 0x4717,  0xbf, 0x81, 0xc2, 0x3d, 0x54, 0x2e, 0xaa, 0x4d );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_QI_MMR,								0x2b6b3a64, 0x2294, 0x4330,  0x96, 0x64, 0xee, 0x4e, 0x62, 0x1c, 0xa0, 0x19 );	// mass fraction of cloud ice water
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_QL_MMR,								0xf44ea1c0, 0x0f2c, 0x4bff,  0x8d, 0x60, 0x67, 0x14, 0x7d, 0xef, 0xbe, 0x90 );	// mass fraction of cloud liquid water
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_QV,									0x73380fa8, 0x66c8, 0x4ffa,  0x9f, 0x40, 0x6e, 0xdb, 0x93, 0x9a, 0xe5, 0x09 );	// specific humidity
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_RH,									0x7272ba90, 0xa96a, 0x4d7c,  0xb3, 0x2b, 0xf0, 0xd1, 0xca, 0xc5, 0x2c, 0x1a );	// Relative Humidity
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_ROOH_CM3,								0xc48fb55f, 0xa2a6, 0x4caf,  0xa7, 0xd8, 0xd4, 0x92, 0x4a, 0xe6, 0x15, 0xba );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_ROO_CM3,								0xab78f29d, 0xa2be, 0x4a5c,  0x98, 0x30, 0x92, 0x7c, 0xc2, 0x60, 0xac, 0xf9 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_SALA,									0xbb72b76f, 0xf7ba, 0x4499,  0xa1, 0xd8, 0xb9, 0x66, 0xd2, 0x95, 0x99, 0x1c );	// Sea Salt aerosol Accum
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_SALC,									0x4a2a1403, 0x1b19, 0x4d65,  0xab, 0x2d, 0x74, 0x33, 0x9d, 0x54, 0x24, 0x5a );	// Sea Salt aerosol Coarse
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_SF6_CM3,								0x47d8cc73, 0xeedf, 0x48c2,  0x86, 0xbd, 0xa2, 0xda, 0xe7, 0x10, 0xf1, 0xe6 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_SO2_CM3,								0x1c6a5a6d, 0xbb24, 0x4574,  0xa0, 0x98, 0xa5, 0x96, 0x07, 0x30, 0x41, 0xfe );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_SO2_VMR,								0x7a6eb67e, 0xeaa7, 0x4ebf,  0xa5, 0xe5, 0x7e, 0xb4, 0xd6, 0x2e, 0xe7, 0x4f );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_SO4_CM3,								0x567270fa, 0xcaf4, 0x464f,  0xaa, 0xef, 0x45, 0x38, 0xeb, 0xf7, 0x7a, 0x8c );	// Sea Salt aerosol Coarse
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_SO4_VMR,								0xb552baa7, 0xc92b, 0x47a0,  0x8b, 0x27, 0x80, 0x2c, 0xc0, 0x24, 0x77, 0x9a );	// Sea Salt aerosol Coarse
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_SURFACE_GEOMETRIC_HEIGHT,				0x36920f4c, 0x3519, 0x4917,  0xa1, 0xad, 0x4e, 0xc3, 0x20, 0xaa, 0xbb, 0x2d );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_SURFACE_GEOPOTENTIAL_HEIGHT,			0xf3066173, 0x79c4, 0x4420,  0x92, 0xac, 0x4e, 0xf0, 0xe6, 0x86, 0xa5, 0xb3 );	// Surface Geopotential Height
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_SURFACE_PRESSURE_PA,					0x1cb6731d, 0x5d14, 0x4a06,  0x84, 0xb5, 0x58, 0x35, 0xf3, 0xa9, 0xea, 0xbe );	// Surface Pressure
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_TEMPERATURE_K,							0x3f4333e6, 0x290f, 0x436a,  0xbb, 0x07, 0x86, 0x37, 0x87, 0x04, 0x22, 0x1e );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_UNDEFINED,								0xe42a521f, 0xde5d, 0x4f2c,  0x96, 0xe9, 0x9d, 0x20, 0x41, 0xf7, 0xa0, 0x79 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_XXX_CM3,								0x531bff83, 0x61d9, 0x4f02,  0x8d, 0x5d, 0xc2, 0x0c, 0xbe, 0x29, 0x95, 0xce );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_XXX_VMR,								0x243b8f9c, 0x1589, 0x4783,  0x9d, 0x2c, 0x29, 0x0d, 0x7f, 0x08, 0xa5, 0x0e );

//	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH3BR_CM3,								0x1a567543, 0xfdd4, 0x4aff,  0xa8, 0xe3, 0xdf, 0xed, 0x06, 0x5c, 0x48, 0x0c );
//	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH3CL_CM3,								0x86d50c4f, 0x7415, 0x47ca,  0x9f, 0x50, 0x01, 0xa3, 0x8f, 0xa0, 0x89, 0x2a );
//	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CH4_CM3,								0x734cb84c, 0xb0c5, 0x47e8,  0x94, 0xb3, 0xa1, 0x6a, 0x96, 0x5e, 0x5c, 0x82 );
//	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CO2_VMR,								0x563c613e, 0x5a1d, 0x43ce,  0x9d, 0x13, 0x16, 0x81, 0x53, 0x56, 0x24, 0x83 );
//	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_CO_CM3,								0xa7bb5257, 0x7d28, 0x4128,  0xbf, 0xf7, 0x64, 0x13, 0xb9, 0x33, 0xc8, 0x93 );
//	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_H2O_CM3,								0x9eaa10ee, 0xbdcf, 0x4f01,  0x94, 0xca, 0x2b, 0x34, 0x8c, 0xca, 0xd4, 0xaa );
//	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_NH3_CM3,								0x8f67bc89, 0xf9c0, 0x48cb,  0xa4, 0x6c, 0x3f, 0x69, 0x7d, 0x1f, 0x7d, 0x19 );


	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_0,							0xa223f0cb, 0xbf59, 0x437a,  0x94, 0xd7, 0x98, 0xe1, 0xb0, 0xb3, 0x8d, 0xef );
	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_1,							0x5ba8a2fa, 0x241c, 0x4ffc,  0xaf, 0x3e, 0xc6, 0x27, 0x7a, 0x0b, 0x76, 0x95 );
	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_2,							0xf1dc3e11, 0x5dff, 0x4408,  0xb4, 0x55, 0x34, 0x17, 0xca, 0x95, 0x93, 0x12 );
	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_3,							0xd4ce809d, 0xa5df, 0x4130,  0xab, 0xb0, 0xab, 0x68, 0xb5, 0x9e, 0x67, 0x49 );
	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_4,							0x33e13d04, 0xf76d, 0x44bf,  0x85, 0xc5, 0xe0, 0x1b, 0x5b, 0xb2, 0xc1, 0xf5 );
	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_5,							0x2ec50d19, 0x9a0e, 0x4d44,  0x8e, 0x1d, 0x9a, 0x7d, 0x41, 0xe3, 0xc5, 0x72 );
	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_6,							0xc039d79c, 0x634a, 0x4302,  0xa5, 0x03, 0xb6, 0xea, 0x86, 0xe9, 0x06, 0x21 );
	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_7,							0x69c34827, 0xe23b, 0x454f,  0xac, 0x7b, 0xb4, 0xd8, 0xcd, 0x5f, 0xf0, 0x48 );
	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_8,							0xb38d611e, 0xe334, 0x4260,  0x9f, 0xe3, 0xb5, 0x13, 0xfc, 0x45, 0x92, 0x3d );
	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_9,							0x29c440a8, 0x7a46, 0x4f2e,  0xa8, 0x1c, 0xd4, 0xff, 0x32, 0xf6, 0xfc, 0x8b );
	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_O2,							0xf50fa01a, 0xbf6b, 0x40fe,  0xb5, 0x27, 0x69, 0xdc, 0xba, 0xe0, 0xb0, 0x59 );
	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_OH,							0xeff9bffa, 0xeb56, 0x4757,  0xaf, 0x12, 0xf2, 0xe3, 0x46, 0x0f, 0x00, 0xba );
	ADD_STANDARD_ENTRY(SKEMISSION_PHOTOCHEMICAL_O3,							0x7f0895b7, 0x0d10, 0x4d70,  0x90, 0x96, 0xad, 0xcb, 0x35, 0x05, 0x89, 0x46 );
	ADD_STANDARD_ENTRY(SKEMISSION_THERMAL,									0xbeefd857, 0x5c8d, 0x47ee,  0xb9, 0xc4, 0x8f, 0x5f, 0x45, 0x7b, 0xaf, 0xce );

// -- These values are deprecated. Dont use them in new code as they will be phased out.
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_GAMMA_EFFECTIVERADIUS_MICRONS,			0x6d7ca761, 0xc752, 0x4df0,  0x9c, 0x1c, 0xc5, 0xa6, 0xf9, 0xe9, 0xf4, 0xbd );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_GAMMA_EFFECTIVEVARIANCE_PERMICRON,		0xfbd93d8d, 0xf5f7, 0x4c9e,  0x87, 0xd1, 0x8a, 0xda, 0x14, 0xa4, 0x16, 0x6c );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_AEROSOL_CM3_MODE2,						0xb1fec2b5, 0x9d27, 0x4dfe,  0xb5, 0x45, 0xa4, 0x28, 0xc9, 0x69, 0x09, 0x51 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_AEROSOL_MODERADIUS_MICRONS,			0xd6ff72a0, 0x3491, 0x4193,  0x8d, 0x0c, 0xa9, 0x0c, 0xb6, 0xed, 0x24, 0xa6 );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_AEROSOL_MODEWIDTH,						0x4540b8e6, 0x2056, 0x4ffa,  0x83, 0x43, 0xbb, 0x65, 0xc3, 0xc1, 0x3d, 0x2d );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_ICE_MODERADIUS_MICRONS,				0x75aeb6b8, 0xf18d, 0x4ad6,  0xae, 0x80, 0x57, 0xa0, 0x14, 0xf0, 0xc3, 0xfc );
	ADD_STANDARD_ENTRY(SKCLIMATOLOGY_ICE_MODEWIDTH,							0x60ebd4ea, 0xb7ea, 0x4ec4,  0xbc, 0x18, 0xd2, 0xe9, 0xdf, 0xa4, 0x89, 0xdc );

   return true;
}






/*-----------------------------------------------------------------------------
 *					Global Constants								2005-6-24
 *	These constants define the GUID's for all of the parameters supported
 *	by the climatatology classes. These are still exported by Python
 *---------------------------------------------------------------------------*/


 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_AOA_DAYS;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_AEROSOLDUST_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_AEROSOLH2SO4_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_AEROSOLICE_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_AEROSOLSURFACEAREA_UM2PerCM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_AEROSOLWATER_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_AEROSOL_CM3; 
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_AEROSOL_EXTINCTIONPERKM;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_ALBEDO;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_Ar_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_BCPI;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_BCPO;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_BRCL_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_BRNO3_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_BRO_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_BRO_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_BRX_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_BRX_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_BRY_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_BRY_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_BR_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_C2H2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_C2H4_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_C2H6_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_C3H6O_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_C3H6O_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_C5H8_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_C5H8_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CCL4_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CCL4_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CF2CL2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CF2CL2_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CF4_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CFCL3_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CFCL3_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CH2O_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CH2O_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CH3BR_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CH3BR_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CH3CL_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CH3CL_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CH3CN_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CH3I_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CH3I_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CH3OH_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CH4_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CH4_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CL2O2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CL2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CLNO3_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CLONO2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CLOUD_FRACTION;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CLO_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CLY_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CLY_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CL_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CO2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CO2_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_COF2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CO_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_CO_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_DUST_0p7mu;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_DUST_1p4mu;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_DUST_2p4mu;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_DUST_4p5mu;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_EFFECTIVESIZE_MICRONS;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_EPV;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_GEOMETRIC_HEIGHT;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_GEOPOTENTIAL_HEIGHT;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_H2CO_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_H2O2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_H2O_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_H2O_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_H2S_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_H2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_H2_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HBR_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HCL_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HCN_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HCOOH_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HF_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HI_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HNO2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HNO2_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HNO3_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HNO3_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HNO4_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HO2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HOBR_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_HOCL_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_H_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_He_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_ICE_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_JH2O;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_LOGNORMAL_MODERADIUS_MICRONS; 
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_LOGNORMAL_MODEWIDTH;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_MECL_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_MECL_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_N2O5_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_N2O_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_N2O_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_N2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_N2_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_NH3_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_NH3_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_NITS;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_NO2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_NO2_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_NO3_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_NOPLUS_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_NOY_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_NOY_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_NO_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_NO_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_N_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_O2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_O2_O2_CM6;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_O2_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_O3_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_O3_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_OCLO_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_OCPI;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_OCPO;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_OCS_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_OH_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_O_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_PAN_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_PAN_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_PH3_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_POTENTIAL_TEMPERATURE_K;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_PRESSURE_PA;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_QI_MMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_QL_MMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_QV;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_RH;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_ROOH_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_ROO_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_SALA;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_SALC;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_SF6_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_SO2_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_SO2_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_SO4_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_SO4_VMR;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_SURFACE_GEOMETRIC_HEIGHT;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_SURFACE_GEOPOTENTIAL_HEIGHT;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_SURFACE_PRESSURE_PA;

 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_TEMPERATURE_K;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_UNDEFINED;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_XXX_CM3;
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_XXX_VMR;
 
 // -- definitions added to support pratmo
// CLIMATOLOGY_HANDLE		 SKCLIMATOLOGY_CH3BR_CM3;
// CLIMATOLOGY_HANDLE		 SKCLIMATOLOGY_CH3CL_CM3;
// CLIMATOLOGY_HANDLE		 SKCLIMATOLOGY_CH4_CM3;
// CLIMATOLOGY_HANDLE		 SKCLIMATOLOGY_CO2_VMR;
// CLIMATOLOGY_HANDLE		 SKCLIMATOLOGY_CO_CM3;
// CLIMATOLOGY_HANDLE		 SKCLIMATOLOGY_H2O_CM3;
// CLIMATOLOGY_HANDLE		 SKCLIMATOLOGY_NH3_CM3;


 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_0;
 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_1;
 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_2;
 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_3;
 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_4;
 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_5;
 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_6;
 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_7;
 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_8;
 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_9;
 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_O2;
 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_OH;
 EMISSION_HANDLE			 SKEMISSION_PHOTOCHEMICAL_O3;
 EMISSION_HANDLE			 SKEMISSION_THERMAL;


// -- These values are deprecated. Dont use them in new code as they will be phased out.
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_GAMMA_EFFECTIVERADIUS_MICRONS;		// { static CLIMATOLOGY_HANDLE x = {  0x6d7ca761, 0xc752, 0x4df0, { 0x9c, 0x1c, 0xc5, 0xa6, 0xf9, 0xe9, 0xf4, 0xbd }  }; return x; }
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_GAMMA_EFFECTIVEVARIANCE_PERMICRON;	// { static CLIMATOLOGY_HANDLE x = {  0xfbd93d8d, 0xf5f7, 0x4c9e, { 0x87, 0xd1, 0x8a, 0xda, 0x14, 0xa4, 0x16, 0x6c }  }; return x; }
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_AEROSOL_CM3_MODE2;					// { static CLIMATOLOGY_HANDLE x = {  0xb1fec2b5, 0x9d27, 0x4dfe, { 0xb5, 0x45, 0xa4, 0x28, 0xc9, 0x69, 0x09, 0x51 }  }; return x; }
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_AEROSOL_MODERADIUS_MICRONS;			// { static CLIMATOLOGY_HANDLE x = {  0xd6ff72a0, 0x3491, 0x4193, { 0x8d, 0x0c, 0xa9, 0x0c, 0xb6, 0xed, 0x24, 0xa6 }  }; return x; }
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_AEROSOL_MODEWIDTH;					// { static CLIMATOLOGY_HANDLE x = {  0x4540b8e6, 0x2056, 0x4ffa, { 0x83, 0x43, 0xbb, 0x65, 0xc3, 0xc1, 0x3d, 0x2d }  }; return x; }
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_ICE_MODERADIUS_MICRONS;				// { static CLIMATOLOGY_HANDLE x = {  0x75aeb6b8, 0xf18d, 0x4ad6, { 0xae, 0x80, 0x57, 0xa0, 0x14, 0xf0, 0xc3, 0xfc }  }; return x; }
 CLIMATOLOGY_HANDLE		SKCLIMATOLOGY_ICE_MODEWIDTH;						// { static CLIMATOLOGY_HANDLE x = {  0x60ebd4ea, 0xb7ea, 0x4ec4, { 0xbc, 0x18, 0xd2, 0xe9, 0xdf, 0xa4, 0x89, 0xdc }  }; return x; }


/*-----------------------------------------------------------------------------
 *					operator < for Climatology handles	2014-1-30*/
/** The less than operator for the CLIMATOLOGY HANDLES **/
/*---------------------------------------------------------------------------*/

bool operator < ( const CLIMATOLOGY_HANDLE& a,  const CLIMATOLOGY_HANDLE& b)
{
	bool less;
	size_t	idx;

	less = a.Data1 < b.Data1;								// see if the first word is less
	if (!less) if ( a.Data1 == b.Data1)						// check to see if its equel
	{														// If first byte is equal
		less = a.Data2 < b.Data2;							// then check 2nd word
		if (!less) if ( a.Data2 == b.Data2)					// if not less then see if equal
		{													// if 1st and 2 words equal
			less = a.Data3 < b.Data3;						// then compare the 3rd word
			if (!less) if ( a.Data3 == b.Data3)				// if 1st, 2nd and 3rd word equal
			{												// then 
				idx = 0;
				while (!less && (idx < 8))
				{
					less = a.Data4[idx] < b.Data4[idx];			// check the 4th word
					if (!less)
					{
						if (a.Data4[idx]  > b.Data4[idx]) break;
						idx++;
					}
				}
			}
		}
	}
	return less;
}




