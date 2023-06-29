/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

/* The REGSAM values are taken from the Windows Registry API. These values should only be used in systems that do not have access to the
 * to the Windows API.
*/
/*
enum REGSAM {	KEY_ALL_ACCESS,				// Combination of KEY_QUERY_VALUE, KEY_ENUMERATE_SUB_KEYS, KEY_NOTIFY, KEY_CREATE_SUB_KEY, KEY_CREATE_LINK, and KEY_SET_VALUE access.
				KEY_CREATE_LINK ,			// Permission to create a symbolic link.
				KEY_CREATE_SUB_KEY,			// Permission to create subkeys.
				KEY_ENUMERATE_SUB_KEYS,		// Permission to enumerate subkeys.
				KEY_EXECUTE,				// Permission for read access.
				KEY_NOTIFY,					// Permission for change notification.
				KEY_QUERY_VALUE,			// Permission to query subkey data.
				KEY_READ,					// Combination of KEY_QUERY_VALUE, KEY_ENUMERATE_SUB_KEYS, and KEY_NOTIFY access.
				KEY_SET_VALUE, 				// Permission to set subkey data.
				KEY_WRITE
			};
*/

/*
typedef const char*	HKEY;
extern HKEY	HKEY_LOCAL_MACHINE;
extern HKEY HKEY_CURRENT_USER;



#if defined( NX_UNIX_VER )
*/

typedef HRESULT (*DLLCANUNLOADNOW)  ( void );
typedef HRESULT (*DLLGETCLASSOBJECT)( REFCLSID rclsid, REFIID riid, void **ppv );
typedef HRESULT (*DLLSETREGISTRY)( const char* dirname );
//--------------------------------------------------------------------------
//						nxComDllEntry
//--------------------------------------------------------------------------

class  nxComDllEntry
{
	private:
		DLLCANUNLOADNOW			m_unload;
		DLLGETCLASSOBJECT		m_getclassobject;
		void*					m_handle;
		nxString				m_dllname;
		CLSID					m_class;


	public:
								nxComDllEntry();
								~nxComDllEntry();
		void					Close();
		nxBOOL					Load( REFCLSID cls, const char* dllname);
		HRESULT					GetClassObject( REFCLSID rclsid, REFIID riid, void **ppv );
		HRESULT					UnloadIfUnused();
		REFCLSID				Class() const {return m_class;}
		const nxString&			DllName() const { return m_dllname;}
		bool					operator== (const nxComDllEntry&       ) const { return nxFALSE; }
		bool					operator!= (const nxComDllEntry&       ) const { return nxFALSE; }
		bool					operator<  (const nxComDllEntry&       ) const { return nxFALSE;}
		bool                    operator>  (const nxComDllEntry&       ) const { return nxFALSE;}
};
//#endif

//---------------------------------------------------------------------------
//						class nxUnixCLSID
//	A class that looks up CLSID in the registry.  Used by UNIX apps to emulate
//	the WIN32 System Regsitry when looking up CLSID relationships.
//--------------------------------------------------------------------------

class nxUnixCLSID
{
	private:
		bool						KeyAsString( REFCLSID cls, nxString* str );

	public:
									nxUnixCLSID		();
								   ~nxUnixCLSID		();
		bool						GetDllName      ( REFCLSID cls, nxString* dllname );
		bool 						CreateCLSIDEntry( REFCLSID cls, nxString dllname  );

//		nxComDllEntry*				Find(REFCLSID cls);
//		void						UnloadUnusedDll();
};
