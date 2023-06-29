#if !defined(NXREGISTRY_INCLUDE_H)
#define NXREGISTRY_INCLUDE_H 1

/** \addtogroup nxwin32registry **/
/** @{ */

//#if !defined(NXREGISTRYCONFIG_DEFAULTCOMPANY)

/** Defines the default company used when accessing registry configuration
 *	with class #nxRegistryConfiguration.  You can change this value by defining
 *	it before including nxRegistryKeyWin32.h  (via nxlib.h for example)
 **/

//#define NXREGISTRYCONFIG_DEFAULTCOMPANY "USask-ARG"
//#endif

class nxRegistryKey;


/*-----------------------------------------------------------------------------
 *					nxRegistryKey_RegistryLocation		 2016- 11- 1*/
/** **/
/*---------------------------------------------------------------------------*/

class nxRegistryKey_RegistryLocation
{
	private:
		nxString			m_basedirectory;
		bool				m_usenativeregistry;

	private:
		bool				Default_BaseDirectory( );

	public:
							nxRegistryKey_RegistryLocation( );
							nxRegistryKey_RegistryLocation( const char* dirname ) {Set_BaseDirectory(dirname);}
		bool				Set_BaseDirectory( const char* dirname );
		const nxString&		BaseDirectory    () const {return  m_basedirectory;}
		bool				UseNativeRegistry() const { return m_usenativeregistry;}

};


/*-----------------------------------------------------------------------------
 *					class nxRegistryConfiguration		2005-3-30*/
/** A class for accessing application specific info in the registry at a
 *	higher level than class #nxRegistryKeyWin32. This class assumes you are storing
 *	low volume configuration data in the registry specific to a given application.
 *	The class uses values stored in a subkey based upon "company" name and current
 *	application.  By default the company name is HKCU/Software/USask-ARG and the
 *	application name is derived from the Module currently executing under the Windows API.
 *
 *	\par Saving and Reading Values
 *	The current code allows for \b nxString, \b double, \b int and \b nxBOOL to be stored in the registry
 *	as named values.  Please note that \b int are only stored to 32 bit precision and are not yet compatible
 *	with 64 bit int.
 *
 *	\par Reading Directory values
 *	Since we anticipate that the registry will be used to save the location of directories
 *	that may be used for caching and other purposes we have provided code, #LocateDirectoryFromKey, to automatically browse
 *	and/or create directories when required.
 *
**/
/*---------------------------------------------------------------------------*/

class nxRegistryConfiguration							//! Save/Read application configuration from registry.
{

	public:
		enum				INIRootLocation { GLOBAL_INI, USER_INI, USER_APPL_INI, GLOBAL_APPL_INI };
		enum				INIAccessRights { READ_INI,   WRITE_INI, FULLIO_INI};

	private:
		bool				m_usenativeregistry;
		nxString			m_filekey;					//!< The name of the file that will store keys. This will not change during a program execution.
		nxString			m_companyname;
		INIRootLocation		m_rootlocation;				//!< What part of the registry is accessed, USER_INI		generated HKCU/Software/USask-ARG/
														//                                          USER_APPL_INI   generates HKCU/Software/USask-ARG/ApplicationSettings/
														//											GLOBAL_INI      generates HKLM/Software/USask-ARG/
														//											GLOBAL_APPL_INI generates HKLM/Software/USask-ARG/ApplicationSettings/
	private:
		nxString&			GetApplicationName					();
		bool				CloseKey							(nxRegistryKey*	 key);
		bool				OpenKey								(nxRegistryKey** key, INIAccessRights accessmode);
		nxBOOL				BrowseForDirectory					( nxString* dirname, const char* promptstr );
		nxString			GenerateKeyName						();
		bool				RemoveLeadingAndTrailingDirChars	( nxString& basekey );
		nxString			BaseKeyName							();

	public:
//							nxRegistryConfiguration	();
							nxRegistryConfiguration	(const char* companyname, const char* filekey, nxRegistryConfiguration::INIRootLocation rootlocation, bool usenativeregistry);
		virtual			   ~nxRegistryConfiguration	();
 static bool				FlushRegistry			( INIRootLocation location);
		nxBOOL				SetRootLocation			( nxRegistryConfiguration::INIRootLocation rootlocation);
		bool				SetFileKeyName			( const char* filekey);
		nxBOOL				SetCompanyKeyName		( const char* keyname);
		nxBOOL				SetString				( const char * subkeyname, const char *  str );
		nxBOOL				GetString				( const char * subkeyname, nxString* str );
		nxBOOL				SetPath					( const char * subkeyname, const char *  str );
		nxBOOL				GetPath					( const char * subkeyname, nxString* str );
		nxBOOL				SetDouble				( const char * subkeyname, double	value);
		nxBOOL				GetDouble				( const char * subkeyname, double* value, double defaultvalue);
		nxBOOL				SetInteger				( const char * subkeyname, int	   value);
		nxBOOL				GetInteger				( const char * subkeyname, int*   value, int defaultvalue);
		nxBOOL				SetBool					( const char * subkeyname, nxBOOL value);
		nxBOOL				GetBool					( const char * subkeyname, nxBOOL * value, nxBOOL defaultvalue);
		nxBOOL				LocateDirectoryFromKey	( const char * subkeyname, nxString* dirname, nxBOOL createifnotexist, nxBOOL browseifnokey, const char* browsepromptstring);
		bool				ShallowCopy				( const nxRegistryConfiguration& other);


};

/*-----------------------------------------------------------------------------
 *					nxRegistryKey		2009-6-2*/
/** A base class for the windows or YAML Registry Key functions**/
/*---------------------------------------------------------------------------*/

class nxRegistryKey
{
	private:
		static nxRegistryKey_RegistryLocation		g_registrylocation;
	
	public:
		virtual				   ~nxRegistryKey		(){}
		virtual bool			DestroyKeyHierarchy	()                                         = 0;	//!< Destroys the key object stored on the heap.
		virtual bool			SetString		( const char * valuename, const char *	str  ) = 0;	//!< Set a named value to this string
		virtual bool			GetString		( const char * valuename, nxString*		str  ) = 0;	//!< Get a string from the named value
		virtual bool			SetDouble		( const char * valuename, double		value) = 0;	//!< Set a named value to this double
		virtual bool			GetDouble		( const char * valuename, double*		value) = 0;	//!< Get a double from this named value
		virtual bool			SetInteger		( const char * valuename, int			value) = 0;	//!< Set a named value to this integer
		virtual bool			GetInteger		( const char * valuename, int*			value) = 0;	//!< Get an integer from the named value
		virtual bool			SetBool			( const char * valuename, bool			value) = 0;	//!< Set a named value to this boolean
		virtual bool			GetBool			( const char * valuename, bool*			value) = 0;	//!< Get a boolean from this named value

	public:
		static const nxRegistryKey_RegistryLocation& 	RegistryLocation		( ) { return  g_registrylocation; }
		static       nxRegistryKey_RegistryLocation* 	RegistryLocationVar		( ) { return &g_registrylocation; }



};

/*-----------------------------------------------------------------------------
 *					class nxRegistryKeyWin32								2004-11-23*/
/** Provides basic manipulation of the Windows Registry. Use the static member
 *	nxRegistryKeyWin32::CreateKey to create keys: this is the only way to create
 *	keys as the constructor is a "private" member.
 **/
/*---------------------------------------------------------------------------*/

#if defined(NX_WINDOWS)
class nxRegistryKeyWin32 : public nxRegistryKey								//! Base registry key utilities.
{
	private:
		nxRegistryKeyWin32*				m_parent;
		HKEY							m_key;
		nxString						m_keyname;

	private:
										nxRegistryKeyWin32	( nxRegistryKeyWin32* parent, const char * name, HKEY basekey, REGSAM access );
		nxBOOL							SetValue			( const char* valuename,  DWORD dwType, const BYTE* lpData, DWORD cbData);
		nxBOOL							GetValue			( const char* valuename,  DWORD dwType,  BYTE* lpData, DWORD* cbData, DWORD mandatorysize);

		static HKEY						RootLocation		( nxRegistryConfiguration::INIRootLocation rootlocation);
		static REGSAM					AccessRights		( nxRegistryConfiguration::INIAccessRights accessrights);

	public:
		virtual						   ~nxRegistryKeyWin32	();
		HKEY							Key					() const {return m_key;}
		static nxRegistryKeyWin32*		CreateKey			(  const char * basekeyname, const char* filekey, nxRegistryConfiguration::INIRootLocation rootlocation, nxRegistryConfiguration::INIAccessRights access);

	public:
		virtual bool					DestroyKeyHierarchy		() { delete this; return true;}
		virtual bool					SetString		( const char * valuename, const char *	str  );
		virtual bool					GetString		( const char * valuename, nxString*		str  );
		virtual bool					SetDouble		( const char * valuename, double		value);
		virtual bool					GetDouble		( const char * valuename, double*		value);
		virtual bool					SetInteger		( const char * valuename, int			value);
		virtual bool					GetInteger		( const char * valuename, int*			value);
		virtual bool					SetBool			( const char * valuename, bool			value);
		virtual bool					GetBool			( const char * valuename, bool*			value);
};
#endif

/*-----------------------------------------------------------------------------
 *					nxRegistryValueUnix		2009-6-1*/
/** **/
/*---------------------------------------------------------------------------*/

class  nxRegistryValueUnix
{
	private:
		nxString	m_value;
		nxString	m_name;

	public:
							nxRegistryValueUnix	()											{};
		void				SetName				( const nxString& name)						{ m_name = name; m_name.MakeLower(); }
		void				SetValue			( const nxString& value)					{ m_value = value;}
		const nxString&		Value				() const									{ return m_value;}
		const nxString&		Name				() const									{ return m_name;}
		bool				operator==			( const nxRegistryValueUnix& other) const	{ return (other.m_name==m_name) && (other.m_value == m_value);}
		bool				operator!=			( const nxRegistryValueUnix& other) const	{ return !(operator==(other));}
		bool				operator <			( const nxRegistryValueUnix& other) const	{ return other.m_name < m_name;}
		bool				operator >			( const nxRegistryValueUnix& other) const	{ return other.m_name > m_name;}
};


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyUnix		2009-6-2*/
/** A class that implements a simple registry like functionality using text
 *	based file. This class implements the functionality of a single registry
 *	key using an ".ini" fil.
 *
 *	\par ini File Format
 *	Section Names are case insensitive.
*   Values  Names are case insensitive
*   Values are case sensitive
*   Exclamation marks (even inside strings!!!) mark the start of a comment.
*
*    [1stSectionName]
*	           = avalue				! Optional Default Value
*      Key1    = value2
*	   Key2    = value1
*	   [subkey1]
*        Key1  = value2
*	     Key2  = value1
*	   []
*	 []
*
*	 [2ndSectionName]
*      Key1  = value2
*	   Key2  = value1
*	 []
*
**/
/*---------------------------------------------------------------------------*/

class nxRegistryKeyUnix : public nxRegistryKey
{
	private:
		nxRegistryConfiguration::INIAccessRights					m_accessmode;
		nxString													m_fullpathname;
		nxString													m_sectionname;
		bool														m_isdirty;
		nxRegistryKeyUnix*											m_parent;
		std::list<nxRegistryKeyUnix*>								m_subkeys;
		std::list<nxRegistryValueUnix>								m_values;

	private:
		typedef std::list<nxRegistryKeyUnix*>::iterator				subkeyIterator;
		typedef std::list<nxRegistryKeyUnix*>::const_iterator		const_subkeyIterator;
		typedef std::list<nxRegistryValueUnix>::iterator			valueIterator;
		typedef std::list<nxRegistryValueUnix>::const_iterator		const_valueIterator;

		enum KEYSTATE  { NEWSUBKEY, NEWVALUE, ENDOFKEY,  NOP };


	private:
		class EqualKeys
		{
			private: nxString	m_kname;
			public:
							EqualKeys			(const char *name )					{ m_kname = name;}
				bool		operator()			( const nxRegistryKeyUnix* other)	{ return (other->SectionName() == m_kname);}
		};

 		class Equalvalues
		{
			private: nxString	m_vname;

			public:
							Equalvalues		( const char *name )				{ m_vname = name;}
				bool		operator()		( const nxRegistryValueUnix& other) { return (other.Name() == m_vname);}
		};


	private:
		KEYSTATE							ParseLine			( nxString &line, nxString* name, nxString* value );
		void								erase				();
		void								ReadKey				( nxFile& infile );
		void								ReadFile			( const char * fileName );
		bool								FindKey				( nxString& keyname,   nxRegistryKeyUnix** key ) ;
		bool								FindValue			( nxString& valuename, nxRegistryValueUnix** value ) ;
		nxRegistryKeyUnix*					RootParent			( );
		bool								CheckAndSaveFile	();				// Writes the contents back to disk, required in FULL ACCESS MODE, only callable from root key
		bool								WriteKey			(nxFile& f);	// Writes the contents of this key to disk in text format.
		void								SetDirty			()				{m_isdirty = true;}
		bool								AddNewValue			( const char* name, const char* str );
		bool								WriteEntriesToFile	( nxFile& f );
		bool								WriteFile			();


	public:
											nxRegistryKeyUnix	( nxRegistryKeyUnix* m_parent, nxRegistryConfiguration::INIAccessRights accessmode);
		virtual							   ~nxRegistryKeyUnix	();
		void								SetSectionName		( const nxString& sectionname );
		bool								GetValue			( const char *name, nxString* value ) const;
		const nxString&						SectionName			() const					{ return m_sectionname;}
		bool								FlushWriteIO		() { return CheckAndSaveFile();}
		static nxRegistryKeyUnix*			CreateKey			(  const char * basekeyname, const char* filekey, nxRegistryConfiguration::INIRootLocation rootlocation, nxRegistryConfiguration::INIAccessRights access);

		public:
			virtual bool					DestroyKeyHierarchy	();
			virtual bool					SetString		( const char * valuename, const char *	str  );
			virtual bool					GetString		( const char * valuename, nxString*		str  );
			virtual bool					SetDouble		( const char * valuename, double		value);
			virtual bool					GetDouble		( const char * valuename, double*		value);
			virtual bool					SetInteger		( const char * valuename, int			value);
			virtual bool					GetInteger		( const char * valuename, int*			value);
			virtual bool					SetBool			( const char * valuename, bool			value);
			virtual bool					GetBool			( const char * valuename, bool*			value);


};


/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml		2009-6-2*/
/** A class that implements a simple registry like functionality using a YAML
 *	based file. 
**/
/*---------------------------------------------------------------------------*/

#include "yaml-cpp/yaml.h"

/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYamlRoot		 2016- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

class nxRegistryKeyYamlRoot
{
	private:
		YAML::Node				m_root;
		const std::string		m_filename;
		std::string				m_fullfilename;
		bool					m_isdirty;

	public:
								nxRegistryKeyYamlRoot			( const char* filename);
							   ~nxRegistryKeyYamlRoot			();
		YAML::Node*				Root							()		{ return &m_root;}
		void					SetDirty						()		{ m_isdirty = true;}
		bool					SetYAMLOutputFolder				( const char* dirname );
		bool					CheckYamlLoaded					();
		bool					CheckDirtyAndSaveIfRequired		();

		static nxRegistryKeyYamlRoot*	GlobalApplKey();
		static nxRegistryKeyYamlRoot*	GlobalKey();
		static nxRegistryKeyYamlRoot*	UserKey();
		static nxRegistryKeyYamlRoot*	UserApplKey();

};



/*-----------------------------------------------------------------------------
 *					nxRegistryKeyYaml		 2016- 11- 7*/
/** **/
/*---------------------------------------------------------------------------*/

class nxRegistryKeyYaml : public nxRegistryKey
{
	private:
		nxRegistryConfiguration::INIAccessRights				m_accessmode;
		nxRegistryKeyYamlRoot*									m_basenode;
		YAML::Node												m_node;

	private:
		nxRegistryKeyYamlRoot*				RootParent				() { return m_basenode;}
		bool								SetStringInternal		( const char * valuename, const char* str  );

	public:
											nxRegistryKeyYaml	( YAML::Node& node, nxRegistryKeyYamlRoot* basenode, nxRegistryConfiguration::INIAccessRights accessmode);
		virtual							   ~nxRegistryKeyYaml	();
		static nxRegistryKeyYamlRoot*		FindRootRegistry	(nxRegistryConfiguration::INIRootLocation		rootlocation);
		static nxRegistryKeyYaml*			CreateKey			(  const char * basekeyname, const char* filekey, nxRegistryConfiguration::INIRootLocation rootlocation, nxRegistryConfiguration::INIAccessRights access);
		static bool							FlushRegistry		( nxRegistryConfiguration::INIRootLocation rootlocation );

		public:
			virtual bool					DestroyKeyHierarchy	() override;
			virtual bool					SetString		( const char * valuename, const char *	str  ) override;
			virtual bool					GetString		( const char * valuename, nxString*		str  ) override;
			virtual bool					SetDouble		( const char * valuename, double		value) override;
			virtual bool					GetDouble		( const char * valuename, double*		value) override;
			virtual bool					SetInteger		( const char * valuename, int			value) override;
			virtual bool					GetInteger		( const char * valuename, int*			value) override;
			virtual bool					SetBool			( const char * valuename, bool			value) override;
			virtual bool					GetBool			( const char * valuename, bool*			value) override;

};


#endif
