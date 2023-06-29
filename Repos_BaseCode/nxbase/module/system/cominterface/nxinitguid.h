/*-------------------------------------------------------------
				UNIX version of MICROSOFT INITGUID.H
	Initguid should only be invoked in one file of a project
	otherwise each project may have significant storage space
	used up by the defintion of GUIDs.
---------------------------------------------------------------*/

#if !defined(NX_INITGUID)
  #error You must define NX_INITGUID at the top of the file and just include nxlib.h. Dont explicitly include nxinitguid.h
#endif

#if !defined(NOMINMAX)
#if (_MSC_VER < 1200)
  #define NOMINMAX
#endif
#endif

#if defined( NX_WINDOWS) && !defined(NX_UNIX_VER)
 #include <objbase.h>
 #include <initguid.h>
#else
 #define INITGUID
 #ifdef DEFINE_GUID
	#undef  DEFINE_GUID
 #endif
#endif
