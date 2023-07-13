/*---------------------------------------------------------------------------
			UNIXOLE2.H
This file replaces the functions and definitions in Microsoft
windows <ole2.h>

---------------------------------------------------------------------------*/
#if !defined( UNIXOLE2_H)
#define UNIXOLE2_H	1
#include "module/system/unix/rpcndr.h"

//typedef long						HRESULT;
#define LPCSTR						const char *
#define FAR
#define EXTERN_C					extern "C"
#define LONG						long
#define DWORD						nxDWORD
#define WORD						nxWORD
#define BYTE						nxBYTE
typedef nxDWORD						ULONG;
typedef nxBOOL						BOOL;
typedef nxBOOL						VARIANT_BOOL;
typedef double						DATE;

#if !defined(TRUE)
#define TRUE	1
#endif

#if !defined(FALSE)
#define FALSE   0
#endif


#define interface					struct
#define STDMETHOD(method)			virtual HRESULT  method
#define STDMETHODIMP				HRESULT
#define STDMETHOD_(type,method)		virtual type  method
#define STDMETHODIMP_(type)			type
#define STDAPI						EXTERN_C HRESULT
#define STDAPI_(type)				EXTERN_C type
#define WINOLEAPI					STDAPI
#define WINOLEAPI_(type)			STDAPI_(type)
#define PURE						= 0
#define THIS_
#define THIS									void
#define DECLARE_INTERFACE(iface)				interface iface
#define DECLARE_INTERFACE_(iface, baseiface)    interface iface : public baseiface
#define LPVOID									void*

#define E_UNEXPECTED                     0x8000FFFFL
#define E_NOINTERFACE                    0x80000004L
#define E_OUTOFMEMORY                    0x80000002L
#define CLASS_E_NOAGGREGATION            0x80040110L
#define CLASS_E_CLASSNOTAVAILABLE        0x80040111L
#define E_FAIL							 0x80004005L
struct GUID
{
    DWORD	Data1;
    WORD	Data2;
    WORD	Data3;
    BYTE	Data4[8];
};

#define IID					GUID
#define	CLSID				GUID
#define REFGUID             const GUID &
#define REFIID              const IID &
#define REFCLSID            const CLSID &
//#define REFFMTID            const FMTID &

inline BOOL operator== (REFGUID a, REFGUID b) { return  (      (a.Data1 == b.Data1)
															&& (a.Data2 == b.Data2)
															&& (a.Data3 == b.Data3)
															&& (a.Data4[0] == b.Data4[0])
															&& (a.Data4[1] == b.Data4[1])
															&& (a.Data4[2] == b.Data4[2])
															&& (a.Data4[3] == b.Data4[3])
															&& (a.Data4[4] == b.Data4[4])
															&& (a.Data4[5] == b.Data4[5])
															&& (a.Data4[6] == b.Data4[6])
															&& (a.Data4[7] == b.Data4[7])
													     );}

typedef unsigned short VARTYPE;

enum VARENUM
  {	VT_EMPTY	= 0,
	VT_NULL	= 1,
	VT_I2	= 2,
	VT_I4	= 3,
	VT_R4	= 4,
	VT_R8	= 5,
	VT_CY	= 6,
	VT_DATE	= 7,
	VT_BSTR	= 8,
	VT_DISPATCH	= 9,
	VT_ERROR	= 10,
	VT_BOOL	= 11,
	VT_VARIANT	= 12,
	VT_UNKNOWN	= 13,
	VT_DECIMAL	= 14,
	VT_I1	= 16,
	VT_UI1	= 17,
	VT_UI2	= 18,
	VT_UI4	= 19,
	VT_I8	= 20,
	VT_UI8	= 21,
	VT_INT	= 22,
	VT_UINT	= 23,
	VT_VOID	= 24,
	VT_HRESULT	= 25,
	VT_PTR	= 26,
	VT_SAFEARRAY	= 27,
	VT_CARRAY	= 28,
	VT_USERDEFINED	= 29,
	VT_LPSTR	= 30,
	VT_LPWSTR	= 31,
	VT_FILETIME	= 64,
	VT_BLOB	= 65,
	VT_STREAM	= 66,
	VT_STORAGE	= 67,
	VT_STREAMED_OBJECT	= 68,
	VT_STORED_OBJECT	= 69,
	VT_BLOB_OBJECT	= 70,
	VT_CF	= 71,
	VT_CLSID	= 72,
	VT_VECTOR	= 0x1000,
	VT_ARRAY	= 0x2000,
	VT_BYREF	= 0x4000,
	VT_RESERVED	= 0x8000,
	VT_ILLEGAL	= 0xffff,
	VT_ILLEGALMASKED	= 0xfff,
	VT_TYPEMASK	= 0xfff
   };


typedef struct tagVARIANT  {
    VARTYPE vt;
    unsigned short wReserved1;
    unsigned short wReserved2;
    unsigned short wReserved3;
    union {
        unsigned char        bVal;                        // VT_UI1.
        short                    iVal;                        // VT_I2    .
        long                    lVal;                        // VT_I4    .
        float                    fltVal;                    // VT_R4    .
        double                dblVal;                    // VT_R8    .
        VARIANT_BOOL        boolVal;                        // VT_BOOL.
		DATE				date;							// VT_DATE
        unsigned char        FAR* pbVal;                // VT_BYREF|VT_UI1.
        short                    FAR* piVal;                // VT_BYREF|VT_I2.
        long                    FAR* plVal;                // VT_BYREF|VT_I4.
        float                    FAR* pfltVal;            // VT_BYREF|VT_R4.
        double                FAR* pdblVal;            // VT_BYREF|VT_R8.
        VARIANT_BOOL        FAR* pboolVal;                // VT_BYREF|VT_BOOL.
    };
} VARIANT;

enum CLSCTX
    {	CLSCTX_INPROC_SERVER	= 0x1,
	CLSCTX_INPROC_HANDLER	= 0x2,
	CLSCTX_LOCAL_SERVER	= 0x4,
	CLSCTX_INPROC_SERVER16	= 0x8,
	CLSCTX_REMOTE_SERVER	= 0x10,
	CLSCTX_INPROC_HANDLER16	= 0x20,
	CLSCTX_INPROC_SERVERX86	= 0x40,
	CLSCTX_INPROC_HANDLERX86	= 0x80
    };

#define SUCCEEDED(Status) ((HRESULT)(Status) >= 0)
#define S_OK                                   ((HRESULT)0x00000000L)
#define S_FALSE                                ((HRESULT)0x00000001L)

interface IUnknown
{
	virtual 	~IUnknown(){}
	STDMETHOD (QueryInterface)(THIS_ REFIID riid, LPVOID FAR* ppvObj) PURE;
    STDMETHOD_(ULONG, AddRef)(THIS) PURE;
    STDMETHOD_(ULONG, Release)(THIS) PURE;
};



interface IClassFactory : public IUnknown
{
	virtual 	~IClassFactory(){}

	STDMETHOD	(CreateInstance)( THIS_ IUnknown * pUnkOuter, REFIID riid, void ** ppvObject ) PURE;
    STDMETHOD   (LockServer)    ( THIS_ BOOL fLock) PURE;
};

interface IDispatch : public IUnknown
{
		virtual 	~IDispatch(){}

	//	define the Dispatch but dont try to make it do anything on Unix systems
};

typedef IUnknown* LPUNKNOWN;


/*--------------------------------------------------------------
					DEFINE_GUID
---------------------------------------------------------------*/

#ifndef NX_INITGUID
	#define NXDEFINE_GUID(name, l, w1, w2, b1, b2, b3, b4, b5, b6, b7, b8) EXTERN_C const GUID name
#else
	#define NXDEFINE_GUID(name, l, w1, w2, b1, b2, b3, b4, b5, b6, b7, b8) EXTERN_C const GUID name = { l, w1, w2, { b1, b2,  b3,  b4,  b5,  b6,  b7,  b8 } }

#endif

EXTERN_C const GUID IID_IUnknown; //  0x00000000L, 0x0000, 0x0000, 0xc0,0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);
EXTERN_C const GUID IID_IClassFactory; //, 0x00000001L, 0x0000, 0x0000, 0xc0,0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x46);


extern "C" HRESULT CoInitialize( LPVOID );
extern "C" void    CoUninitialize();
extern "C" HRESULT CoGetClassObject( REFCLSID rclsid, DWORD clsctx,	LPVOID serverinfo, REFIID riid, LPVOID * ppv);
extern "C" HRESULT CoCreateInstance( REFCLSID rclsid, LPVOID aggregate, DWORD clsctx,  REFIID riid,	LPVOID * ppv);
extern "C" void    CoFreeUnusedLibraries();

extern "C" void VariantInit( VARIANT*  pvarg  );
extern "C" HRESULT VariantClear( VARIANT *  pvarg );
#endif
