/*------------------------------------------------------------------
	Dummy:- rpcndr.h

	13-May-1999	NDL

	This file is included to replicate the types defined by the
	Microsoft MIDL compiler.  It purpose is to just make the code
	compile properly.

	Its intended only for use on UNIX systems and was kluged together by
	Nick Lloyd.

	The file may occasional maintenance if the MIDL compiler is upgraded.

*/



#if !defined(__RPC_FAR)

#if !defined(__cplusplus)
#define __cplusplus
#endif
#define __RPC_FAR
#define __RPC_USER
#define RPC_IF_HANDLE int
#define MIDL_INTERFACE(x) struct
#define STDMETHODCALLTYPE
#define HRESULT long
#define UINT    long
#define interface struct
#define __RPC_STUB
#define IRpcStubBuffer void
#define IRpcChannelBuffer void
#define PRPC_MESSAGE  int
#define boolean int
#define BSTR char*
/*
enum VARENUM {
    VT_EMPTY = 0,
    VT_NULL = 1,
    VT_I2 = 2,
    VT_I4 = 3,
    VT_R4 = 4,
    VT_R8 = 5,
    VT_CY = 6,
    VT_DATE = 7,
    VT_BSTR = 8,
    VT_DISPATCH = 9,
    VT_ERROR = 10,
    VT_BOOL = 11,
    VT_VARIANT = 12,
    VT_UNKNOWN = 13,
    VT_DECIMAL = 14,
    VT_I1 = 16,
    VT_UI1 = 17,
    VT_UI2 = 18,
    VT_UI4 = 19,
    VT_I8 = 20,
    VT_UI8 = 21,
    VT_INT = 22,
    VT_UINT = 23,
    VT_VOID = 24,
    VT_HRESULT  = 25,
    VT_PTR = 26,
    VT_SAFEARRAY = 27,
    VT_CARRAY = 28,
    VT_USERDEFINED = 29,
    VT_LPSTR = 30,
    VT_LPWSTR = 31,
    VT_FILETIME = 64,
    VT_BLOB = 65,
    VT_STREAM = 66,
    VT_STORAGE = 67,
    VT_STREAMED_OBJECT = 68,
    VT_STORED_OBJECT = 69,
    VT_BLOB_OBJECT = 70,
    VT_CF = 71,
    VT_CLSID = 72,
    VT_VECTOR = 0x1000,
    VT_ARRAY = 0x2000,
    VT_BYREF = 0x4000,
    VT_RESERVED = 0x8000,
};
*/

#define __RPCNDR_H_VERSION__
#define COM_NO_WINDOWS_H

#endif
