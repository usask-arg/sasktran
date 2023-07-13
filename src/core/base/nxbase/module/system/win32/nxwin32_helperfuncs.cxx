#include "nxbase_core.h"

//---------------------------------------------------------------------------
//						nxGetClsID
//	Get the CLSID given the ProgID
//---------------------------------------------------------------------------

nxBOOL nxGetClsid( CLSID* cls, const char *pgmname )
{
	WCHAR*			unicode;
	HRESULT			Status;

	int len = (int)strlen(pgmname);
	unicode = new WCHAR[len+2];

	len    = MultiByteToWideChar( CP_ACP, MB_PRECOMPOSED, pgmname,-1, unicode, len+2 );   
	Status = CLSIDFromProgID( unicode, cls);
	delete [] unicode;
	return SUCCEEDED(Status);
}


void  Wait( int n )
{
	nxTimeStamp	tend;
	nxTimeStamp	tnow;

	tend.FromSystem();
	tend = tend + (n/1000.0)*nxTimeStamp::ONESECOND;
	do
	{
		YieldToSystem( 10 );
		tnow.FromSystem();
	} while (tnow < tend);
}

//---------------------------------------------------------------------------
//			YieldToSystem()
//	Yield control to the operating system for approximately N PeekMessages.
//---------------------------------------------------------------------------	

void YieldToSystem( int N )
{
   MSG       msg;

   if (N < 2) N = 2;
   for (int i = 0; i < N; i++)
   {
      while (PeekMessage( &msg, NULL, 0, 0, PM_REMOVE))
      {
			TranslateMessage( &msg  );
			DispatchMessage ( & msg );
      }
   }
}


//---------------------------------------------------------------------------
//						nxStringFromUnicode
//---------------------------------------------------------------------------

void nxStringToBSTR( nxString& str, BSTR* bstr )
{
	*bstr = ::SysAllocStringByteLen( str.DangerousTypecast(), str.GetLength()+1 );
	::MultiByteToWideChar( CP_ACP, 0, str, str.GetLength()+1, *bstr, str.GetLength()+1);
}


