/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_core.h"


void nxTrace(const char* lpszFormat, ...)
{
	va_list args;
	va_start(args, lpszFormat);

#if defined(NX_WINDOWS) || defined(_WIN32)
	int nBuf;
	char szBuffer[1024];

	nBuf = _vsnprintf(szBuffer, sizeof(szBuffer), lpszFormat, args);
	szBuffer[sizeof(szBuffer)-1] = '\0';
	OutputDebugStringA(szBuffer);
#else
	//vprintf( lpszFormat, args );
#endif
	va_end(args);
}
