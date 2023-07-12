#include "nxbase_math.h"
#include "sklapack.h"
#include <boost/thread.hpp>

#ifdef NX_WINDOWS
	#pragma comment( lib, "netlib")
	#pragma message("Automatically linking netlib.lib and netlib.dll into this project")
#else
	#define DLGAMA dlgama_
#endif

F77EXTERNC f77doublereal F77CALL DLGAMA( f77doublereal *x);

static double dlgama_helper( double *x)
{
	static boost::mutex	threadlock;

	boost::lock_guard<boost::mutex> lock(threadlock);
	double y =  DLGAMA(x);
	return y;
}

namespace nxnetlib
{
	double (* dlgama)( double *x) = &dlgama_helper;
};


