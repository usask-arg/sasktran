#if defined(_DEBUG) && defined(SWIG_PYTHON_INTERPRETER_NO_DEBUG)
/* Use debug wrappers with the Python release dll */
# undef _DEBUG
# include <Python.h>
# define _DEBUG 1
#else
# include <Python.h>
#endif
#include "nxlib.h"
#include "nxlogpython.h"


//---------------------------------------------------------------------------
//						nxLogMATLAB::nxLogMATLAB
//---------------------------------------------------------------------------

nxLogPython::nxLogPython()
{
	CheckDefaultLogger();
	this->AddRef();
}


//----------------------------------------------------------------------------
//			nxLogPython::DisplayEntry()
//	A class for writing Log information to the screen.
//----------------------------------------------------------------------------

void nxLogPython::DisplayEntry( const nxLogEntry& entry )
{
	nxString	TheLine;

	TheLine = entry.StatusString( nxTRUE );
	TheLine = TheLine + entry.Message();
	PySys_WriteStdout("%s\n",(const char*)TheLine );
}

