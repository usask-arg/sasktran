#include "Python.h"			// MATLAB include file
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

