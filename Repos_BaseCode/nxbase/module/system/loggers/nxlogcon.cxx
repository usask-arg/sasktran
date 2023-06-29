/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_core.h"
#include <iostream>
#include <fstream>

nxLogConsole::nxLogConsole()
{
	AddRef();
	m_showUT = nxFALSE;
	CheckDefaultLogger();
}

void nxLogConsole::ShowUT( nxBOOL show )
{
	m_showUT = show;
}

//----------------------------------------------------------------------------
//			nxLogConsole::DisplayEntry()
//	A class for writing Log information to the screen.
//----------------------------------------------------------------------------

void nxLogConsole::DisplayEntry( const nxLogEntry& entry )
{
	nxString	TheLine;

	TheLine = entry.StatusString( nxTRUE ) + entry.Message();
	if (m_showUT) TheLine += " [" + entry.Mjd()+"]";
	std::cout.flush();
	fprintf( stdout, "%s\n", (const char*)TheLine);
	fflush(stdout);
}

