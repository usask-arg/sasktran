/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
//#include "stdafx.h"
#include "nxbase_core.h"


//----------------------------------------------------------------------------
//			nxFileSpec::Constructor
//----------------------------------------------------------------------------

nxFileSpec::nxFileSpec( const char * filenam)
{
   *this = filenam;
}

//----------------------------------------------------------------------------
//			nxFileSpec::operator=
//----------------------------------------------------------------------------

const char *nxFileSpec::operator= (const char *lpnewname)
{
	nxString	 newname( lpnewname );

	if (newname.IsEmpty())
	{
		m_directory.Empty(nxFALSE);
		m_drive.Empty(nxFALSE);
		m_name.Empty(nxFALSE);
		m_extension.Empty(nxFALSE);
	}
	else
	{
		int	   idx;
		bool	noextension;

		newname.MakeReverse();
		idx = newname.FindAnyOf(".\\/");									// find the extension character or a directory delimiter
		noextension = (idx < 0) || newname[idx] != '.';						// no extension if we have not '.' or the '.' is in a directory specification
		if (noextension)
		{
			m_extension.Empty(nxFALSE);									// if not found then clear the extension
		}
		else																// otherwise
		{																	// extract the
			m_extension = newname.Left( idx+1 );								// and now remove the extension from the string
			m_extension.MakeReverse();
			newname	    = newname.Right( newname.GetLength() - idx-1 );		// extension from the end of the string including the '.'
		}																	// and that is that
		newname.MakeReverse();

	// This is stuff just for windows
		idx = newname.Find(':');											// find the drive specification on the left of the string
		if (idx < 0)														// if we didnt find ':'
		{
			idx = newname.Find("\\\\");
			if (idx < 0)
			{
				m_drive.Empty(nxFALSE);												// if we did not find it then clear the drive spec
			}
			else
			{
				m_drive = newname.Left(idx+2);
				newname	= newname.Right( newname.GetLength() - idx - 2);
				idx = newname.Find("\\");
				if (idx < 0)
				{
					m_drive += newname + "\\";
					newname.Empty(nxFALSE);
				}
				else
				{
					m_drive += newname.Left(idx+1);
					newname = newname.Right( newname.GetLength() - idx - 1 );
				}
			}
		}
		else																// otherwise
		{																	// get the
			m_drive = newname.Left( idx + 1 );								// drive specification
			newname = newname.Right( newname.GetLength() - idx -1 );		// and remove the drive spec from the string
		}


		newname.MakeReverse();												// reverse the remaining string
		idx  = newname.FindAnyOf("/\\");									// find a directory character
		if (idx >= 0 )														// if we found  a directory symbol
		{																	// then
			m_name      = newname.Left( idx );								// get the reversed name
			m_directory = newname.Right( newname.GetLength() - idx );		// get the reversed directory
			m_directory.MakeReverse();										// reverse the directory
		}																	// and that is that
		else																// otherwise
		{																	// no directory found
			m_directory.Empty(nxFALSE);											// so clear the directory structure
 			m_name = newname;												// copy the reverse name
		}																	// and thats that
		m_name.MakeReverse();												// reverse the name
	}
	return lpnewname;														// and return the pointer passed in.
}




