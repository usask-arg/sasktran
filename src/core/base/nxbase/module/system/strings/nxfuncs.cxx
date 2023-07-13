/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_core.h"

//---------------------------------------------------------------------------
//			function nxStrtok
//	parses a const string and writes the result to a nxStringArray.
//	The string passed in can either be a nxString or a const char *
//---------------------------------------------------------------------------

int nxStrtok( const char *utstring, nxStringArray *List, const char *separator ) 
{
   char *aline;
   char *token;
   static const char defsep[] = " ,\t\n";
   const char *sep;
   nxString StringToken;

   aline = new char [strlen(utstring) + 1];		// copy the string into a temporary buffer
   strcpy( aline, utstring );					// so we can pass it to strtok without problem
   List->RemoveAll();							// Remove all the elements from the output array
   if (separator == NULL) sep = defsep; else sep = separator;

   token = strtok( aline, sep);					// now get the first token in the list
   while (token != NULL )						// if we got a token
   {											// then
      StringToken = token;						// copy the token to a string
      List->Add(StringToken);					// and add the string to our array
      token = strtok( NULL, sep);				// now look fo the next token
   }											// repeat for all of the tokens
   delete [] aline;								// delete the buffer space
   return List->GetSize();						// return the number of tokens
}

