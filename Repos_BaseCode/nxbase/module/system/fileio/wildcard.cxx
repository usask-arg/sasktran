/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_core.h"


//----------------------------------------------------------------------------
//			CheckRestOfMatch
//	Check that the wildspec string occurs at the first character of string
//	Remaining.  If it does then return the pointer to the next free char in
//	remaining, otherwise return NULL.  Code allows wildcard '?' to mean
//	anything (except '\0').
//----------------------------------------------------------------------------

static const char *CheckRestOfMatch( const char *wildspec, const char *Remaining )
{
   if (Remaining == NULL) return NULL;
   if (wildspec  == NULL) return NULL;
   
	nxBOOL matchok = nxTRUE;
	char c,b;
	do
	{
		c = *Remaining++;
		b = *wildspec++;
		if (b != '\0')
		{
			matchok = (c != '\0') && ( b == '?' || (b == c) );
		}
	} while ( matchok && (b != '\0'));
	if (!matchok) Remaining = NULL; else Remaining--;
	return Remaining;
}
 
//----------------------------------------------------------------------------
//			FindFirstMatch
//	Find the first occurence of character Buffer in string Remaining.
//	If a match is found then return the pointer to the start of the
//	match, otherwise return NULL.  The wildcard '?' matches anything.
//----------------------------------------------------------------------------

static const char *FindFirstMatch( const char Buffer, const char *Remaining )
{
	if ((Buffer == '\0') || (Remaining == NULL)) return NULL;   // Watchout for null strings etc.
	if (Buffer == '?') return Remaining;						// The '?' wildcard matches anything, so found first match

	char c;														// Not '?' so search for first occurence.
	nxBOOL found; // = nxFALSE;       	                			// default to not found
	do															// so scan along					
	{															// the string
		c = *Remaining++;										// get the character from the string
		found = (c == Buffer);									// see if its the same
	}															// and
	while (( c != '\0') && !found);								// loop until end of string or match found
	if (found) Remaining--; else Remaining = NULL;				// if match then adjust pointer else set to NULL
	return Remaining;											// return pointer to match (or NULL)
}
 
//----------------------------------------------------------------------------
//			BufferIsSomeWhere
//	returns nxTRUE if the string in wildspec is somewhere in string
//	RestOfThis.  '?' wildcards only are allowed. 
//
//	If Floating is nxTRUE then serach will scan from the string for the
//	wildcard string.  If floating is nxFALSE then the search will start
//	at the first character of the string.
//
//	If a successful match is found then the RestOfThis pointer is updated
//	to point to the last character beyond the match.
//----------------------------------------------------------------------------

static const char *BufferIsSomeWhere( const char * wildspec, const char * Remaining, nxBOOL Floating)
{
   const char *endofmatch = NULL;
   nxBOOL  Found      = nxFALSE;										// and default to not found

	if (Floating)														// if I'm allowed to look around for a match
	{																	// then 
		const char *firstmatch = Remaining;								// start at the beginning of the string
		do 																// then repeatedly
		{             													// look for
			firstmatch = FindFirstMatch( *wildspec, firstmatch );		// where the first match occurs if at all
			if (firstmatch)												// if first character matches
			{															// then
				endofmatch = CheckRestOfMatch( wildspec, firstmatch );	// check the rest of the match
				Found = (endofmatch != NULL);          					// did we find a match.
				firstmatch++;                                			// and move along from the last first match
			}															// and repeat until 
		} while (!Found && (firstmatch != NULL) );						// we get full match or no match at all
	}           														// and that i sthe end of floating search                                 			
	else																// otherwise if fixed search was requested
	{																	// then simply
		endofmatch = CheckRestOfMatch( wildspec, Remaining );			// check match from first character
	}																	// and that is that
	return endofmatch;													// return pointer if good else return NULL.
}

//----------------------------------------------------------------------------
//			nxWildcard::constructor
//----------------------------------------------------------------------------

nxWildcard::nxWildcard( const char *spec)
{
   wildspec = NULL;
   SetSpec( spec );
}

//----------------------------------------------------------------------------
//			nxWildcard::SetSpec
//----------------------------------------------------------------------------

void nxWildcard::SetSpec( const char *spec )
{
   int n;
   if (spec == NULL)
   {
      wildspec = new char [2];
      strcpy(wildspec, "" );
   }
   else
   {
      n = (int)strlen(spec);
      wildspec = new char [n+1];
      strcpy( wildspec, spec );
   }

}

//----------------------------------------------------------------------------
//			nxWildcard::destructor
//----------------------------------------------------------------------------

nxWildcard::~nxWildcard()
{
   if( wildspec != NULL ) delete [] wildspec;
}

//----------------------------------------------------------------------------
//			nxWildcard::Matches
//	compares ThisString to the wildspec and returns nxTRUE if ThisString
//	is compatible with the wildcard. Otherwise it returns nxFALSE
//----------------------------------------------------------------------------


nxBOOL nxWildcard::Match( const char *ThisString)
{

	char       *buffer;
	const char *RestOfThis;
	int         bufferidx = 0;
	int         i = 0;

	if (ThisString == NULL) return nxFALSE;
	int   maxw = (int)strlen(wildspec);											// get the length of the  wildcard
	if    (maxw < 1) return nxFALSE;										// if its less than 1 byte then no mathc found

	nxBOOL  ok       = nxTRUE;                                              // used to flag that the comparison is ok.
	nxBOOL  Floating = nxFALSE;												// signifies a '*' has been encountered and enables a floating scan for sub-strings.

	buffer     = new char[ maxw + 1];
	*buffer     = '\0';
	RestOfThis = ThisString;

	while ( (i < maxw) && ok) 												// Scan through the wildcard
	{																		// and find
		switch (wildspec[i])												// the places where there are
		{																	// '*' characters.
			case '*' :  Floating = nxTRUE;									// nxWildcards make the search float
         				if (bufferidx !=0 )									// if this is the second '*', then check out the buffered data.
         				{													// if this is somewhere in rest of string
         					RestOfThis = BufferIsSomeWhere( buffer, RestOfThis, Floating);	// Sees if buffer is in rest of strinng
         					ok         = RestOfThis != NULL;
         				}
         				bufferidx  = 0;										// reset the buffer index to 0
         				*buffer    = '\0';									// and terminate the buffer string
         				break;
         
			default :	buffer[bufferidx++] = wildspec[i];					// all other characters in the wildcard
						buffer[bufferidx]   = '\0';							// are simply added to the buffer.
						break;												// for later comparisons.
		}
		i++;																// point to next character in string.
	}
   
	if (bufferidx != 0)														// if the buffer has characters in it
	{																		// then
		RestOfThis = BufferIsSomeWhere( buffer, RestOfThis, Floating ); 	// see if the buffered wildcard is in the string
		ok = (RestOfThis != NULL) && (*RestOfThis == '\0' );				// If it is make sure its at the end of the string
	}																		// else if bufferidx == 0, then last char was a '*', so it matches.
   
	delete [] buffer;
	return ok;
}

/*
int main( int argc, char *argv[] )
{

  nxBOOL check;
  WildCard  spec( "*345.*?0");
  
  check = spec.Matches( "P012345.TN0");
  cout << "check = " << check << "\n";
  return check;
}
   
*/
