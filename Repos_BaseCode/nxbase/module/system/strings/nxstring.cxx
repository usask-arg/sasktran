/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_core.h"

static const char*  const NULLSTR = "\0";
#include <ctype.h>


static char *nxstrlwr( char *string )
{
	char *ptr = string;
	if (ptr == NULL) return NULL;
	while (*ptr != '\0')
	{
		*ptr = (char)tolower( *ptr );
		ptr++;
	}
	return string;
}

static char *nxstrupr( char *string )
{
	char *ptr = string;
	if (ptr == NULL) return NULL;
	while (*ptr != '\0')
	{
		*ptr = (char)toupper( *ptr );
		ptr++;
	}
	return string;
}

static char *nxstrrev( char *string )
{
	if (string == NULL) return NULL;
	size_t n      = strlen(string);
	char * source = string + n-1;
	char * ptr    = string;
	char c;

	while (source > ptr)
	{
		c = *ptr;
		*ptr++    = *source;
		*source-- = c;
	}
	return string;
}


//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::Default blank constructor
//	Allocate a blank NULL string.
//-
//---------------------------------------------------------------------------

nxString::nxString()
{
	m_localstorage[0] = '\0';
	m_length     = 0;
	m_str        = &m_localstorage[0];
	m_allocation = N_ELEMENTS(m_localstorage);;
}

nxString::nxString( const char *str )
{
	m_localstorage[0] = '\0';
	m_str        = &m_localstorage[0];
	m_length     = 0;
	m_allocation = N_ELEMENTS(m_localstorage);
	*this = str;
}
 
nxString::nxString( const nxString& otherstring )
{
	m_localstorage[0] = '\0';
	m_str        = &m_localstorage[0];
	m_length     = 0;
	m_allocation = N_ELEMENTS(m_localstorage);
	*this = otherstring;
}
//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::Destructor
//	Ensure that the memory is properly deallocated.
//-
//---------------------------------------------------------------------------

nxString::~nxString()
{
	Empty(nxTRUE);
}


/*-----------------------------------------------------------------------------
 *					nxString::SetToNullStr		2004-12-1*/
/** **/
/*---------------------------------------------------------------------------*/

void nxString::SetToNullStr()
{
	if (m_str != NULL) m_str[0] = '\0';
	m_length = 0;
}

//---------------------------------------------------------------------------
//						nxString::Empty						
//---------------------------------------------------------------------------

void nxString::Empty( nxBOOL clearcache )
{
	if (clearcache)
	{
		if ((m_str != &m_localstorage[0]) && (m_str != NULL))
		{
			delete [] m_str;
		}
		m_str        = &m_localstorage[0];
		m_allocation = N_ELEMENTS(m_localstorage);
	}
	SetToNullStr();
}

//---------------------------------------------------------------------------
//						nxString::CheckAllocation
//	Checks the allocation of the memory buffer.  If the user requests that
//	the current string be saved then it is done.
//---------------------------------------------------------------------------

nxBOOL nxString::CheckAllocation( size_t nbytes, nxBOOL savestring )
{
	nxString	temp;
	nxBOOL		ok;
	size_t		oldallocation;

	ok = (nbytes <= m_allocation);
	if (!ok) 										// if the requested size is greater than this allocation
	{												// then
		oldallocation = m_allocation + 100;			// Allocate a decent chunk of memory (1 byte increment is too inefficient).
		if (savestring) temp = *this;				// if the user requests the string be saved then make a temp copy
		Empty(nxTRUE);								// clear the current string
		if (nbytes < oldallocation) nbytes = oldallocation;				// may as well allocate decent size string (enough to hold a filename perhaps?)
		m_str = new char [nbytes];					// and create the new segment
		ok    = (m_str != NULL);					// check that the string was ok.
		if (ok)										// if it was						
		{											// then
			m_allocation = nbytes;					// allocate the string
			if (savestring) *this = temp;			// and copy back if saving of source required.
		}
		else
		{
			m_allocation = N_ELEMENTS(m_localstorage);
			m_str        = &m_localstorage[0];   
			SetToNullStr();
			NXTHROW(("nxString memory allocation error"));
		}
	}
	return ok;
}
//---------------------------------------------------------------------------
//						nxString::CopyString
//---------------------------------------------------------------------------

nxBOOL nxString::CopyString( const char *str, size_t nbytes )
{
	nxBOOL	ok;

	if (str == NULL)
	{
		SetToNullStr();
		ok       = nxFALSE;
	}
	else
	{
		if (nbytes == 0) nbytes = strlen(str) + 1;
		ok = CheckAllocation( nbytes );
		if (ok)
		{
			strcpy( m_str, str );
		}
		m_length = nbytes-1;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxString::ReserveAndErase		2005-8-2*/
/** Reserve space for nchars charcaters and erase the original string
 **/
/*---------------------------------------------------------------------------*/

nxBOOL nxString::ReserveAndErase( int nchars )
{
	nxBOOL ok;
	
	ok =  CheckAllocation( (size_t)nchars);
	m_length = 0;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					nxString::UpdateLengthAfterExternalWrite		2005-8-2*/
/** This function must be called if an external module has just written
 *	string into our buffer
 **/
/*---------------------------------------------------------------------------*/

size_t nxString::UpdateLengthAfterExternalWrite()
{
	m_length = strlen(m_str);
	if (m_length > m_allocation)
	{
		NXTHROW(( "nxString::UpdateLengthAfterExternalWrite, The external write exceeded the allocated buffer"));
	}
	return m_length;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::operator=
//	Copies a nxString from one string to another.
//-
//---------------------------------------------------------------------------

nxString&	nxString::operator=( const nxString &otherstring)
{
	int		nbytes;

	nbytes = otherstring.GetLength() + 1;
	CopyString( otherstring, nbytes );
	return *this;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::operator=
//	Copies a char * string to a nxString
//-
//---------------------------------------------------------------------------

nxString&	nxString::operator=( const char * otherstring)
{
	CopyString( otherstring );
	return *this;
}


/*---------------------------------------------------------------------------
 *'					nxString::operator=                             2003-6-11
 *-------------------------------------------------------------------------*/

nxString& nxString::operator=(  const WCHAR* widestring )
{
	nxStringFromUnicode( this, widestring );
	return *this;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::operator=
//	Copies a char * string to a nxString
//-
//---------------------------------------------------------------------------

nxString&	nxString::operator=( char  achar)
{
	char dummy[2];

	dummy[0] = achar;
	dummy[1] = '\0';
	CopyString( dummy, 2 );
	return *this;
}
//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::MakeLower
//	Makes the string lower case.
//-
//---------------------------------------------------------------------------

void nxString::MakeLower()
{
	if (m_length > 0) nxstrlwr( m_str );
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::MakeReverse
//	Reverses the order of this string
//-
//---------------------------------------------------------------------------

void nxString::MakeReverse()
{
	if (m_length > 1)  nxstrrev( m_str );
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::MakeUpper
//	converts the entire string to uppercase
//-
//---------------------------------------------------------------------------

void nxString::MakeUpper()
{
	if (m_length > 0) nxstrupr( m_str );
}


//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::GetAt
//	Returns the character at a specific location in the string
//-
//---------------------------------------------------------------------------

char nxString::GetAt( size_t nindex ) const
{
	char c;

	if ((nindex >= 0) && ( nindex < m_length )) c = m_str[nindex]; else c = '\0';
	return c;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::operator []
//	returns the characters at location indexed by array operator
//-
//---------------------------------------------------------------------------

char nxString::operator [] ( int nindex) const
{
	return GetAt(nindex);
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::SetAt
//	Sets the string at a specific location to the character "ch"
//	if ch is the null terminator then the length of the string
//	is adjusted accordingly.
//-
//---------------------------------------------------------------------------

void nxString::SetAt( size_t nindex, char ch )
{
	if ((nindex >= 0) && ( nindex < m_length)) m_str[nindex] = ch;
	if (ch == '\0') m_length = nindex;
}

//---------------------------------------------------------------------------
//						nxString::concatenate
//		make the string by concatenating two string pointers.
//---------------------------------------------------------------------------

void nxString::concatenate( const char * str1, size_t l1, const char *str2, size_t l2 )
{
	nxBOOL	ok;
	char*	outstr;

	ok = CheckAllocation( l1 + l2 + 1 );						// get the total allocation required
	if (ok)														// if it worked
	{															// then
		strcpy( m_str, str1 );									// copy the first string
		outstr = m_str + l1;									// get the start of the second string
		strcpy( outstr, str2 );									// and add the second string
		m_length = l1 + l2;										// and set the length of the string
	}
}		

//---------------------------------------------------------------------------
//+
//NAME:
//						CString operator + 
//-
//---------------------------------------------------------------------------

 nxString operator + (const nxString& str1, const nxString & str2 )
{
	nxString s;
	s.concatenate( str1, str1.GetLength(), str2, str2.GetLength());
	return s;
}

 nxString operator + (const nxString& str1, char ch )
{
	nxString s;
	char	 dummy[2];

	dummy[1] = '\0';
	dummy[0] = ch;
	s.concatenate( str1, str1.GetLength(), dummy, (int)strlen(dummy));
	return s;
}

 nxString operator + (char ch, const nxString& str1 )
{
	nxString s;
	char	 dummy[2];

	dummy[1] = '\0';
	dummy[0] = ch;
	s.concatenate( dummy, (int)strlen(dummy), str1, str1.GetLength());
	return s;
}

 nxString operator + ( const nxString& str1,  const char* str2 )
{
	nxString s;

	if (str2 == NULL) str2 = NULLSTR;
	s.concatenate( str1, str1.GetLength(), str2, (int)strlen(str2) );
	return s;
}

 nxString operator + (const char* str1, const nxString& str2 )
{
	nxString s;

	if ( str1 == NULL) str1 = NULLSTR;
	s.concatenate( str1, (int)strlen(str1), str2, str2.GetLength());
	return s;
}

//---------------------------------------------------------------------------
//					nxString::AddAString
//---------------------------------------------------------------------------

void nxString::AddAString( const char* str, size_t l1 )  
{
	nxBOOL	ok;

	ok = CheckAllocation( l1 + GetLength() + 1, nxTRUE );			// get the total allocation required amd save the orginal string
	if (ok)														// if it worked
	{															// then
		strcat( m_str, str );									// copy the first string
		m_length = l1 + m_length;								// and set the length of the string
	}
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::operator +=
//-
//---------------------------------------------------------------------------

nxString & nxString::operator += ( const nxString& str )
{
	AddAString( str, str.GetLength() );
	return *this;
}

nxString & nxString::operator += ( char ch )
{
	char dummy[2];

	dummy[0] = ch;
	dummy[1] = '\0';
	AddAString( dummy, (int)strlen(dummy) );
	return *this;
}

nxString & nxString::operator += ( const char* str )
{
	if (str != NULL) AddAString( str, (int)strlen(str));
	return *this;
}


static int safestrcmp( const char* s1, const char * s2 )
{
	if (s1 == NULL) s1 = NULLSTR;
	if (s2 == NULL) s2 = NULLSTR;
	return strcmp(s1,s2);
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::operator ==
//-
//---------------------------------------------------------------------------
 nxBOOL operator == ( const nxString& s1, const nxString& s2 )
{
	return safestrcmp( s1, s2 ) == 0;
}
		
 nxBOOL operator == ( const nxString& s1, const char *  s2 )
{
	return safestrcmp( s1, s2 ) == 0;
}

 nxBOOL operator == ( const char * s1, const nxString& s2 )
{
	return safestrcmp( s1, s2 ) == 0;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::operator !=
//-
//---------------------------------------------------------------------------
 nxBOOL operator != ( const nxString& s1, const nxString& s2 )
{
	return safestrcmp( s1, s2 ) != 0;
}
		
 nxBOOL operator != ( const nxString& s1, const char *  s2 )
{
	return safestrcmp( s1, s2 ) != 0;
}

 nxBOOL operator != ( const char * s1, const nxString& s2 )
{
	return safestrcmp( s1, s2 ) != 0;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::operator <
//-
//---------------------------------------------------------------------------

 nxBOOL operator < ( const nxString& s1, const nxString& s2 )
{
	return safestrcmp( s1, s2 ) < 0;
}
		
 nxBOOL operator < ( const nxString& s1, const char *  s2 )
{
	return safestrcmp( s1, s2 ) < 0;
}

 nxBOOL operator < ( const char * s1, const nxString& s2 )
{
	return safestrcmp( s1, s2 ) < 0;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::operator >
//-
//---------------------------------------------------------------------------

 nxBOOL operator > ( const nxString& s1, const nxString& s2 )
{
	return safestrcmp( s1, s2 ) > 0;
}
		
 nxBOOL operator > ( const nxString& s1, const char *  s2 )
{
	return safestrcmp( s1, s2 ) > 0;
}

 nxBOOL operator > ( const char * s1, const nxString& s2 )
{
	return safestrcmp( s1, s2 ) > 0;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::operator <=
//-
//---------------------------------------------------------------------------

 nxBOOL operator <= ( const nxString& s1, const nxString& s2 )
{
	return safestrcmp( s1, s2 ) <= 0;
}
		
 nxBOOL operator <= ( const nxString& s1, const char *  s2 )
{
	return safestrcmp( s1, s2 ) <= 0;
}

 nxBOOL operator <= ( const char * s1, const nxString& s2 )
{
	return safestrcmp( s1, s2 ) <= 0;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::operator >=
//-
//---------------------------------------------------------------------------

 nxBOOL operator >= ( const nxString& s1, const nxString& s2 )
{
	return safestrcmp( s1, s2 ) >= 0;
}
		
 nxBOOL operator >= ( const nxString& s1, const char *  s2 )
{
	return safestrcmp( s1, s2 ) >= 0;
}

 nxBOOL operator >= ( const char * s1, const nxString& s2 )
{
	return safestrcmp( s1, s2 ) >= 0;
}


//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::sprintf
//	Emulates the sprintf function, but allocates a large (1024 byte)
//	buffer to ensure that memory does not overflow.
//-
//---------------------------------------------------------------------------

void nxString::sprintf( const char *format, ... )
{
	va_list ArgPtr;
	char	   buffer[1024] = "";
	int n = -9999;
	
	va_start(ArgPtr, format);
   
	try
	{
#if defined(NX_WINDOWS)
		n = ::vsprintf_s(buffer, N_ELEMENTS(buffer)-1, format, ArgPtr );
#else
		n = ::vsnprintf(buffer, N_ELEMENTS(buffer)-1, format, ArgPtr );
#endif
	}
	catch (...)
	{
		n = -9998;
	}
	if ( (n < 0) || (n >=N_ELEMENTS(buffer) ))
	{
		if (format == nullptr) fprintf( stderr, "nxString:sprintf, exception executing vsprintf, error code = %d. format parameter is NULL\n", (int)n);
		else                   fprintf( stderr, "nxString:sprintf, exception executing vsprintf, error code = %d. format paramater =[%s]\n", (int) n, (const char*)format);
	}
   va_end(ArgPtr);
   *this = buffer;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::vsprintf
//	Emulates the vsprintf function, but allocates a large (1024 byte)
//	buffer to ensure that memory does not overflow
//-
//---------------------------------------------------------------------------

void nxString::vsprintf( const char* format, va_list argptr)
{
	char buffer[1024];

	::vsprintf( buffer, format, argptr );
	*this = buffer;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::Find
//	Find the first occurence of string charstr or character c in this
//	nxString object.  Returns the index of the first character that
//	mathes.  Returns -1 if no macth is found.
//-
//---------------------------------------------------------------------------

int nxString::Find( char c )
{
	char	*ptr;
	int		idx;
 	ptr = strchr( m_str, (int)c); 
 	if (ptr == NULL)
 	{
 		idx = -1;
 	}
 	else
 	{
 		idx = (int) (ptr-m_str);
 	}
 	return idx;
}


/*-----------------------------------------------------------------------------
 *					nxString::FindAnyOf		2009-6-4*/
/** Find the first occurrence of any of the characters in cstr
 **/
/*---------------------------------------------------------------------------*/

int nxString::FindAnyOf( const char* cstr )
{
	int		idx = -1;
	int		idx1;

	while ( (*cstr != '\0') && (idx != 0) )					// While we have more characters to search 
	{														// and we have not found a character in the first character of the string
 		idx1 = Find( *cstr++ );								// the look for this character
		if (idx1 >= 0)										// if we found it 
		{													// then
			if (idx < 0) idx = idx1;						// if we have not yet found any characters then simply assign
			else         idx = (idx1 < idx) ? idx1 : idx;	// otherwise get the minimum of idx and idx1
		}
	}
	return idx;
}

/*---------------------------------------------------------------------------
 *'					nxString::Find}                                 2003-6-27
 *-------------------------------------------------------------------------*/

int nxString::Find( const char* charstr)
{
	char	*ptr;
	int		idx;
 	ptr = strstr( m_str, charstr); 
 	if (ptr == NULL)
 	{
 		idx = -1;
 	}
 	else
 	{
		idx = (int)(ptr-m_str);
 	}
 	return idx;
}


//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::TruncateAt
//	Truncates the string at the first occurence of any character
//	in string charset.  Uses the strcspn function.  Very useful
//	for eliminating comments from the ends of lines etc.
//	Returns the length of the string
//-
//---------------------------------------------------------------------------

size_t nxString::TruncateAt( const char* charset)
{
	size_t idx;
	
	if (m_length > 0)
	{
		idx        = strcspn( m_str, charset );
		m_str[idx] = '\0';
		m_length   = idx;
	}
	return m_length;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::RemoveWhiteSpace
//	Removes leading and trailing white spaces from the text.
//	White space is any non-printable char (<= SPACE char)
//	returns the length of the new string
//-
//---------------------------------------------------------------------------

size_t nxString::RemoveWhiteSpace()
{
	size_t i;
	char*	newstr;
	char*	outstr;

	if (m_length > 0)
	{
		for (i = m_length;  i > 0; )			// scan backwards through the file
		{											// until we find
			--i;
			if (m_str[i] > ' ')						// a non-printable char
			{										// so now
				m_str[i+1] = '\0';					// terminate the string at this point
				break;								// and break out of the loop
			}										// otherwise loop until a printable char
		}											// is found or we hit the end 
		if ( i == 0)								// IF we dont terminate until the last character
		{											// then
			if ( m_str[0] <= ' ') m_str[0] = '\0';	// see if its a white space character.
		}											// and terminate if it is.
	}
	newstr = m_str;								// get the start of the string
	while (*newstr != '\0')						// while we have more characters
	{                            				// then
		if (*newstr > ' ') break;				// hunt for non-printables
		newstr++;								// increment through loop until found;
	}

	outstr = m_str;								// now copy from the original string
	while (*newstr != '\0')						// to the beginning of the original
	{											// I have avoided strcpy as I'm
		*outstr++ = *newstr++;					// not sure how it handles overlapped			
	}											// strings on all platforms
	*outstr = '\0';								// terminate the copied string
	m_length = (outstr - m_str);					// and reset the length of this string.
	return m_length;
}												// and that is that.

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::Replace
//	replace all occurences of character "oldchar" with character "newchar"
//	If newchar is the null terminator then the string length is readjusted.
//-
//--------------------------------------------------------------------------- 

size_t nxString::Replace( char oldchar, char newchar )
{
	char*	str;
	size_t	i = 0;
	
	str = m_str;
	while (*str != '\0' )
	{
		if (*str == oldchar)
		{
			i++;
			*str = newchar;
		}
		str++;
	}
	if (newchar == '\0') m_length = (int)strlen(m_str);
	return i;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::Mid
//	Similar to the BASIC MID$ function returns a sub string from
//	another string.  The function starts at character offset nfirst and
//	generates a new string upto ncount bytes long and returns the value
//	as a new nxString
//-
//---------------------------------------------------------------------------

nxString nxString::Mid( size_t nfirst, size_t ncount )
{
	char*	instr;
	char*	outstr;
	char*	outs;
	size_t	 i;
	nxString result;
	
	if (m_length == 0 || ncount < 1 || nfirst >= m_length) return "";
	
	outstr = new char[ncount+1];
	if (outstr == NULL) NXTHROW(("nxString::Mid memory allocation error"));
	instr  = m_str + nfirst;
	i      = 0;
	outs   = outstr;
	
	while ((*instr != '\0') && ( i < ncount))
	{
		*outs++ = *instr++;
		i++;
	}
	*outs = '\0';
	result = outstr;
	delete [] outstr;
	return result;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::Left
//	Similar to the BASIC Left function returns a sub string from
//	another string.  The function returns the left most characters upt ncount
//-
//---------------------------------------------------------------------------

nxString nxString::Left( size_t ncount )
{
	char*	instr;
	char*	outstr;
	char*	outs;
	size_t	 i;
	nxString result;
	
	if (m_length == 0 || ncount < 1) return "";
	if (ncount > m_length) ncount = m_length;
	outstr = new char[ncount+1];
	if (outstr == NULL) NXTHROW(("nxString::Left memory allocation error"));
	outs   = outstr;
	instr  = m_str;
	for (i = 0; i < ncount; i++) *outs++ = *instr++;
	*outs = '\0';
	result = outstr;
	delete [] outstr;
	return result;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxString::Right
//	Similar to the BASIC Right function returns a sub string from
//	another string.  The function returns the right most characters
//	upto ncount
//-
//---------------------------------------------------------------------------

nxString nxString::Right( size_t ncount )
{
	char*		instr;
	char*		outstr;
	char*		outs;
	size_t		i;
	nxString	result;
	
	if (m_length == 0 || ncount < 1) return "";
	if (ncount > m_length) ncount = m_length;
	outstr = new char[ncount+1];
	if (outstr == NULL) NXTHROW(("nxString::Right memory allocation error"));
	outs   = outstr;
	instr  = m_str + (m_length - ncount);
	for (i = 0; i < ncount; i++) *outs++ = *instr++;
	*outs = '\0';
	result = outstr;
	delete [] outstr;
	return result;
}


const char * nxString::AsFormattedNumber( int value )
{
#if defined(NX_WINDOWS)
	
	char      chrbuf[400];
	nxString  str;

	NUMBERFMT	fmt;
	fmt.NumDigits = 0;
	fmt.LeadingZero = 0;
	fmt.Grouping    = 3;
	fmt.lpDecimalSep = ".";
	fmt.lpThousandSep = ",";
	fmt.NegativeOrder = 0; 
	str.sprintf("%d",(int)value );
	::GetNumberFormat( LOCALE_USER_DEFAULT, 0, (const char *)str, &fmt, chrbuf, sizeof(chrbuf)-1 );  
	*this = chrbuf;
#else
	this->sprintf("%d",(int)value);
#endif
	return  (const char *)(*this);
}


/*---------------------------------------------------------------------------
 *						nxString::EnsureLastCharIsDirectoryChar
 *	Ensures that the last character is a directory char.  If the string is empty
 *	then the string is left empty (corresponds to current directory)
 *-------------------------------------------------------------------------*/

void nxString::EnsureLastCharIsDirectoryChar()
{
	if (!IsEmpty())
	{
		if (m_str[m_length-1] != DIRECTORY_CHAR)
		{
			this->operator +=(DIRECTORY_CHAR);
		}
	}
}


/*-----------------------------------------------------------------------------
 *					nxString::EnsureLastCharIsNotDirectoryChar		2009-6-4*/
/** **/
/*---------------------------------------------------------------------------*/

void nxString::EnsureLastCharIsNotDirectoryChar()
{
	char	c;

	if (!IsEmpty())
	{
		c = m_str[m_length-1];
		if (( c == '\\') || (c == '/'))
		{
			m_str[m_length-1] = '\0';
			--m_length;
		}
	}
}
/*---------------------------------------------------------------------------
 *'					nxString::ConvertToWide	                    2003-6-11
 *-------------------------------------------------------------------------*/

nxStringw nxString::ConvertToWide()
{
	nxStringw str;

	if (m_length  > 0) nxStringwFromChar( &str, m_str );
	return str;
}
	

/*-----------------------------------------------------------------------------
 *					nxString::TrimToComments		2007-4-16*/
/** **/
/*---------------------------------------------------------------------------*/

size_t nxString::TrimToComments( const char* commentchars)
{
	size_t	index;

	index = strcspn( m_str, commentchars);
	m_str[index] = '\0';
	m_length     = index;
	RemoveWhiteSpace();
	return m_length;
}


/*-----------------------------------------------------------------------------
 *					nxString::CountNonWhiteFieldsOnLine		2007-4-16*/
/** **/
/*---------------------------------------------------------------------------*/

size_t nxString::CountNonWhiteFieldsOnLine()
{
	size_t			fieldcount = 0;
	bool			inwhite    = true;
	char			c;
	const char*	white = " \t\n\r";

	for ( size_t i = 0; i < m_length; i++)
	{
		c = m_str[i];
		if (inwhite)							/* If we are in white space */							
		{										/* then */
			if ( strchr( white, c) == NULL)		/* if we have not found white space */
			{									/* then */
				inwhite = 0;					/* Flag that we are not in white space */
				fieldcount++;					/* and increment the number of fields */
			}									/* and that */
		}										/* is that */
		else									/* otherwise we are not in white space */
		{											/* so */
			inwhite = ( strchr( white, c) != NULL);	/* if we have found white space */
		}
	}
	return fieldcount;
}


/*-----------------------------------------------------------------------------
 *					nxString::InputTextLine		2007-4-16*/
/** **/
/*---------------------------------------------------------------------------*/

bool nxString::InputTextLine( std::istream& instr)
{
	char	buffer[256];
	bool	incomplete;		

	instr.clear();
	instr.getline( m_str, (std::streamsize)m_allocation );
	UpdateLengthAfterExternalWrite();

	incomplete  = instr.fail();
	while (incomplete && !instr.eof())
	{
		instr.clear();
		instr.getline( buffer, N_ELEMENTS(buffer), '\n' );
		*this += buffer;
		incomplete  = instr.fail();
	}
	return (!incomplete || instr.eof());
}


/*-----------------------------------------------------------------------------
 *					nxString::MakeDirectorySeparatorsOSConsistent		2013-2-28*/
/** Makes all occurrences of directory separarors OS Consistent
**/
/*---------------------------------------------------------------------------*/

bool nxString::MakeDirectorySeparatorsOSConsistent( char oschar)
{
	char			c;

	for ( size_t i = 0; i < m_length; i++)
	{
		c = m_str[i];
		if (( c == '\\') || ( c=='/'))
		{
			m_str[i] = oschar;
		}
	}
	return true;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringArray::Constructor
//-
//--------------------------------------------------------------------------- 

nxStringArray::nxStringArray()
{
	m_array.reserve(64);		// Reserve space for 64 strings (makes things quick)
}


/*-----------------------------------------------------------------------------
 *					nxStringArray::Reserve		2005-3-22*/
/** **/
/*---------------------------------------------------------------------------*/

void nxStringArray::Reserve( int n)
{
	m_array.reserve(n);		// Reserve space for 64 strings (makes things quick)
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringArray::RemoveAll
//-
//---------------------------------------------------------------------------
void nxStringArray::RemoveAll()
{
	m_array.erase( m_array.begin(), m_array.end() );
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringArray::Destructor
//-
//---------------------------------------------------------------------------
nxStringArray::~nxStringArray()
{
	RemoveAll();
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringArray::GetAt
//-
//---------------------------------------------------------------------------

nxString& nxStringArray::GetAt( int nindex)
{
	static nxString blank;
	if ((nindex >= 0) && (nindex < (int)m_array.size())) return (m_array[nindex]);
	nxLog::Record(NXLOG_WARNING, "nxStringArray::GetAt, Index (%d) is out of bounds (0 to %d)", (int)nindex, (int)m_array.size()-1);
	return blank;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringArray::Add
//-
//---------------------------------------------------------------------------

int nxStringArray::Add( const char * str )
{
	m_array.push_back( nxString(str) );
	return (int)m_array.size();											// return the number of points in the list.
}

//---------------------------------------------------------------------------
//			function nxStrtok
//	parses a const string and writes the result to a nxStringArray.
//	The string passed in can either be a nxString or a const char *
//---------------------------------------------------------------------------

// ---- 2016-09-12 ndl303.
// ---- define a macro to wrap strtok for GNU GCC and Vsual Studio. Visual studio implementation of strtok is thread safe but the GNU version is not and they recommend using strtok_r which is thread safe
// ---- THis issue was uncovered when running the pratmo model code on linux and we got registry warning popping up as the multiple threads were corrupting each other.
//

#if defined (_MSC_VER)
#define STRTOK_R( token, sep, saveptr)	strtok( token, sep)
#else
#define STRTOK_R( token, sep, saveptr)	strtok_r( token, sep, saveptr)
#endif

int nxStringArray::Strtok( const char *utstring, const char *separator ) 
{
   char *token;
   char*	saveptr;
   static const char defsep[] = " ,\t\n";
   const char *sep;
   nxString StringToken;
   nxString	aline(utstring);

//   aline = new char [strlen(utstring) + 1];		// copy the string into a temporary buffer
// if (aline == NULL) NXTHROW (("nxStringArray::Strtok, memory allocation error"));
//  strcpy( aline, utstring );					// so we can pass it to strtok without problem
   RemoveAll();									// Remove all the elements from the nxStringArray
   if (separator == NULL) sep = defsep; else sep = separator;

   saveptr = aline.DangerousTypecast();
   token = STRTOK_R( saveptr, sep, &saveptr);					// now get the first token in the list
   while (token != NULL )						// if we got a token
   {											// then
      StringToken = token;						// copy the token to a string
      Add(StringToken);							// and add the string to our array
      token = STRTOK_R( NULL, sep, &saveptr);	// now look fo the next token
   }											// repeat for all of the tokens
   return GetSize();							// return the number of tokens
}



