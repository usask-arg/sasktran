/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_core.h"
#include <wchar.h>
#include <stdarg.h>

static wchar_t NULLSTRBUFFER[]  = {0,0,0};
static wchar_t*  const NULLSTR  = &NULLSTRBUFFER[1];

//---------------------------------------------------------------------------
//						nxStringFromUnicode
//---------------------------------------------------------------------------

void nxStringFromUnicode( nxString* str, const wchar_t* unicodestr  )
{
	char	quick[1024];
	char*	buffer;
	nxBOOL	needsalloc;
	size_t	n;

#if defined(NX_WINDOWS)
	n = ::WideCharToMultiByte( CP_ACP, 0, unicodestr, -1, NULL, 0, NULL, NULL );		// Get the length of buffer required
#else
	n = wcslen( unicodestr );
#endif

	needsalloc = (n >= N_ELEMENTS(quick));												// Can we use a stack buffer
	if (!needsalloc) buffer = &quick[0];												// Yup so do it
	else             buffer = new char[n+1];											// otherwise allocate our buffer
	NXASSERT(buffer != NULL);															// make sure we dont have NULL

#if defined(NX_WINDOWS)
	n = ::WideCharToMultiByte( CP_ACP, 0, unicodestr, -1, buffer, (int)(n+1), NULL, NULL );	// do the conversion
#else
	for (size_t i=0; i < n; i++)
	{
		buffer[i] = (char)unicodestr[i];
	};
	buffer[n] = '\0';
#endif

	(*str) = (char *)buffer;
	if (needsalloc) delete [] buffer;
}

/*---------------------------------------------------------------------------
 *'					nxStringwFromChar                               2003-6-11
 *-------------------------------------------------------------------------*/

void nxStringwFromChar( nxStringw* str, const char* bytestr  )
{
	wchar_t	quick[512];
	wchar_t*	buffer;
	nxBOOL	needsalloc;
	size_t	n;

#if defined(NX_WINDOWS)
	n = ::MultiByteToWideChar( CP_ACP, 0, bytestr, -1, NULL, 0 );		// Get the length of buffer required
#else
	n = strlen( bytestr);
#endif
	needsalloc = (n >= N_ELEMENTS(quick));													// Can we use a stack buffer
	if (!needsalloc) buffer = &quick[0];												// Yup so do it
	else             buffer = new wchar_t[n+1];												// otherwise allocate our buffer
	NXASSERT(buffer != NULL);															// make sure we dont have NULL
#if defined(NX_WINDOWS)
	n = ::MultiByteToWideChar( CP_ACP, 0, bytestr, -1, buffer, (int)(n+1));	// do the conversion
#else
	for (size_t i=0; i < n; i++)
	{
		buffer[i] = bytestr[i];
	};
	buffer[n] = L'\0';
#endif
	(*str) = buffer;
	if (needsalloc) delete [] buffer;
}

static wchar_t *nxstrlwr( wchar_t *string )
{
	wchar_t *ptr = string;
	if (ptr == NULL) return NULL;
	while (*ptr != L'\0')
	{
		*ptr = towlower( *ptr );
		ptr++;
	}
	return string;
}

static wchar_t *nxstrupr( wchar_t *string )
{
	wchar_t *ptr = string;
	if (ptr == NULL) return NULL;
	while (*ptr != L'\0')
	{
		*ptr = towupper( *ptr );
		ptr++;
	}
	return string;
}

static wchar_t *nxstrrev( wchar_t *string )
{
	if (string == NULL) return NULL;
	size_t n = wcslen(string);
	wchar_t * source = string + n-1;
	wchar_t * ptr    = string;
	wchar_t c;

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
//						nxStringw::Default blank constructor
//	Allocate a blank NULL string.
//-
//---------------------------------------------------------------------------

nxStringw::nxStringw()
{
	m_internalstorage[0] = L'\0';
	m_length     = 0;
	m_str        = &m_internalstorage[0];
	m_allocation = N_ELEMENTS(m_internalstorage);
}

nxStringw::nxStringw( const wchar_t *str )
{
	m_internalstorage[0] = L'\0';
	m_str                = &m_internalstorage[0];
	m_length             = 0;
	m_allocation         = N_ELEMENTS(m_internalstorage);
	*this = str;
}
 
nxStringw::nxStringw( const nxStringw& otherstring )
{
	m_internalstorage[0] = L'\0';
	m_str                = &m_internalstorage[0];
	m_length             = 0;
	m_allocation         = N_ELEMENTS(m_internalstorage);
	*this                = otherstring;
}
//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::Destructor
//	Ensure that the memory is properly deallocated.
//-
//---------------------------------------------------------------------------

nxStringw::~nxStringw()
{
	Empty(nxTRUE);
}

void nxStringw::SetToNullStr()
{
	if (m_str != NULL) m_str[0] = L'\0';
	m_length = 0;
}

//---------------------------------------------------------------------------
//						nxStringw::Empty						
//---------------------------------------------------------------------------

void nxStringw::Empty( nxBOOL clearcache )
{
	if (clearcache)
	{
		if ((m_str != &m_internalstorage[0]) && (m_str != NULL))
		{
			delete [] m_str;				
		}
		m_str        = &m_internalstorage[0];
		m_allocation = N_ELEMENTS(m_internalstorage);
	}
	SetToNullStr();
}

//---------------------------------------------------------------------------
//						nxStringw::CheckAllocation
//	Checks the allocation of the memory buffer.  If the user requests that
//	the current string be saved then it is done.
//---------------------------------------------------------------------------

nxBOOL nxStringw::CheckAllocation( size_t nbytes, nxBOOL savestring )
{
	nxStringw	temp;
	nxBOOL		ok;
	size_t		oldallocation;

	ok = (nbytes <= m_allocation);
	if (!ok) 											// if the requested size is greater than this allocation
	{													// then
		oldallocation = m_allocation + 100;
		if (savestring) temp = *this;						// if the user requests the string be saved then make a temp copy
		Empty(nxTRUE);										// clear the current string
		if (nbytes < oldallocation) nbytes = oldallocation;					// may as well allocate decent size string (enough to hold a filename perhaps?)
		m_str = new wchar_t [nbytes];					// and create the new segment
		ok    = (m_str != NULL);					// check that the string was ok.
		if (ok)										// if it was						
		{											// then
			m_allocation = nbytes;					// allocate the string
			if (savestring) *this = temp;			// and copy back if saving of source required.
		}
		else
		{
			m_str        = &m_internalstorage[0];
			m_allocation = N_ELEMENTS(m_internalstorage);
			SetToNullStr();
			NXTHROW(("nxStringw memory allocation error"));
		}
	}
	return ok;
}
//---------------------------------------------------------------------------
//						nxStringw::CopyString
//---------------------------------------------------------------------------

nxBOOL nxStringw::CopyString( const wchar_t *str, size_t nbytes )
{
	nxBOOL	ok;

	if (str == NULL)
	{
		SetToNullStr();
		ok       = nxFALSE;
	}
	else
	{
		if (nbytes == 0) nbytes = wcslen(str) + 1;
		ok = CheckAllocation( nbytes );
		if (ok)
		{
			wcscpy( m_str, str );
		}
		m_length = nbytes-1;
	}
	return ok;
}
//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::operator=
//	Copies a nxStringw from one string to another.
//-
//---------------------------------------------------------------------------
const nxStringw&	nxStringw::operator=( const nxStringw &otherstring)
{
	size_t		nbytes;

	nbytes = otherstring.GetLength() + 1;
	CopyString( otherstring, nbytes );
	return *this;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::operator=
//	Copies a wchar_t * string to a nxStringw
//-
//---------------------------------------------------------------------------
const nxStringw&	nxStringw::operator=( const wchar_t * otherstring)
{
	CopyString( otherstring );
	return *this;
}

/*---------------------------------------------------------------------------
 *'					nxStringw::operator=                                       2003-6-11
 *-------------------------------------------------------------------------*/

const nxStringw& nxStringw::operator=( const char*      psz)
{
	nxStringwFromChar( this, psz);
	return *this;
}


//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::operator=
//	Copies a wchar_t * string to a nxStringw
//-
//---------------------------------------------------------------------------

const nxStringw&	nxStringw::operator=( wchar_t  achar)
{
	wchar_t dummy[2];

	dummy[0] = achar;
	dummy[1] = L'\0';
	CopyString( dummy, 2 );
	return *this;
}
//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::MakeLower
//	Makes the string lower case.
//-
//---------------------------------------------------------------------------

void nxStringw::MakeLower()
{
	if (m_length > 0) nxstrlwr( m_str );
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::MakeReverse
//	Reverses the order of this string
//-
//---------------------------------------------------------------------------

void nxStringw::MakeReverse()
{
	if (m_length > 1)  nxstrrev( m_str );
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::MakeUpper
//	converts the entire string to uppercase
//-
//---------------------------------------------------------------------------

void nxStringw::MakeUpper()
{
	if (m_length > 0) nxstrupr( m_str );
}


//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::GetAt
//	Returns the character at a specific location in the string
//-
//---------------------------------------------------------------------------

wchar_t nxStringw::GetAt( size_t nindex ) const
{
	wchar_t c;

	if ((nindex >= 0) && ( nindex < m_length )) c = m_str[nindex]; else c = L'\0';
	return c;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::operator []
//	returns the characters at location indexed by array operator
//-
//---------------------------------------------------------------------------

wchar_t nxStringw::operator [] ( size_t nindex) const
{
	return GetAt(nindex);
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::SetAt
//	Sets the string at a specific location to the character "ch"
//	if ch is the null terminator then the length of the string
//	is adjusted accordingly.
//-
//---------------------------------------------------------------------------

void nxStringw::SetAt( size_t nindex, wchar_t ch )
{
	if ((nindex >= 0) && ( nindex < m_length)) m_str[nindex] = ch;
	if (ch == L'\0') m_length = nindex;
}

//---------------------------------------------------------------------------
//						nxStringw::concatenate
//		make the string by concatenating two string pointers.
//---------------------------------------------------------------------------

void nxStringw::concatenate( const wchar_t * str1, size_t l1, const wchar_t *str2, size_t l2 )
{
	nxBOOL	ok;
	wchar_t*	outstr;

	ok = CheckAllocation( l1 + l2 + 1 );						// get the total allocation required
	if (ok)														// if it worked
	{															// then
		wcscpy( m_str, str1 );									// copy the first string
		outstr = m_str + l1;									// get the start of the second string
		wcscpy( outstr, str2 );									// and add the second string
		m_length = l1 + l2;										// and set the length of the string
	}
}		

//---------------------------------------------------------------------------
//+
//NAME:
//						CString operator + 
//-
//---------------------------------------------------------------------------

 nxStringw operator + (const nxStringw& str1, const nxStringw & str2 )
{
	nxStringw s;
	s.concatenate( str1, str1.GetLength(), str2, str2.GetLength());
	return s;
}

 nxStringw operator + (const nxStringw& str1, wchar_t ch )
{
	nxStringw s;
	wchar_t	 dummy[2];

	dummy[1] = L'\0';
	dummy[0] = ch;
	s.concatenate( str1, str1.GetLength(), dummy, wcslen(dummy));
	return s;
}

 nxStringw operator + (wchar_t ch, const nxStringw& str1 )
{
	nxStringw	s;
	wchar_t		dummy[2];

	dummy[1] = L'\0';
	dummy[0] = ch;
	s.concatenate( dummy, wcslen(dummy), str1, str1.GetLength());
	return s;
}

 nxStringw operator + ( const nxStringw& str1,  const wchar_t* str2 )
{
	nxStringw s;

	if (str2 == NULL) str2 = NULLSTR;
	s.concatenate( str1, str1.GetLength(), str2, wcslen(str2) );
	return s;
}

 nxStringw operator + (const wchar_t* str1, const nxStringw& str2 )
{
	nxStringw s;

	if ( str1 == NULL) str1 = NULLSTR;
	s.concatenate( str1, wcslen(str1), str2, str2.GetLength());
	return s;
}

//---------------------------------------------------------------------------
//					nxStringw::AddAString
//---------------------------------------------------------------------------

void nxStringw::AddAString( const wchar_t* str, size_t l1 )  
{
	nxBOOL	ok;

	ok = CheckAllocation( l1 + GetLength() + 1, nxTRUE );			// get the total allocation required amd save the orginal string
	if (ok)														// if it worked
	{															// then
		wcscat( m_str, str );									// copy the first string
		m_length = l1 + m_length;								// and set the length of the string
	}
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::operator +=
//-
//---------------------------------------------------------------------------

const nxStringw & nxStringw::operator += ( const nxStringw& str )
{
	AddAString( str, str.GetLength() );
	return *this;
}

const nxStringw & nxStringw::operator += ( wchar_t ch )
{
	wchar_t dummy[2];

	dummy[0] = ch;
	dummy[1] = L'\0';
	AddAString( dummy, wcslen(dummy) );
	return *this;
}

const nxStringw & nxStringw::operator += ( const wchar_t* str )
{
	if (str != NULL) AddAString( str, wcslen(str));
	return *this;
}


static int safestrcmp( const wchar_t* s1, const wchar_t * s2 )
{
	if (s1 == NULL) s1 = NULLSTR;
	if (s2 == NULL) s2 = NULLSTR;
	return wcscmp(s1,s2);
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::operator ==
//-
//---------------------------------------------------------------------------
 nxBOOL operator == ( const nxStringw& s1, const nxStringw& s2 )
{
	return safestrcmp( s1, s2 ) == 0;
}
		
 nxBOOL operator == ( const nxStringw& s1, const wchar_t *  s2 )
{
	return safestrcmp( s1, s2 ) == 0;
}

 nxBOOL operator == ( const wchar_t * s1, const nxStringw& s2 )
{
	return safestrcmp( s1, s2 ) == 0;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::operator !=
//-
//---------------------------------------------------------------------------
 nxBOOL operator != ( const nxStringw& s1, const nxStringw& s2 )
{
	return safestrcmp( s1, s2 ) != 0;
}
		
 nxBOOL operator != ( const nxStringw& s1, const wchar_t *  s2 )
{
	return safestrcmp( s1, s2 ) != 0;
}

 nxBOOL operator != ( const wchar_t * s1, const nxStringw& s2 )
{
	return safestrcmp( s1, s2 ) != 0;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::operator <
//-
//---------------------------------------------------------------------------

 nxBOOL operator < ( const nxStringw& s1, const nxStringw& s2 )
{
	return safestrcmp( s1, s2 ) < 0;
}
		
 nxBOOL operator < ( const nxStringw& s1, const wchar_t *  s2 )
{
	return safestrcmp( s1, s2 ) < 0;
}

 nxBOOL operator < ( const wchar_t * s1, const nxStringw& s2 )
{
	return safestrcmp( s1, s2 ) < 0;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::operator >
//-
//---------------------------------------------------------------------------

 nxBOOL operator > ( const nxStringw& s1, const nxStringw& s2 )
{
	return safestrcmp( s1, s2 ) > 0;
}
		
 nxBOOL operator > ( const nxStringw& s1, const wchar_t *  s2 )
{
	return safestrcmp( s1, s2 ) > 0;
}

 nxBOOL operator > ( const wchar_t * s1, const nxStringw& s2 )
{
	return safestrcmp( s1, s2 ) > 0;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::operator <=
//-
//---------------------------------------------------------------------------

 nxBOOL operator <= ( const nxStringw& s1, const nxStringw& s2 )
{
	return safestrcmp( s1, s2 ) <= 0;
}
		
 nxBOOL operator <= ( const nxStringw& s1, const wchar_t *  s2 )
{
	return safestrcmp( s1, s2 ) <= 0;
}

 nxBOOL operator <= ( const wchar_t * s1, const nxStringw& s2 )
{
	return safestrcmp( s1, s2 ) <= 0;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::operator >=
//-
//---------------------------------------------------------------------------

 nxBOOL operator >= ( const nxStringw& s1, const nxStringw& s2 )
{
	return safestrcmp( s1, s2 ) >= 0;
}
		
 nxBOOL operator >= ( const nxStringw& s1, const wchar_t *  s2 )
{
	return safestrcmp( s1, s2 ) >= 0;
}

 nxBOOL operator >= ( const wchar_t * s1, const nxStringw& s2 )
{
	return safestrcmp( s1, s2 ) >= 0;
}


//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::sprintf
//	Emulates the sprintf function, but allocates a large (1024 byte)
//	buffer to ensure that memory does not overflow.
//-
//---------------------------------------------------------------------------

void nxStringw::sprintf( const wchar_t *format, ... )
{
   va_list	   ArgPtr;
   wchar_t	   buffer[1024];

   va_start(ArgPtr, format);
#if defined(NX_WINDOWS)
	_vsnwprintf( buffer, N_ELEMENTS(buffer), format, ArgPtr );
#else
    vswprintf(buffer, N_ELEMENTS(buffer), format, ArgPtr );
#endif
   va_end(ArgPtr);
   *this = buffer;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::vsprintf
//	Emulates the vsprintf function, but allocates a large (1024 byte)
//	buffer to ensure that memory does not overflow
//-
//---------------------------------------------------------------------------

void nxStringw::vsprintf( const wchar_t* format, va_list argptr)
{
	wchar_t buffer[1024];

#if defined(NX_WINDOWS)
	_vsnwprintf( buffer, N_ELEMENTS(buffer), format, argptr );
#else
    vswprintf(buffer, N_ELEMENTS(buffer), format, argptr );
#endif
	*this = buffer;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::Find
//	Find the first occurence of string charstr or character c in this
//	nxStringw object.  Returns the index of the first character that
//	mathes.  Returns -1 if no macth is found.
//-
//---------------------------------------------------------------------------

int nxStringw::Find( wchar_t c )
{
	wchar_t	*ptr;
	int		idx;

 	ptr = wcschr( m_str, (int)c); 
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

int nxStringw::Find( const wchar_t* charstr)
{
	wchar_t	*ptr;
	int		idx;
 	ptr = wcsstr( m_str, charstr); 
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
//						nxStringw::TruncateAt
//	Truncates the string at the first occurence of any character
//	in string charset.  Uses the strcspn function.  Very useful
//	for eliminating comments from the ends of lines etc.
//	Returns the length of the string
//-
//---------------------------------------------------------------------------

size_t nxStringw::TruncateAt( const wchar_t* charset)
{
	size_t idx;
	
	if (m_length > 0)
	{
		idx = wcscspn( m_str, charset );
		m_str[idx] = L'\0';
		m_length   = idx;
	}
	return m_length;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::RemoveWhiteSpace
//	Removes leading and trailing white spaces from the text.
//	White space is any non-printable wchar_t (<= SPACE wchar_t)
//	returns the length of the new string
//-
//---------------------------------------------------------------------------
size_t nxStringw::RemoveWhiteSpace()
{
	size_t i;
	wchar_t*	newstr;
	wchar_t*	outstr;
	
	if (m_length > 0)
	{
		for (i = m_length-1;  i >= 0; i--)			// scan backwards through the file
		{											// until we find
			if (m_str[i] > L' ')					// a non-printable wchar_t
			{										// so now
				m_str[i+1] = L'\0';					// terminate the string at this point
				break;								// and break out of the loop
			}										// otherwise loop until a printable wchar_t
		}											// is found or we hit the end 
	}
	newstr = m_str;								// get the start of the string
	while (*newstr != L'\0')					// while we have more characters
	{                            				// then
		if (*newstr > L' ') break;				// hunt for non-printables
		newstr++;								// increment through loop until found;
	}

	outstr = m_str;								// now copy from the original string
	while (*newstr != L'\0')					// to the beginning of the original
	{											// I have avoided strcpy as I'm
		*outstr++ = *newstr++;					// not sure how it handles overlapped			
	}											// strings on all platforms
	*outstr = L'\0';							// terminate the copied string
	m_length = (size_t)(outstr - m_str);		// and reset the length of this string.
	return m_length;
}												// and that is that.

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::Replace
//	replace all occurences of character "oldchar" with character "newchar"
//	If newchar is the null terminator then the string length is readjusted.
//-
//--------------------------------------------------------------------------- 

int	nxStringw::Replace( wchar_t oldchar, wchar_t newchar )
{
	wchar_t *str;
	int	  i = 0;
	
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
	if (newchar == '\0') m_length = wcslen(m_str);
	return i;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::Mid
//	Similar to the BASIC MID$ function returns a sub string from
//	another string.  The function starts at character offset nfirst and
//	generates a new string upto ncount bytes long and returns the value
//	as a new nxStringw
//-
//---------------------------------------------------------------------------

nxStringw nxStringw::Mid( size_t nfirst, size_t ncount )
{
	wchar_t*	instr;
	wchar_t*	outstr;
	wchar_t*	outs;
	size_t	i;
	nxStringw result;
	
	if (m_length == 0 || ncount < 1 || nfirst >= m_length) return L"";
	
	outstr = new wchar_t[ncount+1];
	if (outstr == NULL) NXTHROW(("nxStringw::Mid memory allocation error"));
	instr  = m_str + nfirst;
	i      = 0;
	outs   = outstr;
	
	while ((*instr != L'\0') && ( i < ncount))
	{
		*outs++ = *instr++;
		i++;
	}
	*outs = L'\0';
	result = outstr;
	delete [] outstr;
	return result;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::Left
//	Similar to the BASIC Left function returns a sub string from
//	another string.  The function returns the left most characters upt ncount
//-
//---------------------------------------------------------------------------

nxStringw nxStringw::Left( size_t ncount )
{
	wchar_t*	instr;
	wchar_t*	outstr;
	wchar_t*	outs;
	size_t		i;
	nxStringw result;
	
	if (m_length == 0 || ncount < 1) return L"";
	if (ncount > m_length) ncount = m_length;
	outstr = new wchar_t[ncount+1];
	if (outstr == NULL) NXTHROW(("nxStringw::Left memory allocation error"));
	outs   = outstr;
	instr  = m_str;
	for (i = 0; i < ncount; i++) *outs++ = *instr++;
	*outs = L'\0';
	result = outstr;
	delete [] outstr;
	return result;
}

//---------------------------------------------------------------------------
//+
//NAME:
//						nxStringw::Right
//	Similar to the BASIC Right function returns a sub string from
//	another string.  The function returns the right most characters
//	upto ncount
//-
//---------------------------------------------------------------------------

nxStringw nxStringw::Right( size_t ncount )
{
	wchar_t*	instr;
	wchar_t*	outstr;
	wchar_t*   outs;
	size_t		i;
	nxStringw result;
	
	if (m_length == 0 || ncount < 1) return L"";
	if (ncount > m_length) ncount = m_length;
	outstr = new wchar_t[ncount+1];
	if (outstr == NULL) NXTHROW(("nxStringw::Right memory allocation error"));
	outs   = outstr;
	instr  = m_str + (m_length - ncount);
	for (i = 0; i < ncount; i++) *outs++ = *instr++;
	*outs = L'\0';
	result = outstr;
	delete [] outstr;
	return result;
}

/*

const wchar_t * nxStringw::AsFormattedNumber( int value )
{
#if defined(NX_WINDOWS)
	
	wchar_t      chrbuf[400];
	nxStringw  str;

	NUMBERFMT	fmt;
	fmt.NumDigits = 0;
	fmt.LeadingZero = 0;
	fmt.Grouping    = 3;
	fmt.lpDecimalSep = ".";
	fmt.lpThousandSep = ",";
	fmt.NegativeOrder = 0; 
	str.sprintf(L"%d",(int)value );
	::GetNumberFormat( LOCALE_USER_DEFAULT, 0, (const wchar_t *)str, &fmt, chrbuf, sizeof(chrbuf)-1 );  
	*this = chrbuf;
#else
	this->sprintf("%d",(int)value);
#endif
	return  (const wchar_t *)(*this);
}
*/

/*---------------------------------------------------------------------------
 *						nxStringw::EnsureLastCharIsDirectoryChar
 *	Ensures that the last character is a directory wchar_t.  If the string is empty
 *	then the string is left empty (corresponds to current directory)
 *-------------------------------------------------------------------------*/

void nxStringw::EnsureLastCharIsDirectoryChar()
{
	if (!IsEmpty())
	{
		if (m_str[m_length-1] != DIRECTORY_CHAR)
		{
			this->operator +=(DIRECTORY_CHAR);
		}
	}
}

/*---------------------------------------------------------------------------
 *'					nxStringw::ConvertToChar	                    2003-6-11
 *-------------------------------------------------------------------------*/

nxString nxStringw::ConvertToChar()
{
	nxString str;

	if (m_length  > 0) nxStringFromUnicode( &str, m_str );
	return str;
}





