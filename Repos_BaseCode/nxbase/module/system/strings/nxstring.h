#if !defined(NXBASE_NXSTRING_H)
#define NXBASE_NXSTRING_H 1

#include <vector>
/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/

class nxStringw;

/*----------------------------------------------------------------------------
 *						class nxString										*/
/**	\ingroup system_strings
 *	A class that handles \c char string manipulation.
**/
/*--------------------------------------------------------------------------*/

class  nxString
{
	private:
		char			m_localstorage[100];
		size_t			m_length;
		size_t			m_allocation;
		char*			m_str;


	private:
		nxBOOL			CheckAllocation	( size_t newsize, nxBOOL savestring = nxFALSE );
		nxBOOL			CopyString		( const char* str,  size_t nbytes = 0 );
		void			concatenate		( const char* str1, size_t l1, const char * str2, size_t l2 );
		void 			AddAString		( const char* str,  size_t l1 )  ;
		void			SetToNullStr	();



	public:
						nxString();
						nxString( const char *astring );
						nxString( const nxString& otherstring);
					   ~nxString();
		nxBOOL			ReserveAndErase		( int nchars );
		size_t			UpdateLengthAfterExternalWrite();

		void			Empty				( nxBOOL clearcache);
		int				GetLength			() const	{ return (int)m_length;}
		nxBOOL			IsEmpty				() const	{ return (m_length < 1);}
		char *			DangerousTypecast	()			{ return m_str;}
		void			MakeLower			();
		void			MakeReverse			();
		void			MakeUpper			();
		char			GetAt				( size_t nindex ) const;
		void			SetAt				( size_t nIndex, char ch );
		void			sprintf				( const char* format, ... );
		void			vsprintf			( const char* format, va_list argptr);
		int				Find				( const char* charstr);
		int				Find				( char c );
		int				FindAnyOf			( const char* cstr );
		size_t			TruncateAt			( const char * charset);
		size_t			RemoveWhiteSpace	();
		size_t			TrimToComments		( const char* commentchars);
		size_t			CountNonWhiteFieldsOnLine();
		size_t			Replace				( char oldchar, char newchar);
		nxString		Mid					( size_t nfirst, size_t ncount );
		nxString		Left				( size_t ncount );
		nxString		Right				( size_t ncount );
		const char*		AsFormattedNumber	( int value );
		void			EnsureLastCharIsDirectoryChar		();
		void			EnsureLastCharIsNotDirectoryChar	();
		bool			MakeDirectorySeparatorsOSConsistent	( char oschar = DIRECTORY_CHAR);
		nxStringw		ConvertToWide();
		bool			InputTextLine		( std::istream& instr);




		nxString&		operator=( const nxString &otherstring);
		nxString&		operator=( const char*     psz);
		nxString&		operator=( const wchar_t*    psz);
		nxString&		operator=( char ch );
						operator const char*() const {return m_str;}
		nxString&		operator += ( const nxString& str );
		nxString&		operator += ( char ch );
		nxString&		operator += ( const char* str );
		char			operator [] (int nindex) const;

friend	 nxString operator + (const nxString& str1, const nxString & str2 );
friend	 nxString operator + (const nxString& str1, char ch );
friend	 nxString operator + (char ch, const nxString& str1 );
friend	 nxString operator + ( const nxString& str1,  const char* str2 );
friend	 nxString operator + (const char* str1, const nxString& str2 );

};

extern  nxBOOL	operator == ( const nxString& s1, const nxString& str2 );
extern  nxBOOL	operator == ( const nxString& s1, const char *    str2 );
extern  nxBOOL	operator == ( const char *    s1, const nxString& str2 );
extern  nxBOOL	operator != ( const nxString& s1, const nxString& str2 );
extern  nxBOOL	operator != ( const nxString& s1, const char *    str2 );
extern  nxBOOL	operator != ( const char *    s1, const nxString& str2 );
extern  nxBOOL	operator <  ( const nxString& s1, const nxString& str2 );
extern  nxBOOL	operator <  ( const nxString& s1, const char *    str2 );
extern  nxBOOL	operator <  ( const char *    s1, const nxString& str2 );
extern  nxBOOL	operator >  ( const nxString& s1, const nxString& str2 );
extern  nxBOOL	operator >  ( const nxString& s1, const char *    str2 );
extern  nxBOOL	operator >  ( const char *    s1, const nxString& str2 );
extern  nxBOOL	operator >= ( const nxString& s1, const nxString& str2 );
extern  nxBOOL	operator >= ( const nxString& s1, const char *    str2 );
extern  nxBOOL	operator >= ( const char *    s1, const nxString& str2 );
extern  nxBOOL	operator <= ( const nxString& s1, const nxString& str2 );
extern  nxBOOL	operator <= ( const nxString& s1, const char *    str2 );
extern  nxBOOL	operator <= ( const char *    s1, const nxString& str2 );
extern  nxString operator + (const nxString& str1, const nxString & str2 );
extern  nxString operator + (const nxString& str1, char ch );
extern  nxString operator + (char ch, const nxString& str1 );
extern  nxString operator + ( const nxString& str1,  const char* str2 );
extern  nxString operator + (const char* str1, const nxString& str2 );

/*---------------------------------------------------------------------------
 *						class nxStringw										*/
/**	\ingroup system_strings
 *	A class that handles Unicode string manipulation.
**/
/*--------------------------------------------------------------------------*/

class  nxStringw
{
	private:
		wchar_t			m_internalstorage[100];
		size_t			m_length;
		size_t			m_allocation;
		wchar_t*		m_str;


	private:
		nxBOOL			CheckAllocation( size_t newsize, nxBOOL savestring = nxFALSE );
		nxBOOL			CopyString  ( const wchar_t *str, size_t nbytes = 0 );
		void			concatenate ( const wchar_t* str1, size_t l1, const wchar_t * str2, size_t l2 );
		void 			AddAString  ( const wchar_t* str, size_t l1 )  ;
		void			SetToNullStr();



	public:
						nxStringw();
						nxStringw( const wchar_t *astring );
						nxStringw( const nxStringw& otherstring);
					   ~nxStringw();
		void			Empty( nxBOOL clearcache);
		size_t			GetLength() const { return m_length;}
		nxBOOL			IsEmpty()   const { return m_length < 1;}
		void			MakeLower();
		void			MakeReverse();
		void			MakeUpper();
		wchar_t			GetAt( size_t nindex ) const;
		void			SetAt( size_t nIndex, wchar_t ch );
		void			sprintf( const wchar_t* format, ... );
		void			vsprintf( const wchar_t* format, va_list argptr);
		int				Find( const wchar_t* charstr);
		int				Find( wchar_t c );
		size_t			TruncateAt( const wchar_t * charset);
		size_t			RemoveWhiteSpace();
		int				Replace( wchar_t oldchar, wchar_t newchar);
		nxStringw		Mid   (size_t nfirst, size_t ncount );
		nxStringw		Left  (size_t ncount );
		nxStringw		Right (size_t ncount );
		wchar_t *			DangerousTypecast() { return m_str;}
		const wchar_t*	AsFormattedNumber( int value );
		void			EnsureLastCharIsDirectoryChar();
		nxString		ConvertToChar();



		const nxStringw&	operator=( const nxStringw &otherstring);
		const nxStringw&	operator=( const wchar_t*     psz);
		const nxStringw&	operator=( const char*      psz);
		const nxStringw&	operator=( wchar_t ch );
							operator const wchar_t*() const {return m_str;}
		const nxStringw&	operator += ( const nxStringw& str );
		const nxStringw&	operator += ( wchar_t ch );
		const nxStringw&	operator += ( const wchar_t* str );
		wchar_t				operator [] (size_t nindex) const;

friend	 nxStringw operator + (const nxStringw& str1, const nxStringw & str2 );
friend	 nxStringw operator + (const nxStringw& str1, wchar_t ch );
friend	 nxStringw operator + (wchar_t ch, const nxStringw& str1 );
friend	 nxStringw operator + ( const nxStringw& str1,  const wchar_t* str2 );
friend	 nxStringw operator + (const wchar_t* str1, const nxStringw& str2 );

};

extern  nxBOOL	  operator == ( const nxStringw& s1, const nxStringw& str2 );
extern  nxBOOL	  operator == ( const nxStringw& s1, const wchar_t *    str2 );
extern  nxBOOL	  operator == ( const wchar_t *    s1, const nxStringw& str2 );
extern  nxBOOL	  operator != ( const nxStringw& s1, const nxStringw& str2 );
extern  nxBOOL	  operator != ( const nxStringw& s1, const wchar_t *    str2 );
extern  nxBOOL	  operator != ( const wchar_t *    s1, const nxStringw& str2 );
extern  nxBOOL	  operator <  ( const nxStringw& s1, const nxStringw& str2 );
extern  nxBOOL	  operator <  ( const nxStringw& s1, const wchar_t *    str2 );
extern  nxBOOL	  operator <  ( const wchar_t *    s1, const nxStringw& str2 );
extern  nxBOOL	  operator >  ( const nxStringw& s1, const nxStringw& str2 );
extern  nxBOOL	  operator >  ( const nxStringw& s1, const wchar_t *    str2 );
extern  nxBOOL	  operator >  ( const wchar_t *    s1, const nxStringw& str2 );
extern  nxBOOL	  operator >= ( const nxStringw& s1, const nxStringw& str2 );
extern  nxBOOL	  operator >= ( const nxStringw& s1, const wchar_t *    str2 );
extern  nxBOOL	  operator >= ( const wchar_t *    s1, const nxStringw& str2 );
extern  nxBOOL	  operator <= ( const nxStringw& s1, const nxStringw& str2 );
extern  nxBOOL	  operator <= ( const nxStringw& s1, const wchar_t *    str2 );
extern  nxBOOL	  operator <= ( const wchar_t *    s1, const nxStringw& str2 );
extern  nxStringw operator +  ( const nxStringw& str1, const nxStringw & str2 );
extern  nxStringw operator +  ( const nxStringw& str1, wchar_t ch );
extern  nxStringw operator +  ( wchar_t ch, const nxStringw& str1 );
extern  nxStringw operator +  ( const nxStringw& str1,  const wchar_t* str2 );
extern  nxStringw operator +  ( const wchar_t* str1, const nxStringw& str2 );


/*----------------------------------------------------------------------------
 *						class nxStringArray									*/
/**	\ingroup system_strings
 * A class that supports an array of strings.  Frequently used for Strtok
 * parsing and for directory listings.
**/
/*---------------------------------------------------------------------------*/

class  nxStringArray
{
	private:
		std::vector<nxString>	m_array;		// pointer to an array of pointers

	public:
						nxStringArray();
					   ~nxStringArray();
		void			Reserve( int n);

		int				GetSize() const {return (int)m_array.size();}
		int				Add( const char *	str );
		nxString&		GetAt( int nindex);
		void			RemoveAll();
		int				Strtok( const char *utstring, const char *separator = NULL  );
		nxString&		operator [] (int nindex){return GetAt(nindex);}

};



/*--------------------------------------------------------------------------
 *					nxStringFromUnicode										*/
/**	\ingroup system_strings
 * Converts a Unicode string to an nxString object. Only designed where the UNICODE string
 * can be properly represented as an ASCII string.
 *
 *	\param str
 *		pointer to nxString object that will contain the converted unicode string	
 *
 *	\param unicodestr
 *		The null terminated Unicode string.
**/
/*------------------------------------------------------------------------*/

extern void nxStringFromUnicode( nxString*  str, const wchar_t* unicodestr  );

/*--------------------------------------------------------------------------
 *					nxStringwFromChar										*/
/**	\ingroup system_strings
 * Converts a char string to a Unicode nxStringw object.
 *
 *	\param str
 *		pointer to nxStringw object that will contain the converted unicode string	
 *
 *	\param bytestr
 *		The null terminated char string.
**/
/*------------------------------------------------------------------------*/

extern void nxStringwFromChar  ( nxStringw* str, const char*  bytestr     );


#if defined(NX_WINDOWS)
/*--------------------------------------------------------------------------
 *					nxStringToBSTR											*/
/**	\ingroup system_strings
 * Converts an nxString to a Unicode, BSTR used in IDispatch automation
 *	
 *	\param str
 *		The string to be converted.
 *			
 *`	\param bstr
 *		Returns the BSTR representation
 **/
/*------------------------------------------------------------------------*/
void nxStringToBSTR( nxString& str, BSTR* bstr );
#endif

#include "nxfuncs.h"

#endif
