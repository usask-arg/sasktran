
/*--------------------------------------------------------------------------
 *					class nxBase36                                2002-9-30*/
/**	\ingroup system_strings
 *	A class for encoding numbers using base 36. ie the uppercase alpha
 *	numeric characters.
**/
/*------------------------------------------------------------------------*/

class  nxBase36
{
	private:
		int					m_value;

	private:
		int					NumDigits() const;

	public:
		static int			FromChar( char c );
		static char			AsChar  ( int value );

	public:
							nxBase36()						{ m_value = 0;}
							nxBase36(int val)				{ m_value = val;}
		int					DecodeFromStr( const char *str);
		nxString			EncodeAsStr( int mindigits ) const;
							operator int () const			{ return m_value;}
							operator nxString () const		{ return EncodeAsStr(0);}
		int					operator = (int val)			{ m_value = val; return m_value;}
		int					operator = (const char * str)	{ return DecodeFromStr(str);}
};


