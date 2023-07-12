#include "nxbase_math.h"

//---------------------------------------------------------------------------
//						nxBase36::FromChar
//---------------------------------------------------------------------------

int nxBase36::FromChar( char c )
{
	int value = -1;

	if ( c >= '0' && c <= '9')
	{
		value = c - '0';
	}
	else
	{
		if ( c >= 'A' && c <= 'Z')
		{
			value = c-'A'+10;
		}
		else
		{
			if (c >= 'a' && c <= 'z' ) value = c - 'a' + 10;
		}
	}
	return value;
}

//---------------------------------------------------------------------------
//						nxBase36::AsChar
//---------------------------------------------------------------------------

char nxBase36::AsChar( int value )
{
	if ((value >= 36) || (value < 0)) return '#';
	if (value < 10) return (char)('0' + value);
	return (char)('a' + value-10);
}

//--------------------------------------------------------------------------
//					
int nxBase36::NumDigits() const
{
	int ndigits;
	if (m_value < 36)
	{
		ndigits = 1;
	}
	else
	{
		ndigits = (int)( log((double)m_value)/log(36.0) ) + 1;
	}
	return ndigits;
}

//---------------------------------------------------------------------------
//						EncodeAsBase36
//---------------------------------------------------------------------------

nxString nxBase36::EncodeAsStr( int mindigits ) const
{
	nxString	answer;
	if (m_value < 0) return answer;

	int			i;
	int         thisvalue;
	int			totaldigits;
	int			numdigits = NumDigits();							// get the number of digits to encode
	char *		str;
	int			value = m_value;
	int			power = 36;
	
	if (numdigits > mindigits) totaldigits = numdigits;
	else                       totaldigits = mindigits;
	str = new char [totaldigits+2];
	if (str == NULL) NXTHROW( ("nxBase36::EncodeAsStr, Memory allocation error"));

	str[totaldigits] = '\0';
	for (i=totaldigits-1; i >= 0; i-- )
	{
		thisvalue = value%power;
		str[i] = AsChar(thisvalue);
		value -= thisvalue;
		power *= 36;
	}
	answer = str;
	delete [] str;
	return answer;
}




int nxBase36::DecodeFromStr( const char * astr )
{
	nxString	str(astr);
	int value = 0;
	int multiplier = 1;
	int digit = 0;

	for (int i =0; i < str.GetLength(); i++ )
	{
		digit = FromChar( str[i] );
		if (digit < 0)break;
		value = value*multiplier + digit;
		if (multiplier == 1) multiplier = 36;
	}
	return value;
}


