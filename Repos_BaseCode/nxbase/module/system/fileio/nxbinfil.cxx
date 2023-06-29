/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
#include "nxbase_core.h"

nxBinaryFile& nxBinaryFile::operator<<( const char * str )
{
	int nbytes;

	if (str != NULL)
	{
		nbytes = (int)strlen(str);
		*this << nbytes;
		if (nbytes > 0)
		{
			m_ok = m_ok && (Write( str, sizeof(char), nbytes) == (size_t)nbytes);
		}
	}
	else
	{
		nbytes = -1;
		*this << nbytes;
	}
	return *this;
}

nxBinaryFile& nxBinaryFile::operator>>( nxString& str )
{
	int nbytes = 0;

	m_ok = m_ok && (Read( &nbytes, sizeof(int), 1) == 1);
	if (m_ok)
	{
		if (nbytes < 1)
		{
			str.Empty(nxFALSE);
		}
		else
		{
			char *buffer = new char[nbytes+1];
			if (buffer == NULL) NXTHROW (("nxBinaryFile memory allocation error"));
			m_ok = m_ok && (Read( buffer, sizeof(char), nbytes) == (size_t)nbytes);
			buffer[nbytes] = '\0';
			str = buffer;
			delete [] buffer;
		}
	}
	else
	{
		str.Empty(nxFALSE);
	}
	return *this;
}

nxBinaryFile& nxBinaryFile::operator<<( const double val )
{
	m_ok = m_ok && (Write( &val, sizeof(double), 1 ) == 1);
	return *this;
}

nxBinaryFile& nxBinaryFile::operator>>( double& val )
{
	m_ok = m_ok && (Read( &val, sizeof(double), 1 ) == 1);
	return *this;
}

nxBinaryFile& nxBinaryFile::operator<<( const nxBOOL bval )
{
	int val = bval ? 1 : 0;
	m_ok = m_ok && (Write( &val, sizeof(val), 1 ) == 1);
	return *this;
}
nxBinaryFile& nxBinaryFile::operator>>( nxBOOL& bval )
{
	int val = 0;
	m_ok = m_ok && (Read( &val, sizeof(val), 1 ) == 1);
	bval = (val != 0);
	return *this;
}
nxBinaryFile& nxBinaryFile::operator<<( const float val )
{
	m_ok = m_ok && (Write( &val, sizeof(float), 1) == 1);
	return *this;
}
nxBinaryFile& nxBinaryFile::operator>>( float& val )
{
	m_ok = m_ok && (Read( &val, sizeof(float), 1) == 1);
	return *this;
}
nxBinaryFile& nxBinaryFile::operator<<( const nxBYTE val )
{
	m_ok = m_ok && (Write( &val, sizeof(nxBYTE), 1) == 1);
	return *this;
}
nxBinaryFile& nxBinaryFile::operator>>( nxBYTE& val )
{
	m_ok = m_ok && (Read( &val, sizeof(nxBYTE), 1) == 1);
	return *this;
}

nxBinaryFile& nxBinaryFile::operator<<( const nxWORD val )
{
	m_ok = m_ok && (Write( &val, sizeof(nxWORD), 1) == 1);
	return *this;
}
nxBinaryFile& nxBinaryFile::operator>>( nxWORD& val )
{
	m_ok = m_ok && (Read( &val, sizeof(nxWORD), 1) == 1);
	return *this;
}

nxBinaryFile& nxBinaryFile::operator<<( const nxDWORD val )
{
	m_ok = m_ok && (Write( &val, sizeof(nxDWORD), 1) == 1);
	return *this;
}
nxBinaryFile& nxBinaryFile::operator>>( nxDWORD& val )
{
	m_ok = m_ok && (Read( &val, sizeof(nxDWORD), 1) == 1);
	return *this;
}

nxBinaryFile& nxBinaryFile::operator<<( const int val )
{
	m_ok = m_ok && (Write( &val, sizeof(int), 1) == 1);
	return *this;
}

nxBinaryFile& nxBinaryFile::operator>>( int& val )
{
	m_ok = m_ok && (Read( &val, sizeof(int), 1) == 1);
	return *this;
}

