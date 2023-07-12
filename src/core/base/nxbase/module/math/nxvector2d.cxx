#include  "nxbase_math.h"

nxVector2D	operator*( double mult,     const nxVector2D& other )
{
	return other*mult;
}
nxVector2D	operator+( double constant, const nxVector2D& other )
{
	return other+constant;
}

nxVector2D	operator-( double constant, const nxVector2D& other )
{
	return other-constant;
}

