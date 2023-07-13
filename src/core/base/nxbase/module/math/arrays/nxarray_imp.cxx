#include "nxbase_math.h"
//#include "module/math/arrays/nxarray.h"
//#include "module/math/arrays/nx2darr.h"

/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/


nxBOOL TransferElement( nxString* array, int , nxBinaryFile* outfile, nxBOOL input)
{
//	nxBinaryFile file( outfile );

	if (input) (*outfile) >> (*array);
	else       (*outfile) << (*array);
	return outfile->IsOk();
}

nxString NullValue( const nxString&  )
{
	return "";
}

