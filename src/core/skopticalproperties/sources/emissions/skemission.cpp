#include <skopticalproperties21.h>


/*-----------------------------------------------------------------------------
 *					skEmission::IsotropicEmissionArray		 2015- 3- 11*/
/** Default implementation of IsotropicEmissionArray implements a simple
 *	single threaded loop.
 **/
/*---------------------------------------------------------------------------*/

bool skEmission::IsotropicEmissionArray( const std::vector<double>& wavenumber, std::vector<double>* isotropicradiance)
{
	bool	ok = true;
	bool	ok1;

	isotropicradiance->resize( wavenumber.size() );
	for (size_t i = 0; i < wavenumber.size(); i++)
	{
		ok1 = IsotropicEmission( wavenumber.at(i), &isotropicradiance->at(i) );
		ok = ok && ok1;
	}
	return ok;
}
