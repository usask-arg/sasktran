#include "../sktran_common.h"

bool SKTRAN_Engine_Base::CalculateRadiance( std::vector<SKTRAN_StokesScalar>*		losradiance,
											double									wavelen,
											size_t								    numordersofscatter,
											SKTRAN_AtmosphericOpticalState_V21*	    opticalstate,
											std::vector<skRTStokesVector>*          losvector, 
											bool									updateclimatology,
											SKTRAN_DiagnosticInterface*			    diag)
{
	nxLog::Record(NXLOG_ERROR,"SKTRAN_Engine_Base::CalculateRadiance, You should not be calling the base version. This is an implementation error");
	return false;
}


/*---------------------------------------------------------------------------
 *      SKTRAN_Engine_Base::CalculateMultiWavelengthRadiance      2020-06-07 */
/** Implements a default implementation o fthe multi-wavelength call to
 *	sasktran that simply calls the single wavelength version multiple times
**/
/*---------------------------------------------------------------------------*/

bool SKTRAN_Engine_Base::CalculateMultiWavelengthRadiance(  std::vector< std::vector<SKTRAN_StokesScalar> > *		losradiance,
															const std::vector<double>&								wavelen,
															size_t													numordersofscatter,
															SKTRAN_AtmosphericOpticalState_V21*						opticalstate,
															std::vector< std::vector<skRTStokesVector> > *          losvector, 
															bool													updateclimatology,
															SKTRAN_DiagnosticInterface*								diag)
{
	nxLog::Record(NXLOG_ERROR,"SKTRAN_Engine_Base::CalculateMultiWavelengthRadiance, You should not be calling the base version. This is an implementation error");
	return false;

}
/*
	size_t											numwave = wavelen.size();
	bool											ok = true;
	bool											ok1;
	std::vector<SKTRAN_StokesScalar>				results;
	std::vector<skRTStokesVector>					losvectorarray;
	std::vector<skRTStokesVector>*					losvectorptr = &losvectorarray;

	losradiance->clear();
	losradiance->reserve(numwave);
	if (losvector != nullptr)
	{
		losvector->clear();
		losvector->reserve(numwave);
		losvectorptr = nullptr;
	}

	for (size_t iw = 0; iw < numwave; iw++)
	{
		results.clear();
		if (losvectorptr != nullptr) losvectorptr->clear();

		ok1 = CalculateRadiance( &results, wavelen[iw], numordersofscatter, opticalstate, losvectorptr, (updateclimatology && (iw ==0 )), diag);
		if (losvector != nullptr) losvector->push_back( losvectorarray );
		losradiance->push_back( results );
		ok = ok && ok1;
	}
	return ok;
}
*/


