#include "include/sktran_hr_internals.h"

bool SKTRAN_HR_Diffuse_Table_SZA::SZAWeights( double cossza, SKTRAN_HR_WEIGHT_TYPE* szaweights, size_t* szaindex, size_t& numindex ) const
{
	bool ok = true;

	if (m_cosszas.size() == 1)
	{
		szaindex[0] = 0;
		szaweights[0] = 1;
		numindex = 1;
	}

	std::vector<double>::const_iterator	low,up;
	up = std::upper_bound( m_cosszas.begin(), m_cosszas.end(), cossza );
	if( up == m_cosszas.end() )
	{
		up = m_cosszas.end() - 1;
	}
	if( up == m_cosszas.begin() )
	{
		low = m_cosszas.begin();
	}
	else
	{
		low = up - 1;
	}
	if( fabs(*up - *low) < 1E-8 )
	{
		numindex = 1;
		szaweights[0] = 1;
		szaindex[0] = low - m_cosszas.begin();
	}
	else if ( cossza > *up )
	{
		numindex = 1;
		szaweights[0] = 1;
		szaindex[0] = up - m_cosszas.begin();
	}
	else if ( cossza < *low )
	{
		numindex = 1;
		szaweights[0] = 1;
		szaindex[0] = 0;
	}
	else
	{
		numindex = 2;
		szaweights[0] = SKTRAN_HR_DBL_TO_WEIGHT((cossza - *low) / (*up - *low));
		szaweights[1] = SKTRAN_HR_DBL_TO_WEIGHT((*up - cossza) / ( *up - *low));
		szaindex[0] = up - m_cosszas.begin();
		szaindex[1] = low - m_cosszas.begin();
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Table_SZA::ChooseDiffusePoints		 2017- 2- 8*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Diffuse_Table_SZA::ChooseDiffusePoints( const HELIODETIC_POINT&				pt, 
													   size_t*								diffuseindex, 
													   SKTRAN_HR_WEIGHT_TYPE*				diffuseweights, 
													   size_t*								usernumpoints ) const
{
	bool								ok;
	bool								ok1;
	SKTRAN_HR_WEIGHT_TYPE				altweights[2];
	SKTRAN_HR_WEIGHT_TYPE				szaweights[2];
	size_t								altindex[2];
	size_t								szaindex[2];
	size_t								numpoints;
	size_t								numalt, numsza;

	ok = true;
	numpoints = 0;
	ok = ok && SZAWeights( pt.CosSZA(), szaweights, szaindex, numsza );
	for( size_t profileidx = 0; profileidx < numsza; profileidx++ )
	{
		numalt = N_ELEMENTS(altweights);
		ok1 = AltWeightsForProfile( pt.Altitude(), profileidx, altweights, altindex, &numalt );
		ok1 = ok1 && (numpoints < *usernumpoints);
		if (ok1)
		{
			for( size_t altidx = 0; (altidx < numalt) && ok; altidx++ )
			{
				diffuseindex[numpoints] =  (	ProfileStartIndex()[szaindex[profileidx]] + altindex[altidx] );
				diffuseweights[numpoints] = ( szaweights[profileidx] * altweights[altidx] );
				++numpoints;
				ok1 = ok1 && (numpoints <= *usernumpoints);
			}
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"SKTRAN_HR_Diffuse_Table_SZA::ChooseDiffusePoints, Insufficient space passed in by user to store the interpolated SZA diffuse points");
			}
		}
		ok = ok && ok1;
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"SKTRAN_HR_Diffuse_Table_SZA::ChooseDiffusePoints, There were errors choosing the sza point to interpolate a point. Thats not good.");
	}
	*usernumpoints = ok ? numpoints : 0;
	return ok;
}
bool SKTRAN_HR_Diffuse_Table_SZA::ChooseGroundPoints( const HELIODETIC_POINT&				pt,
													  size_t*								diffuseindex,
													  SKTRAN_HR_WEIGHT_TYPE*				diffuseweights,
													  size_t&								numpoints ) const
{
	bool ok = true;

	SKTRAN_HR_WEIGHT_TYPE				szaweights[2];
	size_t								szaindex[2];
	size_t								numsza;
	numpoints = 0;

	ok = ok && SZAWeights( pt.CosSZA(), szaweights, szaindex, numsza );
	for( size_t profileidx = 0; profileidx < numsza; profileidx++ )
	{
		diffuseindex[numpoints] = (	GroundStartIndex() + szaindex[profileidx] );
		diffuseweights[numpoints] = ( szaweights[profileidx] );
		++numpoints;
	}

	return ok;
}


void SKTRAN_HR_Diffuse_Table_SZA::SetNumProfileInterp( size_t numinterp )
{
	m_cosszas.resize( ProfileStartIndex().size() - 1 );

	for( size_t szaidx = 0; szaidx < m_cosszas.size(); szaidx++ )
	{
		m_cosszas[szaidx] = DiffusePoints()[ProfileStartIndex()[szaidx]].Location().CosSZA();
	}
}