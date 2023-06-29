#include <skopticalproperties21.h>
#include <limits>
/*-----------------------------------------------------------------------------
 *					skEmission_Tabulated_HeightWavelength::skEmission_Tabulated_HeightWavelength		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

skEmission_Tabulated_HeightWavelength::skEmission_Tabulated_HeightWavelength()
{
	m_idxh1      = 0;
	m_idxh2      = 0;
	m_h1         = 0.0;
	m_h2         = 0.0;
	m_isground = false;
}


/*-----------------------------------------------------------------------------
 *					skEmission_Tabulated_HeightWavelength::LookupIndicesAndWeights		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_Tabulated_HeightWavelength::LookupIndicesAndWeights( const nx1dArray<double>& h, double value, double* w1, size_t* idx1, double* w2, size_t* idx2 ) const
{
	nxArrayIter<double>	start;
	nxArrayIter<double>	finish;
	nxArrayIter<double>	iter2;
	bool				ok;
	double				h1;
	double				h2;
	double				f;


	ok     = (h.size() > 0);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skEmission_Tabulated_HeightWavelength::LookupIndicesAndWeights, Cannot lookup up height indices as the height array is empty");
		*idx1 = 0;
		*w1   = 0.0;
		*idx2 = 0;
		*w2   = 0;
	}
	else
	{
		start  = h.begin();
		finish = h.end();
		iter2 = std::upper_bound( start, finish, value);
			
		if (iter2 == start )										// If we are befor ethe first element
		{															// then truncate to the first element
			*idx1 = 0;
			*w1   = 1.0;
			*idx2 = 0;
			*w2   = 0;
		}
		else if (iter2 == finish)									// if we are at the end of the array
		{															// then truncate to the last element
			*idx1 = h.size() -1 ;
			*w1   = 1.0;
			*idx2 = 0;
			*w2   = 0.0;
		}
		else
		{
			*idx2   = (iter2 - start);
			*idx1   = *idx2 - 1;
		
			h2      = h.At(*idx2); 
			h1      = h.At(*idx1);
			f       = (value-h1)/(h2-h1);
			*w1     = (1-f);
			*w2     = f;
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skEmission_Tabulated_HeightWavelength::SetAtmosphericState		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_Tabulated_HeightWavelength::UpdateLocation( const GEODETIC_INSTANT& pt, bool isground)
{
	bool	ok;

	m_isground = isground;
	ok = m_isground;
	if (ok)
	{
		m_h1 = 0;
		m_h2 = 0;
	}
	else
	{
		ok = LookupIndicesAndWeights( m_heights,pt.heightm, &m_h1, &m_idxh1, &m_h2, &m_idxh2 );
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skEmission_Tabulated_HeightWavelength::CalculateCrossSections		2014-2-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skEmission_Tabulated_HeightWavelength::IsotropicEmission( double wavenumber, double* emission)
{
	double		wavelen;
	double		w1;
	double		w2;
	size_t		idx1;
	size_t		idx2;
	double		exth1;
	double		exth2;
	double		ext;
	bool		ok;

	ok = m_isground;
	if (ok)
	{
		*emission = 0.0;
	}
	else
	{
		wavelen = 1.0E7/wavenumber;
		ok         = LookupIndicesAndWeights( m_wavelennm, wavelen, &w1, &idx1, &w2, &idx2 );
		exth1      = w1*m_emission.At(idx1, m_idxh1) + w2*m_emission.At(idx2,m_idxh1);			// Get the extinction at required wavelength at lower altitude
		exth2      = w1*m_emission.At(idx1, m_idxh2) + w2*m_emission.At(idx2,m_idxh2);			// Get the extinction at required wavelength at upper altitude
		ext        = m_h1*exth1 + m_h2*exth2;
		*emission  = ext;
	}
	return ok;
}




/*-----------------------------------------------------------------------------
 *					ReleaseResources		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

void skEmission_Tabulated_HeightWavelength:: ReleaseResources()
{
	m_idxh1      = 0;
	m_idxh2      = 0;
	m_h1         = 0.0;
	m_h2         = 0.0;
	m_heights.erase();
	m_wavelennm.erase();
	m_emission.erase();
}

/*-----------------------------------------------------------------------------
 *					skEmission_Tabulated_HeightWavelength::CreateClone		2009-5-19*/
/** **/
/*---------------------------------------------------------------------------*/

/*
bool skEmission_Tabulated_HeightWavelength::CreateClone( skOpticalProperties ** userclone ) const
{
	bool													ok;
	skEmission_Tabulated_HeightWavelength*	clone;

	clone = new skEmission_Tabulated_HeightWavelength;
	ok = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skEmission_Tabulated_HeightWavelength::CreateClone, Error creating allocating memory for clone");
	}
	else
	{
		clone->AddRef();
		clone->m_isdirty    = true;
		clone->m_isabsorber = m_isabsorber;
		clone->m_idxh1      = m_idxh1;
		clone->m_idxh2      = m_idxh2;
		clone->m_h1         = m_h1;
		clone->m_h2         = m_h2;
		ok =       clone->skOpticalProperties::DeepCopy( *this );
		ok = ok && clone->m_heights.DeepCopy	( m_heights);
		ok = ok && clone->m_wavelennm.DeepCopy	( m_wavelennm );
		ok = ok && clone->m_emission.DeepCopy	( m_emission );
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skEmission_Tabulated_HeightWavelength::CreateClone, Error making clone, I'm going to erase the clone to blanks");
		clone->ReleaseResources();
	}
	*userclone = clone;
	return ok;
};
*/

/*-----------------------------------------------------------------------------
 *					skEmission_Tabulated_HeightWavelength::SetExtinctionTable		2009-5-19*/
/** Sets up the table that will be used for extinction (or absorption) as a
*	function of altitude and wavelength. The user passes in a 2-D array of extinction
*	values and two 1D arrays that define the wavelength and height "axes" of the grid.
*
*	\param extinction
*	A 2D array of atmopsheric isotropic emission. The units are photons/nm/sec/ster
*	The array is double( numwavelens, numheights ) where numwavelens is the number of
*	wavelengths in parameter wavelens and numheights is the number of heights in parameter
*	heights_meets. The code linearly interpolates the table in both height and wavelength.
*	It truncates the interpolation at the end points of the array.
*
*	\param wavelens
*	The wavelengths at which the table is specified. The units are nm. Must be 1 or more. All elements must
*	be provided in ascending order.
*
*	\param heights_meters
*	The heights at which the table is specified. The units are in meters. Muts be 1 or more values. All elements
*	must be provided in ascending order.
**/
/*---------------------------------------------------------------------------*/

bool skEmission_Tabulated_HeightWavelength::SetEmissionTable( const nx2dArray<double>& emission, const nx1dArray<double>& wavelens, const nx1dArray<double>& heights_meters)
{
	const size_t*		dims;
	bool				ok;

	dims = emission.ArrayRankSpecs()->Dims();
	ok   =	   ( dims[0] == wavelens.size()) 
			&& ( dims[1] == heights_meters.size())
			&& ( dims[0] > 0 )
			&& ( dims[1] > 0 );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skEmission_Tabulated_HeightWavelength::SetExtinctionTable, Error setting table as array sizes are incompatible with each other. Did you mix up the heights and wavelengths. Wavelength varies fastest in array and height varies slowest");
		ReleaseResources();
	}
	else
	{
		ok =       m_heights.DeepCopy	( heights_meters);
		ok = ok && m_wavelennm.DeepCopy	( wavelens );
		ok = ok && m_emission.DeepCopy( emission );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skEmission_Tabulated_HeightWavelength::SetExtinctionTable, Error copying arrays over. I'm clearing array to zero");
		}
	}
	if (!ok)
	{
		ReleaseResources();
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skEmission_Tabulated_HeightWavelength::LoadHeightWavelengthProfileFromFile		2009-6-26*/
/** Loads in the tabulated extinction from a file.
 *
 *	Line 1:					numwave	numheights			number of wavelengths and number of heights
 *	Line 2:					wavelengths(numwavelen)		expressed in nanometers.
 *	Line 3:					heights(numheights)			expressed in meters.
 *	Line 4:					extinction(numwave)			extinction for wavelengths at height (0)
 *	Line 5:					extinction(numwave)			extinction for wavelengths at height (1)
 *	..
 *	Line (numheights+3):	extinction(numwave)			extinction for wavelengths at height (numheights-1)
**/
/*---------------------------------------------------------------------------*/

bool skEmission_Tabulated_HeightWavelength::LoadHeightWavelengthProfileFromFile( const char* filename )
{
	std::ifstream			strm;
	size_t					numh;
	size_t					numwave;
	nx1dArray<double>		wavelens;
	nx1dArray<double>		heights;
	nx2dArray<double>		extinction;
	bool					ok;

	strm.open(filename, std::ios_base::in );
	strm  >> numwave >> numh;							// Load in the number of heights and the number of wavelengths
	ok = (numh > 0) && (numwave > 0) && !strm.fail();	// Make sure it is good
	ok = ok && wavelens.SetSize(numwave);
	ok = ok && heights.SetSize (numh );
	ok = ok && extinction.SetSize( numwave, numh );
	if (ok)
	{
		strm >> wavelens; 
		strm >> heights;
		strm >> extinction;
		ok = ok && !strm.fail();
	}
	strm.close();
	ok = ok && SetEmissionTable( extinction, wavelens, heights);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skEmission_Tabulated_HeightWavelength::LoadHeightWavelengthProfileFromFile, error loading and setting the height extinction table from file <%s>", (const char*) filename);
		ReleaseResources();
	}
	return ok;
}


