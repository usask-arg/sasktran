#include <skopticalproperties21.h>
#include <iomanip>

// array of the avg IWC for each of the 18 Baum mean PSDs (g/m3)
static double baum_iwc[18] =
{
	0.0059,
	0.0028,
	0.0019,
	0.0012,
	0.0015,
	0.0277,
	0.0814,
	0.0462,
	0.2574,
	0.0876,
	0.1409,
	0.0614,
	0.1085,
	0.1032,
	0.1214,
	0.1036,
	0.0803,
	0.1401
};
// fractions, 'f', scattered into forward direction using standard criteria used in this class
static double trunc_fraction[18] =
{
	0.384,
	0.420,
	0.470,
	0.464,
	0.519,
	0.553,
	0.573,
	0.596,
	0.615,
	0.635,
	0.642,
	0.656,
	0.650,
	0.648,
	0.643,
	0.640,
	0.641,
	0.629
};

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB		2013-1-9*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB()
{
	m_isabsorber		= false;
	m_numfilecolumns	= 13;
//	m_crossection		= NULL;
//	m_deltafraction		= NULL;
//	m_heights			= NULL;
//	m_wavelennm			= NULL;			// table wavelengths in nm
//	m_scatterangles_db	= NULL;
//	m_truncangles_db	= NULL;
	SetDBArrays( );
	SetDefaultTruncationAngles();
}
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::~skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB		2013-1-9*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::~skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB()
{
	ReleaseResources();
}
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::ReleaseResources		2013-1-9*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::ReleaseResources()
{
	m_idxh1     = 0;
	m_idxh2     = 0;
	m_f1		= 0.0;
	m_f2		= 0.0;
	m_lowerheight = NULL;
	m_upperheight = NULL;
	m_heights.erase();
	m_wavelennm.erase();
	m_wavelennm_db.erase();
	m_effectivesize_db.erase();
	m_crossection.erase();
	m_deltafraction.erase();
	m_heightentries.clear();
	m_truncangles_db.erase();
	m_scatterangles_db.erase();
}

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetLocalDirectory( const char *configdir )
{
	m_configdir.sprintf( "%s",configdir );
	return true;
}

/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetDBArrays		2010-1-18
 *	Define the default effective sizes, wavelengths, and scattering angles. 
 *	For scatter angles, a list of angular sections with different resoutions in 
 *	each section is supplied to sample the forward scatter peak near at high 
 *	resolution while the other angles use a much coarser resolution.
 *---------------------------------------------------------------------------*/
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetDBArrays()
{
	bool			ok = true;
	double			lastend = 0.0;
	size_t			i,j,idx,numangles,nsteps,numrange;

	nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetDBArrays, Effective Size arrays, wavlength and phase matric angles should be acquired from files rather than hard coded into code."
		                         "I'm pretty sure the current wavelength code is wrong. It generates the wrong numbers as there is a gap in Baum between 1.00 and 1.20 microns");

	if ( m_effectivesize_db.size()==0 )											// If not already defined then 
	{																			// defie the baum database.
		m_effectivesize_db.SetSize( 18 );										// The baum database comes as 18 different sizes of De
		for(idx=0; idx < m_effectivesize_db.size(); idx++)
		{
			m_effectivesize_db.At(idx) = 10.0 + idx*10.0;
		}
	}
	if ( m_wavelennm_db.size()==0 )												// If not already defined
	{																			// then
		m_wavelennm_db.SetSize( 144 );											// There are 144 standard wavelengths
		for(idx=0; idx < m_wavelennm_db.size(); idx++)							// that step from 400nm to 1.
		{
			m_wavelennm_db.At(idx) = 1e3*( 0.4 + idx*0.01 );
		}
	}
	numrange	= 6;
	double	startangle[6] = {0.0,  2.0,  5.0,  10.0,  15.0, 176.0};
	double  endangle[6]   = {2.0,  5.0, 10.0,  15.0, 176.0, 180.0};
	double  resolution[6] = {0.01, 0.05, 0.1,   0.5,   1.0,  0.25};

	numangles = 0;
	for (i=0; i < numrange && ok; i++)
	{
		ok = (startangle[i] >= lastend) && (endangle[i] <= 180.0);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::ConfigureScatterAngles, The startangle (%g) must be less than lastangle (%g) and endangle (%g) must be <= 90 degrees", (double)startangle[i], (double)lastend, (double)endangle[i]);
		}
		else
		{
			nsteps     = (size_t)((endangle[i] - startangle[i])/resolution[i] + 0.5);
			numangles += nsteps;
			lastend    = startangle[i] + (nsteps-1)*resolution[i];;
		}
	}
	if (ok)
	{
		numangles++;
		ok = m_scatterangles_db.SetSize(numangles);			
		idx	= 0;
		for (i=0;i<numrange && ok; i++)
		{
			nsteps   = (size_t)((endangle[i] - startangle[i])/resolution[i] + 0.5);
			for (j=0;j<nsteps;j++)
			{
				m_scatterangles_db.At(idx++) =  startangle[i] + j*resolution[i];
			}
		}
		m_scatterangles_db.At(idx) =  180.0;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetDefaultTruncationAngles		2013-1-9*/
/**	Set the default truncation angles for all DB effective sizes.  Selections are 
 *	based on minimizing the absolute geometric mean (over inc. directions) of the 
 *	photon-conservation multipliers.
 *	For reference on these selections, see Wiensz et al, JQSRT, 2012. **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetDefaultTruncationAngles( )
{
	size_t		idx;
	if ( m_truncangles_db.size()==0 )
	{
		m_truncangles_db.SetSize( 18 );
		for(idx=0;idx<3;idx++)				// 10..40 um: use 5 degrees
		{
			m_truncangles_db.At(idx) = 5.0;
		}
		for(idx=3;idx<18;idx++)				// 50..180 um: use 2 degrees
		{
			m_truncangles_db.At(idx) = 2.0;
		}
	}

	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadDatabaseDirectoryFromRegistry		2013-1-9*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadDatabaseDirectoryFromRegistry( )
{
	nxRegistryConfiguration		config( "USask-ARG", "skOpticalProperties/Baum_IceCrystals/",nxRegistryConfiguration::GLOBAL_INI, true);
	nxString					filename;
	bool						ok;

	ok = config.LocateDirectoryFromKey( "DatabaseDirectory", &filename, true, true, "Enter the directory of the Baum in-situ ice crystal parameter database");
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadDatabaseDirectoryFromRegistry, error loading aersol database directory name from registry, using TEMP");
		filename = getenv("TEMP");
	}
	return filename;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadStorageDirectoryFromRegistry		2013-1-9*/
/** **/
/*---------------------------------------------------------------------------*/

nxString skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadStorageDirectoryFromRegistry( )
{
	nxRegistryConfiguration		config( "USask-ARG", "skOpticalProperties/Baum_IceCrystals/",nxRegistryConfiguration::GLOBAL_INI, true);
	nxString					filename;
	bool						ok;

	ok = config.LocateDirectoryFromKey( "CacheDirectory", &filename, true, true, "Enter the directory to store the Baum in-situ ice crystal cached file");
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadStorageDirectoryFromRegistry, error loading aersol cache directory name from registry, using TEMP");
		filename = getenv("TEMP");
	}
	return filename;
}

/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LookupIndicesAndWeights			2011-4-18
 *	Find the indices and weights for interpolation - either on a wavelength or 
 *	on a height basis.
 *---------------------------------------------------------------------------*/
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LookupIndicesAndWeights( const nx1dArray<double>& h, double value, double* w1, size_t* idx1, double* w2, size_t* idx2 ) const
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
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LookupIndicesAndWeights, Cannot llokup up height indices as the height array is empty");
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
		else														// otherwise, linearly interpolate
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
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetAtmosphericState			2011-4-18
 * 
 *---------------------------------------------------------------------------*/
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetAtmosphericState( skClimatology* /*neutralatmosphere*/, const GEODETIC_INSTANT& pt, bool* crosssectionschanged  )
{
	bool	ok;

	ok =		LookupIndicesAndWeights( m_heights,pt.heightm, &m_f1, &m_idxh1, &m_f2, &m_idxh2 );
	ok = ok &&	ConfigureHeightInterpolation( pt.heightm );
	m_isdirty = true;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetAtmosphericState, There was an error setting the atmospheric state for the specified point. That might be a subtle problem as it disables this objects influence on any calculation");
	}
	if (crosssectionschanged != NULL) *crosssectionschanged = m_isdirty;					// Return flag stating whether cross-sections have changed
	return ok;
}


/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::CalculateCrossSections			2011-4-18
 *	Interpolate the cross sections by height and wavelength. Do the same for 
 *	delta-function fraction 'f'
 *---------------------------------------------------------------------------*/
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs, double* forwardscatterfrac  ) const
{
	double		wavelen;
	double		w1;
	double		w2;
	size_t		idx1;
	size_t		idx2;
	double		sigmah1;
	double		sigmah2;
	double		sigma;
	double		f1,f2;
	bool		ok;

	wavelen = 1.0E7/wavenumber;
	ok      = LookupIndicesAndWeights( m_wavelennm, wavelen, &w1, &idx1, &w2, &idx2 );
	sigmah1 = w1*m_crossection.At(idx1, m_idxh1) + w2*m_crossection.At(idx2,m_idxh1);			// Get extinction at required wavelength at lower altitude
	sigmah2 = w1*m_crossection.At(idx1, m_idxh2) + w2*m_crossection.At(idx2,m_idxh2);			// Get extinction at required wavelength at upper altitude
	sigma   = m_f1*sigmah1 + m_f2*sigmah2;
	f1		= w1*m_deltafraction.At( idx1,m_idxh1) + w2*m_deltafraction.At(idx2,m_idxh1);	// Get delta fraction 'f' at lower alt
	f2		= w1*m_deltafraction.At( idx1,m_idxh2) + w2*m_deltafraction.At(idx2,m_idxh2);	// Get delta fraction 'f' at upper alt
	*forwardscatterfrac = m_f1*f1 + m_f2*f2;	

	if (m_isabsorber)
	{
		*absxs = sigma;
		*extxs = sigma;
		*scattxs = 0.0; 
	}
	else
	{
		*absxs   = 0.0;
		*extxs   = sigma;
		*scattxs = sigma; 
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::CalculateCrossSections		2014-2-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::CalculateCrossSections( double wavenumber, double* absxs, double* extxs, double* scattxs )
{
	double	forwardscatterfrac;
	bool	ok;
	
	ok = CalculateCrossSections( wavenumber, absxs, extxs, scattxs, &forwardscatterfrac  );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::CalculatePhaseMatrix		2014-2-21*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::CalculatePhaseMatrix( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix) const
{
	bool															ok;
	double															wavelennm;
	const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry*		upperentry = nullptr;
	const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry*		lowerentry = nullptr;
	double															upperp = 0.0;
	double															lowerp = 0.0;
	

	wavelennm = 1.0E7/wavenumber;
	phasematrix->SetTo(0.0);

	ok =       (m_upperheight != NULL) && (m_lowerheight != NULL);
	ok = ok && m_upperheight->FindWavelengthEntry(wavelennm, &upperentry);
	ok = ok && m_lowerheight->FindWavelengthEntry(wavelennm, &lowerentry);
	ok = ok && upperentry->GetPhaseMatrix( cosscatterangle, &upperp );
	ok = ok && lowerentry->GetPhaseMatrix( cosscatterangle, &lowerp );	// was upperentry-> ...
	
	phasematrix->At(1,1) = (SKRTFLOAT)(m_f1*lowerp + m_f2*upperp);
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::CalculatePhaseMatrix		2013-1-9*/
/**	Return the value of the phase matrix for this wavenum and (cos)scatter angle.
 *	NOTE: This function works with database phase functions that have already been
 *	truncated, and these will be normalized according to 
 *		f=\int[P_{smooth}(\cos\Theta)/2]d\cos\Theta = (1-f),
 *	where f is the fraction of light scattering into the forward direction.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::CalculatePhaseMatrix( double wavenumber , double cosscatterangle, skRTPhaseMatrix* phasematrix)
{
	return CalculatePhaseMatrix( wavenumber , cosscatterangle, phasematrix);
}

/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::ConfigureHeightInterpolation		2011-5-07*/
/** Finds the entries that bound the specified height and configures the linear combination.
 *	The algorithm sets phase matrices to zero for heights outside the upper and lower bounds
 *	of the set of heights. This is important if the phase matrix is only defined 
 *	for tropospheric altitudes as we do not want it being used all the way up.
**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::ConfigureHeightInterpolation( double heightm )
{
	skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry		dummy(heightm);
	iterator											iter;
	bool												ok;
	double												h1;
	double												h2;
	double												dh;

	iter = std::lower_bound( m_heightentries.begin(), m_heightentries.end(), dummy);		// Find the entry that is "bigger" than this entry
	ok   = !(iter == m_heightentries.begin()) && !(iter == m_heightentries.end());			// ok if we are inside the range of altitudes
	if (ok )
	{
		m_upperheight = &(*iter);
		--iter;
		m_lowerheight = &(*iter);
		h2            = m_upperheight->Height();
		h1            = m_lowerheight->Height();
		dh            = h2-h1;
		m_f1          = (h2-heightm)/dh;
		m_f2          = 1-m_f1;
	}
	else
	{
		m_upperheight = NULL;
		m_lowerheight = NULL;
		m_f1          = 0.0;
		m_f2          = 0.0;
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::CreateClone			2011-4-18
 *	
 *---------------------------------------------------------------------------*/
/*
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::CreateClone( skOpticalProperties ** userclone ) const
{
	skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB*			clone=NULL;
	bool																ok;
	bool																ok1;
	const_iterator														iter;
	skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry		newentry;
	skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry*		entry;

	clone = new skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB;
	ok    = (clone != NULL);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::CreateClone, Memory allocation error creating clone. Thats a problem");
	}
	else
	{
		clone->AddRef();
		clone->m_isdirty			= true;
		clone->m_isabsorber			= m_isabsorber;
		clone->m_configdir			= m_configdir;
		clone->m_idxh1				= m_idxh1;
		clone->m_idxh2				= m_idxh2;
		clone->m_f1					= m_f1;
		clone->m_f2					= m_f2;
		clone->m_numfilecolumns		= m_numfilecolumns;	
		ok =       clone->skOpticalProperties::DeepCopy( *this );
		ok = ok && clone->m_heights.DeepCopy			( m_heights);
		ok = ok && clone->m_effectivesize_db.DeepCopy	( m_effectivesize_db );
		ok = ok && clone->m_wavelennm_db.DeepCopy		( m_wavelennm_db );
		ok = ok && clone->m_crossection.DeepCopy		( m_crossection );
		ok = ok && clone->m_deltafraction.DeepCopy		( m_deltafraction );
		ok = ok && clone->m_truncangles_db.DeepCopy		( m_truncangles_db );
		ok = ok && clone->m_scatterangles_db.DeepCopy	( m_scatterangles_db );

		for (iter = m_heightentries.begin(); !(m_heightentries.end() == iter); ++iter)
		{
			clone->m_heightentries.push_back(newentry);
			entry = &(clone->m_heightentries.back());
			ok1 = entry->DeepCopy( *iter );
			ok = ok && ok1;
		}
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::CreateClone, error making copy. Thats a problem.");
			clone->ReleaseResources();
		}
	}
	*userclone = clone;
	return ok;
}
*/

/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::GetXSectsFromDBAndSetMembers			2011-4-18
 *	Given the wavelengths to be modeled, as well as the effective sizes and number
 *	densities (both as a function of height), get the cross sections from the Baum DB,
 *	convert to extinctions, and populate the member extinction array.  
 *	
 *---------------------------------------------------------------------------*/
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::GetXSectsFromDBAndSetMembers( const nx1dArray<double>& heightsm,const nx1dArray<double>& reff,const std::vector<double>& lambdasnm,nx2dArray<double> *cloudcrosssection )
{	
	bool				ok;
	nxString			db_xsectfile,local_crosssectfile;
	size_t				idxl,idxh,idxx;
	size_t				numLambdas,numHeights;
	double				ksca_background(1e-30),detectvectlambda(750);
	nx2dArray<double>	crosssection;
	nx2dArray<double>	db_xsectdata;
	std::vector<double>	db_lambdas;

	local_crosssectfile.sprintf( "%s/IceCrossSections.dat",(const char *)m_configdir );
	db_lambdas	= m_wavelennm_db.STLVector();
	numLambdas	= lambdasnm.size();
	numHeights	= reff.size();
	m_heights.SetSize( numHeights );
	ksca_background		= 1e-30;
	detectvectlambda	= 750.0;		// cloud detection vector w/l: write scatt. extinct to array for this wavelen
	ok =		crosssection.SetSize( lambdasnm.size(),heightsm.size() );
	ok = ok &&	cloudcrosssection->SetSize( heightsm.size(),2 );
	for (idxh=0;idxh<numHeights;idxh++)		
	{
		if ( IsValidDatabaseEffectiveSize(reff.At(idxh)) )	// if reff is listed size, get cross sect from DB
		{
			db_xsectfile.sprintf( "%s/Solar_mix1_d%03u.dat",(const char *)LoadDatabaseDirectoryFromRegistry(),(unsigned int)(reff.At(idxh)) );
			ok		= db_xsectdata.InputColumnMajorText( (const char *)db_xsectfile,m_numfilecolumns );
			if (ok)
			{
				for (idxl=0;idxl<numLambdas;idxl++)			// for each model wavelength, find
				{											// DB equivalent lambda
					idxx	= (size_t)(std::lower_bound( db_lambdas.begin(),db_lambdas.end(),lambdasnm.at(idxl) )-db_lambdas.begin());
					if (idxx >= 0)							// if this wavelength is valid,
					{										// get scat. x-sect
						crosssection.At(idxl,idxh) = db_xsectdata.At(idxx,m_numfilecolumns-2)*1.0e-8;	// convert x-sect in um^2 to cm^2
					}
					else									// otherwise 
					{
						crosssection.At(idxl,idxh) = ksca_background;	// use insignif. value
					}
					if (fabs(lambdasnm.at(idxl)-detectvectlambda)<0.75)	// ndl303, 2013-01-09, Changed abs to fabs. temp: save x-sect profile (to determine tau)
					{
						cloudcrosssection->At( idxh,0 )	= heightsm.At( idxh);
						cloudcrosssection->At( idxh,1 )	= crosssection.At(idxl,idxh);
					}
				}
				db_xsectdata.erase();
			}
			else
			{
				nxLog::Record(NXLOG_ERROR,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::GetXSectsFromDBAndSetMembers, Error reading ice crystal from file <%s>", (const char*)db_xsectfile );
			}
		}
		else		
		{
			for (idxl=0;idxl<numLambdas;idxl++)			//TODO: extrapolate decent value for small particles.  for now, if reff outside range, set crosssection v small
			{
				crosssection.At(idxl,idxh) = ksca_background;		// TODO: fix this. 
			}
		}
	}
	ok = ok && SetCrossSectionTable( crosssection,lambdasnm,heightsm );		// set members
	ok = ok && WriteHeightWavelengthProfileToFile( (const char *)local_crosssectfile,crosssection,heightsm,lambdasnm );
	return ok;
}

/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::WriteHeightWavelengthProfileToFile		2010-3-4
 *
 *---------------------------------------------------------------------------*/
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::WriteHeightWavelengthProfileToFile( const char *filename,nx2dArray<double> crosssection, nx1dArray<double> heightsm,std::vector<double> lambdas ) const
{
	std::ofstream		strm;
	bool				ok;
	size_t				idxl,idxh,numLambdas,numHeights;

	numLambdas	= lambdas.size();
	numHeights	= heightsm.size();
	strm.open(filename, std::ios_base::out );
	strm << lambdas.size() << "\t" << heightsm.size() << std::endl;		// write the number of heights and the number of wavelengths
	for(idxl=0;idxl<numLambdas;idxl++)
	{
		strm << lambdas.at(idxl) << "\t"; 
	}
	strm << std::endl;
	for(idxh=0;idxh<numHeights;idxh++)
	{
		strm << heightsm.At(idxh) << "\t";
	}
	strm << std::endl;
	for(idxh=0;idxh<numHeights;idxh++)
	{
		for(idxl=0;idxl<numLambdas;idxl++)
		{
			strm << crosssection.At(idxl,idxh) << "\t";
		}
		strm << std::endl;
	}
	ok = !strm.fail();
	strm.close();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::WriteHeightWavelengthProfileToFile, error writing the height crosssection table to file <%s>", (const char*) filename);
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadHeightWavelengthProfileFromFile		2011-4-11*/
/** Loads in the tabulated crosssection from a file.
 *
 *	Line 1:					numwave	numheights			number of wavelengths and number of heights
 *	Line 2:					wavelengths(numwavelen)		expressed in nanometers.
 *	Line 3:					heights(numheights)			expressed in meters.
 *	Line 4:					cross section(numwave)		crosssection for wavelengths at height (0)
 *	Line 5:					cross section(numwave)		crosssection for wavelengths at height (1)
 *	..
 *	Line (numheights+3):	cross section(numwave)		crosssection for wavelengths at height (numheights-1)
**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadHeightWavelengthProfileFromFile( const char* filename )
{
	std::ifstream			strm;
	size_t					numh;
	size_t					numwave;
	nx1dArray<double>		wavelens;
	nx1dArray<double>		heights;
	nx2dArray<double>		crosssection;
	bool					ok;

	strm.open(filename, std::ios_base::in );
	strm  >> numwave >> numh;							// Load in the number of heights and the number of wavelengths
	ok = (numh > 0) && (numwave > 0) && !strm.fail();	// Make sure it is good
	ok = ok && wavelens.SetSize(numwave);
	ok = ok && heights.SetSize (numh );
	ok = ok && crosssection.SetSize( numwave, numh );
	if (ok)
	{
		strm >> wavelens; 
		strm >> heights;
		strm >> crosssection;
		ok = ok && !strm.fail();
	}
	strm.close();
	ok = ok && SetCrossSectionTable( crosssection, wavelens.STLVector(), heights);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadHeightWavelengthProfileFromFile, error loading and setting the height crosssection table from file <%s>", (const char*) filename);
		ReleaseResources();
	}
	return ok;
}
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetCrossSectionTable		2013-1-9*/
/** Set up the table that will be used for crosssection (or absorption) as a
 *	function of altitude and wavelength. The user passes in a 2-D array of crosssection
 *	values and two 1D arrays that define the wavelength and height "axes" of the grid.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetCrossSectionTable( const nx2dArray<double>& crosssection, const std::vector<double>& wavelen_nm, const nx1dArray<double>& heights_meters)
{
	const size_t*		dims;
	bool				ok;

	dims = crosssection.ArrayRankSpecs()->Dims();
	ok   =	   ( dims[0] == wavelen_nm.size()) 
			&& ( dims[1] == heights_meters.size())
			&& ( dims[0] > 0 )
			&& ( dims[1] > 0 );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetCrossSectionTable, Error setting table as array sizes are incompatible with each other");
		ReleaseResources();
	}
	else
	{
		ok =       m_heights.DeepCopy	( heights_meters);
		m_wavelennm.SetSize( wavelen_nm.size() );				// changed wavelen arrays passed to STL vects to work with ret. code
		for(size_t idx=0;idx<wavelen_nm.size();idx++)	
			m_wavelennm.At(idx)	= wavelen_nm.at(idx);
		//ok = ok && m_wavelennm.DeepCopy	( wavelen_nm );
		ok = ok && m_crossection.DeepCopy( crosssection );
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetCrossSectionTable, Error copying arrays over. I'm clearing array to zero");
		}
	}
	m_isdirty = true;
	if (!ok)
	{
		ReleaseResources();
	}
	return ok;
}
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetupCloudPropertiesGaussGeneric		2013-1-9*/
/***	Set up a cloud with a semi-Gaussian profile of cloud particle number density 
 *	and constant effective particle size.  'cloud thickness' is defined to be
 *	the FWHM and 'cloud top' is the upper half-way point.
 *	List the cloud properties in ascending altitude order.
 *	NOTE: optionally, make opt thicker by having (x-mu)^4, etc.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetupCloudPropertiesGaussGeneric( nx2dArray<double> *icedata,double hCloudTopKm,double cloudthickm,double particlesize,double maxnumdensity,double opticaldelta )
{
	size_t				idxh,numcheights,numsigma;
	nxString			datafile;
	double				sigmah,muh,h,num;

	muh			= hCloudTopKm*1e3 - cloudthickm/2.0;	// mean of num. density distrib
	sigmah		= cloudthickm/2.3548200;				// stdev of num density distrib (in m)
	numsigma	= 3;									// how many sigmas above/below define the 'cloud' region
	numcheights	= size_t(floor(double(2*numsigma)*sigmah/opticaldelta))+1;

	icedata->SetSize( numcheights,3 );
	for ( idxh=0;idxh<numcheights;idxh++ )
	{
		h		= opticaldelta*nxmath::round(((muh-double(numsigma)*sigmah) + double(idxh)*opticaldelta)/opticaldelta);
		icedata->At( idxh,0 ) = h;
		icedata->At( idxh,1 ) = particlesize;	
		num		= maxnumdensity*std::exp( -pow((h-muh),2.0)/(2*sigmah*sigmah) );
		icedata->At( idxh,2 ) = num;
	}
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetIceRadiusAndCrossSectionProfile		2013-1-9*/
/**	The primary setup function for this class.
 *	Configure the height- and wavelength-dependent cross sections and phase matrices 
 *	for the class using the specified particle size, and build the cross section profile.
 *	'configdir' specifies the location in which several needed working files 
 *	are stored: the master file (maps size(phase function file) to heights), 
 *	phase functions (spec at height, as funct. of wavelength and scatter angle),
 *	and the crosssection file (as a function of height and wavelength)
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetIceRadiusAndCrossSectionProfile( nx2dArray<double> icesizeprofile,std::vector<double> lambdasnm,nx2dArray<double> *cloudscatextinct,const char *configdir="./" )
{
	bool				ok;
	nxString			masterfilename;
	nx1dArray<double>	heightsm,reff,buffer;

	heightsm	= icesizeprofile.XSlice(0, &buffer);
	if ( heightsm.At(heightsm.size()-1)<1000 )
	{
		nxLog::Record(NXLOG_INFO,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetIceRadiusAndCrossSectionProfile, It looks like heights are in km. I am changing them to metres." );
		heightsm *= 1e3;
	}
	reff		= icesizeprofile.XSlice(1, &buffer);

	ok =		SetLocalDirectory		( configdir );								// set up the local dir. for tabulated class
	ok = ok &&	GetXSectsFromDBAndSetMembers		( heightsm,reff,lambdasnm,cloudscatextinct );
	ok = ok &&	GetPhaseFunctsFromDBAndSetMembers	( heightsm,reff,lambdasnm,masterfilename );	// get phase functions from DB, truncate, calc. forward-scat fraction
	ok = ok &&	LoadHeightWavelengthProfileFromMasterFile( masterfilename );		//write to indiv. files, write file list
	
	if (!ok)
	{
		nxLog::Record(NXLOG_ERROR,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetIceRadiusAndCrossSectionProfile, error setting up ice optical properties from database" );
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetIceRadiusAndNumDensityProfileGivenTau		2013-1-9*/
/**	The primary setup function for class IF the user desires to treat the cross-sections 
 *	as extinctions (i.e. fold num density into cross-section and use a 'one' climatology 
 *	for the species.
 *	Configure the height- and wavelength-dependent extinctions and phase matrices 
 *	for the class using the specified particle size, and build the extinction profile 
 *	such that the vertical optical depth at wavelength 'lambdatau' is 'tau' for 
 *	the given optical property spacing 'opticaldelta'.
 *	'configdir' specifies the location in which several needed working files 
 *	are stored: the master file (maps size(phase function file) to heights), 
 *	phase functions (spec at height, as funct. of wavelength and scatter angle),
 *	and the extinction file (as a function of height and wavelength)
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetIceRadiusAndNumDensityProfileGivenTau( 
						std::vector<double> lambdasnm,double hCloudTopKm,double cloudthickM,double opticaldelta,
						double particlesize,double lambdatau,double tau,
						nx2dArray<double> *iceprofile,nx2dArray<double> *cloudscatextinct,std::vector<double> *raytracealts,double minraytracedelta,double *iwp,const char *configdir="./" )
{
	bool				ok(true);
	nxString			masterfilename;
	nx1dArray<double>	heightsm,reff,buffer,numden;
	double				maxnumdensity;
	size_t				idx;

	// set up the effective-size and number denstiy profiles
	ok =		SetupCloudPropertiesGaussGeneric( iceprofile,hCloudTopKm,cloudthickM,particlesize,maxnumdensity=1.0,opticaldelta );
	heightsm	= iceprofile->XSlice(0, &buffer);
	if ( heightsm.At(heightsm.size()-1)<1000 )
	{
		nxLog::Record(NXLOG_INFO,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetIceRadiusAndNumDensityProfileGivenTau, It looks like heights are in km. I am changing them to metres." );
		heightsm *= 1e3;
	}
	reff		= iceprofile->XSlice(1, &buffer);				//
	numden		= iceprofile->XSlice(2, &buffer);				// generic gaussian #den profile

	ok = ok &&	SetLocalDirectory					( configdir );				// set up the local dir. for tabulated class
	ok = ok &&	GetXSectsFromDBAndSetMembersGivenTau( heightsm,reff,lambdasnm,&numden,cloudscatextinct,lambdatau,tau,raytracealts,minraytracedelta,iwp );
	for(idx=0;idx<numden.size();idx++)
	{
		iceprofile->At( idx,2 )	= numden.At( idx );								// copy the numdens that make 'tau'
	}
	if (ok)
	{
		ok = ok &&	GetPhaseFunctsFromDBAndSetMembers	( heightsm,reff,lambdasnm,masterfilename );	// get phase functions from DB, truncate, calc. forward-scat fraction
		ok = ok &&	LoadHeightWavelengthProfileFromMasterFile( masterfilename );		//write to indiv. files, write file list
	}
	else
	{
		nxLog::Record(NXLOG_ERROR,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::SetIceRadiusAndNumDensityProfileGivenTau, error setting up ice optical properties from database" );
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *	skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::GetXSectsFromDBAndSetMembersGivenTau	2011-11-23
 *	Given the wavelengths and effective sizes and number density profiles, get 
 *	the scattering cross sections from the Baum DB, convert to extinctions, 
 *	and populate the member extinction array.  
 *	The number density profile (and consequently extinction @ the retrieval wavelength)
 *	is shifted so that the cloud layer's total optical thickness equals the specified value. 
 *	If the specified optical depth can't be done (too thick for 10m minimum 
 *	cell spacing), an error is returned and the computation is not done.
 *---------------------------------------------------------------------------*/
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::GetXSectsFromDBAndSetMembersGivenTau( const nx1dArray<double>& heightsm,const nx1dArray<double>& reff,const std::vector<double>& lambdasnm,nx1dArray<double> *numden,nx2dArray<double> *cloudscatextinct,double lambdatau,double tau,std::vector<double> *raytracealts,double minraytracedelta,double *iwp )
{	
	bool				ok;
	nxString			db_xsectfile,local_extinctfile;
	size_t				idxl,idxh,idxx;
	size_t				idxtau =0;
	size_t				numLambdas,numHeights;
	double				ksca_background(1e-30),opticaldelta,maxtau,f,tauthresh;
	nx2dArray<double>	extinction;
	nx2dArray<double>	db_xsectdata;
	std::vector<double>	db_lambdas;
	nx1dArray<double>	fraction(18,trunc_fraction);
	nx1dArray<double>	iwc(18,baum_iwc);

	local_extinctfile.sprintf( "%s/IceCrossSections.dat",(const char *)m_configdir );
	db_lambdas	= m_wavelennm_db.STLVector();
	numLambdas	= lambdasnm.size();
	numHeights	= reff.size();
	m_heights.SetSize( numHeights );
	ksca_background		= 1e-30;						// dummy value for lambdas outside database
	tauthresh			= 1.4;							// max cell optical depth threshold
	opticaldelta	= heightsm.At(1)-heightsm.At(0);	// optical props delta specified
	ok =		extinction.SetSize( lambdasnm.size(),heightsm.size() );
	ok = ok &&	cloudscatextinct->SetSize( heightsm.size(),4 );
	for (idxh=0;idxh<numHeights;idxh++)		
	{
		if ( IsValidDatabaseEffectiveSize(reff.At(idxh)) )	// if reff is listed size, get cross sect from DB
		{
			db_xsectfile.sprintf( "%s/Solar_mix1_d%03u.dat",(const char *)LoadDatabaseDirectoryFromRegistry(),(unsigned int)(reff.At(idxh)) );
			ok		= db_xsectdata.InputColumnMajorText( (const char *)db_xsectfile,m_numfilecolumns );
			if (ok)
			{
				for (idxl=0;idxl<numLambdas;idxl++)			// for each model wavelength, find
				{											// DB equivalent lambda
					idxx	= (size_t)(std::lower_bound( db_lambdas.begin(),db_lambdas.end(),lambdasnm.at(idxl) )-db_lambdas.begin());
					if (idxx >= 0)							// if wavelength valid, get scat. x-sect
					{										// and scale by generic num. density
						extinction.At(idxl,idxh) = db_xsectdata.At(idxx,m_numfilecolumns-2)*1.0e-8*numden->At(idxh);// convert x-sect in um^2 to cm^2
					}
					else									// otherwise 
					{
						extinction.At(idxl,idxh) = ksca_background;	// use insignif. value
					}
					if (fabs(lambdasnm.at(idxl)-lambdatau)<0.75)	// ndl303, 2013-01-09, changed abs to fabs. temporary: save extinct profile to file (for tau)
					{
						idxtau	= idxl;
					}
				}
				db_xsectdata.erase();
			}
			else
			{
				nxLog::Record(NXLOG_ERROR,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::GetXSectsFromDBAndSetMembersGivenTau, Error reading ice crystal from file <%s>", (const char*)db_xsectfile );
			}
		}
		else											// if invalid wavelength (<400)
		{
			for (idxl=0;idxl<numLambdas;idxl++)			//TODO: extrapolate decent value for small particles.  for now, if reff outside range, set extinction v small
			{
				extinction.At(idxl,idxh) = ksca_background;		// TODO: fix this. 
			}
		}
	}
	ok = ok && AdjustNumAndExtinctionToFitTau( heightsm,&extinction,numden,idxtau,tau,opticaldelta );// adjust ext,numden to specified tau
	maxtau	= 0;
	*iwp	= 0;
	for (idxh=0;idxh<numHeights;idxh++)		
	{
		f		 = fraction.At(size_t(reff.At( idxh )/10.0)-1);					// get (default) trunc. fraction
		*iwp	+= numden->At(idxh)*iwc.At(size_t(reff.At( idxh )/10.0)-1)*(opticaldelta);	// integrate IWC -> IWP (g/m2)
		cloudscatextinct->At( idxh,0 )	= heightsm.At( idxh);
		cloudscatextinct->At( idxh,1 )	= extinction.At(idxtau,idxh)*1e5;		// ext. per km
		cloudscatextinct->At( idxh,2 )	= (1-f)*extinction.At(idxtau,idxh)*1e2*(opticaldelta);// (vert) segment opt depth
		cloudscatextinct->At( idxh,3 )	= (1-f)*extinction.At(idxtau,idxh)*1e2*sqrt((opticaldelta)*(2.0*(6378.0e3+heightsm.At(idxh))+opticaldelta) );//max along-path segment opt depth
		maxtau	= nxmax( maxtau,cloudscatextinct->At( idxh,3 ) );					// max along-path tau per shell:
	}																			// NOTE: this uses 'opticaldelta' in place of 'raytracedelta'
	ok = (maxtau < tauthresh);						// ensure scattering opt depth doesn't get too large
	if (ok)
	{
		ok = ok && ConfigureRayTracingByExtinctionProfile( cloudscatextinct,minraytracedelta,raytracealts );
		ok = ok && SetCrossSectionTable( extinction,lambdasnm,heightsm );		// set members
		ok = ok && WriteHeightWavelengthProfileToFile( (const char *)local_extinctfile,extinction,heightsm,lambdasnm );
	}
	else
	{
		nxLog::Record(NXLOG_ERROR,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::GetXSectsFromDBAndSetMembersGivenTau, optical depth %.4f too large for this cloud thickness.",tau );
	}
	return ok;
}
/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::ConfigureRayTracingByExtinctionProfile		2013-1-9*/
/**	Set up the SASKTRAN configuration altitudes (ray tracing alts and diffuse 
 *	points) according to the specified extinction profile, usually assuming a 
 *	Gaussian numden profile.  the goal here is to maintain a constant cell optical 
 *	depth.
 *	The maximum cell optical depth is calculated from the point of the max. number
 *	density of the profile, and this per-cell depth is used to scale the ray-tracing
 *	shells outward from this maximum.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::ConfigureRayTracingByExtinctionProfile( nx2dArray<double> *cloudscatextinct,double minraytracedelta,std::vector<double> *raytracealts )
{
	bool					ok(true);
	size_t					idx,maxidx,ctr,idxh;
	double					tau0,deltah,gridDelta,gridCloudBottTop;
	std::vector<double>		shellaltsc,hSTL,shell0,shellf;
	nx1dArray<double>		shellalts0,shellaltsf;
	nx1dArray<double>		h,extinct,buffer;
	double					htop,hbott,outsidedelta,minshellabove,maxshellbelow,minshell,maxshell,maxdelta;
	size_t					numshellsbelow,numshellsabove;
	
	h		= cloudscatextinct->XSlice(0, &buffer);
	extinct	= cloudscatextinct->XSlice(1, &buffer);	

	gridDelta		= 5.0;										// shell spacings must differ by multiple of this
	gridCloudBottTop= 100.0;									// cloud top and bottom 'snapped' to 100m accuracy
	outsidedelta	= 1000.0;									// nominal 1km shell spacing outside

	hSTL		= h.STLVector();
	maxidx		= size_t(floor((h.size()-1)/2.0));				// idx of max extintion (gauss profile)
	tau0		= extinct.At( maxidx )*minraytracedelta;		// opt. depth of max cell: stay ~const over profile
	shellaltsc.resize( h.size()*size_t(minraytracedelta/(h.At(1)-h.At(0))) );
	ctr = 0;
	shellaltsc.at(ctr)	= h.At( maxidx );						// altitute of max extinction
	deltah				= gridDelta*nxmath::round(tau0/extinct.At( maxidx-1)/gridDelta );	// round these to nearest 'delta'
	shellaltsc.at(ctr+1)= shellaltsc.at(ctr) - deltah;			//

	ctr++;
	maxdelta			= deltah;
	while(shellaltsc.at(ctr)>h.At(0) )							// go until get to bottom of optical alts
	{
		idxh		= size_t(std::lower_bound( hSTL.begin(),hSTL.end(),shellaltsc.at(ctr)-deltah )- hSTL.begin());	// index in extinct array
		deltah		= gridDelta*nxmath::round( tau0/extinct.At( idxh )/gridDelta );	// maintain extinct/km ratio
		if (deltah >=10.0 && deltah < 1000.0)					// chack that spacing does not get too large
		{
			shellaltsc.at(ctr+1)	= shellaltsc.at(ctr)-deltah;	// round to nearest grid
		}
		else														// if out of bounds (is too large)
		{															// use prev. spacing
			shellaltsc.at(ctr+1)	= shellaltsc.at(ctr)-maxdelta;
		}
		maxdelta	= nxmax(maxdelta,deltah);
		ctr++;
	}
	ctr++;
	shellaltsc.resize( 2*ctr-1 );
	for( idx=0;idx<ctr-1;idx++)					// reflect heights above max extinct alt
	{
		shellaltsc.at(idx+ctr)	= 2*shellaltsc.at(0) - shellaltsc.at(idx+1);
	}
	std::sort( shellaltsc.begin(),shellaltsc.end() );
	std::copy(shellaltsc.begin(), unique(shellaltsc.begin(),shellaltsc.end() ),shellaltsc.begin() );

	htop		= gridCloudBottTop*nxmath::round( shellaltsc.back()/gridCloudBottTop );	// cloud top and bottom
	hbott		= gridCloudBottTop*nxmath::round( shellaltsc.front()/gridCloudBottTop );// to nearest 100m
	// ray-tracing shells config
	if ((size_t(htop)%size_t(outsidedelta))==0)			// if cloud top on shell boundary,
	{
		minshellabove	= htop + outsidedelta;			// shells start one up from here
	}
	else												// otherwise
	{
		minshellabove	= outsidedelta*ceil(htop/outsidedelta);	// shells start at next 'delta'
	}
	if ((size_t(hbott)%size_t(outsidedelta))==0)					// ditto with cloud bottom
	{
		maxshellbelow	= hbott - outsidedelta;			
	}
	else
	{
		maxshellbelow	= outsidedelta*floor(hbott/outsidedelta);
	}
	minshell	= 0.0;
	maxshell	= 100000.0;
	numshellsbelow		= size_t((maxshellbelow-minshell)/outsidedelta) + 1;	// shells below, including boundary
	numshellsabove		= size_t((maxshell-minshellabove)/outsidedelta) + 1;	// 
	shellalts0.SetSize(numshellsbelow);									// 
	shellalts0.Indgen(numshellsbelow)	*= outsidedelta;				
	(shellalts0)						+= minshell;
	shellaltsf.SetSize(numshellsabove);									
	shellaltsf.Indgen(numshellsabove)	*= outsidedelta;				
	(shellaltsf)						+= minshellabove;
	shell0	= shellalts0.STLVector();
	shellf	= shellaltsf.STLVector();
	raytracealts->clear();
	raytracealts->reserve( shell0.size()+shellaltsc.size()+shellf.size() );
	raytracealts->insert( raytracealts->end(),shell0.begin(),shell0.end() );
	raytracealts->insert( raytracealts->end(),shellaltsc.begin(),shellaltsc.end() );
	raytracealts->insert( raytracealts->end(),shellf.begin(),shellf.end() );
	return ok;
}

/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::AdjustNumAndExtinctionToFitTau		2011-11-23
 *	Adjust the number denstiy and extinction array to make the vertical optical depth 
 *	at the specified wavelength (idxtau) match the specified value (tau)
 *---------------------------------------------------------------------------*/
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::AdjustNumAndExtinctionToFitTau( nx1dArray<double> heightsM,nx2dArray<double> *extinction,nx1dArray<double> *numden,size_t idxtau,double tau,double /*opticaldelta*/ )
{
	bool					ok(true);
	double					tau1;
	nx1dArray<double>		extinct,buffer,integral;
	nxSpline				spline;

	extinct		= (*extinction).YSlice( int(idxtau), &buffer);	// get extinction profile @ ret. wavelen,
	spline.Configure( heightsM,extinct );
	tau1		= spline.Integrate( &integral )*(1e2);			// convert extinct to m-1, and integrate it

	if (fabs(tau1) > 1e-30)
	{
		(*numden)		*= (tau/tau1);							// adjust extinction and numden
		(*extinction)	*= (tau/tau1);							// to adjust opt depth equal to 'tau'
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::GetPhaseFunctsFromDBAndSetMembers		2010-1-18
 *	For given reff(h) and lambdas, get the phase matrix data from the Baum cirrus
 *	properties database.
 *	This function then truncates the phase functions, writes them to file, and 
 *	calls the 'master file'-writing function.
 *---------------------------------------------------------------------------*/
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::GetPhaseFunctsFromDBAndSetMembers( const nx1dArray<double>& heightsm,const nx1dArray<double>& reff,const std::vector<double>& lambdasnm,nxString &masterfilename )
{
	bool					ok(true);
	std::ofstream			strm;
	nx1dArray<nxString>		filenames;
	nx1dArray<double>		cosscatangle,p;
	nx2dArray<double>		db_phasematrixdata;
	nx3dArray<double>		phasefunct;
	std::vector<double>		effectivesizes;
	nxString				db_phasematrixfname;
	size_t					idxr,idxl,idxt,idxh,idxcol,numheights,numlambdas,numangles;

	masterfilename.sprintf( "%s/MasterFile.dat",(const char *)m_configdir );
	numlambdas		= lambdasnm.size();
	numheights		= reff.size();
	numangles		= m_scatterangles_db.size();
	effectivesizes	= m_effectivesize_db.STLVector();
	ok =		filenames.SetSize( numheights );
	ok = ok &&	cosscatangle.SetSize( numangles );
	ok = ok &&	phasefunct.SetSize( numlambdas,numangles,numheights );
	ok = ok &&	m_deltafraction.SetSize( numlambdas,numheights );
	m_deltafraction.SetTo( 0.0 );
	for (idxt=0;idxt<numangles;idxt++)			// copy cos(scatterangle)
	{
		cosscatangle.At(idxt)		= nxmath::cosd( m_scatterangles_db.At(idxt) );
	}
	for (idxl=0;idxl<numlambdas;idxl++)					// for every lambda, get P11 for each effective size
	{
		if (IsValidDatabaseWavelength( lambdasnm.at(idxl)) )	// if lambda is listed size, get phase funct from DB
		{
			db_phasematrixfname.sprintf( "%s/Solar_mix1_spf_col_%1.2f.dat", (const char *)LoadDatabaseDirectoryFromRegistry(),lambdasnm.at(idxl)/1000.0 );
			ok = db_phasematrixdata.InputColumnMajorText( (const char *)db_phasematrixfname,2*m_effectivesize_db.size()+1 );
			if (ok)
			{
				for (idxh=0;idxh<numheights;idxh++)			// for every height, get reff.  get phase function for this lambda, reff
				{
					if (IsValidDatabaseEffectiveSize(reff.At(idxh)))
					{
						idxr	= (size_t)(std::lower_bound( effectivesizes.begin(),effectivesizes.end(),reff.At(idxh))-effectivesizes.begin());
						for (idxt=0;idxt<numangles;idxt++)				
						{
							idxcol							= 2*idxr+1;					// DB includes mean and stdev: offset to correct
							phasefunct.At(idxl,idxt,idxh)	= db_phasematrixdata.At(idxt,idxcol);
						}
						// get phase function and truncate it for given scattering angle
						phasefunct.Slice( (int)idxl,(int)idxl,NXARRAY_STARSELECT,NXARRAY_STARSELECT,(int)idxh,(int)idxh,&p );	// get this phase function
						m_deltafraction.At( idxl,idxh )	= TruncateAndComputeDeltaFraction( &p,m_truncangles_db.At(idxr) );	// replace with truncated version and calculate 'f'
					}
					else														// small size (used for bracketing heights): set isotropic
					{
						for (idxt=0;idxt<numangles;idxt++)			//TODO: make this better.  for now, get phase funct for min effective size
						{
							phasefunct.At(idxl,idxt,idxh)	= 1.0;
						}
					}
				}
			}
			else
			{
				nxLog::Record(NXLOG_ERROR,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::Get_PhaseFunctFromDB, Error reading from file <%s>", (const char*)db_phasematrixfname );
			}
		}
		else
		{
			nxLog::Record(NXLOG_INFO,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::Get_PhaseFunctFromDB, database phase functions not defined for wavelength <%3f> nm", lambdasnm.at(idxl) );
			for (idxt=0;idxt<numangles;idxt++)	
			{
				for (idxh=0;idxh<numheights;idxh++)			
				{
					phasefunct.At(idxl,idxt,idxh)	= 1.0;
				}
			}
		}
	}
	ok	= ok && WritePhaseFunctionFiles( lambdasnm,reff,cosscatangle,phasefunct,filenames );
	ok	= ok && WriteMasterFile( heightsm,filenames,masterfilename );

	db_phasematrixdata.erase();
	cosscatangle.erase();				
	phasefunct.erase();
	filenames.erase();
	return ok;
}

/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::TruncateAndComputeDeltaFraction		2011-7-20
 *	Truncate phase function using a 'cutoff angle' criterion and return the directly
 *	forward-scattered fraction, 'f'.
 *	NOTE: this function overwrites the original phase function with its truncated
 *	version.
 *---------------------------------------------------------------------------*/
double skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::TruncateAndComputeDeltaFraction( nx1dArray<double> *phasefunct,double thetaC )
{
	double				deltafraction;		// fraction of incoming light not forward-scattered
//	bool				ok(true);
	size_t				idx,idxr,idxcut;
	nx1dArray<double>	xmu,pdelta,integral;
	std::vector<double>	thetas_stl;
	nxSpline			spline;

	thetas_stl		= m_scatterangles_db.STLVector();			// get index at which 'thetaC' occurs
	idxcut			= size_t(std::lower_bound( thetas_stl.begin(),thetas_stl.end(),thetaC )-thetas_stl.begin());
	xmu.SetSize( idxcut );
	pdelta.SetSize( idxcut );
	for(idx=0;idx<idxcut;idx++)							// get cos(theta) and delta-funct for angles in fwd peak
	{													// put in order of incr. cos(theat)
		idxr				= idxcut-idx-1;
		xmu.At(idx)			= nxmath::cosd( m_scatterangles_db.At( idxr ) );
		pdelta.At(idx)		= phasefunct->At( idxr ) - phasefunct->At( idxcut );	// delta-funct component
		phasefunct->At(idxr)= phasefunct->At( idxcut );	// replace orig. with truncated version
	}
	spline.Configure( xmu,pdelta );
	deltafraction	= 0.5*spline.Integrate( &integral );
	return deltafraction;
}

/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::ComputeAsymmetryParameter		2011-9-7
 *	Compute the asymmetry parameter g = <cos(Theta)> for the indicated phase function.
 *---------------------------------------------------------------------------*/
double skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::ComputeAsymmetryParameter( nx1dArray<double> *phasefunct,double /*thetaC*/ )
{
	double				asymmetry;
	nxSpline			spline;
	size_t				idx,idxr,numangles;
	nx1dArray<double>	xmu,function,integral;
	numangles	= phasefunct->size();
	function.SetSize( numangles );
	xmu.SetSize( numangles );
	for(idx=0;idx<numangles;idx++)
	{									
		idxr				= numangles-idx-1;
		xmu.At(idx)			= nxmath::cosd( m_scatterangles_db.At( idxr ) );
		function.At(idx)	= phasefunct->At( idxr )*xmu.At(idx);	
	}
	spline.Configure( xmu,function );
	asymmetry	= 0.5*spline.Integrate( &integral );
	return asymmetry;
}



/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::WritePhaseFunctionFile		2010-1-6
 *	Write the file ...
 *	for each height, write P11(mu) into format for skOpticalProperties_TabulatedPhaseMatrix_HeightEntry
 *
 *---------------------------------------------------------------------------*/
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::WritePhaseFunctionFiles(const std::vector<double>& lambdasnm,const nx1dArray<double>& reff,const nx1dArray<double>& cosscatangle,const nx3dArray<double>& phasefunct,nx1dArray<nxString>& filenames ) const
{
	bool				ok = true;
	std::ofstream		strm;
	size_t				idxl,idxh,idxt,numAngles;
	nxString			cachefilename,fullcachefilename;

	numAngles	= cosscatangle.size();

	for (idxh=0;idxh<filenames.size();idxh++)				
	{
		cachefilename.sprintf( "PhaseFunct_%03.1f.dat",reff.At(idxh) );	//TODO: make better caching string format for these?
		fullcachefilename.sprintf( "%s/%s",(const char *)m_configdir,(const char *)cachefilename );
		filenames.At(idxh)	= cachefilename;		// enter full path to each file (used in master file)

		strm.open((const char *)fullcachefilename, std::ios_base::out );	ok = strm.is_open();
		strm << lambdasnm.size() << "\t" << cosscatangle.size() << std::endl;	ok = !strm.fail();				// write number of wavelengths, make sure the stream is ok
		for (idxt=0;idxt<cosscatangle.size();idxt++)			// invert order so that mu=cos(theta) is in asc order
		{
			strm << std::fixed << std::setprecision(12) << cosscatangle.At(numAngles-1-idxt) << "\t";			
			ok = ok && !strm.fail();		// write in the cos(scatter angle) from the file
		}
		strm << std::endl;
		for (idxl=0;idxl<lambdasnm.size();idxl++ )			// invert order so that P11(cos(theta)) in asc order of mu
		{
			strm << lambdasnm.at(idxl) << "\t";
			for (idxt=0;idxt<cosscatangle.size();idxt++)
			{
				strm << phasefunct.At(idxl,numAngles-1-idxt,idxh) << "\t";
			}
			strm << std::endl;
		}
		strm.close();
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::WritePhaseFunctionFiles, Error writing phase function file <%s>", (const char*)fullcachefilename );
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::WriteMasterFile			2010-1-6
 * Write the file ...
 *	write master file for skOpticalProperties_TabulatedPhaseMatrix_HeightWavelength
 *
 *---------------------------------------------------------------------------*/
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::WriteMasterFile( const nx1dArray<double>& heightsm,const nx1dArray<nxString>& filenames,const nxString& masterfilename ) const
{
	bool				ok;
	std::ofstream		strm;

	strm.open((const char *)masterfilename, std::ios_base::out );
	ok =		strm.is_open();
	if (ok)
	{
		ok = (heightsm.size() == filenames.size());
		if (ok)
		{
			for(size_t idxh=0;idxh<filenames.size();idxh++)
			{
				strm << heightsm.At(idxh) << "\t" << filenames.At(idxh) << std::endl;							
				ok = ok && !strm.fail();								// MAke sure the stream is ok
			}
			strm.close();
		}
		else
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::WriteMasterFile, number of heights is different from number of phase function files." );
		}
	}
	else
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::WriteMasterFile, Error writing master file <%s>", (const char*)masterfilename );
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *		skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadHeightWavelengthProfileFromMasterFile		2011-4-11*/
/** Loads in the set of phasematrices asa function of height, wavelength and scattering angle.
 *	This method is responsible for managing the entire load and it controls
 *	loading the height profiles. This consists of reading 
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadHeightWavelengthProfileFromMasterFile( const char* masterfilename )
{
	std::ifstream										strm;
	double												h;
	double												lasth = -9999999.0;
	char												buffer[1000];
	nxString											filename;
	nxString											fullfilename;
	nxFileSpec											spec(masterfilename);
	bool												ok=false;
	skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry		dummy;
	skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry*	entry;

	//ReleaseResources();		//TODO: 

	strm.open(masterfilename, std::ios_base::in );
	while ((!strm.eof() && !strm.fail()) )
	{
		strm >> h;										// Input the height
		if (!strm.eof() && !strm.fail())				// IF we have not hit the end of file or had an error
		{												// then
			ok = (lasth < h);								
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "LoadHeightWavelengthProfileFromMasterFile, the height profiles are not specified in ascending order, this is going to mess up future height indexing and interpolation");
			}
			else
			{
				lasth = h;
				strm.getline( buffer, 1000);
				filename = buffer;
				filename.RemoveWhiteSpace();
				fullfilename = spec.FullDirSpec() + filename;

				#if defined(NX_WINDOWS)
					char fullexpandedname[_MAX_PATH+1];
					if( _fullpath( fullexpandedname,fullfilename, _MAX_PATH ) != NULL ) fullfilename = fullexpandedname;
				#else
					NXTRACE_ONCEONLY(firsttime, ("skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadHeightWavelengthProfileFromMasterFile, Got to check that filenames are properly resolved in the Lunix code\n"));
				#endif

				ok = nxDirectory::FileExists(fullfilename);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_BaumDB::LoadHeightWavelengthProfileFromMasterFile, Cannot load in phase matrix data for height %g as the file <%s> does not exist", (double)h, (const char*)fullfilename);
				}
				else
				{
					m_heightentries.push_back(dummy);
					entry = &m_heightentries.back();
					ok = entry->LoadWavelengthEntriesFromHeightFile( h, fullfilename );
					if (!ok) m_heightentries.pop_back();
				}
			}
		}
	}
	strm.close();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry		2011-4-16*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry()
{
	m_roundedwavelen = 0;
	//m_roundedwavelen = 0.0;
	m_phasematrix.SetStatic();
	m_cosscatterangle.SetStatic();
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::ConfigureEntry		2011-4-17*/
/** Configures this entry. Sets up the wavelength by rounding. Copies over the 
 *	phasematrix and cos(scattering angle) arrays. Checks the cos(scattering angle) 
 *	to make sure it is in proper range and order.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::ConfigureEntry( double wavelennm, const double* phasematrix, const double* cosangle, size_t numpoints)
{
	bool		ok;
	size_t		idx;

	m_roundedwavelen = RoundTheWavelength(wavelennm);

	ok = (cosangle[0] >= -1.00001) && (cosangle[0] <= 1.00001);						// Make sure the cosine of scattering angle is in proper range. Check first element
	for (idx = 1; idx < numpoints; idx++ )											// also check that it is in ascending order (this is easy for the user to get wrong)
	{																				// as scattering anlgle 0 to 180 is opposite to cos(scatter angle)
		ok = ok && (cosangle[idx] >= -1.00001) && (cosangle[idx] <= 1.00001) && (cosangle[idx] > cosangle[idx-1]);
	}
	ok = ok && m_phasematrix.CopyGridArray    ( phasematrix, numpoints );			// Create the phase matrix grid
	ok = ok &&  m_cosscatterangle.CopyGridArray( cosangle,    numpoints );			// and the scattering angle grid. They must be the same size.

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::ConfigureEntry, The entry is invalid as the cos(scattering angles) are not in ascending order (-1 to +1) or are out of the range (-1 to +1) or there was a memory allocation error.");
		ReleaseResources();
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::GetPhaseMatrix		2011-4-16*/
/**	Returns the scalar phasematrix using linear interpolation of the
 *	scattering angle grid. The 
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::GetPhaseMatrix( double cosangle, double* phasematrix ) const
{
	SKTRAN_GridIndex	lowercell;
	SKTRAN_GridIndex	uppercell;
	double				lowerweight;
	double				upperweight;
	bool				ok;
	
	*phasematrix = 0;
	ok = (cosangle >= -1.00001) && (cosangle <= 1.00001);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::GetPhaseMatrix, The cosine(scatteringing) is not in the range -1 to 1. Make sure you are using cosine(scattering angle). Do OT use radians or degree!");
	}
	else
	{
		ok = m_cosscatterangle.FindBoundingIndices( cosangle, SKTRAN_GridDefBase_V2::OUTOFBOUND_ERROR, &lowercell, &lowerweight, &uppercell, &upperweight);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::PhaseMatrix, There was an error looking up phase matrix for cos(scatter angle)  = %g, setting value to 0.0", (double)cosangle);
		}
		else
		{
			*phasematrix = lowerweight* m_phasematrix.At(lowercell) + upperweight*m_phasematrix.At(uppercell);
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::ReleaseResources		2011-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::ReleaseResources()
{
	m_phasematrix.erase();
	m_cosscatterangle.erase();
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::DeepCopy		2011-4-17*/
/** Make a cloned copy of the other object.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::DeepCopy( const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry& other )
{
	bool	ok;

	m_roundedwavelen = other.m_roundedwavelen;
	ok =       m_phasematrix.DeepCopy( other.m_phasematrix );
	ok = ok && m_cosscatterangle.DeepCopy( other.m_cosscatterangle );
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::DeepCopy, Error copying wavelength entries, thats a problem ");
		ReleaseResources();
	}
	return ok;
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry		2011-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry	()
{
	m_heightmeters = -999999.0;
	m_numwavelens = 0;
	m_maxwavelens = 0;
	m_wavelenentries = NULL;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::skOpticalProperties_TabulatedPhaseMatrix_HeightEntry		2011-4-26*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry( double heightm )
{
	m_heightmeters = heightm;
	m_numwavelens = 0;
	m_maxwavelens = 0;
	m_wavelenentries = NULL;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::~skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry		2011-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::~skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedPhaseMatrix_HeightEntry::ReleaseResources		2011-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::ReleaseResources()
{
	if (m_wavelenentries != NULL)
	{
		delete [] m_wavelenentries;
		m_wavelenentries = NULL;
	}
	m_numwavelens = 0;
	m_maxwavelens =0;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::AllocateEntries		2011-4-17*/
/** Allocates N entries. Uses the small optimation os using the same memory if the
 *	the current emory allocation is sufficient.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::AllocateEntries	( size_t maxentries )
{
	bool	ok;

	ok = (maxentries <= m_maxwavelens);
	if (!ok)
	{
		ReleaseResources();
		m_wavelenentries = new skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry[maxentries];
		ok = (m_wavelenentries != NULL);
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry, Error allocating memeory for %u wavelength entries", (unsigned int)maxentries );
		}
		else
		{
			m_maxwavelens  = maxentries;
		}
	}
	if (ok) m_numwavelens = maxentries;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::LoadWavelengthEntriesFromHeightFile		2011-4-17*/
/** Reads in the phase matrix entries for a set of wavelengths from a
 *	single text file associated with one altitude. The description of the file
 *	format can be found in the class description.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::LoadWavelengthEntriesFromHeightFile( double heightm, const char* filename )
{
	std::ifstream		strm;
	bool				ok;
	bool				ok1;
	size_t				numwavelengths;
	size_t				numangles;
	double				wavelen_nm;
	nx1dArray<double>	phasematrix;
	nx1dArray<double>	cosangles;
	size_t				i;
	

	m_heightmeters = heightm;
	strm.open(filename, std::ios_base::in );
	ok = strm.is_open();

	strm >> numwavelengths >> numangles;						// Read in the number of wavelengths
	ok = !strm.fail();											// MAke sure the stream is ok
	ok = ok && AllocateEntries( numwavelengths );				// Allocate space for the wavelength entries
	ok = ok && phasematrix.SetSize (numangles);					// Allocate a temporary buffer for holding the phasematrix
	ok = ok && cosangles.SetSize(numangles );					// and another for holding the cos(scatter angle)
	
	strm >> cosangles; ok = ok && !strm.fail();					// read in the cos(scatter angle) from the file

	for (i = 0; i < numwavelengths; i++ )
	{
		strm >> wavelen_nm;
		strm >> phasematrix;
		ok1 = !strm.fail();
		ok1 = ok1 && m_wavelenentries[i].ConfigureEntry( wavelen_nm, phasematrix.UnsafeArrayBasePtr(), cosangles.UnsafeArrayBasePtr(), numangles );
		if (!ok1)
		{
			nxLog::Record(NXLOG_WARNING, "skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::LoadWavelengthEntriesFromHeightFile, There was an error loading phase-matrix entries for the wavelength entry %u, wavelength = %g", (unsigned int)i, (double)wavelen_nm);
		}
		ok = ok && ok1;
	}
	strm.close();
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::LoadWavelengthEntriesFromHeightFile, There was an error loading the pahse matrix tables from file <%s>, I'm deleting the  whole structuire form memory");
		ReleaseResources();
	}
	return ok;

}

/*-----------------------------------------------------------------------------
 *			skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::DeepCopy		2011-4-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::DeepCopy( const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry& other )
{
	bool	ok;

	m_heightmeters = other.Height();
	ok = AllocateEntries( other.m_numwavelens );
	for (size_t idx = 0; idx < other.m_numwavelens; idx++)
	{
		ok = ok && m_wavelenentries[idx].DeepCopy( other.m_wavelenentries[idx] );
	}
	m_numwavelens = other.m_numwavelens;
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::DeepCopy, Error copying one intsance to another. Thats a problem");
		ReleaseResources();
		m_heightmeters = -9999999.0;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					class entries_equal		2013-1-9*/
/** \internal 
 *	A quick and dirty comparator Predicate class used by the stl find_if in
 *	method FindWavelengthEntry
 **/
/*---------------------------------------------------------------------------*/

class entries_equal
{
	private:
		//double	m_wavelen;
		size_t	m_wavelen;

	public:
					entries_equal	(size_t roundedwavelen)														{ m_wavelen = roundedwavelen;}
		bool		operator		() (const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry& entry) const { return m_wavelen == entry.WavelengthIndex();}
};


/*-----------------------------------------------------------------------------
 *			skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::FindWavelengthEntry		2011-5-7*/
/** Find teh wavelength entry at this altitudde**/
/*---------------------------------------------------------------------------*/
bool skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::FindWavelengthEntry( double wavelennm, const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry** entry ) const
{
	//double														roundednm;
	size_t														roundednm;
	const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry*	start;
	const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry*	iter;
	const skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry*	finish;	

	start  = m_wavelenentries;
	iter   = start;
	finish = start + m_numwavelens;

	roundednm = skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_WavelengthEntry::RoundTheWavelength( wavelennm );
	entries_equal	comparator(roundednm);

	iter = std::find_if( start, finish, comparator );
	if (iter == finish)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_TabulatedCrossSectionAndPhaseMatrix_HeightEntry::FindWavelengthEntry, could not find a match for wavelength %g", (double)wavelennm);
		iter = NULL;
	}
	*entry = iter;
	return (iter != NULL);
}
