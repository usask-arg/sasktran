#include <skopticalproperties21.h>
#include <nxbase_threads.h>

/*-----------------------------------------------------------------------------
 *					skconvolvedabsorbtionfuncptr::skconvolvedabsorbtionfuncptr		2012-4-30*/
/** **/
/*---------------------------------------------------------------------------*/

skconvolvedabsorbtionfuncptr::skconvolvedabsorbtionfuncptr()
{
	m_entry = NULL;
	m_w0    = 0.0;
	m_sd    = 0.0;
	m_xsid  = 1;
}


/*-----------------------------------------------------------------------------
 *					skconvolvedabsorbtionfuncptr::operator()		2012-4-30*/
/** This is a class that evaluates the convolution term at a given wavelength.
  *	It is called by the quadrarture that integrates over the point spread
  * function
 **/
/*---------------------------------------------------------------------------*/

double skconvolvedabsorbtionfuncptr::operator () ( double wavelennm)
{
	double					wavenum;
	double					x;
	bool					ok;
	skOpticalProperties*	optprop;
	double					v=0.0;
	double					absxs,extxs,scattxs;

	ok = (m_entry != NULL);
	if (ok)
	{
		optprop    = m_entry->HighResOpticalProperties();
		wavenum   = 1.0E7/wavelennm;
		ok		   = optprop->CalculateCrossSections(wavenum, &absxs, &extxs, &scattxs);
		if (ok)
		{
			switch ( m_xsid)
			{
				case 1: v = absxs; break;
				case 2: v = scattxs; break;
				case 3: v = extxs; break;
				case 4: v = 1; break;
				default: v = 0.0; break;
			};
			x   = (wavelennm - m_w0)/(m_sd);
			v  *= exp( -0.5*x*x );
		}
	}
	if (!ok) v = 0.0;
	return v;
}


/*-----------------------------------------------------------------------------
 *					skconvolvedabsorbtionfuncptr::Configure		2012-4-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool skconvolvedabsorbtionfuncptr::Configure( skOpticalProperties_ConvolvedDiscreteWavelenEntry* entry, double w0, double sd)
{
	m_entry = entry;
	m_w0 = w0;
	m_sd = sd;
//	m_threadindex = threadindex;
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperty_AdditionalStateInfoKey::skOpticalProperty_AdditionalStateInfoKey		2013-11-19*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperty_AdditionalStateInfoKey::skOpticalProperty_AdditionalStateInfoKey()
{
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperty_AdditionalStateInfoKey::~skOpticalProperty_AdditionalStateInfo
Key		2013-11-19*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperty_AdditionalStateInfoKey::~skOpticalProperty_AdditionalStateInfoKey()
{
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperty_AdditionalStateInfoKey::Clear		2013-11-19*/
/** **/
/*---------------------------------------------------------------------------*/

void skOpticalProperty_AdditionalStateInfoKey::Clear( )
{
	m_numstateparams     = -1;
	for( size_t i = 0; i < N_ELEMENTS(m_stateparameters); i++)
	{
		m_stateparameters[i] = std::numeric_limits<double>::quiet_NaN();
	}
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperty_AdditionalStateInfoKey::SetKeyStateParameters		2013-11-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperty_AdditionalStateInfoKey::SetKeyStateParameters( const double params[], size_t numparams)
{
	bool	ok;

	ok = (numparams <= N_ELEMENTS(m_stateparameters));
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperty_AdditionalStateInfoKey::SetKeyStateParameters, requested number of parameters(%d) bigger than number available (%d), truncating to max allowed", (int)numparams, (int)N_ELEMENTS(m_stateparameters) ) ;
		Clear();
	}
	else
	{
		for (size_t i = 0; i < numparams; i++)
		{
			m_stateparameters[i] = params[i];
		}
		m_numstateparams = (int)numparams;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperty_AdditionalStateInfoKey::operator==		2013-11-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperty_AdditionalStateInfoKey::operator == (const skOpticalProperty_AdditionalStateInfoKey& other ) const
{
	bool	isequal;

	isequal = (m_numstateparams == other.m_numstateparams) && (m_numstateparams>= 0);
	for (int i = 0; isequal && (i < m_numstateparams); i++)
	{
		isequal = isequal && (m_stateparameters[i] == other.m_stateparameters[i]);
	}
	return isequal;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperty_AdditionalStateInfoKey::operator<		2013-11-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperty_AdditionalStateInfoKey::operator <  (const skOpticalProperty_AdditionalStateInfoKey& other ) const
{
	bool	lessthan;
	bool	isequal;
	int		i;

	lessthan = (m_numstateparams <  other.m_numstateparams) && ( m_numstateparams >= 0);	// This index is less if it uses less parameters. Greater than if undefined
	isequal  = (m_numstateparams == other.m_numstateparams) && ( m_numstateparams >= 0);	// They may be equal if they same number of state params;
	i        = 0;
	while (isequal && (i < m_numstateparams) )
	{
		lessthan = (m_stateparameters[i] <  other.m_stateparameters[i]);	// This index is less if this parameter is les sthan the other
		isequal  = (m_stateparameters[i] == other.m_stateparameters[i]);	// only continue the loop if the parameters are equal, otherwise we are greater than or less than
		++i;																// loop over parameters until definitely less than
	}
	return lessthan;														
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenEntry::init		2012-4-30*/
/** The default empty constructor**/
/*---------------------------------------------------------------------------*/

void skOpticalProperties_ConvolvedDiscreteWavelenEntry::init()
{
	m_isdirty             = true;
	m_wavenumber          = 0;
	m_psf_fwhm_nm         = 0;
	m_highres_stepsize_nm = 0;
	m_convolvedabsorption = 0;
	m_convolvedextinction = 0;
	m_convolvedscattering = 0;
	m_highresopticalprops = NULL;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenEntry::skOpticalProperties_ConvolvedDiscreteWavelenEntry		2012-4-30*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_ConvolvedDiscreteWavelenEntry::skOpticalProperties_ConvolvedDiscreteWavelenEntry()
{
	init();
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenEntry::~skOpticalProperties_ConvolvedDiscreteWavelenEntry		2012-4-30*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_ConvolvedDiscreteWavelenEntry::~skOpticalProperties_ConvolvedDiscreteWavelenEntry()
{
	if (m_highresopticalprops != NULL) m_highresopticalprops->Release();
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenEntry::skOpticalProperties_ConvolvedDiscreteWavelenEntry		2012-4-30*/
/** Constructs the wavelngth entry with teh specified wavenumber. This
 *	is a convenience function to help the container indexing code.**/
/*---------------------------------------------------------------------------*/

skOpticalProperties_ConvolvedDiscreteWavelenEntry::skOpticalProperties_ConvolvedDiscreteWavelenEntry(double wavenum)
{
	init();
	m_wavenumber = wavenum;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenEntry::skOpticalProperties_ConvolvedDiscreteWavelenEntry		2012-4-30*/
/** **/
/*---------------------------------------------------------------------------*/

skOpticalProperties_ConvolvedDiscreteWavelenEntry::skOpticalProperties_ConvolvedDiscreteWavelenEntry( const skOpticalProperties_ConvolvedDiscreteWavelenEntry& other )
{
	m_isdirty             = other.m_isdirty;
	m_wavenumber          = other.m_wavenumber;
	m_psf_fwhm_nm         = other.m_psf_fwhm_nm;
	m_highres_stepsize_nm = other.m_highres_stepsize_nm;
	m_convolvedabsorption = other.m_convolvedabsorption;
	m_convolvedextinction = other.m_convolvedextinction;
	m_convolvedscattering = other.m_convolvedscattering;
	m_highresopticalprops = other.m_highresopticalprops;
	if (m_highresopticalprops != NULL) m_highresopticalprops->AddRef();
}



/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenEntry::SetHighResOpticalProperties		2012-4-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenEntry::SetHighResOpticalProperties(skOpticalProperties* optprop, double nm_resolution)
{
	if ( optprop                != NULL) optprop->AddRef();
	if ( m_highresopticalprops  != NULL) m_highresopticalprops->Release();
	m_highresopticalprops = optprop;
	m_highres_stepsize_nm = nm_resolution;
	SetDirty();
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenEntry::SetWavenumberAndFWHM		2012-4-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenEntry::SetWavenumberAndFWHM(double wavenumber, double fwhm_nm)
{
	m_wavenumber          = wavenumber;
	m_psf_fwhm_nm         = fwhm_nm;
	SetDirty();
	return true;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenEntry::GenerateConvolvedCrossSections		2012-4-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenEntry::GenerateConvolvedCrossSections( )
{
	skconvolvedabsorbtionfuncptr								funcptr;
	nxTrapezoidalQuadrature< skconvolvedabsorbtionfuncptr >		trapez;
	double														w0;
	double														wl;
	double														wh;
	int															npts;
	double														sigma;
	double														norm;
	double														normalcheck;	// then
	bool														ok;


	w0     = 1.E7/m_wavenumber;												// Get the central wavelength of the convolution
	sigma  = 0.42466090014400952136075141705144 * m_psf_fwhm_nm;			// Get the standard deviation from the full-wdith half max.	
	norm   = 0.39894228040143267793994605993438/sigma;						// 1/( sqrt(2.pi) * sigma )
	wl     = floor((w0 - 5.0*sigma)/m_highres_stepsize_nm + 0.5);			// get the starting point in units of step size
	wh     = floor((w0 + 5.0*sigma)/m_highres_stepsize_nm + 0.5);			// Get the end      point in units of step size 
	npts   = (int)(wh-wl + 1);												// Get the number of points
	trapez.SetRange	( wl*m_highres_stepsize_nm, wh*m_highres_stepsize_nm);
	trapez.SetOrder	( npts );

	funcptr.Configure  ( this, w0, sigma);											// Configure the function object that will evaluate cross-section and guassian product at eachw avelength

	funcptr.SetVariable( skconvolvedabsorbtionfuncptr::XSID_Normalize()  );			// configure function object to use unity
	normalcheck = trapez.Integrate( funcptr );										// do the integral
	NXASSERT( (normalcheck*norm > 0.995) && (normalcheck*norm < 1.05) );			// and make sure we are close to unity
	ok = (normalcheck*norm > 0.995) && (normalcheck*norm < 1.05);			// and make sure we are close to unity
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenEntry::GenerateConvolvedCrossSections, The convolution normalization is not very accurate (%e, it should be 1.00), you may be experiencing accurcay issues", (double)normalcheck);
	}

	if ( m_highresopticalprops->IsAbsorber() )										// If the high res object is an absorber
	{																				// then 
		funcptr.SetVariable( skconvolvedabsorbtionfuncptr::XSID_Absorption() );		// configure function object to use absorption cross-section
		m_convolvedabsorption = trapez.Integrate( funcptr )/normalcheck;			// and evaluate the convolution and normalize
	}
	else m_convolvedabsorption = 0;													// If the high res properties is not an absorber then the cross-section is zero
	

	if ( m_highresopticalprops->IsScatterer() )										// If the high res object is a scatterer
	{																				// then
		funcptr.SetVariable( skconvolvedabsorbtionfuncptr::XSID_Scatter()    );		// configure the function object to use the scattering cross-section
		m_convolvedscattering = trapez.Integrate( funcptr )/normalcheck;			// and do the integral
	}																				// otherwise
	else m_convolvedscattering = 0;													// not a scatterer so set cross-section to zero

	funcptr.SetVariable( skconvolvedabsorbtionfuncptr::XSID_Extinction() );			// always configure function object for extinction calculation
	m_convolvedextinction = trapez.Integrate( funcptr )/normalcheck;				// do the integral and save the cross-section

	return ok;
}

/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenEntry::CheckDirtyAndUpdate		2012-4-30*/
/** **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenEntry::CheckDirtyAndUpdate( )
{
	bool	ok;


	ok = !m_isdirty;
	if (!ok)
	{
		if (m_highresopticalprops == NULL)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenEntry::CheckDirtyAndUpdate, Cannot get convolved optical propertuies as no underlying optical properties are defined");
		}
		else if (m_wavenumber == 0.0)
		{
			nxLog::Record(NXLOG_WARNING,"skOpticalProperties_ConvolvedDiscreteWavelenEntry::CheckDirtyAndUpdate, Cannot get convolved optical propertuies as wavenumber is undefinedas zero");
		}
		else
		{
			{
				std::unique_lock<std::mutex>	lock(m_mutex);					// Lock this entry before we do the update. MOst multi-thread will not be accessing the same entry at the same time but you never know.
				ok = !m_isdirty;													// Make sure that the entry is still dirty. An earlier thread may have fixed the problem
				if (!ok)
				{
					ok = GenerateConvolvedCrossSections( );
					m_isdirty = !ok;
				}
			}
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING, "skOpticalProperties_ConvolvedDiscreteWavelenEntry::CheckDirtyAndUpdate, There was an error generating the convolved cross-sections");
			}
		}
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *			skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable::skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable		2012-5-3*/
/** Constructor for the table of convolved entries.
**/
/*---------------------------------------------------------------------------*/

skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable::skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable()
{
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable::FindEntry		2012-5-3*/
/** Find the convolved entry for the given wavenumber.
**/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable::FindEntry( double wavenumber, skOpticalProperties_ConvolvedDiscreteWavelenEntry** entry)
{
//	skOpticalProperties_ConvolvedDiscreteWavelenEntry		dummy(wavenumber);
	iterator						iter;
	bool							ok;

	iter   = m_entries.find( wavenumber );
	ok     = !(iter == m_entries.end());
	if (entry != NULL) *entry = ok ? &((*iter).second) : NULL;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable::SetDirty		2012-5-3*/
/** Flag the entire table as dirty. This will force all entries to to re-evaluate
 *	their internal cross-section data.
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable::SetDirty()
{
	iterator iter;

	for ( iter = m_entries.begin(); !(iter == m_entries.end()); ++iter)
	{
		(*iter).second.SetDirty();
	}
	return true;
}


/*-----------------------------------------------------------------------------
 *					skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable::AddEntry		2012-5-4*/
/** Add a new entry for the wavenumber to the table of convolved entries.
 *	If the entry already exists then modify the entry
 **/
/*---------------------------------------------------------------------------*/

bool skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable::AddEntry ( double wavenumber, skOpticalProperties* highresproperties, double nm_resolution, double fwhm_nm )
{
	skOpticalProperties_ConvolvedDiscreteWavelenEntry				dummy;
	skOpticalProperties_ConvolvedDiscreteWavelenEntry*				updateentry;
	bool															ok;
	std::pair <iterator, bool>										pr;
	value_type														entry(wavenumber,dummy);


	ok = !FindEntry(wavenumber,&updateentry);									// Find an existing entry. OK if it does not exist
	if (!ok)																	// if the entry does exist then modify
	{
		ok =       updateentry->SetHighResOpticalProperties( highresproperties, nm_resolution );
		ok = ok && updateentry->SetWavenumberAndFWHM       ( wavenumber, fwhm_nm );
	}
	else																		// if the entry does not exist then create
	{
		entry.second.SetHighResOpticalProperties( highresproperties, nm_resolution );
		entry.second.SetWavenumberAndFWHM       ( wavenumber, fwhm_nm );
		pr = m_entries.insert(  entry );
		ok = pr.second;
	}
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING, "skOpticalProperties_ConvolvedDiscreteWavelenEntriesTable::AddEntry Error adding entry for wavenumber %f, (wavelength %f nm)", (double)wavenumber, (double)(1.0E7/wavenumber) );
	}
	return ok;
}


