#include <skopticalproperties21.h>
#include <omp.h>

size_t	skSpectralLine::g_numinstances = 0;

/*-----------------------------------------------------------------------------
 *					skSpectralLine::skSpectralLine		2013-3-19*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLine::skSpectralLine( ) 
{
	g_numinstances++;
	m_parentmolecule         = NULL;
	m_currentT               = -9999.0;
	m_current_line_intensity = 0.0;
}	

/*-----------------------------------------------------------------------------
 *					skSpectralLine::~skSpectralLine		2013-3-19*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLine::~skSpectralLine( )
{
	g_numinstances--;
};


/*-----------------------------------------------------------------------------
 *					skSpectralLine::SetParentMolecule		2013-3-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLine::SetParentMolecule( const skSpectralLineCollection* parent )
{
	m_currentT               = -9999.0;
	m_current_line_intensity = 0.0;
	m_parentmolecule         = parent; 
	NXASSERT(parent != NULL); 
	return (parent != NULL);
}

/*-----------------------------------------------------------------------------
 *					skSpectralLine::LineIntensity		2013-3-14*/
/** Follows equation A11 for linestrength in the Hitran 1996 paper.
 *	\f[ S_{\eta \eta'}(T)=S_{\eta \eta'}(T_{ref})\frac{Q(T_{ref})}{Q(T)}\frac{\exp(-c_2E_\eta/T)}{\exp(-c_2E_\eta/T_{ref})}\frac{[1-\exp(-c_2\nu_{\eta\eta'})/T]}{[[1-\exp(-c_2\nu_{\eta\eta'})/T_{ref}]}
 *  \f]
**/
/*---------------------------------------------------------------------------*/

bool skSpectralLine::CalculateLineIntensity( double T )
{
	static const double c2	= 1.438775220950165353239701610176;	// cm*K ** 2nd radiation constant hc/k
	double	A;
	double	B;
	double  C;
	double	D;
	double	Qref ;
	double  Qt;
	double	TRef;
	double  nu;
	
	NXASSERT((m_parentmolecule != NULL));
	if (T != m_currentT)
	{
		nu    = Nu();
		TRef  = Tref();
		Qref  = m_parentmolecule->QPartition( TRef );
		Qt    = m_parentmolecule->QPartition( T );
		A     = Qref/Qt;
		B     =     exp(  c2 * ELower() * ( 1.0/TRef - 1.0/T ) );
		C     = 1 - exp( -c2 * nu / T );
		D     = 1 - exp( -c2 * nu / TRef );
		SetCurrentLineIntensity(Snm() * A * B * C / D, T );
	}
	return true;
}

size_t skSpectralLineCollection::g_numinstances = 0;

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::skSpectralLineCollection		2013-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineCollection::skSpectralLineCollection(double microwindow_min, double microwindow_max)
{
	g_numinstances++;

	m_lineshape = NULL;
	m_maxlinestrength = 0.0;
	m_microwindow_minwavenum = microwindow_min;
	m_microwindow_maxwavenum = microwindow_max;

}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::~skSpectralLineCollection		2013-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineCollection::~skSpectralLineCollection()
{
	ClearLines(0);
	if (m_lineshape != NULL)  m_lineshape->Release();
	g_numinstances--;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::SetLineShapeObject		2013-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection::SetLineShapeObject	( skSpectralLineShape* lineshapeobject )
{
	iterator	iter;
	bool		ok = true;
	bool		ok1;

	if (lineshapeobject != NULL) lineshapeobject->AddRef();
	if( m_lineshape     != NULL) m_lineshape->Release();
	m_lineshape = lineshapeobject;

	for (iter =m_lines.begin(); !(iter == m_lines.end()); ++iter)
	{
		ok1 = (*iter)->SetLineShapeObject(m_lineshape);
		ok = ok && ok1;
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *              skSpectralLineCollection::ClearLines              2019-11-07 */
/** **/
/*---------------------------------------------------------------------------*/

void skSpectralLineCollection::ClearLines(size_t reservesize)
{
	iterator	iter;
	m_maxlinestrength = 0.0;
	for (iter= m_lines.begin(); iter != m_lines.end(); ++iter)
	{
		(*iter)->Release();
	}
	m_lines.clear();
	m_linesforthreadaccess.clear();
	m_threadstorage.clear();
	m_lines.reserve( reservesize );
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::AddEntry		2014-2-26*/
/** Add the new spectral line and keep track of the maximum line strength **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection::AddEntry( skSpectralLineEntry* entry )
{
	entry->AddRef();
	m_lines.push_back(entry); 
	m_linesforthreadaccess.resize(0);
	return true;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::SetAtmosphericState		2013-3-19*/
/** Set the atmospheric state given a location and a climatology.
 *	Note that temperature and pressure are passed in as an optimization.
 *	Temperature and pressure will normally be retrieved from the same
 *	geodetic point and climatology by the caller. The purpose is to
 *	avoid each and every spectral line requesting this info from
 *	the atmosphere object.

 *	Each spectral line will calculate and cache its line strength at this point
 *	and each line shape object will work out its cached parameters.
**/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection::UpdateLocation	( double temperature, double pressure, const GEODETIC_INSTANT& geopt, skClimatology* atmosphere )
{
	iterator				iter;
	skSpectralLineEntry*	entry;
	bool					ok = true;
	bool					ok1;

	m_maxlinestrength = 0.0;
	for (iter = begin(); !(iter == end()); ++iter)
	{
		entry = *(iter);
		ok1 = entry->ConfigureLineParameters(	temperature, pressure, geopt, atmosphere);
		if (ok1) 
		{
			double	nu = entry->SpectralLine()->Nu();
			if ( ( nu >= m_microwindow_minwavenum) && (nu < m_microwindow_maxwavenum) )						// See if this line is inside the current micro-window
			{																								// If it is then
				m_maxlinestrength = nxmax( m_maxlinestrength, entry->SpectralLine()->LineIntensity() );		// Keep track of the maximum line strength loaded in.
			}																								// Ignore lines outside the current micro-window
		}
		ok = ok && ok1;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::AbsorptionCrossSection		2013-3-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection::AbsorptionCrossSectionOrEmission(double nu, double* absxsec ) const
{
	const_iterator				iter;
	const skSpectralLineEntry*	entry;
	double					xsec;
	bool					ok = true;
	bool					ok1;

	*absxsec = 0.0;
	for (iter = begin(); !(iter == end()); ++iter)
	{
		entry     = *(iter);
		ok1       = entry->AbsorptionCrossSectionOrEmission(nu, &xsec);
		*absxsec += xsec;
		ok        = ok && ok1;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::CheckNumThreads		 2014- 10- 23*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection::CheckNumThreads()
{
	bool ok;

	ok = m_threadstorage.size() > 0;
	if (!ok)
	{
		ok = SetNumThreads( omp_get_num_procs() );
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::AddAbsorptionCrossSectionArray		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection::AddAbsorptionCrossSectionOrEmissionArrayMultiThreads( const std::vector<double>& nu, std::vector<double>* absxsec )
{
	iterator					iter;
	bool						ok = true;
	int							numlines;
	size_t						thrdindex;

	if (m_linesforthreadaccess.size() != m_lines.size())						// Make sure the vector version  of m_lines is up to date
	{																			// It should be the same size. Adding a line clears thr array
		m_linesforthreadaccess.resize(0);										// Clear the current array
		m_linesforthreadaccess.reserve(m_lines.size());							// reserve space for new storage
		for (iter = begin(); !(iter == end()); ++iter)							// Now copy over the lines
		{																		// Into our thread specific coopy
			m_linesforthreadaccess.push_back( (*iter) );						// We only need the pointers
		}
	}

	for ( thrdindex = 0; thrdindex < m_threadstorage.size(); thrdindex++)				// for each thread allocated by call to SetNumThreads
	{																					// so
		m_threadstorage.at(thrdindex).assign( nu.size(), 0.0 );							// resize the storage
	}

	numlines = (int) m_lines.size();													// Get the number of spectral lines to process

	skSpectralLineEntry*		entry;
	size_t						threadindex;
	#pragma omp parallel for schedule(dynamic) private (entry, threadindex)	num_threads ( (int)m_threadstorage.size())
	for (int iline = 0; iline < numlines; ++iline)
	{
		threadindex = omp_get_thread_num();												// between SetNumThreads and omp_set_num_threads
		entry       = m_linesforthreadaccess.at(iline);									// Get the line
		entry->AddAbsorptionCrossSectionOrEmissionArray(nu, &m_threadstorage.at(threadindex) );	// Calculate the cross-section
	}

	for ( size_t ithrdindex = 0; ithrdindex < m_threadstorage.size(); ithrdindex++)		// Now copy the cross-sections from the multiple threads
	{																					// back into the single copy.
																						
		for ( size_t iw = 0; iw < nu.size(); iw++)
		{
			absxsec->at(iw) += m_threadstorage.at(ithrdindex).at(iw);
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::AddAbsorptionCrossSectionArray		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection::AddAbsorptionCrossSectionOrEmissionArraySingleThreads( const std::vector<double>& nu, std::vector<double>* absxsec )
{
	iterator				iter;
	skSpectralLineEntry*	entry;
	bool					ok = true;
	bool					ok1;

	for (iter = begin(); !(iter == end()); ++iter)
	{
		entry     = *(iter);
		ok1       = entry->AddAbsorptionCrossSectionOrEmissionArray(nu, absxsec);
		ok        = ok && ok1;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::AddAbsorptionCrossSectionArray		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection::AddAbsorptionCrossSectionOrEmissionArray(  const std::vector<double>& nu, std::vector<double>* absxsec) 
{
	CheckNumThreads();

	return  (m_threadstorage.size() > 0) ? AddAbsorptionCrossSectionOrEmissionArrayMultiThreads ( nu, absxsec)
		                                 : AddAbsorptionCrossSectionOrEmissionArraySingleThreads( nu, absxsec);
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::SetNumThreads		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool  skSpectralLineCollection::SetNumThreads( size_t numthreads)
{
	m_threadstorage.resize(numthreads);
	return true;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::SetLineLimitsfromMaxLineStrength		2014-2-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection::SetLineLimitsfromMaxLineStrength( double parentmaxlinestrength)
{
	iterator					iter;
	skSpectralLineEntry*		entry;
	bool						ok = true;
	bool						ok1;

	for (iter = begin(); !(iter == end()); ++iter)
	{
		entry     = *(iter);
		ok1       = entry->SetLineLimitsfromMaxLineStrength(parentmaxlinestrength);
		ok        = ok && ok1;
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineCollection::SetTolerance			2014-2-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineCollection::SetLineTolerance( double tolerance )
{
	iterator					iter;
	skSpectralLineEntry*		entry;
	bool						ok = true;
	bool						ok1;

	for (iter = begin(); !(iter == end()); ++iter)
	{
		entry     = *(iter);
		ok1       = entry->SetTolerance( tolerance);
		ok        = ok && ok1;
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineEntry::skSpectralLineEntry		2013-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineEntry::skSpectralLineEntry()
{
	m_spectralline      = NULL;
	m_lineshapeobject   = NULL;
	m_storagebuffer     = NULL;

}

/*-----------------------------------------------------------------------------
 *					skSpectralLineEntry::skSpectralLineEntry		2013-3-8*/
/** The copy constructor copies the spectral line and line shape object
	but will create a new storage buffer object for this instance.
 */
/*---------------------------------------------------------------------------*/

skSpectralLineEntry::skSpectralLineEntry( const skSpectralLineEntry& other )
{
	m_spectralline      = NULL;
	m_lineshapeobject   = NULL;
	m_storagebuffer     = NULL;
	SetSpectralLine   ( other.m_spectralline );
	SetLineShapeObject( other.m_lineshapeobject );
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineEntry::~skSpectralLineEntry		2013-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralLineEntry::~skSpectralLineEntry()
{
	if (m_storagebuffer     != NULL) m_storagebuffer->Release();
	if (m_spectralline      != NULL) m_spectralline->Release();
	if (m_lineshapeobject   != NULL) m_lineshapeobject->Release();
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineEntry::SetSpectralLine		2013-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineEntry::SetSpectralLine( skSpectralLine* spectralline )
{
	if (spectralline   != NULL) spectralline->AddRef();
	if (m_spectralline != NULL) m_spectralline->Release();
	m_spectralline = spectralline;
	return true;
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineEntry::SetLineShapeObject		2013-3-7*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineEntry::SetLineShapeObject( skSpectralLineShape* lineshapeobject )
{
	bool	ok;

	if (lineshapeobject   != NULL) lineshapeobject->AddRef();
	if (m_lineshapeobject != NULL) m_lineshapeobject->Release();
	if (m_storagebuffer   != NULL) m_storagebuffer->Release();
	m_lineshapeobject = lineshapeobject;
	m_storagebuffer   = NULL;
	ok = (lineshapeobject == NULL);
	if (!ok)
	{
		ok = m_lineshapeobject->CreateStorageBuffer( &m_storagebuffer );
	}
	return ok;

}

/*-----------------------------------------------------------------------------
 *					skSpectralLineEntry::SetLineLimitsfromMaxLineStrength		2014-2-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineEntry::SetLineLimitsfromMaxLineStrength	( double parentmaxlinestrength) 
{ 
	return m_lineshapeobject->SetParentMaxLineStrength( parentmaxlinestrength, m_spectralline, m_storagebuffer);
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineEntry::SetTolerance		2014-2-26*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineEntry::SetTolerance( double tolerance)
{
	return m_lineshapeobject->SetTolerance( tolerance, m_spectralline, m_storagebuffer);
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineEntry::ConfigureLineParameters		2013-3-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineEntry::ConfigureLineParameters(	double temperature, double pressure, const GEODETIC_INSTANT& geopt, skClimatology* atmosphericstate)
{
	bool	ok;

	ok = (m_lineshapeobject != NULL) && (m_spectralline != NULL) && ( m_storagebuffer != NULL);
	ok = ok && m_spectralline->CalculateLineIntensity    (temperature);
	ok = ok && m_lineshapeobject->ConfigureLineParameters( m_spectralline, temperature, pressure, geopt, atmosphericstate, m_storagebuffer);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralLineEntry::ConfigureLineParameters, There were errors configuring the line parameters for this entry");
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineEntry::AbsorptionCrossSection		2013-3-19*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineEntry::AbsorptionCrossSectionOrEmission(	double nu,  double* absxsec ) const
{
	double	f;
//	double	Snm;
	bool	ok;

	ok       = m_lineshapeobject->LineShapeFunction(nu, &f, m_spectralline, m_storagebuffer);
//	Snm      = m_spectralline->LineIntensity(); (//LineShapeFunction includes the line strength
	*absxsec = ok ? f : 0.0;
	return ok;
}


/*-----------------------------------------------------------------------------
 *					skSpectralLineEntry::AbsorptionCrossSectionArray		2014-2-27*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineEntry::AddAbsorptionCrossSectionOrEmissionArray( const std::vector<double>& nu, std::vector<double>* absxsec )
{
	return  m_lineshapeobject->AddLineShapeFunctionArray(nu, absxsec, m_spectralline, m_storagebuffer);
}

/*-----------------------------------------------------------------------------
 *					skSpectralLineShape::LineShapeFunctionArray		2014-2-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralLineShape::AddLineShapeFunctionArray( const std::vector<double>&			nu, 
													 std::vector<double>*				uservalue, 
													 const skSpectralLine*				spectralline, 
													 skSpectralLineShapeStorageBuffer*	storagebuffer )
{
	size_t	numwave = nu.size();
	bool	ok = true;
	bool	ok1;
	size_t	i;

	NXASSERT(( nu.size() == uservalue->size() ));
//	uservalue->resize( numwave);
	for (i =0; i < numwave; i++)
	{
		ok1 = LineShapeFunction( nu.at(i), &uservalue->at(i), spectralline, storagebuffer);
		ok = ok && ok1;
	}
	return ok;
}
