#include <skopticalproperties21.h>


/*-----------------------------------------------------------------------------
 *					skSpectralEmissionCollection_HitranIsotope_OLD::skSpectralEmissionCollection_HitranIsotope_OLD		2013-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralEmissionCollection_HitranIsotope::skSpectralEmissionCollection_HitranIsotope( const skSpectralLineCollection_HitranIsotope& isotope, skSpectralLineShape* lineshapeobject  )
										  : skSpectralLineCollection_HitranIsotope( isotope.ParentChemical(), isotope.MolNum(), isotope.IsoNum(), isotope.MoleculeInfo() ) 
{
	SetLineShapeObject( lineshapeobject );
	CreateSpectralLinesAndUniqueUpperStates( isotope );
	
}


/*-----------------------------------------------------------------------------
 *					skSpectralEmissionCollection_HitranIsotope_OLD::~skSpectralEmissionCollection_HitranIsotope_OLD		2013-3-14*/
/** **/
/*---------------------------------------------------------------------------*/

skSpectralEmissionCollection_HitranIsotope::~skSpectralEmissionCollection_HitranIsotope()
{
}


/*---------------------------------------------------------------------------
 *   skSpectralEmissionCollection_HitranIsotope::UpdateLocation   2020-08-20 */
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralEmissionCollection_HitranIsotope::UpdateLocation( double temperature, double pressure, const GEODETIC_INSTANT& geopt, skClimatology* atmosphere )
{
	bool ok;

	m_upperstates.UpdatePartitions(temperature);
	ok = skSpectralLineCollection::UpdateLocation( temperature, pressure, geopt, atmosphere);
	return ok;
}

/*---------------------------------------------------------------------------
 * skSpectralEmissionCollection_HitranIsotope::CreateSpectralLinesAndUniqueUpperStates  2020-07-27 */
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralEmissionCollection_HitranIsotope::CreateSpectralLinesAndUniqueUpperStates(const skSpectralLineCollection_HitranIsotope& isotope)
{
	std::list< const skSpectralLine_HitranLine*> lines;
	std::list<skSpectralLine_HitranEmission*>	 upperstate_emission_lines;	// Holds the original instance of each Hitran_UpperState_SpectralLine. Pointers to these objects are used elsewhere in this class
	skSpectralLineEntry_HitranEmission*			 entry;

	for ( auto line = isotope.begin(); line != isotope.end(); ++line )
	{
		const skSpectralLineEntry*			entry        = *line;
		const skSpectralLine*				spectralline = entry->SpectralLine();	
		const skSpectralLine_HitranLine*	hitranline   = dynamic_cast<const skSpectralLine_HitranLine*>( spectralline); 
		lines.push_back( hitranline );
	}
	
	m_upper_state_selector.make_list_of_unique_upperstates( lines, &m_upperstates, &upperstate_emission_lines);		// Search through the HITRAN lines and get a list of unique upper states and emission lines.

	ClearLines(upperstate_emission_lines.size());
	for ( auto iter = upperstate_emission_lines.begin(); iter != upperstate_emission_lines.end(); ++iter)	// Copy the emission lines over to this object as a collection of skSpectralLine objects.
	{																										// so
		entry = new skSpectralLineEntry_HitranEmission( *iter, LineShape() ) ;								// The objects in  upperstate_emission_lines are now passed over to the "entry" skSpectralLine object to control lifetime
		AddEntry( entry  );																					// And the entry object "lifettime" is now managed by "this" object
	}
	return true;

}


size_t skSpectralEmissionCollection_HitranChemical::g_numinstances = 0;

/*---------------------------------------------------------------------------
 * skSpectralEmissionCollection_HitranChemical::skSpectralEmissionCollection_HitranChemical2020-08-18 */
/** **/
/*---------------------------------------------------------------------------*/

skSpectralEmissionCollection_HitranChemical::skSpectralEmissionCollection_HitranChemical( const char* chemicalname, double lowerwavenumber, double upperwavenumber, bool hapicompliant, int isotopefilterid, const char* lowerstateglobalquantafilter, const char* upperstateglobalquantafilter)
	                                       : m_listoflines( chemicalname, lowerwavenumber,upperwavenumber, 0.0, hapicompliant, isotopefilterid, lowerstateglobalquantafilter, upperstateglobalquantafilter )
{
	m_lineshape                = nullptr;
	g_numinstances++;
}


/*---------------------------------------------------------------------------
 * skSpectralEmissionCollection_HitranChemical::~skSpectralEmissionCollection_HitranChemical2020-08-18 */
/** **/
/*---------------------------------------------------------------------------*/

skSpectralEmissionCollection_HitranChemical::~skSpectralEmissionCollection_HitranChemical()
{
	ClearIsotopeEmissions();
	if ( m_lineshape != nullptr)				m_lineshape->Release();
	g_numinstances--;
}


/*---------------------------------------------------------------------------
 *   skSpectralEmissionCollection_HitranChemical::ClearIsotopes   2020-08-18 */
/** **/
/*---------------------------------------------------------------------------*/

void skSpectralEmissionCollection_HitranChemical::ClearIsotopeEmissions()
{
	for ( iterator iter = m_isotope_emissions.begin(); iter != m_isotope_emissions.end(); ++iter)
	{
		skSpectralEmissionCollection_HitranIsotope* entry = *iter;
		delete entry;
	}
	m_isotope_emissions.clear();
}


/*---------------------------------------------------------------------------
 * skSpectralEmissionCollection_HitranChemical::SetLineShapeObject2020-08-18 */
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralEmissionCollection_HitranChemical::SetLineShapeObject( skSpectralLineShape* lineshapeobject )
{
	bool ok;

	ok = (lineshapeobject != nullptr) && ( lineshapeobject == m_lineshape);				// Is the new line shape object identical to the existing? Ignore if it is
	if (!ok)																			// Otherwise
	{																					// We destroy the current lisy of upper state emission isotopes 
		if (lineshapeobject != nullptr) lineshapeobject->AddRef();
		if (m_lineshape     != nullptr) m_lineshape->Release();
		m_lineshape = lineshapeobject;
		ClearIsotopeEmissions();
		CreateListOfIsotopeEmissions();
	}
	return true;
}


/*---------------------------------------------------------------------------
 * skSpectralEmissionCollection_HitranChemical::MakeListOfIsotopes2020-08-18 */
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralEmissionCollection_HitranChemical::CreateListOfIsotopeEmissions()
{

	bool	ok;

	ok = (m_isotope_emissions.size() > 0);													// Do we have our list of isotope emission lines already created
	if (!ok)																				// nope. 
	{
		ok = m_lineshape != nullptr;
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"skSpectralEmissionCollection_HitranChemical::CheckListOfIsotopeEmissions(), Cannot create list of isotope upper state emissions as no lineshape obejct is defined");
		}
		else
		{

			const std::map< size_t, skSpectralLineCollection_HitranIsotope>& isotopes = m_listoflines.Isotopes();
			size_t n = isotopes.size();

			for ( auto it = isotopes.begin(); it != isotopes.end(); it++ )
			{
				const skSpectralLineCollection_HitranIsotope&   isotope   = it->second;
				skSpectralEmissionCollection_HitranIsotope*		isotope_emission = new skSpectralEmissionCollection_HitranIsotope( isotope, m_lineshape);
				m_isotope_emissions.push_back( isotope_emission );
			}
		}
	}
	return ok;
}


/*---------------------------------------------------------------------------
 *   skSpectralEmissionCollection_HitranChemical::EmissionArray   2020-08-18 */
/** **/
/*---------------------------------------------------------------------------*/

bool  skSpectralEmissionCollection_HitranChemical::EmissionArray( const std::vector<double>& wavenum, std::vector<double>* emission)
{
	bool	ok = true;
	bool	ok1;
	double	abundance;
	std::vector<double>		signal;

	ok = (m_lineshape != nullptr);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralEmissionCollection_HitranChemical::Emission, you must set the the Lineshape object before calculating the emission");
	}
	else
	{
		signal.resize(wavenum.size());
		for ( iterator iter = m_isotope_emissions.begin(); iter != m_isotope_emissions.end(); ++iter)
		{
			std::fill(signal.begin(), signal.end(), 0.0);
			abundance = (*iter)->MoleculeInfo()->m_abundance;														// Get the isotopic abundance of this isotope
			ok1 = (*iter)->AddAbsorptionCrossSectionOrEmissionArray( wavenum, &signal);								// 
			for (size_t i = 0; i < signal.size(); i++)
			{
				(*emission)[i] += signal[i]*abundance;
			}
			ok = ok && ok1;
		}
	}
	if (!ok) nxLog::Record(NXLOG_WARNING,"skSpectralEmissionCollection_HitranChemical::EmissionArray, There were errors calculating the upper state emission. Thats not good");
	return ok;
}

/*---------------------------------------------------------------------------
 *     skSpectralEmissionCollection_HitranChemical::Emission      2020-08-18 */
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralEmissionCollection_HitranChemical::Emission( double nu, double* emission  ) const 
{
	bool	ok = true;
	bool	ok1;
	double	abundance;
	double	signal;
	
	ok = (m_lineshape != nullptr);
	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"skSpectralEmissionCollection_HitranChemical::Emission, you must set the the Lineshape object before calculating the emission");
	}
	else
	{
		for ( const_iterator iter = m_isotope_emissions.begin(); iter != m_isotope_emissions.end(); ++iter)
		{
			abundance = (*iter)->MoleculeInfo()->m_abundance;														// Get the isotopic abundance of this isotope
			ok1 = (*iter)->AbsorptionCrossSectionOrEmission( nu, &signal);								// 
			(*emission)+= signal*abundance;
			ok = ok && ok1;
		}
	}
	if (!ok) nxLog::Record(NXLOG_WARNING,"skSpectralEmissionCollection_HitranChemical::Emission, There were errors calculating the upper state emission. Thats not good");
	return ok;


}


/*---------------------------------------------------------------------------
 *  skSpectralEmissionCollection_HitranChemical::UpdateLocation   2020-08-20 */
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralEmissionCollection_HitranChemical::UpdateLocation( const GEODETIC_INSTANT& geopt, skClimatology* atmosphere )
{
	iterator	iter;
	bool		ok = true;
	bool		ok1;
	double		temperature;
	double		pressure = 0.0;


	ok =       atmosphere->GetParameter( SKCLIMATOLOGY_TEMPERATURE_K, geopt, &temperature, false );
	ok = ok && atmosphere->GetParameter( SKCLIMATOLOGY_PRESSURE_PA,   geopt, &pressure,    false );

	for (iter =m_isotope_emissions.begin(); !(iter == m_isotope_emissions.end()); ++iter)
	{
		ok1 = (*iter)->UpdateLocation( temperature, pressure, geopt, atmosphere );
		ok = ok && ok1;
	}
	return ok;
}

/*---------------------------------------------------------------------------
 *    skSpectralEmissionCollection_HitranChemical::UpdateCache    2020-08-20 */
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralEmissionCollection_HitranChemical::UpdateCache( const GEODETIC_INSTANT& geopt)
{
	bool ok;
	
	ok = m_listoflines.SelfBroadeningClimatology_UpdateCache(geopt);
	return ok;

}

/*---------------------------------------------------------------------------
 * skSpectralEmissionCollection_HitranChemical::SetInternalNumberDensityClimatology 2020-08-18 */
/** **/
/*---------------------------------------------------------------------------*/

bool skSpectralEmissionCollection_HitranChemical::SetSelfBroadeningClimatology( const CLIMATOLOGY_HANDLE& parameterguid, skClimatology* numberdensityclimatology )
{
	return m_listoflines.SetSelfBroadeningClimatology( parameterguid, numberdensityclimatology);
}

