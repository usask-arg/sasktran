#pragma once

/*---------------------------------------------------------------------------
 *                  Class HitranLine_Upperstate                   2020-07-24 */
/** **/
/*---------------------------------------------------------------------------*/

class HitranLine_UpperState : public nxUnknown
{
	private:
		double		m_E;
		double		m_degeneracy;
		double		m_fractional_occupation;

	public:
					HitranLine_UpperState	( double E, double degeneracy)	{ m_E = E; m_degeneracy = degeneracy; m_fractional_occupation = 0.0;}
		void		SetEnergyLevel			( double E)						{ m_E = E;}
		void		SetFractionalOccupation	( double frac)					{ m_fractional_occupation = frac;}
		double		FractionalOccupation	() const						{ return m_fractional_occupation;}
		double		EnergyLevel				() const						{ return m_E;}
		double		Degeneracy				() const						{ return m_degeneracy;} 
		double		BoltzmannExponent		( double T) const;
};


/*---------------------------------------------------------------------------
 *                  Class HitranLine_UpperStates                  2020-07-24 */
/** **/
/*---------------------------------------------------------------------------*/

class HitranLine_UpperStates
{
	private:
		std::list< HitranLine_UpperState*>	m_upperstates;
		double								m_current_temperature;


	public:
					HitranLine_UpperStates	(){m_current_temperature= -9999.0;}
				   ~HitranLine_UpperStates	();
		bool		CreateNewUpperState		( HitranLine_UpperState** newupperstate, double mu, double degeneracy );
		bool		UpdatePartitions		( double T) ;
};


/*---------------------------------------------------------------------------
 *              Class skSpectralLine_HitranEmission               2020-07-28 */
/** **/
/*---------------------------------------------------------------------------*/

class skSpectralLine_HitranEmission : public skSpectralLine_HitranLine
{
	private:
		const HitranLine_UpperState*			m_upperstate;

	public:
		static size_t							g_numinstances;
	private:
												skSpectralLine_HitranEmission	( const skSpectralLine_HitranEmission& other);		// Dont expose copy constructor
	public:
												skSpectralLine_HitranEmission	( const skSpectralLine_HitranLine* hitranline, const HitranLine_UpperState* upperstate	);
											   ~skSpectralLine_HitranEmission	();
		const HitranLine_UpperState*			UpperState						() const { return m_upperstate; }
		virtual	bool							CalculateLineIntensity			(  double T) override;
		virtual double							Snm								( )	const	override {return 0.0;}
};

/*---------------------------------------------------------------------------
 *          Class HitranLine_SelectUpperStatesFromLines           2020-07-24 */
/** **/
/*---------------------------------------------------------------------------*/

class HitranLine_SelectUpperStatesFromLines
{
	private:
		double		m_delta_E;							// Upper states levels closer than this with the same degeneracy are considered part of the same upper state level

	private:
		double		median_value							( std::vector<double>& sorted_mu_array  );

	public:
					HitranLine_SelectUpperStatesFromLines	() { m_delta_E = 0.0005; }
		bool		make_list_of_unique_upperstates			(	std::list< const skSpectralLine_HitranLine*>&		all_spectral_lines,
																HitranLine_UpperStates*								upperstates,
																std::list<skSpectralLine_HitranEmission*>*			spectrallines);
};



/*-----------------------------------------------------------------------------
 *					class skSpectralLineEntry_HitranEmission						2013-3-7*/
/** A small lightweight class used to hold all the information used to
 *	calculate spectra from line spectra. The class holds a pointer to a
 *	spectral line entry (from Hitran for instance), a pointer to a line shape
 *	object, Voigt for example and an optional work/storage buffer. The line entry
 *	is one of several (or indeed many) lines held by a molecule. The entry, becuase
 *	it is lightweight can be easily passed around.
 *
 *	The entry is especially lightweight if the lineshape object is set to NULL as there
 *	is no requirement to make a storage buffer. Thus higher level classes that load
 *	all the lines for a molecule, which can be many thousands, should not define line shape objects
 *	(unless they dont mind the memory allocation hit). Line shape objects should be set
 *	after lines have been selected into smaller groups of micro-windows.
 */
/*---------------------------------------------------------------------------*/

class skSpectralLineEntry_HitranEmission : public skSpectralLineEntry
{
	public:
												skSpectralLineEntry_HitranEmission	( );
												skSpectralLineEntry_HitranEmission	( skSpectralLine_HitranEmission* line, skSpectralLineShape* m_lineshapeobject 	);
												skSpectralLineEntry_HitranEmission	( const skSpectralLineEntry_HitranEmission& other );
											   ~skSpectralLineEntry_HitranEmission	( );
};

/*---------------------------------------------------------------------------
 *               Class skSpectralEmissionCollection_HitranIsotope 2020-07-27 */
/** A class used to generate a collection of emission lines for a given isotope.
 *	This class inherits the collection of hitran lines for the given isotope
 *	and then provides extra code to determine the list of unique upper states
 *	from the list of Hitran lines.
 **/
/*---------------------------------------------------------------------------*/

class skSpectralEmissionCollection_HitranIsotope : public skSpectralLineCollection_HitranIsotope
{
	private:
		HitranLine_UpperStates					m_upperstates;					// The collection of unique upper states.
		HitranLine_SelectUpperStatesFromLines	m_upper_state_selector;			// An object that decides the unique list of upper states. This can be replaced if the criteria change 

	private:
		bool									CreateSpectralLinesAndUniqueUpperStates		( const skSpectralLineCollection_HitranIsotope&    isotope);
												skSpectralEmissionCollection_HitranIsotope  ( const skSpectralEmissionCollection_HitranIsotope& other);

	public:
												skSpectralEmissionCollection_HitranIsotope	( const skSpectralLineCollection_HitranIsotope& isotope, skSpectralLineShape* lineshapeobject );
		virtual								   ~skSpectralEmissionCollection_HitranIsotope  ();
		virtual bool							UpdateLocation							( double temperature, double pressure, const GEODETIC_INSTANT& geopt, skClimatology* atmosphere ) override;

};



/*---------------------------------------------------------------------------
 *       Class skSpectralEmissionCollection_HitranChemical        2020-08-18 */
/** **/
/*---------------------------------------------------------------------------*/

class skSpectralEmissionCollection_HitranChemical
{
	private:
				std::list< skSpectralEmissionCollection_HitranIsotope*>						m_isotope_emissions;			// The array of isotopes associated with this chemical
		typedef std::list< skSpectralEmissionCollection_HitranIsotope*>::value_type			value_type;
		typedef std::list< skSpectralEmissionCollection_HitranIsotope*>::iterator			iterator;
		typedef std::list< skSpectralEmissionCollection_HitranIsotope*>::const_iterator		const_iterator;

	private:
		skSpectralLineCollection_HitranChemical		m_listoflines;							// The complete list of lines from the Hitran database for this chemical
		skSpectralLineShape*						m_lineshape;


	public:
		static size_t								g_numinstances;

	private:
		bool										CreateListOfIsotopeEmissions				();
		void										ClearIsotopeEmissions						();

	public:
													skSpectralEmissionCollection_HitranChemical	( const char* chemicalname, double lowerwavenumber, double upperwavenumber, bool hapicompliant,  int iostopefilterid, const char* lowerstateglobalquantafilter, const char* upperstateglobalquantafilter);
												   ~skSpectralEmissionCollection_HitranChemical	();
		bool										SetLineShapeObject							( skSpectralLineShape* lineshapeobject );
		bool										SetSelfBroadeningClimatology				( const CLIMATOLOGY_HANDLE& parameterguid, skClimatology* numberdensityclimatology );
		bool										UpdateLocation								( const GEODETIC_INSTANT& geopt, skClimatology* atmosphere );
		bool										UpdateCache									( const GEODETIC_INSTANT& geopt);
		bool										Emission									( double nu, double* absxsec  ) const;
		bool										EmissionArray								( const std::vector<double>& wavenum, std::vector<double>* absxs);
};

class skEmission_HitranChemical;

/*---------------------------------------------------------------------------
 *                 Class Hitran_Emission_Cache                  2019-11-06 */
/** This class is used to try and help speed up engines such as HR which
 *	unable to make use of the speed optimized skEmission_HitranChemical::IsotropicEmissionArray.
 *  This class internally calls IsotropicEmissionArray whenever any new location
 *	in the atmosphere is encountered.
 optimize **/
/*---------------------------------------------------------------------------*/

class Hitran_Emission_Cache
{
	private:
		skEmission_HitranChemical*								m_parent;
		std::vector<double>										m_wavenum;				// Cached Wavenumbers in ascending order.
		std::vector<double>*									m_current_emission;		// Current cached cross-sections 
		std::vector<double>										m_blank_entry;			// A blank entry, so we dont have to use null ptrs when there are errors.

	private:
		        std::map< hitran_geodetic_point, std::vector<double> >				m_cached_entries;			// Cached cross sections from calls to SetLocation
		typedef std::map< hitran_geodetic_point, std::vector<double> >::iterator	iterator;
		typedef std::map< hitran_geodetic_point, std::vector<double> >::value_type	value_type;

	private:
		bool													CreateNewEntry				( const GEODETIC_INSTANT& geo_pt, iterator* iter);

	public:
																Hitran_Emission_Cache		( skEmission_HitranChemical* parent);
		bool													SetLocation					( const GEODETIC_INSTANT& geo_pt );
		bool													CalculateEmissions			( double wavenumber, double *signal );
		bool													SetCachedWavenumbers		( const std::vector<double>& wavenumbers );									
};



/*---------------------------------------------------------------------------
 *                Class skEmission_HitranChemical                 2020-08-18 */
/** **/
/*---------------------------------------------------------------------------*/

class skEmission_HitranChemical : public skEmission
{
	private:
		skClimatology*									m_atmospheric_state;				// The climatology used to calculate atmospheric state parameters such as pressure and temperature
		skClimatology*									m_upperstate_numberdensity;			// The climatology used to calculate the number density of excited upper state molecules. 
		CLIMATOLOGY_HANDLE								m_upperstate_numberdensity_handle;	// The GUID of the species for the number density of the excited upper state molecules.
		double											m_upperstate_numberdensity_value;	// The current number of excited upper state molecules. 
		skClimatology*									m_selfbroadeningclimatology;		// The climatology to calculate the number density of this hitran molecule in molecules/cm3. Used for partial pressure in line self-broadening calcs.
		CLIMATOLOGY_HANDLE								m_selfbroadeningclimatologyhandle;	// The GUID of the species for the number density of this hitran molecule.
		skSpectralEmissionCollection_HitranChemical*	m_hitranchemical;					// The object that calculates the  HITRAN upper state emissions. Not constructed until required
		nxString										m_chemicalname;						// The name of this hitran chmenical. This  is used when constructing m_hitranchemical
		Hitran_Emission_Cache*							m_xs_optimizer;						// An object to optimize cross-sections calculation for engines such as HR which cannot call CalculateCrossSectionsArray
		skSpectralLineShape*							m_lineshapeobject;					// The line shape object, typicall voigt-kuntz.
		bool											m_isdirty;
		double											m_lowwavenum;
		double											m_hihwavenum;
		bool											m_hapicompliant;					// uses files generated from the HAPI compliant interfaces
		int												m_isotopefilterid;
		std::string										m_lowerstateglobalquantafilter;
		std::string										m_upperstateglobalquantafilter;

	private:
		void											SetDirty							();
		bool											CheckDirtyAndUpdate					(const GEODETIC_INSTANT& pt);
		void											init								();
		bool											CalculateEmissionInternal			( double wavenumber, double* signal ) const;
		bool											CheckWavenumberIsAscending			( const std::vector<double>&	wavenumber) const;
		double											FactorToConvertTo_PhotonsPerCM2PerSecPerSteradianPerNM( double wavenumber ) const;

	private:
														skEmission_HitranChemical			( const skEmission_HitranChemical& other );			// Dont allow copy constructor
		skEmission_HitranChemical&						operator =							( const skEmission_HitranChemical& other );			// Dont allow assignment oeprator

	public:
														skEmission_HitranChemical			();
														skEmission_HitranChemical			(const char* chemicalname, double lowerwavenumber, double upperwavenumber);
		virtual										   ~skEmission_HitranChemical			();
		bool											SetChemicalName						( const char*  chemicalname);
		bool											SetWavenumberRange					( double lowwavenumber, double highwavenumber );
		bool											SetUpperStateNumberDensity			( skClimatology* upperstate_numberdensity );
		bool											SetUpperStateNumberDensityHandle	( const CLIMATOLOGY_HANDLE& parameterguid );
		bool											SetSelfBroadeningClimatology		( skClimatology* numberdensityclimatology );
		bool											SetSelfBroadeningClimatologyHandle	( const CLIMATOLOGY_HANDLE& parameterguid );
		bool											SetLineShapeObject					( skSpectralLineShape* lineshapeobject );
		bool											SetAtmosphericState					( skClimatology* climate );
		bool											SetIsotopeIdFilter					( int isotopefilterid );
		bool											SetLowerStateGlobalQuantaFilter		( std::string lowerstateglobalquantfilter);
		bool											SetUpperStateGlobalQuantaFilter		( std::string upperstateglobalquantfilter);

		skSpectralEmissionCollection_HitranChemical*	HitranChemical						() { return m_hitranchemical;}
		skClimatology*									AtmosphericStateClimatology			() { return m_atmospheric_state; }
		bool											EnableCachedEmissions				( double* wavenumbers, size_t numwave );

	public:
		virtual bool									UpdateLocation				( const GEODETIC_INSTANT& pt, bool isground ) override;
		virtual bool									UpdateCache					( const GEODETIC_INSTANT& pt) override;														// Allows the emission object to update any internal caches to this time and location
		virtual bool									IsotropicEmission			( double wavenumber, double* isotropicradiance)  override;									// Calculate the isotropic emission as a radiance at the specified wave-number at the location specified in last call to SetAtmosphericState
		virtual bool									IsotropicEmissionArray		( const std::vector<double>& wavenumber, std::vector<double>* isotropicradiance) override;
};






