

/*-----------------------------------------------------------------------------
 *					class skClimatology_UserDefinedTable		2006-6-16*/
/**  \ingroup skClimmisc
 *	A class built to support climatologies derived from user defined tables and arrays.
 *
 *	\par Supported Species
 *	This model supports any species the user care to choose:-
 *
*	
 **/
/*---------------------------------------------------------------------------*/

class skClimatology_UserDefinedTable : public skClimatology
{
	private:
		nx2dArray<double>		m_profile;
		CLIMATOLOGY_HANDLE*		m_climatetypes;
		size_t					m_numclimate;

	private:
		void					ReleaseResources				();
		bool					AllocateClimatologyTable		( size_t numclimates );
		bool					InterpolateValueAtAltitude		( double* value, double altitude, size_t columnidx );

	public:
								skClimatology_UserDefinedTable	( );
								skClimatology_UserDefinedTable	( CLIMATOLOGY_HANDLE id, const char* filename );
		virtual				   ~skClimatology_UserDefinedTable	( );
		bool					DeepCopy						( const skClimatology_UserDefinedTable& other );

	public:
		virtual	bool			UpdateCache						( const GEODETIC_INSTANT& placeandtime ) override;
		virtual	bool			GetParameter					( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache) override;
		virtual bool			IsSupportedSpecies				( const CLIMATOLOGY_HANDLE& species ) override;
//		virtual bool			CreateClone						(skClimatology** clone)	const override;
		bool					LoadProfileFromTextFile			( const CLIMATOLOGY_HANDLE* species, size_t numspecies, const char* file );
		bool							LoadProfileFrom2dArray			( const CLIMATOLOGY_HANDLE* species, size_t numspecies, const nx2dArray<double>& profile ); 
		const nx2dArray< double >&		Return2dArrayProfile			( ) const { return m_profile; }
		nx2dArray< double >*			Return2dArrayProfileVar			( ) { return &m_profile; }
};

/*-----------------------------------------------------------------------------
 *					class UserTableSplineEntry		2014-1-30*/
/** **/
/*---------------------------------------------------------------------------*/

class UserTableSplineEntry
{
	private:
		nxSpline2				m_spline;						//!< The  
		double					m_badvalue;
		std::vector<double>		m_heights;						// Cached version of heights and profile
		std::vector<double>		m_profile;						// use for piecewise linear interpolation
		bool					m_dolog;						//!< If true then interpolate the logs of the signal
		bool					m_dopiecewiselinear;			//!< If true then do piecewise linear interpolation else do spline interpolation

	private:
		bool					CheckHeightsAreAscending	(  const std::vector<double>& h_meters ) const;
	public:
								UserTableSplineEntry	();
							   ~UserTableSplineEntry	();
		bool					CreateProfile			( const std::vector<double>& h_meters, const std::vector<double>& profile, bool dologinterpolation, bool dopiecewiselinear, double badvalue);
		double					Interpolate				( double h_meters );
};


/*-----------------------------------------------------------------------------
 *					class skClimatology_UserTableSpline		2006-6-16*/
/**  A class built to support spline interpolation, either linear or log, 
 *	of user defined height profiles. The class implements an unlimited number of species
 *	which are uniquely identified by their climatology handle. The class
 *	does not read in tables from files yet but that would not be that
 *	hard to implement.
 **/
/*---------------------------------------------------------------------------*/

class skClimatology_UserTableSpline : public skClimatology
{
	private:
		        std::map< CLIMATOLOGY_HANDLE, UserTableSplineEntry>					m_species;
		typedef std::map< CLIMATOLOGY_HANDLE, UserTableSplineEntry>::value_type		value_type;
		typedef std::map< CLIMATOLOGY_HANDLE, UserTableSplineEntry>::iterator		iterator;

	private:
		void					ReleaseResources				();
		bool					AllocateClimatologyTable		( size_t numclimates );
		bool					InterpolateValueAtAltitude		( double* value, double altitude, size_t columnidx );

	public:
								skClimatology_UserTableSpline	( );
		virtual				   ~skClimatology_UserTableSpline	( );
		void					ClearProfiles					( ) {m_species.clear();}

	public:
		virtual	bool			UpdateCache						( const GEODETIC_INSTANT& placeandtime );
		virtual	bool			GetParameter					( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache);
		virtual bool			IsSupportedSpecies				( const CLIMATOLOGY_HANDLE& species );
		bool					AddProfile						( const CLIMATOLOGY_HANDLE& species, const std::vector<double>& h_meters, const std::vector<double>& profile, bool dologinterpolation, bool dopiecewiselinear, double badvalue = std::numeric_limits<double>::quiet_NaN() );
};


