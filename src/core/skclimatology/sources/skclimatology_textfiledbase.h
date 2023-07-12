/*-----------------------------------------------------------------------------
 *					class skClimatology_TextFileDatabase		2006-6-16*/
/**  \ingroup skClimmisc
 *	\deprecated 
 *	A class built to support Sasktran V1. The class is only for deprecated
 *	usage. Users developing new code should use skClimatology_UserDefinedTable.
 *	This class provides support for loading in the five text based
 *	height profile climatologies used in Sasktran V1.  The heights in each text file
 *	can be specified in meters or in kilometers.
 **/
/*---------------------------------------------------------------------------*/

class skClimatology_TextFileDatabase : public skClimatology
{
	private:
		nxString				m_path;
		nx2dArray<double>		m_ozonedensity;
		nx2dArray<double>		m_neutraldensity;
		nx2dArray<double>		m_aerosoldensity;
		nx2dArray<double>		m_temperatureprofile;
		nx2dArray<double>		m_NO2density;

		bool					LoadProfileFromTextFile		( const char* file, nx2dArray<double>* profile);
		bool					LoadProfileFromArray		( int numalts, double* alts, double* speciesarray, nx2dArray<double>* profile ); 
		bool					InterpolateValueAtAltitude	( double* value, double altitude, nx2dArray<double>& profile, const char* filename_of_cache, const char* speciesname );

	public:
								skClimatology_TextFileDatabase	( );
		virtual				   ~skClimatology_TextFileDatabase	( );

	public:
		bool					SetPath				( const char* path );
//		bool					DeepCopy			( const skClimatology_TextFileDatabase& other );
		virtual	bool			UpdateCache			( const GEODETIC_INSTANT& placeandtime ) override;
		virtual	bool			GetParameter		( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache) override;
		virtual bool			IsSupportedSpecies  ( const CLIMATOLOGY_HANDLE& species ) override;
//		virtual bool			CreateClone			( skClimatology** userclone ) const override;
		bool					SetSpeciesProfile	( const CLIMATOLOGY_HANDLE& species, int numalts, double* alts_km, double* profile );
};

