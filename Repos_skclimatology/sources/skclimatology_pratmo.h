/*-----------------------------------------------------------------------------
 *					class Pratmo				2008-2-11*/
/**	\ingroup no2skClimatology
 *	A climatology that returns the NO2 number density from the Pratmo
 *	database.
 **/
/*---------------------------------------------------------------------------*/

class skClimatology_Pratmo : public skClimatology
{
	private:
		double						m_cachemjd;
		double						m_cachelatitude;
		double						m_cachelongitude;
		nxString					m_path;
		nx2dArray<double>			m_NO2density;
		size_t						m_timeIndex;					// The index of the current time loaded in m_NO2Density
		size_t						m_angleIndex;					// The index o fthe latitude loaded into m_NO2Density
		double						m_sza;							// The solar zenith angle of the data loaded in m_NO2density.
		bool						m_isPM;
		nxRegistryConfiguration		m_config;						//!< The object used to access the registry.

	private:
		void					ClearCurrentIndex			( );
		bool					LoadProfileFromBinaryFile	( const GEODETIC_INSTANT& placeandtime );
		bool					InterpolateValueAtAltitude	( const GEODETIC_INSTANT& placeandtime, double* value );
		bool					InitializeMemory			();
		int 					TimeIndex					( const GEODETIC_INSTANT& geopt );
		int						LatitudeIndex				( const GEODETIC_INSTANT& geopt );
		bool					IsPM						( const GEODETIC_INSTANT& geopt );
		double					SZA							( const GEODETIC_INSTANT& geopt );
		bool					FetchBaseDirectory			( );
		bool					DeepCopy					( const skClimatology_Pratmo& other);
		bool					GoToStartOfTimeIndex		( size_t timeIndex, std::ifstream& nfile);
		bool					GoToStartOfLatitudeIndex	( size_t latitudeIndex, bool isPM, std::ifstream& nfile, size_t* numszaprofiles );
		bool					InterpolateSZAprofiles		( size_t numszaprofiles, double sza, std::ifstream& nfile );




	public:
								skClimatology_Pratmo		( );
		virtual				   ~skClimatology_Pratmo		( );
		bool					SetPath						( const char* path );

	public:
		virtual	bool			UpdateCache					( const GEODETIC_INSTANT& placeandtime ) override;
		virtual	bool			GetParameter				( const CLIMATOLOGY_HANDLE& species, const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache) override;
		virtual bool			IsSupportedSpecies			( const CLIMATOLOGY_HANDLE& species ) override;
//		virtual bool			CreateClone					( skClimatology** userclone ) const override;
};
