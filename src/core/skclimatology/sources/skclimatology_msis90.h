
/*-----------------------------------------------------------------------------
 *					class nxmsis90									2004-11-23*/
/** \internal
 *	A C++ class interface to the MSIS-90 Fortran model.  It provides a
 *	convenient neutral density model for the stratosphere, mesosphere and
 *	thermospshere.  This class will always invoke MSIS whewnever any of the requested
 *	coordinates are changed.  Most users will normally use class #nxmsis90cached which
 *	caches the height profile of the various species for a given latitude, and longitude.
 *
 *	\par Libraries
 *	To use this class you must link your code with the fortran code
 *	stored in library \b msis90lib.lib and you must have \b msis90lib.dll
 *	somewhere on your path.
 **/
/*---------------------------------------------------------------------------*/

class nxmsis90
{
	private:								// input parameters for MSIS
		nxTimeStamp		m_mjd;
		double			m_altitude;			// altitude in kms
		double			m_latitude;			// Geodetic latitude
		double			m_longitude;		// Geodetic longitude
		double			m_f107a;
		double			m_f107;
		double			m_ap[7];
		nxBOOL			m_isdirty;


	private:								// output parameters from MSIS

		double			m_HE;
		double			m_O;
		double			m_N2;
		double			m_O2;
		double			m_AR;
		double			m_TOTALMASS;
		double			m_H;
		double			m_N;
		double			m_exosphere_temp;
		double			m_altitude_temp;

	private:
		bool			IsDirty			() { return m_isdirty;}

	protected:
		nxBOOL			InvokeMsis90	();
		double			CurrentAltitude	() { return m_altitude;}
		void			SetDirty		() { m_isdirty = true; }
	void				SetF107			( double f107);						//!< Set the value of F10.7
	void				SetF107Avg		( double f107avg);					//!< Set the value of F10.7
	bool				SetAP			( const double* ap, int n);					//!< Set the value of F10.7

	public:
						nxmsis90		();
		virtual		   ~nxmsis90		(){}
		bool			DeepCopy		( const nxmsis90& other );

	double				Mjd				()		{ return m_mjd.MJD();}		//!< Return the MJD of the current model run
	double				Height			()		{ return m_altitude;}		//!< Return the altitude of the current model run
	double				Latitude		()		{ return m_latitude;}		//!< Return the latitude of the current model run
	double				Longitude		()		{ return m_longitude;}		//!< Return the longitude of the current model run
	double				F107Avg			()		{ return m_f107a;}			//!< return the average value of F10.7 used in the model run
	double				F107			()		{ return m_f107;}			//!< return the curent value of F10.7 used in the model run
	virtual	double		MeanMolecularWeight();								//!< Get the Mean Molecular weight at the current location
	virtual	double		MeanNumberDensity();								//!< Get the Mean number density ( particles /cm3 ?) at the current location
	virtual	double		T				();									//!< Get the temperature (in Kelvins?) at the current location
	virtual	double		P				();									//!< Get the pressure (Pascals) at the current location
	virtual	double		T_Exosphere		();									//!< Get the temperature of the exospshere from the current model run
	virtual	double		TotalMass		();									//!< Get the total mass density at the current location
	virtual	double		O2				();									//!< Get the molecular O2 number density at current location
	virtual double      O2_O2			();									//!< Get the square of O2 number density. Suitable for O2/O2 collisional induced absorption (CIA)
	virtual	double		O				();									//!< Get the atomic oxygen number density at current location
	virtual	double		N2				();									//!< Get the molecular N2 number density at the current location
	virtual	double		AR				();									//!< Get the AR number density at the current location
	virtual	double		N				();									//!< Get the atomic N number density at the current location
	virtual	double		H				();									//!< Get the atomic H number density at the current location
	virtual	double		HE				();									//!< Get the HE number density at the current location
	virtual void		UpdateCachedProfiles( double mjd, double latitude, double longitude );	//!< Run the model for this latitude, longitude and time
	void				SetLongitude	( double longi );					//!< Set the longitude
	void				SetLatitude		( double lati );					//!< Set the latitude
	void				SetHeight		( double h );						//!< Set the height
	void				SetMjd			( double mjd );						//!< Set the MJD
};



/*---------------------------------------------------------------------------
 *					class nxmsis90cached						   2002-10-23*/
 /** \internal
 */
/*-------------------------------------------------------------------------*/

class nxmsis90cached : public nxmsis90
{
	private:
				std::map< CLIMATOLOGY_HANDLE, nxSpline2>				m_profiles;
		typedef std::map< CLIMATOLOGY_HANDLE, nxSpline2>::iterator		iterator;
		typedef std::map< CLIMATOLOGY_HANDLE, nxSpline2>::value_type	value_type;

		bool				m_cacheisdirty;
		double				m_minH;						// Minimum altitude in KM
		double				m_maxH;						// maximum altitude in KM
		double				m_deltah;					// altitude resolution in KMS

	private:
		void			ConfigureSpline				( const CLIMATOLOGY_HANDLE& key, const nx1dArray<double>&	ht, nx1dArray<double>&	profile );
		void			SetDirtyCache				() { m_cacheisdirty = true; }

	public:
						nxmsis90cached			();
		virtual		   ~nxmsis90cached			(){}
		void			SetCacheAltitudeRange	( double minh, double maxh, double deltah = 1.0);
		bool			CacheIsDirty			() const { return m_cacheisdirty; }
		void			UpdateCachedProfiles	( double mjd, double latitude, double longitude );
		bool			InterpolateToHeight		( const CLIMATOLOGY_HANDLE& key, double h, double* value);
		bool			DeepCopy				( const nxmsis90cached& other );

	public:
		bool			SetHeightSpacingKMS		( double deltah);
		bool			SetMaxHeightKMS			( double maxh );
		bool			SetF10p7				( double f107);
		bool			SetF10p7avg				( double f107a);
		bool			SetAp					( const double* ap, int numpts);
		bool			AddSpecies		        ( const CLIMATOLOGY_HANDLE& key);
		bool			IsSupportedSpecies		( const CLIMATOLOGY_HANDLE& species );
};

/*-----------------------------------------------------------------------------
 *					skClimatology_MSIS90											2005-7-27*/
/** \ingroup skClimAir
 *	This is a class that implements the skClimatology interface for the msis90
 *	model. The MSIS90 database is a fortran subroutine program that is really
 *	intended for mesospheric climatology but has a convenient typical climatology
 *	that extends to the ground.
 *
 *	\par Supported Species
 *	MSIS has the capability to model many mesopsheric constituents.  To date this model only
 *	supports the following species:-
 *	\n
 *	-# CLIMATOLOGY_AIRNUMBERDENSITY_CM3
 *	-# CLIMATOLOGY_PRESSURE_PA
 *	-# CLIMATOLOGY_TEMPERATURE_K
 *	-# CLIMATOLOGY_O2_CM3
 *	\n
 *
 *	\par Range of Valid model output
 *	The model should be valid over the following ranges:
 *	\n
 *	-# All times and dates
 *	-# All latitudes and longitudes
 *	-# Ranges from 0 km to 120 km.
**/
/*---------------------------------------------------------------------------*/

class skClimatology_MSIS90 : public skClimatology
{
	private:
		nxmsis90cached			m_msis;

	public:
								skClimatology_MSIS90			();
		virtual				   ~skClimatology_MSIS90			();
		bool					DeepCopy						( const skClimatology_MSIS90& other); 
		nxmsis90cached&			CachedMSIS						() { return m_msis;}

	public:
		virtual	bool			UpdateCache			( const GEODETIC_INSTANT& placeandtime) override;
		virtual	bool			GetParameter		( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache) override;
		virtual bool			IsSupportedSpecies  ( const CLIMATOLOGY_HANDLE& species ) override;
//		virtual bool			CreateClone			 (skClimatology** clone)	const override;
};


