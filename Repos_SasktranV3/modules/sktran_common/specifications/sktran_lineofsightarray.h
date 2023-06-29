
/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightEntry_V2		2008-4-17*/
/** @ingroup linesofsight
 *	A class that defines the observer and look direction in geographic
 *	coordinates. These values are typically the location and look direction of
 *	an instrument.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_LineOfSightEntry_V2 
{
	public:	// Dont change the numbers below as some code implicitly depends upon the order 
		enum VIEWING_TYPE  {    SKTRAN_VIEWING_TYPE_UNDEFINED              = 0,
				                SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATSPACE = 1,
				                SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATLIMB  = 2,
				                SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATNADIR = 3,
				                SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATSPACE = 4,
				                SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATLIMB  = 5,
				                SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATNADIR = 6,
								SKTRAN_VIEWING_TYPE_NEARGROUND             = 7,
								SKTRAN_VIEWING_TYPE_ENDOFLIST              = 8
							};

	private:
		VIEWING_TYPE	m_viewingtype;
		nxVector		m_observer;				//!< The observer's location
		nxVector		m_look;					//!< The look direction
		double			m_mjd;						//!< Optional MJD set by the user (0 if not used)

	public:
		void				Configure			( const nxVector& observer, const nxVector& lookunit, double mjd );
		const nxVector&		Observer			() const { return m_observer;}					//!< The observer's location
		const nxVector&		Look				() const { return m_look;}						//!< The look direction
		double				Mjd					() const { return m_mjd;}						//!< Optional MJD set by the user (0 if not used)
		double				ZenithAngleOfLook	() const { return m_look.AngleTo( m_observer);}
		VIEWING_TYPE		ViewingType			() const { return m_viewingtype;}
		VIEWING_TYPE		DefaultViewingType	( nxGeodetic& geoid, double toa_altitude ) const;
		static size_t		NumViewingTypes     ()			{ return (size_t)SKTRAN_VIEWING_TYPE_ENDOFLIST - (size_t)SKTRAN_VIEWING_TYPE_UNDEFINED;}
		static VIEWING_TYPE IntegerToViewType	( size_t idx);
		static size_t		ViewTypeToInteger	( VIEWING_TYPE viewtype);
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_LineOfSightArray_V21						2008-4-17*/
/** @ingroup linesofsight
 *	A class that stores a set of observer locations and look directions. It is implicitly
 *	assumed that all of the observations and look directions are in geographic
 *	coordinates. These values are typically the location and look direction of
 *	an instrument. Sasktran must correct the oberver position for the osculating
 *	sphere offset.
 *
 *	It is implicitly assumed that the array of measurements can be adequately
 *	characterized by one sun location and one instant in time.  Given that the sun
 *	is 0.5 degrees wide and moves around the Earth at 1 degree every 4 minutes then
 *	it is reasonable to assume all measurements collected within 2 minutes can be reasonably
 *	approximated by a single sun and single instant in time. 
 *
 *	The sun can be explicitly set by the user or the class can determine the sun from the
 *	current value of mjd.
 *
 *	The mjd of the set of rays can be set manually by the user or it can be determined from the
 *	average of the rays within the valid height range and with valid mjd values (mjd >1000).  
 *
 *	The range of solar zenith angles is calculated from the sza at the end and middle points
 *	of each ray. This is only approximate, especially when you have rays close to
 *	solar zenith angles of zero or 180. In addition it does not account for the
 *	solar zenith angle of the diffuse rays emanating away from the primary rays.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_LineOfSightArray_V21
{
	private:
		        std::vector<SKTRAN_LineOfSightEntry_V2>				m_linesofsight;			//!< The lines of sight. Set by the user
		typedef std::vector<SKTRAN_LineOfSightEntry_V2>::iterator	iterator;				//   Iterator definition


	public:
											SKTRAN_LineOfSightArray_V21					();
		virtual							   ~SKTRAN_LineOfSightArray_V21					();

		bool								DeepCopy									( const SKTRAN_LineOfSightArray_V21& other );
		bool								AddLineOfSight								( const nxVector& observer, const nxVector& lookunit, double mjd = 0.0);

		bool								AddLineOfSightFromTangentPoint				( const GEODETIC_INSTANT& tangentpoint, double geographicbearingofobserver_degrees, double heightofobserver_meters, nxGeodetic* geoid);

		bool								SetRaysFromTangentHeightArray				( double  mjd, double lat, double lng,
																						  double  sza, double saa,
																						  double  rayazi,
																						  double* tangentheight_meters, int numheights,
																						  double  nominal_observerheight_meters,
																						  nxVector*	sunvector);
		bool								AddEquatorialLineOfSight					(double sza, double saa, double tanheight, double satAltitude, double mjd);

		const SKTRAN_LineOfSightEntry_V2*	Entry										(size_t index) const { return &(m_linesofsight[index]);}
		void								Clear										();
//		bool								SetMeanMJD									( double mjd );								//!< Sets MJD of rays, overides  average of ray entry mjd's
		size_t								NumRays										( ) const { return m_linesofsight.size();}
		bool								GetRay										( size_t idx, const SKTRAN_LineOfSightEntry_V2** entry ) const;
		bool								GetRayVar									( size_t idx, SKTRAN_LineOfSightEntry_V2** entry );
		double								MeanMJD										( ) const;
};
