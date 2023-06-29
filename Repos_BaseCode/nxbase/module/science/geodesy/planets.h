/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/


/*-----------------------------------------------------------------------------
 *					PlanetaryBody		2004-11-23*/
/** \ingroup Geodesy
 *  This is the base class for calculating the positions of astronomical entities such
 *	stars, planets and satellites. Users will normally instantiate the specific derived class 
 *	for their purposes.  This base class provides the purely virtual function #UpdateECIPosition
 *	which, when overloaded by derived classes, calculates the ECI position of the object at the
 *	specified time.  Several derived classes will also work out the ECI velocity of the object at
 *	the same time. This code was originally developed in the 1990's and at that time we developed 
 *	polynomial approximation implementations for the Sun, Moon and artificial satellites. In addition
 *	there is also a Keplerian orbit class. This class works in an ECI (true equator, true equinox)
 *	coordinate system and uses meters for objects within the solar system and unit vectors for the
 *	distant stars.
 *
 *	\par Novas
 *	The US Astronomical Observatory Novas code has not yet been integrated into this library. This is still
 *	work to deb done.
 *
 *	\par JPL Ephemerides
 *	The JPL ephemerides have not been integrated into this library.  This is still work to be done
 *
 *	\par Rise/Set calculations
 *  The base class provides a few functions for working out when the object rises and sets as
 *	well as the planetary body's postion in topocentric (observer based) coordinates.
 *
 	\par Example

	\code
	PlanetSun							sun;
	nxTimeStamp							mjd;
	nxVector							geographic;
	nxVector							eci;

	mjd.SetUTC("2001-11-15 16:45:59");							// Get the position of the sun at this time
	sun.UpdateECILocation( mjd );								// Update ECI location of the Sun
	eci        = sun.Location();								// Get the ECI position of the sun
	geographic = eci.EquatorialToGeographic( mjd );				// Convert ECI to geographic uns

	\endcode

	\par Libraries
	nxbaseX.lib*
**/
/*---------------------------------------------------------------------------*/

class  PlanetaryBody
{
   protected:													// @access Protected data
      nxBOOL			m_DistantObject;						// !< Flags whether the object is distant or not.
      nxTimeStamp	        m_time;								// !< The time when object was at "location"
      nxVector 	        m_location;  							// !< Current location of object.

   protected:									// Keep these members public so derived classes can access them
      void 				ADDTHE( double C1, double S1, double C2, double S2,
		   						double &C, double &S);
      double 			FRAC( double X);
      double 			SINE ( double PHI);
      double 			SN( double X );
      double 			CS( double X );
      double 			TN( double X );
      double 			ASN( double X );
      double 			ACS(double X);
      double 			ATN2( double Y, double X );
      double 			ATN( double X );



	protected:	// @access Protected members
		void 	        	ConvertToEquatorialCoords ( nxBOOL UseTrueEclipticOfDate = nxFALSE);	//!< Convert ecliptic to ECI coords
		void 	        	ConvertToEclipticCoords   ( nxBOOL UseTrueEclipticOfDate = nxFALSE);	//!< Convert ECI to ecliptic coords
		nxBOOL				TimeOfContactChange       ( nxTimeStamp TStart, nxTimeStamp Tend, nxTimeStamp *ChangeTime, const nxGeodetic &Site, double minelevation, double stepsize ); //!< Time when contact changes

	public:	
		
		
static	void 	        	NutateEquatorialCoords( nxVector* location, nxTimeStamp tnow );		//!< Nutate the ECI coords
static	double				Ecliptic( const nxTimeStamp& tdt, nxBOOL UseTrueEclipticOfDate);	//!< Get the obliquity of the ecliptic at True Dynamical Time tdt
static	void				Nutation( const nxTimeStamp& tdt, double *DPSI, double *DEPS );		//!< Get the nutation corrections to RA and DEC  at True Dynamical Time tdt.
static	nxTimeStamp			TDT	    ( const nxTimeStamp& utc);									//!< Convert UT time Tnow to True Dynamical Time.
static	double				GAST    (       nxTimeStamp  utc);									//!< Get the Greenwhich Apparent Sidereal time at given UT time

	public:	 

         					PlanetaryBody();																								//!< Default Constructor
		virtual			   ~PlanetaryBody(){}																								//!< Destructor
		nxVector			ApparentECIPosition ( const nxGeodetic &Site);																						//!< Get the apparent position of object as seen from a specific location
		nxVector			Topocentric         ( const nxGeodetic &Site );																						//!< Convert the coordinates to TopoCentric
		nxBOOL				RiseSet             ( const nxTimeStamp &Tnow, nxTimeStamp* rise, nxTimeStamp* sets, const nxGeodetic& Site, double ZenithAngle);	//!< Get the time at which the object rises and sets
		nxBOOL				BelowHorizon        ( const nxTimeStamp &Tnow, const nxGeodetic& Site, double zenith);												//!< Determine if this body is below the observer's horizon
		nxBOOL				InContactWithGround ( nxTimeStamp &Tnow, const nxGeodetic&	Site, double minelevation);												//!< Return True if the object is in contact with the ground
		nxTimeStamp 		StartOfContact      ( const nxGeodetic& Site, double minelevation, nxTimeStamp Tnow, nxTimeStamp TEnd, nxTimeStamp *EndOfContact);	//!< calculate the start and end time of contact with ground

	public:
		const nxVector&		Location() const {return m_location;}	  							//!< Returns the location of the PlanetaryBody in ECI coordinates (metres for clsoe objects, unit vector for infinite objects)

	virtual void    		UpdateECIPosition   ( const nxTimeStamp &TNOW ) = 0;				//!< Purely virtual, Updates position (and possibly velocity) of the PlanetaryBody in ECI coordinates (metres for solar system objects, unit vector for infinite objects).

};

/*-----------------------------------------------------------------------------
 *					class  PlanetSun								2004-11-23*/
/** \ingroup Geodesy
 *	Class calculates the ECI position of the Sun with an accuracy better than
 *	1 arc minute.  In addition to the methods provided by #PlanetaryBody it also
 *	provides a mthod to get the Sun's Ecliptic coordinates and can return true 
 *	if a location is in astronomical dark. See class #nxPlanetaryBdy for an
 *	example on usage.
**/
/*---------------------------------------------------------------------------*/

class  PlanetSun: public PlanetaryBody
{
   	private:
     	double 		C3[9],S3[9];		// ARRAY [-1..7] OF double;
      	double 		C[9],S[9];			// ARRAY [-8..0] OF double;
      	double 		M2,M3,M4,M5,M6;
      	double 		D,A,UU;
      	double 		U,V,DL,DR,DB;
      	double 		T;					// Time in Julian centuries since JD2000

	private:
      	void 		TERM(  int  I1,  int  I,   int  IT,
		 			       double DLC, double DLS, double DRC,
		 			       double DRS, double DBC, double DBS);


      void 			PERTVEN();
      void 			PERTMAR();
      void 			PERTSAT();
      void 			PERTJUP();
      void 			PERTMOO();


   public:															// @access Public Methods
               		PlanetSun();									// !< Default Constructor
      void     		EclipticCoords   ( const nxTimeStamp  Tnow);		// !< Calculate the ecliptic coordinates
      nxBOOL 	    AstronomicalDark ( const nxTimeStamp &Tnow, const nxGeodetic &Site );	// !< Determine if it is astronomical dark or not at a given location
	  double		AUDistance		 ( double mjd );
   virtual  void	UpdateECIPosition( const nxTimeStamp &now );		// !< Calculate the Position of the SUN in ECI coords
};


/*-----------------------------------------------------------------------------
 *					class PlanetMoon								2004-11-23*/
/** \ingroup Geodesy
 *	Class calculates the ECI position of the Moon using a polynomial
 *  approximation with an accuracy better than
 *	1 arc minute.  In addition to the methods provided by #PlanetaryBody it also
 *	provides a mthod to get the Sun's Ecliptic coordinates and can return true 
 *	if a location is in astronomical dark. See class #nxPlanetaryBdy for an
 *	example on usage using #PlanetSun.
 **/
/*---------------------------------------------------------------------------*/

class  PlanetMoon: public PlanetaryBody
{
   	private:
      	double 		DGAM,FAC;
      	double 		DLAM,N,GAM1C,SINPI;
      	double 		L0, L, LS, F, D ,S;
      	double 		DL0,DL,DLS,DF,DD,DS;
      	double 		T;
      	double 		CO[13][4], SI[13][4];				// ARRAY[-6..6,1..4] OF double;

   	private:
      	void		LONG_PERIODIC (  double T, double &DL0, double &DL, double &DLS, double &DF, double &DD, double &DGAM);
      	void 		INIT();
      	void 		TERM( long P, long Q, long R, long S, double &X, double &Y);
      	void 		ADDSOL( double COEFFL, double COEFFS, double COEFFG, double COEFFP, long P, long Q, long R, long S);
      	void 		SOLAR1();
      	void 		SOLAR2();
      	void 		SOLAR3();
      	void 		ADDN( double COEFFN, long P, long Q, long R, long S, double &N);
      	void 		SOLARN( double &N);
      	void 		PLANETARY( double &DLAM);

	
   	public:																		//@access Public Members
                   		PlanetMoon();											//!< Default Constructor
      	void     		EclipticCoords( const nxTimeStamp  &Tnow);				//!< Calculate the ecliptic coords of the Moon
      	double	   		Phase( const nxTimeStamp &Tnow, const nxGeodetic &site);	//!< Get the phase of the moon as seen from a given location
      	virtual void 	UpdateECIPosition( const nxTimeStamp &now );				//!< Get the ECI coordinates of the moon.
};


/*-----------------------------------------------------------------------------
 *					class nxSatelliteBase									2004-11-23*/
/** \ingroup Geodesy
 *	nxSatelliteBase is a base class for artificial satellites orbiting the Earth.
 *  The class depends upon derived classes to provide the calculation of position and velocity
 *	at any instant in time. The class uses Earth Centred Intertial coordinates using the true equator and true
 *	equinox.  Individual satellite predictors should ensure they are consistent
 *	with true equator and true equator and equinox in all calculations. See class #nxPlanetaryBdy for an
 *	example on usage using #PlanetSun.
 **/
/*---------------------------------------------------------------------------*/

class nxSatelliteBase: public PlanetaryBody
{

	protected:																	// @access protected data
		nxTimeStamp		m_epoch;												// !< Start Epoch of this satellite (usually from 2-line elements).
		double			m_DaysPerRev;											// !< Period of the satellit in days
		long    		m_StartOrbitNumber;										// !< Orbit number at time "StartofOrbit"
		nxTimeStamp		m_StartOfOrbit;											// !< Start time of orbitnumber.
		nxVector		m_velocity;												// !< Velocity of the satellite.

	private:

		class EvalZcomponent
		{
			private:
				nxSatelliteBase*	object;

			public:
							EvalZcomponent( nxSatelliteBase* sat ) { object = sat;}
				double		operator() ( double t )	 { return object->ZComponent(t);}
		};

	public:																					//@access Public Methods
						nxSatelliteBase();															//!< Default Constructor
		double        	DepthOfEclipse		( nxTimeStamp Tnow );								//!< Angular distance of satellite from terminator 
		nxBOOL 			InEclipse			( nxTimeStamp & Tnow);								//!< Calculate if satellite is in Eclipse
		nxBOOL 			TimeOfTerminator	( nxTimeStamp TStart,nxTimeStamp Tend, nxTimeStamp &Terminator, double stepsize ); //!< Calculate instant when satellite crosses terminator
		double			Period				();													//!< Get the period of this satellite
		double			ZComponent			( const double mjd );								//!< Get satellite Z component at time mjd.
		nxTimeStamp		EquatorCrossing		( const nxTimeStamp & Tnow);						//!< Get time of last equator crossing before time Tnow
		long			OrbitNumber			( const nxTimeStamp &Tnow, nxTimeStamp* Start );	//!< Get the orbit number and orbit start time 
		long			OrbitNumber			( const nxTimeStamp &Tnow );						//!< Get the orbit number

	public:
 		nxTimeStamp		Epoch() const	 {return m_epoch;}									//!< Returns the Epoch for this satellite (this is normally the start time of the elements)
		nxVector		Velocity () const { return m_velocity;}							   	//!< Returns the velocity of the satellite in m/s


	public:
		virtual void	UpdateECIPosition( const nxTimeStamp &TNOW )=0;						//!< Update position and velocity of satellite.

};


/*<---------------------------------------------------------------------------
 *'			@class KeplerOrbit |
 *  Implements an elliptical Keplerian orbit.  By default it is setup for
 *	Earth orbits
 *
 *''Base Class
 *	public: nxSatelliteBase-->PlanetaryBody
 *>---------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
 *					class KeplerOrbit								2004-11-23*/
/** \ingroup Geodesy
 *  Implements an elliptical Keplerian orbit.  By default it is setup for
 *	Earth orbits
**/
/*---------------------------------------------------------------------------*/

class KeplerOrbit : public nxSatelliteBase
{
	private:
		double      m_mu;// = 3.98601210E14;	//!< The gravitational parameter for Earth (or other massive body)
		double		m_N0;						//!< Mean motion
		double		m_M0;						//!< Mean Anomaly	at the specified epoch
		double		m_e;						//!< eccentricity
		double		m_a;						//!< semi-major axis
		double		m_b;						//!< semi minor axis
		double		m_h;						//!< Magnitude of the angular velocity;
		nxVector	m_xunit;					//!< Get the unit vector pointing to the perigree (i.e. same direction as eccentricy vector)
		nxVector	m_yunit;					//!< Get the unit vector pointing along minor axis of ellips, (circa 90 degrees mean anomaly).
		nxVector	m_zunit;					//!< Get the unit vector perpendicular to the plane (ie in the direction of the angular momentum.
	
	private:								  								//@access private members
		double		EccentricAnomaly( double M, double ecc, double eps ); 	//!< Calculate eccentric anomaly

	public:										 						//@access public members
					KeplerOrbit			();															//!< Default Constructor
		double		OrbitalPeriod		( double apogee);											//!< Returns the orbital period for a given apogee
		void		FromStateVector		( double mjd, int orbitnumber, nxVector r, nxVector v);		//!<	Calculates internal elements from an ECI state vector
		void		FromElements		( double mjd, int orbitnumber, double I0, double RAAN, double W0, double E0, double N0,  double M0 ); //!<  Calculates the internal elements from Keplerian orbital elements.
virtual void    	UpdateECIPosition   ( const nxTimeStamp &mjd );  								//!<	Calculates the position of the satellite.

};



/*-----------------------------------------------------------------------------
 *					class  NoradSatellite		2004-11-23*/
/** \ingroup Geodesy
 *	Class NoradSatellite, implements satellite prediction that is common to
 *	all of the NORAD prediction models.  This class adds support for handling
 *	two line elements while derived classes implement the specific Norad models
 *	(eg SGP8). This class basically provides the methods to read two line elements
 *
**/
/*---------------------------------------------------------------------------*/

class  NoradSatellite: public nxSatelliteBase
{
	protected:
		char    		m_catnr[6];					//!< nxSatelliteBase catalogue number.
		long    		m_elset;					//!< Element set number.
		double  		m_XMO;						//!< The mean anomaly at epoch 			 (stored as radians).
		double			m_XNO;						//!< Mean motion,					 (stored as radians per minute).
		double  		m_XNODEO;          			//!< The mean longitude of ascending node at epoch (stored as radians)
		double  		m_OMEGAO;          			//!< The mean argument of perigree at epoch        (stored as radians)
		double  		m_EO;              			//!< The mean eccentricity` at epoch.
		double  		m_XINCL;           			//!< The mean inclination at epoch          	 (stored as radians).
		double  		m_XNDT20;          			//!< First time derivative of the Mean motion or Ballistic Coefficient. (stored as change/minute)
		double  		m_XNDD60;          			//!< Second time derivative (stored as change/minute)
		double  		m_BSTAR;          	  		//!< the SGP4 type drag coefficient.
		int				m_ideep;

	public:
						NoradSatellite();													//!< Default Constructor
		nxBOOL			GetElementsFrom( const char *Filename, const char *SatelliteName);	//!< Get the two line elements from Filename for SatelliteName
		void			TwoLine( const char* line1, const char * line2);   					//!< Decode two line elements.
		virtual void	Model_Init() = 0;													//!< Initialise the specific satellite model.
};


/*-----------------------------------------------------------------------------
 *					class  SGP8										2004-11-23*/
/** \ingroup Geodesy
 *	The NORAD SGP8 satellite predictor. The current implementation only
 *	handles near earth objects as I ported the orginal SGP8 fortran code
 *	to C++ and did not port the deep space model.  This is work that needs to
 *	be done.
 *
 * This class implements the NORAD SGP8 orbital prediction model.
 * It is suitable for all satellites whose mean period is less than
 * 225 minutes.  Satellites with periods longer than 225 minutes should
 * use the NORAD SDP8 "deep" space model.
 *
 * \par Quality Control.
 * I have tested this code against the NORAD FORTRAN code and against
 * published NORAD Test cases.  In both cases answers have been identical
 * to approximately 8 significant digits.  I.E. The same code on two
 * different computers gives the same answer to approximately 0.1 metre
 * resolution. 
 *
 * I have also found that the NORAD Test cases do not fully test the code.
 * In particular it did not detect a quadrant error I had in a ACTAN
 * subroutine.
 *
 * The code is implemeted as a SGP8 Class that inherits from a nxSatelliteBase
 * class.  The two methods Model_Init and UpdateECIPosition are both
 * overloaded.  As other models are developed eg SDP8 they should also
 * overload Model_Init and UpdateECIPosition.
 *
 * The UpdateECIPosition updates,
 *  - The m_time of the Update.
 *  - The ECI vector m_location of the satellite in metres
 *  -  The ECI velocity of the satellite in m/s.
 *
*/
/*---------------------------------------------------------------------------*/

class  SGP8 : public NoradSatellite
{
    private:
       double A3COF;
       double COSIO2;
       double ED;
       double EDOT;
       double GAMMA;
       int    ISIMP;
       double OVGPP;
       double OMGDT;
       double PP;
       double QQ;
       double SINI;
       double SINIO2;
       double UNM5TH;
       double UNMTH2;
       double XGDT1;
       double XHDT1;
       double XLLDOT;
       double XMDT1;
       double XNODOT;
       double XNDT;
       double XND;
       double COSI;
       double THETA2;
       double TTHMUN;
       double XNODP;
    public:
			 		SGP8();
	virtual void	Model_Init(void);
	virtual void    UpdateECIPosition( const nxTimeStamp &TNOW );
};


