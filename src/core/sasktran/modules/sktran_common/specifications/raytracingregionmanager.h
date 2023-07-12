
class SKTRAN_RayTracingRegionManager;


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingReferencePointEstimator		2014-4-4*/
/** @ingroup rays
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayTracingReferencePointEstimator
{
	private:
		double										m_targetaltitude;							//!< A rough guide indicating the target altitude of interest, eg, 25000.0, used to set the reference point
		double										m_targetrange;								//!< A rough range indicating the range of interesting height eg 15000, ie. desired range is from 10000m to 40000m (25000-15000 to 25000+15000)
		double										m_targetquadraticweights[2];
		SKTRAN_RayTracingRegionManager*				m_regionmanager;
		SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE	m_primaryscantype;							// Stores the code which is the primary scan classification of the lines of sight

	private:
		double									LimbAltitudeWeight							( double altitude_meters );
		bool									AverageObserverLocation						(const SKTRAN_LineOfSightArray_V21& linesofsight, nxVector* referencepoint, double* mjd);
		bool									AverageWeightedLimbLocation					(const SKTRAN_LineOfSightArray_V21& linesofsight, nxVector* referencepoint, double* mjd);
		bool									AverageTargetHeightLocation					(const SKTRAN_LineOfSightArray_V21& linesofsight, nxVector* referencepoint, double* mjd);

	public:
												SKTRAN_RayTracingReferencePointEstimator	( SKTRAN_RayTracingRegionManager* manager, double targetaltitude, double targetrange);
											   ~SKTRAN_RayTracingReferencePointEstimator	();
		bool									ConfigureTargetAltitude						( double targetaltitude_meters, double targetrange_meters, double maxweightattragetaltitude);
		bool									GuessViewingGeometryMode					(const SKTRAN_LineOfSightArray_V21& linesofsight, SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE& primaryscantype);
		bool									EstimateReferencePoint						(const SKTRAN_LineOfSightArray_V21& linesofsight, nxVector* referencepoint, double* mjd);
};

/*-----------------------------------------------------------------------------
 *					class SKTRAN_RayTracingRegionManager		2010-4-27*/
/** @ingroup rays
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayTracingRegionManager
{
	private:
		nxVector								m_referencepoint;							//!< The reference point for Sasktran. 
		nxVector								m_sun;										//!< The sun location set by the user.
		double									m_sza;										//!< The solar zenith angle at the reference point
		double									m_minsza;									//!< The minimum solar zenith angle determined from all of the rays
		double									m_maxsza;									//!< The maximum solar zenith angle determined from all of the rays
		double									m_fixedearthradius;							//!< If not NaN the use this for the Earth radius (in meters).

	private:
		bool									m_isdirty;
		double									m_referencepoint_targetaltitude;			//!< The target altitude when determining the reference altitude ie either limb, nadir or zenith measurements, default 25000 meters
		double									m_referencepoint_targetrange;				//!< The altitude range above and below the traget altitude for enhanced weighting of limb viewing lines of sight, default 15000 meters.
		double									m_maxalt;									//!< The maximum altitude of the atmosphere in meters to consider when determining reference point
		double									m_minalt;									//!< the minium altitude of the atmosphere in meters used when considering reference point
		bool									m_nadir_referencepoint_onground;			//!< True if the reference point should be the average ground intersection for nadir lines of sight; false for average observer position.
		double									m_groundalt;								//!< The ground altitude in meters. Used to determine rays that hit ground
		double									m_mjd;										//!< The mjd of the measurements, set by the user
		nxGeodetic								m_geoid;									//!< The geoid used for calculations
		SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE	m_primaryscantype;						//!< Stores the code which is the primary scan classification of the lines of sight

	private:
		bool									UpdateReferencePointAndMJD					( const SKTRAN_LineOfSightArray_V21& linesofsight );
		bool									UpdateSun									();
		bool									UpdateReferencePointSZA						();
		bool									UpdateMinMaxSZA								( const SKTRAN_LineOfSightArray_V21& linesofsight);
		bool									CheckParameters								() const;
		nxVector								MoveInsideObserverOutside					(const nxVector& observer, const nxVector& look);
		void									SetDirty									() { m_isdirty = true;}
		bool									IsDirty										() const { return  m_isdirty;}
		bool									GetRayEndpointsObserverInside				( const nxVector& observer, const nxVector& look, nxVector* startpt, nxVector* endpt);
		virtual bool							GetRayEndpointsObserverOutside				( const nxVector& observer, const nxVector& look, nxVector* startpt, nxVector* endpt);

	public:
		double									GroundAlt									() const { return m_groundalt; }
		double									MaxAlt										() const { return m_maxalt; }

	public:
												SKTRAN_RayTracingRegionManager				();
		virtual								   ~SKTRAN_RayTracingRegionManager				();
		bool									MakeCoordinateSystem						( std::shared_ptr< const SKTRAN_CoordinateTransform_V2>* usercoordinates, double groundaltitude, double toa_altitude   ) const;
		nxGeodetic&								Geoid										() { return m_geoid;}
		double									MaxAltitude									() { return m_maxalt;}
		void									Clear										();
		bool									SetSun										( const nxVector& sun );					//!< Set Sun location, overrides using PlanetSun and default mjd
		bool									SetSunFromTangentPointZenAzi				( GEODETIC_INSTANT tangentpt, double solar_zenith, double solar_azimuth );

		bool									SetEarthRadius								( double earthradius_meters);												//!< Manually Set the Earth radius
		bool									SetReferencePoint							( double latitude, double longitude, double height_meters, double mjd );		//!< Sets the reference point, overrides the array of rays.
		bool									SetSZARange									( double sza, double minsza, double maxsza );
		bool									SetGroundAltitude							( double groundheight_meters );				//!< Height used for the ground (default is 0.0)
		bool									SetUpperBoundAltitude						( double upperheight_meters );				//!< Height used for upper bound of reference point calculation, default (100,000)
		bool									SetLowerBoundAltitude						( double lowerheight_meters );				//!< Height used for lower bound of reference point calculation, default (0.0)
		bool									SetNadirReferencePointOnGround				( bool onground );
		bool									SetReferencePointTargetAltitude				( double targetaltitude_meters);
		bool									SetReferencePointTargetRange				( double range_meters);
		bool									GetReferencePoint							( double* latitude, double* longitude) const;
		bool									GetSZA										( double* refptsza, double* minsza, double* maxsza) const;
		bool									GetSun										( nxVector* sun) const;
		bool									GetMJD										( double* mjd) const;
		bool									GetNadirReferencePointOnGround				( bool* onground) const;
		bool									GetPrimaryScanType							( SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE* primaryscantype) const;
		bool									ConfigureGeographicRegionManually			( double latitude, double longitude, double mjd, nxVector sun, double deltaszaDegrees);
		virtual bool							UpdateUndefinedParametersFromLinesOfSight	( const SKTRAN_LineOfSightArray_V21& linesofsight );
		bool									IsProperlyDefined							() const;
		bool									GetRayEndpoints								( const nxVector& observer, const nxVector& look, nxVector* startpt, nxVector* endpt);
};
