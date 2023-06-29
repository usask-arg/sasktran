/*-----------------------------------------------------------------------------
 *					GEODETIC_INSTANT								2005-6-24*/
/** \ingroup Climatology
 *	A structure used to represent a 3-D location is space and an instant in
 *	time. The coordinates are the geodetic corodinates in degrees and meters.
 *	This structure uses latitude and longitude in degrees,
 *  height in meters above the geoid and time is UTC expressed as a modified Julian day.  Note that longitude can be expressed as either
 *  -180 to 180 or 0 to 360 or any other circular combination. Latitude should always be in the range -90 to 90.
 **/
/*---------------------------------------------------------------------------*/

//#if !defined(GEODETIC_INSTANT_DEFINED)
class GEODETIC_INSTANT
{
	public:
		double	latitude;			//!< Geodetic latitude in degrees (-90 to +90 )
		double	longitude;			//!< Geodetic longitude in degrees (+ve east, 0 is Greenwich, -180 to 180 or 0 to 360)
		double	heightm;			//!< Height in meters.
		double	mjd;				//!< An instant in time.

	public:
		GEODETIC_INSTANT()
		{
			latitude  = -99999.0;
			longitude = -99999.0;
			heightm   = -99999.0;
			mjd       = -99999.0;
		}

		
		GEODETIC_INSTANT( double alatitude, double alongitude, double aheight, double amjd)
		{
			latitude  = alatitude; 
			longitude = alongitude;
			heightm   = aheight;
			mjd       = amjd;
		}


		bool operator== ( const GEODETIC_INSTANT& other ) const
		{
			return     ( mjd       == other.mjd)
				    && ( heightm   == other.heightm)
					&& ( latitude  == other.latitude)
					&& ( longitude == other.longitude);
		}

		void FromSequence( const double fixedarray[4])
		{
			latitude  = fixedarray[0]; 
			longitude = fixedarray[1];
			heightm   = fixedarray[2];
			mjd       = fixedarray[3];
		}

		GEODETIC_INSTANT	AsSequence() const { return *this;}
};
#define GEODETIC_INSTANT_DEFINED
//#endif
