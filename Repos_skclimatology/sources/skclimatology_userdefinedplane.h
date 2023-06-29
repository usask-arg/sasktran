#include <array>
/*-----------------------------------------------------------------------------
 *					skClimatology_UserDefinedPlane				2014-7-7*/
/** A two dimensional user defined climatology which varies in altitude
 *  and angle in a specified plane.  Interpolation is done bilinearly.  For points
 *  not in the plane, the angle is determined by projecting the point onto the plane.
 **/
/*---------------------------------------------------------------------------*/

class skClimatology_UserDefinedPlane : public skClimatology
{
	private:
		typedef std::vector< std::pair< CLIMATOLOGY_HANDLE, nx2dArray<double> > >::iterator Iterator;
		std::vector<std::pair< CLIMATOLOGY_HANDLE, nx2dArray<double> > >	m_profiles;
		std::vector<double>													m_heights;		//!< Heights in m
		std::vector<double>													m_angles;		//!< grid of angles measured from m_reference in plane
		std::vector<bool>													m_dologinterp;

		nxVector															m_normal;
		nxVector															m_reference;
	private:
		Iterator				IteratorToProfile				( const CLIMATOLOGY_HANDLE& species );
		double					ProjectedAngle					( const GEODETIC_INSTANT& placeandtime );
		double					InterpolateProfile				( double angle, double height, const nx2dArray<double>& profile, bool dolog );
		bool					LinearInterpIndexAndWeight		( double val, const std::vector<double>& table, std::array<size_t,2>& index, std::array<double,2>& weight, bool zeropastboundary );
		bool					IsInValidState					();
	public:
								skClimatology_UserDefinedPlane	( );
		virtual				   ~skClimatology_UserDefinedPlane	( );
		virtual	bool			UpdateCache						( const GEODETIC_INSTANT& placeandtime );
		virtual	bool			GetParameter					( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache);
		virtual bool			IsSupportedSpecies				( const CLIMATOLOGY_HANDLE& species );

		/** Sets the altitude grid used for the climatology.  This must be called before
		*   adding any species to the climatology.
		*
		*   \param heights
		*		Vector of altitudes in [m].
		*/
		bool					SetHeightGrid					( const std::vector<double>& heights );
		/** Sets the angular grid used for the climatology.  This must be called before
		 *  adding any species to the climatology.
		 *
		 *  \param angles
		 *		Angles in degrees.  An angle of 0 corresponds to the reference vector of the plane
		 *      and angles increase towards the direction of normal cross reference.
		 */
		bool					SetAngleGrid					( const std::vector<double>& angles );
		/** Defines the plane used in the climatology.  This must be called before attempting
		 *  to access any of the climatogy paramaters.
		 *
		 * \param normal
		 *		The normal vector for the plane, will be converted to a unit vector internally.
		 *
		 * \param referenceinplane
		 *		A reference vector in the plane where angles are calculated relative to,
		 *		i.e. the x-axis of the plane.
		*/ 
		bool					SetPlane						( const nxVector& normal, const nxVector& referenceinplane );
		/** Adds a species to the climatology.  This can only be called after the height and angular grids
		 *  have been specified.
		 *
		 *	\param profile
		 *		Two dimensional profile (heights,angles) for the species.
		 *
		 *	\param species
		 *		Specifier for the species to add.
		 */
		bool					AddSpecies						( const nx2dArray<double>& profile, const CLIMATOLOGY_HANDLE& species, bool dolog=false );
};