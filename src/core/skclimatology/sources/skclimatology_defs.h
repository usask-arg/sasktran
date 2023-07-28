#include <nxbase_core.h>
#include <nxbase_geodesy.h>


/** \ingroup skClimatology
  * \name Outline
	\par Introduction
    This is a C++ library that provides a standard software interface to 
	geophysical climatologies. The intent of the class is to provide other software
	packages with a standard interface to geophysical parameters and is especially
	focussed on providing a height profile of parameters at a single location.

	When we originally developed our Sasktran Radiative Transfer Model we looked
	at ways to inject climatological data into the mode. When we looked around we saw that the
	majority of other models specified atmospheric parameters as a series of regular
	grid array tables. This is ok and it defintely works but does leave the end-user 
	with the burden of reading and decoding all of the different data formats. All
	of the reading and decoding eventually becomes so annoying, extended and compilcated
	that it overwhelms the radiative transfer model development. We developed the skClimatlogy
	interface to at least try and simplify the radiative transfer models interface to
	the data. As it happens the resultant library is also very useful in other applications.

	Over the years we have built up a selection of climatology classes that meet
	our needs and we make these freely available to anyone who chooses to use them. However
	the real advantage for our group has been the simplification of software design for our
	radiative transfer models. In that light, we now view the climatology problem as being not
	one of looking for existing code libraries but rather development of a new C++ class derived 
	from skClimatology which, for some strange reason, never seems that overwhelming. We think
	our skClimatology interface provides a simple yet robust interface which is well suited
	for our needs as experimental atmospheric physicists.
	

*/

/*-----------------------------------------------------------------------------
 *					class skClimatology								2005-6-24*/
/** \ingroup skclimatology
 *	The skClimatology base interface
 **/
/*---------------------------------------------------------------------------*/

class skClimatology : public nxUnknown
{
	private:
		bool					m_cacheisloaded;
		double					m_missingvalue;

	public:
								skClimatology	();
		virtual				   ~skClimatology	() {}
//		bool					DeepCopy		( const skClimatology& other );
		void					SetCacheIsLoaded( bool isloaded)			{ m_cacheisloaded = isloaded;}
		bool					CacheIsLoaded   ( ) const 					{ return m_cacheisloaded;}

		/** Get the value used for data that missing in the model.  This is the value returned by
		 *	#GetParameter and #GetHeightProfile when there is an error retrieving model values.
		 */
		double					MissingValue	() const					{ return m_missingvalue;}

		/** Set the value returned by #GetParameter and #GetHeightProfile for missing data
		 */

		bool					SetMissingValue	( double missingvalue)		{ m_missingvalue = missingvalue; return nxTRUE;}

		/** Checks to see if the cache is valid. The cache is always valid after a successful
		 *	call to UpdateCache
	     *
		 *	\return
		 *	Returns true if the cache was appropriately updated for the requested location.
		 *
		 */

		bool					CheckCache		();

	public:

		/** Update the models cache at the specified location and instant. All models must update
		 *	their cache before they can retrieve values.  This function can be automatically invoked
		 *	by #GetParameter and GetHeightProfile.
		 *
		 *	\return
		 *	Returns true if the cache was appropriately updated for the requested location.
		 *
		 */

		virtual	bool			UpdateCache			 ( const GEODETIC_INSTANT& placeandtime) = 0;

		/** Get the requested species at the specified location. If the model cannot retrieve a valid value for the
		 *	species at the specified location it will return false and set the return value to the missing value.\n
		 *
		 * \param species
		 *	The species to be retrieved
		 *
		 * \param placeandtime
		 *	The place and time at which the species is required
		 *
		 * \param value
		 *	Write the value of the species into the location poointed to by this variable.  It will return
		 *	\i MissingValue if there is a problem
		 *
		 *	\param updatecache
		 *	If true then call #UpdateCache with this placeandlocation before evaluating the species.
		 *
		 *  \return
		 *	Returns true if the species was successfully retrieved.
		*/

		virtual	bool			GetParameter		 ( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, double* value, bool updatecache) = 0;

		/** Get the a height profile of the requested species at the specified location. This is similar to #GetParameter
		 *	except it fetches a height profile at the user specified resolution.  In addition it returns the number of 
		 *	altitudes that are bad (set to \i MissingValue ).  A return value of 0 means all altitudes
		 *	were properly evaluated.  Dysfunctional conditions (like the number of heights equal to 0) will return a value of 0.
		 *	\n
		 *
		 *	\param species
		 *	The species to be retrieved
		 *
		 *	\param placeandtime
		 *	The place and time at which the species is required.  The heightm field of this structure
		 *	is used if the cach is updated but is ignored for actual evaluation of species.
		 *
		 *	\param geometricheightm
		 *	An array of heights (in meters) at which the requested species is required.
		 *
		 *	\param value
		 *	A pointer to a 1d array.  Upon return the array will contain the value of the species at the
		 *	specified heights.  Any heights that could not be evaluated by the model will set the corresponding
		 *	entry in value to \i MissingValue()
		 *
		 *	\param updatecache
		 *	If true then call #UpdateCache with this placeandlocation before evaluating the species.
		 *
		 *  \return
		 *	Returns the number of species that were not successfully retrieved. I.E. a value of 0 means
		 *	all altitudes were successfully retrieved.

		*/
		virtual	bool			GetHeightProfile	 ( const CLIMATOLOGY_HANDLE& species,  const GEODETIC_INSTANT& placeandtime, const double* geometricheightm, int numheights, double* value, bool updatecache, size_t* numbad);

		/** Query the model to see if it supports a specific species.
		 *\return
			Returns true if the model supports the requested species.
		*/

		virtual bool			IsSupportedSpecies	 ( const CLIMATOLOGY_HANDLE& species ) = 0;
		
		/** Create a clone of this climatology.  This is a DEPRECATED function. It is only used
		 *	for legacy support of SASKTRAN V2. This is the only way to support skClimatology
		 *	in a multithreaded environment: create a clone of a base instance and pass the clone
		 *	to the multithreaded environment. The user must call Release() on the climatology
		 *	passed back (from any thread) once they have finished with the object. Note that the
		 *	clone is a completely separate instance of skClimatology and is not (generally) linked or synced
		 *	with the parent after its creation.
		*/

//		virtual bool			CreateClone			 (skClimatology** clone)	const { throw "You should not be calling the SkClimatology base class CreateClone method"; return false;} 
};




#include "skclimatology_userdefinedtable.h"
#include "skclimatology_userdefined_latlon_table.h"
#include "skclimatology_userdefinedplane.h"
#include "skclimatology_textfiledbase.h"
#include "skclimatology_zero.h"
#include "skclimatology_one.h"
#include "skclimatology_msis90.h"
#include "skclimatology_labowozone.h"
#include "skclimatology_pratmo.h"
#include "skclimatology_linearcombo.h"
