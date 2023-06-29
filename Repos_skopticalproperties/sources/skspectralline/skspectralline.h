
#include <complex>

/*
#if defined(_MSC_VER) && !defined(NXSKSPECTRALLINE_CORE_LIB)
	#pragma message("Automatically linking the skSpectralLine library into this project")
	#if defined (NXDEBUG)
		#pragma comment( lib, "skSpectralLine_Debug")
	#else
		#pragma comment( lib, "skSpectralLine_Release")
	#endif
#else
	#pragma message("NOT Automatically linking the skSpectralLine version 2.1 library into this project")
#endif
*/

class skSpectralLineCollection;
class skSpectralLineShape;

/*-----------------------------------------------------------------------------
 *					skSpectralLine		2013-3-6*/
/** \ingroup spectralline
 *	A base class that reprsents a single transition line in a atom or molecule.
 *	Derived classed may use measured databases to calculate the transition line
 *	parameters, e.g. HITRAN, Cologne, JPL. Alternatively derived classes may
 *	use theoretical calculations to provide line parameters.
 *
 **/
/*---------------------------------------------------------------------------*/

class skSpectralLine : public nxUnknown
{
	private:
		const skSpectralLineCollection*				m_parentmolecule;
		double										m_currentT;
		double										m_current_line_intensity;

	public:
		static size_t								g_numinstances;

	protected:
		void										SetCurrentLineIntensity(  double line_intensity, double T) { m_current_line_intensity = line_intensity; m_currentT = T;}

	public:
		bool										SetParentMolecule		( const skSpectralLineCollection* parent );
		const skSpectralLineCollection*				ParentMolecule			( ) const				{ return m_parentmolecule;}
		double										LineIntensity			( ) const				{ return m_current_line_intensity;}

	public:
													skSpectralLine			( );
		virtual									   ~skSpectralLine			( );
		virtual bool								CalculateLineIntensity	( double T );
		virtual double								Snm						( )	const	= 0;					//!< Spectral Line Intensity from level n to m, at reference temperature. Same value as Hitran database [cm-1/(molecule cm-2)].
		virtual double								Nu						( )	const	= 0;					//!< The spectral line transition frequency [cm-1] from level n to m (in vacuum)  
		virtual double								GammaAir				( )	const	= 0;					//!< Air  Broadened half width HWHM cm-1/atm at reference temperature 
		virtual double								GammaSelf				( )	const	= 0;					//!< Self Broadened half width HWHM cm-1/atm at reference temperature
		virtual double								EinsteinA				( )	const	= 0;					//!< Einstein A coefficient 
		virtual double								ELower					( )	const	= 0;					//!< Lower State energy in cm-1
		virtual double								Nair					( )	const	= 0;					//!< Coefficient of temperature dependence of air broadened half width.
		virtual double								Deltaair				( )	const	= 0;					//!< Air Broadened pressure shift of line transition  in cm-1/atm at reference temeprature
		virtual double								Tref					( )	const	= 0;					//!< Reference Temperature in K for the spectral line parameters (often 296K).
		virtual double								EUpper					( ) const { return ELower() + Nu();}
};

/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer		2013-3-6*/
/** \ingroup spectrallineinternals
 *	Each spectral line shape is given the opportunity to save some storage
 *	space for each spectral line it is expected to process. This can be useful
 *	if we have to calculate the Voigt profile at a large number of
 *	wave numbers for a large number of lines as we can cache many of the Voigt
 *	parameters for each spectral line.
 *
 *	The trick to making code legible is to create the base class with almost no
 *	functionality, it is simply a placeholder for derived classes. We then make
 *	a templated derived class skSpectralLineShapeStorageBuffer_Type<T> where
 *	type T is the buffer object. This technique provides a lot of flexibility
 *	each line shape object has two components, 1) the interface class that
 *	most users see and a worker class, T, that holds all of the algorithm details.
 *	
 **/
/*---------------------------------------------------------------------------*/

class skSpectralLineShapeStorageBuffer : public nxUnknown
{
	public:
		virtual			~skSpectralLineShapeStorageBuffer(){};
};


/*-----------------------------------------------------------------------------
 *					class skSpectralLineEntry						2013-3-7*/
/** \ingroup spectrallineinternals
 *	A small lightweight class used to hold all the information used to
 *	calculate spectra from line spectra. The class holds a pointer to a
 *	spectral line entry (from Hitran for instance), a pointer to a line shape
 *	object, Voigt for example and an optional work/storage buffer. The line entry
 *	is one of several (or indeed many) lines held by a molecule. The entry, becuase
 *	it is lightweight can be easily passed around.
 *
 *	The entry is especially lightweight if the lineshape object is set to NULL as there
 *	is no requirement to make a storage buffer. Thus higher level classes that load
 *	all the lines for a molecule, which can be many thousands, should not define line shape objects
 *	(unless they dont mind the memory allocation hit). Line shape objects should be set
 *	after lines have been selected into smaller groups of micro-windows.
 */
/*---------------------------------------------------------------------------*/

class skSpectralLineEntry : public nxUnknown
{
	private:
		skSpectralLine*						m_spectralline;			// Pointer to the spectral line of interest, eg pointer to Hitran entry, Cologne entry, JPL entry, theoretical entry
		skSpectralLineShape*				m_lineshapeobject;		// Pointer to the spectral line shape object (eg Voigt, Doppler, Lorenz, user defined etc.)
		skSpectralLineShapeStorageBuffer*	m_storagebuffer;		// Pointer to any storage area the line shape obejct needs (may be NULL)

	public:
											skSpectralLineEntry					( );
											skSpectralLineEntry					( const skSpectralLineEntry& other );
										   ~skSpectralLineEntry					( );
		bool								SetSpectralLine						( skSpectralLine* spectralline   );
		const skSpectralLine*				SpectralLine						() const { return m_spectralline;}
		bool								SetLineShapeObject					( skSpectralLineShape* lineshapeobject );		
		bool								ConfigureLineParameters				( double temperature, double pressure, const GEODETIC_INSTANT& geopt, skClimatology*	 atmopshericstate);
		double								Nu00								( )	const									{ return (m_spectralline != NULL) ? m_spectralline->Nu() : 0.0;}
		bool								operator <							( const skSpectralLineEntry& other ) const  { return Nu00() < other.Nu00();}

	public:
		virtual bool						SetLineLimitsfromMaxLineStrength			( double parentmaxlinestrength);
		virtual bool						SetTolerance								( double tolerance);
		virtual bool						AbsorptionCrossSectionOrEmission			( double nu, double *absxsec ) const;
		virtual bool						AddAbsorptionCrossSectionOrEmissionArray	( const std::vector<double>& nu, std::vector<double>* absxsec );

};

/*-----------------------------------------------------------------------------
 *					skSpectralLine_ParentMolecule					2013-3-6*/
/**	\ingroup spectralline
 *	A class used to represent all the lines (of interest) in a molecule.
 *	Different instances of molecule should be created for different isotopes
 *	as they may have different masses, partition functions etc.
 *
 *	The Molecule is really just a collection of spectral lines used to
 *	generate a spectrum in a region of interest. We have built the molecule class
 *	so the spectral lines use a collection of spectral line entry. Each entry
 *	is a lightweight structure consisting of 3 pointers. This means we should
*	not experience significant problem moving spectral line entries from one
*	molecule instance to another molecule instance, this will be useful for
*	micro-window regions of interest from larger sets of data.
 */
/*---------------------------------------------------------------------------*/

class skSpectralLineCollection
{
	private:
		skSpectralLineShape*									m_lineshape;				// The current line shape object. 
		double													m_maxlinestrength;			// The maximum line strength in this collection of lines
		double													m_microwindow_minwavenum;
		double													m_microwindow_maxwavenum;
		std::vector< skSpectralLineEntry* >						m_lines;					// An array of spectral lines.
		std::vector< skSpectralLineEntry* >						m_linesforthreadaccess;		// Storage for the threads. This makes the OMP loop much easier to code up.
		std::vector< std::vector<double>  >						m_threadstorage;			// Storage for the threads. Stores the cross-section from each thread

	public:
		typedef std::vector< skSpectralLineEntry* >::iterator		iterator;
		typedef std::vector< skSpectralLineEntry* >::const_iterator const_iterator;
		iterator												begin()				{ return m_lines.begin();}
		const_iterator											begin() const		{ return m_lines.begin();}
		iterator												end  ()				{ return m_lines.end();}
		const_iterator											end  ()	const 		{ return m_lines.end();}
		size_t													size () const		{ return m_lines.size();}

		static size_t											g_numinstances;

	private:
		bool													AddAbsorptionCrossSectionOrEmissionArraySingleThreads ( const std::vector<double>& nu, std::vector<double>* absxsec ); // Single threaded version
		bool													AddAbsorptionCrossSectionOrEmissionArrayMultiThreads  ( const std::vector<double>& nu, std::vector<double>* absxsec ); // Single threaded version
		bool													CheckNumThreads								();

	public:
																skSpectralLineCollection( double microwindow_min, double microwindow_max );
		virtual												   ~skSpectralLineCollection();
		bool													SetLineShapeObject						( skSpectralLineShape* lineshapeobject );
		double													MicroWindow_MinWavenum					() const { return m_microwindow_minwavenum;}
		double													MicroWindow_MaxWavenum					() const { return m_microwindow_maxwavenum;}
		skSpectralLineShape*									LineShape								() { return m_lineshape;}
		const std::vector< skSpectralLineEntry*  >&				SpectralLines							() const { return m_lines;}
		void													ClearLines								(size_t reservesize);
		bool													AddEntry								( skSpectralLineEntry* entry );
		bool													AbsorptionCrossSectionOrEmission		( double nu, double* absxsec ) const;
		bool													AddAbsorptionCrossSectionOrEmissionArray( const std::vector<double>& nu, std::vector<double>* absxsec );
		bool													SetNumThreads							( size_t numthreads);

	public:
		virtual double											MaxLineStrength							( ) const {return m_maxlinestrength;}			// get the max line strength from this set of lines
		virtual bool											SetLineLimitsfromMaxLineStrength		( double parentmaxlinestrength);
		virtual bool											SetLineTolerance						( double tolerance );
		virtual bool											UpdateLocation							( double temperature, double pressure, const GEODETIC_INSTANT& geopt, skClimatology* atmosphere );

	public:
		
		/** Return partition function Q(T) for the lower level of the spectral line passed in. 
			By design each spectral line collection will only have one internal partition function.
		 */
		virtual double											QPartition			( double T ) const  = 0;	

		/** Return the mass of the molecule in AMU. By design, each spectral line collection
		 *	is only allowed one mass. Thus different isotopes will be implemented as different
		 *	collections.	
		 */
		virtual double											MassAMU				( ) const	 = 0;	

		/** Return the partial pressure of the "molecule" that represents the collection
		 *	This might be used when calculating slef-broadening line width parameters
		 *	in the line shape object. The default implementation returns 0.0 but other implementations
		 *	might query the climatology passed into SetAtmosphericState
		 */

		virtual double											PartialPressure		( const GEODETIC_INSTANT& geopt, double airpressure, double airtemp) const  = 0;	// Return the partial pressure of this molecule, you are given a location, airpressure (pascals) and temperature (K) to assist with calculation. 


};


/*-----------------------------------------------------------------------------
 *					skSpectralLineShapeStorageBuffer_Type		2013-3-6*/
/** \ingroup spectrallineinternals
 *	A small templated class to fetch the lineshape object's storage buffer information
 *	from a base-class skSpectralLineShapeStorageBuffer object. This
 *	avoids excessive use of void* to pass storage pointers aroudn and keeps
 *	all type checking on the up and up. 
 *
 *	The line shape class will typically allocate a new instance  of
 *	skSpectralLineShapeStorageBuffer_Type using the storage class or structure of
 *	their  choice as the template. This will automatically create an instance
 *	inside the m_cache member. This typically occurs when the line-shape
 *	object is asked to Create a storage buffer.
 *
 *	The line shape object can retrieve a pointer to their cached storage class or
 *	structure by calling skSpectralLineShapeStorageBuffer_Type<T>::Cache() passing in
 *	the pointer to the base class. Note the function is static and does not require
 *	an instance to be used.
 **/
/*---------------------------------------------------------------------------*/

template <class T>
class skSpectralLineShapeStorageBuffer_Type : public skSpectralLineShapeStorageBuffer
{
	public:
		T				m_cache;

	public:
		virtual		   ~skSpectralLineShapeStorageBuffer_Type() {}

		static T*		Cache( skSpectralLineShapeStorageBuffer* baseclass )
						{
							skSpectralLineShapeStorageBuffer_Type<T>* derivedclass;

							derivedclass = dynamic_cast<skSpectralLineShapeStorageBuffer_Type<T>* >( baseclass);
							return &derivedclass->m_cache;
						}
};

/*-----------------------------------------------------------------------------
 *					skSpectralLineShape		2013-3-6*/
/** \ingroup spectralline
 *	A base class used to implement line shape functions. Many of the
 *	derived classes will implement Voigt, Lorenz or Gauss but many other
 *	options, like user defined shapes can be added at later dates.
 *
 *	\par Multi-threading
 *	This class will often form the work horse of HITRAN calculations and it
 *	will be critical that speed is fully optimized in the class. Part of that equation
 *	is that the LineShapeFunction method must be thread safe so we can have multiple
 *	threads running simultaneously. To help develop multi-threading code the skSpectralLineShape
 *	class implements the line-shape function as const, implying it is intrinsically thread safe.
 *	In reality this is not quite true as the method request that the user pass in the class
 *	specific storage buffer created by CreateStorageBuffer. This buffer will have all of the
 *	information that can be cached one we know atmospheric state and it will also hold any
 *	internal structures (eg lookup tables) that must be changed during processing of method
 *	LineShapeFunction. The buffer code must ensure that changes to internal structures within this
 *	buffer are properly protected with mutexes etc.
 **/
/*---------------------------------------------------------------------------*/

class skSpectralLineShape : public nxUnknown
{
	public:

		/** [THREAD-SAFE] Calculate the Line Shape Function for a given spectral line. This
		 *	interface is for single wave-number processing and is suitable for the Monte-Carlo 
		 *	and Successive Order engines. The controlling software must successfully called 
		 *	ConfigureAtmopshericState and ConfigureStorageBuffer before calling LineShapeFunction.
		 *	Typically, atmospheric state information, which is common to all of the
		 *	spectral lines processed, is stored within the derived SpectralLineShape object while
		 *	line specific information for the given atmospheric state is stored inside the
		 *	storage buffer. This function must be implemented so it is thread safe as multiple threads
		 *	may be calculating LineShape for different wavenumbers or spectral lines simultaneaously
		 *
		 *	\param nu
		 *		The wavenumber at which to calculate the lineshape. We dont specify vacuum or air but it should be the same as the
		 *		optical properties. HITRAN cross-sections are specified in vacuum. Many other cross-sections, especially in the visible, are specified at STP.
		 *
		 *	\param uservalue
		 *		returns the line shape value at this wavenumber. The line shape is normalized so its total area is unity.
		 *
		 *	\param spectralline
		 *		The spectral line (eg HITRAN line) object that is "using" this LineShape object. Note that one instance of skSpectralLineShape may
		 *		be shared between many different spectral lines.
		 *
		 *	\param storagebuffer
		 *		An object created by this spectral line shape object that is used to cache information about the specific
		 *		spectral line being processed. Access to the cached information must be thread safe using mutexes etc as it is
		 *		possible that multiple threads processing multiple wavenumbers are accessing the cached object simultaneously.
		 */

		virtual bool		LineShapeFunction			(	double								nu,						//!< Calculate line shape function at this wavenumber
															double*								uservalue,				//!< returns the line shape value in this variable
															const skSpectralLine*				spectralline,			//!< Spectral line for line shape function
															skSpectralLineShapeStorageBuffer*	storagebuffer ) = 0;	//!< Storage buffer from an earlier call to ConfigureStorageBuffer

		/** [THREAD-SAFE] Allows the line shape object to update the useful wavenumber limits for
		 *	future Voigt/lineshape function evaluations. This can be a signifcant processing
		 *	optimization. This function is usually invoked after the user has moved to a different point
		 *	in the atmosphere and [pressure and temperature hve changed.
		 *
		 *	\param parentmaxlinestrength
		 *		The maximum line strength considered for this micro-window calculation. The parent classes have passed in the "brightest"
		 *		line in the microwindow. Given this brightness we can limit the range of wavenumbers that we need to evaluate this spectral line
		 *		and still keep the calculation within  a certain threshold/tolerance.
		 *
		 *	\param spectralline
		 *		The spectral line (eg HITRAN line) object that is "using" this LineShape object. Note that one instance of skSpectralLineShape may
		 *		be shared between many different spectral lines.
		 *
		 *	\param storagebuffer
		 *		An object created by this spectral line shape object that is used to cache information about the specific
		 *		spectral line being processed. Access to the cached information must be thread safe using mutexes etc as it is
		 *		possible that multiple threads processing multiple wavenumbers are accessing the cached object simultaneously.
		 */
		virtual bool		SetParentMaxLineStrength	(	double parentmaxlinestrength,
															const skSpectralLine*				spectralline,			//!< Spectral line for line shape function
															skSpectralLineShapeStorageBuffer*	storagebuffer ) const = 0;	//!< Storage buffer from an earlier call to ConfigureStorageBuffer

		/** [THREAD-SAFE] Allows the user to manually set the refection tolerance for this line shape object to set update the useful wavenumber limits for
		 *	future Voigt/lineshape function evaluations. This can be a signifcant processing optimization. This function is usually invoked before the user 
		 *	moves to a new point in the atmosphere and pressure and temperature have changed.
		 *
		 *	\param tolerance
		 *		The tolerance of the minimum accepted line compared to the Max Strength Line. The default tolernace is typically 1.0E-09
		 *
		 *	\param spectralline
		 *		The spectral line (eg HITRAN line) object that is "using" this LineShape object. Note that one instance of skSpectralLineShape may
		 *		be shared between many different spectral lines.
		 *
		 *	\param storagebuffer
		 *		An object created by this spectral line shape object that is used to cache information about the specific
		 *		spectral line being processed. Access to the cached information must be thread safe using mutexes etc as it is
		 *		possible that multiple threads processing multiple wavenumbers are accessing the cached object simultaneously.
		 */
		virtual bool		SetTolerance				( double tolerance, 
														  const skSpectralLine* spectralline, 
														  skSpectralLineShapeStorageBuffer* storagebuffer )  = 0;


		/** [THREAD-SAFE] Calculate the Line Shape Function for an array of wavenumbers. Returns the
		 *	lineshape in array uservalue. This function replicates function #LineShapeFunction except this
		 *	method provides a derived class the ability to optimize the line shape code across an array of
		 *	wavenumbers. This can be a significant saving when finding the X,Y regions of the Voigt function
		 *	for example. Keep in mind that one line shape object is shared between many spectral lines and that the code
		 *	will be possibly running in multiple threads processing muliple spectral lines simultaneously. A default
		 *	implementation is provided which simply loops and calls LineShapeFunction.
		 *
		 *	\param nu
		 *		The array of wavenumbers at which to calculate the lineshape. We dont specify vacuum or air but it should be the same as the
		 *		optical properties. HITRAN cross-sections are specified in vacuum. Many other cross-sections, especially in the visible, are specified at STP.
		 *
		 *	\param uservalue
		 *		returns the array of line shape value. It will be same size as "nu". Note we only use the resize method to avoid unnecessary memory allocations.
		 *		The line shape is normalized so its total area is unity.
		 *
		 *	\param spectralline
		 *		The spectral line (eg HITRAN line) object that is "using" this LineShape object. Note that one instance of skSpectralLineShape may
		 *		be shared between many different spectral lines.
		 *
		 *	\param storagebuffer
		 *		An object created by this spectral line shape object that is used to cache information about the specific
		 *		spectral line being processed. Access to the cached information must be thread safe using mutexes etc as it is
		 *		possible that multiple threads processing multiple wavenumbers are accessing the cached object simultaneously.
		 */

		virtual bool		AddLineShapeFunctionArray	(	const std::vector<double>&			nu,							//!< Calculate line shape function at this wavenumber
															std::vector<double>*				uservalue,					//!< returns the line shape value in this variable
															const skSpectralLine*				spectralline,				//!< Spectral line for line shape function
															skSpectralLineShapeStorageBuffer*	storagebuffer );			//!< Storage buffer from an earlier call to ConfigureStorageBuffer

		/** 
		 *  [NOT THREAD-SAFE] Extract all of the necessary line parameters from the atmospheric state and spectral line.
		 *	If necessary the information can be stored in the storage buffer and this storage buffer
		 *	instance will be passed to future calls to LineShapeFuction.
		 *	We would expect that just about all derived classes will extract temperature, SKCLIMATOLOGY_TEMPERATURE_K,
		 *	and most will extract pressure, SKCLIMATOLOGY_PRESSURE_PA from atmosphericstate. Derived classes will probably
		 *	extract relevant spectral line info from spectral line and possibly generate intermediate results which
		 *	are then stored in the storage buffer. This function does not need to be thread safe in the sense that
		 *	
		 */

		virtual bool		ConfigureLineParameters		(	const skSpectralLine*				spectralline,			//!< The line parameters
															double								temperature,			//!< The temperature at geopt (as an optimization)
															double								pressure,				//!< The pressure at geopt (as an optimization).
															const GEODETIC_INSTANT&				geopt,					//!< Get atmopsheric state at this point
															skClimatology*						atmopshericstate,		//!< Use this climatology to get atmopsheric state
															skSpectralLineShapeStorageBuffer*	storagebuffer ) = 0;	//!< and store it in this buffer.

		/** [NOT THREAD-SAFE] Function to create the line shape storage buffer for the spectral line. All storage buffers must
		 *	derive from skSpectralLineShapeStorageBuffer, we provide templatec class skSpectralLineShapeStorageBuffer_Type
		 *	for this purpose but users can provide their own solution if desired.
		 *	If the creation is valid then the pointer to the storage buffer object is returned in storage buffer
		 *	with a reference count of 1. The calling code must release the storage buffer with a call to "Release".
		 *	This function does not need to be thread-safe.
		*/

		virtual bool		CreateStorageBuffer			(	skSpectralLineShapeStorageBuffer**	storagebuffer ) = 0;	//!< and store it in this buffer.
};

#include "skspectrallineshape_voigtkuntz.h"
#include "skspectrallineshape_humlicekwells.h"
#include "skspectrallineshape_voigttabulated.h"
#include "hitran/hitran_xs_cache.h"
#include "hitran/hitran_partition_cache.h"
#include "hitran/skspectralline_hitran.h"
#include "hitran/hitran_spectrallineio.h"

