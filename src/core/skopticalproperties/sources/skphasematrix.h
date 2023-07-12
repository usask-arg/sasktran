

typedef double SKRTFLOAT;

/*---------------------------------------------------------------------------
 *					class skRTStokesVector					2003-12-8  */
/**	\ingroup Radiative_Transfer
 *	\par Overview
 *	A light, low-overhead, class that represents the 4 element, Stokes
 *	vector [I,Q,U,V] using a 4 element array of #SKRTFLOAT.  It is really
 *	a few indexing methods and overloaded operators acting on the 4 element storage.
 *  We have purposely kept the class minimalistic as radiative transfer codes
 *	generate literally millions of these matrices and total memory burden
 *	becomes an issue.

 *  \par Indexing
 *	The class implements 1 based indexing so comparisons with mathematical
 *	theory is much easier. In summary \code
 *
 *  I = stokesvector.At(1);
 *  Q = stokesvector.At(2);
 *  U = stokesvector.At(3);
 *  V = stokesvector.At(4);
 *
 *  \endcode
 *  Note that offsets from #ArrayBasePtr() are zero-based.
 *
 *	\par Storage and Iteration
 *	The class uses local storage for the array and should be quite speed efficient. Each element
 *	of the array is set to 0.0 during construction. The class provides \e begin and \e end methods
 *	for iterating over the contents.
 *	simplify
 *
 */
/*-------------------------------------------------------------------------*/

class skRTStokesVector
{
	private:
		SKRTFLOAT				m_phasematrixstorage[4];											//!< Local Storage for the m_phasematrix

	public:
		typedef	SKRTFLOAT*		iterator;															//!< Used to iterate over the Stokes vector

	public:
								skRTStokesVector	();												//!< Default constructor creates a blank stokes vector
								skRTStokesVector	( const skRTStokesVector& other );				//!< Copy constructor. Copies the other array over this one
		skRTStokesVector&		operator =			( const skRTStokesVector& other );				//!< Assigns the other stokes vector to this one
		void					SetTo				( float  val );									//!< Sets the stokes vector to value \e val
		void					SetTo				( double val );									//!< Sets the stokes vector to value \e val
		void					SetTo				( double I, double Q, double U, double V);		//!< Sets the stokes vector to value \e val
		void					SetTo				( double* IQUV);								//!< Sets the Stokes vector to the 4 values I,Q, U, V.  Pointer must be at least 4 elements
		SKRTFLOAT*				ArrayBasePtr		( )	{ return &m_phasematrixstorage[0];}			//!< Returns a pointer to the start of the 4 elements of this Stokes vector
		iterator				begin				( )	{ return &m_phasematrixstorage[0]; }		//!< Returns an iterator to the start of this vector
		iterator				end					( )	{ return &m_phasematrixstorage[N_ELEMENTS(m_phasematrixstorage)];}	//!< Gets an iterator to the point just beyond our Stokes vector
		skRTStokesVector&		operator *=			( double val);									//!< Multiplies the stokes vector by this value
		skRTStokesVector&		operator *=			( float  val);									//!< Multiplies the stokes vector by this value
		skRTStokesVector&		operator +=			(const skRTStokesVector& other);				//!< adds the other vector to this one on an element by element basis
		skRTStokesVector		operator +			(const skRTStokesVector& other) const;				//!< Adds two Stokes vectors together. Convenient but not optimal.
		skRTStokesVector		operator -			(const skRTStokesVector& other) const;				//!< Subtracts two Stokes vectors. Convenient but not optimal.
		skRTStokesVector		operator *			(double value) const;									//!< returns this Stokes vector multiplied by a double precision number. Convenient but not optimal.
		skRTStokesVector		operator *			(float value) const;									//!< Returns this Stokes vector multiplies by a floating point number. Convenient but not optimal.
		void                    RotatePolarPlaneThru( double cosEta, double sinEta );                       //!< Rotate the basis for the electric field components used in this vector through angle eta (Tony's thesis page 44)
#if defined(NXDEBUG)
		SKRTFLOAT&				At					(int idx);										//!< Returns a reference to the element at index \e idx (1-4)
		SKRTFLOAT				At					(int idx) const;								//!< Returns a copy of the element at index \e idx (1-4)
#else
		inline SKRTFLOAT&		At					(int idx)		{ return m_phasematrixstorage[ --idx ];}	//!< Returns a reference to the element at index \e idx (1-4)
		inline SKRTFLOAT		At					(int idx) const { return m_phasematrixstorage[ --idx ];}	//!< Returns a copy of the element at index \e idx (1-4)
#endif

	template<typename radtype>
	static void SetToZero(radtype& a);

	template<typename radtype>
	static bool IsNegative(const radtype& a);

	template<typename radtype>
	static void SetNegativesToZero(radtype& a);

	template<typename radtype>
	static void SetNansToZero(radtype& a);
};



/*--------------------------------------------------------------------------
 *					class skRTPhaseMatrix								  */
/**	\ingroup skphasemat
 *	A light, low-overhead, class that represents a 4x4 scattering phase matrix.  It is really
 *	a  a few indexing methods and overloaded operators acting on the 4x4 matrix storage.
 *  We have purposely kept the class minimalistic as the radiative transfer codes
 *	generate literally millions of these matrices and total memory burden becomes an issue.  The
 *	memory burden is why the class is not derived from #nx2dArray which has about a 100 byte overhead.
 *
 *	\par Storage and Layout
 *	The class uses local storage for the array of 16 elements of #SKRTFLOAT.  I
 *	recommend leaving #SKRTFLOAT as type \e float to keep the memory burden to a minimum
 *	in radiative transfer codes. The class uses column major storage, similar to fortran,
 *	and the #At methods access the array using 1 based indexing identical to the notation
 *	in the mathematical theory. In summary \code
 *
 *  P11 = phasematrix.At(1,1); ...  P14 = phasematrix.At(1,4);
 *  P21 = phasematrix.At(2,1); ...  P24 = phasematrix.At(2,4);
 *  P31 = phasematrix.At(3,1); ...  P34 = phasematrix.At(3,4);
 *  P41 = phasematrix.At(4,1); ...  P44 = phasematrix.At(4,4);
 *
 *	\endcode
 *
 *	\par Phase Matrix Normalization
 *	We use the same normalization as Hansen and Travis: the phase
 *	matrix when integrafted over \f$4\pi\f$ solid angle gives \f$4\pi\f$.
 *
 *	\par Iteration
 *	The class provides \e begin and \e end methods for iterating over the
 *	contents of the phase matrix as well as a couple of operator overloads to
 *	simplify
 */
/*-------------------------------------------------------------------------*/


class skRTPhaseMatrix
{
	private:
		SKRTFLOAT				m_phasematrixstorage[4*4];												/*!< Local Storage for the m_phasematrix */

	public:
		typedef	      SKRTFLOAT*	iterator;															/*!< Used to iterate over the phase matrix. Very fast and efficient */
		typedef const SKRTFLOAT*	const_iterator;														/*!< Used to iterate over the phase matrix and guaranteed not to modify it. Very fast and efficient */

	public:
							skRTPhaseMatrix();															/*!< Create a blank instance */
							skRTPhaseMatrix( const skRTPhaseMatrix& other );							/*!< Copy the other phase matri to this instance */
		SKRTFLOAT&			At				( int row, int col);										/*!< Returns a reference to the element index by <em>(row,col)</em> using a 1 based nomenclature*/
		SKRTFLOAT			At				( int row, int col) const;									/*!< Returns a reference to the element index by <em>(row,col)</em> using a 1 based nomenclature*/
		SKRTFLOAT*			ArrayBasePtr	()					{return &m_phasematrixstorage[0];}		/*!< Returns a pointer to the start of the 4x4 matrix */
		iterator			begin			()					{return &m_phasematrixstorage[0];}	 	/*!< Get the iterator pointing to the first element */
		const_iterator      begin           () const            {return &m_phasematrixstorage[0];}		/*!< Get a const version of the iterator pointing to the first element */
		iterator			end			    ()					{return &m_phasematrixstorage[N_ELEMENTS(m_phasematrixstorage)];} 	/*!< Get the iterator pointing to the element just beyond this phase matrix*/
		const_iterator      end             () const            {return &m_phasematrixstorage[N_ELEMENTS(m_phasematrixstorage)];}	/*!< Get a const version of the iterator pointing to the element just beyond this phase matrix*/
		void				SetTo			( double val );												/*!< Set all of the elements of the phase matix to \e val */
		void				SetTo			( float val );												/*!< Set all of the elements of the phase matix to \e val */
		skRTPhaseMatrix&    RMultBy         ( const skRTPhaseMatrix& other );                           /*!< This becomes the result of right multiplying itself by #other */
		skRTPhaseMatrix&    LMultBy         ( const skRTPhaseMatrix& other );                           /*!< This becomes the result of left multiplying itself by #other */
		skRTPhaseMatrix&	operator +=     ( const skRTPhaseMatrix& other);							/*!< Add the other phase matrix to this one */
		skRTPhaseMatrix&	operator *=     ( double value);											/*!< Multiply the phase matrix by \e val */
		skRTPhaseMatrix&	operator *=     ( float  value);											/*!< Multiply the phase matrix by \e val */
		skRTPhaseMatrix&	operator =		( const skRTPhaseMatrix& other);							/*!< Assign (copy) the other phase matrix to this one */
		skRTStokesVector	operator *		( const skRTStokesVector& stokes ) const;   				/*!< Multiply this phase matrix by a Stokes Vector to get a Stokes vector */
		skRTPhaseMatrix		operator *		( double value)                    const;					/*!< Multiply this phase matrix by a constant and return the answer */
		skRTPhaseMatrix		operator *		( float  value)                    const;
		skRTPhaseMatrix		operator -		( const skRTPhaseMatrix& other)    const;
		skRTPhaseMatrix		operator +		( const skRTPhaseMatrix& other)    const;


	public:
		bool					ApplyStokesRotation( double mu, double muprime, double dphi, skRTPhaseMatrix* rotatedmatrix );
		inline static  double	GetScatteringAngle ( double mu, double muprime, double dphi ) { return (mu*muprime + sqrt( (1-mu*mu)*(1-muprime*muprime))*cos(dphi));}
};


/*-----------------------------------------------------------------------------
 *					class skOpticalProperties					 */
/**	\ingroup skopticalprop
 *	\par Overview
 *	The skOpticalProperties class is a base class that provides an interface for calculating the
 *	scattering, absorption and extinction properties of individual atoms, molecules or particles
 *	in the atmosphere. The class is one of the foundations of the radiative transfer engines. A user
 *	will typically combine the cross-sections of individual atoms, molecules or particles provided by this set
 *	of classes with appropriate number densities to calculate extinction, scattering and extinction per cm.
 *	Users will first set the location of the calculation in the atmosphere using method #SetAtmosphericState 
 *	or #InternalClimatology_UpdateCache and then calculate cross-sections and phase matrices using #CalculateCrossSections
 *	and #CalculatePhaseMatrix. A helper function, #GetRotatedPhaseMatrix, is provided to rotate phase
 *	matrices from one reference frame to another
 *
 *  \par Atmospheric State
 *	Most atoms, molecules or particles have absorption, scattering or extinction cross-sections that depend
 *	upon atmospheric condition and this detail must be modelled by this class. For example, many cross-sections
 *	depend upon temperature, pressure or both. The skOpticalProperties class uses a climatology class, see #SetAtmosphericState,
 *	derived from class skClimatology to provide information about atmospheric state. Most users will normally
 *	call #SetAtmosphericState once to initialize the skOpticalProperties class and then use #InternalClimatology_UpdateCache to
 *	"move" aroudn teh atmosphere. The difference is that SetAtmosphericState will refresh the climatologies internal cache, which
 *	can be computationally expensive, while InternalClimatology_UpdateCache simply interpolates within the cached model, which is normally much quicker.
 *	
 *	\par Multi-Threading in Wavelength Only
 *	The optical property code was upgraded in early 2014 to be thread safe for multiple wavelength calculations at one geophysical location. This means that 
 *	methods #CalculateCrossSections and #CalculatePhaseMatrix, which calculate properties at the current location, are thread safe but methods #SetAtmosphericState 
 *	and #InternalClimatology_UpdateCache are not thread safe. Users will typically arranage their code so they multi-thread over wavelength calculations at one location and then
 *	update the location within a single thread.
 *
 *	\par Delta Function Forward scatter
 *	An optimization, SetDeltaFunctionForwardScatterFraction(), is provided for derived aerosol classes that have strong forward scatter which 
 *	can be approximated by a delta function. This forward scattering fraction is approximated as "unscattered" light that is unmodified by the species presence.
 *
 *	\par Convolution of Optical Properties with Instrument Point Spread Functions
 *	We have developed a set of interfaces that convolve optical properties with (instrument) point spread functions. This functionality
 *	relies on two interface classes, skWavelengthToPSF_Table which describes the instrument that measured the cross-section info
 *	and skOpticalProperty_AdditionalStateInfo which generates unique storage keys from the atmospheric state for caching convolved cross-sections.
 *	These two interfaces, or derived variations, are then inherited by the cross-section classes, for example see class skOpticalProperties_O3_SciaBogumilV4.
 *	The actual convolution is performed by class #skOpticalProperties_ConvolvedDiscreteWavelenCachedState and its derived classes which takes all of the above
 *	elements and behaves as a replacement for the source cross-section object except it is now convolved to the desired resolution.
 */
/*-------------------------------------------------------------------------*/

/* 
NOTE THAT skOpticalProperties crosses DLL boundaries so none of the base class methods should  use standard template structures like std::vector
as these do not behave well when mixing debug and release versions of code.
*/

class skOpticalProperties : public nxUnknown
{
	private:
		skOpticalProperties&				operator =							( const skOpticalProperties& other );  // =delete; Visual Studio 2012 does not like this yet			// Dont allow assignment
											skOpticalProperties					( const skOpticalProperties& other );  // =delete;			// Dont allow copy constructor

	public:
											skOpticalProperties					();
		virtual							   ~skOpticalProperties					(){}
		bool								DeepCopy							(const skOpticalProperties& other );
		bool								GetRotatedPhaseMatrix				( double wavenum, double mu, double muprime, double dphi, skRTPhaseMatrix* rotatedmatrix);						//!< Returns the rotated phase matrix necessary for scattering
		void								CheckCosineRange					( double *mu);
		bool								IsDeltaFunctionForwardScatter		() const	{ return DeltaFunctionForwardScatterFraction() > 1e-6;}													//!< Returns true if this cross-section has a delta function forward scatter element.

	public:
		virtual bool						SetAtmosphericState					( skClimatology* neutralatmosphere)                        = 0;			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						SetLocation							( const GEODETIC_INSTANT& pt, bool* crosssectionschanged ) = 0;			//!< Sets the atmospheric state and location for calculating cross-sections, usually temperature, pressure and position
		virtual bool						InternalClimatology_UpdateCache		( const GEODETIC_INSTANT& pt) = 0;
		virtual bool						CalculateCrossSections				( double wavenumber, double* absxs, double* extxs, double* scattxs )      = 0;								//!< Calculate cross-sections at the specified wave-number.
		virtual bool						CalculateCrossSectionsArray			( const double * wavenumber, int numwavenumber, double *absxs,  double* extxs, double* scattxs);			//!< Calculate cross-sections at the specified array of wave-numbers
		virtual bool						CalculatePhaseMatrix				( double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix);//!< Calculate the phase matrix at the specified wavenumber and scattering angle
		virtual bool						CalculateP11                        ( double wavenumber, std::pair<double, size_t> cosscatterandindex, double& p11 );
		virtual bool						IsScatterer							() const                                                         = 0;			//!< Returns true if this particles scatters radiation
		virtual bool						IsAbsorber							() const                                                         = 0;			//!< Returns true if this particles absorbs radiation radiation
		virtual bool						IsHeightDependent					() const    { return true; }
		virtual bool						IsInelasticScatterer				() const	{ return false;}
		virtual double						DeltaFunctionForwardScatterFraction () const	{ return 0.0;}

		// Functions involving Legendre polynomials
		virtual bool						LegendreCoefficientsP11(double wavenumber, double* coeff, int usermaxcoeff, int& opticalmaxcoeff);
        virtual bool                        LegendreCoefficientsPolarized(double wavenumber, double* a1, double* a2, double* a3, double* a4, double* b1, double* b2, int usermaxcoeff, int& opticalmaxcoeff);

		// Internal speedup
		virtual bool						PhaseGridHint(const std::vector<double>& cosscatterangles) { return true; }
};


/*-----------------------------------------------------------------------------
	*					SetAtmosphericState								2011-8-9*/
/** 
	*	\fn  virtual bool skOpticalProperties::SetAtmosphericState	( skClimatology* neutralatmosphere, const GEODETIC_INSTANT& pt, bool* crosssectionschanged ) = 0;
	*	Method that allows the optical properties class to update its cross-sections in response to changing any standard atmospheric parameters that affect cross-section,
	*	the external climatology should support pressure and temperature. The method is called so the optical properties class can respond to changes in the 
	*	geophysical location.The incoming climatology is external to the optical properties class and is provided by the user. The user (eg the radiative transfer model) will
	*	usually call SetAtmosphericState before calling CalculateCrossSections. The method is intended to be called many times and does not directly call sKClimatology->UpdateCache. The external
	*	user is responsible for calling skClimatology->UpdateCache. The incoming climatology is only meant to handle standard atmospheric state parameters such as pressure temperature and density
	*	as this allows it to be easily derived from a suite of climatology classes such as ECMWF, MISIS etc. 
	*
	*	Some optical property objects, eg ice and aerosols, need additional climatologies about non-standard atmospheric properties. For example, mie aerosols need climatologies
	*	of mode radius and mode width. These specialized climatologies are stored internally within the optical properties class and are automatically queried by the specific 
	*	implementation of SetAtmosphericState
	*
	*	\param neutralatmosphere
	*		A climatology passed in by the user. The optical properties class can use this class if they wish to lookup
	*		parameters of interest. The climatology will usually support, pressure, temperature and density.
	*
	*	\param pt
	*		The time and location at which the next set of cross-sections will be needed. Note that the height field in point
	*		will be valid (as this is often used to get the pressure and temperature).
	*
	*	\param crosssectionschanged
	*	returns true if the new geophysical location will cause a change in cross-section. Several species (eg Rayleigh scattering)
	*	are completely unaffected by changes in location.
	*
	**/
/*---------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
	*					skOpticalProperties::InternalClimatology_UpdateCache						2011-8-9*/
/** \fn virtual bool skOpticalProperties::InternalClimatology_UpdateCache ( const GEODETIC_INSTANT& pt) = 0
	*	Method that allows the optical properties class to update any internal
	*	climatologies it may be using. This is useful for optical properties that use
	*	internal climatologies to store various parameters. For example aerosols and ice particles store internal
	*	climatologies of the geographical/height distributions of particle size distribution paramaters, e.g. mode
	*	radius and mode width. The call to InternalClimatology_UpdateCache only reloads the caches of internal climatologies 
	*	and should be followed by subsequent calls to SetAtmosphericState.
	*
	*	\param pt
	*		The time and location for which the internal climatologies (if any) will be cached . Note that the height
	*		field is usually ignored as all internal climatologies are expected to cache height profiles of their respective parameters.
	*
	*	\return
	*		True if success otherwise false.
	**/
/*---------------------------------------------------------------------------*/


/*-----------------------------------------------------------------------------
	*					skOpticalProperties::CalculateCrossSections		2011-8-9*/
/** \fn virtual bool skOpticalProperties::CalculateCrossSections( double wavenumber )= 0;
	*	The main function of the class. The user calls this function so the optical properties object can
	*	calculate and update its three cross-sections: absorption, scattering and extinction.
	*	The user may fetch the calculated cross-sections  by calling #AbsorptionCrossSection, #ScatteringCrossSection
	*  and #ExtinctionCrossSection. This function should be designed to be thread safe.
	*
	*	There is no requirement that the class generate physically meaningful cross-sections
	*	where the sxtinction is the sum of the absorption and the scattering (although this will normally be the case)
	*
	*	Well-written classes will check that the user is not calling for a calculation under identical conditions.
	*
	*	\param wavenumber
	*		The wavenumber (cm-1) at which the calculations are required.
	*
	*	\param absxs
	*	Returns the absorption cross-section in cm-2
	*
	*	\param extxs
	*	Returns the extinction cross-section in cm-2
	*
	*	\param scattxs
	*	Returns the scattering cross-section in cm-2
	*
	*	\return
	*		True if success otherwise false.
	**/
/*---------------------------------------------------------------------------*/

class skOpticalProperties_Inelastic : public skOpticalProperties
{
	public:
						skOpticalProperties_Inelastic() {}
		virtual		   ~skOpticalProperties_Inelastic() {}

	public:
		virtual size_t	NumInelasticLines() const = 0;

		virtual bool	CalculateInelasticCrossSections			( double wavenumin, double* inelxs ) = 0;
		virtual bool	CalculateInelasticCrossSections			( double wavenumin, size_t lineidx, double* wavenumout, double* inelxs) = 0;

		virtual bool	CalculateInelasticCrossSections_Reverse	( double wavenumout, double* inelxs ) = 0;
		virtual bool	CalculateInelasticCrossSections_Reverse	( double wavenumout, size_t lineidx, double* wavenumin, double* inelxs) = 0;

		virtual bool	CalculateInelasticPhaseMatrix			( double wavenumber, double cosscatterangle, skRTPhaseMatrix* phasematrix ) = 0;

		virtual bool	IsInelasticScatterer() const override { return true; }


};
