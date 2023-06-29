

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayFactory_Base	 2016- 7- 12*/
/** @ingroup rays  
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayFactory_Base //: public nxUnknown
{
	private:
		std::shared_ptr< const SKTRAN_CoordinateTransform_V2>		m_coords;
		std::weak_ptr < const SKTRAN_RayFactory_Base>				m_thisfactory;

	protected:
		bool														SetThisFactoryInternal	(std::shared_ptr< const SKTRAN_RayFactory_Base> thisfactory ) { NXASSERT(( thisfactory.get() == this)); m_thisfactory = thisfactory; return true;}
	public:
		const SKTRAN_CoordinateTransform_V2*						CoordsPtr		() const		{ return m_coords.get();}
		std::shared_ptr< const SKTRAN_CoordinateTransform_V2>		CoordsObject	() const		{ return m_coords;}
		const std::weak_ptr <const SKTRAN_RayFactory_Base>&			ThisFactory		() const		{ return m_thisfactory;}

	public:
												SKTRAN_RayFactory_Base  (std::shared_ptr< const SKTRAN_CoordinateTransform_V2>& coords) { m_coords = coords;}
		virtual								   ~SKTRAN_RayFactory_Base	(){}
		virtual const SKTRAN_RayTracer_Base*	RayTracerBase			() const = 0;
		virtual bool							CreateRayObject			( std::unique_ptr<SKTRAN_RayOptical_Base>* rayptr) const = 0;
		virtual bool							ConfigureOptical        ( SKTRAN_AtmosphericOpticalState_V21 *opticalstate, double wavelen_nm, GEODETIC_INSTANT referencepoint) = 0;
//		virtual bool							TraceRay				( SKTRAN_RayOptical_Base* rayptr) const = 0;
//		virtual bool							SetThisFactory			(std::shared_ptr< const SKTRAN_RayFactory_Base> thisfactory ) = 0;
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayFactory								 2014- 11- 7*/
/** @ingroup rays
 *  A templated class used to construct rays for the various radiative
 *	transfer engines. This class coordinates two attributes associated with a ray:
 *
 *	 -(i) the type of ray, 
 *	 -(ii) the type of ray tracer
 *
 *	USers will generate a templated Rayfactory class and this is used to create and trace rays. The
 *	templatallows us to support several ray tracing  algorithms (straight lines, curved etc.) and the different
 *	different storage requirements for different ray type and engines.
 *
 *	Template RAY_TYPE can be any class derived from SKTRAN_RayOptical_Base and typical values
 *	are SKTRAN_RayOptical_Straight and SKTRAN_RayOptical_Curved_Piecewise. This specifies the
 *	type of ray required.
 *
 *	Template RAYTRACER_TYPE is any class derived SKTRAN_Ray_Tracer_Base. The RAYTRACER class must provide 
 *	function called TraceRay that supports the RAY_TYPE as its first argument. Note that we do intentionally
 *	do not virtualize the method TraceRay as we want to make sure that ill-formed RayFactory templates will generate compiler
 *	errors and we also want to try and detect in real-time when straight rays are accidentally passed to curved ray tracers etc.
 *
 *	Template STORAGE_TYPE is the class that is used to store information for the ray. This replaces the MinimumContainer concept
 *	introduced in earlier versions of SasktranV3
 *
 *	
 **/
/*---------------------------------------------------------------------------*/

template <class RAY_TYPE, class RAYTRACER_TYPE, class RAYSTORAGE_TYPE>

class SKTRAN_RayFactory : public SKTRAN_RayFactory_Base
{
	private:
		std::shared_ptr< RAYTRACER_TYPE>					m_raytracer;

	public:
						SKTRAN_RayFactory (std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords)
			                             : SKTRAN_RayFactory_Base( coords),
			                               m_raytracer (new RAYTRACER_TYPE(coords) )
		{
		}

		virtual	~SKTRAN_RayFactory ()
		{
		}

		/*-----------------------------------------------------------------------------
		 *					CreateRayObject		 2014- 11- 20*/
		/** **/
		/*---------------------------------------------------------------------------*/

		virtual bool CreateRayObject( std::unique_ptr<SKTRAN_RayOptical_Base>* rayptr) const override
		{

			rayptr->reset( new RAY_TYPE( std::unique_ptr<RAYSTORAGE_TYPE>( new RAYSTORAGE_TYPE(m_raytracer->CoordsObject())), m_raytracer) );
			return (rayptr->get() != nullptr);
		}

		virtual bool ConfigureOptical(SKTRAN_AtmosphericOpticalState_V21 *opticalstate, double wavelen_nm, GEODETIC_INSTANT referencepoint)
		{
			return m_raytracer->ConfigureOptical(opticalstate, wavelen_nm, referencepoint);
		}




		/*-----------------------------------------------------------------------------
		 *					RayTracerBase		 2014- 11- 20*/
		/** **/
		/*---------------------------------------------------------------------------*/

		const SKTRAN_RayTracer_Base*	RayTracerBase () const
		{ 
			return m_raytracer.get();
		}


		/*-----------------------------------------------------------------------------
		 *					RayTracer		 2014- 11- 20*/
		/** **/
		/*---------------------------------------------------------------------------*/

		RAYTRACER_TYPE*	 RayTracer () 
		{ 
			return m_raytracer.get();
		}
};

