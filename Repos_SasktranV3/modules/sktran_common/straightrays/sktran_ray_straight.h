

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayGeometry_Straight		2014-1-25*/
/** @ingroup rays
 *	This class implements the geometry for a straight line ray that intersects
 *	a user defined grid of shells.
**/
/*---------------------------------------------------------------------------*/
//
//class SKTRAN_RayGeometry_Straight : public SKTRAN_RayGeometry_Base
//{
////	public:	// These are the methods used to access the m_quaddistances.
////		bool												ExtendReserve_Points				( size_t numpoints);
//
//	public:
//															SKTRAN_RayGeometry_Straight			( SKTRAN_RayStorage_Straight*	trajectorystorage);
//		virtual											   ~SKTRAN_RayGeometry_Straight			();
////		bool												LocationAlongRayAsVector			( const double& distancealongray, HELIODETIC_VECTOR* pt ) const;
////		bool												GetQuadratureInterpParams_startPoint( size_t quadraturepoint, size_t cellidx, double* r0, double* t0, double* rt, HELIODETIC_POINT* startpoint ) const;
////		bool												GetQuadratureInterpParams			( size_t cellidx, double* r0, double* r1, double* t0, double* t1, double* rt,
////																								  HELIODETIC_POINT* startpoint, HELIODETIC_POINT* endpoint ) const;
//
////	public: //--- virtuals inherited and overloaded by this class
////		virtual void										NotifyDerived_RayInvalid			()            override;
//
//};

class SKTRAN_RayTracer_Shells;

// Virtual methods required to use the optpropintegrator classes that do straight or piecewise straight approximate quadrature
class SKTRAN_RayOptical_StraightQuadrature_Base : public SKTRAN_RayOptical_Base
{
	public:
		virtual bool												GetQuadratureInterpParams_startPoint(size_t quadraturepoint, size_t cellidx, double* r0, double* t0, double* rt, HELIODETIC_POINT* startpoint) const = 0;
		virtual bool												GetQuadratureInterpParams(size_t cellidx, double* r0, double* r1, double* t0, double* t1, double* rt,
			HELIODETIC_POINT* startpoint, HELIODETIC_POINT* endpoint) const = 0;
};

/*-----------------------------------------------------------------------------
 *					class SKTRAN_RayOptical_Straight				2014-1-24*/
/** @ingroup rays
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayOptical_Straight : public SKTRAN_RayOptical_StraightQuadrature_Base
{
	private:
		std::unique_ptr<SKTRAN_RayStorage_Straight>			m_trajectorystorageobject;		// Note that m_trajectorystorage must apear BEFORE m_geometrystraight in the class definition.
		SKTRAN_RayStorage_Straight*							m_trajectorystorage;		// Note that m_trajectorystorage must apear BEFORE m_geometrystraight in the class definition.
		std::shared_ptr<SKTRAN_RayTracer_Shells>			m_raytracer;


	public:
															SKTRAN_RayOptical_Straight			( std::unique_ptr<SKTRAN_RayStorage_Straight> trajectorystorage, std::shared_ptr< SKTRAN_RayTracer_Shells> raytracer);
		virtual											   ~SKTRAN_RayOptical_Straight			();
//		SKTRAN_RayGeometry_Straight*						GRAY								()		 { return &m_geometrystraight;}
//		const SKTRAN_RayGeometry_Straight*					GRAY								() const { return &m_geometrystraight;}
		SKTRAN_RayStorage_Straight*							StraightStorageVar					( )		 { return m_trajectorystorage;}			//!< Straight rays need a pointer to staright ray storage but the base class holds the unique_ptr

		bool												LocationAlongRayAsVector			( const double& distancealongray, HELIODETIC_VECTOR* pt ) const;
		bool												GetQuadratureInterpParams_startPoint( size_t quadraturepoint, size_t cellidx, double* r0, double* t0, double* rt, HELIODETIC_POINT* startpoint ) const;
		bool												GetQuadratureInterpParams			( size_t cellidx, double* r0, double* r1, double* t0, double* t1, double* rt,
																								  HELIODETIC_POINT* startpoint, HELIODETIC_POINT* endpoint ) const;

		virtual bool										TraceRay_NewMethod					() override;
		virtual void										NotifyDerived_RayInvalid			() override;


};

