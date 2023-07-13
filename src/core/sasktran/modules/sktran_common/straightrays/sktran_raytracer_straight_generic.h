


/*-----------------------------------------------------------------------------
 *					SKTRAN_GeometryObject		 2016- 7- 12*/
/** @ingroup grids
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_GeometryObject
{
	private:
		virtual bool EqualTo( const SKTRAN_GeometryObject& rhs ) const = 0;

	public:
		virtual ~SKTRAN_GeometryObject() { };
		
		// Many intersection calculations are bottlenecked by the sqrt call in position.Magnitude()
		// So we offer two variants of find intersections
		virtual std::array<double, 2> FindIntersections( const nxVector& look, const nxVector& position ) const { return FindIntersections( look, position, position.Magnitude() ); };
		virtual std::array<double, 2> FindIntersections( const nxVector& look, const nxVector& position, double positionmagnitude ) const = 0;

		friend 
		bool operator==( SKTRAN_GeometryObject& lhs, SKTRAN_GeometryObject& rhs ) 
		{
			return lhs.EqualTo( rhs );
		}
};



/*-----------------------------------------------------------------------------
 *					SKTRAN_GeometryObject_Sphere		 2016- 7- 12*/
/** @ingroup grids 
*/
/*---------------------------------------------------------------------------*/

class SKTRAN_GeometryObject_Sphere : public SKTRAN_GeometryObject
{
	private:
		double				m_radius;

	private:
		bool				EqualTo( const SKTRAN_GeometryObject& rhs ) const override;
	public:
		SKTRAN_GeometryObject_Sphere( double radius ) { m_radius = radius; }
		SKTRAN_GeometryObject_Sphere() { m_radius = 0.0; }
		double&	RadiusVar() { return m_radius; }
		double Radius()	const	{ return m_radius; }
		virtual std::array<double, 2> FindIntersections( const nxVector& look, const nxVector& position, double positionmagnitude ) const override;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_GeometryObject_Plane		 2016- 7- 12*/
/**  @ingroup grids
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_GeometryObject_Plane : public SKTRAN_GeometryObject
{
	private:
		nxVector			m_normal;			//!< normal vector to the plane
	private:
		bool				EqualTo( const SKTRAN_GeometryObject& rhs ) const override;
	public:
		SKTRAN_GeometryObject_Plane( const nxVector& normal ) { m_normal = normal; }

		virtual std::array<double,2> FindIntersections( const nxVector& look, const nxVector& position, double positionmagnitude ) const override;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_GeometryObject_Cylinder		 2016- 7- 12*/
/**  @ingroup grids
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_GeometryObject_Cylinder : public SKTRAN_GeometryObject
{
	private:
		nxVector			m_z;
		double				m_radius;
	private:
		bool				EqualTo( const SKTRAN_GeometryObject& rhs ) const override;
	public:
		SKTRAN_GeometryObject_Cylinder( const nxVector& z, double radius ) { m_z = z; m_radius = radius; }
		virtual std::array<double,2> FindIntersections( const nxVector& look, const nxVector& position, double positionmagnitude ) const override;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_GeometryObject_Cone		 2016- 7- 12*/
/** @ingroup grids
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_GeometryObject_Cone : public SKTRAN_GeometryObject
{
	private:
		nxVector			m_dir;					//!< unit vector of cone
		double				m_halfapexangle;		//!< in degrees
	private:
		bool				EqualTo( const SKTRAN_GeometryObject& rhs ) const override;
	public:
		SKTRAN_GeometryObject_Cone( const nxVector& dir, double halfapexangle ) { m_dir = dir; m_halfapexangle = halfapexangle; }
		virtual std::array<double,2> FindIntersections( const nxVector& look, const nxVector& position, double positionmagnitude ) const override;
};

/*-----------------------------------------------------------------------------
 *					class SKTRAN_RayTracer_Straight_Generic		 2014- 11- 6*/
/** @ingroup rays
 * This is the ray tracing class  predominantly used by Dan's HR engineclass pred
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayTracer_Straight_Generic : public SKTRAN_RayTracer_Shells
{
	private:
		nxThreadStorageMap< std::vector<double> >				 m_threadtrajectory;
		std::vector< std::unique_ptr<SKTRAN_GeometryObject>  >   m_container;
		SKTRAN_GeometryObject_Sphere							 m_earth;
		SKTRAN_GeometryObject_Sphere							 m_upperatmo;

	private:
		bool										PushBack_Distance					( double s, SKTRAN_RayStorage_Straight* storage) const;

	public:
													SKTRAN_RayTracer_Straight_Generic	(std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords);
		virtual									   ~SKTRAN_RayTracer_Straight_Generic	() { };
		bool										AddGeometryObject					( std::unique_ptr<SKTRAN_GeometryObject> obj );
		bool										SetEarthRadius						( double radius );
		void										SetUpperAtmoRadius					( double radius ) { m_upperatmo.RadiusVar() = radius; }

		virtual bool								TraceStraightRay					( SKTRAN_RayOptical_Straight*	aray) const override;

};