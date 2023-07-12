

/*-----------------------------------------------------------------------------
 *					class SKTRAN_UnitSphere_V2		2007-12-14*/
/** @ingroup unitsphere
 *	The base class defining the unit sphere used at each point in the
 *	diffuse profile table to define the outbound directions at which the
 *	(scattered) outbound radiances will be calculated.  This class expects 
 *	derived class to define the vertices.
 *
 *	\par Triangulation
 *	This class does do interpolation of the unit sphere using three point
 *	linear interpolation, method #Triangulation.  At the current time the code
 *	cannot handle arbitrary arrays of vertices as it simply assumes that the 3
 *	closest points are the best points for interpolation. This is true for the
 *	buckyball derived class  but may not be true for other geometries.  A more
 *	general solution for the interpolation would require some sort of Delaunay
 *	interpolation.  I also think the linear interpolation used at the moment
 *	is not a proper triangular interpolation and may need future work.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_UnitSphere_V2 : public nxUnknown
{

	private:
		nxVector*									m_unitvectors;			//!< The array of unit vector vertices on a unit sphere;
		size_t										m_numunitvectors;		//!< The number of unit vectos 
		nx1dArray<double>							m_weights;				//!< the array of cubature weights for the array of unit vertices


	private:
		void										ReleaseVertices				();

	public: 
		virtual bool								InterpolateTriangle			( const nxVector& P, size_t* unit_indexptr, double* unit_weightptr ) const; // nx: This is really private but I make it public for testing purposes. srd: This can give negative weights, should be replaced with the Delaunay version 

	protected:
		bool										FindThreeClosestVertices	( const nxVector& unit, size_t* idxmax1, size_t* idxmax2, size_t* idxmax3, size_t* idxmax4  ) const;
		bool										AllocateVertices			( size_t numunitvectors);
		nxVector&									UnitVectorAtVar				( size_t idx );
		double&										CubatureWeightAtVar			( size_t idx );

	public:
		const nxVector&								UnitVectorAt				( size_t idx )	const;
		size_t										NumUnitVectors				()				const { return m_numunitvectors;}
		double										CubatureWeightAt			( size_t idx )	const;
		bool										LocalLookToAziZen			( const nxVector& locallook, double* azi, double*zen ) const;

	public:
													SKTRAN_UnitSphere_V2		(); 
		virtual									   ~SKTRAN_UnitSphere_V2		();
		virtual size_t								MaxNumInterpIndex () const  { return 3; }
		virtual bool								Triangulate					( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const;
		virtual bool								Triangulate					( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices, size_t& speedHelper) const;

};



/*-----------------------------------------------------------------------------
 *					class SKTRAN_UnitSphere_WithLookupTable_V2		2007-12-14*/
/** @ingroup unitsphere
 *	The base class defining the unit sphere used at each point in the
 *	diffuse profile table to define the outbound directions at which the
 *	(scattered) outbound radiances will be calculated.  This class expects 
 *	derived class to define the vertices.
 *
 *	\par Triangulation
 *	This class does do interpolation of the unit sphere using three point
 *	linear interpolation, method #Triangulation.  At the current time the code
 *	cannot handle arbitrary arrays of vertices as it simply assumes that the 3
 *	closest points are the best points for interpolation. This is true for the
 *	buckyball derived class  but may not be true for other geometries.  A more
 *	general solution for the interpolation would require some sort of Delaunay
 *	interpolation.  I also think the linear interpolation used at the moment
 *	is not a proper triangular interpolation and may need future work.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_UnitSphere_WithLookupTable_V2 : public SKTRAN_UnitSphere_V2
{
	protected:
	struct SKTRAN_UnitSphereInterpolationEntry
	{
	bool		isvalid;
	size_t		unit_indexptr[3];
	double		unit_weightptr[3];
	};

	private:
		double										m_deltazen;				//!< Lookup table zenith resolution in radians
		double										m_deltaazi;				//!< Lookup table azimuth resolution in radians
		size_t										m_numazimuth;			//!< Num Azimuths resolution for lookup table;
		size_t										m_numzenith;		    //!< Num Zeniths resolution for lookup table;
		size_t										m_numlookup;			//!< The number of lookup entries for fast lookup of indexing.
		SKTRAN_UnitSphereInterpolationEntry*		m_lookupentries;		//!< The array  of fast lookup entries;

	private:
		bool										AllocateAziZenTable			();
		bool										GenerateLookupEntry			(  SKTRAN_UnitSphere_WithLookupTable_V2::SKTRAN_UnitSphereInterpolationEntry* entry,  size_t zenidx, size_t aziidx );
		void										ReleaseLookupTable			();

	protected:
		bool										InitializeLookupTable			();
		bool										AllocateLookupTable					( size_t numlookup);
		SKTRAN_UnitSphereInterpolationEntry*		InterpolationEntryAtVar				( size_t lookupindex );
		const SKTRAN_UnitSphereInterpolationEntry*	InterpolationEntryAt				( size_t lookupindex ) const;
		size_t										NumLookupEntries					() const				{ return m_numlookup;}

	public:
													SKTRAN_UnitSphere_WithLookupTable_V2	(); 
		virtual									   ~SKTRAN_UnitSphere_WithLookupTable_V2	();
		virtual bool								Triangulate							( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const;
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereBucky_V2		2007-12-17*/
/** @ingroup unitsphere
 *	A class used to define the outbound vertices of a bucky-ball. The
 *	bucky-ball consists of 92 vertices that cover the sphere relatively
 *	smoothly.
 **/
/*---------------------------------------------------------------------------*/

class  SKTRAN_UnitSphereBucky_V2: public SKTRAN_UnitSphere_WithLookupTable_V2
{
	private:
		double								m_phi;

	public:
											SKTRAN_UnitSphereBucky_V2	();
		virtual							   ~SKTRAN_UnitSphereBucky_V2	();
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_GridDefOutboundUnitSphereME100	2008-2-10*/
/** @ingroup unitsphere
 *	A class that specifies the distribution of points over a sphere. I have taken a table
 *	of points using a minimum entropy algorithm. The table came from work by 
 *  Rob Womersley, http://web.maths.unsw.edu.au/~rsw/Sphere/. He has produced distributions
 *	of points over spheres.  In this instance I have the 100 point table.  Others could be coded
 *	into the system.
 **/
/*---------------------------------------------------------------------------*/

class  SKTRAN_UnitSphereME: public SKTRAN_UnitSphere_WithLookupTable_V2
{

	public:
											SKTRAN_UnitSphereME	(size_t numvertexrequested );
		virtual							   ~SKTRAN_UnitSphereME	();
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphereLatLonGrid		2008-3-12*/
/** @ingroup unitsphere
 *	Creates vertices distributed on a Lat/Lon grid. The class takes
 *	a distribution of zenith angles and a distribution of azimuth angles and
 *	places the points over the unit sphere. The class forces a single unit vector
 *	to be placed at the upward and downward poles even if the user has not
 *	requested them. 
 *
 *	This class covers the entire sphere in zenith angle and azimuth (both in radians) using a triangle fan
 *	around the poles and triangle strips at all other latitudes between. The vertices of the 
 *	triangle strips (and fans) map to the unit vectors placed on the sphere. However more than one
 *	triangle vertex can refer to the same unit vector. This is used to provide proper interpolation over the
 *	0 and 2 pi boundaries
 * 
 *	The triangular strips are built from the one azimuth distribution that is offset as we move from one row to
 *	next, these are referred to as the odd and even rows. The idea is that the staggered "odd" and "even"
 *	distributions will help with interpolation. The odd azimuths are offset to the middle of the even distributions.
 *
 *	I have tried to provide a picture that shows the triangle fans and triangle strips below. Its tricky to do with
 *	just text characters but hopefully you get the gist.
 *
 *
 *            0 rads                            2 pi rads
 *            |   Vertex grid between 0 an 2 pi   |  Wrap around zone
 *            |                                   |
 *       o----|-o----o-----o-----o-----o-----o----|-o--  TOP ROW IS POLE, EVEN AZIMUTHS, BUT ONLY ONE UNIT VECTOR 
 *        \   |/\    /\    /\    /\    /\    /\   |/\   
 *         \  /	 \  /  \  /  \  /  \  /  \  /  \  /  \  
 *          \/|	  \/    \/    \/    \/	  \/    \/|   \ 
 *       ---0-|---1-----2-----3-----4-----5-----6-|---7	 ODD ROW AZIMUTHS, 8 AZIMUTHS BUT ONLY 6 UNIT VECTOR BETWEEN 0 and 2 pi
 *         / \|   /\    /\    /\    /\    /\    /\|   /
 *        /   \  /  \  /  \  /  \  /  \  /  \  /  \  / 
 *       /    |\/    \/    \/    \/   \/    \/    |\/   
 *       0----|-1----2-----3-----4-----5-----6----|-7--  EVEN ROW AZIMUTHS, 8 AZIMUTHS BUT ONLY 6 UNIT VECTORS BETWEEN 0 and 2 pi
 *        \   |/\    /\    /\    /\    /\    /\   |/\  
 *         \  /	 \  /  \  /  \  /  \  /  \  /  \  /  \ 
 *          \/|	  \/    \/    \/    \/	  \/    \/|   \
 *       ---0-|---1-----2-----3-----4-----5-----6-|---7	 ODD ROW AZIMUTHS
 *         / \|   /\    /\    /\    /\    /\    /\|   /
 *        /   \  /  \  /  \  /  \  /  \  /  \  /  \  / 
 *       /    |\/    \/    \/    \/   \/    \/    |\/   
 *       0----|-1----2-----3-----4-----5-----6----|-7--	 EVEN ROW AZIMUTHS
 *        \   |/\    /\    /\    /\    /\    /\   |/\   
 *         \  /	 \  /  \  /  \  /  \  /  \  /  \  /  \  
 *          \/|	  \/    \/    \/    \/	  \/    \/|   \ 
 *       ---o-|---o-----o-----o-----o-----o-----o-|---o  BOTTOM ROW IS POLE, EVEN OR ODD AZIMUTHS, BUT ONLY ONE VERTEX (UNLESS ITS GROUND)
 *
 **/
/*---------------------------------------------------------------------------*/

class  SKTRAN_UnitSphereLatLonGrid: public SKTRAN_UnitSphere_V2
{
	private:
		SKTRAN_GridDefDiffuseIncomingZenith_V21*		m_zenith;					// The array of zenith angles in rdaiancs, stored in ascending  order. Must include upper pole and ground/horizon or bottom pole.
		SKTRAN_GridDefDiffuseIncomingAzimuth_V21*	m_azieven;					// Holds the grid of azimuth angles same as the users angles (in radians), add on extra points at beginning and end to manage wraparound.
		SKTRAN_GridDefDiffuseIncomingAzimuth_V21*	m_aziodd;					// Holds the grid of azimuth angles halfway between the even angles, adds on extra points at beginning and end to manage wraparound

	private:
		struct trianglevertexstruct
		{
			size_t										zen;
			size_t										azi;
			SKTRAN_GridDefDiffuseIncomingAzimuth_V21*	azigrid;
		};

	private:
		bool								ConfigureSolidAngles		();
		bool								AllocateInternalZenithGrid	(const SKTRAN_GridDefDiffuseIncomingZenith_V21* userzenith);
		bool								AllocateInternalAzimuthGrid (const SKTRAN_GridDefDiffuseIncomingAzimuth_V21* userazi);
		double								AziInRange0To360			( double val );
		bool								IsPole						( double zenithangle_radians ) const;
		void								ReleaseGrids				();
		bool								CreateVertexGrid			();
		bool								AssignCubatureWeights		();
		size_t								ZenAziIndexToVertexIndex	( size_t zenidx, size_t aziidx ) const;
		size_t								ZenithIndexToVertexIndex	( size_t zenidx ) const;
		size_t								AziIndexToVertexOffset		( size_t aziidx ) const;
		bool								IsGroundPoint				() const;
		bool								IsInsideTriangle			( double zen, double azi, const trianglevertexstruct* zenazivertices, trianglevertexstruct* thetriangle, double* weights) const;
		bool								FindInsideTriangle			( double zen, double azi, const trianglevertexstruct* triangleptr[3], trianglevertexstruct* thegoodtriangle, double* weights) const;
		bool								InterpEvenUpperTriangleStrip( double zen, double azi, size_t zenidx, size_t aziidx, trianglevertexstruct* thegoodtriangle, double* weights) const;
		bool								InterpOddUpperTriangleStrip ( double zen, double azi, size_t zenidx, size_t aziidx, trianglevertexstruct* thegoodtriangle, double* weights) const;
		bool								InterpTriangleFan			( double zen, double azi, size_t zenidx1, SKTRAN_GridDefDiffuseIncomingAzimuth_V21* azigrid, trianglevertexstruct* thegoodtriangle, double* weights) const;

	public:
											SKTRAN_UnitSphereLatLonGrid	();
		virtual							   ~SKTRAN_UnitSphereLatLonGrid	();
		bool								DefineGrid					( const SKTRAN_GridDefDiffuseIncomingZenith_V21* userzenith, const SKTRAN_GridDefDiffuseIncomingAzimuth_V21* userazimuth );
		virtual bool						Triangulate					( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const override;
		
	public:
		size_t								NumZenith					() const;	
		size_t								NumAzimuthVertex			() const { return (m_azieven->NumAngles() -2);}
		bool								GetZenithVertexIndexAndNumVertexAzimuth	( size_t zenidx, size_t* firstvertexidx, size_t* numazi ) const;

};

/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_Delaunay		2013-09-18*/
/** @ingroup unitsphere
 * Constructs the Delaunay triangulation of the vector of unit vectors
 * passed in at ConstructTriangulation. Optionally, the unit
 * vector openAxis can be passed in -- the triangulation will
 * be calculated to include this point, but it will be given
 * zero weight in any lookups. These vectors should be well-
 * seperated enough that the cosine between any two of them is greater
 * than ~10e-5. Also, it must be possible to construct a tetrahedron
 * from unitVecs (and openAxis) that contains the origin. For the 
 * triangulation to look "good" for a set of unit vectors clustered
 * mostly on one part of the sphere openAxis should be chosen to 
 * be opposite this cluster on the sphere. E.g., if unit vectors
 * are clustered around the tangent point, openAxis should be the
 * reflection of the tangent point across the origin. If points
 * are clustered and openAxis is not used the triangulation will
 * not look good, and very likely the triangulation method will fail. 
 *
 * The Delaunay triangulation is unique, but it may not always 
 * intuitive which edges will be joined by the triangulation. It
 * would be good for users to check what triangulations look like
 * for the distribution of point sets they're using either by
 * using the class' PrintTriangulation function (which prints a row
 * of indices for each triangle in the triangulation, indexing vertices
 * in the order they were passed in to CreateTriangulation, openAxis
 * having index {numUnitVecsPassedIn}+1) or by calculating the triangulation
 * in some other program (e.g. using MATLAB method convhull).
 **/
/*---------------------------------------------------------------------------*/
class SKTRAN_UnitSphere_Delaunay : public SKTRAN_UnitSphere_V2
{

	protected:
        template<typename T> struct tuple2 { T d[2];};
	    template<typename T> struct tuple3 { T d[3];};
	    template<typename T> struct tuple4 { T d[4];};

        bool            m_hasDummyPoint;
        nxVector        m_openAxis;
        tuple3<size_t>* m_faces;
        tuple3<size_t>* m_neigs;
		std::vector<nxVector> m_faceNormals;
        size_t          m_numFaces;

	private:
		void										ReleaseResources				();
    protected:
        bool                                        CopyVerticesToInternal          ( const nxVector* unitVecs_xyz, size_t numTriplets, const nxVector* openAxis = NULL );
        bool                                        ConstructSimplex                ( tuple4<size_t>& simplexVerts );
        bool                                        ConstructSimplex_bruteForce     ( tuple4<size_t>& simplexVerts );
        virtual bool                                ConstructLookupObjects          ( ) = 0;
        bool                                        ConstructTriangulation          ( );
        bool                                        CreateDummyPoint                ( );
        bool                                        GetHyperplane                   ( tuple3<size_t>& verts, tuple4<double>& target) const;
        bool                                        TestPointUnderPlane             ( tuple4<double> n4, size_t vertIndex, double& distanceUnder) const;
		virtual bool                                FindThreeNearestVertices        ( const nxVector& unit, size_t* unit_indexptr, size_t maxvertices, size_t& faceIndex) const = 0;
		virtual bool                                FindInterpolationWeights        ( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t faceidx) const = 0;

	public:
		virtual bool								InterpolateTriangle				( const nxVector& P, size_t* unit_indexptr, double* unit_weightptr, size_t faceidx ) const;	// This should be protected, but it's public for debugging in the base class

	public:
													SKTRAN_UnitSphere_Delaunay		(); 
		virtual									   ~SKTRAN_UnitSphere_Delaunay		();
        bool                                        CreateTriangulation             ( const nxVector* unitVecs, size_t numUnitVecs, const nxVector* openAxis = NULL );		
		virtual bool								Triangulate						( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const;
		virtual bool								Triangulate					    ( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices, size_t& speedHelper) const;
		void                                        PrintTriangulation              ( std::string fileToPrintTo ) const;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_Delaunay_nonTabledLookup		2013-09-18*/
/** @ingroup unitsphere
 * Performs triangulation of query points without using a lookup table. May be slower
 * than classes using a lookup table for cases where true random access is needed
 * (e.g. triangulating unsorted, randomly distributed query points would be quite
 * slow), but is very fast if points ordered as described below. 
 * 
 * The class is fast if query point (n) is usually bound by either the same 
 * triangulation element as was point (n-1), or by one of that element's neighbours.
 * It doesn't matter if the set of query points spans the whole sphere, so long as, 
 * as we follow the set of points, the path we take moves from one triangulation 
 * element to its neighbour, and then to that elements neighbour, etc., jumping 
 * more than one neighbour as infrequently as possible. It's not terribly expensive
 * for the path to loop back on itself (e.g. say from the north pole, slowly down to 
 * (lat,lon)=(0,90), back to the north pole, then to (0,-90)) -- it's much more
 * important that sequential points be in the same neighbourhood. Examples of good
 * paths to follow:
 *  - Spiralling down the sphere from around a point
 *  - Following a processing great circle (e.g. OSIRIS orbit)
 *  - Raster pattern
 *
 * If query points follow such a pattern it is best to call Triangulate with argument 
 * speedHelper, which must be unique to each tread using the sphere. 
 * SKTRAN_UnitSphere_Delaunay_nonTabledLookup probably out-performs any other 
 * possible Delaunay lookup algorithm if called in this way. The value in speedHelper
 * must be preserved between function calls for there to be a performance improvement. 
 * 
 * If there is no pattern to the query points but they're usually in the same neighbourhood,
 * (e.g. often sampling near a reference point, not so much far from that point), calling 
 * OptimizeForLookupInNeighbourhoodOf with the tangent point or reference point as the argument
 * and then calling Triangulate *without* speedHelper will perform better. Lookups are much 
 * slower (by a factor of ten is reasonable) if query points can't be given some simple 
 * pattern where point (n) is close to point (n-1). 
 **/
/*---------------------------------------------------------------------------*/
class SKTRAN_UnitSphere_Delaunay_nonTabledLookup : public SKTRAN_UnitSphere_Delaunay
{

    private:
		std::vector< tuple3<nxVector> >             m_edgenorms;
		size_t										m_startFace;

	private:
		void										ReleaseResources				     ( );
        virtual bool                                ConstructLookupObjects               ( );
		virtual bool                                FindThreeNearestVertices             ( const nxVector& unit, size_t* unit_indexptr, size_t maxvertices, size_t& startface) const;
		virtual bool								FindThreeNearestVertices_bruteforce	 ( const nxVector& unit, size_t* unit_indexptr, size_t maxvertices ) const;
		virtual bool                                FindThreeNearestVertices_directed    ( const nxVector& unit, size_t* unit_indexptr, size_t maxvertices, size_t& startface) const;
		virtual bool                                FindThreeNearestVertices_undirected  ( const nxVector& unit, size_t* unit_indexptr, size_t maxvertices) const;
		virtual bool                                FindInterpolationWeights             ( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t faceidx) const;

	public:
                                                    SKTRAN_UnitSphere_Delaunay_nonTabledLookup  ( );
		virtual									   ~SKTRAN_UnitSphere_Delaunay_nonTabledLookup  ( );
		bool										OptimizeForLookupInNeighbourhoodOf      (const nxVector& unit);
		
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_Delaunay_binaryLookup		2013-09-18*/
/** @ingroup unitsphere
 * Skeleton of unfinished class. 
 * May be useful if there is a lot of random access needed on the sphere. The
 * idea was:
 *  - Imagine a plane parallel to the xy-plane that starts at z=-1, moves up
 *    to z=1
 *  - A circle of varying radius is defined by the intersection of the sphere
 *    with this plane. Different intervals of the circle are bound by different
 *    elements of the triangulation. 
 *  - At every z where the sequence of these intervals changes (e.g. the plane
 *    moves up beyond the region bound by a certain element, or it moves up
 *    and now part of it is in an element that wasn't intersected before)
 *    z is pushed onto the zTable
 *  - For every z on the zTable, the edge that bounds each interval on 
 *    the circle (say from the counterclockwise direction) is pushed onto
 *    a vector along with the index of the triangulation element associated
 *    with that interval; this vector is sorted so that iterating through the 
 *    vector corresponds to moving counterclockwise around the circle
 * -  To find which triangulation element bounds query vector u:
 *    - Take the z-component of u; look this up in the zTable (lower_bound)
 *      to find which vector should be searched in the next step
 *    - Start at the beginning of the vector; there is a plane formed by 
 *      the origin--{low-z vector element}--{high-z vector element}. The
 *      normal to this plane can be calculated; if u.Dot(normal) is 
 *      positive then u is counterclockwise from where the edge intersects
 *      the sphere (that is, ccw in rotation around z-axis) OR there
 *      are large faces in the triangulation so that there's a hemisphere
 *      without any unitVectors in it
 *    - You can do a sort of binary search -- check the middle vector element,
 *      1/4-element, etc., looking for where u.Dot(normal) changes sign. 
 *      u is in the triangulation element corresponding to the first vector
 *      element for which u.Dot(normal) is positive. 
 * 
 * The hard part is figuring out what z need to be in the zTable -- this
 * isn't necessarily obvious. Imagine a line on the Earth joining Iqualuit
 * to Anchorage so there are points on the line with z greater than z at 
 * the endpoints -- it might be quite a bit of work to figure out what should
 * be in the zTable if there are many lines like this on the triangulation. 
 * 
 * SKTRAN_UnitSphere_Delaunay_nonTabledLookup is good if the query points can
 * be ordered as described in the comments -- it's probably easier to do that
 * ordering than it is to impliment this class. 
 **/
/*---------------------------------------------------------------------------*/
class SKTRAN_UnitSphere_Delaunay_binaryLookup : SKTRAN_UnitSphere_Delaunay
{
    private:	
    class lutType {
        public:
            nxVector normal;
            size_t fidx;
            lutType(const nxVector& n, size_t fidx) : normal(n), fidx(fidx) {}
    };
    std::vector<double>  m_zTable;
    std::vector<lutType> m_lookUpTable;

    void										GetFaceIndicesOrder				( size_t fidx, tuple3<size_t>& order) const;
    double										GetFaceSignedness				( size_t fidx, const tuple3<size_t>& order) const;
    virtual bool                                ConstructLookupObjects          ( );

	public:
        SKTRAN_UnitSphere_Delaunay_binaryLookup();
       ~SKTRAN_UnitSphere_Delaunay_binaryLookup();

};


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_V2		2014-03-24*/
/** @ingroup unitsphere
 *	Defines a plane in which points are projected on 
 *
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_UnitSphere_Plane : public SKTRAN_UnitSphere_V2
{
	private:
		nxVector						m_normal;
		nxVector						m_x;
		nxVector						m_y;
		nxVector						m_reference;
		SKTRAN_GridDefAngular_V21		m_angles;

	private:
		double							ProjectedAngle( const nxVector& unit ) const;
		nxVector						UnitVectorFromAngle( double th ) const;
		bool							CheckForUniformAngles();
	public:
		SKTRAN_UnitSphere_Plane						();
		virtual ~SKTRAN_UnitSphere_Plane					() { };
		bool							ConstructPlane	( const std::vector<nxVector>& unitvecs, size_t refindex = 0 );
		bool							ConstructPlane  ( const std::vector<double>& angles, const nxVector& reference, const nxVector& normal );
		size_t							MaxNumInterpIndex() const override { return 2; };
		virtual bool					Triangulate		( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const;
		virtual bool					Triangulate		( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices, size_t& speedHelper) const;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_UnitSphere_Dummy		2015-01-29*/
/**  @ingroup unitsphere
 *	 Dummy unitsphere class which always returns back 1, i.e. no lookup
 *   This allows classes which use a unitsphere for 3d lookup to act as 1d tables
 *   Must be constructed with a reference location (usually the reference point)
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_UnitSphere_Dummy : public SKTRAN_UnitSphere_V2
{
	public:
		SKTRAN_UnitSphere_Dummy			( const nxVector& loc );
		virtual ~SKTRAN_UnitSphere_Dummy		() { };
		size_t							MaxNumInterpIndex () const override { return 1; }
		virtual bool					Triangulate		( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const;
		virtual bool					Triangulate		( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices, size_t& speedHelper) const;
};