//#include "sktran_hr_internals.h"

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Index		2014-10-30*/
/** Each incoming ray to a diffuse point gets radiance contributions 
 *  from the outgoing rays of other diffuse points.  This class represents
 *  one of those contributions, an index corresponding to an outgoing ray
 *  of a diffuse point and its weight.  All of the outgoing radiances is stored
 *  as one vector in the DiffuseTable, and index corresponds to the outgoing 
 *  radiance's location in that vector.
 **/
/*---------------------------------------------------------------------------*/

struct SKTRAN_HR_Diffuse_Index
{
	unsigned int			index;
	SKTRAN_HR_WEIGHT_TYPE	weight;
    public:
    inline bool operator< (const SKTRAN_HR_Diffuse_Index& other) const { return index<other.index;}
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Index_Array		2014-10-30*/
/** The collection of all SKTRAN_HR_diffuse_Index's for a single incoming
 *  ray.  Summing over all of these indicies will give the incoming radiance
 *  for that ray.  All of the incoming radiances are stored as one vector
 *  in the diffuse table, and matrixidx corresponds to the incoming radiances
 *  location in that vector.
 **/
/*---------------------------------------------------------------------------*/

struct SKTRAN_HR_Diffuse_Index_Array
{
	std::vector<SKTRAN_HR_Diffuse_Index> diffindex;
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Diffuse_Point		2014-10-30*/
/** A single diffuse point within the diffuse table.  The point stores the incoming and outgoing
 *	unit spheres which are used to perform atmopsheric scattering for points in the atmosphere
 *	and surface scattering for points stored on the ground. 
 *
 *  The diffuse, higher order scattering from the diffuse points is coordinated by the global diffuse table and is typically
 *	broken into three components, following the 3 formulae,
 *
 *  - \f$ A_v(i) = k_{scatt}(\theta).d\Omega \f$
 *  - \f$ I_{out}(k) = \sum_{j}A_v(j).I_{in}(j) \f$
 *  - \f$ I_{in}(j) = \sum_{i} w(i).I_{out}[index(i)] \f$
 *
 *	-# The first equation is performed once, early in the RT multiple scatter calculation and is stored in the global \f$ A_v \f$ matrix 
 *	   stored in the diffuse table, see #SKTRAN_HR_Diffuse_Table_CPU. The calculation is performed in #Avals_ScalarStore::CalcScatteringMatrixPoint.
 *	-# The second equation is performed during the scattering phase of each higher order scatter., see #RadStore_Scalar::ScatterPoint.
 *	-# The third equation ifs performed during the line of sight intergration phase of each higher order of scatter. see #RadStore_Scalar::ComputeNextOrderPoint
 *
 *	\par Ground Scattering
 *	Ground scattering uses the same framework as atmospheric scattering but adjustments are made to ground points to replace the atmospheric scattering terms
 *	with surface scattering terms. In this implementation surface scattering uses the same 3 equations as atmospheric scattering except \f$ k_{scatt}\f$ is replaced 
 *	with the BRDF function and \f$I_{out}\f$ is replaced with downward flux. 
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_Diffuse_Point
{
	private:
		const SKTRAN_UnitSphere_V2*							m_incomingsphere;			//!< Specifies the directions of rays coming into this sphere. Note that the rays are actually "looking away" from the point.
		const SKTRAN_HR_OutgoingSphereObject_Base*          m_outgoingsphereobj;		//!< Specifies the directions of "diffuse" rays shooting out of this point used to cerate the diffuse radiation field in the diffuse table.
		std::vector<HELIODETIC_UNITVECTOR>					m_rotated_incomingunitvectors;		// !< The incoming unit vectors from the incoming sphere, rotated to the global HELIODETIC coordinate system with SUN along Z axis

		HELIODETIC_POINT									m_location;					//!< The location of this point in the Heliodetic systsm
		std::vector<double>									m_firstorderradiances;		//!< Stores the first oder radiances. Its for info only. Its not part of the calculation
		std::vector<SKTRAN_HR_Diffuse_Index_Array>			m_diffuseindexes;			//!< The table of indices and weights used to create the radiances coming into this point. Its a linear combination of outbound radiances in the diffuse table
		bool												m_isgroundpoint;			//!< True if this is a groundpoint. Groundpoints treat scattering as surface scattering using BRDF rather than scattering off the atmospheric particles/molecules
        size_t												m_pointidx;					//!< Unique identifier for this point in the DiffuseTable
		size_t												m_inidx;					//!< The start of incoming ray indices for this diffuse point in the Diffuse Table storage variables
		size_t												m_outidx;					//!< The start of outgoing ray indices for this diffuse point in the Diffuse Table storage variables
		size_t												m_scatvalidx;				//!< The start of entries for the scattering matrix of this diffuse point in the Diffuse Table storage variables

	private:
		virtual bool										CreateRotated_GlobalUnitVectors();
		bool												ReleaseResources();

	public:
															SKTRAN_HR_Diffuse_Point();
		virtual											   ~SKTRAN_HR_Diffuse_Point();
		nxVector											Rotate_Heliodetic_to_Local( const HELIODETIC_UNITVECTOR& helio) const;
		bool												ConfigureSpheres					( const SKTRAN_UnitSphere_V2* incomingsphere, const SKTRAN_HR_OutgoingSphereObject_Base*  outgoingsphereobj, const HELIODETIC_POINT& location, bool isground);
		const std::vector<SKTRAN_HR_Diffuse_Index_Array>&	DiffuseIndices						() const { return m_diffuseindexes;}
		std::vector<SKTRAN_HR_Diffuse_Index_Array>&			DiffuseIndicesVar					()		 { return m_diffuseindexes;}
		bool												CleanDiffuseIndex					();
		const SKTRAN_UnitSphere_V2*							IncomingUnitSphere					() const { return m_incomingsphere;}
		const HELIODETIC_POINT&								Location							() const { return m_location; } 
		size_t												NumIncomingRays						() const { return m_incomingsphere->NumUnitVectors(); }
		size_t												NumOutGoingRays						() const;
		size_t												NumGroundDownwardFlux				() const {  NXASSERT( (m_isgroundpoint)  ); return NumIncomingRays();}
        size_t												UniquePointIdentifier				() const { return m_pointidx;}
		bool												SetPointIndices						( size_t pointidx, size_t incomingcounter, size_t outgoingcounter, size_t scattvalcounter, size_t numdiffuseindices);
		bool												DumpFirstOrderRadiance				( std::string filename );
		void												SetFirstOrderRadiance				( size_t idx, const double& radiance ) { m_firstorderradiances[idx] = radiance; }
		const HELIODETIC_UNITVECTOR&						IncomingRayGlobalCoords				( size_t idx ) const { return m_rotated_incomingunitvectors[idx]; }
		nxVector											IncomingUnitRayLocalCoords			( size_t idx ) const { return m_incomingsphere->UnitVectorAt( idx ); }
		void												OutgoingRayLocalCoords				( size_t idx, nxVector& outray ) const;
		double												InCubatureWeight					( size_t idx ) const { return m_incomingsphere->CubatureWeightAt( idx ); }
		bool												IsGroundPoint						() const { return m_isgroundpoint;}
		double												OutgoingCubatureWeight				( size_t outidx ) const;
		bool												CorrectForIndexDuplicates			( size_t inrayidx );
		bool												CorrectForIndexDuplicates_Sorted	( size_t sortedray_inrayidx );
        double												CosScatteringAngle					( size_t outidx, size_t inidx ) const { nxVector outray; OutgoingRayLocalCoords( outidx, outray); return -(m_incomingsphere->UnitVectorAt(inidx) & outray); }

	public:
		/*virtual*/ size_t									NumUniqueScatterIncoming    () const;
		/*virtual*/ size_t									NumUniqueScatterOutgoing    () const;
		/*virtual*/ size_t									UniqueScatterIncoming       ( size_t inidx ) const;
		/*virtual*/ size_t									UniqueScatterOutgoing       ( size_t outidx) const;
		/*virtual*/ size_t									IncomingRadianceIdx         ( size_t inidx ) const { return m_inidx  + inidx;}
		/*virtual*/ size_t									OutgoingRadianceIdx         ( size_t outidx) const { NXASSERT( (!m_isgroundpoint) ); return m_outidx + outidx;}
		/*virtual*/ size_t									GroundDownwardFluxIdx       ( size_t outidx) const { NXASSERT( (m_isgroundpoint)  ); return m_outidx + outidx;}	
		/*virtual*/ size_t									ScatterPropertyIdx          ( size_t outidx, size_t inidx ) const { NXASSERT( (!m_isgroundpoint)); return (m_scatvalidx + NumOutGoingRays()*inidx + outidx); }	// The analogous index for ground points is GroundDownwardFluxFactorsIdx -- This assert is left in to make sure developers know how this implementation works 
		/*virtual*/ size_t									GroundDownwardFluxFactorsIdx( size_t outidx) const { NXASSERT( (m_isgroundpoint)  ); return m_scatvalidx + outidx;}	
		/*virtual*/ bool									TriangulateOnOutgoing       ( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const;
};

