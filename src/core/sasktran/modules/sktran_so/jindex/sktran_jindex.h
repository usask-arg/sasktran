class SKTRAN_JIndexTableFactory_V2;
//class SKTRAN_EngineGeometry_V2;
//Added 2009/08/25 by Tony Bathgate
class SKTRANSO_JindexTableBase;

/*-----------------------------------------------------------------------------
 *					SKTRANSO_JIndex								2007-12-21*/
/** A class used to hold a descriptor used in the linear interpolation
 *	of the diffuse profile outbound from a unit sphere.  The class
 *	holds the index of the angle and radius of a diffuse point from the
 *	diffuse profile table as well as the index of the directional ray. The
 *	weight to be applied to the radiance is also stored
 **/
/*---------------------------------------------------------------------------*/

#pragma pack(push)
#pragma pack(2)
class SKTRANSO_JIndex
{
	private:
		double						m_weight;					//!< 8  (4) byte weight to be applied to the radiances from the three vertices on the unit sphere.
		uint32_t					m_radiigridindex;			//!< 2  (2) byte index for the radius in the grid. Indexes the diffuse radii in the diffuse profile table.
		uint16_t					m_szagridindex;				//!< 2  (1) byte index for the diffuse profiles in SZA in the diffuse profile table
		uint16_t					m_vertexindex;				//!< 2  (1) byte index for the vertex the 3 unit vectors on the unit sphere used for linear interpolation

	public:
		void						ConfigureLOSSolarTransmissionIndex		( SKTRAN_GridIndex quadpointindex, double weight );
		void						ConfigureDiffuseTableIndex				( SKTRAN_GridIndex positionindex,   SKTRAN_GridIndex heightindex, double szaweight, double radiiweight, size_t vertexindex, double vertexweight);
		void						ConfigureSolarTransmissionTableIndex	( SKTRAN_GridIndex positionindex,   SKTRAN_GridIndex heightindex, double szaweight, double radiiweight, size_t vertexindex, double vertexweight );
		void						ConfigureScatterMatrixTableIndex		( SKTRAN_GridIndex positionindex,   SKTRAN_GridIndex heightindex, double szaweight, double radiiweight, size_t vertexindex, double vertexweight );
		void						ConfigureGroundPointTableIndex			( SKTRAN_GridIndex groundpointindex, double pointweight, size_t vertexidx, double weight );					// get the point indexing, weights and the vertex index points to the incoming ray
		void						ConfigureEmissionTableIndex				( SKTRAN_GridIndex positionindex,   SKTRAN_GridIndex heightindex, double szaweight, double altweight, bool isground);

		SKTRAN_GridIndex			PositionIndex	()	const { return m_szagridindex;  }
		SKTRAN_GridIndex			HeightIndex		()	const { return m_radiigridindex;}
		SKTRAN_GridIndex			VertexIndex		()	const { return m_vertexindex;}
		double						VertexWeight	()	const { return m_weight;}
		bool						IsEmissionGround()  const;
		void						MultiplyWeightBy( double val ) { m_weight = (SKTRAN_GridIndexFloat)( m_weight*val);} 
};
#pragma pack(pop)


/*-----------------------------------------------------------------------------
 *					enum ENUM_SKTRAN_JSOURCE		2008-2-5*/
/** Specifies the tables that will be used for interpolation. 
 **/
/*---------------------------------------------------------------------------*/


/*
	enum	ENUM_SKTRAN_INTERPTABLE {	SKTRAN_TUNDEFINED,				//!< Reserved to place enumeration in "known" undefined state.
									SKTRAN_TDIFFUSEINBOUND,			//!< The inbound rays in the Diffuse Points table
									SKTRAN_TDIFFUSEOUTBOUND,		//!< The outbound source functions in the Diffuse Points table
									SKTRAN_TSOLARTRANS,				//!< The Solar transmission Table 
									SKTRAN_TDIFFUSEGROUNDPOINTS};	//!< The Diffuse Ground Points Table
*/

/*-----------------------------------------------------------------------------
 *					class SKTRAN_QuadratureJIndexTable_V2		2008-2-2*/
/** This is a class that stores indexes and weights to all of the diffuse
 *	radiances from various SASKTRAN tables that contribute to the incoming radiance
 *	seen by the origin of the ray.
 *
 *	The ray consists of NumCells cells and (NumCells+1) shells. Evaluation of the signal
 *	seen at the ray origin requires a quadrature across each cell.  The quadrature
 *	object will decide how many quadrature points at a distance "s" from the user are
 *	selected for each cell.  Typically, SASKTRAN only uses one point per cell
 *	but there is no reason this could not be increased to 3 or 5 points.
 *
 *	The Source function in the ray direction must be evaluated for each of the
 *	quadrature points must be evaluated by linear interpolation of the diffuse
 *	outbound radiances from the last iteration. This could use a linear combination
 *	of anywhere between 0 and 12 (or higher) individual outgoing radiances.  The
 *	code must keep track of how many individual radiances contribute to the evaluation
 *	of a single J and it must also keep track of how many J contribute to the quadrature
 *	for any given cell.
 *
 *	I have implemented this code so all of the index entries for a given ray
 *	are stored in contiguous memory.
 *
 **/
/*---------------------------------------------------------------------------*/

class SKTRANSO_JIndexArray
{
	private:
		uint32_t				m_numcells;			//!< The number of ray cells represented in this array
		uint32_t				m_numQ;				//!< The number of quadrature points represented in this array. One linear combo J per quadrature point
		uint32_t				m_numindices;		//!< The size of m_Jarray;
		uint32_t*				m_Qindex;			//!< Array [m_numQ     + 1] indicating where each Q starts in the list of m_Jarray;
		uint32_t*				m_firstQinCell;		//!< Array [m_numcells + 1] indicating the first Q (indexes m_Qindex ) in each cell
		SKTRANSO_JIndex*		m_Jarray;			//!< Array [m_numindices] of source function descriptors for this ray.

	#ifdef NXDEBUG
		size_t					m_maxcells;
		size_t					m_maxQ;
		size_t					m_maxindices;
	#endif

	private:
		void					ReleaseResources();

	public:
								SKTRANSO_JIndexArray();
		virtual				   ~SKTRANSO_JIndexArray();
		bool					DeepCopy					( const SKTRANSO_JIndexArray& other );
		bool					AllocateMaximumStorage		( size_t numcells, size_t numQ, size_t numindices);
		bool					ResetCounters				();
		bool					InsertQuadraturePointEntries( const SKTRANSO_JIndex* Jlinearsum, size_t numentries);
		bool					FinishCellEntries			();
		void					Clear						()		 { ReleaseResources();}
		size_t					NumCells					() const { return m_numcells;}
		size_t					NumQ						() const { return m_numQ;}
		size_t					NumIndices					() const { return m_numindices;}
		const SKTRANSO_JIndex*	At							( size_t idx) const { NXASSERT(( idx < m_numindices )); return m_Jarray + idx;}
		bool					IsEmpty						() const { return (m_numQ == 0);}
		bool					JIndicesAtQuadraturePoint	( size_t quadraturepointidx, size_t*                  Jindex, size_t* numentries) const;
		bool					JEntriesAtQuadraturePoint	( size_t quadraturepointidx, const SKTRANSO_JIndex** Jstart, size_t* numentries) const;
		bool					QuadraturePointsInCell		( size_t cellidx, size_t* startptindex, size_t* numquadraturepoints) const;
		void					DumpTable					() const;
		void					MultiplyAllWeightsBy		( double value );


};
/*-----------------------------------------------------------------------------
 *					SKTRAN_JValueTable_V21								2008-2-3*/
/** This class holds the optical processing loops counterpart for the
 *	SKTRAN_QuadratureJindexTable_V2. It stores pointers to all of the
 *	radiances that contribute to the source function terms at each quadrature
 *	point alomg each ray.
 *
 *	2012-09-14, ndl303, Changed the m_radiances and m_weights arrays to use
 *	std::vector. This provides easier debugging and better index checking.
 *	Code does rely upon the fact that std::vector produces contiguous elements.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_JValueTable_V21
{
	private:
		const SKTRANSO_JIndexArray*				m_jindex;		//!< The JIndex object, this allows us to index individual quadrature points	
		size_t										m_reservesize;	//!< The maximum number of elements pre-allocated
		size_t										m_numnonzero;	//!, The number of non-zero elements, sensibly set after calling AdjustWeights.
		std::vector<const SKTRAN_StokesScalar*>		m_radiances;	//!< Array[m_numindices] of pointers to the contributing radiances for this ray
		std::vector<SKTRAN_StokesScalar>			m_weights;		//!< Array[m_numindices] of weights for each radiance

	private:
		void							ReleaseResources();

	public:
										SKTRAN_JValueTable_V21();
		virtual						   ~SKTRAN_JValueTable_V21();
		bool							AttachToGeometry			 ( const SKTRANSO_JIndexArray& jindex );
		bool							ReserveStorage				 ( size_t maxpoints );
		const SKTRANSO_JIndexArray*	JIndexTable					 () const		{ return m_jindex;}
//		bool							JWeightsAtQuadraturePoint	 ( size_t quadraturepointidx, SKTRAN_StokesScalar** weights, size_t* numentries);
		// SetRadiancePtr and SetWeight are currently unused because of the changes I maded to  ConvertJIndexTable and adding SetWeightsAndRadiancePtrs - Tony B.
//		bool							SetRadiancePtr				 ( size_t idx, const SKTRAN_StokesScalar* radiance) { NXASSERT(( idx < m_jindex->NumIndices() )); m_radiances[idx] = radiance; return (radiance != NULL); }
//		bool							SetWeight					 ( size_t idx, SKTRAN_StokesScalar weight )         { NXASSERT(( idx < m_jindex->NumIndices() )); m_weights[idx]   = weight; return true;}
		bool							AdjustWeightsAndTrim		 ( const double* factor, size_t numquadraturepoints );
		bool							AdjustWeightsByConstantFactorAndTrim( double factor );
		void							Clear						();
		SKTRAN_StokesScalar				Evaluate					() const;
		void							DumpTable					() const;
		size_t							ReserveSize					() const { return m_reservesize;}
		size_t							size						() const { return m_numnonzero;}
		//Added 2008/08/25 -Tony Bathgate
		bool							SetWeightsAndRadiancePtrs	( const SKTRANSO_JindexTableBase* table, ENUM_SKTRAN_JSOURCE tablesource );
};


/*-----------------------------------------------------------------------------
 *					class SKTRANSO_JindexTableBase		2010-4-5*/
/** Defines a base class for interpolation by the Jindex formulation.**/
/*---------------------------------------------------------------------------*/

class SKTRANSO_JindexTableBase: public nxUnknown
{
	private:
		std::shared_ptr<const SKTRAN_RayFactory_Base>			m_rayfactory;		

	protected:
		std::shared_ptr<const SKTRAN_RayFactory_Base>			RayFactory					() const { NXASSERT(( m_rayfactory)); return m_rayfactory;}

	public:
																SKTRANSO_JindexTableBase		() {}
		virtual												   ~SKTRANSO_JindexTableBase		() {}
		bool													SetRayFactory				( std::shared_ptr<const SKTRAN_RayFactory_Base> rayfactory) { m_rayfactory = rayfactory; return true;}

		virtual bool											InterpolateTable				( const HELIODETIC_POINT&       location,				// Location where interpolation required
																								  const HELIODETIC_UNITVECTOR&	look,					// Look direction where interpolation required
																								  SKTRANSO_JIndex*				vertexdescriptortable, 
																								  size_t						maxpts, 
																								  size_t*						npts,
																								  double						weight ) const = 0;

		virtual const SKTRAN_StokesScalar*						ConvertJIndexToRadiancePtr		( const SKTRANSO_JIndex*		entry,
																								  ENUM_SKTRAN_JSOURCE			jsource
																								) const = 0;

};





