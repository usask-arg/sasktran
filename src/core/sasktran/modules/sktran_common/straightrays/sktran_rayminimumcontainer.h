


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Base		2014-1-24*/
/** @ingroup rays
 * Caching for each point in the ray is done inside this container. 
 * The container is configured to cache only what is needed to
 * perform a radiative transfer calculation -- e.g. if we're tracing
 * through shells in a 1D height atmosphere we only need to cache
 * the point radius and distance from the tangent point; if the 
 * atmosphere is 3D we might also need to cache a unit vector to
 * each point. 
 *
 * Right now I basically have vectors for everything I *might* want
 * to store, and these are filled or left empty depending on which
 * "fill" functions are set at configuration time. Likewise for accessors.
 *
 * Inheritance from this will be messy, so try to just use function pointers
 * or functors. 
 * 
 * Consider using std::bind in the future. This was brought into the standard
 * with C++11, but I'm not using it because the Westgrid c++ compiler is ancient. 
 * Apparently this allows delegates to be implemented much more nicely than what
 * I've done below. 
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayStorage_Base
{
	private:
		bool													m_groundishit;
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	m_coords;
		
	public:
																SKTRAN_RayStorage_Base				( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords);
		virtual												   ~SKTRAN_RayStorage_Base				();
		const SKTRAN_CoordinateTransform_V2*					GetCoordsPtr						() const			{ return m_coords.get();}
		std::shared_ptr<const SKTRAN_CoordinateTransform_V2>	GetCoordsObject						() const			{ return m_coords;}
		bool													GroundIsHit							() const			{ return m_groundishit;}
		void													SetGroundIsHit						( bool groundishit) { m_groundishit = groundishit;}

	public:
		virtual	bool											InitializeObserver					( const HELIODETIC_VECTOR&    observer, const HELIODETIC_UNITVECTOR& look ) = 0;
		virtual void											ClearStorage						() = 0;
		virtual bool											Reserve								( size_t numquadraturepoints ) = 0;
		virtual bool											Resize								( size_t numquadraturepoints ) = 0;
		virtual bool											SplitCell                           ( size_t cellindex ) = 0;
		virtual void											TruncateToNumElements				( size_t numels ) = 0;

		virtual size_t											NumCells							() const = 0;
		virtual size_t											NumQuadraturePoints					() const = 0;
		virtual double											RadiusOfPoint						( size_t quadraturepoint_index ) const = 0;
		virtual double											AltitudeOfPoint						( size_t quadraturepoint_index ) const = 0;							//{ return m_coords->RadiusToAltitude(m_radii[quadraturepoint_index]);}
		virtual bool											LocationOfPoint						( size_t quadraturepoint_index, HELIODETIC_POINT*  pt)	const = 0;
		virtual double											DistanceOfPointFromOrigin			( size_t quadraturepoint_index ) const = 0;
		virtual double											DistanceOfPointFromCellTangentPoint	( size_t quadraturepoint_index, size_t cellindex ) const  = 0;
		virtual double											RadiusOfCellTangentPoint			( size_t cellindex ) const = 0;
		virtual HELIODETIC_UNITVECTOR							AverageLookVectorAwayFromObserver	( size_t cellindex ) const = 0;
		virtual HELIODETIC_UNITVECTOR							AverageLookVectorTowardsObserver	( size_t cellindex ) const = 0;
		virtual double											CellLength							( size_t cellindex ) const = 0;
		virtual bool											CellMidPoint						( size_t cellindex, HELIODETIC_POINT* pt ) const = 0;

		// Ratio of actual cell length to geometric line length
		virtual double											CellCurvature						( size_t cellindex ) const = 0;

		// these are not required to implement 
		virtual double											ExtinctionAtCellStart				( size_t cellindex ) const { return -1; }
		virtual void											SetExtinction						( size_t cellindex, double ext ) const { return; }

		// these are not required to implement HR,SO, MC; used in the TIR engine in the calculation of analytic weighting functions
		virtual double											WFAtPoint							( const CLIMATOLOGY_HANDLE& wf_species, size_t quadraturepoint_index ) const { return -1; }
		virtual void											SetWFAtPoint						( const CLIMATOLOGY_HANDLE& wf_species, size_t quadraturepoint, double value ) const { return; }
		virtual void											AddWFSpecies						( const CLIMATOLOGY_HANDLE& wf_species ) const { return; }
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight		 2014- 11- 20 */
/** @ingroup rays
**/
/* --------------------------------------------------------------------------- */

class SKTRAN_RayStorage_Straight : public SKTRAN_RayStorage_Base
{
	private:
		HELIODETIC_VECTOR										m_observer;					//!< The location of the obserevr
		HELIODETIC_UNITVECTOR									m_lookaway;					//!< The observers look direction away from observer
		HELIODETIC_UNITVECTOR									m_looktowards;				//!< The observers look direction towards observer
		std::vector<SKTRAN_Distance>							m_distancefromorigin;			//!< The array[m_numshellintercepts] of shell intercepts, expressed as distance from oberver, sorted in ascending order.
		std::vector<SKTRAN_Distance>							m_distFromTan;
		std::vector<SKTRAN_Distance>							m_radii;
		double													m_Rt;						//!< The radius of the tangent point.
		double													m_Robs;						//!< The radius of the observer.
		double													m_Tobs;						//!< The distance of tangent point from observer, +ve parallel to line of sight
		double													m_geoidradius;

	private:
		void													ReleaseResources					();
		bool													CalculateTangentPointDetails		();

	public:
																SKTRAN_RayStorage_Straight			( std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords);
		virtual												   ~SKTRAN_RayStorage_Straight			();
		const HELIODETIC_VECTOR&								GetObserver							() const { return m_observer;}					//!< The location of the obserevr
		const HELIODETIC_UNITVECTOR&							LookVector							() const { return m_lookaway;}					//!< The observers look direction
		double													Rt									() const { return m_Rt;}
		double													Robs								() const { return m_Robs;}
		double													Tobs								() const { return m_Tobs;}
		virtual bool											PushBack							( SKTRAN_Distance r, SKTRAN_Distance distFromTan, SKTRAN_Distance s) ;
		virtual bool											Insert								( SKTRAN_Distance r, SKTRAN_Distance distFromTan, SKTRAN_Distance s, size_t index);

	public:
		virtual	bool											InitializeObserver					( const HELIODETIC_VECTOR&    observer, const HELIODETIC_UNITVECTOR& look ) override;
		virtual void											ClearStorage						() override;
		virtual bool											Reserve								( size_t numquadraturepoints ) override;
		virtual bool											Resize								( size_t numquadraturepoints ) override;
		virtual bool											SplitCell                           ( size_t cellindex ) override;
		virtual void											TruncateToNumElements				( size_t numels ) override;
		virtual HELIODETIC_UNITVECTOR							AverageLookVectorAwayFromObserver	( size_t raysegment_index ) const override { return m_lookaway;}
		virtual HELIODETIC_UNITVECTOR							AverageLookVectorTowardsObserver	( size_t raysegment_index ) const override { return m_looktowards;}
		virtual double											DistanceOfPointFromCellTangentPoint	( size_t quadraturepoint_index, size_t cellindex ) const  override;	//{ return m_distFromTan[quadraturepoint_index];}
		virtual double											DistanceOfPointFromOrigin			( size_t quadraturepoint_index ) const  override;	//{ return m_distFromTan[quadraturepoint_index];}
		virtual double											RadiusOfPoint						( size_t quadraturepoint_index ) const override;							//{ return m_radii[quadraturepoint_index];
		virtual double											AltitudeOfPoint						( size_t quadraturepoint_index ) const override;							//{ return m_coords->RadiusToAltitude(m_radii[quadraturepoint_index]);}
		virtual bool											LocationOfPoint						( size_t quadraturepoint_index, HELIODETIC_POINT*  pt)	const override;
		virtual double											RadiusOfCellTangentPoint			( size_t cellindex ) const  override { return m_Rt;}//{ return m_distFromTan[quadraturepoint_index];}
		virtual size_t											NumQuadraturePoints					() const override { return m_distancefromorigin.size();}
		virtual size_t											NumCells							() const override { return (m_distancefromorigin.size() > 0) ?  m_distancefromorigin.size() - 1 : 0;}
		virtual double											CellLength							( size_t cellindex ) const override;
		virtual bool											CellMidPoint						( size_t cellindex, HELIODETIC_POINT* pt) const override;
		virtual double											CellCurvature						( size_t cellindex ) const override { return 1.0; }


};

/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_CurvedPiecewise		 2014- 11- 20*/
/** @ingroup rays
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayStorage_CurvedPiecewise : public SKTRAN_RayStorage_Base
{
	private:
		std::vector<HELIODETIC_POINT>							m_locations;
		std::vector<HELIODETIC_UNITVECTOR>						m_lookvectors;
		std::vector<SKTRAN_Distance>							m_celllengths;
		std::vector<SKTRAN_Distance>							m_tangentradii;
		std::vector<SKTRAN_Distance>							m_disttotangent;

		HELIODETIC_VECTOR										m_observer;
		HELIODETIC_UNITVECTOR									m_observerlookaway;
	private:
		void													CellTangentParams(const HELIODETIC_POINT& loc, const HELIODETIC_UNITVECTOR& look, double* tangentradius, double* disttotangent) const;

	public:
																SKTRAN_RayStorage_CurvedPiecewise(  std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords) : SKTRAN_RayStorage_Base(coords) {}
		virtual												   ~SKTRAN_RayStorage_CurvedPiecewise(){}

		virtual bool											PushBack							( HELIODETIC_UNITVECTOR* uv, HELIODETIC_POINT* pt, SKTRAN_Distance celllength );
		virtual bool											Insert								( HELIODETIC_UNITVECTOR* uv, HELIODETIC_POINT* pt, SKTRAN_Distance celllength, size_t index);

		virtual	bool											InitializeObserver(const HELIODETIC_VECTOR&    observer, const HELIODETIC_UNITVECTOR& look);
		virtual void											ClearStorage();
		virtual bool											Reserve(size_t numquadraturepoints);
		virtual bool											Resize(size_t numquadraturepoints);
		virtual bool											SplitCell(size_t cellindex) override;
		virtual void											TruncateToNumElements(size_t numels);

		virtual size_t											NumCells() const;
		virtual size_t											NumQuadraturePoints() const;
		virtual double											RadiusOfPoint(size_t quadraturepoint_index) const;
		virtual double											AltitudeOfPoint(size_t quadraturepoint_index) const;							//{ return m_coords->RadiusToAltitude(m_radii[quadraturepoint_index]);}
		virtual bool											LocationOfPoint(size_t quadraturepoint_index, HELIODETIC_POINT*  pt)	const;
		virtual double											DistanceOfPointFromOrigin(size_t quadraturepoint_index) const;
		virtual double											DistanceOfPointFromCellTangentPoint(size_t quadraturepoint_index, size_t cellindex) const;
		virtual double											RadiusOfCellTangentPoint(size_t cellindex) const;
		virtual HELIODETIC_UNITVECTOR							AverageLookVectorAwayFromObserver(size_t cellindex) const;
		virtual HELIODETIC_UNITVECTOR							AverageLookVectorTowardsObserver(size_t cellindex) const;
		virtual double											CellLength(size_t cellindex) const;
		virtual bool											CellMidPoint(size_t cellindex, HELIODETIC_POINT* pt) const;
		virtual double											CellCurvature(size_t cellindex) const override;


		// these are not required to implement 
		virtual double											ExtinctionAtCellStart(size_t cellindex) const { return -1; }
		virtual void											SetExtinction(size_t cellindex, double ext) const { return; }
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight_SO		 2014- 11- 20*/
/** @ingroup rays
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayStorage_Straight_SO : public SKTRAN_RayStorage_Straight
{
	public:
		SKTRAN_RayStorage_Straight_SO(  std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords) : SKTRAN_RayStorage_Straight(coords) {}
};




/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight_HR		 2014- 11- 20*/
/** @ingroup rays
 *	In HR rays are only held by one object at a time ( the thread ) and thus
 *  the mutable inside on m_extinction is 'safe' but ugly
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayStorage_Straight_HR : public SKTRAN_RayStorage_Straight
{
	private:
		mutable std::vector<double>								m_extinction;
	public:
		SKTRAN_RayStorage_Straight_HR( std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords) : SKTRAN_RayStorage_Straight(coords) {}
		virtual bool											Reserve								( size_t numquadraturepoints ) override
		{
			m_extinction.reserve( numquadraturepoints );
			return SKTRAN_RayStorage_Straight::Reserve( numquadraturepoints );
		}
		virtual bool											Resize								( size_t numquadraturepoints ) override
		{
			m_extinction.resize( numquadraturepoints );
			return SKTRAN_RayStorage_Straight::Resize( numquadraturepoints );
		}
		virtual void											TruncateToNumElements				( size_t numels ) override
		{
			//m_extinction.resize(     std::min( numels, m_extinction.size() ) );
			m_extinction.resize( numels );
			return SKTRAN_RayStorage_Straight::TruncateToNumElements( numels );
		}
		virtual void											ClearStorage						() override
		{
			m_extinction.clear();
			return SKTRAN_RayStorage_Straight::ClearStorage();
		}
		virtual bool											PushBack							( SKTRAN_Distance r, SKTRAN_Distance distFromTan, SKTRAN_Distance s) override
		{
			m_extinction.push_back(-1);
			return SKTRAN_RayStorage_Straight::PushBack( r, distFromTan, s );
		}
		virtual bool											Insert								( SKTRAN_Distance r, SKTRAN_Distance distFromTan, SKTRAN_Distance s, size_t index) override
		{
			m_extinction.insert( std::begin(m_extinction) + index, -1);
			return SKTRAN_RayStorage_Straight::Insert( r, distFromTan, s, index );
		}
		virtual double											ExtinctionAtCellStart				( size_t cellindex ) const { return m_extinction[cellindex]; }
		virtual void											SetExtinction						( size_t cellindex, double ext ) const { m_extinction[cellindex] = ext; }
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_CurvedPiecewise_HR		 2014- 11- 20*/
/** @ingroup rays
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayStorage_CurvedPiecewise_HR : public SKTRAN_RayStorage_CurvedPiecewise
{
private:
	mutable std::vector<double>								m_extinction;
public:
	SKTRAN_RayStorage_CurvedPiecewise_HR(std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords) : SKTRAN_RayStorage_CurvedPiecewise(coords) {}
	virtual bool											Reserve(size_t numquadraturepoints) override
	{
		m_extinction.reserve(numquadraturepoints);
		return SKTRAN_RayStorage_CurvedPiecewise::Reserve(numquadraturepoints);
	}
	virtual bool											Resize(size_t numquadraturepoints) override
	{
		m_extinction.resize(numquadraturepoints);
		return SKTRAN_RayStorage_CurvedPiecewise::Resize(numquadraturepoints);
	}
	virtual void											TruncateToNumElements(size_t numels) override
	{
		//m_extinction.resize(     std::min( numels, m_extinction.size() ) );
		m_extinction.resize(numels);
		return SKTRAN_RayStorage_CurvedPiecewise::TruncateToNumElements(numels);
	}
	virtual void											ClearStorage() override
	{
		m_extinction.clear();
		return SKTRAN_RayStorage_CurvedPiecewise::ClearStorage();
	}
	virtual bool											PushBack(HELIODETIC_UNITVECTOR* uv, HELIODETIC_POINT* pt, SKTRAN_Distance celllength) override
	{
		m_extinction.push_back(-1);
		return SKTRAN_RayStorage_CurvedPiecewise::PushBack(uv, pt, celllength);
	}
	virtual bool											Insert(HELIODETIC_UNITVECTOR* uv, HELIODETIC_POINT* pt, SKTRAN_Distance celllength, size_t index) override
	{
		m_extinction.insert(std::begin(m_extinction) + index, -1);
		return SKTRAN_RayStorage_CurvedPiecewise::Insert(uv, pt, celllength, index);
	}
	virtual double											ExtinctionAtCellStart(size_t cellindex) const { return m_extinction[cellindex]; }
	virtual void											SetExtinction(size_t cellindex, double ext) const { m_extinction[cellindex] = ext; }
};



/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_Straight_MC		 2014- 11- 20*/
/** @ingroup rays
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayStorage_Straight_MC : public SKTRAN_RayStorage_Straight
{
	public:
											SKTRAN_RayStorage_Straight_MC( std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords) : SKTRAN_RayStorage_Straight(coords) {}
};


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayStorage_CurvedPiecewise_MC		 2014- 11- 20*/
/** @ingroup rays
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayStorage_CurvedPiecewise_MC : public SKTRAN_RayStorage_CurvedPiecewise
{
	public:
											SKTRAN_RayStorage_CurvedPiecewise_MC( std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords);
										   ~SKTRAN_RayStorage_CurvedPiecewise_MC(){}

	// -- implement methods inherited from SKTRAN_RayStorage_CurvedPiecewise

	public:
		virtual bool						PushBack							( HELIODETIC_UNITVECTOR* uv, HELIODETIC_POINT* pt, SKTRAN_Distance celllength) override;
		virtual bool						Insert								( HELIODETIC_UNITVECTOR* uv, HELIODETIC_POINT* pt, SKTRAN_Distance celllength, size_t index ) override;
	
	// -- implement methods inherited from SKTRAN_RayStorage_Base

	public:
		virtual	bool						InitializeObserver					( const HELIODETIC_VECTOR&    observer, const HELIODETIC_UNITVECTOR& look ) override;
		virtual void						ClearStorage						() override;
		virtual bool						Reserve								( size_t numquadraturepoints ) override;
		virtual bool						Resize								( size_t numquadraturepoints ) override;
		virtual HELIODETIC_UNITVECTOR		AverageLookVectorAwayFromObserver	( size_t raysegment_index ) const override;
		virtual HELIODETIC_UNITVECTOR		AverageLookVectorTowardsObserver	( size_t raysegment_index ) const override;
		virtual void						TruncateToNumElements				( size_t numels ) override;
		virtual double						DistanceOfPointFromOrigin			( size_t quadraturepoint_index) const override;
		virtual double						DistanceOfPointFromCellTangentPoint	( size_t quadraturepoint_index, size_t cellindex ) const  override;
		virtual double						AltitudeOfPoint						( size_t quadraturepoint_index ) const override;
		virtual double						RadiusOfPoint						( size_t quadraturepoint_index ) const override;
		virtual bool						LocationOfPoint						( size_t quadraturepoint_index, HELIODETIC_POINT*  pt)	const override;
		virtual double						RadiusOfCellTangentPoint			( size_t cellindex ) const override;							//{ return m_radii[quadraturepoint_index];}
		virtual size_t						NumCells							() const;
		virtual size_t						NumQuadraturePoints					() const;
		virtual double						CellLength							( size_t cellindex ) const override;
		virtual bool						CellMidPoint						( size_t cellindex, HELIODETIC_POINT* pt ) const override;

};

