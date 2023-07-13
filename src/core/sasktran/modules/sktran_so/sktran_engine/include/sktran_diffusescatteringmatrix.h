
class SKTRAN_ScatteringMatrixPointGeometry_V21;

/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointOptical_V21		2008-1-10*/
/**	This is a class involved in scattering the light from the incoming 
 *	directions to the outgoing directions at one of the Diffuse Points
 *	in the diffuse table. This class has potential to rapidly exhaust memory as it
 *	stores an internal scattering array of (numincoming x numoutgoing) and can
 *	quickly get to the 10s or 100s of megabyte size.
 *
 *	The class is closely coupled to the incoming lines of sight into a  Diffuse
 *	Point and the number of outbound points on the unit sphere.
 *
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_ScatteringMatrixPointOptical_V21
{
	private:
		const SKTRAN_ScatteringMatrixPointGeometry_V21*			m_scatmatrixgeometry;	//!< The geometry aspects of this scattering matrix
		size_t													m_numincoming;			//!< The number of incoming rays (sames as the value in the geometry)
//		size_t													m_numoutgoing;			//!< The number of outgoing directions
		SKTRAN_PhaseMatrixScalar*								m_scatterarray;			//!< Array [numincoming x numoutgoing] Storage for the scattering of each incoming ray to each outgoing direction

	private:
		void											ReleaseResources();
		bool											ComputeMultipliersAndAdjustScatterArray( const SKTRAN_TableOpticalProperties_V21* opticalprops );
		SKTRAN_StokesScalar								ComputeMeanPhaseFunctionAndAdjustScatterArray( size_t inidx,const SKTRAN_TableOpticalProperties_V21* opticalprops );

	public:
														SKTRAN_ScatteringMatrixPointOptical_V21();
		virtual										   ~SKTRAN_ScatteringMatrixPointOptical_V21();
		void											SetGeometry			( const SKTRAN_ScatteringMatrixPointGeometry_V21* geometry);
		bool											AttachToGeometry	( );
		bool											ConfigureOptical	( const SKTRAN_TableOpticalProperties_V21* opticalprops );
		bool											ScatterIncomingRays	( SKTRAN_DiffusePointOptical_V21* point );
		size_t											NumIncoming			() const;
		size_t											NumOutgoing			() const;
};

/*-----------------------------------------------------------------------------
 *					SKTRAN_ScatteringMatrixPointGeometry_V21		2008-1-10*/
/** This is a class involved in caching the geometry associated with
 *	the scattering of light from incoming directions to outgoing directions for
 *	points in the diffuse table.
 *
 *	This class caches the incoming azimuth and zenith angle grids as well as the
 *	outbound unit sphere grid. In addition it stores the differential solid
 *	angles required when integrating over the incoming rays. The class is
 *	closely coupled to its optically dependent counterpart
 *	#SKTRAN_ScatteringMatrixPointOptical_V21.
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_ScatteringMatrixPointGeometry_V21
{
	private:
		HELIODETIC_POINT									m_location;				//!< The radius in meters of this shell, used to look up optical properties
		const SKTRAN_UnitSphere_V2*							m_unitsphere;			//!< The unit sphere defining all of the outbound radiance directions.
		const SKTRAN_UnitSphereLatLonGrid*					m_incomingunitsphere;	//!< The unit sphere definiing all of the inbound lines of sight (note the unit vector point outwards)

	private:
		void												ReleaseResources		();

	public:
															SKTRAN_ScatteringMatrixPointGeometry_V21();
		virtual											   ~SKTRAN_ScatteringMatrixPointGeometry_V21();

		bool												ConfigureGeometry_Stage1	( const HELIODETIC_POINT&                             location,
																						  const SKTRAN_UnitSphereLatLonGrid*				  incomingunitsphere,
																						  const SKTRAN_UnitSphere_V2*						  outboundunitsphere );

		bool												ConfigureGeometry_Stage2MT	( );
		const SKTRAN_UnitSphere_V2*							OutboundUnitSphere			()  const { return m_unitsphere;}
		const SKTRAN_UnitSphereLatLonGrid*					InboundUnitSphere			()  const { return m_incomingunitsphere;}
		const HELIODETIC_POINT&								Location					()  const { return m_location;} 
};
