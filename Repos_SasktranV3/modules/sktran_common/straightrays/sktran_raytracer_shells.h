//#pragma once
//
//#include "sktran_common_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracer_Shells		2013-05-28*/
/** @ingroup rays
**/
/*---------------------------------------------------------------------------*/

class SKTRAN_RayTracer_Shells : public SKTRAN_RayTracer_Base
{
	private:
		std::shared_ptr< const SKTRAN_GridDefRayTracingShells_V21>  m_parentgrid;								//!< Ray tracing shell altitude grid

	private:
		bool										TraceRayInternal							( SKTRAN_RayOptical_Straight* ray ) const; // This needs to be virtual for curved ray tracer to work
		bool										AllocatePathElements						( size_t numelements, SKTRAN_RayOptical_Straight* straightray ) const;
		bool										TraceObserverINSIDE_LookingUp				( SKTRAN_RayOptical_Straight* ray ) const;
		bool										TraceObserverABOVE_LOSHitsGround			( SKTRAN_RayOptical_Straight* ray ) const;
		bool										TraceObserverINSIDE_LookingDownHitsGround	( SKTRAN_RayOptical_Straight* ray ) const;
		bool										TraceObserverINSIDE_LookingDownPassesThru	( SKTRAN_RayOptical_Straight* ray ) const;
		bool										TraceObserverABOVE_LOSPassesThrough			( SKTRAN_RayOptical_Straight* ray ) const;
		double										DistanceToTangentPoint_fromTrig				( double radius, double tangentradiussquared) const;

	public:
																SKTRAN_RayTracer_Shells			(std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords);
		virtual												   ~SKTRAN_RayTracer_Shells			();
		bool													Initialize						(std::shared_ptr< const SKTRAN_GridDefRayTracingShells_V21> raytracinggrid );
		const SKTRAN_GridDefRayTracingShells_V21*				RayTracingShellsPtr				() const {return m_parentgrid.get();}	// Only in here temporarily while I switch MC over to the new ray tracing style


	public:
		virtual bool								TraceStraightRay							( SKTRAN_RayOptical_Straight*			ray) const;

};

