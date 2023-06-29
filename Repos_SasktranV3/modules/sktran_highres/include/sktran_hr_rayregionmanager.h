//#include "sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_RayTracingRegionManager		2014-10-30*/
/** Modifications to the RayTracingRegion manager specific to the HR engine.
 *  Additional points are calculated, in and out reference points, which correspond
 *  to the average entry to the atmosphere and the average exit of the atmosphere
 *  for all lines of sight
 **/
/*---------------------------------------------------------------------------*/

class SKTRAN_HR_RayTracingRegionManager : public SKTRAN_RayTracingRegionManager
{
	private:
		nxVector							m_inreferencepoint;
		nxVector							m_outreferencepoint;

		double 								m_losminssa;
		double 								m_losmaxssa;

	private:
		virtual bool						UpdateBoundingReferences ( const SKTRAN_LineOfSightArray_V21& linesofsight );
		virtual bool 						UpdateLOSScatteringAngles( const SKTRAN_LineOfSightArray_V21& linesofsight );
		virtual bool						GetRayEndpointsObserverOutside				( const nxVector& observer, const nxVector& look, nxVector* startpt, nxVector* endpt);
	public:
		bool								GetBoundingReferences( nxVector& in, nxVector& out ) const;
		bool 								GetBoundingLOSScatteringAngles( double& minssa, double& maxssa ) const;
		virtual bool						UpdateUndefinedParametersFromLinesOfSight	( const SKTRAN_LineOfSightArray_V21& linesofsight );
};


