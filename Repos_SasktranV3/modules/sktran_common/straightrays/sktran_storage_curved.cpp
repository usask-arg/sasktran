
#include "../sktran_common.h"



SKTRAN_RayStorage_CurvedPiecewise_MC::SKTRAN_RayStorage_CurvedPiecewise_MC( std::shared_ptr<const SKTRAN_CoordinateTransform_V2> coords) : SKTRAN_RayStorage_CurvedPiecewise(coords) 
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC not yet implemented");
}

bool SKTRAN_RayStorage_CurvedPiecewise_MC::PushBack( HELIODETIC_UNITVECTOR* uv, HELIODETIC_POINT* pt, SKTRAN_Distance celllength)
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::PushBack not yet implemented");
	return false;
}

bool SKTRAN_RayStorage_CurvedPiecewise_MC::Insert( HELIODETIC_UNITVECTOR* uv, HELIODETIC_POINT* pt, SKTRAN_Distance celllength, size_t index)
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::Insert not yet implemented");
	return false;
}

bool SKTRAN_RayStorage_CurvedPiecewise_MC::InitializeObserver( const HELIODETIC_VECTOR&    observer, const HELIODETIC_UNITVECTOR& look )
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::InitializeObserver not yet implemented");
	return false;
}

void SKTRAN_RayStorage_CurvedPiecewise_MC::ClearStorage()
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::ClearStorage not yet implemented");
}

bool SKTRAN_RayStorage_CurvedPiecewise_MC::Reserve( size_t numquadraturepoints )
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::Reserve not yet implemented");
	return false;
}

//bool SKTRAN_RayStorage_CurvedPiecewise_MC::LocationAlongRayAsPoint( double shellintercept, HELIODETIC_POINT*  pt    )const
//{
//	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::LocationAlongRayAsPoint not yet implemented");
//	return false;
//}
//
//bool SKTRAN_RayStorage_CurvedPiecewise_MC::LocationAlongRayAsVector( double shellintercept, HELIODETIC_VECTOR* pt    )	const
//{
//	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::LocationAlongRayAsVector not yet implemented");
//	return false;
//}
//
//bool SKTRAN_RayStorage_CurvedPiecewise_MC::GetCellIntercepts( size_t raysegment_idx, SKTRAN_Distance* startintercept, SKTRAN_Distance*exitintercept ) const
//{
//	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::GetCellIntercepts not yet implemented");
//	return false;
//
//}
//
//HELIODETIC_VECTOR SKTRAN_RayStorage_CurvedPiecewise_MC::EndPoint() const
//{
//	HELIODETIC_VECTOR	l;
//	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::EndPoint not yet implemented");
//	return l;
//}
//
//SKTRAN_Distance SKTRAN_RayStorage_CurvedPiecewise_MC::EndIntercept() const
//{
//	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::EndIntercept not yet implemented");
//	return std::numeric_limits<double>::quiet_NaN();
//
//}

HELIODETIC_UNITVECTOR SKTRAN_RayStorage_CurvedPiecewise_MC::AverageLookVectorAwayFromObserver( size_t raysegment_index ) const
{
	HELIODETIC_UNITVECTOR	l;
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::AverageLookVectorAwayFromObserver not yet implemented");
	return l;
}

HELIODETIC_UNITVECTOR SKTRAN_RayStorage_CurvedPiecewise_MC::AverageLookVectorTowardsObserver( size_t raysegment_index ) const
{
	HELIODETIC_UNITVECTOR	l;
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::AverageLookVectorTowardsObserver not yet implemented");
	return l;
}

bool SKTRAN_RayStorage_CurvedPiecewise_MC::Resize( size_t numquadraturepoints )
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::Resize not yet implemented");
	return false;
}

double SKTRAN_RayStorage_CurvedPiecewise_MC::DistanceOfPointFromOrigin( size_t quadraturepoint_index) const
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::DistanceOfPointFromOrigin not yet implemented");
	return false;
}

double SKTRAN_RayStorage_CurvedPiecewise_MC::CellLength( size_t quadraturepoint_index) const
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::CellLength not yet implemented");
	return 0.0;
}
bool SKTRAN_RayStorage_CurvedPiecewise_MC::CellMidPoint( size_t cell, HELIODETIC_POINT* pt ) const
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::CellMidPoint not yet implemented");
	pt->Clear();
	return false;
}

bool SKTRAN_RayStorage_CurvedPiecewise_MC::LocationOfPoint( size_t quadraturepointindex, HELIODETIC_POINT* pt ) const
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::LocationOfPoint not yet implemented");
	pt->Clear();
	return false;
}

void SKTRAN_RayStorage_CurvedPiecewise_MC::TruncateToNumElements( size_t numels )
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::TruncateToNumElements not yet implemented");

}

double SKTRAN_RayStorage_CurvedPiecewise_MC::DistanceOfPointFromCellTangentPoint( size_t quadraturepoint_index, size_t raysegment_index ) const
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::DistanceFromCellTangentPoint not yet implemented");
	return std::numeric_limits<double>::quiet_NaN();
}

double SKTRAN_RayStorage_CurvedPiecewise_MC::AltitudeOfPoint( size_t quadraturepoint_index ) const
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::AltitudeOfQuadraturePoint not yet implemented");
	return std::numeric_limits<double>::quiet_NaN();
}

double SKTRAN_RayStorage_CurvedPiecewise_MC::RadiusOfPoint( size_t quadraturepoint_index ) const 
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::RadiusOfQuadraturePoint not yet implemented");
	return std::numeric_limits<double>::quiet_NaN();
}
double SKTRAN_RayStorage_CurvedPiecewise_MC::RadiusOfCellTangentPoint( size_t raysegment_index ) const
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_HR::RadiusOfCellTangentPoint not yet implemented");
	return std::numeric_limits<double>::quiet_NaN();
}

size_t SKTRAN_RayStorage_CurvedPiecewise_MC::NumCells() const
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::NumCells() not yet implemented");
	return 0;
}

size_t SKTRAN_RayStorage_CurvedPiecewise_MC::NumQuadraturePoints() const
{
	nxLog::Record( NXLOG_INFO, "class SKTRAN_RayStorage_CurvedPiecewise_MC::NumQuadraturePoints() not yet implemented");
	return 0;
}
