#include "../sktran_common.h"

SKTRAN_RayOptical_Curved::SKTRAN_RayOptical_Curved(std::unique_ptr<SKTRAN_RayStorage_CurvedPiecewise> trajectorystorage, std::shared_ptr<SKTRAN_RayTracer_Curved_Shells> raytracer)
{
	m_trajectorystorageobject = std::move(trajectorystorage);

	m_trajectorystorage = m_trajectorystorageobject.get();
	InitializeStorage(m_trajectorystorage);

	m_raytracer = raytracer;
}

SKTRAN_RayOptical_Curved::~SKTRAN_RayOptical_Curved()
{
}

bool SKTRAN_RayOptical_Curved::LocationAlongRayAsVector(const double& distancealongray, HELIODETIC_VECTOR * pt) const
{
	nxLog::Record(NXLOG_ERROR, "SKTRAN_RayOptical_Curved::LocationAlongRayAsVector Is not yet implemented");
	return false;
}

bool SKTRAN_RayOptical_Curved::GetQuadratureInterpParams_startPoint(size_t quadraturepoint, size_t cellidx, double * r0, double * t0, double * rt, HELIODETIC_POINT * startpoint) const
{
	bool ok = true;

	ok = m_trajectorystorage->LocationOfPoint(quadraturepoint, startpoint);
	*r0 = startpoint->Radius();
	*t0 = m_trajectorystorage->DistanceOfPointFromCellTangentPoint(quadraturepoint, cellidx);
	*rt = m_trajectorystorage->RadiusOfCellTangentPoint(cellidx);

	return ok;
}

bool SKTRAN_RayOptical_Curved::GetQuadratureInterpParams(size_t cellidx, double * r0, double * r1, double * t0, double * t1, double * rt, HELIODETIC_POINT * startpoint, HELIODETIC_POINT * endpoint) const
{
	bool ok = true;
	// do the end point first (equal to start point of next cell)
	ok = ok && GetQuadratureInterpParams_startPoint(cellidx + 1, cellidx, r1, t1, rt, endpoint);

	// now the start point
	ok = ok && GetQuadratureInterpParams_startPoint(cellidx, cellidx, r0, t0, rt, startpoint);

	return ok;
}

bool SKTRAN_RayOptical_Curved::TraceRay_NewMethod()
{
	return m_raytracer->TraceRay(this);
}

void SKTRAN_RayOptical_Curved::NotifyDerived_RayInvalid()
{
}
