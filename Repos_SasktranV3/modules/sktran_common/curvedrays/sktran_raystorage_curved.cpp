#include "../sktran_common.h"


void SKTRAN_RayStorage_CurvedPiecewise::CellTangentParams(const HELIODETIC_POINT& loc, const HELIODETIC_UNITVECTOR& look, double* tangentradius, double* disttotangent) const
{
	double coszenith = loc.CosZenithAngle(look);

	*tangentradius = loc.Radius() * sqrt(1 - coszenith * coszenith);

	// Not really the correct distance to the TP, but this is only used as a parameter in the extinction calculation which requires it
	// to be calculated this way
	*disttotangent = -1.0 * loc.Radius() * coszenith;
}

bool SKTRAN_RayStorage_CurvedPiecewise::PushBack(HELIODETIC_UNITVECTOR* uv, HELIODETIC_POINT* pt, SKTRAN_Distance celllength)
{
	m_locations.push_back(*pt);
	m_lookvectors.push_back(*uv);
	m_celllengths.push_back(celllength);

	double rt, dist;
	CellTangentParams(*pt, *uv, &rt, &dist);

	m_tangentradii.push_back(rt);
	m_disttotangent.push_back(dist);

	return true;
}

bool SKTRAN_RayStorage_CurvedPiecewise::Insert(HELIODETIC_UNITVECTOR* uv, HELIODETIC_POINT* pt, SKTRAN_Distance celllength, size_t index)
{
	m_locations.insert(m_locations.begin() + index, *pt);
	m_lookvectors.insert(m_lookvectors.begin() + index, *uv);
	m_celllengths.insert(m_celllengths.begin() + index, celllength);

	double rt, dist;
	CellTangentParams(*pt, *uv, &rt, &dist);

	m_tangentradii.insert(m_tangentradii.begin() + index, rt);
	m_disttotangent.insert(m_disttotangent.begin() + index, dist);

	return true;
}

bool SKTRAN_RayStorage_CurvedPiecewise::InitializeObserver(const HELIODETIC_VECTOR& observer, const HELIODETIC_UNITVECTOR& look)
{
	m_observer = observer;
	m_observerlookaway = look;
	
	return true;
}

void SKTRAN_RayStorage_CurvedPiecewise::ClearStorage()
{
	m_locations.clear();
	m_lookvectors.clear();
	m_celllengths.clear();
	m_tangentradii.clear();
	m_disttotangent.clear();
}

bool SKTRAN_RayStorage_CurvedPiecewise::Reserve(size_t numquadraturepoints)
{
	m_locations.reserve(numquadraturepoints);
	m_lookvectors.reserve(numquadraturepoints);
	m_celllengths.reserve(numquadraturepoints);
	m_tangentradii.reserve(numquadraturepoints);
	m_disttotangent.reserve(numquadraturepoints);

	return true;
}

bool SKTRAN_RayStorage_CurvedPiecewise::Resize(size_t numquadraturepoints)
{
	m_locations.resize(numquadraturepoints);
	m_lookvectors.resize(numquadraturepoints);
	m_celllengths.resize(numquadraturepoints);
	m_tangentradii.resize(numquadraturepoints);
	m_disttotangent.resize(numquadraturepoints);

	return true;
}

bool SKTRAN_RayStorage_CurvedPiecewise::SplitCell(size_t cellindex)
{
	HELIODETIC_POINT mid;

	CellMidPoint(cellindex, &mid);

	double newds = m_celllengths[cellindex] / 2.0;

	Insert(&m_lookvectors[cellindex], &mid, newds, cellindex);

	m_celllengths[cellindex + 1] = newds;

	return true;
}

void SKTRAN_RayStorage_CurvedPiecewise::TruncateToNumElements(size_t numels)
{
	Resize(numels);
}

size_t SKTRAN_RayStorage_CurvedPiecewise::NumCells() const
{
	if (m_locations.size() == 0)
	{
		return 0;
	}
	else
	{
		return m_locations.size() - 1;
	}
}

size_t SKTRAN_RayStorage_CurvedPiecewise::NumQuadraturePoints() const
{
	return m_locations.size();
}

double SKTRAN_RayStorage_CurvedPiecewise::RadiusOfPoint(size_t quadraturepoint_index) const
{
	return m_locations[quadraturepoint_index].Radius();
}

double SKTRAN_RayStorage_CurvedPiecewise::AltitudeOfPoint(size_t quadraturepoint_index) const
{
	return m_locations[quadraturepoint_index].Altitude();
}

bool SKTRAN_RayStorage_CurvedPiecewise::LocationOfPoint(size_t quadraturepoint_index, HELIODETIC_POINT * pt) const
{
	*pt = m_locations[quadraturepoint_index];

	return true;
}

double SKTRAN_RayStorage_CurvedPiecewise::DistanceOfPointFromOrigin(size_t quadraturepoint_index) const
{
	double distance = 0.0;

	for (size_t idx = 0; idx < quadraturepoint_index; idx++)
	{
		distance += m_celllengths[idx];
	}

	return distance;
}

double SKTRAN_RayStorage_CurvedPiecewise::DistanceOfPointFromCellTangentPoint(size_t quadraturepoint_index, size_t cellindex) const
{
	// TODO: Might need to adjust based on cell index for curved rays?
	double rt, T;

	CellTangentParams(m_locations[quadraturepoint_index], m_lookvectors[cellindex], &rt, &T);

	return T;
}

double SKTRAN_RayStorage_CurvedPiecewise::RadiusOfCellTangentPoint(size_t cellindex) const
{
	return m_tangentradii[cellindex];
}

HELIODETIC_UNITVECTOR SKTRAN_RayStorage_CurvedPiecewise::AverageLookVectorAwayFromObserver(size_t cellindex) const
{
	return m_lookvectors[cellindex];
}

HELIODETIC_UNITVECTOR SKTRAN_RayStorage_CurvedPiecewise::AverageLookVectorTowardsObserver(size_t cellindex) const
{
	HELIODETIC_UNITVECTOR out = m_lookvectors[cellindex];

	out.Negate();

	return out;
}

double SKTRAN_RayStorage_CurvedPiecewise::CellLength(size_t cellindex) const
{
	return m_celllengths[cellindex];
}

bool SKTRAN_RayStorage_CurvedPiecewise::CellMidPoint(size_t cellindex, HELIODETIC_POINT * pt) const
{
	// TODO: More accurate is to account for the cell curvature when calculating the midpoint

	// This assumes a piecewise line
	HELIODETIC_VECTOR start = m_locations[cellindex].Vector();
	HELIODETIC_VECTOR end = m_locations[cellindex + 1].Vector();
	HELIODETIC_VECTOR middle = (start + end) * 0.5;

	pt->FromVector(middle, GetCoordsObject().get());

	return true;
}

double SKTRAN_RayStorage_CurvedPiecewise::CellCurvature(size_t cellindex) const
{
	HELIODETIC_VECTOR diff = (m_locations[cellindex + 1].Vector() - m_locations[cellindex].Vector());
	double fraction = m_celllengths[cellindex] / diff.Magnitude();
	return fraction;
}

