#include "../sktran_common.h"

SKTRAN_UnitSphere_Plane::SKTRAN_UnitSphere_Plane()
{
	
}

bool SKTRAN_UnitSphere_Plane::ConstructPlane  ( const std::vector<nxVector>& unitvecs, size_t refindex ) 
{
	bool ok = true;
	std::vector< std::pair<double, size_t> > sortedangles;
	// copy the unit vectors
	ok = ok && AllocateVertices( unitvecs.size() );
	m_reference = unitvecs[refindex];

	// use the first and last point to make the plane
	nxVector look = (unitvecs[unitvecs.size()-1] - unitvecs[0]).UnitVector();

	m_normal = look.Cross( m_reference ).UnitVector();
	m_x = m_reference;
	m_y = m_reference.Cross(m_normal);
	// convert to angles
	m_angles.AllocateGridArray( unitvecs.size() );
	sortedangles.resize( unitvecs.size() );
	for( size_t idx = 0; idx < unitvecs.size(); idx++ )
	{
		double anglenotproj;
		double angleproj;
		anglenotproj = nxmath::atan2d( unitvecs[idx] & m_y, unitvecs[idx] & m_x );
		angleproj    = ProjectedAngle( unitvecs[idx] );
		if( abs( (anglenotproj - angleproj ) ) > 1E-5 )
		{
			nxLog::Record(NXLOG_WARNING, "Warning, entered points in SKTRAN_UnitSphere_Plane::ConstructPlane are not coplanar");
		}
		sortedangles[idx].first = angleproj;
		sortedangles[idx].second = idx;
	}
	std::sort( std::begin( sortedangles ), std::end( sortedangles ), [](std::pair< double, size_t > p1, std::pair< double, size_t > p2) { return p1.first < p2.first; } );

	for( size_t idx = 0; idx < unitvecs.size(); idx++ )
	{
		UnitVectorAtVar(idx) = unitvecs[sortedangles[idx].second];
		m_angles.AtVar(idx)  = sortedangles[idx].first;
	}

	ok = ok && CheckForUniformAngles();
	return ok;
}

bool SKTRAN_UnitSphere_Plane::ConstructPlane( const std::vector<double>& angles, const nxVector& reference, const nxVector& normal )
{
	bool ok = true;

	m_reference = reference;
	m_normal = normal;
	m_x = m_reference;
	m_y = m_reference.Cross(m_normal);
	std::vector<double> sortedangles = angles;
	std::sort( std::begin( sortedangles ), std::end( sortedangles ) );

	m_angles.AllocateGridArray( sortedangles.size() );
	for( size_t idx = 0; idx < sortedangles.size(); idx++ )
	{
		m_angles.AtVar( idx ) = sortedangles[idx];
		// calculate the unit vector of the point
		UnitVectorAtVar(idx) = UnitVectorFromAngle( m_angles.AtVar( idx ) );
	}

	ok = ok && CheckForUniformAngles();

	return ok;
}

nxVector SKTRAN_UnitSphere_Plane::UnitVectorFromAngle( double th ) const
{
	nxVector ret;
	double sinth = nxmath::sind(th);
	double costh = nxmath::cosd(th);
	double ux = m_normal.X();
	double uy = m_normal.Y();
	double uz = m_normal.Z();
	double rx = m_reference.X();
	double ry = m_reference.Y();
	double rz = m_reference.Z();

	ret.SetCoords( (costh + ux*ux*(1-costh))*rx + (ux*uy*(1-costh) - uz*sinth)*ry + (ux*uz*(1-costh) + uy*sinth)*rz,
				   (uy*ux*(1-costh) + uz*sinth)*rx + (costh + uy*uy*(1-costh))*ry + (uy*uz*(1-costh) + ux*sinth)*rz,
				   (uz*ux*(1-costh) - uy*sinth)*rx + (uz*uy*(1-costh) + ux*sinth)*ry + (costh + uz*uz*(1-costh))*rz );
	return ret.UnitVector();
}

bool SKTRAN_UnitSphere_Plane::CheckForUniformAngles()
{
	if( m_angles.NumAngles() >= 2 )
	{
		double anglediff = m_angles.At(1) - m_angles.At(0);
		bool isuniform = true;
		for( size_t idx = 0; idx < m_angles.NumAngles()-1; idx++ )
		{
			isuniform = isuniform && (abs(m_angles.At(idx+1) - m_angles.At(idx) - anglediff ) < 1E-10);
			anglediff = m_angles.At(idx+1) - m_angles.At(idx);
		}
		if( isuniform )
		{
			//m_angles.SetGridSearchMode( SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM );
		}
	}
	return true;
}

double SKTRAN_UnitSphere_Plane::ProjectedAngle( const nxVector& unit ) const
{
	nxVector projvector = unit - unit.Dot( m_normal ) * unit;
	return nxmath::atan2d( m_y & projvector, m_x & projvector );
}

bool SKTRAN_UnitSphere_Plane::Triangulate	( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices) const
{
	double angle = ProjectedAngle( unit );

	bool ok = m_angles.FindBoundingIndices( angle, SKTRAN_GridDefBase_V2::OUTOFBOUND_TRUNCATE, &unit_indexptr[0], &unit_weightptr[0], &unit_indexptr[1], &unit_weightptr[1]);
	unit_indexptr[2] = 0;
	unit_weightptr[2] = 0;

	return ok;
}

bool SKTRAN_UnitSphere_Plane::Triangulate	( const nxVector& unit, size_t* unit_indexptr, double* unit_weightptr, size_t maxvertices, size_t& speedHelper) const
{
	return Triangulate( unit, unit_indexptr, unit_weightptr, maxvertices );
}