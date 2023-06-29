#include "include/sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::SKTRAN_HR_Specs_Internal_OpticalPropertiesTable		2013-06-14*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::SKTRAN_HR_Specs_Internal_OpticalPropertiesTable()
{
	m_optproptype        = SKTRAN_HR_OpticalPropertiesTableType::SKTRAN_HR_OpticalPropertiesTableType_1d;
	m_maxPolarizationOrder = 0;
	m_polHigherOrderBehaviour = SKTRAN_HR_PolHOType::unpolarized;
	m_scatgridresolution = 0.0;
	m_heightspacing      = 0.0;
	m_numprofiles        = 0;
	m_coneanglesep       = 0.0;
	m_numcones           = 0;
	m_profilepercone     = 0;
	m_forcecacheupdates  = false;
	m_raymanager         = nullptr;
	m_atmosphereHasDelta = SKTRAN_HR_AtmosphereHasDelta::unset;
	m_norefractionsinglescatter = false;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::~SKTRAN_HR_Specs_Internal_OpticalPropertiesTable		2013-06-14*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::~SKTRAN_HR_Specs_Internal_OpticalPropertiesTable()
{
	
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::Configure		2013-06-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::Configure( const SKTRAN_HR_Specs_User& specs, const SKTRAN_HR_RayTracingRegionManager* raymanager )
{
	bool ok = true;

	const SKTRAN_HR_Specs_User_OpticalPropertiesTable& optspecs = specs.OpticalPropertiesSpecsConst();

	m_numprofiles    = optspecs.m_numprofiles;
	m_optproptype    = optspecs.m_optproptype;
	m_maxPolarizationOrder = optspecs.m_maxPolarizationOrder;
	m_polHigherOrderBehaviour = optspecs.m_polHigherOrderBehaviour;
	m_polarizationHigherOrderFraction = optspecs.m_polarizationHigherOrderFraction;

	m_heightspacing = optspecs.m_heightres;
	m_scatgridresolution = optspecs.m_scatres;
	m_raymanager = raymanager;
	m_coneanglesep = optspecs.m_coneanglesep;
	m_profilepercone = optspecs.m_profilepercone;
	m_numcones = optspecs.m_numcones;
	m_precachewavel = optspecs.m_precachewavel;
	m_atmosphereHasDelta = optspecs.m_atmosphereHasDelta;
	m_forcecacheupdates = optspecs.m_forcecacheupdates;


	if( optspecs.m_normal.IsValid() && optspecs.m_reference.IsValid() )
	{
		// user specified normal and reference
		m_normal = optspecs.m_normal.UnitVector();  // edit luf542 2022-01-28
		m_reference = optspecs.m_reference.UnitVector();  // edit luf542 2022-01-28
	}
	else
	{
		// Calculate normal and reference from the LOS
		MakeNormalAndReferenceFromLOS();
	}

	if( optspecs.m_anglegrid.size() != 0 )
	{
		// user specified angle grid
		m_anglegrid = optspecs.m_anglegrid;
	}
	else
	{
		// create default angle grid
		MakeDefaultAngleGrid();
	}
	m_heightgrid = optspecs.m_heightgrid;
	//ok = ok && ConfigureDefaults();

	return ok;
}

SKTRAN_GridDefSLON_V21 SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::MakeSLONGrid( const SKTRAN_CoordinateTransform_V2& coords ) const
{
	SKTRAN_GridDefSLON_V21 ret;

	// Note that if we are in SZA table mode, In3dMode() returns false
	if( !In3dMode() )
	{
		// just need one slon at the TP
		ret.AllocateGridArray( 1 );
		ret.AtVar( 0 ) = 0.0;
		return ret;
	}

	if( m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_3D_UnitSphere )
	{
		// slon is going to vary in each cone
		ret.AllocateGridArray( m_numcones*2 + 1 );
		ret.AtVar( 0 ) = m_numcones * m_coneanglesep * (-1) * nxmath::Pi / 180;
		for( size_t idx = 1; idx < m_numcones*2 + 1; idx++ )
		{
			ret.AtVar( idx ) = ret.AtVar( idx-1 ) + m_coneanglesep * nxmath::Pi / 180;
		}
		return ret;
	}

	if( m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_LOSPlane || m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_SZA )
	{
		// TODO: calculate actual slon grid
		// worst case is when the angle grid is the slon direction
		ret.AllocateGridArray( m_anglegrid.size() );
		for( size_t idx = 0; idx < m_anglegrid.size(); idx++ )
		{
			ret.AtVar( idx ) = m_anglegrid[idx] * nxmath::Pi / 180;
		}
		return ret;
	}

	nxLog::Record( NXLOG_WARNING, "SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::MakeSLONGRID SHOUDLNT BE HERE");
	return ret;

}

void SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::MakeNormalAndReferenceFromLOS()
{
	double lat,lon;
	m_raymanager->GetReferencePoint( &lat, &lon );

	nxGeodetic geo = nxGeodetic(lat, lon);
	m_reference = geo.Location();
	m_reference = m_reference.UnitVector();

	if( SKTRAN_HR_OpticalPropertiesTableType_SZA == m_optproptype )
	{
		// normal vector is reference cross sun
		nxVector sun(0,0,1);
		m_normal = m_reference.Cross( sun ).UnitVector();
		return;
	}
	/* else */


	nxVector in,out;
	m_raymanager->GetBoundingReferences(in, out);

	m_normal = in.Cross(m_reference);
	m_normal = -1.0*m_normal.UnitVector();
}

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::AddOpticalInformationToRayTracer( SKTRAN_RayTracer_Straight_Generic& raytracer, const SKTRAN_CoordinateTransform_V2* coords ) const
{
	if( m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_1d )
	{
		return true;
	}
	else if ( m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_LOSPlane || m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_SZA )
	{
		return AddPlaneInformationToRayTracer( raytracer, *coords );
	}
	else if ( m_optproptype == SKTRAN_HR_OpticalPropertiesTableType_3D_UnitSphere )
	{
		return AddDelaunayInformationToRayTracer( raytracer, *coords );
	}
	return false;
}

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::AddPlaneInformationToRayTracer( SKTRAN_RayTracer_Straight_Generic& raytracer, const SKTRAN_CoordinateTransform_V2& coords ) const
{
	nxVector x,y;
	HELIODETIC_VECTOR xh,yh;
	xh = coords.GeographicToHelio( m_reference );
	yh = coords.GeographicToHelio( m_reference.Cross(m_normal));
	x.SetCoords( xh.X(), xh.Y(), xh.Z() );
	y.SetCoords( yh.X(), yh.Y(), yh.Z() );
	for( size_t idx = 0; idx < m_anglegrid.size(); idx++ )
	{
		nxVector normal;
		normal = nxmath::cosd( m_anglegrid[idx] + 90.0 ) * x + nxmath::sind( m_anglegrid[idx] + 90.0 ) * y;
		raytracer.AddGeometryObject( std::unique_ptr<SKTRAN_GeometryObject_Plane> ( new SKTRAN_GeometryObject_Plane(normal)) );
	}
	return true;
}

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::AddDelaunayInformationToRayTracer( SKTRAN_RayTracer_Straight_Generic& raytracer, const SKTRAN_CoordinateTransform_V2& coords ) const
{
	// First add the cones
	// series of cones centered on the reference point, spaced evenly with m_coneanglesep 
	HELIODETIC_VECTOR refhelio = coords.ReferencePoint(0.0).Vector();
	nxVector coneunitvec = nxVector(refhelio.X(), refhelio.Y(), refhelio.Z()).UnitVector();

	for( size_t idx = 1; idx < m_numcones; idx++ )
	{
		raytracer.AddGeometryObject( std::unique_ptr<SKTRAN_GeometryObject_Cone> ( new SKTRAN_GeometryObject_Cone( coneunitvec, idx*m_coneanglesep ) ) );
	}

	// Now add the planes
	// Define each plane by two vectors, the reference, and a point on a circle 10 degrees off
	// since we stagger the spacing every circle we need twice as many planes
	double zenith = 10 * nxmath::Pi/180;
	double azimuth;
	nxVector secondvec;
	nxVector normal;
	for( size_t idx = 0; idx < m_profilepercone*2; idx++ )
	{
		azimuth = 2*nxmath::Pi * (idx) / (m_profilepercone*2);
		secondvec = CalcRotatedVector( coneunitvec, zenith, azimuth );
		normal = coneunitvec.Cross( secondvec ).UnitVector();
		raytracer.AddGeometryObject( std::unique_ptr<SKTRAN_GeometryObject_Plane> ( new SKTRAN_GeometryObject_Plane( normal ) ) );
	}
	return true;
}

void SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::MakeDefaultAngleGrid()
{
	// default is 15 degrees on either side of the tangent point
	double spacing = 30.0 / double(m_numprofiles-1);
	m_anglegrid.resize(m_numprofiles);
	for( size_t idx = 0; idx < m_numprofiles; idx++ )
	{
		m_anglegrid[idx] = -15.0 + spacing*idx;
	}
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::ConfigureDefaults		2013-06-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::ConfigureDefaults()
{
	bool ok = true;
	
	m_optproptype = SKTRAN_HR_OpticalPropertiesTableType_3D_UnitSphere;
	m_maxPolarizationOrder = 0;
	m_polHigherOrderBehaviour = SKTRAN_HR_PolHOType::unpolarized;
	m_atmosphereHasDelta = SKTRAN_HR_AtmosphereHasDelta::no;
	//m_optproptype = SKTRAN_HR_OpticalPropertiesTableType_1d;
	m_scatgridresolution = 0.5;
	m_heightspacing = 500;

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::Create1dTable		2013-06-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::Create1dTable( OpticalTablePtr& table, const SKTRAN_CoordinateTransform_V2& coords, double toaHeight )
{
	bool ok = true;
	std::unique_ptr<SKTRAN_TableOpticalProperties_3D_UnitSphere> opttable;
	if (SKTRAN_HR_OpticalPropertiesTableType_1d == m_optproptype)
	{
		opttable.reset(new SKTRAN_TableOpticalProperties_3D_UnitSphere);
	}
	else if (SKTRAN_HR_OpticalPropertiesTableType_1d_ConstantLayers == m_optproptype)
	{
		opttable.reset(new SKTRAN_TableOpticalProperties_3D_UnitSphere_Constant);
	}
	else
	{
		ok = false;
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::Create1dTable, invalid table type.");
	}
	SKTRAN_GridDefScatterAngle_V21* scattergrid = new SKTRAN_GridDefScatterAngle_V21;
	SKTRAN_GridDefOpticalPropertiesRadii_V21* heightgrid = new SKTRAN_GridDefOpticalPropertiesRadii_V21;
	std::unique_ptr< SKTRAN_PolarizationProperties_Base> polobj(nullptr);

	nxVector referenceHELIO;
	referenceHELIO.SetCoords( coords.ReferencePoint(0.0).UnitVector().X(),coords.ReferencePoint(0.0).UnitVector().Y(),coords.ReferencePoint(0.0).UnitVector().Z() );

	SKTRAN_UnitSphere_V2* unitsphere = new SKTRAN_UnitSphere_Dummy( referenceHELIO );

	ok = ok && MakeScatterAngleGrid( *scattergrid );
	ok = ok && MakeHeightGrid( *heightgrid, toaHeight );

	ok = ok && CreatePolarizationObject( polobj );
	opttable->SetPolarizationProperties( polobj );
	opttable->SetAltitudes( *heightgrid );
	opttable->SetScatterGrid( *scattergrid );
	opttable->SetUnitSphere( *unitsphere );
	opttable->SetWavelengths( m_precachewavel );
	opttable->ConfigureGeometry( nullptr );
	table = std::move( opttable );
	table->AddRef();

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::MakeScatterAngleGrid		2013-06-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::MakeScatterAngleGrid( SKTRAN_GridDefScatterAngle_V21& scattergrid ) const
{
	bool ok = true;

	if(!m_norefractionsinglescatter)
	{
		ok = ok && scattergrid.Configure( m_scatgridresolution );
		scattergrid.SetGridSearchMode( SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM );
	}
	else
	{
		double minssa, maxssa;
		m_raymanager->GetBoundingLOSScatteringAngles(minssa, maxssa);

		double mintablessa = std::max(0.0, minssa - m_scatgridresolution);
		double maxtablessa = std::min(180.0, maxssa + m_scatgridresolution);

		ok = ok && scattergrid.Configure(m_scatgridresolution, mintablessa, maxtablessa);
	}
	
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::MakeUniformHeightGrid		2013-06-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::MakeHeightGrid( SKTRAN_GridDefOpticalPropertiesRadii_V21& heightgrid, double toaHeight )
{
	bool ok = true;
	
	if( m_heightgrid.size() != 0 )
	{
		return heightgrid.ConfigureAltitudes( &m_heightgrid[0], m_heightgrid.size() );
	}

	size_t numalts = static_cast<size_t>(ceil( (toaHeight) / m_heightspacing)) + 1;

	std::vector<double> altgrid;
	altgrid.resize(numalts);
	for( size_t i = 0; i < numalts; i++ )
	{
		altgrid[i] = i*m_heightspacing;
	}
	ok = ok && heightgrid.ConfigureAltitudes( &altgrid[0], numalts );
	heightgrid.SetGridSearchMode( SKTRAN_GridDefBase_V2::GRIDSEARCH_UNIFORM );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::MakeUnitSphere		2013-09-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::MakeUnitSphere( SKTRAN_UnitSphere_V2** unitsphere )
{
	bool ok = true;

	SKTRAN_UnitSphereME* spherelocal = new SKTRAN_UnitSphereME( m_numprofiles );
	*unitsphere = spherelocal;

	return ok;
}



bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::MakeDelaunaySphere( SKTRAN_UnitSphere_V2** unitsphere, const SKTRAN_CoordinateTransform_V2& coords )
{
	size_t numcircles = m_numcones;				// number of circles around the tangent point
	size_t pointspercircle = m_profilepercone;
	double circlesep = m_coneanglesep;
	double zenith;						// the angle off of the tangent
	double azimuth;						// angle around the cone
	HELIODETIC_POINT temp;
	HELIODETIC_VECTOR tempv;
	nxVector			tangentpoint;

	size_t numpoints = numcircles * pointspercircle + 1;	// extra 1 is for tangent point

	std::vector< nxVector > unitvecs;
	unitvecs.resize( numpoints );

	bool hascreatedsphere = false;
	SKTRAN_UnitSphere_Delaunay_nonTabledLookup* spherelocal;

	temp = coords.ReferencePoint( 0.0 );
	tempv.SetCoords( temp.Vector().X(), 0, temp.Vector().Z() );
	tangentpoint.SetCoords( tempv.UnitVector().X(),
						   tempv.UnitVector().Y(),
						   tempv.UnitVector().Z() ); // tangent point;
	unitvecs[0] = tangentpoint;

	// The delaunay triangulation code fails frequently, if it does fail nudge the tangent
	// point very slightly off and retry
	size_t numnudges = 0;
	while ( !hascreatedsphere && numnudges < 100 )
	{
		for( size_t circleidx = 0; circleidx < numcircles; circleidx++ )
		{
			for( size_t pointidx = 0; pointidx < pointspercircle; pointidx++ )
			{ 
				zenith = (circleidx+1) * circlesep;
				azimuth = 2*nxmath::Pi * (pointidx + (circleidx % 2)*0.5) / (pointspercircle);
				unitvecs[pointidx + circleidx * pointspercircle + 1] = CalcRotatedVector( unitvecs[0], zenith, azimuth );
			}
		}

		nxVector tangentopposite;
		tangentopposite.SetCoords( -1.0*unitvecs[0].X(),
								   -1.0*unitvecs[0].Y(),
								   -1.0*unitvecs[0].Z() );
		spherelocal = new SKTRAN_UnitSphere_Delaunay_nonTabledLookup;
		hascreatedsphere = spherelocal->CreateTriangulation( &unitvecs[0], numpoints, &tangentopposite );
		if( !hascreatedsphere )
		{
			delete spherelocal;
			// nudge the reference point off slightly
			unitvecs[0] = CalcRotatedVector( tangentpoint, 0.001, numnudges*2*nxmath::Pi/100 );
			++numnudges;
		}
	}
	if( hascreatedsphere )
	{
		*unitsphere = spherelocal;
	}
	else
	{
		nxLog::Record( NXLOG_WARNING, "Delaunay triangulation failed even after 100 nudges" );
	}


	return true;

}

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::MakeLOSPlaneSphere( SKTRAN_UnitSphere_V2** sphere, const SKTRAN_CoordinateTransform_V2& coords )
{
	bool ok = true;
	nxVector look;
	std::vector<nxVector> locations;

	SKTRAN_UnitSphere_Plane* spherelocal = new SKTRAN_UnitSphere_Plane;

	nxVector x,y;
	HELIODETIC_VECTOR xh, yh;
	xh = coords.GeographicToHelio( m_reference );
	yh = coords.GeographicToHelio( m_reference.Cross(m_normal));
	x.SetCoords( xh.X(), xh.Y(), xh.Z() );
	y.SetCoords( yh.X(), yh.Y(), yh.Z() );

	for( const auto& th : m_anglegrid )
	{
		nxVector loc = nxmath::cosd(th) * x + nxmath::sind(th) * y;
		locations.emplace_back( loc );
	}

	spherelocal->ConstructPlane( locations );
	*sphere = spherelocal;

	return ok;
}

nxVector SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::CalcRotatedVector( const nxVector& tangent, double zenith, double azimuth ) const
{
	// note that tangent_y is 0 since we force the tangent to be at slon 0
	// so rotate the zenith in 2 dimensions first

	nxVector result;

	result.SetCoords( nxmath::cosd(zenith) * tangent.X() - nxmath::sind(zenith) * tangent.Z(),
					  0,
					  nxmath::sind(zenith) * tangent.X() + nxmath::cosd(zenith) * tangent.Z() );

	// now we need to rotate about the tangent point, but still remembering
	// tangent_y = result_y = 0
	double cosazi = cos(azimuth);
	double sinazi = sin(azimuth);
	double ux = tangent.X();
	double uz = tangent.Z();
	double rx = result.X();
	double rz = result.Z();

	result.SetCoords( ( cosazi + ux*ux*(1-cosazi) ) * rx + (ux*uz*(1-cosazi)) * rz,
		              uz*sinazi*rx + -1.0*ux*sinazi*rz,
					  uz*ux*(1-cosazi)*rx + (cosazi + uz*uz*(1-cosazi))*rz );

	return result;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::CreateOpticalTable		2013-06-14*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::CreateOpticalTable( OpticalTablePtr& table, const SKTRAN_CoordinateTransform_V2& coords, double toaHeight )
{
	bool ok = true;

	if( SKTRAN_HR_OpticalPropertiesTableType_1d == m_optproptype ||
		SKTRAN_HR_OpticalPropertiesTableType_1d_ConstantLayers == m_optproptype )
	{
		ok = ok && Create1dTable( table, coords, toaHeight );
	}
	else if( SKTRAN_HR_OpticalPropertiesTableType_3D_UnitSphere == m_optproptype || 
		SKTRAN_HR_OpticalPropertiesTableType_LOSPlane == m_optproptype ||
		SKTRAN_HR_OpticalPropertiesTableType_SZA == m_optproptype )
	{
		ok = ok && Created3DUnitSphereTable( table, coords, toaHeight );
	}
	else
	{
		ok = false;
		nxLog::Record(NXLOG_WARNING, "Error Should not be here, SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::CreateOpticalTable");
	}
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::CreateOpticalTable, Could not create optical properties table.");

	return ok;
}



/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::Created3DUnitSphereTable		2013-09-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::CreatePolarizationObject( std::unique_ptr< SKTRAN_PolarizationProperties_Base >& polobject ) const
{
	bool ok = true;
	
	if(ok){
		if(0==m_maxPolarizationOrder){
			polobject = std::unique_ptr< SKTRAN_PolarizationProperties_Base > (new SKTRAN_PolarizationProperties_NoPolarization);
		} else{
			switch ( m_atmosphereHasDelta ){
			case SKTRAN_HR_AtmosphereHasDelta::yes:
				polobject = std::unique_ptr< SKTRAN_PolarizationProperties_Base > (new SKTRAN_PolarizationProperties_Polarized_Eddington);
				break;
			case SKTRAN_HR_AtmosphereHasDelta::no:
				polobject = std::unique_ptr< SKTRAN_PolarizationProperties_Base > (new SKTRAN_PolarizationProperties_Polarized);
				break;
			case SKTRAN_HR_AtmosphereHasDelta::unset:
				polobject = std::unique_ptr< SKTRAN_PolarizationProperties_Base > (nullptr);
				nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::CreatePolarizationObject, Simulation is polarized, but atmosphere was not checked for delta.");
				ok = false;
				break;
			}
		}
	}

	ok = ok && nullptr!=polobject;
	if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::CreatePolarizationObject, Could not create polarization object.");
	
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::Created3DUnitSphereTable		2013-09-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_OpticalPropertiesTable::Created3DUnitSphereTable( OpticalTablePtr& table, const SKTRAN_CoordinateTransform_V2& coords, double toaHeight )
{
	bool ok = true;

	std::unique_ptr<SKTRAN_TableOpticalProperties_3D_UnitSphere> opttable ( new SKTRAN_TableOpticalProperties_3D_UnitSphere );
	SKTRAN_GridDefScatterAngle_V21* scattergrid = new SKTRAN_GridDefScatterAngle_V21;
	SKTRAN_GridDefOpticalPropertiesRadii_V21* heightgrid = new SKTRAN_GridDefOpticalPropertiesRadii_V21;
	std::unique_ptr< SKTRAN_PolarizationProperties_Base > polobj;
	SKTRAN_UnitSphere_V2* unitsphere;

	ok = ok && MakeScatterAngleGrid( *scattergrid );
	ok = ok && MakeHeightGrid( *heightgrid, toaHeight );
	//ok = ok && MakeUnitSphere( &unitsphere );
	if( SKTRAN_HR_OpticalPropertiesTableType_LOSPlane == m_optproptype || SKTRAN_HR_OpticalPropertiesTableType_SZA == m_optproptype )
	{
		ok = ok && MakeLOSPlaneSphere(&unitsphere, coords );
	}
	else if ( SKTRAN_HR_OpticalPropertiesTableType_3D_UnitSphere == m_optproptype )
	{
		ok = ok && MakeDelaunaySphere(&unitsphere, coords );
	}
	else
	{
		nxLog::Record(NXLOG_WARNING, "Invalid type in Created3DUnitSphereTable");
	}

	opttable->SetAltitudes( *heightgrid );
	opttable->SetScatterGrid( *scattergrid );
	opttable->SetUnitSphere( *unitsphere );
	opttable->SetWavelengths( m_precachewavel );
	opttable->SetForceCacheUpdates( m_forcecacheupdates );

	ok = ok && CreatePolarizationObject( polobj );
	opttable->SetPolarizationProperties( polobj );
	opttable->ConfigureGeometry( NULL );

	table = std::move( opttable );
	table->AddRef();

	return ok;
}
