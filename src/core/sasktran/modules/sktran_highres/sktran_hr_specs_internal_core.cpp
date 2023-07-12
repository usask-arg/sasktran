#include "include/sktran_hr_internals.h"


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Core::SKTRAN_HR_Specs_Internal_Core		2013-06-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Specs_Internal_Core::SKTRAN_HR_Specs_Internal_Core()
{
	m_numprofiles  = 0;
	m_numofflook   = 0;
	m_numinterp    = 0;
	m_angleofflook = 0.0;
	m_diffuseheightres = 0.0;
	m_diffusemaxheight = std::numeric_limits<double>::quiet_NaN();
	m_linesofsight = nullptr;
	m_in3dmode = false;
	m_calcwf   = false;
	m_track = false;
	m_userdiffuseplacementtype = SKTRAN_HR_DiffuseProfilePlacementType::SKTRAN_HR_DiffuseProfilePlacement_LinearLOS;
	m_internaldiffuseplacementtype = InternalDiffusePlacementType::Undefined;
	m_maxPolarizationOrder = 0;
	m_polHigherOrderBehaviour = SKTRAN_HR_PolHOType::unpolarized;
	//	SKTRAN_HR_Specs_Internal_wf									m_wfspecs;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Core::~SKTRAN_HR_Specs_Internal_Core		2013-06-17*/
/** **/
/*---------------------------------------------------------------------------*/

SKTRAN_HR_Specs_Internal_Core::~SKTRAN_HR_Specs_Internal_Core()
{
	ReleaseResources();
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Core::Configure		2013-06-17*/
/** Configure is called when we decide to reconcile all of the various user
 *	settings.  It is at this time that we can decide if the user has set sensible
 *	settings or not.
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_Core::Configure( const SKTRAN_Specifications_Base& specs, const SKTRAN_LineOfSightArray_V21& linesofsight )
{
	bool ok = true;
	const SKTRAN_HR_Specs_User		userspecs = dynamic_cast<const SKTRAN_HR_Specs_User&>(specs);
	m_linesofsight = &linesofsight;

	ok = ok && m_raytracerspecs.Configure( userspecs, (this) );
	ok = ok && m_integratorspecs.Configure( userspecs );
	ok = ok && m_opttablespecs.Configure( userspecs, &m_raymanager );
	ok = ok && m_diffusespecs.Configure( userspecs );
	ok = ok && m_wfspecs.Configure( userspecs );

	m_numprofiles          = userspecs.DiffuseSpecsConst().GetNumProfiles();
	m_numofflook           = userspecs.DiffuseSpecsConst().GetNumOffLook();
	m_angleofflook         = userspecs.DiffuseSpecsConst().GetAngleOffLook();
	m_diffuseheightres     = userspecs.DiffuseSpecsConst().GetHeightRes();
	m_diffusemaxheight     = userspecs.DiffuseSpecsConst().GetMaxDiffuseHeight();			// Note this value is usually Nan
	m_manualdiffuseheights = userspecs.DiffuseSpecsConst().GetManualDiffuseHeights();
	m_userdiffuseplacementtype     = userspecs.DiffuseSpecsConst().m_diffuseplacement;

	m_scatterorder = userspecs.ScatterOrder();

	m_maxPolarizationOrder    = userspecs.OpticalPropertiesSpecsConst().GetMaxPolarizationOrder();
	m_polHigherOrderBehaviour = userspecs.OpticalPropertiesSpecsConst().GetPolarizationHigherOrderBehaviour();
	m_polarizationHigherOrderFraction = userspecs.OpticalPropertiesSpecsConst().GetPolarizationHigherOrderFraction();

	m_diffuseprofilelocations = userspecs.DiffuseSpecsConst().m_manualdiffuselocations;
	m_diffuseprofilelatlons = userspecs.DiffuseSpecsConst().m_manualdiffuselatlons;
	m_diffuseprofileszas = userspecs.DiffuseSpecsConst().m_manualdiffuseszas;
	m_diffuseprofilelospositions = userspecs.DiffuseSpecsConst().m_manualdiffuselospositions;
	m_diffuseplanereference = userspecs.DiffuseSpecsConst().m_diffuseplanereference.UnitVector();
	m_diffuseplanenormal = userspecs.DiffuseSpecsConst().m_diffuseplanenormal.UnitVector();
	m_diffuseplaneangles = userspecs.DiffuseSpecsConst().m_diffuseplaneangles;

	m_trackedscatords = userspecs.DiffuseSpecsConst().m_diagnosticscatorders;
	m_trackeddiffprofs = userspecs.DiffuseSpecsConst().m_diagnosticdiffprofs;

	if (m_trackeddiffprofs.size() > 0 || m_trackedscatords.size() > 0)
	{
		m_track = true;
		if (m_trackeddiffprofs.size() == 0)
		{
			m_trackeddiffprofs.resize(1);
			m_trackeddiffprofs[0] = 0;
		}
		if (m_trackeddiffprofs.size() > m_numprofiles)
		{
			m_trackeddiffprofs.resize(m_numprofiles);
			for (size_t count = 0; count < m_numprofiles; count++)
			{
				m_trackeddiffprofs[count] = count;
			}
		}
		if (m_trackedscatords.size() == 0)
		{
			m_trackedscatords.resize(1);
			m_trackedscatords[0] = 1;
		}
	}

	m_in3dmode = m_opttablespecs.In3dMode();
	m_numinterp = 2;

	if (m_scatterorder == 1 && !m_raytracerspecs.LineOfSightRefractionEnabled())
	{
		m_opttablespecs.SetNoRefractionSingleScatter(true);
	}

	return ok;
}

const char* SKTRAN_HR_Specs_Internal_Core::DiffuseLocationTypeChar(InternalDiffusePlacementType type) const
{
	switch (type)
	{
		case InternalDiffusePlacementType::ManualLocation:
			return "manual locations";
		case InternalDiffusePlacementType::ManualLatLon:
			return "manual lat/lon pairs";
		case InternalDiffusePlacementType::ManualSZA:
			return "manual SZAs";
		case InternalDiffusePlacementType::ManualLOS:
			return "manual LOS positions";
		case InternalDiffusePlacementType::ManualPlane:
			return "manual plane angles";
		case InternalDiffusePlacementType::LinearSZA:
			return "linear SZAs";
		case InternalDiffusePlacementType::LinearSZAForceTP:
			return "linear SZAs with forced tangent point";
		case InternalDiffusePlacementType::SmartSZA:
			return "smart SZAs";
		default:
			return "unknown internal diffuse location type";
	}
}

bool SKTRAN_HR_Specs_Internal_Core::DiffusePlacementWarning(InternalDiffusePlacementType selected, InternalDiffusePlacementType overridden) const
{
	nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_Internal_Core::CreateDiffuseLocations, Multiple diffuse profile placement options were selected. The %s have been overriden by %s.", DiffuseLocationTypeChar(overridden), DiffuseLocationTypeChar(selected));
	return true;
}

bool SKTRAN_HR_Specs_Internal_Core::CreateDiffuseTableType( std::unique_ptr<SKTRAN_HR_Diffuse_Table_CPU>& diffusetable, bool usecache )
{
	bool ok = true;
	SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE primaryscantype;

	if( usecache )
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_Internal_Core::Cached diffuse table is not currently supported");
		return false;
	}

	m_internaldiffuseplacementtype = InternalDiffusePlacementType::Undefined;
	if (m_diffuseprofilelocations.size() > 0)
	{
		m_internaldiffuseplacementtype = InternalDiffusePlacementType::ManualLocation;
	}
	if (m_diffuseprofilelatlons.size() > 0)
	{
		if (m_internaldiffuseplacementtype == InternalDiffusePlacementType::Undefined) m_internaldiffuseplacementtype = InternalDiffusePlacementType::ManualLatLon;
		else DiffusePlacementWarning(m_internaldiffuseplacementtype, InternalDiffusePlacementType::ManualLatLon);
	}
	if (m_diffuseprofileszas.size() > 0)
	{
		if (m_internaldiffuseplacementtype == InternalDiffusePlacementType::Undefined) m_internaldiffuseplacementtype = InternalDiffusePlacementType::ManualSZA;
		else DiffusePlacementWarning(m_internaldiffuseplacementtype, InternalDiffusePlacementType::ManualSZA);
	}
	if (m_diffuseprofilelospositions.size() > 0)
	{
		if (m_internaldiffuseplacementtype == InternalDiffusePlacementType::Undefined) m_internaldiffuseplacementtype = InternalDiffusePlacementType::ManualLOS;
		else DiffusePlacementWarning(m_internaldiffuseplacementtype, InternalDiffusePlacementType::ManualLOS);
	}
	if (m_diffuseplanenormal.IsValid() && !m_diffuseplanenormal.IsZero() && m_diffuseplanereference.IsValid() && !m_diffuseplanereference.IsZero() && m_diffuseplaneangles.size() > 0)
	{
		if (m_internaldiffuseplacementtype == InternalDiffusePlacementType::Undefined) m_internaldiffuseplacementtype = InternalDiffusePlacementType::ManualPlane;
		else DiffusePlacementWarning(m_internaldiffuseplacementtype, InternalDiffusePlacementType::ManualPlane);
	}
	if (m_userdiffuseplacementtype == SKTRAN_HR_DiffuseProfilePlacement_LinearSZA)
	{
		if (m_internaldiffuseplacementtype == InternalDiffusePlacementType::Undefined) m_internaldiffuseplacementtype = InternalDiffusePlacementType::LinearSZA;
		else DiffusePlacementWarning(m_internaldiffuseplacementtype, InternalDiffusePlacementType::LinearSZA); 
	}
	if (m_userdiffuseplacementtype == SKTRAN_HR_DiffuseProfilePlacement_LinearSZAForceTP)
	{
		if (m_internaldiffuseplacementtype == InternalDiffusePlacementType::Undefined) m_internaldiffuseplacementtype = InternalDiffusePlacementType::LinearSZAForceTP;
		else DiffusePlacementWarning(m_internaldiffuseplacementtype, InternalDiffusePlacementType::LinearSZAForceTP);
	}
	if (m_userdiffuseplacementtype == SKTRAN_HR_DiffuseProfilePlacement_SmartSZA)
	{
		if (m_internaldiffuseplacementtype == InternalDiffusePlacementType::Undefined) m_internaldiffuseplacementtype = InternalDiffusePlacementType::SmartSZA; 
		else DiffusePlacementWarning(m_internaldiffuseplacementtype, InternalDiffusePlacementType::SmartSZA);
	}
	if (m_userdiffuseplacementtype == SKTRAN_HR_DiffuseProfilePlacement_LinearLOS)
	{
		if (m_internaldiffuseplacementtype == InternalDiffusePlacementType::Undefined)
		{
			primaryscantype = SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE::SKTRAN_VIEWING_TYPE_UNDEFINED;
			ok = ok && m_raymanager.GetPrimaryScanType(&primaryscantype);

			switch (primaryscantype)
			{
				case SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATLIMB:
				case SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE::SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATLIMB:
				{ // omit diffuse profiles at entry and exit points, tangent point replaced with reference point (usually the same)
					if (m_numofflook % 2 == 0)
						nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_Internal_Core::CreateDiffuseLocations, Number of off-look diffuse planes must be odd; %d will be reduced to %d", m_numofflook, m_numofflook - 1);
					if (m_numofflook > 2)
						m_internaldiffuseplacementtype = InternalDiffusePlacementType::LinearLOS2D;  
					else
						m_internaldiffuseplacementtype = InternalDiffusePlacementType::LinearLOS1D;
					break;
				}
				case SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATNADIR:
				{ // omit entry point, include exit point, replace one point with reference point if it is nearby
					if (m_numofflook > 1)
						nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_Internal_Core::CreateDiffuseLocations, Multiple diffuse planes not supported for nadir-pointing lines of sight, one plane will be used.");
					m_internaldiffuseplacementtype = InternalDiffusePlacementType::LinearLOS1DExit;
					break;
				}
				case SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE::SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATNADIR:
				{ // include entry/exit points, replace one point with reference point if it is nearby
					if (m_numofflook > 1)
						nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_Internal_Core::CreateDiffuseLocations, Multiple diffuse planes not supported for nadir-pointing lines of sight, one plane will be used.");
					m_internaldiffuseplacementtype = InternalDiffusePlacementType::LinearLOS1DEntryExit; 
					break;
				}
				case SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE::SKTRAN_VIEWING_TYPE_INSPACE_LOOKINGATSPACE:
				case SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE::SKTRAN_VIEWING_TYPE_INATMOS_LOOKINGATSPACE:
				case SKTRAN_LineOfSightEntry_V2::VIEWING_TYPE::SKTRAN_VIEWING_TYPE_NEARGROUND:
				{ // include entry point, omit exit point, replace one point with reference point if it is nearby
					if (m_numofflook > 2)
						nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_Internal_Core::CreateDiffuseLocations, multiple diffuse planes not supported for space-pointing lines of sight, one plane will be used.");
					m_internaldiffuseplacementtype = InternalDiffusePlacementType::LinearLOS1DEntry;
					break;
				}
				default:
					ok = false;
			}
		}
		// no warning in this case because SKTRAN_HR_DiffuseProfilePlacement_LinearLOS is the default
	}

	switch (m_internaldiffuseplacementtype)
	{
		case InternalDiffusePlacementType::LinearSZA:
		case InternalDiffusePlacementType::LinearSZAForceTP:
		case InternalDiffusePlacementType::SmartSZA:
		case InternalDiffusePlacementType::LinearLOS1D:
		case InternalDiffusePlacementType::LinearLOS1DEntry:
		case InternalDiffusePlacementType::LinearLOS1DExit:
		case InternalDiffusePlacementType::LinearLOS1DEntryExit: 
		{
			std::unique_ptr<SKTRAN_HR_Diffuse_Table_SZA> diffuselocal(new SKTRAN_HR_Diffuse_Table_SZA);
			diffusetable = std::move(diffuselocal);
			break;
		}
		default:
		{
			std::unique_ptr<SKTRAN_HR_Diffuse_Table_CPU> diffuselocal(new SKTRAN_HR_Diffuse_Table_CPU);
			diffusetable = std::move(diffuselocal);
			break;
		}
	}
	
	ok = ok && m_internaldiffuseplacementtype != InternalDiffusePlacementType::Undefined;
	ok = ok && CreateDiffuseTableAValues(diffusetable);
    
	return ok;
}


bool SKTRAN_HR_Specs_Internal_Core::CreateDiffuseTableAValues ( std::unique_ptr< SKTRAN_HR_Diffuse_Table_CPU >& diffusetable )
{
	//diffusetable->SetPolarization( m_polOrder, m_polHigherOrderBehaviour );
	bool ok = true; 
	std::unique_ptr< RadStore_Base > radStorage; 
	std::unique_ptr< Avals_Base >    avals;

	if( 0==m_maxPolarizationOrder ){
		radStorage = std::unique_ptr< RadStore_Base > ( new RadStore_Scalar );
		avals      = std::unique_ptr< Avals_Base >    ( new Avals_ScalarStore );
	} else if( 1==m_maxPolarizationOrder ){
		radStorage = std::unique_ptr< RadStore_Base > ( new RadStore_PV1 );
		avals      = std::unique_ptr< Avals_Base >    ( new Avals_ScalarStore );
	} else if( 1<m_maxPolarizationOrder ){
		bool useConstHOPolarization = SKTRAN_HR_PolHOType::constOut==m_polHigherOrderBehaviour;
		 auto tempradstore= std::unique_ptr< RadStore_Polarized > ( new RadStore_Polarized );
		tempradstore->ConfigurePolarizedScattering(m_maxPolarizationOrder-1,useConstHOPolarization, m_polarizationHigherOrderFraction);
        radStorage = std::move( tempradstore );

        switch( DiffuseSpecs( ).GetDiffuseMatrixStorageMethod( ) ){
            case SKTRAN_HR_DiffuseMatrixStorageMethod::scalar: 
		        avals = std::unique_ptr< Avals_Base > ( new Avals_ScalarStore );
                break;
            case SKTRAN_HR_DiffuseMatrixStorageMethod::scatter_cache: 
		        avals = std::unique_ptr< Avals_Base > ( new Avals_MatrixStore< SKTRAN_ScatMat_MIMSNC > );
                break;
            case SKTRAN_HR_DiffuseMatrixStorageMethod::scatter_interpolate: 
                avals = std::unique_ptr< Avals_Base > ( new Avals_MatrixTable );
                break;
            case SKTRAN_HR_DiffuseMatrixStorageMethod::phase:
		        avals = std::unique_ptr< Avals_Base > ( new Avals_MatrixStore< SKTRAN_PhaseMat_MIMSNC > );
                break;
            default:
                nxLog::Record( NXLOG_ERROR, "SKTRAN_HR_Specs_Internal_Core::CreateDiffuseTableAValues, Diffuse matrix storage method not recognized.");
                ok = false;
                break;
        }
	} else{
		ok = false;
	}

	ok = ok && nullptr!=radStorage;
	ok = ok && nullptr!=avals;

	if(ok){
		diffusetable->ConfigureStorage( radStorage, avals );

        //diffusetable->m_radStorage = std::move( radStorage );
        //diffusetable->m_Avals      = std::move( avals );
        //diffusetable->m_Avals->SetOpticalTable( diffusetable->m_opticaltable );
	} else{
		nxLog::Record(NXLOG_ERROR, "SKTRAN_HR_Specs_Internal_Core::CreateDiffuseTableAValues, Could not create diffuse helper objects.");
	}

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Core::CreateCoordinates		2013-06-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_Core::CreateCoordinates( std::shared_ptr< const SKTRAN_CoordinateTransform_V2>*	coords,
													   nxVector&												sun,
													   const SKTRAN_LineOfSightArray_V21&						linesofsight, 
													   double													surfaceHeight,
                                                       double                                                   toaHeight,
													   bool														nadirReferencePointOnGround )
{
	bool ok = true;

	if( sun.IsValid() )
	{
		m_raymanager.SetSun( sun );
	}
	ok = ok && m_raymanager.SetNadirReferencePointOnGround(nadirReferencePointOnGround);
	ok = ok && m_raymanager.SetGroundAltitude( surfaceHeight );
    ok = ok && m_raymanager.SetUpperBoundAltitude( toaHeight );
    ok = ok && m_raymanager.UpdateUndefinedParametersFromLinesOfSight( linesofsight );
	//nxLog::Record(NXLOG_INFO,"SKTRAN_HR_Specs_Internal_Core::CreateCoordinates ***** TODO **** Set the min and max height in the coordinates object. Currently set to 0 and 100000");
    ok = ok && m_raymanager.MakeCoordinateSystem( coords, 0.0, toaHeight ); // m_rayTracerSpecs is not yet configured 

	m_coords = *coords;

	ok = ok && m_diffusespecs.Initialize( *coords );

	ok = ok && m_raymanager.GetSun( &sun );

    if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_HR_Specs_Internal_Core::CreateCoordinates, Error creating coordinates." );
	return ok;
}

/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Core::CreateDiffuseTable		2013-06-17*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_Core::CreateDiffuseTable( std::unique_ptr<SKTRAN_HR_Diffuse_Table_CPU>& diffusetable )
{
	bool ok = true;

	std::unique_ptr<SKTRAN_HR_Diffuse_Table_CPU>  diffuse_local ( new SKTRAN_HR_Diffuse_Table_CPU );
	ok = ok && CreateDiffuseTableType( diffuse_local, false );
	ok = ok && CreateDiffusePoints( *diffuse_local ); 
	if (ok) diffusetable = std::move ( diffuse_local );
	if (!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_HR_Specs_Internal_Core::CreateDiffuseTable, Error creating diffuse table.");

	return ok;
}

bool SKTRAN_HR_Specs_Internal_Core::ReleaseResources()
{
	return true;
}


bool SKTRAN_HR_Specs_Internal_Core::RotateStartAndEnd( const HELIODETIC_UNITVECTOR& look, const HELIODETIC_UNITVECTOR& in, HELIODETIC_UNITVECTOR& out, double theta )
{
	bool ok = true;
	
	double costheta = nxmath::cosd( theta );
	double sintheta = nxmath::sind( theta );

	double ux = look.X();
	double uy = look.Y();
	double uz = look.Z();

	out.SetCoords( (costheta + ux*ux*(1-costheta)) * in.X() + (ux*uy*(1-costheta) - uz*sintheta)*in.Y() + (ux*uz*(1-costheta) + uy*sintheta)*in.Z(),
				   (uy*ux*(1-costheta) + uz*sintheta) * in.X() + (costheta + uy*uy*(1-costheta))*in.Y() + (uy*uz*(1-costheta) - ux*sintheta)*in.Z(),
				   (uz*ux*(1-costheta) - uy*sintheta) * in.X() + (uz*uy*(1-costheta) + ux*sintheta)*in.Y() + (costheta + uz*uz*(1-costheta))*in.Z());

	return ok;
}

bool SKTRAN_HR_Specs_Internal_Core::CalcReferencePoints( nxVector& in, nxVector& out )
{
	bool ok = true;

	ok = ok && m_raymanager.GetBoundingReferences( in, out );

	return ok;
}

bool SKTRAN_HR_Specs_Internal_Core::CreateDiffusePoints( SKTRAN_HR_Diffuse_Table_CPU& diffusetable )
{
	bool ok = true;

	std::vector<double> heights;
	std::vector<HELIODETIC_POINT> locations;	

	size_t numgroundpoints;
	size_t numprofiles;
	size_t numdiffusepoints;

	size_t count;
	std::vector<size_t> profilestartidx;

	HELIODETIC_POINT pt;
	double surfaceheight;

	ok = ok && CreateDiffuseHeights(&heights);
	ok = ok && CreateDiffuseLocations( locations );

	numgroundpoints = locations.size();
	numprofiles = locations.size();
	numdiffusepoints = locations.size() * heights.size();

	// make the diffuse points
	count = 0;
	diffusetable.AllocateDiffusePoints(numdiffusepoints, numgroundpoints);
	for (size_t profileidx = 0; profileidx < numprofiles; profileidx++)
	{
		profilestartidx.push_back(count);
		for (size_t heightidx = 0; heightidx < heights.size(); heightidx++)
		{
			pt.Initialize(locations[profileidx].Vector().UnitVector(), locations[profileidx].Vector().Magnitude() + heights[heightidx], CoordinatePtr());
			ok = ok && m_diffusespecs.MakeDiffusePoint(diffusetable.DiffusePointAtVar(count), pt, false);
			++count;
		}
	}
	diffusetable.SetGroundStartIdx(count);
	profilestartidx.push_back(count);

	// set diagnostic if applicable
	if (m_track)
	{
		std::vector<size_t> trackstartidx;
		trackstartidx.resize(m_trackeddiffprofs.size());
		for (size_t idx = 0; idx < m_trackeddiffprofs.size(); idx++)
		{
			if (m_trackeddiffprofs[idx] < numprofiles)
			{
				trackstartidx[idx] = profilestartidx[m_trackeddiffprofs[idx]];
			}
		}
		diffusetable.SetDiagnostic(trackstartidx, m_trackedscatords);
	}

	// and the ground points
	surfaceheight = RayTracerSpecs().SurfaceHeight() + 1e-6;
	for (size_t groundidx = 0; groundidx < numgroundpoints; groundidx++)
	{
		pt.Initialize(locations[groundidx].UnitVector(), locations[groundidx].Vector().Magnitude() + surfaceheight, CoordinatePtr());
		ok = ok && m_diffusespecs.MakeDiffusePoint(diffusetable.DiffusePointAtVar(count), pt, true);
		++count;
	}

	ok = ok && diffusetable.ConfigureIndexes();
	ok = ok && diffusetable.SetProfileStartIdx(profilestartidx);

	if (m_numinterp == 2)
	{
		// if we only have one diffuse plane then we can make a plane unit sphere which
		// handles the interpolation for us
		// otherwise just let the difufse table do the default interpolation
		std::vector<nxVector> profileloc;
		profileloc.resize(numprofiles);
		for (size_t idx = 0; idx < profileloc.size(); idx++)
		{
			profileloc[idx] = nxVector(locations[idx].Vector().UnitVector().X(), locations[idx].Vector().UnitVector().Y(), locations[idx].Vector().UnitVector().Z());
		}
		std::unique_ptr<SKTRAN_UnitSphere_Plane> sphere(new SKTRAN_UnitSphere_Plane);
		sphere->ConstructPlane(profileloc);
		diffusetable.SetDiffuseLocationSphere(std::move(sphere));
	}

	// Set the number of profiles used for each horizontal interpolation
	diffusetable.SetNumProfileInterp(m_numinterp);
	

	return ok;
}


bool SKTRAN_HR_Specs_Internal_Core::CreateOffLOSDiffuseLocations( std::vector<HELIODETIC_POINT>& locations )
{
	bool ok = true;

	double					maxrotateangle = m_angleofflook;
	size_t					numprofiles = m_numprofiles;
	size_t					numperp = m_numofflook;
	
	size_t					numalonglos;
	size_t					count = 0;
	HELIODETIC_VECTOR		in,out,temp;
	HELIODETIC_UNITVECTOR	temprotate;
	HELIODETIC_UNITVECTOR	rotateunit;
	HELIODETIC_UNITVECTOR	look;
	HELIODETIC_POINT		pt;
	double					totaldist;
	double					t;
	double					rotangle;
	nxVector				inreference;
	nxVector				outreference;
	size_t					refidx;
	bool					includeEntry, includeExit;
	HELIODETIC_UNITVECTOR	u1, u2;
	double					fraction;

	// get the line of sight endpoints
	ok = ok && CalcReferencePoints( inreference, outreference );


	switch (m_internaldiffuseplacementtype)
	{
		case InternalDiffusePlacementType::LinearLOS1D:
		{
			// calculate entry/exit points without osc sphere translation for consistency with older versions
			in = m_coords->GeographicToHelio(inreference);
			out = m_coords->GeographicToHelio(outreference);

			numperp = 1;
			if (numprofiles % 2 == 0) // enforce odd number of profiles
				--numprofiles;

			includeEntry = false;
			includeExit= false;
			refidx = numprofiles / 2; // replace center point with reference point
			break;
		}
		case InternalDiffusePlacementType::LinearLOS2D:
		{ 			
			// calculate entry/exit points without osc sphere translation for consistency with older versions
			in = m_coords->GeographicToHelio(inreference);
			out = m_coords->GeographicToHelio(outreference); 

			if (numperp % 2 == 0) // enforce odd number of planes
				--numperp;
			if (numprofiles % 2 == 0) // enforce odd number of profiles
				--numprofiles;

			includeEntry = false;
			includeExit = false;
			refidx = numprofiles / 2; // replace center point with reference point
			break;
		}
		case InternalDiffusePlacementType::LinearLOS1DEntry:
		{
			// calculate entry/exit points with osc sphere translation
			in = m_coords->GeographicToHelio(m_coords->TranslateGeoidToOsculatingSphere(inreference)); // may want to make in = observer for the in-atmosphere case (currently it traces the los backwards)
			out = m_coords->GeographicToHelio(m_coords->TranslateGeoidToOsculatingSphere(outreference));

			numperp = 1;
			includeEntry = true;
			includeExit = false;
			refidx = numprofiles; // replace point with reference point if nearby
			break;
		}
		case InternalDiffusePlacementType::LinearLOS1DExit:
		{
			// calculate entry/exit points with osc sphere translation
			in = m_coords->GeographicToHelio(m_coords->TranslateGeoidToOsculatingSphere(inreference));
			out = m_coords->GeographicToHelio(m_coords->TranslateGeoidToOsculatingSphere(outreference));

			numperp = 1;
			includeEntry = false;
			includeExit = true;
			refidx = numprofiles; // replace point with reference point if nearby
			break;
		}
		case InternalDiffusePlacementType::LinearLOS1DEntryExit:
		{
			// calculate entry/exit points with osc sphere translation
			in = m_coords->GeographicToHelio(m_coords->TranslateGeoidToOsculatingSphere(inreference));
			out = m_coords->GeographicToHelio(m_coords->TranslateGeoidToOsculatingSphere(outreference));

			numperp = 1;
			includeEntry = true; // place diffuse profile at entry/obs
			includeExit = true;
			refidx = numprofiles; // replace point with reference point if nearby
			break;
		}
		default:
		{
			nxLog::Record(NXLOG_ERROR, "SKTRAN_HR_Specs_Internal_Core::CreateOffLOSDiffuseLocations, Invalid diffuse placement type.");
			return false;
		}
	}

	// flag if the reference point lines up with the entry/exit point
	numalonglos = numprofiles;
	if (refidx == numprofiles)
	{
		if (includeEntry && (m_coords->ReferencePoint(0.0).UnitVector() & in.UnitVector()) > 0.9999996) // if vectors are within ~0.05 deg
			refidx = 0;
		else if (includeExit && (m_coords->ReferencePoint(0.0).UnitVector() & out.UnitVector()) > 0.9999996) // if vectors are within ~0.05 deg
			refidx = numprofiles - 1;
		else
			numalonglos = numprofiles - 1; // no match; place the reference point independently and place the remaining points along the los
	}

	totaldist = (out - in).Magnitude();
	look = (out - in).UnitVector();
	fraction = 1.0 / max(numalonglos + 1 - includeEntry - includeExit, (size_t)1); // distance between los points

	locations.resize(numprofiles*numperp);
	for (size_t idx = 0; idx < numalonglos; idx++)
	{
		// place points in los plane
		t = (idx + !includeEntry) * fraction * totaldist;
		temp.SetCoords(in.X() + look.X() * t,
			in.Y() + look.Y() * t,
			in.Z() + look.Z() * t);
		if (idx != refidx)
		{
			locations[count].Initialize(temp.UnitVector(), m_coords->AltitudeToRadius(0.0), CoordinatePtr());
		}
		else
		{
			locations[count] = m_coords->ReferencePoint(0.0);
		}
		++count;
		// place points in off-los planes
		for (size_t rotidx = 0; rotidx < numperp - 1; rotidx += 2)
		{
			rotangle = (maxrotateangle * 2) / (numperp - 1) * (rotidx / 2 + 1);
			//rotangle = maxrotateangle / (rotidx/2 + 1 );
			rotateunit = locations[count - 1].Vector().UnitVector();
			RotateStartAndEnd(look, rotateunit, temprotate, rotangle);
			locations[count].Initialize(temprotate, m_coords->AltitudeToRadius(0.0), CoordinatePtr());
			++count;
			RotateStartAndEnd(look, rotateunit, temprotate, -1.0 * rotangle);
			locations[count].Initialize(temprotate, m_coords->AltitudeToRadius(0.0), CoordinatePtr());
			++count;
		}
	}
	if (numalonglos < numprofiles) // place the reference point if it isn't already
	{
		locations[count] = m_coords->ReferencePoint(0.0);
		++count;
	}

	return ok;
}

bool SKTRAN_HR_Specs_Internal_Core::CreateManualDiffuseLocations(std::vector<HELIODETIC_POINT>& locations )
{
	bool ok = true;

	size_t numprofiles = m_diffuseprofilelocations.size();
	locations.resize(numprofiles);

	for( size_t profidx = 0; profidx < numprofiles; profidx++ )
	{
		locations[profidx].FromVector(m_coords->GeographicToHelio(m_diffuseprofilelocations[profidx]), CoordinatePtr());
	}

	// set interpolation parameters
	if (numprofiles == 1)
	{
		m_numinterp = 1;
	}
	else if (m_userdiffuseplacementtype == SKTRAN_HR_DiffuseProfilePlacement_LinearLOS)
	{
		m_numinterp = 2;
	}
	else
	{
		m_numinterp = 3;
	}

	return ok;
}

bool SKTRAN_HR_Specs_Internal_Core::CreateLatLonDiffuseLocations(std::vector<HELIODETIC_POINT>& locations )
{
	bool ok = true;

	std::vector<size_t> profilestartidx;
	size_t numprofiles = m_diffuseprofilelatlons.size() / 2;

	HELIODETIC_VECTOR hloc;
	locations.resize(numprofiles);
	for (size_t profidx = 0; profidx < numprofiles; profidx++)
	{
		// place diffuse points on the osculating sphere at the given lat/lon coordinates
		nxVector nxloc;
		HELIODETIC_POINT pt;
		nxloc.FromLatLong(m_diffuseprofilelatlons[profidx * 2], m_diffuseprofilelatlons[profidx * 2 + 1], m_coords->AltitudeToRadius(0.0));
		pt.FromVector(m_coords->GeographicToHelio(nxloc), CoordinatePtr());
		locations[profidx] = pt;
	}
	
	// set interpolation parameters
	if (numprofiles == 1)
	{
		m_numinterp = 1;
	}
	else if (m_userdiffuseplacementtype == SKTRAN_HR_DiffuseProfilePlacement_LinearLOS)
	{
		m_numinterp = 2;
	}
	else
	{
		m_numinterp = 3;
	}

	return ok;
}

bool SKTRAN_HR_Specs_Internal_Core::CreateManualLOSDiffuseLocations(std::vector<HELIODETIC_POINT>& locations )
{
	bool ok = true;


	size_t					numprofiles;
	size_t					tangentidx;
	size_t					count = 0;
	HELIODETIC_VECTOR		in, out, temp;
	HELIODETIC_UNITVECTOR	look;
	HELIODETIC_POINT		pt;
	double					totaldist;
	double					t;
	double					rotangle;
	nxVector				inreference;
	nxVector				outreference;
	std::vector<size_t>		profilestartidx;
	double					surfaceheight;

	// get the line of sight endpoints
	ok = ok && CalcReferencePoints(inreference, outreference);
	
	numprofiles = m_diffuseprofilelospositions.size();

	// convert to osculating sphere centered coordinates
	in = m_coords->GeographicToHelio(m_coords->TranslateGeoidToOsculatingSphere(inreference));
	out = m_coords->GeographicToHelio(m_coords->TranslateGeoidToOsculatingSphere(outreference));
	totaldist = (out - in).Magnitude();
	look = (out - in).UnitVector();

	locations.resize(numprofiles);
	for (size_t idx = 0; idx < numprofiles; idx++)
	{
		t = totaldist * m_diffuseprofilelospositions[idx];
		temp.SetCoords(in.X() + look.X() * t,
			in.Y() + look.Y() * t,
			in.Z() + look.Z() * t);
		locations[idx].Initialize(temp.UnitVector(), m_coords->AltitudeToRadius(0.0), CoordinatePtr());
	}

	// set interpolation parameters
	if (numprofiles == 1)
	{
		m_numinterp = 1;
	}
	else
	{
		m_numinterp = 2;
	}

	return ok;
}

bool SKTRAN_HR_Specs_Internal_Core::CreatePlaneDiffuseLocations(std::vector<HELIODETIC_POINT>& locations )
{
	bool ok = true;

	HELIODETIC_UNITVECTOR x = m_coords->GeographicToHelioUnitVector(m_diffuseplanereference);
	HELIODETIC_UNITVECTOR z = m_coords->GeographicToHelioUnitVector(m_diffuseplanenormal);
	//HELIODETIC_UNITVECTOR y = m_coords->GeographicToHelioUnitVector(m_diffuseplanenormal.Cross(m_diffuseplanereference));
	HELIODETIC_UNITVECTOR y = m_coords->GeographicToHelioUnitVector(m_diffuseplanereference.Cross(m_diffuseplanenormal));

	size_t numprofiles = m_diffuseplaneangles.size();
	locations.resize(numprofiles);

	for (size_t idx = 0; idx < numprofiles; idx++)
	{
		HELIODETIC_VECTOR u = HELIODETIC_VECTOR(x, nxmath::cosd(m_diffuseplaneangles[idx])) + HELIODETIC_VECTOR(y, nxmath::sind(m_diffuseplaneangles[idx]));
		locations[idx].Initialize(u.UnitVector(), m_coords->AltitudeToRadius(0.0), CoordinatePtr());
	}

	// set interpolation parameters
	if (numprofiles == 1)
	{
		m_numinterp = 1;
	}
	else
	{
		m_numinterp = 2;
	}

	return ok;
}


bool SKTRAN_HR_Specs_Internal_Core::CreateLinearSZADiffuseLocations( std::vector<HELIODETIC_POINT>& locations )
{
	bool ok = true;

	std::vector<double>		szas;
	std::vector<size_t>		profilestartidx;
	size_t numprofiles;

	ok = ok && CreateSZAs( szas );
	numprofiles = szas.size();

	locations.resize(numprofiles);
	for( size_t profileidx = 0; profileidx < numprofiles; profileidx++ )
	{
		HELIODETIC_UNITVECTOR hvec;
		hvec.SetCoords( nxmath::sind( szas[profileidx] ), 0.0, nxmath::cosd( szas[profileidx] ));
		locations[profileidx].Initialize(hvec, m_coords->AltitudeToRadius(0.0), CoordinatePtr());
	}

	// set interpolation parameters
	if (numprofiles == 1)
	{
		m_numinterp = 1;
	}
	else
	{
		m_numinterp = 2;
	}

	return ok;
}

bool SKTRAN_HR_Specs_Internal_Core::CreateSZAs( std::vector<double>& szas ) const
{
	bool ok = true;

	double refsza;
	double minsza;
	double maxsza;

	size_t numsza = m_numprofiles;


	ok = ok && m_raymanager.GetSZA( &refsza, &minsza, &maxsza );
	if ( m_internaldiffuseplacementtype == InternalDiffusePlacementType::ManualSZA )
	{
		szas = m_diffuseprofileszas;
	}
	else if( m_internaldiffuseplacementtype == InternalDiffusePlacementType::LinearSZAForceTP )
	{

		if( numsza % 2 == 0) 
			szas.resize( numsza+1 );
		else
			szas.resize( numsza );
	
		szas[0] = refsza;
		double topdiff = abs(maxsza - refsza);
		double mindiff = abs(minsza - refsza);
		size_t numiter = (szas.size()-1)/2;
		for( size_t szaidx = 0; szaidx < numiter; szaidx++ )
		{
			szas[1+szaidx] = topdiff / (numiter) * (szaidx+1) + refsza;
			szas[szas.size() - 1 - szaidx] = refsza - (mindiff / (numiter) * (szaidx+1));
		}
	}
	else if ( m_internaldiffuseplacementtype == InternalDiffusePlacementType::LinearSZA )
	{
		/* SASKTRANV21 PLACEMENT CODE
		szas.resize(numsza);
		double cosoffset;
		if( szas.size() == 1 )
			cosoffset = 0;
		else
			cosoffset = (nxmath::cosd(minsza) - nxmath::cosd(maxsza)) / (szas.size() - 1);
		int idxoffset = (int)(szas.size()-1)/2;

		for( int szaidx = 0; szaidx < szas.size(); szaidx++ )
		{
			szas[szaidx] = nxmath::acosd( nxmath::cosd(refsza) + ( szaidx - idxoffset ) * cosoffset );
		}
		END SASKTRANV21 PLACEMENT CODE */

		// we place diffuse profiles linearly in sza, including the max and min points
		szas.resize(numsza);
		if( numsza == 1 )
		{
			szas[0] = refsza;
		}
		else
		{
			for( int szaidx = 0; szaidx < (int)szas.size(); szaidx++ )
			{
				szas[szaidx] = minsza + szaidx * (maxsza - minsza) / (numsza-1);
			}
		}
	}
	else if( m_internaldiffuseplacementtype == InternalDiffusePlacementType::SmartSZA )
	{
		szas.resize(numsza);
		if( maxsza < 90.0 || (minsza > 90.0) )
		{
			// either fully before the terminator or fully after, place linearly
			if( numsza == 1 )
			{
				szas[0] = refsza;
			}
			else
			{
				for( int szaidx = 0; szaidx < (int)szas.size(); szaidx++ )
				{
					szas[szaidx] = minsza + szaidx * (maxsza - minsza) / (numsza-1);
				}
			}
		}
		else
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_HR_Specs_Internal_Core::CreateSZA option not currently supported");
			//double pastterminatorfraction = 0.7;		// place 70% of profiles past the terminator and 30% before
			//auto placementF = [](double x) { return x*x; };

			//size_t numpast = ceil(pastterminatorfraction * numsza);
			//size_t numbefore  = numsza - numpast;

			//size_t count = 0;
		}
	}
	else
	{
		nxLog::Record(NXLOG_WARNING, "SKTRAN_HR_Specs_Internal_Core::CreateSZAs error invalid diffuse placement type");
	}
	
	std::sort( szas.begin(), szas.end() );
	std::reverse( szas.begin(), szas.end() );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_Core::CreateDiffuseHeights		 2017- 2- 6*/
/** Creates the diffuse heights. This method needs to create the actual diffuse
	heights used by the model
 **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_Core::CreateDiffuseHeights( std::vector<double>* heights )
{
	bool	ok = true;
	bool	ok1;
	double	h;
	size_t	numdiffuse;
	double  startindex;
	double	endindex;
	double	minheight;
	double  lasth;
	double  minaltitude, maxaltitude;

	if( m_manualdiffuseheights.size() > 0 )
	{
		*heights = m_manualdiffuseheights;

		if( m_track )
		// ---- dump heights to H5 file if applicable
		{
			hid_t H5Fid, H5Sid, H5Did; // H5 id-s for file, dataset and space
			hid_t propList;
			H5Fid = H5Fopen( "DiagnosticData.h5", H5F_ACC_RDWR, H5P_DEFAULT );
			if( H5Fid < 0 )
			{
				nxLog::Record( NXLOG_WARNING, "Could not open h5 diagnostic file. Thats not good" );
				return false;
			}
			hsize_t dims1[1] = { heights->size() };

			H5Sid = H5Screate_simple( 1, dims1, NULL );
			propList = H5Pcreate( H5P_DATASET_CREATE );
			H5Pset_layout( propList, H5D_CHUNKED );
			H5Pset_chunk( propList, 1, dims1 );
			H5Did = H5Dcreate2( H5Fid, "heights", H5T_NATIVE_DOUBLE, H5Sid, H5P_DEFAULT, propList, H5P_DEFAULT );
			H5Dwrite( H5Did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, heights->data() );
			H5Dclose( H5Did );

			H5Pclose( propList );
			H5Sclose( H5Sid );
			H5Fclose( H5Fid );
		}
		return ok;
	}
	else
	{
		minaltitude = RayTracerSpecs().SurfaceHeight();
		maxaltitude = NXFINITE(m_diffusemaxheight )? m_diffusemaxheight : RayTracerSpecs().TOAHeight();
		startindex = floor( minaltitude/m_diffuseheightres - 0.5);	// Calculate the minimum start altitude on a regular spaced grid of heightres 
		endindex   = ceil ( maxaltitude/m_diffuseheightres + 0.5);	// Calculate the maximum end altitude on a regular spaced grid of heightres
		minheight  =  0.01 + minaltitude;								// Force the an additional point just above the ground. This will be our first diffuse point
		numdiffuse = (size_t)(endindex - startindex) + 3;					// Calculate the maximum number of points we need in the heights array.
		heights->resize(0);													// Resize the array to zero
		heights->reserve( numdiffuse);										// Reset the number of points required
		heights->push_back( minheight );									// force an additional point just above the ground makes a slight difference for some cases near 340nm
		for( size_t heightidx = 0; heightidx < numdiffuse; heightidx++ )	// For each possible diffuse point
		{
			h = (startindex + 0.5 + heightidx)*m_diffuseheightres;				// Get the height of diffuse point offset to the middle of the regular spaced grid
			if ( (h > (minheight+1) ) && ( h <= (maxaltitude-1)))				// If this height is in the range we desire then
			{																	// then
				heights->push_back(h);											// add it to our list
			}
		}
		if ( (maxaltitude - heights->back()) > 0.1 ) heights->push_back(maxaltitude);	// Cap of fthe list with the top of the atmopshere point
	}

	// ---- perform a few reality checks to make sure
	ok = ok && (heights->size() > 0);
	if (ok)
	{
		lasth = heights->at(0);
		for (size_t ih=1; (ih < heights->size()) && ok; ++ih)
		{
			h = heights->at(ih);
			ok = ok && (h > lasth);
			if (!ok)
			{
				nxLog::Record(NXLOG_WARNING,"SKTRAN_HR_Specs_Internal_Core::CreateDiffuseHeights, the heights of the diffuse points is not montonically increasing");
			}
			lasth = h;
		}
		if (ok)
		{

			ok1 =       (heights->front() - RayTracerSpecs().SurfaceHeight()) < 0.11;
			ok1 = ok1 && (heights->back () - RayTracerSpecs().TOAHeight()) > -0.11;
			if (!ok1)
			{
				if (!NXFINITE(m_diffusemaxheight ))			// dont print the message if user has used "hidden" property to override the setting
				{
					nxLog::Record(NXLOG_WARNING,"SKTRAN_HR_Specs_Internal_Core::CreateDiffuseHeights, The span of the diffuse heights (%e to %e) does not span the entire range of the atmosphere %e to %e",
					                        (double)heights->front(), (double)heights->back (), (double)RayTracerSpecs().SurfaceHeight(), (double)RayTracerSpecs().TOAHeight()     );
				}
			}
		}
		if (heights->back() < 2000)
		{
			nxLog::Record(NXLOG_WARNING,"SKTRAN_HR_Specs_Internal_Core::CreateDiffuseHeights, It looks like the diffuse heights are specificed in kilometers rather than meters, max diffuse height = %e meters. PLease specify diffuse heights in meters", (double)heights->back());
		}
	}

	if( m_track )
	// ---- dump heights to H5 file if applicable
	{
		hid_t H5Fid, H5Sid, H5Did; // H5 id-s for file, dataset and space
		hid_t propList;
		H5Fid = H5Fopen( "DiagnosticData.h5", H5F_ACC_RDWR, H5P_DEFAULT );
		if( H5Fid < 0 )
		{
			nxLog::Record( NXLOG_WARNING, "Could not open h5 diagnostic file. Thats not good" );
			return false;
		}
		hsize_t dims1[1] = { heights->size() };

		H5Sid = H5Screate_simple( 1, dims1, NULL );
		propList = H5Pcreate( H5P_DATASET_CREATE );
		H5Pset_layout( propList, H5D_CHUNKED );
		H5Pset_chunk( propList, 1, dims1 );
		H5Did = H5Dcreate2( H5Fid, "heights", H5T_NATIVE_DOUBLE, H5Sid, H5P_DEFAULT, propList, H5P_DEFAULT );
		H5Dwrite( H5Did, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, heights->data() );
		H5Dclose( H5Did );

		H5Pclose( propList );
		H5Sclose( H5Sid );
		H5Fclose( H5Fid );
	}
	return ok;
}

bool SKTRAN_HR_Specs_Internal_Core::CreateDiffuseLocations(std::vector<HELIODETIC_POINT>& locations )
{
	bool ok = true;

	switch (m_internaldiffuseplacementtype)
	{
		case InternalDiffusePlacementType::ManualLocation:
			ok = ok && CreateManualDiffuseLocations(locations);
			break;
		case InternalDiffusePlacementType::ManualLatLon:
			ok = ok && CreateLatLonDiffuseLocations(locations);
			break;
		case InternalDiffusePlacementType::ManualSZA:
			ok = ok && CreateLinearSZADiffuseLocations(locations);
			break;
		case InternalDiffusePlacementType::ManualLOS:
			ok = ok && CreateManualLOSDiffuseLocations(locations);
			break;
		case InternalDiffusePlacementType::ManualPlane:
			ok = ok && CreatePlaneDiffuseLocations(locations);
			break;
		case InternalDiffusePlacementType::LinearSZA:
		case InternalDiffusePlacementType::LinearSZAForceTP:
		case InternalDiffusePlacementType::SmartSZA:
			ok = ok && CreateLinearSZADiffuseLocations(locations);
			break;
		case InternalDiffusePlacementType::LinearLOS1D:
		case InternalDiffusePlacementType::LinearLOS1DEntry:
		case InternalDiffusePlacementType::LinearLOS1DExit:
		case InternalDiffusePlacementType::LinearLOS1DEntryExit:
		case InternalDiffusePlacementType::LinearLOS2D:
			ok = ok && CreateOffLOSDiffuseLocations(locations);
			break;
		default:
			ok = false;
	}
	return ok;
}


// currently deprecated, functionality is duplicated by CreateOffLOSDiffuseProfiles
bool SKTRAN_HR_Specs_Internal_Core::CreateLinearDiffuseLocations( std::vector<HELIODETIC_POINT>& locations )
{
	m_numofflook = 1;
	return CreateOffLOSDiffuseLocations( locations );
}
