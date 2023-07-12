#include "include/sktran_hr_internals.h"

bool SKTRAN_HR_Specs_Internal_wf::Configure( const SKTRAN_HR_Specs_User& specs )
{
	m_doextinction = specs.WeightingFunctionSpecsConst().m_doextinction;
	m_doscatextinction = specs.WeightingFunctionSpecsConst().m_doscatextinction;
	m_wfmode = specs.WeightingFunctionSpecsConst().m_wfmode;
	m_manualwfresolution = specs.WeightingFunctionSpecsConst().m_manualwfresolution;
	m_maxwfheight = specs.WeightingFunctionSpecsConst().m_maxwfheight;
	m_wfinterpwidth = specs.WeightingFunctionSpecsConst().m_wfinterpwidth;
	m_wfspecies = specs.WeightingFunctionSpecsConst().m_wfspecies;
	m_wfspeciesmode = specs.WeightingFunctionSpecsConst().m_wfspeciesmode;
	m_wfprecision = specs.WeightingFunctionSpecsConst().m_wfprecision;
	m_wfaerosolmode = specs.WeightingFunctionSpecsConst().m_aerosolmode;
	m_aerosolsizepercentchange = specs.WeightingFunctionSpecsConst().m_aerosolsizepercentchange;

	m_manualwfheights = specs.WeightingFunctionSpecsConst().m_manualwfheights;
	m_manualwfwidth   = specs.WeightingFunctionSpecsConst().m_manualwfwidth;
	m_manualwfwidthleft = specs.WeightingFunctionSpecsConst().m_manualwfwidthleft;
	m_manualwfwidthright = specs.WeightingFunctionSpecsConst().m_manualwfwidthright;


	if( m_manualwfheights.size() == 0 )
	{
		size_t numpert = (size_t) (ceil(m_maxwfheight / m_manualwfresolution));
		m_manualwfheights.resize(numpert);
		for( size_t idx = 0; idx < numpert; idx++ )
		{
			m_manualwfheights[idx] = (idx+0.5)*m_manualwfresolution;
		}
	}

	if( m_manualwfwidth.size() == 0 && (m_manualwfwidthleft.size() == 0) && (m_manualwfwidthright.size() == 0) )
	{
		m_manualwfwidth.resize( m_manualwfheights.size());
		for( size_t idx = 0; idx < m_manualwfwidth.size(); idx++ )
		{
			m_manualwfwidth[idx] = m_wfinterpwidth;
		}
	}


	return true;
}


bool SKTRAN_HR_Specs_Internal_wf::AddWfGeometryToRayTracer( SKTRAN_RayTracer_Straight_Generic& raytracer, std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords ) const
{
	bool ok = true;

	if( m_wfmode == SKTRAN_HR_wf_Mode_2d )
	{
		for( size_t angleidx = 0; angleidx < m_normals.size(); angleidx++ )
		{
			ok = ok && raytracer.AddGeometryObject( std::unique_ptr<SKTRAN_GeometryObject> ( new SKTRAN_GeometryObject_Plane( m_normals[angleidx] ) ) );
		}
		for( size_t altidx = 0; altidx < m_alts.size(); altidx++ )
		{
			ok = ok && raytracer.AddGeometryObject( std::unique_ptr<SKTRAN_GeometryObject> ( new SKTRAN_GeometryObject_Sphere( m_alts[altidx]+coords->AltitudeToRadius(0.0) ) ) );
		}
	}
	else
	{
		m_pertlist.AddGeometryToRayTracer(coords, raytracer);
	}
	return ok;
}



bool SKTRAN_HR_Specs_Internal_wf::MakePerturbationList( const SKTRAN_RayTracingRegionManager&					raymanager,
														 std::shared_ptr< const SKTRAN_CoordinateTransform_V2>& coords,
															  SKTRAN_HR_LinesOfSightTable&						linesofsight,
														const SKTRAN_HR_Specs_Internal_OpticalPropertiesTable&	opttablespecs )
{
	if( m_wfmode == SKTRAN_HR_wf_Mode_1d_los )
	{
		return MakeOneDimLOSPerturbation( linesofsight, coords );
	}
	else if ( m_wfmode == SKTRAN_HR_wf_Mode_1d_uniform )
	{
		return MakeOneDimUniformPerturbations();
	}
	else if ( m_wfmode == SKTRAN_HR_wf_Mode_2d )
	{
		return MakeTwoDimPerturbation( coords, opttablespecs );
	}
	else
	{
		return true;
	}
}



bool SKTRAN_HR_Specs_Internal_wf::MakeTwoDimPerturbation(  std::shared_ptr< const SKTRAN_CoordinateTransform_V2> & coords,
														  const SKTRAN_HR_Specs_Internal_OpticalPropertiesTable& opttablespecs )
{
	size_t numheight = m_manualwfheights.size();
	size_t numangle  = opttablespecs.m_anglegrid.size();

	nxVector x,y;
	HELIODETIC_VECTOR xh,yh;
	xh = coords->GeographicToHelio( opttablespecs.m_reference );
	yh = coords->GeographicToHelio( opttablespecs.m_reference.Cross(opttablespecs.m_normal));
	x.SetCoords( xh.X(), xh.Y(), xh.Z() );
	y.SetCoords( yh.X(), yh.Y(), yh.Z() );
	
	m_normals.resize( numangle );
	nxVector normal1, normal2;

	//const int numpert = numangle * numheight;
	// Allocate a contiguous block of memory for the perturbations
	// These are placed in a unique pointer below so we do not have to worry about freeing the memory later

	for( size_t angleidx = 0; angleidx < numangle; angleidx++ )
	{
		double th1 = opttablespecs.m_anglegrid[angleidx];
		double th2;
		if( angleidx == numangle -1 )
		{
			th2 = opttablespecs.m_anglegrid[angleidx-1];
		}
		else
		{
			th2 = opttablespecs.m_anglegrid[angleidx+1];
		}
		double thsep = abs(th2 - th1);
		normal1 = nxmath::cosd(th1-thsep+90.0) * x + nxmath::sind(th1-thsep+90.0) * y;
		normal2 = nxmath::cosd(th1+thsep+90.0) * x + nxmath::sind(th1+thsep+90.0) * y;
		for( size_t altidx = 0; altidx < numheight; altidx++ )
		{
			SKTRAN_HR_Perturbation_Absorption_Box pert;
			pert.Initialize( m_manualwfheights[altidx], m_manualwfwidth[altidx], normal1, normal2, 1.0 );
			m_pertlist.add_item(pert);
		}
		m_normals[angleidx] = normal1;
	}
	m_normals[m_normals.size()-1] = normal2;
	for( size_t altidx = 0; altidx < numheight; altidx++ )
	{
		m_alts.push_back(m_manualwfheights[altidx] - m_manualwfwidth[altidx]);
		m_alts.push_back(m_manualwfheights[altidx] + m_manualwfwidth[altidx]);
	}
	std::sort( std::begin(m_alts), std::end(m_alts) );
	m_alts.erase( std::unique( std::begin(m_alts), std::end(m_alts), [] (double a, double b) { return abs(a-b) < 1E-10; } ), std::end( m_alts ) );


	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_wf::MakeOneDimUniformPerturbations		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/

bool SKTRAN_HR_Specs_Internal_wf::MakeOneDimUniformPerturbations()
{

	size_t numpert = m_manualwfheights.size();
	for( size_t idx = 0; idx < numpert; idx++ )
	{
		SKTRAN_HR_Perturbation_Absorption_Linear perttest;
		if (m_manualwfwidthleft.size() > 0 && m_manualwfwidthright.size() > 0)
		{
			perttest.Initialize(m_manualwfheights[idx], m_manualwfwidthleft[idx], m_manualwfwidthright[idx], 1.0);
		}
		else
		{
			perttest.Initialize(m_manualwfheights[idx], m_manualwfwidth[idx], m_manualwfwidth[idx], 1.0);
		}
		m_pertlist.add_item(perttest);
	}

	return true;
}


/*-----------------------------------------------------------------------------
 *					SKTRAN_HR_Specs_Internal_wf::MakeOneDimLOSPerturbation		 2014- 11- 13*/
/** **/
/*---------------------------------------------------------------------------*/
bool SKTRAN_HR_Specs_Internal_wf::MakeOneDimLOSPerturbation( SKTRAN_HR_LinesOfSightTable&							 linesofsight,
															 std::shared_ptr< const SKTRAN_CoordinateTransform_V2>&	 coords )
{
	std::unique_ptr<SKTRAN_RayStorage_Straight_HR>			trajectorystorageobject( new SKTRAN_RayStorage_Straight_HR(coords) );
	std::shared_ptr<SKTRAN_RayTracer_Shells>				raytracerobject        ( new SKTRAN_RayTracer_Shells(coords) );



	SKTRAN_RayOptical_Straight ray( std::move(trajectorystorageobject), std::move(raytracerobject) );


	double rtprev = 0;
	double pertwidthabove, pertwidthbelow;

	for( size_t idx = 0; idx < linesofsight.LinesOfSightArray()->NumRays(); idx++ )
	{	
		SKTRAN_HR_Perturbation_Absorption_Linear perttest;

		const SKTRAN_LineOfSightEntry_V2* entry;
		linesofsight.LinesOfSightArray()->GetRay(idx, &entry );
		ray.MoveObserver( coords->GeographicToHelio( entry->Observer() ), coords->GeographicToHelio( entry->Look() ).UnitVector() );
		double robs, tobs, rt;
		ray.CalculateBaseLineTangentPointDetails( 0.0, &robs, &tobs, &rt );

		double alt = coords->RadiusToAltitude(rt);

		if( idx != linesofsight.LinesOfSightArray()->NumRays() - 1 )
		{
			linesofsight.LinesOfSightArray()->GetRay(idx+1,&entry);
			ray.MoveObserver( coords->GeographicToHelio( entry->Observer() ), coords->GeographicToHelio( entry->Look() ).UnitVector() );
			double rtnext;
			ray.CalculateBaseLineTangentPointDetails( 0.0, &robs, &tobs, &rtnext );
			if( idx != 0 )
			{
				pertwidthbelow = pertwidthabove;
			}
			pertwidthabove = abs(rt - rtnext);
			if( idx == 0 )
			{
				pertwidthbelow = 500;
			}
		}
		else
		{
			pertwidthabove = 500;
			pertwidthbelow = abs(rt - rtprev);
		}
		//pertwidth = 1000;
		//if( idx == 0 )
		//{
		//	std::unique_ptr<SKTRAN_HR_Perturbation_Absorption> pert ( new SKTRAN_HR_Perturbation_Absorption );
		//	pert->Initialize( alt + pertwidthabove, 0, 1.0 );
		//	m_pertlist.emplace_back( std::move( pert ) );
		//}
		//else
		//{
			perttest.Initialize( alt, pertwidthbelow, pertwidthabove, 1.0 );
			m_pertlist.add_item(perttest);
		//}
		//perttest->Initialize( alt, 2000, 2000, 1.0);
		rtprev = rt;
	}

	return true;
}

