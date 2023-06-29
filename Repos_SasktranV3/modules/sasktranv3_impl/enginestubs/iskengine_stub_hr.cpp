#include "../dllimplementation/stdafx.h"
#include "modules/sasktranv3_impl/sktranif_impl_helperclasses.h"
#include <boost/regex.hpp>

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::ISKEngine_Stub_HR		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEngine_Stub_HR::ISKEngine_Stub_HR()
{
	m_radiance.SetReuseMemory(true);
	m_numordersofscatter = 50;
	m_numthreads = 0;
	m_updateclimatology  = true;
	m_geometryisconfigured = false;
	m_storeStokes = false;
	m_issetdiagnostics = false;
	m_storeOpticalDepth = false;
	MakeScalarSetFunctions();
	MakeVectorSetFunctions();
	MakeScalarGetFunctions();
	MakeVectorGetFunctions();
	MakeObjectSetFunctions();
	MakeStringSetFunctions();
	//MakeDefaultOpticalState();

	m_specs.IntegratorSpecs().SetMaxOpticalDepth( 100000 );
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::~ISKEngine_Stub_HR		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

ISKEngine_Stub_HR::~ISKEngine_Stub_HR()
{
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::AddLineOfSight		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::AddLineOfSight( double mjd,  const nxVector& obs, const nxVector& look, int* losindex  )
{
	bool		ok;

	ok = m_linesofsight.AddLineOfSight( obs, look, mjd );
	if (ok)
	{
		*losindex = (int)m_linesofsight.NumRays() - 1;
	}
	else
	{
		*losindex = -999999;
	}
	m_geometryisconfigured = false;
	m_radiance.erase();
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::CheckModelNotInitalized		 2015- 11- 12*/
/** **/
/*---------------------------------------------------------------------------*/

bool  ISKEngine_Stub_HR::CheckModelNotInitalized( const char* propertystr) const
{
	if (m_geometryisconfigured)  nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Cannot change property %s after the geometry has been initialized", (const char*)propertystr);
	return !m_geometryisconfigured;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::MakeScalarSetFunctions		 2015- 11- 12*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::MakeScalarSetFunctions()
{
	AddSetScalarFunction("storeopticaldepth",
		[&, this](double d)
		{
			int specifier = (int)ceil(d-0.5);

			if (specifier == 1)
			{
				m_storeOpticalDepth = true;
			}
			else
			{
				m_storeOpticalDepth = false;
			}

			return true;
		}
	);

	AddSetScalarFunction("numordersofscatter",
		[&, this](double d)
		{ 
			bool ok;
			int scatorder = (int)ceil(d-0.5);

			ok = (scatorder >= 1 ) && (scatorder < 1000);
			if (!ok) nxLog::Record(NXLOG_WARNING,"ISKEngine HR, Invalid number of orders of scatter (%d). Use a number between 1 and 1000 inclusive", (int)scatorder);
			else
			{
				m_numordersofscatter = scatorder;
				m_specs.SetScatterOrder(scatorder);
			}
			return ok;
		}
	);

	AddSetScalarFunction("useshellraytracer", 
		[&, this](double d )
		{
			int specifier = (int)ceil(d-0.5);
			bool ok = true;

			if( specifier == 1 )
			{	
				ok =  CheckModelNotInitalized("useshellraytracer");
				if (ok)
				{
					ok = ok && m_specs.RayTracingSpecs().SetLinesOfSightType( SKTRAN_HR_RayTracer_Shells );
					ok = ok && m_specs.RayTracingSpecs().SetSolarType( SKTRAN_HR_RayTracer_Shells );
					ok = ok && m_specs.RayTracingSpecs().SetDiffuseType( SKTRAN_HR_RayTracer_Shells );
					if (!ok)
					{
						nxLog::Record(NXLOG_WARNING,"ISKEngine HR, There were errors configuring the model to useshellraytracer 1");
					}
				}
			}
			else
			{
				ok = false;
				nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Unknown specifier (%d) for property useshellraytracer", (int)specifier);
			}
			return true;
		}
	);

	AddSetScalarFunction( "userefraction",
		[&, this](double d)
		{
			int specifier = (int)ceil(d - 0.5);
			bool ok = true;

			if (specifier == 1)
			{
				ok = CheckModelNotInitalized("userefraction");
				if (ok)
				{
					ok = ok && m_specs.RayTracingSpecs().SetLinesOfSightType(SKTRAN_HR_RayTracer_Curved);
					if (!ok)
					{
						nxLog::Record(NXLOG_WARNING, "ISKEngine HR, There were errors configuring the model to userefraction");
					}
				}
			}
			else if (specifier == 2)
			{
				ok = CheckModelNotInitalized("userefraction");
				if (ok)
				{
					ok = ok && m_specs.RayTracingSpecs().SetLinesOfSightType(SKTRAN_HR_RayTracer_Curved);
					ok = ok && m_specs.RayTracingSpecs().SetSolarType(SKTRAN_HR_RayTracer_Curved);
					if (!ok)
					{
						nxLog::Record(NXLOG_WARNING, "ISKEngine HR, There were errors configuring the model to userefraction");
					}
				}
			}
			else if (specifier == 3)
			{
				ok = CheckModelNotInitalized("userefraction");
				if (ok)
				{
					ok = ok && m_specs.RayTracingSpecs().SetLinesOfSightType(SKTRAN_HR_RayTracer_Curved);
					ok = ok && m_specs.RayTracingSpecs().SetSolarType(SKTRAN_HR_RayTracer_Curved);
					ok = ok && m_specs.RayTracingSpecs().SetDiffuseType(SKTRAN_HR_RayTracer_Curved);
					if (!ok)
					{
						nxLog::Record(NXLOG_WARNING, "ISKEngine HR, There were errors configuring the model to userefraction");
					}
				}
			}
			else
			{
				ok = false;
				nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Unknown specifier (%d) for property userefraction", (int)specifier);
			}
			return ok;
		}
	);

	AddSetScalarFunction( "calcwf",
		[&, this](double d)
		{
			int specifier = (int)ceil(d-0.5);
			bool ok;

			ok =  CheckModelNotInitalized("calcwf");
			if (ok)
			{
				if( specifier == 0 )
				{
					m_specs.WeightingFunctionSpecs().SetWeightingFunctionMode( SKTRAN_HR_wf_Mode_None );
				}
				else if( specifier == 1 )
				{
					m_specs.WeightingFunctionSpecs().SetWeightingFunctionMode( SKTRAN_HR_wf_Mode_1d_los );
				}
				else if ( specifier == 2 )
				{
					m_specs.WeightingFunctionSpecs().SetWeightingFunctionMode( SKTRAN_HR_wf_Mode_1d_uniform );
				}
				else if ( specifier == 3 )
				{
					m_specs.WeightingFunctionSpecs().SetWeightingFunctionMode( SKTRAN_HR_wf_Mode_2d );
				}
				else
				{
					ok = false;
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Unknown specifier (%d) for property calcwf", (int)specifier);
				}
			}
			return ok;
		}
	);

    AddSetScalarFunction( "scattermatrixstoragemethod",
		[&, this](double d)
		{
			int specifier = (int)ceil(d-0.5);
			bool ok = true;

			if( 0==specifier )
			{
				m_specs.DiffuseSpecs().SetScatterMatrixStorageMethod( SKTRAN_HR_DiffuseMatrixStorageMethod::scalar );
			}
			else if( 1==specifier )
			{
				m_specs.DiffuseSpecs().SetScatterMatrixStorageMethod( SKTRAN_HR_DiffuseMatrixStorageMethod::scatter_cache );
			}
			else if ( 2==specifier )
			{
				m_specs.DiffuseSpecs().SetScatterMatrixStorageMethod( SKTRAN_HR_DiffuseMatrixStorageMethod::scatter_interpolate );
			}
			else if ( 3==specifier )
			{
				m_specs.DiffuseSpecs().SetScatterMatrixStorageMethod ( SKTRAN_HR_DiffuseMatrixStorageMethod::phase );
			}
			else
			{
				ok = false;
				nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Unknown specifier (%d) for property scattermatrixstoragemethod", (int)specifier);
			}

			return ok;
		}
	);

	AddSetScalarFunction("aerowfmode",
		[&, this](double d)
		{
			int specifier = (int)ceil(d - 0.5);
			bool ok;

			ok = CheckModelNotInitalized("aerowfmode");
			if (ok)
			{
				if (specifier == 0)
				{
					m_specs.WeightingFunctionSpecs().SetAerosolWeightingFunctionMode(SKTRAN_HR_wf_aerosol_Mode_numberdensity);
				}
				else if (specifier == 1)
				{
					m_specs.WeightingFunctionSpecs().SetAerosolWeightingFunctionMode(SKTRAN_HR_wf_aerosol_Mode_moderadius);
				}
				else if (specifier == 2)
				{
					m_specs.WeightingFunctionSpecs().SetAerosolWeightingFunctionMode(SKTRAN_HR_wf_aerosol_Mode_modewidth);
				}
				else
				{
					ok = false;
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Unknown specifier (%d) for property calcwf", (int)specifier);
				}
			}
			return ok;
		}
	);

	AddSetScalarFunction( "diffuseplacementtype",
		[&, this](double d)
		{
			int specifier = (int)ceil(d-0.5);
			bool ok;

			ok =  CheckModelNotInitalized("diffuseplacementtype");
			if (ok)
			{

				if( specifier == 0 )
				{
					m_specs.DiffuseSpecs().SetDiffusePlacementType( SKTRAN_HR_DiffuseProfilePlacement_LinearLOS );
				}
				else if( specifier == 1 )
				{
					m_specs.DiffuseSpecs().SetDiffusePlacementType( SKTRAN_HR_DiffuseProfilePlacement_LinearSZAForceTP );
				}
				else if( specifier == 2 )
				{
					m_specs.DiffuseSpecs().SetDiffusePlacementType( SKTRAN_HR_DiffuseProfilePlacement_LinearSZA );
				}
				else
				{
					ok = false;
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Unknown specifier (%d) for property diffuseplacementtype", (int)specifier);
				}
			}
			return ok;
		}
	);

	AddSetScalarFunction( "polarizationhigherorderfraction",
		[&, this](double d)
		{
			bool ok = d >= 0.0;

			if (ok)
			{
				m_specs.OpticalPropertiesSpecs().SetPolarizationHigherOrderFraction(d);
			}

			return ok;
		}
	);

	AddSetScalarFunction( "forcev21diffuseincoming",
		[&, this](double d)
		{
			int specifier = (int)ceil(d-0.5);
			bool ok;

			ok =  CheckModelNotInitalized("forcev21diffuseincoming");
			if (ok)
			{
				if( specifier == 0 )
				{
					m_specs.DiffuseSpecs().SetDiffuseIncomingType( SKTRAN_HR_DiffuseIncomingType_Default );
				}
				else if ( specifier == 1 )
				{
					m_specs.DiffuseSpecs().SetDiffuseIncomingType( SKTRAN_HR_DiffuseIncomingType_HardCode );
				}
				else if( specifier == 2 )
				{
					m_specs.DiffuseSpecs().SetDiffuseIncomingType( SKTRAN_HR_DiffuseIncomingType_v21Shifted );
				}
				else
				{
					ok = false;
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Unknown specifier (%d) for property forcev21diffuseincoming", (int)specifier);
				}
			}
			return ok;
		}
	);

	AddSetScalarFunction( "horizonsize",
		[&, this](double d )
		{
			bool ok;

			ok =  CheckModelNotInitalized("horizonsize");
			if (ok)
			{
				m_specs.DiffuseSpecs().SetHorizonSize( d );
			}
			return ok;
		}
	);

	AddSetScalarFunction( "numthreads",
		[&, this](double d)
		{
			int numthreads = (int)ceil(d-0.5);
			bool ok;


			ok =  CheckModelNotInitalized("numthreads");
			if (ok)
			{
				ok = (numthreads >= 0 && numthreads < 1000);
				if (!ok)
				{
					nxLog::Record( NXLOG_WARNING,"ISKEngine HR, Invalid number of threads (%d) for property numthreads. Enter a number between 0 and 1000 inclusive", (int)numthreads);
				}
				else
				{
					m_numthreads = numthreads;
				}
			}
			return ok;
		}
	);

	AddSetScalarFunction( "integrationtechnique",
		[&, this](double d)
		{
			int specifier = (int)ceil(d-0.5);
			bool ok;

			ok =  CheckModelNotInitalized("integrationtechnique");
			if (ok)
			{
				if( specifier == 0 )
				{
					m_specs.IntegratorSpecs().UseLegacySasktran21Technique(true);
				}
				else if( specifier == 1 )
				{
					m_specs.IntegratorSpecs().UseLegacySasktran21Technique(false);
					return true;
				}
				else if ( specifier == 2 )
				{
					m_specs.IntegratorSpecs().SetIntegratorType( SKTRAN_HR_IntegratorType_Constant );
				}
				else
				{
					ok = false;
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Unknown specifier (%d) for property integrationtechnique", (int)specifier);
				}
			}
			return ok;
		}
	);

	AddSetScalarFunction( "solarraytracingshells",
		[&, this](double d)
		{
			bool ok;

			ok =  CheckModelNotInitalized("solarraytracingshells");
			if (ok)
			{
				ok = m_specs.RayTracingSpecs().SetSolarShellSpacing( d );
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"ISKEngine HR, Error setting solar shell spacing to %e", (double)d);
				}
			}
			return ok;
		}
	);

	AddSetScalarFunction( "aerosolwfsizepercentchange",
		[&, this](double d)
		{
			bool ok;

			ok = CheckModelNotInitalized("aerosolwfsizepercentchange");
			if (ok)
			{
				m_specs.WeightingFunctionSpecs().SetAerosolSizePercentChange(d);
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Error setting aerosol wf size percent change to %e", (double)d);
				}
			}
			return ok;
		}
	);


	AddSetScalarFunction("raytracingshells",
		[&, this](double d)
		{
			bool ok;

			ok =  CheckModelNotInitalized("raytracingshells");
			if (ok)
			{
				ok = m_specs.RayTracingSpecs().SetShellSpacing( d );
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"ISKEngine HR, Error setting shell spacing to %e", (double)d);
				}
			}
			return ok;
		}
	);


	AddSetScalarFunction( "maxopticaldepthofcell", 
		[&, this](double d)
		{
			bool ok;

			ok =  CheckModelNotInitalized("maxopticaldepthofcell");
			if (ok)
			{
				ok = ( d > 0 );
				if (ok)
				{
					m_specs.IntegratorSpecs().SetMaxOpticalDepth( d );
				}
				else
				{
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Invalid Maximum Optical Depth of Cell [%e] entered for property maxopticaldepthofcell", d );
					ok = false;
				}
			}
			return ok;
		}
	);

	AddSetScalarFunction("maxadaptiveopticaldepthofray",
		[&, this](double d)
	{
		bool ok;

		ok = CheckModelNotInitalized("maxadaptiveopticaldepthofray");
		if (ok)
		{
			ok = (d > 0);
			if (ok)
			{
				m_specs.IntegratorSpecs().SetMaxRayOpticalDepth(d);
			}
			else
			{
				nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Invalid Maximum Optical Depth of Cell [%e] entered for property maxadaptiveopticaldepthofray", d);
				ok = false;
			}
		}
		return ok;
	}
	);

	AddSetScalarFunction( "minextinctionratioofcell",
		[&, this](double d)
		{
			bool ok;

			ok =  CheckModelNotInitalized("minextinctionratioofcell");
			if (ok)
			{
				ok = ( d > 0 && d <= 1.0 );
				if (ok)
				{
					m_specs.IntegratorSpecs().SetMaxExtinctionGradient( d );
				}
				else
				{
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Invalid minimim extinction ratio of cell [%e] entered", (double)d );
				}
			}
			return ok;
		}
	);

	AddSetScalarFunction( "numdiffuseprofilesinplane", 
		[&, this](double d)
		{
			int num = (int)ceil(d-0.5);
			bool ok;

			ok =  CheckModelNotInitalized("numdiffuseprofilesinplane");
			if (ok)
			{
				ok = m_specs.DiffuseSpecs().SetNumProfiles( num );
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"ISKEngine HR, Error setting the number of diffuse profiles in plane to %d", (int)num);
				}
			}
			return ok;
		}
	);

	AddSetScalarFunction( "forcelineardiffuseplacement", 
		[&, this](double d)
		{
			int specifier = (int)ceil(d-0.5);
			bool ok;

			ok =  CheckModelNotInitalized("forcelineardiffuseplacement");
			if (ok)
			{
				if( specifier == 0 )
				{
					m_specs.DiffuseSpecs().SetForceLinearPlacement(false);
				}
				else if( specifier == 1 )
				{
					m_specs.DiffuseSpecs().SetForceLinearPlacement(true);
				}
				else
				{
					ok = false;
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Unknown specifier (%d) for property forcelineardiffuseplacement", (int)specifier);
				}
			}
			return ok;
		}
	);

	AddSetScalarFunction( "numdiffuseoutgoing",
		[&, this](double d)
		{
			int num = (int)ceil(d-0.5);
			bool ok;

			ok =  CheckModelNotInitalized("numdiffuseoutgoing");
			if (ok)
			{
				m_specs.DiffuseSpecs().SetNumOutgoing(num);
			}
			return ok;
		}
	);

	AddSetScalarFunction( "diffuseheightres",
		[&, this](double d)
		{
			bool ok;

			ok =  CheckModelNotInitalized("diffuseheightres");
			if (ok)
			{
				m_specs.DiffuseSpecs().SetHeightRes( d );
			}
			return ok;
		}
	);

	AddSetScalarFunction( "maxdiffuseheight",
		[&, this](double d)
		{
			bool ok;

			ok =  CheckModelNotInitalized("maxdiffuseheight");
			if (ok)
			{
				m_specs.DiffuseSpecs().SetMaxDiffuseHeight( d );
			}
			return ok;
		}
	);

	AddSetScalarFunction( "numdiffuseplanes",
		[&, this](double d)
		{
			int num = (int)ceil(d-0.5);
			bool ok;

			ok =  CheckModelNotInitalized("numdiffuseplanes");
			if (ok)
			{
				m_specs.DiffuseSpecs().SetNumoffLook( num );
			}
			return ok;
		}
	);

	AddSetScalarFunction( "diffusemaxangleoffplane",
		[&, this](double d)
		{
			bool ok;

			ok =  CheckModelNotInitalized("diffusemaxangleoffplane");
			if (ok)
			{
				m_specs.DiffuseSpecs().SetAngleOffLook( d );
			}
			return ok;
		}
	);

	AddSetScalarFunction( "solartablenumcossza",
		[&, this](double d)
		{
			bool ok;

			ok = CheckModelNotInitalized("solartablenumcossza");
			if (ok)
			{
				int num = (int)ceil(d - 0.5);
				m_engine.SetNumCosSza((size_t)num);
			}
			return ok;
		}
	);

	AddSetScalarFunction( "opticalpropertiesheightres",
		[&, this](double d)
		{
			bool ok;

			ok =  CheckModelNotInitalized("opticalpropertiesheightres");
			if (ok)
			{
				m_specs.OpticalPropertiesSpecs().SetHeightRes( d );
			}
			return ok;
		}
	);
	
	AddSetScalarFunction( "opticaltabletype",
		[&, this](double d )
		{
			int specifier = (int)ceil(d-0.5);
			bool ok;

			ok =  CheckModelNotInitalized("opticaltabletype");
			if (ok)
			{
				if( specifier == 0 )
				{
					m_specs.OpticalPropertiesSpecs().SetOpticalPropertiesType( SKTRAN_HR_OpticalPropertiesTableType_1d );
				}
				else if( specifier == 1 )
				{
					m_specs.OpticalPropertiesSpecs().SetOpticalPropertiesType( SKTRAN_HR_OpticalPropertiesTableType_3D_UnitSphere );
				}
				else if( specifier == 2 )
				{
					m_specs.OpticalPropertiesSpecs().SetOpticalPropertiesType( SKTRAN_HR_OpticalPropertiesTableType_LOSPlane );
				}
				else if( specifier == 3 )
				{
					m_specs.OpticalPropertiesSpecs().SetOpticalPropertiesType( SKTRAN_HR_OpticalPropertiesTableType_SZA );
				}
				else if (specifier == 4)
				{
					m_specs.OpticalPropertiesSpecs().SetOpticalPropertiesType(SKTRAN_HR_OpticalPropertiesTableType_1d_ConstantLayers);
				}
				else
				{
					ok = false;
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Unknown specifier (%d) for property opticaltabletype", (int)specifier);
				}
			}
			return ok;
		}
	);

	AddSetScalarFunction( "wfprecision",
		[&, this](double d )
		{
			int specifier = (int) ceil(d-0.5);
			bool ok;

			ok =  CheckModelNotInitalized("wfprecision");
			if (ok)
			{
				if( specifier == 0 )
				{
					m_specs.WeightingFunctionSpecs().SetWeightingFunctionprecision( SKTRAN_HR_wf_precision_all );
				}
				else
				{
					m_specs.WeightingFunctionSpecs().SetWeightingFunctionprecision( SKTRAN_HR_wf_precision_onlylos );
				}
			}
			return ok;
		}
	);

	AddSetScalarFunction( "forceopticalcacheupdates",
		[&, this](double d)
		{
			int specifier = (int) ceil(d-0.5);
			bool ok;

			ok =  CheckModelNotInitalized("forceopticalcacheupdates");
			if (ok)
			{
				if( specifier ==0 )
				{
					m_specs.OpticalPropertiesSpecs().SetForceCacheUpdates( false );
				}
				else
				{
					m_specs.OpticalPropertiesSpecs().SetForceCacheUpdates( true );
				}
			}
			return ok;
		}
	);

    AddSetScalarFunction( "polarizationtype",
		[&, this](double d )
		{

			bool ok;
			int specifier = (int)ceil(d-0.5);

			ok =  CheckModelNotInitalized("polarizationtype");
			if (ok)
			{
				ok = SetPolarizationMode( specifier );
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Error setting polarization mode (%d)", (int)specifier);
				}
			}

			return ok;

		}
	);

	AddSetScalarFunction( "groundshiftalt",
		[&, this](double d)
		{
			bool ok;

			ok =  CheckModelNotInitalized("groundshiftalt");
			if (ok)
			{
				ok = m_specs.RayTracingSpecs().SetGroundShiftAlt( d );
				if (!ok) nxLog::Record(NXLOG_WARNING,"ISKEngine HR, Error setting property groundshiftalt to %e", (double)d);
			}
			return ok;
		}
	);

	AddSetScalarFunction( "surfaceheight",
		[&, this](double d)
		{
			bool ok;

			ok =  CheckModelNotInitalized("surfaceheight");
			if (ok)
			{
				ok = m_specs.RayTracingSpecs().SetGroundShiftAlt( d );
				if (!ok) nxLog::Record(NXLOG_WARNING,"ISKEngine HR, Error setting property surfaceheight to %e", (double)d);
			}
			return ok;
		}
	);

	AddSetScalarFunction( "toaheight",
		[&, this](double d)
		{
			bool ok;

			ok =  CheckModelNotInitalized("toaheight");
			if (ok)
			{
				ok =       m_specs.RayTracingSpecs().SetTOAHeight( d );
				if (!ok) nxLog::Record(NXLOG_WARNING,"ISKEngine HR, Error setting property toaheight to %e", (double)d);
			}
			return ok;
		}
	);

	AddSetScalarFunction( "earthradius",
		[&, this](double d)
		{
			bool ok;

			ok =  CheckModelNotInitalized("earthradius");
			if (ok)
			{
				ok = m_engine.InternalSpecsVar()->RayManager().SetEarthRadius(d) ;
				if (!ok) nxLog::Record(NXLOG_WARNING,"ISKEngine HR, Error setting property earthradius to %e", (double)d);
			}
			return ok;
		}
	);


    AddSetScalarFunction( "usesolartransmission",
		[&, this](double d )
		{
			bool ok = true;
			int specifier = (int)ceil(d-0.5);

			if (0 == specifier) {
				m_specs.IntegratorSpecs().SetUseSolarTransmission( false );
			}else {
				m_specs.IntegratorSpecs().SetUseSolarTransmission( true );
			}
			return ok;
		}
	);
    
    AddSetScalarFunction( "useemissions",
		[&, this](double d )
		{
			bool ok = true;
			int specifier = (int)ceil(d-0.5);

			if (0 == specifier) {
				m_specs.IntegratorSpecs().SetUseEmissions( false );
			}else {
				m_specs.IntegratorSpecs().SetUseEmissions( true );
			}
			return ok;
		}
	);

	AddSetScalarFunction("disablelosscattering",
		[&, this](double d)
	{
		bool ok = true;
		int specifier = (int)ceil(d - 0.5);

		if (0 == specifier) {
			m_engine.SetUseLOSScattering(true);
		}
		else {
			m_engine.SetUseLOSScattering(false);
		}
		return ok;
	}
	);

	AddSetScalarFunction("prefillsolartransmission",
		[&, this](double d)
		{
			bool ok = true;
			int specifier = (int)ceil(d-0.5);

			if (0 == specifier)
			{
				m_engine.SetPrefillSolarTransmission(false);
			}
			else
			{
				m_engine.SetPrefillSolarTransmission(true);
			}
			return ok;
		}
	);

	AddSetScalarFunction("usesolartableforsinglescattering",
		[&, this](double d)
		{
			bool ok = true;
			int specifier = (int)ceil(d-0.5);

			if (0 == specifier)
			{
				m_engine.SetUseSolarTableForSingleScattering(false);
			}
			else
			{
				m_engine.SetUseSolarTableForSingleScattering(true);
			}
			return ok;
		}
	);

	AddSetScalarFunction("opticalphasetableresolution",
		[&, this](double d)
	{
		bool ok = true;
		m_specs.OpticalPropertiesSpecs().SetScatterResolution(d);
		return ok;
	}
	);

	AddSetScalarFunction("numsolartablesza",
		[&, this](double d)
		{
			bool ok = true;
			int numsza = (int)ceil(d-0.5);
			m_engine.SetNumSolarTableSZA(numsza);
			return ok;
		}
	);

	AddSetScalarFunction("nadirreferencepointonground",
		[&, this](double d)
		{
			bool ok = true;
			int onground = (int)ceil(d - 0.5);
			if (onground == 0)
			{
				m_specs.RayTracingSpecs().SetNadirReferencePointOnGround(false);
			}
			else
			{
				m_specs.RayTracingSpecs().SetNadirReferencePointOnGround(true);
			}
			return ok;
		}
	);


	return true;
}



/*---------------------------------------------------------------------------
 *           ISKEngine_Stub_HR::MakeVectorSetFunctions            2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::MakeVectorSetFunctions()
{
	AddSetVectorFunction( "setsun",
		[&, this](const double* sun, int n) 
		{
			bool ok;

			ok =  CheckModelNotInitalized("setsun");
			if (ok)
			{
				ok = ( n == 3 );
				ok = ok && m_engine.SetSun( nxVector( sun[0], sun[1], sun[2] ) );
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING,"ISKEngine HR, Error setting property setsun");
				}
			}
			return ok;
		}
	);

	AddSetVectorFunction( "diffuseincomingresolution",
		[&, this](const double* reso, int n)
		{
			bool ok;

			ok =  CheckModelNotInitalized("diffuseincomingresolution");
			if (ok)
			{
				ok = ( n == 4 );
				if (ok)
				{
					size_t	numazi = (size_t) ceil(reso[3]-0.5);

					m_specs.DiffuseSpecs().SetNumBeforeHoriz( (size_t) ceil(reso[0]-0.5) );
					m_specs.DiffuseSpecs().SetNumHoriz( (size_t) ceil(reso[1]-0.5) );
					m_specs.DiffuseSpecs().SetNumAfterHoriz( (size_t) ceil(reso[2]-0.5 ));
					if (numazi > 0) m_specs.DiffuseSpecs().SetNumAzi( numazi );
				}
				else
				{
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Error setting property diffuseincomingresolution");
				}
			}
			return ok;
		}
	);
			
	AddSetVectorFunction( "setreferencepoint",
		[&, this](const double* value, int n)
		{
			bool ok;

			ok =  CheckModelNotInitalized("setreferencepoint");
			ok = ok && (n == 4);
			if (ok) 
			{
				ok = m_engine.InternalSpecsVar()->RayManager().SetReferencePoint( value[0], value[1], value[2], value[3] );
			}
			return ok;
		}
	);


	AddSetVectorFunction( "manualdiffuseheights", 
		[&, this](const double* h, int n)
		{
			bool ok;

			ok =  CheckModelNotInitalized("manualdiffuseheights");
			if (ok)
			{
				m_specs.DiffuseSpecs().SetManualDiffuseHeights( std::vector<double>(h, h+n) );
			}
			return ok;
		}
	);

	AddSetVectorFunction( "precachewavel", 
		[&, this](const double* h, int n)
		{
			bool ok;

			ok =  CheckModelNotInitalized("precachewavel");
			if (ok)
			{
				m_specs.OpticalPropertiesSpecs().SetPrecacheWavel( std::vector<double>(h, h+n) );
			}
			return ok;
		}
	);

	AddSetVectorFunction( "manualopticalheights",
		[&, this](const double* h, int n)
		{
			bool ok;

			ok =  CheckModelNotInitalized("manualopticalheights");
			if (ok)
			{
				m_specs.OpticalPropertiesSpecs().SetHeightGrid( std::vector<double>(h,h+n) );
			}
			return ok;
		}
	);

	AddSetVectorFunction( "opticalnormalandreference",
		[&, this](const double* h, int n)
		{
			bool ok;

			ok =  CheckModelNotInitalized("opticalnormalandreference");
			if (ok)
			{
				ok = (n == 6);
				if (ok)
				{
					nxVector normal( h[0], h[1], h[2] );
					nxVector reference( h[3], h[4], h[5] );
					m_specs.OpticalPropertiesSpecs().SetNormalAndReference( normal, reference );
				}
				else
				{
					nxLog::Record(NXLOG_WARNING,"ISKEngine HR, The opticalnormalandreference property requires an array of exactly 6 numbers. We received %d parameters", (int)n);
				}
			}
			return ok;
		}
	);

	AddSetVectorFunction( "opticalanglegrid", 
		[&, this]( const double* th, int n )
		{
			bool ok;

			ok =  CheckModelNotInitalized("opticalanglegrid");
			if (ok)
			{
				m_specs.OpticalPropertiesSpecs().SetAngleGrid( std::vector<double>(th,th+n) );
			}
			return ok;
		}
	);

	AddSetVectorFunction( "threedopticaltableparam",
		[&, this]( const double* h, int n )
		{
			bool ok;

			ok =  CheckModelNotInitalized("threedopticaltableparam");
			if (ok)
			{
				ok = (n == 3);
				if (ok)
				{
					m_specs.OpticalPropertiesSpecs().SetDelaunayParam( h[0], (size_t)ceil(h[1]-0.5), (size_t)ceil(h[2]-0.5) );	
				}
				else
				{
					nxLog::Record(NXLOG_WARNING,"ISKEngine HR, The threedopticaltableparam property requires an array of exactly 3 numbers. We received %d parameters", (int)n);
				}
			}
			return ok;
		}
	);

	AddSetVectorFunction( "manualdiffuselocations", 
		[&, this] ( const double* h, int n )
		{
			bool ok;

			ok =  CheckModelNotInitalized("manualdiffuselocations");
			if (ok)
			{
				size_t numprof = n/3;
				std::vector<nxVector> profloc;
				profloc.resize(numprof);

				ok  = (numprof > 0);
				if (ok)
				{

					for( size_t idx = 0; idx < numprof; idx++ )
					{
						profloc[idx].SetCoords( nxVector(h[3*idx], h[3*idx+1], h[3*idx+2]) );
					}
					m_specs.DiffuseSpecs().SetManualDiffuseLocations( profloc );
				}
				else
				{
					nxLog::Record(NXLOG_WARNING,"ISKEngine HR, The manualdiffuselocations property requires at least 1 three element profile locations. We received %d elements or %d locations", (int)n, (int)numprof);
				}
			}
			return ok;
		}
	);

	AddSetVectorFunction( "manualdiffuselatlons", 
		[&, this] ( const double* h, int n )
		{
			bool ok;

			ok =  CheckModelNotInitalized("manualdiffuselatlons");
			if (ok)
			{
				size_t numpairs = n/2;
				ok = numpairs > 0;
				if (ok)
				{
					std::vector<double> proflatlon(h, h + n);
					m_specs.DiffuseSpecs().SetManualDiffuseLatLons( proflatlon );
				}
				else
				{
					nxLog::Record(NXLOG_WARNING,"ISKEngine HR, The manualdiffuselatlon property requires at least 1 lat/lon pair. We received %d elements or %d pairs", (int)n, (int)numpairs);
				}
			}
			return ok;
		}
	);

	AddSetVectorFunction("manualdiffuseszas",
		[&, this](const double* h, int n)
	{
		bool ok;

		ok = CheckModelNotInitalized("manualdiffuseszas");
		if (ok)
		{
			std::vector<double> szas(h, h + n);

			ok = (n > 0);
			if (ok)
			{
				m_specs.DiffuseSpecs().SetManualDiffuseSZAs(szas);
			}
			else
			{
				nxLog::Record(NXLOG_WARNING, "ISKEngine HR, The manualdiffuseszas property requires at least 1 solar zenith angle. We received %d. ", (int)n );
			}
		}
		return ok;
	}
	);


	AddSetVectorFunction("manualdiffuselospositions",
		[&, this](const double* h, int n)
		{
			bool ok;

			ok = CheckModelNotInitalized("manualdiffuselospositions");
			if (ok)
			{
				std::vector<double> positions(h, h + n);

				ok = (n > 0);
				if (ok)
				{
					m_specs.DiffuseSpecs().SetManualDiffuseLOSPositions(positions);
				}
				else
				{
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, The manualdiffuselospositions property requires at least 1 diffuse profile position. We received %d. ", (int)n);
				}
			}
			return ok;
		}
	);

	AddSetVectorFunction("manualdiffuseplanenormalandreference",
		[&, this](const double* h, int n)
		{
			bool ok;

			ok = CheckModelNotInitalized("manualdiffuseplanenormalandreference");
			if (ok)
			{
				ok = (n == 6);
				if (ok)
				{
					nxVector normal(h[0], h[1], h[2]);
					nxVector reference(h[3], h[4], h[5]);
					m_specs.DiffuseSpecs().SetManualDiffusePlaneNormal(normal);
					m_specs.DiffuseSpecs().SetManualDiffusePlaneReference(reference);
				}
				else
				{
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, The manualdiffuseplanenormalandreference property requires an array of exactly 6 numbers. We received %d parameters", (int)n);
				}
			}
			return ok;
		}
	);

	AddSetVectorFunction("manualdiffuseplaneangles",
		[&, this](const double* h, int n)
		{
			bool ok;

			ok = CheckModelNotInitalized("manualdiffuseplaneangles");
			if (ok)
			{
				std::vector<double> angles(h, h + n);

				ok = (n > 0);
				if (ok)
				{
					m_specs.DiffuseSpecs().SetManualDiffusePlaneAngles(angles);
				}
				else
				{
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, The manualdiffuseplaneangles property requires at least 1 diffuse profile position. We received %d. ", (int)n);
				}
			}
			return ok;
		}
	);

	AddSetVectorFunction( "diagnosticscatterorders",
		[&, this](const double * h, int n)
		{
			bool ok;
		
			ok = CheckModelNotInitalized("diagnosticscatterorders");
			if (ok)
			{
				std::vector<size_t> scatords;
				scatords.resize(n);

				if( n > 0 )
				{
					for( size_t idx = 0; idx < n; idx++ )
					{
						scatords[idx] = (size_t)ceil( h[idx] - 0.5 );
					}
					m_specs.DiffuseSpecs().SetDiagnosticScatterOrders( scatords );
					if (!m_issetdiagnostics)
					{
						hid_t fid;
						fid = H5Fcreate( "DiagnosticData.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
						if( fid < 0 )
						{
							nxLog::Record( NXLOG_WARNING, "Could not create h5 diagnostic file. Thats not good" );
							ok = false;
						}
						H5Fclose( fid );
					}
					m_issetdiagnostics = true;
				}
				else
				{
					m_specs.DiffuseSpecs().SetDiagnosticScatterOrders( scatords );
					m_issetdiagnostics = false;
					nxLog::Record( NXLOG_WARNING, "ISKEngine HR, The diagnosticscatterorders property has %i length; Diagnostics will most likely be disabled", n );
				}
			}
			return ok;
		}
	);

	AddSetVectorFunction( "diagnosticdiffuseprofiles", 
		[&, this](const double * h, int n)
		{
			bool ok;

			ok = CheckModelNotInitalized("diagnosticdiffuseprofiles");
			if (ok)
			{
				std::vector<size_t> profs;
				profs.resize(n);

				if( n > 0 )
				{
					for( size_t idx = 0; idx < n; idx++ )
					{
						profs[idx] = (size_t)ceil( h[idx] - 0.5 );
					}
					m_specs.DiffuseSpecs().SetDiagnosticDiffuseProfiles( profs );
					if (!m_issetdiagnostics)
					{
						hid_t fid;
						fid = H5Fcreate( "DiagnosticData.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
						if( fid < 0 )
						{
							nxLog::Record( NXLOG_WARNING, "Could not create h5 diagnostic file. Thats not good" );
							ok = false;
						}
						H5Fclose( fid );
					}
					m_issetdiagnostics = true;
				}
				else
				{
					m_specs.DiffuseSpecs().SetDiagnosticDiffuseProfiles( profs );
					m_issetdiagnostics = false;
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, The diagnosticdiffuseprofiles property has %i length; Diagnostics will most likely be disabled", n );
				}
			}
			return ok;
		}
	);


	AddSetVectorFunction( "wfheights", 
		[&, this](const double* h, int n )
		{
			bool ok;

			ok =  CheckModelNotInitalized("wfheights");
			if (ok)
			{
				m_specs.WeightingFunctionSpecs().SetWeightingFunctionHeight(std::vector<double>(h,h+n) );
			}
			return ok;
		}
	);

	AddSetVectorFunction( "spectralalbedo", 
		[&, this](const double* h, int n)
		{
			int half = n/2;

			skBRDF_VariableAlbedo* alb = new skBRDF_VariableAlbedo;
			alb->SetAlbedo(h + half, h, half);
			m_opticalstate.SetAlbedoObject( alb );
			return true;
		}
	);

	AddSetVectorFunction( "wfwidths",
		[&, this](const double* h, int n )
		{
			bool ok;

			ok =  CheckModelNotInitalized("wfwidths");
			if (ok)
			{
				m_specs.WeightingFunctionSpecs().SetWeightingFunctionWidth(std::vector<double>(h,h+n));
			}
			return ok;
		}
	);

	AddSetVectorFunction("wfwidthsleft",
		[&, this](const double* h, int n)
		{
			bool ok;

			ok = CheckModelNotInitalized("wfwidthsleft");
			if (ok)
			{
				m_specs.WeightingFunctionSpecs().SetWeightingFunctionLeftWidth(std::vector<double>(h, h + n));
			}
			return ok;
		}
	);

	AddSetVectorFunction("wfwidthsright",
		[&, this](const double* h, int n)
			{
			bool ok;

			ok = CheckModelNotInitalized("wfwidthsright");
			if (ok)
			{
				m_specs.WeightingFunctionSpecs().SetWeightingFunctionRightWidth(std::vector<double>(h, h + n));
			}
			return ok;
		}
	);

	AddSetVectorFunction( "manualraytracingshells",
		[&, this](const double* h, int n)
		{
			bool ok;

			ok = CheckModelNotInitalized("manualraytracingshells");
			if (ok)
			{
				ok = ok && m_specs.RayTracingSpecs().SetManualShells(std::vector<double>(h, h + n));
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Error setting the manualraytracingshells property");
				}
			}
			return ok;
		}
	);

	AddSetVectorFunction( "manualsolarraytracingshells",
		[&, this](const double* h, int n)
		{
			bool ok;

			ok = CheckModelNotInitalized("manualsolarraytracingshells");
			if (ok)
			{
				ok = ok && m_specs.RayTracingSpecs().SetManualSolarShells(std::vector<double>(h, h + n));
				if (!ok)
				{
					nxLog::Record(NXLOG_WARNING, "ISKEngine HR, Error setting the manualsolarraytracingshells property");
				}
			}
			return ok;
		}
	);
	return true;

}


/*---------------------------------------------------------------------------
 *           ISKEngine_Stub_HR::MakeScalarGetFunctions            2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::MakeScalarGetFunctions()
{
	AddGetScalarFunction( "numlinesofsight",
		[&, this]( double* val )
		{
			*val = (double)m_linesofsight.NumRays();
			return true;
		}
	);
	return true;
}

bool ISKEngine_Stub_HR::MakeVectorGetFunctions()
{
	AddGetVectorFunction( "referencepoint",
		[&, this]( int index )
		{
			GEODETIC_INSTANT pt = m_engine.ReferencePoint();
			m_getpropertybuffer.resize(4);
			m_getpropertybuffer[0] = pt.latitude;
			m_getpropertybuffer[1] = pt.longitude;
			m_getpropertybuffer[2] = pt.heightm;
			m_getpropertybuffer[3] = pt.mjd;
			return true;
		}
	);

	AddGetVectorFunction( "sun",
		[&, this]( int index )
		{
			nxVector sun = m_engine.GetSun();
			m_getpropertybuffer.resize(3);
			m_getpropertybuffer[0] = sun.X();
			m_getpropertybuffer[1] = sun.Y();
			m_getpropertybuffer[2] = sun.Z();
			return true;
		}
	);

	AddGetVectorFunction( "wavel",
		[&, this]( int index )
		{
			m_getpropertybuffer = m_wavelen;
			return true;
		}
	);

	AddGetVectorFunction( "brdfwf",
		[&, this](int index)
		{
			size_t size = m_brdfwfbuffer.size();
			m_getpropertybuffer.resize(size);
			for (int idx = 0; idx < size; idx++)
			{
				m_getpropertybuffer[idx] = m_brdfwfbuffer.ArrayBasePtr()[idx];
			}

			return true;
		}
	);

	AddGetVectorFunction( "observer",
		[&, this]( int index )
		{
			bool ok = true;

			ok = ok && (int)m_linesofsight.NumRays() > index;
			ok = ok && index >= 0;
			if( !ok )
			{
				nxLog::Record(NXLOG_WARNING, "ISKENGINE_Stub_HR::GetPropertyArray, The observer array getproperty specifies an index (%d) which is out of the line of sight range [0..%d]", (int)index, (int)m_linesofsight.NumRays());
			}
			/* else */
			const SKTRAN_LineOfSightEntry_V2* entry;
			ok = ok && m_linesofsight.GetRay( index, &entry );
			nxVector pt = entry->Observer();
			m_getpropertybuffer.resize(3);
			m_getpropertybuffer[0] = pt.X();
			m_getpropertybuffer[1] = pt.Y();
			m_getpropertybuffer[2] = pt.Z();
			return ok;
		}
	);
	AddGetVectorFunction( "look", 
		[&, this]( int index )
		{
			bool ok = true;

			ok = ok && (int)m_linesofsight.NumRays() > index;
			ok = ok && index >= 0;
			if( !ok )
			{
				nxLog::Record(NXLOG_WARNING, "ISKENGINE_Stub_HR::GetPropertyArray, The observer array getproperty specifies an index (%d) which is out of the line of sight range [0..%d]", (int)index, (int)m_linesofsight.NumRays());
			}
			/* else */
			const SKTRAN_LineOfSightEntry_V2* entry;
			ok = ok && m_linesofsight.GetRay( index, &entry );
			nxVector pt = entry->Look();
			m_getpropertybuffer.resize(3);
			m_getpropertybuffer[0] = pt.X();
			m_getpropertybuffer[1] = pt.Y();
			m_getpropertybuffer[2] = pt.Z();
			return ok;
		}
	);

	AddGetVectorFunction( "stokesvec",
		[&, this]( int index )
		{
			bool ok = true;

			ok = ok && 0<=index && index < (int)(m_radiancePol.XSize()*m_radiancePol.YSize());
			if(ok){
				int numrays = (int) m_radiance.XSize();
				int numwav  = (int) m_radiance.YSize();
				int xindex = index % numrays;
				int yindex = index / numrays;
				m_getpropertybuffer.resize(4);

				if( m_storeStokes ){ 
					m_getpropertybuffer[0] = m_radiancePol.At(xindex, yindex).At(1);
					m_getpropertybuffer[1] = m_radiancePol.At(xindex, yindex).At(2);
					m_getpropertybuffer[2] = m_radiancePol.At(xindex, yindex).At(3);
					m_getpropertybuffer[3] = m_radiancePol.At(xindex, yindex).At(4);
				} else{
					m_getpropertybuffer[0] = m_radiance.At(xindex, yindex);
					m_getpropertybuffer[1] = 0.0;
					m_getpropertybuffer[2] = 0.0;
					m_getpropertybuffer[3] = 0.0;
				}
			} else{
				nxLog::Record(NXLOG_WARNING, "ISKENGINE_Stub_HR::GetPropertyArray, No stokes vector stored for index %i.", (int)index );
			}
			return ok;
		}
	);

	AddGetVectorFunction( "basis",
		[&, this]( int index )
		{
			bool ok = true;

			ok = ok && 0 <= index && index < (int)m_linesofsight.NumRays();
			if( ok )
			{
				GEOGRAPHIC_BASIS	basis;
				ok = GetBasisGeo( &basis, index);
				m_getpropertybuffer.resize(9);
				m_getpropertybuffer[0] = basis.X().X();
				m_getpropertybuffer[1] = basis.X().Y();
				m_getpropertybuffer[2] = basis.X().Z();
				m_getpropertybuffer[3] = basis.Y().X();
				m_getpropertybuffer[4] = basis.Y().Y();
				m_getpropertybuffer[5] = basis.Y().Z();
				m_getpropertybuffer[6] = basis.Z().X();
				m_getpropertybuffer[7] = basis.Z().Y();
				m_getpropertybuffer[8] = basis.Z().Z();
			}
			else
			{
				nxLog::Record(NXLOG_WARNING, "ISKENGINE_Stub_HR::GetPropertyArray, The observer array getproperty specifies an index (%d) which is out of the line of sight range [0..%d]", (int)index, (int)m_linesofsight.NumRays());
			}
			return ok;
		}
	);

	AddGetVectorFunction( "tangentalts",
		[&]( int index )
		{
			bool ok = true;
			m_getpropertybuffer.resize(m_linesofsight.NumRays());
			double ground = m_engine.Coordinates()->AltitudeToRadius(0.0);
			for( int rayidx = 0; rayidx < m_linesofsight.NumRays(); rayidx++ )
			{
				auto ray = m_engine.LinesOfSight().RayAt(rayidx);
				double Robs, Tobs, Rt;
				ray->CalculateBaseLineTangentPointDetails(0, &Robs, &Tobs, &Rt);

				m_getpropertybuffer[rayidx] = Rt - ground;
			}

			return ok;
		}
	);

	AddGetVectorFunction( "loscellstartopticaldepth",
		[&]( int index )
		{
			bool ok = true;
			auto& opticaldepths = m_cellOpticalDepthBuffer[index];
			size_t numcells = opticaldepths.size();

			m_getpropertybuffer.resize(numcells);

			for( int cellidx = 0; cellidx < numcells; cellidx++ )
			{
				m_getpropertybuffer[cellidx] = opticaldepths[cellidx];
			}

			return ok;
		}
	
	);

	AddGetVectorFunction( "loscellstartdistance",
		[&]( int index )
		{
			bool ok = true;
			auto& distances = m_cellDistanceBuffer[index];
			size_t numcells = distances.size();

			m_getpropertybuffer.resize(numcells);

			for( int cellidx = 0; cellidx < numcells; cellidx++ )
			{
				m_getpropertybuffer[cellidx] = distances[cellidx];
			}

			return ok;
		}
	
	);

	return true;
}


/*---------------------------------------------------------------------------
 *           ISKEngine_Stub_HR::MakeObjectSetFunctions            2019-10-01 */
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::MakeObjectSetFunctions()
{
	AddSetObjectFunction( "albedoclim",
		[&, this]( nxUnknown* obj )
		{
			skClimatology* clim = dynamic_cast<skClimatology*>(obj);
			if( clim == nullptr )
			{
				return false;
			}
			else
			{
				skBRDF_AlbedoPlane* alb = new skBRDF_AlbedoPlane(clim);
				m_opticalstate.SetAlbedoObject( alb );
				return true;
			}
		}
	);

	AddSetObjectFunction( "solarspectrum",
		[&, this]( nxUnknown* obj )
		{
			skSolarSpectrum* solarspectrum = dynamic_cast<skSolarSpectrum*>(obj);
			if ( solarspectrum == nullptr )
			{
				return false;
			}
			else
			{
				m_opticalstate.EmissionObjectVar()->SetSolarSpectrum( solarspectrum );
				return true;
			}
		}
	);

	return true;
}

bool ISKEngine_Stub_HR::MakeStringSetFunctions()
{
	AddSetStringFunction("wfspecies",
		[&, this](const char* cstr)
	{
		std::vector<std::string> wf_handles;
		std::string s(cstr);
		boost::regex e("(\\w+)");
		boost::smatch m;

		while (boost::regex_search(s, m, e)) {
			std::string handle_str = m[0];
			wf_handles.push_back(handle_str);
			s = m.suffix().str();
		}

		m_specs.WeightingFunctionSpecs().SetWeightingFunctionSpeciesString(wf_handles);

		return true;
	}
	);
	return true;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::SetInternalPolarizationMode		 2015- 10- 30*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::SetPolarizationMode( int specifier)
{
	bool	ok = true;

	if( 0==specifier )
	{
		m_specs.OpticalPropertiesSpecs().SetMaxPolarizationOrder( 0 );
		m_storeStokes = false;
	} else if ( 0<specifier ){
		m_specs.OpticalPropertiesSpecs().SetMaxPolarizationOrder( specifier );
		m_storeStokes = true;
	} else{
		nxLog::Record(NXLOG_WARNING, "SKEngine HR, Invalid Specifier [%d] for PolarizationType", specifier );
		ok = false;
	}

	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::GetBasis		 2015- 10- 29*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::GetBasisHelio( GEOGRAPHIC_BASIS* basis, size_t losindex)
{
	bool								ok = true;
	const SKTRAN_LineOfSightEntry_V2*	entry;
			
	ok = ok && nullptr!=basis;
	ok = ok && m_engine.LinesOfSight().LinesOfSightArray()->GetRay( losindex, &entry );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_HR::GetBasis, error fetching ray (%d). Cannot fetch BASIS for this ray", (int)losindex);
	}
	else
	{
		const HELIODETIC_POINT observerh;

		HELIODETIC_VECTOR observervector = m_engine.Coordinates()->GeographicToHelio(entry->Observer());
		HELIODETIC_POINT observerpoint;

		observerpoint.FromVector(observervector, m_engine.Coordinates());

		HELIODETIC_UNITVECTOR lookh;
		lookh.FromVector(m_engine.Coordinates()->GeographicToHelio(entry->Look()));

		HELIODETIC_BASIS heliobasis(observerpoint, lookh);

		const nxVector propagation = m_engine.Coordinates()->HelioUnitVectorToGeographic(heliobasis.X());
		const nxVector solar_theta = m_engine.Coordinates()->HelioUnitVectorToGeographic(heliobasis.Y());
		const nxVector solar_phi = m_engine.Coordinates()->HelioUnitVectorToGeographic(heliobasis.Z());

		const nxVector geo_phi = entry->Observer().UnitVector().Cross(propagation).UnitVector();
		const nxVector geo_theta = entry->Observer();

		*basis = GEOGRAPHIC_BASIS(propagation, solar_theta, solar_phi);
	}
	return ok;
}

bool ISKEngine_Stub_HR::GetBasisGeo( GEOGRAPHIC_BASIS* basis, size_t losindex)
{
	bool								ok = true;
	const SKTRAN_LineOfSightEntry_V2*	entry;
			
	ok = ok && nullptr!=basis;
	ok = ok && m_engine.LinesOfSight().LinesOfSightArray()->GetRay( losindex, &entry );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_HR::GetBasis, error fetching ray (%d). Cannot fetch BASIS for this ray", (int)losindex);
	}
	else
	{
		const nxVector propagation = -1 * entry->Look();
		const nxVector geo_theta = entry->Observer().UnitVector().Cross(propagation).UnitVector();
		const nxVector geo_phi = propagation.Cross(geo_theta);

		*basis = GEOGRAPHIC_BASIS(propagation, geo_theta, geo_phi);
	}
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::ParseCommandAndIndex		 2015- 10- 29*/
/** **/
/*---------------------------------------------------------------------------*/

//bool ISKEngine_Stub_HR::ParseCommandAndIndex( const nxString& input, nxString& cmd, int& index )
//{
//	nxStringArray	tokens;
//	int				numtoken;
//	bool			ok;
//
//	numtoken = nxStrtok( input, &tokens, "([]) ,:;");
//	ok = (numtoken == 2);
//	if (ok)
//	{
//		cmd = tokens.GetAt(0);
//		index = atoi( tokens.GetAt(1) );
//		cmd.MakeLower();
//	}
//


	//// first check if we have an input with an index
	//std::regex e("^([^\\[\\(]*)[\\[\\(](\\d*)[\\]\\)]");			// i feel bad
	//std::smatch m;
	//std::string str = (const char*)input;

	//bool ok = true;
	//ok = ok && std::regex_search( str, m, e );
	//if( ok )
	//{
	//	// we have an input with an index
	//	ok = (m.size() == 3);
	//	if ( ok )
	//	{
	//		cmd = m[1].str().c_str();
	//		index = atoi(m[2].str().c_str());
	//	}
	//}
//	if( !ok )
//	{
//		cmd = input;
//		cmd.MakeLower();
//		index = -1;
//	}
//	return true;
//}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::AddSpecies		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::AddSpecies( const CLIMATOLOGY_HANDLE& species, ISKClimatology_Stub* climatology, ISKOpticalProperty_Stub* opticalproperty)
{
	skClimatology*				climptr;
	nxUnknown*					optbaseptr;
	skOpticalProperties*		optptr;
	bool						ok;

	optbaseptr = (opticalproperty != nullptr)? opticalproperty->RawObjectPointer() : nullptr;
	climptr   = dynamic_cast<skClimatology*>( climatology->RawObjectPointer() );
	optptr    = (optbaseptr != nullptr) ? dynamic_cast<skOpticalProperties*>( optbaseptr) : nullptr;
	ok        = m_opticalstate.AddSpecies( species, climptr, optptr);

	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::AddEmission		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::AddEmission( const EMISSION_HANDLE& species, ISKEmission_Stub* emissionobject)
{
	skEmission*				emission;
	bool						ok;

	emission = dynamic_cast<skEmission*>( emissionobject->RawObjectPointer() );
	ok      = m_opticalstate.AddEmission( species, emission );
	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::SetAlbedo		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::SetAlbedo( double albedo )
{
	return m_opticalstate.SetAlbedo( albedo );
}

/*-----------------------------------------------------------------------------
*					ISKEngine_Stub_HR::SetBRDF		2017-3-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::SetBRDF(ISKBrdf_Stub* brdf)
{
	skBRDF* brdf_cast = dynamic_cast<skBRDF*>(brdf->RawObjectPointer());

	if (!brdf_cast)
	{
		nxLog::Record(NXLOG_WARNING, "Error in SetBRDF, input BRDF is not a valid BRDF object");
		return false;
	}
	else
	{
		return m_opticalstate.SetAlbedoObject(brdf_cast);
	}
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::SetWavelengths		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::SetWavelengths( const double* wavelen, int numwavelen )
{
	m_wavelen.assign( wavelen, wavelen + numwavelen );
	m_radiance.erase();
	return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::InitializeModel		2014-3-13*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::InitializeModel()
{
	bool ok = true;
	if ( !m_geometryisconfigured )
	{
		ok = ok && m_engine.ConfigureModel( m_specs, m_linesofsight, m_numthreads );
		m_geometryisconfigured = ok;
		if (!ok)
		{
			nxLog::Record(NXLOG_WARNING,"ISKEngine HR, Error initializing the model. Thats not good");
		}
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::CalculateRadiance		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::CalculateRadiance(const double** radiance, int* numwavelens , int* numlinesofsight)
{
	std::vector<double>		r;
	bool					ok = true;
	bool					ok1;
	size_t					numrays = m_linesofsight.NumRays();
	size_t					numwave = m_wavelen.size();

	std::vector< skRTStokesVector >* vecptr = nullptr;
	skRTStokesVector	blank;

	ok = ok && InitializeModel();

	r.resize( numrays );
	m_radiance.SetSize( numrays, numwave );
	if( m_storeStokes )
	{
		vecptr = new std::vector< skRTStokesVector >;
		vecptr->resize( numrays );
		m_radiancePol.SetSize( numrays, numwave );
	}

	if ( m_storeOpticalDepth )
	{
		m_cellDistanceBuffer.resize(numrays * numwave);
		m_cellOpticalDepthBuffer.resize(numrays * numwave);
	}

	m_wfbuffer.erase();
	m_brdfwfbuffer.erase();

	int numwf = (int) m_engine.InternalSpecs().WeightingFunctionSpecs().PertList().StoreSize();
	int numwfspecies = (int)m_engine.InternalSpecs().WeightingFunctionSpecs().WFSpecies().size();
	if (numwf != 0)
	{
		m_wfbuffer.SetSize(numwf * numwfspecies, numrays, numwave);
		m_brdfwfbuffer.SetSize(numrays, numwave);
	}
	for (size_t i = 0; i <numwave; i++)
	{
		ok1 = m_engine.CalculateRadiance( &r, m_wavelen.at(i), m_numordersofscatter, &m_opticalstate, vecptr, m_updateclimatology );
		if (!ok1 && m_storeStokes) 
		{
			blank.At(1) = std::numeric_limits<double>::quiet_NaN(); 
			blank.At(2) = std::numeric_limits<double>::quiet_NaN(); 
			blank.At(3) = std::numeric_limits<double>::quiet_NaN(); 
			blank.At(4) = std::numeric_limits<double>::quiet_NaN(); 
		}
		for (size_t l = 0; l < numrays; l++)
		{
			m_radiance.At(l,i )                       = ok1 ? r.at(l)       : std::numeric_limits<double>::quiet_NaN();
			if( m_storeStokes ) m_radiancePol.At(l,i) = ok1 ? vecptr->at(l) : blank;
			if ((numwf != 0) && ok1)
			{
				for (size_t j = 0; j < numwfspecies; j++)
				{
					auto wf = m_engine.WFForRay(j, l);
					for (size_t wfidx = 0; wfidx < numwf; wfidx++)
					{
						m_wfbuffer.At(wfidx + j*numwf, l, i) = wf.WeightAt(wfidx);
					}
				}
				auto wf = m_engine.WFForRay(0, l);
				m_brdfwfbuffer.At(l, i) = wf.brdfwf();
			}

			if ( m_storeOpticalDepth )
			{
				auto& cellDistance = m_cellDistanceBuffer[i*numrays + l];
				auto& cellOpticalDepth = m_cellOpticalDepthBuffer[i*numrays + l];

				auto ray = m_engine.LinesOfSight().RayAt(l);

				size_t numcells = ray->GetNumCells();

				cellDistance.resize(numcells);
				cellOpticalDepth.resize(numcells);

				for( int cellidx = 0; cellidx < numcells; cellidx++ )
				{
					cellDistance[cellidx] = ray->Storage()->DistanceOfPointFromOrigin(cellidx);
					cellOpticalDepth[cellidx] = ray->OpticalDepthArray()[cellidx];
				}
			}
		}
		m_updateclimatology = false;
		ok = ok && ok1;
	}
	*radiance = m_radiance.ArrayBasePtr();
	*numwavelens = (int)numwave;
	*numlinesofsight = (int)numrays;

	if (vecptr != nullptr) delete vecptr;

	return ok;
}



/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::CalculateStokesVector		 2015- 10- 29*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::CalculateStokesVector(const ISKStokesVector** radiance, int* numwavelens , int* numlinesofsight)
{
	std::vector<double>					r;
	bool								ok = true;
	bool								ok1;
	size_t								numrays = m_linesofsight.NumRays();
	size_t								numwave = m_wavelen.size();
	IQUV								iquv;
	GEOGRAPHIC_BASIS					basis_helio, basis_geo;
	ISKBasisDirection					basdir, basdirgeo;
	std::vector< skRTStokesVector >		vecptr;
	int									numwf, numwfspecies;
	const SKTRAN_CoordinateTransform_V2*	coords = m_engine.Coordinates();

	if( !m_geometryisconfigured )
	{
		m_engine.ConfigureModel(m_specs, m_linesofsight, 0);
		m_geometryisconfigured = true;
	}

	r.resize( numrays );
	m_radiance.SetSize( numrays, numwave );
	vecptr.resize( numrays );
	m_radiancePol.SetSize( numrays, numwave );
	m_radiancepolarized.SetSize( numrays, numwave );

	numwf = (int) m_engine.InternalSpecs().WeightingFunctionSpecs().PertList().StoreSize();
	numwfspecies = (int)m_engine.InternalSpecs().WeightingFunctionSpecs().WFSpecies().size();
	m_wfbuffer.erase();
	m_brdfwfbuffer.erase();
	if (numwf != 0)
	{
		m_wfbuffer.SetSize(numwf * numwfspecies*4, numrays, numwave);
		m_brdfwfbuffer.SetSize(numrays, numwave);
	}

    ISKStokesVector temp_wf;

	for (size_t i = 0; i <numwave; i++)
	{
		ok1 = m_engine.CalculateRadiance( &r, m_wavelen.at(i), m_numordersofscatter, &m_opticalstate, &vecptr, m_updateclimatology );
		for (size_t l = 0; l < numrays; l++)
		{
			m_radiance.At(l,i )   = r.at(l);
			m_radiancePol.At(l,i) = vecptr.at(l);
			iquv.I = vecptr.at(l).At(1);
			iquv.Q = vecptr.at(l).At(2);
			iquv.U = vecptr.at(l).At(3);
			iquv.V = vecptr.at(l).At(4);
			ok1 = ok1 && GetBasisHelio( &basis_helio, l);
			ok1 = ok1 && GetBasisGeo(&basis_geo, l);
			basdir.Assign( basis_helio.X(), basis_helio.Y(), basis_helio.Z() );
			m_radiancepolarized.At(l,i).Assign( iquv, basdir);
			basdirgeo.Assign( basis_geo.X(), basis_geo.Y(), basis_geo.Z() );
			m_radiancepolarized.At(l,i).to_new_basis(basdirgeo);
			if( numwf != 0 )
			{
				for (size_t j = 0; j < numwfspecies; j++)
				{
                    if(m_storeStokes){
                        auto wf = m_engine.PolarizedWFForRay(j, l);

                        for (size_t wfidx = 0; wfidx < numwf; wfidx++)
                        {
                            iquv.I = wf.WeightAt(wfidx).I();
                            iquv.Q = wf.WeightAt(wfidx).Q();
                            iquv.U = wf.WeightAt(wfidx).U();
                            iquv.V = wf.WeightAt(wfidx).V();

                            temp_wf.Assign(iquv, basdir);
                            temp_wf.to_new_basis(basdirgeo);

                            m_wfbuffer.At((wfidx + j*numwf)*4, l, i) = temp_wf.I();
                            m_wfbuffer.At((wfidx + j*numwf)*4 + 1, l, i) = temp_wf.Q();
                            m_wfbuffer.At((wfidx + j*numwf)*4 + 2, l, i) = temp_wf.U();
                            m_wfbuffer.At((wfidx + j*numwf)*4 + 3, l, i) = temp_wf.V();
                        }
                    } else {
                        auto wf = m_engine.WFForRay(j, l);

                        for (size_t wfidx = 0; wfidx < numwf; wfidx++)
                        {
                            m_wfbuffer.At((wfidx + j*numwf)*4, l, i) = wf.WeightAt(wfidx);
                            m_wfbuffer.At((wfidx + j*numwf)*4 + 1, l, i) = 0;
                            m_wfbuffer.At((wfidx + j*numwf)*4 + 2, l, i) = 0;
                            m_wfbuffer.At((wfidx + j*numwf)*4 + 3, l, i) = 0;
                        }
                    }

				}
				//auto wf = m_engine.WFForRay(0, l);
				//m_brdfwfbuffer.At(l, i) = wf.brdfwf();
			}
		}
		m_updateclimatology = false;
		ok = ok && ok1;
	}
	*radiance = m_radiancepolarized.ArrayBasePtr();
	*numwavelens = (int)numwave;
	*numlinesofsight = (int)numrays;
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::Radiance		2014-2-8*/
/** **/
/*---------------------------------------------------------------------------*/
//
//double ISKEngine_Stub_HR::Radiance( int userlosindex, int userwavelenindex)
//{
//	bool	ok;
//	size_t	losindex = userlosindex;
//	size_t	windex   = userwavelenindex;
//	double  value;
//
//	ok = (windex < m_radiance.YSize()) && ( losindex < m_radiance.XSize() );
//	value = ok ? m_radiance.At(losindex,windex) : std::numeric_limits<double>::quiet_NaN();
//	return value;
//};

bool ISKEngine_Stub_HR::GetWeightingFunctions(const double** wf, int* numwavel, int* numlos, int* numwf)
{
	int numcalcwf = (int) m_engine.InternalSpecs().WeightingFunctionSpecs().PertList().StoreSize();
	if(numcalcwf == 0)
	{
		*wf = nullptr;
		*numwavel = 0;
		*numlos = 0;
		*numwf = 0;
	}
	else
	{
		*numwavel = (int) m_wfbuffer.ZSize();
		*numlos = (int) m_wfbuffer.YSize();
		*numwf = (int) m_wfbuffer.XSize();


		*wf = m_wfbuffer.ArrayBasePtr();
	}

	return true;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_SO::SetAtmosphericState		2014-02-10*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::SetAtmosphericState( ISKClimatology_Stub* climatology )
{
	bool ok = true;

	ok = ok && m_opticalstate.SetAtmosphericStateModel( dynamic_cast<skClimatology*>(climatology->RawObjectPointer()) );

	return ok;
}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::MakeDefaultOpticalState		2014-3-11*/
/** **/
/*---------------------------------------------------------------------------*/

bool ISKEngine_Stub_HR::MakeDefaultOpticalState()
{
	skOpticalProperties_RayleighDryAir*		rayleigh;					// Optical properties of one air molecule
	skClimatology_MSIS90*					msis90;	
	skClimatology_LabowOzoneVMR*			o3numberdensity;
	skOpticalProperties_O3_OSIRISRes*		o3_opticalprops;			// optical properties of one O3 molecule
	bool									ok;

	rayleigh         = new skOpticalProperties_RayleighDryAir;
	msis90           = new skClimatology_MSIS90;
	o3numberdensity  = new skClimatology_LabowOzoneVMR;
	o3_opticalprops  = new skOpticalProperties_O3_OSIRISRes;

	m_opticalstate.erase();
	ok =       (rayleigh != NULL) && (msis90 != NULL) && (o3_opticalprops != NULL) && (o3numberdensity != NULL);

	ok = ok && m_opticalstate.AddSpecies( SKCLIMATOLOGY_AIRNUMBERDENSITY_CM3, msis90,          rayleigh );
	ok = ok && m_opticalstate.AddSpecies( SKCLIMATOLOGY_PRESSURE_PA, msis90, rayleigh );
	ok = ok && m_opticalstate.AddSpecies( SKCLIMATOLOGY_TEMPERATURE_K, msis90, rayleigh );
	ok = ok && m_opticalstate.AddSpecies( SKCLIMATOLOGY_O3_CM3, o3numberdensity, o3_opticalprops);
	ok = ok && m_opticalstate.SetAlbedo( 1.0 );

	if (!ok)
	{
		nxLog::Record(NXLOG_WARNING,"MakeSpeciesList, there was an error making the species list");
	}
	return ok;
}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::Set		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/
//
//bool ISKEngine_Stub_HR::SetPropertyScalar( const char* propertyname, double value )
//{
//	nxString str(propertyname);
//	str.MakeLower();
//	auto funciterator = m_scalarsetfunctions.find( str );
//	if( funciterator == std::end(m_scalarsetfunctions) )
//	{
//		// 2015-11-20 ndl303, We dont want Set Scalar reporting an error as matlab checks it when passed 1 element arrays. We need it to fail quietly.
//		//nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_HR::Set, this object does not support any scalar properties including [%s]\n", (const char*)propertyname); 
//		return false;
//	}
//	else
//	{
//		return funciterator->second(value);
//	}
//}
//
/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::Set		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/
//
//bool ISKEngine_Stub_HR::SetPropertyArray( const char* propertyname, const double* value, int numpoints )
//{
//	bool	ok;
//	nxString str(propertyname);
//	str.MakeLower();
//
//	auto funciterator = m_vectorsetfunctions.find( str );
//	if( funciterator == std::end( m_vectorsetfunctions ) )
//	{
//		nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_HR::Set, this object does not support array property [%s]\n", (const char*)propertyname);
//		return false;
//	}
//	else
//	{
//		return funciterator->second(value, numpoints);
//	}
//	return ok;
//}


/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::Set		2014-3-6*/
/** **/
/*---------------------------------------------------------------------------*/
//
//bool ISKEngine_Stub_HR::SetPropertyObject( const char* propertyname, nxUnknown* object )
//{
//	nxString str(propertyname);
//	str.MakeLower();
//	auto funciterator = m_objectsetfunctions.find( str );
//	if( funciterator == std::end(m_objectsetfunctions) )
//	{
//		nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_HR::Set, this object does not support any object properties including [%s]\n", (const char*)propertyname);
//		return false;
//	}
//	else
//	{
//		return funciterator->second(object);
//	};
//}
//

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::GetPropertyScalar		 2014- 5- 16*/
/** **/
/*---------------------------------------------------------------------------*/
//
//bool ISKEngine_Stub_HR::GetPropertyScalar( const char* propertyname, double* value )
//{
//	nxString str(propertyname);
//	str.MakeLower();
//	auto funciterator = m_scalargetfunctions.find( str );
//	if( funciterator == std::end(m_scalargetfunctions) )
//	{
//		return false;
//	}
//	else
//	{
//		return funciterator->second(value);
//	}
//}
//

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::GetPropertyArray		 2014- 5- 16*/
/** **/
/*---------------------------------------------------------------------------*/

//bool ISKEngine_Stub_HR::GetProperty( const char* propertyname, const double** value, int* numpoints )
//{
//	nxString name( propertyname);
//	nxString cmd;
//	int index;
//	double	scalarvalue;
//	bool ok;
//
//	ok = GetPropertyScalar( propertyname, &scalarvalue);
//	if (ok)
//	{
//		m_getpropertybuffer.resize(1, scalarvalue);
//		*numpoints = 0;
//		*value = &m_getpropertybuffer[0];
//	}
//	else
//	{
//
//		ok = ParseCommandAndIndex( name, cmd, index );
//		auto funcinterator = m_vectorgetfunctions.find( cmd );
//		if( funcinterator == std::end(m_vectorgetfunctions ) )
//		{
//			nxLog::Record(NXLOG_WARNING,"ISKEngine(HR), The HR engine not support property <%s>", (const char*) cmd );
//			*numpoints = 0;
//			*value     = NULL;
//			ok = false;
//		}
//		else
//		{
//			ok = funcinterator->second(index);
//			*numpoints = (int)m_getpropertybuffer.size();
//			*value = &m_getpropertybuffer[0];
//		}
//	}
//	return ok;
//}

/*-----------------------------------------------------------------------------
 *					ISKEngine_Stub_HR::GetPropertyObject		 2014- 5- 16*/
/** **/
/*---------------------------------------------------------------------------*/
//
//bool ISKEngine_Stub_HR::GetPropertyObject( const char* propertyname, nxUnknown** object )
//{
//	nxString name( propertyname);
//
//	*object = nullptr;
//	nxLog::Record(NXLOG_WARNING,"ISKEngine_Stub_HR::GetPropertyObject, The HR Object Property does not support parameter %s", (const char*) name );
//	return false;
//}
//
