#include "../sktran_common.h"


SKTRAN_EmissionTable_1D::SKTRAN_EmissionTable_1D ( )
{
}


SKTRAN_EmissionTable_1D::~SKTRAN_EmissionTable_1D()
{
}


bool SKTRAN_EmissionTable_1D::SourceTermAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double* source   ) const
{
    bool ok = true;

    size_t loindex;
    double loweight;
    size_t hiindex;
    double hiweight;

    ok = ok && m_altgrid.FindBoundingIndices( qobj.GetPoint().Altitude(), SKTRAN_GridDefBase_V2::OUTOFBOUND_TRUNCATE, &loindex, &loweight, &hiindex, &hiweight );
    *source = m_emission[loindex]*loweight + m_emission[hiindex]*hiweight;
    
    return ok;
}


bool SKTRAN_EmissionTable_1D::GroundSourceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double* source   ) const
{
    bool ok = true;

    *source = m_groundemission;

    return ok;
}


bool SKTRAN_EmissionTable_1D::MonteCarlo_SingleScatteredRadianceAtPoint( const SKTRAN_SourceTermQueryObject_Base& qobj, double& radiance ) const
{
    bool ok = true;

    //ok = ok = SourceTermAtPoint( qobj, radiance );

    size_t loindex;
    double loweight;
    size_t hiindex;
    double hiweight;

    ok = ok && m_altgrid.FindBoundingIndices( qobj.GetPoint().Altitude(), SKTRAN_GridDefBase_V2::OUTOFBOUND_TRUNCATE, &loindex, &loweight, &hiindex, &hiweight );
    radiance = m_emissionSample[loindex]*loweight + m_emissionSample[hiindex]*hiweight;

    return ok;
}


bool SKTRAN_EmissionTable_1D::MonteCarlo_GroundScatteredRadianceAtPoint(const SKTRAN_SourceTermQueryObject_Base& qobj, double& radiance) const
{
    bool ok = true;

    ok = ok && 0*GroundSourceAtPoint( qobj, &radiance );

    return ok;
}


bool  SKTRAN_EmissionTable_1D::InitializeGeometry( std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords,
                                                   const SKTRAN_GridDefRayTracingShells_V21&             altgrid,
                                                   const SKTRAN_GridDefCosSZA_V21&                       cosszagrid,
                                                   const SKTRAN_GridDefSLON_V21&                         slongrid )
{
    bool ok = true; 

    ok = ok && 1==cosszagrid.NumGridPoints();
    ok = ok && 1==slongrid.NumGridPoints();

    if(!ok){
        nxLog::Record(NXLOG_ERROR, "SKTRAN_EmissionTable_1D::InitializeGeometry, Trying to initialize 1D table with multiple geographic locations." );
    } else{
        GEODETIC_INSTANT point( coords->ReferencePtLatitude(), slongrid.At(0), 0.0, coords->ReferencePointMJD() );
        m_point = point;

        ok = ok && m_altgrid.DeepCopy( altgrid );
        m_emission      .resize( m_altgrid.NumShells() );
        m_emissionSample.resize( m_altgrid.NumShells() );
    }

    if(!ok) nxLog::Record(NXLOG_ERROR, "SKTRAN_EmissionTable_1D::InitializeGeometry, Could not initialize table geometry. " );

    return ok;
}


bool SKTRAN_EmissionTable_1D::ConfigureOptical ( double wavelen, SKTRAN_AtmosphericEmission* emissionObject, const SKTRAN_TableOpticalProperties_Base * opticaltable ) 
{
    bool ok = true; 
	bool useThermalEmissions = false;
    
    nxVector          geoVector;
    GEODETIC_INSTANT  geoinstant;
    HELIODETIC_POINT  heliopoint;

	skEmission*			emissionPtr;
	skEmission_Thermal*	thermalEmissionPtr = nullptr;
    
    
	// check list of emissions for thermal emission
	useThermalEmissions = emissionObject->GetSpeciesEmissionObject(SKEMISSION_THERMAL, &emissionPtr);
	if (useThermalEmissions)
	{
		// dynamic cast to ensure that the call to SetAbsorptionPerM only occurs for emission objects of the class skEmission_Thermal
		thermalEmissionPtr = dynamic_cast<skEmission_Thermal*>(emissionPtr);
	}
	
	ok = ok && emissionObject->SetWavelength( wavelen );

    geoinstant.mjd       = opticaltable->Coordinates()->ReferencePointMJD();
    geoinstant.latitude  = opticaltable->Coordinates()->ReferencePtLatitude();
    geoinstant.longitude = opticaltable->Coordinates()->ReferencePtLongitude();

    double kscatpermeter;
	double kabspermeter;
    
    geoinstant.heightm = 0.0;
    ok = ok && emissionObject->SetTimeAndLocation( geoinstant, true, false );
    m_groundemission = emissionObject->IsotropicRadiance();
    heliopoint = opticaltable->Coordinates()->ReferencePoint( 0.0 );
    opticaltable->ScatteringExtinctionPerCM( heliopoint );
    for ( int altidx=0; altidx<m_altgrid.NumShells(); ++altidx )
    {
        geoinstant.heightm = m_altgrid.At( altidx );
        heliopoint = opticaltable->Coordinates()->ReferencePoint( m_altgrid.At(altidx) );
        ok = ok && emissionObject->SetTimeAndLocation( geoinstant, false, false );
		if (thermalEmissionPtr)
		{
			// send absorption per meter at current location to the thermal emission object
			kabspermeter = (opticaltable->TotalExtinctionPerCM(heliopoint) - opticaltable->ScatteringExtinctionPerCM(heliopoint)) * 100.0;
			thermalEmissionPtr->SetAbsorptionPerM(kabspermeter);
		}
        m_emission[altidx]  = emissionObject->IsotropicRadiance();
        kscatpermeter = opticaltable->ScatteringExtinctionPerCM( heliopoint ) * 100.0;
        m_emissionSample[altidx] = m_emission[altidx] * (kscatpermeter>0 ? 1.0/kscatpermeter : 0.0); // kext/albedo = 1/kscat -- this factor sometimes gets multiplied onto the emission and has to be divided out now to make everything work properly without repetetively looking up kscatt 
        //m_emission[altidx] += PlanckBlackbody( temperatureK );
    }

    return ok;
}


