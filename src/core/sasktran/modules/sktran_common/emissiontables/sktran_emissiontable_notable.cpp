#include "../sktran_common.h"

SKTRAN_EmissionTable_NoTable::SKTRAN_EmissionTable_NoTable ( )
{
    m_filename = "C:/ARGsoftware/a-band-spectra/tabulatedEmission_nxFormat.txt";
 //   m_neutralatmosphere = new skClimatology_MSIS90;
 //   m_neutralatmosphere->AddRef( );
    m_solarSpectrum     = new skSolarSpectrum_SAO2010;
    m_solarSpectrum->AddRef( );
    m_emission          = new skEmission_Tabulated_HeightWavelength;
    m_emission->AddRef( );
    m_coords = nullptr;

    m_mjd        = std::numeric_limits<double>::quiet_NaN();;
    m_wavenumber = std::numeric_limits<double>::quiet_NaN();
}


SKTRAN_EmissionTable_NoTable::~SKTRAN_EmissionTable_NoTable()
{
//    m_neutralatmosphere->Release( );
    m_solarSpectrum->Release( );
    m_coords->Release( );
}


void SKTRAN_EmissionTable_NoTable::SetTableFile( std::string& filename )
{
    m_filename = filename;
}


bool SKTRAN_EmissionTable_NoTable::SourceTermAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double* source   ) const
{
    bool ok = true;
    GEODETIC_INSTANT point ( m_coords->PointToGeodetic( qobj.GetPoint(), m_mjd ) );
#pragma omp critical
{
    ok = ok && m_table.SetTimeAndLocation( point, false, false );
    *source = m_table.IsotropicRadiance();
}
    return ok;
}


bool SKTRAN_EmissionTable_NoTable::GroundSourceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double* source   ) const
{
    bool ok = true;
    GEODETIC_INSTANT point ( m_coords->PointToGeodetic( qobj.GetPoint(), m_mjd ) );
#pragma omp critical
{
    ok = ok && m_table.SetTimeAndLocation( point, true, false );
    *source = m_table.IsotropicRadiance();
}
    return ok;
}


bool SKTRAN_EmissionTable_NoTable::MonteCarlo_SingleScatteredRadianceAtPoint ( const SKTRAN_SourceTermQueryObject_Base& qobj, double& radiance ) const
{
    bool ok = false;

    return ok;
}


bool SKTRAN_EmissionTable_NoTable::MonteCarlo_GroundScatteredRadianceAtPoint(const SKTRAN_SourceTermQueryObject_Base& qobj, double& radiance) const
{
    bool ok = false;

    return ok;
}


bool SKTRAN_EmissionTable_NoTable::Initialize ( )
{
    bool ok = true;
    ok = ok && m_emission->LoadHeightWavelengthProfileFromFile( m_filename.c_str() );
    if(!ok) nxLog::Record( NXLOG_WARNING, "SKTRAN_EmissionTable_1D::Initialize, Could not load table from file %s", m_filename.c_str() );
//    ok = ok && m_table.SetAtmosphericStateModel ( m_neutralatmosphere );
    ok = ok && m_table.SetSolarSpectrum         ( m_solarSpectrum );
    ok = ok && m_table.AddEmission( SKEMISSION_PHOTOCHEMICAL_0, m_emission );

    return ok;
}


void SKTRAN_EmissionTable_NoTable::InitializeGeometry ( const SKTRAN_CoordinateTransform_V2* coords ) 
{
    coords->AddRef( );
    if( nullptr!=m_coords ) m_coords->Release( );
    m_coords = coords; 
    m_mjd = m_coords->ReferencePointMJD(); 
}


bool SKTRAN_EmissionTable_NoTable::ConfigureOptical(double wavelen, SKTRAN_AtmosphericEmission*, const SKTRAN_TableOpticalProperties_Base* ) 
{
    bool ok = true;
    
    ok = ok && m_table.SetWavelength ( wavelen );

    return ok;
}


