#include "include/sktran_hr_internals.h"

SKTRAN_HR_Specs_Internal_Integrator::SKTRAN_HR_Specs_Internal_Integrator()
{
	m_integratortype          = SKTRAN_HR_Integrator_Type::SKTRAN_HR_IntegratorType_Straight;
	m_maxadaptiveopticaldepth = 0.0;
	m_maxextinctiongradient   = 0.0;
	m_maxadaptiverayopticaldepth = 0.0;
    m_usesolartransmission    = false;
    m_useemissions            = false;

}

SKTRAN_HR_Specs_Internal_Integrator::~SKTRAN_HR_Specs_Internal_Integrator()
{

}

bool SKTRAN_HR_Specs_Internal_Integrator::ConfigureDefaults()
{
	bool ok = true;

	//m_integratortype = SKTRAN_HR_IntegratorType_Straight;
	//m_integratortype = SKTRAN_HR_IntegratorType_Adaptive;

	return ok;
}

bool SKTRAN_HR_Specs_Internal_Integrator::Configure( const SKTRAN_HR_Specs_User& specs )
{
	bool ok = true;
	m_integratortype          = specs.IntegratorSpecsConst().GetIntegratorType();
	m_maxadaptiveopticaldepth = specs.IntegratorSpecsConst().GetMaxOpticalDepth();
	m_maxadaptiverayopticaldepth = specs.IntegratorSpecsConst().GetMaxRayOpticalDepth();
	m_maxextinctiongradient   = specs.IntegratorSpecsConst().GetMaxExtinctionGradient();
    m_useemissions            = specs.IntegratorSpecsConst().GetUseEmissions();
    m_usesolartransmission    = specs.IntegratorSpecsConst().GetUseSolarTransmission();

	return ok;
}

bool SKTRAN_HR_Specs_Internal_Integrator::CreateIntegrator( const SKTRAN_TableOpticalProperties_Base& optproptable, OptIntegratorPtr& optint, SrcIntegratorPtr& srcint )
{
	bool ok = true;

	if( SKTRAN_HR_IntegratorType_Straight == m_integratortype )
	{
		std::unique_ptr<SKTRAN_OpticalPropertiesIntegrator_Straight> optintegrator_local ( new SKTRAN_OpticalPropertiesIntegrator_Straight );
		std::unique_ptr<SKTRAN_SourceTermIntegrator_Order0>          srcintegrator_local ( new SKTRAN_SourceTermIntegrator_Order0 );

		optintegrator_local->AddRef();
		srcintegrator_local->AddRef();
		ok = ok && optintegrator_local->SetOpticalProps( &optproptable );
		ok = ok && srcintegrator_local->SetOpticalProps( &optproptable );
		optint = std::move( optintegrator_local );
		srcint = std::move( srcintegrator_local );
	}
	else if ( SKTRAN_HR_IntegratorType_Adaptive == m_integratortype )
	{
		
		std::unique_ptr<SKTRAN_OpticalPropertiesIntegrator_Adaptive> optintegrator_local ( new SKTRAN_OpticalPropertiesIntegrator_Adaptive );
		std::unique_ptr<SKTRAN_SourceTermIntegrator_Order2>          srcintegrator_local ( new SKTRAN_SourceTermIntegrator_Order2 );

		optintegrator_local->AddRef();
		srcintegrator_local->AddRef();
		ok = ok && optintegrator_local->SetOpticalProps( &optproptable );
		ok = ok && srcintegrator_local->SetOpticalProps( &optproptable );
		optintegrator_local->SetMaxOpticalDepthOfCell( m_maxadaptiveopticaldepth );
		optintegrator_local->SetMaxExtinctionGradientOfCell( m_maxextinctiongradient );
		optintegrator_local->SetMaxRayOpticalDepthToSplit( m_maxadaptiverayopticaldepth );
		optint = std::move( optintegrator_local );
		srcint = std::move( srcintegrator_local );
	}
	else if ( SKTRAN_HR_IntegratorType_Constant == m_integratortype )
	{
		std::unique_ptr<SKTRAN_OpticalPropertiesIntegrator_ConstantLayers> optintegrator_local(new SKTRAN_OpticalPropertiesIntegrator_ConstantLayers);
		std::unique_ptr<SKTRAN_SourceTermIntegrator_Order2>          srcintegrator_local(new SKTRAN_SourceTermIntegrator_Order2);

		optintegrator_local->AddRef();
		srcintegrator_local->AddRef();
		ok = ok && optintegrator_local->SetOpticalProps(&optproptable);
		ok = ok && srcintegrator_local->SetOpticalProps(&optproptable);
		optint = std::move(optintegrator_local);
		srcint = std::move(srcintegrator_local);
	}
	else
	{
		ok = false;
		nxLog::Record(NXLOG_WARNING, "Error should not be here, SKTRAN_HR_Specs_Internal_Integrator::CreateIntegrator");
	}

	return ok;
}


bool SKTRAN_HR_Specs_Internal_Integrator::GetUseSolarTransmission() const
{
    return m_usesolartransmission;
}


bool SKTRAN_HR_Specs_Internal_Integrator::GetUseEmissions() const
{
    return m_useemissions;
}
