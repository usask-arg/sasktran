/**
 * SASKTRAN TIR Internal Integrator Specifications
 * 2019-04-16
 */

/**
 * SKTRAN_TIR_Specs_Internal_Integrator
 *
 * Internal specifications class for configuring the TIR line of sight
 * integration. Provides methods for reading in the user's integrator specs and
 * creating the integrator object with these settings. See the TIR User
 * Integrator Specs header file for a detailed description of the available
 * settings.
 */
class SKTRAN_TIR_Specs_Internal_Integrator
{
private:
	double								m_maxopticaldepth;
	double								m_minextinctionratio;
	OpticalPropertiesIntegratorTypeTIR	m_opttype;
	SourceTermOrderTIR					m_srctype;
	bool								m_docurvedrays;
	LayerExtinctionTypeTIR				m_extinctiontype;
	SpeciesWFUnitTIR					m_wfunit;

private:
	virtual bool ConfigureDefaults();

public:
	SKTRAN_TIR_Specs_Internal_Integrator();
	virtual ~SKTRAN_TIR_Specs_Internal_Integrator();

	virtual bool Configure(const SKTRAN_TIR_Specs_User& specs);
	virtual bool CreateIntegrator(const SKTRAN_TIR_TableOpticalProperties& optproptable,
								  IntegratorPtrTIR& integrator);
};
