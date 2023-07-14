/**
 * SASKTRAN TIR User Integrator Specifications
 */

/**
 * SKTRAN_TIR_Specs_User_Integrator
 *
 * Allows the user to configure properties of the TIR line-of-sight
 * integration. The available settings are:
 *
 *   IntegratorType - Sets the integration to use a straight or adaptive
 *     method. If the adaptive method is specified, the integrator will split
 *     cells if the optical depth exceeds a given threshold or the ratio
 *     between the extinction at one end of the cell to other end becomes too
 *     small.
 *   MaxOpticalDepth - If the optical depth of a single segment exceeds this
 *     value and the adaptive integration method is enabled, the segment will
 *     be split at its midpoint.
 *   MinExtinctionRatio - The ratio between the extinction at one end of a path
 *     segment (whichever is smaller) to the extinction at the other end of the
 *     segment (whichever is larger). Used if the adaptive integration method
 *     is selected. If the ratio is less than the value of this parameter, the
 *     segment will be split in the middle.
 *   SourceTermType - Sets the source term within a cell to either a 0th or 2nd
 *     order approximation.
 *   LayerExtinctionType - Determines whether the optical depth calculation
 *     will assume a constant extinction value within a cell or compute the
 *     optical depth by assuming that the extinction varies linearly with
 *     height.
 */
class SKTRAN_TIR_Specs_User_Integrator
{
private:
	double								m_maxopticaldepth;
	double								m_minextinctionratio;
	OpticalPropertiesIntegratorTypeTIR	m_opttype;
	SourceTermOrderTIR					m_srctype;
	LayerExtinctionTypeTIR				m_extinctiontype;

private:
	void ConfigureDefaults();

public:
	SKTRAN_TIR_Specs_User_Integrator() { ConfigureDefaults(); }

	void SetMaxOpticalDepth(double opticaldepth) { m_maxopticaldepth = opticaldepth; }
	void SetMinExtinctionRatio(double ratio) { m_minextinctionratio = ratio; }
	void SetOpticalPropertiesType(OpticalPropertiesIntegratorTypeTIR t) { m_opttype = t; }
	void SetSourceTermType(SourceTermOrderTIR t) { m_srctype = t; }
	void SetLayerExtinctionType(LayerExtinctionTypeTIR t) { m_extinctiontype = t; }

	double GetMaxOpticalDepth() const { return m_maxopticaldepth; }
	double GetMinExtinctionRatio() const { return m_minextinctionratio; }
	OpticalPropertiesIntegratorTypeTIR GetOpticalPropertiesType() const { return m_opttype; }
	SourceTermOrderTIR GetSourceTermOrder() const { return m_srctype; }
	LayerExtinctionTypeTIR GetLayerExtinctionType() const { return m_extinctiontype; }
};
