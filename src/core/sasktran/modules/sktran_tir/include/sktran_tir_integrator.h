/**
 * SASKTRAN TIR Integrator
 */

/**
 * SKTRAN_TIR_Integrator
 * 2018-09-13
 * 
 * Base class for performing line of sight integration. Similar to the optical
 * property and source term integration classes found in sktran_common, but due
 * to the simpler integration method in the thermal IR regime, they are
 * combined here. Additionally, the optical properties class specific to the
 * TIR engine is used. This class supports adaptive integration, where the ray
 * is subdivided into smaller segments if the optical depth or the extinction
 * gradient becomes too large. Order 2 source terms are supported for radiance
 * calculations, but currently are not used in the computation of analytic
 * weighting functions due to an unresolved numerical instability.
 */
class SKTRAN_TIR_Integrator : public nxUnknown
{
protected:
	const SKTRAN_TIR_TableOpticalProperties* m_opticalprops;
	double m_maxopticaldepth;
	double m_minextinctionratio;
	OpticalPropertiesIntegratorTypeTIR m_optinttype;
	SourceTermOrderTIR m_srcorder;
	LayerExtinctionTypeTIR m_extinctiontype;
	SpeciesWFUnitTIR m_wfunit;

protected:
	virtual double OpticalDepthOfCell(const SKTRAN_RayOptical_Base* ray,
									  size_t cellidx,
									  size_t wavelidx) const;
	virtual double OpticalDepthOfSegment_LinearWithHeight(size_t cellidx,
														  size_t wavelidx,
														  const SKTRAN_RayOptical_Base* ray) const;
	virtual double OpticalDepthOfSegment_Constant(size_t cellidx,
												  size_t wavelidx,
												  const SKTRAN_RayOptical_Base* ray) const;

public:
	SKTRAN_TIR_Integrator();
	~SKTRAN_TIR_Integrator();

	void ReleaseResources();

	bool SetOpticalProps(const SKTRAN_TIR_TableOpticalProperties* optprop);
	const SKTRAN_TIR_TableOpticalProperties* GetOpticalProps() const { return m_opticalprops; }

	virtual bool IntegrateRay(SKTRAN_RayOptical_Base* baseray,
							  double& radiance,
							  const std::vector<CLIMATOLOGY_HANDLE>& wf_handles,
							  size_t wavelidx) const;

	void SetMaxOpticalDepthOfCell(double opticaldepth) { m_maxopticaldepth = opticaldepth; }
	void SetMinExtinctionRatioOfCell(double extinctionratio) { m_minextinctionratio = extinctionratio; }
	void SetOpticalPropertiesIntegrationType(OpticalPropertiesIntegratorTypeTIR optinttype) { m_optinttype = optinttype; }
	void SetSourceTermOrder(SourceTermOrderTIR srcorder) { m_srcorder = srcorder; }
	void SetLayerExtinctionType(LayerExtinctionTypeTIR extinctiontype) { m_extinctiontype = extinctiontype; }
	void SetSpeciesWFUnit(SpeciesWFUnitTIR unit) { m_wfunit = unit; }

	static bool GetSecondOrderTerms(const double& sourcestart,
									const double& sourcemid,
									const double& sourceend,
									double ds,
									double opticaldepthcell,
									double transmissioncell,
									double rad,
									double& sourcecell,
									double& derirad,
									double& drad_dopt);
	static bool GetQuadraticCoeff(const double& sourcestart,
								  const double& sourcemid,
								  const double& sourceend,
								  double ds,
								  double& a,
								  double& b,
								  double& c);
};
