/**
 * SASKTRAN TIR Internal Optical Properties Specs
 */

/**
 * SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable
 *
 * Internal specifications class for configuring the TIR optical properties
 * table. Provides methods for reading in user specs and creating a table with
 * the appropriate dimensionality and resolution. See the TIR User Optical
 * Properties Table Specs header file for a detailed description of the
 * available settings. Supports 1D (height) and 2D (height and angle measured
 * along line-of-sight). These settings determine which points atmospheric
 * properties will be explicitly calculated at, including thermal emissions.
 */
class SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable
{
private:
	AtmosphereDimensionTIR						m_tabledim;
	size_t										m_numprofiles;
	double										m_heightspacing;

	nxVector									m_normal;
	nxVector									m_reference;
	std::vector<double>							m_anglegrid;
	std::vector<double>							m_heightgrid;

	bool										m_forcecacheupdates;

	const SKTRAN_TIR_RayTracingRegionManager*	m_raymanager;

private:
	virtual bool Create1dTable(OpticalTablePtrTIR& table,
							   const SKTRAN_CoordinateTransform_V2& coords,
							   double toaHeight);
	virtual bool Create3DUnitSphereTable(OpticalTablePtrTIR& table,
										 const SKTRAN_CoordinateTransform_V2& coords,
										 double toaHeight);
	void MakeDefaultAngleGrid();
	void MakeNormalAndReferenceFromLOS();
	virtual bool MakeHeightGrid(SKTRAN_GridDefOpticalPropertiesRadii_V21& heightgrid,
								double toaHeight);
	virtual bool MakeLOSPlaneSphere(SKTRAN_UnitSphere_V2** unitsphere,
									const SKTRAN_CoordinateTransform_V2& coords);

	bool AddPlaneInformationToRayTracer(SKTRAN_RayTracer_Straight_Generic& raytracer,
										const SKTRAN_CoordinateTransform_V2& coords) const;

public:
	SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable();
	virtual ~SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable();

	virtual bool Configure(const SKTRAN_TIR_Specs_User& specs,
						   const SKTRAN_TIR_RayTracingRegionManager* raymanager);
	virtual bool CreateOpticalTable(OpticalTablePtrTIR& table,
									const SKTRAN_CoordinateTransform_V2& coords,
									double toaHeight);
	virtual bool ConfigureDefaults();
	virtual bool In3dMode() const {
		return m_tabledim == AtmosphereDimensionTIR::dim2;
	}
	double GetHeightSpacing() const { return m_heightspacing; }

	bool AddOpticalInformationToRayTracer(SKTRAN_RayTracer_Straight_Generic& raytracer,
										  const SKTRAN_CoordinateTransform_V2* coords) const;

	friend class SKTRAN_TIR_Specs_Internal_wf;
};

