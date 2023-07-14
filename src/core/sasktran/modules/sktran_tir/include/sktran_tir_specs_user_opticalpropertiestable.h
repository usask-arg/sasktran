/**
 * SASKTRAN TIR User Optical Properties Specifications
 */

/**
 * SKTRAN_TIR_Specs_User_OpticalPropertiesTable
 *
 * Allows the user to configure properties of the TIR optical properties table.
 * The available settings are:
 *
 *   OpticalPropertiesTableDimension - 1D and 2D tables are supported. 1D
 *     tables vary in height while 2D tables vary in height and an angle along
 *     the line of sight.
 *   NumProfiles - If a 2D table is specified, this sets the number of
 *     horizontal location points (i.e. a height profile is created at each
 *     point)
 *   HeightRes - Sets the spacing between height points in meters
 *   NormalAndReference - Allows manually defining the optical plane of a 2D
 *     table. Normal is a unit vector normal to the optical table plane, while
 *     reference is the x axis of the plane.
 *   AngleGrid - A vector containing angles about the reference point where the
 *     2D table is defined, in degrees. An angle of 0 deg corresponds to the
 *     reference point.
 *   HeightGrid - A vector of heights, in meters, to define the table height.
 *     Generally this setting is not used and instead HeightRes is set which
 *     creates an evenly spaced table.
 *   DelaunayParam
 *   GroundEmissivity - Constant emissivity of the Earth's surface. Thermal
 *     emissions from the Earth's surface are modelled by blackbody emission
 *     where the temperature is determined from the bottom altitude of the
 *     atmospheric table.
 *   ForceCacheUpdates - If true, forces a cache update everything the time and
 *     place of the optical state object is changed. Can be left at its default
 *     value of false because the engine will perform cache updates when they
 *     are needed.
 */
class SKTRAN_TIR_Specs_User_OpticalPropertiesTable
{
private:
	AtmosphereDimensionTIR				m_tabledim;
	size_t								m_numprofiles;
	double								m_heightres;

	nxVector							m_normal;
	nxVector							m_reference;
	std::vector<double>					m_anglegrid;
	std::vector<double>					m_heightgrid;

	bool								m_forcecacheupdates;

	double								m_groundemissivity;

private:
	bool ConfigureDefaults();

public:
	SKTRAN_TIR_Specs_User_OpticalPropertiesTable() { ConfigureDefaults(); }
	~SKTRAN_TIR_Specs_User_OpticalPropertiesTable() {};
	void SetOpticalPropertiesTableDimension(AtmosphereDimensionTIR type) { m_tabledim = type; }
	void SetNumProfiles(size_t numprofiles) { m_numprofiles = numprofiles; }
	void SetHeightRes(double h) { m_heightres = h; }
	void SetNormalAndReference(const nxVector& normal, const nxVector& reference) { m_normal = normal; m_reference = reference; }
	void SetAngleGrid(const std::vector<double>& anglegrid) { m_anglegrid = anglegrid; }
	void SetHeightGrid(const std::vector<double>& heightgrid) { m_heightgrid = heightgrid; }
	void SetGroundEmissivity(const double emissivity) { m_groundemissivity = emissivity; }
	void SetForceCacheUpdates(bool update) { m_forcecacheupdates = update; }

	const double GetGroundEmissivity() const { return m_groundemissivity; }

	friend class SKTRAN_TIR_Specs_Internal_OpticalPropertiesTable;
};
