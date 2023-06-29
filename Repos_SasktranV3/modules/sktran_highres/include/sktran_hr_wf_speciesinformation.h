#pragma once

class SKTRAN_HR_WF_SpeciesInformationBase
{
private:
	std::unique_ptr<SKTRAN_GridDefScatterAngle_V21> m_phaseanglegrid;
	std::unique_ptr<nx2dArray<SKTRAN_ScatMat_MIMSNC>> m_phasefunction;
	bool m_isscatterer;
	std::unique_ptr<std::vector<double>> m_xs;

protected:

	void Initialize(std::unique_ptr<SKTRAN_GridDefScatterAngle_V21>&& phaseanglegrid, 
		            std::unique_ptr<nx2dArray<SKTRAN_ScatMat_MIMSNC>>&& phasefunction,
		            bool isscatterer,
		            std::unique_ptr<std::vector<double>>&& xs);

public:
	SKTRAN_HR_WF_SpeciesInformationBase();
	const SKTRAN_GridDefScatterAngle_V21& PhaseAngleGrid() const { return *m_phaseanglegrid; }
	const nx2dArray<SKTRAN_ScatMat_MIMSNC>& PhaseFunction() const { return *m_phasefunction; }
	const bool IsScatterer() const { return m_isscatterer; }
	const std::vector<double>& CrossSection() const { return *m_xs; }
};

class SKTRAN_HR_WF_SpeciesInformationStandard : public SKTRAN_HR_WF_SpeciesInformationBase
{
private:
	
public:
	SKTRAN_HR_WF_SpeciesInformationStandard(skClimatology& backgroundstate, skOpticalProperties& opticalproperty, const SKTRAN_HR_WF_Store& wfstore,
		const SKTRAN_CoordinateTransform_V2& coordinates, double wavelength, const SKTRAN_HR_Specs_Internal_OpticalPropertiesTable& opticalspecs);
};

class SKTRAN_HR_WF_SpeciesInformationAerosolLogNormal : public SKTRAN_HR_WF_SpeciesInformationBase
{
private:

public:
	SKTRAN_HR_WF_SpeciesInformationAerosolLogNormal(skClimatology& backgroundstate, skOpticalProperties& opticalproperty, const SKTRAN_HR_WF_Store& wfstore,
		const SKTRAN_CoordinateTransform_V2& coords, double wavelength, SKTRAN_AtmosphericOpticalState_V21& opticalstate, skClimatology* aero_clim, const CLIMATOLOGY_HANDLE& aerohandle, double moderadiuschange, double modewidthchange,
		const SKTRAN_HR_Specs_Internal_OpticalPropertiesTable& opticalspecs);
};